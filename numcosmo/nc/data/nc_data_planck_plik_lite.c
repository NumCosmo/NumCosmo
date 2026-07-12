/***************************************************************************
 *            nc_data_planck_plik_lite.c
 *
 *  Fri June 27 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_planck_plik_lite.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 *
 * numcosmo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * numcosmo is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcDataPlanckPlikLite:
 *
 * Planck plik_lite (CMB-marginalized, foreground-free) high-ℓ likelihood.
 *
 * Native reimplementation of the Planck `plik_cmbonly` likelihood as an
 * #NcmDataGaussCov subclass. The data are binned cross-spectra bandpowers
 * (raw $C_\ell$, in $\mu K^2$) with a fixed covariance; the model bandpower for
 * each bin is the normalized weighted average of the theory $C_\ell$ over the
 * bin, divided by the absolute calibration $A_\mathrm{planck}^2$. Being a
 * #NcmDataGaussCov, it supports resampling (unlike the clik wrapper).
 *
 * Bin operator (per active bin $b$, absolute multipole range
 * $[\ell^{\min}_b, \ell^{\max}_b]$ with weights summing to one):
 * $$ X^{\rm model}_b = \frac{1}{A_\mathrm{planck}^2}
 *    \sum_{\ell=\ell^{\min}_b}^{\ell^{\max}_b} w_{b,\ell}\, C_\ell^{s(b)}, $$
 * where $s(b) \in \{TT, EE, TE\}$ is the bin's spectrum.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/data/nc_data_planck_plik_lite.h"
#include "nc/cmb/nc_planck_fi.h"
#include "nc/cmb/nc_planck_fi_cor_tt.h"
#include "ncm/core/ncm_serialize.h"

/* Per-bin spectrum tag (index into the cached theory Cl vectors). */
typedef enum
{
  NC_DATA_PLIK_LITE_SPEC_TT = 0,
  NC_DATA_PLIK_LITE_SPEC_EE,
  NC_DATA_PLIK_LITE_SPEC_TE,
  NC_DATA_PLIK_LITE_SPEC_LEN,
} NcDataPlikLiteSpec;

struct _NcDataPlanckPlikLite
{
  /*< private >*/
  NcmDataGaussCov parent_instance;
  NcHIPertBoltzmann *pb;
  NcmVector *bin_lmin;    /* absolute ell start, per active bin */
  NcmVector *bin_lmax;    /* absolute ell end,   per active bin */
  NcmVector *bin_weight;  /* flattened weights, concatenated per bin (sum = 1 per bin) */
  NcmVector *spectrum_id; /* NcDataPlikLiteSpec per active bin */
  guint lmax;             /* max multipole needed from the Boltzmann code */
  gchar *calib_name;      /* nuisance parameter name for the absolute calibration */
  /* runtime caches (not serialized) */
  NcmVector *cl_cache[NC_DATA_PLIK_LITE_SPEC_LEN];
  guint calib_pi;
  gboolean calib_pi_set;
};

enum
{
  PROP_0,
  PROP_PERT_BOLTZMANN,
  PROP_BIN_LMIN,
  PROP_BIN_LMAX,
  PROP_BIN_WEIGHT,
  PROP_SPECTRUM_ID,
  PROP_LMAX,
  PROP_CALIB_NAME,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataPlanckPlikLite, nc_data_planck_plik_lite, NCM_TYPE_DATA_GAUSS_COV)

static void
nc_data_planck_plik_lite_init (NcDataPlanckPlikLite *plik)
{
  guint i;

  plik->pb           = NULL;
  plik->bin_lmin     = NULL;
  plik->bin_lmax     = NULL;
  plik->bin_weight   = NULL;
  plik->spectrum_id  = NULL;
  plik->lmax         = 0;
  plik->calib_name   = NULL;
  plik->calib_pi     = 0;
  plik->calib_pi_set = FALSE;

  for (i = 0; i < NC_DATA_PLIK_LITE_SPEC_LEN; i++)
    plik->cl_cache[i] = NULL;
}

static void
nc_data_planck_plik_lite_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataPlanckPlikLite *plik = NC_DATA_PLANCK_PLIK_LITE (object);

  g_return_if_fail (NC_IS_DATA_PLANCK_PLIK_LITE (object));

  switch (prop_id)
  {
    case PROP_PERT_BOLTZMANN:
      nc_data_planck_plik_lite_set_hipert_boltzmann (plik, g_value_get_object (value));
      break;
    case PROP_BIN_LMIN:
      ncm_vector_substitute (&plik->bin_lmin, g_value_get_object (value), TRUE);
      break;
    case PROP_BIN_LMAX:
      ncm_vector_substitute (&plik->bin_lmax, g_value_get_object (value), TRUE);
      break;
    case PROP_BIN_WEIGHT:
      ncm_vector_substitute (&plik->bin_weight, g_value_get_object (value), TRUE);
      break;
    case PROP_SPECTRUM_ID:
      ncm_vector_substitute (&plik->spectrum_id, g_value_get_object (value), TRUE);
      break;
    case PROP_LMAX:
      plik->lmax = g_value_get_uint (value);
      break;
    case PROP_CALIB_NAME:
      g_clear_pointer (&plik->calib_name, g_free);
      plik->calib_name   = g_value_dup_string (value);
      plik->calib_pi_set = FALSE;
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_data_planck_plik_lite_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataPlanckPlikLite *plik = NC_DATA_PLANCK_PLIK_LITE (object);

  g_return_if_fail (NC_IS_DATA_PLANCK_PLIK_LITE (object));

  switch (prop_id)
  {
    case PROP_PERT_BOLTZMANN:
      g_value_set_object (value, plik->pb);
      break;
    case PROP_BIN_LMIN:
      g_value_set_object (value, plik->bin_lmin);
      break;
    case PROP_BIN_LMAX:
      g_value_set_object (value, plik->bin_lmax);
      break;
    case PROP_BIN_WEIGHT:
      g_value_set_object (value, plik->bin_weight);
      break;
    case PROP_SPECTRUM_ID:
      g_value_set_object (value, plik->spectrum_id);
      break;
    case PROP_LMAX:
      g_value_set_uint (value, plik->lmax);
      break;
    case PROP_CALIB_NAME:
      g_value_set_string (value, plik->calib_name);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_data_planck_plik_lite_dispose (GObject *object)
{
  NcDataPlanckPlikLite *plik = NC_DATA_PLANCK_PLIK_LITE (object);
  guint i;

  nc_hipert_boltzmann_clear (&plik->pb);
  ncm_vector_clear (&plik->bin_lmin);
  ncm_vector_clear (&plik->bin_lmax);
  ncm_vector_clear (&plik->bin_weight);
  ncm_vector_clear (&plik->spectrum_id);

  for (i = 0; i < NC_DATA_PLIK_LITE_SPEC_LEN; i++)
    ncm_vector_clear (&plik->cl_cache[i]);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_planck_plik_lite_parent_class)->dispose (object);
}

static void
nc_data_planck_plik_lite_finalize (GObject *object)
{
  NcDataPlanckPlikLite *plik = NC_DATA_PLANCK_PLIK_LITE (object);

  g_clear_pointer (&plik->calib_name, g_free);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_planck_plik_lite_parent_class)->finalize (object);
}

static void _nc_data_planck_plik_lite_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_planck_plik_lite_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp);

static void
nc_data_planck_plik_lite_class_init (NcDataPlanckPlikLiteClass *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class          = NCM_DATA_CLASS (klass);
  NcmDataGaussCovClass *gauss_class = NCM_DATA_GAUSS_COV_CLASS (klass);

  object_class->set_property = &nc_data_planck_plik_lite_set_property;
  object_class->get_property = &nc_data_planck_plik_lite_get_property;
  object_class->dispose      = &nc_data_planck_plik_lite_dispose;
  object_class->finalize     = &nc_data_planck_plik_lite_finalize;

  g_object_class_install_property (object_class,
                                   PROP_PERT_BOLTZMANN,
                                   g_param_spec_object ("hipert-boltzmann",
                                                        NULL,
                                                        "Perturbations (Cls source)",
                                                        NC_TYPE_HIPERT_BOLTZMANN,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_BIN_LMIN,
                                   g_param_spec_object ("bin-lmin",
                                                        NULL,
                                                        "Per-bin lower multipole (absolute ell)",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_BIN_LMAX,
                                   g_param_spec_object ("bin-lmax",
                                                        NULL,
                                                        "Per-bin upper multipole (absolute ell)",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_BIN_WEIGHT,
                                   g_param_spec_object ("bin-weight",
                                                        NULL,
                                                        "Flattened per-bin averaging weights",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SPECTRUM_ID,
                                   g_param_spec_object ("spectrum-id",
                                                        NULL,
                                                        "Per-bin spectrum tag (0=TT,1=EE,2=TE)",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LMAX,
                                   g_param_spec_uint ("lmax",
                                                      NULL,
                                                      "Maximum multipole required from the Boltzmann code",
                                                      0, G_MAXUINT, 2508,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_CALIB_NAME,
                                   g_param_spec_string ("calib-name",
                                                        NULL,
                                                        "Nuisance parameter name for the absolute calibration",
                                                        "A_planck",
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->prepare    = &_nc_data_planck_plik_lite_prepare;
  gauss_class->mean_func = &_nc_data_planck_plik_lite_mean_func;
}

static void
_nc_data_planck_plik_lite_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataPlanckPlikLite *plik = NC_DATA_PLANCK_PLIK_LITE (data);
  NcHICosmo *cosmo           = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));

  if (plik->pb == NULL)
    g_error ("_nc_data_planck_plik_lite_prepare: cannot prepare without a NcHIPertBoltzmann object. "
             "Use nc_data_planck_plik_lite_set_hipert_boltzmann to set it.");

  if (cosmo == NULL)
    g_error ("_nc_data_planck_plik_lite_prepare: cannot prepare without a NcHICosmo object in the NcmMSet.");

  nc_hipert_boltzmann_prepare_if_needed (plik->pb, cosmo);
}

static NcmVector *
_nc_data_planck_plik_lite_get_cls (NcDataPlanckPlikLite *plik, NcDataPlikLiteSpec spec)
{
  if (plik->cl_cache[spec] == NULL)
    plik->cl_cache[spec] = ncm_vector_new (plik->lmax + 1);

  switch (spec)
  {
    case NC_DATA_PLIK_LITE_SPEC_TT:
      nc_hipert_boltzmann_get_TT_Cls (plik->pb, plik->cl_cache[spec]);
      break;
    case NC_DATA_PLIK_LITE_SPEC_EE:
      nc_hipert_boltzmann_get_EE_Cls (plik->pb, plik->cl_cache[spec]);
      break;
    case NC_DATA_PLIK_LITE_SPEC_TE:
      nc_hipert_boltzmann_get_TE_Cls (plik->pb, plik->cl_cache[spec]);
      break;
    default:                                            /* LCOV_EXCL_LINE */
      g_assert_not_reached ();                          /* LCOV_EXCL_LINE */
      break;                                            /* LCOV_EXCL_LINE */
  }

  return plik->cl_cache[spec];
}

static void
_nc_data_planck_plik_lite_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp)
{
  NcDataPlanckPlikLite *plik = NC_DATA_PLANCK_PLIK_LITE (gauss);
  NcPlanckFI *pfi            = NC_PLANCK_FI (ncm_mset_peek (mset, nc_planck_fi_id ()));
  const guint np             = ncm_data_gauss_cov_get_size (gauss);
  gboolean need[NC_DATA_PLIK_LITE_SPEC_LEN] = { FALSE, FALSE, FALSE };
  NcmVector *cls[NC_DATA_PLIK_LITE_SPEC_LEN];
  gdouble A_planck, calib2;
  guint b, woff, s;

  if (pfi == NULL)
    g_error ("_nc_data_planck_plik_lite_mean_func: a NcPlanckFI model is required in the NcmMSet.");

  if (!plik->calib_pi_set)
  {
    const gchar *name = (plik->calib_name != NULL) ? plik->calib_name : "A_planck";

    if (!ncm_model_param_index_from_name (NCM_MODEL (pfi), name, &plik->calib_pi, NULL))
      g_error ("_nc_data_planck_plik_lite_mean_func: cannot find calibration parameter `%s' in model `%s'.",
               name, ncm_model_name (NCM_MODEL (pfi)));

    plik->calib_pi_set = TRUE;
  }

  A_planck = ncm_model_param_get (NCM_MODEL (pfi), plik->calib_pi);
  calib2   = A_planck * A_planck;

  /* Which spectra are needed. */
  for (b = 0; b < np; b++)
  {
    s = (guint) ncm_vector_get (plik->spectrum_id, b);
    g_assert_cmpuint (s, <, NC_DATA_PLIK_LITE_SPEC_LEN);
    need[s] = TRUE;
  }

  for (s = 0; s < NC_DATA_PLIK_LITE_SPEC_LEN; s++)
    cls[s] = need[s] ? _nc_data_planck_plik_lite_get_cls (plik, (NcDataPlikLiteSpec) s) : NULL;

  /* Bin: model bandpower = normalized weighted average of raw Cl / A_planck^2. */
  woff = 0;

  for (b = 0; b < np; b++)
  {
    const guint spec = (guint) ncm_vector_get (plik->spectrum_id, b);
    const guint lo   = (guint) ncm_vector_get (plik->bin_lmin, b);
    const guint hi   = (guint) ncm_vector_get (plik->bin_lmax, b);
    NcmVector *Cl    = cls[spec];
    gdouble acc      = 0.0;
    guint l;

    for (l = lo; l <= hi; l++)
    {
      acc += ncm_vector_get (Cl, l) * ncm_vector_get (plik->bin_weight, woff);
      woff++;
    }

    ncm_vector_set (vp, b, acc / calib2);
  }
}

/**
 * nc_data_planck_plik_lite_new_from_file:
 * @filename: a serialized #NcDataPlanckPlikLite file.
 *
 * Creates a new #NcDataPlanckPlikLite from @filename.
 *
 * Returns: (transfer full): the newly created #NcDataPlanckPlikLite.
 */
NcDataPlanckPlikLite *
nc_data_planck_plik_lite_new_from_file (const gchar *filename)
{
  NcDataPlanckPlikLite *plik = NC_DATA_PLANCK_PLIK_LITE (ncm_serialize_global_from_file (filename));

  g_assert (NC_IS_DATA_PLANCK_PLIK_LITE (plik));

  return plik;
}

/**
 * nc_data_planck_plik_lite_set_hipert_boltzmann:
 * @plik: a #NcDataPlanckPlikLite
 * @pb: a #NcHIPertBoltzmann
 *
 * Sets the perturbations object used as the source of the theory $C_\ell$.
 *
 */
void
nc_data_planck_plik_lite_set_hipert_boltzmann (NcDataPlanckPlikLite *plik, NcHIPertBoltzmann *pb)
{
  nc_hipert_boltzmann_clear (&plik->pb);

  if (pb != NULL)
    plik->pb = nc_hipert_boltzmann_ref (pb);
}

/**
 * nc_data_planck_plik_lite_peek_hipert_boltzmann:
 * @plik: a #NcDataPlanckPlikLite
 *
 * Returns: (transfer none): the perturbations object, or %NULL.
 */
NcHIPertBoltzmann *
nc_data_planck_plik_lite_peek_hipert_boltzmann (NcDataPlanckPlikLite *plik)
{
  return plik->pb;
}
