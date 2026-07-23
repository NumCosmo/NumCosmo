/***************************************************************************
 *            nc_data_planck_lensing.c
 *
 *  Thu July 23 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_planck_lensing.c
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
 * NcDataPlanckLensing:
 *
 * Planck 2018 CMB lensing band-power likelihood, native port (clik itype=4).
 *
 * An #NcmDataGauss over the 9 reconstructed $C_L^{\phi\phi}$ band-powers. The
 * data vector is @pp_hat and the inverse covariance is @siginv (clik already
 * ships the inverse, so #NcmDataGauss uses it as-is). The model band-power
 * $$ b_i = \sum_L B_{iL}\, C_L^{\phi\phi}\,\frac{L^2(L+1)^2}{2\pi} - c^0_i $$
 * is optionally augmented, when @hascl selects CMB spectra, by the estimator
 * renormalization
 * $$ b_i \mathrel{+}= \sum_L R^{\phi}_{iL}\, C_L^{\phi\phi}\,\frac{L^2(L+1)^2}{2\pi}
 *      + \frac{1}{A^2}\sum_X \sum_L R^{X}_{iL}\, C_L^{X}\,\frac{L(L+1)}{2\pi}, $$
 * with $X$ running over the selected CMB spectra (TT, EE, TE, ...) and $A$ the
 * calibration. The CMB-marginalized file has @hascl all zero and @cors %NULL,
 * so $b_i$ depends on $C_L^{\phi\phi}$ alone. Being an #NcmDataGauss it can
 * resample; unlike the clik wrapper it is self-contained.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/data/nc_data_planck_lensing.h"
#include "nc/cmb/nc_planck_fi.h"
#include "nc/data/nc_data_cmb.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_serialize.h"

#include <math.h>

struct _NcDataPlanckLensing
{
  /*< private >*/
  NcmDataGauss parent_instance;
  NcHIPertBoltzmann *pb;
  guint lmax;
  guint nbins;
  gboolean has_calib;
  gchar *calib_name;
  GArray *hascl;   /* 6 flags [TT,EE,BB,TE,TB,EB] */
  NcmMatrix *bins; /* nbins x (lmax+1) */
  NcmVector *cor0; /* nbins */
  NcmMatrix *cors; /* nbins x nlt, or NULL */
  /* runtime (constructed) */
  guint nl;                      /* lmax + 1 */
  guint ncmb;                    /* number of selected CMB spectra */
  guint nlt;                     /* (1 + ncmb) * nl (or nl if ncmb == 0) */
  NcDataCMBDataType cmb_type[6]; /* selected CMB types, in canonical order */
  NcmVector *cl_phi;             /* lmax+1 scratch */
  NcmVector *cl_cmb[6];          /* lmax+1 scratch per selected CMB spectrum */
  gint calib_pi;                 /* resolved calib param index (-1 => unresolved) */
};

enum
{
  PROP_0,
  PROP_PB,
  PROP_LMAX,
  PROP_NBINS,
  PROP_HAS_CALIB,
  PROP_CALIB_NAME,
  PROP_HASCL,
  PROP_BINS,
  PROP_COR0,
  PROP_CORS,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataPlanckLensing, nc_data_planck_lensing, NCM_TYPE_DATA_GAUSS)

static void
nc_data_planck_lensing_init (NcDataPlanckLensing *lens)
{
  guint k;

  lens->pb         = NULL;
  lens->lmax       = 0;
  lens->nbins      = 0;
  lens->has_calib  = FALSE;
  lens->calib_name = NULL;
  lens->hascl      = g_array_new (FALSE, FALSE, sizeof (guint32));
  lens->bins       = NULL;
  lens->cor0       = NULL;
  lens->cors       = NULL;
  lens->nl         = 0;
  lens->ncmb       = 0;
  lens->nlt        = 0;
  lens->cl_phi     = NULL;
  lens->calib_pi   = -1;

  for (k = 0; k < 6; k++)
  {
    lens->cmb_type[k] = 0;
    lens->cl_cmb[k]   = NULL;
  }
}

/* clik hascl order [TT,EE,BB,TE,TB,EB] -> NcDataCMBDataType flag. */
static NcDataCMBDataType
lensing_cli_to_type (guint cli)
{
  switch (cli)
  {
    case 0:
      return NC_DATA_CMB_TYPE_TT;

    case 1:
      return NC_DATA_CMB_TYPE_EE;

    case 2:
      return NC_DATA_CMB_TYPE_BB;

    case 3:
      return NC_DATA_CMB_TYPE_TE;

    case 4:
      return NC_DATA_CMB_TYPE_TB;

    case 5:
      return NC_DATA_CMB_TYPE_EB;

    default:
      g_assert_not_reached ();

      return 0;
  }
}

static void
lensing_get_cls_by_type (NcHIPertBoltzmann *pb, NcDataCMBDataType type, NcmVector *cls)
{
  switch (type)
  {
    case NC_DATA_CMB_TYPE_TT:
      nc_hipert_boltzmann_get_TT_Cls (pb, cls);
      break;
    case NC_DATA_CMB_TYPE_EE:
      nc_hipert_boltzmann_get_EE_Cls (pb, cls);
      break;
    case NC_DATA_CMB_TYPE_BB:
      nc_hipert_boltzmann_get_BB_Cls (pb, cls);
      break;
    case NC_DATA_CMB_TYPE_TE:
      nc_hipert_boltzmann_get_TE_Cls (pb, cls);
      break;
    case NC_DATA_CMB_TYPE_TB:
      nc_hipert_boltzmann_get_TB_Cls (pb, cls);
      break;
    case NC_DATA_CMB_TYPE_EB:
      nc_hipert_boltzmann_get_EB_Cls (pb, cls);
      break;
    default:
      g_assert_not_reached ();
  }
}

static void
nc_data_planck_lensing_constructed (GObject *object)
{
  NcDataPlanckLensing *lens = NC_DATA_PLANCK_LENSING (object);
  guint cli;

  g_assert_cmpuint (lens->lmax, >, 0);
  g_assert_cmpuint (lens->nbins, >, 0);
  g_assert_nonnull (lens->hascl);
  g_assert_cmpuint (lens->hascl->len, ==, 6);

  lens->nl = lens->lmax + 1;

  /* Derive the ordered list of selected CMB spectra from hascl. */
  lens->ncmb = 0;

  for (cli = 0; cli < 6; cli++)
  {
    if (g_array_index (lens->hascl, guint32, cli) != 0)
    {
      lens->cmb_type[lens->ncmb] = lensing_cli_to_type (cli);
      lens->cl_cmb[lens->ncmb]   = ncm_vector_new (lens->nl);
      lens->ncmb++;
    }
  }

  lens->nlt = (lens->ncmb > 0) ? (1 + lens->ncmb) * lens->nl : lens->nl;

  /* Bins and cor0 are always present; cors only when CMB spectra are selected. */
  g_assert_nonnull (lens->bins);
  g_assert_cmpuint (ncm_matrix_nrows (lens->bins), ==, lens->nbins);
  g_assert_cmpuint (ncm_matrix_ncols (lens->bins), ==, lens->nl);
  g_assert_nonnull (lens->cor0);
  g_assert_cmpuint (ncm_vector_len (lens->cor0), ==, lens->nbins);

  if (lens->ncmb > 0)
  {
    g_assert_nonnull (lens->cors);
    g_assert_cmpuint (ncm_matrix_nrows (lens->cors), ==, lens->nbins);
    g_assert_cmpuint (ncm_matrix_ncols (lens->cors), ==, lens->nlt);
  }

  lens->cl_phi = ncm_vector_new (lens->nl);

  /* Size the Gaussian data block; the observed mean and inverse covariance are
   * installed separately (NcmDataGauss "mean"/"inv-cov"). */
  if (ncm_data_gauss_get_size (NCM_DATA_GAUSS (lens)) != lens->nbins)
    ncm_data_gauss_set_size (NCM_DATA_GAUSS (lens), lens->nbins);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_planck_lensing_parent_class)->constructed (object);
}

static void
nc_data_planck_lensing_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataPlanckLensing *lens = NC_DATA_PLANCK_LENSING (object);

  g_return_if_fail (NC_IS_DATA_PLANCK_LENSING (object));

  switch (prop_id)
  {
    case PROP_PB:
      nc_data_planck_lensing_set_hipert_boltzmann (lens, g_value_get_object (value));
      break;
    case PROP_LMAX:
      lens->lmax = g_value_get_uint (value);
      break;
    case PROP_NBINS:
      lens->nbins = g_value_get_uint (value);
      break;
    case PROP_HAS_CALIB:
      lens->has_calib = g_value_get_boolean (value);
      break;
    case PROP_CALIB_NAME:
      g_free (lens->calib_name);
      lens->calib_name = g_value_dup_string (value);
      break;
    case PROP_HASCL:
    {
      GVariant *var = g_value_get_variant (value);

      if (var != NULL)
        ncm_cfg_array_set_variant (lens->hascl, var);

      break;
    }
    case PROP_BINS:
      ncm_matrix_substitute (&lens->bins, g_value_get_object (value), TRUE);
      break;
    case PROP_COR0:
      ncm_vector_substitute (&lens->cor0, g_value_get_object (value), TRUE);
      break;
    case PROP_CORS:
      ncm_matrix_substitute (&lens->cors, g_value_get_object (value), TRUE);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_data_planck_lensing_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataPlanckLensing *lens = NC_DATA_PLANCK_LENSING (object);

  g_return_if_fail (NC_IS_DATA_PLANCK_LENSING (object));

  switch (prop_id)
  {
    case PROP_PB:
      g_value_set_object (value, lens->pb);
      break;
    case PROP_LMAX:
      g_value_set_uint (value, lens->lmax);
      break;
    case PROP_NBINS:
      g_value_set_uint (value, lens->nbins);
      break;
    case PROP_HAS_CALIB:
      g_value_set_boolean (value, lens->has_calib);
      break;
    case PROP_CALIB_NAME:
      g_value_set_string (value, lens->calib_name);
      break;
    case PROP_HASCL:

      if (lens->hascl != NULL)
        g_value_take_variant (value, ncm_cfg_array_to_variant (lens->hascl, G_VARIANT_TYPE ("u")));

      break;
    case PROP_BINS:
      g_value_set_object (value, lens->bins);
      break;
    case PROP_COR0:
      g_value_set_object (value, lens->cor0);
      break;
    case PROP_CORS:
      g_value_set_object (value, lens->cors);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_data_planck_lensing_dispose (GObject *object)
{
  NcDataPlanckLensing *lens = NC_DATA_PLANCK_LENSING (object);
  guint k;

  nc_hipert_boltzmann_clear (&lens->pb);
  ncm_matrix_clear (&lens->bins);
  ncm_vector_clear (&lens->cor0);
  ncm_matrix_clear (&lens->cors);
  ncm_vector_clear (&lens->cl_phi);

  for (k = 0; k < 6; k++)
    ncm_vector_clear (&lens->cl_cmb[k]);

  g_clear_pointer (&lens->hascl, g_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_planck_lensing_parent_class)->dispose (object);
}

static void
nc_data_planck_lensing_finalize (GObject *object)
{
  NcDataPlanckLensing *lens = NC_DATA_PLANCK_LENSING (object);

  g_clear_pointer (&lens->calib_name, g_free);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_planck_lensing_parent_class)->finalize (object);
}

static void _nc_data_planck_lensing_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_planck_lensing_mean_func (NcmDataGauss *gauss, NcmMSet *mset, NcmVector *vp);

static void
nc_data_planck_lensing_class_init (NcDataPlanckLensingClass *klass)
{
  GObjectClass *object_class     = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class       = NCM_DATA_CLASS (klass);
  NcmDataGaussClass *gauss_class = NCM_DATA_GAUSS_CLASS (klass);

  object_class->set_property = &nc_data_planck_lensing_set_property;
  object_class->get_property = &nc_data_planck_lensing_get_property;
  object_class->constructed  = &nc_data_planck_lensing_constructed;
  object_class->dispose      = &nc_data_planck_lensing_dispose;
  object_class->finalize     = &nc_data_planck_lensing_finalize;

  g_object_class_install_property (object_class,
                                   PROP_PB,
                                   g_param_spec_object ("hipert-boltzmann",
                                                        NULL,
                                                        "Perturbations (Cls source)",
                                                        NC_TYPE_HIPERT_BOLTZMANN,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LMAX,
                                   g_param_spec_uint ("lmax",
                                                      NULL,
                                                      "Maximum multipole",
                                                      1, G_MAXUINT, 2500,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NBINS,
                                   g_param_spec_uint ("nbins",
                                                      NULL,
                                                      "Number of band-power bins",
                                                      1, G_MAXUINT, 9,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_HAS_CALIB,
                                   g_param_spec_boolean ("has-calib",
                                                         NULL,
                                                         "Whether a calibration parameter scales the CMB renormalization",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_CALIB_NAME,
                                   g_param_spec_string ("calib-name",
                                                        NULL,
                                                        "Calibration parameter name",
                                                        "A_planck",
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_HASCL,
                                   g_param_spec_variant ("hascl",
                                                         NULL,
                                                         "CMB spectra selection flags [TT,EE,BB,TE,TB,EB]",
                                                         G_VARIANT_TYPE ("au"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_BINS,
                                   g_param_spec_object ("bins",
                                                        NULL,
                                                        "Binning matrix (nbins x lmax+1)",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_COR0,
                                   g_param_spec_object ("cor0",
                                                        NULL,
                                                        "Per-bin constant renormalization offset (nbins)",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_CORS,
                                   g_param_spec_object ("cors",
                                                        NULL,
                                                        "Renormalization response matrix (nbins x nlt), or NULL",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->prepare    = &_nc_data_planck_lensing_prepare;
  gauss_class->mean_func = &_nc_data_planck_lensing_mean_func;
}

static void
_nc_data_planck_lensing_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataPlanckLensing *lens = NC_DATA_PLANCK_LENSING (data);
  NcHICosmo *cosmo          = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  guint k;

  if (lens->pb == NULL)
    g_error ("_nc_data_planck_lensing_prepare: no NcHIPertBoltzmann set.");

  if (cosmo == NULL)
    g_error ("_nc_data_planck_lensing_prepare: no NcHICosmo in mset.");

  /* Ensure the theory targets and their lmax cover what the model needs.
   * Idempotent: append is a bit-or, lmax is only grown. */
  nc_hipert_boltzmann_append_target_Cls (lens->pb, NC_DATA_CMB_TYPE_PHIPHI);

  if (nc_hipert_boltzmann_get_PHIPHI_lmax (lens->pb) < lens->lmax)
    nc_hipert_boltzmann_set_PHIPHI_lmax (lens->pb, lens->lmax);

  for (k = 0; k < lens->ncmb; k++)
  {
    nc_hipert_boltzmann_append_target_Cls (lens->pb, lens->cmb_type[k]);

    switch (lens->cmb_type[k])
    {
      case NC_DATA_CMB_TYPE_TT:

        if (nc_hipert_boltzmann_get_TT_lmax (lens->pb) < lens->lmax)
          nc_hipert_boltzmann_set_TT_lmax (lens->pb, lens->lmax);

        break;
      case NC_DATA_CMB_TYPE_EE:

        if (nc_hipert_boltzmann_get_EE_lmax (lens->pb) < lens->lmax)
          nc_hipert_boltzmann_set_EE_lmax (lens->pb, lens->lmax);

        break;
      case NC_DATA_CMB_TYPE_TE:

        if (nc_hipert_boltzmann_get_TE_lmax (lens->pb) < lens->lmax)
          nc_hipert_boltzmann_set_TE_lmax (lens->pb, lens->lmax);

        break;
      default:
        break;
    }
  }

  nc_hipert_boltzmann_prepare_if_needed (lens->pb, cosmo);
}

/* Resolve the calibration factor 1/A^2 from the NcPlanckFI model in @mset. */
static gdouble
lensing_calib (NcDataPlanckLensing *lens, NcmMSet *mset)
{
  NcPlanckFI *pfi = NC_PLANCK_FI (ncm_mset_peek (mset, nc_planck_fi_id ()));
  gdouble cal;

  if (pfi == NULL)
    g_error ("lensing_calib: calibration '%s' requested but no NcPlanckFI in mset.", lens->calib_name);

  if (lens->calib_pi < 0)
  {
    guint idx;

    if (!ncm_model_param_index_from_name (NCM_MODEL (pfi), lens->calib_name, &idx, NULL))
      g_error ("lensing_calib: calibration parameter '%s' not in model.", lens->calib_name);

    lens->calib_pi = idx;
  }

  cal = ncm_model_param_get (NCM_MODEL (pfi), lens->calib_pi);

  return 1.0 / (cal * cal);
}

static void
_nc_data_planck_lensing_mean_func (NcmDataGauss *gauss, NcmMSet *mset, NcmVector *vp)
{
  NcDataPlanckLensing *lens = NC_DATA_PLANCK_LENSING (gauss);
  const guint lmax          = lens->lmax;
  const gdouble *phi        = NULL;
  gdouble calib             = 1.0;
  guint ib, l, k;

  nc_hipert_boltzmann_get_PHIPHI_Cls (lens->pb, lens->cl_phi);
  phi = ncm_vector_data (lens->cl_phi);

  for (k = 0; k < lens->ncmb; k++)
    lensing_get_cls_by_type (lens->pb, lens->cmb_type[k], lens->cl_cmb[k]);

  if ((lens->cors != NULL) && lens->has_calib)
    calib = lensing_calib (lens, mset);

  for (ib = 0; ib < lens->nbins; ib++)
  {
    gdouble b_ib = 0.0;

    /* Band-power projection: b_i = sum_L bins[i,L] Clpp[L] L^2(L+1)^2/2pi. */
    for (l = 0; l <= lmax; l++)
    {
      const gdouble w = (gdouble) (l * l) * (l + 1.0) * (l + 1.0) / 2.0 / M_PI;

      b_ib += ncm_matrix_get (lens->bins, ib, l) * phi[l] * w;
    }

    b_ib -= ncm_vector_get (lens->cor0, ib);

    if (lens->cors != NULL)
    {
      gdouble bb = 0.0;

      /* phi renormalization block (cors columns 0..lmax). */
      for (l = 0; l <= lmax; l++)
      {
        const gdouble w = (gdouble) (l * l) * (l + 1.0) * (l + 1.0) / 2.0 / M_PI;

        bb += ncm_matrix_get (lens->cors, ib, l) * phi[l] * w;
      }

      /* CMB renormalization blocks, scaled by the calibration. */
      for (k = 0; k < lens->ncmb; k++)
      {
        const guint i0     = (k + 1) * lens->nl;
        const gdouble *clx = ncm_vector_data (lens->cl_cmb[k]);

        for (l = 0; l <= lmax; l++)
        {
          const gdouble w = (l * (l + 1.0)) / 2.0 / M_PI;

          bb += ncm_matrix_get (lens->cors, ib, i0 + l) * clx[l] * w * calib;
        }
      }

      b_ib += bb;
    }

    ncm_vector_set (vp, ib, b_ib);
  }
}

/**
 * nc_data_planck_lensing_new:
 * @lmax: maximum multipole
 * @nbins: number of band-power bins
 * @has_calib: whether a calibration parameter scales the CMB renormalization
 * @hascl: (element-type guint32): 6 CMB spectra flags [TT,EE,BB,TE,TB,EB]
 * @bins: (transfer none): binning matrix (nbins x lmax+1)
 * @cor0: (transfer none): per-bin constant offset (length nbins)
 * @cors: (transfer none) (nullable): renormalization response (nbins x nlt), or %NULL
 *
 * Builds a fully-specified #NcDataPlanckLensing. The observed band-powers and
 * inverse covariance are set through the #NcmDataGauss "mean" and "inv-cov"
 * properties; the theory $C_\ell$ source is set separately.
 *
 * Returns: (transfer full): the new #NcDataPlanckLensing.
 */
NcDataPlanckLensing *
nc_data_planck_lensing_new (guint     lmax,
                            guint     nbins,
                            gboolean  has_calib,
                            GArray    *hascl,
                            NcmMatrix *bins,
                            NcmVector *cor0,
                            NcmMatrix *cors)
{
  GVariant *hascl_var       = ncm_cfg_array_to_variant (hascl, G_VARIANT_TYPE ("u"));
  NcDataPlanckLensing *lens = g_object_new (NC_TYPE_DATA_PLANCK_LENSING,
                                            "lmax", lmax,
                                            "nbins", nbins,
                                            "has-calib", has_calib,
                                            "hascl", hascl_var,
                                            "bins", bins,
                                            "cor0", cor0,
                                            "cors", cors,
                                            NULL);

  g_variant_unref (hascl_var);

  return lens;
}

/**
 * nc_data_planck_lensing_new_from_file:
 * @filename: a serialized #NcDataPlanckLensing file.
 *
 * Returns: (transfer full): the new #NcDataPlanckLensing.
 */
NcDataPlanckLensing *
nc_data_planck_lensing_new_from_file (const gchar *filename)
{
  NcDataPlanckLensing *lens = NC_DATA_PLANCK_LENSING (ncm_serialize_global_from_file (filename));

  g_assert (NC_IS_DATA_PLANCK_LENSING (lens));

  return lens;
}

/**
 * nc_data_planck_lensing_set_hipert_boltzmann:
 * @lens: a #NcDataPlanckLensing
 * @pb: a #NcHIPertBoltzmann
 *
 * Sets the theory $C_\ell$ source.
 */
void
nc_data_planck_lensing_set_hipert_boltzmann (NcDataPlanckLensing *lens, NcHIPertBoltzmann *pb)
{
  nc_hipert_boltzmann_clear (&lens->pb);

  if (pb != NULL)
    lens->pb = nc_hipert_boltzmann_ref (pb);
}

/**
 * nc_data_planck_lensing_peek_hipert_boltzmann:
 * @lens: a #NcDataPlanckLensing
 *
 * Returns: (transfer none): the perturbations object, or %NULL.
 */
NcHIPertBoltzmann *
nc_data_planck_lensing_peek_hipert_boltzmann (NcDataPlanckLensing *lens)
{
  return lens->pb;
}

