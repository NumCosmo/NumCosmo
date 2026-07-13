/***************************************************************************
 *            nc_data_planck_smica.c
 *
 *  Fri June 27 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_planck_smica.c
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
 * NcDataPlanckSmica:
 *
 * Planck plik high-ℓ SMICA likelihood (TT), native reimplementation.
 *
 * An #NcmDataGaussCov whose data vector is the masked upper-triangular entries
 * of the observed cross-frequency bandpower covariances $\hat{R}_q$, with
 * covariance the inverse of the SMICA Gaussian criterion matrix. The model
 * vector (mean_func) assembles, per bin, the $m\times m$ cross-frequency
 * covariance $R_q$ from a CMB term plus the foreground and calibration
 * components (CIB, tSZ, kSZ, point sources, galactic dust, subpixel, beam
 * leakage, inter-calibration, absolute calibration), then extracts the same
 * masked entries. Being an #NcmDataGaussCov it supports resampling.
 *
 * The assembly reproduces the clik `smica` likelihood; the exact formulas are
 * validated against the file `check_value` in the port's Python reference.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/data/nc_data_planck_smica.h"
#include "nc/cmb/nc_planck_fi.h"
#include "nc/cmb/nc_planck_fi_cor_tt.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_serialize.h"

#include <math.h>

/* SMICA nuisance parameters, in the order NcPlanckFICorTT exposes them (== the
 * clik check_param tail order). Indices into the NcPlanckFI model. */
typedef enum _NcDataPlanckSmicaPar
{
  NC_DATA_PLANCK_SMICA_PAR_A_CIB_217 = 0,
  NC_DATA_PLANCK_SMICA_PAR_CIB_INDEX,
  NC_DATA_PLANCK_SMICA_PAR_XI_SZ_CIB,
  NC_DATA_PLANCK_SMICA_PAR_A_SZ,
  NC_DATA_PLANCK_SMICA_PAR_PS_100_100,
  NC_DATA_PLANCK_SMICA_PAR_PS_143_143,
  NC_DATA_PLANCK_SMICA_PAR_PS_143_217,
  NC_DATA_PLANCK_SMICA_PAR_PS_217_217,
  NC_DATA_PLANCK_SMICA_PAR_KSZ_NORM,
  NC_DATA_PLANCK_SMICA_PAR_GAL545_100,
  NC_DATA_PLANCK_SMICA_PAR_GAL545_143,
  NC_DATA_PLANCK_SMICA_PAR_GAL545_143_217,
  NC_DATA_PLANCK_SMICA_PAR_GAL545_217,
  NC_DATA_PLANCK_SMICA_PAR_SBPX_100_100,
  NC_DATA_PLANCK_SMICA_PAR_SBPX_143_143,
  NC_DATA_PLANCK_SMICA_PAR_SBPX_143_217,
  NC_DATA_PLANCK_SMICA_PAR_SBPX_217_217,
  NC_DATA_PLANCK_SMICA_PAR_CALIB_100T,
  NC_DATA_PLANCK_SMICA_PAR_CALIB_217T,
  NC_DATA_PLANCK_SMICA_PAR_A_PLANCK,
  NC_DATA_PLANCK_SMICA_PAR_LEN,
} NcDataPlanckSmicaPar;

struct _NcDataPlanckSmica
{
  /*< private >*/
  NcmDataGaussCov parent_instance;
  NcHIPertBoltzmann *pb;
  guint lmin;
  guint lmax;
  guint nell; /* lmax - lmin + 1 */
  guint m;
  guint m2;
  guint nbins;
  NcmVector *freqs;       /* m: [100,143,217] */
  NcmVector *a_cmb;       /* m: CMB mixing */
  NcmVector *sz_color;    /* m: tSZ colour corrections */
  NcmVector *gcib_conv;   /* m: CIB muK->MJ/sr conversion */
  NcmVector *gibxsz_conv; /* m: gibXsz conversion (100 zeroed) */
  GArray *bin_lmin;       /* nbins: per-bin ell offset (from lmin) */
  GArray *bin_lmax;       /* nbins */
  NcmVector *bin_weight;  /* flattened per-(bin,ell) weights */
  GArray *quad_idx;       /* np: flat indices into (nbins*m*m) R_q to extract */
  /* per-ell templates (converter stores sz/ksz PRE-NORMALIZED at ell=3000) */
  NcmVector *tmpl_gcib;   /* (10001, m0, m0) 4x4 layout, offset 0 */
  NcmVector *tmpl_sz;     /* ell 2.. */
  NcmVector *tmpl_ksz;    /* ell 0.. */
  NcmVector *tmpl_gibxsz; /* ell 0.. */
  NcmVector *tmpl_dust;   /* (3001, 12, 12) */
  NcmVector *tmpl_leak;   /* (3001, 12, 12) */
  NcmVector *tmpl_sbpx;   /* (3001, 12, 12) */
  /* runtime caches (not serialized) */
  NcmVector *cl_TT;
  gdouble *Rq;   /* nbins*m*m */
  gdouble *perl; /* nell*m*m scratch */
};

enum
{
  PROP_0,
  PROP_PB,
  PROP_LMIN,
  PROP_LMAX,
  PROP_M,
  PROP_NBINS,
  PROP_FREQS,
  PROP_A_CMB,
  PROP_SZ_COLOR,
  PROP_GCIB_CONV,
  PROP_GIBXSZ_CONV,
  PROP_BIN_LMIN,
  PROP_BIN_LMAX,
  PROP_BIN_WEIGHT,
  PROP_QUAD_IDX,
  PROP_TMPL_GCIB,
  PROP_TMPL_SZ,
  PROP_TMPL_KSZ,
  PROP_TMPL_GIBXSZ,
  PROP_TMPL_DUST,
  PROP_TMPL_LEAK,
  PROP_TMPL_SBPX,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataPlanckSmica, nc_data_planck_smica, NCM_TYPE_DATA_GAUSS_COV)

#define GCIB_NF 4 /* gcib/cnoise template frequency dimension */
#define CN_NF 12  /* cnoise 12x12 (4 freq x 3 TEB), T uses first 3 */

static void
nc_data_planck_smica_init (NcDataPlanckSmica *smica)
{
  smica->pb          = NULL;
  smica->lmin        = 0;
  smica->lmax        = 0;
  smica->m           = 0;
  smica->nbins       = 0;
  smica->freqs       = NULL;
  smica->a_cmb       = NULL;
  smica->sz_color    = NULL;
  smica->gcib_conv   = NULL;
  smica->gibxsz_conv = NULL;
  smica->bin_lmin    = g_array_new (FALSE, FALSE, sizeof (guint32));
  smica->bin_lmax    = g_array_new (FALSE, FALSE, sizeof (guint32));
  smica->bin_weight  = NULL;
  smica->quad_idx    = g_array_new (FALSE, FALSE, sizeof (guint32));
  smica->tmpl_gcib   = NULL;
  smica->tmpl_sz     = NULL;
  smica->tmpl_ksz    = NULL;
  smica->tmpl_gibxsz = NULL;
  smica->tmpl_dust   = NULL;
  smica->tmpl_leak   = NULL;
  smica->tmpl_sbpx   = NULL;
  smica->cl_TT       = NULL;
  smica->Rq          = NULL;
  smica->perl        = NULL;
}

static void
nc_data_planck_smica_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataPlanckSmica *s = NC_DATA_PLANCK_SMICA (object);

  g_return_if_fail (NC_IS_DATA_PLANCK_SMICA (object));

  switch (prop_id)
  {
    case PROP_PB:
      nc_data_planck_smica_set_hipert_boltzmann (s, g_value_get_object (value));
      break;
    case PROP_LMIN:
      s->lmin = g_value_get_uint (value);
      break;
    case PROP_LMAX:
      s->lmax = g_value_get_uint (value);
      break;
    case PROP_M:
      s->m = g_value_get_uint (value);
      break;
    case PROP_NBINS:
      s->nbins = g_value_get_uint (value);
      break;
    case PROP_FREQS:
      ncm_vector_substitute (&s->freqs, g_value_get_object (value), TRUE);
      break;
    case PROP_A_CMB:
      ncm_vector_substitute (&s->a_cmb, g_value_get_object (value), TRUE);
      break;
    case PROP_SZ_COLOR:
      ncm_vector_substitute (&s->sz_color, g_value_get_object (value), TRUE);
      break;
    case PROP_GCIB_CONV:
      ncm_vector_substitute (&s->gcib_conv, g_value_get_object (value), TRUE);
      break;
    case PROP_GIBXSZ_CONV:
      ncm_vector_substitute (&s->gibxsz_conv, g_value_get_object (value), TRUE);
      break;
    case PROP_BIN_LMIN:
    {
      GVariant *var = g_value_get_variant (value);

      if (var != NULL)
        ncm_cfg_array_set_variant (s->bin_lmin, var);

      break;
    }
    case PROP_BIN_LMAX:
    {
      GVariant *var = g_value_get_variant (value);

      if (var != NULL)
        ncm_cfg_array_set_variant (s->bin_lmax, var);

      break;
    }
    case PROP_BIN_WEIGHT:
      ncm_vector_substitute (&s->bin_weight, g_value_get_object (value), TRUE);
      break;
    case PROP_QUAD_IDX:
    {
      GVariant *var = g_value_get_variant (value);

      if (var != NULL)
        ncm_cfg_array_set_variant (s->quad_idx, var);

      break;
    }
    case PROP_TMPL_GCIB:
      ncm_vector_substitute (&s->tmpl_gcib, g_value_get_object (value), TRUE);
      break;
    case PROP_TMPL_SZ:
      ncm_vector_substitute (&s->tmpl_sz, g_value_get_object (value), TRUE);
      break;
    case PROP_TMPL_KSZ:
      ncm_vector_substitute (&s->tmpl_ksz, g_value_get_object (value), TRUE);
      break;
    case PROP_TMPL_GIBXSZ:
      ncm_vector_substitute (&s->tmpl_gibxsz, g_value_get_object (value), TRUE);
      break;
    case PROP_TMPL_DUST:
      ncm_vector_substitute (&s->tmpl_dust, g_value_get_object (value), TRUE);
      break;
    case PROP_TMPL_LEAK:
      ncm_vector_substitute (&s->tmpl_leak, g_value_get_object (value), TRUE);
      break;
    case PROP_TMPL_SBPX:
      ncm_vector_substitute (&s->tmpl_sbpx, g_value_get_object (value), TRUE);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_data_planck_smica_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataPlanckSmica *s = NC_DATA_PLANCK_SMICA (object);

  g_return_if_fail (NC_IS_DATA_PLANCK_SMICA (object));

  switch (prop_id)
  {
    case PROP_PB:
      g_value_set_object (value, s->pb);
      break;
    case PROP_LMIN:
      g_value_set_uint (value, s->lmin);
      break;
    case PROP_LMAX:
      g_value_set_uint (value, s->lmax);
      break;
    case PROP_M:
      g_value_set_uint (value, s->m);
      break;
    case PROP_NBINS:
      g_value_set_uint (value, s->nbins);
      break;
    case PROP_FREQS:
      g_value_set_object (value, s->freqs);
      break;
    case PROP_A_CMB:
      g_value_set_object (value, s->a_cmb);
      break;
    case PROP_SZ_COLOR:
      g_value_set_object (value, s->sz_color);
      break;
    case PROP_GCIB_CONV:
      g_value_set_object (value, s->gcib_conv);
      break;
    case PROP_GIBXSZ_CONV:
      g_value_set_object (value, s->gibxsz_conv);
      break;
    case PROP_BIN_LMIN:

      if (s->bin_lmin != NULL)
        g_value_take_variant (value, ncm_cfg_array_to_variant (s->bin_lmin, G_VARIANT_TYPE ("u")));

      break;
    case PROP_BIN_LMAX:

      if (s->bin_lmax != NULL)
        g_value_take_variant (value, ncm_cfg_array_to_variant (s->bin_lmax, G_VARIANT_TYPE ("u")));

      break;
    case PROP_BIN_WEIGHT:
      g_value_set_object (value, s->bin_weight);
      break;
    case PROP_QUAD_IDX:

      if (s->quad_idx != NULL)
        g_value_take_variant (value, ncm_cfg_array_to_variant (s->quad_idx, G_VARIANT_TYPE ("u")));

      break;
    case PROP_TMPL_GCIB:
      g_value_set_object (value, s->tmpl_gcib);
      break;
    case PROP_TMPL_SZ:
      g_value_set_object (value, s->tmpl_sz);
      break;
    case PROP_TMPL_KSZ:
      g_value_set_object (value, s->tmpl_ksz);
      break;
    case PROP_TMPL_GIBXSZ:
      g_value_set_object (value, s->tmpl_gibxsz);
      break;
    case PROP_TMPL_DUST:
      g_value_set_object (value, s->tmpl_dust);
      break;
    case PROP_TMPL_LEAK:
      g_value_set_object (value, s->tmpl_leak);
      break;
    case PROP_TMPL_SBPX:
      g_value_set_object (value, s->tmpl_sbpx);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

/* Assert a per-channel config vector is present with length @len. */
#define SMICA_ASSERT_VLEN(v, len)                           \
        G_STMT_START {                                      \
          g_assert_nonnull (v);                             \
          g_assert_cmpuint (ncm_vector_len (v), ==, (len)); \
        } G_STMT_END

static void
nc_data_planck_smica_constructed (GObject *object)
{
  NcDataPlanckSmica *s = NC_DATA_PLANCK_SMICA (object);
  const guint lmax     = s->lmax;
  guint bin_tot        = 0;
  guint b, k;

  g_assert_null (s->cl_TT);
  g_assert_null (s->Rq);
  g_assert_null (s->perl);
  g_assert_cmpuint (s->lmax, >=, s->lmin);
  g_assert_cmpuint (s->nbins, >, 0);
  g_assert_cmpuint (s->m, >, 0);

  s->nell = s->lmax - s->lmin + 1;
  s->m2   = s->m * s->m;

  /* Per-channel config: exactly m entries each. */
  SMICA_ASSERT_VLEN (s->freqs, s->m);
  SMICA_ASSERT_VLEN (s->a_cmb, s->m);
  SMICA_ASSERT_VLEN (s->sz_color, s->m);
  SMICA_ASSERT_VLEN (s->gcib_conv, s->m);
  SMICA_ASSERT_VLEN (s->gibxsz_conv, s->m);

  /* Binning: one [lo,hi] pair per bin, ordered and inside the ell window;
   * bin_weight holds one weight per (bin,ell) pair, so its length is the
   * total window width summed over bins (smica_bin_add walks it linearly). */
  g_assert_nonnull (s->bin_lmin);
  g_assert_nonnull (s->bin_lmax);
  g_assert_cmpuint (s->bin_lmin->len, ==, s->nbins);
  g_assert_cmpuint (s->bin_lmax->len, ==, s->nbins);

  for (b = 0; b < s->nbins; b++)
  {
    const guint lo = g_array_index (s->bin_lmin, guint32, b);
    const guint hi = g_array_index (s->bin_lmax, guint32, b);

    g_assert_cmpuint (lo, <=, hi);
    g_assert_cmpuint (hi, <, s->nell);
    bin_tot += hi - lo + 1;
  }

  SMICA_ASSERT_VLEN (s->bin_weight, bin_tot);

  /* Masked-entry indices: non-empty and in-bounds of the assembled R_q. */
  g_assert_nonnull (s->quad_idx);
  g_assert_cmpuint (s->quad_idx->len, >, 0);

  for (k = 0; k < s->quad_idx->len; k++)
    g_assert_cmpuint (g_array_index (s->quad_idx, guint32, k), <, s->nbins * s->m2);

  /* Templates: present and long enough for the multipoles they are indexed at
   * (gcib/dust reach fixed pivots ell=3000/200 in addition to lmax). */
  g_assert_nonnull (s->tmpl_gcib);
  g_assert_nonnull (s->tmpl_sz);
  g_assert_nonnull (s->tmpl_ksz);
  g_assert_nonnull (s->tmpl_gibxsz);
  g_assert_nonnull (s->tmpl_dust);
  g_assert_nonnull (s->tmpl_leak);
  g_assert_nonnull (s->tmpl_sbpx);
  g_assert_cmpuint (ncm_vector_len (s->tmpl_gcib), >=, (MAX (lmax, 3000) + 1) * GCIB_NF * GCIB_NF);
  g_assert_cmpuint (ncm_vector_len (s->tmpl_sz), >=, lmax - 1);     /* indexed [ell-2] */
  g_assert_cmpuint (ncm_vector_len (s->tmpl_ksz), >=, lmax + 1);    /* indexed [ell]   */
  g_assert_cmpuint (ncm_vector_len (s->tmpl_gibxsz), >=, lmax - 1); /* indexed [ell-2] */
  g_assert_cmpuint (ncm_vector_len (s->tmpl_dust), >=, (lmax + 1) * CN_NF * CN_NF);
  g_assert_cmpuint (ncm_vector_len (s->tmpl_leak), >=, (lmax + 1) * CN_NF * CN_NF);
  g_assert_cmpuint (ncm_vector_len (s->tmpl_sbpx), >=, (lmax + 1) * CN_NF * CN_NF);

  s->cl_TT = ncm_vector_new (s->lmax + 1);
  s->Rq    = g_new0 (gdouble, s->nbins * s->m2);
  s->perl  = g_new0 (gdouble, s->nell * s->m2);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_planck_smica_parent_class)->constructed (object);
}

#undef SMICA_ASSERT_VLEN

static void
nc_data_planck_smica_dispose (GObject *object)
{
  NcDataPlanckSmica *s = NC_DATA_PLANCK_SMICA (object);

  nc_hipert_boltzmann_clear (&s->pb);
  ncm_vector_clear (&s->freqs);
  ncm_vector_clear (&s->a_cmb);
  ncm_vector_clear (&s->sz_color);
  ncm_vector_clear (&s->gcib_conv);
  ncm_vector_clear (&s->gibxsz_conv);
  g_clear_pointer (&s->bin_lmin, g_array_unref);
  g_clear_pointer (&s->bin_lmax, g_array_unref);
  ncm_vector_clear (&s->bin_weight);
  g_clear_pointer (&s->quad_idx, g_array_unref);
  ncm_vector_clear (&s->tmpl_gcib);
  ncm_vector_clear (&s->tmpl_sz);
  ncm_vector_clear (&s->tmpl_ksz);
  ncm_vector_clear (&s->tmpl_gibxsz);
  ncm_vector_clear (&s->tmpl_dust);
  ncm_vector_clear (&s->tmpl_leak);
  ncm_vector_clear (&s->tmpl_sbpx);
  ncm_vector_clear (&s->cl_TT);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_planck_smica_parent_class)->dispose (object);
}

static void
nc_data_planck_smica_finalize (GObject *object)
{
  NcDataPlanckSmica *s = NC_DATA_PLANCK_SMICA (object);

  g_clear_pointer (&s->Rq, g_free);
  g_clear_pointer (&s->perl, g_free);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_planck_smica_parent_class)->finalize (object);
}

static void _nc_data_planck_smica_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_planck_smica_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp);

static void
nc_data_planck_smica_class_init (NcDataPlanckSmicaClass *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class          = NCM_DATA_CLASS (klass);
  NcmDataGaussCovClass *gauss_class = NCM_DATA_GAUSS_COV_CLASS (klass);

  object_class->set_property = &nc_data_planck_smica_set_property;
  object_class->get_property = &nc_data_planck_smica_get_property;
  object_class->constructed  = &nc_data_planck_smica_constructed;
  object_class->dispose      = &nc_data_planck_smica_dispose;
  object_class->finalize     = &nc_data_planck_smica_finalize;

  g_object_class_install_property (object_class,
                                   PROP_PB,
                                   g_param_spec_object ("hipert-boltzmann",
                                                        NULL,
                                                        "Perturbations (Cls source)",
                                                        NC_TYPE_HIPERT_BOLTZMANN,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB)); /* runtime source, mutable */
  g_object_class_install_property (object_class,
                                   PROP_LMIN,
                                   g_param_spec_uint ("lmin",
                                                      NULL,
                                                      "Minimum multipole",
                                                      0, G_MAXUINT, 30,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LMAX,
                                   g_param_spec_uint ("lmax",
                                                      NULL,
                                                      "Maximum multipole",
                                                      0, G_MAXUINT, 2508,
                                                      G_PARAM_READWRITE  | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_M,
                                   g_param_spec_uint ("m-channels",
                                                      NULL,
                                                      "Number of frequency channels",
                                                      0, G_MAXUINT, 3,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NBINS,
                                   g_param_spec_uint ("nbins",
                                                      NULL,
                                                      "Number of bandpower bins",
                                                      0, G_MAXUINT, 215,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_FREQS,
                                   g_param_spec_object ("freqs",
                                                        NULL,
                                                        "Channel frequencies (GHz)",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_A_CMB,
                                   g_param_spec_object ("a-cmb",
                                                        NULL,
                                                        "CMB mixing vector",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SZ_COLOR,
                                   g_param_spec_object ("sz-color",
                                                        NULL,
                                                        "tSZ colour corrections",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_GCIB_CONV,
                                   g_param_spec_object ("gcib-conv",
                                                        NULL,
                                                        "CIB muK->MJ/sr conversion factors",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_GIBXSZ_CONV,
                                   g_param_spec_object ("gibxsz-conv",
                                                        NULL,
                                                        "CIBxtSZ conversion factors",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_BIN_LMIN,
                                   g_param_spec_variant ("bin-lmin",
                                                         NULL,
                                                         "Per-bin lower multipole offset",
                                                         G_VARIANT_TYPE ("au"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_BIN_LMAX,
                                   g_param_spec_variant ("bin-lmax",
                                                         NULL,
                                                         "Per-bin upper multipole offset",
                                                         G_VARIANT_TYPE ("au"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_BIN_WEIGHT,
                                   g_param_spec_object ("bin-weight",
                                                        NULL,
                                                        "Flattened per-(bin,ell) binning weights",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_QUAD_IDX,
                                   g_param_spec_variant ("quad-idx",
                                                         NULL,
                                                         "Flat indices of the masked R_q entries",
                                                         G_VARIANT_TYPE ("au"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TMPL_GCIB,
                                   g_param_spec_object ("tmpl-gcib",
                                                        NULL,
                                                        "Clustered CIB template",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TMPL_SZ,
                                   g_param_spec_object ("tmpl-sz",
                                                        NULL,
                                                        "tSZ template (normalized at ell=3000)",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TMPL_KSZ,
                                   g_param_spec_object ("tmpl-ksz",
                                                        NULL,
                                                        "kSZ template (normalized at ell=3000)",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TMPL_GIBXSZ,
                                   g_param_spec_object ("tmpl-gibxsz",
                                                        NULL,
                                                        "CIBxtSZ correlation template",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TMPL_DUST,
                                   g_param_spec_object ("tmpl-dust",
                                                        NULL,
                                                        "Galactic dust template",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TMPL_LEAK,
                                   g_param_spec_object ("tmpl-leak",
                                                        NULL,
                                                        "Beam leakage template",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TMPL_SBPX,
                                   g_param_spec_object ("tmpl-sbpx",
                                                        NULL,
                                                        "Subpixel effect template",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->prepare    = &_nc_data_planck_smica_prepare;
  gauss_class->mean_func = &_nc_data_planck_smica_mean_func;
}

static void
_nc_data_planck_smica_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataPlanckSmica *s = NC_DATA_PLANCK_SMICA (data);
  NcHICosmo *cosmo     = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));

  if (s->pb == NULL)
    g_error ("_nc_data_planck_smica_prepare: no NcHIPertBoltzmann set.");

  if (cosmo == NULL)
    g_error ("_nc_data_planck_smica_prepare: no NcHICosmo in mset.");

  nc_hipert_boltzmann_prepare_if_needed (s->pb, cosmo);
}

/* Bin a per-ell m*m contribution (indexed il=ell-lmin) into Rq, accumulating. */
static void
smica_bin_add (NcDataPlanckSmica *s)
{
  guint b, il, i, bb = 0;

  for (b = 0; b < s->nbins; b++)
  {
    const guint lo = g_array_index (s->bin_lmin, guint32, b);
    const guint hi = g_array_index (s->bin_lmax, guint32, b);
    gdouble *Rqb   = s->Rq + b * s->m2;

    for (il = lo; il <= hi; il++)
    {
      const gdouble w   = ncm_vector_get (s->bin_weight, bb++);
      const gdouble *pl = s->perl + il * s->m2;

      for (i = 0; i < s->m2; i++)
        Rqb[i] += w * pl[i];
    }
  }
}

static gdouble
sz_spectrum (gdouble nu)
{
  const gdouble nu0 = 143.0;
  const gdouble x   = nu / 56.78;
  const gdouble x0  = nu0 / 56.78;

  return (2.0 - x / 2.0 / tanh (x / 2.0)) / (2.0 - x0 / 2.0 / tanh (x0 / 2.0));
}

static void
_nc_data_planck_smica_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp)
{
  NcDataPlanckSmica *s = NC_DATA_PLANCK_SMICA (gauss);
  NcPlanckFI *pfi      = NC_PLANCK_FI (ncm_mset_peek (mset, nc_planck_fi_id ()));
  const guint np       = ncm_data_gauss_cov_get_size (gauss);
  const gdouble twopi  = 2.0 * M_PI;
  gdouble p[NC_DATA_PLANCK_SMICA_PAR_LEN];
  gdouble fnu[8], fnu_gx[8];
  guint b, il, i, j, k, ell;

  if (pfi == NULL)
    g_error ("_nc_data_planck_smica_mean_func: need a NcPlanckFI model in mset.");

  /* nuisance params (NcPlanckFICorTT order == check_param tail order) */
  for (i = 0; i < NC_DATA_PLANCK_SMICA_PAR_LEN; i++)
    p[i] = ncm_model_param_get (NCM_MODEL (pfi), i);

  nc_hipert_boltzmann_get_TT_Cls (s->pb, s->cl_TT);

  memset (s->Rq, 0, sizeof (gdouble) * s->nbins * s->m2);

  /* --- CMB: P_q = binned raw Cl; Rq += P_q * a_cmb[i]*a_cmb[j] --- */
  for (il = 0; il < s->nell; il++)
  {
    const gdouble cl = ncm_vector_get (s->cl_TT, il + s->lmin);
    gdouble *pl      = s->perl + il * s->m2;

    for (i = 0; i < s->m; i++)
      for (j = 0; j < s->m; j++)
        pl[i * s->m + j] = cl * ncm_vector_get (s->a_cmb, i) * ncm_vector_get (s->a_cmb, j);
  }

  smica_bin_add (s);

  /* --- sz: A_sz * 2pi/(l(l+1)) * tmpl[l-2] * fnu_i*fnu_j --- */
  for (i = 0; i < s->m; i++)
    fnu[i] = sz_spectrum (ncm_vector_get (s->freqs, i)) * ncm_vector_get (s->sz_color, i);

  for (il = 0; il < s->nell; il++)
  {
    ell = il + s->lmin;

    const gdouble c = p[NC_DATA_PLANCK_SMICA_PAR_A_SZ] * twopi / (ell * (ell + 1.0)) *
                      ncm_vector_get (s->tmpl_sz, ell - 2);
    gdouble *pl = s->perl + il * s->m2;

    for (i = 0; i < s->m; i++)
      for (j = 0; j < s->m; j++)
        pl[i * s->m + j] = c * fnu[i] * fnu[j];
  }

  smica_bin_add (s);

  /* --- ksz: ksz_norm * 2pi/(l(l+1)) * tmpl[l] (freq-independent) --- */
  for (il = 0; il < s->nell; il++)
  {
    ell = il + s->lmin;

    const gdouble c = p[NC_DATA_PLANCK_SMICA_PAR_KSZ_NORM] * twopi / (ell * (ell + 1.0)) *
                      ncm_vector_get (s->tmpl_ksz, ell);
    gdouble *pl = s->perl + il * s->m2;

    for (i = 0; i < s->m2; i++)
      pl[i] = c;
  }

  smica_bin_add (s);

  /* --- pointsource: flat, ps_A_{ij}/(3000*3001/2pi); adds to all bins --- */
  {
    const gdouble nrm = 3000.0 * 3001.0 / twopi;
    gdouble ps[9];

    /* clik defaults every ps_A pair to 1; only these four are model params. */
    for (i = 0; i < 9; i++)
      ps[i] = 1.0 / nrm;

    ps[0 * s->m + 0] = p[NC_DATA_PLANCK_SMICA_PAR_PS_100_100] / nrm;
    ps[1 * s->m + 1] = p[NC_DATA_PLANCK_SMICA_PAR_PS_143_143] / nrm;
    ps[1 * s->m + 2] = ps[2 * s->m + 1] = p[NC_DATA_PLANCK_SMICA_PAR_PS_143_217] / nrm;
    ps[2 * s->m + 2] = p[NC_DATA_PLANCK_SMICA_PAR_PS_217_217] / nrm;

    for (b = 0; b < s->nbins; b++)
      for (i = 0; i < s->m2; i++)
        s->Rq[b * s->m2 + i] += ps[i];
  }

  /* --- gcib: rigid CIB (template GCIB_NF x GCIB_NF; irigid=217=index 2) --- */
  {
    const gdouble lp    = 3000.0;
    const gdouble *gt   = ncm_vector_data (s->tmpl_gcib);
    const gdouble conv2 = ncm_vector_get (s->gcib_conv, 2);
    const gdouble tref  = gt[(guint) lp * GCIB_NF * GCIB_NF + 2 * GCIB_NF + 2];
    gdouble A[9];

    for (i = 0; i < s->m; i++)
      for (j = 0; j < s->m; j++)
      {
        A[i * s->m + j] = (p[NC_DATA_PLANCK_SMICA_PAR_A_CIB_217] / tref *
                           ncm_vector_get (s->gcib_conv, i) * ncm_vector_get (s->gcib_conv, j) /
                           (conv2 * conv2)) / lp / (lp + 1.0) * twopi;
      }

    for (il = 0; il < s->nell; il++)
    {
      ell = il + s->lmin;

      const gdouble v = pow (ell / lp, p[NC_DATA_PLANCK_SMICA_PAR_CIB_INDEX] - (-1.3));
      gdouble *pl     = s->perl + il * s->m2;

      for (i = 0; i < s->m; i++)
        for (j = 0; j < s->m; j++)
          pl[i * s->m + j] = v * gt[ell * GCIB_NF * GCIB_NF + i * GCIB_NF + j] * A[i * s->m + j];
    }

    smica_bin_add (s);
  }

  /* --- gibXsz: CIB x tSZ cross --- */
  {
    const gdouble *gxt = ncm_vector_data (s->tmpl_gibxsz);
    gdouble A[9];

    for (i = 0; i < s->m; i++)
      fnu_gx[i] = sz_spectrum (ncm_vector_get (s->freqs, i)) * ncm_vector_get (s->sz_color, i);

    for (i = 0; i < s->m; i++)
      for (j = 0; j < s->m; j++)
        A[i * s->m + j] = -p[NC_DATA_PLANCK_SMICA_PAR_XI_SZ_CIB] * sqrt (p[NC_DATA_PLANCK_SMICA_PAR_A_SZ]) *
                          (sqrt (fnu_gx[i] * p[NC_DATA_PLANCK_SMICA_PAR_A_CIB_217] * ncm_vector_get (s->gibxsz_conv, j)) +
                           sqrt (fnu_gx[j] * p[NC_DATA_PLANCK_SMICA_PAR_A_CIB_217] * ncm_vector_get (s->gibxsz_conv, i)));

    for (il = 0; il < s->nell; il++)
    {
      ell = il + s->lmin;

      const gdouble c = gxt[ell - 2] * twopi / (ell * (ell + 1.0)); /* corr_template ell-2 */
      gdouble *pl     = s->perl + il * s->m2;

      for (i = 0; i < s->m; i++)
        for (j = 0; j < s->m; j++)
          pl[i * s->m + j] = A[i * s->m + j] * c;
    }

    smica_bin_add (s);
  }

  /* --- cnoise (dust abso=1 lpiv=200; leak/sbpx abso=0). template CN_NF^2 --- */
#define CNOISE(tmpl, get_amp, abso, lpiv)                                     \
        G_STMT_START {                                                        \
          const gdouble *ct = ncm_vector_data (tmpl);                         \
          for (il = 0; il < s->nell; il++)                                    \
          {                                                                   \
            ell = il + s->lmin;                                               \
            gdouble *pl = s->perl + il * s->m2;                               \
            for (i = 0; i < s->m; i++)                                        \
            for (j = 0; j < s->m; j++)                                        \
            {                                                                 \
              const gdouble v = get_amp;                                      \
              gdouble a;                                                      \
              if (abso)                                                       \
              a = v / (ct[(lpiv) * CN_NF * CN_NF + i * CN_NF + j] *           \
                       (lpiv) * ((lpiv) + 1.0) / twopi);                      \
              else                                                            \
              a                = v;                                           \
              pl[i * s->m + j] = ct[ell * CN_NF * CN_NF + i * CN_NF + j] * a; \
            }                                                                 \
          }                                                                   \
          smica_bin_add (s);                                                  \
        } G_STMT_END

  /* dust (gal545): pairs 00,11,12,22; others 0 */
  {
    gdouble dg[9];

    memset (dg, 0, sizeof (dg));
    dg[0 * s->m + 0] = p[NC_DATA_PLANCK_SMICA_PAR_GAL545_100];
    dg[1 * s->m + 1] = p[NC_DATA_PLANCK_SMICA_PAR_GAL545_143];
    dg[1 * s->m + 2] = dg[2 * s->m + 1] = p[NC_DATA_PLANCK_SMICA_PAR_GAL545_143_217];
    dg[2 * s->m + 2] = p[NC_DATA_PLANCK_SMICA_PAR_GAL545_217];
    CNOISE (s->tmpl_dust, dg[i * s->m + j], 1, 200);
  }
  /* beam leakage: all pairs = 1, abso=0 */
  CNOISE (s->tmpl_leak, 1.0, 0, 0);

  /* subpixel: pairs 00,11,12,22 = A_sbpx; abso=0 */
  {
    gdouble sg[9];

    memset (sg, 0, sizeof (sg));
    sg[0 * s->m + 0] = p[NC_DATA_PLANCK_SMICA_PAR_SBPX_100_100];
    sg[1 * s->m + 1] = p[NC_DATA_PLANCK_SMICA_PAR_SBPX_143_143];
    sg[1 * s->m + 2] = sg[2 * s->m + 1] = p[NC_DATA_PLANCK_SMICA_PAR_SBPX_143_217];
    sg[2 * s->m + 2] = p[NC_DATA_PLANCK_SMICA_PAR_SBPX_217_217];
    CNOISE (s->tmpl_sbpx, sg[i * s->m + j], 0, 0);
  }
#undef CNOISE

  /* --- icalTP (mult): Rq[i,j] *= 1/sqrt(cal_i)/sqrt(cal_j); 143 = ref --- */
  {
    gdouble cal[3];

    cal[0] = 1.0 / sqrt (p[NC_DATA_PLANCK_SMICA_PAR_CALIB_100T]);
    cal[1] = 1.0;
    cal[2] = 1.0 / sqrt (p[NC_DATA_PLANCK_SMICA_PAR_CALIB_217T]);

    /* --- totcal (mult): Rq *= 1/A_planck^2 --- */
    const gdouble tc = 1.0 / (p[NC_DATA_PLANCK_SMICA_PAR_A_PLANCK] * p[NC_DATA_PLANCK_SMICA_PAR_A_PLANCK]);

    for (b = 0; b < s->nbins; b++)
      for (i = 0; i < s->m; i++)
        for (j = 0; j < s->m; j++)
          s->Rq[b * s->m2 + i * s->m + j] *= cal[i] * cal[j] * tc;
  }

  /* --- extract masked entries via quad_idx --- */
  for (k = 0; k < np; k++)
    ncm_vector_set (vp, k, s->Rq[(guint) g_array_index (s->quad_idx, guint32, k)]);
}

/**
 * nc_data_planck_smica_new:
 * @lmin: lowest multipole of the theory window
 * @lmax: highest multipole of the theory window
 * @m: number of frequency channels
 * @nbins: number of bandpower bins
 * @freqs: (transfer none): channel frequencies (GHz), length @m
 * @a_cmb: (transfer none): CMB mixing vector, length @m
 * @sz_color: (transfer none): tSZ colour corrections, length @m
 * @gcib_conv: (transfer none): CIB muK->MJ/sr conversions, length @m
 * @gibxsz_conv: (transfer none): CIBxtSZ conversions, length @m
 * @bin_lmin: (element-type guint32): per-bin lower multipole offset, length @nbins
 * @bin_lmax: (element-type guint32): per-bin upper multipole offset, length @nbins
 * @bin_weight: (transfer none): flattened per-(bin,ell) binning weights
 * @quad_idx: (element-type guint32): flat indices of the masked R_q entries
 * @tmpl_gcib: (transfer none): clustered CIB template
 * @tmpl_sz: (transfer none): tSZ template (normalized at ell=3000)
 * @tmpl_ksz: (transfer none): kSZ template (normalized at ell=3000)
 * @tmpl_gibxsz: (transfer none): CIBxtSZ correlation template
 * @tmpl_dust: (transfer none): galactic dust template
 * @tmpl_leak: (transfer none): beam leakage template
 * @tmpl_sbpx: (transfer none): subpixel effect template
 *
 * Builds a fully-specified #NcDataPlanckSmica. Every argument maps to a
 * construction-only property; the object validates their mutual consistency
 * (lengths, index bounds) at construction. The theory $C_\ell$ source and the
 * observed mean/covariance are set separately (see
 * nc_data_planck_smica_set_hipert_boltzmann() and the #NcmDataGaussCov API).
 *
 * Returns: (transfer full): the new #NcDataPlanckSmica.
 */
NcDataPlanckSmica *
nc_data_planck_smica_new (guint     lmin,
                          guint     lmax,
                          guint     m,
                          guint     nbins,
                          NcmVector *freqs,
                          NcmVector *a_cmb,
                          NcmVector *sz_color,
                          NcmVector *gcib_conv,
                          NcmVector *gibxsz_conv,
                          GArray    *bin_lmin,
                          GArray    *bin_lmax,
                          NcmVector *bin_weight,
                          GArray    *quad_idx,
                          NcmVector *tmpl_gcib,
                          NcmVector *tmpl_sz,
                          NcmVector *tmpl_ksz,
                          NcmVector *tmpl_gibxsz,
                          NcmVector *tmpl_dust,
                          NcmVector *tmpl_leak,
                          NcmVector *tmpl_sbpx)
{
  GVariant *v_bin_lmin     = ncm_cfg_array_to_variant (bin_lmin, G_VARIANT_TYPE ("u"));
  GVariant *v_bin_lmax     = ncm_cfg_array_to_variant (bin_lmax, G_VARIANT_TYPE ("u"));
  GVariant *v_quad_idx     = ncm_cfg_array_to_variant (quad_idx, G_VARIANT_TYPE ("u"));
  NcDataPlanckSmica *smica = g_object_new (NC_TYPE_DATA_PLANCK_SMICA,
                                           "lmin", lmin,
                                           "lmax", lmax,
                                           "m-channels", m,
                                           "nbins", nbins,
                                           "freqs", freqs,
                                           "a-cmb", a_cmb,
                                           "sz-color", sz_color,
                                           "gcib-conv", gcib_conv,
                                           "gibxsz-conv", gibxsz_conv,
                                           "bin-lmin", v_bin_lmin,
                                           "bin-lmax", v_bin_lmax,
                                           "bin-weight", bin_weight,
                                           "quad-idx", v_quad_idx,
                                           "tmpl-gcib", tmpl_gcib,
                                           "tmpl-sz", tmpl_sz,
                                           "tmpl-ksz", tmpl_ksz,
                                           "tmpl-gibxsz", tmpl_gibxsz,
                                           "tmpl-dust", tmpl_dust,
                                           "tmpl-leak", tmpl_leak,
                                           "tmpl-sbpx", tmpl_sbpx,
                                           NULL);

  g_variant_unref (v_bin_lmin);
  g_variant_unref (v_bin_lmax);
  g_variant_unref (v_quad_idx);

  return smica;
}

/**
 * nc_data_planck_smica_new_from_file:
 * @filename: a serialized #NcDataPlanckSmica file.
 *
 * Returns: (transfer full): the new #NcDataPlanckSmica.
 */
NcDataPlanckSmica *
nc_data_planck_smica_new_from_file (const gchar *filename)
{
  NcDataPlanckSmica *s = NC_DATA_PLANCK_SMICA (ncm_serialize_global_from_file (filename));

  g_assert (NC_IS_DATA_PLANCK_SMICA (s));

  return s;
}

/**
 * nc_data_planck_smica_set_hipert_boltzmann:
 * @smica: a #NcDataPlanckSmica
 * @pb: a #NcHIPertBoltzmann
 *
 * Sets the theory $C_\ell$ source.
 */
void
nc_data_planck_smica_set_hipert_boltzmann (NcDataPlanckSmica *smica, NcHIPertBoltzmann *pb)
{
  nc_hipert_boltzmann_clear (&smica->pb);

  if (pb != NULL)
    smica->pb = nc_hipert_boltzmann_ref (pb);
}

/**
 * nc_data_planck_smica_peek_hipert_boltzmann:
 * @smica: a #NcDataPlanckSmica
 *
 * Returns: (transfer none): the perturbations object, or %NULL.
 */
NcHIPertBoltzmann *
nc_data_planck_smica_peek_hipert_boltzmann (NcDataPlanckSmica *smica)
{
  return smica->pb;
}

