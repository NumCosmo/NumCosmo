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
 * An #NcmDataGauss whose data vector is the masked upper-triangular entries
 * of the observed cross-frequency bandpower covariances $\hat{R}_q$, with
 * inverse covariance the SMICA Gaussian criterion matrix (clik ships the
 * criterion directly, so #NcmDataGauss uses it as-is -- no inversion). The model
 * vector (mean_func) assembles, per bin, the $m\times m$ cross-frequency
 * covariance $R_q$ from a CMB term plus the foreground and calibration
 * components (CIB, tSZ, kSZ, point sources, galactic dust, subpixel, beam
 * leakage, inter-calibration, absolute calibration), then extracts the same
 * masked entries. Being an #NcmDataGauss it supports resampling.
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
  NcmDataGauss parent_instance;
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
  /* polarization (TTTEEE); all NULL/empty for TT (=> has_pol FALSE) */
  GArray *field;       /* m: per-channel field 0=T, 1=E (empty => all T) */
  NcmVector *tmpl_e2e; /* (3001,12,12) EE end-to-end cnoise; NULL for TT */
  GArray *ical_im;     /* icalTP calibrated-map indices */
  NcmVector *ical_w;   /* icalTP mixing weights, length m*m*2 */
  GArray *ical_other;  /* icalTP mixing other-map indices, length m*m*2 */
  /* runtime caches (not serialized) */
  NcmVector *cl_TT;
  NcmVector *cl_EE; /* only allocated when has_pol */
  NcmVector *cl_TE; /* only allocated when has_pol */
  guint mT;         /* number of T channels (== m when no field) */
  gboolean has_pol;
  guint cn[16];    /* channel -> cnoise 12x12 index (teb-major) */
  gint pol_pi[27]; /* resolved NcPlanckFI indices for pol params */
  gboolean pol_pi_set;
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
  PROP_FIELD,
  PROP_TMPL_E2E,
  PROP_ICAL_IM,
  PROP_ICAL_W,
  PROP_ICAL_OTHER,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataPlanckSmica, nc_data_planck_smica, NCM_TYPE_DATA_GAUSS)

#define GCIB_NF 4 /* gcib/cnoise template frequency dimension */
#define CN_NF 12  /* cnoise 12x12 (4 freq x 3 TEB), T uses first 3 */

/* Polarization nuisance params, resolved by name in the NcPlanckFICorTTTEEE
 * model (indices 0..19 are shared with TT and read positionally). Slot order
 * matches the PP_* accessors below. */
static const gchar *nc_data_planck_smica_pol_par[] = {
  "galf_EE_A_100", "galf_EE_A_100_143", "galf_EE_A_100_217",
  "galf_EE_A_143", "galf_EE_A_143_217", "galf_EE_A_217", "galf_EE_index",
  "galf_TE_A_100", "galf_TE_A_100_143", "galf_TE_A_100_217",
  "galf_TE_A_143", "galf_TE_A_143_217", "galf_TE_A_217", "galf_TE_index",
  "A_cnoise_e2e_100_100_EE", "A_cnoise_e2e_143_143_EE", "A_cnoise_e2e_217_217_EE",
  "A_sbpx_100_100_EE", "A_sbpx_100_143_EE", "A_sbpx_100_217_EE",
  "A_sbpx_143_143_EE", "A_sbpx_143_217_EE", "A_sbpx_217_217_EE",
  "calib_100P", "calib_143P", "calib_217P", "A_pol",
};

#define NC_SMICA_POL_NP (G_N_ELEMENTS (nc_data_planck_smica_pol_par))
#define PP_GALF_EE_A(k) (0 + (k)) /* k=0..5 */
#define PP_GALF_EE_INDEX 6
#define PP_GALF_TE_A(k) (7 + (k)) /* k=0..5 */
#define PP_GALF_TE_INDEX 13
#define PP_E2E(k) (14 + (k))     /* k=0..2 */
#define PP_SBPX_EE(k) (17 + (k)) /* k=0..5 */
#define PP_CALIB_P(k) (23 + (k)) /* k=0..2 (100P,143P,217P) */
#define PP_A_POL 26

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
  smica->field       = g_array_new (FALSE, FALSE, sizeof (guint32));
  smica->tmpl_e2e    = NULL;
  smica->ical_im     = g_array_new (FALSE, FALSE, sizeof (guint32));
  smica->ical_w      = NULL;
  smica->ical_other  = g_array_new (FALSE, FALSE, sizeof (guint32));
  smica->cl_TT       = NULL;
  smica->cl_EE       = NULL;
  smica->cl_TE       = NULL;
  smica->mT          = 0;
  smica->has_pol     = FALSE;
  smica->pol_pi_set  = FALSE;
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
    case PROP_FIELD:
    {
      GVariant *var = g_value_get_variant (value);

      if (var != NULL)
        ncm_cfg_array_set_variant (s->field, var);

      break;
    }
    case PROP_TMPL_E2E:
      ncm_vector_substitute (&s->tmpl_e2e, g_value_get_object (value), TRUE);
      break;
    case PROP_ICAL_IM:
    {
      GVariant *var = g_value_get_variant (value);

      if (var != NULL)
        ncm_cfg_array_set_variant (s->ical_im, var);

      break;
    }
    case PROP_ICAL_W:
      ncm_vector_substitute (&s->ical_w, g_value_get_object (value), TRUE);
      break;
    case PROP_ICAL_OTHER:
    {
      GVariant *var = g_value_get_variant (value);

      if (var != NULL)
        ncm_cfg_array_set_variant (s->ical_other, var);

      break;
    }
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
    case PROP_FIELD:

      if (s->field != NULL)
        g_value_take_variant (value, ncm_cfg_array_to_variant (s->field, G_VARIANT_TYPE ("u")));

      break;
    case PROP_TMPL_E2E:
      g_value_set_object (value, s->tmpl_e2e);
      break;
    case PROP_ICAL_IM:

      if (s->ical_im != NULL)
        g_value_take_variant (value, ncm_cfg_array_to_variant (s->ical_im, G_VARIANT_TYPE ("u")));

      break;
    case PROP_ICAL_W:
      g_value_set_object (value, s->ical_w);
      break;
    case PROP_ICAL_OTHER:

      if (s->ical_other != NULL)
        g_value_take_variant (value, ncm_cfg_array_to_variant (s->ical_other, G_VARIANT_TYPE ("u")));

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

  /* Field descriptor (optional): empty => all T (TT). Derive mT / has_pol and
   * the channel -> cnoise 12x12 index map (teb-major: T at 0..3, E at 4..7). */
  g_assert_cmpuint (s->m, <=, G_N_ELEMENTS (s->cn));

  if ((s->field != NULL) && (s->field->len > 0))
  {
    guint nt = 0, ne = 0;

    g_assert_cmpuint (s->field->len, ==, s->m);

    for (k = 0; k < s->m; k++)
    {
      const guint fk = g_array_index (s->field, guint32, k);

      g_assert_cmpuint (fk, <=, 1);

      if (fk == 0)
        s->cn[k] = nt++;
      else
        s->cn[k] = 4 + ne++;
    }

    s->mT      = nt;
    s->has_pol = (nt < s->m);
    /* This port assumes the T channels precede the E channels. */
    g_assert_cmpuint (s->mT + ne, ==, s->m);
  }
  else
  {
    for (k = 0; k < s->m; k++)
      s->cn[k] = k;

    s->mT      = s->m;
    s->has_pol = FALSE;
  }

  if (s->has_pol)
  {
    /* Polarization payload must be complete and consistently sized. */
    g_assert_nonnull (s->tmpl_e2e);
    g_assert_cmpuint (ncm_vector_len (s->tmpl_e2e), >=, (lmax + 1) * CN_NF * CN_NF);
    SMICA_ASSERT_VLEN (s->ical_w, 2 * s->m2);
    g_assert_nonnull (s->ical_im);
    g_assert_cmpuint (s->ical_im->len, >, 0);
    g_assert_nonnull (s->ical_other);
    g_assert_cmpuint (s->ical_other->len, ==, 2 * s->m2);

    s->cl_EE = ncm_vector_new (s->lmax + 1);
    s->cl_TE = ncm_vector_new (s->lmax + 1);
  }

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
  g_clear_pointer (&s->field, g_array_unref);
  ncm_vector_clear (&s->tmpl_e2e);
  g_clear_pointer (&s->ical_im, g_array_unref);
  ncm_vector_clear (&s->ical_w);
  g_clear_pointer (&s->ical_other, g_array_unref);
  ncm_vector_clear (&s->cl_TT);
  ncm_vector_clear (&s->cl_EE);
  ncm_vector_clear (&s->cl_TE);

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
static void _nc_data_planck_smica_mean_func (NcmDataGauss *gauss, NcmMSet *mset, NcmVector *vp);

static void
nc_data_planck_smica_class_init (NcDataPlanckSmicaClass *klass)
{
  GObjectClass *object_class     = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class       = NCM_DATA_CLASS (klass);
  NcmDataGaussClass *gauss_class = NCM_DATA_GAUSS_CLASS (klass);

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
  g_object_class_install_property (object_class,
                                   PROP_FIELD,
                                   g_param_spec_variant ("field",
                                                         NULL,
                                                         "Per-channel field type (0=T, 1=E); empty => all T",
                                                         G_VARIANT_TYPE ("au"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TMPL_E2E,
                                   g_param_spec_object ("tmpl-e2e",
                                                        NULL,
                                                        "EE end-to-end correlated-noise template (polarization)",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ICAL_IM,
                                   g_param_spec_variant ("ical-im",
                                                         NULL,
                                                         "icalTP calibrated-map indices (polarization)",
                                                         G_VARIANT_TYPE ("au"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ICAL_W,
                                   g_param_spec_object ("ical-w",
                                                        NULL,
                                                        "icalTP calibration mixing weights, length m*m*2 (polarization)",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ICAL_OTHER,
                                   g_param_spec_variant ("ical-other",
                                                         NULL,
                                                         "icalTP mixing other-map indices, length m*m*2 (polarization)",
                                                         G_VARIANT_TYPE ("au"), NULL,
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

/* Add a cnoise-family component: per-ell pl[i,j] = tmpl[cn_i,cn_j] * a, where
 * a = amp (abso==0) or amp normalized at l_pivot (abso==1). @amp is an m*m
 * matrix (0 outside the component's pairs). Zeroed entries contribute 0 (no
 * template access), so a pair absent from a component adds nothing. For the
 * pol case perl is zeroed first so untouched blocks don't leak stale values;
 * for TT (mT==m, cn[i]==i) this reproduces the original inline assembly. */
static void
smica_cnoise_add (NcDataPlanckSmica *s, NcmVector *tmpl, const gdouble *amp, gboolean abso, guint lpiv)
{
  const gdouble twopi = 2.0 * M_PI;
  const gdouble *ct   = ncm_vector_data (tmpl);
  guint il, i, j, ell;

  if (s->has_pol)
    memset (s->perl, 0, sizeof (gdouble) * s->nell * s->m2);

  for (il = 0; il < s->nell; il++)
  {
    gdouble *pl = s->perl + il * s->m2;

    ell = il + s->lmin;

    for (i = 0; i < s->m; i++)
      for (j = 0; j < s->m; j++)
      {
        const gdouble v = amp[i * s->m + j];
        gdouble a;

        if (v == 0.0)
        {
          pl[i * s->m + j] = 0.0;
          continue;
        }

        if (abso)
          a = v / (ct[lpiv * CN_NF * CN_NF + s->cn[i] * CN_NF + s->cn[j]] *
                   lpiv * (lpiv + 1.0) / twopi);
        else
          a = v;

        pl[i * s->m + j] = ct[ell * CN_NF * CN_NF + s->cn[i] * CN_NF + s->cn[j]] * a;
      }
  }

  smica_bin_add (s);
}

/* Add a powerlaw-free-emissivity (galactic polarization dust) component:
 * pl[i,j] = (amp[i,j]/nrmit) * (ell/l_pivot)^index, l_pivot=500, l2norm=1 =>
 * nrmit = l_pivot(l_pivot+1)/2pi. Pol-only, so perl is always zeroed first. */
static void
smica_pwfe_add (NcDataPlanckSmica *s, const gdouble *amp, gdouble index)
{
  const gdouble twopi = 2.0 * M_PI;
  const gdouble lpiv  = 500.0;
  const gdouble nrmit = lpiv * (lpiv + 1.0) / twopi;
  guint il, i, j;

  memset (s->perl, 0, sizeof (gdouble) * s->nell * s->m2);

  for (il = 0; il < s->nell; il++)
  {
    const gdouble v = pow ((il + s->lmin) / lpiv, index);
    gdouble *pl     = s->perl + il * s->m2;

    for (i = 0; i < s->m; i++)
      for (j = 0; j < s->m; j++)
        if (amp[i * s->m + j] != 0.0)
          pl[i * s->m + j] = (amp[i * s->m + j] / nrmit) * v;
  }

  smica_bin_add (s);
}

/* Upper-triangular pair index for @a<=@b over @n frequencies (matches the
 * galf_EE/TE A_{f}[_{f2}] parameter ordering). */
static guint
smica_utri (guint a, guint b, guint n)
{
  return a * n - a * (a - 1) / 2 + (b - a);
}

static void
_nc_data_planck_smica_mean_func (NcmDataGauss *gauss, NcmMSet *mset, NcmVector *vp)
{
  NcDataPlanckSmica *s = NC_DATA_PLANCK_SMICA (gauss);
  NcPlanckFI *pfi      = NC_PLANCK_FI (ncm_mset_peek (mset, nc_planck_fi_id ()));
  const guint np       = ncm_data_gauss_get_size (gauss);
  const guint m        = s->m;
  const guint mT       = s->mT;
  const gdouble twopi  = 2.0 * M_PI;
  gdouble p[NC_DATA_PLANCK_SMICA_PAR_LEN];
  gdouble pp[NC_SMICA_POL_NP];
  gdouble fnu[16], fnu_gx[16];
  const guint32 *fld = NULL;
  guint b, il, i, j, k, ell;

  if (pfi == NULL)
    g_error ("_nc_data_planck_smica_mean_func: need a NcPlanckFI model in mset.");

  /* Shared TT-order nuisance params (indices 0..19 identical in TT & TTTEEE). */
  for (i = 0; i < NC_DATA_PLANCK_SMICA_PAR_LEN; i++)
    p[i] = ncm_model_param_get (NCM_MODEL (pfi), i);

  if (s->has_pol)
  {
    fld = (const guint32 *) s->field->data;

    if (!s->pol_pi_set)
    {
      for (k = 0; k < NC_SMICA_POL_NP; k++)
      {
        guint idx;

        if (!ncm_model_param_index_from_name (NCM_MODEL (pfi), nc_data_planck_smica_pol_par[k], &idx, NULL))
          g_error ("_nc_data_planck_smica_mean_func: pol parameter '%s' not in model.",
                   nc_data_planck_smica_pol_par[k]);

        s->pol_pi[k] = idx;
      }

      s->pol_pi_set = TRUE;
    }

    for (k = 0; k < NC_SMICA_POL_NP; k++)
      pp[k] = ncm_model_param_get (NCM_MODEL (pfi), s->pol_pi[k]);
  }

  nc_hipert_boltzmann_get_TT_Cls (s->pb, s->cl_TT);

  if (s->has_pol)
  {
    nc_hipert_boltzmann_get_EE_Cls (s->pb, s->cl_EE);
    nc_hipert_boltzmann_get_TE_Cls (s->pb, s->cl_TE);
  }

  memset (s->Rq, 0, sizeof (gdouble) * s->nbins * s->m2);

  /* --- CMB: Rq += a_cmb_i*a_cmb_j * C_l^{field(i),field(j)} --- */
  for (il = 0; il < s->nell; il++)
  {
    const gdouble clTT = ncm_vector_get (s->cl_TT, il + s->lmin);
    const gdouble clEE = s->has_pol ? ncm_vector_get (s->cl_EE, il + s->lmin) : 0.0;
    const gdouble clTE = s->has_pol ? ncm_vector_get (s->cl_TE, il + s->lmin) : 0.0;
    gdouble *pl        = s->perl + il * s->m2;

    for (i = 0; i < m; i++)
      for (j = 0; j < m; j++)
      {
        gdouble cl = clTT;

        if (s->has_pol)
          cl = (fld[i] == 0) ? ((fld[j] == 0) ? clTT : clTE)
                             : ((fld[j] == 0) ? clTE : clEE);

        pl[i * m + j] = cl * ncm_vector_get (s->a_cmb, i) * ncm_vector_get (s->a_cmb, j);
      }
  }

  smica_bin_add (s);

  /* --- sz: A_sz * 2pi/(l(l+1)) * tmpl[l-2] * fnu_i*fnu_j (temperature) --- */
  for (i = 0; i < mT; i++)
    fnu[i] = sz_spectrum (ncm_vector_get (s->freqs, i)) * ncm_vector_get (s->sz_color, i);

  if (s->has_pol)
    memset (s->perl, 0, sizeof (gdouble) * s->nell * s->m2);

  for (il = 0; il < s->nell; il++)
  {
    ell = il + s->lmin;

    const gdouble c = p[NC_DATA_PLANCK_SMICA_PAR_A_SZ] * twopi / (ell * (ell + 1.0)) *
                      ncm_vector_get (s->tmpl_sz, ell - 2);
    gdouble *pl = s->perl + il * s->m2;

    for (i = 0; i < mT; i++)
      for (j = 0; j < mT; j++)
        pl[i * m + j] = c * fnu[i] * fnu[j];
  }

  smica_bin_add (s);

  /* --- ksz: ksz_norm * 2pi/(l(l+1)) * tmpl[l] (freq-indep, temperature) --- */
  if (s->has_pol)
    memset (s->perl, 0, sizeof (gdouble) * s->nell * s->m2);

  for (il = 0; il < s->nell; il++)
  {
    ell = il + s->lmin;

    const gdouble c = p[NC_DATA_PLANCK_SMICA_PAR_KSZ_NORM] * twopi / (ell * (ell + 1.0)) *
                      ncm_vector_get (s->tmpl_ksz, ell);
    gdouble *pl = s->perl + il * s->m2;

    for (i = 0; i < mT; i++)
      for (j = 0; j < mT; j++)
        pl[i * m + j] = c;
  }

  smica_bin_add (s);

  /* --- pointsource: flat, ps_A_{ij}/(3000*3001/2pi) (temperature) --- */
  {
    const gdouble nrm = 3000.0 * 3001.0 / twopi;
    gdouble ps[36];

    memset (ps, 0, sizeof (ps));

    /* clik defaults every T ps_A pair to 1; only these four are model params. */
    for (i = 0; i < mT; i++)
      for (j = 0; j < mT; j++)
        ps[i * m + j] = 1.0 / nrm;

    ps[0 * m + 0] = p[NC_DATA_PLANCK_SMICA_PAR_PS_100_100] / nrm;
    ps[1 * m + 1] = p[NC_DATA_PLANCK_SMICA_PAR_PS_143_143] / nrm;
    ps[1 * m + 2] = ps[2 * m + 1] = p[NC_DATA_PLANCK_SMICA_PAR_PS_143_217] / nrm;
    ps[2 * m + 2] = p[NC_DATA_PLANCK_SMICA_PAR_PS_217_217] / nrm;

    for (b = 0; b < s->nbins; b++)
      for (i = 0; i < s->m2; i++)
        s->Rq[b * s->m2 + i] += ps[i];
  }

  /* --- gcib: rigid CIB (4x4 template; irigid=217=index 2; temperature) --- */
  {
    const gdouble lp    = 3000.0;
    const gdouble *gt   = ncm_vector_data (s->tmpl_gcib);
    const gdouble conv2 = ncm_vector_get (s->gcib_conv, 2);
    const gdouble tref  = gt[(guint) lp * GCIB_NF * GCIB_NF + 2 * GCIB_NF + 2];
    gdouble A[16];

    for (i = 0; i < mT; i++)
      for (j = 0; j < mT; j++)
        A[i * mT + j] = (p[NC_DATA_PLANCK_SMICA_PAR_A_CIB_217] / tref *
                         ncm_vector_get (s->gcib_conv, i) * ncm_vector_get (s->gcib_conv, j) /
                         (conv2 * conv2)) / lp / (lp + 1.0) * twopi;

    if (s->has_pol)
      memset (s->perl, 0, sizeof (gdouble) * s->nell * s->m2);

    for (il = 0; il < s->nell; il++)
    {
      ell = il + s->lmin;

      const gdouble v = pow (ell / lp, p[NC_DATA_PLANCK_SMICA_PAR_CIB_INDEX] - (-1.3));
      gdouble *pl     = s->perl + il * s->m2;

      for (i = 0; i < mT; i++)
        for (j = 0; j < mT; j++)
          pl[i * m + j] = v * gt[ell * GCIB_NF * GCIB_NF + i * GCIB_NF + j] * A[i * mT + j];
    }

    smica_bin_add (s);
  }

  /* --- gibXsz: CIB x tSZ cross (temperature) --- */
  {
    const gdouble *gxt = ncm_vector_data (s->tmpl_gibxsz);
    gdouble A[16];

    for (i = 0; i < mT; i++)
      fnu_gx[i] = sz_spectrum (ncm_vector_get (s->freqs, i)) * ncm_vector_get (s->sz_color, i);

    for (i = 0; i < mT; i++)
      for (j = 0; j < mT; j++)
        A[i * mT + j] = -p[NC_DATA_PLANCK_SMICA_PAR_XI_SZ_CIB] * sqrt (p[NC_DATA_PLANCK_SMICA_PAR_A_SZ]) *
                        (sqrt (fnu_gx[i] * p[NC_DATA_PLANCK_SMICA_PAR_A_CIB_217] * ncm_vector_get (s->gibxsz_conv, j)) +
                         sqrt (fnu_gx[j] * p[NC_DATA_PLANCK_SMICA_PAR_A_CIB_217] * ncm_vector_get (s->gibxsz_conv, i)));

    if (s->has_pol)
      memset (s->perl, 0, sizeof (gdouble) * s->nell * s->m2);

    for (il = 0; il < s->nell; il++)
    {
      ell = il + s->lmin;

      const gdouble c = gxt[ell - 2] * twopi / (ell * (ell + 1.0)); /* corr_template ell-2 */
      gdouble *pl     = s->perl + il * s->m2;

      for (i = 0; i < mT; i++)
        for (j = 0; j < mT; j++)
          pl[i * m + j] = A[i * mT + j] * c;
    }

    smica_bin_add (s);
  }

  /* --- cnoise: gal545 dust (abso, T), beam leakage (all pairs), subpixel --- */
  {
    gdouble amp[36];

    memset (amp, 0, sizeof (amp));
    amp[0 * m + 0] = p[NC_DATA_PLANCK_SMICA_PAR_GAL545_100];
    amp[1 * m + 1] = p[NC_DATA_PLANCK_SMICA_PAR_GAL545_143];
    amp[1 * m + 2] = amp[2 * m + 1] = p[NC_DATA_PLANCK_SMICA_PAR_GAL545_143_217];
    amp[2 * m + 2] = p[NC_DATA_PLANCK_SMICA_PAR_GAL545_217];
    smica_cnoise_add (s, s->tmpl_dust, amp, TRUE, 200);

    for (i = 0; i < s->m2; i++) /* beam leakage: all pairs = 1 */
      amp[i] = 1.0;

    smica_cnoise_add (s, s->tmpl_leak, amp, FALSE, 0);

    memset (amp, 0, sizeof (amp)); /* subpixel (temperature) */
    amp[0 * m + 0] = p[NC_DATA_PLANCK_SMICA_PAR_SBPX_100_100];
    amp[1 * m + 1] = p[NC_DATA_PLANCK_SMICA_PAR_SBPX_143_143];
    amp[1 * m + 2] = amp[2 * m + 1] = p[NC_DATA_PLANCK_SMICA_PAR_SBPX_143_217];
    amp[2 * m + 2] = p[NC_DATA_PLANCK_SMICA_PAR_SBPX_217_217];
    smica_cnoise_add (s, s->tmpl_sbpx, amp, FALSE, 0);
  }

  /* --- polarization foregrounds (TTTEEE only) --- */
  if (s->has_pol)
  {
    const guint mE = m - mT;
    gdouble amp[36];

    /* galf dust_EE (E x E) */
    memset (amp, 0, sizeof (amp));

    for (i = 0; i < mE; i++)
      for (j = i; j < mE; j++)
      {
        const gdouble v = pp[PP_GALF_EE_A (smica_utri (i, j, mE))];

        amp[(mT + i) * m + (mT + j)] = v;
        amp[(mT + j) * m + (mT + i)] = v;
      }

    smica_pwfe_add (s, amp, pp[PP_GALF_EE_INDEX]);

    /* galf dust_TE (T x E), param keyed by the sorted frequency pair */
    memset (amp, 0, sizeof (amp));

    for (i = 0; i < mT; i++)
      for (j = 0; j < mE; j++)
      {
        const gdouble v = pp[PP_GALF_TE_A (smica_utri (MIN (i, j), MAX (i, j), mE))];

        amp[i * m + (mT + j)] = v;
        amp[(mT + j) * m + i] = v;
      }

    smica_pwfe_add (s, amp, pp[PP_GALF_TE_INDEX]);

    /* EE end-to-end correlated noise (E diagonal) */
    memset (amp, 0, sizeof (amp));

    for (i = 0; i < mE; i++)
      amp[(mT + i) * m + (mT + i)] = pp[PP_E2E (i)];

    smica_cnoise_add (s, s->tmpl_e2e, amp, FALSE, 0);

    /* subpixel (E x E) */
    memset (amp, 0, sizeof (amp));

    for (i = 0; i < mE; i++)
      for (j = i; j < mE; j++)
      {
        const gdouble v = pp[PP_SBPX_EE (smica_utri (i, j, mE))];

        amp[(mT + i) * m + (mT + j)] = v;
        amp[(mT + j) * m + (mT + i)] = v;
      }

    smica_cnoise_add (s, s->tmpl_sbpx, amp, FALSE, 0);
  }

  /* --- multiplicative calibration --- */
  if (!s->has_pol)
  {
    /* icalTP (143 = ref) then totcal: Rq[i,j] *= cal_i*cal_j / A_planck^2 */
    const gdouble tc = 1.0 / (p[NC_DATA_PLANCK_SMICA_PAR_A_PLANCK] * p[NC_DATA_PLANCK_SMICA_PAR_A_PLANCK]);
    gdouble cal[3];

    cal[0] = 1.0 / sqrt (p[NC_DATA_PLANCK_SMICA_PAR_CALIB_100T]);
    cal[1] = 1.0;
    cal[2] = 1.0 / sqrt (p[NC_DATA_PLANCK_SMICA_PAR_CALIB_217T]);

    for (b = 0; b < s->nbins; b++)
      for (i = 0; i < m; i++)
        for (j = 0; j < m; j++)
          s->Rq[b * s->m2 + i * m + j] *= cal[i] * cal[j] * tc;
  }
  else
  {
    /* icalTP two-term mixing, then totcalP (A_pol on E legs), then totcal. */
    const gdouble *w   = ncm_vector_data (s->ical_w);
    const guint32 *oth = (const guint32 *) s->ical_other->data;
    const gdouble apol = pp[PP_A_POL];
    const gdouble tc   = 1.0 / (p[NC_DATA_PLANCK_SMICA_PAR_A_PLANCK] * p[NC_DATA_PLANCK_SMICA_PAR_A_PLANCK]);
    gdouble calvec[16];
    gdouble calval[5];

    g_assert_cmpuint (s->ical_im->len, ==, 5);
    calval[0] = p[NC_DATA_PLANCK_SMICA_PAR_CALIB_100T];
    calval[1] = p[NC_DATA_PLANCK_SMICA_PAR_CALIB_217T];
    calval[2] = pp[PP_CALIB_P (0)];
    calval[3] = pp[PP_CALIB_P (1)];
    calval[4] = pp[PP_CALIB_P (2)];

    for (i = 0; i < m; i++)
      calvec[i] = 1.0;

    for (k = 0; k < s->ical_im->len; k++)
      calvec[g_array_index (s->ical_im, guint32, k)] = 1.0 / sqrt (calval[k]);

    for (b = 0; b < s->nbins; b++)
    {
      gdouble *Rqb = s->Rq + b * s->m2;

      for (i = 0; i < m; i++)
        for (j = i; j < m; j++)
        {
          const guint mpos = (i * m + j) * 2;
          const gdouble w0 = w[mpos];
          const gdouble w1 = w[mpos + 1];

          if (w0 + w1 != 0.0)
          {
            Rqb[i * m + j] *= w0 * calvec[i] * calvec[j] +
                              w1 * calvec[oth[mpos]] * calvec[oth[mpos + 1]];
            Rqb[j * m + i] = Rqb[i * m + j];
          }
        }

      for (i = 0; i < m; i++)
        for (j = 0; j < m; j++)
        {
          gdouble f = tc;

          if (fld[i] == 1)
            f /= apol;

          if (fld[j] == 1)
            f /= apol;

          Rqb[i * m + j] *= f;
        }
    }
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
 * nc_data_planck_smica_set_hipert_boltzmann() and the #NcmDataGauss API).
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

