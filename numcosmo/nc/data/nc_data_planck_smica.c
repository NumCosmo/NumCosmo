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
#include "ncm/core/ncm_serialize.h"

#include <math.h>

/* SMICA nuisance parameters, in the order NcPlanckFICorTT exposes them (== the
 * clik check_param tail order). Indices into the NcPlanckFI model. */
typedef enum
{
  P_A_cib_217 = 0, P_cib_index, P_xi_sz_cib, P_A_sz,
  P_ps_100_100, P_ps_143_143, P_ps_143_217, P_ps_217_217,
  P_ksz_norm, P_gal545_100, P_gal545_143, P_gal545_143_217, P_gal545_217,
  P_sbpx_100_100, P_sbpx_143_143, P_sbpx_143_217, P_sbpx_217_217,
  P_calib_100T, P_calib_217T, P_A_planck, P_LEN,
} SmicaPar;

struct _NcDataPlanckSmica
{
  /*< private >*/
  NcmDataGaussCov parent_instance;
  NcHIPertBoltzmann *pb;
  guint lmin, lmax, m, nbins;
  NcmVector *freqs;      /* m: [100,143,217] */
  NcmVector *a_cmb;      /* m: CMB mixing */
  NcmVector *sz_color;   /* m: tSZ colour corrections */
  NcmVector *gcib_conv;  /* m: CIB muK->MJ/sr conversion */
  NcmVector *gibxsz_conv;/* m: gibXsz conversion (100 zeroed) */
  NcmVector *bin_lmin;   /* nbins: per-bin ell offset (from lmin) */
  NcmVector *bin_lmax;   /* nbins */
  NcmVector *bin_weight; /* flattened per-(bin,ell) weights */
  NcmVector *quad_idx;   /* np: flat indices into (nbins*m*m) R_q to extract */
  /* per-ell templates (converter stores sz/ksz PRE-NORMALIZED at ell=3000) */
  NcmVector *tmpl_gcib;  /* (10001, m0, m0) 4x4 layout, offset 0 */
  NcmVector *tmpl_sz;    /* ell 2.. */
  NcmVector *tmpl_ksz;   /* ell 0.. */
  NcmVector *tmpl_gibxsz;/* ell 0.. */
  NcmVector *tmpl_dust;  /* (3001, 12, 12) */
  NcmVector *tmpl_leak;  /* (3001, 12, 12) */
  NcmVector *tmpl_sbpx;  /* (3001, 12, 12) */
  /* runtime caches (not serialized) */
  NcmVector *cl_TT;
  gdouble *Rq;           /* nbins*m*m */
  gdouble *perl;         /* nell*m*m scratch */
  guint calib_pi;        /* index of A_cib_217 in the NcPlanckFI model */
  gboolean calib_pi_set;
};

enum
{
  PROP_0, PROP_PB, PROP_LMIN, PROP_LMAX, PROP_M, PROP_NBINS,
  PROP_FREQS, PROP_A_CMB, PROP_SZ_COLOR, PROP_GCIB_CONV, PROP_GIBXSZ_CONV,
  PROP_BIN_LMIN, PROP_BIN_LMAX, PROP_BIN_WEIGHT, PROP_QUAD_IDX,
  PROP_TMPL_GCIB, PROP_TMPL_SZ, PROP_TMPL_KSZ, PROP_TMPL_GIBXSZ,
  PROP_TMPL_DUST, PROP_TMPL_LEAK, PROP_TMPL_SBPX, PROP_SIZE,
};

G_DEFINE_TYPE (NcDataPlanckSmica, nc_data_planck_smica, NCM_TYPE_DATA_GAUSS_COV)

#define GCIB_NF 4     /* gcib/cnoise template frequency dimension */
#define CN_NF 12      /* cnoise 12x12 (4 freq x 3 TEB), T uses first 3 */

static void
nc_data_planck_smica_init (NcDataPlanckSmica *smica)
{
  smica->pb = NULL;
  smica->lmin = smica->lmax = smica->m = smica->nbins = 0;
  smica->freqs = smica->a_cmb = smica->sz_color = NULL;
  smica->gcib_conv = smica->gibxsz_conv = NULL;
  smica->bin_lmin = smica->bin_lmax = smica->bin_weight = smica->quad_idx = NULL;
  smica->tmpl_gcib = smica->tmpl_sz = smica->tmpl_ksz = smica->tmpl_gibxsz = NULL;
  smica->tmpl_dust = smica->tmpl_leak = smica->tmpl_sbpx = NULL;
  smica->cl_TT = NULL;
  smica->Rq = NULL;
  smica->perl = NULL;
  smica->calib_pi = 0;
  smica->calib_pi_set = FALSE;
}

static void
nc_data_planck_smica_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataPlanckSmica *s = NC_DATA_PLANCK_SMICA (object);

  switch (prop_id)
  {
    case PROP_PB: nc_data_planck_smica_set_hipert_boltzmann (s, g_value_get_object (value)); break;
    case PROP_LMIN: s->lmin = g_value_get_uint (value); break;
    case PROP_LMAX: s->lmax = g_value_get_uint (value); break;
    case PROP_M: s->m = g_value_get_uint (value); break;
    case PROP_NBINS: s->nbins = g_value_get_uint (value); break;
    case PROP_FREQS: ncm_vector_substitute (&s->freqs, g_value_get_object (value), TRUE); break;
    case PROP_A_CMB: ncm_vector_substitute (&s->a_cmb, g_value_get_object (value), TRUE); break;
    case PROP_SZ_COLOR: ncm_vector_substitute (&s->sz_color, g_value_get_object (value), TRUE); break;
    case PROP_GCIB_CONV: ncm_vector_substitute (&s->gcib_conv, g_value_get_object (value), TRUE); break;
    case PROP_GIBXSZ_CONV: ncm_vector_substitute (&s->gibxsz_conv, g_value_get_object (value), TRUE); break;
    case PROP_BIN_LMIN: ncm_vector_substitute (&s->bin_lmin, g_value_get_object (value), TRUE); break;
    case PROP_BIN_LMAX: ncm_vector_substitute (&s->bin_lmax, g_value_get_object (value), TRUE); break;
    case PROP_BIN_WEIGHT: ncm_vector_substitute (&s->bin_weight, g_value_get_object (value), TRUE); break;
    case PROP_QUAD_IDX: ncm_vector_substitute (&s->quad_idx, g_value_get_object (value), TRUE); break;
    case PROP_TMPL_GCIB: ncm_vector_substitute (&s->tmpl_gcib, g_value_get_object (value), TRUE); break;
    case PROP_TMPL_SZ: ncm_vector_substitute (&s->tmpl_sz, g_value_get_object (value), TRUE); break;
    case PROP_TMPL_KSZ: ncm_vector_substitute (&s->tmpl_ksz, g_value_get_object (value), TRUE); break;
    case PROP_TMPL_GIBXSZ: ncm_vector_substitute (&s->tmpl_gibxsz, g_value_get_object (value), TRUE); break;
    case PROP_TMPL_DUST: ncm_vector_substitute (&s->tmpl_dust, g_value_get_object (value), TRUE); break;
    case PROP_TMPL_LEAK: ncm_vector_substitute (&s->tmpl_leak, g_value_get_object (value), TRUE); break;
    case PROP_TMPL_SBPX: ncm_vector_substitute (&s->tmpl_sbpx, g_value_get_object (value), TRUE); break;
    default: G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); break;
  }
}

static void
nc_data_planck_smica_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataPlanckSmica *s = NC_DATA_PLANCK_SMICA (object);

  switch (prop_id)
  {
    case PROP_PB: g_value_set_object (value, s->pb); break;
    case PROP_LMIN: g_value_set_uint (value, s->lmin); break;
    case PROP_LMAX: g_value_set_uint (value, s->lmax); break;
    case PROP_M: g_value_set_uint (value, s->m); break;
    case PROP_NBINS: g_value_set_uint (value, s->nbins); break;
    case PROP_FREQS: g_value_set_object (value, s->freqs); break;
    case PROP_A_CMB: g_value_set_object (value, s->a_cmb); break;
    case PROP_SZ_COLOR: g_value_set_object (value, s->sz_color); break;
    case PROP_GCIB_CONV: g_value_set_object (value, s->gcib_conv); break;
    case PROP_GIBXSZ_CONV: g_value_set_object (value, s->gibxsz_conv); break;
    case PROP_BIN_LMIN: g_value_set_object (value, s->bin_lmin); break;
    case PROP_BIN_LMAX: g_value_set_object (value, s->bin_lmax); break;
    case PROP_BIN_WEIGHT: g_value_set_object (value, s->bin_weight); break;
    case PROP_QUAD_IDX: g_value_set_object (value, s->quad_idx); break;
    case PROP_TMPL_GCIB: g_value_set_object (value, s->tmpl_gcib); break;
    case PROP_TMPL_SZ: g_value_set_object (value, s->tmpl_sz); break;
    case PROP_TMPL_KSZ: g_value_set_object (value, s->tmpl_ksz); break;
    case PROP_TMPL_GIBXSZ: g_value_set_object (value, s->tmpl_gibxsz); break;
    case PROP_TMPL_DUST: g_value_set_object (value, s->tmpl_dust); break;
    case PROP_TMPL_LEAK: g_value_set_object (value, s->tmpl_leak); break;
    case PROP_TMPL_SBPX: g_value_set_object (value, s->tmpl_sbpx); break;
    default: G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); break;
  }
}

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
  ncm_vector_clear (&s->bin_lmin);
  ncm_vector_clear (&s->bin_lmax);
  ncm_vector_clear (&s->bin_weight);
  ncm_vector_clear (&s->quad_idx);
  ncm_vector_clear (&s->tmpl_gcib);
  ncm_vector_clear (&s->tmpl_sz);
  ncm_vector_clear (&s->tmpl_ksz);
  ncm_vector_clear (&s->tmpl_gibxsz);
  ncm_vector_clear (&s->tmpl_dust);
  ncm_vector_clear (&s->tmpl_leak);
  ncm_vector_clear (&s->tmpl_sbpx);
  ncm_vector_clear (&s->cl_TT);

  G_OBJECT_CLASS (nc_data_planck_smica_parent_class)->dispose (object);
}

static void
nc_data_planck_smica_finalize (GObject *object)
{
  NcDataPlanckSmica *s = NC_DATA_PLANCK_SMICA (object);

  g_clear_pointer (&s->Rq, g_free);
  g_clear_pointer (&s->perl, g_free);

  G_OBJECT_CLASS (nc_data_planck_smica_parent_class)->finalize (object);
}

static void _nc_data_planck_smica_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_planck_smica_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp);

static void
install_vec (GObjectClass *oc, guint id, const gchar *name)
{
  g_object_class_install_property (oc, id,
                                   g_param_spec_object (name, NULL, name, NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

static void
nc_data_planck_smica_class_init (NcDataPlanckSmicaClass *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class          = NCM_DATA_CLASS (klass);
  NcmDataGaussCovClass *gauss_class = NCM_DATA_GAUSS_COV_CLASS (klass);

  object_class->set_property = &nc_data_planck_smica_set_property;
  object_class->get_property = &nc_data_planck_smica_get_property;
  object_class->dispose      = &nc_data_planck_smica_dispose;
  object_class->finalize     = &nc_data_planck_smica_finalize;

  g_object_class_install_property (object_class, PROP_PB,
                                   g_param_spec_object ("hipert-boltzmann", NULL, "Cls source",
                                                        NC_TYPE_HIPERT_BOLTZMANN,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
#define UINT_PROP(id, name, def) g_object_class_install_property (object_class, id, \
    g_param_spec_uint (name, NULL, name, 0, G_MAXUINT, def, \
                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB))
  UINT_PROP (PROP_LMIN, "lmin", 30);
  UINT_PROP (PROP_LMAX, "lmax", 2508);
  UINT_PROP (PROP_M, "m-channels", 3);
  UINT_PROP (PROP_NBINS, "nbins", 215);
#undef UINT_PROP

  install_vec (object_class, PROP_FREQS, "freqs");
  install_vec (object_class, PROP_A_CMB, "a-cmb");
  install_vec (object_class, PROP_SZ_COLOR, "sz-color");
  install_vec (object_class, PROP_GCIB_CONV, "gcib-conv");
  install_vec (object_class, PROP_GIBXSZ_CONV, "gibxsz-conv");
  install_vec (object_class, PROP_BIN_LMIN, "bin-lmin");
  install_vec (object_class, PROP_BIN_LMAX, "bin-lmax");
  install_vec (object_class, PROP_BIN_WEIGHT, "bin-weight");
  install_vec (object_class, PROP_QUAD_IDX, "quad-idx");
  install_vec (object_class, PROP_TMPL_GCIB, "tmpl-gcib");
  install_vec (object_class, PROP_TMPL_SZ, "tmpl-sz");
  install_vec (object_class, PROP_TMPL_KSZ, "tmpl-ksz");
  install_vec (object_class, PROP_TMPL_GIBXSZ, "tmpl-gibxsz");
  install_vec (object_class, PROP_TMPL_DUST, "tmpl-dust");
  install_vec (object_class, PROP_TMPL_LEAK, "tmpl-leak");
  install_vec (object_class, PROP_TMPL_SBPX, "tmpl-sbpx");

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
  const guint m2 = s->m * s->m;
  guint b, il, i, bb = 0;

  for (b = 0; b < s->nbins; b++)
  {
    const guint lo = (guint) ncm_vector_get (s->bin_lmin, b);
    const guint hi = (guint) ncm_vector_get (s->bin_lmax, b);
    gdouble *Rqb   = s->Rq + b * m2;

    for (il = lo; il <= hi; il++)
    {
      const gdouble w = ncm_vector_get (s->bin_weight, bb++);
      const gdouble *pl = s->perl + il * m2;

      for (i = 0; i < m2; i++)
        Rqb[i] += w * pl[i];
    }
  }
}

static gdouble
sz_spectrum (gdouble nu)
{
  const gdouble nu0 = 143.0;
  const gdouble x = nu / 56.78, x0 = nu0 / 56.78;

  return (2.0 - x / 2.0 / tanh (x / 2.0)) / (2.0 - x0 / 2.0 / tanh (x0 / 2.0));
}

static void
_nc_data_planck_smica_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp)
{
  NcDataPlanckSmica *s = NC_DATA_PLANCK_SMICA (gauss);
  NcPlanckFI *pfi      = NC_PLANCK_FI (ncm_mset_peek (mset, nc_planck_fi_id ()));
  const guint m = s->m, m2 = s->m * s->m, nell = s->lmax - s->lmin + 1;
  const guint np = ncm_data_gauss_cov_get_size (gauss);
  const gdouble twopi = 2.0 * M_PI;
  gdouble p[P_LEN];
  gdouble fnu[8], fnu_gx[8];
  guint b, il, i, j, k, ell, l;

  if (pfi == NULL)
    g_error ("_nc_data_planck_smica_mean_func: need a NcPlanckFI model in mset.");

  /* nuisance params (NcPlanckFICorTT order == check_param tail order) */
  for (i = 0; i < P_LEN; i++)
    p[i] = ncm_model_param_get (NCM_MODEL (pfi), i);

  if (s->cl_TT == NULL)
    s->cl_TT = ncm_vector_new (s->lmax + 1);
  if (s->Rq == NULL)
    s->Rq = g_new (gdouble, s->nbins * m2);
  if (s->perl == NULL)
    s->perl = g_new (gdouble, nell * m2);

  nc_hipert_boltzmann_get_TT_Cls (s->pb, s->cl_TT);

  memset (s->Rq, 0, sizeof (gdouble) * s->nbins * m2);

  /* --- CMB: P_q = binned raw Cl; Rq += P_q * a_cmb[i]*a_cmb[j] --- */
  for (il = 0; il < nell; il++)
  {
    const gdouble cl = ncm_vector_get (s->cl_TT, il + s->lmin);
    gdouble *pl = s->perl + il * m2;

    for (i = 0; i < m; i++)
      for (j = 0; j < m; j++)
        pl[i * m + j] = cl * ncm_vector_get (s->a_cmb, i) * ncm_vector_get (s->a_cmb, j);
  }
  smica_bin_add (s);

  /* --- sz: A_sz * 2pi/(l(l+1)) * tmpl[l-2] * fnu_i*fnu_j --- */
  for (i = 0; i < m; i++)
    fnu[i] = sz_spectrum (ncm_vector_get (s->freqs, i)) * ncm_vector_get (s->sz_color, i);
  for (il = 0; il < nell; il++)
  {
    ell = il + s->lmin;
    const gdouble c = p[P_A_sz] * twopi / (ell * (ell + 1.0)) *
                      ncm_vector_get (s->tmpl_sz, ell - 2);
    gdouble *pl = s->perl + il * m2;
    for (i = 0; i < m; i++)
      for (j = 0; j < m; j++)
        pl[i * m + j] = c * fnu[i] * fnu[j];
  }
  smica_bin_add (s);

  /* --- ksz: ksz_norm * 2pi/(l(l+1)) * tmpl[l] (freq-independent) --- */
  for (il = 0; il < nell; il++)
  {
    ell = il + s->lmin;
    const gdouble c = p[P_ksz_norm] * twopi / (ell * (ell + 1.0)) *
                      ncm_vector_get (s->tmpl_ksz, ell);
    gdouble *pl = s->perl + il * m2;
    for (i = 0; i < m2; i++)
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
    ps[0 * m + 0] = p[P_ps_100_100] / nrm;
    ps[1 * m + 1] = p[P_ps_143_143] / nrm;
    ps[1 * m + 2] = ps[2 * m + 1] = p[P_ps_143_217] / nrm;
    ps[2 * m + 2] = p[P_ps_217_217] / nrm;
    for (b = 0; b < s->nbins; b++)
      for (i = 0; i < m2; i++)
        s->Rq[b * m2 + i] += ps[i];
  }

  /* --- gcib: rigid CIB (template GCIB_NF x GCIB_NF; irigid=217=index 2) --- */
  {
    const gdouble lp = 3000.0;
    const gdouble *gt = ncm_vector_data (s->tmpl_gcib);
    const gdouble conv2 = ncm_vector_get (s->gcib_conv, 2);
    const gdouble tref = gt[(guint) lp * GCIB_NF * GCIB_NF + 2 * GCIB_NF + 2];
    gdouble A[9];
    for (i = 0; i < m; i++)
      for (j = 0; j < m; j++)
      {
        A[i * m + j] = (p[P_A_cib_217] / tref *
                        ncm_vector_get (s->gcib_conv, i) * ncm_vector_get (s->gcib_conv, j) /
                        (conv2 * conv2)) / lp / (lp + 1.0) * twopi;
      }
    for (il = 0; il < nell; il++)
    {
      ell = il + s->lmin;
      const gdouble v = pow (ell / lp, p[P_cib_index] - (-1.3));
      gdouble *pl = s->perl + il * m2;
      for (i = 0; i < m; i++)
        for (j = 0; j < m; j++)
          pl[i * m + j] = v * gt[ell * GCIB_NF * GCIB_NF + i * GCIB_NF + j] * A[i * m + j];
    }
    smica_bin_add (s);
  }

  /* --- gibXsz: CIB x tSZ cross --- */
  {
    const gdouble *gxt = ncm_vector_data (s->tmpl_gibxsz);
    gdouble A[9];
    for (i = 0; i < m; i++)
      fnu_gx[i] = sz_spectrum (ncm_vector_get (s->freqs, i)) * ncm_vector_get (s->sz_color, i);
    for (i = 0; i < m; i++)
      for (j = 0; j < m; j++)
        A[i * m + j] = -p[P_xi_sz_cib] * sqrt (p[P_A_sz]) *
                       (sqrt (fnu_gx[i] * p[P_A_cib_217] * ncm_vector_get (s->gibxsz_conv, j)) +
                        sqrt (fnu_gx[j] * p[P_A_cib_217] * ncm_vector_get (s->gibxsz_conv, i)));
    for (il = 0; il < nell; il++)
    {
      ell = il + s->lmin;
      const gdouble c = gxt[ell - 2] * twopi / (ell * (ell + 1.0)); /* corr_template ell-2 */
      gdouble *pl = s->perl + il * m2;
      for (i = 0; i < m; i++)
        for (j = 0; j < m; j++)
          pl[i * m + j] = A[i * m + j] * c;
    }
    smica_bin_add (s);
  }

  /* --- cnoise (dust abso=1 lpiv=200; leak/sbpx abso=0). template CN_NF^2 --- */
#define CNOISE(tmpl, get_amp, abso, lpiv)                                            \
  do {                                                                              \
    const gdouble *ct = ncm_vector_data (tmpl);                                     \
    for (il = 0; il < nell; il++) {                                                 \
      ell = il + s->lmin;                                                           \
      gdouble *pl = s->perl + il * m2;                                              \
      for (i = 0; i < m; i++)                                                       \
        for (j = 0; j < m; j++) {                                                   \
          const gdouble v = get_amp;                                                \
          gdouble a;                                                                \
          if (abso)                                                                 \
            a = v / (ct[(lpiv) * CN_NF * CN_NF + i * CN_NF + j] *                    \
                     (lpiv) * ((lpiv) + 1.0) / twopi);                              \
          else                                                                      \
            a = v;                                                                  \
          pl[i * m + j] = ct[ell * CN_NF * CN_NF + i * CN_NF + j] * a;              \
        }                                                                           \
    }                                                                               \
    smica_bin_add (s);                                                              \
  } while (0)

  /* dust (gal545): pairs 00,11,12,22; others 0 */
  {
    gdouble dg[9];
    memset (dg, 0, sizeof (dg));
    dg[0 * m + 0] = p[P_gal545_100];
    dg[1 * m + 1] = p[P_gal545_143];
    dg[1 * m + 2] = dg[2 * m + 1] = p[P_gal545_143_217];
    dg[2 * m + 2] = p[P_gal545_217];
    CNOISE (s->tmpl_dust, dg[i * m + j], 1, 200);
  }
  /* beam leakage: all pairs = 1, abso=0 */
  CNOISE (s->tmpl_leak, 1.0, 0, 0);
  /* subpixel: pairs 00,11,12,22 = A_sbpx; abso=0 */
  {
    gdouble sg[9];
    memset (sg, 0, sizeof (sg));
    sg[0 * m + 0] = p[P_sbpx_100_100];
    sg[1 * m + 1] = p[P_sbpx_143_143];
    sg[1 * m + 2] = sg[2 * m + 1] = p[P_sbpx_143_217];
    sg[2 * m + 2] = p[P_sbpx_217_217];
    CNOISE (s->tmpl_sbpx, sg[i * m + j], 0, 0);
  }
#undef CNOISE

  /* --- icalTP (mult): Rq[i,j] *= 1/sqrt(cal_i)/sqrt(cal_j); 143 = ref --- */
  {
    gdouble cal[3];
    cal[0] = 1.0 / sqrt (p[P_calib_100T]);
    cal[1] = 1.0;
    cal[2] = 1.0 / sqrt (p[P_calib_217T]);
    /* --- totcal (mult): Rq *= 1/A_planck^2 --- */
    const gdouble tc = 1.0 / (p[P_A_planck] * p[P_A_planck]);
    for (b = 0; b < s->nbins; b++)
      for (i = 0; i < m; i++)
        for (j = 0; j < m; j++)
          s->Rq[b * m2 + i * m + j] *= cal[i] * cal[j] * tc;
  }

  /* --- extract masked entries via quad_idx --- */
  for (k = 0; k < np; k++)
    ncm_vector_set (vp, k, s->Rq[(guint) ncm_vector_get (s->quad_idx, k)]);

  (void) l;
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
