/***************************************************************************
 *            ncm_stats_dist1d_epdf.c
 *
 *  Sat March 14 19:32:06 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_dist1d_epdf.c
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmStatsDist1dEPDF:
 *
 * One dimensional probability distribution based on an EPDF.
 *
 * Reconstruction of an arbitrary one dimensional probability distribution based on a
 * Empirical Probability Distribution Function (EPDF).
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_spline_func.h"
#include "math/ncm_stats_dist1d_epdf.h"
#include "math/ncm_c.h"
#include "math/ncm_cfg.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <complex.h>
#ifdef HAVE_FFTW3
#include <fftw3.h>
#endif /* HAVE_FFTW3 */
#include <gsl/gsl_sort.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_MAX_OBS,
  PROP_NOBS,
  PROP_BANDWIDTH,
  PROP_H_FIXED,
  PROP_SD_MIN_SCALE,
  PROP_OUTLIERS_THRESHOLD,
};


struct _NcmStatsDist1dEPDF
{
  /*< private >*/
  NcmStatsDist1d parent_instance;
  NcmStatsVec *obs_stats;
  guint max_obs;
  NcmStatsDist1dEPDFBw bw;
  gdouble h_fixed;
  gdouble sd_min_scale;
  gdouble outliers_threshold;
  gdouble h;
  guint n_obs;
  guint np_obs;
  gdouble WT;
  GArray *obs;
  gdouble min;
  gdouble max;
  gboolean list_sorted;
  guint fftsize;
  NcmVector *Iv;
  NcmVector *p_data;
  NcmVector *p_tilde;
  NcmVector *p_tilde2;
  NcmVector *p_est;
  NcmVector *xv;
  NcmVector *pv;
  gpointer fft_data_to_tilde;
  gpointer fft_tilde_to_est;
  NcmSpline *ph_spline;
  NcmSpline *p_spline;
  gboolean bw_set;
};

G_DEFINE_TYPE (NcmStatsDist1dEPDF, ncm_stats_dist1d_epdf, NCM_TYPE_STATS_DIST1D)

typedef struct _NcmStatsDist1dEPDFObs
{
  gdouble x;
  gdouble w;
} NcmStatsDist1dEPDFObs;

static void
ncm_stats_dist1d_epdf_init (NcmStatsDist1dEPDF *epdf1d)
{
  epdf1d->obs_stats          = ncm_stats_vec_new (1, NCM_STATS_VEC_VAR, FALSE);
  epdf1d->max_obs            = 0;
  epdf1d->bw                 = NCM_STATS_DIST1D_EPDF_BW_LEN;
  epdf1d->h_fixed            = 0.0;
  epdf1d->sd_min_scale       = 0.0;
  epdf1d->outliers_threshold = 0.0;
  epdf1d->h                  = 0.0;
  epdf1d->n_obs              = 0;
  epdf1d->np_obs             = 0;
  epdf1d->WT                 = 0.0;
  epdf1d->obs                = NULL;
  epdf1d->min                = GSL_POSINF;
  epdf1d->max                = GSL_NEGINF;

  epdf1d->fftsize           = 0;
  epdf1d->Iv                = NULL;
  epdf1d->p_data            = NULL;
  epdf1d->p_est             = NULL;
  epdf1d->p_tilde           = NULL;
  epdf1d->p_tilde2          = NULL;
  epdf1d->xv                = NULL;
  epdf1d->pv                = NULL;
  epdf1d->fft_data_to_tilde = NULL;
  epdf1d->fft_tilde_to_est  = NULL;

  epdf1d->ph_spline = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  epdf1d->p_spline  = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  epdf1d->bw_set    = FALSE;

  ncm_stats_vec_enable_quantile (epdf1d->obs_stats, 0.5);
}

static void
ncm_stats_dist1d_epdf_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_stats_dist1d_epdf_parent_class)->constructed (object);
  {
    /*NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (object);*/
    NcmStatsDist1dEPDF *epdf1d = NCM_STATS_DIST1D_EPDF (object);

    epdf1d->obs = g_array_sized_new (FALSE, FALSE, sizeof (NcmStatsDist1dEPDFObs), epdf1d->max_obs);
  }
}

static void
ncm_stats_dist1d_epdf_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmStatsDist1dEPDF *epdf1d = NCM_STATS_DIST1D_EPDF (object);

  g_return_if_fail (NCM_IS_STATS_DIST1D_EPDF (object));

  switch (prop_id)
  {
    case PROP_MAX_OBS:
      epdf1d->max_obs = g_value_get_uint (value);
      break;
    case PROP_BANDWIDTH:
      ncm_stats_dist1d_epdf_set_bw_type (epdf1d, g_value_get_enum (value));
      break;
    case PROP_H_FIXED:
      epdf1d->h_fixed = g_value_get_double (value);
      break;
    case PROP_SD_MIN_SCALE:
      epdf1d->sd_min_scale = g_value_get_double (value);
      break;
    case PROP_OUTLIERS_THRESHOLD:
      epdf1d->outliers_threshold = g_value_get_double (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_stats_dist1d_epdf_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmStatsDist1dEPDF *epdf1d = NCM_STATS_DIST1D_EPDF (object);

  g_return_if_fail (NCM_IS_STATS_DIST1D_EPDF (object));

  switch (prop_id)
  {
    case PROP_MAX_OBS:
      g_value_set_uint (value, epdf1d->max_obs);
      break;
    case PROP_NOBS:
      g_value_set_uint (value, epdf1d->n_obs);
      break;
    case PROP_BANDWIDTH:
      g_value_set_enum (value, ncm_stats_dist1d_epdf_get_bw_type (epdf1d));
      break;
    case PROP_H_FIXED:
      g_value_set_double (value, epdf1d->h_fixed);
      break;
    case PROP_SD_MIN_SCALE:
      g_value_set_double (value, epdf1d->sd_min_scale);
      break;
    case PROP_OUTLIERS_THRESHOLD:
      g_value_set_double (value, epdf1d->outliers_threshold);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_stats_dist1d_epdf_dispose (GObject *object)
{
  NcmStatsDist1dEPDF *epdf1d = NCM_STATS_DIST1D_EPDF (object);

  ncm_stats_vec_clear (&epdf1d->obs_stats);
  g_clear_pointer (&epdf1d->obs, g_array_unref);

  ncm_vector_clear (&epdf1d->Iv);
  ncm_vector_clear (&epdf1d->p_data);
  ncm_vector_clear (&epdf1d->p_tilde);
  ncm_vector_clear (&epdf1d->p_tilde2);
  ncm_vector_clear (&epdf1d->p_est);
  ncm_vector_clear (&epdf1d->xv);
  ncm_vector_clear (&epdf1d->pv);
  ncm_spline_clear (&epdf1d->ph_spline);
  ncm_spline_clear (&epdf1d->p_spline);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist1d_epdf_parent_class)->dispose (object);
}

static void
ncm_stats_dist1d_epdf_finalize (GObject *object)
{
  NcmStatsDist1dEPDF *epdf1d = NCM_STATS_DIST1D_EPDF (object);

#ifdef HAVE_FFTW3
  g_clear_pointer (&epdf1d->fft_data_to_tilde, fftw_destroy_plan);
  g_clear_pointer (&epdf1d->fft_tilde_to_est, fftw_destroy_plan);
#endif /* HAVE_FFTW3 */

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist1d_epdf_parent_class)->finalize (object);
}

static gdouble _ncm_stats_dist1d_epdf_p (NcmStatsDist1d *sd1, gdouble x);
static gdouble _ncm_stats_dist1d_epdf_m2lnp (NcmStatsDist1d *sd1, gdouble x);
static void _ncm_stats_dist1d_epdf_prepare (NcmStatsDist1d *sd1);
static gdouble _ncm_stats_dist1d_epdf_get_current_h (NcmStatsDist1d *sd1);

static void
ncm_stats_dist1d_epdf_class_init (NcmStatsDist1dEPDFClass *klass)
{
  GObjectClass *object_class     = G_OBJECT_CLASS (klass);
  NcmStatsDist1dClass *sd1_class = NCM_STATS_DIST1D_CLASS (klass);

  object_class->constructed  = ncm_stats_dist1d_epdf_constructed;
  object_class->set_property = ncm_stats_dist1d_epdf_set_property;
  object_class->get_property = ncm_stats_dist1d_epdf_get_property;
  object_class->dispose      = ncm_stats_dist1d_epdf_dispose;
  object_class->finalize     = ncm_stats_dist1d_epdf_finalize;

  g_object_class_install_property (object_class,
                                   PROP_MAX_OBS,
                                   g_param_spec_uint ("max-obs",
                                                      NULL,
                                                      "Maximum observations before compacting",
                                                      10, G_MAXUINT, 100000,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NOBS,
                                   g_param_spec_uint ("n-obs",
                                                      NULL,
                                                      "Number of observations",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_BANDWIDTH,
                                   g_param_spec_enum ("bandwidth",
                                                      NULL,
                                                      "Bandwidth method",
                                                      NCM_TYPE_STATS_DIST1D_EPDF_BW, NCM_STATS_DIST1D_EPDF_BW_AUTO,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_H_FIXED,
                                   g_param_spec_double ("h-fixed",
                                                        NULL,
                                                        "Fixed bandwidth",
                                                        1.0e-5, 1.0e5, 0.1,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SD_MIN_SCALE,
                                   g_param_spec_double ("sd-min-scale",
                                                        NULL,
                                                        "Percentage of the standard deviation to use as minimum distance",
                                                        1.0e-20, 1.0e20, 1.0e-3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_OUTLIERS_THRESHOLD,
                                   g_param_spec_double ("outliers-threshold",
                                                        NULL,
                                                        "How many sigmas to consider an outlier",
                                                        1.0, 1000.0, 20.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  sd1_class->p             = &_ncm_stats_dist1d_epdf_p;
  sd1_class->m2lnp         = &_ncm_stats_dist1d_epdf_m2lnp;
  sd1_class->prepare       = &_ncm_stats_dist1d_epdf_prepare;
  sd1_class->get_current_h = &_ncm_stats_dist1d_epdf_get_current_h;
}

#define _NCM_STATS_DIST1D_HROT(sd, R, n) (pow (4.0 / 3.0, 1.0 / 5.0) * GSL_MIN ((sd), ((R) / 1.34)) * pow (n * 1.0, -1.0 / 5.0))

static gint
_ncm_stats_dist1d_epdf_cmp_double (gconstpointer a,
                                   gconstpointer b)
{
#define A (*((gdouble *) a))
#define B (*((gdouble *) b))

  return (A == B) ? 0.0 : ((A < B) ? -1 : 1);
}

#undef A
#undef B

#define _NCM_STATS_DIST1D_EPDF_OBS_N(obs, epdf1d, sd) \
        0.5 * (erf (((obs)->x - (epdf1d)->min) / (sd)) + erf (((epdf1d)->max - (obs)->x) / (sd)))

static gdouble _ncm_stats_dist1d_epdf_p_gk (NcmStatsDist1dEPDF *epdf1d, gdouble x);

static void
_ncm_stats_dist1d_epdf_compact_obs (NcmStatsDist1dEPDF *epdf1d)
{
  register guint i, j;
  gint obs_len = epdf1d->obs->len;

  if (epdf1d->list_sorted)
    return;

  g_array_sort (epdf1d->obs, _ncm_stats_dist1d_epdf_cmp_double);
  j = 0;

  {
    NcmStatsDist1dEPDFObs *obs_j = &g_array_index (epdf1d->obs, NcmStatsDist1dEPDFObs, j);
    const gdouble sd_e           = ncm_stats_vec_get_sd (epdf1d->obs_stats, 0);
    const gdouble min_dist       = sd_e * epdf1d->sd_min_scale;

    /*const gdouble R_e      = ncm_stats_vec_get_quantile_spread (epdf1d->obs_stats, 0);*/
    /*const gdouble sd       = _NCM_STATS_DIST1D_SROT (sd_e, R_e, epdf1d->n_obs);*/

    for (i = 1; i < epdf1d->obs->len; i++)
    {
      NcmStatsDist1dEPDFObs *obs_i = &g_array_index (epdf1d->obs, NcmStatsDist1dEPDFObs, i);

      if (fabs (obs_j->x - obs_i->x) < min_dist)
      {
        obs_j->x = (obs_j->x * obs_j->w + obs_i->x * obs_i->w) / (obs_i->w + obs_j->w);
        obs_j->w = obs_i->w + obs_j->w;

        obs_len--;
      }
      else
      {
        j++;

        if (i != j)
          g_array_index (epdf1d->obs, NcmStatsDist1dEPDFObs, j) = *obs_i;

        obs_j = &g_array_index (epdf1d->obs, NcmStatsDist1dEPDFObs, j);
      }
    }
  }

  g_array_set_size (epdf1d->obs, obs_len);
  epdf1d->list_sorted = TRUE;
}

static guint
_ncm_stats_dist1d_epdf_bsearch (GArray *obs, const gdouble x, const guint l, const guint u)
{
  if (u > l + 1)
  {
    const guint m = (u + l) / 2;

    if (g_array_index (obs, NcmStatsDist1dEPDFObs, m).x > x)
      return _ncm_stats_dist1d_epdf_bsearch (obs, x, l, m);
    else
      return _ncm_stats_dist1d_epdf_bsearch (obs, x, m, u);
  }
  else
  {
    return l;
  }
}

static gdouble
_ncm_stats_dist1d_epdf_estimate_df2 (NcmVector *p_tilde2, NcmVector *Iv, const guint n, const guint l, const gdouble t)
{
  const gdouble pi2  = M_PI * M_PI;
  const gdouble pi2l = gsl_pow_int (pi2, l);
  gdouble s          = 0.0;
  guint i;

  for (i = 0; i < n; i++)
  {
    const gdouble Ii         = ncm_vector_fast_get (Iv, i);
    const gdouble Ili        = gsl_pow_int (Ii, l);
    const gdouble p_tilde_i2 = ncm_vector_fast_get (p_tilde2, i);

    s += Ili * p_tilde_i2 * exp (-Ii * pi2 * t);
  }

  return 0.5 * pi2l * s;
}

static gdouble
_ncm_stats_dist1d_epdf_estimate_h (NcmVector *p_tilde2, NcmVector *Iv, const guint obs_len, const guint n, const guint l, const gdouble t)
{
  const gdouble df2     = _ncm_stats_dist1d_epdf_estimate_df2 (p_tilde2, Iv, n, l, t);
  const gdouble ln_Ndf2 = log (df2 * obs_len);
  const gdouble lp05    = l + 0.5;
  const gdouble ln_fact = log1p (exp2 (-lp05)) + lp05 * M_LN2  - ncm_c_lnpi () + lgamma (lp05) - ncm_c_ln3 ();
  const gdouble tn      = exp ((ln_fact - ln_Ndf2) / (1.0 + lp05));

  g_assert (l >= 2);

  /*printf ("# l %u | % 20.15g => % 20.15g\n", l, t, tn);*/

  if (l == 2)
  {
    const gdouble df2s = _ncm_stats_dist1d_epdf_estimate_df2 (p_tilde2, Iv, n, l, tn);

    return pow (2.0 * obs_len * ncm_c_pi () * df2s, -2.0 / 5.0);
  }
  else
  {
    return _ncm_stats_dist1d_epdf_estimate_h (p_tilde2, Iv, obs_len, n, l - 1, tn);
  }
}

static void
_ncm_stats_dist1d_epdf_autobw (NcmStatsDist1dEPDF *epdf1d)
{
#ifdef HAVE_FFTW3
  const guint nbins        = exp2 (14.0 /*ceil (log2 (epdf1d->obs->len * 10))*/);
  const gdouble delta_l    = (epdf1d->max - epdf1d->min) * 2.0;
  const gdouble deltax     = delta_l / nbins;
  const gdouble xm         = (epdf1d->max + epdf1d->min) * 0.5;
  const gdouble lb         = xm - delta_l * 0.5;
  gdouble xc               = lb + deltax;
  guint fftw_default_flags = ncm_cfg_get_fftw_default_flag ();
  guint i, j;

  if (epdf1d->fftsize != nbins)
  {
    ncm_vector_clear (&epdf1d->Iv);
    epdf1d->Iv = ncm_vector_new_fftw (nbins);

    ncm_vector_clear (&epdf1d->p_data);
    epdf1d->p_data = ncm_vector_new_fftw (nbins);

    ncm_vector_clear (&epdf1d->p_tilde);
    epdf1d->p_tilde = ncm_vector_new_fftw (nbins);

    ncm_vector_clear (&epdf1d->p_tilde2);
    epdf1d->p_tilde2 = ncm_vector_new_fftw (nbins);

    ncm_vector_clear (&epdf1d->p_est);
    epdf1d->p_est = ncm_vector_new_fftw (nbins);

    ncm_vector_clear (&epdf1d->xv);
    epdf1d->xv = ncm_vector_new_fftw (nbins / 2 + 1);

    ncm_vector_clear (&epdf1d->pv);
    epdf1d->pv = ncm_vector_new_fftw (nbins / 2 + 1);

    ncm_spline_set (epdf1d->ph_spline, epdf1d->xv, epdf1d->pv, FALSE);

    epdf1d->fftsize = nbins;

    {
      G_LOCK_DEFINE_STATIC (prepare_fft_lock);
      G_LOCK (prepare_fft_lock);

      ncm_cfg_load_fftw_wisdom ("ncm_stats_dist1d_wisdown");

      ncm_cfg_lock_plan_fftw ();

      epdf1d->fft_data_to_tilde = fftw_plan_r2r_1d (nbins, ncm_vector_data (epdf1d->p_data), ncm_vector_data (epdf1d->p_tilde),
                                                    FFTW_REDFT10, fftw_default_flags | FFTW_DESTROY_INPUT);
      epdf1d->fft_tilde_to_est = fftw_plan_r2r_1d (nbins, ncm_vector_data (epdf1d->p_tilde), ncm_vector_data (epdf1d->p_est),
                                                   FFTW_REDFT01, fftw_default_flags | FFTW_DESTROY_INPUT);

      ncm_cfg_unlock_plan_fftw ();

      ncm_cfg_save_fftw_wisdom ("ncm_stats_dist1d_wisdown");

      G_UNLOCK (prepare_fft_lock);
    }

    for (i = 0; i < nbins; i++)
      ncm_vector_fast_set (epdf1d->Iv, i, gsl_pow_2 (i + 0.5));
  }

  ncm_vector_set_zero (epdf1d->p_data);

  j = 0;
  {
    const guint obs_len = epdf1d->obs->len;

    for (i = 0; i < obs_len; i++)
    {
      NcmStatsDist1dEPDFObs *obs_i = &g_array_index (epdf1d->obs, NcmStatsDist1dEPDFObs, i);

      while (obs_i->x > xc)
      {
        j++;
        xc += deltax;
      }

      ncm_vector_fast_addto (epdf1d->p_data, j, obs_i->w / epdf1d->WT);
    }
  }

  fftw_execute (epdf1d->fft_data_to_tilde);

  for (i = 0; i < nbins; i++)
  {
    const gdouble p_tilde_i  = ncm_vector_fast_get (epdf1d->p_tilde, i);
    const gdouble p_tilde_i2 = p_tilde_i * p_tilde_i;

    ncm_vector_fast_set (epdf1d->p_tilde2, i, p_tilde_i2);
  }

  {
    gdouble t  = gsl_pow_2 (epdf1d->h / delta_l);
    gdouble tn = 0.0;

    j = 0;

    while (fabs (1.0 - tn / t) > 1.0e-7)
    {
      const gdouble tni = _ncm_stats_dist1d_epdf_estimate_h (epdf1d->p_tilde2, epdf1d->Iv, epdf1d->n_obs /*obs_len*/, nbins, 7, t);

      tn = t;
      t  = tni;
      /*printf ("% 20.15g => % 20.15g\n", tn, t);*/
      j++;

      if (j >= 10000)
        g_error ("_ncm_stats_dist1d_epdf_autobw: too many steps to find bandwidth.");  /* LCOV_EXCL_LINE */
    }

    epdf1d->h = sqrt (t) * delta_l;
  }

  if (FALSE)
  {
    const gdouble kb_exp = gsl_pow_2 (M_PI  * epdf1d->h / delta_l) * 0.5;
    const guint ni       = nbins / 4;
    const guint nf       = ni * 3;

    for (i = 0; i < nbins; i++)
    {
      const gdouble Ii   = ncm_vector_fast_get (epdf1d->Iv, i);
      const gdouble kb_i = exp (-kb_exp * Ii);

      ncm_vector_fast_mulby (epdf1d->p_tilde, i, kb_i);

      if (G_UNLIKELY (kb_i == 0))
        break;
    }

    if (i < nbins)
      memset (ncm_vector_ptr (epdf1d->p_tilde, i), 0, (nbins - i) * sizeof (gdouble));

    fftw_execute (epdf1d->fft_tilde_to_est);

    j = 0;

    for (i = ni; i <= nf; i++)
    {
      const gdouble x = lb + deltax * i;

      ncm_vector_fast_set (epdf1d->xv, j, x);
      ncm_vector_fast_set (epdf1d->pv, j, log (fabs (ncm_vector_fast_get (epdf1d->p_est, i) / (2.0 * delta_l))));
      j++;
    }

    ncm_spline_prepare (epdf1d->ph_spline);
  }

#else
  g_error ("ncm_stats_dist1d_epdf_autobw: FFTW3 not available."); /* LCOV_EXCL_LINE */
#endif /* HAVE_FFTW3 */
}

static void
_ncm_stats_dist1d_epdf_set_bw (NcmStatsDist1dEPDF *epdf1d)
{
  if (epdf1d->bw_set)
  {
    return;
  }
  else
  {
    const gdouble sd_e = ncm_stats_vec_get_sd (epdf1d->obs_stats, 0);
    const gdouble R_e  = ncm_stats_vec_get_quantile_spread (epdf1d->obs_stats, 0);
    const gdouble h    = _NCM_STATS_DIST1D_HROT (sd_e, R_e, epdf1d->n_obs);

    epdf1d->h = h;

    switch (epdf1d->bw)
    {
      case NCM_STATS_DIST1D_EPDF_BW_FIXED:
        epdf1d->h = epdf1d->h_fixed;
        break;
      case NCM_STATS_DIST1D_EPDF_BW_RoT:
        break;
      case NCM_STATS_DIST1D_EPDF_BW_AUTO:
        _ncm_stats_dist1d_epdf_autobw (epdf1d);
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
        break;                   /* LCOV_EXCL_LINE */
    }

    epdf1d->bw_set = TRUE;
  }
}

static gdouble
_ncm_stats_dist1d_epdf_p_gk (NcmStatsDist1dEPDF *epdf1d, gdouble x)
{
  NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (epdf1d);
  gdouble res         = 0.0;

  if ((x < epdf1d->min) || (x > epdf1d->max))
    return 0.0;

  g_assert_cmpuint (epdf1d->obs->len, >, 0);

  _ncm_stats_dist1d_epdf_compact_obs (epdf1d);
  _ncm_stats_dist1d_epdf_set_bw (epdf1d);

  if (FALSE)
  {
    const gdouble bias_corr = 0.5 * (erf ((x - epdf1d->min) / (M_SQRT2 * epdf1d->h)) + erf ((epdf1d->max - x) / (M_SQRT2 * epdf1d->h)));
    const gdouble phat      = exp (ncm_spline_eval (epdf1d->ph_spline, x));

    return phat / bias_corr;
  }

  {
    guint s = _ncm_stats_dist1d_epdf_bsearch (epdf1d->obs, x, 0, epdf1d->obs->len - 1);
    guint i;
    gint j;

    for (i = s; i < epdf1d->obs->len; i++)
    {
      NcmStatsDist1dEPDFObs *obs = &g_array_index (epdf1d->obs, NcmStatsDist1dEPDFObs, i);
      const gdouble x_i          = obs->x;
      const gdouble de_i         = (x - x_i) / epdf1d->h;
      const gdouble de2_i        = de_i * de_i;
      const gdouble wexp_i       = obs->w * exp (-de2_i * 0.5);

      res += wexp_i;

      if (wexp_i / res < GSL_DBL_EPSILON)
        break;
    }

    for (j = s - 1; j >= 0; j--)
    {
      NcmStatsDist1dEPDFObs *obs = &g_array_index (epdf1d->obs, NcmStatsDist1dEPDFObs, j);
      const gdouble x_j          = obs->x;
      const gdouble de_j         = (x - x_j) / epdf1d->h;
      const gdouble de2_j        = de_j * de_j;
      const gdouble wexp_j       = obs->w * exp (-de2_j * 0.5);

      res += wexp_j;

      if (wexp_j / res < GSL_DBL_EPSILON)
        break;
    }
  }

  {
    const gdouble xi        = ncm_stats_dist1d_get_xi (sd1);
    const gdouble xf        = ncm_stats_dist1d_get_xf (sd1);
    const gdouble phat      = (res / (sqrt (2.0 * M_PI) * epdf1d->h) + 1.0 / (xf - xi)) / (epdf1d->WT + 1.0);
    const gdouble bias_corr = 0.5 * (erf ((x - epdf1d->min) / (M_SQRT2 * epdf1d->h)) + erf ((epdf1d->max - x) / (M_SQRT2 * epdf1d->h)));

    return phat / bias_corr;
  }
}

static gdouble
_ncm_stats_dist1d_epdf_p (NcmStatsDist1d *sd1, gdouble x)
{
  NcmStatsDist1dEPDF *epdf1d = NCM_STATS_DIST1D_EPDF (sd1);

  /*return fabs (ncm_spline_eval (epdf1d->p_spline, x));*/
  return _ncm_stats_dist1d_epdf_p_gk (epdf1d, x);
}

static gdouble
_ncm_stats_dist1d_epdf_m2lnp (NcmStatsDist1d *sd1, gdouble x)
{
  /*NcmStatsDist1dEPDF *epdf1d = NCM_STATS_DIST1D_EPDF (sd1);
   *  return -2.0 * log (fabs (ncm_spline_eval (epdf1d->p_spline, x)));*/
  return -2.0 * log (_ncm_stats_dist1d_epdf_p (sd1, x));
}

static void
_ncm_stats_dist1d_epdf_update_limits (NcmStatsDist1dEPDF *epdf1d)
{
  NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (epdf1d);

  ncm_stats_dist1d_set_xi (sd1, epdf1d->min);
  ncm_stats_dist1d_set_xf (sd1, epdf1d->max);

  return;
}

static void
_ncm_stats_dist1d_epdf_prepare (NcmStatsDist1d *sd1)
{
  NcmStatsDist1dEPDF *epdf1d = NCM_STATS_DIST1D_EPDF (sd1);

  _ncm_stats_dist1d_epdf_compact_obs (epdf1d);
  _ncm_stats_dist1d_epdf_set_bw (epdf1d);

  if (G_UNLIKELY (epdf1d->min == epdf1d->max))
  {
    ncm_stats_dist1d_set_xi (sd1, epdf1d->min);
    ncm_stats_dist1d_set_xf (sd1, epdf1d->min);
  }
  else
  {
    _ncm_stats_dist1d_epdf_update_limits (epdf1d);
  }

  return;
}

static gdouble
_ncm_stats_dist1d_epdf_get_current_h (NcmStatsDist1d *sd1)
{
  NcmStatsDist1dEPDF *epdf1d = NCM_STATS_DIST1D_EPDF (sd1);

  return epdf1d->h;
}

/**
 * ncm_stats_dist1d_epdf_new_full:
 * @max_obs: maximum observations before compacting
 * @bw: a #NcmStatsDist1dEPDFBw
 * @h_fixed: fixed bandwidth
 * @sd_min_scale: scale of the minimum distance
 *
 * Creates a new EPDF object, it creates an interpolated
 * PDF from the observations.
 *
 * Returns: (transfer full): a new #NcmStatsDist1dEPDF
 */
NcmStatsDist1dEPDF *
ncm_stats_dist1d_epdf_new_full (guint max_obs, NcmStatsDist1dEPDFBw bw, gdouble h_fixed, gdouble sd_min_scale)
{
  NcmStatsDist1dEPDF *epdf1d = g_object_new (NCM_TYPE_STATS_DIST1D_EPDF,
                                             "max-obs", max_obs,
                                             "bandwidth", bw,
                                             "h-fixed", h_fixed,
                                             "sd-min-scale", sd_min_scale,
                                             NULL);

  return epdf1d;
}

/**
 * ncm_stats_dist1d_epdf_new:
 * @sd_min_scale: scale of the minimum distance
 *
 * Creates a new EPDF object, it creates an interpolated
 * PDF from the observations.
 *
 * Returns: (transfer full): a new #NcmStatsDist1dEPDF
 */
NcmStatsDist1dEPDF *
ncm_stats_dist1d_epdf_new (gdouble sd_min_scale)
{
  NcmStatsDist1dEPDF *epdf1d = g_object_new (NCM_TYPE_STATS_DIST1D_EPDF,
                                             "sd-min-scale", sd_min_scale,
                                             NULL);

  return epdf1d;
}

/**
 * ncm_stats_dist1d_epdf_ref:
 * @epdf1d: a #NcmStatsDist1dEPDF
 *
 * Increases the reference count of @epdf1d by one.
 *
 * Returns: (transfer full): @epdf1d.
 */
NcmStatsDist1dEPDF *
ncm_stats_dist1d_epdf_ref (NcmStatsDist1dEPDF *epdf1d)
{
  return g_object_ref (epdf1d);
}

/**
 * ncm_stats_dist1d_epdf_free:
 * @epdf1d: a #NcmStatsDist1dEPDF
 *
 * Atomically decrements the reference count of @epdf1d by one. If the reference count drops to 0,
 * all memory allocated by @epdf1d is released.
 *
 */
void
ncm_stats_dist1d_epdf_free (NcmStatsDist1dEPDF *epdf1d)
{
  g_object_unref (epdf1d);
}

/**
 * ncm_stats_dist1d_epdf_clear:
 * @epdf1d: a #NcmStatsDist1dEPDF
 *
 * Atomically decrements the reference count of @epdf1d by one. If the reference count drops to 0,
 * all memory allocated by @epdf1d is released. Set the pointer to NULL;
 *
 */
void
ncm_stats_dist1d_epdf_clear (NcmStatsDist1dEPDF **epdf1d)
{
  g_clear_object (epdf1d);
}

/**
 * ncm_stats_dist1d_epdf_set_bw_type:
 * @epdf1d: a #NcmStatsDist1dEPDF
 * @bw: a #NcmStatsDist1dEPDFBw
 *
 * Sets the bandwidth computation type to @bw. The object
 * must be (re)prepared after the call to this method to be used.
 *
 */
void
ncm_stats_dist1d_epdf_set_bw_type (NcmStatsDist1dEPDF *epdf1d, NcmStatsDist1dEPDFBw bw)
{
  epdf1d->bw     = bw;
  epdf1d->bw_set = FALSE;
}

/**
 * ncm_stats_dist1d_epdf_get_bw_type:
 * @epdf1d: a #NcmStatsDist1dEPDF
 *
 * Returns: the current bandwidth computation type #NcmStatsDist1dEPDFBw.
 */
NcmStatsDist1dEPDFBw
ncm_stats_dist1d_epdf_get_bw_type (NcmStatsDist1dEPDF *epdf1d)
{
  return epdf1d->bw;
}

/**
 * ncm_stats_dist1d_epdf_set_h_fixed:
 * @epdf1d: a #NcmStatsDist1dEPDF
 * @h_fixed: fixed bandwidth
 *
 * Sets the fixed bandwidth to @h_fixed. The object
 * must be (re)prepared after the call to this method to be used.
 * This value is used only if the bandwidth computation type is
 * #NCM_STATS_DIST1D_EPDF_BW_FIXED.
 *
 */
void
ncm_stats_dist1d_epdf_set_h_fixed (NcmStatsDist1dEPDF *epdf1d, gdouble h_fixed)
{
  epdf1d->h_fixed = h_fixed;
  epdf1d->bw_set  = FALSE;
}

/**
 * ncm_stats_dist1d_epdf_get_h_fixed:
 * @epdf1d: a #NcmStatsDist1dEPDF
 *
 * Returns: the current fixed bandwidth.
 */
gdouble
ncm_stats_dist1d_epdf_get_h_fixed (NcmStatsDist1dEPDF *epdf1d)
{
  return epdf1d->h_fixed;
}

/**
 * ncm_stats_dist1d_epdf_add_obs_weight:
 * @epdf1d: a #NcmStatsDist1dEPDF
 * @x: an observation
 * @w: observation weight
 *
 * Adds a new observation @x with weight @w to the @epdf1d updating
 * the internal approximation of the EPDF when necessary.
 *
 */
void
ncm_stats_dist1d_epdf_add_obs_weight (NcmStatsDist1dEPDF *epdf1d, const gdouble x, const gdouble w)
{
  NcmStatsDist1dEPDFObs obs = {x, w};
  const gdouble new_max     = GSL_MAX (epdf1d->max, x);
  const gdouble new_min     = GSL_MIN (epdf1d->min, x);

  if (!gsl_finite (x) || (w < 0.0))
  {
    g_warning ("ncm_stats_dist1d_epdf_add_obs_weight: invalid observation %u [x = %g, w = %g], skipping...\n", epdf1d->n_obs, x, w); /* LCOV_EXCL_LINE */

    return;
  }

  if (w == 0.0)
    return;

  epdf1d->bw_set = FALSE;

  ncm_stats_vec_set (epdf1d->obs_stats, 0, x);
  ncm_stats_vec_update_weight (epdf1d->obs_stats, w);

  epdf1d->n_obs++;
  epdf1d->WT += w;

  g_array_append_val (epdf1d->obs, obs);
  epdf1d->list_sorted = FALSE;
  epdf1d->np_obs++;

  if (epdf1d->np_obs > epdf1d->max_obs)
  {
    _ncm_stats_dist1d_epdf_compact_obs (epdf1d);
    epdf1d->max_obs = GSL_MAX (10 * epdf1d->obs->len, epdf1d->max_obs);
    epdf1d->np_obs  = 0;
  }

  epdf1d->max = new_max;
  epdf1d->min = new_min;
}

/**
 * ncm_stats_dist1d_epdf_add_obs:
 * @epdf1d: a #NcmStatsDist1dEPDF
 * @x: an observation
 *
 * Adds a new observation @x (weight 1.0) to the @epdf1d updating
 * the internal approximation of the EPDF when necessary.
 *
 */
void
ncm_stats_dist1d_epdf_add_obs (NcmStatsDist1dEPDF *epdf1d, gdouble x)
{
  ncm_stats_dist1d_epdf_add_obs_weight (epdf1d, x, 1.0);
}

/**
 * ncm_stats_dist1d_epdf_reset:
 * @epdf1d: a #NcmStatsDist1dEPDF
 *
 * Empty the object @epdf1d discarding all observations.
 *
 */
void
ncm_stats_dist1d_epdf_reset (NcmStatsDist1dEPDF *epdf1d)
{
  ncm_stats_vec_reset (epdf1d->obs_stats, TRUE);
  g_array_set_size (epdf1d->obs, 0);

  epdf1d->n_obs  = 0;
  epdf1d->np_obs = 0;
  epdf1d->WT     = 0.0;
  epdf1d->min    = GSL_POSINF;
  epdf1d->max    = GSL_NEGINF;
}

/**
 * ncm_stats_dist1d_epdf_set_min:
 * @epdf1d: a #NcmStatsDist1dEPDF
 * @min: sets min observation value
 *
 * Sets the lower bound for the distribution. It may be updated if a observation
 * with value smaller than @min is added.
 *
 */
void
ncm_stats_dist1d_epdf_set_min (NcmStatsDist1dEPDF *epdf1d, const gdouble min)
{
  epdf1d->min = min;
}

/**
 * ncm_stats_dist1d_epdf_set_max:
 * @epdf1d: a #NcmStatsDist1dEPDF
 * @max: sets max observation value
 *
 * Sets the upper bound for the distribution. It may be updated if a observation
 * with value larger than @max is added.
 *
 */
void
ncm_stats_dist1d_epdf_set_max (NcmStatsDist1dEPDF *epdf1d, const gdouble max)
{
  epdf1d->max = max;
}

/**
 * ncm_stats_dist1d_epdf_get_obs_mean:
 * @epdf1d: a #NcmStatsDist1dEPDF
 *
 * Calculates the mean value of the observations.
 *
 * Returns: the mean value.
 */
gdouble
ncm_stats_dist1d_epdf_get_obs_mean (NcmStatsDist1dEPDF *epdf1d)
{
  return ncm_stats_vec_get_mean (epdf1d->obs_stats, 0);
}

