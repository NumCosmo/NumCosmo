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
 * Kernel density estimate of a one-dimensional distribution from weighted observations.
 *
 * With observations $x_j$ of weights $w_j$, total weight $W$ and bandwidth $h$, the density is
 * $$
 * p(x) = \frac{1}{c(x)(W + 1)}\left[\frac{1}{\sqrt{2\pi}h}\sum_j w_j e^{-(x - x_j)^2/(2h^2)} + \frac{1}{x_f - x_i}\right],
 * $$
 * a Gaussian kernel sum plus a uniform component of weight 1, where
 * $c(x) = [\mathrm{erf}((x - x_i)/(\sqrt{2}h)) + \mathrm{erf}((x_f - x)/(\sqrt{2}h))]/2$ is the
 * fraction of a kernel at $x$ inside the support $[x_i, x_f]$. The support is the range of
 * the observations, or wider when set by ncm_stats_dist1d_epdf_set_min() and
 * ncm_stats_dist1d_epdf_set_max().
 *
 * Observations closer than $\sigma\,s$, with $\sigma$ their standard deviation and $s$
 * #NcmStatsDist1dEPDF:sd-min-scale, are merged into their weighted mean. Merging happens in
 * ncm_stats_dist1d_prepare() and whenever more than #NcmStatsDist1dEPDF:max-obs observations
 * were added since the last merge; the limit then grows to ten times the merged count.
 *
 * The bandwidth follows #NcmStatsDist1dEPDF:bandwidth, see #NcmStatsDist1dEPDFBw. The
 * automatic bandwidth is a diffusion (Botev-type) plug-in selector on a $2^{14}$-bin
 * histogram of the merged observations over $[x_i - R/2, x_f + R/2]$, $R = x_f - x_i$,
 * iterated from the rule-of-thumb value. For Gaussian data it is within 20% of the
 * AMISE-optimal bandwidth for $N$ from $10^4$ to $10^6$; on the claw density of Marron and
 * Wand it is five times the optimum at $N = 10^5$. The rule-of-thumb and automatic bandwidths use the
 * number of observations $N$, not their weights, and the interquartile range ignores the
 * weights.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/stats/ncm_stats_dist1d_epdf.h"
#include "ncm/core/ncm_c.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <complex.h>
#include <fftw3.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_MAX_OBS,
  PROP_NOBS,
  PROP_BANDWIDTH,
  PROP_H_FIXED,
  PROP_SD_MIN_SCALE,
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
  gpointer fft_data_to_tilde;
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
  epdf1d->obs_stats    = ncm_stats_vec_new (1, NCM_STATS_VEC_VAR, FALSE);
  epdf1d->max_obs      = 0;
  epdf1d->bw           = NCM_STATS_DIST1D_EPDF_BW_LEN;
  epdf1d->h_fixed      = 0.0;
  epdf1d->sd_min_scale = 0.0;
  epdf1d->h            = 0.0;
  epdf1d->n_obs        = 0;
  epdf1d->np_obs       = 0;
  epdf1d->WT           = 0.0;
  epdf1d->obs          = NULL;
  epdf1d->min          = GSL_POSINF;
  epdf1d->max          = GSL_NEGINF;

  epdf1d->fftsize           = 0;
  epdf1d->Iv                = NULL;
  epdf1d->p_data            = NULL;
  epdf1d->p_tilde           = NULL;
  epdf1d->p_tilde2          = NULL;
  epdf1d->fft_data_to_tilde = NULL;

  epdf1d->bw_set = FALSE;

  ncm_stats_vec_enable_quantile (epdf1d->obs_stats, 0.5);
}

static void
ncm_stats_dist1d_epdf_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_stats_dist1d_epdf_parent_class)->constructed (object);
  {
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

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist1d_epdf_parent_class)->dispose (object);
}

static void
ncm_stats_dist1d_epdf_finalize (GObject *object)
{
  NcmStatsDist1dEPDF *epdf1d = NCM_STATS_DIST1D_EPDF (object);

  g_clear_pointer (&epdf1d->fft_data_to_tilde, ncm_cfg_fftw_plan_destroy);

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
                                                      "Number of added observations that triggers a merge",
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
                                                        "Merging distance in units of the standard deviation",
                                                        1.0e-20, 1.0e20, 1.0e-3,
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

static gdouble _ncm_stats_dist1d_epdf_p_gk (NcmStatsDist1dEPDF *epdf1d, gdouble x);

static void
_ncm_stats_dist1d_epdf_compact_obs (NcmStatsDist1dEPDF *epdf1d)
{
  guint i, j;
  gint obs_len = epdf1d->obs->len;

  if (epdf1d->list_sorted)
    return;

  g_array_sort (epdf1d->obs, _ncm_stats_dist1d_epdf_cmp_double);
  j = 0;

  {
    NcmStatsDist1dEPDFObs *obs_j = &g_array_index (epdf1d->obs, NcmStatsDist1dEPDFObs, j);
    const gdouble sd_e           = ncm_stats_vec_get_sd (epdf1d->obs_stats, 0);
    const gdouble min_dist       = sd_e * epdf1d->sd_min_scale;


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

/* w[l][i] = I_i^l p_tilde2_i; the sum stops where exp underflows to exactly zero, since the I_i increase */
static gdouble
_ncm_stats_dist1d_epdf_estimate_df2 (gdouble * const *w, NcmVector *Iv, const guint n, const guint l, const gdouble t)
{
  const gdouble pi2  = M_PI * M_PI;
  const gdouble pi2l = gsl_pow_int (pi2, l);
  gdouble s          = 0.0;
  guint i;

  for (i = 0; i < n; i++)
  {
    const gdouble Ii = ncm_vector_fast_get (Iv, i);

    if (Ii * pi2 * t > 746.0)
      break;

    s += w[l][i] * exp (-Ii * pi2 * t);
  }

  return 0.5 * pi2l * s;
}

static gdouble
_ncm_stats_dist1d_epdf_estimate_h (gdouble * const *w, NcmVector *Iv, const guint obs_len, const guint n, const guint l, const gdouble t)
{
  const gdouble df2     = _ncm_stats_dist1d_epdf_estimate_df2 (w, Iv, n, l, t);
  const gdouble ln_Ndf2 = log (df2 * obs_len);
  const gdouble lp05    = l + 0.5;
  const gdouble ln_fact = log1p (exp2 (-lp05)) + lp05 * M_LN2  - ncm_c_lnpi () + lgamma (lp05) - ncm_c_ln3 ();
  const gdouble tn      = exp ((ln_fact - ln_Ndf2) / (1.0 + lp05));

  g_assert (l >= 2);

  if (l == 2)
  {
    const gdouble df2s = _ncm_stats_dist1d_epdf_estimate_df2 (w, Iv, n, l, tn);

    return pow (2.0 * obs_len * ncm_c_pi () * df2s, -2.0 / 5.0);
  }
  else
  {
    return _ncm_stats_dist1d_epdf_estimate_h (w, Iv, obs_len, n, l - 1, tn);
  }
}

static void
_ncm_stats_dist1d_epdf_autobw (NcmStatsDist1dEPDF *epdf1d)
{
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

    epdf1d->fftsize = nbins;

    {
      G_LOCK_DEFINE_STATIC (prepare_fft_lock);

      gboolean first;

      G_LOCK (prepare_fft_lock);

      first = ncm_cfg_fftw_plan_begin ("ncm_stats_dist1d_epdf_redft10_01_%u", nbins);

      epdf1d->fft_data_to_tilde = fftw_plan_r2r_1d (nbins, ncm_vector_data (epdf1d->p_data), ncm_vector_data (epdf1d->p_tilde),
                                                    FFTW_REDFT10, fftw_default_flags | FFTW_DESTROY_INPUT);

      ncm_cfg_fftw_plan_end (first);

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
    gdouble *w[8];
    guint l;

    for (l = 2; l < 8; l++)
    {
      w[l] = g_new (gdouble, nbins);

      for (i = 0; i < nbins; i++)
        w[l][i] = gsl_pow_int (ncm_vector_fast_get (epdf1d->Iv, i), l) * ncm_vector_fast_get (epdf1d->p_tilde2, i);
    }

    w[0] = w[1] = NULL;

    j = 0;

    while (fabs (1.0 - tn / t) > 1.0e-7)
    {
      const gdouble tni = _ncm_stats_dist1d_epdf_estimate_h (w, epdf1d->Iv, epdf1d->n_obs /*obs_len*/, nbins, 7, t);

      tn = t;
      t  = tni;
      j++;

      if (j >= 10000)
        g_error ("_ncm_stats_dist1d_epdf_autobw: too many steps to find bandwidth.");  /* LCOV_EXCL_LINE */
    }

    epdf1d->h = sqrt (t) * delta_l;

    for (l = 2; l < 8; l++)
      g_free (w[l]);
  }
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

  return _ncm_stats_dist1d_epdf_p_gk (epdf1d, x);
}

static gdouble
_ncm_stats_dist1d_epdf_m2lnp (NcmStatsDist1d *sd1, gdouble x)
{
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

  if (epdf1d->n_obs == 0)
    g_error ("_ncm_stats_dist1d_epdf_prepare: no observations.");

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
 * @max_obs: number of added observations that triggers a merge
 * @bw: a #NcmStatsDist1dEPDFBw
 * @h_fixed: bandwidth for #NCM_STATS_DIST1D_EPDF_BW_FIXED
 * @sd_min_scale: merging distance in units of the standard deviation
 *
 * Creates a new #NcmStatsDist1dEPDF, see #NcmStatsDist1dEPDF:max-obs,
 * #NcmStatsDist1dEPDF:bandwidth, #NcmStatsDist1dEPDF:h-fixed and
 * #NcmStatsDist1dEPDF:sd-min-scale.
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
 * @sd_min_scale: merging distance in units of the standard deviation
 *
 * Creates a new #NcmStatsDist1dEPDF with the automatic bandwidth and the default
 * #NcmStatsDist1dEPDF:max-obs.
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
 * Decreases the reference count of @epdf1d by one.
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
 * Decreases the reference count of *@epdf1d by one and sets *@epdf1d to %NULL.
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
 * Sets #NcmStatsDist1dEPDF:bandwidth; takes effect at the next ncm_stats_dist1d_prepare().
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
 * Returns: the bandwidth type.
 */
NcmStatsDist1dEPDFBw
ncm_stats_dist1d_epdf_get_bw_type (NcmStatsDist1dEPDF *epdf1d)
{
  return epdf1d->bw;
}

/**
 * ncm_stats_dist1d_epdf_set_h_fixed:
 * @epdf1d: a #NcmStatsDist1dEPDF
 * @h_fixed: bandwidth
 *
 * Sets #NcmStatsDist1dEPDF:h-fixed, the bandwidth of #NCM_STATS_DIST1D_EPDF_BW_FIXED; takes
 * effect at the next ncm_stats_dist1d_prepare().
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
 * Returns: the bandwidth of #NCM_STATS_DIST1D_EPDF_BW_FIXED.
 */
gdouble
ncm_stats_dist1d_epdf_get_h_fixed (NcmStatsDist1dEPDF *epdf1d)
{
  return epdf1d->h_fixed;
}

/**
 * ncm_stats_dist1d_epdf_add_obs_weight:
 * @epdf1d: a #NcmStatsDist1dEPDF
 * @x: observation
 * @w: weight, non-negative
 *
 * Adds the observation @x with weight @w; it enters the estimate at the next
 * ncm_stats_dist1d_prepare(). A zero weight is ignored; a non-finite @x or a negative @w
 * is skipped with a warning.
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
 * @x: observation
 *
 * Adds the observation @x with weight 1, see ncm_stats_dist1d_epdf_add_obs_weight().
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
 * Discards all observations and the bounds; ncm_stats_dist1d_prepare() aborts until an
 * observation is added.
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
 * @min: lower bound
 *
 * Sets the lower bound of the support; a smaller observation added later lowers it. Takes
 * effect at the next ncm_stats_dist1d_prepare().
 */
void
ncm_stats_dist1d_epdf_set_min (NcmStatsDist1dEPDF *epdf1d, const gdouble min)
{
  epdf1d->min    = min;
  epdf1d->bw_set = FALSE;
}

/**
 * ncm_stats_dist1d_epdf_set_max:
 * @epdf1d: a #NcmStatsDist1dEPDF
 * @max: upper bound
 *
 * Sets the upper bound of the support; a larger observation added later raises it. Takes
 * effect at the next ncm_stats_dist1d_prepare().
 */
void
ncm_stats_dist1d_epdf_set_max (NcmStatsDist1dEPDF *epdf1d, const gdouble max)
{
  epdf1d->max    = max;
  epdf1d->bw_set = FALSE;
}

/**
 * ncm_stats_dist1d_epdf_get_obs_mean:
 * @epdf1d: a #NcmStatsDist1dEPDF
 *
 * Returns: the weighted mean of the observations.
 */
gdouble
ncm_stats_dist1d_epdf_get_obs_mean (NcmStatsDist1dEPDF *epdf1d)
{
  return ncm_stats_vec_get_mean (epdf1d->obs_stats, 0);
}

