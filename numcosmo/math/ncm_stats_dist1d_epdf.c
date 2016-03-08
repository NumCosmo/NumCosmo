/***************************************************************************
 *            ncm_stats_dist1d_epdf.c
 *
 *  Sat March 14 19:32:06 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_stats_dist1d_epdf.c
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_stats_dist1d_epdf
 * @title: NcmStatsDist1dEPDF
 * @short_description: One dimensional probability distribution based on an EPDF.
 * 
 * Empirical Probability Distribution Function (EPDF).
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_stats_dist1d_epdf.h"

enum
{
  PROP_0,
  PROP_MAX_OBS,
  PROP_NOBS,
  PROP_SD_SCALE,
  PROP_SD_MIN_SCALE,
  PROP_OUTLIERS_THRESHOLD,
};

G_DEFINE_TYPE (NcmStatsDist1dEPDF, ncm_stats_dist1d_epdf, NCM_TYPE_STATS_DIST1D);

#define _NCM_STATS_DIST1D_SROT(sd,R,n) (0.9 * GSL_MIN ((sd), ((R) / 1.34)) * pow (n * 1.0, -1.0 / 5.0))

gint 
_ncm_stats_dist1d_epdf_cmp_double (gconstpointer a,
                                   gconstpointer b)
{
#define A (*((gdouble *)a))
#define B (*((gdouble *)b))
  return (A == B) ? 0.0 : ((A < B) ? -1 : 1);  
}
#undef A
#undef B

typedef struct _NcmStatsDist1dEPDFObs
{
  gdouble x;
  gdouble w;
  gdouble n;
} NcmStatsDist1dEPDFObs;

#define _NCM_STATS_DIST1D_EPDF_OBS_N(obs,epdf1d,sd) \
0.5 * (erf (((obs)->x - (epdf1d)->min) / (M_SQRT2 * (sd))) + erf (((epdf1d)->max - (obs)->x) / (M_SQRT2 * (sd))))

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
    const gdouble sd_e     = ncm_stats_vec_get_sd (epdf1d->obs_stats, 0);
    const gdouble min_dist = sd_e * epdf1d->sd_min_scale;
    const gdouble R_e      = ncm_stats_vec_get_quantile_spread (epdf1d->obs_stats, 0);
    const gdouble sd       = _NCM_STATS_DIST1D_SROT (sd_e, R_e, epdf1d->n_obs);

    for (i = 1; i < epdf1d->obs->len; i++)
    {
      NcmStatsDist1dEPDFObs *obs_i = &g_array_index (epdf1d->obs, NcmStatsDist1dEPDFObs, i);
      
      if (fabs (obs_j->x - obs_i->x) < min_dist)
      {
#ifdef _NCM_STATS_DIST1D_EPDF_GROUPING_ALL  
        gdouble wx = (obs_j->x * obs_j->w + obs_i->x * obs_i->w);
        gdouble w  = obs_i->w + obs_j->w;
        guint k;
        obs_len--;
        for (k = i + 1; k < epdf1d->obs->len; k++)
        {
          NcmStatsDist1dEPDFObs *obs_ik = &g_array_index (epdf1d->obs, NcmStatsDist1dEPDFObs, k);
          if (fabs (obs_ik->x - obs_i->x) < min_dist)
          {
            wx += obs_ik->x * obs_ik->w;
            w  += obs_ik->w;
            obs_len--;
            i++;
            obs_i = obs_ik;
          }
          else
            break;
        }
        obs_j->x = wx / w;
        obs_j->w = w;
#else
        obs_j->x = (obs_j->x * obs_j->w + obs_i->x * obs_i->w) / (obs_i->w + obs_j->w);
        obs_j->w = obs_i->w + obs_j->w;

        obs_j->n = _NCM_STATS_DIST1D_EPDF_OBS_N (obs_j, epdf1d, sd);
        
        obs_len--;
#endif /* _NCM_STATS_DIST1D_EPDF_GROUPING_ALL */
      }
      else
      {
        j++;
        
        obs_i->n = _NCM_STATS_DIST1D_EPDF_OBS_N (obs_i, epdf1d, sd);        

        if (i != j)
        {
          g_array_index (epdf1d->obs, NcmStatsDist1dEPDFObs, j) = *obs_i; 
        }
        obs_j = &g_array_index (epdf1d->obs, NcmStatsDist1dEPDFObs, j);
      }
    }
  }

  g_array_set_size (epdf1d->obs, obs_len);
  epdf1d->list_sorted = TRUE;
}

static void
_ncm_stats_dist1d_epdf_trim_outliers (NcmStatsDist1dEPDF *epdf1d)
{
  const gdouble obs_mean = ncm_stats_vec_get_mean (epdf1d->obs_stats, 0);
  const gdouble obs_sd   = ncm_stats_vec_get_sd (epdf1d->obs_stats, 0);
  const gdouble max_dist = obs_sd * epdf1d->outliers_threshold;
  gint obs_len = epdf1d->obs->len;
  gdouble min = GSL_POSINF;
  gdouble max = GSL_NEGINF;
  register guint i, j;

  ncm_stats_vec_reset (epdf1d->obs_stats, TRUE);
  g_array_sort (epdf1d->obs, _ncm_stats_dist1d_epdf_cmp_double);

  j = 0;
  for (i = 0; i < epdf1d->obs->len; i++)
  {
    NcmStatsDist1dEPDFObs *obs_c = &g_array_index (epdf1d->obs, NcmStatsDist1dEPDFObs, i);
    if (fabs (obs_c->x - obs_mean) < max_dist)
    {
      if (i != j)
      {
        g_array_index (epdf1d->obs, NcmStatsDist1dEPDFObs, j) = obs_c[0]; 
      }
      min = GSL_MIN (min, obs_c->x);
      max = GSL_MAX (max, obs_c->x);
      ncm_stats_vec_set (epdf1d->obs_stats, 0, obs_c->x);
      ncm_stats_vec_update_weight (epdf1d->obs_stats, obs_c->w);      
      j++;
    }
    else
    {
      ncm_message ("# discarding outlier % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g...\n", obs_c->x, obs_c->w, obs_mean, obs_sd, max_dist);
      obs_len--;
    }
  }

  epdf1d->min = min;
  epdf1d->max = max;
  g_array_set_size (epdf1d->obs, obs_len);
}



static void
ncm_stats_dist1d_epdf_init (NcmStatsDist1dEPDF *epdf1d)
{
  epdf1d->obs_stats          = ncm_stats_vec_new (1, NCM_STATS_VEC_VAR, FALSE);
  epdf1d->max_obs            = 0;
  epdf1d->min_knots          = 0;
  epdf1d->max_knots          = 0;
  epdf1d->sd_scale           = 0.0;
  epdf1d->sd_min_scale       = 0.0;
  epdf1d->outliers_threshold = 0.0;
  epdf1d->sd                 = 0.0;
  epdf1d->n_obs              = 0;
  epdf1d->np_obs             = 0;
  epdf1d->obs                = NULL;
  epdf1d->min                = GSL_POSINF;
  epdf1d->max                = GSL_NEGINF;

  ncm_stats_vec_enable_quantile (epdf1d->obs_stats, 0.5);
}

static void
ncm_stats_dist1d_epdf_constructed (GObject *object)
{
  /*NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (object);*/
  NcmStatsDist1dEPDF *epdf1d = NCM_STATS_DIST1D_EPDF (object);
  epdf1d->obs = g_array_sized_new (FALSE, FALSE, sizeof (NcmStatsDist1dEPDFObs), epdf1d->max_obs);
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
    case PROP_SD_SCALE:
      epdf1d->sd_scale = g_value_get_double (value);
      break;
    case PROP_SD_MIN_SCALE:
      epdf1d->sd_min_scale = g_value_get_double (value);
      break;
    case PROP_OUTLIERS_THRESHOLD:
      epdf1d->outliers_threshold = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
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
    case PROP_SD_SCALE:
      g_value_set_double (value, epdf1d->sd_scale);
      break;
    case PROP_SD_MIN_SCALE:
      g_value_set_double (value, epdf1d->sd_min_scale);
      break;
    case PROP_OUTLIERS_THRESHOLD:
      g_value_set_double (value, epdf1d->outliers_threshold);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_stats_dist1d_epdf_dispose (GObject *object)
{
  NcmStatsDist1dEPDF *epdf1d = NCM_STATS_DIST1D_EPDF (object);

  ncm_stats_vec_clear (&epdf1d->obs_stats);
  g_clear_pointer (&epdf1d->obs, g_array_unref);
  
  /* Chain up : end */  
  G_OBJECT_CLASS (ncm_stats_dist1d_epdf_parent_class)->dispose (object);
}

static void
ncm_stats_dist1d_epdf_finalize (GObject *object)
{
  /*NcmStatsDist1dEPDF *epdf1d = NCM_STATS_DIST1D_EPDF (object);*/
  

  /* Chain up : end */  
  G_OBJECT_CLASS (ncm_stats_dist1d_epdf_parent_class)->finalize (object);
}

static gdouble ncm_stats_dist1d_epdf_p (NcmStatsDist1d *sd1, gdouble x);
static gdouble ncm_stats_dist1d_epdf_m2lnp (NcmStatsDist1d *sd1, gdouble x);
static void ncm_stats_dist1d_epdf_prepare (NcmStatsDist1d *sd1);

static void
ncm_stats_dist1d_epdf_class_init (NcmStatsDist1dEPDFClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
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
                                                      "Maximum saved observations",
                                                      10, G_MAXUINT, 1000,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NOBS,
                                   g_param_spec_uint ("n-obs",
                                                      NULL,
                                                      "Number of observations",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SD_SCALE,
                                   g_param_spec_double ("sd-scale",
                                                        NULL,
                                                        "Percentage of the standard deviation to use",
                                                        1.0e-5, 1.0e5, 0.1,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SD_MIN_SCALE,
                                   g_param_spec_double ("sd-min-scale",
                                                        NULL,
                                                        "Percentage of the standard deviation to use as minimum distance",
                                                        1.0e-20, 1.0e20, 1.0e-2,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_OUTLIERS_THRESHOLD,
                                   g_param_spec_double ("outliers-threshold",
                                                        NULL,
                                                        "How many sigmas to consider an outlier",
                                                        1.0, 1000.0, 20.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  sd1_class->p       = &ncm_stats_dist1d_epdf_p;
  sd1_class->m2lnp   = &ncm_stats_dist1d_epdf_m2lnp;
  sd1_class->prepare = &ncm_stats_dist1d_epdf_prepare;
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
    return l;
}

static gdouble 
ncm_stats_dist1d_epdf_p_gk (NcmStatsDist1dEPDF *epdf1d, gdouble x)
{
  NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (epdf1d);
  const gdouble sd_e  = ncm_stats_vec_get_sd (epdf1d->obs_stats, 0);
  const gdouble R_e   = ncm_stats_vec_get_quantile_spread (epdf1d->obs_stats, 0);
  const gdouble sd    = _NCM_STATS_DIST1D_SROT (sd_e, R_e, epdf1d->n_obs);
  gdouble res = 0.0;
  gint i;

  _ncm_stats_dist1d_epdf_compact_obs (epdf1d);

#ifndef _NCM_STATS_DIST1D_EPDF_DIRECT 
  {
    gint s = _ncm_stats_dist1d_epdf_bsearch (epdf1d->obs, x, 0, epdf1d->obs->len - 1);

    for (i = s; i < epdf1d->obs->len; i++)
    {
      NcmStatsDist1dEPDFObs *obs = &g_array_index (epdf1d->obs, NcmStatsDist1dEPDFObs, i);
      const gdouble x_i   = obs->x;
      const gdouble de_i  = (x - x_i) / sd;
      const gdouble de2_i = de_i * de_i;
      const gdouble exp_i = obs->w * exp (-0.5 * de2_i);
      res += exp_i / obs->n;
      /*printf ("% 20.15g % 20.15g % 20.15g [% 20.15g % 20.15g] % 20.15g\n", x_i, obs->w, obs->n, epdf1d->min, epdf1d->max, x);*/
      if (exp_i / res < GSL_DBL_EPSILON)
        break;
    }
    for (i = s - 1; i >= 0; i--)
    {
      NcmStatsDist1dEPDFObs *obs = &g_array_index (epdf1d->obs, NcmStatsDist1dEPDFObs, i);
      const gdouble x_i   = obs->x;
      const gdouble de_i  = (x - x_i) / sd;
      const gdouble de2_i = de_i * de_i;
      const gdouble exp_i = obs->w * exp (-0.5 * de2_i);
      res += exp_i / obs->n;
      /*printf ("% 20.15g % 20.15g % 20.15g [% 20.15g % 20.15g] % 20.15g\n", x_i, obs->w, obs->n, epdf1d->min, epdf1d->max, x);*/
      if (exp_i / res < GSL_DBL_EPSILON)
        break;
    }
    /*printf ("#----------------------------------------------------------------------------------------#\n");*/
  }
#else   
  for (i = 0; i < epdf1d->obs->len; i++)
  {
    NcmStatsDist1dEPDFObs *obs = &g_array_index (epdf1d->obs, NcmStatsDist1dEPDFObs, i);
    const gdouble x_i   = obs->x;
    const gdouble de_i  = (x - x_i) / sd;
    const gdouble de2_i = de_i * de_i;
    const gdouble exp_i = obs->w * exp (-0.5 * de2_i);
    res += exp_i / obs->n;
  }
#endif /* _NCM_STATS_DIST1D_EPDF_DIRECT */

  return (res / (sqrt (2.0 * M_PI) * sd) + 1.0 / (sd1->xf - sd1->xi)) / (epdf1d->n_obs + 1.0);
}

static gdouble 
ncm_stats_dist1d_epdf_p (NcmStatsDist1d *sd1, gdouble x)
{
  NcmStatsDist1dEPDF *epdf1d = NCM_STATS_DIST1D_EPDF (sd1);
  return ncm_stats_dist1d_epdf_p_gk (epdf1d, x);
}

static gdouble 
ncm_stats_dist1d_epdf_m2lnp (NcmStatsDist1d *sd1, gdouble x)
{
  return -2.0 * log (ncm_stats_dist1d_epdf_p (sd1, x));
}

static void 
ncm_stats_dist1d_epdf_update_limits (NcmStatsDist1dEPDF *epdf1d)
{
  NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (epdf1d);
  const gdouble sd_e = ncm_stats_vec_get_sd (epdf1d->obs_stats, 0);
  const gdouble R_e  = ncm_stats_vec_get_quantile_spread (epdf1d->obs_stats, 0);
  const gdouble sd = _NCM_STATS_DIST1D_SROT (sd_e, R_e, epdf1d->n_obs);
  const gdouble lnNorma2 = 2.0 * log (1.0 * epdf1d->n_obs) + log (2.0 * M_PI) + 2.0 * log (sd);
  const gdouble exp_b    = -2.0 * NCM_STATS_DIST1D (epdf1d)->reltol - lnNorma2;
  /*const gdouble exp_b    = -2.0 * GSL_LOG_DBL_EPSILON - lnNorma2;*/

  epdf1d->sd = sd;
  if (exp_b < 0.0)
  {
    sd1->xi = epdf1d->min;
    sd1->xf = epdf1d->max;
  }
  else
  {
    sd1->xi = epdf1d->min - sqrt (exp_b) * sd;
    sd1->xf = epdf1d->max + sqrt (exp_b) * sd;
  }
  return;
}

static void 
ncm_stats_dist1d_epdf_prepare (NcmStatsDist1d *sd1)
{
  NcmStatsDist1dEPDF *epdf1d = NCM_STATS_DIST1D_EPDF (sd1);

  if (FALSE)
  {
    gdouble sd = ncm_stats_vec_get_sd (epdf1d->obs_stats, 0);
    gdouble sd_old = 0.0;

    do {
      _ncm_stats_dist1d_epdf_trim_outliers (epdf1d);
      sd_old = sd;
      sd = ncm_stats_vec_get_sd (epdf1d->obs_stats, 0);
      printf ("# sd % 20.15g % 20.15g\n", sd, sd_old);
    } while (sd / sd_old < 0.9);
  }
  
  if (G_UNLIKELY (epdf1d->min == epdf1d->max))
    sd1->xi = sd1->xf = epdf1d->min;
  else
    ncm_stats_dist1d_epdf_update_limits (epdf1d);
  
  /*sd1->abstol = 0.0;*/
  return;
}

/**
 * ncm_stats_dist1d_epdf_new:
 * @max_obs: maximum saved observations
 * @sd_scale: standard deviation scale for basis functions
 * @sd_min_scale: scale of the minimum distance 
 * 
 * Creates a new EPDF object, it creates an interpolated
 * PDF from the observations. 
 * 
 * Returns: a new #NcmStatsDist1dEPDF
 */
NcmStatsDist1dEPDF *
ncm_stats_dist1d_epdf_new (guint max_obs, gdouble sd_scale, gdouble sd_min_scale)
{
  NcmStatsDist1dEPDF *epdf1d = g_object_new (NCM_TYPE_STATS_DIST1D_EPDF,
                                             "max-obs", max_obs,
                                             "sd-scale", sd_scale,
                                             "sd-min-scale", sd_min_scale,
                                             NULL);
  return epdf1d;
}

/**
 * ncm_stats_dist1d_epdf_add_obs:
 * @epdf1d: a #NcmStatsDist1dEPDF
 * @x: an observation
 * 
 * Adds a new observation @x to the @epdf1d updating the internal
 * approximation of the EPDF when necessary.
 * 
 */
void
ncm_stats_dist1d_epdf_add_obs (NcmStatsDist1dEPDF *epdf1d, gdouble x)
{
  const NcmStatsDist1dEPDFObs obs = {x, 1.0, 1.0};
  const gdouble new_max = GSL_MAX (epdf1d->max, x);
  const gdouble new_min = GSL_MIN (epdf1d->min, x);

  if (!gsl_finite (x))
  {
    g_warning ("# non-finite observation number %u [%g], skipping...\n", epdf1d->n_obs, x);
    return;
  }

  ncm_stats_vec_set (epdf1d->obs_stats, 0, x);
  ncm_stats_vec_update (epdf1d->obs_stats);

  g_array_append_val (epdf1d->obs, obs);
  epdf1d->list_sorted = FALSE;

  epdf1d->sd = ncm_stats_vec_get_sd (epdf1d->obs_stats, 0) * epdf1d->sd_scale;
  epdf1d->n_obs++;
  epdf1d->np_obs++;

  if (epdf1d->np_obs > epdf1d->max_obs)
  {    
    _ncm_stats_dist1d_epdf_compact_obs (epdf1d);
    epdf1d->np_obs = 0;
    /*printf ("# new size %u %u\n", epdf1d->obs->len, epdf1d->max_obs);*/
  }
  
  epdf1d->max         = new_max;
  epdf1d->min         = new_min;
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
  epdf1d->n_obs = 0;
  epdf1d->np_obs = 0;
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
