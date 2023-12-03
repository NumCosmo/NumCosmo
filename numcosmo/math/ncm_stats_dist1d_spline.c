/***************************************************************************
 *            ncm_stats_dist1d_spline.c
 *
 *  Thu February 12 16:51:07 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_dist1d_spline.c
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
 * SECTION:ncm_stats_dist1d_spline
 * @title: NcmStatsDist1dSpline
 * @short_description: One dimensional probability distribution based on a spline
 *
 * Reconstruction of an arbitrary one dimensional probability distribution based on a spline.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_stats_dist1d_spline.h"

enum
{
  PROP_0,
  PROP_M2LNP,
  PROP_TAIL_SIGMA,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcmStatsDist1dSpline, ncm_stats_dist1d_spline, NCM_TYPE_STATS_DIST1D)

static void
ncm_stats_dist1d_spline_init (NcmStatsDist1dSpline *sd1s)
{
  sd1s->m2lnp      = NULL;
  sd1s->tail_sigma = 0.0;
}

static void
ncm_stats_dist1d_spline_dispose (GObject *object)
{
  NcmStatsDist1dSpline *sd1s = NCM_STATS_DIST1D_SPLINE (object);
  
  ncm_spline_clear (&sd1s->m2lnp);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist1d_spline_parent_class)->dispose (object);
}

static void
ncm_stats_dist1d_spline_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist1d_spline_parent_class)->finalize (object);
}

static void
ncm_stats_dist1d_spline_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmStatsDist1dSpline *sd1s = NCM_STATS_DIST1D_SPLINE (object);
  
  g_return_if_fail (NCM_IS_STATS_DIST1D_SPLINE (object));
  
  switch (prop_id)
  {
    case PROP_M2LNP:
      ncm_spline_clear (&sd1s->m2lnp);
      sd1s->m2lnp = g_value_dup_object (value);
      break;
    case PROP_TAIL_SIGMA:
      sd1s->tail_sigma = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_stats_dist1d_spline_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmStatsDist1dSpline *sd1s = NCM_STATS_DIST1D_SPLINE (object);
  
  g_return_if_fail (NCM_IS_STATS_DIST1D_SPLINE (object));
  
  switch (prop_id)
  {
    case PROP_M2LNP:
      g_value_set_object (value, sd1s->m2lnp);
      break;
    case PROP_TAIL_SIGMA:
      g_value_set_double (value, sd1s->tail_sigma);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static gdouble ncm_stats_dist1d_spline_p (NcmStatsDist1d *sd1, gdouble x);
static gdouble ncm_stats_dist1d_spline_m2lnp (NcmStatsDist1d *sd1, gdouble x);
static void ncm_stats_dist1d_spline_prepare (NcmStatsDist1d *sd1);

static void
ncm_stats_dist1d_spline_class_init (NcmStatsDist1dSplineClass *klass)
{
  GObjectClass *object_class     = G_OBJECT_CLASS (klass);
  NcmStatsDist1dClass *sd1_class = NCM_STATS_DIST1D_CLASS (klass);
  
  object_class->set_property = ncm_stats_dist1d_spline_set_property;
  object_class->get_property = ncm_stats_dist1d_spline_get_property;
  object_class->dispose      = ncm_stats_dist1d_spline_dispose;
  object_class->finalize     = ncm_stats_dist1d_spline_finalize;
  
  g_object_class_install_property (object_class,
                                   PROP_M2LNP,
                                   g_param_spec_object ("m2lnp",
                                                        NULL,
                                                        "m2lnp",
                                                        NCM_TYPE_SPLINE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  g_object_class_install_property (object_class,
                                   PROP_TAIL_SIGMA,
                                   g_param_spec_double ("tail-sigma",
                                                        NULL,
                                                        "Tail sigma",
                                                        1.0e-100, G_MAXDOUBLE, 1.0e-2,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  sd1_class->p       = &ncm_stats_dist1d_spline_p;
  sd1_class->m2lnp   = &ncm_stats_dist1d_spline_m2lnp;
  sd1_class->prepare = &ncm_stats_dist1d_spline_prepare;
}

static void
ncm_stats_dist1d_spline_tail_init (NcmStatsDist1dSplineTail *tail, gdouble xb, gdouble p, gdouble d1p, gdouble d2p)
{
  tail->xb    = xb;
  tail->a     = p;
  tail->b     = d1p;
  tail->c     = d2p;
}

static gdouble
ncm_stats_dist1d_spline_tail_eval (NcmStatsDist1dSplineTail *tail, gdouble x)
{
  const gdouble xmxb   = x - tail->xb;
  const gdouble xmxb2  = xmxb * xmxb;
  
  return tail->a + tail->b * xmxb + 0.5 * tail->c * xmxb2;
}

static gdouble
ncm_stats_dist1d_spline_m2lnp (NcmStatsDist1d *sd1, gdouble x)
{
  NcmStatsDist1dSpline *sd1s = NCM_STATS_DIST1D_SPLINE (sd1);
  
  if (x < sd1s->left_tail.xb)
    return ncm_stats_dist1d_spline_tail_eval (&sd1s->left_tail, x);
  else if (x > sd1s->right_tail.xb)
    return ncm_stats_dist1d_spline_tail_eval (&sd1s->right_tail, x);
  else
    return ncm_spline_eval (sd1s->m2lnp, x);
}

static gdouble
ncm_stats_dist1d_spline_p (NcmStatsDist1d *sd1, gdouble x)
{
  const gdouble m2lnp = ncm_stats_dist1d_spline_m2lnp (sd1, x);
  
  return exp (-0.5 * m2lnp);
}

static void
ncm_stats_dist1d_spline_prepare (NcmStatsDist1d *sd1)
{
  NcmStatsDist1dSpline *sd1s = NCM_STATS_DIST1D_SPLINE (sd1);
  gdouble m2lnp_lb, m2lnp_ub;
  gdouble d1m2lnp_lb, d1m2lnp_ub;
  gdouble d2m2lnp_lb, d2m2lnp_ub;
  gdouble x_lb = 0.0, x_ub = 0.0;
  
  ncm_spline_prepare (sd1s->m2lnp);
  
  ncm_spline_get_bounds (sd1s->m2lnp, &x_lb, &x_ub);
  
  sd1->xi = x_lb;
  sd1->xf = x_ub;
  
  m2lnp_lb = ncm_spline_eval (sd1s->m2lnp, x_lb);
  m2lnp_ub = ncm_spline_eval (sd1s->m2lnp, x_ub);
  
  d1m2lnp_lb = ncm_spline_eval_deriv (sd1s->m2lnp, x_lb);
  d1m2lnp_ub = ncm_spline_eval_deriv (sd1s->m2lnp, x_ub);
  
  d2m2lnp_lb = ncm_spline_eval_deriv2 (sd1s->m2lnp, x_lb);
  d2m2lnp_ub = ncm_spline_eval_deriv2 (sd1s->m2lnp, x_ub);
  
  ncm_stats_dist1d_spline_tail_init (&sd1s->left_tail,  x_lb, m2lnp_lb, d1m2lnp_lb, d2m2lnp_lb);
  ncm_stats_dist1d_spline_tail_init (&sd1s->right_tail, x_ub, m2lnp_ub, d1m2lnp_ub, d2m2lnp_ub);
}

/**
 * ncm_stats_dist1d_spline_new:
 * @m2lnp: a #NcmSpline
 *
 * Returns a new #NcmStatsDist1dSpline where @m2lnp, $-2\ln(p(x))$, is a #NcmSpline, where $p(x)$
 * is the probability density. The probability density $p(x)$ do not need to be normalized.
 *
 * Returns: a new #NcmStatsDist1dSpline
 */
NcmStatsDist1dSpline *
ncm_stats_dist1d_spline_new (NcmSpline *m2lnp)
{
  NcmStatsDist1dSpline *sd1s = g_object_new (NCM_TYPE_STATS_DIST1D_SPLINE,
                                             "m2lnp", m2lnp,
                                             NULL);
  
  return sd1s;
}

