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
 * NcmStatsDist1dSpline:
 *
 * One-dimensional distribution tabulated by a spline, of $-2\ln p(x)$ or of $p(x)$.
 *
 * Exactly one of #NcmStatsDist1dSpline:m2lnp and #NcmStatsDist1dSpline:density is set. The
 * support is the spline's knot range, which ncm_stats_dist1d_prepare() sets as
 * $[x_i, x_f]$, replacing any value set before. The density need not be normalized.
 *
 * With #NcmStatsDist1dSpline:m2lnp, a spline $m_2(x)$, the density is
 * $p(x) = e^{-[m_2(x) - m_\mathrm{min}]/2}$, with $m_\mathrm{min}$ the smallest knot
 * value, so ncm_stats_dist1d_eval_norma() is relative to $e^{-m_\mathrm{min}/2}$. Outside
 * the knot range, at a distance $\delta$ from the nearest bound $x_b$, $m_2$ continues as
 * $m_2(x_b) \pm m_2^\prime(x_b)\,\delta + c\,\delta^2/2$ with
 * $c = \max[m_2^{\prime\prime}(x_b), 1/(x_f - x_i)^2, s^2/2]$, where the last term enters
 * only when the outward slope $s$ is negative. The continuation matches the value and the
 * slope at $x_b$, and also the second derivative when that is the largest term. It grows
 * to $+\infty$, so the density goes to zero, and $m_2$ dips below $m_2(x_b)$ by at most 1.
 *
 * With #NcmStatsDist1dSpline:density, a spline $P(x)$, the density is $\max[P(x), 0]$ in the
 * knot range and zero outside it.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/stats/ncm_stats_dist1d_spline.h"

enum
{
  PROP_0,
  PROP_M2LNP,
  PROP_DENSITY,
};

typedef struct _NcmStatsDist1dSplineTail
{
  gdouble xb;
  gdouble a;
  gdouble b;
  gdouble c;
} NcmStatsDist1dSplineTail;

struct _NcmStatsDist1dSpline
{
  /*< private >*/
  NcmStatsDist1d parent_instance;
  NcmSpline *m2lnp;
  NcmSpline *p;
  gdouble m2lnp_min;
  gdouble x_lb;
  gdouble x_ub;
  NcmStatsDist1dSplineTail left_tail;
  NcmStatsDist1dSplineTail right_tail;
};

G_DEFINE_TYPE (NcmStatsDist1dSpline, ncm_stats_dist1d_spline, NCM_TYPE_STATS_DIST1D)

static void
ncm_stats_dist1d_spline_init (NcmStatsDist1dSpline *sd1s)
{
  sd1s->m2lnp     = NULL;
  sd1s->p         = NULL;
  sd1s->m2lnp_min = 0.0;
  sd1s->x_lb      = 0.0;
  sd1s->x_ub      = 0.0;
}

static void
ncm_stats_dist1d_spline_dispose (GObject *object)
{
  NcmStatsDist1dSpline *sd1s = NCM_STATS_DIST1D_SPLINE (object);

  ncm_spline_clear (&sd1s->m2lnp);
  ncm_spline_clear (&sd1s->p);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist1d_spline_parent_class)->dispose (object);
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
    case PROP_DENSITY:
      ncm_spline_clear (&sd1s->p);
      sd1s->p = g_value_dup_object (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
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
    case PROP_DENSITY:
      g_value_set_object (value, sd1s->p);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
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

  g_object_class_install_property (object_class,
                                   PROP_M2LNP,
                                   g_param_spec_object ("m2lnp",
                                                        NULL,
                                                        "Spline of -2 ln p",
                                                        NCM_TYPE_SPLINE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_DENSITY,
                                   g_param_spec_object ("density",
                                                        NULL,
                                                        "Spline of the density p",
                                                        NCM_TYPE_SPLINE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  sd1_class->p       = &ncm_stats_dist1d_spline_p;
  sd1_class->m2lnp   = &ncm_stats_dist1d_spline_m2lnp;
  sd1_class->prepare = &ncm_stats_dist1d_spline_prepare;
}

static gdouble
ncm_stats_dist1d_spline_tail_eval (NcmStatsDist1dSplineTail *tail, gdouble x)
{
  const gdouble xmxb = x - tail->xb;

  return tail->a + tail->b * xmxb + 0.5 * tail->c * xmxb * xmxb;
}

static gdouble
ncm_stats_dist1d_spline_m2lnp (NcmStatsDist1d *sd1, gdouble x)
{
  NcmStatsDist1dSpline *sd1s = NCM_STATS_DIST1D_SPLINE (sd1);

  if (sd1s->p != NULL)
  {
    const gdouble p = ncm_stats_dist1d_spline_p (sd1, x);

    return (p > 0.0) ? -2.0 * log (p) : GSL_POSINF;
  }

  if (x < sd1s->x_lb)
    return ncm_stats_dist1d_spline_tail_eval (&sd1s->left_tail, x);
  else if (x > sd1s->x_ub)
    return ncm_stats_dist1d_spline_tail_eval (&sd1s->right_tail, x);
  else
    return ncm_spline_eval (sd1s->m2lnp, x);
}

static gdouble
ncm_stats_dist1d_spline_p (NcmStatsDist1d *sd1, gdouble x)
{
  NcmStatsDist1dSpline *sd1s = NCM_STATS_DIST1D_SPLINE (sd1);

  if (sd1s->p != NULL)
  {
    if ((x < sd1s->x_lb) || (x > sd1s->x_ub))
      return 0.0;

    return GSL_MAX (ncm_spline_eval (sd1s->p, x), 0.0);
  }

  return exp (-0.5 * (ncm_stats_dist1d_spline_m2lnp (sd1, x) - sd1s->m2lnp_min));
}

/* Continuation of m2lnp beyond the bound xb; outward is the sign of the outward direction */
static void
_ncm_stats_dist1d_spline_tail_init (NcmStatsDist1dSplineTail *tail, NcmSpline *m2lnp, gdouble xb, gdouble outward, gdouble L)
{
  const gdouble d1 = ncm_spline_eval_deriv (m2lnp, xb);
  const gdouble d2 = ncm_spline_eval_deriv2 (m2lnp, xb);
  const gdouble s  = outward * d1;
  gdouble c        = GSL_MAX (d2, 1.0 / (L * L));

  if (s < 0.0)
    c = GSL_MAX (c, 0.5 * s * s);

  tail->xb = xb;
  tail->a  = ncm_spline_eval (m2lnp, xb);
  tail->b  = d1;
  tail->c  = c;
}

static void
ncm_stats_dist1d_spline_prepare (NcmStatsDist1d *sd1)
{
  NcmStatsDist1dSpline *sd1s = NCM_STATS_DIST1D_SPLINE (sd1);

  if ((sd1s->m2lnp == NULL) == (sd1s->p == NULL))
    g_error ("ncm_stats_dist1d_spline_prepare: exactly one of the m2lnp and density splines must be set.");

  if (sd1s->p != NULL)
  {
    ncm_spline_prepare (sd1s->p);
    ncm_spline_get_bounds (sd1s->p, &sd1s->x_lb, &sd1s->x_ub);
  }
  else
  {
    NcmVector *yv = ncm_spline_peek_yv (sd1s->m2lnp);
    gdouble L;

    ncm_spline_prepare (sd1s->m2lnp);
    ncm_spline_get_bounds (sd1s->m2lnp, &sd1s->x_lb, &sd1s->x_ub);

    sd1s->m2lnp_min = ncm_vector_get_min (yv);
    L               = sd1s->x_ub - sd1s->x_lb;

    _ncm_stats_dist1d_spline_tail_init (&sd1s->left_tail,  sd1s->m2lnp, sd1s->x_lb, -1.0, L);
    _ncm_stats_dist1d_spline_tail_init (&sd1s->right_tail, sd1s->m2lnp, sd1s->x_ub, +1.0, L);
  }

  ncm_stats_dist1d_set_xi (sd1, sd1s->x_lb);
  ncm_stats_dist1d_set_xf (sd1, sd1s->x_ub);
}

/**
 * ncm_stats_dist1d_spline_new:
 * @m2lnp: a #NcmSpline of $-2\ln p(x)$
 *
 * Creates a new #NcmStatsDist1dSpline with #NcmStatsDist1dSpline:m2lnp set to @m2lnp;
 * ncm_stats_dist1d_prepare() prepares @m2lnp.
 *
 * Returns: (transfer full): a new #NcmStatsDist1dSpline
 */
NcmStatsDist1dSpline *
ncm_stats_dist1d_spline_new (NcmSpline *m2lnp)
{
  NcmStatsDist1dSpline *sd1s = g_object_new (NCM_TYPE_STATS_DIST1D_SPLINE,
                                             "m2lnp", m2lnp,
                                             NULL);

  return sd1s;
}

/**
 * ncm_stats_dist1d_spline_new_from_density:
 * @p: a #NcmSpline of the density $p(x)$
 *
 * Creates a new #NcmStatsDist1dSpline with #NcmStatsDist1dSpline:density set to @p;
 * ncm_stats_dist1d_prepare() prepares @p.
 *
 * Returns: (transfer full): a new #NcmStatsDist1dSpline
 */
NcmStatsDist1dSpline *
ncm_stats_dist1d_spline_new_from_density (NcmSpline *p)
{
  NcmStatsDist1dSpline *sd1s = g_object_new (NCM_TYPE_STATS_DIST1D_SPLINE,
                                             "density", p,
                                             NULL);

  return sd1s;
}

