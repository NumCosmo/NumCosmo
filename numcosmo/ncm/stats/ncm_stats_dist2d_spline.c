/***************************************************************************
 *            ncm_stats_dist2d_spline.c
 *
 *  Sat July 22 22:31:17 2017
 *  Copyright  2017  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * ncm_stats_dist2d_spline.c
 * Copyright (C) 2017 Mariana Penna Lima <pennalima@gmail.com>
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
 * NcmStatsDist2dSpline:
 *
 * Two-dimensional distribution with $-2\ln p(x, y)$ given by a spline.
 *
 * The support is the knot range of the #NcmStatsDist2dSpline:m2lnp spline, and the density
 * $p(x, y) = e^{-m_2(x, y)/2}$ is not normalized. Only ncm_stats_dist2d_eval_m2lnp(),
 * ncm_stats_dist2d_eval_pdf() and the bounds are implemented; the cumulative distribution,
 * the marginals and the conditional quantiles abort.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/stats/ncm_stats_dist2d_spline.h"

enum
{
  PROP_0,
  PROP_M2LNP,
};

struct _NcmStatsDist2dSpline
{
  /*< private >*/
  NcmStatsDist2d parent_instance;
  NcmSpline2d *m2lnp;
};

G_DEFINE_TYPE (NcmStatsDist2dSpline, ncm_stats_dist2d_spline, NCM_TYPE_STATS_DIST2D)

static void
ncm_stats_dist2d_spline_init (NcmStatsDist2dSpline *sd2s)
{
  sd2s->m2lnp = NULL;
}

static void
_ncm_stats_dist2d_spline_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmStatsDist2dSpline *sd2s = NCM_STATS_DIST2D_SPLINE (object);

  g_return_if_fail (NCM_IS_STATS_DIST2D_SPLINE (object));

  switch (prop_id)
  {
    case PROP_M2LNP:
      ncm_spline2d_clear (&sd2s->m2lnp);
      sd2s->m2lnp = g_value_dup_object (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_stats_dist2d_spline_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmStatsDist2dSpline *sd2s = NCM_STATS_DIST2D_SPLINE (object);

  g_return_if_fail (NCM_IS_STATS_DIST2D_SPLINE (object));

  switch (prop_id)
  {
    case PROP_M2LNP:
      g_value_set_object (value, sd2s->m2lnp);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_stats_dist2d_spline_dispose (GObject *object)
{
  NcmStatsDist2dSpline *sd2s = NCM_STATS_DIST2D_SPLINE (object);

  ncm_spline2d_clear (&sd2s->m2lnp);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist2d_spline_parent_class)->dispose (object);
}

static void _ncm_stats_dist2d_spline_xbounds (NcmStatsDist2d *sd2, gdouble *xi, gdouble *xf);
static void _ncm_stats_dist2d_spline_ybounds (NcmStatsDist2d *sd2, gdouble *yi, gdouble *yf);
static gdouble _ncm_stats_dist2d_spline_pdf (NcmStatsDist2d *sd2, const gdouble x, const gdouble y);
static gdouble _ncm_stats_dist2d_spline_m2lnp (NcmStatsDist2d *sd2, const gdouble x, const gdouble y);
static void _ncm_stats_dist2d_spline_prepare (NcmStatsDist2d *sd2);

static void
ncm_stats_dist2d_spline_class_init (NcmStatsDist2dSplineClass *klass)
{
  GObjectClass *object_class     = G_OBJECT_CLASS (klass);
  NcmStatsDist2dClass *sd2_class = NCM_STATS_DIST2D_CLASS (klass);

  object_class->set_property = &_ncm_stats_dist2d_spline_set_property;
  object_class->get_property = &_ncm_stats_dist2d_spline_get_property;
  object_class->dispose      = &_ncm_stats_dist2d_spline_dispose;

  g_object_class_install_property (object_class,
                                   PROP_M2LNP,
                                   g_param_spec_object ("m2lnp",
                                                        NULL,
                                                        "Spline of -2 ln p",
                                                        NCM_TYPE_SPLINE2D,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  sd2_class->xbounds = &_ncm_stats_dist2d_spline_xbounds;
  sd2_class->ybounds = &_ncm_stats_dist2d_spline_ybounds;
  sd2_class->pdf     = &_ncm_stats_dist2d_spline_pdf;
  sd2_class->m2lnp   = &_ncm_stats_dist2d_spline_m2lnp;
  sd2_class->prepare = &_ncm_stats_dist2d_spline_prepare;
}

static void
_ncm_stats_dist2d_spline_xbounds (NcmStatsDist2d *sd2, gdouble *xi, gdouble *xf)
{
  NcmStatsDist2dSpline *sd2s = NCM_STATS_DIST2D_SPLINE (sd2);
  NcmVector *xv              = ncm_spline2d_peek_xv (sd2s->m2lnp);

  *xi = ncm_vector_get (xv, 0);
  *xf = ncm_vector_get (xv, ncm_vector_len (xv) - 1);
}

static void
_ncm_stats_dist2d_spline_ybounds (NcmStatsDist2d *sd2, gdouble *yi, gdouble *yf)
{
  NcmStatsDist2dSpline *sd2s = NCM_STATS_DIST2D_SPLINE (sd2);
  NcmVector *yv              = ncm_spline2d_peek_yv (sd2s->m2lnp);

  *yi = ncm_vector_get (yv, 0);
  *yf = ncm_vector_get (yv, ncm_vector_len (yv) - 1);
}

static gdouble
_ncm_stats_dist2d_spline_m2lnp (NcmStatsDist2d *sd2, const gdouble x, const gdouble y)
{
  NcmStatsDist2dSpline *sd2s = NCM_STATS_DIST2D_SPLINE (sd2);

  return ncm_spline2d_eval (sd2s->m2lnp, x, y);
}

static gdouble
_ncm_stats_dist2d_spline_pdf (NcmStatsDist2d *sd2, const gdouble x, const gdouble y)
{
  return exp (-0.5 * _ncm_stats_dist2d_spline_m2lnp (sd2, x, y));
}

static void
_ncm_stats_dist2d_spline_prepare (NcmStatsDist2d *sd2)
{
  NcmStatsDist2dSpline *sd2s = NCM_STATS_DIST2D_SPLINE (sd2);

  if (sd2s->m2lnp == NULL)
    g_error ("_ncm_stats_dist2d_spline_prepare: no m2lnp spline set.");

  ncm_spline2d_prepare (sd2s->m2lnp);
}

/**
 * ncm_stats_dist2d_spline_new:
 * @m2lnp: a #NcmSpline2d of $-2\ln p(x, y)$
 *
 * Creates a new #NcmStatsDist2dSpline with #NcmStatsDist2dSpline:m2lnp set to @m2lnp;
 * ncm_stats_dist2d_prepare() prepares @m2lnp.
 *
 * Returns: (transfer full): a new #NcmStatsDist2dSpline
 */
NcmStatsDist2dSpline *
ncm_stats_dist2d_spline_new (NcmSpline2d *m2lnp)
{
  NcmStatsDist2dSpline *sd2s = g_object_new (NCM_TYPE_STATS_DIST2D_SPLINE,
                                             "m2lnp", m2lnp,
                                             NULL);

  return sd2s;
}

