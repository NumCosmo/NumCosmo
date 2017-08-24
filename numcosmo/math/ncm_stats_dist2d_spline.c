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
 * SECTION:ncm_stats_dist2d_spline
 * @title: NcmStatsDist2dSpline
 * @short_description: Two-dimensional probability distribution based on a spline
 * 
 * Reconstruction of an arbitrary two-dimensional probability distribution based on a spline.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_stats_dist2d_spline.h"

enum
{
  PROP_0,
  PROP_M2LNL,
  PROP_MARGINAL_X, 
	PROP_SIZE,
};

G_DEFINE_TYPE (NcmStatsDist2dSpline, ncm_stats_dist2d_spline, NCM_TYPE_STATS_DIST2D);

static void
ncm_stats_dist2d_spline_init (NcmStatsDist2dSpline *sd2s)
{
  sd2s->m2lnp      = NULL;
  sd2s->marginal_x = TRUE; 
}

static void
ncm_stats_dist2d_spline_dispose (GObject *object)
{
  NcmStatsDist2dSpline *sd2s = NCM_STATS_DIST2D_SPLINE (object);
  
  ncm_spline2d_clear (&sd2s->m2lnp);
  
  /* Chain up : end */  
  G_OBJECT_CLASS (ncm_stats_dist2d_spline_parent_class)->dispose (object);
}

static void
ncm_stats_dist2d_spline_finalize (GObject *object)
{

  /* Chain up : end */  
  G_OBJECT_CLASS (ncm_stats_dist2d_spline_parent_class)->finalize (object);
}

static void
ncm_stats_dist2d_spline_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmStatsDist2dSpline *sd2s = NCM_STATS_DIST2D_SPLINE (object);
  g_return_if_fail (NCM_IS_STATS_DIST2D_SPLINE (object));

  switch (prop_id)
  {
    case PROP_M2LNL:
      ncm_spline2d_clear (&sd2s->m2lnp);
      sd2s->m2lnp = g_value_dup_object (value);
      break;
    case PROP_MARGINAL_X:
      sd2s->marginal_x = g_value_get_boolean (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_stats_dist2d_spline_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmStatsDist2dSpline *sd2s = NCM_STATS_DIST2D_SPLINE (object);
  g_return_if_fail (NCM_IS_STATS_DIST2D_SPLINE (object));

  switch (prop_id)
  {
    case PROP_M2LNL:
      g_value_set_object (value, sd2s->m2lnp);
      break;
    case PROP_MARGINAL_X:
      g_value_set_boolean (value, sd2s->marginal_x);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void _ncm_stats_dist2d_spline_xbounds (NcmStatsDist2d *sd2, gdouble *xi, gdouble *xf);
static void _ncm_stats_dist2d_spline_ybounds (NcmStatsDist2d *sd2, gdouble *yi, gdouble *yf);
static gdouble _ncm_stats_dist2d_spline_pdf (NcmStatsDist2d *sd2, gdouble x, gdouble y);
static gdouble _ncm_stats_dist2d_spline_m2lnp (NcmStatsDist2d *sd2, gdouble x, gdouble y);
static gdouble _ncm_stats_dist2d_spline_cdf (NcmStatsDist2d *sd2, gdouble x, gdouble y);
static gdouble _ncm_stats_dist2d_spline_marginal (NcmStatsDist2d *sd2, gdouble u);
static gdouble _ncm_stats_dist2d_spline_inv_pdf (NcmStatsDist2d *sd2, gdouble u, gdouble y);
static void _ncm_stats_dist2d_spline_prepare (NcmStatsDist2d *sd2);

static void
ncm_stats_dist2d_spline_class_init (NcmStatsDist2dSplineClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmStatsDist2dClass *sd2_class = NCM_STATS_DIST2D_CLASS (klass);

  object_class->set_property = ncm_stats_dist2d_spline_set_property;
  object_class->get_property = ncm_stats_dist2d_spline_get_property;
  object_class->dispose      = ncm_stats_dist2d_spline_dispose;
  object_class->finalize     = ncm_stats_dist2d_spline_finalize;

  g_object_class_install_property (object_class,
                                   PROP_M2LNL,
                                   g_param_spec_object ("m2lnp",
                                                        NULL,
                                                        "m2lnp",
                                                        NCM_TYPE_SPLINE2D,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_MARGINAL_X,
                                   g_param_spec_boolean ("marginal-x",
                                                         NULL,
                                                         "Compute marginal with respect to x if True, and y if False.",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	sd2_class->xbounds      = &_ncm_stats_dist2d_spline_xbounds;
	sd2_class->ybounds      = &_ncm_stats_dist2d_spline_ybounds;
  sd2_class->pdf          = &_ncm_stats_dist2d_spline_pdf;
  sd2_class->m2lnp        = &_ncm_stats_dist2d_spline_m2lnp;
	sd2_class->cdf          = &_ncm_stats_dist2d_spline_cdf;
	sd2_class->marginal_pdf = &_ncm_stats_dist2d_spline_marginal;
	sd2_class->inv_cond     = &_ncm_stats_dist2d_spline_inv_pdf;
  sd2_class->prepare      = &_ncm_stats_dist2d_spline_prepare;
}

static void 
_ncm_stats_dist2d_spline_xbounds (NcmStatsDist2d *sd2, gdouble *xi, gdouble *xf)
{
	NcmStatsDist2dSpline *sd2s = NCM_STATS_DIST2D_SPLINE (sd2);
	const gint len = ncm_vector_len (sd2s->m2lnp->xv);
	*xi = ncm_vector_get (sd2s->m2lnp->xv, 0);
	*xf = ncm_vector_get (sd2s->m2lnp->xv, len - 1);
}

static void 
_ncm_stats_dist2d_spline_ybounds (NcmStatsDist2d *sd2, gdouble *yi, gdouble *yf)
{
	NcmStatsDist2dSpline *sd2s = NCM_STATS_DIST2D_SPLINE (sd2);
	const gint len = ncm_vector_len (sd2s->m2lnp->yv);
	*yi = ncm_vector_get (sd2s->m2lnp->yv, 0);
	*yf = ncm_vector_get (sd2s->m2lnp->yv, len - 1);
}

static gdouble 
_ncm_stats_dist2d_spline_m2lnp (NcmStatsDist2d *sd2, gdouble x, gdouble y)
{
  NcmStatsDist2dSpline *sd2s = NCM_STATS_DIST2D_SPLINE (sd2);

  return ncm_spline2d_eval (sd2s->m2lnp, x, y);
}

static gdouble 
_ncm_stats_dist2d_spline_pdf (NcmStatsDist2d *sd2, gdouble x, gdouble y)
{
  const gdouble m2lnp = _ncm_stats_dist2d_spline_m2lnp (sd2, x, y);
  return exp (-0.5 * m2lnp);
}

static gdouble 
_ncm_stats_dist2d_spline_cdf (NcmStatsDist2d *sd2, gdouble x, gdouble y)
{
  NcmStatsDist2dSpline *sd2s = NCM_STATS_DIST2D_SPLINE (sd2);
	gdouble xi, yi, xf, yf;
  	
  _ncm_stats_dist2d_spline_xbounds (sd2, &xi, &xf);
	_ncm_stats_dist2d_spline_ybounds (sd2, &yi, &yf);
	
  return ncm_spline2d_integ_dxdy (sd2s->m2lnp, xi, x, yi, y);
}

static gdouble 
_ncm_stats_dist2d_spline_marginal (NcmStatsDist2d *sd2, gdouble u)
{
  NcmStatsDist2dSpline *sd2s = NCM_STATS_DIST2D_SPLINE (sd2);
	NCM_UNUSED (sd2s);	
  return 0.0;
}

static gdouble 
_ncm_stats_dist2d_spline_inv_pdf (NcmStatsDist2d *sd2, gdouble u, gdouble y)
{
  NcmStatsDist2dSpline *sd2s = NCM_STATS_DIST2D_SPLINE (sd2);
	NCM_UNUSED (sd2s);	
  return 0.0;
}

static void 
_ncm_stats_dist2d_spline_prepare (NcmStatsDist2d *sd2)
{
  NcmStatsDist2dSpline *sd2s = NCM_STATS_DIST2D_SPLINE (sd2);

	ncm_spline2d_prepare (sd2s->m2lnp);
}

/**
 * ncm_stats_dist2d_spline_new:
 * @m2lnp: a #NcmSpline2d
 * 
 * Returns a new #NcmStatsDist2dSpline where @m2lnp, $-2\ln(L(x, y))$, is a #NcmSpline2d, where $L(x, y)$ 
 * is the probability density function.
 * 
 * Returns: a new #NcmStatsDist2dSpline
 */
NcmStatsDist2dSpline *
ncm_stats_dist2d_spline_new (NcmSpline2d *m2lnp)
{
  NcmStatsDist2dSpline *sd2s = g_object_new (NCM_TYPE_STATS_DIST2D_SPLINE,
                                             "m2lnp", m2lnp,
                                             NULL);
  return sd2s;
}
