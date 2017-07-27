/***************************************************************************
 *            ncm_stats_dist2d.c
 *
 *  Sat July 22 16:21:25 2017
 *  Copyright  2017  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * ncm_stats_dist2d.c
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
 * SECTION:ncm_stats_dist2d
 * @title: NcmStatsDist2d
 * @short_description: Abstract class for implementing two-dimensional probability distributions
 *
 * Abstract class to reconstruct an arbitrary two-dimensional probability distribution.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_stats_dist2d.h"
#include "math/ncm_spline2d_bicubic.h"

enum
{
  PROP_0,
  PROP_SIZE,
};

G_DEFINE_ABSTRACT_TYPE (NcmStatsDist2d, ncm_stats_dist2d, G_TYPE_OBJECT);

static void
ncm_stats_dist2d_init (NcmStatsDist2d *sd2)
{
	
}

static void
ncm_stats_dist2d_dispose (GObject *object)
{
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist2d_parent_class)->dispose (object);
}

static void
ncm_stats_dist2d_finalize (GObject *object)
{
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist2d_parent_class)->finalize (object);
}

static void
ncm_stats_dist2d_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmStatsDist2d *sd2 = NCM_STATS_DIST2D (object);
  g_return_if_fail (NCM_IS_STATS_DIST2D (object));

	NCM_UNUSED (sd2);
	
  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_stats_dist2d_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmStatsDist2d *sd2 = NCM_STATS_DIST2D (object);
  g_return_if_fail (NCM_IS_STATS_DIST2D (object));

	NCM_UNUSED (sd2);
	
  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_stats_dist2d_class_init (NcmStatsDist2dClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->dispose      = ncm_stats_dist2d_dispose;
  object_class->finalize     = ncm_stats_dist2d_finalize;
  object_class->set_property = ncm_stats_dist2d_set_property;
  object_class->get_property = ncm_stats_dist2d_get_property;

  klass->xbounds  = NULL;
	klass->ybounds  = NULL;
	klass->pdf      = NULL;
	klass->m2lnL    = NULL;
	klass->cdf      = NULL;
  klass->marginal = NULL;
	klass->inv_pdf  = NULL;
	klass->prepare  = NULL;
}

/**
 * ncm_stats_dist2d_ref:
 * @sd2: a #NcmStatsDist2d
 *
 * Increases the reference count of @sd2.
 * 
 * Returns: (transfer full): @sd2.
 */
NcmStatsDist2d *
ncm_stats_dist2d_ref (NcmStatsDist2d *sd2)
{
  return g_object_ref (sd2);
}

/**
 * ncm_stats_dist2d_free:
 * @sd2: a #NcmStatsDist2d
 *
 * Decreases the reference count of @sd2.
 *
 */
void 
ncm_stats_dist2d_free (NcmStatsDist2d *sd2)
{
  g_object_unref (sd2);
}

/**
 * ncm_stats_dist2d_clear:
 * @sd2: a #NcmStatsDist2d
 *
 * Decreases the reference count of *@sd2 and sets the pointer *@sd2 to NULL.
 *
 */
void 
ncm_stats_dist2d_clear (NcmStatsDist2d **sd2)
{
  g_clear_object (sd2);
}

/**
 * ncm_stats_dist2d_prepare:
 * @sd2: a #NcmStatsDist2d
 *
 * Prepares the object for calculations.
 */
void 
ncm_stats_dist2d_prepare (NcmStatsDist2d *sd2)
{
  NcmStatsDist2dClass *sd2_class = NCM_STATS_DIST2D_GET_CLASS (sd2); 

  if (sd2_class->prepare != NULL)
    sd2_class->prepare (sd2);

}

/**
 * ncm_stats_dist2d_eval_pdf:
 * @sd2: a #NcmStatsDist2d
 * @x: random variable value
 * @y: random variable value
 *
 * Calculates the value of the probability density function (PDF) at @x and @y.
 * 
 * Returns: the PDF value at @x and @y
 */
gdouble 
ncm_stats_dist2d_eval_pdf (NcmStatsDist2d *sd2, gdouble x, gdouble y)
{
	return NCM_STATS_DIST2D_GET_CLASS (sd2)->pdf (sd2, x, y);
}

/**
 * ncm_stats_dist2d_eval_m2lnp:
 * @sd2: a #NcmStatsDist2d
 * @x: random variable value
 * @y: random variable value 
 *
 * Calculates the value of the $-2\ln(L(x, y))$ for the probability density function.
 * 
 * Returns: the value of $-2\ln(L(x, y))$.
 */
gdouble 
ncm_stats_dist2d_eval_m2lnp (NcmStatsDist2d *sd2, gdouble x, gdouble y)
{
  return NCM_STATS_DIST2D_GET_CLASS (sd2)->m2lnL (sd2, x, y);
}

/**
 * ncm_stats_dist2d_eval_cdf:
 * @sd2: a #NcmStatsDist2d
 * @x: random variable value
 * @y: random variable value
 *
 * Calculates the value of the cumulative distribution function (CDF) within [x_i, @x] and [y_i, @y].
 * 
 * Returns: the CDF value given the intervals [x_i, @x] and [y_i, @y]
 */
gdouble 
ncm_stats_dist2d_eval_cdf (NcmStatsDist2d *sd2, gdouble x, gdouble y)
{
  return NCM_STATS_DIST2D_GET_CLASS (sd2)->cdf (sd2, x, y);
}

/**
 * ncm_stats_dist2d_eval_inv_pdf:
 * @sd2: a #NcmStatsDist2d
 * @u: a number between [0, 1]
 *
 * FIXME 
 * Calculates the value of the random variable $x$ for which the cumulative
 * distribution satisfy $\int_{x_i}^x\mathrm{d}x^\prime p(x^\prime) = u$.
 * 
 * Returns: the value of y.
 */
gdouble
ncm_stats_dist2d_eval_marginal (NcmStatsDist2d *sd2, gdouble u)
{
	return NCM_STATS_DIST2D_GET_CLASS (sd2)->marginal (sd2, u);
}

/**
 * ncm_stats_dist2d_eval_inv_pdf:
 * @sd2: a #NcmStatsDist2d
 * @u: a number between [0, 1]
 * @y: FIXME
 *
 * FIXME 
 * Calculates the value of the random variable $x$ for which the cumulative
 * distribution satisfy $\int_{x_i}^x\mathrm{d}x^\prime p(x^\prime) = u$.
 * 
 * Returns: the value of x.
 */
gdouble
ncm_stats_dist2d_eval_inv_pdf (NcmStatsDist2d *sd2, gdouble u, gdouble y)
{
	return NCM_STATS_DIST2D_GET_CLASS (sd2)->inv_pdf (sd2, u, y);
}

