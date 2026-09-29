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
 * NcmStatsDist2d:
 *
 * Base class for two-dimensional probability distributions.
 *
 * A subclass implements the methods it supports; the others abort naming the subclass.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/stats/ncm_stats_dist2d.h"

G_DEFINE_ABSTRACT_TYPE (NcmStatsDist2d, ncm_stats_dist2d, G_TYPE_OBJECT)

static void
ncm_stats_dist2d_init (NcmStatsDist2d *sd2)
{
}

#define _NCM_STATS_DIST2D_NOT_IMPLEMENTED(name) \
        g_error ("ncm_stats_dist2d_" name ": `%s' does not implement " name ".", G_OBJECT_TYPE_NAME (sd2))

static void
_ncm_stats_dist2d_xbounds (NcmStatsDist2d *sd2, gdouble *xi, gdouble *xf)
{
  _NCM_STATS_DIST2D_NOT_IMPLEMENTED ("xbounds");
}

static void
_ncm_stats_dist2d_ybounds (NcmStatsDist2d *sd2, gdouble *yi, gdouble *yf)
{
  _NCM_STATS_DIST2D_NOT_IMPLEMENTED ("ybounds");
}

static gdouble
_ncm_stats_dist2d_pdf (NcmStatsDist2d *sd2, const gdouble x, const gdouble y)
{
  _NCM_STATS_DIST2D_NOT_IMPLEMENTED ("pdf");

  return 0.0;
}

static gdouble
_ncm_stats_dist2d_m2lnp (NcmStatsDist2d *sd2, const gdouble x, const gdouble y)
{
  _NCM_STATS_DIST2D_NOT_IMPLEMENTED ("m2lnp");

  return 0.0;
}

static gdouble
_ncm_stats_dist2d_cdf (NcmStatsDist2d *sd2, const gdouble x, const gdouble y)
{
  _NCM_STATS_DIST2D_NOT_IMPLEMENTED ("cdf");

  return 0.0;
}

static gdouble
_ncm_stats_dist2d_marginal_pdf (NcmStatsDist2d *sd2, const gdouble xy)
{
  _NCM_STATS_DIST2D_NOT_IMPLEMENTED ("marginal_pdf");

  return 0.0;
}

static gdouble
_ncm_stats_dist2d_marginal_cdf (NcmStatsDist2d *sd2, const gdouble xy)
{
  _NCM_STATS_DIST2D_NOT_IMPLEMENTED ("marginal_cdf");

  return 0.0;
}

static gdouble
_ncm_stats_dist2d_marginal_inv_cdf (NcmStatsDist2d *sd2, const gdouble u)
{
  _NCM_STATS_DIST2D_NOT_IMPLEMENTED ("marginal_inv_cdf");

  return 0.0;
}

static gdouble
_ncm_stats_dist2d_inv_cond (NcmStatsDist2d *sd2, const gdouble u, const gdouble xy)
{
  _NCM_STATS_DIST2D_NOT_IMPLEMENTED ("inv_cond");

  return 0.0;
}

static void
ncm_stats_dist2d_class_init (NcmStatsDist2dClass *klass)
{
  klass->xbounds          = &_ncm_stats_dist2d_xbounds;
  klass->ybounds          = &_ncm_stats_dist2d_ybounds;
  klass->pdf              = &_ncm_stats_dist2d_pdf;
  klass->m2lnp            = &_ncm_stats_dist2d_m2lnp;
  klass->cdf              = &_ncm_stats_dist2d_cdf;
  klass->marginal_pdf     = &_ncm_stats_dist2d_marginal_pdf;
  klass->marginal_cdf     = &_ncm_stats_dist2d_marginal_cdf;
  klass->marginal_inv_cdf = &_ncm_stats_dist2d_marginal_inv_cdf;
  klass->inv_cond         = &_ncm_stats_dist2d_inv_cond;
  klass->prepare          = NULL;
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
 * Decreases the reference count of *@sd2 and sets the pointer *@sd2 to %NULL.
 */
void
ncm_stats_dist2d_clear (NcmStatsDist2d **sd2)
{
  g_clear_object (sd2);
}

/**
 * ncm_stats_dist2d_prepare: (virtual prepare)
 * @sd2: a #NcmStatsDist2d
 *
 * Calls the subclass prepare, if any; must be called before evaluating @sd2.
 */
void
ncm_stats_dist2d_prepare (NcmStatsDist2d *sd2)
{
  NcmStatsDist2dClass *sd2_class = NCM_STATS_DIST2D_GET_CLASS (sd2);

  if (sd2_class->prepare != NULL)
    sd2_class->prepare (sd2);
}

/**
 * ncm_stats_dist2d_xbounds: (virtual xbounds)
 * @sd2: a #NcmStatsDist2d
 * @xi: (out): lower bound of $x$
 * @xf: (out): upper bound of $x$
 *
 * Gets the range of $x$ of the support.
 */
void
ncm_stats_dist2d_xbounds (NcmStatsDist2d *sd2, gdouble *xi, gdouble *xf)
{
  NCM_STATS_DIST2D_GET_CLASS (sd2)->xbounds (sd2, xi, xf);
}

/**
 * ncm_stats_dist2d_ybounds: (virtual ybounds)
 * @sd2: a #NcmStatsDist2d
 * @yi: (out): lower bound of $y$
 * @yf: (out): upper bound of $y$
 *
 * Gets the range of $y$ of the support.
 */
void
ncm_stats_dist2d_ybounds (NcmStatsDist2d *sd2, gdouble *yi, gdouble *yf)
{
  NCM_STATS_DIST2D_GET_CLASS (sd2)->ybounds (sd2, yi, yf);
}

/**
 * ncm_stats_dist2d_eval_pdf: (virtual pdf)
 * @sd2: a #NcmStatsDist2d
 * @x: first variable
 * @y: second variable
 *
 * Evaluates the density at (@x, @y).
 *
 * Returns: the density $p(x, y)$.
 */
gdouble
ncm_stats_dist2d_eval_pdf (NcmStatsDist2d *sd2, const gdouble x, const gdouble y)
{
  return NCM_STATS_DIST2D_GET_CLASS (sd2)->pdf (sd2, x, y);
}

/**
 * ncm_stats_dist2d_eval_m2lnp: (virtual m2lnp)
 * @sd2: a #NcmStatsDist2d
 * @x: first variable
 * @y: second variable
 *
 * Evaluates $-2\ln p(x, y)$.
 *
 * Returns: $-2\ln p(x, y)$.
 */
gdouble
ncm_stats_dist2d_eval_m2lnp (NcmStatsDist2d *sd2, const gdouble x, const gdouble y)
{
  return NCM_STATS_DIST2D_GET_CLASS (sd2)->m2lnp (sd2, x, y);
}

/**
 * ncm_stats_dist2d_eval_cdf: (virtual cdf)
 * @sd2: a #NcmStatsDist2d
 * @x: first variable
 * @y: second variable
 *
 * Evaluates the probability of $[x_i, x] \times [y_i, y]$.
 *
 * Returns: the cumulative distribution at (@x, @y).
 */
gdouble
ncm_stats_dist2d_eval_cdf (NcmStatsDist2d *sd2, const gdouble x, const gdouble y)
{
  return NCM_STATS_DIST2D_GET_CLASS (sd2)->cdf (sd2, x, y);
}

/**
 * ncm_stats_dist2d_eval_marginal_pdf: (virtual marginal_pdf)
 * @sd2: a #NcmStatsDist2d
 * @xy: value of the marginal's variable
 *
 * Evaluates the marginal density of one variable; which one depends on the subclass.
 *
 * Returns: the marginal density at @xy.
 */
gdouble
ncm_stats_dist2d_eval_marginal_pdf (NcmStatsDist2d *sd2, const gdouble xy)
{
  return NCM_STATS_DIST2D_GET_CLASS (sd2)->marginal_pdf (sd2, xy);
}

/**
 * ncm_stats_dist2d_eval_marginal_cdf: (virtual marginal_cdf)
 * @sd2: a #NcmStatsDist2d
 * @xy: value of the marginal's variable
 *
 * Evaluates the cumulative distribution of the marginal of ncm_stats_dist2d_eval_marginal_pdf().
 *
 * Returns: the marginal cumulative distribution at @xy.
 */
gdouble
ncm_stats_dist2d_eval_marginal_cdf (NcmStatsDist2d *sd2, const gdouble xy)
{
  return NCM_STATS_DIST2D_GET_CLASS (sd2)->marginal_cdf (sd2, xy);
}

/**
 * ncm_stats_dist2d_eval_marginal_inv_cdf: (virtual marginal_inv_cdf)
 * @sd2: a #NcmStatsDist2d
 * @u: probability, in $[0, 1]$
 *
 * Evaluates the inverse of ncm_stats_dist2d_eval_marginal_cdf().
 *
 * Returns: the quantile of the marginal.
 */
gdouble
ncm_stats_dist2d_eval_marginal_inv_cdf (NcmStatsDist2d *sd2, const gdouble u)
{
  return NCM_STATS_DIST2D_GET_CLASS (sd2)->marginal_inv_cdf (sd2, u);
}

/**
 * ncm_stats_dist2d_eval_inv_cond: (virtual inv_cond)
 * @sd2: a #NcmStatsDist2d
 * @u: probability, in $[0, 1]$
 * @xy: value of the marginal's variable
 *
 * Evaluates the quantile @u of the other variable conditional on the marginal's variable
 * being @xy.
 *
 * Returns: the conditional quantile.
 */
gdouble
ncm_stats_dist2d_eval_inv_cond (NcmStatsDist2d *sd2, const gdouble u, const gdouble xy)
{
  return NCM_STATS_DIST2D_GET_CLASS (sd2)->inv_cond (sd2, u, xy);
}

