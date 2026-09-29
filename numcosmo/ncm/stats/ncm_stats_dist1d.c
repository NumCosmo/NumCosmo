/***************************************************************************
 *            ncm_stats_dist1d.c
 *
 *  Thu February 12 15:37:11 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_dist1d.c
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
 * NcmStatsDist1d:
 *
 * Base class for one-dimensional probability distributions on $[x_i, x_f]$.
 *
 * A subclass provides the density $p(x)$, which need not be normalized, and $-2\ln p(x)$.
 * With #NcmStatsDist1d:compute-cdf, ncm_stats_dist1d_prepare() integrates the cumulative
 * distribution as an ODE to #NcmStatsDist1d:reltol, in at least 1000 steps, and in units of
 * $p(x_\mathrm{mode})\,(x_f - x_i)$, so its accuracy depends neither on the scale of $p$ nor
 * on that of $x$. The inverse is a Steffen spline through the same knots with $x$ and the
 * probability swapped; it is monotone and stays in $[x_i, x_f]$. Setting $x_i = x_f$ gives a point mass at
 * $x_i$. The mode search uses an internal minimizer, so one object must not be used from
 * several threads at once.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/stats/ncm_stats_dist1d.h"
#include "ncm/spline/ncm_spline_cubic_notaknot.h"
#include "ncm/spline/ncm_spline_gsl.h"
#include "ncm/core/ncm_cfg.h"

enum
{
  PROP_0,
  PROP_XI,
  PROP_XF,
  PROP_NORMA,
  PROP_RELTOL,
  PROP_ABSTOL,
  PROP_COMPUTE_CDF,
};

typedef struct _NcmStatsDist1dPrivate
{
  gdouble xi;
  gdouble xf;
  gdouble norma;
  gdouble p_mode;
  gdouble cdf_norma;
  gdouble reltol;
  gdouble abstol;
  gboolean compute_cdf;
  NcmSpline *inv_cdf;
  GArray *inv_u;
  GArray *inv_x;
  NcmOdeSpline *pdf;
  gsl_min_fminimizer *fmin;
} NcmStatsDist1dPrivate;

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmStatsDist1d, ncm_stats_dist1d, G_TYPE_OBJECT)

/* The CDF integration steps are at most (xf - xi) / NCM_STATS_DIST1D_MIN_SUBDIVISIONS, so no feature wider than that is stepped over */
#define NCM_STATS_DIST1D_MIN_SUBDIVISIONS 1000

static gdouble
_ncm_stats_dist1d_pdf_dydx (gdouble y, gdouble x, gpointer userdata)
{
  NcmStatsDist1d *sd1         = NCM_STATS_DIST1D (userdata);
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);

  return NCM_STATS_DIST1D_GET_CLASS (sd1)->p (sd1, x) / self->p_mode / (self->xf - self->xi);
}

static void
ncm_stats_dist1d_init (NcmStatsDist1d *sd1)
{
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);
  NcmSpline *s1               = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());

  self->xi          = 0.0;
  self->xf          = 0.0;
  self->norma       = 0.0;
  self->p_mode      = 1.0;
  self->cdf_norma   = 1.0;
  self->reltol      = 0.0;
  self->inv_cdf     = NCM_SPLINE (ncm_spline_gsl_new (gsl_interp_steffen));
  self->inv_u       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->inv_x       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->pdf         = ncm_ode_spline_new (s1, _ncm_stats_dist1d_pdf_dydx);
  self->fmin        = gsl_min_fminimizer_alloc (gsl_min_fminimizer_brent);
  self->compute_cdf = FALSE;

  ncm_ode_spline_set_min_subdivisions (self->pdf, NCM_STATS_DIST1D_MIN_SUBDIVISIONS);

  ncm_spline_free (s1);
}

static void
ncm_stats_dist1d_dispose (GObject *object)
{
  NcmStatsDist1d *sd1         = NCM_STATS_DIST1D (object);
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);

  ncm_spline_clear (&self->inv_cdf);
  ncm_ode_spline_clear (&self->pdf);
  g_clear_pointer (&self->inv_u, g_array_unref);
  g_clear_pointer (&self->inv_x, g_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist1d_parent_class)->dispose (object);
}

static void
ncm_stats_dist1d_finalize (GObject *object)
{
  NcmStatsDist1d *sd1         = NCM_STATS_DIST1D (object);
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);

  gsl_min_fminimizer_free (self->fmin);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist1d_parent_class)->finalize (object);
}

static void
ncm_stats_dist1d_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmStatsDist1d *sd1         = NCM_STATS_DIST1D (object);
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);

  g_return_if_fail (NCM_IS_STATS_DIST1D (object));

  switch (prop_id)
  {
    case PROP_XI:
      self->xi = g_value_get_double (value);
      break;
    case PROP_XF:
      self->xf = g_value_get_double (value);
      break;
    case PROP_RELTOL:
      self->reltol = g_value_get_double (value);
      break;
    case PROP_ABSTOL:
      self->abstol = g_value_get_double (value);
      break;
    case PROP_COMPUTE_CDF:
      ncm_stats_dist1d_set_compute_cdf (sd1, g_value_get_boolean (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_stats_dist1d_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmStatsDist1d *sd1         = NCM_STATS_DIST1D (object);
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);

  g_return_if_fail (NCM_IS_STATS_DIST1D (object));

  switch (prop_id)
  {
    case PROP_XI:
      g_value_set_double (value, ncm_stats_dist1d_get_xi (sd1));
      break;
    case PROP_XF:
      g_value_set_double (value, ncm_stats_dist1d_get_xf (sd1));
      break;
    case PROP_NORMA:
      g_value_set_double (value, self->norma);
      break;
    case PROP_RELTOL:
      g_value_set_double (value, self->reltol);
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, self->abstol);
      break;
    case PROP_COMPUTE_CDF:
      g_value_set_boolean (value, ncm_stats_dist1d_get_compute_cdf (sd1));
      break;

    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static gdouble _ncm_stats_dist1d_p_not_implemented (NcmStatsDist1d *sd1, gdouble x);
static gdouble _ncm_stats_dist1d_m2lnp_not_implemented (NcmStatsDist1d *sd1, gdouble x);
static gdouble _ncm_stats_dist1d_get_current_h_not_implemented (NcmStatsDist1d *sd1);

static void _ncm_stats_dist1d_prepare_inv_cdf (NcmStatsDist1dPrivate *self);

static void
ncm_stats_dist1d_class_init (NcmStatsDist1dClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->dispose      = ncm_stats_dist1d_dispose;
  object_class->finalize     = ncm_stats_dist1d_finalize;
  object_class->set_property = ncm_stats_dist1d_set_property;
  object_class->get_property = ncm_stats_dist1d_get_property;

  g_object_class_install_property (object_class,
                                   PROP_XI,
                                   g_param_spec_double ("xi",
                                                        NULL,
                                                        "x_i",
                                                        -G_MAXDOUBLE, +G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_XF,
                                   g_param_spec_double ("xf",
                                                        NULL,
                                                        "x_f",
                                                        -G_MAXDOUBLE, +G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NORMA,
                                   g_param_spec_double ("norma",
                                                        NULL,
                                                        "Distribution norma",
                                                        0.0, +G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "relative tolerance",
                                                        0.0, 1.0, 1.0e-14,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_ABSTOL,
                                   g_param_spec_double ("abstol",
                                                        NULL,
                                                        "Absolute tolerance on the location of the mode",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_COMPUTE_CDF,
                                   g_param_spec_boolean ("compute-cdf",
                                                         NULL,
                                                         "Whether to compute CDF and inverse CDF",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->p             = &_ncm_stats_dist1d_p_not_implemented;
  klass->m2lnp         = &_ncm_stats_dist1d_m2lnp_not_implemented;
  klass->prepare       = NULL;
  klass->get_current_h = &_ncm_stats_dist1d_get_current_h_not_implemented;
}

static gdouble
_ncm_stats_dist1d_p_not_implemented (NcmStatsDist1d *sd1, gdouble x)
{
  g_error ("_ncm_stats_dist1d_p: `%s' does not implement p.", G_OBJECT_TYPE_NAME (sd1));

  return 0.0;
}

static gdouble
_ncm_stats_dist1d_m2lnp_not_implemented (NcmStatsDist1d *sd1, gdouble x)
{
  g_error ("_ncm_stats_dist1d_m2lnp: `%s' does not implement m2lnp.", G_OBJECT_TYPE_NAME (sd1));

  return 0.0;
}

static gdouble
_ncm_stats_dist1d_get_current_h_not_implemented (NcmStatsDist1d *sd1)
{
  g_error ("_ncm_stats_dist1d_get_current_h: `%s' does not implement get_current_h.", G_OBJECT_TYPE_NAME (sd1));

  return 0.0;
}

/**
 * ncm_stats_dist1d_ref:
 * @sd1: a #NcmStatsDist1d
 *
 * Increases the reference count of @sd1.
 *
 * Returns: (transfer full): @sd1.
 */
NcmStatsDist1d *
ncm_stats_dist1d_ref (NcmStatsDist1d *sd1)
{
  return g_object_ref (sd1);
}

/**
 * ncm_stats_dist1d_free:
 * @sd1: a #NcmStatsDist1d
 *
 * Decreases the reference count of @sd1.
 */
void
ncm_stats_dist1d_free (NcmStatsDist1d *sd1)
{
  g_object_unref (sd1);
}

/**
 * ncm_stats_dist1d_clear:
 * @sd1: a #NcmStatsDist1d
 *
 * Decreases the reference count of *@sd1 and sets the pointer *@sd1 to %NULL.
 */
void
ncm_stats_dist1d_clear (NcmStatsDist1d **sd1)
{
  g_clear_object (sd1);
}

/**
 * ncm_stats_dist1d_prepare:
 * @sd1: a #NcmStatsDist1d
 *
 * Calls the subclass prepare and then, when $x_i \neq x_f$ and #NcmStatsDist1d:compute-cdf is
 * %TRUE, locates the mode, integrates the cumulative distribution and the normalization, and
 * builds the inverse. Must be called after changing $x_i$, $x_f$ or the density. Aborts if
 * $x_f < x_i$ or if the density at the mode is not positive and finite.
 */
void
ncm_stats_dist1d_prepare (NcmStatsDist1d *sd1)
{
  NcmStatsDist1dClass *sd1_class = NCM_STATS_DIST1D_GET_CLASS (sd1);
  NcmStatsDist1dPrivate *self    = ncm_stats_dist1d_get_instance_private (sd1);

  if (sd1_class->prepare != NULL)
    sd1_class->prepare (sd1);

  if (self->xf < self->xi)
    g_error ("ncm_stats_dist1d_prepare: `%s' has xf = % 22.15g below xi = % 22.15g.",
             G_OBJECT_TYPE_NAME (sd1), self->xf, self->xi);

  self->p_mode    = 1.0;
  self->cdf_norma = 1.0;
  self->norma     = 1.0;

  if (G_LIKELY (self->xi != self->xf) && self->compute_cdf)
  {
    /* The CDF is integrated in units of p(mode) (xf - xi), so its tolerances depend neither on the scale of p nor on that of x */
    self->p_mode = sd1_class->p (sd1, ncm_stats_dist1d_eval_mode (sd1));

    if (!(gsl_finite (self->p_mode) && (self->p_mode > 0.0)))
      g_error ("ncm_stats_dist1d_prepare: `%s' has density % 22.15g at its mode.",
               G_OBJECT_TYPE_NAME (sd1), self->p_mode);

    ncm_ode_spline_set_reltol (self->pdf, self->reltol);
    ncm_ode_spline_set_abstol (self->pdf, GSL_DBL_EPSILON * 10.0);
    ncm_ode_spline_set_interval (self->pdf, 0.0, self->xi, self->xf);

    ncm_ode_spline_prepare (self->pdf, sd1);
    self->cdf_norma = ncm_spline_eval (ncm_ode_spline_peek_spline (self->pdf), self->xf);
    self->norma     = self->cdf_norma * self->p_mode * (self->xf - self->xi);

    _ncm_stats_dist1d_prepare_inv_cdf (self);
  }
}

/*
 * The inverse interpolates the CDF knots with x and u swapped. A run of knots with the
 * same u (zero density) keeps its first knot and, at the next representable u, its last,
 * so a u above the run starts at the end of the zero-density interval.
 */
static void
_ncm_stats_dist1d_prepare_inv_cdf (NcmStatsDist1dPrivate *self)
{
  NcmSpline *cdf    = ncm_ode_spline_peek_spline (self->pdf);
  NcmVector *xv     = ncm_spline_peek_xv (cdf);
  NcmVector *yv     = ncm_spline_peek_yv (cdf);
  const guint len   = ncm_vector_len (xv);
  gdouble last_u    = -1.0;
  gdouble run_end_x = 0.0;
  gboolean in_run   = FALSE;
  guint i;

  g_array_set_size (self->inv_u, 0);
  g_array_set_size (self->inv_x, 0);

  for (i = 0; i < len; i++)
  {
    const gdouble u_i = (i + 1 == len) ? 1.0 : GSL_MIN (ncm_vector_get (yv, i) / self->cdf_norma, 1.0);
    const gdouble x_i = ncm_vector_get (xv, i);

    if (u_i > last_u)
    {
      if (in_run)
      {
        const gdouble u_end = nextafter (last_u, 2.0);

        if (u_end < u_i)
        {
          g_array_append_val (self->inv_u, u_end);
          g_array_append_val (self->inv_x, run_end_x);
        }

        in_run = FALSE;
      }

      g_array_append_val (self->inv_u, u_i);
      g_array_append_val (self->inv_x, x_i);
      last_u = u_i;
    }
    else
    {
      run_end_x = x_i;
      in_run    = TRUE;
    }
  }

  ncm_spline_set_array (self->inv_cdf, self->inv_u, self->inv_x, TRUE);
}

/**
 * ncm_stats_dist1d_set_xi:
 * @sd1: a #NcmStatsDist1d
 * @xi: lower bound $x_i$
 *
 * Sets #NcmStatsDist1d:xi.
 */
void
ncm_stats_dist1d_set_xi (NcmStatsDist1d *sd1, gdouble xi)
{
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);

  self->xi = xi;
}

/**
 * ncm_stats_dist1d_set_xf:
 * @sd1: a #NcmStatsDist1d
 * @xf: upper bound $x_f$
 *
 * Sets #NcmStatsDist1d:xf.
 */
void
ncm_stats_dist1d_set_xf (NcmStatsDist1d *sd1, gdouble xf)
{
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);

  self->xf = xf;
}

/**
 * ncm_stats_dist1d_get_xi:
 * @sd1: a #NcmStatsDist1d
 *
 * Returns: the lower bound $x_i$.
 */
gdouble
ncm_stats_dist1d_get_xi (NcmStatsDist1d *sd1)
{
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);

  return self->xi;
}

/**
 * ncm_stats_dist1d_get_xf:
 * @sd1: a #NcmStatsDist1d
 *
 * Returns: the upper bound $x_f$.
 */
gdouble
ncm_stats_dist1d_get_xf (NcmStatsDist1d *sd1)
{
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);

  return self->xf;
}

/**
 * ncm_stats_dist1d_get_current_h: (virtual get_current_h)
 * @sd1: a #NcmStatsDist1d
 *
 * Gets the kernel bandwidth of a kernel density estimate. Aborts for a subclass that
 * does not implement it.
 *
 * Returns: the current bandwidth $h$.
 */
gdouble
ncm_stats_dist1d_get_current_h (NcmStatsDist1d *sd1)
{
  NcmStatsDist1dClass *sd1_class = NCM_STATS_DIST1D_GET_CLASS (sd1);

  return sd1_class->get_current_h (sd1);
}

/**
 * ncm_stats_dist1d_set_compute_cdf:
 * @sd1: a #NcmStatsDist1d
 * @compute_cdf: whether to compute the cumulative distribution
 *
 * Sets #NcmStatsDist1d:compute-cdf. Without it ncm_stats_dist1d_prepare() computes neither
 * the normalization nor the cumulative distribution and its inverse.
 */
void
ncm_stats_dist1d_set_compute_cdf (NcmStatsDist1d *sd1, gboolean compute_cdf)
{
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);

  self->compute_cdf = compute_cdf;
}

/**
 * ncm_stats_dist1d_get_compute_cdf:
 * @sd1: a #NcmStatsDist1d
 *
 * Returns: %TRUE if @sd1 computes the cumulative distribution and its inverse.
 */
gboolean
ncm_stats_dist1d_get_compute_cdf (NcmStatsDist1d *sd1)
{
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);

  return self->compute_cdf;
}

/**
 * ncm_stats_dist1d_eval_p:
 * @sd1: a #NcmStatsDist1d
 * @x: random variable value
 *
 * Evaluates the density at @x divided by the normalization. Without
 * #NcmStatsDist1d:compute-cdf the normalization is 1 and the density is not normalized.
 *
 * Returns: the density $p(x)$.
 */
gdouble
ncm_stats_dist1d_eval_p (NcmStatsDist1d *sd1, gdouble x)
{
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);

  if (G_UNLIKELY (self->xi == self->xf))
    return self->xi == x ? 1.0 : 0.0;

  return NCM_STATS_DIST1D_GET_CLASS (sd1)->p (sd1, x) / self->norma;
}

/**
 * ncm_stats_dist1d_eval_m2lnp:
 * @sd1: a #NcmStatsDist1d
 * @x: random variable value
 *
 * Evaluates $-2\ln p(x)$ of the density as given by the subclass, without the
 * normalization, see ncm_stats_dist1d_eval_norma().
 *
 * Returns: $-2\ln p(x)$.
 */
gdouble
ncm_stats_dist1d_eval_m2lnp (NcmStatsDist1d *sd1, gdouble x)
{
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);

  if (G_UNLIKELY (self->xi == self->xf))
    return self->xi == x ? 0.0 : GSL_POSINF;

  return NCM_STATS_DIST1D_GET_CLASS (sd1)->m2lnp (sd1, x);
}

/**
 * ncm_stats_dist1d_eval_pdf:
 * @sd1: a #NcmStatsDist1d
 * @x: random variable value, in $[x_i, x_f]$
 *
 * Evaluates the cumulative distribution $\int_{x_i}^x p(x^\prime)\,\mathrm{d}x^\prime$.
 * Requires #NcmStatsDist1d:compute-cdf; @x is not checked.
 *
 * Returns: the probability of $[x_i, x]$.
 */
gdouble
ncm_stats_dist1d_eval_pdf (NcmStatsDist1d *sd1, gdouble x)
{
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);

  if (G_UNLIKELY (self->xi == self->xf))
    return self->xi <= x ? 1.0 : 0.0;

  return ncm_spline_eval (ncm_ode_spline_peek_spline (self->pdf), x) / self->cdf_norma;
}

/**
 * ncm_stats_dist1d_eval_norma:
 * @sd1: a #NcmStatsDist1d
 *
 * Gets the integral of the subclass density over $[x_i, x_f]$, computed by
 * ncm_stats_dist1d_prepare(); 1 without #NcmStatsDist1d:compute-cdf.
 *
 * Returns: the normalization.
 */
gdouble
ncm_stats_dist1d_eval_norma (NcmStatsDist1d *sd1)
{
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);

  if (G_UNLIKELY (self->xi == self->xf))
    return 1.0;

  return self->norma;
}

/**
 * ncm_stats_dist1d_eval_inv_pdf:
 * @sd1: a #NcmStatsDist1d
 * @u: probability, in $[0, 1]$
 *
 * Evaluates the inverse of the cumulative distribution, the $x$ with
 * $\int_{x_i}^x p(x^\prime)\,\mathrm{d}x^\prime = u$. Returns $x_i$ for $u \leq 0$ and $x_f$ for
 * $u \geq 1$. Requires #NcmStatsDist1d:compute-cdf.
 *
 * Returns: the quantile $x$.
 */
gdouble
ncm_stats_dist1d_eval_inv_pdf (NcmStatsDist1d *sd1, const gdouble u)
{
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);

  if (G_UNLIKELY (self->xi == self->xf))
    return self->xi;
  else if (G_UNLIKELY (u <= 0.0))
    return self->xi;
  else if (G_UNLIKELY (u >= 1.0))
    return self->xf;
  else
    return ncm_spline_eval (self->inv_cdf, u);
}

/**
 * ncm_stats_dist1d_eval_inv_pdf_tail:
 * @sd1: a #NcmStatsDist1d
 * @v: probability, in $[0, 1]$
 *
 * Evaluates the $x$ with $\int_x^{x_f} p(x^\prime)\,\mathrm{d}x^\prime = v$, that is
 * ncm_stats_dist1d_eval_inv_pdf() at $1 - v$, so @v below $\epsilon$ is not resolved.
 * Returns $x_f$ for $v \leq 0$ and $x_i$ for $v \geq 1$. Requires #NcmStatsDist1d:compute-cdf.
 *
 * Returns: the quantile $x$.
 */
gdouble
ncm_stats_dist1d_eval_inv_pdf_tail (NcmStatsDist1d *sd1, const gdouble v)
{
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);

  if (G_UNLIKELY (self->xi == self->xf))
    return self->xi;
  else if (G_UNLIKELY (v <= 0.0))
    return self->xf;
  else if (G_UNLIKELY (v >= 1.0))
    return self->xi;
  else
    return ncm_spline_eval (self->inv_cdf, 1.0 - v);
}

/**
 * ncm_stats_dist1d_gen:
 * @sd1: a #NcmStatsDist1d
 * @rng: a #NcmRNG
 *
 * Draws a value from the distribution by inverting the cumulative distribution at a
 * uniform deviate. Requires #NcmStatsDist1d:compute-cdf.
 *
 * Returns: the drawn value.
 */
gdouble
ncm_stats_dist1d_gen (NcmStatsDist1d *sd1, NcmRNG *rng)
{
  const gdouble u = ncm_rng_uniform01_gen (rng);

  return ncm_stats_dist1d_eval_inv_pdf (sd1, u);
}

static gdouble
_ncm_stats_dist1d_m2lnp (gdouble x, gpointer p)
{
  NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (p);

  return ncm_stats_dist1d_eval_m2lnp (sd1, x);
}

/**
 * ncm_stats_dist1d_eval_mode:
 * @sd1: a #NcmStatsDist1d
 *
 * Locates the maximum of the density: the minimum of $-2\ln p$ on 1000 equally spaced
 * points, refined by Brent's method between the two neighbouring grid points to a relative
 * tolerance $\sqrt{\mathrm{reltol}}$ and the absolute tolerance #NcmStatsDist1d:abstol.
 * When the two best grid points tie, Brent's method refines between them. The grid point is
 * returned when it is $x_i$ or $x_f$, or when a neighbour has zero density.
 * Warns if the refinement stops before its tolerance.
 *
 * Returns: the mode.
 */
gdouble
ncm_stats_dist1d_eval_mode (NcmStatsDist1d *sd1)
{
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);
  const gdouble reltol        = sqrt (self->reltol);
  const gint max_iter         = 1000000;
  const gint linear_search    = 1000;
  const gdouble dx            = (self->xf - self->xi) / (linear_search - 1.0);
  gdouble x                   = 0.5 * (self->xf + self->xi);
  gint k_min                  = -1;
  gint iter                   = 0;
  gdouble x0, x1, last_x0, last_x1;
  gsl_function F;
  gdouble fmin;
  gint status;
  gint ret;

  if (G_UNLIKELY (self->xi == self->xf))
    return self->xi;

  F.params   = sd1;
  F.function = &_ncm_stats_dist1d_m2lnp;

  fmin = GSL_POSINF;

  for (iter = 0; iter < linear_search; iter++)
  {
    const gdouble x_try = self->xi + dx * iter;
    const gdouble f_try = ncm_stats_dist1d_eval_m2lnp (sd1, x_try);

    if (f_try < fmin)
    {
      fmin  = f_try;
      x     = x_try;
      k_min = iter;
    }
  }

  if ((k_min <= 0) || (k_min >= linear_search - 1))
    return x;

  /* Brent refines inside the grid neighbours of the minimum, and needs f(x) finite and below both */
  x0      = self->xi + dx * (k_min - 1);
  x1      = self->xi + dx * (k_min + 1);
  last_x0 = x0;
  last_x1 = x1;

  {
    const gdouble f_x0 = ncm_stats_dist1d_eval_m2lnp (sd1, x0);
    const gdouble f_x1 = ncm_stats_dist1d_eval_m2lnp (sd1, x1);

    if (!(gsl_finite (f_x0) && gsl_finite (f_x1)))
      return x;

    /* A tie with the right neighbour (a mode between two grid points) brackets between them */
    if (f_x1 == fmin)
    {
      const gdouble x_mid = x + 0.5 * dx;

      if (!(ncm_stats_dist1d_eval_m2lnp (sd1, x_mid) < fmin))
        return x;

      x0 = x;
      x1 = x + dx;
      x  = x_mid;
    }
    else if ((fmin >= f_x0) || (fmin >= f_x1))
    {
      return x;
    }

    last_x0 = x0;
    last_x1 = x1;
  }

  iter = 0;

  ret = gsl_min_fminimizer_set (self->fmin, &F, x, x0, x1);
  NCM_TEST_GSL_RESULT ("ncm_stats_dist1d_eval_mode", ret);

  do {
    iter++;
    status = gsl_min_fminimizer_iterate (self->fmin);

    if (status)
      g_error ("ncm_stats_dist1d_eval_mode: cannot find the minimum (%s)", gsl_strerror (status));  /* LCOV_EXCL_LINE */

    x  = gsl_min_fminimizer_x_minimum (self->fmin);
    x0 = gsl_min_fminimizer_x_lower (self->fmin);
    x1 = gsl_min_fminimizer_x_upper (self->fmin);

    status = gsl_min_test_interval (x0, x1, self->abstol, reltol);

    if ((status == GSL_CONTINUE) && (x0 == last_x0) && (x1 == last_x1))
    {
      g_warning ("ncm_stats_dist1d_eval_mode: minimization not improving, giving up. (% 22.15g) [% 22.15g % 22.15g]", x, x0, x1); /* LCOV_EXCL_LINE */
      break;                                                                                                                      /* LCOV_EXCL_LINE */
    }

    last_x0 = x0;
    last_x1 = x1;
  } while ((status == GSL_CONTINUE) && (iter < max_iter));

  if (status != GSL_SUCCESS)
    g_warning ("ncm_stats_dist1d_eval_mode: minimization tolerance not achieved" /* LCOV_EXCL_LINE */
               " in %d iterations, giving up. (% 22.15g) [% 22.15g % 22.15g]", max_iter, x, x0, x1);  /* LCOV_EXCL_LINE */

  return x;
}

