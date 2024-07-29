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
 * SECTION:ncm_stats_dist1d
 * @title: NcmStatsDist1d
 * @short_description: Abstract class for implementing one dimensional probability distributions
 *
 * Abstract class to reconstruct an arbitrary one dimensional probability distribution.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_stats_dist1d.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_cfg.h"

enum
{
  PROP_0,
  PROP_XI,
  PROP_XF,
  PROP_NORMA,
  PROP_RELTOL,
  PROP_ABSTOL,
  PROP_MAX_PROB,
  PROP_COMPUTE_CDF,
  PROP_SIZE,
};

typedef struct _NcmStatsDist1dPrivate
{
  /*< private >*/
  GObject parent_instance;
  gdouble xi;
  gdouble xf;
  gdouble norma;
  gdouble reltol;
  gdouble abstol;
  gdouble max_prob;
  gboolean compute_cdf;
  NcmOdeSpline *inv_cdf;
  NcmOdeSpline *pdf;
  gsl_min_fminimizer *fmin;
} NcmStatsDist1dPrivate;

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmStatsDist1d, ncm_stats_dist1d, G_TYPE_OBJECT)

static gdouble
_ncm_stats_dist1d_inv_cdf_dydx (gdouble y, gdouble x, gpointer userdata)
{
  NcmStatsDist1d *sd1         = NCM_STATS_DIST1D (userdata);
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);

  return self->norma / (NCM_STATS_DIST1D_GET_CLASS (sd1)->p (sd1, y) * gsl_pow_2 (cosh (x)));
}

static gdouble
_ncm_stats_dist1d_pdf_dydx (gdouble y, gdouble x, gpointer userdata)
{
  NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (userdata);

  /*printf ("B % 20.15g % 20.15g % 20.15g\n", x, y, ncm_stats_dist1d_eval_p (sd1, x));*/
  return NCM_STATS_DIST1D_GET_CLASS (sd1)->p (sd1, x);
}

static void
ncm_stats_dist1d_init (NcmStatsDist1d *sd1)
{
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);
  NcmSpline *s1               = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  NcmSpline *s2               = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  NcmSpline *s3               = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());

  self->xi          = 0.0;
  self->xf          = 0.0;
  self->norma       = 0.0;
  self->reltol      = 0.0;
  self->max_prob    = 0.0;
  self->inv_cdf     = ncm_ode_spline_new (s1, _ncm_stats_dist1d_inv_cdf_dydx);
  self->pdf         = ncm_ode_spline_new (s3, _ncm_stats_dist1d_pdf_dydx);
  self->fmin        = gsl_min_fminimizer_alloc (gsl_min_fminimizer_brent);
  self->compute_cdf = FALSE;

  /* hnil is not an error in this case */
  /*self->inv_pdf->stop_hnil = FALSE;*/

  ncm_spline_free (s1);
  ncm_spline_free (s2);
  ncm_spline_free (s3);
}

static void
ncm_stats_dist1d_dispose (GObject *object)
{
  NcmStatsDist1d *sd1         = NCM_STATS_DIST1D (object);
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);

  ncm_ode_spline_clear (&self->inv_cdf);
  ncm_ode_spline_clear (&self->pdf);

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
    case PROP_MAX_PROB:
      self->max_prob = g_value_get_double (value);
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
    case PROP_MAX_PROB:
      g_value_set_double (value, self->max_prob);
      break;
    case PROP_COMPUTE_CDF:
      g_value_set_boolean (value, ncm_stats_dist1d_get_compute_cdf (sd1));
      break;

    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

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
                                                        "Absolute tolerance on the random variables",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_MAX_PROB,
                                   g_param_spec_double ("max-prob",
                                                        NULL,
                                                        "Maximal probability considered",
                                                        0.0, 1.0, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_COMPUTE_CDF,
                                   g_param_spec_boolean ("compute-cdf",
                                                         NULL,
                                                         "Whether to compute CDF and inverse CDF",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->p             = NULL;
  klass->prepare       = NULL;
  klass->get_current_h = NULL;
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
 *
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
 * Decreases the reference count of *@sd1 and sets the pointer *@sd1 to NULL.
 *
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
 * Prepares the object for calculations.
 */
void
ncm_stats_dist1d_prepare (NcmStatsDist1d *sd1)
{
  NcmStatsDist1dClass *sd1_class = NCM_STATS_DIST1D_GET_CLASS (sd1);
  NcmStatsDist1dPrivate *self    = ncm_stats_dist1d_get_instance_private (sd1);

  if (sd1_class->prepare != NULL)
    sd1_class->prepare (sd1);

  if (G_LIKELY (self->xi != self->xf) && self->compute_cdf)
  {
    ncm_ode_spline_set_reltol (self->inv_cdf, self->reltol);
    ncm_ode_spline_set_reltol (self->pdf, self->reltol);

    ncm_ode_spline_set_abstol (self->inv_cdf, self->abstol);
    ncm_ode_spline_set_abstol (self->pdf, GSL_DBL_EPSILON * 10.0); /* Avoid too much accuracy problems. */

    ncm_ode_spline_set_xi (self->inv_cdf, 0.0);
    ncm_ode_spline_set_yi (self->inv_cdf, self->xi);
    ncm_ode_spline_set_yf (self->inv_cdf, self->xf);

    ncm_ode_spline_set_interval (self->pdf, 0.0, self->xi, self->xf);

    self->norma = 1.0;
    ncm_ode_spline_prepare (self->pdf, sd1);
    self->norma = ncm_spline_eval (ncm_ode_spline_peek_spline (self->pdf), self->xf);

    ncm_ode_spline_prepare (self->inv_cdf, sd1);
  }
  else
  {
    self->norma = 1.0;
  }
}

/**
 * ncm_stats_dist1d_set_xi:
 * @sd1: a #NcmStatsDist1d
 * @xi: a double
 *
 * Sets the lower bound of the distribution $x_i$.
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
 * @xf: a double
 *
 * Sets the upper bound of the distribution $x_f$.
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
 * Returns: the lower bound of the distribution $x_i$.
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
 * Returns: the upper bound of the distribution $x_f$.
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
 * Returns: the current value of the bandwidth h.
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
 * @compute_cdf: a boolean
 *
 * Enable/Disable the computation of the CDF and inverse CDF
 * whenever @compute_cdf is TRUE/FALSE.
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
 * Returns: If the @sd1 is computing the CDF and inverse CDF.
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
 * Calculates the value of the probability density at @x.
 *
 * Returns: the value of the probability density at @x.
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
 * Calculates the value of the $-2\ln(p(x))$ for the probability density.
 * It can be unnormalized, the norma can be retrieved using
 * ncm_stats_dist1d_eval_norma().
 *
 * Returns: the value of $-2\ln(p(x))$.
 */
gdouble
ncm_stats_dist1d_eval_m2lnp (NcmStatsDist1d *sd1, gdouble x)
{
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);

  if (G_UNLIKELY (self->xi == self->xf))
    return self->xi == x ? 0.0 : GSL_NEGINF;

  return NCM_STATS_DIST1D_GET_CLASS (sd1)->m2lnp (sd1, x);
}

/**
 * ncm_stats_dist1d_eval_pdf:
 * @sd1: a #NcmStatsDist1d
 * @x: random variable value
 *
 * Calculates the value of the probability of the interval [x_i, @x].
 *
 * Returns: the value of the probability of the interval [x_i, @x].
 */
gdouble
ncm_stats_dist1d_eval_pdf (NcmStatsDist1d *sd1, gdouble x)
{
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);

  if (G_UNLIKELY (self->xi == self->xf))
    return self->xi <= x ? 1.0 : 0.0;

  return ncm_spline_eval (ncm_ode_spline_peek_spline (self->pdf), x) / self->norma;
}

/**
 * ncm_stats_dist1d_eval_norma:
 * @sd1: a #NcmStatsDist1d
 *
 * Calculates the norma of the distribution. If the probability
 * density is already normalized it will return 1.0.
 *
 * Returns: the value distribution normalization.
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
 * @u: a number between [0, 1]
 *
 * Calculates the value of the random variable $x$ for which the cumulative
 * distribution satisfy $\int_{x_i}^x\mathrm{d}x^\prime p(x^\prime) = u$.
 *
 * Returns: the value of x.
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
    return ncm_spline_eval (ncm_ode_spline_peek_spline (self->inv_cdf), atanh (u));
}

/**
 * ncm_stats_dist1d_eval_inv_pdf_tail:
 * @sd1: a #NcmStatsDist1d
 * @v: a number between [0, 1]
 *
 * Calculates the value of the random variable $x$ for which the cumulative
 * distribution satisfy $\int_{x}^{x_f}\mathrm{d}x^\prime p(x^\prime) = v$.
 *
 * Returns: the value of x.
 */
gdouble
ncm_stats_dist1d_eval_inv_pdf_tail (NcmStatsDist1d *sd1, const gdouble v)
{
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);

  if (G_UNLIKELY (self->xi == self->xf))
  {
    return self->xi;
  }
  else if (G_UNLIKELY (v <= 0.0))
  {
    return self->xf;
  }
  else if (G_UNLIKELY (v >= 1.0))
  {
    return self->xi;
  }
  else
  {
    const gdouble atan_1mv = 0.5 * (M_LN2 + log1p (-0.5 * v) - log (v));

    return ncm_spline_eval (ncm_ode_spline_peek_spline (self->inv_cdf), atan_1mv);
  }
}

/**
 * ncm_stats_dist1d_gen:
 * @sd1: a #NcmStatsDist1d
 * @rng: a #NcmRNG
 *
 * Generates a realization of the probability distribution.
 *
 * Returns: the value of the probability of the interval [x_i, @x].
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
 * Calculates the mode of the distribution.
 *
 * Returns: the mode of the probability distribution.
 */
gdouble
ncm_stats_dist1d_eval_mode (NcmStatsDist1d *sd1)
{
  NcmStatsDist1dPrivate *self = ncm_stats_dist1d_get_instance_private (sd1);
  const gdouble reltol        = sqrt (self->reltol);
  const gint max_iter         = 1000000;
  const gint linear_search    = 1000;
  gdouble x0                  = self->xi;
  gdouble x1                  = self->xf;
  gdouble x                   = 0.5 * (self->xf + self->xi);
  gdouble last_x0             = x0;
  gdouble last_x1             = x1;
  gint iter                   = 0;
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
    const gdouble x_try = self->xi + (self->xf - self->xi) / (linear_search - 1.0) * iter;
    const gdouble f_try = ncm_stats_dist1d_eval_m2lnp (sd1, x_try);

    if (f_try < fmin)
    {
      fmin = f_try;
      x    = x_try;
    }
  }

  iter = 0;

  ret = gsl_min_fminimizer_set (self->fmin, &F, x, x0, x1);
  NCM_TEST_GSL_RESULT ("ncm_stats_dist1d_eval_mode", ret);

  do {
    iter++;
    status = gsl_min_fminimizer_iterate (self->fmin);

    if (status)
      g_error ("ncm_stats_dist1d_mode: Cannot find minimum (%s)", gsl_strerror (status));  /* LCOV_EXCL_LINE */

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

