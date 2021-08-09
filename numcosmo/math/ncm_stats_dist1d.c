/***************************************************************************
 *            ncm_stats_dist1d.c
 *
 *  Thu February 12 15:37:11 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_stats_dist1d.c
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

enum
{
  PROP_0,
  PROP_XI,
  PROP_XF,
  PROP_NORMA,
  PROP_RELTOL,
  PROP_ABSTOL,
  PROP_MAX_PROB,
  PROP_SIZE,
};

G_DEFINE_ABSTRACT_TYPE (NcmStatsDist1d, ncm_stats_dist1d, G_TYPE_OBJECT);

static gdouble
_ncm_stats_dist1d_inv_cdf_dydx (gdouble y, gdouble x, gpointer userdata)
{
  NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (userdata);
  
  return sd1->norma / (NCM_STATS_DIST1D_GET_CLASS (sd1)->p (sd1, y) * gsl_pow_2 (cosh (x)));
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
  NcmSpline *s1 = ncm_spline_cubic_notaknot_new ();
  NcmSpline *s2 = ncm_spline_cubic_notaknot_new ();
  NcmSpline *s3 = ncm_spline_cubic_notaknot_new ();
  
  sd1->xi       = 0.0;
  sd1->xf       = 0.0;
  sd1->norma    = 0.0;
  sd1->reltol   = 0.0;
  sd1->max_prob = 0.0;
  sd1->inv_cdf  = ncm_ode_spline_new (s1, _ncm_stats_dist1d_inv_cdf_dydx);
  sd1->pdf      = ncm_ode_spline_new (s3, _ncm_stats_dist1d_pdf_dydx);
  sd1->fmin     = gsl_min_fminimizer_alloc (gsl_min_fminimizer_brent);
  
  /* hnil is not an error in this case */
  /*sd1->inv_pdf->stop_hnil = FALSE;*/
  
  ncm_spline_free (s1);
  ncm_spline_free (s2);
  ncm_spline_free (s3);
}

static void
ncm_stats_dist1d_dispose (GObject *object)
{
  NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (object);
  
  ncm_ode_spline_clear (&sd1->inv_cdf);
  ncm_ode_spline_clear (&sd1->pdf);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist1d_parent_class)->dispose (object);
}

static void
ncm_stats_dist1d_finalize (GObject *object)
{
  NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (object);
  
  gsl_min_fminimizer_free (sd1->fmin);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist1d_parent_class)->finalize (object);
}

static void
ncm_stats_dist1d_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (object);
  
  g_return_if_fail (NCM_IS_STATS_DIST1D (object));
  
  switch (prop_id)
  {
    case PROP_XI:
      sd1->xi = g_value_get_double (value);
      break;
    case PROP_XF:
      sd1->xf = g_value_get_double (value);
      break;
    case PROP_RELTOL:
      sd1->reltol = g_value_get_double (value);
      break;
    case PROP_ABSTOL:
      sd1->abstol = g_value_get_double (value);
      break;
    case PROP_MAX_PROB:
      sd1->max_prob = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_stats_dist1d_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (object);
  
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
      g_value_set_double (value, sd1->norma);
      break;
    case PROP_RELTOL:
      g_value_set_double (value, sd1->reltol);
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, sd1->abstol);
      break;
    case PROP_MAX_PROB:
      g_value_set_double (value, sd1->max_prob);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
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
  
  klass->p       = NULL;
  klass->prepare = NULL;
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
  
  if (sd1_class->prepare != NULL)
    sd1_class->prepare (sd1);
  
  if (G_LIKELY (sd1->xi != sd1->xf))
  {
    ncm_ode_spline_set_reltol (sd1->inv_cdf, sd1->reltol);
    ncm_ode_spline_set_reltol (sd1->pdf, sd1->reltol);
    
    ncm_ode_spline_set_abstol (sd1->inv_cdf, sd1->abstol);
    ncm_ode_spline_set_abstol (sd1->pdf, GSL_DBL_EPSILON * 10.0); /* Avoid too much accuracy problems. */
    
    ncm_ode_spline_set_xi (sd1->inv_cdf, 0.0);
    ncm_ode_spline_set_yi (sd1->inv_cdf, sd1->xi);
    ncm_ode_spline_set_yf (sd1->inv_cdf, sd1->xf);
    
    ncm_ode_spline_set_interval (sd1->pdf, 0.0, sd1->xi, sd1->xf);
    
    sd1->norma = 1.0;
    ncm_ode_spline_prepare (sd1->pdf, sd1);
    sd1->norma = ncm_spline_eval (ncm_ode_spline_peek_spline (sd1->pdf), sd1->xf);
    
    ncm_ode_spline_prepare (sd1->inv_cdf, sd1);
  }
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
  return sd1->xi;
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
  return sd1->xf;
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
  if (G_UNLIKELY (sd1->xi == sd1->xf))
    return sd1->xi == x ? 1.0 : 0.0;
  
  return NCM_STATS_DIST1D_GET_CLASS (sd1)->p (sd1, x) / sd1->norma;
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
  if (G_UNLIKELY (sd1->xi == sd1->xf))
    return sd1->xi == x ? 0.0 : GSL_NEGINF;
  
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
  if (G_UNLIKELY (sd1->xi == sd1->xf))
    return sd1->xi <= x ? 1.0 : 0.0;
  
  return ncm_spline_eval (ncm_ode_spline_peek_spline (sd1->pdf), x) / sd1->norma;
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
  if (G_UNLIKELY (sd1->xi == sd1->xf))
    return 1.0;
  
  return sd1->norma;
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
  if (G_UNLIKELY (sd1->xi == sd1->xf))
    return sd1->xi;
  else if (G_UNLIKELY (u <= 0.0))
    return sd1->xi;
  else if (G_UNLIKELY (u >= 1.0))
    return sd1->xf;
  else
    return ncm_spline_eval (ncm_ode_spline_peek_spline (sd1->inv_cdf), atanh (u));
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
  if (G_UNLIKELY (sd1->xi == sd1->xf))
  {
    return sd1->xi;
  }
  else if (G_UNLIKELY (v <= 0.0))
  {
    return sd1->xf;
  }
  else if (G_UNLIKELY (v >= 1.0))
  {
    return sd1->xi;
  }
  else
  {
    const gdouble atan_1mv = 0.5 * (M_LN2 + log1p (-0.5 * v) - log (v));
    
    return ncm_spline_eval (ncm_ode_spline_peek_spline (sd1->inv_cdf), atan_1mv);
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
  const gdouble u = gsl_rng_uniform (rng->r);
  
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
  gint status;
  gint iter = 0, max_iter = 1000000;
  gsl_function F;
  gdouble x0      = sd1->xi;
  gdouble x1      = sd1->xf;
  gdouble x       = 0.5 * (sd1->xf + sd1->xi);
  gdouble last_x0 = x0;
  gdouble last_x1 = x1;
  
  if (G_UNLIKELY (sd1->xi == sd1->xf))
    return sd1->xi;
  
  F.params   = sd1;
  F.function = &_ncm_stats_dist1d_m2lnp;
  
  gsl_min_fminimizer_set (sd1->fmin, &F, x, x0, x1);
  
  do {
    iter++;
    status = gsl_min_fminimizer_iterate (sd1->fmin);
    
    if (status)
      g_error ("ncm_stats_dist1d_mode: Cannot find minimum (%s)", gsl_strerror (status));
    
    x  = gsl_min_fminimizer_x_minimum (sd1->fmin);
    x0 = gsl_min_fminimizer_x_lower (sd1->fmin);
    x1 = gsl_min_fminimizer_x_upper (sd1->fmin);
    
    status = gsl_min_test_interval (x0, x1, sd1->abstol, sd1->reltol);
    
    if ((status == GSL_CONTINUE) && (x0 == last_x0) && (x1 == last_x1))
    {
      g_warning ("ncm_stats_dist1d_eval_mode: minimization not improving, giving up. (% 22.15g) [% 22.15g % 22.15g]", x, x0, x1);
      break;
    }
    
    last_x0 = x0;
    last_x1 = x1;
  } while (status == GSL_CONTINUE && iter < max_iter);
  
  return x;
}

