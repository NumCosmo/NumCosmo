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
 * @title: FIXME
 * @short_description: FIXME
 *
 * FIXME
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
};

G_DEFINE_ABSTRACT_TYPE (NcmStatsDist1d, ncm_stats_dist1d, G_TYPE_OBJECT);

static gdouble 
_ncm_stats_dist1d_inv_pdf_dydx (gdouble y, gdouble x, gpointer userdata)
{
  NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (userdata);
  /*printf ("A % 20.15g % 20.15g % 20.15g\n", x, y, sd1->norma / ncm_stats_dist1d_p (sd1, y));*/
  return sd1->norma / ncm_stats_dist1d_p (sd1, y);
}

static gdouble 
_ncm_stats_dist1d_pdf_dydx (gdouble y, gdouble x, gpointer userdata)
{
  NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (userdata);
  /*printf ("B % 20.15g % 20.15g % 20.15g\n", x, y, ncm_stats_dist1d_p (sd1, x));*/
  return ncm_stats_dist1d_p (sd1, x);
}

static void
ncm_stats_dist1d_init (NcmStatsDist1d *sd1)
{
  NcmSpline *s1 = ncm_spline_cubic_notaknot_new ();
  NcmSpline *s2 = ncm_spline_cubic_notaknot_new ();
  sd1->xi       = 0.0;
  sd1->xf       = 0.0;
  sd1->norma    = 0.0;
  sd1->reltol   = 0.0;
  sd1->inv_pdf  = ncm_ode_spline_new (s1, _ncm_stats_dist1d_inv_pdf_dydx);
  sd1->pdf      = ncm_ode_spline_new (s2, _ncm_stats_dist1d_pdf_dydx);
  sd1->fmin     = gsl_min_fminimizer_alloc (gsl_min_fminimizer_brent);

  ncm_spline_free (s1);
  ncm_spline_free (s2);
}

static void
ncm_stats_dist1d_dispose (GObject *object)
{
  NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (object);
  
  ncm_ode_spline_clear (&sd1->inv_pdf);
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
      g_value_set_double (value, sd1->xi);
      break;
    case PROP_XF:
      g_value_set_double (value, sd1->xf);
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
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_stats_dist1d_class_init (NcmStatsDist1dClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

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
                                                        0.0, 1.0, NCM_ODE_SPLINE_DEFAULT_RELTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  g_object_class_install_property (object_class,
                                   PROP_ABSTOL,
                                   g_param_spec_double ("abstol",
                                                        NULL,
                                                        "Absolute tolerance on the random variables",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->p       = NULL;
  klass->prepare = NULL;
}

/**
 * ncm_stats_dist1d_ref:
 * @sd1: a #NcmStatsDist1d.
 *
 * Increase the reference count of @sd1.
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
 * @sd1: a #NcmStatsDist1d.
 *
 * Decrease the reference count of @sd1.
 *
 */
void 
ncm_stats_dist1d_free (NcmStatsDist1d *sd1)
{
  g_object_unref (sd1);
}

/**
 * ncm_stats_dist1d_clear:
 * @sd1: a #NcmStatsDist1d.
 *
 * Decrease the reference count of *@sd1 and sets the pointer *@sd1 to NULL.
 *
 */
void 
ncm_stats_dist1d_clear (NcmStatsDist1d **sd1)
{
  g_clear_object (sd1);
}

/**
 * ncm_stats_dist1d_prepare:
 * @sd1: a #NcmStatsDist1d.
 *
 * Prepare the object for calculations.
 */
void 
ncm_stats_dist1d_prepare (NcmStatsDist1d *sd1)
{
  NcmStatsDist1dClass *sd1_class = NCM_STATS_DIST1D_GET_CLASS (sd1); 

  if (sd1_class->prepare != NULL)
    sd1_class->prepare (sd1);

  ncm_ode_spline_set_reltol (sd1->inv_pdf, sd1->reltol);
  ncm_ode_spline_set_reltol (sd1->pdf, sd1->reltol);
  
  ncm_ode_spline_set_abstol (sd1->pdf, GSL_DBL_EPSILON);
  ncm_ode_spline_set_reltol (sd1->pdf, sd1->abstol);
  
  ncm_ode_spline_set_interval (sd1->inv_pdf, sd1->xi, 0.0, 1.0);
  ncm_ode_spline_set_interval (sd1->pdf, 0.0, sd1->xi, sd1->xf);

  ncm_ode_spline_prepare (sd1->pdf, sd1);
  sd1->norma = 1.0;
  sd1->norma = ncm_spline_eval (sd1->pdf->s, sd1->xf);

  ncm_ode_spline_prepare (sd1->inv_pdf, sd1);
}

/**
 * ncm_stats_dist1d_p:
 * @sd1: a #NcmStatsDist1d.
 * @x: random variable value.
 *
 * Calculates the value of the probability density at @x.
 * It can be unnormalized, the norma can be retrieved using
 * ncm_stats_dist1d_norma().
 * 
 * Returns: the value of the probability density at @x.
 */
gdouble 
ncm_stats_dist1d_p (NcmStatsDist1d *sd1, gdouble x)
{
  return NCM_STATS_DIST1D_GET_CLASS (sd1)->p (sd1, x);
}

/**
 * ncm_stats_dist1d_m2lnp:
 * @sd1: a #NcmStatsDist1d.
 * @x: random variable value.
 *
 * Calculates the value of the $-2\ln(p(x))$ for the probability density.
 * It can be unnormalized, the norma can be retrieved using
 * ncm_stats_dist1d_norma().
 * 
 * Returns: the value of $-2\ln(p(x))$.
 */
gdouble 
ncm_stats_dist1d_m2lnp (NcmStatsDist1d *sd1, gdouble x)
{
  return NCM_STATS_DIST1D_GET_CLASS (sd1)->m2lnp (sd1, x);
}

/**
 * ncm_stats_dist1d_pdf:
 * @sd1: a #NcmStatsDist1d.
 * @x: random variable value.
 *
 * Calculates the value of the probability of the interval [x_i, @x].
 * 
 * Returns: the value of the probability of the interval [x_i, @x].
 */
gdouble 
ncm_stats_dist1d_pdf (NcmStatsDist1d *sd1, gdouble x)
{
  return ncm_spline_eval (sd1->pdf->s, x) / sd1->norma;
}

/**
 * ncm_stats_dist1d_norma:
 * @sd1: a #NcmStatsDist1d.
 *
 * Calculates the norma of the distribution. If the probability
 * density is already normalized it will return 1.0.
 * 
 * Returns: the value distribution normalization.
 */
gdouble 
ncm_stats_dist1d_norma (NcmStatsDist1d *sd1)
{
  return sd1->norma;
}

/**
 * ncm_stats_dist1d_inv_pdf:
 * @sd1: a #NcmStatsDist1d.
 * @u: a number between [0, 1].
 *
 * Calculates the value of the random variable $x$ for which the cumulative
 * distribution satisfy $\int_{x_i}^x\mathrm{d}x^\prime p(x^\prime) = u$.
 * 
 * Returns: the value of x.
 */
gdouble 
ncm_stats_dist1d_inv_pdf (NcmStatsDist1d *sd1, gdouble x)
{
  return ncm_spline_eval (sd1->inv_pdf->s, x);
}

/**
 * ncm_stats_dist1d_gen:
 * @sd1: a #NcmStatsDist1d.
 * @rng: a #NcmRNG.
 * 
 * Generates a realization of the probability distribution.
 * 
 * Returns: the value of the probability of the interval [x_i, @x].
 */
gdouble 
ncm_stats_dist1d_gen (NcmStatsDist1d *sd1, NcmRNG *rng)
{
  const gdouble u = gsl_rng_uniform (rng->r);
  return ncm_stats_dist1d_inv_pdf (sd1, u);
}

static gdouble
_ncm_stats_dist1d_m2lnp (gdouble x, gpointer p)
{
  NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (p);
  return ncm_stats_dist1d_m2lnp (sd1, x);
}

/**
 * ncm_stats_dist1d_mode:
 * @sd1: a #NcmStatsDist1d.
 * 
 * Calculates the mode of the distribution.
 * 
 * Returns: the mode of the probability distribution.
 */
gdouble
ncm_stats_dist1d_mode (NcmStatsDist1d *sd1)
{
  gint status;
  gint iter = 0, max_iter = 1000000;
  gsl_function F;
  gdouble x0 = sd1->xi;
  gdouble x1 = sd1->xf;
  gdouble x  = 0.5 * (sd1->xf + sd1->xi);

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
  } while (status == GSL_CONTINUE && iter < max_iter);

  return x;
}
