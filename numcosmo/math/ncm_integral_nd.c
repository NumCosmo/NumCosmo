/***************************************************************************
 *            ncm_integral_nd.c
 *
 *  Thu July 20 08:39:30 2023
 *  Copyright  2023 Eduardo José Barroso
 *  <eduardo.jsbarroso@uel.br>
 ****************************************************************************/
/*
 * ncm_integral_nd.c
 * Copyright (C) 2023 Eduardo José Barroso <eduardo.jsbarroso@uel.br>
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
 * NcmIntegralND:
 *
 * N-dimensional integration object.
 *
 * This object is used to perform n-dimensional integration of a function using
 * different methods.
 *
 * The integration can be performed using the cubature library. The cubature
 * library is a library for adaptive multidimensional integration. The original
 * code can be found at https://github.com/stevengj/cubature.
 *
 * To use this object, the user must initialize a child object that implements
 * two functions: get_dimensions and integrand. To do so, the user must
 * first define these two functions as the prototypes described
 * in ncm_integral_nd.h. The get dimensions function should return the
 * dimension of the arguments to be integrated and the function dimension.
 * For instance, if the integrand is given by $F(x,y) = x + y$, the get
 * dimensions function must return $(2,1)$, such that it will compute
 * \begin{align}
 * \int \int F(x, y) dxdy
 * ,\end{align}
 * returning a scalar for the integral evaluated in the given intervals.
 *
 * This object can also be used with multi-dimensional functions that return
 * an array instead of a scalar. Considering the integrand
 * $F(x,y,z) = [x^2, y+z ]$, the get dimensions method should return $(3,2)$
 * and the object will compute the integral
 * \begin{align}
 * \int \int \int F(x,y,z) dxdy = [\frac{yzx^3}{3}, frac{x(y^2+z^2)}{2}]
 * \end{align}
 * for the given intervals.
 *
 * Having the functions, the user must instantiate an object of the type
 * #NcmIntegralNDClass defined with these functions. To do so, one must call the macro
 * #NCM_INTEGRAL_ND_DEFINE_TYPE to define the new object type, which
 * will later be instantiable. Examples of how to define the objects
 * containing the integrand can be found in the test folder under
 * test\textunderscore ncm\textunderscore integral\textunderscore nd.c.
 * For an example of the Python implementation of the integrand in a class,
 * check test\textunderscore py\textunderscore integralnd.py in the same folder.
 * This object cannot be used without the child object containing the cited functions.
 *
 * After defining the child class with the necessary functions,
 * the user may use the integration object with the prefered method
 * from the cubature library.
 *
 * The user may provide the input values for: @rel_tol - ncm_integral_nd_set_reltol(), @abs_tol - ncm_integral_nd_set_abstol(),
 * @integ_method - ncm_integral_nd_set_method(), @max_eval - ncm_integral_nd_set_maxeval(), @error - ncm_integral_nd_set_error().
 * If these functions are not called, default parameters are chosen.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_integral_nd.h"
#include "ncm_enum_types.h"
#include "math/ncm_c.h"
#include "math/ncm_cfg.h"

#include "misc/cubature.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcmIntegralNDPrivate
{
  NcmIntegralNDMethod method;
  NcmIntegralNDError error;
  guint maxeval;
  gdouble reltol;
  gdouble abstol;
  NcmVector *x_vec;
  NcmVector *fval_vec;
} NcmIntegralNDPrivate;

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmIntegralND, ncm_integral_nd, G_TYPE_OBJECT)

enum
{
  PROP_0,
  PROP_METHOD,
  PROP_ERROR,
  PROP_MAXEVAL,
  PROP_RELTOL,
  PROP_ABSTOL,
  PROP_SIZE,
};

static void
ncm_integral_nd_init (NcmIntegralND *intnd)
{
  NcmIntegralNDPrivate * const self = ncm_integral_nd_get_instance_private (intnd);

  self->method  = NCM_INTEGRAL_ND_METHOD_LEN;
  self->error   = NCM_INTEGRAL_ND_ERROR_LEN;
  self->maxeval = 0;
  self->reltol  = 0.0;
  self->abstol  = 0.0;

  /* These are dummy vectors to be used in the integrand function */
  self->x_vec    = ncm_vector_new_data_static ((gdouble *) 1, 1, 1);
  self->fval_vec = ncm_vector_new_data_static ((gdouble *) 1, 1, 1);
}

static void
ncm_integral_nd_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmIntegralND *intnd = NCM_INTEGRAL_ND (object);

  g_return_if_fail (NCM_IS_INTEGRAL_ND (object));

  switch (prop_id)
  {
    case PROP_METHOD:
      ncm_integral_nd_set_method (intnd, g_value_get_enum (value));
      break;
    case PROP_ERROR:
      ncm_integral_nd_set_error (intnd, g_value_get_enum (value));
      break;
    case PROP_MAXEVAL:
      ncm_integral_nd_set_maxeval (intnd, g_value_get_uint (value));
      break;
    case PROP_RELTOL:
      ncm_integral_nd_set_reltol (intnd, g_value_get_double (value));
      break;
    case PROP_ABSTOL:
      ncm_integral_nd_set_abstol (intnd, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_integral_nd_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmIntegralND *intnd = NCM_INTEGRAL_ND (object);

  g_return_if_fail (NCM_IS_INTEGRAL_ND (object));

  switch (prop_id)
  {
    case PROP_METHOD:
      g_value_set_enum (value, ncm_integral_nd_get_method (intnd));
      break;
    case PROP_ERROR:
      g_value_set_enum (value, ncm_integral_nd_get_error (intnd));
      break;
    case PROP_MAXEVAL:
      g_value_set_uint (value, ncm_integral_nd_get_maxeval (intnd));
      break;
    case PROP_RELTOL:
      g_value_set_double (value, ncm_integral_nd_get_reltol (intnd));
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, ncm_integral_nd_get_abstol (intnd));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_integral_nd_finalize (GObject *object)
{
  NcmIntegralND *intnd              = NCM_INTEGRAL_ND (object);
  NcmIntegralNDPrivate * const self = ncm_integral_nd_get_instance_private (intnd);

  ncm_vector_clear (&self->x_vec);
  ncm_vector_clear (&self->fval_vec);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_integral_nd_parent_class)->finalize (object);
}

void _ncm_integral_nd_eval (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);
void _ncm_integral_nd_get_dimensions (NcmIntegralND *intnd, guint *dim, guint *fdim);

static void
ncm_integral_nd_class_init (NcmIntegralNDClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &ncm_integral_nd_set_property;
  object_class->get_property = &ncm_integral_nd_get_property;
  object_class->finalize     = &ncm_integral_nd_finalize;

  g_object_class_install_property (object_class,
                                   PROP_METHOD,
                                   g_param_spec_enum ("method",
                                                      NULL,
                                                      "Integration method",
                                                      NCM_TYPE_INTEGRAL_ND_METHOD,
                                                      NCM_INTEGRAL_ND_METHOD_CUBATURE_H,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_ERROR,
                                   g_param_spec_enum ("error",
                                                      NULL,
                                                      "Error measure",
                                                      NCM_TYPE_INTEGRAL_ND_ERROR,
                                                      NCM_INTEGRAL_ND_ERROR_INDIVIDUAL,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_MAXEVAL,
                                   g_param_spec_uint ("maxeval",
                                                      NULL,
                                                      "Maximum number of function evaluations (0 means unlimited)",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Integral relative tolerance",
                                                        0.0, 1.0, NCM_INTEGRAL_ND_DEFAULT_RELTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ABSTOL,
                                   g_param_spec_double ("abstol",
                                                        NULL,
                                                        "Integral absolute tolerance",
                                                        0.0, G_MAXDOUBLE, NCM_INTEGRAL_ND_DEFAULT_ABSTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->integrand      = &_ncm_integral_nd_eval;
  klass->get_dimensions = &_ncm_integral_nd_get_dimensions;
}

/* LCOV_EXCL_START */

void
_ncm_integral_nd_eval (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval)
{
  g_error ("ncm_integral_nd_eval not implemented by subclass: `%s`.", G_OBJECT_TYPE_NAME (intnd));
}

void
_ncm_integral_nd_get_dimensions (NcmIntegralND *intnd, guint *dim, guint *fdim)
{
  g_error ("ncm_integral_nd_get_dimensions not implemented by subclass: `%s`.", G_OBJECT_TYPE_NAME (intnd));
}

/* LCOV_EXCL_STOP */

/**
 * ncm_integral_nd_ref:
 * @intnd: a #NcmIntegralND
 *
 * Increases the reference count of @intnd by one.
 *
 * Returns: (transfer full): @intnd.
 */
NcmIntegralND *
ncm_integral_nd_ref (NcmIntegralND *intnd)
{
  return g_object_ref (intnd);
}

/**
 * ncm_integral_nd_free:
 * @intnd: a #NcmIntegralND
 *
 * Decreases the reference count of @intnd by one.
 *
 */
void
ncm_integral_nd_free (NcmIntegralND *intnd)
{
  g_object_unref (intnd);
}

/**
 * ncm_integral_nd_clear:
 * @intnd: a #NcmIntegralND
 *
 * If *@intnd is different from NULL, decreases the reference
 * count of *@intnd by one and sets *@intnd to NULL.
 *
 */
void
ncm_integral_nd_clear (NcmIntegralND **intnd)
{
  g_clear_object (intnd);
}

/**
 * ncm_integral_nd_set_method:
 * @intnd: a #NcmIntegralND
 * @method: a #NcmIntegralNDMethod
 *
 * Sets the integration method to use.
 *
 */
void
ncm_integral_nd_set_method (NcmIntegralND *intnd, NcmIntegralNDMethod method)
{
  NcmIntegralNDPrivate * const self = ncm_integral_nd_get_instance_private (intnd);

  self->method = method;
}

/**
 * ncm_integral_nd_set_error:
 * @intnd: a #NcmIntegralND
 * @error: a #NcmIntegralNDError
 *
 * Sets the error measure to use.
 *
 */
void
ncm_integral_nd_set_error (NcmIntegralND *intnd, NcmIntegralNDError error)
{
  NcmIntegralNDPrivate * const self = ncm_integral_nd_get_instance_private (intnd);

  self->error = error;
}

/**
 * ncm_integral_nd_set_maxeval:
 * @intnd: a #NcmIntegralND
 * @maxeval: maximum number of function evaluations
 *
 * Sets the maximum number of function evaluations to use.
 * Zero means unlimited.
 */
void
ncm_integral_nd_set_maxeval (NcmIntegralND *intnd, guint maxeval)
{
  NcmIntegralNDPrivate * const self = ncm_integral_nd_get_instance_private (intnd);

  self->maxeval = maxeval;
}

/**
 * ncm_integral_nd_set_reltol:
 * @intnd: a #NcmIntegralND
 * @reltol: relative tolerance
 *
 * Sets the relative tolerance @reltol to use.
 *
 */
void
ncm_integral_nd_set_reltol (NcmIntegralND *intnd, gdouble reltol)
{
  NcmIntegralNDPrivate * const self = ncm_integral_nd_get_instance_private (intnd);

  self->reltol = reltol;
}

/**
 * ncm_integral_nd_set_abstol:
 * @intnd: a #NcmIntegralND
 * @abstol: absolute tolerance
 *
 * Sets the absolute tolerance @reltol to use.
 *
 */
void
ncm_integral_nd_set_abstol (NcmIntegralND *intnd, gdouble abstol)
{
  NcmIntegralNDPrivate * const self = ncm_integral_nd_get_instance_private (intnd);

  self->abstol = abstol;
}

/**
 * ncm_integral_nd_get_method:
 * @intnd: a #NcmIntegralND
 *
 * Returns: the integration method used.
 */
NcmIntegralNDMethod
ncm_integral_nd_get_method (NcmIntegralND *intnd)
{
  NcmIntegralNDPrivate * const self = ncm_integral_nd_get_instance_private (intnd);

  return self->method;
}

/**
 * ncm_integral_nd_get_error:
 * @intnd: a #NcmIntegralND
 *
 * Returns: the error measure used.
 */
NcmIntegralNDError
ncm_integral_nd_get_error (NcmIntegralND *intnd)
{
  NcmIntegralNDPrivate * const self = ncm_integral_nd_get_instance_private (intnd);

  return self->error;
}

/**
 * ncm_integral_nd_get_maxeval:
 * @intnd: a #NcmIntegralND
 *
 * Returns: the maximum number of function evaluations used.
 */
guint
ncm_integral_nd_get_maxeval (NcmIntegralND *intnd)
{
  NcmIntegralNDPrivate * const self = ncm_integral_nd_get_instance_private (intnd);

  return self->maxeval;
}

/**
 * ncm_integral_nd_get_reltol:
 * @intnd: a #NcmIntegralND
 *
 * Returns: the relative tolerance used.
 */
gdouble
ncm_integral_nd_get_reltol (NcmIntegralND *intnd)
{
  NcmIntegralNDPrivate * const self = ncm_integral_nd_get_instance_private (intnd);

  return self->reltol;
}

/**
 * ncm_integral_nd_get_abstol:
 * @intnd: a #NcmIntegralND
 *
 * Returns: the absolute tolerance used.
 */
gdouble
ncm_integral_nd_get_abstol (NcmIntegralND *intnd)
{
  NcmIntegralNDPrivate * const self = ncm_integral_nd_get_instance_private (intnd);

  return self->abstol;
}

static gint
_ncm_integral_nd_cubature_int (unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
  NcmIntegralND *intnd              = NCM_INTEGRAL_ND (fdata);
  NcmIntegralNDPrivate * const self = ncm_integral_nd_get_instance_private (intnd);

  ncm_vector_replace_data_full (self->x_vec, (gdouble *) x, ndim, 1);
  ncm_vector_replace_data_full (self->fval_vec, fval, fdim, 1);

  NCM_INTEGRAL_ND_GET_CLASS (intnd)->integrand (intnd, self->x_vec, ndim, 1, fdim, self->fval_vec);

  return 0;
}

static gint
_ncm_integral_nd_cubature_vint (unsigned ndim, size_t npt, const double *x, void *fdata, unsigned fdim, double *fval)
{
  NcmIntegralND *intnd              = NCM_INTEGRAL_ND (fdata);
  NcmIntegralNDPrivate * const self = ncm_integral_nd_get_instance_private (intnd);

  ncm_vector_replace_data_full (self->x_vec, (gdouble *) x, ndim * npt, 1);
  ncm_vector_replace_data_full (self->fval_vec, fval, fdim * npt, 1);

  NCM_INTEGRAL_ND_GET_CLASS (intnd)->integrand (intnd, self->x_vec, ndim, npt, fdim, self->fval_vec);

  return 0;
}

/**
 * ncm_integral_nd_eval:
 * @intnd: a #NcmIntegralND
 * @xi: a #NcmVector containing the inferior integration limit $x_i$
 * @xf: a #NcmVector containing the superior integration limit $x_f$
 * @res: a #NcmVector containing the result of the integration
 * @err: a #NcmVector containing the error of the integration
 *
 * Evaluated the integral $I_F(x_i, x_f) = \int_{x_i}^{x_f}F(x)\mathrm{d}x$.
 *
 */
void
ncm_integral_nd_eval (NcmIntegralND *intnd, const NcmVector *xi, const NcmVector *xf, NcmVector *res, NcmVector *err)
{
  NcmIntegralNDPrivate * const self = ncm_integral_nd_get_instance_private (intnd);
  gint error                        = 0;

  guint dim, fdim;
  gint ret;


  NCM_INTEGRAL_ND_GET_CLASS (intnd)->get_dimensions (intnd, &dim, &fdim);

  g_assert_cmpuint (ncm_vector_len (xi), ==, dim);
  g_assert_cmpuint (ncm_vector_len (xf), ==, dim);
  g_assert_cmpuint (ncm_vector_len (res), ==, fdim);
  g_assert_cmpuint (ncm_vector_len (err), ==, fdim);

  switch (self->error)
  {
    case NCM_INTEGRAL_ND_ERROR_INDIVIDUAL:
      error = ERROR_INDIVIDUAL;
      break;
    case NCM_INTEGRAL_ND_ERROR_PAIRWISE:
      error = ERROR_PAIRED;
      break;
    case NCM_INTEGRAL_ND_ERROR_L2:
      error = ERROR_L2;
      break;
    case NCM_INTEGRAL_ND_ERROR_L1:
      error = ERROR_L1;
      break;
    case NCM_INTEGRAL_ND_ERROR_LINF:
      error = ERROR_LINF;
      break;
    default:
      g_error ("ncm_integral_nd_eval: invalid error measure: `%d`.", self->error);
      break;
  }


  switch (self->method)
  {
    case NCM_INTEGRAL_ND_METHOD_CUBATURE_H:
      ret = hcubature (
        fdim,
        _ncm_integral_nd_cubature_int,
        intnd,
        dim,
        ncm_vector_const_data (xi),
        ncm_vector_const_data (xf),
        self->maxeval,
        self->abstol,
        self->reltol,
        error,
        ncm_vector_data (res),
        ncm_vector_data (err)
                      );
      break;
    case NCM_INTEGRAL_ND_METHOD_CUBATURE_P:
      ret = pcubature (
        fdim,
        _ncm_integral_nd_cubature_int,
        intnd,
        dim,
        ncm_vector_const_data (xi),
        ncm_vector_const_data (xf),
        self->maxeval,
        self->abstol,
        self->reltol,
        error,
        ncm_vector_data (res),
        ncm_vector_data (err)
                      );
      break;
    case NCM_INTEGRAL_ND_METHOD_CUBATURE_H_V:
      ret = hcubature_v (
        fdim,
        _ncm_integral_nd_cubature_vint,
        intnd,
        dim,
        ncm_vector_const_data (xi),
        ncm_vector_const_data (xf),
        self->maxeval,
        self->abstol,
        self->reltol,
        error,
        ncm_vector_data (res),
        ncm_vector_data (err)
                        );
      break;
    case NCM_INTEGRAL_ND_METHOD_CUBATURE_P_V:
      ret = pcubature_v (
        fdim,
        _ncm_integral_nd_cubature_vint,
        intnd,
        dim,
        ncm_vector_const_data (xi),
        ncm_vector_const_data (xf),
        self->maxeval,
        self->abstol,
        self->reltol,
        error,
        ncm_vector_data (res),
        ncm_vector_data (err)
                        );
      break;
    default:
      g_error ("ncm_integral_nd_eval: invalid method: `%d`.", self->method);
      break;
  }

  g_assert_cmpint (ret, ==, 0);
}

