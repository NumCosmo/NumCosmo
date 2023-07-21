/***************************************************************************
 *            ncm_integralnd.c
 *
 *  Thu July 20 08:39:30 2023
 *  Copyright  2023 Eduardo José Barroso
 *  <eduardo.jsbarroso@uel.br>
 ****************************************************************************/
/*
 * ncm_integralnd.c
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
 * SECTION:ncm_integralnd
 * @title: NcmIntegralnd
 * @short_description: N-dimensional integration object.
 * @stability: Stable
 * @include: numcosmo/math/ncm_integralnd.h
 *
 * This object is used to perform n-dimensional integration of a function
 * using different methods.
 *
 * The integration can be performed using the cubature library. The cubature
 * library is a library for adaptive multidimensional integration.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_integralnd.h"
#include "ncm_enum_types.h"
#include "math/ncm_c.h"
#include "math/ncm_cfg.h"

#include "misc/cubature.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmIntegralndPrivate
{
  NcmIntegralndMethod method;
  NcmIntegralndError error;
  guint maxeval;
  gdouble reltol;
  gdouble abstol;
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmIntegralnd, ncm_integralnd, G_TYPE_OBJECT);

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
ncm_integralnd_init (NcmIntegralnd *intnd)
{
  NcmIntegralndPrivate * const self = intnd->priv = ncm_integralnd_get_instance_private (intnd);

  self->method  = NCM_INTEGRALND_METHOD_LEN;
  self->error   = NCM_INTEGRALND_ERROR_LEN;
  self->maxeval = 0;
  self->reltol  = 0.0;
  self->abstol  = 0.0;
}

static void
ncm_integralnd_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmIntegralnd *intnd = NCM_INTEGRALND (object);

  g_return_if_fail (NCM_IS_INTEGRALND (object));

  switch (prop_id)
  {
    case PROP_METHOD:
      ncm_integralnd_set_method (intnd, g_value_get_enum (value));
      break;
    case PROP_ERROR:
      ncm_integralnd_set_error (intnd, g_value_get_enum (value));
      break;
    case PROP_MAXEVAL:
      ncm_integralnd_set_maxeval (intnd, g_value_get_uint (value));
      break;
    case PROP_RELTOL:
      ncm_integralnd_set_reltol (intnd, g_value_get_double (value));
      break;
    case PROP_ABSTOL:
      ncm_integralnd_set_abstol (intnd, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_integralnd_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmIntegralnd *intnd = NCM_INTEGRALND (object);

  g_return_if_fail (NCM_IS_INTEGRALND (object));

  switch (prop_id)
  {
    case PROP_METHOD:
      g_value_set_enum (value, ncm_integralnd_get_method (intnd));
      break;
    case PROP_ERROR:
      g_value_set_enum (value, ncm_integralnd_get_error (intnd));
      break;
    case PROP_MAXEVAL:
      g_value_set_uint (value, ncm_integralnd_get_maxeval (intnd));
      break;
    case PROP_RELTOL:
      g_value_set_double (value, ncm_integralnd_get_reltol (intnd));
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, ncm_integralnd_get_abstol (intnd));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_integralnd_finalize (GObject *object)
{
  /* NcmIntegralnd *intnd              = NCM_INTEGRALND (object); */
  /* NcmIntegralndPrivate * const self = intnd->priv; */

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_integralnd_parent_class)->finalize (object);
}

void _ncm_integralnd_eval (NcmIntegralnd *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);
void _ncm_integralnd_get_dimensions (NcmIntegralnd *intnd, guint *dim, guint *fdim);

static void
ncm_integralnd_class_init (NcmIntegralndClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &ncm_integralnd_set_property;
  object_class->get_property = &ncm_integralnd_get_property;
  object_class->finalize     = &ncm_integralnd_finalize;

  g_object_class_install_property (object_class,
                                   PROP_METHOD,
                                   g_param_spec_enum ("method",
                                                      NULL,
                                                      "Integration method",
                                                      NCM_TYPE_INTEGRALND_METHOD,
                                                      NCM_INTEGRALND_METHOD_CUBATURE_H,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_ERROR,
                                   g_param_spec_enum ("error",
                                                      NULL,
                                                      "Error measure",
                                                      NCM_TYPE_INTEGRALND_ERROR,
                                                      NCM_INTEGRALND_ERROR_INDIVIDUAL,
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
                                                        0.0, 1.0, NCM_INTEGRALND_DEFAULT_RELTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ABSTOL,
                                   g_param_spec_double ("abstol",
                                                        NULL,
                                                        "Integral absolute tolerance",
                                                        0.0, G_MAXDOUBLE, NCM_INTEGRALND_DEFAULT_ABSTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->integrand      = &_ncm_integralnd_eval;
  klass->get_dimensions = &_ncm_integralnd_get_dimensions;
}

void
_ncm_integralnd_eval (NcmIntegralnd *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval)
{
  g_error ("ncm_integralnd_eval not implemented by subclass: `%s`.", G_OBJECT_TYPE_NAME (intnd));
}

void
_ncm_integralnd_get_dimensions (NcmIntegralnd *intnd, guint *dim, guint *fdim)
{
  g_error ("ncm_integralnd_get_dimensions not implemented by subclass: `%s`.", G_OBJECT_TYPE_NAME (intnd));
}

/**
 * ncm_integralnd_ref:
 * @intnd: a #NcmIntegralnd
 *
 * Increases the reference count of @intnd by one.
 *
 * Returns: (transfer full): @intnd.
 */
NcmIntegralnd *
ncm_integralnd_ref (NcmIntegralnd *intnd)
{
  return g_object_ref (intnd);
}

/**
 * ncm_integralnd_free:
 * @intnd: a #NcmIntegralnd
 *
 * Decreases the reference count of @intnd by one.
 *
 */
void
ncm_integralnd_free (NcmIntegralnd *intnd)
{
  g_object_unref (intnd);
}

/**
 * ncm_integralnd_clear:
 * @intnd: a #NcmIntegralnd
 *
 * If *@intnd is different from NULL, decreases the reference
 * count of *@intnd by one and sets *@intnd to NULL.
 *
 */
void
ncm_integralnd_clear (NcmIntegralnd **intnd)
{
  g_clear_object (intnd);
}

/**
 * ncm_integralnd_set_method:
 * @intnd: a #NcmIntegralnd
 * @method: a #NcmIntegralndMethod
 *
 * Sets the integration method to use.
 *
 */
void
ncm_integralnd_set_method (NcmIntegralnd *intnd, NcmIntegralndMethod method)
{
  NcmIntegralndPrivate * const self = intnd->priv;

  self->method = method;
}

/**
 * ncm_integralnd_set_error:
 * @intnd: a #NcmIntegralnd
 * @error: a #NcmIntegralndError
 *
 * Sets the error measure to use.
 *
 */
void
ncm_integralnd_set_error (NcmIntegralnd *intnd, NcmIntegralndError error)
{
  NcmIntegralndPrivate * const self = intnd->priv;

  self->error = error;
}

/**
 * ncm_integralnd_set_maxeval:
 * @intnd: a #NcmIntegralnd
 * @maxeval: maximum number of function evaluations
 *
 * Sets the maximum number of function evaluations to use.
 * Zero means unlimited.
 */
void
ncm_integralnd_set_maxeval (NcmIntegralnd *intnd, guint maxeval)
{
  NcmIntegralndPrivate * const self = intnd->priv;

  self->maxeval = maxeval;
}

/**
 * ncm_integralnd_set_reltol:
 * @intnd: a #NcmIntegralnd
 * @reltol: relative tolerance
 *
 * Sets the relative tolerance @reltol to use.
 *
 */
void
ncm_integralnd_set_reltol (NcmIntegralnd *intnd, gdouble reltol)
{
  NcmIntegralndPrivate * const self = intnd->priv;

  self->reltol = reltol;
}

/**
 * ncm_integralnd_set_abstol:
 * @intnd: a #NcmIntegralnd
 * @abstol: absolute tolerance
 *
 * Sets the absolute tolerance @reltol to use.
 *
 */
void
ncm_integralnd_set_abstol (NcmIntegralnd *intnd, gdouble abstol)
{
  NcmIntegralndPrivate * const self = intnd->priv;

  self->abstol = abstol;
}

/**
 * ncm_integralnd_get_method:
 * @intnd: a #NcmIntegralnd
 *
 * Returns: the integration method used.
 */
NcmIntegralndMethod
ncm_integralnd_get_method (NcmIntegralnd *intnd)
{
  NcmIntegralndPrivate * const self = intnd->priv;

  return self->method;
}

/**
 * ncm_integralnd_get_error:
 * @intnd: a #NcmIntegralnd
 *
 * Returns: the error measure used.
 */
NcmIntegralndError
ncm_integralnd_get_error (NcmIntegralnd *intnd)
{
  NcmIntegralndPrivate * const self = intnd->priv;

  return self->error;
}

/**
 * ncm_integralnd_get_maxeval:
 * @intnd: a #NcmIntegralnd
 *
 * Returns: the maximum number of function evaluations used.
 */
guint
ncm_integralnd_get_maxeval (NcmIntegralnd *intnd)
{
  NcmIntegralndPrivate * const self = intnd->priv;

  return self->maxeval;
}

/**
 * ncm_integralnd_get_reltol:
 * @intnd: a #NcmIntegralnd
 *
 * Returns: the relative tolerance used.
 */
gdouble
ncm_integralnd_get_reltol (NcmIntegralnd *intnd)
{
  NcmIntegralndPrivate * const self = intnd->priv;

  return self->reltol;
}

/**
 * ncm_integralnd_get_abstol:
 * @intnd: a #NcmIntegralnd
 *
 * Returns: the absolute tolerance used.
 */
gdouble
ncm_integralnd_get_abstol (NcmIntegralnd *intnd)
{
  NcmIntegralndPrivate * const self = intnd->priv;

  return self->abstol;
}

static gint
_ncm_integralnd_cubature_int (unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
  NcmIntegralnd *intnd = NCM_INTEGRALND (fdata);
  NcmVector *x_vec     = ncm_vector_new_data_static ((gdouble *) x, ndim, 1);
  NcmVector *fval_vec  = ncm_vector_new_data_static (fval, fdim, 1);


  NCM_INTEGRALND_GET_CLASS (intnd)->integrand (intnd, x_vec, ndim, 1, fdim, fval_vec);

  ncm_vector_free (x_vec);
  ncm_vector_free (fval_vec);

  return 0;
}

static gint
_ncm_integralnd_cubature_vint (unsigned ndim, size_t npt, const double *x, void *fdata, unsigned fdim, double *fval)
{
  NcmIntegralnd *intnd = NCM_INTEGRALND (fdata);
  NcmVector *x_vec     = ncm_vector_new_data_static ((gdouble *) x, ndim * npt, 1);
  NcmVector *fval_vec  = ncm_vector_new_data_static (fval, fdim * npt, 1);

  NCM_INTEGRALND_GET_CLASS (intnd)->integrand (intnd, x_vec, ndim, npt, fdim, fval_vec);

  ncm_vector_free (x_vec);
  ncm_vector_free (fval_vec);

  return 0;
}

/**
 * ncm_integralnd_eval:
 * @intnd: a #NcmIntegralnd
 * @xi: a #NcmVector containing the inferior integration limit $x_i$
 * @xf: a #NcmVector containing the superior integration limit $x_f$
 * @res: a #NcmVector containing the result of the integration
 * @err: a #NcmVector containing the error of the integration
 *
 * Evaluated the integral $I_F(x_i, x_f) = \int_{x_i}^{x_f}F(x)\mathrm{d}x$.
 *
 */
void
ncm_integralnd_eval (NcmIntegralnd *intnd, const NcmVector *xi, const NcmVector *xf, NcmVector *res, NcmVector *err)
{
  NcmIntegralndPrivate * const self = intnd->priv;
  gint error                        = 0;

  guint dim, fdim;
  gint ret;


  NCM_INTEGRALND_GET_CLASS (intnd)->get_dimensions (intnd, &dim, &fdim);

  g_assert_cmpuint (ncm_vector_len (xi), ==, dim);
  g_assert_cmpuint (ncm_vector_len (xf), ==, dim);
  g_assert_cmpuint (ncm_vector_len (res), ==, fdim);
  g_assert_cmpuint (ncm_vector_len (err), ==, fdim);

  switch (self->error)
  {
    case NCM_INTEGRALND_ERROR_INDIVIDUAL:
      error = ERROR_INDIVIDUAL;
      break;
    case NCM_INTEGRALND_ERROR_PAIRWISE:
      error = ERROR_PAIRED;
      break;
    case NCM_INTEGRALND_ERROR_L2:
      error = ERROR_L2;
      break;
    case NCM_INTEGRALND_ERROR_L1:
      error = ERROR_L1;
      break;
    case NCM_INTEGRALND_ERROR_LINF:
      error = ERROR_LINF;
      break;
    default:
      g_error ("ncm_integralnd_eval: invalid error measure: `%d`.", self->error);
      break;
  }


  switch (self->method)
  {
    case NCM_INTEGRALND_METHOD_CUBATURE_H:
      ret = hcubature (
        fdim,
        _ncm_integralnd_cubature_int,
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
    case NCM_INTEGRALND_METHOD_CUBATURE_P:
      ret = pcubature (
        fdim,
        _ncm_integralnd_cubature_int,
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
    case NCM_INTEGRALND_METHOD_CUBATURE_H_V:
      ret = hcubature_v (
        fdim,
        _ncm_integralnd_cubature_vint,
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
    case NCM_INTEGRALND_METHOD_CUBATURE_P_V:
      ret = pcubature_v (
        fdim,
        _ncm_integralnd_cubature_vint,
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
      g_error ("ncm_integralnd_eval: invalid method: `%d`.", self->method);
      break;
  }

  g_assert (ret == 0);
}

