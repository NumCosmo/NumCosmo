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
 * Abstract class for integrals of vector-valued functions over hyperrectangles.
 *
 * Computes $\int F(\vec{x})\,\mathrm{d}^n x$ for $F: \mathbb{R}^n \to \mathbb{R}^m$ over the
 * box between two corners, with the adaptive cubature library of S. G. Johnson
 * (https://github.com/stevengj/cubature). A subclass implements the integrand, a
 * #NcmIntegralNDF, and the dimensions $(n, m)$, a #NcmIntegralNDGetDimensions; the macros
 * %NCM_INTEGRAL_ND_DEFINE_TYPE and %NCM_INTEGRAL_ND_DEFINE_TYPE_WITH_FREE define such a
 * subclass with a user-data field.
 *
 * #NcmIntegralND:method selects the h-adaptive or the p-adaptive algorithm, each with a
 * scalar or a vectorized integrand, and #NcmIntegralND:error how the error of a
 * vector-valued integral is measured. The integration stops when the error meets
 * #NcmIntegralND:reltol or #NcmIntegralND:abstol; reaching #NcmIntegralND:maxeval
 * evaluations before that aborts. When the p-adaptive algorithm fails it is retried with the h-adaptive one,
 * with at most %NCM_INTEGRAL_ND_RETRY_MAXEVAL evaluations when #NcmIntegralND:maxeval is
 * zero; a failure of that retry, or of the h-adaptive algorithm, aborts. The buffers
 * passed to the integrand belong to the object, so evaluation is not reentrant.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/integration/ncm_integral_nd.h"
#include "ncm_enum_types.h"
#include "ncm/core/ncm_c.h"
#include "ncm/core/ncm_cfg.h"

#include "external/misc/cubature.h"

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

  /* Views over the cubature buffers, repointed at each integrand call */
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

  /**
   * NcmIntegralND:method:
   *
   * The cubature algorithm, see #NcmIntegralNDMethod.
   */
  g_object_class_install_property (object_class,
                                   PROP_METHOD,
                                   g_param_spec_enum ("method",
                                                      NULL,
                                                      "Integration method",
                                                      NCM_TYPE_INTEGRAL_ND_METHOD,
                                                      NCM_INTEGRAL_ND_METHOD_CUBATURE_H,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmIntegralND:error:
   *
   * How the error of a vector-valued integral is measured, see #NcmIntegralNDError.
   */
  g_object_class_install_property (object_class,
                                   PROP_ERROR,
                                   g_param_spec_enum ("error",
                                                      NULL,
                                                      "Error measure",
                                                      NCM_TYPE_INTEGRAL_ND_ERROR,
                                                      NCM_INTEGRAL_ND_ERROR_INDIVIDUAL,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmIntegralND:maxeval:
   *
   * The maximum number of integrand evaluations, or zero for no limit.
   */
  g_object_class_install_property (object_class,
                                   PROP_MAXEVAL,
                                   g_param_spec_uint ("maxeval",
                                                      NULL,
                                                      "Maximum number of function evaluations (0 means unlimited)",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmIntegralND:reltol:
   *
   * The relative tolerance.
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Integral relative tolerance",
                                                        0.0, 1.0, NCM_INTEGRAL_ND_DEFAULT_RELTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmIntegralND:abstol:
   *
   * The absolute tolerance.
   */
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
 * If *@intnd is not %NULL, decreases its reference count by one and sets *@intnd to %NULL.
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
 * Sets #NcmIntegralND:method.
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
 * Sets #NcmIntegralND:error.
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
 * @maxeval: the maximum number of integrand evaluations
 *
 * Sets #NcmIntegralND:maxeval.
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
 * @reltol: the relative tolerance
 *
 * Sets #NcmIntegralND:reltol.
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
 * @abstol: the absolute tolerance
 *
 * Sets #NcmIntegralND:abstol.
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
 * Gets #NcmIntegralND:method.
 *
 * Returns: the cubature algorithm.
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
 * Gets #NcmIntegralND:error.
 *
 * Returns: the error measure.
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
 * Gets #NcmIntegralND:maxeval.
 *
 * Returns: the maximum number of integrand evaluations.
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
 * Gets #NcmIntegralND:reltol.
 *
 * Returns: the relative tolerance.
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
 * Gets #NcmIntegralND:abstol.
 *
 * Returns: the absolute tolerance.
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

static const gchar *
_ncm_integral_nd_method_name (NcmIntegralNDMethod method)
{
  switch (method)
  {
    case NCM_INTEGRAL_ND_METHOD_CUBATURE_H:
      return "hcubature";

    case NCM_INTEGRAL_ND_METHOD_CUBATURE_P:
      return "pcubature";

    case NCM_INTEGRAL_ND_METHOD_CUBATURE_H_V:
      return "hcubature_v";

    case NCM_INTEGRAL_ND_METHOD_CUBATURE_P_V:
      return "pcubature_v";

    default:                   /* LCOV_EXCL_LINE */
      return "unknown method"; /* LCOV_EXCL_LINE */
  }
}

/* Evaluations of the h-adaptive retry when maxeval is zero: unlimited, h-adaptive
 * subdivision never reports failure and grows until memory is exhausted. */
#define NCM_INTEGRAL_ND_RETRY_MAXEVAL (10000000)

/* The convergence test of cubature, applied to the final result: cubature returns
 * success also when it stops at maxeval without converging */
static gboolean
_ncm_integral_nd_converged (unsigned fdim, const double *val_v, const double *err_v, double reqAbsError, double reqRelError, error_norm norm)
#define ERR(j) err_v[j]
#define VAL(j) val_v[j]
#include "external/misc/converged.h"
#undef ERR
#undef VAL

static gboolean
_ncm_integral_nd_method_is_p (NcmIntegralNDMethod method)
{
  return (method == NCM_INTEGRAL_ND_METHOD_CUBATURE_P) ||
         (method == NCM_INTEGRAL_ND_METHOD_CUBATURE_P_V);
}

static NcmIntegralNDMethod
_ncm_integral_nd_method_h_of_p (NcmIntegralNDMethod method)
{
  return (method == NCM_INTEGRAL_ND_METHOD_CUBATURE_P) ?
         NCM_INTEGRAL_ND_METHOD_CUBATURE_H :
         NCM_INTEGRAL_ND_METHOD_CUBATURE_H_V;
}

/* Runs one cubature method */
static gint
_ncm_integral_nd_run (NcmIntegralND *intnd, NcmIntegralNDMethod method, guint maxeval, guint dim, guint fdim, gint error, const NcmVector *xi, const NcmVector *xf, NcmVector *res, NcmVector *err)
{
  NcmIntegralNDPrivate * const self = ncm_integral_nd_get_instance_private (intnd);
  gint ret                          = 0;

  switch (method)
  {
    case NCM_INTEGRAL_ND_METHOD_CUBATURE_H:
      ret = hcubature (
        fdim,
        _ncm_integral_nd_cubature_int,
        intnd,
        dim,
        ncm_vector_const_data (xi),
        ncm_vector_const_data (xf),
        maxeval,
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
        maxeval,
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
        maxeval,
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
        maxeval,
        self->abstol,
        self->reltol,
        error,
        ncm_vector_data (res),
        ncm_vector_data (err)
      );
      break;
    default:                                                           /* LCOV_EXCL_LINE */
      g_error ("ncm_integral_nd_eval: invalid method: `%d`.", method); /* LCOV_EXCL_LINE */
      break;                                                           /* LCOV_EXCL_LINE */
  }

  return ret;
}

/**
 * ncm_integral_nd_eval:
 * @intnd: a #NcmIntegralND
 * @xi: the lower corner, of length $n$
 * @xf: the upper corner, of length $n$
 * @res: the output integral, of length $m$
 * @err: the output error estimate, of length $m$
 *
 * Computes $\int_{\vec{x}_i}^{\vec{x}_f} F(\vec{x})\,\mathrm{d}^n x$ into @res, see
 * #NcmIntegralND.
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


  ret = _ncm_integral_nd_run (intnd, self->method, self->maxeval, dim, fdim, error, xi, xf, res, err);

  /* A p-adaptive failure (Clenshaw-Curtis levels exhausted) is retried h-adaptively,
   * which converges on integrands smooth only in parts of the domain, with a finite
   * budget, see NCM_INTEGRAL_ND_RETRY_MAXEVAL. */
  if ((ret != 0) && _ncm_integral_nd_method_is_p (self->method))
  {
    const NcmIntegralNDMethod fallback = _ncm_integral_nd_method_h_of_p (self->method);
    const guint retry_maxeval          = (self->maxeval == 0) ? NCM_INTEGRAL_ND_RETRY_MAXEVAL : self->maxeval;

    g_debug ("ncm_integral_nd_eval: %s failed (%d) on %s, retrying with %s (maxeval %u).",
             _ncm_integral_nd_method_name (self->method), ret,
             G_OBJECT_TYPE_NAME (intnd), _ncm_integral_nd_method_name (fallback),
             retry_maxeval);

    ret = _ncm_integral_nd_run (intnd, fallback, retry_maxeval, dim, fdim, error, xi, xf, res, err);
  }

  if ((ret == 0) && !_ncm_integral_nd_converged (fdim, ncm_vector_data (res), ncm_vector_data (err),
                                                 self->abstol, self->reltol, (fdim <= 1) ? ERROR_INDIVIDUAL : error))
    g_error ("ncm_integral_nd_eval: %s on %s stopped at maxeval %u without reaching reltol %.17g "
             "or abstol %.17g.",
             _ncm_integral_nd_method_name (self->method), G_OBJECT_TYPE_NAME (intnd),
             self->maxeval, self->reltol, self->abstol);

  if (ret != 0)
  {
    GString *bounds = g_string_new (NULL);
    guint i;

    for (i = 0; i < dim; i++)
      g_string_append_printf (bounds, "%s[% 22.15g, % 22.15g]", (i > 0) ? ", " : "",
                              ncm_vector_get (xi, i), ncm_vector_get (xf, i));

    g_error ("ncm_integral_nd_eval: %s failed (%d) on %s integrating %u dimension(s) "
             "over %s to %u component(s), reltol %.17g, abstol %.17g, maxeval %u. "
             "%s",
             _ncm_integral_nd_method_name (self->method), ret,
             G_OBJECT_TYPE_NAME (intnd), dim, bounds->str, fdim,
             self->reltol, self->abstol, self->maxeval,
             _ncm_integral_nd_method_is_p (self->method) ?
             "The p-adaptive method ran out of Clenshaw-Curtis levels before reaching "
             "the requested tolerance, which is what happens when the integrand is not "
             "accurate or smooth enough to support the tolerance asked of it, and the "
             "h-adaptive retry did not converge either. Loosen reltol to match the "
             "accuracy the integrand actually carries." : "");

    g_string_free (bounds, TRUE);
  }
}

