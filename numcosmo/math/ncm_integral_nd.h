/***************************************************************************
 *            ncm_integral_nd.h
 *
 *  Thu July 20 08:39:30 2023
 *  Copyright  2023 Eduardo José Barroso
 *  <eduardo.jsbarroso@uel.br>
 ****************************************************************************/
/*
 * ncm_integral_nd.h
 * Copyright (C) 2023 Eduardo José Barroso <eduardo.jsbarroso@uel.br>*
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

#ifndef _NCM_INTEGRAL_ND_H_
#define _NCM_INTEGRAL_ND_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

#include <numcosmo/math/ncm_vector.h>

G_BEGIN_DECLS

#define NCM_TYPE_INTEGRAL_ND (ncm_integral_nd_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmIntegralND, ncm_integral_nd, NCM, INTEGRAL_ND, GObject)

/**
 * NcmIntegralNDF:
 * @intnd: a #NcmIntegralND
 * @x: a #NcmVector containing the value of the variable of integration
 * @dim: the dimension of the integral argument
 * @npoints: the number of points in the array @x
 * @fdim: the dimension of the function to be integrated
 * @fval: a #NcmVector containing the @fdim values of the integrand at all points in @x
 *
 * The type of the function that must be implemented by a subclass of #NcmIntegralND.
 *
 * This function receives @npoints points in the array @x (size @dim * @npoints), and
 * returns an array (size @fdim * @npoints) of @npoints values of the integrand at all
 * points in @x. The @x is an array of doubles in row-major order
 * (i.e. the first @dim elements of @x are the coordinates of the first point, the next
 * @dim elements are the coordinates of the second point, and so on). The return value
 * is an array of @fdim values of the integrand at all points in @x (e.g. the first
 * @fdim elements are the values of the integrand at the first point, the next @fdim
 * elements are the values of the integrand at the second point, and so on).
 *
 */
typedef void (*NcmIntegralNDF) (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);

/**
 * NcmIntegralNDGetDimensions:
 * @intnd: a #NcmIntegralND
 * @dim: (out): the dimension of the integral argument
 * @fdim: (out): the dimension of the function to be integrated
 *
 * The type of the function that must be implemented by a subclass of #NcmIntegralND.
 *
 * This function returns the dimension of the integral argument and the dimension of
 * the function to be integrated.
 *
 */
typedef void (*NcmIntegralNDGetDimensions) (NcmIntegralND *intnd, guint *dim, guint *fdim);

struct _NcmIntegralNDClass
{
  /*< private >*/
  GObjectClass parent_class;

  NcmIntegralNDF integrand;
  NcmIntegralNDGetDimensions get_dimensions;
  gpointer padding[16];
};

/**
 * NcmIntegralNDMethod:
 * @NCM_INTEGRAL_ND_METHOD_CUBATURE_H: adapative integration by partitioning the integration domain ("h-adaptive")
 *  and using the same fixed-degree quadrature in each subdomain, recursively, until convergence is achieved.
 * @NCM_INTEGRAL_ND_METHOD_CUBATURE_P: adaptive integration by increasing the degree of (tensor-product
 *  Clenshaw-Curtis) quadrature rules ("p-adaptive"), rather than subdividing the domain ("h-adaptive").
 *  Possibly better for smooth integrands in low dimensions.
 * @NCM_INTEGRAL_ND_METHOD_CUBATURE_H_V: same as @NCM_INTEGRAL_ND_METHOD_CUBATURE_H with vectorized integrand
 * @NCM_INTEGRAL_ND_METHOD_CUBATURE_P_V: same as @NCM_INTEGRAL_ND_METHOD_CUBATURE_P with vectorized integrand
 *
 * The type of the method used to perform the integral.
 *
 */
typedef enum _NcmIntegralNDMethod
{
  NCM_INTEGRAL_ND_METHOD_CUBATURE_H,
  NCM_INTEGRAL_ND_METHOD_CUBATURE_P,
  NCM_INTEGRAL_ND_METHOD_CUBATURE_H_V,
  NCM_INTEGRAL_ND_METHOD_CUBATURE_P_V,
  /* < private > */
  NCM_INTEGRAL_ND_METHOD_LEN, /*< skip >*/
} NcmIntegralNDMethod;

/**
 * NcmIntegralNDError:
 * @NCM_INTEGRAL_ND_ERROR_INDIVIDUAL: error is estimated for each integrand separately
 * @NCM_INTEGRAL_ND_ERROR_PAIRWISE: error is estimated for each pair of integrands
 * @NCM_INTEGRAL_ND_ERROR_L2: error is estimated for the L2 norm of the vector of integrands
 * @NCM_INTEGRAL_ND_ERROR_L1: error is estimated for the L1 norm of the vector of integrands
 * @NCM_INTEGRAL_ND_ERROR_LINF: error is estimated for the L-infinity norm of the vector of integrands
 *
 * The type of the error estimation used to perform the integral.
 *
 */
typedef enum _NcmIntegralNDError
{
  NCM_INTEGRAL_ND_ERROR_INDIVIDUAL,
  NCM_INTEGRAL_ND_ERROR_PAIRWISE,
  NCM_INTEGRAL_ND_ERROR_L2,
  NCM_INTEGRAL_ND_ERROR_L1,
  NCM_INTEGRAL_ND_ERROR_LINF,
  /* < private > */
  NCM_INTEGRAL_ND_ERROR_LEN, /*< skip >*/
} NcmIntegralNDError;

NcmIntegralND *ncm_integral_nd_ref (NcmIntegralND *intnd);
void ncm_integral_nd_free (NcmIntegralND *intnd);
void ncm_integral_nd_clear (NcmIntegralND **intnd);

void ncm_integral_nd_set_method (NcmIntegralND *intnd, NcmIntegralNDMethod method);
void ncm_integral_nd_set_error (NcmIntegralND *intnd, NcmIntegralNDError error);
void ncm_integral_nd_set_maxeval (NcmIntegralND *intnd, guint maxeval);
void ncm_integral_nd_set_reltol (NcmIntegralND *intnd, gdouble reltol);
void ncm_integral_nd_set_abstol (NcmIntegralND *intnd, gdouble abstol);

NcmIntegralNDMethod ncm_integral_nd_get_method (NcmIntegralND *intnd);
NcmIntegralNDError ncm_integral_nd_get_error (NcmIntegralND *intnd);
guint ncm_integral_nd_get_maxeval (NcmIntegralND *intnd);
gdouble ncm_integral_nd_get_reltol (NcmIntegralND *intnd);
gdouble ncm_integral_nd_get_abstol (NcmIntegralND *intnd);

void ncm_integral_nd_eval (NcmIntegralND *intnd, const NcmVector *xi, const NcmVector *xf, NcmVector *res, NcmVector *err);

#define NCM_INTEGRAL_ND_DEFAULT_RELTOL 1e-7
#define NCM_INTEGRAL_ND_DEFAULT_ABSTOL 0.0

#define NCM_INTEGRAL_ND_DEFINE_TYPE_WITH_FREE(MODULE, OBJ_NAME, ModuleObjName, module_obj_name, method_get_dimensions, method_integrand, user_data, user_data_free) \
  G_DECLARE_FINAL_TYPE (ModuleObjName, module_obj_name, MODULE, OBJ_NAME, NcmIntegralND) \
  struct _ ## ModuleObjName { NcmIntegralND parent_instance; user_data data; }; \
  G_DEFINE_TYPE (ModuleObjName, module_obj_name, NCM_TYPE_INTEGRAL_ND) \
  static void \
  module_obj_name ## _init (ModuleObjName * intnd) \
  { \
  } \
  static void \
  module_obj_name ## _finalize (GObject * object) \
  { \
    ModuleObjName *intnd = MODULE ## _ ## OBJ_NAME (object); \
    user_data_free (&intnd->data); \
    G_OBJECT_CLASS (module_obj_name ## _parent_class)->finalize (object); \
  } \
  static void module_obj_name ## _class_init (ModuleObjName ## Class * klass) \
  { \
    NcmIntegralNDClass *intnd_class = NCM_INTEGRAL_ND_CLASS (klass); \
    GObjectClass *gobject_class     = G_OBJECT_CLASS (klass); \
    gobject_class->finalize     = &module_obj_name ## _finalize; \
    intnd_class->get_dimensions = &method_get_dimensions; \
    intnd_class->integrand      = &method_integrand; \
  } \

#define NCM_INTEGRAL_ND_DEFINE_TYPE(MODULE, OBJ_NAME, ModuleObjName, module_obj_name, method_get_dimensions, method_integrand, user_data) \
  NCM_INTEGRAL_ND_DEFINE_TYPE_WITH_FREE (MODULE, OBJ_NAME, ModuleObjName, module_obj_name, method_get_dimensions, method_integrand, user_data, (void)) \

G_END_DECLS

#endif /* _NCM_INTEGRAL_ND_H_ */

