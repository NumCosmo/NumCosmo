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

#include <numcosmo/ncm/algebra/ncm_vector.h>

G_BEGIN_DECLS

#define NCM_TYPE_INTEGRAL_ND (ncm_integral_nd_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmIntegralND, ncm_integral_nd, NCM, INTEGRAL_ND, GObject)

/**
 * NcmIntegralNDF:
 * @intnd: a #NcmIntegralND
 * @x: the points, @npoints groups of @dim coordinates
 * @dim: the number of variables $n$
 * @npoints: the number of points
 * @fdim: the number of components $m$
 * @fval: the output values, @npoints groups of @fdim components
 *
 * The integrand of #NcmIntegralND, computing $F$ at each point of @x into @fval, point
 * after point. @npoints is one except for the vectorized methods.
 */
typedef void (*NcmIntegralNDF) (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);

/**
 * NcmIntegralNDGetDimensions:
 * @intnd: a #NcmIntegralND
 * @dim: (out): the number of variables $n$
 * @fdim: (out): the number of components $m$
 *
 * Gets the dimensions of the integrand of #NcmIntegralND.
 */
typedef void (*NcmIntegralNDGetDimensions) (NcmIntegralND *intnd, guint *dim, guint *fdim);

struct _NcmIntegralNDClass
{
  /*< private >*/
  GObjectClass parent_class;

  NcmIntegralNDF integrand;
  NcmIntegralNDGetDimensions get_dimensions;

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[16];
};

/**
 * NcmIntegralNDMethod:
 * @NCM_INTEGRAL_ND_METHOD_CUBATURE_H: h-adaptive, subdividing the domain with a
 *   fixed-degree rule in each part
 * @NCM_INTEGRAL_ND_METHOD_CUBATURE_P: p-adaptive, raising the degree of a tensor-product
 *   Clenshaw-Curtis rule on the whole domain, suited to smooth integrands in few dimensions
 * @NCM_INTEGRAL_ND_METHOD_CUBATURE_H_V: h-adaptive with a vectorized integrand
 * @NCM_INTEGRAL_ND_METHOD_CUBATURE_P_V: p-adaptive with a vectorized integrand
 *
 * The cubature algorithms of #NcmIntegralND.
 */
typedef enum _NcmIntegralNDMethod /*< prefix=NCM_INTEGRAL_ND_METHOD_CUBATURE >*/
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
 * @NCM_INTEGRAL_ND_ERROR_INDIVIDUAL: each component meets the tolerance
 * @NCM_INTEGRAL_ND_ERROR_PAIRWISE: each pair of consecutive components, taken as the real
 *   and imaginary parts of a complex value, meets the tolerance in modulus
 * @NCM_INTEGRAL_ND_ERROR_L2: the $L^2$ norm of the error vector meets the tolerance
 *   relative to the $L^2$ norm of the integral
 * @NCM_INTEGRAL_ND_ERROR_L1: as %NCM_INTEGRAL_ND_ERROR_L2 with the $L^1$ norm
 * @NCM_INTEGRAL_ND_ERROR_LINF: as %NCM_INTEGRAL_ND_ERROR_L2 with the maximum norm
 *
 * How #NcmIntegralND measures the error of a vector-valued integral.
 */
typedef enum _NcmIntegralNDError /*< prefix=NCM_INTEGRAL_ND_ERROR >*/
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

/**
 * NCM_INTEGRAL_ND_DEFAULT_RELTOL:
 *
 * Default #NcmIntegralND:reltol.
 */
#define NCM_INTEGRAL_ND_DEFAULT_RELTOL 1e-7

/**
 * NCM_INTEGRAL_ND_DEFAULT_ABSTOL:
 *
 * Default #NcmIntegralND:abstol.
 */
#define NCM_INTEGRAL_ND_DEFAULT_ABSTOL 0.0

/**
 * NCM_INTEGRAL_ND_DEFINE_TYPE_WITH_FREE:
 * @MODULE: the module prefix, in upper case
 * @OBJ_NAME: the type name without the prefix, in upper case
 * @ModuleObjName: the type name, in camel case
 * @module_obj_name: the type name, in snake case
 * @method_get_dimensions: the #NcmIntegralNDGetDimensions
 * @method_integrand: the #NcmIntegralNDF
 * @user_data: the type of the user-data field `data`
 * @user_data_free: the function called with a pointer to `data` on finalization
 *
 * Defines a final subclass of #NcmIntegralND whose instance holds a field `data` of type
 * @user_data.
 */
#define NCM_INTEGRAL_ND_DEFINE_TYPE_WITH_FREE(MODULE, OBJ_NAME, ModuleObjName, module_obj_name, method_get_dimensions, method_integrand, user_data, user_data_free) \
        G_DECLARE_FINAL_TYPE (ModuleObjName, module_obj_name, MODULE, OBJ_NAME, NcmIntegralND)                                                                      \
        struct _ ## ModuleObjName { NcmIntegralND parent_instance; user_data data; };                                                                               \
        G_DEFINE_TYPE (ModuleObjName, module_obj_name, NCM_TYPE_INTEGRAL_ND)                                                                                        \
        static void                                                                                                                                                 \
        module_obj_name ## _init (ModuleObjName * intnd)                                                                                                            \
        {                                                                                                                                                           \
        }                                                                                                                                                           \
        static void                                                                                                                                                 \
        module_obj_name ## _finalize (GObject * object)                                                                                                             \
        {                                                                                                                                                           \
          ModuleObjName *intnd = MODULE ## _ ## OBJ_NAME (object);                                                                                                  \
          user_data_free (&intnd->data);                                                                                                                            \
          G_OBJECT_CLASS (module_obj_name ## _parent_class)->finalize (object);                                                                                     \
        }                                                                                                                                                           \
        static void module_obj_name ## _class_init (ModuleObjName ## Class * klass)                                                                                 \
        {                                                                                                                                                           \
          NcmIntegralNDClass *intnd_class = NCM_INTEGRAL_ND_CLASS (klass);                                                                                          \
          GObjectClass *gobject_class     = G_OBJECT_CLASS (klass);                                                                                                 \
          gobject_class->finalize     = &module_obj_name ## _finalize;                                                                                              \
          intnd_class->get_dimensions = &method_get_dimensions;                                                                                                     \
          intnd_class->integrand      = &method_integrand;                                                                                                          \
        }                                                                                                                                                           \


/**
 * NCM_INTEGRAL_ND_DEFINE_TYPE:
 * @MODULE: the module prefix, in upper case
 * @OBJ_NAME: the type name without the prefix, in upper case
 * @ModuleObjName: the type name, in camel case
 * @module_obj_name: the type name, in snake case
 * @method_get_dimensions: the #NcmIntegralNDGetDimensions
 * @method_integrand: the #NcmIntegralNDF
 * @user_data: the type of the user-data field `data`
 *
 * Same as %NCM_INTEGRAL_ND_DEFINE_TYPE_WITH_FREE with nothing to free.
 */
#define NCM_INTEGRAL_ND_DEFINE_TYPE(MODULE, OBJ_NAME, ModuleObjName, module_obj_name, method_get_dimensions, method_integrand, user_data)                    \
        NCM_INTEGRAL_ND_DEFINE_TYPE_WITH_FREE (MODULE, OBJ_NAME, ModuleObjName, module_obj_name, method_get_dimensions, method_integrand, user_data, (void)) \

G_END_DECLS

#endif /* _NCM_INTEGRAL_ND_H_ */

