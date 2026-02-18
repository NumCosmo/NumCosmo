/***************************************************************************
 *            nc_xcor_kernel_component.h
 *
 *  Wed February 12 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_xcor_kernel_component.h
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_XCOR_KERNEL_COMPONENT_H_
#define _NC_XCOR_KERNEL_COMPONENT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/nc_hicosmo.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_KERNEL_COMPONENT (nc_xcor_kernel_component_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcXcorKernelComponent, nc_xcor_kernel_component, NC, XCOR_KERNEL_COMPONENT, GObject)

/**
 * NcXcorKernelComponentEvalKernel:
 * @comp: a #NcXcorKernelComponent
 * @cosmo: a #NcHICosmo
 * @xi: comoving distance
 * @k: wave number
 *
 * Evaluates the kernel function K(k, xi) at the given comoving distance and wave number.
 *
 * Returns: the value of K(k, xi)
 */
typedef gdouble (*NcXcorKernelComponentEvalKernel) (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble xi, gdouble k);

/**
 * NcXcorKernelComponentEvalPrefactor:
 * @comp: a #NcXcorKernelComponent
 * @cosmo: a #NcHICosmo
 * @k: wave number
 * @l: multipole
 *
 * Evaluates the prefactor that may depend on k and ℓ.
 *
 * Returns: the prefactor value
 */
typedef gdouble (*NcXcorKernelComponentEvalPrefactor) (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble k, gint l);

/**
 * NcXcorKernelComponentGetLimits:
 * @comp: a #NcXcorKernelComponent
 * @cosmo: a #NcHICosmo
 * @xi_min: (out): minimum comoving distance
 * @xi_max: (out): maximum comoving distance
 * @k_min: (out): minimum wave number
 * @k_max: (out): maximum wave number
 *
 * Gets the valid integration ranges for this component.
 */
typedef void (*NcXcorKernelComponentGetLimits) (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble *xi_min, gdouble *xi_max, gdouble *k_min, gdouble *k_max);

struct _NcXcorKernelComponentClass
{
  /*< private >*/
  GObjectClass parent_class;

  NcXcorKernelComponentEvalKernel eval_kernel;
  NcXcorKernelComponentEvalPrefactor eval_prefactor;
  NcXcorKernelComponentGetLimits get_limits;

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[15];
};

NcXcorKernelComponent *nc_xcor_kernel_component_ref (NcXcorKernelComponent *comp);
void nc_xcor_kernel_component_free (NcXcorKernelComponent *comp);
void nc_xcor_kernel_component_clear (NcXcorKernelComponent **comp);

void nc_xcor_kernel_component_set_epsilon (NcXcorKernelComponent *comp, gdouble epsilon);
gdouble nc_xcor_kernel_component_get_epsilon (NcXcorKernelComponent *comp);

void nc_xcor_kernel_component_set_ny (NcXcorKernelComponent *comp, guint ny);
guint nc_xcor_kernel_component_get_ny (NcXcorKernelComponent *comp);

void nc_xcor_kernel_component_set_max_iter (NcXcorKernelComponent *comp, guint max_iter);
guint nc_xcor_kernel_component_get_max_iter (NcXcorKernelComponent *comp);

void nc_xcor_kernel_component_set_tol (NcXcorKernelComponent *comp, gdouble tol);
gdouble nc_xcor_kernel_component_get_tol (NcXcorKernelComponent *comp);

void nc_xcor_kernel_component_prepare (NcXcorKernelComponent *comp, NcHICosmo *cosmo);

gdouble nc_xcor_kernel_component_eval_k_max (NcXcorKernelComponent *comp, gdouble y);
gdouble nc_xcor_kernel_component_eval_K_max (NcXcorKernelComponent *comp, gdouble y);
gdouble nc_xcor_kernel_component_eval_k_epsilon (NcXcorKernelComponent *comp, gdouble y);

gdouble nc_xcor_kernel_component_eval_kernel (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble xi, gdouble k);
gdouble nc_xcor_kernel_component_eval_prefactor (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble k, gint l);
void nc_xcor_kernel_component_get_limits (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble *xi_min, gdouble *xi_max, gdouble *k_min, gdouble *k_max);

#define NC_XCOR_KERNEL_COMPONENT_DEFAULT_EPSILON 1.0e-12

/**
 * NC_XCOR_KERNEL_COMPONENT_DEFINE_TYPE:
 * @MODULE: the name of the module defining the type, all capitalized
 * @OBJ_NAME: the name of the type to define, all capitalized
 * @ModuleObjName: the name of the type to define, camel case
 * @module_obj_name: the name of the type to define, snake case
 * @method_eval_kernel: the name of the method that evaluates K(k, xi)
 * @method_eval_prefactor: the name of the method that evaluates the prefactor
 * @method_get_limits: the name of the method that returns the valid ranges
 * @user_data: the type of the user data
 * @user_data_free: the name of the method that frees the user data
 *
 * A convenience macro to define a subclass of #NcXcorKernelComponent with a custom user data type.
 */
#define NC_XCOR_KERNEL_COMPONENT_DEFINE_TYPE(MODULE, OBJ_NAME, ModuleObjName, module_obj_name, method_eval_kernel, method_eval_prefactor, method_get_limits, user_data, user_data_free) \
        G_DECLARE_FINAL_TYPE (ModuleObjName, module_obj_name, MODULE, OBJ_NAME, NcXcorKernelComponent)                                                                                  \
        struct _ ## ModuleObjName { NcXcorKernelComponent parent_instance; user_data data; };                                                                                           \
        G_DEFINE_TYPE (ModuleObjName, module_obj_name, NC_TYPE_XCOR_KERNEL_COMPONENT)                                                                                                   \
        static void                                                                                                                                                                     \
        module_obj_name ## _init (ModuleObjName * comp)                                                                                                                                 \
        {                                                                                                                                                                               \
        }                                                                                                                                                                               \
        static void                                                                                                                                                                     \
        module_obj_name ## _finalize (GObject * object)                                                                                                                                 \
        {                                                                                                                                                                               \
          ModuleObjName *comp = MODULE ## _ ## OBJ_NAME (object);                                                                                                                       \
          user_data_free (&comp->data);                                                                                                                                                 \
          G_OBJECT_CLASS (module_obj_name ## _parent_class)->finalize (object);                                                                                                         \
        }                                                                                                                                                                               \
        static void module_obj_name ## _class_init (ModuleObjName ## Class * klass)                                                                                                     \
        {                                                                                                                                                                               \
          NcXcorKernelComponentClass *comp_class = NC_XCOR_KERNEL_COMPONENT_CLASS (klass);                                                                                              \
          GObjectClass *gobject_class            = G_OBJECT_CLASS (klass);                                                                                                              \
          gobject_class->finalize    = &module_obj_name ## _finalize;                                                                                                                   \
          comp_class->eval_kernel    = &method_eval_kernel;                                                                                                                             \
          comp_class->eval_prefactor = &method_eval_prefactor;                                                                                                                          \
          comp_class->get_limits     = &method_get_limits;                                                                                                                              \
        }                                                                                                                                                                               \

G_END_DECLS

#endif /* _NC_XCOR_KERNEL_COMPONENT_H_ */

