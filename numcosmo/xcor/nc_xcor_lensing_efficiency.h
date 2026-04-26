/***************************************************************************
 *            nc_xcor_lensing_efficiency.h
 *
 *  Thu January 02 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_xcor_lensing_efficiency.h
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

#ifndef _NC_XCOR_LENSING_EFFICIENCY_H_
#define _NC_XCOR_LENSING_EFFICIENCY_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/math/ncm_model_ctrl.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_LENSING_EFFICIENCY (nc_xcor_lensing_efficiency_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcXcorLensingEfficiency, nc_xcor_lensing_efficiency, NC, XCOR_LENSING_EFFICIENCY, GObject)

/**
 * NcXcorLensingEfficiencyEvalSource:
 * @lens_eff: a #NcXcorLensingEfficiency
 * @z: redshift
 *
 * Evaluates the source weight function $W_{\mathrm{src}}(z)$ for the lensing
 * efficiency computation. This is used in the integral:
 * \begin{equation}
 * g(z) = \int_z^{z_{\max}} dz' \left(1 - \frac{\chi(z)}{\chi(z')}\right) W_{\mathrm{src}}(z')
 * \end{equation}
 *
 * Returns: $W_{\mathrm{src}}(z)$
 */
typedef gdouble (*NcXcorLensingEfficiencyEvalSource) (NcXcorLensingEfficiency *lens_eff, gdouble z);

/**
 * NcXcorLensingEfficiencyGetZRange:
 * @lens_eff: a #NcXcorLensingEfficiency
 * @zmin: (out): minimum source redshift
 * @zmax: (out): maximum source redshift
 *
 * Returns the redshift range for the source distribution used in lensing efficiency
 * computation.
 *
 */
typedef void (*NcXcorLensingEfficiencyGetZRange) (NcXcorLensingEfficiency *lens_eff, gdouble *zmin, gdouble *zmax);

struct _NcXcorLensingEfficiencyClass
{
  /*< private >*/
  GObjectClass parent_class;

  NcXcorLensingEfficiencyEvalSource eval_source;
  NcXcorLensingEfficiencyGetZRange get_z_range;

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[16];
};

NcXcorLensingEfficiency *nc_xcor_lensing_efficiency_ref (NcXcorLensingEfficiency *lens_eff);
void nc_xcor_lensing_efficiency_free (NcXcorLensingEfficiency *lens_eff);
void nc_xcor_lensing_efficiency_clear (NcXcorLensingEfficiency **lens_eff);

void nc_xcor_lensing_efficiency_set_distance (NcXcorLensingEfficiency *lens_eff, NcDistance *dist);
NcDistance *nc_xcor_lensing_efficiency_peek_distance (NcXcorLensingEfficiency *lens_eff);

void nc_xcor_lensing_efficiency_set_reltol (NcXcorLensingEfficiency *lens_eff, gdouble reltol);
void nc_xcor_lensing_efficiency_set_abstol (NcXcorLensingEfficiency *lens_eff, gdouble abstol);

gdouble nc_xcor_lensing_efficiency_get_reltol (NcXcorLensingEfficiency *lens_eff);
gdouble nc_xcor_lensing_efficiency_get_abstol (NcXcorLensingEfficiency *lens_eff);

void nc_xcor_lensing_efficiency_prepare (NcXcorLensingEfficiency *lens_eff, NcHICosmo *cosmo);
gdouble nc_xcor_lensing_efficiency_eval (NcXcorLensingEfficiency *lens_eff, gdouble z);

#define NC_XCOR_LENSING_EFFICIENCY_DEFAULT_RELTOL (1.0e-13)
#define NC_XCOR_LENSING_EFFICIENCY_DEFAULT_ABSTOL (1.0e-50)

/**
 * NC_XCOR_LENSING_EFFICIENCY_DEFINE_TYPE:
 * @MODULE: the name of the module defining the type, all capitalized
 * @OBJ_NAME: the name of the type to define, all capitalized
 * @ModuleObjName: the name of the type to define, camel case
 * @module_obj_name: the name of the type to define, snake case
 * @method_eval_source: the name of the method that evaluates $W_{\mathrm{src}}(z)$
 * @method_get_z_range: the name of the method that returns the source redshift range
 * @user_data: the type of the user data
 *
 * A convenience macro to define a subclass of #NcXcorLensingEfficiency with a custom
 * user data type. This follows the same pattern as #NcmIntegralND for inline type
 * definition.
 *
 * Example:
 * |[<!-- language="C" -->
 * typedef struct _MyLensEffData { MyKernel *kernel; } MyLensEffData;
 *
 * static gdouble
 * _my_lens_eff_eval_source (NcXcorLensingEfficiency *lens_eff, gdouble z)
 * {
 *   MyLensEff *my_eff = MY_LENS_EFF (lens_eff);
 *   return my_eff->data.kernel->dndz (z);
 * }
 *
 * static void
 * _my_lens_eff_get_z_range (NcXcorLensingEfficiency *lens_eff, gdouble *zmin, gdouble *zmax)
 * {
 *   MyLensEff *my_eff = MY_LENS_EFF (lens_eff);
 *   *zmin = my_eff->data.kernel->zmin;
 *   *zmax = my_eff->data.kernel->zmax;
 * }
 *
 * NC_XCOR_LENSING_EFFICIENCY_DEFINE_TYPE (MY, LENS_EFF, MyLensEff, my_lens_eff,
 *                                         _my_lens_eff_eval_source,
 *                                         _my_lens_eff_get_z_range,
 *                                         MyLensEffData)
 * ]|
 */
#define NC_XCOR_LENSING_EFFICIENCY_DEFINE_TYPE(MODULE, OBJ_NAME, ModuleObjName, module_obj_name, method_eval_source, method_get_z_range, user_data) \
        G_DECLARE_FINAL_TYPE (ModuleObjName, module_obj_name, MODULE, OBJ_NAME, NcXcorLensingEfficiency)                                            \
        struct _ ## ModuleObjName { NcXcorLensingEfficiency parent_instance; user_data data; };                                                     \
        G_DEFINE_TYPE (ModuleObjName, module_obj_name, NC_TYPE_XCOR_LENSING_EFFICIENCY)                                                             \
        static void                                                                                                                                 \
        module_obj_name ## _init (ModuleObjName * lens_eff)                                                                                         \
        {                                                                                                                                           \
        }                                                                                                                                           \
        static void                                                                                                                                 \
        module_obj_name ## _finalize (GObject * object)                                                                                             \
        {                                                                                                                                           \
          G_OBJECT_CLASS (module_obj_name ## _parent_class)->finalize (object);                                                                     \
        }                                                                                                                                           \
        static void                                                                                                                                 \
        module_obj_name ## _class_init (ModuleObjName ## Class * klass)                                                                             \
        {                                                                                                                                           \
          NcXcorLensingEfficiencyClass *lens_eff_class = NC_XCOR_LENSING_EFFICIENCY_CLASS (klass);                                                  \
          GObjectClass *gobject_class                  = G_OBJECT_CLASS (klass);                                                                    \
          gobject_class->finalize     = &module_obj_name ## _finalize;                                                                              \
          lens_eff_class->eval_source = &method_eval_source;                                                                                        \
          lens_eff_class->get_z_range = &method_get_z_range;                                                                                        \
        }

G_END_DECLS

#endif /* _NC_XCOR_LENSING_EFFICIENCY_H_ */

