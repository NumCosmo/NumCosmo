/***************************************************************************
 *            nc_mass_function.h
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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

#ifndef _NC_MASS_FUNCTION_H_
#define _NC_MASS_FUNCTION_H_

#include <glib-object.h>

G_BEGIN_DECLS
#define NC_TYPE_MASS_FUNCTION             (nc_mass_function_get_type ())
#define NC_MASS_FUNCTION(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_MASS_FUNCTION, NcMassFunction))
#define NC_MASS_FUNCTION_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_MASS_FUNCTION, NcMassFunctionClass))
#define NC_IS_MASS_FUNCTION(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_MASS_FUNCTION))
#define NC_IS_MASS_FUNCTION_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_MASS_FUNCTION))
#define NC_MASS_FUNCTION_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_MASS_FUNCTION, NcMassFunctionClass))
typedef struct _NcMassFunctionClass NcMassFunctionClass;
typedef struct _NcMassFunction NcMassFunction;

struct _NcMassFunctionClass
{
  /*< private > */
  GObjectClass parent_class;
};

struct _NcMassFunction
{
  /*< private > */
  GObject parent_instance;
  NcDistance *dist;
  NcMatterVar *vp;
  NcGrowthFunc *gf;
  NcMultiplicityFunc *mulf;
  gdouble area_survey;
  gdouble N_sigma;
  gdouble growth;		/* Internal use only */
  NcmSpline2d *d2NdzdlnM;
  gdouble lnMi;
  gdouble lnMf;
  gdouble zi;
  gdouble zf;
  NcmModelCtrl *ctrl;
};

/**
 * NcMassFunctionSplineOptimize:
 * @NC_MASS_FUNCTION_SPLINE_NONE: FIXME
 * @NC_MASS_FUNCTION_SPLINE_LNM: FIXME
 * @NC_MASS_FUNCTION_SPLINE_Z: FIXME
 *
 * FIXME
 *
 */
typedef enum _NcMassFunctionSplineOptimize
{
  NC_MASS_FUNCTION_SPLINE_NONE = 0,
  NC_MASS_FUNCTION_SPLINE_LNM,
  NC_MASS_FUNCTION_SPLINE_Z,
} NcMassFunctionSplineOptimize;

GType
nc_mass_function_get_type (void)
G_GNUC_CONST;


NcMassFunction *nc_mass_function_new (NcDistance * dist, NcMatterVar * vp, NcGrowthFunc * gf, NcMultiplicityFunc * mulf);
NcMassFunction *nc_mass_function_copy (NcMassFunction * mfp);
void nc_mass_function_free (NcMassFunction * mfp);
gdouble nc_mass_function (NcMassFunction * mfp, NcHICosmo * model, gdouble lnM, gdouble z);
void nc_mass_function_sigma (NcMassFunction * mfp, NcHICosmo * model, gdouble lnM, gdouble z,
                             gdouble * dn_dlnM_ptr, gdouble * sigma_ptr);
void nc_mass_function_alpha_eff (NcMatterVar * vp, NcHICosmo * model, gdouble lnM, gdouble * a_eff_ptr);
gdouble nc_mass_function_dn_M_to_inf_dv (NcMassFunction * mfp, NcHICosmo * model, gdouble M, gdouble z);
gdouble nc_mass_function_dn_M1_to_M2_dv (NcMassFunction * mfp, NcHICosmo * model, gdouble M1, gdouble M2,
                                         gdouble z);
gdouble nc_mass_function_dcomoving_volume_dzdomega (NcMassFunction * mfp, NcHICosmo * model, gdouble z);
gdouble nc_mass_function_d2NdzdlnM (NcMassFunction *mfp, NcHICosmo *model, gdouble lnM, gdouble z);
gdouble nc_mass_function_dNdz (NcMassFunction *mfp, NcHICosmo *model, gdouble lnMl, gdouble lnMu, gdouble z,
                               gboolean spline);
gdouble nc_mass_function_N (NcMassFunction * mfp, NcHICosmo * model, gdouble lnMl, gdouble lnMu, gdouble zl,
                            gdouble zu, NcMassFunctionSplineOptimize spline);
void nc_mass_function_d2NdzdlnM_optimize (NcMassFunction * mfp, NcHICosmo * model, gdouble lnMi, gdouble lnMf,
                                          gdouble zi, gdouble zf);

G_END_DECLS
#endif /* _NC_MASS_FUNCTION_H_ */
