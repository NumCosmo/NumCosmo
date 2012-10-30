/***************************************************************************
 *            nc_matter_var.h
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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

#ifndef _NC_MATTER_VAR_H_
#define _NC_MATTER_VAR_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/math/ncm_model_ctrl.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/lss/nc_window.h>
#include <numcosmo/lss/nc_transfer_func.h>

G_BEGIN_DECLS

#define NC_TYPE_MATTER_VAR             (nc_matter_var_get_type ())
#define NC_MATTER_VAR(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_MATTER_VAR, NcMatterVar))
#define NC_MATTER_VAR_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_MATTER_VAR, NcMatterVarClass))
#define NC_IS_MATTER_VAR(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_MATTER_VAR))
#define NC_IS_MATTER_VAR_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_MATTER_VAR))
#define NC_MATTER_VAR_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_MATTER_VAR, NcMatterVarClass))

typedef struct _NcMatterVarClass NcMatterVarClass;
typedef struct _NcMatterVar NcMatterVar;

/**
 * NcMatterVarStrategy:
 * @NC_MATTER_VAR_NUMINT: Compute variance with numerical integration.
 * @NC_MATTER_VAR_SPLINEINT: Compute variance with spline.
 * @NC_MATTER_VAR_FFT: Compute using fft.
 *
 * FIXME	 
 */	
typedef enum _NcMatterVarStrategy
{
  NC_MATTER_VAR_NUMINT,
  NC_MATTER_VAR_SPLINEINT,
  NC_MATTER_VAR_FFT,
} NcMatterVarStrategy;

struct _NcMatterVarClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcMatterVar
{
  /*< private >*/
  GObject parent_instance;
  NcWindow *wp;            
  NcTransferFunc *tf;  
  NcMatterVarStrategy vs;   
  gdouble r;                
  gint spline_init;          
  gint n_points_spline;      
  NcmSpline *integrand_overw2_spline; 
  NcmSpline *sigma2_over_growth; 
  NcmSpline *deriv_sigma2_over_growth;
  NcmModelCtrl *ctrl;
};

GType nc_matter_var_get_type (void) G_GNUC_CONST;

NcMatterVar *nc_matter_var_new (NcMatterVarStrategy vs, NcWindow *wp, NcTransferFunc *tf);
NcMatterVar *nc_matter_var_copy (NcMatterVar *vp);
void nc_matter_var_free (NcMatterVar *vp);
void nc_matter_var_prepare (NcMatterVar *vp, NcHICosmo *model);
gdouble nc_matter_var_var0 (NcMatterVar *vp, NcHICosmo *model, gdouble lnR);
gdouble nc_matter_var_dlnvar0_dR (NcMatterVar *vp, NcHICosmo *model, gdouble lnR);
gdouble nc_matter_var_dlnvar0_dlnR (NcMatterVar *vp, NcHICosmo *model, gdouble lnR); 
gdouble nc_matter_var_mass_to_R (NcMatterVar *vp, NcHICosmo *model, gdouble M);
gdouble nc_matter_var_R_to_mass (NcMatterVar *vp, NcHICosmo *model, gdouble R);
gdouble nc_matter_var_lnM_to_lnR (NcMatterVar *vp, NcHICosmo *model, gdouble lnM);
gdouble nc_matter_var_lnR_to_lnM (NcMatterVar *vp, NcHICosmo *model, gdouble lnR);
gdouble nc_matter_var_integrand_over_window2 (NcMatterVar *vp, NcHICosmo *model, gdouble k);

gdouble nc_matter_var_spectral_moment_over_growth2 (NcMatterVar *vp, NcHICosmo *model, gint n);
gdouble nc_matter_var_spectral_moment_over_growth2_tophat (NcMatterVar *vp, NcHICosmo *model, gint n);
gdouble nc_matter_var_spectral_moment_over_growth2_gaussian (NcMatterVar *vp, NcHICosmo *model, gint n);

gdouble nc_matter_var_dsigma0_dR (NcMatterVar *vp, NcHICosmo *model, gdouble lnR);
gdouble nc_matter_var_sigma8_sqrtvar0 (NcMatterVar *vp, NcHICosmo *model);
/*
NcmFunc *nc_matter_var_new_variance_over_growth2 (NcMatterVar *vp);
NcmFunc *nc_matter_var_new_dvariance_over_growth2_dR (NcMatterVar *vp);
NcmFuncConst *nc_matter_var_new_variance_normalization (NcMatterVar *vp);
*/

G_END_DECLS

#endif /* _NC_MATTER_VAR_H_ */
