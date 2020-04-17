/***************************************************************************
 *            nc_halo_mass_function.h
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

#ifndef _NC_HALO_MASS_FUNCTION_H_
#define _NC_HALO_MASS_FUNCTION_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/lss/nc_multiplicity_func.h>
#include <numcosmo/math/ncm_spline2d.h>
#include <numcosmo/math/ncm_powspec_filter.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_MASS_FUNCTION             (nc_halo_mass_function_get_type ())
#define NC_HALO_MASS_FUNCTION(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HALO_MASS_FUNCTION, NcHaloMassFunction))
#define NC_HALO_MASS_FUNCTION_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HALO_MASS_FUNCTION, NcHaloMassFunctionClass))
#define NC_IS_HALO_MASS_FUNCTION(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HALO_MASS_FUNCTION))
#define NC_IS_HALO_MASS_FUNCTION_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HALO_MASS_FUNCTION))
#define NC_HALO_MASS_FUNCTION_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HALO_MASS_FUNCTION, NcHaloMassFunctionClass))
typedef struct _NcHaloMassFunctionClass NcHaloMassFunctionClass;
typedef struct _NcHaloMassFunction NcHaloMassFunction;

struct _NcHaloMassFunctionClass
{
  /*< private > */
  GObjectClass parent_class;
};

struct _NcHaloMassFunction
{
  /*< private > */
  GObject parent_instance;
  NcDistance *dist;
  NcMultiplicityFunc *mulf;
  NcmPowspecFilter *psf;
  gdouble area_survey;
  NcmSpline2d *d2NdzdlnM;
  gdouble lnMi;
  gdouble lnMf;
  gdouble zi;
  gdouble zf;
  gdouble prec;
  NcmModelCtrl *ctrl_cosmo;
  NcmModelCtrl *ctrl_reion;
};

/**
 * NcHaloMassFunctionSplineOptimize:
 * @NC_HALO_MASS_FUNCTION_SPLINE_NONE: FIXME
 * @NC_HALO_MASS_FUNCTION_SPLINE_LNM: FIXME
 * @NC_HALO_MASS_FUNCTION_SPLINE_Z: FIXME
 *
 * FIXME
 *
 */
typedef enum _NcHaloMassFunctionSplineOptimize
{
  NC_HALO_MASS_FUNCTION_SPLINE_NONE = 0,
  NC_HALO_MASS_FUNCTION_SPLINE_LNM,
  NC_HALO_MASS_FUNCTION_SPLINE_Z,
} NcHaloMassFunctionSplineOptimize;

GType nc_halo_mass_function_get_type (void) G_GNUC_CONST;

NcHaloMassFunction *nc_halo_mass_function_new (NcDistance *dist, NcmPowspecFilter *psf, NcMultiplicityFunc *mulf);
void nc_halo_mass_function_free (NcHaloMassFunction *mfp);
void nc_halo_mass_function_clear (NcHaloMassFunction **mfp);

void nc_halo_mass_function_set_area (NcHaloMassFunction *mfp, gdouble area);
void nc_halo_mass_function_set_prec (NcHaloMassFunction *mfp, gdouble prec);
void nc_halo_mass_function_set_area_sd (NcHaloMassFunction *mfp, gdouble area_sd);
void nc_halo_mass_function_set_eval_limits (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble lnMi, gdouble lnMf, gdouble zi, gdouble zf);
void nc_halo_mass_function_prepare (NcHaloMassFunction *mfp, NcHICosmo *cosmo);
NCM_INLINE void nc_halo_mass_function_prepare_if_needed (NcHaloMassFunction *mfp, NcHICosmo *cosmo);

gdouble nc_halo_mass_function_lnM_to_lnR (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble lnM);
gdouble nc_halo_mass_function_lnR_to_lnM (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble lnR);

void nc_halo_mass_function_dn_dlnR_sigma (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble lnR, gdouble z, gdouble *sigma_ptr, gdouble *dn_dlnR_ptr);
void nc_halo_mass_function_dn_dlnM_sigma (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *sigma_ptr, gdouble *dn_dlnM_ptr);

gdouble nc_halo_mass_function_dn_dlnR (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble lnR, gdouble z);
gdouble nc_halo_mass_function_dn_dlnM (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble lnM, gdouble z);

gdouble nc_halo_mass_function_dv_dzdomega (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble z);
NCM_INLINE gdouble nc_halo_mass_function_d2n_dzdlnM (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble lnM, gdouble z);
gdouble nc_halo_mass_function_dn_dz (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble lnMl, gdouble lnMu, gdouble z, gboolean spline);
gdouble nc_halo_mass_function_n (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble lnMl, gdouble lnMu, gdouble zl, gdouble zu, NcHaloMassFunctionSplineOptimize spline);

G_END_DECLS

#endif /* _NC_HALO_MASS_FUNCTION_H_ */

#ifndef _NC_HALO_MASS_FUNCTION_INLINE_H_
#define _NC_HALO_MASS_FUNCTION_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE void
nc_halo_mass_function_prepare_if_needed (NcHaloMassFunction *mfp, NcHICosmo *cosmo)
{
  gboolean cosmo_up = ncm_model_ctrl_update (mfp->ctrl_cosmo, NCM_MODEL (cosmo));
  
  if (cosmo_up)
	  nc_halo_mass_function_prepare (mfp, cosmo);
}

NCM_INLINE gdouble
nc_halo_mass_function_d2n_dzdlnM (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble lnM, gdouble z)
{
  return ncm_spline2d_eval (mfp->d2NdzdlnM, lnM, z);
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_HALO_MASS_FUNCTION_INLINE_H_ */
