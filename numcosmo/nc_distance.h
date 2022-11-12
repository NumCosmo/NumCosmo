/***************************************************************************
 *            nc_distance.h
 *
 *  Tue Mar 18 13:36:26 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/

/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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

#ifndef _NC_DISTANCE_H_
#define _NC_DISTANCE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/nc_recomb.h>
#include <numcosmo/math/ncm_ode_spline.h>
#include <numcosmo/math/ncm_model_ctrl.h>
#include <numcosmo/math/ncm_function_cache.h>

G_BEGIN_DECLS

#define NC_TYPE_DISTANCE             (nc_distance_get_type ())
#define NC_DISTANCE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DISTANCE, NcDistance))
#define NC_DISTANCE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DISTANCE, NcDistanceClass))
#define NC_IS_DISTANCE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DISTANCE))
#define NC_IS_DISTANCE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DISTANCE))
#define NC_DISTANCE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DISTANCE, NcDistanceClass))

typedef struct _NcDistanceClass NcDistanceClass;
typedef struct _NcDistance NcDistance;
typedef gdouble (*NcDistanceFunc0) (NcDistance *dist, NcHICosmo *cosmo);
typedef gdouble (*NcDistanceFunc1) (NcDistance *dist, NcHICosmo *cosmo, gdouble z);

struct _NcDistanceClass
{
  /*< private >*/
  GObjectClass parent_class;
};

/**
 * NcDistanceComovingMethod:
 * @NC_DISTANCE_COMOVING_METHOD_INT_E: performs a numerical evaluation.
 * @NC_DISTANCE_COMOVING_METHOD_FROM_MODEL: uses the method defined by the implementation of #NcHICosmo.
 *
 * Enumeration to define which method to be applied in order to compute the cosmological distances.
 * 
 */
typedef enum _NcDistanceComovingMethod
{
  NC_DISTANCE_COMOVING_METHOD_INT_E = 0,
  NC_DISTANCE_COMOVING_METHOD_FROM_MODEL,
  /* < private > */
  NC_DISTANCE_COMOVING_METHOD_LEN,   /*< skip >*/  
} NcDistanceComovingMethod;

struct _NcDistance
{
  /*< private >*/
  GObject parent_instance;
  NcmOdeSpline *comoving_distance_spline;
  NcmSpline *inv_comoving_dist;
  NcmFunctionCache *comoving_distance_cache;
	NcmFunctionCache *comoving_infinity;
  NcmFunctionCache *time_cache;
  NcmFunctionCache *lookback_time_cache;
  NcmFunctionCache *conformal_time_cache;
  NcmFunctionCache *sound_horizon_cache;
  NcmModelCtrl *ctrl;
  gdouble zf;
  gboolean use_cache;
  gboolean cpu_inv_comoving;
  NcRecomb *recomb;
  NcDistanceComovingMethod cmethod;
};

typedef struct _NcDistanceFunc
{
  const gchar *name;
  const gchar *desc;
  NcDistanceFunc0 f;
  NcHICosmoImpl impl;
} NcDistanceFunc;

typedef struct _NcDistanceFuncZ
{
  const gchar *name;

  
  const gchar *desc;
  NcDistanceFunc1 f;
  NcHICosmoImpl impl;
} NcDistanceFuncZ;

GType nc_distance_get_type (void) G_GNUC_CONST;

NcDistance *nc_distance_new (gdouble zf);
NcDistance *nc_distance_ref (NcDistance *dist);

void nc_distance_require_zf (NcDistance *dist, const gdouble zf);
void nc_distance_set_recomb (NcDistance *dist, NcRecomb *recomb);
void nc_distance_compute_inv_comoving (NcDistance *dist, gboolean cpu_inv_xi);

void nc_distance_prepare (NcDistance *dist, NcHICosmo *cosmo);
NCM_INLINE void nc_distance_prepare_if_needed (NcDistance *dist, NcHICosmo *cosmo);

void nc_distance_free (NcDistance *dist);
void nc_distance_clear (NcDistance **dist);

/***************************************************************************
 * Redshift independent 'distances'
 ****************************************************************************/

gdouble nc_distance_hubble (NcDistance *dist, NcHICosmo *cosmo);
gdouble nc_distance_decoupling_redshift (NcDistance *dist, NcHICosmo *cosmo);
gdouble nc_distance_drag_redshift (NcDistance *dist, NcHICosmo *cosmo);
gdouble nc_distance_shift_parameter_lss (NcDistance *dist, NcHICosmo *cosmo);
gdouble nc_distance_comoving_lss (NcDistance *dist, NcHICosmo *cosmo);
gdouble nc_distance_acoustic_scale (NcDistance *dist, NcHICosmo *cosmo);
gdouble nc_distance_theta100CMB (NcDistance *dist, NcHICosmo *cosmo);
gdouble nc_distance_angular_diameter_curvature_scale (NcDistance *dist, NcHICosmo *cosmo);
gdouble nc_distance_r_zd (NcDistance *dist, NcHICosmo *cosmo);
gdouble nc_distance_r_zd_Mpc (NcDistance *dist, NcHICosmo *cosmo);

/***************************************************************************
 * Redshift dependent 'distances'
 ****************************************************************************/

gdouble nc_distance_comoving (NcDistance *dist, NcHICosmo *cosmo, const gdouble z);
gdouble nc_distance_transverse (NcDistance *dist, NcHICosmo *cosmo, const gdouble z);
gdouble nc_distance_dtransverse_dz (NcDistance *dist, NcHICosmo *cosmo, const gdouble z);
gdouble nc_distance_luminosity (NcDistance *dist, NcHICosmo *cosmo, const gdouble z);
gdouble nc_distance_angular_diameter (NcDistance *dist, NcHICosmo *cosmo, const gdouble z);
gdouble nc_distance_dmodulus (NcDistance *dist, NcHICosmo *cosmo, const gdouble z);
gdouble nc_distance_luminosity_hef (NcDistance *dist, NcHICosmo *cosmo, const gdouble z_he, const gdouble z_cmb);
gdouble nc_distance_dmodulus_hef (NcDistance *dist, NcHICosmo *cosmo, const gdouble z_he, const gdouble z_cmb);
gdouble nc_distance_shift_parameter (NcDistance *dist, NcHICosmo *cosmo, const gdouble z);
gdouble nc_distance_dilation_scale (NcDistance *dist, NcHICosmo *cosmo, const gdouble z);
gdouble nc_distance_bao_A_scale (NcDistance *dist, NcHICosmo *cosmo, const gdouble z);
gdouble nc_distance_sound_horizon (NcDistance *dist, NcHICosmo *cosmo, const gdouble z);
gdouble nc_distance_dsound_horizon_dz (NcDistance *dist, NcHICosmo *cosmo, const gdouble z);
gdouble nc_distance_bao_r_Dv (NcDistance *dist, NcHICosmo *cosmo, const gdouble z);
gdouble nc_distance_DH_r (NcDistance *dist, NcHICosmo *cosmo, const gdouble z);
gdouble nc_distance_DA_r (NcDistance *dist, NcHICosmo *cosmo, const gdouble z);
gdouble nc_distance_Dt_r (NcDistance *dist, NcHICosmo *cosmo, const gdouble z);
gdouble nc_distance_comoving_z_to_infinity (NcDistance *dist, NcHICosmo *cosmo, const gdouble z);
gdouble nc_distance_transverse_z_to_infinity (NcDistance *dist, NcHICosmo *cosmo, const gdouble z);

gdouble nc_distance_transverse_z1_z2 (NcDistance *dist, NcHICosmo *cosmo, const gdouble z1, const gdouble z2);
gdouble nc_distance_angular_diameter_z1_z2 (NcDistance *dist, NcHICosmo *cosmo, const gdouble z1, const gdouble z2);

/***************************************************************************
 * Inverse 'distances'
 ****************************************************************************/
gdouble nc_distance_inv_comoving (NcDistance *dist, NcHICosmo *cosmo, gdouble xi);

/***************************************************************************
 *            cosmic_time.h
 *
 *  Wed Nov 12 17:06:43 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/

gdouble nc_distance_cosmic_time (NcDistance *dist, NcHICosmo *cosmo, const gdouble z);
gdouble nc_distance_lookback_time (NcDistance *dist, NcHICosmo *cosmo, const gdouble z);
gdouble nc_distance_conformal_time (NcDistance *dist, NcHICosmo *cosmo, const gdouble z);
gdouble nc_distance_conformal_lookback_time (NcDistance *dist, NcHICosmo *cosmo, const gdouble z);

G_END_DECLS

#endif /* _NC_DISTANCE_INLINE_H_ */

#ifndef _NC_DISTANCE_INLINE_H_
#define _NC_DISTANCE_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE void
nc_distance_prepare_if_needed (NcDistance *dist, NcHICosmo *cosmo)
{
  if (ncm_model_ctrl_update (dist->ctrl, NCM_MODEL (cosmo)))
    nc_distance_prepare (dist, cosmo);
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_DISTANCE_INLINE_H_ */

