/***************************************************************************
 *            nc_distance.h
 *
 *  Tue Mar 18 13:36:26 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
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

#include <glib-object.h>

G_BEGIN_DECLS

#define NC_TYPE_DISTANCE             (nc_distance_get_type ())
#define NC_DISTANCE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DISTANCE, NcDistance))
#define NC_DISTANCE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DISTANCE, NcDistanceClass))
#define NC_IS_DISTANCE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DISTANCE))
#define NC_IS_DISTANCE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DISTANCE))
#define NC_DISTANCE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DISTANCE, NcDistanceClass))

typedef struct _NcDistanceClass NcDistanceClass;
typedef struct _NcDistance NcDistance;
typedef gdouble (*NcDistanceFunc0) (NcDistance *dist, NcHICosmo *model);
typedef gdouble (*NcDistanceFunc1) (NcDistance *dist, NcHICosmo *model, gdouble z);

struct _NcDistanceClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcDistance
{
  /*< private >*/
  GObject parent_instance;
  gpointer cache;
  NcmOdeSpline *comoving_distance_spline;
  NcFunctionCache *comoving_distance_cache;
  NcFunctionCache *time_cache;
  NcFunctionCache *lookback_time_cache;
  NcFunctionCache *conformal_time_cache;
  NcFunctionCache *sound_horizon_cache;
  GStaticMutex cache_lock;
  GStaticMutex update_lock;
  NcmModelCtrl *ctrl;
  gdouble z_f;
  gboolean use_cache;
};

GType nc_distance_get_type (void) G_GNUC_CONST;

NcDistance *nc_distance_new (gdouble z_f);
NcDistance *nc_distance_ref (NcDistance *dist);
void nc_distance_prepare (NcDistance *dist, NcHICosmo *model);
G_INLINE_FUNC void nc_distance_prepare_if_needed (NcDistance *dist, NcHICosmo *model);
void nc_distance_free (NcDistance *dist);

NcmMSetFunc *nc_distance_func0_new (NcDistance *dist, NcDistanceFunc0 f0);
NcmMSetFunc *nc_distance_func1_new (NcDistance *dist, NcDistanceFunc1 f1);

/***************************************************************************
 * Redshift independent 'distances'
 ****************************************************************************/

gdouble nc_distance_hubble (NcDistance *dist, NcHICosmo *model);
gdouble nc_distance_decoupling_redshift (NcDistance *dist, NcHICosmo *model);
gdouble nc_distance_drag_redshift (NcDistance *dist, NcHICosmo *model);
gdouble nc_distance_shift_parameter_lss (NcDistance *dist, NcHICosmo *model);
gdouble nc_distance_comoving_lss (NcDistance *dist, NcHICosmo *model);
gdouble nc_distance_comoving_a0_lss (NcDistance *dist, NcHICosmo *model);
gdouble nc_distance_acoustic_scale (NcDistance *dist, NcHICosmo *model);
gdouble nc_distance_dilation_scale_ratio (NcDistance *dist, NcHICosmo *model);
gdouble nc_distance_Omega_k (NcDistance *dist, NcHICosmo *model);

/***************************************************************************
 * Redshift dependent 'distances'
 ****************************************************************************/

gdouble nc_distance_comoving (NcDistance *dist, NcHICosmo *model, gdouble z);
gdouble nc_distance_comoving_a0 (NcDistance *dist, NcHICosmo *model, gdouble z);
gdouble nc_distance_curvature (NcDistance *dist, NcHICosmo *model, gdouble z);
gdouble nc_distance_transverse (NcDistance *dist, NcHICosmo *model, gdouble z);
gdouble nc_distance_luminosity (NcDistance *dist, NcHICosmo *model, gdouble z);
gdouble nc_distance_modulo (NcDistance *dist, NcHICosmo *model, gdouble z);
gdouble nc_distance_shift_parameter (NcDistance *dist, NcHICosmo *model, gdouble z);
gdouble nc_distance_dilation_scale (NcDistance *dist, NcHICosmo *model, gdouble z);
gdouble nc_distance_bao_A_scale (NcDistance *dist, NcHICosmo *model, gdouble z);
gdouble nc_distance_sound_horizon (NcDistance *dist, NcHICosmo *model, gdouble z);
gdouble nc_distance_bao_r_Dv (NcDistance *dist, NcHICosmo *model, gdouble z);

/***************************************************************************
 * Distances transformations
 ****************************************************************************/

gdouble nc_distance_comoving2luminosity (gint k, gdouble sqrt_Omega_k, gdouble hubble_dist, gdouble z, gdouble cd);
gdouble nc_distance_luminosity2modulo (gint k, gdouble sqrt_Omega_k, gdouble hubble_dist, gdouble z, gdouble ld);
gdouble nc_distance_comoving2modulo (gint k, gdouble sqrt_Omega_k, gdouble hubble_dist, gdouble z, gdouble cd);

gdouble nc_distance_luminosity2comoving (gint k, gdouble sqrt_Omega_k, gdouble hubble_dist, gdouble z, gdouble ld);
gdouble nc_distance_modulo2luminosity (gint k, gdouble sqrt_Omega_k, gdouble hubble_dist, gdouble z, gdouble dm);
gdouble nc_distance_modulo2comoving (gint k, gdouble sqrt_Omega_k, gdouble hubble_dist, gdouble z, gdouble dm);

gdouble nc_distance_dmodulo (gint k, gdouble Omega_k, gdouble sqrt_Omega_k, gdouble cd, gdouble dcd, gdouble d_omega_t);


/***************************************************************************
 *            cosmic_time.h
 *
 *  Wed Nov 12 17:06:43 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/

gdouble nc_distance_cosmic_time (NcDistance *dist, NcHICosmo *model, gdouble z);
gdouble nc_distance_lookback_time (NcDistance *dist, NcHICosmo *model, gdouble z);
gdouble nc_distance_conformal_time (NcDistance *dist, NcHICosmo *model, gdouble z);
gdouble nc_distance_conformal_lookback_time (NcDistance *dist, NcHICosmo *model, gdouble z);

gdouble nc_distance_cosmic_time_mks_scale (NcDistance *dist, NcHICosmo *model);
gdouble nc_distance_conformal_time_mks_scale (NcDistance *dist, NcHICosmo *model);

G_END_DECLS

#endif /* _NC_DISTANCE_INLINE_H_ */

#ifndef _NC_DISTANCE_INLINE_H_
#define _NC_DISTANCE_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

G_INLINE_FUNC void
nc_distance_prepare_if_needed (NcDistance *dist, NcHICosmo *model)
{
  if (ncm_model_ctrl_update (dist->ctrl, NCM_MODEL (model)))
	nc_distance_prepare (dist, model);
}

G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_DISTANCE_INLINE_H_ */
