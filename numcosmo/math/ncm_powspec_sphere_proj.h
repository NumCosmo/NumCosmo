/***************************************************************************
 *            ncm_powspec_sphere_proj.h
 *
 *  Mon April 01 11:07:10 2019
 *  Copyright  2019  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_powspec_sphere_proj.h
 * Copyright (C) 2019 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_POWSPEC_SPHERE_PROJ_H_
#define _NCM_POWSPEC_SPHERE_PROJ_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_powspec.h>
#include <numcosmo/math/ncm_fftlog.h>
#include <numcosmo/math/ncm_spline2d.h>

G_BEGIN_DECLS

#define NCM_TYPE_POWSPEC_SPHERE_PROJ             (ncm_powspec_sphere_proj_get_type ())
#define NCM_POWSPEC_SPHERE_PROJ(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_POWSPEC_SPHERE_PROJ, NcmPowspecSphereProj))
#define NCM_POWSPEC_SPHERE_PROJ_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_POWSPEC_SPHERE_PROJ, NcmPowspecSphereProjClass))
#define NCM_IS_POWSPEC_SPHERE_PROJ(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_POWSPEC_SPHERE_PROJ))
#define NCM_IS_POWSPEC_SPHERE_PROJ_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_POWSPEC_SPHERE_PROJ))
#define NCM_POWSPEC_SPHERE_PROJ_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_POWSPEC_SPHERE_PROJ, NcmPowspecSphereProjClass))

typedef struct _NcmPowspecSphereProjClass NcmPowspecSphereProjClass;
typedef struct _NcmPowspecSphereProj NcmPowspecSphereProj;
typedef struct _NcmPowspecSphereProjPrivate NcmPowspecSphereProjPrivate;

struct _NcmPowspecSphereProjClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcmPowspecSphereProj
{
  /*< private >*/
  GObject parent_instance;
  NcmPowspecSphereProjPrivate *priv;
};

GType ncm_powspec_sphere_proj_get_type (void) G_GNUC_CONST;

NcmPowspecSphereProj *ncm_powspec_sphere_proj_new (NcmPowspec *ps, guint ell_min, guint ell_max);
NcmPowspecSphereProj *ncm_powspec_sphere_proj_ref (NcmPowspecSphereProj *psp);
void ncm_powspec_sphere_proj_free (NcmPowspecSphereProj *psp);
void ncm_powspec_sphere_proj_clear (NcmPowspecSphereProj **psp);

void ncm_powspec_sphere_proj_prepare (NcmPowspecSphereProj *psp, NcmModel *model);
void ncm_powspec_sphere_proj_prepare_if_needed (NcmPowspecSphereProj *psp, NcmModel *model);

void ncm_powspec_sphere_proj_set_xi_i (NcmPowspecSphereProj *psp, gdouble xi_i);
void ncm_powspec_sphere_proj_set_xi_f (NcmPowspecSphereProj *psp, gdouble xi_f);
void ncm_powspec_sphere_proj_set_k_pivot (NcmPowspecSphereProj *psp, gdouble k_pivot);

gdouble ncm_powspec_sphere_proj_get_r_min (NcmPowspecSphereProj *psp);
gdouble ncm_powspec_sphere_proj_get_r_max (NcmPowspecSphereProj *psp);

void ncm_powspec_sphere_proj_set_ell_min (NcmPowspecSphereProj *psp, guint ell_min);
void ncm_powspec_sphere_proj_set_ell_max (NcmPowspecSphereProj *psp, guint ell_max);
guint ncm_powspec_sphere_proj_get_ell_min (NcmPowspecSphereProj *psp);
guint ncm_powspec_sphere_proj_get_ell_max (NcmPowspecSphereProj *psp);

gdouble ncm_powspec_sphere_proj_eval_lnvar_lnr (NcmPowspecSphereProj *psp, const gdouble z, const gdouble lnr);

gdouble ncm_powspec_sphere_proj_get_w (NcmPowspecSphereProj *psp, const guint w_i);
void ncm_powspec_sphere_proj_get_ell (NcmPowspecSphereProj *psp, const guint w_i, const gint ell, GArray **lnxi, GArray **Cell);
NcmSpline *ncm_powspec_sphere_proj_peek_ell_spline (NcmPowspecSphereProj *psp, const guint w_i, const gint ell);
gdouble ncm_powspec_sphere_proj_eval_Cell_xi1_xi2 (NcmPowspecSphereProj *psp, NcmModel *model, const gint ell, const gdouble z1, const gdouble z2, const gdouble xi1, const gdouble xi2);

#define NCM_POWSPEC_SPHERE_PROJ_DEFAULT_SIZE (200)

G_END_DECLS

#endif /* _NCM_POWSPEC_SPHERE_PROJ_H_ */
