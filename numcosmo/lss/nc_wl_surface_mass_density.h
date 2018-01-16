/***************************************************************************
 *            nc_wl_surface_mass_density.h
 *
 *  Tue Aug 15 17:22:45 2017
 *  Copyright  2017  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2017 <pennalima@gmail.com>
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

#ifndef _NC_WL_SURFACE_MASS_DENSITY_H_
#define _NC_WL_SURFACE_MASS_DENSITY_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/math/ncm_ode_spline.h>
#include <numcosmo/math/ncm_model_ctrl.h>

G_BEGIN_DECLS

#define NC_TYPE_WL_SURFACE_MASS_DENSITY             (nc_wl_surface_mass_density_get_type ())
#define NC_WL_SURFACE_MASS_DENSITY(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_WL_SURFACE_MASS_DENSITY, NcWLSurfaceMassDensity))
#define NC_WL_SURFACE_MASS_DENSITY_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_WL_SURFACE_MASS_DENSITY, NcWLSurfaceMassDensityClass))
#define NC_IS_WL_SURFACE_MASS_DENSITY(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_WL_SURFACE_MASS_DENSITY))
#define NC_IS_WL_SURFACE_MASS_DENSITY_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_WL_SURFACE_MASS_DENSITY))
#define NC_WL_SURFACE_MASS_DENSITY_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_WL_SURFACE_MASS_DENSITY, NcWLSurfaceMassDensityClass))

typedef struct _NcWLSurfaceMassDensityClass NcWLSurfaceMassDensityClass;
typedef struct _NcWLSurfaceMassDensity NcWLSurfaceMassDensity;

struct _NcWLSurfaceMassDensityClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcWLSurfaceMassDensity
{
  /*< private >*/
  GObject parent_instance;
  NcmModelCtrl *ctrl;
  gboolean misc;
  gboolean twoh;
  gboolean cg;
  gboolean nw;
};

GType nc_wl_surface_mass_density_get_type (void) G_GNUC_CONST;

NcWLSurfaceMassDensity *nc_wl_surface_mass_density_new (gdouble zf);
NcWLSurfaceMassDensity *nc_wl_surface_mass_density_ref (NcWLSurfaceMassDensity *smd);

void nc_wl_surface_mass_density_require_zf (NcWLSurfaceMassDensity *smd, const gdouble zf);

void nc_wl_surface_mass_density_prepare (NcWLSurfaceMassDensity *smd, NcHICosmo *cosmo);

void nc_wl_surface_mass_density_free (NcWLSurfaceMassDensity *smd);
void nc_wl_surface_mass_density_clear (NcWLSurfaceMassDensity **smd);


G_END_DECLS

#endif /* _NC_WL_SURFACE_MASS_DENSITY_INLINE_H_ */

#ifndef _NC_WL_SURFACE_MASS_DENSITY_INLINE_H_
#define _NC_WL_SURFACE_MASS_DENSITY_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

G_INLINE_FUNC void
nc_wl_surface_mass_density_prepare_if_needed (NcWLSurfaceMassDensity *smd, NcHICosmo *cosmo)
{
  if (ncm_model_ctrl_update (smd->ctrl, NCM_MODEL (cosmo)))
    nc_wl_surface_mass_density_prepare (smd, cosmo);
}

G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_WL_SURFACE_MASS_DENSITY_INLINE_H_ */
