/***************************************************************************
 *            nc_wl_cluster_lens.h
 *
 *  Fri November 27 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_wl_cluster_lens.h
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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

#ifndef _NC_WL_CLUSTER_LENS_H_
#define _NC_WL_CLUSTER_LENS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/background/nc_hicosmo.h>
#include <numcosmo/nc/background/nc_distance.h>
#include <numcosmo/nc/lss/halo/nc_halo_position.h>
#include <numcosmo/nc/lss/halo/nc_halo_density_profile.h>
#include <numcosmo/ncm/algebra/ncm_spectral_panels.h>

G_BEGIN_DECLS

#define NC_TYPE_WL_CLUSTER_LENS (nc_wl_cluster_lens_get_type ())

G_DECLARE_FINAL_TYPE (NcWLClusterLens, nc_wl_cluster_lens, NC, WL_CLUSTER_LENS, GObject)

NcWLClusterLens *nc_wl_cluster_lens_new (NcDistance * dist);
NcWLClusterLens *nc_wl_cluster_lens_ref (NcWLClusterLens *lens);
void nc_wl_cluster_lens_free (NcWLClusterLens *lens);
void nc_wl_cluster_lens_clear (NcWLClusterLens **lens);

NcDistance *nc_wl_cluster_lens_peek_distance (NcWLClusterLens *lens);

void nc_wl_cluster_lens_prepare (NcWLClusterLens *lens, NcHICosmo *cosmo, NcHaloPosition *hp, NcHaloDensityProfile *dp);
void nc_wl_cluster_lens_prepare_sigma_crit_inv (NcWLClusterLens *lens, NcHICosmo *cosmo, const gdouble z_max);
gdouble nc_wl_cluster_lens_get_z_cl (NcWLClusterLens *lens);
guint64 nc_wl_cluster_lens_get_geometry_hash (NcWLClusterLens *lens);
guint64 nc_wl_cluster_lens_get_profile_hash (NcWLClusterLens *lens);

void nc_wl_cluster_lens_set_sigma_crit_inv_reltol (NcWLClusterLens *lens, const gdouble reltol);
gdouble nc_wl_cluster_lens_get_sigma_crit_inv_reltol (NcWLClusterLens *lens);
NcmSpectralPanels *nc_wl_cluster_lens_peek_sigma_crit_inv (NcWLClusterLens *lens);

gdouble nc_wl_cluster_lens_sigma_crit (NcWLClusterLens *lens, NcHICosmo *cosmo, const gdouble z_s);
gdouble nc_wl_cluster_lens_sigma_crit_inv (NcWLClusterLens *lens, NcHICosmo *cosmo, const gdouble z_s);
void nc_wl_cluster_lens_sigma (NcWLClusterLens *lens, NcHaloDensityProfile *dp, const gdouble R, gdouble *sigma, gdouble *mean_sigma);

NCM_INLINE gdouble nc_wl_cluster_lens_reduced_shear (const gdouble sigma, const gdouble mean_sigma, const gdouble sigma_crit_inv);

gdouble nc_wl_cluster_lens_reduced_shear_at (NcWLClusterLens *lens, NcHICosmo *cosmo, const gdouble sigma, const gdouble mean_sigma, const gdouble z_s);

G_END_DECLS

#endif /* _NC_WL_CLUSTER_LENS_H_ */

#ifndef _NC_WL_CLUSTER_LENS_INLINE_H_
#define _NC_WL_CLUSTER_LENS_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE gdouble
nc_wl_cluster_lens_reduced_shear (const gdouble sigma, const gdouble mean_sigma, const gdouble sigma_crit_inv)
{
  return (mean_sigma - sigma) * sigma_crit_inv / (1.0 - sigma * sigma_crit_inv);
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_WL_CLUSTER_LENS_INLINE_H_ */

