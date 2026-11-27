/***************************************************************************
 *            nc_galaxy_wl_binned_scatter.h
 *
 *  Thu November 26 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_wl_binned_scatter.h
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

#ifndef _NC_GALAXY_WL_BINNED_SCATTER_H_
#define _NC_GALAXY_WL_BINNED_SCATTER_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/model/ncm_model.h>
#include <numcosmo/ncm/model/ncm_mset.h>
#include <numcosmo/ncm/core/ncm_rng.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_wl_obs_binned.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_WL_BINNED_SCATTER (nc_galaxy_wl_binned_scatter_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcGalaxyWLBinnedScatter, nc_galaxy_wl_binned_scatter, NC, GALAXY_WL_BINNED_SCATTER, NcmModel)

struct _NcGalaxyWLBinnedScatterClass
{
  /*< private >*/
  NcmModelClass parent_class;

  gdouble (*m2lnL) (NcGalaxyWLBinnedScatter *scatter, NcGalaxyWLObsBinned *binned, NcmVector *model);
  void (*gen) (NcGalaxyWLBinnedScatter *scatter, NcGalaxyWLObsBinned *binned, NcmVector *model, NcmRNG *rng, NcmVector *data);

  /* Padding to allow adding up to 8 more virtual functions without breaking ABI. */
  gpointer padding[8];
};

NCM_MSET_MODEL_DECLARE_ID (nc_galaxy_wl_binned_scatter);

NcGalaxyWLBinnedScatter *nc_galaxy_wl_binned_scatter_ref (NcGalaxyWLBinnedScatter *scatter);
void nc_galaxy_wl_binned_scatter_free (NcGalaxyWLBinnedScatter *scatter);
void nc_galaxy_wl_binned_scatter_clear (NcGalaxyWLBinnedScatter **scatter);

gdouble nc_galaxy_wl_binned_scatter_m2lnL (NcGalaxyWLBinnedScatter *scatter, NcGalaxyWLObsBinned *binned, NcmVector *model);
void nc_galaxy_wl_binned_scatter_gen (NcGalaxyWLBinnedScatter *scatter, NcGalaxyWLObsBinned *binned, NcmVector *model, NcmRNG *rng, NcmVector *data);

G_END_DECLS

#endif /* _NC_GALAXY_WL_BINNED_SCATTER_H_ */

