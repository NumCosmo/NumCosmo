/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_wl_likelihood.h
 *
 *  Mon May 08 18:13:24 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_wl_likelihood.h
 * Copyright (C) 2023 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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

#ifndef _NC_GALAXY_WL_LIKELIHOOD_H_
#define _NC_GALAXY_WL_LIKELIHOOD_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_stats_dist_kde.h>
#include <numcosmo/math/ncm_stats_dist_kernel_gauss.h>
#include <numcosmo/galaxy/nc_galaxy_sd_shape.h>
#include <numcosmo/galaxy/nc_galaxy_sd_z_proxy.h>
#include <numcosmo/galaxy/nc_galaxy_sd_position.h>


G_BEGIN_DECLS

#define NC_TYPE_GALAXY_WL_LIKELIHOOD             (nc_galaxy_wl_likelihood_get_type ())
#define NC_GALAXY_WL_LIKELIHOOD(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_GALAXY_WL_LIKELIHOOD, NcGalaxyWLLikelihood))
#define NC_GALAXY_WL_LIKELIHOOD_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_GALAXY_WL_LIKELIHOOD, NcGalaxyWLLikelihoodClass))
#define NC_IS_GALAXY_WL_LIKELIHOOD(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_GALAXY_WL_LIKELIHOOD))
#define NC_IS_GALAXY_WL_LIKELIHOOD_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_GALAXY_WL_LIKELIHOOD))
#define NC_GALAXY_WL_LIKELIHOOD_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_GALAXY_WL_LIKELIHOOD, NcGalaxyWLLikelihoodClass))

typedef struct _NcGalaxyWLLikelihoodClass NcGalaxyWLLikelihoodClass;
typedef struct _NcGalaxyWLLikelihood NcGalaxyWLLikelihood;
typedef struct _NcGalaxyWLLikelihoodPrivate NcGalaxyWLLikelihoodPrivate;

struct _NcGalaxyWLLikelihoodClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcGalaxyWLLikelihood
{
  /*< private >*/
  GObject parent_instance;
  NcGalaxyWLLikelihoodPrivate *priv;
};

GType nc_galaxy_wl_likelihood_get_type (void) G_GNUC_CONST;

NcGalaxyWLLikelihood *nc_galaxy_wl_likelihood_new (NcGalaxySDShape *s_dist, NcGalaxySDZProxy *zp_dist, NcGalaxySDPosition *rz_dist);
NcGalaxyWLLikelihood *nc_galaxy_wl_likelihood_ref (NcGalaxyWLLikelihood *gwl);

void nc_galaxy_wl_likelihood_free (NcGalaxyWLLikelihood *gwl);
void nc_galaxy_wl_likelihood_clear (NcGalaxyWLLikelihood **gwl);

void nc_galaxy_wl_likelihood_set_obs (NcGalaxyWLLikelihood *gwl, NcmMatrix *obs);
NcmMatrix *nc_galaxy_wl_likelihood_peek_obs (NcGalaxyWLLikelihood *gwl);
NcmStatsDistKDE *nc_galaxy_wl_likelihood_peek_kde (NcGalaxyWLLikelihood *gwl);
void nc_galaxy_wl_likelihood_prepare (NcGalaxyWLLikelihood *gwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster);
gdouble nc_galaxy_wl_likelihood_eval_m2lnP (NcGalaxyWLLikelihood *gwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster);
gdouble nc_galaxy_wl_likelihood_kde_eval_m2lnP (NcGalaxyWLLikelihood *gwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster);
guint nc_galaxy_wl_likelihood_len (NcGalaxyWLLikelihood *gwll);
void nc_galaxy_wl_likelihood_set_cut (NcGalaxyWLLikelihood *gwl, const gdouble zp_min, const gdouble zp_max, const gdouble r_min, const gdouble r_max);
void nc_galaxy_wl_likelihood_set_ndata (NcGalaxyWLLikelihood *gwl, gdouble ndata);

G_END_DECLS

#endif /* _NC_GALAXY_WL_LIKELIHOOD_H_ */

