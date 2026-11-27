/***************************************************************************
 *            nc_galaxy_wl_binned_estimator_private.h
 *
 *  Thu November 26 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_wl_binned_estimator_private.h
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

#ifndef _NC_GALAXY_WL_BINNED_ESTIMATOR_PRIVATE_H_
#define _NC_GALAXY_WL_BINNED_ESTIMATOR_PRIVATE_H_

#include <glib.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_wl_obs.h>

G_BEGIN_DECLS

/*
 * The galaxies of a per-galaxy catalog that fall in one of the bins, with the quantities
 * every classic estimator needs. Built by nc_galaxy_wl_binned_estimator_build() and passed
 * to the build_bins() vfunc, which reads it only.
 *
 * - len: the number of galaxies (all arrays have this length).
 * - n_bins: the number of non-empty bins; bin[i] < n_bins.
 * - row: the row of the galaxy in the catalog.
 * - bin: the index of its bin, counting only the non-empty bins.
 * - radius: projected radius at the cluster redshift (proper Mpc).
 * - e_t, e_x: tangential and cross components of the observed ellipticity, in the
 *   convention of the catalog (ellip_conv).
 * - c_t, c_x: the same rotation of the additive bias (c1, c2).
 * - c1, c2: the additive bias in the frame of the catalog.
 * - rot_re, rot_im: cos (2 phi), sin (2 phi), the rotation to the tangential frame, so that
 *   e_t + i e_x = (e1 + i e2) (rot_re - i rot_im).
 * - m: multiplicative bias.
 * - std_noise: per-component shape measurement noise.
 * - e_rms: per-component intrinsic ellipticity rms.
 * - inv_sigma_crit: <Sigma_cr^-1> over the redshift density, zero below the cluster.
 * - inv_sigma_crit2: <Sigma_cr^-2> over the redshift density.
 * - sigma_crit_inf: Sigma_cr for a source at infinity.
 * - shear_scale: the factor between the ellipticity and the shear, before responsivity:
 *   2 for the distortion (TRACE) and 1 for the ellipticity (TRACE_DET).
 * - use_responsivity: whether the convention has a responsivity (TRACE).
 */
struct _NcGalaxyWLBinnedEstimatorSources
{
  guint len;
  guint n_bins;
  GArray *row;
  GArray *bin;
  GArray *radius;
  GArray *e_t;
  GArray *e_x;
  GArray *c_t;
  GArray *c_x;
  GArray *c1;
  GArray *c2;
  GArray *rot_re;
  GArray *rot_im;
  GArray *m;
  GArray *std_noise;
  GArray *e_rms;
  GArray *inv_sigma_crit;
  GArray *inv_sigma_crit2;
  gdouble sigma_crit_inf;
  gdouble shear_scale;
  gboolean use_responsivity;
};

#define NC_WLBE_SRC(src, name, i) (g_array_index ((src)->name, gdouble, (i)))
#define NC_WLBE_SRC_BIN(src, i) (g_array_index ((src)->bin, guint, (i)))

G_END_DECLS

#endif /* _NC_GALAXY_WL_BINNED_ESTIMATOR_PRIVATE_H_ */

