/***************************************************************************
 *            nc_halo_density_profile.h
 *
 *  Sat June 07 19:45:55 2014
 *  Copyright  2014
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_halo_density_profile.h
 * Copyright (C) 2014 Mariana Penna Lima <pennalima@gmail.com>
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

#ifndef _NC_HALO_DENSITY_PROFILE_H_
#define _NC_HALO_DENSITY_PROFILE_H_

#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/lss/nc_halo_mass_summary.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_DENSITY_PROFILE (nc_halo_density_profile_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcHaloDensityProfile, nc_halo_density_profile, NC, HALO_DENSITY_PROFILE, NcmModel)

struct _NcHaloDensityProfileClass
{
  /*< private >*/
  NcmModelClass parent_class;

  gdouble (*eval_dl_density) (NcHaloDensityProfile *dp, const gdouble x);
  gdouble (*eval_dl_spher_mass) (NcHaloDensityProfile *dp, const gdouble x);
  gdouble (*eval_dl_2d_density) (NcHaloDensityProfile *dp, const gdouble X);
  gdouble (*eval_dl_cyl_mass) (NcHaloDensityProfile *dp, const gdouble X);
};

NCM_MSET_MODEL_DECLARE_ID (nc_halo_density_profile);

NcHaloDensityProfile *nc_halo_density_profile_ref (NcHaloDensityProfile *dp);

void nc_halo_density_profile_free (NcHaloDensityProfile *dp);
void nc_halo_density_profile_clear (NcHaloDensityProfile **dp);

void nc_halo_density_profile_set_reltol (NcHaloDensityProfile *dp, const gdouble reltol);
void nc_halo_density_profile_set_lnXi (NcHaloDensityProfile *dp, const gdouble lnXi);
void nc_halo_density_profile_set_lnXf (NcHaloDensityProfile *dp, const gdouble lnXf);

gdouble nc_halo_density_profile_get_reltol (NcHaloDensityProfile *dp);
gdouble nc_halo_density_profile_get_lnXi (NcHaloDensityProfile *dp);
gdouble nc_halo_density_profile_get_lnXf (NcHaloDensityProfile *dp);

void nc_halo_density_profile_get_phys_limts (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z, gdouble *Ri, gdouble *Rf);
NcHaloMassSummary *nc_halo_density_profile_peek_mass_summary (NcHaloDensityProfile *dp);

gdouble nc_halo_density_profile_eval_dl_density (NcHaloDensityProfile *dp, const gdouble x);
gdouble nc_halo_density_profile_eval_dl_spher_mass (NcHaloDensityProfile *dp, const gdouble x);
gdouble nc_halo_density_profile_eval_dl_2d_density (NcHaloDensityProfile *dp, const gdouble X);
gdouble nc_halo_density_profile_eval_dl_cyl_mass (NcHaloDensityProfile *dp, const gdouble X);

gdouble nc_halo_density_profile_rho_s (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z);
gdouble nc_halo_density_profile_r_s (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z);
void nc_halo_density_profile_r_s_rho_s (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z, gdouble *r_s, gdouble *rho_s);

gdouble nc_halo_density_profile_eval_density (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble r, const gdouble z);
gdouble nc_halo_density_profile_eval_spher_mass (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble r, const gdouble z);
gdouble nc_halo_density_profile_eval_spher_mass_delta (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z);
gdouble nc_halo_density_profile_eval_2d_density (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z);
gdouble nc_halo_density_profile_eval_cyl_mass (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z);

GArray *nc_halo_density_profile_eval_density_array (NcHaloDensityProfile *dp, NcHICosmo *cosmo, GArray *r, gdouble fin, gdouble fout, const gdouble z);
GArray *nc_halo_density_profile_eval_2d_density_array (NcHaloDensityProfile *dp, NcHICosmo *cosmo, GArray *R, gdouble fin, gdouble fout, const gdouble z);
GArray *nc_halo_density_profile_eval_cyl_mass_array (NcHaloDensityProfile *dp, NcHICosmo *cosmo, GArray *R, gdouble fin, gdouble fout, const gdouble z);

gdouble nc_halo_density_profile_eval_numint_dl_spher_mass (NcHaloDensityProfile *dp, const gdouble X);
gdouble nc_halo_density_profile_eval_numint_dl_2d_density (NcHaloDensityProfile *dp, const gdouble X);
gdouble nc_halo_density_profile_eval_numint_dl_cyl_mass (NcHaloDensityProfile *dp, const gdouble X);

void nc_halo_density_profile_get_numint_splines (NcHaloDensityProfile *dp, NcmSpline **spher_mass, NcmSpline **twod_density, NcmSpline **cyl_mass);

G_END_DECLS

#endif /* _NC_HALO_DENSITY_PROFILE_H_ */

