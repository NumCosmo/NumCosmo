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

G_BEGIN_DECLS

#define NC_TYPE_HALO_DENSITY_PROFILE             (nc_halo_density_profile_get_type ())
#define NC_HALO_DENSITY_PROFILE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HALO_DENSITY_PROFILE, NcHaloDensityProfile))
#define NC_HALO_DENSITY_PROFILE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HALO_DENSITY_PROFILE, NcHaloDensityProfileClass))
#define NC_IS_HALO_DENSITY_PROFILE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HALO_DENSITY_PROFILE))
#define NC_IS_HALO_DENSITY_PROFILE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HALO_DENSITY_PROFILE))
#define NC_HALO_DENSITY_PROFILE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HALO_DENSITY_PROFILE, NcHaloDensityProfileClass))

typedef struct _NcHaloDensityProfileClass NcHaloDensityProfileClass;
typedef struct _NcHaloDensityProfile NcHaloDensityProfile;

struct _NcHaloDensityProfileClass
{
  /*< private >*/
  NcmModelClass parent_class;
  gdouble (*eval_density) (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble r, const gdouble z);
  gdouble (*integral_density_los) (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z);
  gdouble (*integral_density_2d) (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z);
  gdouble (*eval_fourier) (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble k, const gdouble M, const gdouble z);
	gdouble (*scale_radius) (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z);
	gdouble (*central_density) (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z);
};

/**
 * NcHaloDensityProfileMassDef:
 * @NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN: halo mass defined in terms of the mean density
 * @NC_HALO_DENSITY_PROFILE_MASS_DEF_CRITICAL: halo mass defined in terms of the critical density
 * @NC_HALO_DENSITY_PROFILE_MASS_DEF_VIRIAL: halo mass defined in terms of virial overdensity times the critical density 
 * 
 * Spherical overdensity halo mass: $$M = \frac{4\pi}{3} \Delta \rho R^3,$$ 
 * where $\rho$ is the mean density of the universe at redshift z, $\rho_m (z)$, or the critical density at z, 
 * $\rho_c (z)$, or $\Delta_{\text{vir}}$ times $\rho_c (z)$.
 * 
 */
typedef enum _NcHaloDensityProfileMassDef
{
	NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN = 0,
	NC_HALO_DENSITY_PROFILE_MASS_DEF_CRITICAL,
	NC_HALO_DENSITY_PROFILE_MASS_DEF_VIRIAL,
	/* < private > */
  NC_HALO_DENSITY_PROFILE_MASS_DEF_LEN, /*< skip >*/	
} NcHaloDensityProfileMassDef; 

/**
 * NcHaloDensityProfileSParams:
 * @NC_HALO_DENSITY_PROFILE_C: concentration parameter
 * @NC_HALO_DENSITY_PROFILE_M_DELTA: halo mass
 *
 * FIXME
 */
typedef enum /*< enum,underscore_name=NC_HALO_DENSITY_PROFILE_SPARAMS >*/
{
  NC_HALO_DENSITY_PROFILE_C = 0,
  NC_HALO_DENSITY_PROFILE_M_DELTA, 
  /* < private > */
  NC_HALO_DENSITY_PROFILE_SPARAM_LEN, /*< skip >*/
} NcHaloDensityProfileSParams;

struct _NcHaloDensityProfile
{
  /*< private >*/
  NcmModel parent_instance;
	NcHaloDensityProfileMassDef mdef; /* mass definition*/
	gdouble z;                    /* redshift */
	gdouble oDelta;               /* overdensity constant */
};

GType nc_halo_density_profile_get_type (void) G_GNUC_CONST;

NCM_MSET_MODEL_DECLARE_ID (nc_halo_density_profile);

NcHaloDensityProfile *nc_halo_density_profile_new_from_name (gchar *density_profile_name);
NcHaloDensityProfile *nc_halo_density_profile_ref (NcHaloDensityProfile *dp);
void nc_halo_density_profile_free (NcHaloDensityProfile *dp);
void nc_halo_density_profile_clear (NcHaloDensityProfile **dp);

gdouble nc_halo_density_profile_eval_density (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble r, const gdouble z);
gdouble nc_halo_density_profile_integral_density_los (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z);
gdouble nc_halo_density_profile_integral_density_2d (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z);
gdouble nc_halo_density_profile_eval_fourier (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble k, const gdouble M, const gdouble z);
gdouble nc_halo_density_profile_scale_radius (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z);
gdouble nc_halo_density_profile_Delta (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z);
gdouble nc_halo_density_profile_mass_density (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z);
gdouble nc_halo_density_profile_mass_density_threshold (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z);
gdouble nc_halo_density_profile_central_density (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z); 

#define NC_HALO_DENSITY_PROFILE_DEFAULT_C        (4.0)
#define NC_HALO_DENSITY_PROFILE_DEFAULT_M_DELTA  (2.0e14)
#define NC_HALO_DENSITY_PROFILE_DEFAULT_PARAMS_ABSTOL (0.0)

G_END_DECLS

#endif /* _NC_HALO_DENSITY_PROFILE_H_ */

