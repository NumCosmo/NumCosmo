/***************************************************************************
 *            nc_density_profile.h
 *
 *  Sat June 07 19:45:55 2014
 *  Copyright  2014  
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_density_profile.h
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

#ifndef _NC_DENSITY_PROFILE_H_
#define _NC_DENSITY_PROFILE_H_

#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/nc_hicosmo.h>

G_BEGIN_DECLS

#define NC_TYPE_DENSITY_PROFILE             (nc_density_profile_get_type ())
#define NC_DENSITY_PROFILE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DENSITY_PROFILE, NcDensityProfile))
#define NC_DENSITY_PROFILE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DENSITY_PROFILE, NcDensityProfileClass))
#define NC_IS_DENSITY_PROFILE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DENSITY_PROFILE))
#define NC_IS_DENSITY_PROFILE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DENSITY_PROFILE))
#define NC_DENSITY_PROFILE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DENSITY_PROFILE, NcDensityProfileClass))

typedef struct _NcDensityProfileClass NcDensityProfileClass;
typedef struct _NcDensityProfile NcDensityProfile;

struct _NcDensityProfileClass
{
  /*< private >*/
  NcmModelClass parent_class;
  gdouble (*eval_density) (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble r, const gdouble z);
  gdouble (*integral_density_los) (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z);
  gdouble (*integral_density_2d) (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z);
  gdouble (*eval_fourier) (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble k, const gdouble M, const gdouble z);
	gdouble (*scale_radius) (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble z);
	gdouble (*central_density) (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z);
};

/**
 * NcDensityProfileMassDef:
 * @NC_DENSITY_PROFILE_MEAN: halo mass defined in terms of the mean density
 * @NC_DENSITY_PROFILE_CRITICAL: halo mass defined in terms of the critical density
 * @NC_DENSITY_PROFILE_VIRIAL: halo mass defined in terms of virial overdensity times the critical density 
 * 
 * Spherical overdensity halo mass: $$M = \frac{4\pi}{3} \Delta \rho R^3,$$ 
 * where $\rho$ is the mean density of the universe at redshift z, $\rho_m (z)$, or the critical density at z, 
 * $\rho_c (z)$, or $\Delta_{\text{vir}}$ times $\rho_c (z)$.
 * 
 */
typedef enum _NcDensityProfileMassDef
{
	NC_DENSITY_PROFILE_MASS_DEF_MEAN = 0,
	NC_DENSITY_PROFILE_MASS_DEF_CRITICAL,
	NC_DENSITY_PROFILE_MASS_DEF_VIRIAL,
	/* < private > */
  NC_DENSITY_PROFILE_MASS_DEF_LEN, /*< skip >*/	
} NcDensityProfileMassDef; 

/**
 * NcDensityProfileParams:
 * @NC_DENSITY_PROFILE_C: concentration parameter
 * @NC_DENSITY_PROFILE_M_DELTA: halo mass
 *
 * FIXME
 */
typedef enum /*< enum,underscore_name=NC_DENSITY_PROFILE_SPARAMS >*/
{
  NC_DENSITY_PROFILE_C = 0,
  NC_DENSITY_PROFILE_M_DELTA, 
  /* < private > */
  NC_DENSITY_PROFILE_SPARAM_LEN, /*< skip >*/
} NcDensityProfileSParams;

struct _NcDensityProfile
{
  /*< private >*/
  NcmModel parent_instance;
	NcDensityProfileMassDef mdef; /* mass definition*/
	gdouble z;                    /* redshift */
	gdouble oDelta;               /* overdensity constant */
};

GType nc_density_profile_get_type (void) G_GNUC_CONST;

NCM_MSET_MODEL_DECLARE_ID (nc_density_profile);

NcDensityProfile *nc_density_profile_new_from_name (gchar *density_profile_name);
NcDensityProfile *nc_density_profile_ref (NcDensityProfile *dp);
void nc_density_profile_free (NcDensityProfile *dp);
void nc_density_profile_clear (NcDensityProfile **dp);

gdouble nc_density_profile_eval_density (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble r, const gdouble z);
gdouble nc_density_profile_integral_density_los (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z);
gdouble nc_density_profile_integral_density_2d (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z);
gdouble nc_density_profile_eval_fourier (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble k, const gdouble M, const gdouble z);
gdouble nc_density_profile_scale_radius (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble z);
gdouble nc_density_profile_Delta (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble z);
gdouble nc_density_profile_mass_density (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble z);
gdouble nc_density_profile_mass_density_threshold (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble z);
gdouble nc_density_profile_central_density (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z); 

#define NC_DENSITY_PROFILE_DEFAULT_C        (4.0)
#define NC_DENSITY_PROFILE_DEFAULT_M_DELTA  (2.0e14)
#define NC_DENSITY_PROFILE_DEFAULT_PARAMS_ABSTOL (0.0)

G_END_DECLS

#endif /* _NC_DENSITY_PROFILE_H_ */

