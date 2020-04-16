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
typedef struct _NcHaloDensityProfilePrivate NcHaloDensityProfilePrivate;

struct _NcHaloDensityProfileClass
{
  /*< private >*/
  NcmModelClass parent_class;
  
  gdouble (*eval_dl_density) (NcHaloDensityProfile *dp, const gdouble x);
  gdouble (*eval_dl_spher_mass) (NcHaloDensityProfile *dp, const gdouble x);
  gdouble (*eval_dl_2d_density) (NcHaloDensityProfile *dp, const gdouble X);
  gdouble (*eval_dl_cyl_mass) (NcHaloDensityProfile *dp, const gdouble X);
};

/**
 * NcHaloDensityProfileMassDef:
 * @NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN: halo mass defined in terms of the mean density $\rho_\mathrm{bg} = \rho_m(z)$
 * @NC_HALO_DENSITY_PROFILE_MASS_DEF_CRITICAL: halo mass defined in terms of the critical density $\rho_\mathrm{bg} = \rho_\mathrm{crit}(z)$
 * @NC_HALO_DENSITY_PROFILE_MASS_DEF_VIRIAL: halo mass defined in terms of virial overdensity times the critical density $\rho_\mathrm{bg} = \rho_\mathrm{crit}(z)$
 *
 * Spherical overdensity halo mass: $$M_\Delta = \frac{4\pi}{3} \Delta \rho_\mathrm{bg} r_\Delta^3,$$
 * where $\rho_\mathrm{bg}$ is the background density of the universe at redshift z, $\rho_\mathrm{bg} (z)$.
 * For @NC_HALO_DENSITY_PROFILE_MASS_DEF_VIRIAL, the parameter #NcHaloDensityProfile:MDelta is ignored and 
 * \begin{equation}\label{def:DVir}
 * \Delta_\mathrm{Vir} = 18 \pi^2 + 82 x - 39 x^2, \quad x \equiv \Omega_m(z) - 1.
 * \end{equation}
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
 * @NC_HALO_DENSITY_PROFILE_C_DELTA: concentration parameter $r_\Delta$
 * @NC_HALO_DENSITY_PROFILE_M_DELTA: halo mass $M_\Delta$
 *
 * Fundamental parametrization of the profile $\rho(r)$,
 * any additional parameter must be included in the implementation
 * of this class.
 * 
 */
typedef enum /*< enum,underscore_name=NC_HALO_DENSITY_PROFILE_SPARAMS >*/
{
  NC_HALO_DENSITY_PROFILE_C_DELTA = 0,
  NC_HALO_DENSITY_PROFILE_M_DELTA,
  /* < private > */
  NC_HALO_DENSITY_PROFILE_SPARAM_LEN, /*< skip >*/
} NcHaloDensityProfileSParams;

struct _NcHaloDensityProfile
{
  /*< private >*/
  NcmModel parent_instance;
  NcHaloDensityProfilePrivate *priv;
};

GType nc_halo_density_profile_get_type (void) G_GNUC_CONST;

NCM_MSET_MODEL_DECLARE_ID (nc_halo_density_profile);

NcHaloDensityProfile *nc_halo_density_profile_new_from_name (gchar *density_profile_name);
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

gdouble nc_halo_density_profile_eval_dl_density (NcHaloDensityProfile *dp, const gdouble x);
gdouble nc_halo_density_profile_eval_dl_spher_mass (NcHaloDensityProfile *dp);
gdouble nc_halo_density_profile_eval_dl_2d_density (NcHaloDensityProfile *dp, const gdouble X);
gdouble nc_halo_density_profile_eval_dl_cyl_mass (NcHaloDensityProfile *dp, const gdouble X);

gdouble nc_halo_density_profile_Delta (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z);
gdouble nc_halo_density_profile_rho_bg (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z);
gdouble nc_halo_density_profile_Delta_rho_bg (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z);

gdouble nc_halo_density_profile_rho_s (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z);
gdouble nc_halo_density_profile_r_s (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z);
void nc_halo_density_profile_r_s_rho_s (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z, gdouble *r_s, gdouble *rho_s);

gdouble nc_halo_density_profile_eval_density (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble r, const gdouble z);
gdouble nc_halo_density_profile_eval_spher_mass (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z);
gdouble nc_halo_density_profile_eval_2d_density (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z);
gdouble nc_halo_density_profile_eval_cyl_mass (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z);

GArray *nc_halo_density_profile_eval_density_array (NcHaloDensityProfile *dp, NcHICosmo *cosmo, GArray *r, gdouble fin, gdouble fout, const gdouble z);
GArray *nc_halo_density_profile_eval_2d_density_array (NcHaloDensityProfile *dp, NcHICosmo *cosmo, GArray *R, gdouble fin, gdouble fout, const gdouble z);
GArray *nc_halo_density_profile_eval_cyl_mass_array (NcHaloDensityProfile *dp, NcHICosmo *cosmo, GArray *R, gdouble fin, gdouble fout, const gdouble z);

#define NC_HALO_DENSITY_PROFILE_DEFAULT_C_DELTA (4.0   )
#define NC_HALO_DENSITY_PROFILE_DEFAULT_M_DELTA (2.0e14)

#define NC_HALO_DENSITY_PROFILE_DEFAULT_PARAMS_ABSTOL (0.0)

G_END_DECLS

#endif /* _NC_HALO_DENSITY_PROFILE_H_ */
