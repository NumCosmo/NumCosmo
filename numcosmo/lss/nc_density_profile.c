/***************************************************************************
 *            nc_density_profile.c
 *
 *  Sat June 07 19:46:31 2014
 *  Copyright  2014  
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_density_profile.c
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

/**
 * SECTION:nc_density_profile
 * @title: NcDensityProfile
 * @short_description: Abstract class for density profile functions. 
 * 
 * This module comprises the set of functions to compute the matter density profile in both real 
 * and Fourier spaces. 
 * 
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_enum_types.h"
#include "lss/nc_density_profile.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"

G_DEFINE_TYPE (NcDensityProfile, nc_density_profile, NCM_TYPE_MODEL);

#define VECTOR (NCM_MODEL (dp)->params)
#define M_DELTA (ncm_vector_get (VECTOR, NC_DENSITY_PROFILE_M_DELTA))
#define C_DELTA   (ncm_vector_get (VECTOR, NC_DENSITY_PROFILE_C))

enum
{
  PROP_0,
	PROP_MDEF,
  PROP_DELTA,
  PROP_SIZE,
};

static void
nc_density_profile_init (NcDensityProfile *dp)
{
  dp->mdef   = NC_DENSITY_PROFILE_MASS_DEF_LEN;
	dp->z      = 0.0;
	dp->oDelta = 0.0;
}

static void
_nc_density_profile_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcDensityProfile *dp = NC_DENSITY_PROFILE (object);
  g_return_if_fail (NC_IS_DENSITY_PROFILE (object));

  switch (prop_id)
  {
		case PROP_MDEF:
			dp->mdef = g_value_get_enum (value); 
			break;
    case PROP_DELTA:
      dp->oDelta = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_density_profile_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDensityProfile *dp = NC_DENSITY_PROFILE (object);
  g_return_if_fail (NC_IS_DENSITY_PROFILE (object));

  switch (prop_id)
  {
		case PROP_MDEF:
      g_value_set_enum (value, dp->mdef);
      break;
    case PROP_DELTA:
      g_value_set_double (value, dp->oDelta);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_density_profile_dispose (GObject *object)
{
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_density_profile_parent_class)->dispose (object);
}

static void
_nc_density_profile_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_density_profile_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_density_profile, NC_TYPE_DENSITY_PROFILE);

static gdouble _nc_density_profile_eval_density (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble r, const gdouble z) { g_error ("Method eval_density not implemented by `%s' class.", G_OBJECT_CLASS_NAME (dp)); return 0.0; }
static gdouble _nc_density_profile_integral_density_los (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z) { g_error ("Method integral_density_los not implemented by `%s' class.", G_OBJECT_CLASS_NAME (dp)); return 0.0; }
static gdouble _nc_density_profile_integral_density_2d (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z) { g_error ("Method integral_density_2d not implemented by `%s' class.", G_OBJECT_CLASS_NAME (dp)); return 0.0; }
static gdouble _nc_density_profile_eval_fourier (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble k, const gdouble M, const gdouble z) { g_error ("Method eval_fourier not implemented by `%s' class.", G_OBJECT_CLASS_NAME (dp)); return 0.0; }
static gdouble _nc_density_profile_scale_radius (NcDensityProfile *dp, NcHICosmo *cosmo, gdouble z);
static gdouble _nc_density_profile_central_density (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z);

static void
nc_density_profile_class_init (NcDensityProfileClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmModelClass* model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_density_profile_set_property;
  model_class->get_property = &_nc_density_profile_get_property;
  object_class->dispose     = &_nc_density_profile_dispose;
  object_class->finalize    = &_nc_density_profile_finalize;

  ncm_mset_model_register_id (model_class,
                              "NcDensityProfile",
                              "NcDensityProfile.",
                              NULL,
                              FALSE,
                              NCM_MSET_MODEL_MAIN);

  ncm_model_class_set_name_nick (model_class, "Matter Density Profile", "DensityProfile");
  ncm_model_class_add_params (model_class, NC_DENSITY_PROFILE_SPARAM_LEN, 0, PROP_SIZE);

  /**
   * NcDensityProfile:c:
   * 
   * Concentration parameter.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_DENSITY_PROFILE_C, "c_{\\Delta}", "c",
                              0.5,  10.0, 1.0e-1,
                              NC_DENSITY_PROFILE_DEFAULT_PARAMS_ABSTOL, NC_DENSITY_PROFILE_DEFAULT_C,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcDensityProfile:MDelta:
   * 
   * Cluster mass within $R_\Delta$, where $\Delta$ is the overdensity.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_DENSITY_PROFILE_M_DELTA, "M_{\\Delta}", "MDelta",
                              1.0e13,  1.0e16, 1.0e13,
                              NC_DENSITY_PROFILE_DEFAULT_PARAMS_ABSTOL, NC_DENSITY_PROFILE_DEFAULT_M_DELTA,
                              NCM_PARAM_TYPE_FIXED);
	/**
   * NcDensityProfile:mdef:
   *
   * Mass definition.
   */
	g_object_class_install_property (object_class,
	                                 PROP_MDEF,
	                                 g_param_spec_enum ("mass-def",
	                                                    NULL,
	                                                    "Mass definition",
	                                                    NC_TYPE_DENSITY_PROFILE_MASS_DEF, NC_DENSITY_PROFILE_MASS_DEF_MEAN,
	                                                    G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));	
  /**
   * NcDensityProfile:Delta:
   *
   * Constant that indicates the overdensity with respect to the critical density. 
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_DELTA,
                                   g_param_spec_double ("Delta",
                                                        NULL,
                                                        "Overdensity constant",
                                                        200, 1500, 200,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	
	klass->eval_density         = &_nc_density_profile_eval_density; 
  klass->integral_density_los = &_nc_density_profile_integral_density_los; 
	klass->integral_density_2d  = &_nc_density_profile_integral_density_2d;
  klass->eval_fourier         = &_nc_density_profile_eval_fourier;
	klass->scale_radius         = &_nc_density_profile_scale_radius;
	klass->central_density      = &_nc_density_profile_central_density;
}

/**
 * nc_density_profile_new_from_name:
 * @density_profile_name: Density profile's name
 * 
 * This function returns a new #NcDensityProfile whose type is defined by @density_profile_name string.
 * 
 * Returns: A new #NcDensityProfile.
 */
NcDensityProfile *
nc_density_profile_new_from_name (gchar *density_profile_name)
{
  GObject *obj = ncm_serialize_global_from_string (density_profile_name);
  GType density_profile_type = G_OBJECT_TYPE (obj);
  if (!g_type_is_a (density_profile_type, NC_TYPE_DENSITY_PROFILE))
    g_error ("nc_density_profile_new_from_name: NcDensityProfile %s do not descend from %s.", density_profile_name, 
             g_type_name (NC_TYPE_DENSITY_PROFILE));

	return NC_DENSITY_PROFILE (obj);
}

/**
 * nc_density_profile_ref:
 * @dp: a #NcDensityProfile
 *
 * Increases the reference count of @dp by one.
 *
 * Returns: (transfer full): @dp
 */
NcDensityProfile *
nc_density_profile_ref (NcDensityProfile *dp)
{
  return g_object_ref (dp);
}

/**
 * nc_density_profile_free: 
 * @dp: a #NcDensityProfile.
 *
 * Atomically decrements the reference count of @dp by one. If the reference count drops to 0, 
 * all memory allocated by @dp is released.
 *
 */
void 
nc_density_profile_free (NcDensityProfile *dp)
{
  g_object_unref (dp);
}

/**
 * nc_density_profile_clear: 
 * @dp: a #NcDensityProfile.
 *
 * Atomically decrements the reference count of @dp by one. If the reference count drops to 0, 
 * all memory allocated by @dp is released. Set the pointer to NULL;
 *
 */
void 
nc_density_profile_clear (NcDensityProfile **dp)
{
  g_clear_object (dp);
}

/**
 * nc_density_profile_eval_density:
 * @dp: a #NcDensityProfile
 * @cosmo: a #NcHICosmo
 * @r: radius 
 * @z: redshift
 *  
 * This function computes the density profile in real space. 
 *
 * Returns: The value of the density profile in real space.
 */
gdouble 
nc_density_profile_eval_density (NcDensityProfile *dp, NcHICosmo *cosmo, gdouble r, gdouble z) 
{ 
  return NC_DENSITY_PROFILE_GET_CLASS (dp)->eval_density (dp, cosmo, r, z);
}

/**
 * nc_density_profile_integral_density_los:
 * @dp: a #NcDensityProfile
 * @cosmo: a #NcHICosmo
 * @R: radius 
 * @z: redshift
 *  
 * This function computes the integral of density profile along the line of sight. 
 *
 * Returns: The value of the integral of the density profile along the line of sight.
 */
gdouble 
nc_density_profile_integral_density_los (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z) 
{ 
  return NC_DENSITY_PROFILE_GET_CLASS (dp)->integral_density_los (dp, cosmo, R, z);
}

/**
 * nc_density_profile_integral_density_2d:
 * @dp: a #NcDensityProfile
 * @cosmo: a #NcHICosmo
 * @R: radius 
 * @z: redshift
 *  
 * This function computes the 2D integral of density profile along the line of sight and the cluster radius ($r < R$). 
 *
 * Returns: The value of the 2D integral of the density profile along the line of sight and the cluster radius ($r < R$).
 */
gdouble 
nc_density_profile_integral_density_2d (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z) 
{ 
  return NC_DENSITY_PROFILE_GET_CLASS (dp)->integral_density_2d (dp, cosmo, R, z);
}

/**
 * nc_density_profile_eval_fourier:
 * @dp: a #NcDensityProfile
 * @cosmo: a #NcHICosmo
 * @k: mode 
 * @M: mass
 * @z: redshift
 *  
 * This function computes the density profile in the Fourier space. 
 *
 * Returns: The value of the density profile in the Fourier space.
 */
gdouble 
nc_density_profile_eval_fourier (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble k, const gdouble M, const gdouble z) 
{ 
  return NC_DENSITY_PROFILE_GET_CLASS (dp)->eval_fourier (dp, cosmo, k, M, z);
}

static gdouble 
_nc_density_profile_scale_radius (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble z)
{ 
	const gdouble M   = M_DELTA;
	const gdouble c   = C_DELTA;
	const gdouble rho = nc_density_profile_mass_density_threshold (dp, cosmo, z);
	const gdouble R3  = (3.0 * M) / (4.0 * M_PI * rho);
	const gdouble R   = cbrt (R3);

	return R / c;                          
}

/**
 * nc_density_profile_scale_radius:
 * @dp: a #NcDensityProfile
 * @cosmo: a #NcHICosmo
 * @z: redshift
 *  
 * This function computes the scale radius $r_s$ in Mpc. 
 *
 * Returns: The value of the scale radius $r_s$.
 */
gdouble
nc_density_profile_scale_radius (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble z) 
{ 
	return NC_DENSITY_PROFILE_GET_CLASS (dp)->scale_radius (dp, cosmo, z);
}

/**
 * nc_density_profile_Delta:
 * @dp: a #NcDensityProfile
 * @cosmo: a #NcHICosmo
 * @z: redshift
 *  
 * This function computes the overdensity with respect to the mass density $\Delta$. 
 *
 * The virial overdensity in units of the critical density.
 * Following Colossus code (Diemer 2018) INCLUIR REF!
 * This function uses the fitting formula of Bryan & Norman 1998 INCLUIR REF!
 *
 * Returns: $\Delta$.
 */
gdouble 
nc_density_profile_Delta (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble z)
{
	switch (dp->mdef)
	{
		case NC_DENSITY_PROFILE_MASS_DEF_MEAN:
		case NC_DENSITY_PROFILE_MASS_DEF_CRITICAL:
			return dp->oDelta;
			break;
		case NC_DENSITY_PROFILE_MASS_DEF_VIRIAL:
		{
  		const gdouble x = nc_hicosmo_E2Omega_m (cosmo, z) / nc_hicosmo_E2 (cosmo, z) - 1.0;
	  	return 18.0 * M_PI * M_PI + 82.0 * x - 39.0 * x * x;
			break;
	  }
  	default:
			g_assert_not_reached ();
			return 0.0;
			break;
	}
}

/**
 * nc_density_profile_mass_density:
 * @dp: a #NcDensityProfile
 * @cosmo: a #NcHICosmo
 * @z: redshift
 *  
 * This function computes the mass density in $M_\odot\mathrm{Mpc}^{-3}$. 
 *
 * Returns: The value of the mass density in $M_\odot\mathrm{Mpc}^{-3}$.
 */
gdouble 
nc_density_profile_mass_density (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble z)
{
	switch (dp->mdef)
	{
		case NC_DENSITY_PROFILE_MASS_DEF_MEAN:
			return ncm_c_crit_mass_density_h2_solar_mass_Mpc3 () * nc_hicosmo_h2 (cosmo) * nc_hicosmo_E2Omega_m (cosmo, z);
			break;
		case NC_DENSITY_PROFILE_MASS_DEF_CRITICAL:
		case NC_DENSITY_PROFILE_MASS_DEF_VIRIAL:
			return ncm_c_crit_mass_density_h2_solar_mass_Mpc3 () * nc_hicosmo_h2 (cosmo) * nc_hicosmo_E2 (cosmo, z);
			break;
  	default:
			g_assert_not_reached ();
			return 0.0;
			break;
	}
}

/**
 * nc_density_profile_mass_density_threshold:
 * @dp: a #NcDensityProfile
 * @cosmo: a #NcHICosmo
 * @z: redshift
 *  
 * This function computes the mass density threshold in $M_\odot\mathrm{Mpc}^{-3}$. 
 *
 * Returns: The value of the mass density threshold in $M_\odot\mathrm{Mpc}^{-3}$.
 */
gdouble 
nc_density_profile_mass_density_threshold (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble z)
{
	const gdouble rhoM  = nc_density_profile_mass_density (dp, cosmo, z);
	const gdouble Delta = nc_density_profile_Delta (dp, cosmo, z);

	return Delta * rhoM;
}

static gdouble _nc_density_profile_central_density (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble r, const gdouble z)
{
	const gdouble M   = M_DELTA;
	const gdouble rs   = nc_density_profile_scale_radius (dp, cosmo, z);
	const gdouble rs3  = rs * rs * rs;
	const gdouble x    = r/rs;
	const gdouble mu   = log (1.0 + x) - x / (1.0 + x);  
	const gdouble rhos = M / (4.0 * M_PI * rs3 * mu); /* as defined in Colossus */
	return rhos;
}

/**
 * nc_density_profile_central_density:
 * @dp: a #NcDensityProfile
 * @cosmo: a #NcHICosmo
 * @R: FIXME  
 * @z: redshift
 *  
 * This function computes the central density $\rho_s$ in $M_\odot / Mpc^3$ . 
 *
 * Returns: The value of the central density $\rho_s$.
 */
gdouble
nc_density_profile_central_density (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z) 
{ 
	return NC_DENSITY_PROFILE_GET_CLASS (dp)->central_density (dp, cosmo, R, z);
}
