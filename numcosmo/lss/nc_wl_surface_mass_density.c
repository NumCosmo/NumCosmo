/***************************************************************************
 *            nc_wl_surface_mass_density.c
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

/**
 * SECTION:nc_wl_surface_mass_density
 * @title: NcWLSurfaceMassDensity
 * @short_description: Weak lensing surface mass density
 *
 * FIXME
 * 
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_wl_surface_mass_density.h"
#include "nc_distance.h"
#include "lss/nc_density_profile.h"
#include "math/integral.h"
#include "math/ncm_c.h"
#include "math/ncm_cfg.h"
#include "math/ncm_spline_cubic_notaknot.h"

G_DEFINE_TYPE (NcWLSurfaceMassDensity, nc_wl_surface_mass_density, NCM_TYPE_MODEL);

#define VECTOR (NCM_MODEL (smd)->params)
#define PCC   (ncm_vector_get (VECTOR, NC_WL_SURFACE_MASS_DENSITY_PCC))
#define ROFF   (ncm_vector_get (VECTOR, NC_WL_SURFACE_MASS_DENSITY_ROFF))

enum
{
  PROP_0,
	PROP_DENSITY_PROFILE,
	PROP_DISTANCE,
	PROP_ZSOURCE,
	PROP_ZLENS,
	PROP_ZCLUSTER,
  PROP_SIZE,
};


static void
nc_wl_surface_mass_density_init (NcWLSurfaceMassDensity *smd)
{
  smd->dp         = NULL;
	smd->dist       = NULL;
  smd->zsource    = 0.0;
  smd->zlens      = 0.0;
	smd->zcluster   = 0.0;
	smd->ctrl_cosmo = ncm_model_ctrl_new (NULL);
  
}

static void
_nc_wl_surface_mass_density_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcWLSurfaceMassDensity *smd = NC_WL_SURFACE_MASS_DENSITY (object);
  g_return_if_fail (NC_IS_WL_SURFACE_MASS_DENSITY (object));

  switch (prop_id)
  {
    case PROP_DENSITY_PROFILE:
      smd->dp = g_value_dup_object (value);
      break;
		case PROP_DISTANCE:
      smd->dist = g_value_dup_object (value);
      break;
		case PROP_ZSOURCE:
			smd->zsource = g_value_get_double (value);
			break;
		case PROP_ZLENS:
			smd->zlens = g_value_get_double (value);
			break;
		case PROP_ZCLUSTER:
			smd->zcluster = g_value_get_double (value);
			break;	
		default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  } 
}

static void
_nc_wl_surface_mass_density_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcWLSurfaceMassDensity *smd = NC_WL_SURFACE_MASS_DENSITY (object);
  g_return_if_fail (NC_IS_WL_SURFACE_MASS_DENSITY (object));

  switch (prop_id)
  {
    case PROP_DISTANCE:
      g_value_set_object (value, smd->dist);
      break;
		case PROP_DENSITY_PROFILE:
		  g_value_set_object (value, smd->dp);
      break;
		case PROP_ZSOURCE:
			g_value_set_double (value, smd->zsource);
			break;
		case PROP_ZLENS:
			g_value_set_double (value, smd->zlens);
			break;
		case PROP_ZCLUSTER:
			g_value_set_double (value, smd->zcluster);
			break;	
		default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_wl_surface_mass_density_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_wl_surface_mass_density_parent_class)->constructed (object);
  {
     NcWLSurfaceMassDensity *smd = NC_WL_SURFACE_MASS_DENSITY (object);

		g_assert_cmpfloat (smd->zlens, <, smd->zsource);
  }
}

static void
_nc_wl_surface_mass_density_dispose (GObject *object)
{
  NcWLSurfaceMassDensity *smd = NC_WL_SURFACE_MASS_DENSITY (object);

	nc_density_profile_clear (&smd->dp);
	nc_distance_clear (&smd->dist);
    
  ncm_model_ctrl_clear (&smd->ctrl_cosmo);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_wl_surface_mass_density_parent_class)->dispose (object);
}

static void
_nc_wl_surface_mass_density_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_wl_surface_mass_density_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_wl_surface_mass_density, NC_TYPE_WL_SURFACE_MASS_DENSITY);

static void
nc_wl_surface_mass_density_class_init (NcWLSurfaceMassDensityClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
	NcmModelClass* model_class = NCM_MODEL_CLASS (klass);

  object_class->constructed  = &_nc_wl_surface_mass_density_constructed;
  model_class->set_property = &_nc_wl_surface_mass_density_set_property;
  model_class->get_property = &_nc_wl_surface_mass_density_get_property;
  object_class->dispose      = &_nc_wl_surface_mass_density_dispose;
  object_class->finalize     = &_nc_wl_surface_mass_density_finalize;

  ncm_model_class_set_name_nick (model_class, "WL Surface Mass Density", "WLSMD");
  ncm_model_class_add_params (model_class, NC_WL_SURFACE_MASS_DENSITY_SPARAM_LEN, 0, PROP_SIZE);
	
  	/**
   * NcWLSurfaceMassDensity:pcc:
   * 
   * Percentage of correctly centered clusters.
   * Interval: [0.0, 1.0]
   */
  ncm_model_class_set_sparam (model_class, NC_WL_SURFACE_MASS_DENSITY_PCC, "p_{cc}", "pcc",
                              0.0,  1.0, 5e-2,
                              NC_WL_SURFACE_MASS_DENSITY_DEFAULT_PARAMS_ABSTOL, NC_WL_SURFACE_MASS_DENSITY_DEFAULT_PCC,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcWLSurfaceMassDensity:Roff:
   * 
   * Scale length of the miscentering probability distribution. 
   * FIXME Set correct values (limits) Units: Mpc/h
   */
  ncm_model_class_set_sparam (model_class, NC_WL_SURFACE_MASS_DENSITY_ROFF, "R_{off}", "Roff",
                              0.0,  2.0, 1e-1,
                              NC_WL_SURFACE_MASS_DENSITY_DEFAULT_PARAMS_ABSTOL, NC_WL_SURFACE_MASS_DENSITY_DEFAULT_ROFF,
                              NCM_PARAM_TYPE_FIXED);

	/**
   * NcWLSurfaceMassDensity:density-profile:
   *
   * This property keeps the object NcDensityProfile.
   */
  g_object_class_install_property (object_class,
                                   PROP_DISTANCE,
                                   g_param_spec_object ("density-profile",
                                                        NULL,
                                                        "Density Profile",
                                                        NC_TYPE_DENSITY_PROFILE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME
                                                        | G_PARAM_STATIC_BLURB));
	
	 /**
   * NcWLSurfaceMassDensity:distance:
   *
   * This property keeps the object NcDistance.
   */
  g_object_class_install_property (object_class,
                                   PROP_DISTANCE,
                                   g_param_spec_object ("distance",
                                                        NULL,
                                                        "Distance",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME
                                                        | G_PARAM_STATIC_BLURB));
	/**
   * NcWLSurfaceMassDensity:zsource:
   *
   * This property sets the source's redshift.
   * 
   */
  g_object_class_install_property (object_class,
                                   PROP_ZSOURCE,
                                   g_param_spec_double ("zsource",
                                                        NULL,
                                                        "z_{source}",
                                                        0.0, 6.0, 1.5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	/**
   * NcWLSurfaceMassDensity:zlens:
   *
   * This property sets the lens' redshift.
   * Usually zlens = zcluster, but we define these two properties 
	 * in order to handle cases where shear signal has been rescaled to a different cluster redshift (following D. Applegate's code.). 
	 *
   */
  g_object_class_install_property (object_class,
                                   PROP_ZLENS,
                                   g_param_spec_double ("zlens",
                                                        NULL,
                                                        "z_{lens}",
                                                        0.0, 5.5, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	/**
   * NcWLSurfaceMassDensity:zcluster:
   *
   * This property sets the cluster's redshift.
	 * 
   */
  g_object_class_install_property (object_class,
                                   PROP_ZCLUSTER,
                                   g_param_spec_double ("zcluster",
                                                        NULL,
                                                        "z_{cluster}",
                                                        0.0, 5.5, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  
}

/**
 * nc_wl_surface_mass_density_new:
 * @dp: a #NcDensityProfile
 * @dist: a #NcDistance
 * 
 * This function allocates memory for a new #NcWLSurfaceMassDensity object and sets its properties to the values from
 * the input arguments.
 *
 * Returns: a new #NcWLSurfaceMassDensity
 */
NcWLSurfaceMassDensity *
nc_wl_surface_mass_density_new (NcDensityProfile *dp, NcDistance *dist)
{
	NcWLSurfaceMassDensity *smd = g_object_new (NC_TYPE_WL_SURFACE_MASS_DENSITY, 
                       "density-profile", dp,
                       "distance", dist, 
                       NULL);
  return smd;
}

/**
 * nc_wl_surface_mass_density_ref:
 * @smd: a #NcWLSurfaceMassDensity
 *
 * Increases the reference count of @smd by one.
 *
 * Returns: (transfer full): @smd
 */
NcWLSurfaceMassDensity *
nc_wl_surface_mass_density_ref (NcWLSurfaceMassDensity *smd)
{
  return g_object_ref (smd);
}

/**
 * nc_wl_surface_mass_density_free:
 * @smd: a #NcWLSurfaceMassDensity
 *
 * Atomically decrements the reference count of @smd by one. If the reference count drops to 0,
 * all memory allocated by @smd is released.
 *
 */
void
nc_wl_surface_mass_density_free (NcWLSurfaceMassDensity *smd)
{
  g_object_unref (smd);
}

/**
 * nc_wl_surface_mass_density_clear:
 * @smd: a #NcWLSurfaceMassDensity
 *
 * Atomically decrements the reference count of @smd by one. If the reference count drops to 0,
 * all memory allocated by @smd is released. Set pointer to NULL.
 *
 */
void
nc_wl_surface_mass_density_clear (NcWLSurfaceMassDensity **smd)
{
  g_clear_object (smd);
}


/**
 * nc_wl_surface_mass_density_prepare:
 * @smd: a #NcWLSurfaceMassDensity
 * @cosmo: a #NcHICosmo
 *
 * FIXME
 *
 */
void
nc_wl_surface_mass_density_prepare (NcWLSurfaceMassDensity *smd, NcHICosmo *cosmo)
{
  /*IMPLEMENT FIXME */
  ncm_model_ctrl_update (smd->ctrl_cosmo, NCM_MODEL (cosmo));

  return;
}

/**
 * nc_wl_surface_mass_density_sigma_critical:
 * @smd: a #NcWLSurfaceMassDensity
 * @cosmo: a #NcHICosmo
 *
 * Computes the critical surface density, 
 * \begin{equation}\label{eq:def:SigmaC}
 * \Sigma_c = \frac{c^2}{4\pi G} \frac{D_s}{D_l D_{ls})},
 * \end{equation}
 * where $c^2$ is the speed of light squared [ncm_c_c2 ()], $G$ is the gravitational constant in units of $m^3/s^2 M_\odot^{-1}$ [ncm_c_G_mass_solar()],
 * $D_s$ ($D_l$) is the angular diameter distance from the observer to the source (lens), and $D_{sl}$ is the angular diameter distance between 
 * the lens and the source.
 *
 * Returns: the critical surface density $\Sigma_c$ in units of $M_\odot / Mpc^2$
 */
gdouble nc_wl_surface_mass_density_sigma_critical (NcWLSurfaceMassDensity *smd, NcHICosmo *cosmo)
{
	//g_assert_cmpfloat (nc_hicosmo_Omega_k0 (cosmo), >=, 0.0);
	const gdouble a = ncm_c_c2 () / (4.0 * M_PI * ncm_c_G_mass_solar ()) * ncm_c_Mpc (); /* [ M_solar / Mpc ] */
	gdouble Ds = nc_distance_angular_diameter (smd->dist, cosmo, smd->zsource);
	gdouble Dl = nc_distance_angular_diameter (smd->dist, cosmo, smd->zlens);
	/* Dls below is only valid for Omega_k >= 0 */
	gdouble Dls = (nc_distance_transverse (smd->dist, cosmo, smd->zsource) - nc_distance_transverse (smd->dist, cosmo, smd->zlens))/ (1.0 + smd->zsource);
	gdouble RH_Mpc = nc_hicosmo_RH_Mpc (cosmo);
  //gdouble beta = Dls / Ds; 

  //printf ("a = %.5g beta = %.5g\n", a, beta);
  //printf ("Ds = %.10g Dl = %.10g Dls = %.10g\n", Ds * RH_Mpc, Dl * RH_Mpc, Dls* RH_Mpc);
	
  return a * Ds / (Dl * Dls * RH_Mpc);
}

/**
 * nc_wl_surface_mass_density_sigma:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcDensityProfile  
 * @cosmo: a #NcHICosmo
 * @R: FIXME
 *
 * FIXME
 *
 * Returns: FIXME 
 */
gdouble
nc_wl_surface_mass_density_sigma (NcWLSurfaceMassDensity *smd, NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R)
{
	smd->dp = dp;
	gdouble sigma_2 = nc_density_profile_integral_density_los (smd->dp, cosmo, R, smd->zcluster);
		
	return 2.0 * sigma_2;
}

/**
 * nc_wl_surface_mass_density_sigma_mean:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcDensityProfile 
 * @cosmo: a #NcHICosmo
 * @R: FIXME
 *
 * FIXME
 *
 * Returns: FIXME 
 */
gdouble 
nc_wl_surface_mass_density_sigma_mean (NcWLSurfaceMassDensity *smd, NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R)
{
	gdouble rs = nc_density_profile_scale_radius (dp, cosmo, smd->zcluster);
	gdouble x  = R / rs;
  gdouble x2 = x * x;
	
	gdouble mean_sigma_2x2 = nc_density_profile_integral_density_2d (dp, cosmo, R, smd->zcluster);
		
	return (2.0 / x2) * mean_sigma_2x2;
}

/* Correction term: Central point mass */

static gdouble 
_nc_wl_surface_mass_density_cg (NcWLSurfaceMassDensity *smd, gdouble M0, gdouble R)
{
	gdouble R2 = R * R;
	
	return M0 / (M_PI * R2);
}

/* Correction term: no-weak shear */

static gdouble 
_nc_wl_surface_mass_density_nw (NcWLSurfaceMassDensity *smd, gdouble M0, gdouble R)
{
	gdouble R2 = R * R;
	
	return 0.0;
}

/**
 * nc_wl_surface_mass_density_convergence:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcDensityProfile  
 * @cosmo: a #NcHICosmo
 * @R: FIXME
 *
 * Computes the convergence 
 *
 * Returns: $\kappa(R)$
 */
gdouble 
nc_wl_surface_mass_density_convergence (NcWLSurfaceMassDensity *smd, NcDensityProfile *dp, NcHICosmo *cosmo, gdouble R)
{
	smd->dp = dp;
  gdouble sigma      = nc_wl_surface_mass_density_sigma (smd, smd->dp, cosmo, R);
	gdouble sigma_crit = nc_wl_surface_mass_density_sigma_critical (smd, cosmo);

	return sigma / sigma_crit;
}

/**
 * nc_wl_surface_mass_density_shear:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcDensityProfile  
 * @cosmo: a #NcHICosmo
 * @R: FIXME
 *
 * FIXME
 *
 * Returns: $\gamma(R)$ 
 */
gdouble 
nc_wl_surface_mass_density_shear (NcWLSurfaceMassDensity *smd, NcDensityProfile *dp, NcHICosmo *cosmo, gdouble R)
{
  gdouble sigma      = nc_wl_surface_mass_density_sigma (smd, dp, cosmo, R);
	gdouble mean_sigma = nc_wl_surface_mass_density_sigma_mean (smd, dp, cosmo, R);
	gdouble sigma_crit = nc_wl_surface_mass_density_sigma_critical (smd, cosmo);

	return (mean_sigma - sigma) / sigma_crit;
}

/**
 * nc_wl_surface_mass_density_reduced_shear:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcDensityProfile  
 * @cosmo: a #NcHICosmo
 * @R: FIXME
 *
 * FIXME 
 * $$ g(R) = \frac{\gamma(R)}{1 - \kappa(R)}$$ 
 *
 * Returns: $g(R)$
 */
gdouble 
nc_wl_surface_mass_density_reduced_shear (NcWLSurfaceMassDensity *smd, NcDensityProfile *dp, NcHICosmo *cosmo, gdouble R)
{
	/* Optimize it to compute sigma and sigma_c just once*/
  gdouble convergence   = nc_wl_surface_mass_density_convergence (smd, dp, cosmo, R);
	gdouble shear         = nc_wl_surface_mass_density_shear (smd, dp, cosmo, R);
	gdouble reduced_shear = shear / (1.0 - convergence); 

	return reduced_shear;
}

/**
 * nc_wl_surface_mass_density_reduced_shear_infinity:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcDensityProfile  
 * @cosmo: a #NcHICosmo
 * @R: FIXME
 *
 * FIXME 
 *
 * Returns: $g_{\infity}(R)$
 */
gdouble 
nc_wl_surface_mass_density_reduced_shear_infinity (NcWLSurfaceMassDensity *smd, NcDensityProfile *dp, NcHICosmo *cosmo, gdouble R)
{
	/* Optimize it to compute sigma and sigma_c just once, and distances at inf */
  gdouble Ds              = nc_distance_angular_diameter (smd->dist, cosmo, smd->zsource);
	gdouble Dls             = (nc_distance_transverse (smd->dist, cosmo, smd->zsource) - nc_distance_transverse (smd->dist, cosmo, smd->zlens))/ (1.0 + smd->zsource);
  gdouble Dinf            = nc_distance_angular_diameter (smd->dist, cosmo, 1.0e6);
	gdouble Dlinf           = (nc_distance_transverse (smd->dist, cosmo, 1.0e6) - nc_distance_transverse (smd->dist, cosmo, smd->zlens))/ (1.0 + 1.0e6);
	gdouble beta_s          = (Dls / Ds) * (Dinf / Dlinf);
  gdouble convergence_inf = beta_s * nc_wl_surface_mass_density_convergence (smd, dp, cosmo, R);
	gdouble shear_inf       = beta_s * nc_wl_surface_mass_density_shear (smd, dp, cosmo, R);
	
	return shear_inf / (1.0 - convergence_inf);
}
