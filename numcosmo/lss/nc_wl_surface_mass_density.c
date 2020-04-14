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
 * This object implements the projected surface mass density and related quantities, such as the convergence and tangential shear.
 * 
 * The projected surface mass density is [nc_wl_surface_mass_density_sigma()]  
 * \begin{equation}\label{eq:sigma}
 * \Sigma (R) = \int \mathrm{d}\chi \, \rho\left(\sqrt{R^2 + \chi^2} \right), 
 * \end{equation} 
 * where $\rho(r)$ is the three-dimensional mass density profile (#NcDensityProfile), $r^2 = R^2 + \chi^2$ is a three-dimensional vector in space, $R$ is a 
 * two-dimensional vector from the halo center. In particular, we consider a projection $\Sigma (R)$ onto the lens plane. 
 * $\chi$ is the distance along the line of sight. 
 * 
 * The mean surface mass density within a circular aperture of radius $R$ is, [nc_wl_surface_mass_density_sigma_mean()]
 * \begin{equation}\label{eq:sigma_mean}
 * \overline{\Sigma} (<R) = \frac{2}{R^2} \int_0^R \mathrm{d}R^\prime \, \Sigma (R^\prime).
 * \end{equation}
 * 
 * The convergence $\kappa (R)$ [nc_wl_surface_mass_density_convergence()] and the shear $\gamma(R)$ [nc_wl_surface_mass_density_shear()] 
 * are given by, respectively,
 * \begin{equation}\label{eq:convergence}
 * \kappa (R) = \frac{\Sigma (R)}{\Sigma_{crit}},
 * \end{equation}
 * * \begin{equation}\label{eq:shear}
 * \gamma (R) = \frac{\Delta\Sigma (R)}{\Sigma_{crit}} = \frac{\overline{\Sigma} (<R) - \Sigma (R)}{\Sigma_{crit}},
 * \end{equation}
 * where $\Sigma_{crit}$ is the critical surface density [nc_wl_surface_mass_density_sigma_critical()],
 * \begin{equation}\label{eq:sigma_critical}
 * \Sigma_{crit} = \frac{c^2}{4\pi G} \frac{D_s}{D_l D_{ls}}.
 * \end{equation}
 * where $c^2$ is the speed of light squared [ncm_c_c2()], $G$ is the gravitational constant [ncm_c_G()], $D_s$ and $Dl$ are the angular diameter distances 
 * to the source and lens, respectively, and $D_{ls}$ is the angular diameter distance between the lens and source.  
 * 
 * See, e.g., [Mandelbaum (2006)][XMandelbaum2006], [Umetsu (2012)][XUmetsu2012], [Applegate (2014)][XApplegate2014], 
 * [Melchior (2017)][XMelchior2017], [Parroni (2017)][XParroni2017].
 * 
 * Usually $z_{lens} = z_{cluster}, but we define these as two different arguments in order to handle cases where shear signal has been 
 * rescaled to a different cluster redshift (following D. Applegate's code.).
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
	PROP_DISTANCE,
	PROP_ZSOURCE,
	PROP_ZLENS,
	PROP_ZCLUSTER,
  PROP_SIZE,
};


static void
nc_wl_surface_mass_density_init (NcWLSurfaceMassDensity *smd)
{
	smd->dist       = NULL;
	smd->ctrl_cosmo = ncm_model_ctrl_new (NULL);
	smd->ctrl_dp    = ncm_model_ctrl_new (NULL);
}

static void
_nc_wl_surface_mass_density_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcWLSurfaceMassDensity *smd = NC_WL_SURFACE_MASS_DENSITY (object);
  g_return_if_fail (NC_IS_WL_SURFACE_MASS_DENSITY (object));

  switch (prop_id)
  {
    case PROP_DISTANCE:
      smd->dist = g_value_dup_object (value);
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
		/*NcWLSurfaceMassDensity *smd = NC_WL_SURFACE_MASS_DENSITY (object);*/

		
  }
}

static void
_nc_wl_surface_mass_density_dispose (GObject *object)
{
  NcWLSurfaceMassDensity *smd = NC_WL_SURFACE_MASS_DENSITY (object);

	nc_distance_clear (&smd->dist);

	ncm_model_ctrl_clear (&smd->ctrl_dp);
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

  object_class->constructed = &_nc_wl_surface_mass_density_constructed;
  model_class->set_property = &_nc_wl_surface_mass_density_set_property;
  model_class->get_property = &_nc_wl_surface_mass_density_get_property;
  object_class->dispose     = &_nc_wl_surface_mass_density_dispose;
  object_class->finalize    = &_nc_wl_surface_mass_density_finalize;

  ncm_mset_model_register_id (model_class,
                              "NcWLSurfaceMassDensity",
                              "NcWLSurfaceMassDensity.",
                              NULL,
                              FALSE,
                              NCM_MSET_MODEL_MAIN);
	
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
                              NC_WL_SURFACE_MASS_DENSITY_DEFAULT_PARAMS_ABSTOL, NC_WL_SURFACE_MASS_DENSITY_DEFAULT_PCC, NCM_PARAM_TYPE_FIXED);

  /**
   * NcWLSurfaceMassDensity:Roff:
   * 
   * Scale length of the miscentering probability distribution. 
   * FIXME Set correct values (limits) Units: Mpc
   */
  ncm_model_class_set_sparam (model_class, NC_WL_SURFACE_MASS_DENSITY_ROFF, "R_{off}", "Roff",
                              0.0,  2.0, 1e-1,
                              NC_WL_SURFACE_MASS_DENSITY_DEFAULT_PARAMS_ABSTOL, NC_WL_SURFACE_MASS_DENSITY_DEFAULT_ROFF, NCM_PARAM_TYPE_FIXED);

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
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));  
}

/**
 * nc_wl_surface_mass_density_new:
 * @dist: a #NcDistance
 * 
 * This function allocates memory for a new #NcWLSurfaceMassDensity object and sets its properties to the values from
 * the input arguments.
 *
 * Returns: a new #NcWLSurfaceMassDensity
 */
NcWLSurfaceMassDensity *
nc_wl_surface_mass_density_new (NcDistance *dist)
{
	NcWLSurfaceMassDensity *smd = g_object_new (NC_TYPE_WL_SURFACE_MASS_DENSITY, 
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
 * Prepares the #NcWLSurfaceMassDensity object @smd for computation.
 *
 */
void 
nc_wl_surface_mass_density_prepare (NcWLSurfaceMassDensity *smd, NcHICosmo *cosmo)
{
  nc_distance_prepare_if_needed (smd->dist, cosmo);
}

/**
 * nc_wl_surface_mass_density_prepare_if_needed:
 * @smd: a #NcWLSurfaceMassDensity
 * @cosmo: a #NcHICosmo
 * 
 * Prepares the #NcWLSurfaceMassDensity object @smd for computation if necessary.
 *
 */
void 
nc_wl_surface_mass_density_prepare_if_needed (NcWLSurfaceMassDensity *smd, NcHICosmo *cosmo)
{
  if (ncm_model_ctrl_update (smd->ctrl_cosmo, NCM_MODEL (cosmo)))
    nc_wl_surface_mass_density_prepare (smd, cosmo);
}

/**
 * nc_wl_surface_mass_density_sigma_critical:
 * @smd: a #NcWLSurfaceMassDensity
 * @cosmo: a #NcHICosmo
 * @zs: source redshift $z_\mathrm{source}$  
 * @zl: lens redshift $z_\mathrm{lens}$
 * @zc: cluster redshift $z_\mathrm{cluster}$
 *
 * Computes the critical surface density, 
 * \begin{equation}\label{eq:def:SigmaC}
 * \Sigma_c = \frac{c^2}{4\pi G} \frac{D_s}{D_l D_{ls}},
 * \end{equation}
 * where $c^2$ is the speed of light squared [ncm_c_c2 ()], $G$ is the gravitational constant in units of $m^3/s^2 M_\odot^{-1}$ [ncm_c_G_mass_solar()],
 * $D_s$ ($D_l$) is the angular diameter distance from the observer to the source (lens), and $D_{ls}$ is the angular diameter distance between 
 * the lens and the source.
 *
 * Returns: the critical surface density $\Sigma_c$ in units of $M_\odot / Mpc^2$
 */
gdouble 
nc_wl_surface_mass_density_sigma_critical (NcWLSurfaceMassDensity *smd, NcHICosmo *cosmo, const gdouble zs, const gdouble zl, const gdouble zc)
{
	//g_assert_cmpfloat (nc_hicosmo_Omega_k0 (cosmo), >=, 0.0);
	const gdouble a      = ncm_c_c2 () / (4.0 * M_PI * ncm_c_G_mass_solar ()) * ncm_c_Mpc (); /* [ M_solar / Mpc ] */
	const gdouble Ds     = nc_distance_angular_diameter (smd->dist, cosmo, zs);
	const gdouble Dl     = nc_distance_angular_diameter (smd->dist, cosmo, zl);
	/* Dls below is only valid for Omega_k >= 0 */
	const gdouble Dls    = (nc_distance_transverse (smd->dist, cosmo, zs) - nc_distance_transverse (smd->dist, cosmo, zl)) / (1.0 + zs);
	const gdouble RH_Mpc = nc_hicosmo_RH_Mpc (cosmo);
  
  return a * Ds / (Dl * Dls * RH_Mpc);
}

/**
 * nc_wl_surface_mass_density_sigma_critical_infinity:
 * @smd: a #NcWLSurfaceMassDensity
 * @cosmo: a #NcHICosmo
 * @zl: lens redshift $z_\mathrm{lens}$
 * @zc: cluster redshift $z_\mathrm{cluster}$
 *
 * Computes the critical surface density, 
 * \begin{equation}\label{eq:def:SigmaC}
 * \Sigma_c = \frac{c^2}{4\pi G} \frac{D_\infty}{D_l D_{l\infty}},
 * \end{equation}
 * where $c^2$ is the speed of light squared [ncm_c_c2 ()], $G$ is the gravitational constant in units of $m^3/s^2 M_\odot^{-1}$ [ncm_c_G_mass_solar()],
 * $D_\infty$ ($D_l$) is the angular diameter distance from the observer to the source at infinite redshift (lens), and $D_{l\infty}$ is the angular diameter distance between 
 * the lens and the source.
 *
 * Returns: the critical surface density $\Sigma_c$ in units of $M_\odot / Mpc^2$
 */
gdouble 
nc_wl_surface_mass_density_sigma_critical_infinity (NcWLSurfaceMassDensity *smd, NcHICosmo *cosmo, const gdouble zl, const gdouble zc)
{
	//g_assert_cmpfloat (nc_hicosmo_Omega_k0 (cosmo), >=, 0.0);
	const gdouble a      = ncm_c_c2 () / (4.0 * M_PI * ncm_c_G_mass_solar ()) * ncm_c_Mpc (); /* [ M_solar / Mpc ] */
	const gdouble Dinf     = nc_distance_transverse_z_to_infinity (smd->dist, cosmo, 0.0);
	const gdouble Dl     = nc_distance_angular_diameter (smd->dist, cosmo, zl);
	const gdouble Dlinf = nc_distance_transverse_z_to_infinity (smd->dist, cosmo, zl);
	
	const gdouble RH_Mpc = nc_hicosmo_RH_Mpc (cosmo);
  
  return a * Dinf / (Dl * Dlinf * RH_Mpc);
}

/**
 * nc_wl_surface_mass_density_sigma:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcDensityProfile  
 * @cosmo: a #NcHICosmo
 * @R: projected radius with respect to the center of the lens / halo
 * @zc: cluster redshift $z_\mathrm{cluster}$
 *
 * Computes the surface mass density at @R, see Eq. $\eqref{eq:sigma}$.
 *
 * Returns: $\Sigma (R)$
 */
gdouble
nc_wl_surface_mass_density_sigma (NcWLSurfaceMassDensity *smd, NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble zc)
{
	gdouble sigma_2 = nc_density_profile_integral_density_los (dp, cosmo, R, zc);
		
	return 2.0 * sigma_2;
}

/**
 * nc_wl_surface_mass_density_sigma_mean:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcDensityProfile 
 * @cosmo: a #NcHICosmo
 * @R: projected radius with respect to the center of the lens / halo
 * @zc: cluster redshift $z_\mathrm{cluster}$
 *
 * Computes the mean surface mass density inside the circle with radius @R, Eq. $\eqref{eq:sigma_mean}$.
 *
 * Returns: $\overline{\Sigma} (<R)$
 */
gdouble 
nc_wl_surface_mass_density_sigma_mean (NcWLSurfaceMassDensity *smd, NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble zc)
{
	const gdouble rs             = nc_density_profile_scale_radius (dp, cosmo, zc);
	const gdouble x              = R / rs;
  const gdouble x2             = x * x;
	const gdouble mean_sigma_2x2 = nc_density_profile_integral_density_2d (dp, cosmo, R, zc);

	return (2.0 / x2) * mean_sigma_2x2;
}

/* Correction term: Central point mass */

/*
static gdouble 
_nc_wl_surface_mass_density_cg (NcWLSurfaceMassDensity *smd, gdouble M0, const gdouble R)
{
	gdouble R2 = R * R;
	
	return M0 / (M_PI * R2);
}
*/

/**
 * nc_wl_surface_mass_density_convergence:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcDensityProfile  
 * @cosmo: a #NcHICosmo
 * @R: projected radius with respect to the center of the lens / halo
 * @zs: source redshift $z_\mathrm{source}$  
 * @zl: lens redshift $z_\mathrm{lens}$
 * @zc: cluster redshift $z_\mathrm{cluster}$
 *
 * Computes the convergence $\kappa(R)$ at @R, see Eq $\eqref{eq:convergence}$.
 *
 * Returns: $\kappa(R)$
 */
gdouble 
nc_wl_surface_mass_density_convergence (NcWLSurfaceMassDensity *smd, NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble zs, const gdouble zl, const gdouble zc)
{
  gdouble sigma      = nc_wl_surface_mass_density_sigma (smd, dp, cosmo, R, zc);
	gdouble sigma_crit = nc_wl_surface_mass_density_sigma_critical (smd, cosmo, zs, zl, zc);

	return sigma / sigma_crit;
}

/**
 * nc_wl_surface_mass_density_convergence_infinity:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcDensityProfile  
 * @cosmo: a #NcHICosmo
 * @R: projected radius with respect to the center of the lens / halo
 * @zl: lens redshift $z_\mathrm{lens}$
 * @zc: cluster redshift $z_\mathrm{cluster}$
 *
 * Computes the convergence $\kappa_\infty(R)$ at @R, see Eq $\eqref{eq:convergence}$, and sources at infinite redshift.
 *
 * Returns: $\kappa_\infty(R)$
 */
gdouble 
nc_wl_surface_mass_density_convergence_infinity (NcWLSurfaceMassDensity *smd, NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble zl, const gdouble zc)
{
  gdouble sigma      = nc_wl_surface_mass_density_sigma (smd, dp, cosmo, R, zc);
	gdouble sigma_crit_inf = nc_wl_surface_mass_density_sigma_critical_infinity (smd, cosmo, zl, zc);

	return sigma / sigma_crit_inf;
}

/**
 * nc_wl_surface_mass_density_shear:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcDensityProfile  
 * @cosmo: a #NcHICosmo
 * @R: projected radius with respect to the center of the lens / halo
 * @zs: source redshift $z_\mathrm{source}$  
 * @zl: lens redshift $z_\mathrm{lens}$
 * @zc: cluster redshift $z_\mathrm{cluster}$
 *
 * Computes the shear $\gamma(R)$ at @R, see Eq $\eqref{eq:shear}$. 
 *
 * Returns: $\gamma(R)$ 
 */
gdouble 
nc_wl_surface_mass_density_shear (NcWLSurfaceMassDensity *smd, NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble zs, const gdouble zl, const gdouble zc)
{
  const gdouble sigma      = nc_wl_surface_mass_density_sigma (smd, dp, cosmo, R, zc);
	const gdouble mean_sigma = nc_wl_surface_mass_density_sigma_mean (smd, dp, cosmo, R, zc);
	const gdouble sigma_crit = nc_wl_surface_mass_density_sigma_critical (smd, cosmo, zs, zl, zc);

	return (mean_sigma - sigma) / sigma_crit;
}

/**
 * nc_wl_surface_mass_density_shear_infinity:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcDensityProfile  
 * @cosmo: a #NcHICosmo
 * @R: projected radius with respect to the center of the lens / halo 
 * @zl: lens redshift $z_\mathrm{lens}$
 * @zc: cluster redshift $z_\mathrm{cluster}$
 *
 * Computes the shear $\gamma_\infty (R)$ at @R, see Eq $\eqref{eq:shear}$, and source at infinite redshift. 
 *
 * Returns: $\gamma_\infty (R)$ 
 */
gdouble 
nc_wl_surface_mass_density_shear_infinity (NcWLSurfaceMassDensity *smd, NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble zl, const gdouble zc)
{
  gdouble sigma      = nc_wl_surface_mass_density_sigma (smd, dp, cosmo, R, zc);
	gdouble mean_sigma = nc_wl_surface_mass_density_sigma_mean (smd, dp, cosmo, R, zc);
	gdouble sigma_crit_inf = nc_wl_surface_mass_density_sigma_critical_infinity (smd, cosmo, zl, zc);

	return (mean_sigma - sigma) / sigma_crit_inf;
}

/**
 * nc_wl_surface_mass_density_reduced_shear:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcDensityProfile  
 * @cosmo: a #NcHICosmo
 * @R: projected radius with respect to the center of the lens / halo
 * @zs: source redshift $z_\mathrm{source}$  
 * @zl: lens redshift $z_\mathrm{lens}$
 * @zc: cluster redshift $z_\mathrm{cluster}$
 *
 * Computes the reduced shear:
 * $$ g(R) = \frac{\gamma(R)}{1 - \kappa(R)},$$ 
 * where $\gamma(R)$ is the shear [nc_wl_surface_mass_density_shear()] and $\kappa(R)$ is the convergence
 * [nc_wl_surface_mass_density_convergence()].	 
 * 
 * Returns: $g(R)$
 */
gdouble 
nc_wl_surface_mass_density_reduced_shear (NcWLSurfaceMassDensity *smd, NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble zs, const gdouble zl, const gdouble zc)
{
	/* Optimize it to compute sigma and sigma_c just once*/
  gdouble convergence   = nc_wl_surface_mass_density_convergence (smd, dp, cosmo, R, zs, zl, zc);
	gdouble shear         = nc_wl_surface_mass_density_shear (smd, dp, cosmo, R, zs, zl, zc);
	gdouble reduced_shear = shear / (1.0 - convergence); 

	return reduced_shear;
}

/**
 * nc_wl_surface_mass_density_reduced_shear_infinity:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcDensityProfile  
 * @cosmo: a #NcHICosmo
 * @R: projected radius with respect to the center of the lens / halo
 * @zs: source redshift $z_\mathrm{source}$  
 * @zl: lens redshift $z_\mathrm{lens}$
 * @zc: cluster redshift $z_\mathrm{cluster}$
 *
 * Computes the reduced shear assuming a lensed source at infinite redshift:
 * $$ g(R) = \frac{\beta_s(zb)\gamma(R)}{1 - \beta_s(zb) \kappa(R)}, $$
 * where $\gamma(R)$ is the shear [nc_wl_surface_mass_density_shear()], $\kappa(R)$ is the convergence
 * [nc_wl_surface_mass_density_convergence()], $z_b$ is the background-galaxy redshift and 
 * $$\beta_s = \frac{D_s}{D_l D_{ls}} \frac{D_\infty}{D_l D_{l\infty}}.$$ 
 * See [Applegate (2014)][XApplegate2014] 
 * 
 * Returns: $g(R)$, source at $z = \infty$
 */
gdouble 
nc_wl_surface_mass_density_reduced_shear_infinity (NcWLSurfaceMassDensity *smd, NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble zs, const gdouble zl, const gdouble zc)
{
	/* Optimize it to compute sigma and sigma_c just once, and distances at inf */
  gdouble Dinf, betainf, beta_s;
	gdouble convergence_inf, shear_inf, g;
	gdouble Ds              = nc_distance_angular_diameter (smd->dist, cosmo, zs);
	gdouble Dls             = (nc_distance_transverse (smd->dist, cosmo, zs) - nc_distance_transverse (smd->dist, cosmo, zl))/ (1.0 + zs);
  gdouble beta            = Dls / Ds;
	
 	Dinf            =  nc_distance_transverse_z_to_infinity (smd->dist, cosmo, 0.0);
	/* Angular diamater distance lens to infinity over angular diameter distance obs to infinity. It's equivalent to compute the ratio of the transverse distances. Source at infinity */
	betainf         = nc_distance_transverse_z_to_infinity (smd->dist, cosmo, zl) / Dinf;
	beta_s          = beta / betainf;

	convergence_inf = beta_s * nc_wl_surface_mass_density_convergence_infinity (smd, dp, cosmo, R, zl, zc);
	shear_inf       = beta_s * nc_wl_surface_mass_density_shear_infinity (smd, dp, cosmo, R, zl, zc);
	
	g = shear_inf / (1.0 - convergence_inf);
	
	return g;
}
