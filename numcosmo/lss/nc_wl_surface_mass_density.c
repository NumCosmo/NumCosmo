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
 * where $\rho(r)$ is the three-dimensional mass density profile (#NcHaloDensityProfile), $r^2 = R^2 + \chi^2$ is a three-dimensional vector in space, $R$ is a
 * two-dimensional vector from the halo center. In particular, we consider a projection $\Sigma (R)$ onto the lens plane.
 * $\chi$ is the distance along the line of sight.
 *
 * The mean surface mass density within a circular aperture of radius $R$ is, [nc_wl_surface_mass_density_sigma_mean()]
 * \begin{equation}\label{eq:sigma_mean}
 * \overline{\Sigma} (<R) = \frac{2}{R^2} \int_0^R \mathrm{d}R^\prime \, R^\prime \Sigma (R^\prime).
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
#include "lss/nc_halo_density_profile.h"
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
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);
  
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
  if (zs < zc)
    return GSL_POSINF;
  else
  {
    const gdouble a   = ncm_c_c2 () / (4.0 * M_PI * ncm_c_G_mass_solar ()) * ncm_c_Mpc (); /* [ M_solar / Mpc ] */
    const gdouble Ds  = nc_distance_angular_diameter (smd->dist, cosmo, zs);
    const gdouble Dl  = nc_distance_angular_diameter (smd->dist, cosmo, zl);
    const gdouble Dls = nc_distance_angular_diameter_z1_z2 (smd->dist, cosmo, zl, zs);

    const gdouble RH_Mpc = nc_hicosmo_RH_Mpc (cosmo);

    return a * Ds / (Dl * Dls * RH_Mpc);
  }
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
  /*g_assert_cmpfloat (nc_hicosmo_Omega_k0 (cosmo), >=, 0.0); */
  const gdouble a     = ncm_c_c2 () / (4.0 * M_PI * ncm_c_G_mass_solar ()) * ncm_c_Mpc (); /* [ M_solar / Mpc ] */
  const gdouble Dinf  = nc_distance_transverse_z_to_infinity (smd->dist, cosmo, 0.0);
  const gdouble Dl    = nc_distance_angular_diameter (smd->dist, cosmo, zl);
  const gdouble Dlinf = nc_distance_transverse_z_to_infinity (smd->dist, cosmo, zl);
  
  const gdouble RH_Mpc = nc_hicosmo_RH_Mpc (cosmo);
  
  return a * Dinf / (Dl * Dlinf * RH_Mpc);
}

/**
 * nc_wl_surface_mass_density_sigma:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcHaloDensityProfile
 * @cosmo: a #NcHICosmo
 * @R: projected radius with respect to the center of the lens / halo
 * @zc: cluster redshift $z_\mathrm{cluster}$
 *
 * Computes the surface mass density at @R, see Eq. $\eqref{eq:sigma}$.
 *
 * Returns: $\Sigma (R)$
 */
gdouble
nc_wl_surface_mass_density_sigma (NcWLSurfaceMassDensity *smd, NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble zc)
{
  const gdouble Sigma = nc_halo_density_profile_eval_2d_density (dp, cosmo, R, zc);
  
  return Sigma;
}

/**
 * nc_wl_surface_mass_density_sigma_mean:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcHaloDensityProfile
 * @cosmo: a #NcHICosmo
 * @R: projected radius with respect to the center of the lens / halo
 * @zc: cluster redshift $z_\mathrm{cluster}$
 *
 * Computes the mean surface mass density inside the circle with radius @R, Eq. $\eqref{eq:sigma_mean}$.
 *
 * Returns: $\overline{\Sigma} (<R)$
 */
gdouble
nc_wl_surface_mass_density_sigma_mean (NcWLSurfaceMassDensity *smd, NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble zc)
{
  const gdouble mean_sigma_2x2 = nc_halo_density_profile_eval_cyl_mass (dp, cosmo, R, zc) / (ncm_c_pi () * R * R);

  return mean_sigma_2x2;
}

/**
 * nc_wl_surface_mass_density_sigma_excess:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcHaloDensityProfile
 * @cosmo: a #NcHICosmo
 * @R: projected radius with respect to the center of the lens / halo
 * @zc: cluster redshift $z_\mathrm{cluster}$
 *
 * Computes difference between the mean surface mass density inside the circle with radius @R (Eq. $\eqref{eq:sigma_mean}$) and the surface mass density at @R (Eq. $\eqref{eq:sigma}$).
 *
 * Returns: $\overline{\Sigma} (<R) - \Sigma (R)$
 */
gdouble
nc_wl_surface_mass_density_sigma_excess (NcWLSurfaceMassDensity *smd, NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble zc)
{
  const gdouble mean_sigma_2x2 = nc_wl_surface_mass_density_sigma_mean (smd, dp, cosmo, R, zc);
  const gdouble sigma          = nc_wl_surface_mass_density_sigma (smd, dp, cosmo, R, zc);
  
  return (mean_sigma_2x2 - sigma);
}

/* Correction term: Central point mass */

/*
 *  static gdouble
 *  _nc_wl_surface_mass_density_cg (NcWLSurfaceMassDensity *smd, gdouble M0, const gdouble R)
 *  {
 *  gdouble R2 = R * R;
 *
 *  return M0 / (M_PI * R2);
 *  }
 */

/**
 * nc_wl_surface_mass_density_convergence:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcHaloDensityProfile
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
nc_wl_surface_mass_density_convergence (NcWLSurfaceMassDensity *smd, NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble zs, const gdouble zl, const gdouble zc)
{
  if (zs < zc)
    return 0.0;
  else
  {
    gdouble sigma      = nc_wl_surface_mass_density_sigma (smd, dp, cosmo, R, zc);
    gdouble sigma_crit = nc_wl_surface_mass_density_sigma_critical (smd, cosmo, zs, zl, zc);
  
    return sigma / sigma_crit;
  }
}

/**
 * nc_wl_surface_mass_density_convergence_infinity:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcHaloDensityProfile
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
nc_wl_surface_mass_density_convergence_infinity (NcWLSurfaceMassDensity *smd, NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble zl, const gdouble zc)
{
  gdouble sigma          = nc_wl_surface_mass_density_sigma (smd, dp, cosmo, R, zc);
  gdouble sigma_crit_inf = nc_wl_surface_mass_density_sigma_critical_infinity (smd, cosmo, zl, zc);
  
  return sigma / sigma_crit_inf;
}

/**
 * nc_wl_surface_mass_density_shear:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcHaloDensityProfile
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
nc_wl_surface_mass_density_shear (NcWLSurfaceMassDensity *smd, NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble zs, const gdouble zl, const gdouble zc)
{
  if (zs < zc)
    return 0.0;
  else
  {
    const gdouble sigma      = nc_wl_surface_mass_density_sigma (smd, dp, cosmo, R, zc);
    const gdouble mean_sigma = nc_wl_surface_mass_density_sigma_mean (smd, dp, cosmo, R, zc);
    const gdouble sigma_crit = nc_wl_surface_mass_density_sigma_critical (smd, cosmo, zs, zl, zc);
  
    return (mean_sigma - sigma) / sigma_crit;
  }
}

/**
 * nc_wl_surface_mass_density_shear_infinity:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcHaloDensityProfile
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
nc_wl_surface_mass_density_shear_infinity (NcWLSurfaceMassDensity *smd, NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble zl, const gdouble zc)
{
  gdouble sigma          = nc_wl_surface_mass_density_sigma (smd, dp, cosmo, R, zc);
  gdouble mean_sigma     = nc_wl_surface_mass_density_sigma_mean (smd, dp, cosmo, R, zc);
  gdouble sigma_crit_inf = nc_wl_surface_mass_density_sigma_critical_infinity (smd, cosmo, zl, zc);
  
  return (mean_sigma - sigma) / sigma_crit_inf;
}

/**
 * nc_wl_surface_mass_density_reduced_shear:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcHaloDensityProfile
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
nc_wl_surface_mass_density_reduced_shear (NcWLSurfaceMassDensity *smd, NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble zs, const gdouble zl, const gdouble zc)
{
  gdouble r_s, rho_s;

  nc_halo_density_profile_r_s_rho_s (dp, cosmo, zc, &r_s, &rho_s);

  {
    const gdouble a             = ncm_c_c2 () / (4.0 * M_PI * ncm_c_G_mass_solar ()) * ncm_c_Mpc () / nc_hicosmo_RH_Mpc (cosmo); /* [ M_solar / Mpc^2 ] */
    const gdouble Omega_k0      = nc_hicosmo_Omega_k0 (cosmo);
    const gdouble sqrt_Omega_k0 = sqrt (fabs (Omega_k0));
    const gint k                = fabs (Omega_k0) < NCM_ZERO_LIMIT ? 0 : (Omega_k0 > 0.0 ? -1 : 1);

    switch (k)
    {
      case -1:
        {
          const gdouble X          = R / r_s;
          const gdouble mean_sigma = (2.0 * nc_halo_density_profile_eval_dl_cyl_mass (dp, X) / (X * X)) * rho_s * r_s;
          const gdouble sigma      = (nc_halo_density_profile_eval_dl_2d_density (dp, X)) * rho_s * r_s;
          const gdouble dl         = nc_distance_comoving (smd->dist, cosmo, zl);
          const gdouble Dl         = sinh (sqrt_Omega_k0 * dl) / ((1.0 + zl) * sqrt_Omega_k0);
          if (zs > zl)
          {
            const gdouble x_i           = 1.0 + zs;
            const gdouble ds            = nc_distance_comoving (smd->dist, cosmo, zs);
            const gdouble Ds            = sinh (sqrt_Omega_k0 * ds)        / (x_i * sqrt_Omega_k0);
            const gdouble Dls           = sinh (sqrt_Omega_k0 * (ds - dl)) / (x_i * sqrt_Omega_k0);
            const gdouble sigma_crit    = a * Ds / (Dl * Dls);
            const gdouble mu            = sigma / sigma_crit;
            const gdouble shear         = (mean_sigma - sigma) / sigma_crit;
            const gdouble reduced_shear = shear / (1.0 - mu);

            return reduced_shear;
          }
          else
            return 0.0;
        }
        break;
      case 0:
        {
          const gdouble X          = R / r_s;
          const gdouble mean_sigma = (2.0 * nc_halo_density_profile_eval_dl_cyl_mass (dp, X) / (X * X)) * rho_s * r_s;
          const gdouble sigma      = (nc_halo_density_profile_eval_dl_2d_density (dp, X)) * rho_s * r_s;
          const gdouble dl         = nc_distance_comoving (smd->dist, cosmo, zl);
          const gdouble Dl         = dl / (1.0 + zl);
          if (zs > zl)
          {
            const gdouble x_i           = 1.0 + zs;
            const gdouble ds            = nc_distance_comoving (smd->dist, cosmo, zs);
            const gdouble Ds            = ds        / x_i;
            const gdouble Dls           = (ds - dl) / x_i;
            const gdouble sigma_crit    = a * Ds / (Dl * Dls);
            const gdouble mu            = sigma / sigma_crit;
            const gdouble shear         = (mean_sigma - sigma) / sigma_crit;
            const gdouble reduced_shear = shear / (1.0 - mu);

            return reduced_shear;
          }
          else
            return 0.0;
        }
        break;
      case 1:
        {
          const gdouble X          = R / r_s;
          const gdouble mean_sigma = (2.0 * nc_halo_density_profile_eval_dl_cyl_mass (dp, X) / (X * X)) * rho_s * r_s;
          const gdouble sigma      = (nc_halo_density_profile_eval_dl_2d_density (dp, X)) * rho_s * r_s;
          const gdouble dl         = nc_distance_comoving (smd->dist, cosmo, zl);
          const gdouble Dl         = sin (sqrt_Omega_k0 * dl) / ((1.0 + zl) * sqrt_Omega_k0);
          if (zs > zl)
          {
            const gdouble x_i           = 1.0 + zs;
            const gdouble ds            = nc_distance_comoving (smd->dist, cosmo, zs);
            const gdouble Ds            = sin (sqrt_Omega_k0 * ds)        / (x_i * sqrt_Omega_k0);
            const gdouble Dls           = sin (sqrt_Omega_k0 * (ds - dl)) / (x_i * sqrt_Omega_k0);
            const gdouble sigma_crit    = a * Ds / (Dl * Dls);
            const gdouble mu            = sigma / sigma_crit;
            const gdouble shear         = (mean_sigma - sigma) / sigma_crit;
            const gdouble reduced_shear = shear / (1.0 - mu);

            return reduced_shear;
          }
          else
            return 0.0;
        }
        break;
      default:
        g_assert_not_reached ();
        break;
    }
  }

  /* Old un-optimized code */
  if (zs < zc)
	return 0.0;
  else
  {
    gdouble convergence   = nc_wl_surface_mass_density_convergence (smd, dp, cosmo, R, zs, zl, zc);
    gdouble shear         = nc_wl_surface_mass_density_shear (smd, dp, cosmo, R, zs, zl, zc);
    gdouble reduced_shear = shear / (1.0 - convergence);
  
    return reduced_shear;
  }
}

/**
 * nc_wl_surface_mass_density_reduced_shear_optzs_prep: (skip)
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcHaloDensityProfile
 * @cosmo: a #NcHICosmo
 * @R: projected radius with respect to the center of the lens / halo
 * @zl: lens redshift $z_\mathrm{lens}$
 * @zc: cluster redshift $z_\mathrm{cluster}$
 * @optzs: a #NcWLSurfaceMassDensityOptzs
 *
 * Computes the reduced shear:
 * $$ g(R) = \frac{\gamma(R)}{1 - \kappa(R)},$$
 * where $\gamma(R)$ is the shear [nc_wl_surface_mass_density_shear()] and $\kappa(R)$ is the convergence
 * [nc_wl_surface_mass_density_convergence()].
 *
 * FIXME
 *
 */
void
nc_wl_surface_mass_density_reduced_shear_optzs_prep (NcWLSurfaceMassDensity *smd, NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble zl, const gdouble zc, NcWLSurfaceMassDensityOptzs *optzs)
{
  gdouble r_s, rho_s;

  nc_halo_density_profile_r_s_rho_s (dp, cosmo, zc, &r_s, &rho_s);

  {
    const gdouble a             = ncm_c_c2 () / (4.0 * M_PI * ncm_c_G_mass_solar ()) * ncm_c_Mpc () / nc_hicosmo_RH_Mpc (cosmo); /* [ M_solar / Mpc^2 ] */
    const gdouble Omega_k0      = nc_hicosmo_Omega_k0 (cosmo);

    optzs->sqrt_Omega_k0 = sqrt (fabs (Omega_k0));

    optzs->k  = fabs (Omega_k0) < NCM_ZERO_LIMIT ? 0 : (Omega_k0 > 0.0 ? -1 : 1);
    optzs->dl = nc_distance_comoving (smd->dist, cosmo, zl);

    switch (optzs->k)
    {
      case -1:
        {
          const gdouble X          = R / r_s;
          const gdouble Dl         = sinh (optzs->sqrt_Omega_k0 * optzs->dl) / ((1.0 + zl) * optzs->sqrt_Omega_k0);

          optzs->mean_sigma = (2.0 * nc_halo_density_profile_eval_dl_cyl_mass (dp, X) / (X * X)) * rho_s * r_s;
          optzs->sigma      = (nc_halo_density_profile_eval_dl_2d_density (dp, X)) * rho_s * r_s;
          optzs->sc_Dls_Ds  = a / Dl;
        }
        break;
      case 0:
        {
          const gdouble X          = R / r_s;
          const gdouble Dl         = optzs->dl / (1.0 + zl);

          optzs->mean_sigma = (2.0 * nc_halo_density_profile_eval_dl_cyl_mass (dp, X) / (X * X)) * rho_s * r_s;
          optzs->sigma      = (nc_halo_density_profile_eval_dl_2d_density (dp, X)) * rho_s * r_s;
          optzs->sc_Dls_Ds  = a / Dl;
        }
        break;
      case 1:
        {
          const gdouble X          = R / r_s;
          const gdouble Dl         = sin (optzs->sqrt_Omega_k0 * optzs->dl) / ((1.0 + zl) * optzs->sqrt_Omega_k0);

          optzs->mean_sigma = (2.0 * nc_halo_density_profile_eval_dl_cyl_mass (dp, X) / (X * X)) * rho_s * r_s;
          optzs->sigma      = (nc_halo_density_profile_eval_dl_2d_density (dp, X)) * rho_s * r_s;
          optzs->sc_Dls_Ds  = a / Dl;
        }
        break;
      default:
        g_assert_not_reached ();
        break;
    }
  }
}

/**
 * nc_wl_surface_mass_density_reduced_shear_optzs:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcHaloDensityProfile
 * @cosmo: a #NcHICosmo
 * @zs: source redshift $z_\mathrm{source}$
 * @optzs: a #NcWLSurfaceMassDensityOptzs
 *
 * Computes the reduced shear:
 * $$ g(R) = \frac{\gamma(R)}{1 - \kappa(R)},$$
 * where $\gamma(R)$ is the shear [nc_wl_surface_mass_density_shear()] and $\kappa(R)$ is the convergence
 * [nc_wl_surface_mass_density_convergence()].
 *
 * FIXME
 *
 * Returns: $g(R)$
 */
gdouble
nc_wl_surface_mass_density_reduced_shear_optzs (NcWLSurfaceMassDensity *smd, NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble zs, const gdouble zl, NcWLSurfaceMassDensityOptzs *optzs)
{
  if (zs < zl)
    return 0.0;

  switch (optzs->k)
  {
    case -1:
    {
      const gdouble x_i           = 1.0 + zs;
      const gdouble ds            = nc_distance_comoving (smd->dist, cosmo, zs);
      const gdouble Ds            = sinh (optzs->sqrt_Omega_k0 * ds)               / (x_i * optzs->sqrt_Omega_k0);
      const gdouble Dls           = sinh (optzs->sqrt_Omega_k0 * (ds - optzs->dl)) / (x_i * optzs->sqrt_Omega_k0);
      const gdouble sigma_crit    = optzs->sc_Dls_Ds * Ds / Dls;
      const gdouble mu            = optzs->sigma / sigma_crit;
      const gdouble shear         = (optzs->mean_sigma - optzs->sigma) / sigma_crit;
      const gdouble reduced_shear = shear / (1.0 - mu);

      return reduced_shear;
    }
    break;
    case 0:
    {
      const gdouble x_i           = 1.0 + zs;
      const gdouble ds            = nc_distance_comoving (smd->dist, cosmo, zs);
      const gdouble Ds            = ds               / x_i;
      const gdouble Dls           = (ds - optzs->dl) / x_i;
      const gdouble sigma_crit    = optzs->sc_Dls_Ds * Ds / Dls;
      const gdouble mu            = optzs->sigma / sigma_crit;
      const gdouble shear         = (optzs->mean_sigma - optzs->sigma) / sigma_crit;
      const gdouble reduced_shear = shear / (1.0 - mu);

      return reduced_shear;
    }
    break;
    case 1:
    {
      const gdouble x_i           = 1.0 + zs;
      const gdouble ds            = nc_distance_comoving (smd->dist, cosmo, zs);
      const gdouble Ds            = sin (optzs->sqrt_Omega_k0 * ds)               / (x_i * optzs->sqrt_Omega_k0);
      const gdouble Dls           = sin (optzs->sqrt_Omega_k0 * (ds - optzs->dl)) / (x_i * optzs->sqrt_Omega_k0);
      const gdouble sigma_crit    = optzs->sc_Dls_Ds * Ds / Dls;
      const gdouble mu            = optzs->sigma / sigma_crit;
      const gdouble shear         = (optzs->mean_sigma - optzs->sigma) / sigma_crit;
      const gdouble reduced_shear = shear / (1.0 - mu);

      return reduced_shear;
    }
    break;
    default:
      g_assert_not_reached ();
      break;
  }
}

/**
 * nc_wl_surface_mass_density_reduced_shear_infinity:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcHaloDensityProfile
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
nc_wl_surface_mass_density_reduced_shear_infinity (NcWLSurfaceMassDensity *smd, NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble zs, const gdouble zl, const gdouble zc)
{
  if (zs < zc)
    return 0.0;
  else
  {
    gdouble Dinf, betainf, beta_s;
    gdouble convergence_inf, shear_inf, g;
    gdouble Ds   = nc_distance_angular_diameter (smd->dist, cosmo, zs);
    gdouble Dls  = (nc_distance_transverse (smd->dist, cosmo, zs) - nc_distance_transverse (smd->dist, cosmo, zl)) / (1.0 + zs);
    gdouble beta = Dls / Ds;

    Dinf =  nc_distance_transverse_z_to_infinity (smd->dist, cosmo, 0.0);
    /* Angular diamater distance lens to infinity over angular diameter distance obs to infinity. It's equivalent to compute the ratio of the transverse distances. Source at infinity */
    betainf = nc_distance_transverse_z_to_infinity (smd->dist, cosmo, zl) / Dinf;
    beta_s  = beta / betainf;

    convergence_inf = beta_s * nc_wl_surface_mass_density_convergence_infinity (smd, dp, cosmo, R, zl, zc);
    shear_inf       = beta_s * nc_wl_surface_mass_density_shear_infinity (smd, dp, cosmo, R, zl, zc);

    g = shear_inf / (1.0 - convergence_inf);

    return g;
  }
}

/**
 * nc_wl_surface_mass_density_magnification:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcHaloDensityProfile
 * @cosmo: a #NcHICosmo
 * @R: projected radius with respect to the center of the lens / halo
 * @zs: source redshift $z_\mathrm{source}$
 * @zl: lens redshift $z_\mathrm{lens}$
 * @zc: cluster redshift $z_\mathrm{cluster}$
 *
 * Computes the reduced shear:
 * $$ \mu(R) = \frac{1}{(1 - \kappa(R))^2 - \vert\gamma^2(R) \vert},$$
 * where $\gamma(R)$ is the shear [nc_wl_surface_mass_density_shear()] and $\kappa(R)$ is the convergence
 * [nc_wl_surface_mass_density_convergence()].
 *
 * Returns: $\mu(R)$
 */
gdouble
nc_wl_surface_mass_density_magnification (NcWLSurfaceMassDensity *smd, NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble zs, const gdouble zl, const gdouble zc)
{
  if (zs < zc)
	return 1.0;
  else
  {
    gdouble convergence   = nc_wl_surface_mass_density_convergence (smd, dp, cosmo, R, zs, zl, zc);
    gdouble shear         = nc_wl_surface_mass_density_shear (smd, dp, cosmo, R, zs, zl, zc);
    gdouble one_m_conv    = 1.0 - convergence;
    gdouble shear2        = shear * shear;
    gdouble magnification = 1.0 / (one_m_conv * one_m_conv - shear2);
  
    return magnification;
  }
}

/**
 * nc_wl_surface_mass_density_sigma_array:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcHaloDensityProfile
 * @cosmo: a #NcHICosmo
 * @R: (element-type gdouble): projected radius with respect to the center of the lens / halo
 * @fin: factor to multiply $R$, it should be $1$ or the appropriated unit conversion
 * @fout: factor to multiply $\rho(R)$, it should be $1$ or the appropriated unit conversion
 * @zc: cluster redshift $z_\mathrm{cluster}$
 *
 * Computes the surface mass density at @R, see Eq. $\eqref{eq:sigma}$.
 *
 * Returns: (transfer full) (element-type gdouble): $\Sigma (R)$
 */
GArray *
nc_wl_surface_mass_density_sigma_array (NcWLSurfaceMassDensity *smd, NcHaloDensityProfile *dp, NcHICosmo *cosmo, GArray *R, gdouble fin, gdouble fout, const gdouble zc)
{
  return nc_halo_density_profile_eval_2d_density_array (dp, cosmo, R, fin, fout, zc);
}

/**
 * nc_wl_surface_mass_density_sigma_excess_array:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcHaloDensityProfile
 * @cosmo: a #NcHICosmo
 * @R: (element-type gdouble): projected radius with respect to the center of the lens / halo
 * @fin: factor to multiply $R$, it should be $1$ or the appropriated unit conversion
 * @fout: factor to multiply $\rho(R)$, it should be $1$ or the appropriated unit conversion
 * @zc: cluster redshift $z_\mathrm{cluster}$
 *
 * Computes difference between the mean surface mass density inside the circle with radius @R (Eq. $\eqref{eq:sigma_mean}$) and the surface mass density at @R (Eq. $\eqref{eq:sigma}$).
 *
 * Returns: (transfer full) (element-type gdouble): $\overline{\Sigma} (<R) - \Sigma (R)$
 */
GArray *
nc_wl_surface_mass_density_sigma_excess_array (NcWLSurfaceMassDensity *smd, NcHaloDensityProfile *dp, NcHICosmo *cosmo, GArray *R, gdouble fin, gdouble fout, const gdouble zc)
{
  gdouble r_s, rho_s;
  
  g_assert_cmpint (R->len, >, 0);
  
  nc_halo_density_profile_r_s_rho_s (dp, cosmo, zc, &r_s, &rho_s);
  
  fin  = fin / r_s;
  fout = fout * rho_s * r_s;
  
  {
    GArray *res = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), R->len);
    guint i;
    
    g_array_set_size (res, R->len);
    
    for (i = 0; i < R->len; i++)
    {
      const gdouble X        = g_array_index (R, gdouble, i) * fin;
      const gdouble barSigma = 2.0 * nc_halo_density_profile_eval_dl_cyl_mass (dp, X) / (X * X);
      const gdouble sigma    = nc_halo_density_profile_eval_dl_2d_density (dp, X);
      
      g_array_index (res, gdouble, i) = fout * (barSigma - sigma);
    }
    
    return res;
  }
}

/**
 * nc_wl_surface_mass_density_reduced_shear_array:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcHaloDensityProfile
 * @cosmo: a #NcHICosmo
 * @R: (element-type gdouble): projected radius with respect to the center of the lens / halo
 * @fin: factor to multiply $R$, it should be $1$ or the appropriated unit conversion
 * @fout: factor to multiply $g(R)$, it should be $1$ or the appropriated unit conversion
 * @zs: (element-type gdouble): source redshift $z_\mathrm{source}$
 * @zl: lens redshift $z_\mathrm{lens}$
 * @zc: cluster redshift $z_\mathrm{cluster}$
 *
 * Computes the reduced shear:
 * $$ g(R) = \frac{\gamma(R)}{1 - \kappa(R)},$$
 * where $\gamma(R)$ is the shear [nc_wl_surface_mass_density_shear()] and $\kappa(R)$ is the convergence
 * [nc_wl_surface_mass_density_convergence()].
 *
 * Returns: (transfer full) (element-type gdouble): $g(R)$
 */
GArray *
nc_wl_surface_mass_density_reduced_shear_array (NcWLSurfaceMassDensity *smd, NcHaloDensityProfile *dp, NcHICosmo *cosmo, GArray *R, gdouble fin, gdouble fout, GArray *zs, const gdouble zl, const gdouble zc)
{
  gdouble r_s, rho_s;

  g_assert_cmpint (R->len, >, 0);
  g_assert_cmpint (zs->len, >, 0);

  nc_halo_density_profile_r_s_rho_s (dp, cosmo, zc, &r_s, &rho_s);

  fin  = fin / r_s;

  {
    GArray *res                 = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), R->len * zs->len);
    const gdouble a             = ncm_c_c2 () / (4.0 * M_PI * ncm_c_G_mass_solar ()) * ncm_c_Mpc () / nc_hicosmo_RH_Mpc (cosmo); /* [ M_solar / Mpc^2 ] */
    const gdouble Omega_k0      = nc_hicosmo_Omega_k0 (cosmo);
    const gdouble sqrt_Omega_k0 = sqrt (fabs (Omega_k0));
    const gint k                = fabs (Omega_k0) < NCM_ZERO_LIMIT ? 0 : (Omega_k0 > 0.0 ? -1 : 1);
    gint i,j;

    g_array_set_size (res, R->len * zs->len);

    switch (k)
    {
      case -1:
        for (i = 0; i < R->len; i++)
        {
          const gdouble X          = g_array_index (R, gdouble, i) * fin;
          const gdouble mean_sigma = (2.0 * nc_halo_density_profile_eval_dl_cyl_mass (dp, X) / (X * X)) * rho_s * r_s;
          const gdouble sigma      = (nc_halo_density_profile_eval_dl_2d_density (dp, X)) * rho_s * r_s;
          const gdouble dl         = nc_distance_comoving (smd->dist, cosmo, zl);
          const gdouble Dl         = sinh (sqrt_Omega_k0 * dl) / ((1.0 + zl) * sqrt_Omega_k0);
          const gint pad           = zs->len * i;

          for (j = 0; j < zs->len; j++)
          {
            const gdouble zs_j          = g_array_index (zs, gdouble, j);
            if (zs_j > zl)
            {
              const gdouble x_j           = 1.0 + zs_j;
              const gdouble ds            = nc_distance_comoving (smd->dist, cosmo, zs_j);
              const gdouble Ds            = sinh (sqrt_Omega_k0 * ds)        / (x_j * sqrt_Omega_k0);
              const gdouble Dls           = sinh (sqrt_Omega_k0 * (ds - dl)) / (x_j * sqrt_Omega_k0);
              const gdouble sigma_crit    = a * Ds / (Dl * Dls);
              const gdouble mu            = sigma / sigma_crit;
              const gdouble shear         = (mean_sigma - sigma) / sigma_crit;
              const gdouble reduced_shear = shear / (1.0 - mu);

              g_array_index (res, gdouble, pad + j) = reduced_shear * fout;
            }
            else
              g_array_index (res, gdouble, pad + j) = 0.0;
          }
        }
        break;
      case 0:
        for (i = 0; i < R->len; i++)
        {
          const gdouble X          = g_array_index (R, gdouble, i) * fin;
          const gdouble mean_sigma = (2.0 * nc_halo_density_profile_eval_dl_cyl_mass (dp, X) / (X * X)) * rho_s * r_s;
          const gdouble sigma      = (nc_halo_density_profile_eval_dl_2d_density (dp, X)) * rho_s * r_s;
          const gdouble dl         = nc_distance_comoving (smd->dist, cosmo, zl);
          const gdouble Dl         = dl / (1.0 + zl);
          const gint pad           = zs->len * i;

          for (j = 0; j < zs->len; j++)
          {
            const gdouble zs_j          = g_array_index (zs, gdouble, j);
            if (zs_j > zl)
            {
              const gdouble x_j           = 1.0 + zs_j;
              const gdouble ds            = nc_distance_comoving (smd->dist, cosmo, zs_j);
              const gdouble Ds            = ds        / x_j;
              const gdouble Dls           = (ds - dl) / x_j;
              const gdouble sigma_crit    = a * Ds / (Dl * Dls);
              const gdouble mu            = sigma / sigma_crit;
              const gdouble shear         = (mean_sigma - sigma) / sigma_crit;
              const gdouble reduced_shear = shear / (1.0 - mu);

              g_array_index (res, gdouble, pad + j) = reduced_shear * fout;
            }
            else
              g_array_index (res, gdouble, pad + j) = 0.0;
          }
        }
        break;
      case 1:
        for (i = 0; i < R->len; i++)
        {
          const gdouble X          = g_array_index (R, gdouble, i) * fin;
          const gdouble mean_sigma = (2.0 * nc_halo_density_profile_eval_dl_cyl_mass (dp, X) / (X * X)) * rho_s * r_s;
          const gdouble sigma      = (nc_halo_density_profile_eval_dl_2d_density (dp, X)) * rho_s * r_s;
          const gdouble dl         = nc_distance_comoving (smd->dist, cosmo, zl);
          const gdouble Dl         = sin (sqrt_Omega_k0 * dl) / ((1.0 + zl) * sqrt_Omega_k0);
          const gint pad           = zs->len * i;

          for (j = 0; j < zs->len; j++)
          {
            const gdouble zs_j          = g_array_index (zs, gdouble, j);
            if (zs_j > zl)
            {
              const gdouble x_j           = 1.0 + zs_j;
              const gdouble ds            = nc_distance_comoving (smd->dist, cosmo, zs_j);
              const gdouble Ds            = sin (sqrt_Omega_k0 * ds)        / (x_j * sqrt_Omega_k0);
              const gdouble Dls           = sin (sqrt_Omega_k0 * (ds - dl)) / (x_j * sqrt_Omega_k0);
              const gdouble sigma_crit    = a * Ds / (Dl * Dls);
              const gdouble mu            = sigma / sigma_crit;
              const gdouble shear         = (mean_sigma - sigma) / sigma_crit;
              const gdouble reduced_shear = shear / (1.0 - mu);

              g_array_index (res, gdouble, pad + j) = reduced_shear * fout;
            }
            else
              g_array_index (res, gdouble, pad + j) = 0.0;
          }
        }
        break;
      default:
        g_assert_not_reached ();
        break;
    }

    return res;
  }
}

/**
 * nc_wl_surface_mass_density_reduced_shear_array_equal:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcHaloDensityProfile
 * @cosmo: a #NcHICosmo
 * @R: (element-type gdouble): projected radius with respect to the center of the lens / halo
 * @fin: factor to multiply $R$, it should be $1$ or the appropriated unit conversion
 * @fout: factor to multiply $g(R)$, it should be $1$ or the appropriated unit conversion
 * @zs: (element-type gdouble): source redshift $z_\mathrm{source}$
 * @zl: lens redshift $z_\mathrm{lens}$
 * @zc: cluster redshift $z_\mathrm{cluster}$
 *
 * Computes the reduced shear:
 * $$ g(R) = \frac{\gamma(R)}{1 - \kappa(R)},$$
 * where $\gamma(R)$ is the shear [nc_wl_surface_mass_density_shear()] and $\kappa(R)$ is the convergence
 * [nc_wl_surface_mass_density_convergence()].
 *
 * Returns: (transfer full) (element-type gdouble): $g(R)$
 */
GArray *
nc_wl_surface_mass_density_reduced_shear_array_equal (NcWLSurfaceMassDensity *smd, NcHaloDensityProfile *dp, NcHICosmo *cosmo, GArray *R, gdouble fin, gdouble fout, GArray *zs, const gdouble zl, const gdouble zc)
{
  gdouble r_s, rho_s;

  g_assert_cmpint (R->len, ==, zs->len);
  g_assert_cmpint (zs->len, >, 0);

  nc_halo_density_profile_r_s_rho_s (dp, cosmo, zc, &r_s, &rho_s);

  fin  = fin / r_s;

  {
    GArray *res                 = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), R->len);
    const gdouble a             = ncm_c_c2 () / (4.0 * M_PI * ncm_c_G_mass_solar ()) * ncm_c_Mpc () / nc_hicosmo_RH_Mpc (cosmo); /* [ M_solar / Mpc^2 ] */
    const gdouble Omega_k0      = nc_hicosmo_Omega_k0 (cosmo);
    const gdouble sqrt_Omega_k0 = sqrt (fabs (Omega_k0));
    const gint k                = fabs (Omega_k0) < NCM_ZERO_LIMIT ? 0 : (Omega_k0 > 0.0 ? -1 : 1);
    gint i;

    g_array_set_size (res, R->len);

    switch (k)
    {
      case -1:
        for (i = 0; i < R->len; i++)
        {
          const gdouble X          = g_array_index (R, gdouble, i) * fin;
          const gdouble mean_sigma = (2.0 * nc_halo_density_profile_eval_dl_cyl_mass (dp, X) / (X * X)) * rho_s * r_s;
          const gdouble sigma      = (nc_halo_density_profile_eval_dl_2d_density (dp, X)) * rho_s * r_s;
          const gdouble dl         = nc_distance_comoving (smd->dist, cosmo, zl);
          const gdouble Dl         = sinh (sqrt_Omega_k0 * dl) / ((1.0 + zl) * sqrt_Omega_k0);
          const gdouble zs_i       = g_array_index (zs, gdouble, i);
          if (zs_i > zl)
          {
            const gdouble x_i           = 1.0 + zs_i;
            const gdouble ds            = nc_distance_comoving (smd->dist, cosmo, zs_i);
            const gdouble Ds            = sinh (sqrt_Omega_k0 * ds)        / (x_i * sqrt_Omega_k0);
            const gdouble Dls           = sinh (sqrt_Omega_k0 * (ds - dl)) / (x_i * sqrt_Omega_k0);
            const gdouble sigma_crit    = a * Ds / (Dl * Dls);
            const gdouble mu            = sigma / sigma_crit;
            const gdouble shear         = (mean_sigma - sigma) / sigma_crit;
            const gdouble reduced_shear = shear / (1.0 - mu);

            g_array_index (res, gdouble, i) = reduced_shear * fout;
          }
          else
            g_array_index (res, gdouble, i) = 0.0;
        }
        break;
      case 0:
        for (i = 0; i < R->len; i++)
        {
          const gdouble X          = g_array_index (R, gdouble, i) * fin;
          const gdouble mean_sigma = (2.0 * nc_halo_density_profile_eval_dl_cyl_mass (dp, X) / (X * X)) * rho_s * r_s;
          const gdouble sigma      = (nc_halo_density_profile_eval_dl_2d_density (dp, X)) * rho_s * r_s;
          const gdouble dl         = nc_distance_comoving (smd->dist, cosmo, zl);
          const gdouble Dl         = dl / (1.0 + zl);
          const gdouble zs_i       = g_array_index (zs, gdouble, i);
          if (zs_i > zl)
          {
            const gdouble x_i           = 1.0 + zs_i;
            const gdouble ds            = nc_distance_comoving (smd->dist, cosmo, zs_i);
            const gdouble Ds            = ds        / x_i;
            const gdouble Dls           = (ds - dl) / x_i;
            const gdouble sigma_crit    = a * Ds / (Dl * Dls);
            const gdouble mu            = sigma / sigma_crit;
            const gdouble shear         = (mean_sigma - sigma) / sigma_crit;
            const gdouble reduced_shear = shear / (1.0 - mu);

            g_array_index (res, gdouble, i) = reduced_shear * fout;
          }
          else
            g_array_index (res, gdouble, i) = 0.0;
        }
        break;
      case 1:
        for (i = 0; i < R->len; i++)
        {
          const gdouble X          = g_array_index (R, gdouble, i) * fin;
          const gdouble mean_sigma = (2.0 * nc_halo_density_profile_eval_dl_cyl_mass (dp, X) / (X * X)) * rho_s * r_s;
          const gdouble sigma      = (nc_halo_density_profile_eval_dl_2d_density (dp, X)) * rho_s * r_s;
          const gdouble dl         = nc_distance_comoving (smd->dist, cosmo, zl);
          const gdouble Dl         = sin (sqrt_Omega_k0 * dl) / ((1.0 + zl) * sqrt_Omega_k0);
          const gdouble zs_i       = g_array_index (zs, gdouble, i);
          if (zs_i > zl)
          {
            const gdouble x_i           = 1.0 + zs_i;
            const gdouble ds            = nc_distance_comoving (smd->dist, cosmo, zs_i);
            const gdouble Ds            = sin (sqrt_Omega_k0 * ds)        / (x_i * sqrt_Omega_k0);
            const gdouble Dls           = sin (sqrt_Omega_k0 * (ds - dl)) / (x_i * sqrt_Omega_k0);
            const gdouble sigma_crit    = a * Ds / (Dl * Dls);
            const gdouble mu            = sigma / sigma_crit;
            const gdouble shear         = (mean_sigma - sigma) / sigma_crit;
            const gdouble reduced_shear = shear / (1.0 - mu);

            g_array_index (res, gdouble, i) = reduced_shear * fout;
          }
          else
            g_array_index (res, gdouble, i) = 0.0;
        }
        break;
      default:
        g_assert_not_reached ();
        break;
    }

    return res;
  }
}

