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
 * NcWLSurfaceMassDensity:
 *
 * Weak lensing surface mass density.
 *
 * This object implements the projected surface mass density and related
 * weak-lensing observables — convergence and tangential shear — from a halo
 * density profile (#NcHaloDensityProfile). The projected surface mass density is
 * \begin{equation*}
 * \Sigma(R) = \int \mathrm{d}\chi\,\rho\!\left(\sqrt{R^2 + \chi^2}\right),
 * \end{equation*}
 * and the convergence and shear follow as $\kappa = \Sigma/\Sigma_\mathrm{crit}$
 * and $\gamma = \Delta\Sigma/\Sigma_\mathrm{crit}$, with $\Sigma_\mathrm{crit}$
 * the critical surface density set by the lens/source geometry (#NcDistance).
 *
 * For the full definitions, the mean surface density, and the critical surface
 * density, see the theoretical background page:
 * <a href="../../theory/wl_surface_mass_density.html">Weak-Lensing Surface Mass Density</a>.
 *
 * Usually $z_\mathrm{lens} = z_\mathrm{cluster}$, but these are kept as separate
 * arguments to handle cases where the shear signal has been rescaled to a
 * different cluster redshift (following D. Applegate's code).
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/wl/nc_wl_surface_mass_density.h"
#include "nc/background/nc_distance.h"
#include "nc/lss/halo/nc_halo_density_profile.h"
#include "ncm/integration/ncm_integrate.h"
#include "ncm/core/ncm_c.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/spline/ncm_spline_cubic_notaknot.h"

G_DEFINE_TYPE (NcWLSurfaceMassDensity, nc_wl_surface_mass_density, NCM_TYPE_MODEL)

#define VECTOR (NCM_MODEL (smd))
#define PCC    (ncm_model_orig_param_get (VECTOR, NC_WL_SURFACE_MASS_DENSITY_PCC))
#define ROFF   (ncm_model_orig_param_get (VECTOR, NC_WL_SURFACE_MASS_DENSITY_ROFF))

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
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
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
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
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
   * Scale length of the miscentering probability distribution (units: Mpc).
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
  return nc_distance_sigma_critical (smd->dist, cosmo, zs, zl);
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
  return nc_distance_sigma_critical_infinity (smd->dist, cosmo, zl);
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
  {
    return 0.0;
  }
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
  {
    return 0.0;
  }
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
        {
          return 0.0;
        }
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
        {
          return 0.0;
        }
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
        {
          return 0.0;
        }
      }
      break;
      default:
        g_assert_not_reached ();
        break;
    }
  }

  /* Old un-optimized code */
  if (zs < zc)
  {
    return 0.0;
  }
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
 *
 */
void
nc_wl_surface_mass_density_reduced_shear_optzs_prep (NcWLSurfaceMassDensity *smd, NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble zl, const gdouble zc, NcWLSurfaceMassDensityOptzs *optzs)
{
  NcWLSurfaceMassDensityLensCtx ctx;
  gdouble r_s, rho_s;

  nc_wl_surface_mass_density_lens_ctx_prep (&ctx, smd, cosmo, zl);
  nc_halo_density_profile_r_s_rho_s (dp, cosmo, zc, &r_s, &rho_s);
  nc_wl_surface_mass_density_reduced_shear_optzs_prep_with_lens_ctx (dp, &ctx, R, r_s, rho_s, optzs);
}

/**
 * nc_wl_surface_mass_density_reduced_shear_optzs_prep_with_lens_ctx: (skip)
 * @dp: a #NcHaloDensityProfile
 * @ctx: a #NcWLSurfaceMassDensityLensCtx already prepared via #nc_wl_surface_mass_density_lens_ctx_prep()
 * @R: projected radius with respect to the center of the lens / halo
 * @r_s: scale radius from #nc_halo_density_profile_r_s_rho_s()
 * @rho_s: scale density from #nc_halo_density_profile_r_s_rho_s()
 * @optzs: (out): a #NcWLSurfaceMassDensityOptzs
 *
 * Variant of #nc_wl_surface_mass_density_reduced_shear_optzs_prep() that reuses
 * a pre-computed lens-side context and pre-computed (r_s, rho_s). Only the per-R
 * surface-mass-density profile evaluations are computed per call.
 *
 */
void
nc_wl_surface_mass_density_reduced_shear_optzs_prep_with_lens_ctx (NcHaloDensityProfile *dp, const NcWLSurfaceMassDensityLensCtx *ctx, const gdouble R, const gdouble r_s, const gdouble rho_s, NcWLSurfaceMassDensityOptzs *optzs)
{
  const gdouble X = R / r_s;

  optzs->k             = ctx->k;
  optzs->sqrt_Omega_k0 = ctx->sqrt_Omega_k0;
  optzs->dl            = ctx->dl;
  optzs->sc_Dls_Ds     = ctx->sc_over_Dl;
  optzs->mean_sigma    = (2.0 * nc_halo_density_profile_eval_dl_cyl_mass (dp, X) / (X * X)) * rho_s * r_s;
  optzs->sigma         = (nc_halo_density_profile_eval_dl_2d_density (dp, X)) * rho_s * r_s;
}

/**
 * nc_wl_surface_mass_density_reduced_shear_optzs:
 * @smd: a #NcWLSurfaceMassDensity
 * @dp: a #NcHaloDensityProfile
 * @cosmo: a #NcHICosmo
 * @zs: source redshift $z_\mathrm{source}$
 * @zl: lens redshift $z_\mathrm{lens}$
 * @optzs: a #NcWLSurfaceMassDensityOptzs
 *
 * Computes the reduced shear:
 * $$ g(R) = \frac{\gamma(R)}{1 - \kappa(R)},$$
 * where $\gamma(R)$ is the shear [nc_wl_surface_mass_density_shear()] and $\kappa(R)$ is the convergence
 * [nc_wl_surface_mass_density_convergence()].
 *
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
 * nc_wl_surface_mass_density_reduced_shear_crit_cache_prep: (skip)
 * @smd: a #NcWLSurfaceMassDensity
 * @cosmo: a #NcHICosmo
 * @zl: lens redshift $z_\mathrm{lens}$
 * @zc: cluster redshift $z_\mathrm{cluster}$
 * @zs: source redshift $z_\mathrm{source}$
 * @crit_cache: (out): a #NcWLSurfaceMassDensityCritCache
 *
 * Computes and stores the critical surface density.
 *
 */
void
nc_wl_surface_mass_density_reduced_shear_crit_cache_prep (NcWLSurfaceMassDensity *smd, NcHICosmo *cosmo, const gdouble zl, const gdouble zc, const gdouble zs, NcWLSurfaceMassDensityCritCache *crit_cache)
{
  NCM_UNUSED (zc);
  NcWLSurfaceMassDensityLensCtx ctx;

  nc_wl_surface_mass_density_lens_ctx_prep (&ctx, smd, cosmo, zl);
  nc_wl_surface_mass_density_reduced_shear_crit_cache_prep_with_lens_ctx (smd, cosmo, &ctx, zs, crit_cache);
}

/**
 * nc_wl_surface_mass_density_lens_ctx_prep: (skip)
 * @ctx: a #NcWLSurfaceMassDensityLensCtx
 * @smd: a #NcWLSurfaceMassDensity
 * @cosmo: a #NcHICosmo
 * @zl: lens redshift $z_\mathrm{lens}$
 *
 * Precomputes the lens-side constants (curvature, comoving distance, sigma_crit
 * prefactor) used by repeated #nc_wl_surface_mass_density_reduced_shear_crit_cache_prep_with_lens_ctx()
 * calls that share the same @cosmo and @zl. Hoists cosmology getters and one
 * comoving-distance evaluation out of an inner-loop call site.
 *
 */
void
nc_wl_surface_mass_density_lens_ctx_prep (NcWLSurfaceMassDensityLensCtx *ctx, NcWLSurfaceMassDensity *smd, NcHICosmo *cosmo, const gdouble zl)
{
  const gdouble a             = ncm_c_c2 () / (4.0 * M_PI * ncm_c_G_mass_solar ()) * ncm_c_Mpc () / nc_hicosmo_RH_Mpc (cosmo);
  const gdouble Omega_k0      = nc_hicosmo_Omega_k0 (cosmo);
  const gdouble sqrt_Omega_k0 = sqrt (fabs (Omega_k0));
  const gint k                = fabs (Omega_k0) < NCM_ZERO_LIMIT ? 0 : (Omega_k0 > 0.0 ? -1 : 1);
  const gdouble dl            = nc_distance_comoving (smd->dist, cosmo, zl);
  gdouble Dl;

  switch (k)
  {
    case -1:
      Dl = sinh (sqrt_Omega_k0 * dl) / ((1.0 + zl) * sqrt_Omega_k0);
      break;
    case 0:
      Dl = dl / (1.0 + zl);
      break;
    case 1:
      Dl = sin (sqrt_Omega_k0 * dl) / ((1.0 + zl) * sqrt_Omega_k0);
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  ctx->k             = k;
  ctx->sqrt_Omega_k0 = sqrt_Omega_k0;
  ctx->zl            = zl;
  ctx->dl            = dl;
  ctx->sc_over_Dl    = a / Dl;
}

/**
 * nc_wl_surface_mass_density_reduced_shear_crit_cache_prep_with_lens_ctx: (skip)
 * @smd: a #NcWLSurfaceMassDensity
 * @cosmo: a #NcHICosmo
 * @ctx: a #NcWLSurfaceMassDensityLensCtx already prepared via #nc_wl_surface_mass_density_lens_ctx_prep()
 * @zs: source redshift $z_\mathrm{source}$
 * @crit_cache: (out): a #NcWLSurfaceMassDensityCritCache
 *
 * Variant of #nc_wl_surface_mass_density_reduced_shear_crit_cache_prep() that
 * reuses a pre-computed lens-side context. Only the source-side comoving
 * distance and curvature transform are computed per call.
 *
 */
void
nc_wl_surface_mass_density_reduced_shear_crit_cache_prep_with_lens_ctx (NcWLSurfaceMassDensity *smd, NcHICosmo *cosmo, const NcWLSurfaceMassDensityLensCtx *ctx, const gdouble zs, NcWLSurfaceMassDensityCritCache *crit_cache)
{
  if (zs < ctx->zl)
  {
    crit_cache->is_source  = FALSE;
    crit_cache->sigma_crit = GSL_NAN;
    return;
  }

  {
    const gdouble x_i = 1.0 + zs;
    const gdouble ds  = nc_distance_comoving (smd->dist, cosmo, zs);
    gdouble Ds, Dls;

    switch (ctx->k)
    {
      case -1:
        Ds  = sinh (ctx->sqrt_Omega_k0 * ds) / (x_i * ctx->sqrt_Omega_k0);
        Dls = sinh (ctx->sqrt_Omega_k0 * (ds - ctx->dl)) / (x_i * ctx->sqrt_Omega_k0);
        break;
      case 0:
        Ds  = ds / x_i;
        Dls = (ds - ctx->dl) / x_i;
        break;
      case 1:
        Ds  = sin (ctx->sqrt_Omega_k0 * ds) / (x_i * ctx->sqrt_Omega_k0);
        Dls = sin (ctx->sqrt_Omega_k0 * (ds - ctx->dl)) / (x_i * ctx->sqrt_Omega_k0);
        break;
      default:
        g_assert_not_reached ();
        break;
    }

    crit_cache->is_source  = TRUE;
    crit_cache->sigma_crit = ctx->sc_over_Dl * Ds / Dls;
  }
}

/**
 * nc_wl_surface_mass_density_reduced_shear_sigma_cache_prep: (skip)
 * @dp: a #NcHaloDensityProfile
 * @cosmo: a #NcHICosmo
 * @R: projected radius with respect to the center of the lens / halo
 * @zl: lens redshift $z_\mathrm{lens}$
 * @zc: cluster redshift $z_\mathrm{cluster}$
 * @sigma_cache: (out): a #NcWLSurfaceMassDensitySigmaCache
 *
 * Computes and stores the surface mass density $\Sigma(R)$ and the mean surface mass
 * density $\overline{\Sigma}(<R)$, see Eqs. $\eqref{eq:sigma}$ and
 * $\eqref{eq:sigma_mean}$.
 *
 */
void
nc_wl_surface_mass_density_reduced_shear_sigma_cache_prep (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble zl, const gdouble zc, NcWLSurfaceMassDensitySigmaCache *sigma_cache)
{
  NCM_UNUSED (zl);

  gdouble r_s, rho_s;

  nc_halo_density_profile_r_s_rho_s (dp, cosmo, zc, &r_s, &rho_s);

  const gdouble X          = R / r_s;
  const gdouble mean_sigma = (2.0 * nc_halo_density_profile_eval_dl_cyl_mass (dp, X) / (X * X)) * rho_s * r_s;
  const gdouble sigma      = (nc_halo_density_profile_eval_dl_2d_density (dp, X)) * rho_s * r_s;

  sigma_cache->sigma      = sigma;
  sigma_cache->mean_sigma = mean_sigma;
}

/**
 * nc_wl_surface_mass_density_reduced_shear_cache:
 * @crit_cache: a #NcWLSurfaceMassDensityCritCache
 * @sigma_cache: a #NcWLSurfaceMassDensitySigmaCache
 *
 * Computes the reduced shear:
 * $$ g(R) = \frac{\gamma(R)}{1 - \kappa(R)},$$
 * where $\gamma(R)$ is the shear [nc_wl_surface_mass_density_shear()] and $\kappa(R)$
 * is the convergence [nc_wl_surface_mass_density_convergence()], using the
 * pre-computed critical surface density and surface mass density.
 *
 * Returns: $g(R)$
 *
 */
gdouble
nc_wl_surface_mass_density_reduced_shear_cache (NcWLSurfaceMassDensityCritCache *crit_cache, NcWLSurfaceMassDensitySigmaCache *sigma_cache)
{
  if (!crit_cache->is_source)
  {
    return 0.0;
  }
  else
  {
    const gdouble mu            = sigma_cache->sigma / crit_cache->sigma_crit;
    const gdouble shear         = (sigma_cache->mean_sigma - sigma_cache->sigma) / crit_cache->sigma_crit;
    const gdouble reduced_shear = shear / (1.0 - mu);

    return reduced_shear;
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
 * See [Applegate (2014)](https://arxiv.org/abs/1208.0605)
 *
 * Returns: $g(R)$, source at $z = \infty$
 */
gdouble
nc_wl_surface_mass_density_reduced_shear_infinity (NcWLSurfaceMassDensity *smd, NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble zs, const gdouble zl, const gdouble zc)
{
  if (zs < zc)
  {
    return 0.0;
  }
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
  {
    return 1.0;
  }
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

  fin = fin / r_s;

  {
    GArray *res                 = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), R->len * zs->len);
    const gdouble a             = ncm_c_c2 () / (4.0 * M_PI * ncm_c_G_mass_solar ()) * ncm_c_Mpc () / nc_hicosmo_RH_Mpc (cosmo); /* [ M_solar / Mpc^2 ] */
    const gdouble Omega_k0      = nc_hicosmo_Omega_k0 (cosmo);
    const gdouble sqrt_Omega_k0 = sqrt (fabs (Omega_k0));
    const gint k                = fabs (Omega_k0) < NCM_ZERO_LIMIT ? 0 : (Omega_k0 > 0.0 ? -1 : 1);
    guint i, j;

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
            const gdouble zs_j = g_array_index (zs, gdouble, j);

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
            {
              g_array_index (res, gdouble, pad + j) = 0.0;
            }
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
            const gdouble zs_j = g_array_index (zs, gdouble, j);

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
            {
              g_array_index (res, gdouble, pad + j) = 0.0;
            }
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
            const gdouble zs_j = g_array_index (zs, gdouble, j);

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
            {
              g_array_index (res, gdouble, pad + j) = 0.0;
            }
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

  fin = fin / r_s;

  {
    GArray *res                 = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), R->len);
    const gdouble a             = ncm_c_c2 () / (4.0 * M_PI * ncm_c_G_mass_solar ()) * ncm_c_Mpc () / nc_hicosmo_RH_Mpc (cosmo); /* [ M_solar / Mpc^2 ] */
    const gdouble Omega_k0      = nc_hicosmo_Omega_k0 (cosmo);
    const gdouble sqrt_Omega_k0 = sqrt (fabs (Omega_k0));
    const gint k                = fabs (Omega_k0) < NCM_ZERO_LIMIT ? 0 : (Omega_k0 > 0.0 ? -1 : 1);
    guint i;

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
          {
            g_array_index (res, gdouble, i) = 0.0;
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
          {
            g_array_index (res, gdouble, i) = 0.0;
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
          {
            g_array_index (res, gdouble, i) = 0.0;
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

