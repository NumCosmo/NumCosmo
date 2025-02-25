/***************************************************************************
 *            nc_halo_density_profile.c
 *
 *  Sat June 07 19:46:31 2014
 *  Copyright  2014
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_halo_density_profile.c
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
 * NcHaloDensityProfile:
 *
 * Abstract class for density profile functions.
 *
 * This abstract class describes the radial matter density profile in real space. Each
 * implementation must provide, at least, the dimensionless 3D density:
 * \begin{equation}\label{def:dlrho}
 * \hat\rho(x) \equiv \frac{\rho(x r_s)}{\rho_s}, \quad \rho(r) = \rho_s \hat\rho\left(\frac{r}{r_s}\right),
 * \end{equation}
 * where $\rho(r)$ is the actual density profile, $\rho_s$ is the profile scale and
 * $r_s$ the scale radius. This function corresponds to the virtual function
 * nc_halo_density_profile_eval_dl_density().
 *
 * # Parametrization
 *
 * The two parameters $\rho_s$ and $r_s$ are described by the fundamental
 * parametrization in terms of $M_\Delta$ and the concentration $c_\Delta$ given a mass
 * defined by $\Delta$ and a background density $\rho_\mathrm{bg}$. These quantities are
 * provided by the halo mass summary object, which is passed to the profile object at
 * construction time.
 *
 * ## Computing $r_s$
 *
 * The mass-radius relation defined in terms of the background density is
 * \begin{equation}\label{eq:mrr}
 * M_\Delta = \frac{4\pi}{3}r_\Delta^3\Delta\,\rho_\mathrm{bg},
 * \end{equation}
 * which implicitly defines $r_\Delta$. The concentration $c_\Delta$ is then
 * defined by
 * \begin{equation}\label{def:cDelta}
 * c_\Delta \equiv \frac{r_\Delta}{r_s}.
 * \end{equation}
 * Consequently, the scale radius $r_s$ can be computed from $M_\Delta$ and $c_\Delta$
 * using
 * \begin{equation}\label{def:r_s}
 * r_s = \frac{1}{c_\Delta}\left(\frac{3M_\Delta}{4\pi\Delta\,\rho_\mathrm{bg}}\right)^{1/3}
 * = \frac{r_{s0}}{(\Delta\,\rho_\mathrm{bg})^{1/3}}, \qquad
 * r_{s0} \equiv \frac{1}{c_\Delta}\left(\frac{3M_\Delta}{4\pi}\right)^{1/3}.
 * \end{equation}
 * We split the expression of $r_s$ in a constant part $r_{s0}$ and a redshift
 * dependent (time-depedent) part $(\Delta\,\rho_\mathrm{bg}(z))^{-1/3}$.
 *
 * Note that, the parameter $r_s$ can be computed directly from $(M_\Delta,\; c_\Delta)$,
 * given the mass definition, without refering to $\hat\rho(x)$.
 *
 * ## Computing $\rho_s$
 *
 * Now, applying the mass definition $M_\Delta$ in terms of the radius $r_\Delta$
 * to our profile results in
 * \begin{equation}\label{eq:def:Mr}
 * M_\Delta = \int_0^{r_\Delta}4\pi r^2\rho(r)\mathrm{d}r
 *          = 4\pi r_s^3 \rho_s \int_0^{c_\Delta}x^2\hat\rho(x)\mathrm{d}x
 *          = 4\pi r_s^3 \rho_s I_{x^2\hat\rho}(c_\Delta),
 * \end{equation}
 * where we defined
 * \begin{equation}\label{def:Ix2_dld}
 * I_{x^2\hat\rho}(c_\Delta) \equiv \int_0^{c_\Delta}x^2\hat\rho(x)\mathrm{d}x.
 * \end{equation}
 * This integral can be implemented through the virtual method
 * nc_halo_density_profile_eval_dl_spher_mass(), otherwise it will
 * be computed numerically using nc_halo_density_profile_eval_dl_density().
 * This same mass can be obtained from the background density using
 * mass-radius relation \eqref{eq:mrr}, consequently
 * \begin{equation}\label{def:rho_s}
 * \rho_s = \frac{c_\Delta^3\Delta\,\rho_\mathrm{bg}}{3I_{x^2\hat\rho}(c_\Delta)}.
 * \end{equation}
 * The only redshift dependency (time-dependency) here comes from the value of
 * $\rho_\mathrm{bg}(z)$, for this reason it is convenient to define a constant
 * quantity
 * \begin{equation}\label{def:rho_s0}
 * \rho_{s0} \equiv \frac{\rho_s}{\Delta\,\rho_\mathrm{bg}} = \frac{c_\Delta^3}{3I_{x^2\hat\rho}(c_\Delta)}.
 * \end{equation}
 *
 * # 2D projection
 *
 * The surface density obtained from the projection of the density
 * profile along the line-of-sight is given by
 * \begin{align}
 * \Sigma(R) &= \int_{-\infty}^\infty\rho(\sqrt{R^2+z^2})\mathrm{d}z, \\\\
 *           &= 2\rho_s\int_{0}^\infty\hat\rho(\sqrt{R^2/r_s^2 + z^2/r_s^2})\mathrm{d}z, \\\\ \label{eq:def:hatSigma}
 *           &= r_s\rho_s\hat\Sigma(R / r_s), & \hat\Sigma(X) &\equiv 2\int_{0}^\infty\hat\rho(\sqrt{X^2 + u^2})\mathrm{d}u.
 * \end{align}
 * In the equation above we obtain the 2D projection $\Sigma(R)$ in terms
 * of its dimensionless version $\hat\Sigma(X)$, where $X = R / r_s$.
 * The user can implement the method nc_halo_density_profile_eval_dl_2d_density()
 * providing $\hat\Sigma(X)$ directly or rely on the numerical implementation.
 *
 * ## Mass on the cylinder of radius $R$
 *
 * Using the 2D projection $\Sigma(R)$ one computes the total mass
 * inside an infinite cylinder of radius $R$ using
 * \begin{align}
 * \overline{M}(R) &= \int_0^R\Sigma(R')2\pi R'\mathrm{d}R' = 2\pi r_s^3\rho_s \hat{\overline{M}}(<R/r_s), \\\\ \label{eq:def:cylmass}
 * \hat{\overline{M}}(X) &\equiv \int_0^X\hat\Sigma(X')X'\mathrm{d}X'.
 * \end{align}
 * Here it is possible to implement the function $\hat{\overline{M}}(X)$
 * through the method nc_halo_density_profile_eval_dl_cyl_mass() or to use
 * the default numerical implementation.
 *
 * # Numerical computation
 *
 * If the implementation (i.e., a particular radial profile implementation of this abstract class) does not provide any of the functions:
 * nc_halo_density_profile_eval_dl_spher_mass(),
 * nc_halo_density_profile_eval_dl_2d_density(),
 * nc_halo_density_profile_eval_dl_cyl_mass(),
 * they will be computed numerically integrating the density
 * $\hat{\rho}$ (nc_halo_density_profile_eval_dl_density()).
 * These functions will be prepared to be computed inside the
 * interval $(X_i,\,X_f)$ defined by #NcHaloDensityProfile:lnXi
 * and #NcHaloDensityProfile:lnXf and using the relative
 * tolerance #NcHaloDensityProfile:reltol. See the following
 * functions to control this behavior:
 *  nc_halo_density_profile_set_reltol(),
 *
 *  nc_halo_density_profile_set_lnXi(),
 *
 *  nc_halo_density_profile_set_lnXf().
 *
 * # Units
 *
 * Distance: $[r] = \mathrm{Mpc}$; Mass: $[M_\Delta] = \mathrm{M}_\odot$; Density:
 * $[\rho] = \mathrm{M}_\odot \, \mathrm{Mpc}^{-3}$; Surface mass density:
 * $[\Sigma] = \mathrm{M}_\odot \, \mathrm{Mpc}^{-2}$.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_enum_types.h"
#include "lss/nc_halo_density_profile.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include "math/ncm_integrate.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_ode_spline.h"
#include "math/ncm_memory_pool.h"
#include "lss/nc_halo_mass_summary.h"

typedef struct _NcHaloDensityProfilePrivate
{
  NcHaloMassSummary *hms;
  NcmModelCtrl *hms_ctrl;
  gdouble z;
  gdouble reltol;
  gdouble lnXi;
  gdouble lnXf;
  NcmOdeSpline *dl_spher_mass_ode_s;
  NcmSpline *dl_spher_mass_s;
  NcmSpline *dl_2d_density_s;
  NcmSpline *dl_cyl_mass_s;
  gdouble rho_s0;
  gdouble r_s0;
} NcHaloDensityProfilePrivate;

enum
{
  PROP_0,
  PROP_RELTOL,
  PROP_LNXI,
  PROP_LNXF,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHaloDensityProfile, nc_halo_density_profile, NCM_TYPE_MODEL)

#define VECTOR       (NCM_MODEL (dp))
/*#define LOG10M_DELTA (ncm_model_orig_param_get (VECTOR, NC_HALO_DENSITY_PROFILE_LOG10M_DELTA)) */
/*#define C_DELTA      (ncm_model_orig_param_get (VECTOR, NC_HALO_DENSITY_PROFILE_C_DELTA)) */

static gdouble _nc_halo_density_profile_prepare_dl_spher_mass_int (gdouble y, gdouble x, gpointer userdata);

static void
nc_halo_density_profile_init (NcHaloDensityProfile *dp)
{
  NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);
  NcmSpline *s                             = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());

  self->hms                 = NULL;
  self->hms_ctrl            = ncm_model_ctrl_new (NULL);
  self->z                   = 0.0;
  self->reltol              = 0.0;
  self->lnXi                = 0.0;
  self->lnXf                = 0.0;
  self->dl_spher_mass_ode_s = NCM_ODE_SPLINE (ncm_ode_spline_new_full (s, &_nc_halo_density_profile_prepare_dl_spher_mass_int, 0.0, 0.0, 1.0));
  self->dl_spher_mass_s     = NULL;
  self->dl_2d_density_s     = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  self->dl_cyl_mass_s       = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  self->rho_s0              = 0.0;
  self->r_s0                = 0.0;

  ncm_spline_free (s);
}

static void
_nc_halo_density_profile_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHaloDensityProfile *dp = NC_HALO_DENSITY_PROFILE (object);

  g_return_if_fail (NC_IS_HALO_DENSITY_PROFILE (object));

  switch (prop_id)
  {
    case PROP_RELTOL:
      nc_halo_density_profile_set_reltol (dp, g_value_get_double (value));
      break;
    case PROP_LNXI:
      nc_halo_density_profile_set_lnXi (dp, g_value_get_double (value));
      break;
    case PROP_LNXF:
      nc_halo_density_profile_set_lnXf (dp, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_halo_density_profile_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHaloDensityProfile *dp = NC_HALO_DENSITY_PROFILE (object);

  g_return_if_fail (NC_IS_HALO_DENSITY_PROFILE (object));

  switch (prop_id)
  {
    case PROP_RELTOL:
      g_value_set_double (value, nc_halo_density_profile_get_reltol (dp));
      break;
    case PROP_LNXI:
      g_value_set_double (value, nc_halo_density_profile_get_lnXi (dp));
      break;
    case PROP_LNXF:
      g_value_set_double (value, nc_halo_density_profile_get_lnXf (dp));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_halo_density_profile_dispose (GObject *object)
{
  NcHaloDensityProfile *dp                 = NC_HALO_DENSITY_PROFILE (object);
  NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);

  ncm_ode_spline_clear (&self->dl_spher_mass_ode_s);
  self->dl_spher_mass_s = NULL;
  ncm_spline_clear (&self->dl_2d_density_s);
  ncm_spline_clear (&self->dl_cyl_mass_s);
  ncm_model_ctrl_clear (&self->hms_ctrl);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_density_profile_parent_class)->dispose (object);
}

static void
_nc_halo_density_profile_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_density_profile_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_halo_density_profile, NC_TYPE_HALO_DENSITY_PROFILE);

static void _nc_halo_density_profile_add_submodel (NcmModel *model, NcmModel *submodel);

static gdouble _nc_halo_density_profile_eval_dl_density (NcHaloDensityProfile *dp, const gdouble x);

static void _nc_halo_density_profile_prepare_ctes (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z);
static void _nc_halo_density_profile_prepare_dl_spher_mass (NcHaloDensityProfile *dp);
static void _nc_halo_density_profile_prepare_dl_2d_density (NcHaloDensityProfile *dp);
static void _nc_halo_density_profile_prepare_dl_cyl_mass (NcHaloDensityProfile *dp);

static void
nc_halo_density_profile_class_init (NcHaloDensityProfileClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_halo_density_profile_set_property;
  model_class->get_property = &_nc_halo_density_profile_get_property;
  object_class->dispose     = &_nc_halo_density_profile_dispose;
  object_class->finalize    = &_nc_halo_density_profile_finalize;

  ncm_mset_model_register_id (model_class,
                              "NcHaloDensityProfile",
                              "NcHaloDensityProfile.",
                              NULL,
                              FALSE,
                              NCM_MSET_MODEL_MAIN);

  ncm_model_class_set_name_nick (model_class, "Matter Density Profile", "DensityProfile");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  /**
   * NcHaloDensityProfile:reltol:
   *
   * Relative tolerance used in the numerical computations.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance",
                                                        GSL_DBL_EPSILON, 1.0, 1.0e-7,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcHaloDensityProfile:lnXi:
   *
   * Logarithm of the lower limit of the interval where
   * the projected densities are computed $\ln(X_i)$.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_LNXI,
                                   g_param_spec_double ("lnXi",
                                                        NULL,
                                                        "Computation interval lower limit",
                                                        -G_MAXDOUBLE, +G_MAXDOUBLE, log (1.0e-4),
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcHaloDensityProfile:lnXf:
   *
   * Logarithm of the upper limit of the interval where
   * the projected densities are computed $\ln(X_f)$.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_LNXF,
                                   g_param_spec_double ("lnXf",
                                                        NULL,
                                                        "Computation interval upper limit",
                                                        -G_MAXDOUBLE, +G_MAXDOUBLE, log (1.0e+4),
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  model_class->add_submodel = &_nc_halo_density_profile_add_submodel;

  klass->eval_dl_density    = &_nc_halo_density_profile_eval_dl_density;
  klass->eval_dl_spher_mass = &nc_halo_density_profile_eval_numint_dl_spher_mass;
  klass->eval_dl_2d_density = &nc_halo_density_profile_eval_numint_dl_2d_density;
  klass->eval_dl_cyl_mass   = &nc_halo_density_profile_eval_numint_dl_cyl_mass;
}

static void
_nc_halo_density_profile_add_submodel (NcmModel *model, NcmModel *submodel)
{
  /* Chain up : start */
  NCM_MODEL_CLASS (nc_halo_density_profile_parent_class)->add_submodel (model, submodel);
  {
    NcHaloDensityProfile *dp                 = NC_HALO_DENSITY_PROFILE (model);
    NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);

    if (ncm_model_id (submodel) == nc_halo_mass_summary_id ())
    {
      if (self->hms != NULL)
        g_error ("Halo mass summary was already defined in `%s' class.", G_OBJECT_CLASS_NAME (dp));

      self->hms = NC_HALO_MASS_SUMMARY (submodel);
    }
  }
}

static gdouble
_nc_halo_density_profile_eval_dl_density (NcHaloDensityProfile *dp, const gdouble x)
{
  g_error ("Required method eval_dl_density not implemented by `%s' class.", G_OBJECT_CLASS_NAME (dp));

  return 0.0;
}

enum
{
  PREPARE_CTES = 0,
  PREPARE_DL_SPHER_MASS,
  PREPARE_DL_2D_DENSITY,
  PREPARE_DL_CYL_MASS,
};

static void
_nc_halo_density_profile_hms_is_update (NcHaloDensityProfile *dp)
{
  NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);
  gboolean hms_updated                     = ncm_model_ctrl_update (self->hms_ctrl, NCM_MODEL (self->hms));

  if (hms_updated)
    ncm_model_state_mark_outdated (NCM_MODEL (dp));
}

static void
_nc_halo_density_profile_prepare_ctes (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z)
{
  _nc_halo_density_profile_hms_is_update (dp);

  if (!ncm_model_lstate_is_update (NCM_MODEL (dp), PREPARE_CTES))
  {
    NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);
    const gdouble cDelta                     = nc_halo_mass_summary_concentration (self->hms, cosmo, z);
    const gdouble MDelta                     = nc_halo_mass_summary_mass (self->hms);

    self->r_s0   = cbrt (3.0 * MDelta / (4.0 * ncm_c_pi ())) / cDelta;
    self->rho_s0 = gsl_pow_3 (cDelta) / (3.0 * nc_halo_density_profile_eval_dl_spher_mass (dp, cDelta));

    ncm_model_lstate_set_update (NCM_MODEL (dp), PREPARE_CTES);
  }
}

static gdouble
_nc_halo_density_profile_prepare_dl_spher_mass_int (gdouble y, gdouble x, gpointer userdata)
{
  NcHaloDensityProfile *dp = NC_HALO_DENSITY_PROFILE (userdata);

  if (x == 0.0)
    return 0.0;

  return x * x * nc_halo_density_profile_eval_dl_density (dp, x);
}

static void
_nc_halo_density_profile_prepare_dl_spher_mass (NcHaloDensityProfile *dp)
{
  _nc_halo_density_profile_hms_is_update (dp);

  if (!ncm_model_lstate_is_update (NCM_MODEL (dp), PREPARE_DL_SPHER_MASS))
  {
    NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);

    ncm_ode_spline_set_xf (self->dl_spher_mass_ode_s, exp (self->lnXf));
    ncm_ode_spline_set_reltol (self->dl_spher_mass_ode_s, self->reltol * 1.0e-1);
    ncm_ode_spline_set_abstol (self->dl_spher_mass_ode_s, 1.0e-30);
    ncm_ode_spline_prepare (self->dl_spher_mass_ode_s, dp);
    self->dl_spher_mass_s = ncm_ode_spline_peek_spline (self->dl_spher_mass_ode_s);

    ncm_model_lstate_set_update (NCM_MODEL (dp), PREPARE_DL_SPHER_MASS);
  }
}

typedef struct _NcHaloDensityProfile2D
{
  NcHaloDensityProfile *dp;
  gdouble X;
  gsl_function *F;
  gsl_function *F2;
  gsl_integration_workspace *w;
} NcHaloDensityProfile2D;

static gdouble
_nc_halo_density_profile_prepare_dl_2d_density_X_u (gdouble lnu, gpointer userdata)
{
  NcHaloDensityProfile2D *dp2D = (NcHaloDensityProfile2D *) userdata;
  const gdouble u              = exp (lnu);

  return nc_halo_density_profile_eval_dl_density (dp2D->dp, hypot (dp2D->X, u)) * u;
}

static gdouble
_nc_halo_density_profile_prepare_dl_2d_density_X (gdouble lnX, gpointer userdata)
{
  NcHaloDensityProfile2D *dp2D             = (NcHaloDensityProfile2D *) userdata;
  NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp2D->dp);
  gdouble err, dl_2d_density_dX, dl_2d_density_X, lnx0, lnx1;
  const gdouble dlnx = 2.0 * M_LN10;
  gdouble abstol     = 0.0;
  gint key           = 6;

  dp2D->X = exp (lnX);

  dl_2d_density_X = 0.0;
  lnx0            = lnX + GSL_LOG_DBL_EPSILON;
  lnx1            = lnx0 + dlnx;

  do {
    gsl_integration_qag (dp2D->F, lnx0, lnx1, abstol, self->reltol, NCM_INTEGRAL_PARTITION, key, dp2D->w, &dl_2d_density_dX, &err);

    dl_2d_density_X += dl_2d_density_dX;
    abstol           = dl_2d_density_X * self->reltol;

    lnx0 = lnx1;
    lnx1 = lnx0 + dlnx;
  } while (fabs (dl_2d_density_dX / dl_2d_density_X) > self->reltol);

  return log (2.0 * dl_2d_density_X);
}

static void
_nc_halo_density_profile_prepare_dl_2d_density (NcHaloDensityProfile *dp)
{
  _nc_halo_density_profile_hms_is_update (dp);

  if (!ncm_model_lstate_is_update (NCM_MODEL (dp), PREPARE_DL_2D_DENSITY))
  {
    NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);
    gsl_integration_workspace **w            = ncm_integral_get_workspace ();
    NcHaloDensityProfile2D dp2d;
    gsl_function F1, F2;

    dp2d.dp = dp;
    dp2d.X  = 0.0;
    dp2d.F  = &F1;
    dp2d.w  = *w;

    F1.function = &_nc_halo_density_profile_prepare_dl_2d_density_X_u;
    F1.params   = &dp2d;

    F2.function = &_nc_halo_density_profile_prepare_dl_2d_density_X;
    F2.params   = &dp2d;

    ncm_spline_set_func (self->dl_2d_density_s, NCM_SPLINE_FUNCTION_SPLINE, &F2, self->lnXi, self->lnXf, 0, self->reltol);

    ncm_memory_pool_return (w);
    ncm_model_lstate_set_update (NCM_MODEL (dp), PREPARE_DL_2D_DENSITY);
  }
}

static gdouble
_nc_halo_density_profile_prepare_dl_cyl_mass_X_x_1 (gdouble lnx, gpointer userdata)
{
  NcHaloDensityProfile2D *dp2D = (NcHaloDensityProfile2D *) userdata;
  const gdouble x              = exp (lnx);
  const gdouble x2             = x * x;

  return x2 * x * nc_halo_density_profile_eval_dl_density (dp2D->dp, x);
}

static gdouble
_nc_halo_density_profile_prepare_dl_cyl_mass_X_x_2 (gdouble lnmu, gpointer userdata)
{
  NcHaloDensityProfile2D *dp2D = (NcHaloDensityProfile2D *) userdata;
  const gdouble fa             = sqrt (fabs (expm1 (-2.0 * lnmu)));
  const gdouble mu             = exp (lnmu);

  return nc_halo_density_profile_eval_dl_density (dp2D->dp, dp2D->X * mu) * mu / (1.0 + fa);
}

static gdouble
_nc_halo_density_profile_prepare_dl_cyl_mass_X (gdouble lnX, gpointer userdata)
{
  NcHaloDensityProfile2D *dp2D             = (NcHaloDensityProfile2D *) userdata;
  NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp2D->dp);

  gdouble dl_cyl_mass_X = 0.0;
  gdouble err, dl_cyl_mass_X_i, X3;
  gdouble abstol = 0.0;
  gint key       = 6;

  dp2D->X = exp (lnX);
  X3      = gsl_pow_3 (dp2D->X);

  gsl_integration_qag (dp2D->F, 2.0 * GSL_LOG_DBL_EPSILON, lnX, abstol, self->reltol, NCM_INTEGRAL_PARTITION, key, dp2D->w, &dl_cyl_mass_X_i, &err);
  dl_cyl_mass_X += dl_cyl_mass_X_i;

  gsl_integration_qag (dp2D->F2, 0.0, 1.0, abstol, self->reltol, NCM_INTEGRAL_PARTITION, key, dp2D->w, &dl_cyl_mass_X_i, &err);
  dl_cyl_mass_X += X3 * dl_cyl_mass_X_i;

  abstol += dl_cyl_mass_X_i * self->reltol;

  gsl_integration_qag (dp2D->F2, 1.0, 10.0, abstol, self->reltol, NCM_INTEGRAL_PARTITION, key, dp2D->w, &dl_cyl_mass_X_i, &err);
  dl_cyl_mass_X += X3 * dl_cyl_mass_X_i;

  abstol += dl_cyl_mass_X_i * self->reltol;

  gsl_integration_qag (dp2D->F2, 10.0, GSL_LOG_DBL_MAX / 4.0, abstol, self->reltol, NCM_INTEGRAL_PARTITION, key, dp2D->w, &dl_cyl_mass_X_i, &err);
  dl_cyl_mass_X += X3 * dl_cyl_mass_X_i;

/*printf ("% 22.15e\n", dl_cyl_mass_X_i); */

  return log (2.0 * dl_cyl_mass_X);
}

static void
_nc_halo_density_profile_prepare_dl_cyl_mass (NcHaloDensityProfile *dp)
{
  _nc_halo_density_profile_hms_is_update (dp);

  if (!ncm_model_lstate_is_update (NCM_MODEL (dp), PREPARE_DL_CYL_MASS))
  {
    NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);
    gsl_integration_workspace **w            = ncm_integral_get_workspace ();
    NcHaloDensityProfile2D dp2d;
    gsl_function F1, F12, F2;

    dp2d.dp = dp;
    dp2d.X  = 0.0;
    dp2d.F  = &F1;
    dp2d.F2 = &F12;
    dp2d.w  = *w;

    F1.function = &_nc_halo_density_profile_prepare_dl_cyl_mass_X_x_1;
    F1.params   = &dp2d;

    F12.function = &_nc_halo_density_profile_prepare_dl_cyl_mass_X_x_2;
    F12.params   = &dp2d;

    F2.function = &_nc_halo_density_profile_prepare_dl_cyl_mass_X;
    F2.params   = &dp2d;

    ncm_spline_set_func (self->dl_cyl_mass_s, NCM_SPLINE_FUNCTION_SPLINE, &F2, self->lnXi, self->lnXf, 0, self->reltol);

    ncm_memory_pool_return (w);
    ncm_model_lstate_set_update (NCM_MODEL (dp), PREPARE_DL_CYL_MASS);
  }
}

/**
 * nc_halo_density_profile_ref:
 * @dp: a #NcHaloDensityProfile
 *
 * Increases the reference count of @dp by one.
 *
 * Returns: (transfer full): @dp
 */
NcHaloDensityProfile *
nc_halo_density_profile_ref (NcHaloDensityProfile *dp)
{
  return g_object_ref (dp);
}

/**
 * nc_halo_density_profile_free:
 * @dp: a #NcHaloDensityProfile
 *
 * Atomically decrements the reference count of @dp by one. If the reference count drops to 0,
 * all memory allocated by @dp is released.
 *
 */
void
nc_halo_density_profile_free (NcHaloDensityProfile *dp)
{
  g_object_unref (dp);
}

/**
 * nc_halo_density_profile_clear:
 * @dp: a #NcHaloDensityProfile
 *
 * Atomically decrements the reference count of @dp by one. If the reference count drops to 0,
 * all memory allocated by @dp is released. Set the pointer to NULL;
 *
 */
void
nc_halo_density_profile_clear (NcHaloDensityProfile **dp)
{
  g_clear_object (dp);
}

/**
 * nc_halo_density_profile_set_reltol:
 * @dp: a #NcHaloDensityProfile
 * @reltol: relative tolerance
 *
 * Sets the relative tolerance used in the numerical computations.
 *
 */
void
nc_halo_density_profile_set_reltol (NcHaloDensityProfile *dp, const gdouble reltol)
{
  NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);

  g_assert_cmpfloat (reltol, >, GSL_DBL_EPSILON);
  g_assert_cmpfloat (reltol, <, 1.0);

  self->reltol = reltol;
  ncm_model_state_mark_outdated (NCM_MODEL (dp));
}

/**
 * nc_halo_density_profile_set_lnXi:
 * @dp: a #NcHaloDensityProfile
 * @lnXi: interval lower limit $\ln(X_i)$
 *
 * Sets the numerical computation lower limit.
 *
 */
void
nc_halo_density_profile_set_lnXi (NcHaloDensityProfile *dp, const gdouble lnXi)
{
  NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);

  self->lnXi = lnXi;
  ncm_model_state_mark_outdated (NCM_MODEL (dp));
}

/**
 * nc_halo_density_profile_set_lnXf:
 * @dp: a #NcHaloDensityProfile
 * @lnXf: interval upper limit $\ln(X_f)$
 *
 * Sets the numerical computation upper limit.
 *
 */
void
nc_halo_density_profile_set_lnXf (NcHaloDensityProfile *dp, const gdouble lnXf)
{
  NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);

  self->lnXf = lnXf;
  ncm_model_state_mark_outdated (NCM_MODEL (dp));
}

/**
 * nc_halo_density_profile_get_reltol:
 * @dp: a #NcHaloDensityProfile
 *
 * Gets the current relative tolerance.
 *
 * Returns: reltol.
 */
gdouble
nc_halo_density_profile_get_reltol (NcHaloDensityProfile *dp)
{
  NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);

  return self->reltol;
}

/**
 * nc_halo_density_profile_get_lnXi:
 * @dp: a #NcHaloDensityProfile
 *
 * Gets the computation interval lower limit $\ln(X_i)$.
 *
 * Returns: $\ln(X_i)$.
 */
gdouble
nc_halo_density_profile_get_lnXi (NcHaloDensityProfile *dp)
{
  NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);

  return self->lnXi;
}

/**
 * nc_halo_density_profile_get_lnXf:
 * @dp: a #NcHaloDensityProfile
 *
 * Gets the computation interval upper limit $\ln(X_f)$.
 *
 * Returns: $\ln(X_f)$.
 */
gdouble
nc_halo_density_profile_get_lnXf (NcHaloDensityProfile *dp)
{
  NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);

  return self->lnXf;
}

/**
 * nc_halo_density_profile_get_phys_limts:
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 * @dp: a #NcHaloDensityProfile
 * @Ri: (out): lower limit $R_i\;\left[\mathrm{Mpc}\right]$
 * @Rf: (out): lower limit $R_f\;\left[\mathrm{Mpc}\right]$
 *
 * Gets the physical computation interval $(R_i,\, R_f)$.
 * This interval is relevant only if the object relies on
 * the numerical computation of the functions:
 * - nc_halo_density_profile_eval_dl_2d_density()
 * - nc_halo_density_profile_eval_dl_cyl_mass()
 *
 */
void
nc_halo_density_profile_get_phys_limts (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z, gdouble *Ri, gdouble *Rf)
{
  NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);
  const gdouble r_s                        = nc_halo_density_profile_r_s (dp, cosmo, z);

  Ri[0] = exp (self->lnXi) * r_s;
  Rf[0] = exp (self->lnXf) * r_s;
}

/**
 * nc_halo_density_profile_peek_mass_summary:
 * @dp: a #NcHaloDensityProfile
 *
 * Peeks the #NcHaloMassSummary submodel.
 *
 * Returns: (transfer none): the #NcHaloMassSummary submodel.
 */
NcHaloMassSummary *
nc_halo_density_profile_peek_mass_summary (NcHaloDensityProfile *dp)
{
  NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);

  return self->hms;
}

/**
 * nc_halo_density_profile_eval_dl_density: (virtual eval_dl_density)
 * @dp: a #NcHaloDensityProfile
 * @x: dimensionless radius $x = r / r_s$
 *
 * This function computes the dimensionless density profile,
 * see Eq. \eqref{def:dlrho}.
 *
 * Returns: the value of the dimensionless density profile $\hat\rho(x)$.
 */
gdouble
nc_halo_density_profile_eval_dl_density (NcHaloDensityProfile *dp, const gdouble x)
{
  return NC_HALO_DENSITY_PROFILE_GET_CLASS (dp)->eval_dl_density (dp, x);
}

/**
 * nc_halo_density_profile_eval_dl_spher_mass: (virtual eval_dl_spher_mass)
 * @dp: a #NcHaloDensityProfile
 * @x: dimensionless radius $x = r / r_s$
 *
 * This function computes the 2d projection of the dimensionless density
 * profile as described in Eq. \eqref{def:Ix2_dld}.
 *
 * Returns: the value of the integral $I_{x^2\hat\rho}(x)$.
 */
gdouble
nc_halo_density_profile_eval_dl_spher_mass (NcHaloDensityProfile *dp, const gdouble x)
{
  /*NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp); */

  return NC_HALO_DENSITY_PROFILE_GET_CLASS (dp)->eval_dl_spher_mass (dp, x);
}

/**
 * nc_halo_density_profile_eval_dl_2d_density: (virtual eval_dl_2d_density)
 * @dp: a #NcHaloDensityProfile
 * @X: dimensionless 2D radius $X = R / r_s$
 *
 * This function computes the dimensionless 2D density profile,
 * see Eq. \eqref{eq:def:hatSigma}.
 *
 * Returns: the value of the dimensionless 2D density profile $\hat\Sigma(X)$.
 */
gdouble
nc_halo_density_profile_eval_dl_2d_density (NcHaloDensityProfile *dp, const gdouble X)
{
  return NC_HALO_DENSITY_PROFILE_GET_CLASS (dp)->eval_dl_2d_density (dp, X);
}

/**
 * nc_halo_density_profile_eval_dl_cyl_mass: (virtual eval_dl_cyl_mass)
 * @dp: a #NcHaloDensityProfile
 * @X: dimensionless 2D radius $X = R / r_s$
 *
 * This function computes the dimensionless cylinder mass,
 * see Eq. \eqref{eq:def:cylmass}.
 *
 * Returns: the value of the dimensionless cylinder mass $\hat{\overline{\Sigma}}(X)$.
 */
gdouble
nc_halo_density_profile_eval_dl_cyl_mass (NcHaloDensityProfile *dp, const gdouble X)
{
  return NC_HALO_DENSITY_PROFILE_GET_CLASS (dp)->eval_dl_cyl_mass (dp, X);
}

/**
 * nc_halo_density_profile_rho_s:
 * @dp: a #NcHaloDensityProfile
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * This function computes the $\rho_s$ parameter as described in
 * Eqs. \eqref{def:rho_s} and \eqref{def:rho_s0}.
 *
 * Returns: the value of $\rho_s(z)\left[M_\odot\times\mathrm{Mpc}^{-3}\right]$.
 */
gdouble
nc_halo_density_profile_rho_s (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z)
{
  NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);

  _nc_halo_density_profile_prepare_ctes (dp, cosmo, z);

  return self->rho_s0 * nc_halo_mass_summary_Delta_rho_bg (self->hms, cosmo, z);
}

/**
 * nc_halo_density_profile_r_s:
 * @dp: a #NcHaloDensityProfile
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * This function computes the $r_s$ parameter as described in
 * Eq. \eqref{def:r_s}.
 *
 * Returns: the value of $r_s(z)\;\left[\mathrm{Mpc}\right]$.
 */
gdouble
nc_halo_density_profile_r_s (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z)
{
  NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);

  _nc_halo_density_profile_prepare_ctes (dp, cosmo, z);

  return self->r_s0 / cbrt (nc_halo_mass_summary_Delta_rho_bg (self->hms, cosmo, z));
}

/**
 * nc_halo_density_profile_r_s_rho_s:
 * @dp: a #NcHaloDensityProfile
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 * @r_s: (out): $r_s\;\left[\mathrm{Mpc}\right]$
 * @rho_s: (out): $\rho_s\;\left[M_\odot\times\mathrm{Mpc}^{-3}\right]$
 *
 * This function computes $r_s$ and $\rho_s$ parameters as described in
 * Eqs. \eqref{def:r_s}, \eqref{def:rho_s} and \eqref{def:rho_s0}.
 *
 */
void
nc_halo_density_profile_r_s_rho_s (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z, gdouble *r_s, gdouble *rho_s)
{
  NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);
  const gdouble Delta_rho_bg               = nc_halo_mass_summary_Delta_rho_bg (self->hms, cosmo, z);

  _nc_halo_density_profile_prepare_ctes (dp, cosmo, z);

  r_s[0]   = self->r_s0 / cbrt (Delta_rho_bg);
  rho_s[0] = self->rho_s0 * Delta_rho_bg;
}

/**
 * nc_halo_density_profile_eval_density:
 * @dp: a #NcHaloDensityProfile
 * @cosmo: a #NcHICosmo
 * @r: radius $r\;\left[\mathrm{Mpc}\right]$
 * @z: redshift $z$
 *
 * This function computes the density profile in real space.
 *
 * Returns: the value of the density profile $\rho(r)\;\left[M_\odot\times\mathrm{Mpc}^{-3}\right]$.
 */
gdouble
nc_halo_density_profile_eval_density (NcHaloDensityProfile *dp, NcHICosmo *cosmo, gdouble r, gdouble z)
{
  gdouble r_s, rho_s;

  nc_halo_density_profile_r_s_rho_s (dp, cosmo, z, &r_s, &rho_s);

  return rho_s * nc_halo_density_profile_eval_dl_density (dp, r / r_s);
}

/**
 * nc_halo_density_profile_eval_spher_mass:
 * @dp: a #NcHaloDensityProfile
 * @cosmo: a #NcHICosmo
 * @r: radius in Mpc
 * @z: redshift $z$
 *
 * This function computes the total mass enclosed in the
 * sphere of radius $r$, see Eq. \eqref{eq:def:Mr}.
 *
 * Returns: the total spherical mass $M(< r)\;\left[M_\odot\right]$.
 */
gdouble
nc_halo_density_profile_eval_spher_mass (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble r, const gdouble z)
{
  /*NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp); */
  gdouble r_s, rho_s, sVol, x;

  nc_halo_density_profile_r_s_rho_s (dp, cosmo, z, &r_s, &rho_s);

  sVol = 4.0 * ncm_c_pi () * gsl_pow_3 (r_s) * rho_s;
  x    = r / r_s;

  return sVol * nc_halo_density_profile_eval_dl_spher_mass (dp, x);
}

/**
 * nc_halo_density_profile_eval_spher_mass_delta:
 * @dp: a #NcHaloDensityProfile
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * This function computes the total mass enclosed in the
 * sphere of radius $r_\Delta$, see Eq. \eqref{eq:def:Mr}.
 *
 * Returns: the total spherical mass $M(< r_\Delta)\;\left[M_\odot\right]$.
 */
gdouble
nc_halo_density_profile_eval_spher_mass_delta (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble z)
{
  NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);
  gdouble concentration                    = nc_halo_mass_summary_concentration (self->hms, cosmo, z);
  gdouble r_s, rho_s, sVol;

  nc_halo_density_profile_r_s_rho_s (dp, cosmo, z, &r_s, &rho_s);

  sVol = 4.0 * ncm_c_pi () * gsl_pow_3 (r_s) * rho_s;

  return sVol * nc_halo_density_profile_eval_dl_spher_mass (dp, concentration);
}

/**
 * nc_halo_density_profile_eval_2d_density:
 * @dp: a #NcHaloDensityProfile
 * @cosmo: a #NcHICosmo
 * @R: radius $R$ in Mpc
 * @z: redshift $z$
 *
 * This function computes the 2D projection of the density profile
 * at radius $R$ and redshift $z$, see Eq. \eqref{}.
 *
 * Returns: the value of $\Sigma(R)\left[M_\odot\times\mathrm{Mpc}^{-2}\right]$.
 */
gdouble
nc_halo_density_profile_eval_2d_density (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z)
{
  gdouble r_s, rho_s;

  nc_halo_density_profile_r_s_rho_s (dp, cosmo, z, &r_s, &rho_s);

  return r_s * rho_s * nc_halo_density_profile_eval_dl_2d_density (dp, R / r_s);
}

/**
 * nc_halo_density_profile_eval_cyl_mass:
 * @dp: a #NcHaloDensityProfile
 * @cosmo: a #NcHICosmo
 * @R: radius $R\left[\mathrm{Mpc}\right]$
 * @z: redshift $z$
 *
 * This function computes the total mass enclose in the
 * cylinder of radius $R$, see Eq. \eqref{}.
 *
 * Returns: the value of $\overline{\Sigma}(R)\left[M_\odot\right]$.
 */
gdouble
nc_halo_density_profile_eval_cyl_mass (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z)
{
  gdouble r_s, rho_s, sVol;

  nc_halo_density_profile_r_s_rho_s (dp, cosmo, z, &r_s, &rho_s);

  sVol = 2.0 * ncm_c_pi () * gsl_pow_3 (r_s) * rho_s;

  return sVol * nc_halo_density_profile_eval_dl_cyl_mass (dp, R / r_s);
}

/**
 * nc_halo_density_profile_eval_density_array:
 * @dp: a #NcHaloDensityProfile
 * @cosmo: a #NcHICosmo
 * @r: (in) (element-type gdouble): radius $r\;\left[\mathrm{Mpc}\right]$
 * @fin: input array factor
 * @fout: output array factor
 * @z: redshift $z$
 *
 * This function computes the density profile in real space.
 *
 * Returns: (transfer full) (element-type gdouble): the value of the density profile $\rho(r)\;\left[M_\odot\times\mathrm{Mpc}^{-3}\right]$.
 */
GArray *
nc_halo_density_profile_eval_density_array (NcHaloDensityProfile *dp, NcHICosmo *cosmo, GArray *r, gdouble fin, gdouble fout, const gdouble z)
{
  gdouble r_s, rho_s;

  g_assert_cmpint (r->len, >, 0);

  nc_halo_density_profile_r_s_rho_s (dp, cosmo, z, &r_s, &rho_s);

  fin  = fin / r_s;
  fout = fout * rho_s;

  {
    GArray *res = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), r->len);
    guint i;

    g_array_set_size (res, r->len);

    for (i = 0; i < r->len; i++)
    {
      g_array_index (res, gdouble, i) = fout * nc_halo_density_profile_eval_dl_density (dp, g_array_index (r, gdouble, i) * fin);
    }

    return res;
  }
}

/**
 * nc_halo_density_profile_eval_2d_density_array:
 * @dp: a #NcHaloDensityProfile
 * @cosmo: a #NcHICosmo
 * @R: (in) (element-type gdouble): radius $r\;\left[\mathrm{Mpc}\right]$
 * @fin: input array factor
 * @fout: output array factor
 * @z: redshift $z$
 *
 * This function computes 2D projection of the density profile
 * at radius $R$ and redshift $z$, see Eq. \eqref{}.
 *
 * Returns: (transfer full) (element-type gdouble): the value of $\Sigma(R)\left[M_\odot\times\mathrm{Mpc}^{-2}\right]$.
 */
GArray *
nc_halo_density_profile_eval_2d_density_array (NcHaloDensityProfile *dp, NcHICosmo *cosmo, GArray *R, gdouble fin, gdouble fout, const gdouble z)
{
  gdouble r_s, rho_s;

  g_assert_cmpint (R->len, >, 0);

  nc_halo_density_profile_r_s_rho_s (dp, cosmo, z, &r_s, &rho_s);

  fin  = fin / r_s;
  fout = fout * r_s * rho_s;

  {
    GArray *res = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), R->len);
    guint i;

    g_array_set_size (res, R->len);

    for (i = 0; i < R->len; i++)
    {
      g_array_index (res, gdouble, i) = fout * nc_halo_density_profile_eval_dl_2d_density (dp, g_array_index (R, gdouble, i) * fin);
    }

    return res;
  }
}

/**
 * nc_halo_density_profile_eval_cyl_mass_array:
 * @dp: a #NcHaloDensityProfile
 * @cosmo: a #NcHICosmo
 * @R: (in) (element-type gdouble): radius $r\;\left[\mathrm{Mpc}\right]$
 * @fin: input array factor
 * @fout: output array factor
 * @z: redshift $z$
 *
 * This function computes the total mass enclose in the
 * cylinder of radius $R$, see Eq. \eqref{}.
 *
 * Returns: (transfer full) (element-type gdouble): the value of $\overline{\Sigma}(R)\left[M_\odot\right]$.
 */
GArray *
nc_halo_density_profile_eval_cyl_mass_array (NcHaloDensityProfile *dp, NcHICosmo *cosmo, GArray *R, gdouble fin, gdouble fout, const gdouble z)
{
  gdouble r_s, rho_s;

  g_assert_cmpint (R->len, >, 0);

  nc_halo_density_profile_r_s_rho_s (dp, cosmo, z, &r_s, &rho_s);

  fin  = fin / r_s;
  fout = fout * 2.0 * ncm_c_pi () * rho_s * gsl_pow_3 (r_s);

  {
    GArray *res = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), R->len);
    guint i;

    g_array_set_size (res, R->len);

    for (i = 0; i < R->len; i++)
    {
      g_array_index (res, gdouble, i) = fout * nc_halo_density_profile_eval_dl_cyl_mass (dp, g_array_index (R, gdouble, i) * fin);
    }

    return res;
  }
}

/**
 * nc_halo_density_profile_eval_numint_dl_spher_mass:
 * @dp: a #NcHaloDensityProfile
 * @X: dimensionless radius $X = r /r_s$
 *
 * This function computes the enclosed mass in a sphere within radius $X$, i.e.,
 * the 2d projection of the dimensionless density
 * profile as described in Eq. \eqref{def:Ix2_dld}. This is the default
 * implementation that will be used unless the child object provides one.
 * This interface is present for testing purpose.
 *
 * Returns: the value of the integral $I_{x^2\hat\rho}(c_\Delta)$.
 */
gdouble
nc_halo_density_profile_eval_numint_dl_spher_mass (NcHaloDensityProfile *dp, const gdouble X)
{
  NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);

  _nc_halo_density_profile_prepare_dl_spher_mass (dp);

  return ncm_spline_eval (self->dl_spher_mass_s, X);
}

/**
 * nc_halo_density_profile_eval_numint_dl_2d_density:
 * @dp: a #NcHaloDensityProfile
 * @X: dimensionless 2D radius $X = R / r_s$
 *
 * This function computes the dimensionless 2D density profile,
 * see Eq. \eqref{eq:def:hatSigma}. This is the default
 * implementation that will be used unless the child object provides one.
 * This interface is present for testing purpose.
 *
 * Returns: the value of the dimensionless 2D density profile $\hat\Sigma(X)$.
 */
gdouble
nc_halo_density_profile_eval_numint_dl_2d_density (NcHaloDensityProfile *dp, const gdouble X)
{
  NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);

  _nc_halo_density_profile_prepare_dl_2d_density (dp);

  return exp (ncm_spline_eval (self->dl_2d_density_s, log (X)));
}

/**
 * nc_halo_density_profile_eval_numint_dl_cyl_mass:
 * @dp: a #NcHaloDensityProfile
 * @X: dimensionless 2D radius $X = R / r_s$
 *
 * This function computes the dimensionless cylinder mass,
 * see Eq. \eqref{eq:def:cylmass}. This is the default
 * implementation that will be used unless the child object
 * provides one. This interface is present for testing purpose.
 *
 * Returns: the value of the dimensionless cylinder mass $\hat{\overline{\Sigma}}(X)$.
 */
gdouble
nc_halo_density_profile_eval_numint_dl_cyl_mass (NcHaloDensityProfile *dp, const gdouble X)
{
  NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);

  _nc_halo_density_profile_prepare_dl_cyl_mass (dp);

  return exp (ncm_spline_eval (self->dl_cyl_mass_s, log (X)));
}

/**
 * nc_halo_density_profile_get_numint_splines:
 * @dp: a #NcHaloDensityProfile
 * @spher_mass: (out) (transfer full): size of the dl_spher_mass_s spline
 * @twod_density: (out) (transfer full): size of the dl_2d_density_s spline
 * @cyl_mass: (out) (transfer full): size of the dl_cyl_mass_s spline
 *
 * This function provides the splines computed in this object, #NcHaloDensityProfile:
 * dl_2d_density_s and dl_cyl_mass_s.
 *
 */
void
nc_halo_density_profile_get_numint_splines (NcHaloDensityProfile *dp, NcmSpline **spher_mass, NcmSpline **twod_density, NcmSpline **cyl_mass)
{
  NcHaloDensityProfilePrivate * const self = nc_halo_density_profile_get_instance_private (dp);

  _nc_halo_density_profile_prepare_dl_spher_mass (dp);
  _nc_halo_density_profile_prepare_dl_cyl_mass (dp);
  _nc_halo_density_profile_prepare_dl_2d_density (dp);

  *spher_mass   = ncm_spline_ref (self->dl_spher_mass_s);
  *twod_density = ncm_spline_ref (self->dl_2d_density_s);
  *cyl_mass     = ncm_spline_ref (self->dl_cyl_mass_s);
}

