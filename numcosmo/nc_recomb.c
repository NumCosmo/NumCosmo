/***************************************************************************
 *            nc_recomb.c
 *
 *  Sun Oct  5 20:40:30 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_recomb
 * @title: NcRecomb
 * @short_description: Abstract class for cosmic recombination.
 * @include: numcosmo/nc_recomb.h
 *
 * $
 *  \newcommand{\He}{\text{He}}
 *  \newcommand{\HeI}{\text{HeI}}
 *  \newcommand{\HeII}{\text{HeII}}
 *  \newcommand{\HeIII}{\text{HeIII}}
 *  \newcommand{\Hy}{\text{H}}
 *  \newcommand{\HyI}{\text{HI}}
 *  \newcommand{\HyII}{\text{HII}}
 *  \newcommand{\e}{{\text{e}^-}}
 * $
 *
 * The #NcRecomb abstract object describe a general recombination process.
 * To describe the functions we have the following definitions,
 * more details in [Weinberg (2008)][XWeinberg2008].
 *
 * We refer to the total number of hydrogen nucleus (ionized or not) as
 * $n_{\Hy}$, the neutral hydrogen atoms as $n_{\HyI}$ and ionized hydrogen as
 * $n_{\HyII}$ and, therefore, $n_{\HyI} + n_{\HyII} = n_{\Hy}$. Similarly,
 * the number of helium nuclei is $n_{\He}$ and the neutral,
 * single and double ionized are $n_{\HeI}$, $n_{\HeII}$ and $n_{\HeIII}$,
 * respectively.
 *
 * We also define the helium primordial abundance as the ratio of the helium
 * mass to the total baryonic mass, i.e.,
 * \begin{align}\label{def:Yp}
 * Y_p = \frac{n_{\He} m_{\He}}{(n_{\He} m_{\He} + n_{\Hy} m_{\Hy})},
 * \end{align}
 * where $m_{\Hy}$ and $m_{\He}$ are the hydrogen and helium masses.
 *
 * The element abundances are defined as the ratio of the element by the total
 * number of free protons $n_p \equiv n_\Hy$:
 * \begin{align}
 * X_{f} = \frac{n_{f}}{n_p},
 * \end{align}
 * where $f$ is any one of the elements described above and $\e$ represent
 * the number of free electrons.
 *
 * These fractions have the following properties:
 * \begin{align}\label{eq:Hy:add}
 * X_\HyI + X_{\HyII} &= 1, \\\\ \label{eq:He:add}
 * X_{\HeI} + X_{\HeII} + X_{\HeIII} &= X_{\He}, \\\\
 * X_{\He} &\equiv \frac{m_p}{m_{\He}}\frac{Y_p}{1-Y_p}.
 * \end{align}
 * We also define the number of free electrons as $n_\e$. Assuming a neutral
 * universe we have
 * \begin{equation}\label{eq:def:Xe}
 * X_\e = X_{\HyII} + X_{\HeII} + 2X_{\HeIII}.
 * \end{equation}
 *
 * <bridgehead>Equilibrium fractions</bridgehead>
 *
 * Equilibrium ratio $X_{\HyII}X_\e / X_\Hy$ through Saha equation, i.e.,
 * \begin{equation}\label{eq:saha:HyI}
 * \frac{X_{\HyII}X_\e}{X_\HyI} = \frac{e^{-\HyI_{1s}/(k_BT)}}{n_{\Hy}\lambda_{\e}^3},
 * \end{equation}
 * where $\HyI_{1s}$ is the hydrogen $1s$ binding energy ncm_c_HI_ion_E_1s_2S0_5(),
 * $\lambda_{\e}$ is the electron thermal wavelength, i.e.,
 * \begin{equation}
 * \lambda_{\e} = \sqrt{\frac{2\pi\hbar^2}{m_\e{}k_BT}},
 * \end{equation}
 * $k_K$ is the Boltzmann constant ncm_c_kb (), $m_\e$ the electron mass and
 * $\hbar$ is the Planck constant ncm_c_hbar ().
 *
 * This calculation is done using the Saha equation as in
 * [Weinberg (2008)][XWeinberg2008].
 *
 * The equilibrium single/non-ionized helium ratio $X_{\HeII}X_\e/X_{\HeI}$
 * through Saha equation, i.e.,
 * \begin{equation}\label{eq:saha:HeI}
 * \frac{X_{\HeII}X_\e}{X_{\HeI}} = \frac{e^{-\HeI_{1s}/(k_BT)}}{4n_{\Hy}\lambda_{\e}^3},
 * \end{equation}
 * where $\HeI_{1s}$ is the helium I $1s$ binding energy ncm_c_HeI_ion_E_1s_1S0().
 * This calculation is done using the Saha equation as in
 * [Seager (2000)][XSeager2000].
 *
 * The equilibrium double/single-ionized helium ratio $X_{\HeIII}X_\e/X_{\HeII}$
 * through Saha equation, i.e.,
 * \begin{equation}\label{eq:saha:HeII}
 * \frac{X_{\HeIII}X_\e}{X_{\HeII}} = \frac{e^{-\HeII_{1s}/(k_BT)}}{4n_{\Hy}\lambda_{\e}^3},
 * \end{equation}
 * where $\HeII_{1s}$ is the helium II $1s$ binding energy ncm_c_HeII_ion_E_1s_2S0_5().
 * This calculation is done using the Saha equation as in
 * [Seager (2000)][XSeager2000].
 *
 * The default value of the helium primordial abundance
 * is given by nc_hicosmo_Yp_4He().
 * The primordial Helium fraction is define by nc_hicosmo_XHe().
 *
 * <bridgehead>Optical depth and visibility function</bridgehead>
 *
 * The derivative of the optical depth $\tau$ with respect to the redshift
 * time $\lambda \equiv -\ln(x) = -\ln(1 + z)$ is
 * \begin{equation}\label{eq:def:dtaudlambda}
 * \frac{d\tau}{d\lambda} = -\frac{c\sigma_Tn_BX_\e}{H},
 * \end{equation}
 * where $c$ is the speed of light [ncm_c_c ()], $\sigma_T$ is the Thomson
 * cross section [ncm_c_thomson_cs ()], $n_B$ is the number density of baryons
 * and $H$ the Hubble function. We define the optical depth $\tau$ integrating
 * from the present time, i.e.,
 * \begin{equation}\label{eq:def:tau}
 * \tau = \int_0^\lambda\frac{d\tau}{d\lambda}.
 * \end{equation}
 * Using the equations above we define the visibility function $v_\tau$ as
 * \begin{equation}\label{eq:def:vtau}
 * v_\tau = \frac{d\tau}{d\lambda}e^{-\tau}.
 * \end{equation}
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_recomb.h"
#include "math/integral.h"
#include "math/ncm_util.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"
#include "math/ncm_spline_func.h"
#include "math/ncm_spline_cubic_notaknot.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_lambert.h>
#include <stdlib.h>
#endif /* NUMCOSMO_GIR_SCAN */

#define X_TODAY 1.0

G_DEFINE_ABSTRACT_TYPE (NcRecomb, nc_recomb, G_TYPE_OBJECT);

enum
{
  PROP_0,
  PROP_ZI,
  PROP_INIT_FRAC,
  PROP_PREC,
  PROP_SIZE,
};

static gdouble _nc_recomb_mdtau_dlambda (gdouble y, gdouble x, gpointer userdata);
static gdouble _nc_recomb_mdtau_drag_dlambda (gdouble y, gdouble x, gpointer userdata);

static void
nc_recomb_init (NcRecomb *recomb)
{
  NcmSpline *tau_s       = ncm_spline_cubic_notaknot_new ();
  NcmSpline *tau_drag_s  = ncm_spline_cubic_notaknot_new ();

  recomb->zi             = 0.0;
  recomb->prec           = 0.0;
  recomb->lambdai        = 0.0;
  recomb->lambdaf        = 0.0;
  recomb->init_frac      = 0.0;
  recomb->fsol           = gsl_root_fsolver_alloc (gsl_root_fsolver_brent);
  recomb->fmin           = gsl_min_fminimizer_alloc (gsl_min_fminimizer_brent);
  recomb->dtau_dlambda_s = NULL;
  recomb->tau_s          = NULL;
  recomb->tau_ode_s      = ncm_ode_spline_new (tau_s, &_nc_recomb_mdtau_dlambda);
  recomb->tau_drag_ode_s = ncm_ode_spline_new (tau_drag_s, &_nc_recomb_mdtau_drag_dlambda);
  recomb->ctrl_cosmo     = ncm_model_ctrl_new (NULL);
  recomb->ctrl_reion     = ncm_model_ctrl_new (NULL);

  recomb->v_tau_max_z    = 0.0;
  recomb->tau_z          = 0.0;
  recomb->tau_drag_z     = 0.0;
  recomb->tau_cutoff_z   = 0.0;
  
  ncm_spline_free (tau_s);
  ncm_spline_free (tau_drag_s);
}

static void
nc_recomb_dispose (GObject *object)
{
  NcRecomb *recomb = NC_RECOMB (object);

  ncm_spline_clear (&recomb->dtau_dlambda_s);
  ncm_spline_clear (&recomb->tau_s);

  ncm_ode_spline_clear (&recomb->tau_ode_s);
  ncm_ode_spline_clear (&recomb->tau_drag_ode_s);

  ncm_model_ctrl_clear (&recomb->ctrl_cosmo);
  ncm_model_ctrl_clear (&recomb->ctrl_reion);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_recomb_parent_class)->dispose (object);
}

static void
nc_recomb_finalize (GObject *object)
{
  NcRecomb *recomb = NC_RECOMB (object);

  gsl_root_fsolver_free (recomb->fsol);
  gsl_min_fminimizer_free (recomb->fmin);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_recomb_parent_class)->finalize (object);
}

static void
nc_recomb_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcRecomb *recomb = NC_RECOMB (object);
  g_return_if_fail (NC_IS_RECOMB (object));

  switch (prop_id)
  {
    case PROP_INIT_FRAC:
      recomb->init_frac = g_value_get_double (value);
      ncm_model_ctrl_force_update (recomb->ctrl_cosmo);
      break;
    case PROP_ZI:
      nc_recomb_set_zi (recomb, g_value_get_double (value));
      break;
    case PROP_PREC:
      recomb->prec = g_value_get_double (value);
      ncm_model_ctrl_force_update (recomb->ctrl_cosmo);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_recomb_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcRecomb *recomb = NC_RECOMB (object);
  g_return_if_fail (NC_IS_RECOMB (object));

  switch (prop_id)
  {
    case PROP_ZI:
      g_value_set_double (value, nc_recomb_get_zi (recomb));
      break;
    case PROP_INIT_FRAC:
      g_value_set_double (value, recomb->init_frac);
      break;
    case PROP_PREC:
      g_value_set_double (value, recomb->prec);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static gdouble _nc_recomb_Xe (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
static gdouble _nc_recomb_XHII (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
static gdouble _nc_recomb_XHeII (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda); 

static void
nc_recomb_class_init (NcRecombClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &nc_recomb_set_property;
  object_class->get_property = &nc_recomb_get_property;
  object_class->dispose      = &nc_recomb_dispose;
  object_class->finalize     = &nc_recomb_finalize;

  /**
   * NcRecomb:init_frac:
   *
   * $\HeIII/\He$ fraction to consider $\HeIII$ completly combined.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_INIT_FRAC,
                                   g_param_spec_double ("init-frac",
                                                        NULL,
                                                        "Initial fraction to start numerical integration",
                                                        0.0, 1.0, 1.0e-11,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcRecomb:zi:
   *
   * Initial redshift to prepare the recombination functions.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_ZI,
                                   g_param_spec_double ("zi",
                                                        NULL,
                                                        "Initial redshift for recombination calculations",
                                                        0.0, G_MAXDOUBLE, NC_RECOMB_STARTING_X,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcRecomb:prec:
   *
   * The precision used in the calculations.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_PREC,
                                   g_param_spec_double ("prec",
                                                        NULL,
                                                        "Precision for recombination calculations",
                                                        0.0, 1.0, 1e-7,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->prepare = NULL;
  klass->Xe      = &_nc_recomb_Xe;
  klass->XHII    = &_nc_recomb_XHII;
  klass->XHeII   = &_nc_recomb_XHeII;
}

static gdouble 
_nc_recomb_Xe (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda)
{
  g_error ("_nc_recomb_Xe: not implemented by `%s'", G_OBJECT_TYPE_NAME (recomb));
  return 0.0;
}

static gdouble 
_nc_recomb_XHII (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda)
{
  g_error ("_nc_recomb_XHII: not implemented by `%s'", G_OBJECT_TYPE_NAME (recomb));
  return 0.0;
}

static gdouble 
_nc_recomb_XHeII (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda)
{
  g_error ("_nc_recomb_XHeII: not implemented by `%s'", G_OBJECT_TYPE_NAME (recomb));
  return 0.0;
}

/**
 * nc_recomb_new_from_name:
 * @recomb_name: a string representing a #NcRecomb object
 *
 * This function instantiates a new object of type #NcRecomb given @recomb_name.
 *
 * Returns: a new #NcRecomb
 */
NcRecomb *
nc_recomb_new_from_name (const gchar *recomb_name)
{
  GObject *obj = ncm_serialize_global_from_string (recomb_name);
  GType multiplicity_type = G_OBJECT_TYPE (obj);
  if (!g_type_is_a (multiplicity_type, NC_TYPE_RECOMB))
    g_error ("nc_recomb_new_from_name: NcRecomb %s do not descend from %s.", recomb_name, g_type_name (NC_TYPE_RECOMB));

  return NC_RECOMB (obj);
}

/**
 * nc_recomb_ref:
 * @recomb: a #NcRecomb
 *
 * Increases the reference count of @recomb.
 *
 * Returns: (transfer full): @recomb
 */
NcRecomb *
nc_recomb_ref (NcRecomb *recomb)
{
  return NC_RECOMB (g_object_ref (recomb));
}

/**
 * nc_recomb_free:
 * @recomb: a #NcRecomb
 *
 * Decreases the reference count of @recomb
 *
 */
void
nc_recomb_free (NcRecomb *recomb)
{
  g_object_unref (recomb);
}

/**
 * nc_recomb_clear:
 * @recomb: a #NcRecomb
 *
 * Decreases the reference count of *@recomb if
 * *@recomb is not NULL, then sets *@recomb to NULL.
 *
 */
void
nc_recomb_clear (NcRecomb **recomb)
{
  g_clear_object (recomb);
}

/**
 * nc_recomb_prepare:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 *
 * Prepares the object using the model @cosmo.
 *
 */
void
nc_recomb_prepare (NcRecomb *recomb, NcHICosmo *cosmo)
{
  NC_RECOMB_GET_CLASS (recomb)->prepare (recomb, cosmo);

  ncm_model_ctrl_update (recomb->ctrl_cosmo, NCM_MODEL (cosmo));
}

/**
 * nc_recomb_prepare_if_needed:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 *
 * Prepares the object using the model @cosmo if it was changed
 * since last preparation.
 *
 */

/**
 * nc_recomb_set_zi:
 * @recomb: a #NcRecomb
 * @zi: initial redshift to prepare the recombination functions $z_i$
 * 
 * Sets the initial redshift of the integration $[0, z_i]$.
 * 
 */
void 
nc_recomb_set_zi (NcRecomb *recomb, const gdouble zi)
{
  if (recomb->zi != zi)
  {
    recomb->zi      = zi;
    recomb->lambdai = -log (1.0 + recomb->zi);
    ncm_model_ctrl_force_update (recomb->ctrl_cosmo);
  }
}

/**
 * nc_recomb_require_zi:
 * @recomb: a #NcRecomb
 * @zi: initial redshift to prepare the recombination functions $z_i$
 *
 * Requires the initial redshift of at least $z_i$ = @zi.
 *
 */
void 
nc_recomb_require_zi (NcRecomb *recomb, const gdouble zi)
{
  if (zi > recomb->zi)
  {
    recomb->zi = zi;
    ncm_model_ctrl_force_update (recomb->ctrl_cosmo);
  }
}

/**
 * nc_recomb_get_zi:
 * @recomb: a #NcRecomb
 *
 * Returns: the current initial redshift $z_i$.
 */
gdouble 
nc_recomb_get_zi (NcRecomb *recomb)
{
  return recomb->zi;
}

static gdouble _nc_recomb_XHeII (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda); 

static gdouble _nc_recomb_root (NcRecomb *recomb, gsl_function *F, gdouble x0, gdouble x1);
static gdouble _nc_recomb_min (NcRecomb *recomb, gsl_function *F, gdouble x0, gdouble x1, gdouble x);

/**
 * nc_recomb_HI_ion_saha:
 * @cosmo: a #NcHICosmo
 * @x: redshift factor, $x = 1+z$
 *
 * Calculates the equilibrium ionized/non-ionized hydrogen abundance ratio
 * $X_{\HyII}X_\e / X_{\HyI}$. See Eq. \eqref{eq:saha:HyI}.
 *
 * Returns: the abundance ratio $X_{\HyII}X_\e / X_{\HyI}$
 */
gdouble
nc_recomb_HI_ion_saha (NcHICosmo *cosmo, gdouble x)
{
  const gdouble T0 = nc_hicosmo_T_gamma0 (cosmo);
  const gdouble T = T0 * x;
  const gdouble kbT = ncm_c_kb () * T;
  const gdouble x3 = gsl_pow_3 (x);
  const gdouble lambda_e3 = gsl_pow_3 (ncm_c_thermal_wl_e () / sqrt (T));
  const gdouble n_H0 = nc_hicosmo_H_number_density (cosmo);
  const gdouble n_H  = n_H0 * x3;

  return gsl_sf_exp_mult (-ncm_c_HI_ion_E_1s_2S0_5 () / kbT, 1.0 / (n_H * lambda_e3));
}

/**
 * nc_recomb_HeI_ion_saha:
 * @cosmo: a #NcHICosmo
 * @x: redshift factor, $x = 1+z$
 *
 * Calculates the equilibrium single/non-ionized helium ratio
 * $X_{\HeII}X_\e/X_{\HeI}$. See Eq. \eqref{eq:saha:HeI}.
 *
 * Returns: the ratio $X_{\HeII}X_\e/X_{\HeI}$
 */
gdouble
nc_recomb_HeI_ion_saha (NcHICosmo *cosmo, gdouble x)
{
  const gdouble T0 = nc_hicosmo_T_gamma0 (cosmo);
  const gdouble T = T0 * x;
  const gdouble kbT = ncm_c_kb () * (T);
  const gdouble x3 = gsl_pow_3 (x);
  const gdouble lambda_e3 = gsl_pow_3 (ncm_c_thermal_wl_e () / sqrt(T));
  const gdouble n_H0 = nc_hicosmo_H_number_density (cosmo);
  const gdouble n_H  = n_H0 * x3;

  return gsl_sf_exp_mult (-ncm_c_HeI_ion_E_1s_1S0 () / kbT, 4.0 / (n_H * lambda_e3));
}

/**
 * nc_recomb_HeII_ion_saha:
 * @cosmo: a #NcHICosmo
 * @x: redshift factor, $x = 1+z$
 *
 * Calculates the equilibrium double/single ionized helium ratio
 * ($X_{\HeIII}X_\e/X_{\HeII}$). See Eq. \eqref{eq:saha:HeII}.
 *
 * Returns: the ratio $X_{\HeIII}X_\e/X_{\HeII}$
 */
gdouble
nc_recomb_HeII_ion_saha (NcHICosmo *cosmo, gdouble x)
{
  const gdouble T0 = nc_hicosmo_T_gamma0 (cosmo);
  const gdouble T = T0 * x;
  const gdouble kbT = ncm_c_kb () * T;
  const gdouble x3 = gsl_pow_3 (x);
  const gdouble lambda_e3 = gsl_pow_3 (ncm_c_thermal_wl_e () / sqrt(T0)) / sqrt(x3);
  const gdouble n_H0 = nc_hicosmo_H_number_density (cosmo);
  const gdouble n_H  = n_H0 * x3;

  return gsl_sf_exp_mult (-ncm_c_HeII_ion_E_1s_2S0_5 () / kbT, 1.0 / (n_H * lambda_e3));
}

/**
 * nc_recomb_HeII_ion_saha_x:
 * @cosmo: a #NcHICosmo
 * @f: $X_{\HeIII}X_\e/X_{\HeII}$
 *
 * Calculates the redshift where the ratio $$X_{\HeIII}X_\e/X_{\HeII} = f.$$
 *
 * This calculation is done by finding the value of $x$ where
 * $$\frac{e^{-\HeII_{1s}/(k_BT)}}{4n_{\Hy}\lambda_{\e}^3} = f.$$
 *
 * Returns: the value of $x$ where the ratio @f occur
 */
gdouble
nc_recomb_HeII_ion_saha_x (NcHICosmo *cosmo, gdouble f)
{
  const gdouble T0 = nc_hicosmo_T_gamma0 (cosmo);
  const gdouble kbT0 = ncm_c_kb () * T0;
  const gdouble mE_kbT0 = -ncm_c_HeII_ion_E_1s_2S0_5 () / kbT0;
  const gdouble lambda_e3_0 = gsl_pow_3 (ncm_c_thermal_wl_e () / sqrt (T0));
  const gdouble n_H0 = nc_hicosmo_H_number_density (cosmo);

  const gdouble y = 3.0 / 2.0 * gsl_sf_lambert_Wm1 (2.0 / 3.0 * mE_kbT0 * cbrt (gsl_pow_2 (lambda_e3_0 * n_H0 * f))) / mE_kbT0;

  return (1.0 / y);
}

/**
 * nc_recomb_HeII_ion_saha_x_by_HeIII_He:
 * @cosmo: a #NcHICosmo
 * @f: value of $X_{\HeIII}/X_{\He}$
 *
 * Calculates the redshift where the ratio $$X_{\HeIII}/X_{\He} = f.$$
 *
 * This calculation is done assuming that hydrogen and helium are
 * fully ionized, i.e., $\HyI = 0 = \HeI$. In this case
 * $$\frac{X_{\HeIII}X_\e}{X_{\HeII}} =
 * \frac{f}{1-f}\left[1 + X_\He(1+f)\right].$$
 *
 * Returns: the value of $x$ where the ratio @f occur
 *
 */
gdouble
nc_recomb_HeII_ion_saha_x_by_HeIII_He (NcHICosmo *cosmo, gdouble f)
{
  const gdouble XHe = nc_hicosmo_XHe (cosmo);
  const gdouble ratio = f * (1.0 + XHe * (1.0 + f)) / (1.0 - f);
  return nc_recomb_HeII_ion_saha_x (cosmo, ratio);
}

/**
 * nc_recomb_He_fully_ionized_Xe:
 * @cosmo: a #NcHICosmo
 * @x: redshift factor, $x = 1+z$
 *
 * Assuming that all helium is single or double ionized and all hydrogen is
 * ionized, we have $$X_\e = 1 + X_\HeII + 2X_\HeIII,\quad X_\He = X_\HeII +
 * X_\HeIII,$$ thus, $$X_\HeIII = X_\e-X_\He-1,\quad X_\HeII = 1 + 2X_\He -
 * X_\e.$$ Using nc_recomb_HeII_ion_saha() and nc_hicosmo_XHe() we obtain
 * $X_\e$.
 *
 * Returns: $X_\e$
 *
 */
gdouble
nc_recomb_He_fully_ionized_Xe (NcHICosmo *cosmo, gdouble x)
{
  const gdouble XHe = nc_hicosmo_XHe (cosmo);
  const gdouble XHeIIIXe_XHeII = nc_recomb_HeII_ion_saha (cosmo, x);
  const gdouble arg = XHe * (XHe + (2.0 + 6.0 * XHeIIIXe_XHeII)) / ((1.0 + XHeIIIXe_XHeII) * (1.0 + XHeIIIXe_XHeII));
  const gdouble Xe = (2.0 + XHe + (1.0 + XHeIIIXe_XHeII) * ncm_util_sqrt1px_m1 (arg)) / 2.0;
  return Xe;
}

typedef struct __nc_recomb_XHe_params
{
  gdouble fHI;
  gdouble fHeI;
  gdouble fHeII;
  gdouble XHe;
  NcHICosmo *cosmo;
} _nc_recomb_XHe_params;

static gdouble
_nc_recomb_XHe_XHII (gdouble XHII, gpointer p)
{
  _nc_recomb_XHe_params *XHe_p = (_nc_recomb_XHe_params *) p;
  const gdouble fHI   = XHe_p->fHI;
  const gdouble fHeI  = XHe_p->fHeI;
  const gdouble fHeII = XHe_p->fHeII;
  const gdouble XHe   = ((fHI * (-1.0 + XHII) + gsl_pow_2 (XHII)) * (gsl_pow_2 (fHI) * gsl_pow_2 (-1.0 + XHII) -
   fHeI * fHI * (-1.0 + XHII) * XHII +
   fHeI * fHeII * gsl_pow_2 (XHII))) / (fHeI * gsl_pow_2 (XHII) * (fHI * (-1.0 + XHII) - 2.0 * fHeII * XHII));

  return XHe - XHe_p->XHe;
}

static gdouble
_nc_recomb_XHe_XHI (gdouble XHI, gpointer p)
{
  _nc_recomb_XHe_params *XHe_p = (_nc_recomb_XHe_params *) p;
  const gdouble fHI   = XHe_p->fHI;
  const gdouble fHeI  = XHe_p->fHeI;
  const gdouble fHeII = XHe_p->fHeII;
  const gdouble XHe   = ((-fHI * XHI + gsl_pow_2 (1.0 - XHI)) * (gsl_pow_2 (fHI) * gsl_pow_2 (XHI) +
   fHeI * fHI * XHI * (1.0 - XHI) +
   fHeI * fHeII * gsl_pow_2 (1.0 - XHI))) / (fHeI * gsl_pow_2 (1.0 - XHI) * (-fHI * XHI - 2.0 * fHeII * (1.0 - XHI)));

  return XHe - XHe_p->XHe;
}


/**
 * nc_recomb_equilibrium_Xe:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 * @x: redshift factor, $x = 1+z$
 *
 * Calculates the ionization fraction $X_\e$ assuming equilibrium at all times.
 * It solves the system containing all Saha's equations Eqs \eqref{eq:saha:HyI},
 * \eqref{eq:saha:HeI} and \eqref{eq:saha:HeII} and the constraints Eq
 * \eqref{eq:Hy:add}, \eqref{eq:He:add} and \eqref{eq:def:Xe}.
 *
 * Returns: $X_\e$
 */
gdouble
nc_recomb_equilibrium_Xe (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble x)
{
   _nc_recomb_XHe_params XHe_p;
  gsl_function F;

  F.params = &XHe_p;

  g_assert_cmpfloat (x, >=, 0.0);

  XHe_p.cosmo  = cosmo;
  XHe_p.fHI    = nc_recomb_HI_ion_saha (cosmo, x);
  XHe_p.fHeI   = nc_recomb_HeI_ion_saha (cosmo, x);
  XHe_p.fHeII  = nc_recomb_HeII_ion_saha (cosmo, x);
  XHe_p.XHe    = nc_hicosmo_XHe (cosmo);

  if (XHe_p.fHeII < 1e-30)
  {
    return 0.5 * (sqrt (XHe_p.fHI * (4.0 + XHe_p.fHI)) - XHe_p.fHI);
  }

  if (XHe_p.fHI < 1.0e3)
  {
    F.function = &_nc_recomb_XHe_XHII;
    {
      const gdouble XHII = _nc_recomb_root (recomb, &F, 1e-30, 1.0 - GSL_DBL_EPSILON);
      return XHe_p.fHI * (-1.0 + 1.0 / XHII);
    }
  }
  else
  {
    F.function = &_nc_recomb_XHe_XHI;
    {
      const gdouble XHI = _nc_recomb_root (recomb, &F, 1e-30, 1.0 - GSL_DBL_EPSILON);
      return XHe_p.fHI * XHI / (1.0 - XHI);
    }
  }

}

/**
 * nc_recomb_equilibrium_XHI:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 * @x: redshift factor, $x = 1+z$
 *
 * Calculates the hydrogen-I fraction $X_\HyI$ assuming equilibrium at all times.
 * It solves the system containing all Saha's equations Eqs \eqref{eq:saha:HyI},
 * \eqref{eq:saha:HeI} and \eqref{eq:saha:HeII} and the constraints Eq
 * \eqref{eq:Hy:add}, \eqref{eq:He:add} and \eqref{eq:def:Xe}.
 *
 * Returns: $X_\HyI$
 */
gdouble
nc_recomb_equilibrium_XHI (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble x)
{
  const gdouble saha_XHII_Xe_XHI = nc_recomb_HI_ion_saha (cosmo, x);
  if (saha_XHII_Xe_XHI == 0.0)
  {
    return 1.0;
  }
  else
  {
    const gdouble Xe = nc_recomb_equilibrium_Xe (recomb, cosmo, x);
    const gdouble f  = saha_XHII_Xe_XHI / Xe;
    return 1.0 / (1.0 + f);
  }
}

/**
 * nc_recomb_equilibrium_XHII:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 * @x: redshift factor, $x = 1+z$
 *
 * Calculates the hydrogen-II fraction $X_\HyII$ assuming equilibrium at all times.
 * It solves the system containing all Saha's equations Eqs \eqref{eq:saha:HyI},
 * \eqref{eq:saha:HeI} and \eqref{eq:saha:HeII} and the constraints Eq
 * \eqref{eq:Hy:add}, \eqref{eq:He:add} and \eqref{eq:def:Xe}.
 *
 * Returns: $X_\HyII$
 */
gdouble
nc_recomb_equilibrium_XHII (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble x)
{
  const gdouble saha_XHII_Xe_XHI = nc_recomb_HI_ion_saha (cosmo, x);
  if (saha_XHII_Xe_XHI == 0.0)
  {
    return 0.0;
  }
  else
  {
    const gdouble Xe = nc_recomb_equilibrium_Xe (recomb, cosmo, x);
    const gdouble f  = saha_XHII_Xe_XHI / Xe;
    return f / (1.0 + f);
  }
}

/**
 * nc_recomb_equilibrium_XHeI:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 * @x: redshift factor, $x = 1+z$
 *
 * Calculates the helium-I fraction $X_\HeI$ assuming equilibrium at all times.
 * It solves the system containing all Saha's equations Eqs \eqref{eq:saha:HyI},
 * \eqref{eq:saha:HeI} and \eqref{eq:saha:HeII} and the constraints Eq
 * \eqref{eq:Hy:add}, \eqref{eq:He:add} and \eqref{eq:def:Xe}.
 *
 * Returns: $X_\HyII$
 */
gdouble
nc_recomb_equilibrium_XHeI (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble x)
{
  const gdouble saha_XHeIII_Xe_XHeII = nc_recomb_HeII_ion_saha (cosmo, x);
  const gdouble saha_XHeII_Xe_XHeI   = nc_recomb_HeI_ion_saha (cosmo, x);
  const gdouble Xe                   = nc_recomb_equilibrium_Xe (recomb, cosmo, x);
  const gdouble f1                   = saha_XHeIII_Xe_XHeII / Xe;
  const gdouble f2                   = saha_XHeII_Xe_XHeI / Xe;

  return nc_hicosmo_XHe (cosmo) / (1.0 + f2 + f1 * f2);
}

/**
 * nc_recomb_equilibrium_XHeII:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 * @x: redshift factor, $x = 1+z$
 *
 * Calculates the helium-II fraction $X_\HeII$ assuming equilibrium at all times.
 * It solves the system containing all Saha's equations Eqs \eqref{eq:saha:HyI},
 * \eqref{eq:saha:HeI} and \eqref{eq:saha:HeII} and the constraints Eq
 * \eqref{eq:Hy:add}, \eqref{eq:He:add} and \eqref{eq:def:Xe}.
 *
 * Returns: $X_\HyII$
 */
gdouble
nc_recomb_equilibrium_XHeII (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble x)
{
  const gdouble saha_XHeIII_Xe_XHeII = nc_recomb_HeII_ion_saha (cosmo, x);
  const gdouble saha_XHeII_Xe_XHeI   = nc_recomb_HeI_ion_saha (cosmo, x);
  const gdouble Xe                   = nc_recomb_equilibrium_Xe (recomb, cosmo, x);
  const gdouble f1                   = saha_XHeIII_Xe_XHeII / Xe;
  const gdouble f2                   = saha_XHeII_Xe_XHeI / Xe;

  return f2 * nc_hicosmo_XHe (cosmo) / (1.0 + f2 + f1 * f2);
}

/**
 * nc_recomb_equilibrium_XHeIII:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 * @x: redshift factor, $x = 1+z$
 *
 * Calculates the helium-III fraction $X_\HeIII$ assuming equilibrium at all times.
 * It solves the system containing all Saha's equations Eqs \eqref{eq:saha:HyI},
 * \eqref{eq:saha:HeI} and \eqref{eq:saha:HeII} and the constraints Eq
 * \eqref{eq:Hy:add}, \eqref{eq:He:add} and \eqref{eq:def:Xe}.
 *
 * Returns: $X_\HyII$
 */
gdouble
nc_recomb_equilibrium_XHeIII (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble x)
{
  const gdouble saha_XHeIII_Xe_XHeII = nc_recomb_HeII_ion_saha (cosmo, x);
  const gdouble saha_XHeII_Xe_XHeI   = nc_recomb_HeI_ion_saha (cosmo, x);
  const gdouble Xe                   = nc_recomb_equilibrium_Xe (recomb, cosmo, x);
  const gdouble f1                   = saha_XHeIII_Xe_XHeII / Xe;
  const gdouble f2                   = saha_XHeII_Xe_XHeI / Xe;

  return f1 * f2 * nc_hicosmo_XHe (cosmo) / (1.0 + f2 + f1 * f2);
}

/**
 * nc_recomb_dtau_dlambda_Xe:
 * @cosmo: a #NcHICosmo
 * @lambda: $\lambda = - \ln(1+z)$
 *
 * The derivative of the optical depth [Eq. \eqref{eq:def:dtaudlambda}] over
 * the ionization fraction $X_\e$ [Eq. \eqref{eq:def:Xe}].
 *
 * Returns: $X_\e^{-1}d\tau/d\lambda$
 */
/**
 * nc_recomb_He_fully_ionized_dtau_dlambda:
 * @cosmo: a #NcHICosmo
 * @lambda: $\lambda = - \ln(1+z)$
 *
 * The derivative of the optical depth [Eq. \eqref{eq:def:dtaudlambda}], considering
 * fully ionized helium and hydrogen [nc_recomb_He_fully_ionized_Xe ()].
 *
 * Returns: $d\tau/d\lambda$
 */


/**
 * nc_recomb_dtau_dx:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 * @lambda: $\lambda = - \ln(1+z)$
 *
 * Computes the derivative of the optical depth with respect to $x = 1 + z$, 
 * $\frac{\mathrm{d}\tau}{\mathrm{d}x}$, at @lambda.
 *
 * Returns: $\left.\frac{\mathrm{d}\tau}{\mathrm{d}x}\right|_{\lambda}$
 */
gdouble
nc_recomb_dtau_dx (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda)
{
  NCM_UNUSED (cosmo);
  return gsl_sf_exp_mult (lambda, - ncm_spline_eval (recomb->dtau_dlambda_s, lambda));
}

/**
 * nc_recomb_dtau_dlambda:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 * @lambda: $\lambda = - \ln(1+z)$
 *
 * Computes the derivative of the optical depth with respect to $\lambda$, 
 * $\frac{\mathrm{d}\tau}{\mathrm{d}\lambda}$, at @lambda.
 *
 * Returns: $\left.\frac{\mathrm{d}\tau}{\mathrm{d}\lambda}\right|_{\lambda}$
 */
gdouble
nc_recomb_dtau_dlambda (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda)
{
  NCM_UNUSED (cosmo);
  return ncm_spline_eval (recomb->dtau_dlambda_s, lambda);
}

/**
 * nc_recomb_d2tau_dlambda2:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 * @lambda: $\lambda = - \ln(1+z)$
 *
 * Computes the second derivative of the optical depth with respect to $\lambda$, 
 * $\frac{\mathrm{d}^2\tau}{\mathrm{d}\lambda^2}$, at @lambda.
 *
 * Returns: $\left.\frac{\mathrm{d}^2\tau}{\mathrm{d}\lambda^2}\right|_{\lambda}$
 */
gdouble
nc_recomb_d2tau_dlambda2 (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda)
{
  NCM_UNUSED (cosmo);
  return ncm_spline_eval_deriv (recomb->dtau_dlambda_s, lambda);
}

/**
 * nc_recomb_d3tau_dlambda3:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 * @lambda: $\lambda = - \ln(1+z)$
 *
 * Computes the third derivative of the optical depth with respect to $\lambda$, 
 * $\frac{\mathrm{d}^3\tau}{\mathrm{d}\lambda^3}$, at @lambda.
 *
 * Returns: $\left.\frac{\mathrm{d}^3\tau}{\mathrm{d}\lambda^3}\right|_{\lambda}$
 */
gdouble
nc_recomb_d3tau_dlambda3 (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda)
{
  NCM_UNUSED (cosmo);
  return ncm_spline_eval_deriv2 (recomb->dtau_dlambda_s, lambda);
}

/**
 * nc_recomb_tau:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 * @lambda: $\lambda = - \ln(1+z)$
 *
 * Computes the optical depth [Eq. \eqref{eq:def:tau}] at @lambda.
 *
 * Returns: $\tau$
 */
gdouble
nc_recomb_tau (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda)
{
  NCM_UNUSED (cosmo);
  return ncm_spline_eval (recomb->tau_s, -lambda);
}

/**
 * nc_recomb_tau_drag:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 * @lambda: $\lambda = - \ln(1+z)$
 *
 * Computes the drag depth,
 * $$\tau_d \equiv \int_0^\lambda \mathrm{d}\lambda R^{-1} \frac{n_e \sigma_T}{1+z},$$ 
 * where $R = \frac{3\bar{\rho}_b}{4\bar{\rho}_\gamma}$.
 * 
 * Returns: $\tau_\mathrm{drag}$
 */
gdouble
nc_recomb_tau_drag (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda)
{
  NCM_UNUSED (cosmo);
  return ncm_spline_eval (ncm_ode_spline_peek_spline (recomb->tau_drag_ode_s), -lambda);
}

/**
 * nc_recomb_tau_lambda0_lambda1:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 * @lambda0: $\lambda_0$
 * @lambda1: $\lambda_1$
 *
 * Computes the optical depth between @lambda1 and @lambda0. 
 *
 * Returns: $\tau(\lambda_1) - \tau(\lambda_0)$
 */
gdouble
nc_recomb_tau_lambda0_lambda1 (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda0, const gdouble lambda1)
{
  NCM_UNUSED (cosmo);
  return ncm_spline_eval (recomb->tau_s, -lambda1) -
    ncm_spline_eval (recomb->tau_s, -lambda0);
}

/**
 * nc_recomb_log_v_tau:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 * @lambda: $\lambda = - \ln(1+z)$
 *
 * Computes the logarithm base e of the visibility function [Eq. \eqref{eq:def:vtau}] at @lambda.
 *
 * Returns: $\log(v_\tau)$
 */
gdouble
nc_recomb_log_v_tau (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda)
{
  const gdouble tau          = nc_recomb_tau (recomb, cosmo, lambda);
  const gdouble dtau_dlambda = nc_recomb_dtau_dlambda (recomb, cosmo, lambda);

  return - tau + log (fabs (dtau_dlambda));
}

/**
 * nc_recomb_v_tau:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 * @lambda: $\lambda = - \ln(1+z)$
 *
 * Computes the visibility function [Eq. \eqref{eq:def:vtau}] at @lambda.
 *
 * Returns: $v_\tau$
 */
gdouble
nc_recomb_v_tau (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda)
{
  const gdouble tau          = nc_recomb_tau (recomb, cosmo, lambda);
  const gdouble dtau_dlambda = nc_recomb_dtau_dlambda (recomb, cosmo, lambda);

  return gsl_sf_exp_mult (-tau, dtau_dlambda);
}

/**
 * nc_recomb_dv_tau_dlambda:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 * @lambda: $\lambda = - \ln(1+z)$
 *
 * Computes the derivative of the visibility function [Eq. \eqref{eq:def:vtau}] at @lambda.
 *
 * Returns: $\left. \frac{\mathrm{d}v_\tau}{\mathrm{d}\lambda}\right|_\lambda$
 */
gdouble
nc_recomb_dv_tau_dlambda (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda)
{
  const gdouble tau = nc_recomb_tau (recomb, cosmo, lambda);
  const gdouble dtau_dlambda = nc_recomb_dtau_dlambda (recomb, cosmo, lambda);
  const gdouble d2tau_dlambda2 = nc_recomb_d2tau_dlambda2 (recomb, cosmo, lambda);
  return gsl_sf_exp_mult (-tau, (d2tau_dlambda2 - dtau_dlambda * dtau_dlambda));
}

/**
 * nc_recomb_d2v_tau_dlambda2:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 * @lambda: $\lambda = - \ln(1+z)$
 *
 * Computes the second derivative of the visibility function [Eq. \eqref{eq:def:vtau}] at @lambda.
 *
 * Returns: $\left. \frac{\mathrm{d}^2v_\tau}{\mathrm{d}\lambda^2}\right|_\lambda$
 */
gdouble
nc_recomb_d2v_tau_dlambda2 (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda)
{
  const gdouble tau = nc_recomb_tau (recomb, cosmo, lambda);
  const gdouble dtau_dlambda = nc_recomb_dtau_dlambda (recomb, cosmo, lambda);
  const gdouble d2tau_dlambda2 = nc_recomb_d2tau_dlambda2 (recomb, cosmo, lambda);
  const gdouble d3tau_dlambda3 = nc_recomb_d3tau_dlambda3 (recomb, cosmo, lambda);
  return gsl_sf_exp_mult (-tau, (gsl_pow_3 (dtau_dlambda) + d3tau_dlambda3
                                 - 3.0 * dtau_dlambda * d2tau_dlambda2));
}

/**
 * nc_recomb_Xe: (virtual Xe)
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 * @lambda: $\lambda = - \ln(1+z)$
 *
 * Calculates the value of $X_\e$ at @lambda.
 *
 * Returns: $X_\e$
 */
/**
 * nc_recomb_XHII: (virtual XHII)
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 * @lambda: $\lambda = - \ln(1+z)$
 *
 * Calculates the value of $X_\HyII$ at @lambda.
 *
 * Returns: $X_\HyII$
 */
/**
 * nc_recomb_XHeII: (virtual XHeII)
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 * @lambda: $\lambda = - \ln(1+z)$
 *
 * Calculates the value of $X_\HeII$ at @lambda.
 *
 * Returns: $X_\HeII$
 */

static gdouble
_nc_recomb_root (NcRecomb *recomb, gsl_function *F, gdouble x0, gdouble x1)
{
  gint status;
  gint iter = 0, max_iter = 1000000;
  gdouble x = 0.0;
  const gdouble prec = GSL_MIN (recomb->prec, 1e-1);

  gsl_root_fsolver_set (recomb->fsol, F, x0, x1);
  do {
    iter++;
    status = gsl_root_fsolver_iterate (recomb->fsol);
    if (status)
      g_error ("_nc_recomb_root_brent: Cannot find root (%s)", gsl_strerror (status));
    x = gsl_root_fsolver_root (recomb->fsol);
    x0 = gsl_root_fsolver_x_lower (recomb->fsol);
    x1 = gsl_root_fsolver_x_upper (recomb->fsol);

    status = gsl_root_test_interval (x0, x1, 0, prec);
  } while (status == GSL_CONTINUE && iter < max_iter);

  return x;
}

static gdouble
_nc_recomb_min (NcRecomb *recomb, gsl_function *F, gdouble x0, gdouble x1, gdouble x)
{
  gint status;
  gint iter = 0, max_iter = 1000000;
  const gdouble prec = GSL_MIN (recomb->prec, 1e-1);
  gdouble lx = 0.0, lx0 = 0.0, lx1 = 0.0;

  gsl_min_fminimizer_set (recomb->fmin, F, x, x0, x1);

  do {
    iter++;
    status = gsl_min_fminimizer_iterate (recomb->fmin);
    if (status)
      g_error ("_nc_recomb_min: Cannot find minimum (%s)", gsl_strerror (status));

    x  = gsl_min_fminimizer_x_minimum (recomb->fmin);
    x0 = gsl_min_fminimizer_x_lower (recomb->fmin);
    x1 = gsl_min_fminimizer_x_upper (recomb->fmin);

    status = gsl_min_test_interval (x0, x1, 0.0, prec);

    if (x == lx && x0 == lx0 && x1 == lx1)
    {
      const gdouble fl = gsl_min_fminimizer_f_lower (recomb->fmin);
      const gdouble fu = gsl_min_fminimizer_f_upper (recomb->fmin);
      if (fabs (2.0 * (fu - fl) / (fu + fl)) < prec)
        status = GSL_SUCCESS;
      else
        g_error ("_nc_recomb_min: Cannot find minimum.");
    }

    lx = x; lx0 = x0; lx1 = x1;
  } while (status == GSL_CONTINUE && iter < max_iter);

  return x;
}

typedef struct __nc_recomb_func
{
  NcRecomb *recomb;
  NcHICosmo *cosmo;
  gdouble ref;
} _nc_recomb_func;

static gdouble
_nc_recomb_v_tau_features (gdouble lambda, gpointer params)
{
  _nc_recomb_func *func = (_nc_recomb_func *) params;
  return nc_recomb_log_v_tau (func->recomb, func->cosmo, lambda) - func->ref;
}

static gdouble
_nc_recomb_v_tau_min (gdouble lambda, gpointer params)
{
  _nc_recomb_func *func = (_nc_recomb_func *) params;
  return -nc_recomb_log_v_tau (func->recomb, func->cosmo, lambda);
}

/**
 * nc_recomb_v_tau_lambda_features:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 * @logref: the logarithm of the reference scale
 * @lambda_max: (out): $\lambda_\text{max}$
 * @lambda_l: (out): $\lambda_l$
 * @lambda_u: (out): $\lambda_u$
 *
 * Calculates the maximum of the visibility function [Eq \eqref{eq:def:vtau}], i.e,
 * the value of $\lambda_\text{max}$ where
 * $dv_\tau(\lambda_\text{max})/d\lambda = 0$, and the values where the
 * visibility drops to $v_\tau(\lambda_\text{max})e^{-\text{logref}}$ to the left
 * $\lambda_l$ and to the right $\lambda_u$ of $\lambda_\text{max}$.
 *
 */
void
nc_recomb_v_tau_lambda_features (NcRecomb *recomb, NcHICosmo *cosmo, gdouble logref, gdouble *lambda_max, gdouble *lambda_l, gdouble *lambda_u)
{
  _nc_recomb_func func;
  gsl_function F;
  gdouble log_v_tau_max;
  gdouble log_v_tau_f;
  const gdouble cmb_ref_z  = 1090.0;
  const gdouble lambda_ref = -log (10.0 * cmb_ref_z + 10.0);
  const gdouble lambda_try = -log (cmb_ref_z + 1.0);

  F.function = &_nc_recomb_v_tau_min;
  F.params   = &func;

  func.recomb = recomb;
  func.cosmo  = cosmo;

  *lambda_max = _nc_recomb_min (recomb, &F, lambda_ref, recomb->lambdaf, lambda_try);

  log_v_tau_max = -_nc_recomb_v_tau_min (*lambda_max, &func);
  log_v_tau_f   = -_nc_recomb_v_tau_min (recomb->lambdaf, &func);
  func.ref      = log_v_tau_max - logref;

  func.ref      = GSL_MAX (func.ref, (log_v_tau_max + log_v_tau_f) * 0.5);

  F.function    = &_nc_recomb_v_tau_features;

  *lambda_l     = _nc_recomb_root (recomb, &F, lambda_ref, *lambda_max);
  *lambda_u     = _nc_recomb_root (recomb, &F, *lambda_max, recomb->lambdaf);
}

static gdouble
_nc_recomb_tau_ref (gdouble lambda, gpointer params)
{
  _nc_recomb_func *func = (_nc_recomb_func *) params;
  return nc_recomb_tau (func->recomb, func->cosmo, lambda) - func->ref;
}

static gdouble
_nc_recomb_tau_drag_ref (gdouble lambda, gpointer params)
{
  _nc_recomb_func *func = (_nc_recomb_func *) params;
  return nc_recomb_tau_drag (func->recomb, func->cosmo, lambda) - func->ref;
}

/**
 * nc_recomb_get_v_tau_max_lambda:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 *
 * Calculates the maximum of the visibility function [Eq. \eqref{eq:def:vtau}],
 * the value of $\lambda_\text{max}$ where
 * $dv_\tau(\lambda_\text{max})/d\lambda = 0$.
 *
 * Returns: $\lambda_\text{max}$
 */
/**
 * nc_recomb_get_tau_lambda:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 *
 * Calculates the value of $\lambda$ where the optical depth
 * [Eq \eqref{eq:def:tau}] is equal to one, i.e., $\tau(\lambda^\star) = 1$.
 *
 * Returns: $\lambda^\star$
 */
/**
 * nc_recomb_get_tau_drag_lambda:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 *
 * Calculates the value of $\lambda$ where the optical depth times $R$
 * [Eq. \eqref{eq:def:tau}] is equal to one, i.e., $\tau_\mathrm{drag}(\lambda^\star) = 1$.
 *
 * Returns: $\lambda^\star$
 */
/**
 * nc_recomb_get_tau_cutoff_lambda:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 *
 * Calculates the value of $\lambda$ where the optical depth
 * [Eq \eqref{eq:def:tau}] attains a value such that
 * $e^{-\tau(\lambda_\text{cutoff})} = \epsilon_\text{double}$, i.e.,
 * equal to the smallest value of a double which add to one.
 *
 * Returns: $\lambda_\text{cutoff}$
 */

/**
 * nc_recomb_get_v_tau_max_z:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 *
 * Calculates the maximum of the visibility function [Eq. \eqref{eq:def:vtau}],
 * the value of $z(\lambda_\text{max})$ where
 * $dv_\tau(\lambda_\text{max})/d\lambda = 0$.
 *
 * Returns: $z(\lambda_\text{max})$
 */
/**
 * nc_recomb_get_tau_z:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 *
 * Calculates the value of $z(\lambda)$ where the optical depth
 * [Eq \eqref{eq:def:tau}] is equal to one, i.e., $\tau(\lambda^\star) = 1$.
 *
 * Returns: $z(\lambda^\star)$
 */
/**
 * nc_recomb_get_tau_drag_z:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 *
 * Calculates the value of $z(\lambda)$ where the optical depth times $R$ (drag depth)
 * [Eq. \eqref{eq:def:tau}] is equal to one, i.e., $\tau_\mathrm{drag}(\lambda^\star) = 1$.
 *
 * Returns: $z(\lambda^\star)$
 */
/**
 * nc_recomb_get_tau_cutoff_z:
 * @recomb: a #NcRecomb
 * @cosmo: a #NcHICosmo
 *
 * Calculates the value of $z(\lambda)$ where the optical depth
 * [Eq \eqref{eq:def:tau}] attains a value such that
 * $e^{-\tau(\lambda_\text{cutoff})} = \epsilon_\text{double}$, i.e.,
 * equal to the smallest value of a double which add to one.
 *
 * Returns: $z(\lambda_\text{cutoff})$
 */

static gdouble
_nc_recomb_dtau_dlambda (gdouble lambda, gpointer p)
{
  _nc_recomb_func *func = (_nc_recomb_func *) p;
  const gdouble Xe      = nc_recomb_Xe (func->recomb, func->cosmo, lambda);
  return Xe * nc_recomb_dtau_dlambda_Xe (func->cosmo, lambda);
}

static gdouble 
_nc_recomb_mdtau_dlambda (gdouble y, gdouble x, gpointer userdata)
{
  const gdouble lambda = - x;
  return - _nc_recomb_dtau_dlambda (lambda, userdata);
}

static gdouble 
_nc_recomb_mdtau_drag_dlambda (gdouble y, gdouble x, gpointer userdata)
{
  _nc_recomb_func *func = (_nc_recomb_func *) userdata;
  const gdouble lambda = - x;
  const gdouble z      = expm1 (x);
  const gdouble R      = 4.0 / 3.0 * nc_hicosmo_E2Omega_g (func->cosmo, z) / nc_hicosmo_E2Omega_b (func->cosmo, z);
  
  return -R * _nc_recomb_dtau_dlambda (lambda, userdata);
}

void
_nc_recomb_prepare_tau_splines (NcRecomb *recomb, NcHICosmo *cosmo)
{
  _nc_recomb_func func;
  gsl_function F;

  if (recomb->tau_s == NULL)
    recomb->tau_s = ncm_spline_cubic_notaknot_new ();

  if (recomb->dtau_dlambda_s == NULL)
    recomb->dtau_dlambda_s = ncm_spline_cubic_notaknot_new ();

  func.recomb = recomb;
  func.cosmo  = cosmo;

  F.params   = &func;

  F.function = &_nc_recomb_dtau_dlambda;
  ncm_spline_set_func (recomb->dtau_dlambda_s, NCM_SPLINE_FUNCTION_SPLINE,
                       &F, recomb->lambdai, recomb->lambdaf, 0, recomb->prec);

  ncm_ode_spline_auto_abstol (recomb->tau_ode_s, TRUE);
  ncm_ode_spline_set_interval (recomb->tau_ode_s, 0.0, -recomb->lambdaf, -recomb->lambdai);
  ncm_ode_spline_prepare (recomb->tau_ode_s, &func);

  ncm_ode_spline_auto_abstol (recomb->tau_drag_ode_s, TRUE);
  ncm_ode_spline_set_interval (recomb->tau_drag_ode_s, 0.0, -recomb->lambdaf, -recomb->lambdai);
  ncm_ode_spline_prepare (recomb->tau_drag_ode_s, &func);
  
  ncm_spline_clear (&recomb->tau_s);
  recomb->tau_s = ncm_spline_ref (ncm_ode_spline_peek_spline (recomb->tau_ode_s));
}

static gdouble
_nc_recomb_v_tau_max_lambda (NcRecomb *recomb, NcHICosmo *cosmo)
{
  _nc_recomb_func func;
  gsl_function F;
  const gdouble cmb_ref_z  = 1090.0;
  const gdouble lambda_ref = -log (10.0 * cmb_ref_z + 10.0);
  const gdouble lambda_try = -log (cmb_ref_z + 1.0);

  F.function = &_nc_recomb_v_tau_min;
  F.params   = &func;

  func.recomb = recomb;
  func.cosmo  = cosmo;

  return _nc_recomb_min (recomb, &F, lambda_ref, recomb->lambdaf, lambda_try);
}

static gdouble
_nc_recomb_tau_lambda (NcRecomb *recomb, NcHICosmo *cosmo)
{
  _nc_recomb_func func;
  gsl_function F;
  const gdouble cmb_ref_z = 1090.0;

  F.function = &_nc_recomb_tau_ref;
  F.params   = &func;

  func.recomb = recomb;
  func.cosmo  = cosmo;
  func.ref    = 1.0;

  return _nc_recomb_root (recomb, &F,
                          -log (10.0 * cmb_ref_z + 10.0),
                          recomb->lambdaf);
}

static gdouble
_nc_recomb_tau_drag_lambda (NcRecomb *recomb, NcHICosmo *cosmo)
{
  _nc_recomb_func func;
  gsl_function F;
  const gdouble cmb_ref_z = 1060.0;

  F.function = &_nc_recomb_tau_drag_ref;
  F.params   = &func;

  func.recomb = recomb;
  func.cosmo  = cosmo;
  func.ref    = 1.0;

  return _nc_recomb_root (recomb, &F,
                          -log (10.0 * cmb_ref_z + 10.0),
                          recomb->lambdaf);
}

static gdouble
_nc_recomb_tau_cutoff_lambda (NcRecomb *recomb, NcHICosmo *cosmo)
{
  _nc_recomb_func func;
  gsl_function F;

  F.function = &_nc_recomb_tau_ref;
  F.params = &func;

  func.recomb = recomb;
  func.cosmo  = cosmo;
  func.ref    = -GSL_LOG_DBL_EPSILON;

  return _nc_recomb_root (recomb, &F, 
                          recomb->lambdai,
                          recomb->lambdaf);
}

void
_nc_recomb_prepare_redshifts (NcRecomb *recomb, NcHICosmo *cosmo)
{
  recomb->v_tau_max_lambda  = _nc_recomb_v_tau_max_lambda (recomb, cosmo);
  recomb->tau_lambda        = _nc_recomb_tau_lambda (recomb, cosmo);
  recomb->tau_drag_lambda   = _nc_recomb_tau_drag_lambda (recomb, cosmo);
  recomb->tau_cutoff_lambda = _nc_recomb_tau_cutoff_lambda (recomb, cosmo);

  recomb->v_tau_max_z       = expm1 (-recomb->v_tau_max_lambda);
  recomb->tau_z             = expm1 (-recomb->tau_lambda);
  recomb->tau_drag_z        = expm1 (-recomb->tau_drag_lambda);
  recomb->tau_cutoff_z      = expm1 (-recomb->tau_cutoff_lambda);
}
