/***************************************************************************
 *            recomb.c
 *
 *  Sun Oct  5 20:40:30 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * SECTION:recomb
 * @title: Recombination
 * @short_description: Cosmic recombination object
 *
 * The default value of the <link linkend="def_Y_p">helium primordial abundance</link>
 * is given the macro #NC_C_PRIM_HE_Y_P.
 * The primordial helium fraction is define by #NC_C_PRIM_HE_XHe.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "thermodyn/recomb.h"
#include "math/integral.h"
#include "math/cvode_util.h"
#include "math/util.h"
#include "perturbations/linear.h"

#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>
#include <nvector/nvector_serial.h>

#define X_TODAY 1.0

/**
 * nc_thermodyn_H_ionization_saha:
 * @model: a #NcmModel, cosmological parameters for a given model
 * @x: normalized scale factor inverse [; x = 1 + z = a_0/a ;]
 *
 * @brief Calculate the equilibrium ionized/non-ionized hydrogen abundance ratio [; X_{H^+}X_e / X_{H} ;].
 *
 * This calculation is done using the saha equation as in <link linkend="XWeinberg2008">Weinberg's cosmology book</link>
 * page XXX. Note that this function is valid only in the equilibrium. See also <link linkend="numcosmo-Recombination">Recombination</link>.
 *
 * Returns: the abundance ratio [; X_{H^+}X_e / X_{H} ;]
 */
gdouble
nc_thermodyn_H_ionization_saha (NcHICosmo *model, gdouble x)
{
  const gdouble T0 = nc_hicosmo_T_gamma0 (model);
  const gdouble T = T0 * x;
  const gdouble kbT = NC_C_kb * (T);
  const gdouble x3 = gsl_pow_3 (x);
  const gdouble h2 = nc_hicosmo_h2 (model);
  const gdouble Omega_b = nc_hicosmo_Omega_b (model);
  const gdouble lambda_e3 = gsl_pow_3 (NC_C_THERMAL_WAVELENGTH_e / sqrt(T0)) / sqrt(x3);
  const gdouble n_H = NC_C_PRIM_H_FRAC * Omega_b * x3 * (NC_C_CRIT_DENSITY * h2) / NC_C_ENERGY_p;

  return gsl_sf_exp_mult (-NC_C_HYDROGEN_BINDING_ENERGY_1s / kbT, 1.0 / (n_H * lambda_e3)); /* FIXME */
}

/**
 * nc_thermodyn_HeI_ionization_saha:
 * @model: a #NcmModel, cosmological parameters for a given model
 * @x: normalized scale factor inverse [; x = 1 + z = a_0/a ;]
 *
 * @brief Calculate the equilibrium single/non-ionized helium ratio ([; X_{He^{+}}X_e/X_{He} ;]).
 *
 * This calculation is done using the saha equation as in \ref seager1999 .
 * Note that this function is valid only in the equilibrium. See also \ref recomb_sec_var.
 *
 * Returns: the ratio [; X_{He^{+}}X_e/X_{He} ;]
 */
gdouble
nc_thermodyn_HeI_ionization_saha (NcHICosmo *model, gdouble x)
{
  const gdouble T0 = nc_hicosmo_T_gamma0 (model);
  const gdouble T = T0 * x;
  const gdouble kbT = NC_C_kb*(T);
  const gdouble x3 = gsl_pow_3 (x);
  const gdouble h2 = nc_hicosmo_h2 (model);
  const gdouble Omega_b = nc_hicosmo_Omega_b (model);
  const gdouble lambda_e3 = gsl_pow_3 (NC_C_THERMAL_WAVELENGTH_e / sqrt(T0)) / sqrt(x3);
  const gdouble n_H = NC_C_PRIM_H_FRAC * Omega_b * x3 * (NC_C_CRIT_DENSITY * h2) / NC_C_ENERGY_p;

  return gsl_sf_exp_mult (-NC_C_HeI_BINDING_ENERGY_1s / kbT, 4.0 / (n_H * lambda_e3));
}

/**
 * nc_thermodyn_HeII_ionization_saha:
 * @model: a #NcmModel, cosmological parameters for a given model
 * @x: normalized scale factor inverse [; x = 1 + z = a_0/a ;]
 *
 * @brief Calculate the equilibrium double/single ionized helium ratio ([; X_{He^{++}}X_e/X_{He^+} ;]).
 *
 * This calculation is done using the saha equation as in \ref seager1999 .
 * Note that this function is valid only in the equilibrium. See also \ref recomb_sec_var.
 *
 * Returns: the ratio [; X_{He^{++}}X_e/X_{He^+} ;]
 */
gdouble
nc_thermodyn_HeII_ionization_saha (NcHICosmo *model, gdouble x)
{
  const gdouble T0 = nc_hicosmo_T_gamma0 (model);
  const gdouble T = T0 * x;
  const gdouble kbT = NC_C_kb * T;
  const gdouble x3 = gsl_pow_3 (x);
  const gdouble h2 = nc_hicosmo_h2 (model);
  const gdouble Omega_b = nc_hicosmo_Omega_b (model);
  const gdouble lambda_e3 = gsl_pow_3 (NC_C_THERMAL_WAVELENGTH_e / sqrt(T0)) / sqrt(x3);
  const gdouble n_H = NC_C_PRIM_H_FRAC * Omega_b * x3 * (NC_C_CRIT_DENSITY * h2) / NC_C_ENERGY_p;

  return gsl_sf_exp_mult (-NC_C_HeII_BINDING_ENERGY_1s / kbT, 1.0 / (n_H * lambda_e3));
}

/**
 * nc_thermodyn_HeII_ionization_saha_x:
 * @model: a #NcmModel, cosmological parameters for a given model
 * @frac: [; X_{He^{++}}X_e/X_{He^+} ;]
 *
 * @brief Calculate the redshift where the ratio #frac occur.
 *
 * This calculation is done using the saha equation as in \ref seager1999 .
 * Note that this function is valid only in the equilibrium. See also \ref recomb_sec_var.
 *
 * Returns: the value of [; x ;] where the ratio #frac occur.
 */
gdouble
nc_thermodyn_HeII_ionization_saha_x (NcHICosmo *model, gdouble frac)
{
  const gdouble T0 = nc_hicosmo_T_gamma0 (model);
  const gdouble kbT0 = NC_C_kb * T0;
  const gdouble mE_kbT0 = -NC_C_HeII_BINDING_ENERGY_1s / kbT0;
  const gdouble lambda_e3_0 = gsl_pow_3 (NC_C_THERMAL_WAVELENGTH_e/sqrt(T0));
  const gdouble h2 = nc_hicosmo_h2 (model);
  const gdouble Omega_b = nc_hicosmo_Omega_b (model);
  const gdouble n_H0 = NC_C_PRIM_H_FRAC * Omega_b * (NC_C_CRIT_DENSITY * h2) / NC_C_ENERGY_p;
  const gdouble y = 3.0 / 2.0 * gsl_sf_lambert_Wm1 (2.0 / 3.0 * mE_kbT0 * cbrt (gsl_pow_2 (lambda_e3_0 * n_H0 * frac))) / mE_kbT0;

  return (1.0 / y);
}

/**
 * nc_thermodyn_recomb_H_case_B:
 * @model: a #NcmModel, cosmological parameters for a given model
 * @Tm: the matter (baryons) temperature [; T_b ;]
 *
 * @brief The case B [; H^+ ;] recombination coefficient.
 *
 * The fitting formula of the case B recombination coefficient for [; H^+ ;] as
 * in \ref pequignot1991 "Pequignot".
 *
 * Returns: the value of the case B recombination coefficient for [; H^+ ;] [; \alpha_H ;] .
 */
gdouble
nc_thermodyn_recomb_H_case_B (NcHICosmo *model, gdouble Tm)
{
  const gdouble F =  1.14;   /* fudge factor */
  const gdouble G =  1e-19;
  const gdouble a =  4.309;
  const gdouble b = -0.6166;
  const gdouble c =  0.6703;
  const gdouble d =  0.5300;
  const gdouble t = Tm * 1e-4;
  const gdouble res = F * G * a * pow (t, b) / (1.0 + c * pow (t, d));
  return res;
}

/**
 * nc_thermodyn_recomb_H_case_B_dTm:
 * @model: a #NcmModel, cosmological parameters for a given model
 * @Tm: the matter (baryons) temperature [; T_b ;]
 *
 * @brief The case B [; H^+ ;] recombination coefficient derivative with respect to Tm.
 *
 * The derivative of the fitting formula of the case B recombination coefficient for [; H^+ ;]
 * #nc_thermodyn_recomb_H_case_B.
 *
 * Returns: the value of the case B recombination coefficient for [; H^+ ;] [; d\alpha_H/dT_b ;].
 */
gdouble
nc_thermodyn_recomb_H_case_B_dTm (NcHICosmo *model, gdouble Tm)
{
  const gdouble F =  1.14;   /* fudge factor */
  const gdouble G =  1e-19;
  const gdouble a =  4.309;
  const gdouble b = -0.6166;
  const gdouble c =  0.6703;
  const gdouble d =  0.5300;
  const gdouble t = Tm * 1e-4;
  const gdouble t_b = pow (t, b);
  const gdouble t_d = pow (t, d);
  const gdouble res = a * F * G * (b + c * (b - d) * t_d) * t_b / (Tm * gsl_pow_2 (1.0 + c * t_d));
  return res;
}

/**
 * nc_thermodyn_recomb_HeI_case_B:
 * @model: a #NcmModel, cosmological parameters for a given model
 * @Tm: the matter (baryons) temperature [; T_b ;]
 *
 * @brief The case B [; He^+ ;] recombination coefficient.
 *
 * The fitting formula of the case B recombination coefficient for [; He^+ ;] as
 * in \ref hummer1998 "Hummer and Storey".
 *
 * Returns: the value of the case B recombination coefficient for [; He^+ ;] [; \alpha_H ;] .
 */
gdouble
nc_thermodyn_recomb_HeI_case_B (NcHICosmo *model, gdouble Tm)
{
  const gdouble sqrt_Tm = sqrt(Tm);
  const gdouble sqrt_T1 = pow (10.0, 5.114 / 2.0);
  const gdouble sqrt_T2 = sqrt (3.0);
  const gdouble p = 0.711;
  const gdouble q = pow (10.0, -16.744);
  const gdouble sqrt_Tm_T2 = sqrt_Tm / sqrt_T2;
  const gdouble sqrt_Tm_T1 = sqrt_Tm / sqrt_T1;
  const gdouble res = q / (sqrt_Tm_T2 * pow (1.0 + sqrt_Tm_T2, 1.0 - p) * pow (1.0 + sqrt_Tm_T1, 1.0 + p));
  return res;
}

/**
 * nc_thermodyn_recomb_HeI_case_B_dTm:
 * @model: a #NcmModel, cosmological parameters for a given model
 * @Tm: the matter (baryons) temperature [; T_b ;]
 *
 * @brief The case B [; He^+ ;] recombination coefficient derivative with respect to Tm.
 *
 * The derivative of the fitting formula of the case B recombination coefficient for [; He^+ ;]
 * #nc_thermodyn_recomb_HeI_case_B.
 *
 * Returns: the value of the case B recombination coefficient for [; He^+ ;] [; d\alpha_H/dT_b ;].
 */
gdouble
nc_thermodyn_recomb_HeI_case_B_dTm (NcHICosmo *model, gdouble Tm)
{
  const gdouble sqrt_Tm = sqrt(Tm);
  const gdouble q = pow (10.0, -16.744);
  const gdouble T1 = pow (10.0, 5.114);
  const gdouble sqrt_T1 = sqrt (T1);
  const gdouble sqrt_T2 = sqrt (3.0);
  const gdouble p = 0.711;
  const gdouble T1_2 = T1 * T1;
  const gdouble sqrt_Tm_T2 = sqrt_Tm / sqrt_T2;
  const gdouble sqrt_Tm_T1 = sqrt_Tm / sqrt_T1;
  const gdouble Tm_T1_3_2 = gsl_pow_3 (sqrt_Tm_T1);
  const gdouble res = -q *
    (
      Tm * (2.0 + p + 3.0 * sqrt_Tm_T2) +
      T1 * sqrt_Tm_T1 * (1.0 + (2.0 - p) * sqrt_Tm_T2)
    ) /
    (
      2.0 * T1_2 * Tm_T1_3_2 * sqrt_Tm_T2 *
      pow (1.0 + sqrt_Tm_T2, 2.0 - p) *
      pow (1.0 + sqrt_Tm_T1, 2.0 + p)
      );
  return res;
}

/**
 * nc_thermodyn_H_ionization_rate_old:
 * @model: a #NcmModel, cosmological parameters for a given model
 * @Xp: FIXME
 * @Tm: FIXME
 * @XHeII: FIXME
 * @x: normalized scale factor inverse [; x = 1 + z = a_0/a ;]
 *
 * dX_e_dx implemented using Weinbergs book
 *
 * Returns: FIXME
 */
gdouble
nc_thermodyn_H_ionization_rate_old (NcHICosmo *model, gdouble Xp, gdouble Tm, gdouble XHeII, gdouble x)
{
  const gdouble Xe = Xp + XHeII;
  const gdouble alpha = nc_thermodyn_recomb_H_case_B (model, Tm);
  const gdouble h2 = nc_hicosmo_h2 (model);
  const gdouble Omega_b = nc_hicosmo_Omega_b (model);
  const gdouble n_b0 = Omega_b * NC_C_CRIT_NUMBER_DENSITY_p * h2;
  const gdouble n_0 = (1.0 - NC_C_PRIM_HE_Y_P) * n_b0;
  const gdouble H = nc_hicosmo_H (model, x - 1.0) / (NC_C_PARSEC * 1e3);
  const gdouble f1 = alpha * n_0 * x * x / H;

  const gdouble x3 = x*x*x;
  const gdouble Tm3_2 = sqrt(gsl_pow_3 (Tm));
  const gdouble f2na = NC_C_DECAY_H_2S_1S * (1.0 - Xp);
  const gdouble f2nb = H / x3 / (n_0 * NC_C_HYDROGEN_LYMAN_2p_WL3_8PI);
  const gdouble f2n = f2na + f2nb;
  const gdouble f2da = NC_BOLTZMAN_FACTOR_H_2s(Tm) * Tm3_2 * alpha * (1.0 - Xp);
  const gdouble f2d = f2na + f2nb + f2da;
  const gdouble f2 = f2n / f2d;

  const gdouble S = NC_C_BOLTZMAN_FACTOR_H_1s (Tm) * Tm3_2 / (n_0 * x3);
  const gdouble f3 = ( Xp * Xe - (1.0 - Xp) * S);

  return f1 * f2 * f3;
}

/**
 * nc_thermodyn_H_ionization_rate:
 * @model: a #NcmModel, cosmological parameters for a given model
 * @XH: FIXME
 * @Tm: FIXME
 * @XHeII: FIXME
 * @x: normalized scale factor inverse [; x = 1 + z = a_0/a ;]
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_thermodyn_H_ionization_rate (NcHICosmo *model, gdouble XH, gdouble Tm, gdouble XHeII, gdouble x)
{
  const gdouble Xe = XH + XHeII;
  const gdouble kbTm = NC_C_kb * Tm;
  const gdouble x3 = gsl_pow_3 (x);
  const gdouble h2 = nc_hicosmo_h2 (model);
  const gdouble Omega_b = nc_hicosmo_Omega_b (model);
  const gdouble Tm3_2 = sqrt(gsl_pow_3 (Tm));
  const gdouble lambda_e3_T3_2 = gsl_pow_3 (NC_C_THERMAL_WAVELENGTH_e);
  const gdouble n_H = NC_C_PRIM_H_FRAC * Omega_b * x3 * (NC_C_CRIT_DENSITY * h2) / NC_C_ENERGY_p;
  const gdouble alpha_H = nc_thermodyn_recomb_H_case_B (model, Tm);
  const gdouble beta_H = gsl_sf_exp_mult (-NC_C_HYDROGEN_BINDING_ENERGY_2s / kbTm, alpha_H * Tm3_2 / lambda_e3_T3_2);
  const gdouble beta_H_exp_mE_kbT = gsl_sf_exp_mult (-NC_C_HYDROGEN_LYMAN_2s / kbTm, beta_H);
  const gdouble f1 = (Xe * XH * n_H * alpha_H - beta_H_exp_mE_kbT * (1.0 - XH));

  const gdouble H = nc_hicosmo_H (model, x - 1.0) / NC_C_KILO_PARSEC;
  const gdouble KH = NC_C_HYDROGEN_LYMAN_2p_WL3_8PI / H;
  const gdouble Lambda_H = 8.22458;
  const gdouble f2 = 1.0 + KH * Lambda_H * n_H * (1.0 - XH);

  const gdouble f3 = H * x * (1.0 + KH * (Lambda_H + beta_H) * n_H * (1.0 - XH));

  return f1 * f2 / f3;
}

/**
 * nc_thermodyn_H_ionization_rate_grad: (skip)
 * @model: a #NcmModel, cosmological parameters for a given model
 * @Xp: FIXME
 * @Tm: FIXME
 * @XHeII: FIXME
 * @x: normalized scale factor inverse [; x = 1 + z = a_0/a ;]
 * @grad: FIXME
 *
 * FIXME
 * dX_e_dx implemented using Weinbergs book
 *
 */
void
nc_thermodyn_H_ionization_rate_grad (NcHICosmo *model, gdouble Xp, gdouble Tm, gdouble XHeII, gdouble x, gsl_vector *grad)
{
  const gdouble x3 = x*x*x;
  const gdouble Xe = Xp + XHeII;
  const gdouble Tm2 = Tm * Tm;
  const gdouble Tm3_2 = sqrt(Tm * Tm2);
  const gdouble alpha = nc_thermodyn_recomb_H_case_B (model, Tm);
  const gdouble beta = NC_BOLTZMAN_FACTOR_H_2s(Tm) * Tm3_2 * alpha;

  const gdouble h2 = nc_hicosmo_h2 (model);
  const gdouble Omega_b = nc_hicosmo_Omega_b (model);
  const gdouble n_b0 = Omega_b * NC_C_CRIT_NUMBER_DENSITY_p * h2;
  const gdouble n_0 = 0.76 * n_b0;
  const gdouble H = nc_hicosmo_H (model, x - 1.0) / (NC_C_PARSEC * 1e3);

  const gdouble f1 = alpha * n_0 * x * x / H;

  const gdouble f2na = NC_C_DECAY_H_2S_1S * (1.0 - Xp);
  const gdouble f2nb = H / x3 / (n_0 * NC_C_HYDROGEN_LYMAN_2p_WL3_8PI);
  const gdouble f2n = f2na + f2nb;

  const gdouble f2da = beta * (1.0 - Xp);
  const gdouble f2d = f2na + f2nb + f2da;

  const gdouble f2 = f2n / f2d;

  const gdouble S = NC_C_BOLTZMAN_FACTOR_H_1s(Tm) * Tm3_2 / (n_0 * x3);
  const gdouble f3 = ( Xp * Xe - (1.0 - Xp) * S);

  const gdouble ddX = f1 * (-NC_C_DECAY_H_2S_1S / f2d) * f3 +
    f1 * ((NC_C_DECAY_H_2S_1S + beta) * f2n / gsl_pow_2 (f2d)) * f3 +
    f1 * f2n / f2d * (Xe + Xp + S);

  const gdouble dalpha = nc_thermodyn_recomb_H_case_B_dTm (model, Tm);
  const gdouble df1 = f1 * dalpha / alpha;
  const gdouble dbeta =
    (3.0 / 2.0 / Tm + NC_C_HYDROGEN_BINDING_ENERGY_1s / NC_C_kb / Tm2 + dalpha / alpha) * beta;
  const gdouble df2 = - f2n * (dbeta * (1.0 - Xp)) / gsl_pow_2 (f2d);
  const gdouble dS =
    (3.0 / 2.0 / Tm + NC_C_HYDROGEN_BINDING_ENERGY_2s / NC_C_kb / Tm2) * S;
  const gdouble df3 = -(1.0 - Xp) * dS;
  const gdouble ddTm = df1 * f2 * f3 + f1 * df2 * f3 + f1 * f2 * df3;

  const gdouble ddXHeII = f1 * f2 * Xp;

  gsl_vector_set (grad, 0, ddX);
  gsl_vector_set (grad, 1, ddTm);
  gsl_vector_set (grad, 2, ddXHeII);

  return;
}

/**
 * nc_thermodyn_HeII_ionization_rate:
 * @model: a #NcmModel, cosmological parameters for a given model
 * @Xp: FIXME
 * @Tm: FIXME
 * @XHeII: FIXME
 * @x: normalized scale factor inverse [; x = 1 + z = a_0/a ;]
 *
 * dX_e_dx implemented using Weinbergs book
 *
 * Returns: FIXME
 */
gdouble
nc_thermodyn_HeII_ionization_rate (NcHICosmo *model, gdouble Xp, gdouble Tm, gdouble XHeII, gdouble x)
{
  const gdouble Xe = Xp + XHeII;
  const gdouble alpha = nc_thermodyn_recomb_HeI_case_B (model, Tm);
  const gdouble h2 = nc_hicosmo_h2 (model);
  const gdouble Omega_b = nc_hicosmo_Omega_b (model);
  const gdouble n_b0 = Omega_b * NC_C_CRIT_NUMBER_DENSITY_p * h2;
  const gdouble n_0 = 0.76 * n_b0;
  const gdouble H = nc_hicosmo_H (model, x - 1.0) / (NC_C_PARSEC * 1e3);
  const gdouble f1 = alpha * n_0 * x * x / H;

  const gdouble x3 = x*x*x;
  const gdouble Tm3_2 = sqrt(gsl_pow_3 (Tm));
  const gdouble f2na = NC_C_DECAY_He_2S_1S * (NC_C_PRIM_HE_XHe - XHeII);
  const gdouble f2nb = gsl_sf_exp_mult (-NC_C_HeI_2s_m_2p_Kb / Tm, H / x3 / (n_0 * NC_C_HeI_LYMAN_2p_WL3_8PI));
  const gdouble f2n = f2na + f2nb;
  const gdouble f2da = 4.0 * NC_BOLTZMAN_FACTOR_HeI_2s(Tm) * Tm3_2 * alpha * (NC_C_PRIM_HE_XHe - XHeII);
  const gdouble f2d = f2na + f2nb + f2da;
  const gdouble f2 = f2n / f2d;

  const gdouble S = 4.0 * NC_BOLTZMAN_FACTOR_HeI_1s(Tm) * Tm3_2 / (n_0 * x3);
  const gdouble f3 = ( XHeII * Xe - (NC_C_PRIM_HE_XHe - XHeII) * S);

  return f1 * f2 * f3;
}

/**
 * nc_thermodyn_HeII_ionization_rate_grad: (skip)
 * @model: a #NcmModel, cosmological parameters for a given model
 * @Xp: FIXME
 * @Tm: FIXME
 * @XHeII: FIXME
 * @x: normalized scale factor inverse [; x = 1 + z = a_0/a ;]
 * @grad: FIXME
 *
 * dX_e_dx implemented using Weinbergs book
 *
 */
void
nc_thermodyn_HeII_ionization_rate_grad (NcHICosmo *model, gdouble Xp, gdouble Tm, gdouble XHeII, gdouble x, gsl_vector *grad)
{
  const gdouble Xe = Xp + XHeII;
  const gdouble Tm2 = Tm * Tm;
  const gdouble alpha = nc_thermodyn_recomb_HeI_case_B (model, Tm);
  const gdouble h2 = nc_hicosmo_h2 (model);
  const gdouble Omega_b = nc_hicosmo_Omega_b (model);
  const gdouble n_b0 = Omega_b * NC_C_CRIT_NUMBER_DENSITY_p * h2;
  const gdouble n_0 = 0.76 * n_b0;
  const gdouble H = nc_hicosmo_H (model, x - 1.0) / (NC_C_PARSEC * 1e3);
  const gdouble f1 = alpha * n_0 * x * x / H;

  const gdouble x3 = x*x*x;
  const gdouble Tm3_2 = sqrt(gsl_pow_3 (Tm));
  const gdouble f2na = NC_C_DECAY_He_2S_1S * (NC_C_PRIM_HE_XHe - XHeII);
  const gdouble f2nb = gsl_sf_exp_mult (-NC_C_HeI_2s_m_2p_Kb / Tm, H / x3 / (n_0 * NC_C_HeI_LYMAN_2p_WL3_8PI));
  const gdouble f2n = f2na + f2nb;
  const gdouble beta = 4.0 * NC_BOLTZMAN_FACTOR_HeI_2s(Tm) * Tm3_2 * alpha;
  const gdouble f2da = beta * (NC_C_PRIM_HE_XHe - XHeII);
  const gdouble f2d = f2na + f2nb + f2da;
  const gdouble f2 = f2n / f2d;

  const gdouble S = 4.0 * NC_BOLTZMAN_FACTOR_HeI_1s(Tm) * Tm3_2 / (n_0 * x3);
  const gdouble f3 = ( XHeII * Xe - (NC_C_PRIM_HE_XHe - XHeII) * S);

  const gdouble ddXHeII = f1 * (-NC_C_DECAY_He_2S_1S / f2d) * f3 +
    f1 * ((NC_C_DECAY_He_2S_1S + beta) * f2n / gsl_pow_2 (f2d)) * f3 +
    f1 * f2n / f2d * (Xe + XHeII + S);

  const gdouble dalpha = nc_thermodyn_recomb_HeI_case_B_dTm (model, Tm);
  const gdouble df1 = f1 * dalpha / alpha;
  const gdouble dbeta = (3.0 / 2.0 / Tm + NC_C_HeI_BINDING_ENERGY_1s / NC_C_kb / Tm2 + dalpha / alpha) * beta;
  const gdouble df2nb = NC_C_HeI_2s_m_2p_Kb / Tm2 * f2nb;
  const gdouble df2 = df2nb / f2d - f2n * (dbeta * (NC_C_PRIM_HE_XHe - XHeII) + df2nb) / gsl_pow_2 (f2d);
  const gdouble dS =
    (3.0 / 2.0 / Tm + NC_C_HeI_BINDING_ENERGY_2s / NC_C_kb / Tm2) * S;
  const gdouble df3 = -(NC_C_PRIM_HE_XHe - XHeII) * dS;
  const gdouble ddTm = df1 * f2 * f3 + f1 * df2 * f3 + f1 * f2 * df3;

  const gdouble ddX = f1 * f2 * XHeII;

  gsl_vector_set (grad, 0, ddX);
  gsl_vector_set (grad, 1, ddTm);
  gsl_vector_set (grad, 2, ddXHeII);

}

/**
 * nc_thermodyn_matter_temperature_dx_old:
 * @model: a #NcmModel, cosmological parameters for a given model
 * @Xp: FIXME
 * @Tm: FIXME
 * @XHeII: FIXME
 * @x: normalized scale factor inverse [; x = 1 + z = a_0/a ;]
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_thermodyn_matter_temperature_dx_old (NcHICosmo *model, gdouble Xp, gdouble Tm, gdouble XHeII, gdouble x)
{
  const gdouble T0 = nc_hicosmo_T_gamma0 (model);
  const gdouble T = T0 * x;
  const gdouble H = nc_hicosmo_H (model, x - 1.0) / (NC_C_PARSEC * 1e3);
  const gdouble f1 = (8.0  * NC_C_THOMPSON_CS * NC_C_AR *
    T0 * T0 * T0 * T0 /
    (3.0 * NC_C_c * NC_C_MASS_e)) *
    gsl_pow_3 (x) / H;

  const gdouble Xe = Xp + XHeII;
  const gdouble f2 = Xe * (Tm - T) / (1.0 + NC_C_PRIM_HE_XHe + Xe);

  const gdouble f3 = 2.0 * Tm / x;

  return f1 * f2 + f3;
}

/**
 * nc_thermodyn_matter_temperature_dx:
 * @model: a #NcmModel, cosmological parameters for a given model
 * @Xp: FIXME
 * @Tm: FIXME
 * @XHeII: FIXME
 * normalized scale factor inverse [; x = 1 + z = a_0/a ;]
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_thermodyn_matter_temperature_dx (NcHICosmo *model, gdouble Xp, gdouble Tm, gdouble XHeII, gdouble x)
{
  const gdouble T0 = nc_hicosmo_T_gamma0 (model);
  const gdouble Tr = T0 * x;
  const gdouble H = nc_hicosmo_H (model, x - 1.0) / NC_C_KILO_PARSEC;

  const gdouble f1 = (8.0 * NC_C_THOMPSON_CS * NC_C_AR *
    T0 * T0 * T0 * T0 /
    (3.0 * NC_C_c * NC_C_MASS_e)) *
    gsl_pow_3 (x) / H;

  const gdouble Xe = Xp + XHeII;
  const gdouble f2 = Xe * (Tm - Tr) / (1.0 + NC_C_PRIM_HE_XHe + Xe);

  const gdouble f3 = 2.0 * Tm / x;

  return f1 * f2 + f3;
}

/**
 * nc_thermodyn_matter_temperature_dx_grad: (skip)
 * @model: a #NcmModel, cosmological parameters for a given model
 * @Xp: FIXME
 * @Tm: FIXME
 * @XHeII: FIXME
 * @x: normalized scale factor inverse [; x = 1 + z = a_0/a ;]
 * @grad: FIXME
 *
 * FIXME
 */
void
nc_thermodyn_matter_temperature_dx_grad (NcHICosmo *model, gdouble Xp, gdouble Tm, gdouble XHeII, gdouble x, gsl_vector *grad)
{
  const gdouble T0 = nc_hicosmo_T_gamma0 (model);
  const gdouble T = T0 * x;
  const gdouble H = nc_hicosmo_H (model, x - 1.0) / (NC_C_PARSEC * 1e3);
  const gdouble f1 = (8.0  * NC_C_THOMPSON_CS * NC_C_AR *
    T0 * T0 * T0 * T0 /
    (3.0 * NC_C_c * NC_C_MASS_e)) *
    gsl_pow_3 (x) / H;

  const gdouble Xe = Xp + XHeII;
  const gdouble f2 = Xe * (Tm - T) / (1.0 + NC_C_PRIM_HE_XHe + Xe);

  const gdouble df2_dX = (1.0 / Xe - 1.0 / (1.0 + NC_C_PRIM_HE_XHe + Xe)) * f2;
  const gdouble ddXp = f1 * df2_dX;

  const gdouble ddTm = f1 * Xe / (1.0 + NC_C_PRIM_HE_XHe + Xe) + 2.0 / x;

  const gdouble ddXHeII = ddXp;

  gsl_vector_set (grad, 0, ddXp);
  gsl_vector_set (grad, 1, ddTm);
  gsl_vector_set (grad, 2, ddXHeII);

  return;
}

static gint
H_ion_single_f (realtype x, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHICosmo *model =  NC_HICOSMO (f_data);
  const gdouble T0 = nc_hicosmo_T_gamma0 (model);
  const gdouble T = T0 * x;
  const gdouble XHeIIi = 0.0;

  NV_Ith_S(ydot, 0) = nc_thermodyn_H_ionization_rate (model, NV_Ith_S(y,0), T, XHeIIi, x);
  return 0;
}

static gint
H_ion_single_J (_NCM_SUNDIALS_INT_TYPE N, realtype x, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHICosmo *model = NC_HICOSMO (jac_data);
  const gdouble T0 = nc_hicosmo_T_gamma0 (model);
  const gdouble T = T0 * x;
  const gdouble XHeIIi = 0.0;
  gdouble grad_data[3];

  gsl_vector_view grad_view = gsl_vector_view_array (grad_data, 3);
  gsl_vector *grad = &grad_view.vector;
  nc_thermodyn_H_ionization_rate_grad (model, NV_Ith_S(y,0), T, XHeIIi, x, grad);

  DENSE_ELEM (J, 0, 0) = grad_data[0];

  return(0);
}

static gint
H_ion_single_Tm_f (realtype x, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHICosmo *model = NC_HICOSMO(f_data);
  gdouble Xp = NV_Ith_S(y,0);
  gdouble Tm = NV_Ith_S(y,1);
  const gdouble XHeIIi = 0.0;

  NV_Ith_S(ydot, 0) = nc_thermodyn_H_ionization_rate (model, Xp, Tm, XHeIIi, x);
  NV_Ith_S(ydot, 1) = nc_thermodyn_matter_temperature_dx (model, Xp, Tm, XHeIIi, x);

  return GSL_SUCCESS;
}

static gint
H_ion_single_Tm_J (_NCM_SUNDIALS_INT_TYPE N, realtype x, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHICosmo *model = NC_HICOSMO(jac_data);
  gdouble Xp = NV_Ith_S(y,0);
  gdouble Tm = NV_Ith_S(y,1);
  gdouble grad_data[3];
  const gdouble XHeIIi = 0.0;
  gsl_vector_view grad_view = gsl_vector_view_array (grad_data, 3);
  gsl_vector *grad = &grad_view.vector;

  nc_thermodyn_H_ionization_rate_grad (model, Xp, Tm, XHeIIi, x, grad);
  DENSE_ELEM (J,0,0) = grad_data[0];
  DENSE_ELEM (J,0,1) = grad_data[1];

  nc_thermodyn_matter_temperature_dx_grad (model, Xp, Tm, XHeIIi, x, grad);
  DENSE_ELEM (J,1,0) = grad_data[0];
  DENSE_ELEM (J,1,1) = grad_data[1];

  return(0);
}

static gint
H_ion_full_f (realtype x, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHICosmo *model = NC_HICOSMO(f_data);
  const gdouble Xp = NV_Ith_S(y, 0);
  const gdouble Tm = NV_Ith_S(y, 1);
  const gdouble XHeII = NV_Ith_S(y, 2);
  NV_Ith_S(ydot, 0) = nc_thermodyn_H_ionization_rate (model, Xp, Tm, XHeII, x);
  NV_Ith_S(ydot, 1) = nc_thermodyn_matter_temperature_dx (model, Xp, Tm, XHeII, x);
  NV_Ith_S(ydot, 2) = nc_thermodyn_HeII_ionization_rate (model, Xp, Tm, XHeII, x);

  return GSL_SUCCESS;
}

static gint
H_ion_full_J (_NCM_SUNDIALS_INT_TYPE N, realtype x, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHICosmo *model = NC_HICOSMO(jac_data);
  gdouble Xp = NV_Ith_S(y,0);
  gdouble Tm = NV_Ith_S(y,1);
  gdouble XHeII = NV_Ith_S(y,2);
  gdouble grad_data[3];
  gsl_vector_view grad_view = gsl_vector_view_array (grad_data, 3);
  gsl_vector *grad = &grad_view.vector;

  nc_thermodyn_H_ionization_rate_grad (model, Xp, Tm, XHeII, x, grad);
  DENSE_ELEM (J,0,0) = grad_data[0];
  DENSE_ELEM (J,0,1) = grad_data[1];
  DENSE_ELEM (J,0,2) = grad_data[2];

  nc_thermodyn_matter_temperature_dx_grad (model, Xp, Tm, XHeII, x, grad);
  DENSE_ELEM (J,1,0) = grad_data[0];
  DENSE_ELEM (J,1,1) = grad_data[1];
  DENSE_ELEM (J,1,2) = grad_data[2];

  nc_thermodyn_HeII_ionization_rate_grad (model, Xp, Tm, XHeII, x, grad);
  DENSE_ELEM (J,2,0) = grad_data[0];
  DENSE_ELEM (J,2,1) = grad_data[1];
  DENSE_ELEM (J,2,2) = grad_data[2];

  return(0);
}

/**
 * nc_thermodyn_recomb_new: (skip)
 * @model: a #NcmModel, cosmological parameters for a given model
 * @type: a #NcThermodynRecombType
 *
 * FIXME
 */
NcThermodynRecomb *
nc_thermodyn_recomb_new (NcHICosmo *model, NcThermodynRecombType type)
{
  NcThermodynRecomb *recomb = g_slice_new (NcThermodynRecomb);
  gint flag;

  recomb->model = model;

  recomb->dtau_dx_accel = gsl_interp_accel_alloc ();
  recomb->tau_accel = gsl_interp_accel_alloc ();
  recomb->dtau_dx_spline = NULL;
  recomb->dtau_dxR_spline = NULL;
  recomb->init_spline = TRUE;

  recomb->H_ion_fraction = nc_function_cache_new (3, NC_INT_ABS_ERROR, NC_INT_ERROR);
  recomb->optical_depth = nc_function_cache_new (1, NC_INT_ABS_ERROR, NC_INT_ERROR);

  switch (type)
  {
    case NC_THERMODYN_RECOMB_SINGLE:
      recomb->n = 1;
      recomb->ion = &H_ion_single_f;
      recomb->ion_J = &H_ion_single_J;
      break;
    case NC_THERMODYN_RECOMB_SINGLE_Tm:
      recomb->n = 2;
      recomb->ion = &H_ion_single_Tm_f;
      recomb->ion_J = &H_ion_single_Tm_J;
      break;
    case NC_THERMODYN_RECOMB_FULL:
      recomb->n = 3;
      recomb->ion = &H_ion_full_f;
      recomb->ion_J = &H_ion_full_J;
      break;
  }

  recomb->y0 = N_VNew_Serial(recomb->n);
  recomb->y = N_VNew_Serial(recomb->n);

  {
    /* Assuming XHI = 1 and XHe_n = 0 obtain the ratio [; X_{He^{++}}X_e/X_{He^+} ;] */
    const gdouble frac  = 1e-11;
    const gdouble ratio = frac * (1.0 + NC_C_PRIM_HE_XHe * (1.0 + frac)) / ((1.0 - frac));
    recomb->He_single_ionized_x = nc_thermodyn_HeII_ionization_saha_x (model, ratio);
  }

  recomb->reltol = 1e-11;
  recomb->abstol = 1e-20;

  recomb->cvode = CVodeCreate (CV_BDF, CV_NEWTON);
  CVODE_CHECK((void *)recomb->cvode, "CVodeCreate", 0, NULL);

  {
    const gdouble T0 = nc_hicosmo_T_gamma0 (model);
    gdouble XH, XHeII, Tm, XeXHeII_XHeI;
    gdouble x0 = recomb->He_single_ionized_x;

    Tm = T0 * x0;
    XH = 1; /* Hydrogen is completly ionized */
    XeXHeII_XHeI = nc_thermodyn_HeI_ionization_saha (model, x0);
    XHeII = (XeXHeII_XHeI * ncm_sqrt1px_m1 ( (1.0 + 2.0 * XeXHeII_XHeI + 4.0 * NC_C_PRIM_HE_XHe * XeXHeII_XHeI) / (XeXHeII_XHeI * XeXHeII_XHeI) ) - 1.0) / 2.0;

    recomb->x0 = x0;
    NV_Ith_S(recomb->y0, 0) = XH;
    if (recomb->n > 1)
      NV_Ith_S(recomb->y0, 1) = Tm;
    if (recomb->n > 2)
      NV_Ith_S(recomb->y0, 2) = XHeII;
  }

  flag = CVodeInit (recomb->cvode, recomb->ion, recomb->x0, recomb->y0);
  CVODE_CHECK(&flag, "CVodeMalloc", 1, FALSE);
  recomb->init = TRUE;

  nc_thermodyn_recomb_reset (recomb);

  return recomb;
}

/**
 * nc_thermodyn_recomb_free:
 * @recomb: a #NcThermodynRecomb
 *
 * FIXME
*/
void
nc_thermodyn_recomb_free (NcThermodynRecomb *recomb)
{
	/* FIXME LEAK!!! */
	nc_function_cache_free (recomb->H_ion_fraction);
	nc_function_cache_free (recomb->optical_depth);
}

/**
 * nc_thermodyn_recomb_reinit:
 * @recomb: a #NcThermodynRecomb
 *
 * FIXME
 *
 * Returns: FIXME
*/
gboolean
nc_thermodyn_recomb_reinit (NcThermodynRecomb *recomb)
{
  gint flag;

  flag = CVodeReInit (recomb->cvode, recomb->x0, recomb->y0);
  CVODE_CHECK(&flag, "CVodeReInit", 1, FALSE);

  nc_thermodyn_recomb_reset (recomb);
  return TRUE;
}

/**
 * nc_thermodyn_recomb_reset:
 * @recomb: a #NcThermodynRecomb
 *
 * FIXME
 *
 * Returns: FIXME
*/
gboolean
nc_thermodyn_recomb_reset (NcThermodynRecomb *recomb)
{
  gint flag;

  recomb->H_ion_fraction->clear = TRUE;

  flag = CVodeSStolerances (recomb->cvode, recomb->reltol, recomb->abstol);
  CVODE_CHECK(&flag, "CVodeSStolerances", 1, FALSE);

  flag = CVodeSetUserData (recomb->cvode , recomb->model);
  CVODE_CHECK(&flag, "CVodeSetUserData", 1, FALSE);

  flag = CVodeSetMaxNumSteps (recomb->cvode, 50000);
  CVODE_CHECK(&flag, "CVodeSetMaxNumSteps", 1, FALSE);

  flag = CVDense (recomb->cvode, recomb->n);
  CVODE_CHECK(&flag, "CVDense", 1, FALSE);

  flag = CVDlsSetDenseJacFn (recomb->cvode, recomb->ion_J);
  CVODE_CHECK(&flag, "CVDlsSetDenseJacFn", 1, FALSE);

  flag = CVodeSetStopTime (recomb->cvode, X_TODAY);
  CVODE_CHECK(&flag, "CVodeSetStopTime", 1, FALSE);

  return TRUE;
}

/**
 * nc_thermodyn_recomb_evolve:
 * @recomb: a #NcThermodynRecomb
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gboolean
nc_thermodyn_recomb_evolve (NcThermodynRecomb *recomb, gdouble x)
{
  gsl_vector *pf;
  gdouble x_found;
  gdouble xi;
  gint flag;

  if (x == recomb->He_single_ionized_x)
    return TRUE;

  if (nc_function_cache_get_near (recomb->H_ion_fraction, x, &x_found, &pf, NC_FUNCTION_CACHE_SEARCH_GT))
  {
    NV_Ith_S(recomb->y, 0) = gsl_vector_get(pf, 0);
    NV_Ith_S(recomb->y, 1) = gsl_vector_get(pf, 1);
    NV_Ith_S(recomb->y, 2) = gsl_vector_get(pf, 2);
    flag = CVodeReInit (recomb->cvode, x_found, recomb->y);
    CVODE_CHECK(&flag, "CVodeReInit", 1, FALSE);
    nc_thermodyn_recomb_reset (recomb);
    if (x_found == x)
      return TRUE;
  }
  else
    nc_thermodyn_recomb_reinit (recomb);

  pf = gsl_vector_alloc (3);
  flag = CVode(recomb->cvode, x, recomb->y, &xi, CV_NORMAL);
  CVODE_CHECK(&flag, "CVode", 1, FALSE);

  recomb->He_single_ionized_x = xi;

  gsl_vector_set (pf, 0, NV_Ith_S(recomb->y,0));
  gsl_vector_set (pf, 1, NV_Ith_S(recomb->y,1));
  gsl_vector_set (pf, 2, NV_Ith_S(recomb->y,2));
  nc_function_cache_insert_vector (recomb->H_ion_fraction, x, pf);

  return TRUE;
}

/**
 * nc_thermodyn_recomb_get_Xe:
 * @recomb: a #NcThermodynRecomb
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_thermodyn_recomb_get_Xe (NcThermodynRecomb *recomb)
{
  return NV_Ith_S (recomb->y, 0) + NV_Ith_S (recomb->y, 2);
}

/**
 * nc_thermodyn_recomb_get_Xe_at:
 * @recomb: a #NcThermodynRecomb
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_thermodyn_recomb_get_Xe_at (NcThermodynRecomb *recomb, gdouble x)
{
  if (x < recomb->x0)
  {
    nc_thermodyn_recomb_evolve (recomb, x);
    return NV_Ith_S (recomb->y, 0) + NV_Ith_S (recomb->y, 2);
  }
  else
  {
    /*
     * Above this value recomb->x0 all helium is single or double ionized,
     * and all hydrogen is ionized.
     */
    const gdouble XHeIIIXe_XHeII = nc_thermodyn_HeII_ionization_saha (recomb->model, x);
    const gdouble arg = NC_C_PRIM_HE_XHe * (NC_C_PRIM_HE_XHe + (2.0 + 6.0 * XHeIIIXe_XHeII)) / ((1.0 + XHeIIIXe_XHeII) * (1.0 + XHeIIIXe_XHeII));
    const gdouble Xe = (2.0 + NC_C_PRIM_HE_XHe + (1.0 + XHeIIIXe_XHeII) * ncm_sqrt1px_m1 (arg)) / 2.0;
    return Xe;
  }
}

static gdouble
dopt_depth (gdouble x, gpointer p)
{
  NcThermodynRecomb *recomb = (NcThermodynRecomb *)p;
  NcHICosmo *model = recomb->model;
  const gdouble x2 = x * x;
  const gdouble h2 = nc_hicosmo_h2 (model);
  const gdouble Omega_b = nc_hicosmo_Omega_b (model);
  const gdouble n_b0 = Omega_b * NC_C_CRIT_NUMBER_DENSITY_p * h2;
  const gdouble n_0 = 0.76 * n_b0;
  const gdouble H = nc_hicosmo_H (model, x - 1.0) / (NC_C_KILO_PARSEC);
  const gdouble Xe = nc_thermodyn_recomb_get_Xe_at (recomb, x);
//printf ("# %.15e %.15e\n", x, Xe);
  return NC_C_c * NC_C_THOMPSON_CS * n_0 * Xe * x2 / H;
}

static gdouble
dopt_depth_over_Xe (NcHICosmo *model, gdouble x)
{
  const gdouble x2 = x * x;
  const gdouble h2 = nc_hicosmo_h2 (model);
  const gdouble Omega_b = nc_hicosmo_Omega_b (model);
  const gdouble n_b0 = Omega_b * NC_C_CRIT_NUMBER_DENSITY_p * h2;
  const gdouble n_0 = 0.76 * n_b0;
  const gdouble H = nc_hicosmo_H (model, x - 1.0) / (NC_C_PARSEC * 1e3);

  return NC_C_c * NC_C_THOMPSON_CS * n_0 * x2 / H;
}

static void _nc_thermodyn_recomb_calc_peak_width (NcThermodynRecomb *recomb);

/**
 * nc_thermodyn_recomb_dtau_dx_init_spline:
 * @recomb: a #NcThermodynRecomb
 *
 * FIXME
 *
*/
void
nc_thermodyn_recomb_dtau_dx_init_spline (NcThermodynRecomb *recomb)
{
#define _NC_PRERECOMB_PARTITION 2000.0
  if (recomb->init_spline)
  {
    gdouble xi, mx, xini = 1.1 * NC_PERTURBATION_START_X;
    gdouble Xe, dtau_dx, dtau_dxR;
    GArray *x = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), 2000);
    GArray *y = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), 2000);
    GArray *yR = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), 2000);
    gdouble exp_arg = log(recomb->x0/xini) / _NC_PRERECOMB_PARTITION;
    const gdouble Omega_r = nc_hicosmo_Omega_r (NC_HICOSMO (recomb->model));
    const gdouble Omega_b = nc_hicosmo_Omega_b (NC_HICOSMO (recomb->model));
    const gdouble R0 = 4.0 * Omega_r / (3.0 * Omega_b);
    gint i;

    for (i = 0; i < _NC_PRERECOMB_PARTITION; i++)
    {
      mx = -xini * exp (i * exp_arg);
      dtau_dx = dopt_depth (-mx, recomb);
      dtau_dxR = -dtau_dx * R0 * mx;
//printf ("%.15g %.15g\n", -1/mx, nc_thermodyn_recomb_get_Xe_at (recomb, -mx));
      g_array_append_val (x, mx);
      g_array_append_val (y, dtau_dx);
      g_array_append_val (yR, dtau_dxR);
    }

//    printf ("# Recomb: Changing to ode at %.15g\n", recomb->x0);
    Xe = NV_Ith_S (recomb->y0, 0) + NV_Ith_S (recomb->y0, 2);
    dtau_dx = dopt_depth_over_Xe (recomb->model, recomb->x0) * Xe;
    mx = -recomb->x0;
    dtau_dxR = -dtau_dx * R0 * mx;

//printf ("%.15g %.15g\n", -1/mx, nc_thermodyn_recomb_get_Xe_at (recomb, -mx));

    g_array_append_val (x, mx);
    g_array_append_val (y, dtau_dx);
    g_array_append_val (yR, dtau_dxR);

    nc_thermodyn_recomb_reinit (recomb);

    while (1)
    {
      CVode (recomb->cvode, 1.0, recomb->y, &xi, CV_ONE_STEP);
      Xe = NV_Ith_S (recomb->y, 0) + NV_Ith_S (recomb->y, 2);
//printf ("%.15g %.15g %.15g %.15g\n", -1/mx, NV_Ith_S (recomb->y, 0), NV_Ith_S (recomb->y, 1), NV_Ith_S (recomb->y, 2));
      dtau_dx = dopt_depth_over_Xe (recomb->model, xi) * Xe;
      if (fabs((mx+xi)/mx) > 1e-7)
      {
        mx = -xi;
        dtau_dxR = dtau_dx * R0 * xi;
        g_array_append_val (x, mx);
        g_array_append_val (y, dtau_dx);
        g_array_append_val (yR, dtau_dxR);
      }
      mx = -xi;
      if (xi == 1.0)
        break;
    }

    recomb->tau_spline = gsl_spline_alloc (gsl_interp_cspline, x->len);
    recomb->log_g_spline = gsl_spline_alloc (gsl_interp_cspline, x->len);
    recomb->dtau_dx_spline = gsl_spline_alloc (gsl_interp_cspline, x->len);
    recomb->dtau_dxR_spline = gsl_spline_alloc (gsl_interp_cspline, x->len);

    gsl_spline_init (recomb->dtau_dx_spline,
                     &(g_array_index(x, gdouble, 0)),
                     &(g_array_index(y, gdouble, 0)),
                     x->len);

    gsl_spline_init (recomb->dtau_dxR_spline,
                     &(g_array_index(x, gdouble, 0)),
                     &(g_array_index(yR, gdouble, 0)),
                     x->len);

    {
      gdouble tau = 0.0;
      gdouble last_x = g_array_index(x, gdouble, x->len - 1);

      g_array_index(y, gdouble, x->len - 1) = 0.0;

      for (i = x->len - 2; i >= 0; i--)
      {
        const gdouble E = sqrt(nc_hicosmo_E2 (NC_HICOSMO (recomb->model), -last_x - 1.0));
        g_array_index (yR, gdouble, i + 1) = -tau + log (fabs(E * gsl_spline_eval(recomb->dtau_dx_spline, last_x, recomb->dtau_dx_accel)));
        tau += gsl_spline_eval_integ (recomb->dtau_dx_spline, g_array_index(x, gdouble, i), last_x, recomb->dtau_dx_accel);
        g_array_index (y, gdouble, i) = tau;
        last_x = g_array_index(x, gdouble, i);
      }
      g_array_index (yR, gdouble, 0) = -tau + log (fabs(gsl_spline_eval(recomb->dtau_dx_spline, last_x, recomb->dtau_dx_accel)));
    }

    gsl_spline_init (recomb->tau_spline,
                     &(g_array_index(x, gdouble, 0)),
                     &(g_array_index(y, gdouble, 0)),
                     x->len);

    gsl_spline_init (recomb->log_g_spline,
                     &(g_array_index(x, gdouble, 0)),
                     &(g_array_index(yR, gdouble, 0)),
                     x->len);

    _nc_thermodyn_recomb_calc_peak_width (recomb);

    g_array_free (x, TRUE);
    g_array_free (y, TRUE);
    g_array_free (yR, TRUE);

    recomb->init_spline = FALSE;
  }
}

typedef struct __nc_thermodyn_recomb_log_g_data
{
  NcThermodynRecomb *recomb;
  gdouble ref;
} _nc_thermodyn_recomb_log_g_data;

static gdouble
_nc_thermodyn_recomb_log_g (gdouble x, gpointer params)
{
  _nc_thermodyn_recomb_log_g_data *log_g_data = (_nc_thermodyn_recomb_log_g_data *) params;
  NcThermodynRecomb *recomb = log_g_data->recomb;
  const gdouble log_g = -gsl_spline_eval (recomb->log_g_spline, -x, recomb->tau_accel);
  return log_g;
}

static gdouble
_nc_thermodyn_recomb_log_g_g_ref_2 (gdouble x, gpointer params)
{
  _nc_thermodyn_recomb_log_g_data *log_g_data = (_nc_thermodyn_recomb_log_g_data *) params;
  NcThermodynRecomb *recomb = log_g_data->recomb;
  const gdouble log_g = -gsl_spline_eval (recomb->log_g_spline, -x, recomb->tau_accel);
  const gdouble log_g_g_ref = log_g - log_g_data->ref;
  return log_g_g_ref * log_g_g_ref;
}

static gdouble
_nc_thermodyn_recomb_tau_m_tauref_2 (gdouble x, gpointer params)
{
  _nc_thermodyn_recomb_log_g_data *log_g_data = (_nc_thermodyn_recomb_log_g_data *) params;
  NcThermodynRecomb *recomb = log_g_data->recomb;
  const gdouble tau = gsl_spline_eval (recomb->tau_spline, -x, recomb->tau_accel);
  const gdouble tau_m_tau_ref = tau - log_g_data->ref;
  return tau_m_tau_ref * tau_m_tau_ref;
}

static void
_nc_thermodyn_recomb_calc_peak_width_iterate (gsl_min_fminimizer *s, gsl_function *F, glong max_iter, gdouble *xm, gdouble *xi, gdouble *xf)
{
  glong iter = 0;
  gint status;
  gsl_min_fminimizer_set (s, F, *xm, *xi, *xf);
  do
  {
    iter++;
    gsl_min_fminimizer_iterate (s);
    *xm = gsl_min_fminimizer_x_minimum (s);
    *xi = gsl_min_fminimizer_x_lower (s);
    *xf = gsl_min_fminimizer_x_upper (s);
    status = gsl_min_test_interval (*xi, *xf, 0.0, 1e-7);
  }
  while (status == GSL_CONTINUE && iter < max_iter);
}


static void
_nc_thermodyn_recomb_calc_peak_width (NcThermodynRecomb *recomb)
{
  const gsl_min_fminimizer_type *T = gsl_min_fminimizer_brent;
  gsl_min_fminimizer *s = gsl_min_fminimizer_alloc (T);
  gdouble xi = 1.0;
  gdouble xf = recomb->x0;
  gdouble xm = 1090.0; /* Just a starting value */
  gdouble min_log_g, log_gi;
  guint max_iter = 1000;
  _nc_thermodyn_recomb_log_g_data log_g_data;
  gsl_function F;

  log_g_data.recomb = recomb;
  log_g_data.ref = 0.0;

  F.function = &_nc_thermodyn_recomb_log_g;
  F.params = &log_g_data;

  _nc_thermodyn_recomb_calc_peak_width_iterate (s, &F, max_iter, &xm, &xi, &xf);
  recomb->x_rec = xm;
  min_log_g = (gsl_min_fminimizer_f_minimum (s));
  log_gi = (_nc_thermodyn_recomb_log_g (1.0, &log_g_data));
  /* log_gf = (_nc_thermodyn_recomb_log_g (recomb->x0, &log_g_data)); */

  F.function = &_nc_thermodyn_recomb_log_g_g_ref_2;
  log_g_data.ref = GSL_MIN (min_log_g + 2.0 * M_LN10, (log_gi + min_log_g) * 0.5);
  xm = (1.0 + recomb->x_rec) / 2.0;
  xi = 1.0;
  xf = recomb->x_rec;
  _nc_thermodyn_recomb_calc_peak_width_iterate (s, &F, max_iter, &xm, &xi, &xf);
  recomb->x_rec_10m2_max[0] = xm;

  xm = (recomb->x0 + recomb->x_rec) / 2.0;
  xi = recomb->x_rec;
  xf = recomb->x0;
  _nc_thermodyn_recomb_calc_peak_width_iterate (s, &F, max_iter, &xm, &xi, &xf);
  recomb->x_rec_10m2_max[1] = xm;

  F.function = &_nc_thermodyn_recomb_tau_m_tauref_2;
  log_g_data.ref = -GSL_LOG_DBL_EPSILON;
  xm = (recomb->x0 + recomb->x_rec) / 2.0;
  xi = recomb->x_rec;
  xf = recomb->x0;
  _nc_thermodyn_recomb_calc_peak_width_iterate (s, &F, max_iter, &xm, &xi, &xf);
  recomb->x_opt_cutoff = xm;

  gsl_min_fminimizer_free (s);
}

/**
 * nc_thermodyn_recomb_dtau_dx:
 * @recomb: a #NcThermodynRecomb
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_thermodyn_recomb_dtau_dx (NcThermodynRecomb *recomb, gdouble x)
{
  nc_thermodyn_recomb_dtau_dx_init_spline (recomb);
  return gsl_spline_eval (recomb->dtau_dx_spline, -x, recomb->dtau_dx_accel);
}

/**
 * nc_thermodyn_recomb_taubar:
 * @recomb: a #NcThermodynRecomb
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_thermodyn_recomb_taubar (NcThermodynRecomb *recomb, gdouble x)
{
  nc_thermodyn_recomb_dtau_dx_init_spline (recomb);
  return -x * gsl_spline_eval (recomb->dtau_dx_spline, -x, recomb->dtau_dx_accel);
}

/**
 * nc_thermodyn_recomb_taubbar:
 * @recomb: a #NcThermodynRecomb
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_thermodyn_recomb_taubbar (NcThermodynRecomb *recomb, gdouble x)
{
  nc_thermodyn_recomb_dtau_dx_init_spline (recomb);
  {
    const gdouble dtau_dx = nc_thermodyn_recomb_dtau_dx (recomb, x);
    const gdouble d2tau_dx2 = -gsl_spline_eval_deriv (recomb->dtau_dx_spline, -x, recomb->dtau_dx_accel);

    return x * dtau_dx + x * x * d2tau_dx2;
  }
}

/**
 * nc_thermodyn_recomb_taubbar_val:
 * @recomb: a #NcThermodynRecomb
 * @x: FIXME
 * @taubar: FIXME
 * @taubbar: FIXME
 *
 * FIXME
 *
*/
void
nc_thermodyn_recomb_taubbar_val (NcThermodynRecomb *recomb, gdouble x, gdouble *taubar, gdouble *taubbar)
{
  nc_thermodyn_recomb_dtau_dx_init_spline (recomb);
  {
    const gdouble dtau_dx = nc_thermodyn_recomb_dtau_dx (recomb, x);
    const gdouble d2tau_dx2 = -gsl_spline_eval_deriv (recomb->dtau_dx_spline, -x, recomb->dtau_dx_accel);

    *taubar = -x * dtau_dx;
    *taubbar = x * dtau_dx + x * x * d2tau_dx2;
  }
}

/**
 * nc_thermodyn_recomb_taubbbar:
 * @recomb: a #NcThermodynRecomb
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_thermodyn_recomb_taubbbar (NcThermodynRecomb *recomb, gdouble x)
{
  nc_thermodyn_recomb_dtau_dx_init_spline (recomb);
  {
    const gdouble dtau_dx = nc_thermodyn_recomb_dtau_dx (recomb, x);
    const gdouble d2tau_dx2 = -gsl_spline_eval_deriv (recomb->dtau_dx_spline, -x, recomb->dtau_dx_accel);
    const gdouble d3tau_dx3 = gsl_spline_eval_deriv2 (recomb->dtau_dx_spline, -x, recomb->dtau_dx_accel);
    const gdouble x2 = x * x;
    const gdouble x3 = x2 * x;

    return -(x * dtau_dx + 3.0 * x2 * d2tau_dx2 + x3 * d3tau_dx3);
  }
}

/**
 * nc_thermodyn_recomb_optical_depth:
 * @recomb: a #NcThermodynRecomb
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_thermodyn_recomb_optical_depth (NcThermodynRecomb *recomb, gdouble x)
{
  nc_thermodyn_recomb_dtau_dx_init_spline (recomb);

//  return gsl_spline_eval_integ (recomb->dtau_dx_spline, -x, -1.0, recomb->dtau_dx_accel);
  return gsl_spline_eval (recomb->tau_spline, -x, recomb->dtau_dx_accel);
}

/**
 * nc_thermodyn_recomb_optical_depth_x0_x1:
 * @recomb: a #NcThermodynRecomb
 * @x0: FIXME
 * @x1: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_thermodyn_recomb_optical_depth_x0_x1 (NcThermodynRecomb *recomb, gdouble x0, gdouble x1)
{
  gdouble result;
  gint ret;
  g_assert (x1 >= x0);
  nc_thermodyn_recomb_dtau_dx_init_spline (recomb);
  ret = gsl_spline_eval_integ_e (recomb->dtau_dx_spline, -x1, -x0, recomb->dtau_dx_accel, &result);
  if (ret != GSL_SUCCESS)
    g_error ("%s", gsl_strerror (ret));

  return -result;
}

/**
 * nc_thermodyn_recomb_optical_depth_R_x0_x1:
 * @recomb: a #NcThermodynRecomb
 * @x0: FIXME
 * @x1: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_thermodyn_recomb_optical_depth_R_x0_x1 (NcThermodynRecomb *recomb, gdouble x0, gdouble x1)
{
  gdouble result;
  gint ret;
  g_assert (x1 >= x0);
  nc_thermodyn_recomb_dtau_dx_init_spline (recomb);
  ret = gsl_spline_eval_integ_e (recomb->dtau_dxR_spline, -x1, -x0, recomb->dtau_dx_accel, &result);
  if (ret != GSL_SUCCESS)
    g_error ("%s", gsl_strerror (ret));

  return -result;
}

/**
 * nc_thermodyn_recomb_log_g:
 * @recomb: a #NcThermodynRecomb
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_thermodyn_recomb_log_g (NcThermodynRecomb *recomb, gdouble x)
{
  return gsl_spline_eval (recomb->log_g_spline, -x, recomb->tau_accel);
}

/**
 * nc_thermodyn_recomb_g:
 * @recomb: a #NcThermodynRecomb
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_thermodyn_recomb_g (NcThermodynRecomb *recomb, gdouble x)
{
  gdouble tau = nc_thermodyn_recomb_optical_depth (recomb, x);
  gdouble dtau_dx = nc_thermodyn_recomb_dtau_dx (recomb, x);

  return x * dtau_dx * exp(-tau);
}

/**
 * nc_thermodyn_recomb_gbar:
 * @recomb: a #NcThermodynRecomb
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_thermodyn_recomb_gbar (NcThermodynRecomb *recomb, gdouble x)
{
  const gdouble tau = nc_thermodyn_recomb_optical_depth (recomb, x);
  const gdouble dtau_dx = nc_thermodyn_recomb_dtau_dx (recomb, x);
  const gdouble d2tau_dx2 = -gsl_spline_eval_deriv (recomb->dtau_dx_spline, -x, recomb->dtau_dx_accel);
  const gdouble x2 = x * x;

  return (x2 * gsl_pow_2 (dtau_dx) - (x * dtau_dx + x2 * d2tau_dx2)) * exp(-tau);
}

/**
 * nc_thermodyn_recomb_gbbar:
 * @recomb: a #NcThermodynRecomb
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_thermodyn_recomb_gbbar (NcThermodynRecomb *recomb, gdouble x)
{
  const gdouble tau = nc_thermodyn_recomb_optical_depth (recomb, x);
  const gdouble dtau_dx = nc_thermodyn_recomb_dtau_dx (recomb, x);
  const gdouble d2tau_dx2 = -gsl_spline_eval_deriv (recomb->dtau_dx_spline, -x, recomb->dtau_dx_accel);
  const gdouble d3tau_dx3 = gsl_spline_eval_deriv2 (recomb->dtau_dx_spline, -x, recomb->dtau_dx_accel);
  const gdouble x2 = x * x;
  const gdouble x3 = x2 * x;

  return (
    x3 * gsl_pow_3 (dtau_dx)
    + (x * dtau_dx + 3.0 * x2 * d2tau_dx2 + x3 * d3tau_dx3 )
    - 3.0 * dtau_dx * (x2 * dtau_dx + x3 * d2tau_dx2)
    ) * exp(-tau);
}

static gdouble
gbar_exp_tau_f (gdouble x, gpointer params)
{
  NcThermodynRecomb *recomb = (NcThermodynRecomb *) params;
  const gdouble dtau_dx = nc_thermodyn_recomb_dtau_dx (recomb, x);
  const gdouble d2tau_dx2 = -gsl_spline_eval_deriv (recomb->dtau_dx_spline, -x, recomb->dtau_dx_accel);
  const gdouble x2 = x * x;

  return (x2 * gsl_pow_2 (dtau_dx) - (x * dtau_dx + x2 * d2tau_dx2));
}

static gdouble
gbar_exp_tau_df (gdouble x, gpointer params)
{
  NcThermodynRecomb *recomb = (NcThermodynRecomb *) params;
  const gdouble dtau_dx = nc_thermodyn_recomb_dtau_dx (recomb, x);
  const gdouble d2tau_dx2 = -gsl_spline_eval_deriv (recomb->dtau_dx_spline, -x, recomb->dtau_dx_accel);
  const gdouble d3tau_dx3 = -gsl_spline_eval_deriv2 (recomb->dtau_dx_spline, -x, recomb->dtau_dx_accel);
  const gdouble x2 = x * x;

  return 2.0 * x * dtau_dx * (dtau_dx + x * d2tau_dx2) - dtau_dx - 3.0 * x * d2tau_dx2 - x2 * d3tau_dx3;
}

static void
gbar_exp_tau_fdf (gdouble x, gpointer params, gdouble *y, gdouble *dy)
{
  NcThermodynRecomb *recomb = (NcThermodynRecomb *) params;
  const gdouble dtau_dx = nc_thermodyn_recomb_dtau_dx (recomb, x);
  const gdouble d2tau_dx2 = -gsl_spline_eval_deriv (recomb->dtau_dx_spline, -x, recomb->dtau_dx_accel);
  const gdouble d3tau_dx3 = -gsl_spline_eval_deriv2 (recomb->dtau_dx_spline, -x, recomb->dtau_dx_accel);
  const gdouble x2 = x * x;

  *y = (x2 * gsl_pow_2 (dtau_dx) - (x * dtau_dx + x2 * d2tau_dx2));
  *dy = 2.0 * x * dtau_dx * (dtau_dx + x * d2tau_dx2) - dtau_dx - 3.0 * x * d2tau_dx2 - x2 * d3tau_dx3;
}

/**
 * nc_thermodyn_recomb_peak_x:
 * @recomb: a #NcThermodynRecomb
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_thermodyn_recomb_peak_x (NcThermodynRecomb *recomb)
{
  nc_thermodyn_recomb_dtau_dx_init_spline (recomb);
  {
    gint iter = 0, status, max_iter = 100;
    const gsl_root_fdfsolver_type *T;
    gsl_root_fdfsolver *s;
    double x0, x = 1000.0;
    gsl_function_fdf FDF;

    FDF.f = &gbar_exp_tau_f;
    FDF.df = &gbar_exp_tau_df;
    FDF.fdf = &gbar_exp_tau_fdf;
    FDF.params = recomb;

//    T = gsl_root_fdfsolver_newton;
    T = gsl_root_fdfsolver_steffenson;
    s = gsl_root_fdfsolver_alloc (T);
    gsl_root_fdfsolver_set (s, &FDF, x);

    do
    {
      iter++;
      gsl_root_fdfsolver_iterate (s);
      x0 = x;
      x = gsl_root_fdfsolver_root (s);
      status = gsl_root_test_delta (x, x0, 0, 1e-13);
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fdfsolver_free (s);
    return x;
  }
}
