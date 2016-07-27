/***************************************************************************
 *            nc_hipert_two_fluids.c
 *
 *  Tue June 10 19:14:57 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_two_fluids.c
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_hipert_two_fluids
 * @title: NcHIPertTwoFluids
 * @short_description: Perturbation object for a two fluids system. 
 *
 * This object provides the computation of the two fluid system of cosmological 
 * perturbations. This problem is decribed by two fluids with energy density and pressure
 * given respectively by $\bar{\rho}_i$ and $\bar{p}_i$ for $i = 1,2$.
 *
 * The system is written in terms of the gauge invariant variable
 * $$\zeta \equiv \Psi - \frac{2\bar{K}}{\kappa(\bar{\rho} + \bar{p})} + E\mathcal{V},$$
 * and the entropy mode $$S = \frac{\kappa\varpi}{x^3 H}(\mathcal{U}_1 - \mathcal{U}_2),$$ 
 * where $\mathcal{U}_i \equiv \psi + E\mathcal{V}_i$ and
 * $$\varpi \equiv \frac{(\bar{\rho}_1+\bar{p}_1)(\bar{\rho}_2+\bar{p}_2)}{\bar{\rho}+\bar{p}}.$$
 * 
 * Their momentum are 
 * \begin{split}
 * P_\zeta &= \frac{2\bar{D}^2_\bar{K}\Psi}{x^3E}, \\\\
 * P_S &= \frac{\delta\rho_2}{\bar{\rho}_2+\bar{p}_2} - \frac{\delta\rho_1}{\bar{\rho}_1 + \bar{p}_1}.
 * \end{split}
 * 
 * The equations of motion in their first order form are
 * \begin{align}
 * \zeta^\prime &= \frac{P_\zeta}{m_\zeta} + Y S, \\\\
 * P_\zeta^\prime &= -m_\zeta\mu_\zeta^2\zeta, \\\\
 * S^\prime &= \frac{P_S}{m_S} + Y \zeta, \\\\
 * P_S^\prime &= -m_S\mu_S^2S.
 * \end{align}
 * The mass $m_\zeta$ and the frequency $\mu_\zeta$ are defined by
 * \begin{align}
 * m_\zeta     &= \frac{3\Delta_\bar{K}(\bar{\rho} + \bar{p})}{\rho_\text{crit0} N x^3 c_s^2 E^2}, \\\\
 * \mu_\zeta^2 &= x^2N^2c_s^2k^2, \\\\
 * m_S         &= \frac{x^3}{c_m^2\varpi N}, \\\\
 * \mu_S^2     &= x^2N^2c_m^2k^2, \\\
 * Y           &= \frac{c_n^2}{c_s^2c_m^2}\frac{1}{m_\zeta m_S \Delta_\bar{K} N E}.
 * \end{align}
 * where $\bar{\rho} + \bar{p}$ is the background total energy density plus pressure,
 * $E^2 = H^2/H_0^2$ is the dimensionless Hubble function squared (nc_hicosmo_E2()), $c_s^2$ the speed of sound,
 * $N$ is the lapse function that in this case (using $\alpha$ as time variable) is $N \equiv \vert{}E\vert^{-1}$, $\rho_\text{crit0}$
 * is the critical density today defined by $\rho_\text{crit0} \equiv 3H_0^2/\kappa$ and $$\Delta_\bar{K} \equiv \frac{k^2}{k^2 + \Omega_{k0}}.$$
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline_cubic_notaknot.h"
#include "perturbations/nc_hipert_two_fluids.h"
#include "perturbations/nc_hipert_itwo_fluids.h"

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_band.h>
#include <nvector/nvector_serial.h> 
#include <gsl/gsl_roots.h>

G_DEFINE_TYPE (NcHIPertTwoFluids, nc_hipert_two_fluids, NC_TYPE_HIPERT);

typedef struct _NcHIPertTwoFluidsArg
{
  NcHICosmo *cosmo;
  NcHIPertTwoFluids *ptf;
} NcHIPertTwoFluidsArg;

static void
nc_hipert_two_fluids_init (NcHIPertTwoFluids *ptf)
{
  ptf->wkb_zeta = NULL;
  ptf->wkb_S    = NULL;
}

static void
nc_hipert_two_fluids_dispose (GObject *object)
{
  NcHIPertTwoFluids *ptf = NC_HIPERT_TWO_FLUIDS (object);

  nc_hipert_wkb_clear (&ptf->wkb_zeta);
  nc_hipert_wkb_clear (&ptf->wkb_S);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_two_fluids_parent_class)->dispose (object);
}


static void
nc_hipert_two_fluids_finalize (GObject *object)
{
  /*NcHIPertTwoFluids *ptf = NC_HIPERT_TWO_FLUIDS (object);*/
    
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_two_fluids_parent_class)->finalize (object);
}

static void _nc_hipert_two_fluids_set_mode_k (NcHIPert *pert, gdouble k);
static void _nc_hipert_two_fluids_set_abstol (NcHIPert *pert, gdouble abstol);
static void _nc_hipert_two_fluids_set_reltol (NcHIPert *pert, gdouble reltol); 

static void
nc_hipert_two_fluids_class_init (NcHIPertTwoFluidsClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->dispose  = nc_hipert_two_fluids_dispose;
  object_class->finalize = nc_hipert_two_fluids_finalize;

  NC_HIPERT_CLASS (klass)->set_mode_k = &_nc_hipert_two_fluids_set_mode_k;
  NC_HIPERT_CLASS (klass)->set_abstol = &_nc_hipert_two_fluids_set_abstol;
  NC_HIPERT_CLASS (klass)->set_reltol = &_nc_hipert_two_fluids_set_reltol;  
}

static void 
_nc_hipert_two_fluids_set_mode_k (NcHIPert *pert, gdouble k) 
{
  NC_HIPERT_CLASS (nc_hipert_two_fluids_parent_class)->set_mode_k (pert, k);
  /* Chain up : start */
  if (!pert->prepared)
  {
/*
    NcHIPertTwoFluids *ptf = NC_HIPERT_TWO_FLUIDS (pert);
    nc_hipert_set_mode_k (NC_HIPERT (ptf->wkb_zeta), k);
    nc_hipert_set_mode_k (NC_HIPERT (ptf->wkb_S), k);
    */
  }
}

static void 
_nc_hipert_two_fluids_set_abstol (NcHIPert *pert, gdouble abstol) 
{
  NC_HIPERT_CLASS (nc_hipert_two_fluids_parent_class)->set_abstol (pert, abstol);
  /* Chain up : start */
  if (!pert->prepared)
  {
/*
    NcHIPertTwoFluids *ptf = NC_HIPERT_TWO_FLUIDS (pert);
    nc_hipert_set_abstol (NC_HIPERT (ptf->wkb_zeta), abstol);
    nc_hipert_set_abstol (NC_HIPERT (ptf->wkb_S), abstol);
    */
  }
}

static void 
_nc_hipert_two_fluids_set_reltol (NcHIPert *pert, gdouble reltol) 
{
  NC_HIPERT_CLASS (nc_hipert_two_fluids_parent_class)->set_reltol (pert, reltol);
  /* Chain up : start */
  if (!pert->prepared)
  {
/*
    NcHIPertTwoFluids *ptf = NC_HIPERT_TWO_FLUIDS (pert);
    nc_hipert_set_reltol (NC_HIPERT (ptf->wkb_zeta), reltol);
    nc_hipert_set_reltol (NC_HIPERT (ptf->wkb_S), reltol);
    */
  }
}

/**
 * nc_hipert_two_fluids_new:
 * 
 * Creates a new #NcHIPertTwoFluids object.
 * 
 * Returns: (transfer full): a new #NcHIPertTwoFluids.
 */
NcHIPertTwoFluids *
nc_hipert_two_fluids_new (void)
{
  NcHIPertTwoFluids *ptf = g_object_new (NC_TYPE_HIPERT_TWO_FLUIDS,
                                         "sys-size", 2 * NC_HIPERT_ITWO_FLUIDS_VARS_LEN,
                                         NULL);

  return ptf;
}

/**
 * nc_hipert_two_fluids_ref:
 * @ptf: a #NcHIPertTwoFluids.
 * 
 * Increases the reference count of @ptf.
 * 
 * Returns: (transfer full): @ptf. 
 */
NcHIPertTwoFluids *
nc_hipert_two_fluids_ref (NcHIPertTwoFluids *ptf)
{
  return g_object_ref (ptf);
}

/**
 * nc_hipert_two_fluids_free:
 * @ptf: a #NcHIPertTwoFluids.
 * 
 * Decreases the reference count of @ptf.
 * 
 */
void 
nc_hipert_two_fluids_free (NcHIPertTwoFluids *ptf)
{
  g_object_unref (ptf);
}

/**
 * nc_hipert_two_fluids_clear:
 * @ptf: a #NcHIPertTwoFluids.
 * 
 * Decreases the reference count of *@ptf and sets *@ptf to NULL.
 * 
 */
void 
nc_hipert_two_fluids_clear (NcHIPertTwoFluids **ptf)
{
  g_clear_object (ptf);
}

/**
 * nc_hipert_two_fluids_prepare_wkb_zeta:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @prec: Required precision.
 * @alpha_i: initial log-redshift time.
 * @alpha_f: final log-redshift time.
 * 
 * Prepare the zeta component of the object for WKB calculations using the cosmology @cosmo.
 * 
 */
void 
nc_hipert_two_fluids_prepare_wkb_zeta (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble prec, gdouble alpha_i, gdouble alpha_f)
{
  nc_hipert_wkb_prepare (ptf->wkb_zeta, G_OBJECT (cosmo), prec, alpha_i, alpha_f);
}

/**
 * nc_hipert_two_fluids_prepare_wkb_S:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @prec: Required precision.
 * @alpha_i: initial log-redshift time.
 * @alpha_f: final log-redshift time.
 * 
 * Prepare the zeta component of the object for WKB calculations using the cosmology @cosmo.
 * 
 */
void 
nc_hipert_two_fluids_prepare_wkb_S (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble prec, gdouble alpha_i, gdouble alpha_f)
{
  nc_hipert_wkb_prepare (ptf->wkb_S, G_OBJECT (cosmo), prec, alpha_i, alpha_f);
}

/**
 * nc_hipert_two_fluids_nuA:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha: log-redshift time.
 * 
 * FIXME
 * 
 * Return: FIXME
 */
gdouble 
nc_hipert_two_fluids_nuA (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha)
{
  return nc_hipert_wkb_nuA (ptf->wkb_zeta, G_OBJECT (cosmo), alpha);
}

/**
 * nc_hipert_two_fluids_nuB:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha: log-redshift time.
 * 
 * FIXME
 * 
 * Return: FIXME
 */
gdouble 
nc_hipert_two_fluids_nuB (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha)
{
  return nc_hipert_wkb_nuA (ptf->wkb_S, G_OBJECT (cosmo), alpha);
}

/**
 * nc_hipert_two_fluids_eom:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha: log-redshift time.
 * @eom: (out callee-allocates) (transfer none): Equation of motion variables.
 * 
 * FIXME
 * 
 */
void 
nc_hipert_two_fluids_eom (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, NcHIPertITwoFluidsEOM **eom)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  return;
}

/**
 * nc_hipert_two_fluids_get_init_cond:
 * @ptf: a #NcHIPertTwoFluids
 * @cosmo: a #NcHICosmo
 * @alpha: the log-redshift time
 * @main_mode: main mode
 * @beta_R: mode $R$ initial phase
 * @state: (out caller-allocates) (array fixed-size=8):
 * 
 * FIXME 
 * 
 */
void 
nc_hipert_two_fluids_get_init_cond (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, guint main_mode, const gdouble beta_R, gdouble *state)
{
  NcHIPert *pert             = NC_HIPERT (ptf);
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  const gdouble beta_I       = beta_R + 0.5 * M_PI; 

  switch (main_mode)
  {
    case 1:
    {
      const gdouble dsigma1     = eom->nu1 - 0.5 * eom->gammabar11;
      const gdouble dsigma2     = eom->nu2 - 0.5 * eom->gammabar22;
      const gdouble Deltadsigma = dsigma1 - dsigma2;

      const complex double a_R_1  = cexp (+ I * beta_R) / M_SQRT2;
      const complex double ac_R_1 = cexp (- I * beta_R) / M_SQRT2;
      const complex double a_I_1  = cexp (+ I * beta_I) / M_SQRT2;
      const complex double ac_I_1 = cexp (- I * beta_I) / M_SQRT2;

      const complex double A_R11  = a_R_1 - 0.25 * eom->gammabar11 * ac_R_1 / dsigma1;
      const complex double A_R12  = 0.5 * (a_R_1 * (- eom->gammabar12 + I * eom->taubar) - ac_R_1 * eom->gammabar12) / Deltadsigma;

      const complex double A_I11  = a_I_1 - 0.25 * eom->gammabar11 * ac_I_1 / dsigma1;
      const complex double A_I12  = 0.5 * (a_I_1 * (- eom->gammabar12 + I * eom->taubar) - ac_I_1 * eom->gammabar12) / Deltadsigma;

      state[NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1] = NC_HIPERT_TWO_FLUIDS_A2Q (A_R11);
      state[NC_HIPERT_ITWO_FLUIDS_VARS_P_R1] = NC_HIPERT_TWO_FLUIDS_A2P (A_R11);

      state[NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2] = NC_HIPERT_TWO_FLUIDS_A2Q (A_R12);
      state[NC_HIPERT_ITWO_FLUIDS_VARS_P_R2] = NC_HIPERT_TWO_FLUIDS_A2P (A_R12);
    
      state[NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1] = NC_HIPERT_TWO_FLUIDS_A2Q (A_I11);
      state[NC_HIPERT_ITWO_FLUIDS_VARS_P_I1] = NC_HIPERT_TWO_FLUIDS_A2P (A_I11);

      state[NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2] = NC_HIPERT_TWO_FLUIDS_A2Q (A_I12);
      state[NC_HIPERT_ITWO_FLUIDS_VARS_P_I2] = NC_HIPERT_TWO_FLUIDS_A2P (A_I12);
    
      break;
    }
    case 2:
    {
      const gdouble dsigma1     = eom->nu1 - 0.5 * eom->gammabar11;
      const gdouble dsigma2     = eom->nu2 - 0.5 * eom->gammabar22;
      const gdouble Deltadsigma = dsigma1 - dsigma2;

      const complex double a_R_2  = cexp (+ I * beta_R) / M_SQRT2;
      const complex double ac_R_2 = cexp (- I * beta_R) / M_SQRT2;
      const complex double a_I_2  = cexp (+ I * beta_I) / M_SQRT2;
      const complex double ac_I_2 = cexp (- I * beta_I) / M_SQRT2;
    
      const complex double A_R22  = a_R_2 - 0.25 * eom->gammabar22 * ac_R_2 / dsigma2;
      const complex double A_R21  = 0.5 * (a_R_2 * (eom->gammabar12 + I * eom->taubar) - ac_R_2 * eom->gammabar12) / Deltadsigma;

      const complex double A_I22  = a_I_2 - 0.25 * eom->gammabar22 * ac_I_2 / dsigma2;
      const complex double A_I21  = 0.5 * (a_I_2 * (eom->gammabar12 + I * eom->taubar) - ac_I_2 * eom->gammabar12) / Deltadsigma;
      
      state[NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1] = NC_HIPERT_TWO_FLUIDS_A2Q (A_R21);
      state[NC_HIPERT_ITWO_FLUIDS_VARS_P_R1] = NC_HIPERT_TWO_FLUIDS_A2P (A_R21);

      state[NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2] = NC_HIPERT_TWO_FLUIDS_A2Q (A_R22);
      state[NC_HIPERT_ITWO_FLUIDS_VARS_P_R2] = NC_HIPERT_TWO_FLUIDS_A2P (A_R22);

      state[NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1] = NC_HIPERT_TWO_FLUIDS_A2Q (A_I21);
      state[NC_HIPERT_ITWO_FLUIDS_VARS_P_I1] = NC_HIPERT_TWO_FLUIDS_A2P (A_I21);

      state[NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2] = NC_HIPERT_TWO_FLUIDS_A2Q (A_I22);
      state[NC_HIPERT_ITWO_FLUIDS_VARS_P_I2] = NC_HIPERT_TWO_FLUIDS_A2P (A_I22);

      break;
    }
    default:
       g_error ("nc_hipert_two_fluids_get_init_cond: Unknown mode %u.", main_mode);
      break;
  } 
}

/**
 * nc_hipert_two_fluids_to_zeta_s:
 * @ptf: a #NcHIPertTwoFluids
 * @cosmo: a #NcHICosmo
 * @alpha: the log-redshift time
 * @state: (inout) (array fixed-size=8):
 * 
 * FIXME 
 * 
 */
void 
nc_hipert_two_fluids_to_zeta_s (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, gdouble *state)
{
  NcHIPert *pert           = NC_HIPERT (ptf);
  NcHIPertITwoFluidsTV *tv = nc_hipert_itwo_fluids_tv_eval (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  gdouble zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_LEN] = {
    0.0, 0.0, 0.0, 0.0, 
    0.0, 0.0, 0.0, 0.0
  };
  const guint syssize = NC_HIPERT_ITWO_FLUIDS_VARS_LEN / 2;
  gint i;

  for (i = 0; i < syssize; i++)
  {
    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R]  += tv->zeta[i]  * state[i];
    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_S_R]     += tv->s[i]     * state[i];
    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R] += tv->Pzeta[i] * state[i];
    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_PS_R]    += tv->Ps[i]    * state[i];

    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I]  += tv->zeta[i]  * state[syssize + i];
    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_S_I]     += tv->s[i]     * state[syssize + i];
    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I] += tv->Pzeta[i] * state[syssize + i];
    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_PS_I]    += tv->Ps[i]    * state[syssize + i];
    
  }

  for (i = 0; i < NC_HIPERT_ITWO_FLUIDS_VARS_LEN; i++)
    state[i] = zeta_s[i];
}

static gint
_nc_hipert_two_fluids_f1 (realtype alpha, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHIPertTwoFluidsArg *arg = (NcHIPertTwoFluidsArg *) f_data;
  const gdouble k = NC_HIPERT (arg->ptf)->k;
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);

  const gdouble Q_R1        = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1);
  const gdouble P_R1        = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1);
  const gdouble Q_R2        = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2);
  const gdouble P_R2        = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2);

  const gdouble Q_I1        = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1);
  const gdouble P_I1        = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1);
  const gdouble Q_I2        = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2);
  const gdouble P_I2        = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2);

  const complex double A_R1 = NC_HIPERT_TWO_FLUIDS_QP2A (Q_R1, P_R1);
  const complex double A_R2 = NC_HIPERT_TWO_FLUIDS_QP2A (Q_R2, P_R2);
  
  const complex double A_I1 = NC_HIPERT_TWO_FLUIDS_QP2A (Q_I1, P_I1);
  const complex double A_I2 = NC_HIPERT_TWO_FLUIDS_QP2A (Q_I2, P_I2);

  const gdouble gammabar11  = eom->gammabar11;
  const gdouble gammabar22  = eom->gammabar22;
  const gdouble gammabar12  = eom->gammabar12;
  const gdouble taubar12    = eom->taubar;

  const complex double dA_R1 = I * eom->nu1 * A_R1 + gammabar11 * Q_R1 + gammabar12 * Q_R2 + 0.5 * taubar12 * A_R2;
  const complex double dA_R2 = I * eom->nu2 * A_R2 + gammabar22 * Q_R2 + gammabar12 * Q_R1 - 0.5 * taubar12 * A_R1;

  const complex double dA_I1 = I * eom->nu1 * A_I1 + gammabar11 * Q_I1 + gammabar12 * Q_I2 + 0.5 * taubar12 * A_I2;
  const complex double dA_I2 = I * eom->nu2 * A_I2 + gammabar22 * Q_I2 + gammabar12 * Q_I1 - 0.5 * taubar12 * A_I1;

  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1) = NC_HIPERT_TWO_FLUIDS_A2Q (dA_R1);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1) = NC_HIPERT_TWO_FLUIDS_A2P (dA_R1);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2) = NC_HIPERT_TWO_FLUIDS_A2Q (dA_R2);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2) = NC_HIPERT_TWO_FLUIDS_A2P (dA_R2);

  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1) = NC_HIPERT_TWO_FLUIDS_A2Q (dA_I1);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1) = NC_HIPERT_TWO_FLUIDS_A2P (dA_I1);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2) = NC_HIPERT_TWO_FLUIDS_A2Q (dA_I2);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2) = NC_HIPERT_TWO_FLUIDS_A2P (dA_I2);
  
  NV_Ith_S (ydot, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = NV_Ith_S (y, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) / eom->m_zeta + eom->y * NV_Ith_S (y, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_R);
  NV_Ith_S (ydot, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = NV_Ith_S (y, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    / eom->m_s    + eom->y * NV_Ith_S (y, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R);
  NV_Ith_S (ydot, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = - eom->mnu2_zeta * NV_Ith_S (y, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R);
  NV_Ith_S (ydot, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = - eom->mnu2_s    * NV_Ith_S (y, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_S_R);

  NV_Ith_S (ydot, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = NV_Ith_S (y, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) / eom->m_zeta + eom->y * NV_Ith_S (y, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_I);
  NV_Ith_S (ydot, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = NV_Ith_S (y, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    / eom->m_s    + eom->y * NV_Ith_S (y, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I);
  NV_Ith_S (ydot, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = - eom->mnu2_zeta * NV_Ith_S (y, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I);
  NV_Ith_S (ydot, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = - eom->mnu2_s    * NV_Ith_S (y, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_S_I);
  
#ifdef NHACA
  if (last_alpha + 1.0e-2 < alpha)
  {
    /*printf ("% 21.15g |  F % 21.15g % 21.15g % 21.15g % 21.15g\n", alpha, creal (A1), cimag (A1), creal (A2), cimag (A2));*/
    /*printf ("% 21.15g | DF % 21.15g % 21.15g % 21.15g % 21.15g\n", alpha, creal (dA1), cimag (dA1), creal (dA2), cimag (dA2));*/
  
    printf ("% 21.15f | [ % 21.15e % 21.15e % 21.15e % 21.15e = % 21.15e <=> % 21.15e ]||[ % 21.15e % 21.15e = % 21.15e <=> % 21.15e ]||[ % 21.15e % 21.15e % 21.15e % 21.15e = % 21.15e <=> % 21.15e ]||[ % 21.15e % 21.15e = % 21.15e <=> % 21.15e ]\n", 
            alpha, 
            - eom->nu1 * Q1, gammabar11 * Q1, gammabar12 * Q2,  0.5 * taubar12 * P2, - eom->nu1 * Q1 + gammabar11 * Q1 + gammabar12 * Q2 + 0.5 * taubar12 * P2, NV_Ith_S (y, 0), 
              eom->nu1 * P1,   0.5 * taubar12 * Q2, eom->nu1 * P1 +  0.5 * taubar12 * Q2, NV_Ith_S (y, 1),
            - eom->nu2 * Q2, gammabar22 * Q2, gammabar12 * Q1, -0.5 * taubar12 * P1, - eom->nu2 * Q2 + gammabar22 * Q2 + gammabar12 * Q1 - 0.5 * taubar12 * P1, NV_Ith_S (y, 2),
              eom->nu2 * P2, - 0.5 * taubar12 * Q1, eom->nu2 * P2 - 0.5 * taubar12 * Q1, NV_Ith_S (y, 3)
            );
   
    //last_alpha = alpha;
  }
#endif
/*
  const gdouble sin_2theta  = sin (2.0 * theta1);
  const gdouble cos_2theta  = cos (2.0 * theta1);
  const complex double F2   = ReF2 + I * ImF2;
  
  NV_Ith_S (ydot, 0) = dsigma1 - 0.5 * (ReF2 * gammabar12 - ImF2 * taubar12) + 0.5 * (gammabar11 + gammabar12 * ReF2) * cos_2theta - 0.5 * ImF2 * sin_2theta * gammabar12;
  NV_Ith_S (ydot, 1) = 0.5 * (sin_2theta * (gammabar11 + gammabar12 * ReF2) - cos_2theta * gammabar12 * ImF2 + gammabar12 * ImF2 + taubar12 * ReF2);

  {
    const complex double expm2Itheta = cos_2theta - I * sin_2theta;
    const complex double F2star      = conj (F2);  
    const complex double dF2         = -I * Deltadsigma * F2 - 0.5 * (taubar12 + I * gammabar12) 
      + 0.5 * expm2Itheta * ( I * gammabar22 * F2star + I * gammabar12 - I * gammabar11 * F2)
      - 0.5 * ((taubar12 - I * gammabar12) * F2 * F2 + I * gammabar12 * F2 * F2star * expm2Itheta);

    NV_Ith_S (ydot, 2) = creal (dF2);
    NV_Ith_S (ydot, 3) = cimag (dF2);
  }
  */
  return 0;
}

static gint
_nc_hipert_two_fluids_J (_NCM_SUNDIALS_INT_TYPE N, realtype alpha, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHIPertTwoFluidsArg *arg = (NcHIPertTwoFluidsArg *) jac_data;
  const gdouble k = NC_HIPERT (arg->ptf)->k;
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);
  
  const gdouble gammabar11  = eom->gammabar11;
  const gdouble gammabar22  = eom->gammabar22;
  const gdouble gammabar12  = eom->gammabar12;
  const gdouble taubar12    = eom->taubar;

  /************************************************************************************************************/
  /* Solution R */
  /************************************************************************************************************/

  /* eom->nu1 * P1 +  0.5 * taubar12 * Q2 */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R1) = eom->nu1;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2) = 0.5 * taubar12;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R2) = 0.0;

  /* - eom->nu1 * Q1 + gammabar11 * Q1 + gammabar12 * Q2 + 0.5 * taubar12 * P2 */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1) = - eom->nu1 + gammabar11;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2) = gammabar12;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R2) = 0.5 * taubar12;

  /* eom->nu2 * P2 - 0.5 * taubar12 * Q1 */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1) = - 0.5 * taubar12;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R2) = eom->nu2;

  /* - eom->nu2 * Q2 + gammabar22 * Q2 + gammabar12 * Q1 - 0.5 * taubar12 * P1 */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1) = gammabar12;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R1) = - 0.5 * taubar12;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2) = - eom->nu2 + gammabar22;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R2) = 0.0;

  /************************************************************************************************************/
  /* Solution I */
  /************************************************************************************************************/

  /* eom->nu1 * P1 +  0.5 * taubar12 * Q2 */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I1) = eom->nu1;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2) = 0.5 * taubar12;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I2) = 0.0;

  /* - eom->nu1 * Q1 + gammabar11 * Q1 + gammabar12 * Q2 + 0.5 * taubar12 * P2 */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1) = - eom->nu1 + gammabar11;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2) = gammabar12;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I2) = 0.5 * taubar12;

  /* eom->nu2 * P2 - 0.5 * taubar12 * Q1 */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1) = - 0.5 * taubar12;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I2) = eom->nu2;

  /* - eom->nu2 * Q2 + gammabar22 * Q2 + gammabar12 * Q1 - 0.5 * taubar12 * P1 */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1) = gammabar12;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I1) = - 0.5 * taubar12;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2) = - eom->nu2 + gammabar22;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I2) = 0.0;

  /************************************************************************************************************/
  /* Solution R */
  /************************************************************************************************************/

  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R,  8 + NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = 0.0;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R,  8 + NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = 0.0;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R,  8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = 1.0 / eom->m_zeta;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R,  8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = eom->y;

  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_S_R,     8 + NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = 0.0;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_S_R,     8 + NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = 0.0;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_S_R,     8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = eom->y;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_S_R,     8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = 1.0 / eom->m_s;

  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = - eom->mnu2_zeta;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = 0.0;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = 0.0;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = 0.0;

  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_R,    8 + NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = 0.0;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_R,    8 + NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = - eom->mnu2_s;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_R,    8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = 0.0;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_R,    8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = 0.0;
  
  /************************************************************************************************************/
  /* Solution I */
  /************************************************************************************************************/

  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I,  8 + NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = 0.0;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I,  8 + NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = 0.0;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I,  8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = 1.0 / eom->m_zeta;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I,  8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = eom->y;

  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_S_I,     8 + NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = 0.0;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_S_I,     8 + NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = 0.0;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_S_I,     8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = eom->y;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_S_I,     8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = 1.0 / eom->m_s;

  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = - eom->mnu2_zeta;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = 0.0;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = 0.0;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = 0.0;

  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_I,    8 + NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = 0.0;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_I,    8 + NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = - eom->mnu2_s;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_I,    8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = 0.0;
  DENSE_ELEM (J, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_I,    8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = 0.0;
  
  return 0;
}

/**
 * nc_hipert_two_fluids_set_init_cond:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alphai: the log-redshift time.
 * @vars: (in) (array fixed-size=8) (element-type double): Perturbations variables conforming to #NcHIPertTwoFluidsVars.
 * 
 * Sets the initial conditions for the two fluids system evolution. 
 * 
 */
void 
nc_hipert_two_fluids_set_init_cond (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alphai, gdouble *vars)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  gint flag;
  guint i;

  pert->alpha0 = alphai;

  for (i = 0; i < 8; i++)
  {
    NV_Ith_S (pert->y, i) = vars[i];
  }

  nc_hipert_two_fluids_to_zeta_s (ptf, cosmo, alphai, vars);
  
  for (i = 0; i < 8; i++)
  {
    NV_Ith_S (pert->y, 8 + i) = vars[i];
  }
 
  if (!pert->cvode_init)
  {
    flag = CVodeInit (pert->cvode, &_nc_hipert_two_fluids_f1, alphai, pert->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
    pert->cvode_init = TRUE;
  }
  else
  {
    flag = CVodeReInit (pert->cvode, alphai, pert->y);
    NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );
  }

  flag = CVodeSetMaxStep (pert->cvode, 1.0);
  NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1,);
    
  flag = CVodeSStolerances (pert->cvode, pert->reltol, 0.0);
  NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1,);

  flag = CVodeSetMaxNumSteps (pert->cvode, 100000000);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

  flag = CVDense (pert->cvode, 2 * NC_HIPERT_ITWO_FLUIDS_VARS_LEN);
  NCM_CVODE_CHECK (&flag, "CVDense", 1, );

  flag = CVDlsSetDenseJacFn (pert->cvode, &_nc_hipert_two_fluids_J);
  NCM_CVODE_CHECK (&flag, "CVDlsSetDenseJacFn", 1, );  
}

/**
 * nc_hipert_two_fluids_evolve:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alphaf: the final log-redshift time.
 * 
 * Evolve the system until @alphaf.
 * 
 */
void 
nc_hipert_two_fluids_evolve (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alphaf)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  NcHIPertTwoFluidsArg arg;
  gdouble alpha_i = 0.0;
  gdouble alpha_t;
  const gdouble delta_alpha = 1.0e-3;
  gint flag;
  
  arg.cosmo = cosmo;
  arg.ptf   = ptf;

  flag = CVodeSetUserData (pert->cvode, &arg);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

  alpha_t = pert->alpha0 + delta_alpha;
  
  while (TRUE)
  {
    flag = CVode (pert->cvode, alpha_t, pert->y, &alpha_i, CV_NORMAL);
    NCM_CVODE_CHECK (&flag, "CVode[nc_hipert_two_fluids_evolve]", 1, );

    {
      const gdouble Q_R1        = NV_Ith_S (pert->y, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1);
      const gdouble Q_R2        = NV_Ith_S (pert->y, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2);
      const gdouble P_R1        = NV_Ith_S (pert->y, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1);
      const gdouble P_R2        = NV_Ith_S (pert->y, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2);
      const gdouble Q_I1        = NV_Ith_S (pert->y, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1);
      const gdouble Q_I2        = NV_Ith_S (pert->y, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2);
      const gdouble P_I1        = NV_Ith_S (pert->y, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1);
      const gdouble P_I2        = NV_Ith_S (pert->y, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2);
      /*const complex double A_R1 = NC_HIPERT_TWO_FLUIDS_QP2A (Q_R1, P_R1);*/
      /*const complex double A_I1 = NC_HIPERT_TWO_FLUIDS_QP2A (Q_I1, P_I1);*/
      /*const complex double A_R2 = NC_HIPERT_TWO_FLUIDS_QP2A (Q_R2, P_R2);*/
      /*const complex double A_I2 = NC_HIPERT_TWO_FLUIDS_QP2A (Q_I2, P_I2);*/
      
      NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha_i, pert->k);
      gdouble state[8];

      /*nc_hipert_two_fluids_get_init_cond (ptf, cosmo, alpha_i, 1, carg (A_R1), state);*/

      state[NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1] = Q_R1; 
      state[NC_HIPERT_ITWO_FLUIDS_VARS_P_R1] = P_R1; 
      state[NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2] = Q_R2;
      state[NC_HIPERT_ITWO_FLUIDS_VARS_P_R2] = P_R2; 

      state[NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1] = Q_I1; 
      state[NC_HIPERT_ITWO_FLUIDS_VARS_P_I1] = P_I1; 
      state[NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2] = Q_I2;
      state[NC_HIPERT_ITWO_FLUIDS_VARS_P_I2] = P_I2; 
      
      nc_hipert_two_fluids_to_zeta_s (ptf, cosmo, alpha_i, state);

#ifndef CUCUCU
      printf ("% 21.15f % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e\n", 
              alpha_i, nc_hicosmo_x_alpha (cosmo, alpha_i),
    /* 03 */  Q_R1, 
    /* 04 */  Q_R2, 
    /* 05 */  P_R1, 
    /* 06 */  P_R2, 
    /* 07 */  Q_I1, 
    /* 08 */  Q_I2, 
    /* 09 */  P_I1, 
    /* 10 */  P_I2, 
    /* 11 */  state[NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R],
    /* 12 */  state[NC_HIPERT_ITWO_FLUIDS_VARS_S_R],
    /* 13 */  state[NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R],
    /* 14 */  state[NC_HIPERT_ITWO_FLUIDS_VARS_PS_R],
    /* 15 */  state[NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I],
    /* 16 */  state[NC_HIPERT_ITWO_FLUIDS_VARS_S_I],
    /* 17 */  state[NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I],
    /* 18 */  state[NC_HIPERT_ITWO_FLUIDS_VARS_PS_I],
    /* 19 */  NV_Ith_S (pert->y, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R),
    /* 20 */  NV_Ith_S (pert->y, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_S_R),
    /* 21 */  NV_Ith_S (pert->y, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R),
    /* 22 */  NV_Ith_S (pert->y, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_R),
    /* 23 */  NV_Ith_S (pert->y, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I),
    /* 24 */  NV_Ith_S (pert->y, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_S_I),
    /* 25 */  NV_Ith_S (pert->y, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I),
    /* 26 */  NV_Ith_S (pert->y, 8 + NC_HIPERT_ITWO_FLUIDS_VARS_PS_I),
    /* 27 */  eom->nu1,
    /* 28 */  eom->nu2,
    /* 29 */  eom->gammabar11,
    /* 30 */  eom->gammabar22,
    /* 31 */  eom->gammabar12,
    /* 32 */  eom->taubar
              );
#endif
/*
      printf ("% 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e\n", 
              alpha_i, nc_hicosmo_x_alpha (cosmo, alpha_i),
              Q1, Q2, P1, P2,
              fabs (state[NC_HIPERT_ITWO_FLUIDS_VARS_Q1] / Q1 - 1.0), 
              fabs (state[NC_HIPERT_ITWO_FLUIDS_VARS_Q2] / Q2 - 1.0), 
              fabs (state[NC_HIPERT_ITWO_FLUIDS_VARS_P1] / P1 - 1.0), 
              fabs (state[NC_HIPERT_ITWO_FLUIDS_VARS_P2] / P2 - 1.0),
              eom->nu1, eom->gammabar11, eom->nu2, eom->gammabar22, eom->gammabar12, eom->taubar
              );
*/      
    }

    if (alpha_i >= alphaf)
      break;
    else
    {
      alpha_t += delta_alpha;
    }
  }
  
  pert->alpha0 = alpha_i;
}

/**
 * nc_hipert_two_fluids_get_values:
 * @ptf: a #NcHIPertTwoFluids.
 * @alphai: (out caller-allocates): Current time.
 * @vars: (inout) (array fixed-size=8) (element-type double): Perturbations variables conforming to #NcHIPertTwoFluidsVars.
 * 
 * Get the current time and values of the numerical solution.
 * 
 */
void
nc_hipert_two_fluids_get_values (NcHIPertTwoFluids *ptf, gdouble *alphai, gdouble **vars)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  guint i;

  *alphai = pert->alpha0;
  for (i = 0; i < NC_HIPERT_TWO_FLUIDS_LEN; i++)
  {
    vars[0][i] = NV_Ith_S (pert->y, i);
  }
}
