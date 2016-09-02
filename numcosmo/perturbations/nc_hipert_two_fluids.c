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
#include <gsl/gsl_odeiv2.h>

#include <arkode/arkode.h>
#include <arkode/arkode_dense.h>


G_DEFINE_TYPE (NcHIPertTwoFluids, nc_hipert_two_fluids, NC_TYPE_HIPERT);

typedef struct _NcHIPertTwoFluidsArg
{
  NcHICosmo *cosmo;
  NcHIPertTwoFluids *ptf;
  gdouble prec;
} NcHIPertTwoFluidsArg;

static void
nc_hipert_two_fluids_init (NcHIPertTwoFluids *ptf)
{
  ptf->wkb_zeta = NULL;
  ptf->wkb_S    = NULL;
  ptf->abstol   = NULL;
  ptf->useQP    = FALSE;
  ptf->state    = ncm_vector_new (NC_HIPERT_ITWO_FLUIDS_VARS_LEN);

  ptf->arg = g_new0 (NcHIPertTwoFluidsArg, 1);

  ptf->arkode = ARKodeCreate ();
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
  NcHIPertTwoFluids *ptf = NC_HIPERT_TWO_FLUIDS (object);

  if (ptf->arkode != NULL)
  {
    ARKodeFree (&ptf->arkode);
    ptf->arkode = NULL;
  }
  
  g_free (ptf->arg);
    
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
                                         "sys-size", NC_HIPERT_ITWO_FLUIDS_VARS_LEN,
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
  /*nc_hipert_wkb_prepare (ptf->wkb_zeta, G_OBJECT (cosmo), prec, alpha_i, alpha_f);*/
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
  /*nc_hipert_wkb_prepare (ptf->wkb_S, G_OBJECT (cosmo), prec, alpha_i, alpha_f);*/
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
 * nc_hipert_two_fluids_get_init_cond_QP:
 * @ptf: a #NcHIPertTwoFluids
 * @cosmo: a #NcHICosmo
 * @alpha: the log-redshift time
 * @main_mode: main mode
 * @beta_R: mode $R$ initial phase
 * @init_cond: a #NcmVector (size >= 8) where to put the initial conditions
 * 
 * Calculates the initial condition for the $(Q,\,P)$ system with initial phase
 * for the R solution $\beta_R = $  @beta_R. The variable @main_mode chooses
 * which mode is excited (1 or 2).
 * 
 */
void 
nc_hipert_two_fluids_get_init_cond_QP (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, guint main_mode, const gdouble beta_R, NcmVector *init_cond)
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

      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1, NC_HIPERT_TWO_FLUIDS_A2Q (A_R11));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1, NC_HIPERT_TWO_FLUIDS_A2P (A_R11));

      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2, NC_HIPERT_TWO_FLUIDS_A2Q (A_R12));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2, NC_HIPERT_TWO_FLUIDS_A2P (A_R12));
    
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1, NC_HIPERT_TWO_FLUIDS_A2Q (A_I11));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1, NC_HIPERT_TWO_FLUIDS_A2P (A_I11));

      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2, NC_HIPERT_TWO_FLUIDS_A2Q (A_I12));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2, NC_HIPERT_TWO_FLUIDS_A2P (A_I12));
    
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
      const complex double A_R21  = 0.5 * (a_R_2 * (I * eom->taubar + eom->gammabar12) - ac_R_2 * eom->gammabar12) / Deltadsigma;

      const complex double A_I22  = a_I_2 - 0.25 * eom->gammabar22 * ac_I_2 / dsigma2;
      const complex double A_I21  = 0.5 * (a_I_2 * (I * eom->taubar + eom->gammabar12) - ac_I_2 * eom->gammabar12) / Deltadsigma;
      
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1, NC_HIPERT_TWO_FLUIDS_A2Q (A_R21));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1, NC_HIPERT_TWO_FLUIDS_A2P (A_R21));

      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2, NC_HIPERT_TWO_FLUIDS_A2Q (A_R22));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2, NC_HIPERT_TWO_FLUIDS_A2P (A_R22));

      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1, NC_HIPERT_TWO_FLUIDS_A2Q (A_I21));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1, NC_HIPERT_TWO_FLUIDS_A2P (A_I21));

      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2, NC_HIPERT_TWO_FLUIDS_A2Q (A_I22));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2, NC_HIPERT_TWO_FLUIDS_A2P (A_I22));

      break;
    }
    default:
       g_error ("nc_hipert_two_fluids_get_init_cond: Unknown main mode %u.", main_mode);
      break;
  } 
}

/**
 * nc_hipert_two_fluids_get_init_cond_zetaS:
 * @ptf: a #NcHIPertTwoFluids
 * @cosmo: a #NcHICosmo
 * @alpha: the log-redshift time
 * @main_mode: main mode
 * @beta_R: mode $R$ initial phase
 * @init_cond: a #NcmVector (size >= 8) where to put the initial conditions
 * 
 * Calculates the initial condition for the $\zeta{}S$ system with initial phase
 * for the R solution $\beta_R = $ @beta_R. The variable @main_mode chooses
 * which mode is excited (1 or 2).
 * 
 */
void 
nc_hipert_two_fluids_get_init_cond_zetaS (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, guint main_mode, const gdouble beta_R, NcmVector *init_cond)
{
  nc_hipert_two_fluids_get_init_cond_QP (ptf, cosmo, alpha, main_mode, beta_R, init_cond);
  nc_hipert_two_fluids_to_zeta_s (ptf, cosmo, alpha, init_cond);
}

/**
 * nc_hipert_two_fluids_to_zeta_s:
 * @ptf: a #NcHIPertTwoFluids
 * @cosmo: a #NcHICosmo
 * @alpha: the log-redshift time
 * @state: a #NcmVector (size >= 8) current state in $(Q,\,P)$ variables
 * 
 * Transform in-place the variables @init_cond from $(Q,\,P)$ to $(\zeta,\,S)$, assuming
 * they are calculated at $\alpha$ = @alpha.
 * 
 */
void 
nc_hipert_two_fluids_to_zeta_s (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, NcmVector *state)
{
  NcHIPert *pert           = NC_HIPERT (ptf);
  NcHIPertITwoFluidsTV *tv = nc_hipert_itwo_fluids_tv_eval (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  const guint syssize      = NC_HIPERT_ITWO_FLUIDS_VARS_LEN / 2;
  gdouble zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_LEN] = {
    0.0, 0.0, 0.0, 0.0, 
    0.0, 0.0, 0.0, 0.0
  };
  gint i;

  for (i = 0; i < syssize; i++)
  {
    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R]  += tv->zeta[i]  * ncm_vector_get (state, i);
    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_S_R]     += tv->s[i]     * ncm_vector_get (state, i);
    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R] += tv->Pzeta[i] * ncm_vector_get (state, i);
    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_PS_R]    += tv->Ps[i]    * ncm_vector_get (state, i);

    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I]  += tv->zeta[i]  * ncm_vector_get (state, syssize + i);
    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_S_I]     += tv->s[i]     * ncm_vector_get (state, syssize + i);
    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I] += tv->Pzeta[i] * ncm_vector_get (state, syssize + i);
    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_PS_I]    += tv->Ps[i]    * ncm_vector_get (state, syssize + i);
  }

  for (i = 0; i < NC_HIPERT_ITWO_FLUIDS_VARS_LEN; i++)
    ncm_vector_set (state, i, zeta_s[i]);
}

static gint
_nc_hipert_two_fluids_f_QP (realtype alpha, N_Vector y, N_Vector ydot, gpointer f_data)
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

  return 0;
}

static gint
_nc_hipert_two_fluids_f_QP_mode1 (realtype alpha, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHIPertTwoFluidsArg *arg = (NcHIPertTwoFluidsArg *) f_data;
  const gdouble k = NC_HIPERT (arg->ptf)->k;
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);

  const gdouble Q_R1        = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1);
  const gdouble P_R1        = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1);

  const gdouble Q_I1        = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1);
  const gdouble P_I1        = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1);

  const complex double A_R1 = NC_HIPERT_TWO_FLUIDS_QP2A (Q_R1, P_R1);  
  const complex double A_I1 = NC_HIPERT_TWO_FLUIDS_QP2A (Q_I1, P_I1);

  const gdouble gammabar11  = eom->gammabar11;
  const gdouble gammabar12  = eom->gammabar12;
  const gdouble taubar12    = eom->taubar;

  const complex double dA_R1 = I * eom->nu1 * A_R1 + gammabar11 * Q_R1;
  const complex double dA_R2 = gammabar12 * Q_R1 - 0.5 * taubar12 * A_R1;

  const complex double dA_I1 = I * eom->nu1 * A_I1 + gammabar11 * Q_I1;
  const complex double dA_I2 = gammabar12 * Q_I1 - 0.5 * taubar12 * A_I1;

  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1) = NC_HIPERT_TWO_FLUIDS_A2Q (dA_R1);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1) = NC_HIPERT_TWO_FLUIDS_A2P (dA_R1);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2) = NC_HIPERT_TWO_FLUIDS_A2Q (dA_R2);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2) = NC_HIPERT_TWO_FLUIDS_A2P (dA_R2);

  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1) = NC_HIPERT_TWO_FLUIDS_A2Q (dA_I1);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1) = NC_HIPERT_TWO_FLUIDS_A2P (dA_I1);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2) = NC_HIPERT_TWO_FLUIDS_A2Q (dA_I2);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2) = NC_HIPERT_TWO_FLUIDS_A2P (dA_I2);

  return 0;
}

static gint
_nc_hipert_two_fluids_f_QP_mode2 (realtype alpha, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHIPertTwoFluidsArg *arg = (NcHIPertTwoFluidsArg *) f_data;
  const gdouble k = NC_HIPERT (arg->ptf)->k;
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);

  const gdouble Q_R2        = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2);
  const gdouble P_R2        = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2);

  const gdouble Q_I2        = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2);
  const gdouble P_I2        = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2);

  const complex double A_R2 = NC_HIPERT_TWO_FLUIDS_QP2A (Q_R2, P_R2);
  
  const complex double A_I2 = NC_HIPERT_TWO_FLUIDS_QP2A (Q_I2, P_I2);

  const gdouble gammabar22  = eom->gammabar22;
  const gdouble gammabar12  = eom->gammabar12;
  const gdouble taubar12    = eom->taubar;

  const complex double dA_R1 = gammabar12 * Q_R2 + 0.5 * taubar12 * A_R2;
  const complex double dA_R2 = I * eom->nu2 * A_R2 + gammabar22 * Q_R2;

  const complex double dA_I1 = gammabar12 * Q_I2 + 0.5 * taubar12 * A_I2;
  const complex double dA_I2 = I * eom->nu2 * A_I2 + gammabar22 * Q_I2;

  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1) = NC_HIPERT_TWO_FLUIDS_A2Q (dA_R1);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1) = NC_HIPERT_TWO_FLUIDS_A2P (dA_R1);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2) = NC_HIPERT_TWO_FLUIDS_A2Q (dA_R2);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2) = NC_HIPERT_TWO_FLUIDS_A2P (dA_R2);

  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1) = NC_HIPERT_TWO_FLUIDS_A2Q (dA_I1);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1) = NC_HIPERT_TWO_FLUIDS_A2P (dA_I1);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2) = NC_HIPERT_TWO_FLUIDS_A2Q (dA_I2);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2) = NC_HIPERT_TWO_FLUIDS_A2P (dA_I2);

  return 0;
}

static gint
_nc_hipert_two_fluids_f_QP_mode1sub (realtype alpha, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHIPertTwoFluidsArg *arg = (NcHIPertTwoFluidsArg *) f_data;
  const gdouble k = NC_HIPERT (arg->ptf)->k;
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);

  const gdouble Q_R2        = NV_Ith_S (y, 0);
  const gdouble P_R2        = NV_Ith_S (y, 1);

  const gdouble Q_I2        = NV_Ith_S (y, 2);
  const gdouble P_I2        = NV_Ith_S (y, 3);

  const complex double A_R2 = NC_HIPERT_TWO_FLUIDS_QP2A (Q_R2, P_R2);
  const complex double A_I2 = NC_HIPERT_TWO_FLUIDS_QP2A (Q_I2, P_I2);
  
  const gdouble gammabar11  = eom->gammabar11;
  const gdouble gammabar22  = eom->gammabar22;
  const gdouble gammabar12  = eom->gammabar12;
  const gdouble taubar12    = eom->taubar;
  const gdouble dsigma1     = eom->nu1 - 0.5 * gammabar11;
  const gdouble dsigma2     = eom->nu2 - 0.5 * gammabar22;

  const complex double Lambda = gammabar12 + I * taubar12;
  
  const complex double A_R1   = (
                                 + A_R2        * ( Lambda     / (2.0 * (dsigma1 - dsigma2)) + gammabar11    * gammabar12 / (8.0 * dsigma1 * (dsigma1 + dsigma2)) ) 
                                 - conj (A_R2) * ( gammabar12 / (2.0 * (dsigma1 + dsigma2)) + conj (Lambda) * gammabar11 / (8.0 * dsigma1 * (dsigma1 - dsigma2)) )
                                 ) / (1.0 - gsl_pow_2 (gammabar11 / (4.0 * dsigma1)));
  const complex double A_I1   = (
                                 + A_I2        * ( Lambda     / (2.0 * (dsigma1 - dsigma2)) + gammabar11    * gammabar12 / (8.0 * dsigma1 * (dsigma1 + dsigma2)) ) 
                                 - conj (A_I2) * ( gammabar12 / (2.0 * (dsigma1 + dsigma2)) + conj (Lambda) * gammabar11 / (8.0 * dsigma1 * (dsigma1 - dsigma2)) )
                                 ) / (1.0 - gsl_pow_2 (gammabar11 / (4.0 * dsigma1)));

  const complex double dA_R2 = I * eom->nu2 * A_R2 + gammabar22 * Q_R2 + gammabar12 * cimag (A_R1) - 0.5 * taubar12 * A_R1;
  const complex double dA_I2 = I * eom->nu2 * A_I2 + gammabar22 * Q_I2 + gammabar12 * cimag (A_I1) - 0.5 * taubar12 * A_I1;

  NV_Ith_S (ydot, 0) = NC_HIPERT_TWO_FLUIDS_A2Q (dA_R2);
  NV_Ith_S (ydot, 1) = NC_HIPERT_TWO_FLUIDS_A2P (dA_R2);

  NV_Ith_S (ydot, 2) = NC_HIPERT_TWO_FLUIDS_A2Q (dA_I2);
  NV_Ith_S (ydot, 3) = NC_HIPERT_TWO_FLUIDS_A2P (dA_I2);

  return 0;
}

static gint
_nc_hipert_two_fluids_f_zetaS (realtype alpha, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHIPertTwoFluidsArg *arg = (NcHIPertTwoFluidsArg *) f_data;
  const gdouble k = NC_HIPERT (arg->ptf)->k;
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);

  const gdouble zeta_R      = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R);
  const gdouble S_R         = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_S_R);
  const gdouble Pzeta_R     = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R);
  const gdouble PS_R        = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R);

  const gdouble zeta_I      = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I);
  const gdouble S_I         = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_S_I);
  const gdouble Pzeta_I     = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I);
  const gdouble PS_I        = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I);
  
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = Pzeta_R / eom->m_zeta + eom->y * PS_R;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = PS_R    / eom->m_s    + eom->y * Pzeta_R;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = - eom->mnu2_zeta * zeta_R;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = - eom->mnu2_s    * S_R;
  
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = Pzeta_I / eom->m_zeta + eom->y * PS_I;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = PS_I    / eom->m_s    + eom->y * Pzeta_I;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = - eom->mnu2_zeta * zeta_I;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = - eom->mnu2_s    * S_I;
  
  return 0;
}

static gint
_nc_hipert_two_fluids_f_zetaS_zeta (realtype alpha, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHIPertTwoFluidsArg *arg = (NcHIPertTwoFluidsArg *) f_data;
  const gdouble k = NC_HIPERT (arg->ptf)->k;
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);

  const gdouble zeta_R      = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R);
  const gdouble Pzeta_R     = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R);

  const gdouble zeta_I      = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I);
  const gdouble Pzeta_I     = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I);
  
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = Pzeta_R / eom->m_zeta;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = eom->y * Pzeta_R;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = - eom->mnu2_zeta * zeta_R;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = 0.0;
  
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = Pzeta_I / eom->m_zeta;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = eom->y * Pzeta_I;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = - eom->mnu2_zeta * zeta_I;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = 0.0;
  
  return 0;
}

static gint
_nc_hipert_two_fluids_f_zetaS_S (realtype alpha, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHIPertTwoFluidsArg *arg = (NcHIPertTwoFluidsArg *) f_data;
  const gdouble k = NC_HIPERT (arg->ptf)->k;
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);

  const gdouble S_R         = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_S_R);
  const gdouble PS_R        = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R);

  const gdouble S_I         = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_S_I);
  const gdouble PS_I        = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I);
  
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = eom->y * PS_R;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = PS_R / eom->m_s;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = 0.0;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = - eom->mnu2_s    * S_R;
  
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = eom->y * PS_I;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = PS_I / eom->m_s;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = 0.0;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = - eom->mnu2_s    * S_I;

  return 0;
}

static gint
_nc_hipert_two_fluids_J_QP (_NCM_SUNDIALS_INT_TYPE N, realtype alpha, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
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
  
  return 0;
}

static gint
_nc_hipert_two_fluids_J_QP_mode1 (_NCM_SUNDIALS_INT_TYPE N, realtype alpha, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHIPertTwoFluidsArg *arg = (NcHIPertTwoFluidsArg *) jac_data;
  const gdouble k = NC_HIPERT (arg->ptf)->k;
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);
  
  const gdouble gammabar11  = eom->gammabar11;
  const gdouble gammabar12  = eom->gammabar12;
  const gdouble taubar12    = eom->taubar;

  /************************************************************************************************************/
  /* Solution R */
  /************************************************************************************************************/

  /* eom->nu1 * P1 */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R1) = eom->nu1;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R2) = 0.0;

  /* - eom->nu1 * Q1 + gammabar11 * Q1  */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1) = - eom->nu1 + gammabar11;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R2) = 0.0;

  /* - 0.5 * taubar12 * Q1 */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1) = - 0.5 * taubar12;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R2) = 0.0;

  /* gammabar12 * Q1 - 0.5 * taubar12 * P1 */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1) = gammabar12;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R1) = - 0.5 * taubar12;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R2) = 0.0;

  /************************************************************************************************************/
  /* Solution I */
  /************************************************************************************************************/

  /* eom->nu1 * P1 */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I1) = eom->nu1;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I2) = 0.0;

  /* - eom->nu1 * Q1 + gammabar11 * Q1 */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1) = - eom->nu1 + gammabar11;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I2) = 0.0;

  /* - 0.5 * taubar12 * Q1 */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1) = - 0.5 * taubar12;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I2) = 0.0;

  /* gammabar12 * Q1 - 0.5 * taubar12 * P1 */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1) = gammabar12;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I1) = - 0.5 * taubar12;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I2) = 0.0;
  
  return 0;
}

static gint
_nc_hipert_two_fluids_J_QP_mode2 (_NCM_SUNDIALS_INT_TYPE N, realtype alpha, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHIPertTwoFluidsArg *arg = (NcHIPertTwoFluidsArg *) jac_data;
  const gdouble k = NC_HIPERT (arg->ptf)->k;
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);

  const gdouble gammabar22  = eom->gammabar22;
  const gdouble gammabar12  = eom->gammabar12;
  const gdouble taubar12    = eom->taubar;

  /************************************************************************************************************/
  /* Solution R */
  /************************************************************************************************************/

  /* 0.5 * taubar12 * Q2 */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2) = 0.5 * taubar12;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R2) = 0.0;

  /* gammabar12 * Q2 + 0.5 * taubar12 * P2 */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2) = gammabar12;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R2) = 0.5 * taubar12;

  /* eom->nu2 * P2 */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R2) = eom->nu2;

  /* - eom->nu2 * Q2 + gammabar22 * Q2 */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2) = - eom->nu2 + gammabar22;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R2) = 0.0;

  /************************************************************************************************************/
  /* Solution I */
  /************************************************************************************************************/

  /* 0.5 * taubar12 * Q2 */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2) = 0.5 * taubar12;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I2) = 0.0;

  /* gammabar12 * Q2 + 0.5 * taubar12 * P2 */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2) = gammabar12;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I2) = 0.5 * taubar12;

  /* eom->nu2 * P2 */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I2) = eom->nu2;

  /* - eom->nu2 * Q2 + gammabar22 * Q2 */
  
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I1) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2) = - eom->nu2 + gammabar22;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I2) = 0.0;
  
  return 0;
}

static gint
_nc_hipert_two_fluids_J_zetaS (_NCM_SUNDIALS_INT_TYPE N, realtype alpha, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHIPertTwoFluidsArg *arg = (NcHIPertTwoFluidsArg *) jac_data;
  const gdouble k = NC_HIPERT (arg->ptf)->k;
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);

  /************************************************************************************************************/
  /* Solution R */
  /************************************************************************************************************/

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R,  NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R,  NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R,  NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = 1.0 / eom->m_zeta;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R,  NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = eom->y;

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_R,     NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_R,     NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_R,     NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = eom->y;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_R,     NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = 1.0 / eom->m_s;

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = - eom->mnu2_zeta;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R, NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = 0.0;

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R,    NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R,    NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = - eom->mnu2_s;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R,    NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R,    NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = 0.0;
  
  /************************************************************************************************************/
  /* Solution I */
  /************************************************************************************************************/

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I,  NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I,  NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I,  NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = 1.0 / eom->m_zeta;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I,  NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = eom->y;

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_I,     NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_I,     NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_I,     NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = eom->y;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_I,     NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = 1.0 / eom->m_s;

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = - eom->mnu2_zeta;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I, NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = 0.0;

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I,    NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I,    NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = - eom->mnu2_s;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I,    NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I,    NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = 0.0;
  
  return 0;
}

static gint
_nc_hipert_two_fluids_J_zetaS_zeta (_NCM_SUNDIALS_INT_TYPE N, realtype alpha, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHIPertTwoFluidsArg *arg = (NcHIPertTwoFluidsArg *) jac_data;
  const gdouble k = NC_HIPERT (arg->ptf)->k;
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);

  /************************************************************************************************************/
  /* Solution R */
  /************************************************************************************************************/

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R,  NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R,  NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R,  NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = 1.0 / eom->m_zeta;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R,  NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = 0.0;

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_R,     NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_R,     NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_R,     NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = eom->y;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_R,     NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = 0.0;

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = - eom->mnu2_zeta;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R, NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = 0.0;

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R,    NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R,    NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R,    NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R,    NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = 0.0;
  
  /************************************************************************************************************/
  /* Solution I */
  /************************************************************************************************************/

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I,  NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I,  NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I,  NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = 1.0 / eom->m_zeta;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I,  NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = 0.0;

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_I,     NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_I,     NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_I,     NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = eom->y;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_I,     NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = 0.0;

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = - eom->mnu2_zeta;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I, NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = 0.0;

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I,    NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I,    NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I,    NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I,    NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = 0.0;
  
  return 0;
}

static gint
_nc_hipert_two_fluids_J_zetaS_S (_NCM_SUNDIALS_INT_TYPE N, realtype alpha, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHIPertTwoFluidsArg *arg = (NcHIPertTwoFluidsArg *) jac_data;
  const gdouble k = NC_HIPERT (arg->ptf)->k;
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);

  /************************************************************************************************************/
  /* Solution R */
  /************************************************************************************************************/

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R,  NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R,  NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R,  NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R,  NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = eom->y;

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_R,     NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_R,     NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_R,     NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_R,     NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = 1.0 / eom->m_s;

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R, NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = 0.0;

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R,    NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R,    NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = - eom->mnu2_s;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R,    NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R,    NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = 0.0;
  
  /************************************************************************************************************/
  /* Solution I */
  /************************************************************************************************************/

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I,  NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I,  NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I,  NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I,  NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = eom->y;

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_I,     NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_I,     NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_I,     NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_I,     NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = 1.0 / eom->m_s;

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I, NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = 0.0;

  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I,    NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I,    NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = - eom->mnu2_s;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I,    NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = 0.0;
  DENSE_ELEM (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I,    NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = 0.0;
  
  return 0;
}

/**
 * nc_hipert_two_fluids_set_init_cond:
 * @ptf: a #NcHIPertTwoFluids
 * @cosmo: a #NcHICosmo
 * @alpha: the log-redshift time
 * @main_mode: main mode
 * @useQP: whether to use the $(Q,\,P)$ system
 * @init_cond: a #NcmVector (size >= 8) containing the initial conditions 
 * 
 * Sets the initial conditions for the two fluids system evolution. 
 * 
 */
void 
nc_hipert_two_fluids_set_init_cond (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, guint main_mode, gboolean useQP, NcmVector *init_cond)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  gint vtype     = useQP ? 1 : 0;
  gint c_vtype   = ptf->useQP ? 1 : 0;
  gint flag;
  guint i;
  ARKRhsFn fE, fI;
  ARKDlsDenseJacFn dfI_dy;

  if (vtype != c_vtype)
    nc_hipert_reset_solver (pert);

  nc_hipert_set_sys_size (pert, NC_HIPERT_ITWO_FLUIDS_VARS_LEN);
  
  pert->alpha0 = alpha;

  if (useQP)
  {
    switch (main_mode)
    {
      case 1:
        fE     = _nc_hipert_two_fluids_f_QP_mode1;
        fI     = _nc_hipert_two_fluids_f_QP_mode2;
        dfI_dy = _nc_hipert_two_fluids_J_QP_mode2;
        break;
      case 2:
        fE     = _nc_hipert_two_fluids_f_QP_mode2;
        fI     = _nc_hipert_two_fluids_f_QP_mode1;
        dfI_dy = _nc_hipert_two_fluids_J_QP_mode1;
        break;
      case 3:
        fE     = _nc_hipert_two_fluids_f_QP;
        fI     = NULL;
        dfI_dy = NULL;
        break;
      case 4:
        fE     = NULL;
        fI     = _nc_hipert_two_fluids_f_QP;
        dfI_dy = _nc_hipert_two_fluids_J_QP;
        break;
      default:
        g_error ("nc_hipert_two_fluids_set_init_cond: Unknown main mode %u.", main_mode);
        break;
    }
  }
  else
  {
    switch (main_mode)
    {
      case 1:
        fE     = _nc_hipert_two_fluids_f_zetaS_zeta;
        fI     = _nc_hipert_two_fluids_f_zetaS_S;
        dfI_dy = _nc_hipert_two_fluids_J_zetaS_S;
        break;
      case 2:
        fE     = _nc_hipert_two_fluids_f_zetaS_S;
        fI     = _nc_hipert_two_fluids_f_zetaS_zeta;
        dfI_dy = _nc_hipert_two_fluids_J_zetaS_zeta;
        break;
      case 3:
        fE     = _nc_hipert_two_fluids_f_zetaS;
        fI     = NULL;
        dfI_dy = NULL;        
        break;
      case 4:
        fE     = NULL;
        fI     = _nc_hipert_two_fluids_f_zetaS;
        dfI_dy = _nc_hipert_two_fluids_J_zetaS;        
        break;
      default:
        g_error ("nc_hipert_two_fluids_set_init_cond: Unknown main mode %u.", main_mode);
        break;
    }
  }

  for (i = 0; i < NC_HIPERT_ITWO_FLUIDS_VARS_LEN; i++)
  {
    NV_Ith_S (pert->y, i) = ncm_vector_get (init_cond, i);
  }
 
  if (!pert->cvode_init)
  {
    flag = ARKodeInit (ptf->arkode, fE, fI, alpha, pert->y);
    NCM_CVODE_CHECK (&flag, "ARKodeInit", 1, );
    
    if (useQP)
    {
      flag = CVodeInit (pert->cvode, &_nc_hipert_two_fluids_f_QP, alpha, pert->y);
      NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
      
      pert->cvode_init = TRUE;      
    }
    else
    {
      flag = CVodeInit (pert->cvode, &_nc_hipert_two_fluids_f_zetaS, alpha, pert->y);
      NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

      pert->cvode_init = TRUE;
    }
  }
  else
  {
    flag = CVodeReInit (pert->cvode, alpha, pert->y);
    NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );

    flag = ARKodeReInit (ptf->arkode, fE, fI, alpha, pert->y);
    NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );
  }
/*
  flag = CVodeSetMaxStep (pert->cvode, 1.0);
  NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1,);
*/
  
  flag = CVodeSStolerances (pert->cvode, pert->reltol, 0.0);
  NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1,);

  flag = ARKodeSStolerances (ptf->arkode, pert->reltol, 0.0);
  NCM_CVODE_CHECK (&flag, "ARKodeSStolerances", 1,);
  
  flag = CVodeSetMaxNumSteps (pert->cvode, G_MAXUINT32);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

  flag = ARKodeSetMaxNumSteps (ptf->arkode, G_MAXUINT32);
  NCM_CVODE_CHECK (&flag, "ARKodeSetMaxNumSteps", 1, );

  flag = CVDense (pert->cvode, NC_HIPERT_ITWO_FLUIDS_VARS_LEN);
  NCM_CVODE_CHECK (&flag, "CVDense", 1, );

  flag = ARKDense (ptf->arkode, NC_HIPERT_ITWO_FLUIDS_VARS_LEN);
  NCM_CVODE_CHECK (&flag, "ARKDense", 1, );

  flag = ARKodeSetLinear (ptf->arkode, 1);
  NCM_CVODE_CHECK (&flag, "ARKodeSetLinear", 1, );
  
  flag = ARKDlsSetDenseJacFn (ptf->arkode, dfI_dy);
  NCM_CVODE_CHECK (&flag, "ARKDlsSetDenseJacFn", 1, );
  
  if (useQP)
  {
    flag = CVDlsSetDenseJacFn (pert->cvode, &_nc_hipert_two_fluids_J_QP);
    NCM_CVODE_CHECK (&flag, "CVDlsSetDenseJacFn", 1, );
  }
  else
  {
    flag = CVDlsSetDenseJacFn (pert->cvode, &_nc_hipert_two_fluids_J_zetaS);
    NCM_CVODE_CHECK (&flag, "CVDlsSetDenseJacFn", 1, );
  }

  switch (main_mode)
  {
    case 1:
    case 2:
      flag = ARKodeSetImEx (ptf->arkode);
      NCM_CVODE_CHECK (&flag, "ARKodeSetImEx", 1, );

      flag = ARKodeSetOrder (ptf->arkode, 5);
      NCM_CVODE_CHECK (&flag, "ARKodeSetOrder", 1, );  

      //flag = ARKodeSetARKTableNum (ptf->arkode, 22, 9);
      //NCM_CVODE_CHECK (&flag, "ARKodeSetIRKTableNum", 1, );  
      break;
    case 3:
      flag = ARKodeSetExplicit (ptf->arkode);
      NCM_CVODE_CHECK (&flag, "ARKodeSetExplicit", 1, );

      flag = ARKodeSetOrder (ptf->arkode, 6);
      NCM_CVODE_CHECK (&flag, "ARKodeSetOrder", 1, );  

      flag = ARKodeSetERKTableNum (ptf->arkode, 10);
      NCM_CVODE_CHECK (&flag, "ARKodeSetERKTableNum", 1, );  
      break;
    case 4:
      flag = ARKodeSetImplicit (ptf->arkode);
      NCM_CVODE_CHECK (&flag, "ARKodeSetImplicit", 1, );

      flag = ARKodeSetIRKTableNum (ptf->arkode, 11 + 11);
      NCM_CVODE_CHECK (&flag, "ARKodeSetIRKTableNum", 1, );  
      break;
    default:
      g_error ("nc_hipert_two_fluids_set_init_cond: Unknown main mode %u.", main_mode);
      break;
  }
  
//  ARKodeSetAdaptivityMethod (ptf->arkode, 5, 1, 0, NULL);
//  NCM_CVODE_CHECK (&flag, "ARKodeSetAdaptivityMethod", 1, );
  
  ptf->useQP = useQP;
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
  NcHIPertTwoFluidsArg *arg = ptf->arg;
  gdouble *yi = NV_DATA_S (pert->y);
  gdouble alpha_i = 0.0;
  gint flag;
  
  arg->cosmo = cosmo;
  arg->ptf   = ptf;

  if (NV_LENGTH_S (pert->y) != NC_HIPERT_ITWO_FLUIDS_VARS_LEN)
    g_error ("nc_hipert_two_fluids_evolve: cannot evolve subsidiary approximated system, use the appropriated evolve function.");
  
  if (TRUE)
  {
    flag = ARKodeSetUserData (ptf->arkode, arg);
    NCM_CVODE_CHECK (&flag, "ARKodeSetUserData", 1, );

    //ARKodeSetDiagnostics (ptf->arkode, stderr);

    flag = ARKode (ptf->arkode, alphaf, pert->y, &alpha_i, CV_NORMAL);
    NCM_CVODE_CHECK (&flag, "CVode[nc_hipert_two_fluids_evolve]", 1, );

    pert->alpha0 = alpha_i;
  }
  else
  {
    flag = CVodeSetUserData (pert->cvode, arg);
    NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

    flag = CVode (pert->cvode, alphaf, pert->y, &alpha_i, CV_NORMAL);
    NCM_CVODE_CHECK (&flag, "CVode[nc_hipert_two_fluids_evolve]", 1, );

    pert->alpha0 = alpha_i;
  }

  if (FALSE)
  {
    gdouble cmp_data[8];
    NcmVector *cmp = ncm_vector_new_data_static (cmp_data, 8, 1);
    NcmVector *ci  = ncm_vector_new_data_static (yi, 8, 1);

    ncm_vector_memcpy (cmp, ci);

    nc_hipert_two_fluids_get_init_cond_QP (ptf, cosmo, pert->alpha0, 2, 
                                           carg (NC_HIPERT_TWO_FLUIDS_QP2A (yi[NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2], yi[NC_HIPERT_ITWO_FLUIDS_VARS_P_R2])), 
                                           cmp);
    ncm_vector_sub (cmp, ci);
    ncm_vector_div (cmp, ci);
    
    printf ("%.5f % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 8.2e\n", 
            pert->alpha0, 
            cmp_data[0], cmp_data[1], cmp_data[2], cmp_data[3], cmp_data[4], cmp_data[5], cmp_data[6], cmp_data[7],
            (1.0 - nc_hipert_two_fluids_get_state_mod (ptf)));

    ncm_vector_free (ci);
    ncm_vector_free (cmp);
  }
}

/**
 * nc_hipert_two_fluids_peek_state:
 * @ptf: a #NcHIPertTwoFluids
 * @cosmo: a #NcHICosmo
 * @alpha: (out): current time
 *
 * Get the current time and values of the numerical solution.
 * 
 * Returns: (transfer none): current solution state.
 */
NcmVector *
nc_hipert_two_fluids_peek_state (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble *alpha)
{
  NcHIPert *pert  = NC_HIPERT (ptf);
  const guint len = NV_LENGTH_S (pert->y);
  guint i;

  if (len == NC_HIPERT_ITWO_FLUIDS_VARS_LEN)
  {
    for (i = 0; i < len; i++)
    {
      ncm_vector_set (ptf->state, i, NV_Ith_S (pert->y, i));
    }
  }
  else if (len == 4)
  {
    const gdouble k = pert->k;
    NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (cosmo), pert->alpha0, k);

    const gdouble Q_R2        = NV_Ith_S (pert->y, 0);
    const gdouble P_R2        = NV_Ith_S (pert->y, 1);

    const gdouble Q_I2        = NV_Ith_S (pert->y, 2);
    const gdouble P_I2        = NV_Ith_S (pert->y, 3);

    const complex double A_R2 = NC_HIPERT_TWO_FLUIDS_QP2A (Q_R2, P_R2);
    const complex double A_I2 = NC_HIPERT_TWO_FLUIDS_QP2A (Q_I2, P_I2);

    const gdouble gammabar11  = eom->gammabar11;
    const gdouble gammabar22  = eom->gammabar22;
    const gdouble gammabar12  = eom->gammabar12;
    const gdouble taubar12    = eom->taubar;
    const gdouble dsigma1     = eom->nu1 - 0.5 * gammabar11;
    const gdouble dsigma2     = eom->nu2 - 0.5 * gammabar22;

    const complex double Lambda = gammabar12 + I * taubar12;

    const complex double A_R1 = (
                                 + A_R2        * ( Lambda     / (2.0 * (dsigma1 - dsigma2)) + gammabar11    * gammabar12 / (8.0 * dsigma1 * (dsigma1 + dsigma2)) ) 
                                 - conj (A_R2) * ( gammabar12 / (2.0 * (dsigma1 + dsigma2)) + conj (Lambda) * gammabar11 / (8.0 * dsigma1 * (dsigma1 - dsigma2)) )
                                 ) / (1.0 - gsl_pow_2 (gammabar11 / (4.0 * dsigma1)));
    const complex double A_I1 = (
                                 + A_I2        * ( Lambda     / (2.0 * (dsigma1 - dsigma2)) + gammabar11    * gammabar12 / (8.0 * dsigma1 * (dsigma1 + dsigma2)) ) 
                                 - conj (A_I2) * ( gammabar12 / (2.0 * (dsigma1 + dsigma2)) + conj (Lambda) * gammabar11 / (8.0 * dsigma1 * (dsigma1 - dsigma2)) )
                                 ) / (1.0 - gsl_pow_2 (gammabar11 / (4.0 * dsigma1)));

    const gdouble Q_R1        = NC_HIPERT_TWO_FLUIDS_A2Q (A_R1);
    const gdouble P_R1        = NC_HIPERT_TWO_FLUIDS_A2P (A_R1);

    const gdouble Q_I1        = NC_HIPERT_TWO_FLUIDS_A2Q (A_I1);
    const gdouble P_I1        = NC_HIPERT_TWO_FLUIDS_A2P (A_I1);

    ncm_vector_set (ptf->state, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1, Q_R1);
    ncm_vector_set (ptf->state, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1, P_R1);

    ncm_vector_set (ptf->state, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2, Q_R2);
    ncm_vector_set (ptf->state, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2, P_R2);

    ncm_vector_set (ptf->state, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1, Q_I1);
    ncm_vector_set (ptf->state, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1, P_I1);

    ncm_vector_set (ptf->state, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2, Q_I2);
    ncm_vector_set (ptf->state, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2, P_I2);
  }
  else
  {
    g_assert_not_reached ();
  }
    

  alpha[0] = pert->alpha0;

  return ptf->state;
}

/**
 * nc_hipert_two_fluids_set_init_cond_mode1sub:
 * @ptf: a #NcHIPertTwoFluids
 * @cosmo: a #NcHICosmo
 * @alpha: the log-redshift time
 * @init_cond: a #NcmVector (size >= 4) containing the initial conditions 
 * 
 * Sets the initial conditions for the two fluids system evolution. 
 * 
 */
void 
nc_hipert_two_fluids_set_init_cond_mode1sub (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, NcmVector *init_cond)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  gint flag;
  const guint sys_size = 4;

  nc_hipert_set_sys_size (pert, sys_size);

  pert->alpha0 = alpha;

  NV_Ith_S (pert->y, 0) = ncm_vector_get (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2);
  NV_Ith_S (pert->y, 1) = ncm_vector_get (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2);

  NV_Ith_S (pert->y, 2) = ncm_vector_get (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2);
  NV_Ith_S (pert->y, 3) = ncm_vector_get (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2);

  if (!pert->cvode_init)
  {
    flag = CVodeInit (pert->cvode, &_nc_hipert_two_fluids_f_QP_mode1sub, alpha, pert->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
    pert->cvode_init = TRUE;
  }
  else
  {
    flag = CVodeReInit (pert->cvode, alpha, pert->y);
    NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );
  }

  flag = CVodeSetMaxStep (pert->cvode, 1.0);
  NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1,);
    
  flag = CVodeSStolerances (pert->cvode, pert->reltol, 0.0);
  NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1,);

  flag = CVodeSetMaxNumSteps (pert->cvode, G_MAXUINT32);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

  flag = CVDense (pert->cvode, sys_size);
  NCM_CVODE_CHECK (&flag, "CVDense", 1, );

/*
  flag = CVDlsSetDenseJacFn (pert->cvode, &_nc_hipert_two_fluids_J_QP);
  NCM_CVODE_CHECK (&flag, "CVDlsSetDenseJacFn", 1, );
*/  
}

/**
 * nc_hipert_two_fluids_evolve_mode1sub:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alphaf: the final log-redshift time.
 * 
 * Evolve the system until @alphaf.
 * 
 */
void 
nc_hipert_two_fluids_evolve_mode1sub (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alphaf)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  NcHIPertTwoFluidsArg *arg = ptf->arg;
  gdouble alpha_i = 0.0;
  gint flag;
  
  arg->cosmo = cosmo;
  arg->ptf   = ptf;

  if (NV_LENGTH_S (pert->y) != 4)
    g_error ("nc_hipert_two_fluids_evolve: cannot evolve full system, use the appropriated evolve function.");

  flag = CVodeSetUserData (pert->cvode, arg);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

  flag = CVode (pert->cvode, alphaf, pert->y, &alpha_i, CV_NORMAL);
  NCM_CVODE_CHECK (&flag, "CVode[nc_hipert_two_fluids_evolve]", 1, );

  pert->alpha0 = alpha_i;
}

/**
 * nc_hipert_two_fluids_get_state_mod:
 * @ptf: a #NcHIPertTwoFluids
 * 
 * Get the current module for the solution.
 * 
 * Returns: state module.
 */
gdouble
nc_hipert_two_fluids_get_state_mod (NcHIPertTwoFluids *ptf)
{
  NcHIPert *pert = NC_HIPERT (ptf);

  const complex double zeta  = NV_Ith_S (pert->y, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  + I * NV_Ith_S (pert->y, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I);
  const complex double S     = NV_Ith_S (pert->y, NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     + I * NV_Ith_S (pert->y, NC_HIPERT_ITWO_FLUIDS_VARS_S_I);
  const complex double Pzeta = NV_Ith_S (pert->y, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) + I * NV_Ith_S (pert->y, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I);
  const complex double PS    = NV_Ith_S (pert->y, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    + I * NV_Ith_S (pert->y, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I);

/*
   printf ("% 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g\n",
          cabs (zeta), cabs (Pzeta), carg (zeta) / M_PI, carg (Pzeta) / M_PI,
          cabs (S), cabs (PS), carg (S) / M_PI, carg (PS) / M_PI
          );
*/  
  return 
    + 2.0 * cabs (zeta) * cabs (Pzeta) * sin (carg (zeta) - carg (Pzeta)) 
    + 2.0 * cabs (S) * cabs (PS) * sin (carg (S) - carg (PS));
}

static gdouble 
_nc_hipert_two_fluids_cross_time_mode1main_root (gdouble alpha, gpointer userdata)
{
  NcHIPertTwoFluidsArg *arg = (NcHIPertTwoFluidsArg *) userdata;
  
  const gdouble k = NC_HIPERT (arg->ptf)->k;
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);

  const gdouble p = fabs (eom->gammabar11 / eom->nu1);

  return log (p / arg->prec);
}

static gdouble 
_nc_hipert_two_fluids_cross_time_mode2main_root (gdouble alpha, gpointer userdata)
{
  NcHIPertTwoFluidsArg *arg = (NcHIPertTwoFluidsArg *) userdata;
  
  const gdouble k = NC_HIPERT (arg->ptf)->k;
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);

  const gdouble p = fabs (eom->gammabar22 / eom->nu2);

  return log (p / arg->prec);
}

static gdouble 
_nc_hipert_two_fluids_cross_time_mode1sub_root (gdouble alpha, gpointer userdata)
{
  NcHIPertTwoFluidsArg *arg = (NcHIPertTwoFluidsArg *) userdata;
  
  const gdouble k = NC_HIPERT (arg->ptf)->k;
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);

  gdouble p = 0.0;

  p = sqrt (
            gsl_pow_2 (eom->gammabar11 / eom->nu1) +
            gsl_pow_2 (eom->gammabar22 / eom->nu1) +
            gsl_pow_2 (eom->gammabar12 / eom->nu1) +
            gsl_pow_2 (eom->taubar     / eom->nu1)
            );

  return log (p / arg->prec);
}

static gdouble 
_nc_hipert_two_fluids_cross_time_mode2sub_root (gdouble alpha, gpointer userdata)
{
  NcHIPertTwoFluidsArg *arg = (NcHIPertTwoFluidsArg *) userdata;
  
  const gdouble k = NC_HIPERT (arg->ptf)->k;
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);

  gdouble p = 0.0;

  p = sqrt (
            gsl_pow_2 (eom->gammabar11 / eom->nu2) +
            gsl_pow_2 (eom->gammabar22 / eom->nu2) +
            gsl_pow_2 (eom->gammabar12 / eom->nu2) +
            gsl_pow_2 (eom->taubar     / eom->nu2)
            );

  return log (p / arg->prec);
}

/**
 * nc_hipert_two_fluids_get_cross_time:
 * @ptf: a #NcHIPertTwoFluids
 * @cosmo: a #NcHICosmo
 * @cross: a #NcHIPertTwoFluidsCross
 * @alpha_i: initial try
 * @prec: precision
 * 
 * Get the initial time where the approximate solution is valid 
 * within precision @prec.
 * 
 * Returns: initial time $\alpha_i$.
 */
gdouble
nc_hipert_two_fluids_get_cross_time (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, NcHIPertTwoFluidsCross cross, gdouble alpha_i, gdouble prec)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  NcHIPertTwoFluidsArg *arg = ptf->arg;
  gint status;
  gint iter = 0, max_iter = 1000000;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;
  gdouble alpha0, alpha1, alpha = alpha_i;

  switch (cross)
  {
    case NC_HIPERT_TWO_FLUIDS_CROSS_MODE1MAIN:
      F.function = &_nc_hipert_two_fluids_cross_time_mode1main_root;
      break;
    case NC_HIPERT_TWO_FLUIDS_CROSS_MODE2MAIN:
      F.function = &_nc_hipert_two_fluids_cross_time_mode2main_root;
      break;
    case NC_HIPERT_TWO_FLUIDS_CROSS_MODE1SUB:
      F.function = &_nc_hipert_two_fluids_cross_time_mode1sub_root;
      break;
    case NC_HIPERT_TWO_FLUIDS_CROSS_MODE2SUB:
      F.function = &_nc_hipert_two_fluids_cross_time_mode2sub_root;
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  F.params   = arg;

  arg->cosmo = cosmo;
  arg->ptf   = ptf;
  arg->prec  = prec;

  alpha0 = alpha_i;

  while (GSL_FN_EVAL (&F, alpha0) > 0.0)
  {
    alpha0 -= 10.0;
  }

  alpha1 = alpha_i;
  while (GSL_FN_EVAL (&F, alpha1) < 0.0)
  {
    alpha1 += 10.0;
  }

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, alpha0, alpha1);
  
  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    if (status)
    {
      g_warning ("%s", gsl_strerror (status));
      break;
    }

    alpha  = gsl_root_fsolver_root (s);
    alpha0 = gsl_root_fsolver_x_lower (s);
    alpha1 = gsl_root_fsolver_x_upper (s);

    status = gsl_root_test_residual (GSL_FN_EVAL (&F, alpha), pert->reltol);
    
    if (status == GSL_CONTINUE) 
      status = gsl_root_test_interval (alpha0, alpha1, 0, prec);
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  if (status)
    g_warning ("%s", gsl_strerror (status));
  
  if (iter >= max_iter)
    g_warning ("nc_hipert_two_fluids_get_cross_time: maximum number of interations reached.");

  gsl_root_fsolver_free (s);

  return alpha;
}

