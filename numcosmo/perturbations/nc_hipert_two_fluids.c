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
 * @title: Two fluids Perturbation Object
 * @short_description: Perturbation object for a two fluids system 
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
 * where $\bar{\rho} + \bar{p}$ is the background total energy density plus pressure defined by nc_hicosmo_rhopp(),
 * $E^2 = H^2/H_0^2$ is the adimensional Hubble function squared (nc_hicosmo_E2()), $c_s^2$ the speed of sound (nc_hicosmo_cs2()),
 * $N$ is the lapse function that in this case (using $\alpha$ as time variable) is $N \equiv \vert{}E\vert^{-1}$, $\rho_\text{crit0}$
 * is the critical density today defined by $\rho_\text{crit0} \equiv 3H_0^2/\kappa$ and $$\Delta_\bar{K} \equiv \frac{k^2}{k^2 + \Omega_k}.$$
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline_cubic_notaknot.h"
#include "perturbations/nc_hipert_two_fluids.h"
#include "perturbations/nc_hipert_iadiab.h"
#include "perturbations/nc_hipert_itwo_fluids.h"

#include "math/cvode_util.h"

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_band.h>
#include <nvector/nvector_serial.h> 
#include <gsl/gsl_roots.h>

G_DEFINE_TYPE (NcHIPertTwoFluids, nc_hipert_two_fluids, NC_TYPE_HIPERT);

static gdouble _nc_hipert_two_fluids_wkb_zeta_phase (gdouble y, gdouble x, gpointer userdata);
static gdouble _nc_hipert_two_fluids_wkb_S_phase (gdouble y, gdouble x, gpointer userdata);

static void
nc_hipert_two_fluids_init (NcHIPertTwoFluids *ptf)
{
  NcmSpline *s1      = ncm_spline_cubic_notaknot_new ();
  NcmSpline *s2      = ncm_spline_cubic_notaknot_new ();

  ptf->wkb_zeta_phase = ncm_ode_spline_new (s1, &_nc_hipert_two_fluids_wkb_zeta_phase);
  ptf->wkb_S_phase    = ncm_ode_spline_new (s2, &_nc_hipert_two_fluids_wkb_S_phase);

  ptf->zeta_wkb_prepared = FALSE;
  ptf->S_wkb_prepared    = FALSE;
  
  ncm_ode_spline_set_abstol (ptf->wkb_zeta_phase, 0.0);
  ncm_ode_spline_set_abstol (ptf->wkb_S_phase, 0.0);
  
  ncm_spline_free (s1);
  ncm_spline_free (s2);
}

static void
nc_hipert_two_fluids_finalize (GObject *object)
{
  NcHIPertTwoFluids *ptf = NC_HIPERT_TWO_FLUIDS (object);
  
  ncm_ode_spline_clear (&ptf->wkb_zeta_phase);
  ncm_ode_spline_clear (&ptf->wkb_S_phase);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_two_fluids_parent_class)->finalize (object);
}

static void
nc_hipert_two_fluids_class_init (NcHIPertTwoFluidsClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);


  object_class->finalize = nc_hipert_two_fluids_finalize;
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

typedef struct _NcHIPertTwoFluidsWKBPhaseArg
{
  NcHICosmo *cosmo;
  gdouble k;
} NcHIPertTwoFluidsWKBPhaseArg;

static gdouble 
_nc_hipert_two_fluids_wkb_zeta_phase (gdouble y, gdouble x, gpointer userdata)
{
  NcHIPertTwoFluidsWKBPhaseArg *arg = (NcHIPertTwoFluidsWKBPhaseArg *) userdata;
  const gdouble nuA2 = nc_hipert_iadiab_nuA2 (NC_HIPERT_IADIAB (arg->cosmo), x, arg->k);
  const gdouble nuA = sqrt (nuA2);
  return nuA;
}

static gdouble 
_nc_hipert_two_fluids_wkb_S_phase (gdouble y, gdouble x, gpointer userdata)
{
  NcHIPertTwoFluidsWKBPhaseArg *arg = (NcHIPertTwoFluidsWKBPhaseArg *) userdata;
  
  const gdouble nuB2 = nc_hipert_itwo_fluids_nuB2 (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), x, arg->k);
  const gdouble nuB = sqrt (nuB2);
  return nuB;
}

/**
 * nc_hipert_two_fluids_prepare_wkb_zeta:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha_i: initial log-redshift time.
 * @alpha_f: final log-redshift time.
 * 
 * Prepare the zeta component of the object for WKB calculations using the cosmology @cosmo.
 * 
 */
void 
nc_hipert_two_fluids_prepare_wkb_zeta (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha_i, gdouble alpha_f)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  if (!pert->prepared)
  {
    ptf->zeta_wkb_prepared = FALSE;
    ptf->S_wkb_prepared = FALSE;
    pert->prepared = TRUE;
  }

  if (!ptf->zeta_wkb_prepared)
  {
    NcHIPertTwoFluidsWKBPhaseArg arg;

    arg.cosmo = cosmo;
    arg.k     = pert->k;

    ncm_ode_spline_set_interval (ptf->wkb_zeta_phase, 0.25 * M_PI, alpha_i, alpha_f);
    ncm_ode_spline_prepare (ptf->wkb_zeta_phase, &arg);

    ptf->zeta_wkb_prepared = TRUE;
  }
}

/**
 * nc_hipert_two_fluids_prepare_wkb_S:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha_i: initial log-redshift time.
 * @alpha_f: final log-redshift time.
 * 
 * Prepare the zeta component of the object for WKB calculations using the cosmology @cosmo.
 * 
 */
void 
nc_hipert_two_fluids_prepare_wkb_S (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha_i, gdouble alpha_f)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  if (!pert->prepared)
  {
    ptf->zeta_wkb_prepared = FALSE;
    ptf->S_wkb_prepared = FALSE;
    pert->prepared = TRUE;
  }

  if (!ptf->S_wkb_prepared)
  {
    NcHIPertTwoFluidsWKBPhaseArg arg;

    arg.cosmo = cosmo;
    arg.k     = pert->k;

    ncm_ode_spline_set_interval (ptf->wkb_S_phase, 0.25 * M_PI, alpha_i, alpha_f);
    ncm_ode_spline_prepare (ptf->wkb_S_phase, &arg);
    
    ptf->S_wkb_prepared = TRUE;
  }
}

/**
 * nc_hipert_two_fluids_wkb_zeta:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha: the log-redshift time.
 * @Re_zeta: (out caller-allocates): Real part of the wkb solution.
 * @Im_zeta: (out caller-allocates): Imaginary part of the wkb solution.
 * 
 * Computes the WKB solution $\zeta_\text{WKB}$ for the mode $k$ at the time $\alpha$. 
 * 
 */
void
nc_hipert_two_fluids_wkb_zeta (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, gdouble *Re_zeta, gdouble *Im_zeta)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  complex double zeta;
  gdouble int_nuA;
  const gdouble nuA2 = nc_hipert_iadiab_nuA2 (NC_HIPERT_IADIAB (cosmo), alpha, pert->k);
  const gdouble nuA = sqrt (nuA2); 

  g_assert (pert->prepared && ptf->zeta_wkb_prepared);

  int_nuA = ncm_spline_eval (ptf->wkb_zeta_phase->s, alpha); 
  
  zeta = cexp (-I * int_nuA) / sqrt (2.0 * eom->mzeta * nuA);

  *Re_zeta = creal (zeta);
  *Im_zeta = cimag (zeta);
}

/**
 * nc_hipert_two_fluids_wkb_zeta_Pzeta:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha: the log-redshift time.
 * @Re_zeta: (out caller-allocates): Real part of the wkb solution.
 * @Im_zeta: (out caller-allocates): Imaginary part of the wkb solution.
 * @Re_Pzeta: (out caller-allocates): Real part of the wkb solution momentum.
 * @Im_Pzeta: (out caller-allocates): Imaginary part of the wkb solution momentum.
 * 
 * Computes the WKB solution $\zeta_\text{WKB}$ and its momentum for the mode $k$ at the time $\alpha$. 
 * 
 */
void
nc_hipert_two_fluids_wkb_zeta_Pzeta (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, gdouble *Re_zeta, gdouble *Im_zeta, gdouble *Re_Pzeta, gdouble *Im_Pzeta)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  complex double zeta, Pzeta;
  gdouble int_nuA;
  gdouble dmzetanuA_nuA;
  const gdouble nuA2 = nc_hipert_iadiab_nuA2 (NC_HIPERT_IADIAB (cosmo), alpha, pert->k);
  const gdouble nuA = sqrt (nuA2); 

  g_assert (pert->prepared && ptf->zeta_wkb_prepared);

  int_nuA = ncm_spline_eval (ptf->wkb_zeta_phase->s, alpha);
  dmzetanuA_nuA = nc_hipert_iadiab_dmzetanuA_nuA (NC_HIPERT_IADIAB (cosmo), alpha, pert->k);

  zeta = cexp (-I * int_nuA) / sqrt (2.0 * eom->mzeta * nuA);

  Pzeta = -I * cexp (-I * int_nuA) * sqrt (0.5 * eom->mzeta * nuA) - 0.5 * dmzetanuA_nuA * zeta;

  *Re_zeta = creal (zeta);
  *Im_zeta = cimag (zeta);

  *Re_Pzeta = creal (Pzeta);
  *Im_Pzeta = cimag (Pzeta);  
}

/**
 * nc_hipert_two_fluids_wkb_Q:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha: the log-redshift time.
 * @Re_Q: (out caller-allocates): Real part of the wkb solution.
 * @Im_Q: (out caller-allocates): Imaginary part of the wkb solution.
 * 
 * Computes the WKB solution $Q_\text{WKB}$ for the mode $k$ at the time $\alpha$. 
 * 
 */
void
nc_hipert_two_fluids_wkb_Q (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, gdouble *Re_Q, gdouble *Im_Q)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  complex double Q;
  gdouble int_nuB;
  const gdouble nuB2 = nc_hipert_itwo_fluids_nuB2 (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  const gdouble nuB = sqrt (nuB2); 

  g_assert (pert->prepared && ptf->S_wkb_prepared);

  int_nuB = ncm_spline_eval (ptf->wkb_S_phase->s, alpha); 
  
  Q = cexp (-I * int_nuB) / sqrt (2.0 * eom->mS * nuB);

  *Re_Q = creal (Q);
  *Im_Q = cimag (Q);
}

/**
 * nc_hipert_two_fluids_wkb_Q_PQ:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha: the log-redshift time.
 * @Re_Q: (out caller-allocates): Real part of the wkb solution.
 * @Im_Q: (out caller-allocates): Imaginary part of the wkb solution.
 * @Re_PQ: (out caller-allocates): Real part of the wkb solution momentum.
 * @Im_PQ: (out caller-allocates): Imaginary part of the wkb solution momentum.
 * 
 * Computes the WKB solution $Q_\text{WKB}$ and its momentum for the mode $k$ at the time $\alpha$. 
 * 
 */
void
nc_hipert_two_fluids_wkb_Q_PQ (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, gdouble *Re_Q, gdouble *Im_Q, gdouble *Re_PQ, gdouble *Im_PQ)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  complex double Q, PQ;
  gdouble int_nuB;
  gdouble dmSnuB_nuB;
  const gdouble nuB2 = nc_hipert_itwo_fluids_nuB2 (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  const gdouble nuB = sqrt (nuB2); 

  g_assert (pert->prepared && ptf->S_wkb_prepared);

  int_nuB = ncm_spline_eval (ptf->wkb_S_phase->s, alpha);

  dmSnuB_nuB = nc_hipert_itwo_fluids_dmSnuB_nuB (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);

  Q = cexp (-I * int_nuB) / sqrt (2.0 * eom->mS * nuB);

  PQ = -I * cexp (-I * int_nuB) * sqrt (0.5 * eom->mS * nuB) - 0.5 * dmSnuB_nuB * Q;

  *Re_Q = creal (Q);
  *Im_Q = cimag (Q);

  *Re_PQ = creal (PQ);
  *Im_PQ = cimag (PQ);  
}

/**
 * nc_hipert_two_fluids_wkb_full_zeta:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha: the log-redshift time.
 * @vars: (inout) (array fixed-size=8) (element-type double): Perturbations variables conforming to #NcHIPertTwoFluidsVars.
 * 
 * Computes the complete WKB solution their momenta for the mode $k$ at the time $\alpha$ for the principal mode $\zeta$. 
 * 
 */
void
nc_hipert_two_fluids_wkb_full_zeta (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, gdouble **vars)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  complex double zeta, Pzeta, Q, PQ, C;
  gdouble int_nuA;
  const gdouble nuA2 = nc_hipert_iadiab_nuA2 (NC_HIPERT_IADIAB (cosmo), alpha, pert->k);
  const gdouble nuB2 = nc_hipert_itwo_fluids_nuB2 (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  const gdouble nuA = sqrt (nuA2);
  const gdouble nuB = sqrt (nuB2);
  const gdouble k = pert->k;
  const gdouble dmSnuB_nuB    = nc_hipert_itwo_fluids_dmSnuB_nuB (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  const gdouble dmzetanuA_nuA = nc_hipert_iadiab_dmzetanuA_nuA (NC_HIPERT_IADIAB (cosmo), alpha, pert->k);
  const gdouble LA = 0.5 * k * dmzetanuA_nuA / (eom->mzeta * nuA);
  const gdouble LB = 0.5 * k * dmSnuB_nuB / (eom->mS * nuB);
  const gdouble YAB = eom->Y * sqrt (nuA * nuB * eom->mzeta * eom->mS) / k;
  const gdouble k2 = k * k;
  
  g_assert (pert->prepared && ptf->S_wkb_prepared && ptf->zeta_wkb_prepared);

  int_nuA = ncm_spline_eval (ptf->wkb_zeta_phase->s, alpha);

  zeta = cexp (-I * int_nuA) / sqrt (2.0 * eom->mzeta * nuA);

  Pzeta = -I * cexp (-I * int_nuA) * sqrt (0.5 * eom->mzeta * nuA) - 0.5 * dmzetanuA_nuA * zeta;

  C = YAB * (k - I * LA) * cexp (-I * int_nuA) / (sqrt (2.0 * k) * k * (nuA * nuA - nuB * nuB));
  
  Q = sqrt (k / (eom->mS * nuB)) * (k * nuA + I * LB * nuB) * C;
  PQ = - I * sqrt ((eom->mS * nuB) / k) * (k2 + LB * LB) * nuB * C;

  vars[0][NC_HIPERT_TWO_FLUIDS_RE_ZETA] = creal (zeta);
  vars[0][NC_HIPERT_TWO_FLUIDS_IM_ZETA] = cimag (zeta);

  vars[0][NC_HIPERT_TWO_FLUIDS_RE_PZETA] = creal (Pzeta);
  vars[0][NC_HIPERT_TWO_FLUIDS_IM_PZETA] = cimag (Pzeta);

  vars[0][NC_HIPERT_TWO_FLUIDS_RE_Q] = creal (Q);
  vars[0][NC_HIPERT_TWO_FLUIDS_IM_Q] = cimag (Q);

  vars[0][NC_HIPERT_TWO_FLUIDS_RE_PQ] = creal (PQ);
  vars[0][NC_HIPERT_TWO_FLUIDS_IM_PQ] = cimag (PQ);  
}

/**
 * nc_hipert_two_fluids_wkb_full_Q:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha: the log-redshift time.
 * @vars: (inout) (array fixed-size=8) (element-type double): Perturbations variables conforming to #NcHIPertTwoFluidsVars.
 * 
 * Computes the complete WKB solution their momenta for the mode $k$ at the time $\alpha$ for the principal mode $Q$. 
 * 
 */
void
nc_hipert_two_fluids_wkb_full_Q (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, gdouble **vars)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  complex double zeta, Pzeta, Q, PQ, C;
  gdouble int_nuB;
  const gdouble nuA2 = nc_hipert_iadiab_nuA2 (NC_HIPERT_IADIAB (cosmo), alpha, pert->k);
  const gdouble nuB2 = nc_hipert_itwo_fluids_nuB2 (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  const gdouble nuA  = sqrt (nuA2);
  const gdouble nuB  = sqrt (nuB2);
  const gdouble k    = pert->k;
  const gdouble dmSnuB_nuB    = nc_hipert_itwo_fluids_dmSnuB_nuB (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  const gdouble dmzetanuA_nuA = nc_hipert_iadiab_dmzetanuA_nuA (NC_HIPERT_IADIAB (cosmo), alpha, pert->k);
  const gdouble LA   = 0.5 * k * dmzetanuA_nuA / (eom->mzeta * nuA);
  const gdouble LB   = 0.5 * k * dmSnuB_nuB / (eom->mS * nuB);
  const gdouble YAB  = eom->Y * sqrt (nuA * nuB * eom->mzeta * eom->mS) / k;
  const gdouble k2   = k * k;
  
  g_assert (pert->prepared && ptf->S_wkb_prepared && ptf->zeta_wkb_prepared);

  int_nuB = ncm_spline_eval (ptf->wkb_S_phase->s, alpha);

  Q = cexp (-I * int_nuB) / sqrt (2.0 * eom->mS * nuB);

  PQ = -I * cexp (-I * int_nuB) * sqrt (0.5 * eom->mS * nuB) - 0.5 * dmSnuB_nuB * Q;

  C = YAB * (k - I * LB) * cexp (-I * int_nuB) / (sqrt (2.0 * k) * k * (nuB * nuB - nuA * nuA));

  zeta  = sqrt (k / (eom->mzeta * nuA)) * (k * nuB + I * LA * nuA) * C;
  Pzeta = - I * sqrt ((eom->mzeta * nuA) / k) * (k2 + LA * LA) * nuA * C;

/*  printf ("# WKB % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", alpha, eom->mzeta, nuA2, eom->mS, nuB2, eom->Y);*/

  vars[0][NC_HIPERT_TWO_FLUIDS_RE_ZETA]  = creal (zeta);
  vars[0][NC_HIPERT_TWO_FLUIDS_IM_ZETA]  = cimag (zeta);

  vars[0][NC_HIPERT_TWO_FLUIDS_RE_PZETA] = creal (Pzeta);
  vars[0][NC_HIPERT_TWO_FLUIDS_IM_PZETA] = cimag (Pzeta);

  vars[0][NC_HIPERT_TWO_FLUIDS_RE_Q]     = creal (Q);
  vars[0][NC_HIPERT_TWO_FLUIDS_IM_Q]     = cimag (Q);

  vars[0][NC_HIPERT_TWO_FLUIDS_RE_PQ]    = creal (PQ);
  vars[0][NC_HIPERT_TWO_FLUIDS_IM_PQ]    = cimag (PQ);  
}

static gdouble 
_nc_hipert_two_fluids_wkb_nuA2 (gdouble x, gpointer userdata)
{
  NcHIPertTwoFluidsWKBPhaseArg *arg = (NcHIPertTwoFluidsWKBPhaseArg *) userdata;
  const gdouble nuA2 = nc_hipert_iadiab_nuA2 (NC_HIPERT_IADIAB (arg->cosmo), x, arg->k);
  return nuA2;
}

static gdouble 
_nc_hipert_two_fluids_wkb_nuB2 (gdouble x, gpointer userdata)
{
  NcHIPertTwoFluidsWKBPhaseArg *arg = (NcHIPertTwoFluidsWKBPhaseArg *) userdata;
  const gdouble nuB2 = nc_hipert_itwo_fluids_nuB2 (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), x, arg->k);
  return nuB2;
}

static gint
_nc_hipert_two_fluids_wkb_nu2_root (NcHIPertTwoFluids *ptf, gdouble (*f) (gdouble, gpointer), NcHICosmo *cosmo, gdouble alpha0, gdouble *alpha)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  gint status;
  gint iter = 0, max_iter = 1000000;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;
  gdouble prec = pert->reltol, alpha1 = *alpha;
  NcHIPertTwoFluidsWKBPhaseArg arg;

  arg.cosmo = cosmo;
  arg.k     = pert->k;

  F.function = f;
  F.params   = &arg;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, alpha0, alpha1);

  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    if (status)
    {
      gsl_root_fsolver_free (s);
      return status;
    }

    *alpha = gsl_root_fsolver_root (s);
    alpha0 = gsl_root_fsolver_x_lower (s);
    alpha1 = gsl_root_fsolver_x_upper (s);

    status = gsl_root_test_interval (alpha0, alpha1, 0, prec);
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);

  return 0;
}

/**
 * nc_hipert_two_fluids_wkb_zeta_maxtime:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha0: the initial log-redshift time.
 * @alpha1: the final log-redshift time.
 * 
 * Search for the root of $\nu_A^2$ between $\alpha_0$ and $\alpha_1$. 
 * 
 * Returns: the root of $\nu_A^2$ between $\alpha_0$ and $\alpha_1$ or NaN if not found.
 */
gdouble 
nc_hipert_two_fluids_wkb_zeta_maxtime (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha0, gdouble alpha1)
{
  gdouble alpha = alpha1;
  guint ret = _nc_hipert_two_fluids_wkb_nu2_root (ptf, &_nc_hipert_two_fluids_wkb_nuA2, cosmo, alpha0, &alpha);

  if (ret == 0)
    return alpha;
  else
    return GSL_NAN;
}

/**
 * nc_hipert_two_fluids_wkb_S_maxtime:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha0: the initial log-redshift time.
 * @alpha1: the final log-redshift time.
 * 
 * Search for the root of $\nu_B^2$ between $\alpha_0$ and $\alpha_1$. 
 * 
 * Returns: the root of $\nu_B^2$ between $\alpha_0$ and $\alpha_1$ or NaN if not found.
 */
gdouble 
nc_hipert_two_fluids_wkb_S_maxtime (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha0, gdouble alpha1)
{
  gdouble alpha = alpha1;
  guint ret = _nc_hipert_two_fluids_wkb_nu2_root (ptf, &_nc_hipert_two_fluids_wkb_nuB2, cosmo, alpha0, &alpha);

  if (ret == 0)
    return alpha;
  else
    return GSL_NAN;
}

typedef struct _NcHIPertTwoFluidsArg
{
  NcHICosmo *cosmo;
  NcHIPertTwoFluids *ptf;
} NcHIPertTwoFluidsArg;

static gint
_nc_hipert_two_fluids_f (realtype tau, N_Vector y, N_Vector ydot, gpointer f_data)
{
  const gdouble alpha = NC_HIPERT_TWO_FLUIDS_TAU2ALPHA (tau);
  NcHIPertTwoFluidsArg *arg = (NcHIPertTwoFluidsArg *) f_data;
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, NC_HIPERT (arg->ptf)->k);
  const gdouble dadt = NC_HIPERT_TWO_FLUIDS_DALPHADTAU (alpha); 
/*
  guint i;
  printf ("% 20.15g", alpha);
  for (i = 0; i < 8; i++)
    printf ("% 15.8g ", NV_Ith_S (y, i));
  printf ("\n");
*/
/*
  printf ("% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", x, 
          NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_RE_PZETA) / eom->mzeta, 
          eom->Y * NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_RE_PQ),
          NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_RE_PQ) / eom->mS, 
          eom->Y * NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_RE_PZETA));
*/  
  NV_Ith_S (ydot, NC_HIPERT_TWO_FLUIDS_RE_ZETA)  = (NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_RE_PZETA) / eom->mzeta + eom->Y * NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_RE_PQ)) * dadt;
  NV_Ith_S (ydot, NC_HIPERT_TWO_FLUIDS_IM_ZETA)  = (NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_IM_PZETA) / eom->mzeta + eom->Y * NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_IM_PQ)) * dadt;

  NV_Ith_S (ydot, NC_HIPERT_TWO_FLUIDS_RE_Q)     = (NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_RE_PQ) / eom->mS + eom->Y * NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_RE_PZETA)) * dadt;
  NV_Ith_S (ydot, NC_HIPERT_TWO_FLUIDS_IM_Q)     = (NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_IM_PQ) / eom->mS + eom->Y * NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_IM_PZETA)) * dadt;

  NV_Ith_S (ydot, NC_HIPERT_TWO_FLUIDS_RE_PZETA) = (- eom->mzeta * eom->nuzeta2 * NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_RE_ZETA)) * dadt;
  NV_Ith_S (ydot, NC_HIPERT_TWO_FLUIDS_IM_PZETA) = (- eom->mzeta * eom->nuzeta2 * NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_IM_ZETA)) * dadt;

  NV_Ith_S (ydot, NC_HIPERT_TWO_FLUIDS_RE_PQ)    = (- eom->mS * eom->nuS2 * NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_RE_Q)) * dadt;
  NV_Ith_S (ydot, NC_HIPERT_TWO_FLUIDS_IM_PQ)    = (- eom->mS * eom->nuS2 * NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_IM_Q)) * dadt;

  return 0;
}

static gint
_nc_hipert_two_fluids_J (_NCM_SUNDIALS_INT_TYPE N, realtype tau, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  const gdouble alpha = NC_HIPERT_TWO_FLUIDS_TAU2ALPHA (tau);
  NcHIPertTwoFluidsArg *arg = (NcHIPertTwoFluidsArg *) jac_data;
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, NC_HIPERT (arg->ptf)->k);
  const gdouble dadt = NC_HIPERT_TWO_FLUIDS_DALPHADTAU (alpha);

  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_RE_ZETA,  NC_HIPERT_TWO_FLUIDS_RE_PZETA) = (1.0 / eom->mzeta) * dadt;
  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_IM_ZETA,  NC_HIPERT_TWO_FLUIDS_IM_PZETA) = (1.0 / eom->mzeta) * dadt;
  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_RE_PZETA, NC_HIPERT_TWO_FLUIDS_RE_ZETA) = (- eom->mzeta * eom->nuzeta2) * dadt;
  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_IM_PZETA, NC_HIPERT_TWO_FLUIDS_IM_ZETA) = (- eom->mzeta * eom->nuzeta2) * dadt;

  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_RE_Q,  NC_HIPERT_TWO_FLUIDS_RE_PQ) = (1.0 / eom->mS) * dadt;
  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_IM_Q,  NC_HIPERT_TWO_FLUIDS_IM_PQ) = (1.0 / eom->mS) * dadt;
  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_RE_PQ, NC_HIPERT_TWO_FLUIDS_RE_Q) = (- eom->mS * eom->nuS2) * dadt;
  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_IM_PQ, NC_HIPERT_TWO_FLUIDS_IM_Q) = (- eom->mS * eom->nuS2) * dadt;

  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_RE_ZETA, NC_HIPERT_TWO_FLUIDS_RE_PQ) = eom->Y * dadt;
  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_IM_ZETA, NC_HIPERT_TWO_FLUIDS_IM_PQ) = eom->Y * dadt;

  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_RE_Q, NC_HIPERT_TWO_FLUIDS_RE_PZETA) = eom->Y * dadt;
  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_IM_Q, NC_HIPERT_TWO_FLUIDS_IM_PZETA) = eom->Y * dadt;
  
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
  gdouble tau_i = NC_HIPERT_TWO_FLUIDS_ALPHA2TAU (alphai);

  pert->alpha0 = alphai;

  if (pert->y == NULL)
    pert->y = N_VNew_Serial (NC_HIPERT_TWO_FLUIDS_LEN);

  for (i = 0; i < NC_HIPERT_TWO_FLUIDS_LEN; i++)
  {
    NV_Ith_S (pert->y, i) = vars[i];
  }

  if (!pert->cvode_init)
  {
    flag = CVodeInit (pert->cvode, &_nc_hipert_two_fluids_f, tau_i, pert->y);
    CVODE_CHECK (&flag, "CVodeInit", 1, );
    pert->cvode_init = TRUE;
  }
  else
  {
    flag = CVodeReInit (pert->cvode, tau_i, pert->y);
    CVODE_CHECK (&flag, "CVodeReInit", 1, );
  }

  flag = CVodeSStolerances (pert->cvode, pert->reltol, pert->abstol);
  CVODE_CHECK(&flag, "CVodeSStolerances", 1,);

  flag = CVodeSetMaxNumSteps (pert->cvode, 100000000);
  CVODE_CHECK(&flag, "CVodeSetMaxNumSteps", 1, );

  flag = CVDense (pert->cvode, NC_HIPERT_TWO_FLUIDS_LEN);
  CVODE_CHECK(&flag, "CVDense", 1, );

  flag = CVDlsSetDenseJacFn (pert->cvode, &_nc_hipert_two_fluids_J);
  CVODE_CHECK(&flag, "CVDlsSetDenseJacFn", 1, );  
}

/**
 * nc_hipert_two_fluids_set_init_cond_wkb_zeta:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alphai: the log-redshift time.
 * 
 * Sets the initial conditions for the system evolution using the value of the WKB solution at @alphai with main mode $\zeta$. 
 * 
 */
void 
nc_hipert_two_fluids_set_init_cond_wkb_zeta (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alphai)
{
  gdouble vars[8] = {0.0, };
  gdouble *vars_ptr = vars;

  nc_hipert_two_fluids_wkb_full_zeta (ptf, cosmo, alphai, &vars_ptr);
  nc_hipert_two_fluids_set_init_cond (ptf, cosmo, alphai, vars);
}

/**
 * nc_hipert_two_fluids_set_init_cond_wkb_Q:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alphai: the log-redshift time.
 * 
 * Sets the initial conditions for the system evolution using the value of the WKB solution at @alphai with main mode $Q$. 
 * 
 */
void 
nc_hipert_two_fluids_set_init_cond_wkb_Q (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alphai)
{
  gdouble vars[8] = {0.0, };
  gdouble *vars_ptr = vars;

  nc_hipert_two_fluids_wkb_full_Q (ptf, cosmo, alphai, &vars_ptr);
  nc_hipert_two_fluids_set_init_cond (ptf, cosmo, alphai, vars);
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
  gdouble tau_i = 0.0;
  gint flag;
  gdouble tau_f = NC_HIPERT_TWO_FLUIDS_ALPHA2TAU (alphaf);
  
  arg.cosmo = cosmo;
  arg.ptf = ptf;

  flag = CVodeSetUserData (pert->cvode, &arg);
  CVODE_CHECK (&flag, "CVodeSetFdata", 1, );

  flag = CVode (pert->cvode, tau_f, pert->y, &tau_i, CV_NORMAL);
  CVODE_CHECK (&flag, "CVode[nc_hipert_two_fluids_evolve]", 1, );

  pert->alpha0 = NC_HIPERT_TWO_FLUIDS_TAU2ALPHA (tau_i);
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
