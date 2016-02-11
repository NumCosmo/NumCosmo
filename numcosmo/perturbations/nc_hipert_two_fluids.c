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
  ptf->wkb_zeta = nc_hipert_wkb_new (NC_TYPE_HIPERT_ITWO_FLUIDS, 
                                     (NcHIPertWKBFunc) &nc_hipert_itwo_fluids_nuA2,
                                     (NcHIPertWKBFunc) &nc_hipert_itwo_fluids_VA,
                                     (NcHIPertWKBFunc) &nc_hipert_itwo_fluids_dmzetanuA_nuA,
                                     (NcHIPertWKBEom) &nc_hipert_itwo_fluids_wkb_zeta_eom);

  ptf->wkb_S = nc_hipert_wkb_new (NC_TYPE_HIPERT_ITWO_FLUIDS, 
                                  (NcHIPertWKBFunc) &nc_hipert_itwo_fluids_nuB2,
                                  (NcHIPertWKBFunc) &nc_hipert_itwo_fluids_VB,
                                  (NcHIPertWKBFunc) &nc_hipert_itwo_fluids_dmSnuB_nuB,
                                  (NcHIPertWKBEom) &nc_hipert_itwo_fluids_wkb_S_eom);
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
    NcHIPertTwoFluids *ptf = NC_HIPERT_TWO_FLUIDS (pert);
    nc_hipert_set_mode_k (NC_HIPERT (ptf->wkb_zeta), k);
    nc_hipert_set_mode_k (NC_HIPERT (ptf->wkb_S), k);
  }
}

static void 
_nc_hipert_two_fluids_set_abstol (NcHIPert *pert, gdouble abstol) 
{
  NC_HIPERT_CLASS (nc_hipert_two_fluids_parent_class)->set_abstol (pert, abstol);
  /* Chain up : start */
  if (!pert->prepared)
  {
    NcHIPertTwoFluids *ptf = NC_HIPERT_TWO_FLUIDS (pert);
    nc_hipert_set_abstol (NC_HIPERT (ptf->wkb_zeta), abstol);
    nc_hipert_set_abstol (NC_HIPERT (ptf->wkb_S), abstol);
  }
}

static void 
_nc_hipert_two_fluids_set_reltol (NcHIPert *pert, gdouble reltol) 
{
  NC_HIPERT_CLASS (nc_hipert_two_fluids_parent_class)->set_reltol (pert, reltol);
  /* Chain up : start */
  if (!pert->prepared)
  {
    NcHIPertTwoFluids *ptf = NC_HIPERT_TWO_FLUIDS (pert);
    nc_hipert_set_reltol (NC_HIPERT (ptf->wkb_zeta), reltol);
    nc_hipert_set_reltol (NC_HIPERT (ptf->wkb_S), reltol);
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
                                         "sys-size", /*NC_HIPERT_TWO_FLUIDS_LEN*/ 4,
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
 * nc_hipert_two_fluids_YAB:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha: log-redshift time.
 * 
 * FIXME
 * 
 * Return: FIXME
 */
gdouble 
nc_hipert_two_fluids_YAB (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  const gdouble YAB = eom->Y * sqrt (eom->mzeta * eom->mS);
  return YAB;
}

/**
 * nc_hipert_two_fluids_LA:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha: log-redshift time.
 * 
 * FIXME
 * 
 * Return: FIXME
 */
gdouble 
nc_hipert_two_fluids_LA (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  const gdouble nuA = nc_hipert_wkb_nuA (ptf->wkb_zeta, G_OBJECT (cosmo), alpha);
  const gdouble dmzetanuA_nuA = nc_hipert_itwo_fluids_dmzetanuA_nuA (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  const gdouble LA = 0.5 * dmzetanuA_nuA / (eom->mzeta * nuA);

  return LA;
}

/**
 * nc_hipert_two_fluids_LB:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha: log-redshift time.
 * 
 * FIXME
 * 
 * Return: FIXME
 */
gdouble 
nc_hipert_two_fluids_LB (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  const gdouble nuB = nc_hipert_wkb_nuA (ptf->wkb_S, G_OBJECT (cosmo), alpha);
  const gdouble dmSnuB_nuB    = nc_hipert_itwo_fluids_dmSnuB_nuB (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  const gdouble LB = 0.5 * dmSnuB_nuB / (eom->mS * nuB);
  return LB;
}

/**
 * nc_hipert_two_fluids_eom_full:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha: log-redshift time.
 * @eom_full: (out) (transfer none): Equation of motion variables.
 * 
 * FIXME
 * 
 */
void 
nc_hipert_two_fluids_eom_full (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, NcHIPertITwoFluidsEOM **eom_full)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  *eom_full = nc_hipert_itwo_fluids_eom_full (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  return;
}


/**
 * nc_hipert_two_fluids_nuzeta:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha: log-redshift time.
 * 
 * FIXME
 * 
 * Return: FIXME
 */
gdouble 
nc_hipert_two_fluids_nuzeta2 (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  return eom->nuzeta2;
}

/**
 * nc_hipert_two_fluids_nuS2:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha: log-redshift time.
 * 
 * FIXME
 * 
 * Return: FIXME
 */
gdouble 
nc_hipert_two_fluids_nuS2 (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  return eom->nuS2;
}

/**
 * nc_hipert_two_fluids_laq2:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha: log-redshift time.
 * 
 * FIXME
 * 
 * Return: FIXME
 */
gdouble 
nc_hipert_two_fluids_laq2 (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  const gdouble y = nc_hipert_two_fluids_YAB (ptf, cosmo, alpha);
  const gdouble y2 = y * y;
  const gdouble abs_nuzeta2_m_nuS2 = fabs (eom->nuzeta2 - eom->nuS2);
  const gdouble F = ncm_sqrt1px_m1 (4.0 * eom->nuzeta2 * eom->nuS2 * y2 / (abs_nuzeta2_m_nuS2 * abs_nuzeta2_m_nuS2));
    
  if (eom->nuzeta2 > eom->nuS2)
    return abs_nuzeta2_m_nuS2 * (2.0 + F); 
  else
    return abs_nuzeta2_m_nuS2 * F;
}

/**
 * nc_hipert_two_fluids_lazeta2:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha: log-redshift time.
 * 
 * FIXME
 * 
 * Return: FIXME
 */
gdouble 
nc_hipert_two_fluids_lazeta2 (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  const gdouble y = nc_hipert_two_fluids_YAB (ptf, cosmo, alpha);
  const gdouble y2 = y * y;
  const gdouble abs_nuzeta2_m_nuS2 = fabs (eom->nuzeta2 - eom->nuS2);
  const gdouble F = ncm_sqrt1px_m1 (4.0 * eom->nuzeta2 * eom->nuS2 * y2 / (abs_nuzeta2_m_nuS2 * abs_nuzeta2_m_nuS2));
    
  if (eom->nuzeta2 < eom->nuS2)
    return abs_nuzeta2_m_nuS2 * (2.0 + F); 
  else
    return abs_nuzeta2_m_nuS2 * F;
}

/**
 * nc_hipert_two_fluids_nuplus2:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha: log-redshift time.
 * 
 * FIXME
 * 
 * Return: FIXME
 */
gdouble 
nc_hipert_two_fluids_nuplus2 (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  const gdouble y = nc_hipert_two_fluids_YAB (ptf, cosmo, alpha);
  const gdouble y2 = y * y;
  const gdouble abs_nuzeta2_m_nuS2 = fabs (eom->nuzeta2 - eom->nuS2);
  const gdouble F = ncm_sqrt1px_m1 (4.0 * eom->nuzeta2 * eom->nuS2 * y2 / (abs_nuzeta2_m_nuS2 * abs_nuzeta2_m_nuS2));
    
  if (eom->nuzeta2 > eom->nuS2)
    return eom->nuzeta2 + abs_nuzeta2_m_nuS2 * F * 0.5; 
  else
    return eom->nuS2 + abs_nuzeta2_m_nuS2 * F * 0.5;
}

/**
 * nc_hipert_two_fluids_numinus2:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha: log-redshift time.
 * 
 * FIXME
 * 
 * Return: FIXME
 */
gdouble 
nc_hipert_two_fluids_numinus2 (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  const gdouble y = nc_hipert_two_fluids_YAB (ptf, cosmo, alpha);
  const gdouble y2 = y * y;
  const gdouble abs_nuzeta2_m_nuS2 = fabs (eom->nuzeta2 - eom->nuS2);
  const gdouble F = ncm_sqrt1px_m1 (4.0 * eom->nuzeta2 * eom->nuS2 * y2 / (abs_nuzeta2_m_nuS2 * abs_nuzeta2_m_nuS2));
    
  if (eom->nuzeta2 > eom->nuS2)
    return eom->nuS2 - abs_nuzeta2_m_nuS2 * F * 0.5; 
  else
    return eom->nuzeta2 - abs_nuzeta2_m_nuS2 * F * 0.5;
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
  nc_hipert_wkb_q (ptf->wkb_zeta, G_OBJECT (cosmo), alpha, Re_zeta, Im_zeta);  
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
  nc_hipert_wkb_q_p (ptf->wkb_zeta, G_OBJECT (cosmo), alpha, Re_zeta, Im_zeta, Re_Pzeta, Im_Pzeta);  
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
  nc_hipert_wkb_q (ptf->wkb_S, G_OBJECT (cosmo), alpha, Re_Q, Im_Q);  
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
  nc_hipert_wkb_q_p (ptf->wkb_S, G_OBJECT (cosmo), alpha, Re_Q, Im_Q, Re_PQ, Im_PQ);  
}

/**
 * nc_hipert_two_fluids_wkb_full_zeta:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @alpha: the log-redshift time.
 * @vars: (inout) (array fixed-size=8) (element-type double): Perturbations variables conforming to #NcHIPertTwoFluidsVars.
 * 
 * Computes the approximated solution for the mode $Q$ and its momentum for the mode $k$ 
 * at the time $\alpha$ using the solution for the mode $\zeta$ present in @vars.
 * 
 */
void
nc_hipert_two_fluids_wkb_full_zeta (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, gdouble **vars)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  complex double zeta, Q, PQ,C;
  /*gdouble Re_zeta, Im_zeta, Re_Pzeta, Im_Pzeta;*/
  const gdouble nuA = nc_hipert_wkb_nuA (ptf->wkb_zeta, G_OBJECT (cosmo), alpha);
  const gdouble nuB = nc_hipert_wkb_nuA (ptf->wkb_S, G_OBJECT (cosmo), alpha);
  const gdouble k = pert->k;
  const gdouble dmSnuB_nuB    = nc_hipert_itwo_fluids_dmSnuB_nuB (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  const gdouble dmzetanuA_nuA = nc_hipert_itwo_fluids_dmzetanuA_nuA (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  const gdouble LA = 0.5 * k * dmzetanuA_nuA / (eom->mzeta * nuA);
  const gdouble LB = 0.5 * k * dmSnuB_nuB / (eom->mS * nuB);
  const gdouble YAB = eom->Y * sqrt (nuA * nuB * eom->mzeta * eom->mS) / k;
  const gdouble k2 = k * k;

  /*nc_hipert_wkb_q_p (ptf->wkb_zeta, G_OBJECT (cosmo), alpha, &Re_zeta, &Im_zeta, &Re_Pzeta, &Im_Pzeta);*/
  
  zeta  = vars[0][NC_HIPERT_TWO_FLUIDS_RE_ZETA]  + I * vars[0][NC_HIPERT_TWO_FLUIDS_IM_ZETA];
  /*Pzeta = vars[0][NC_HIPERT_TWO_FLUIDS_RE_PZETA] + I * vars[0][NC_HIPERT_TWO_FLUIDS_IM_PZETA];*/
  
  C = YAB * (k - I * LA) * zeta * sqrt (2.0 * eom->mzeta * nuA) / (sqrt (2.0 * k) * k * (nuA * nuA - nuB * nuB));

  /*  printf ("% 20.15g % 20.15g % 20.15e % 20.15e % 20.15e\n", alpha, nc_hicosmo_x_alpha (cosmo, alpha), fabs (alpha / nuA), fabs (alpha / nuB), fabs (YAB));*/
  
  Q = sqrt (k / (eom->mS * nuB)) * (k * nuA + I * LB * nuB) * C;
  PQ = - I * sqrt ((eom->mS * nuB) / k) * (k2 + LB * LB) * nuB * C;
/*
  vars[0][NC_HIPERT_TWO_FLUIDS_RE_ZETA] = creal (zeta);
  vars[0][NC_HIPERT_TWO_FLUIDS_IM_ZETA] = cimag (zeta);

  vars[0][NC_HIPERT_TWO_FLUIDS_RE_PZETA] = creal (Pzeta);
  vars[0][NC_HIPERT_TWO_FLUIDS_IM_PZETA] = cimag (Pzeta);
*/
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
 * Computes the approximated solution for the mode $\zeta$ and its momentum for the mode $k$ 
 * at the time $\alpha$ using the solution for the mode $Q$ present in @vars. 
 * 
 */
void
nc_hipert_two_fluids_wkb_full_Q (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, gdouble **vars)
{
  NcHIPert *pert = NC_HIPERT (ptf);
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  complex double zeta, Pzeta, Q, C;
  /*gdouble Re_Q, Im_Q, Re_PQ, Im_PQ;*/
  const gdouble nuA = nc_hipert_wkb_nuA (ptf->wkb_zeta, G_OBJECT (cosmo), alpha);
  const gdouble nuB = nc_hipert_wkb_nuA (ptf->wkb_S, G_OBJECT (cosmo), alpha);
  const gdouble k    = pert->k;
  const gdouble dmSnuB_nuB    = nc_hipert_itwo_fluids_dmSnuB_nuB (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  const gdouble dmzetanuA_nuA = nc_hipert_itwo_fluids_dmzetanuA_nuA (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, pert->k);
  const gdouble LA   = 0.5 * k * dmzetanuA_nuA / (eom->mzeta * nuA);
  const gdouble LB   = 0.5 * k * dmSnuB_nuB / (eom->mS * nuB);
  const gdouble YAB  = eom->Y * sqrt (nuA * nuB * eom->mzeta * eom->mS) / k;
  const gdouble k2   = k * k;
  
  /*nc_hipert_wkb_q_p (ptf->wkb_S, G_OBJECT (cosmo), alpha, &Re_Q, &Im_Q, &Re_PQ, &Im_PQ);*/
  
  Q  = vars[0][NC_HIPERT_TWO_FLUIDS_RE_Q]  + I * vars[0][NC_HIPERT_TWO_FLUIDS_IM_Q];
  /*PQ = vars[0][NC_HIPERT_TWO_FLUIDS_RE_PQ] + I * vars[0][NC_HIPERT_TWO_FLUIDS_IM_PQ];*/

  C = YAB * (k - I * LB) * sqrt (2.0 * eom->mS * nuB) * Q / (sqrt (2.0 * k) * k * (nuB * nuB - nuA * nuA));

  zeta  = sqrt (k / (eom->mzeta * nuA)) * (k * nuB + I * LA * nuA) * C;
  Pzeta = - I * sqrt ((eom->mzeta * nuA) / k) * (k2 + LA * LA) * nuA * C;

  /*  printf ("% 20.15g % 20.15g % 20.15e % 20.15e % 20.15e\n", alpha, nc_hicosmo_x_alpha (cosmo, alpha), fabs (alpha / nuA), fabs (alpha / nuB), fabs (YAB));  */
  /*  printf ("# WKB % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", alpha, eom->mzeta, nuA2, eom->mS, nuB2, eom->Y);  */

  vars[0][NC_HIPERT_TWO_FLUIDS_RE_ZETA]  = creal (zeta);
  vars[0][NC_HIPERT_TWO_FLUIDS_IM_ZETA]  = cimag (zeta);

  vars[0][NC_HIPERT_TWO_FLUIDS_RE_PZETA] = creal (Pzeta);
  vars[0][NC_HIPERT_TWO_FLUIDS_IM_PZETA] = cimag (Pzeta);

/*  
  vars[0][NC_HIPERT_TWO_FLUIDS_RE_Q]     = creal (Q);
  vars[0][NC_HIPERT_TWO_FLUIDS_IM_Q]     = cimag (Q);

  vars[0][NC_HIPERT_TWO_FLUIDS_RE_PQ]    = creal (PQ);
  vars[0][NC_HIPERT_TWO_FLUIDS_IM_PQ]    = cimag (PQ);  
*/
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
  return nc_hipert_wkb_maxtime (ptf->wkb_zeta, G_OBJECT (cosmo), alpha0, alpha1);  
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
  return nc_hipert_wkb_maxtime (ptf->wkb_S, G_OBJECT (cosmo), alpha0, alpha1);  
}

/**
 * nc_hipert_two_fluids_wkb_zeta_maxtime_prec:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @cmp: Comparison type.
 * @prec: Required precision.
 * @alpha0: the initial log-redshift time.
 * @alpha1: the final log-redshift time.
 * 
 * Search for the instant at which the WKB approximation starts to fails within the asked precision. 
 * 
 * Returns: the instant $\alpha$ between $\alpha_0$ and $\alpha_1$ or NaN if not found.
 */
gdouble 
nc_hipert_two_fluids_wkb_zeta_maxtime_prec (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, NcHIPertWKBCmp cmp, gdouble prec, gdouble alpha0, gdouble alpha1)
{
  return nc_hipert_wkb_maxtime_prec (ptf->wkb_zeta, G_OBJECT (cosmo), cmp, prec, alpha0, alpha1);  
}

/**
 * nc_hipert_two_fluids_wkb_S_maxtime_prec:
 * @ptf: a #NcHIPertTwoFluids.
 * @cosmo: a #NcHICosmo.
 * @cmp: Comparison type.
 * @prec: Required precision.
 * @alpha0: the initial log-redshift time.
 * @alpha1: the final log-redshift time.
 * 
 * Search for the instant at which the WKB approximation starts to fails within the asked precision. 
 * 
 * Returns: the instant $\alpha$ between $\alpha_0$ and $\alpha_1$ or NaN if not found.
 */
gdouble 
nc_hipert_two_fluids_wkb_S_maxtime_prec (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, NcHIPertWKBCmp cmp, gdouble prec, gdouble alpha0, gdouble alpha1)
{
  return nc_hipert_wkb_maxtime_prec (ptf->wkb_S, G_OBJECT (cosmo), cmp, prec, alpha0, alpha1);
}

static gint
_nc_hipert_two_fluids_f (realtype alpha, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHIPertTwoFluidsArg *arg = (NcHIPertTwoFluidsArg *) f_data;
  const gdouble k = NC_HIPERT (arg->ptf)->k;
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_full (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);
/*  
  NV_Ith_S (ydot, NC_HIPERT_TWO_FLUIDS_RE_ZETA)  = NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_RE_PZETA) / eom->mzeta + eom->Y * NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_RE_PQ);
  NV_Ith_S (ydot, NC_HIPERT_TWO_FLUIDS_IM_ZETA)  = NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_IM_PZETA) / eom->mzeta + eom->Y * NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_IM_PQ);

  NV_Ith_S (ydot, NC_HIPERT_TWO_FLUIDS_RE_Q)     = NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_RE_PQ) / eom->mS + eom->Y * NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_RE_PZETA);
  NV_Ith_S (ydot, NC_HIPERT_TWO_FLUIDS_IM_Q)     = NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_IM_PQ) / eom->mS + eom->Y * NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_IM_PZETA);

  NV_Ith_S (ydot, NC_HIPERT_TWO_FLUIDS_RE_PZETA) = - eom->mzeta * eom->nuzeta2 * NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_RE_ZETA);
  NV_Ith_S (ydot, NC_HIPERT_TWO_FLUIDS_IM_PZETA) = - eom->mzeta * eom->nuzeta2 * NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_IM_ZETA);

  NV_Ith_S (ydot, NC_HIPERT_TWO_FLUIDS_RE_PQ)    = - eom->mS * eom->nuS2 * NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_RE_Q);
  NV_Ith_S (ydot, NC_HIPERT_TWO_FLUIDS_IM_PQ)    = - eom->mS * eom->nuS2 * NV_Ith_S (y, NC_HIPERT_TWO_FLUIDS_IM_Q);
*/
  const gdouble Ip   = NV_Ith_S (y, 0);
  const gdouble Im   = NV_Ith_S (y, 1);
  const gdouble phip = NV_Ith_S (y, 2);
  const gdouble phim = NV_Ith_S (y, 3);
  
  NV_Ith_S (ydot, 0) = - eom->US[0] * Ip * sin (2.0 * phip) / eom->nu_plus + 2.0 * sqrt (Im * Ip) * 
    (-eom->US[2] * cos (phip) * sin (phim) + eom->UA * (eom->nu_plus * sin (phip) * sin (phim) + eom->nu_minus * cos (phip) * cos (phim))) / k;

  NV_Ith_S (ydot, 1) = - eom->US[1] * Im * sin (2.0 * phim) / eom->nu_minus + 2.0 * sqrt (Im * Ip) * 
    (-eom->US[2] * cos (phim) * sin (phip) - eom->UA * (eom->nu_plus * cos (phip) * cos (phim) + eom->nu_minus * sin (phip) * sin (phim))) / k;

  NV_Ith_S (ydot, 2) = (eom->nu_plus - 0.5 * eom->US[0] / eom->nu_plus) - 0.5 * eom->US[0] * cos (2.0 * phip) / eom->nu_plus 
    + sqrt (Im / Ip) * (eom->US[2] * sin (phip) * sin (phim) + eom->UA * (eom->nu_plus * cos (phip) * sin (phim) - eom->nu_minus * sin (phip) * cos (phim))) / k;

  NV_Ith_S (ydot, 3) = (eom->nu_minus - 0.5 * eom->US[1] / eom->nu_minus) - 0.5 * eom->US[1] * cos (2.0 * phim) / eom->nu_minus 
    + sqrt (Ip / Im) * (eom->US[2] * sin (phip) * sin (phim) + eom->UA * (eom->nu_plus * cos (phip) * sin (phim) - eom->nu_minus * sin (phip) * cos (phim))) / k;

/*  
  NV_Ith_S (ydot, 0) = - eom->Uplus  * Ip * cos (2.0 * phip) + Im * (eom->Wplus * sin (phip) * sin (phim) - eom->Wminus * cos (phip) * cos (phim));
  NV_Ith_S (ydot, 1) = - eom->Uminus * Im * cos (2.0 * phim) - Ip * (eom->Wplus * cos (phip) * cos (phim) - eom->Wminus * sin (phip) * sin (phim));
  NV_Ith_S (ydot, 2) = eom->nu_plus  + eom->Uplus  * sin (2.0 * phip) + Im / Ip * (eom->Wplus * cos (phip) * sin (phim) + eom->Wminus * cos (phim) * sin (phip));
  NV_Ith_S (ydot, 3) = eom->nu_minus + eom->Uminus * sin (2.0 * phim) + Ip / Im * (eom->Wplus * cos (phip) * sin (phim) + eom->Wminus * cos (phim) * sin (phip));
*/
/*
printf ("#-------------------------------------------------------------------------------------------------------\n");
printf ("% 20.15g % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e\n", alpha, eom->nu_plus, eom->nu_minus, eom->Uplus, eom->Uminus, eom->Wplus, eom->Wminus);
printf ("% 20.15e | % 20.15e % 20.15e % 20.15e = % 20.15e\n", Ip, - eom->Uplus  * Ip * cos (2.0 * phip), Im * eom->Wplus * sin (phip) * sin (phim), - Im * eom->Wminus * cos (phip) * cos (phim), NV_Ith_S (ydot, 0));
printf ("% 20.15e | % 20.15e % 20.15e % 20.15e = % 20.15e\n", Im, - eom->Uminus * Im * cos (2.0 * phim), - Ip * eom->Wplus * cos (phip) * cos (phim), + Ip * eom->Wminus * sin (phip) * sin (phim), NV_Ith_S (ydot, 1));
printf ("% 20.15e | % 20.15e % 20.15e % 20.15e % 20.15e = % 20.15e\n", phip, eom->nu_plus, eom->Uplus  * sin (2.0 * phip), + Im / Ip * eom->Wplus * cos (phip) * sin (phim), + Im / Ip * eom->Wminus * cos (phim) * sin (phip), NV_Ith_S (ydot, 2));
printf ("% 20.15e | % 20.15e % 20.15e % 20.15e % 20.15e = % 20.15e\n", phim, eom->nu_minus, eom->Uminus * sin (2.0 * phim), + Ip / Im * eom->Wplus * cos (phip) * sin (phim), + Ip / Im *  eom->Wminus * cos (phim) * sin (phip), NV_Ith_S (ydot, 3));
*/

printf ("#-----------------------------------------------------------------------------\n");
printf ("% 20.15g % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e\n", alpha, eom->nu_plus, eom->nu_minus, eom->US[0] / gsl_pow_2 (eom->nu_plus), eom->US[1] / gsl_pow_2 (eom->nu_minus), eom->US[2], eom->UA);
printf ("% 20.15g % 20.15e % 20.15e % 20.15e % 20.15e\n", alpha, Ip, Im, phip, phim);
printf ("% 20.15g % 20.15e % 20.15e % 20.15e % 20.15e\n", alpha, NV_Ith_S (ydot, 0), NV_Ith_S (ydot, 1), NV_Ith_S (ydot, 2), NV_Ith_S (ydot, 3));


  return 0;
}

static gint
_nc_hipert_two_fluids_J (_NCM_SUNDIALS_INT_TYPE N, realtype alpha, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHIPertTwoFluidsArg *arg = (NcHIPertTwoFluidsArg *) jac_data;
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_full (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, NC_HIPERT (arg->ptf)->k);
/*
  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_RE_ZETA,  NC_HIPERT_TWO_FLUIDS_RE_PZETA) = 1.0 / eom->mzeta;
  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_IM_ZETA,  NC_HIPERT_TWO_FLUIDS_IM_PZETA) = 1.0 / eom->mzeta;
  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_RE_PZETA, NC_HIPERT_TWO_FLUIDS_RE_ZETA) = - eom->mzeta * eom->nuzeta2;
  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_IM_PZETA, NC_HIPERT_TWO_FLUIDS_IM_ZETA) = - eom->mzeta * eom->nuzeta2;

  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_RE_Q,  NC_HIPERT_TWO_FLUIDS_RE_PQ) = 1.0 / eom->mS;
  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_IM_Q,  NC_HIPERT_TWO_FLUIDS_IM_PQ) = 1.0 / eom->mS;
  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_RE_PQ, NC_HIPERT_TWO_FLUIDS_RE_Q) = - eom->mS * eom->nuS2;
  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_IM_PQ, NC_HIPERT_TWO_FLUIDS_IM_Q) = - eom->mS * eom->nuS2;

  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_RE_ZETA, NC_HIPERT_TWO_FLUIDS_RE_PQ) = eom->Y;
  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_IM_ZETA, NC_HIPERT_TWO_FLUIDS_IM_PQ) = eom->Y;

  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_RE_Q, NC_HIPERT_TWO_FLUIDS_RE_PZETA) = eom->Y;
  DENSE_ELEM (J, NC_HIPERT_TWO_FLUIDS_IM_Q, NC_HIPERT_TWO_FLUIDS_IM_PZETA) = eom->Y;
  */
  const gdouble Ip   = NV_Ith_S (y, 0);
  const gdouble Im   = NV_Ith_S (y, 1);
  const gdouble phip = NV_Ith_S (y, 2);
  const gdouble phim = NV_Ith_S (y, 3);

  /*- eom->Uplus * Ip * cos (2.0 * phip) + Im * (eom->Wplus * sin (phip) * sin (phim) - eom->Wminus * cos (phip) * cos (phim))*/
  DENSE_ELEM (J, 0, 0) = - eom->Uplus * cos (2.0 * phip);
  DENSE_ELEM (J, 0, 1) = + (eom->Wplus * sin (phip) * sin (phim) - eom->Wminus * cos (phip) * cos (phim));
  DENSE_ELEM (J, 0, 2) = + 2.0 * eom->Uplus * Ip * sin (2.0 * phip) + Im * (eom->Wplus * cos (phip) * sin (phim) + eom->Wminus * sin (phip) * cos (phim));
  DENSE_ELEM (J, 0, 3) =                                            + Im * (eom->Wplus * sin (phip) * cos (phim) + eom->Wminus * cos (phip) * sin (phim));

  /*- eom->Uminus * Im * cos (2.0 * phim) - Ip * (eom->Wplus * cos (phip) * cos (phim) - eom->Wminus * sin (phip) * sin (phim))*/
  DENSE_ELEM (J, 1, 0) = - (eom->Wplus * cos (phip) * cos (phim) - eom->Wminus * sin (phip) * sin (phim));
  DENSE_ELEM (J, 1, 1) = - eom->Uminus * cos (2.0 * phim);
  DENSE_ELEM (J, 1, 2) =                                             + Ip * (eom->Wplus * sin (phip) * cos (phim) + eom->Wminus * cos (phip) * sin (phim));
  DENSE_ELEM (J, 1, 3) = + 2.0 * eom->Uminus * Im * sin (2.0 * phim) + Ip * (eom->Wplus * cos (phip) * sin (phim) + eom->Wminus * sin (phip) * cos (phim));

  /*eom->nu_plus + eom->Uplus  * sin (2.0 * phip) + Im / Ip * (eom->Wplus * cos (phip) * sin (phim) + eom->Wminus * cos (phim) * sin (phip))*/
  DENSE_ELEM (J, 2, 0) = - Im / (Ip * Ip) * (eom->Wplus * cos (phip) * sin (phim) + eom->Wminus * cos (phim) * sin (phip));
  DENSE_ELEM (J, 2, 1) = + 1.0 / Ip       * (eom->Wplus * cos (phip) * sin (phim) + eom->Wminus * cos (phim) * sin (phip));
  DENSE_ELEM (J, 2, 2) = + 2.0 * eom->Uplus  * cos (2.0 * phip) - Im / Ip * (eom->Wplus * sin (phip) * sin (phim) - eom->Wminus * cos (phim) * cos (phip));
  DENSE_ELEM (J, 2, 3) =                                        + Im / Ip * (eom->Wplus * cos (phip) * cos (phim) - eom->Wminus * sin (phim) * sin (phip));

  /*eom->nu_minus + eom->Uminus * sin (2.0 * phim) + Ip / Im * (eom->Wplus * cos (phip) * sin (phim) + eom->Wminus * cos (phim) * sin (phip))*/
  DENSE_ELEM (J, 3, 0) = + 1.0 / Im       * (eom->Wplus * cos (phip) * sin (phim) + eom->Wminus * cos (phim) * sin (phip));
  DENSE_ELEM (J, 3, 1) = - Ip / (Im * Im) * (eom->Wplus * cos (phip) * sin (phim) + eom->Wminus * cos (phim) * sin (phip));
  DENSE_ELEM (J, 3, 2) =                                        - Ip / Im * (eom->Wplus * sin (phip) * sin (phim) - eom->Wminus * cos (phim) * cos (phip));
  DENSE_ELEM (J, 3, 3) = + 2.0 * eom->Uminus * cos (2.0 * phim) + Ip / Im * (eom->Wplus * cos (phip) * cos (phim) - eom->Wminus * sin (phim) * sin (phip));

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

  if (ptf->abstol == NULL)
  {
    ptf->abstol = N_VNew_Serial (/*NC_HIPERT_TWO_FLUIDS_LEN*/ 4);
    NV_Ith_S (ptf->abstol, 0) = 0.0;
    NV_Ith_S (ptf->abstol, 1) = 0.0;
    NV_Ith_S (ptf->abstol, 2) = pert->reltol;
    NV_Ith_S (ptf->abstol, 3) = pert->reltol;
  }

  /*printf ("# Initial conditions\n% 20.15g ", alphai);*/  
  for (i = 0; i < /*NC_HIPERT_TWO_FLUIDS_LEN*/ 4; i++)
  {
    NV_Ith_S (pert->y, i) = vars[i];
    /*printf ("% 20.15g ", vars[i]);*/
  }
  /*printf ("\n");*/

  if (!pert->cvode_init)
  {
    flag = CVodeInit (pert->cvode, &_nc_hipert_two_fluids_f, alphai, pert->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
    pert->cvode_init = TRUE;
  }
  else
  {
    flag = CVodeReInit (pert->cvode, alphai, pert->y);
    NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );
  }

  flag = CVodeSVtolerances (pert->cvode, pert->reltol, ptf->abstol);
  NCM_CVODE_CHECK (&flag, "CVodeSVtolerances", 1,);

  flag = CVodeSetMaxNumSteps (pert->cvode, 100000000);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

  flag = CVDense (pert->cvode, /*NC_HIPERT_TWO_FLUIDS_LEN*/ 4);
  NCM_CVODE_CHECK (&flag, "CVDense", 1, );

  flag = CVDlsSetDenseJacFn (pert->cvode, &_nc_hipert_two_fluids_J);
  NCM_CVODE_CHECK (&flag, "CVDlsSetDenseJacFn", 1, );  
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

  nc_hipert_two_fluids_wkb_zeta_Pzeta (ptf, cosmo, alphai,
                                       &vars[NC_HIPERT_TWO_FLUIDS_RE_ZETA],  &vars[NC_HIPERT_TWO_FLUIDS_IM_ZETA],
                                       &vars[NC_HIPERT_TWO_FLUIDS_RE_PZETA], &vars[NC_HIPERT_TWO_FLUIDS_IM_PZETA]
                                       );
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

  nc_hipert_two_fluids_wkb_Q_PQ (ptf, cosmo, alphai,
                                 &vars[NC_HIPERT_TWO_FLUIDS_RE_Q],  &vars[NC_HIPERT_TWO_FLUIDS_IM_Q],
                                 &vars[NC_HIPERT_TWO_FLUIDS_RE_PQ], &vars[NC_HIPERT_TWO_FLUIDS_IM_PQ]
                                 );
  
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
  gdouble alpha_i = 0.0;
  gint flag;
  
  arg.cosmo = cosmo;
  arg.ptf = ptf;

  flag = CVodeSetUserData (pert->cvode, &arg);
  NCM_CVODE_CHECK (&flag, "CVodeSetFdata", 1, );

  while (TRUE)
  {
    flag = CVode (pert->cvode, alphaf, pert->y, &alpha_i, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "CVode[nc_hipert_two_fluids_evolve]", 1, );

    if (NV_Ith_S (pert->y, 2) > 1.0e3 || NV_Ith_S (pert->y, 3) > 1.0e3)
    {
      NV_Ith_S (pert->y, 2) = fmod (NV_Ith_S (pert->y, 2), 2.0 * M_PI);
      NV_Ith_S (pert->y, 3) = fmod (NV_Ith_S (pert->y, 3), 2.0 * M_PI);

      flag = CVodeReInit (pert->cvode, alpha_i, pert->y);
      NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );

      flag = CVodeSVtolerances (pert->cvode, pert->reltol, ptf->abstol);
      NCM_CVODE_CHECK (&flag, "CVodeSVtolerances", 1,);

      flag = CVodeSetMaxNumSteps (pert->cvode, 100000000);
      NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

      flag = CVDense (pert->cvode, /*NC_HIPERT_TWO_FLUIDS_LEN*/ 4);
      NCM_CVODE_CHECK (&flag, "CVDense", 1, );

      flag = CVDlsSetDenseJacFn (pert->cvode, &_nc_hipert_two_fluids_J);
      NCM_CVODE_CHECK (&flag, "CVDlsSetDenseJacFn", 1, );
    }

    if (alpha_i >= alphaf)
      break;
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
