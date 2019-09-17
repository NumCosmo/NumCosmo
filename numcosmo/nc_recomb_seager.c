/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_recomb_seager.c
 *
 *  Mon November 05 18:28:23 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
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
 * SECTION:nc_recomb_seager
 * @title: NcRecombSeager
 * @short_description: Cosmic recombination implementing Seager (1999).
 * @include: numcosmo/nc_recomb_seager.h
 *
 * Cosmic recobination as initally describe in [Seager (1999)][XSeager1999]
 * and [Seager (2000)][XSeager2000]. The code includes now all modifications
 * as in [recfast 1.5.2](http://www.astro.ubc.ca/people/scott/recfast.html),
 * which includes the modifications discussed in [Wong (2008)][XWong2008].
 * Nonetheless, we do not include the modification for the matter temperature
 * evolution as describe in [Scott (2009)][XScott2009]. Since we use a more
 * robust integration method such modification for the temperature evolution
 * is simply unnecessary.
 *
 * See [NcRecomb][NcRecomb.description] for symbol definitions.
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
 * This code solves the system of equations for the singly ionized hydrogen $X_\HyII$
 * and helium $X_\HeII$ as well as for the baryon temperature $T_m$.
 *
 * The equations are:
 * \begin{align}
 * \frac{\mathrm{d}X_\HyII}{\mathrm{d}x} &= \frac{X_\HyII X_\e n_\Hy - X_\HyI B_{\HyI, 1s\,{}^2\\!S_{1/2}}(T_m)}{H x}\left[\alpha_\Hy(T_m)\frac{n_\Hy K_{\HyI} X_\HyI \Lambda_\Hy + 1}{n_\Hy K_{\HyI} X_\HyI \left[\Lambda_\Hy + B_{\HyI, 2s\,{}^2\\!S_{1/2}}(T_m) \alpha_\Hy(T_m)\right] + 1}\right], \\\\
 * \frac{\mathrm{d}T_m}{\mathrm{d}x} &= \frac{c}{Hx}\frac{8\sigma_\mathrm{T}a_\mathrm{R} T^4_r}{3m_\mathrm{e}c^2} \frac{X_\e(T_m - T_r)}{1 + X_\He + X_\e} + \frac{2T_m}{x}, \\\\
 * \frac{\mathrm{d}X_\HeII}{\mathrm{d}x} &= \frac{X_\HeII X_\e n_\Hy - X_\HeI B_{\HeI, 1s\,{}^1\\!S_{0}}(T_m)}{H x}\Bigg\\{\left[\alpha_\He(T_m)\frac{n_\Hy K_{\HeI} X_\HeI \Lambda_\He + B^{\HeI, 2s\,{}^1\\!S_{0}}_{\HeI, 2p\,{}^1\\!P_{1}}(T_m)}{n_H K_{\HeI} X_\HeI \left[\Lambda_\He + B_{\HeI, 2s\,{}^1\\!S_{0}}(T_m) \alpha_\He(T_m)\right] + B^{\HeI, 2s\,{}^1\\!S_{0}}_{\HeI, 2p\,{}^1\\!P_{1}}(T_m)}\right] \nonumber\\\\
 * &+ \alpha_\He^\mathrm{t}(T_m)\frac{1}{n_H K_{\HeI}^\mathrm{t} X_\HeI B_{\HeI, 2p\,{}^3\\!P_\mathrm{mean}}(T_m) \alpha_\He^\mathrm{t}(T_m) + 1} \Bigg\\},
 * \end{align}
 *
 * The Boltzmann factor for hydrogen levels are given by
 * $$B_{\HyI, l}(T_m) = k_\mathrm{e}^3(T_m)\,\exp\left[-E_{\HyI, l} / (k_\mathrm{B}T_m)\right],$$
 * for $l = 1s\,{}^2\\!S_{1/2}, 2s\,{}^2\\!S_{1/2}$, see ncm_c_boltzmann_factor_HI_1s_2S0_5(),
 * ncm_c_boltzmann_factor_HI_2s_2S0_5() and ncm_c_thermal_wn_e() for the definition
 * of the electron thermal wavenumber $k_\mathrm{e}$.
 *
 * For the helium-I levels the Boltzmann factors are
 * $$B_{\HeI, l}(T_m) = 4 k_\mathrm{e}^3(T_m)\,\exp\left[-E_{\HeI, l} / (k_\mathrm{B}T_m)\right],$$
 * where the levels $l$ used are $1s\,{}^1\\!S_{0}$, $2s\,{}^1\\!S_{0}$ and $2p\,{}^1\\!P_{1}$,
 * see ncm_c_HeI_ion_wn_1s_1S0(), ncm_c_HeI_ion_wn_2s_1S0(),
 * ncm_c_HeI_ion_wn_2p_1P1(). The symbol $B^{\HeI, 2s\,{}^1\\!S_{0}}_{\HeI, 2p\,{}^1\\!P_{1}}(T_m)$
 * represents the ratio of two Boltzmann factors, i.e.,
 * $$B^{\HeI, 2s\,{}^1\\!S_{0}}_{\HeI, 2p\,{}^1\\!P_{1}}(T_m) =
 * \exp\left[-(E_{\HeI, 2s\,{}^1\\!S_{0}} - E_{\HeI, 2p\,{}^1\\!P_{1}})/ (k_\mathrm{B}T_m)\right].$$
 *
 * The two photon decaying rates for $\Hy$ is $\Lambda_\Hy$ and is given by ncm_c_decay_H_rate_2s_1s(),
 * and for helium-I is $\Lambda_\He$ given by ncm_c_decay_He_rate_2s_1s().
 *
 * The Case B coefficient for $\Hy$, $\alpha_\Hy$ is calculated by nc_recomb_seager_pequignot_HI_case_B(),
 * while the Case B coefficient for $\He$, $\alpha_\He$ and $\alpha_\He^\mathrm{t}$ are given
 * respectively by nc_recomb_seager_hummer_HeI_case_B() and nc_recomb_seager_hummer_HeI_case_B_trip().
 *
 * The flags in #NcRecombSeagerOpt define which are the $K$ factors to be used and whether to include
 * triplets factors in the $\He$ rate.
 *
 * The code uses nc_recomb_HeII_ion_saha_x_by_HeIII_He() to obtain the value of $\lambda$
 * where the numerical integration will start. The integration is performed considering all
 * components without any switching or approximation.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_spline_func.h"
#include "math/ncm_util.h"
#include "nc_enum_types.h"
#include "nc_recomb.h"
#include "nc_recomb_seager.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_hyperg.h>

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#define SUN_DENSE_ACCESS SM_ELEMENT_D
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcRecombSeagerPrivate
{
  NcRecombSeagerOpt opts;
  NcRecombSeagerKHI2p2Pmean K_HI_2p_2Pmean;
  NcRecombSeagerKHeI2p KX_HeI_2p_1P1;
  NcRecombSeagerKHeI2p KX_HeI_2p_3Pmean;
  NcRecombSeagerKHeI2pGrad KX_HeI_2p_1P1_grad;
  NcRecombSeagerKHeI2pGrad KX_HeI_2p_3Pmean_grad;
  gdouble H_fudge;
  gdouble AGauss1, AGauss2;
  gdouble zGauss1, zGauss2;
  gdouble wGauss1, wGauss2;
  gdouble A2P_s;
  gdouble A2P_t;
  gdouble sigma_He_2P_s;
  gdouble sigma_He_2P_t;
  gdouble Pb, Pb_t;
  gdouble Qb, Qb_t;
  gpointer cvode;
  gboolean init;
  N_Vector y0;
  N_Vector y;
  N_Vector abstol;
  SUNMatrix A;
  SUNLinearSolver LS;
  guint n;
  NcmSpline *Xe_s;
  NcmSpline *Xe_reion_s;
  NcmSpline *Xe_recomb_s;
  NcmSpline *XHII_s;
  NcmSpline *XHeII_s;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcRecombSeager, nc_recomb_seager, NC_TYPE_RECOMB);

enum
{
  PROP_0,
  PROP_OPTS,
  PROP_SIZE,
};

static gint H_ion_full_f (realtype lambda, N_Vector y, N_Vector ydot, gpointer f_data);

static gint H_ion_full_J (realtype lambda, N_Vector y, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static gdouble _nc_recomb_seager_K_HI_2p_2Pmean (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble H);
static gdouble _nc_recomb_seager_K_HI_2p_2Pmean_gcor (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble H);

static gdouble _nc_recomb_seager_KX_HeI_2p_1P1 (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble XHI, const gdouble Tm, const gdouble XHeI, const gdouble H, const gdouble n_H);
static gdouble _nc_recomb_seager_KX_HeI_2p_1P1_sobolev (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble XHI, const gdouble Tm, const gdouble XHeI, const gdouble H, const gdouble n_H);
static gdouble _nc_recomb_seager_KX_HeI_2p_1P1_sobolev_cont (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble XHI, const gdouble Tm, const gdouble XHeI, const gdouble H, const gdouble n_H);

static void _nc_recomb_seager_KX_HeI_2p_1P1_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble XHI, const gdouble Tm, const gdouble XHeI, const gdouble H, const gdouble n_H, gdouble grad[3]);
static void _nc_recomb_seager_KX_HeI_2p_1P1_sobolev_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble XHI, const gdouble Tm, const gdouble XHeI, const gdouble H, const gdouble n_H, gdouble grad[3]);
static void _nc_recomb_seager_KX_HeI_2p_1P1_sobolev_cont_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble XHI, const gdouble Tm, const gdouble XHeI, const gdouble H, const gdouble n_H, gdouble grad[3]);

static gdouble _nc_recomb_seager_KX_HeI_2p_3Pmean_sobolev (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble XHI, const gdouble Tm, const gdouble XHeI, const gdouble H, const gdouble n_H);
static gdouble _nc_recomb_seager_KX_HeI_2p_3Pmean_sobolev_cont (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble XHI, const gdouble Tm, const gdouble XHeI, const gdouble H, const gdouble n_H);

static void _nc_recomb_seager_KX_HeI_2p_3Pmean_sobolev_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble XHI, const gdouble Tm, const gdouble XHeI, const gdouble H, const gdouble n_H, gdouble grad[3]);
static void _nc_recomb_seager_KX_HeI_2p_3Pmean_sobolev_cont_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble XHI, const gdouble Tm, const gdouble XHeI, const gdouble H, const gdouble n_H, gdouble grad[3]);

static void
nc_recomb_seager_init (NcRecombSeager *recomb_seager)
{
  NcRecombSeagerPrivate * const self = recomb_seager->priv = nc_recomb_seager_get_instance_private (recomb_seager);

  self->cvode = CVodeCreate (CV_BDF);
  NCM_CVODE_CHECK ((void*)self->cvode, "CVodeCreate", 0, );

  self->init                  = FALSE;
  self->opts                  = 0;
  self->H_fudge               = 0.0;
  self->AGauss1               = -0.140;
  self->AGauss2               = 0.079;
  self->zGauss1               = 7.280;
  self->zGauss2               = 6.730;
  self->wGauss1               = 0.180;
  self->wGauss2               = 0.330;

  self->A2P_s                 = 1.798287e9; /* [s^-1] from recfast.for 1.5.2 (Morton, Wu & Drake (2006)) */
  self->A2P_t                 = 177.58; /* [s^-1] from recfast.for 1.5.2 (Lach & Pachuski (2001)) */
  self->sigma_He_2P_s         = 1.436289e-22; /* [m^2]  from recfast.for 1.5.2 (Hummer & Storey (1998)) */
  self->sigma_He_2P_t         = 1.484872e-22; /* [m^2]  from recfast.for 1.5.2 (Hummer & Storey (1998)) */
  self->Pb                    = 0.36; /* value from KIV (2007) */
  self->Qb                    = 0.86; /* b_He fudge factor     */
  self->Pb_t                  = 0.66; /* value from KIV (2007) */
  self->Qb_t                  = 0.90; /* value from KIV (2007) */

  self->K_HI_2p_2Pmean        = &_nc_recomb_seager_K_HI_2p_2Pmean;
  self->KX_HeI_2p_1P1         = &_nc_recomb_seager_KX_HeI_2p_1P1;

  self->KX_HeI_2p_1P1_grad    = &_nc_recomb_seager_KX_HeI_2p_1P1_grad;

  self->KX_HeI_2p_3Pmean      = NULL;
  self->KX_HeI_2p_3Pmean_grad = NULL;

  self->n                     = 3;

  self->y0                    = N_VNew_Serial (self->n);
  self->y                     = N_VNew_Serial (self->n);
  self->abstol                = N_VNew_Serial (self->n);

  self->A                     = SUNDenseMatrix (self->n, self->n);
  self->LS                    = SUNDenseLinearSolver (self->y, self->A);
  
  NCM_CVODE_CHECK ((gpointer)self->A, "SUNDenseMatrix", 0, );
  NCM_CVODE_CHECK ((gpointer)self->LS, "SUNDenseLinearSolver", 0, );
  
  self->Xe_s                  = ncm_spline_cubic_notaknot_new ();
  self->Xe_reion_s            = ncm_spline_cubic_notaknot_new ();
  self->Xe_recomb_s           = ncm_spline_cubic_notaknot_new ();
  self->XHII_s                = ncm_spline_cubic_notaknot_new ();
  self->XHeII_s               = ncm_spline_cubic_notaknot_new ();
}

static void
_nc_recomb_seager_set_property (GObject* object, guint prop_id, const GValue* value, GParamSpec* pspec)
{
  NcRecombSeager *recomb_seager = NC_RECOMB_SEAGER (object);
  g_return_if_fail (NC_IS_RECOMB_SEAGER (object));

  switch (prop_id)
  {
    case PROP_OPTS:
      nc_recomb_seager_set_options (recomb_seager, g_value_get_flags (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_recomb_seager_get_property (GObject* object, guint prop_id, GValue* value, GParamSpec* pspec)
{
  NcRecombSeager *recomb_seager = NC_RECOMB_SEAGER (object);
  g_return_if_fail (NC_IS_RECOMB_SEAGER (object));

  switch (prop_id)
  {
    case PROP_OPTS:
      g_value_set_flags (value, nc_recomb_seager_get_options (recomb_seager));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_recomb_seager_constructed (GObject* object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_recomb_seager_parent_class)->constructed (object);
  {
  }
}

static void
_nc_recomb_seager_dispose (GObject* object)
{
  NcRecombSeager *recomb_seager = NC_RECOMB_SEAGER (object);
  NcRecombSeagerPrivate * const self = recomb_seager->priv;

  ncm_spline_clear (&self->Xe_s);
  ncm_spline_clear (&self->Xe_reion_s);
  ncm_spline_clear (&self->Xe_recomb_s);
  ncm_spline_clear (&self->XHII_s);
  ncm_spline_clear (&self->XHeII_s);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_recomb_seager_parent_class)->dispose (object);
}

static void
_nc_recomb_seager_finalize (GObject* object)
{
  NcRecombSeager *recomb_seager = NC_RECOMB_SEAGER (object);
  NcRecombSeagerPrivate * const self = recomb_seager->priv;

  N_VDestroy (self->y);
  N_VDestroy (self->y0);
  N_VDestroy (self->abstol);

  if (self->A != NULL)
  {
    SUNMatDestroy (self->A);
    self->A = NULL;
  }
  
  if (self->LS != NULL)
  {
    SUNLinSolFree (self->LS);
    self->LS = NULL;
  }

  CVodeFree (&self->cvode);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_recomb_seager_parent_class)->finalize (object);
}

static void _nc_recomb_seager_prepare (NcRecomb *recomb, NcHICosmo *cosmo);
static gdouble _nc_recomb_seager_Xe (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
static gdouble _nc_recomb_seager_XHII (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
static gdouble _nc_recomb_seager_XHeII (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);

static void
nc_recomb_seager_class_init (NcRecombSeagerClass *klass)
{
  GObjectClass* object_class  = G_OBJECT_CLASS (klass);
  NcRecombClass* recomb_class = NC_RECOMB_CLASS (klass);

  object_class->set_property  = &_nc_recomb_seager_set_property;
  object_class->get_property  = &_nc_recomb_seager_get_property;
  object_class->constructed   = &_nc_recomb_seager_constructed;
  object_class->finalize      = &_nc_recomb_seager_finalize;
  object_class->dispose       = &_nc_recomb_seager_dispose;

  /**
   * NcRecombSeager:options:
   *
   * Integration options.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_OPTS,
                                   g_param_spec_flags ("options",
                                                       NULL,
                                                       "Integration options",
                                                       NC_TYPE_RECOMB_SEAGER_OPT, NC_RECOM_SEAGER_OPT_ALL,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  recomb_class->prepare = &_nc_recomb_seager_prepare;
  recomb_class->Xe      = &_nc_recomb_seager_Xe;
  recomb_class->XHII    = &_nc_recomb_seager_XHII;
  recomb_class->XHeII   = &_nc_recomb_seager_XHeII;  
}

static gdouble
_nc_recomb_seager_Xe (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda)
{
  NcRecombSeager *recomb_seager = NC_RECOMB_SEAGER (recomb);
  NcRecombSeagerPrivate * const self = recomb_seager->priv;
  NCM_UNUSED (cosmo);
  return ncm_spline_eval (self->Xe_s, lambda);
}

static gdouble
_nc_recomb_seager_XHII (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda)
{
  NcRecombSeager *recomb_seager = NC_RECOMB_SEAGER (recomb);
  NcRecombSeagerPrivate * const self = recomb_seager->priv;
  NCM_UNUSED (cosmo);
  return ncm_spline_eval (self->XHII_s, lambda);
}

static gdouble
_nc_recomb_seager_XHeII (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda)
{
  NcRecombSeager *recomb_seager = NC_RECOMB_SEAGER (recomb);
  NcRecombSeagerPrivate * const self = recomb_seager->priv;
  NCM_UNUSED (cosmo);
  return ncm_spline_eval (self->XHeII_s, lambda);
}

static gdouble nc_recomb_seager_HII_ion_rate (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHI, const gdouble XHII, const gdouble Tm, const gdouble XHeI, const gdouble XHeII, const gdouble x);
static gdouble nc_recomb_seager_HeII_ion_rate (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHI, const gdouble XHII, const gdouble Tm, const gdouble XHeI, const gdouble XHeII, const gdouble x);
static gdouble nc_recomb_seager_Tm_dx (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHI, const gdouble XHII, const gdouble Tm, const gdouble XHeI, const gdouble XHeII, const gdouble x);

static void nc_recomb_seager_HII_ion_rate_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHI, const gdouble XHII, const gdouble Tm, const gdouble XHeI, const gdouble XHeII, const gdouble x, gdouble* grad);
static void nc_recomb_seager_HeII_ion_rate_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHI, const gdouble XHII, const gdouble Tm, const gdouble XHeI, const gdouble XHeII, const gdouble x, gdouble* grad);
static void nc_recomb_seager_Tm_dx_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHI, const gdouble XHII, const gdouble Tm, const gdouble XHeI, const gdouble XHeII, const gdouble x, gdouble* grad);

static gdouble
_nc_recomb_seager_K_HI_2p_2Pmean (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble H)
{
  const gdouble K_HI_2p_2Pmean = ncm_c_HI_Lyman_wl3_8pi_2p_2Pmean () / H;

  return K_HI_2p_2Pmean;
}

static gdouble
_nc_recomb_seager_K_HI_2p_2Pmean_gcor (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble H)
{
  NcRecombSeagerPrivate * const self = recomb_seager->priv;
  
  const gdouble lnx = log (x);
  const gdouble K_gcor = 1.0 +
    self->AGauss1 * exp (-gsl_pow_2 ((lnx - self->zGauss1) / self->wGauss1)) +
    self->AGauss2 * exp (-gsl_pow_2 ((lnx - self->zGauss2) / self->wGauss2));
  const gdouble K_HI_2p_2Pmean = ncm_c_HI_Lyman_wl3_8pi_2p_2Pmean () * K_gcor / H;

  return K_HI_2p_2Pmean;
}

static gdouble
_nc_recomb_seager_KX_HeI_2p_1P1 (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble XHI, const gdouble Tm, const gdouble XHeI, const gdouble H, const gdouble n_H)
{
  const gdouble K_HeI = ncm_c_HeI_Lyman_wl3_8pi_2p_1P1 () * XHeI / H;

  return K_HeI;
}

static void
_nc_recomb_seager_KX_HeI_2p_1P1_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble XHI, const gdouble Tm, const gdouble XHeI, const gdouble H, const gdouble n_H, gdouble grad[3])
{
  grad[0] = 0.0;
  grad[1] = 0.0;
  grad[2] = -ncm_c_HeI_Lyman_wl3_8pi_2p_1P1 () / H;

  return;
}

static gdouble
_nc_recomb_seager_KX_HeI_2p_1P1_sobolev (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble XHI, const gdouble Tm, const gdouble XHeI, const gdouble H, const gdouble n_H)
{
  NcRecombSeagerPrivate * const self = recomb_seager->priv;
  const gdouble K_HeI    = ncm_c_HeI_Lyman_wl3_8pi_2p_1P1 () / H;
  const gdouble tau_XHeI = K_HeI * 3.0 * self->A2P_s * n_H;
  const gdouble tau      = tau_XHeI * XHeI;
  const gdouble P        = gsl_sf_exprel (-tau);

  return K_HeI / (P * tau_XHeI);
}

static void
_nc_recomb_seager_KX_HeI_2p_1P1_sobolev_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble XHI, const gdouble Tm, const gdouble XHeI, const gdouble H, const gdouble n_H, gdouble grad[3])
{
  NcRecombSeagerPrivate * const self = recomb_seager->priv;
  const gdouble K_HeI    = ncm_c_HeI_Lyman_wl3_8pi_2p_1P1 () / H;
  const gdouble tau_XHeI = K_HeI * 3.0 * self->A2P_s * n_H;
  const gdouble tau      = tau_XHeI * XHeI;
  const gdouble P        = gsl_sf_exprel (-tau);
  const gdouble dP_dtau  = -0.5 * gsl_sf_hyperg_1F1_int (2, 3, -tau);

  grad[0] = 0.0;
  grad[1] = 0.0;
  grad[2] = K_HeI * dP_dtau / (P * P);

  return;
}

static gdouble
_nc_recomb_seager_KX_HeI_2p_1P1_sobolev_cont (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble XHI, const gdouble Tm, const gdouble XHeI, const gdouble H, const gdouble n_H)
{
  NcRecombSeagerPrivate * const self = recomb_seager->priv;
  const gdouble K_HeI     = ncm_c_HeI_Lyman_wl3_8pi_2p_1P1 () / H;
  const gdouble tau_XHeI  = K_HeI * 3.0 * self->A2P_s * n_H;
  const gdouble tau       = tau_XHeI * XHeI;
  const gdouble P         = gsl_sf_exprel (-tau);

  const gdouble Doppler   = sqrt (2.0 * ncm_c_pi () * ncm_c_kb () * Tm / ncm_c_rest_energy_4He ());
  const gdouble gamma2P_s = H * tau / (n_H * ncm_c_c () * Doppler * self->sigma_He_2P_s * fabs (XHI));
  const gdouble A_Hcon    = 1.0 / (1.0 + self->Pb * pow (gamma2P_s, self->Qb));

  return K_HeI / ((P + A_Hcon) * tau_XHeI);
}

static void
_nc_recomb_seager_KX_HeI_2p_1P1_sobolev_cont_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble XHI, const gdouble Tm, const gdouble XHeI, const gdouble H, const gdouble n_H, gdouble grad[3])
{
  NcRecombSeagerPrivate * const self = recomb_seager->priv;
  const gdouble K_HeI           = ncm_c_HeI_Lyman_wl3_8pi_2p_1P1 () / H;
  const gdouble tau_XHeI        = K_HeI * 3.0 * self->A2P_s * n_H;
  const gdouble tau             = tau_XHeI * XHeI;
  const gdouble P               = gsl_sf_exprel (-tau);
  const gdouble dP_dtau         = -0.5 * gsl_sf_hyperg_1F1_int (2, 3, -tau);

  const gdouble Doppler         = sqrt (2.0 * ncm_c_pi () * ncm_c_kb () * Tm / ncm_c_rest_energy_4He ());
  const gdouble gamma2P_s_XHeI  = H * tau_XHeI / (n_H * ncm_c_c () * Doppler * self->sigma_He_2P_s * fabs (XHI));
  const gdouble gamma2P_s       = gamma2P_s_XHeI * XHeI;
  const gdouble dgamma2P_s_dTm  = -0.5 * gamma2P_s / Tm;
  const gdouble dgamma2P_s_dXHI = -gamma2P_s / XHI;
  const gdouble A_Hcon_f        = self->Pb * pow (gamma2P_s, self->Qb);
  const gdouble A_Hcon          = 1.0 / (1.0 + A_Hcon_f);
  const gdouble dA_Hcon_dgamma  = -A_Hcon * A_Hcon * self->Qb * A_Hcon_f / gamma2P_s;

  const gdouble denom           = gsl_pow_2 (P + A_Hcon) * tau_XHeI;

  grad[0] = K_HeI * dA_Hcon_dgamma * dgamma2P_s_dXHI / denom;
  grad[1] = -K_HeI * dA_Hcon_dgamma * dgamma2P_s_dTm / denom;
  grad[2] = K_HeI * (dP_dtau * tau_XHeI + dA_Hcon_dgamma * gamma2P_s_XHeI) / denom;

  return;
}

static gdouble
_nc_recomb_seager_KX_HeI_2p_3Pmean_sobolev (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble XHI, const gdouble Tm, const gdouble XHeI, const gdouble H, const gdouble n_H)
{
  NcRecombSeagerPrivate * const self = recomb_seager->priv;
  const gdouble K_HeI    = ncm_c_HeI_Lyman_wl3_8pi_2p_3Pmean () / H;
  const gdouble tau_XHeI = K_HeI * 3.0 * self->A2P_t * n_H;
  const gdouble tau      = tau_XHeI * XHeI;
  const gdouble P        = gsl_sf_exprel (-tau);

  return K_HeI / (P * tau_XHeI);
}

static void
_nc_recomb_seager_KX_HeI_2p_3Pmean_sobolev_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble XHI, const gdouble Tm, const gdouble XHeI, const gdouble H, const gdouble n_H, gdouble grad[3])
{
  NcRecombSeagerPrivate * const self = recomb_seager->priv;
  const gdouble K_HeI    = ncm_c_HeI_Lyman_wl3_8pi_2p_3Pmean () / H;
  const gdouble tau_XHeI = K_HeI * 3.0 * self->A2P_t * n_H;
  const gdouble tau      = tau_XHeI * XHeI;
  const gdouble P        = gsl_sf_exprel (-tau);
  const gdouble dP_dtau  = -0.5 * gsl_sf_hyperg_1F1_int (2, 3, -tau);

  grad[0] = 0.0;
  grad[1] = 0.0;
  grad[2] = K_HeI * dP_dtau / (P * P);

  return;
}

static gdouble
_nc_recomb_seager_KX_HeI_2p_3Pmean_sobolev_cont (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble XHI, const gdouble Tm, const gdouble XHeI, const gdouble H, const gdouble n_H)
{
  NcRecombSeagerPrivate * const self = recomb_seager->priv;
  const gdouble K_HeI     = ncm_c_HeI_Lyman_wl3_8pi_2p_3Pmean () / H;
  const gdouble tau_XHeI  = K_HeI * 3.0 * self->A2P_t * n_H;
  const gdouble tau       = tau_XHeI * XHeI;
  const gdouble P         = gsl_sf_exprel (-tau);
  const gdouble one_3     = 1.0 / 3.0;

  const gdouble Doppler   = sqrt (2.0 * ncm_c_pi () * ncm_c_kb () * Tm / ncm_c_rest_energy_4He ());
  const gdouble gamma2P_t = H * tau / (n_H * ncm_c_c () * Doppler * self->sigma_He_2P_t * fabs (XHI));
  const gdouble A_Hcon_t  = one_3 / (1.0 + self->Pb_t * pow (gamma2P_t, self->Qb_t));

  return K_HeI / ((P + A_Hcon_t) * tau_XHeI);
}

static void
_nc_recomb_seager_KX_HeI_2p_3Pmean_sobolev_cont_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble XHI, const gdouble Tm, const gdouble XHeI, const gdouble H, const gdouble n_H, gdouble grad[3])
{
  NcRecombSeagerPrivate * const self = recomb_seager->priv;
  const gdouble K_HeI           = ncm_c_HeI_Lyman_wl3_8pi_2p_3Pmean () / H;
  const gdouble tau_XHeI        = K_HeI * 3.0 * self->A2P_t * n_H;
  const gdouble tau             = tau_XHeI * XHeI;
  const gdouble P               = gsl_sf_exprel (-tau);
  const gdouble one_3           = 1.0 / 3.0;
  const gdouble dP_dtau         = -0.5 * gsl_sf_hyperg_1F1_int (2, 3, -tau);

  const gdouble Doppler         = sqrt (2.0 * ncm_c_pi () * ncm_c_kb () * Tm / ncm_c_rest_energy_4He ());
  const gdouble gamma2P_t_XHeI  = H * tau_XHeI / (n_H * ncm_c_c () * Doppler * self->sigma_He_2P_t * fabs (XHI));
  const gdouble gamma2P_t       = gamma2P_t_XHeI * XHeI;
  const gdouble dgamma2P_t_dTm  = -0.5 * gamma2P_t / Tm;
  const gdouble dgamma2P_t_dXHI = -gamma2P_t / XHI;
  const gdouble A_Hcon_f        = self->Pb_t * pow (gamma2P_t, self->Qb_t);
  const gdouble A_Hcon          = one_3 / (1.0 + A_Hcon_f);
  const gdouble dA_Hcon_dgamma  = -3.0 * A_Hcon * A_Hcon * self->Qb_t * A_Hcon_f / gamma2P_t;

  const gdouble denom           = gsl_pow_2 (P + A_Hcon) * tau_XHeI;

  grad[0] = K_HeI * dA_Hcon_dgamma * dgamma2P_t_dXHI / denom;
  grad[1] = -K_HeI * dA_Hcon_dgamma * dgamma2P_t_dTm / denom;
  grad[2] = K_HeI * (dP_dtau * tau_XHeI + dA_Hcon_dgamma * gamma2P_t_XHeI) / denom;

  return;
}

typedef struct _NcRecombSeagerParams
{
  NcRecombSeager *recomb_seager;
  NcHICosmo *cosmo;
} NcRecombSeagerParams;

static gint
H_ion_full_f (realtype lambda, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcRecombSeagerParams* rsp = (NcRecombSeagerParams*)f_data;
  const gdouble x     = exp (-lambda);
  const gdouble XHII  = NV_Ith_S (y, 0);
  const gdouble Tm    = NV_Ith_S (y, 1);
  const gdouble XHe   = nc_hicosmo_XHe (rsp->cosmo);
  const gdouble XHeII = NV_Ith_S (y, 2);
  const gdouble XHI   = (1.0 - XHII);
  const gdouble XHeI  = (XHe - XHeII);

  NV_Ith_S (ydot, 0) = -x * nc_recomb_seager_HII_ion_rate (rsp->recomb_seager, rsp->cosmo, XHI, XHII, Tm, XHeI, XHeII, x);
  NV_Ith_S (ydot, 1) = -x * nc_recomb_seager_Tm_dx (rsp->recomb_seager, rsp->cosmo, XHI, XHII, Tm, XHeI, XHeII, x);
  NV_Ith_S (ydot, 2) = -x * nc_recomb_seager_HeII_ion_rate (rsp->recomb_seager, rsp->cosmo, XHI, XHII, Tm, XHeI, XHeII, x);

  return GSL_SUCCESS;
}

static gint
H_ion_full_J (realtype lambda, N_Vector y, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcRecombSeagerParams* rsp = (NcRecombSeagerParams*)jac_data;
  const gdouble x     = exp (-lambda);
  const gdouble XHII  = NV_Ith_S (y, 0);
  const gdouble Tm    = NV_Ith_S (y, 1);
  const gdouble XHe   = nc_hicosmo_XHe (rsp->cosmo);
  const gdouble XHeII = NV_Ith_S (y, 2);
  const gdouble XHI   = (1.0 - XHII);
  const gdouble XHeI  = (XHe - XHeII);
  gdouble grad[3];

  NCM_UNUSED (fy);
  NCM_UNUSED (tmp1);
  NCM_UNUSED (tmp2);
  NCM_UNUSED (tmp3);

  nc_recomb_seager_HII_ion_rate_grad (rsp->recomb_seager, rsp->cosmo, XHI, XHII, Tm, XHeI, XHeII, x, grad);
  SUN_DENSE_ACCESS (J, 0, 0) = -x * grad[0];
  SUN_DENSE_ACCESS (J, 0, 1) = -x * grad[1];
  SUN_DENSE_ACCESS (J, 0, 2) = -x * grad[2];

  nc_recomb_seager_Tm_dx_grad (rsp->recomb_seager, rsp->cosmo, XHI, XHII, Tm, XHeI, XHeII, x, grad);
  SUN_DENSE_ACCESS (J, 1, 0) = -x * grad[0];
  SUN_DENSE_ACCESS (J, 1, 1) = -x * grad[1];
  SUN_DENSE_ACCESS (J, 1, 2) = -x * grad[2];

  nc_recomb_seager_HeII_ion_rate_grad (rsp->recomb_seager, rsp->cosmo, XHI, XHII, Tm, XHeI, XHeII, x, grad);
  SUN_DENSE_ACCESS (J, 2, 0) = -x * grad[0];
  SUN_DENSE_ACCESS (J, 2, 1) = -x * grad[1];
  SUN_DENSE_ACCESS (J, 2, 2) = -x * grad[2];

  return 0;
}

static gdouble
_nc_recomb_He_fully_ionized_Xe (gdouble lambda, gpointer p)
{
  NcHICosmo *cosmo = NC_HICOSMO (p);
  return nc_recomb_He_fully_ionized_Xe (cosmo, exp (-lambda));
}

static gdouble
_nc_recomb_Xe_reion (gdouble lambda, gpointer p)
{
  gpointer *m = (gpointer *) p;
  NcRecombSeager *recomb_seager = NC_RECOMB_SEAGER (m[0]);
  NcRecombSeagerPrivate * const self = recomb_seager->priv;
  NcHICosmo *cosmo              = NC_HICOSMO (m[1]);
  NcHIReion *reion              = NC_HIREION (m[2]);
  
  const gdouble Xe_recomb = ncm_spline_eval (self->Xe_recomb_s, lambda);
  const gdouble Xe_reion  = nc_hireion_get_Xe (reion, cosmo, lambda, Xe_recomb);

  return Xe_reion;
}

static void
_nc_recomb_seager_prepare (NcRecomb *recomb, NcHICosmo *cosmo)
{
  NcRecombSeager *recomb_seager = NC_RECOMB_SEAGER (recomb);
  NcRecombSeagerPrivate * const self = recomb_seager->priv;
  
  const gdouble XHe          = nc_hicosmo_XHe (cosmo);
  const gdouble x_HeIII      = nc_recomb_HeII_ion_saha_x_by_HeIII_He (cosmo, recomb->init_frac);
  const gdouble lambdai      = recomb->lambdai;
  const gdouble lambda_HeIII = -log (x_HeIII);
  const gdouble lambdaf      = -log (1.0);
  NcHIReion *reion           = NC_HIREION (ncm_model_peek_submodel_by_mid (NCM_MODEL (cosmo), nc_hireion_id ()));
  NcRecombSeagerParams pparams = { recomb_seager, cosmo };
  gsl_function F;

  F.function = &_nc_recomb_He_fully_ionized_Xe;
  F.params   = cosmo;

  ncm_spline_set_func (self->Xe_recomb_s, NCM_SPLINE_FUNCTION_SPLINE,
                       &F, lambdai, lambda_HeIII, 0, recomb->prec);

  /*****************************************************************************
   * Assuming hydrogen is completly ionized and no more double ionized helium
   * i.e., $X_\HeIII = 0$.
   ****************************************************************************/
  {
    const gdouble T0   = nc_hicosmo_T_gamma0 (cosmo);
    const gdouble XHII = 1.0;
    gdouble XHeII, Tm, XeXHeII_XHeI;

    Tm           = T0 * x_HeIII;
    XeXHeII_XHeI = nc_recomb_HeI_ion_saha (cosmo, x_HeIII);
    XHeII        = (XeXHeII_XHeI + 1.0) * ncm_util_sqrt1px_m1 (4.0 * XHe * XeXHeII_XHeI / gsl_pow_2 (XeXHeII_XHeI + 1.0)) / 2.0;

    NV_Ith_S (self->y0, 0) = XHII;
    NV_Ith_S (self->y0, 1) = Tm;
    NV_Ith_S (self->y0, 2) = XHeII;
  }

  if (!self->init)
  {
    gint flag = CVodeInit (self->cvode, &H_ion_full_f, lambda_HeIII, self->y0);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
    self->init = TRUE;
  }
  else
  {
    gint flag = CVodeReInit (self->cvode, lambda_HeIII, self->y0);
    NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );
  }

  {
    const gdouble reltol = GSL_MIN (recomb->prec, 1e-11);
    gint flag;

    NV_Ith_S (self->abstol, 0) = GSL_MIN (reltol, recomb->prec * 1e-5);
    NV_Ith_S (self->abstol, 2) = GSL_MIN (reltol, recomb->prec * 1e-5);
    NV_Ith_S (self->abstol, 1) = 0.0;

    flag = CVodeSVtolerances (self->cvode, reltol, self->abstol);
    NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );

    flag = CVodeSetUserData (self->cvode, &pparams);
    NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

    flag = CVodeSetMaxNumSteps (self->cvode, 100000);
    NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

    flag = CVodeSetLinearSolver (self->cvode, self->LS, self->A);
    NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, );

    flag = CVodeSetJacFn (self->cvode, &H_ion_full_J);
    NCM_CVODE_CHECK (&flag, "CVodeSetJacFn", 1, );

    flag = CVodeSetStopTime (self->cvode, lambdaf);
    NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

    flag = CVodeSetMaxErrTestFails (self->cvode, 14);
    NCM_CVODE_CHECK (&flag, "CVodeSetMaxErrTestFails", 1, );
  }

  {
    NcmVector *lambda_v = ncm_spline_get_xv (self->Xe_recomb_s);
    NcmVector *Xe_v     = ncm_spline_get_yv (self->Xe_recomb_s);
    GArray *lambda_a    = ncm_vector_get_array (lambda_v);
    GArray *Xe_a        = ncm_vector_get_array (Xe_v);
    gdouble lambda_last = lambda_HeIII;
    GArray *XHII_a, *XHeII_a;

    ncm_vector_free (lambda_v);
    ncm_vector_free (Xe_v);
    
    if (ncm_spline_is_empty (self->XHII_s))
    {
      XHII_a = g_array_new (FALSE, FALSE, sizeof (gdouble));
    }
    else
    {
      NcmVector *XHII_v = ncm_spline_get_yv (self->XHII_s);

      XHII_a = ncm_vector_get_array (XHII_v);
      g_array_set_size (XHII_a, 0);
      
      ncm_vector_free (XHII_v);
    }

    if (ncm_spline_is_empty (self->XHeII_s))
    {
      XHeII_a = g_array_new (FALSE, FALSE, sizeof (gdouble));
    }
    else
    {
      NcmVector *XHeII_v = ncm_spline_get_yv (self->XHeII_s);
      
      XHeII_a = ncm_vector_get_array (XHeII_v);
      g_array_set_size (XHeII_a, 0);

      ncm_vector_free (XHeII_v);
    }
    
    {
      const guint lenXe = ncm_spline_get_len (self->Xe_recomb_s);
      guint i;

      for (i = 0; i < lenXe; i++)
      {
        const gdouble lambda_i = ncm_vector_get (lambda_v, i);
        const gdouble x        = exp (-lambda_i);

        const gdouble XHII     = nc_recomb_equilibrium_XHII (recomb, cosmo, x);
        const gdouble XHeII    = nc_recomb_equilibrium_XHeII (recomb, cosmo, x);

        g_array_append_val (XHII_a, XHII);
        g_array_append_val (XHeII_a, XHeII);
      }
    }

    while (TRUE)
    {
      gdouble lambda_i;
      gint flag = CVode (self->cvode, lambdaf, self->y, &lambda_i, CV_ONE_STEP);
      NCM_CVODE_CHECK (&flag, "CVode111", 1, );
      
      {
        const gdouble XHII  = NV_Ith_S (self->y, 0);
        const gdouble XHeII = NV_Ith_S (self->y, 2);
        const gdouble Xe    = XHII + XHeII;

        if (fabs ((lambda_last - lambda_i) / lambda_last) > 1.0e-7)
        {
          g_array_append_val (lambda_a, lambda_i);
          g_array_append_val (Xe_a, Xe);

          g_array_append_val (XHII_a, XHII);
          g_array_append_val (XHeII_a, XHeII);

          lambda_last = lambda_i;
        }

        if (lambda_i == lambdaf)
          break;
      }
    }

    ncm_spline_set_array (self->Xe_recomb_s, lambda_a, Xe_a, TRUE);
    ncm_spline_set_array (self->XHII_s,      lambda_a, XHII_a, TRUE);
    ncm_spline_set_array (self->XHeII_s,     lambda_a, XHeII_a, TRUE);

    if (reion != NULL)
    {
      gpointer m[3] = {recomb, cosmo, reion};

      F.function = &_nc_recomb_Xe_reion;
      F.params   = m;
            
      ncm_spline_set_func (self->Xe_reion_s, NCM_SPLINE_FUNCTION_SPLINE,
                           &F, lambdai, lambdaf, 0, recomb->prec);

      ncm_spline_clear (&self->Xe_s);
      self->Xe_s = ncm_spline_ref (self->Xe_reion_s);
    }
    else
    {
      ncm_spline_clear (&self->Xe_s);
      self->Xe_s = ncm_spline_ref (self->Xe_recomb_s);
    }  
    
    g_array_unref (lambda_a);
    g_array_unref (Xe_a);
    g_array_unref (XHII_a);
    g_array_unref (XHeII_a);
  }

  recomb->tau_s          = ncm_spline_copy_empty (self->Xe_s);
  recomb->dtau_dlambda_s = ncm_spline_copy_empty (self->Xe_s);

  _nc_recomb_prepare_tau_splines (recomb, cosmo);
  _nc_recomb_prepare_redshifts (recomb, cosmo);
}

static gdouble
nc_recomb_seager_HII_ion_rate (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHI, const gdouble XHII, const gdouble Tm, const gdouble XHeI, const gdouble XHeII, const gdouble x)
{
  NcRecombSeagerPrivate * const self = recomb_seager->priv;
  const gdouble Xe            = XHII + XHeII;
  const gdouble x2            = x * x;
  const gdouble x3            = x2 * x;

  const gdouble alpha_H       = nc_recomb_seager_pequignot_HI_case_B (recomb_seager, cosmo, Tm);
  const gdouble n_H0          = nc_hicosmo_H_number_density (cosmo);
  const gdouble n_H           = n_H0 * x3;
  const gdouble H             = nc_hicosmo_H (cosmo, x - 1.0) / ncm_c_kpc ();
  const gdouble Tm3_2         = sqrt (gsl_pow_3 (Tm));
  const gdouble K_HI          = self->K_HI_2p_2Pmean (recomb_seager, cosmo, x, H);
  const gdouble Lambda_H      = ncm_c_decay_H_rate_2s_1s ();
  const gdouble nKX           = n_H * K_HI * XHI;
  const gdouble B_HI_1s_2S0_5 = ncm_c_boltzmann_factor_HI_1s_2S0_5 (Tm) * Tm3_2;
  const gdouble B_HI_2s_2S0_5 = ncm_c_boltzmann_factor_HI_2s_2S0_5 (Tm) * Tm3_2;

  const gdouble saha_factor   = (XHII * Xe * n_H - XHI * B_HI_1s_2S0_5) / (H * x);

  const gdouble R_n           = (nKX * Lambda_H + 1.0) * alpha_H;
  const gdouble R_d           = nKX * (Lambda_H + B_HI_2s_2S0_5 * alpha_H) + 1.0;

  const gdouble R             = R_n / R_d;

  return saha_factor * R;
}

static void
nc_recomb_seager_HII_ion_rate_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHI, const gdouble XHII, const gdouble Tm, const gdouble XHeI, const gdouble XHeII, const gdouble x, gdouble* grad)
{
  NcRecombSeagerPrivate * const self = recomb_seager->priv;
  
  const gdouble Xe                  = XHII + XHeII;
  const gdouble x2                  = x * x;
  const gdouble x3                  = x2 * x;

  const gdouble alpha_H             = nc_recomb_seager_pequignot_HI_case_B (recomb_seager, cosmo, Tm);
  const gdouble n_H0                = nc_hicosmo_H_number_density (cosmo);
  const gdouble n_H                 = n_H0 * x3;
  const gdouble H                   = nc_hicosmo_H (cosmo, x - 1.0) / ncm_c_kpc ();
  const gdouble Tm3_2               = sqrt (gsl_pow_3 (Tm));
  const gdouble K_HI                = self->K_HI_2p_2Pmean (recomb_seager, cosmo, x, H);
  const gdouble Lambda_H            = ncm_c_decay_H_rate_2s_1s ();
  const gdouble nKX                 = n_H * K_HI * XHI;
  const gdouble B_HI_1s_2S0_5       = ncm_c_boltzmann_factor_HI_1s_2S0_5 (Tm) * Tm3_2;
  const gdouble B_HI_2s_2S0_5       = ncm_c_boltzmann_factor_HI_2s_2S0_5 (Tm) * Tm3_2;
  const gdouble Hx                  = H * x;

  const gdouble saha_factor         = (XHII * Xe * n_H - XHI * B_HI_1s_2S0_5) / Hx;

  const gdouble R_n                 = (nKX * Lambda_H + 1.0) * alpha_H;
  const gdouble R_d                 = nKX * (Lambda_H + B_HI_2s_2S0_5 * alpha_H) + 1.0;
  const gdouble R                   = R_n / R_d;

  const gdouble R_d_pow_2           = R_d * R_d;
  const gdouble alpha_H2            = alpha_H * alpha_H;
  const gdouble Tm_pow_2            = Tm * Tm;
  const gdouble dalpha_H_dTm        = nc_recomb_seager_pequignot_HI_case_B_dTm (recomb_seager, cosmo, Tm);
  const gdouble dB_HI_1s_2S0_5_dTm  = B_HI_1s_2S0_5 * (ncm_c_HI_ion_E_1s_2S0_5 () / (ncm_c_kb () * Tm_pow_2) + 1.5 / Tm);

  const gdouble dsaha_factor_dXHII  = ((XHII + Xe) * n_H + B_HI_1s_2S0_5) / Hx;
  const gdouble dsaha_factor_dTm    = -XHI * dB_HI_1s_2S0_5_dTm / Hx;
  const gdouble dsaha_factor_dXHeII = XHII * n_H / Hx;

  const gdouble dR_dKX              = -B_HI_2s_2S0_5 * n_H * alpha_H2 / R_d_pow_2;
  const gdouble dR_dBalpha          = -R * nKX / R_d;
  const gdouble dBalpha_Tm          = B_HI_2s_2S0_5 * alpha_H * (ncm_c_HI_ion_E_2s_2S0_5 () / (ncm_c_kb () * Tm_pow_2) + 1.5 / Tm + dalpha_H_dTm / alpha_H);

  const gdouble dR_dXHII            = -dR_dKX * K_HI;
  const gdouble dR_dTm              = dR_dBalpha * dBalpha_Tm + R * dalpha_H_dTm / alpha_H;

  grad[0] = R * dsaha_factor_dXHII + saha_factor * dR_dXHII;
  grad[1] = R * dsaha_factor_dTm + saha_factor * dR_dTm;
  grad[2] = R * dsaha_factor_dXHeII;
}

static gdouble
nc_recomb_seager_HeII_ion_rate (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHI, const gdouble XHII, const gdouble Tm, const gdouble XHeI, const gdouble XHeII, const gdouble x)
{
  NcRecombSeagerPrivate * const self = recomb_seager->priv;
  
  const gdouble Xe              = XHII + XHeII;
  const gdouble x2              = x * x;
  const gdouble x3              = x2 * x;

  const gdouble alpha           = nc_recomb_seager_hummer_HeI_case_B (recomb_seager, cosmo, Tm);
  const gdouble n_H0            = nc_hicosmo_H_number_density (cosmo);
  const gdouble n_H             = n_H0 * x3;
  const gdouble H               = nc_hicosmo_H (cosmo, x - 1.0) / ncm_c_kpc ();
  const gdouble Tm3_2           = sqrt (gsl_pow_3 (Tm));
  const gdouble KX_HeI          = self->KX_HeI_2p_1P1 (recomb_seager, cosmo, x, XHI, Tm, XHeI, H, n_H);
  const gdouble Lambda_He       = ncm_c_decay_He_rate_2s_1s ();
  const gdouble nKX             = n_H * KX_HeI;
  const gdouble exp_E_HeI_2s_2p = exp (-ncm_c_HeI_Balmer_E_kb_2p_1P1_2s_1S0 () / Tm);
  const gdouble B_HeI_1s_1S0    = 4.0 * ncm_c_boltzmann_factor_HeI_1s_1S0 (Tm) * Tm3_2;
  const gdouble B_HeI_2s_1S0    = 4.0 * ncm_c_boltzmann_factor_HeI_2s_1S0 (Tm) * Tm3_2;

  const gdouble saha_factor     = (XHeII * Xe * n_H - XHeI * B_HeI_1s_1S0) / (H * x);

  const gdouble R_sd_n          = (nKX * Lambda_He + exp_E_HeI_2s_2p) * alpha;
  const gdouble R_sd_d          = nKX * (Lambda_He + B_HeI_2s_1S0 * alpha) + exp_E_HeI_2s_2p;

  const gdouble R_sd            = R_sd_n / R_sd_d;

  if (self->KX_HeI_2p_3Pmean == NULL)
  {
    return saha_factor * R_sd;
  }
  else
  {
    const gdouble alpha_t         = nc_recomb_seager_hummer_HeI_case_B_trip (recomb_seager, cosmo, Tm);
    const gdouble KX_HeI_t        = self->KX_HeI_2p_3Pmean (recomb_seager, cosmo, x, XHI, Tm, XHeI, H, n_H);
    const gdouble nKX_t           = n_H * KX_HeI_t;
    const gdouble B_HeI_2p_3Pmean = 4.0 * ncm_c_boltzmann_factor_HeI_2p_3Pmean (Tm) * Tm3_2;

    const gdouble R_t_n           = alpha_t;
    const gdouble R_t_d           = nKX_t * B_HeI_2p_3Pmean * alpha_t + 1.0;
    const gdouble R_t             = R_t_n / R_t_d;
    const gdouble R               = R_sd + R_t;

    return saha_factor * R;
  }
}

static void
nc_recomb_seager_HeII_ion_rate_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHI, const gdouble XHII, const gdouble Tm, const gdouble XHeI, const gdouble XHeII, const gdouble x, gdouble* grad)
{
  NcRecombSeagerPrivate * const self = recomb_seager->priv;
  
  const gdouble Xe                  = XHII + XHeII;
  const gdouble x2                  = x * x;
  const gdouble x3                  = x2 * x;

  const gdouble alpha               = nc_recomb_seager_hummer_HeI_case_B (recomb_seager, cosmo, Tm);
  const gdouble n_H0                = nc_hicosmo_H_number_density (cosmo);
  const gdouble n_H                 = n_H0 * x3;
  const gdouble H                   = nc_hicosmo_H (cosmo, x - 1.0) / ncm_c_kpc ();
  const gdouble Tm3_2               = sqrt (gsl_pow_3 (Tm));
  const gdouble KX_HeI              = self->KX_HeI_2p_1P1 (recomb_seager, cosmo, x, XHI, Tm, XHeI, H, n_H);
  const gdouble Lambda_He           = ncm_c_decay_He_rate_2s_1s ();
  const gdouble nKX                 = n_H * KX_HeI;
  const gdouble exp_E_HeI_2s_2p     = exp (-ncm_c_HeI_Balmer_E_kb_2p_1P1_2s_1S0 () / Tm);
  const gdouble B_HeI_1s_1S0        = 4.0 * ncm_c_boltzmann_factor_HeI_1s_1S0 (Tm) * Tm3_2;
  const gdouble B_HeI_2s_1S0        = 4.0 * ncm_c_boltzmann_factor_HeI_2s_1S0 (Tm) * Tm3_2;
  const gdouble Hx                  = H * x;

  const gdouble saha_factor         = (XHeII * Xe * n_H - XHeI * B_HeI_1s_1S0) / Hx;

  const gdouble R_sd_n              = (nKX * Lambda_He + exp_E_HeI_2s_2p) * alpha;
  const gdouble R_sd_d              = nKX * (Lambda_He + B_HeI_2s_1S0 * alpha) + exp_E_HeI_2s_2p;

  const gdouble R_sd                = R_sd_n / R_sd_d;

  const gdouble R_sd_d_pow_2        = R_sd_d * R_sd_d;
  const gdouble alpha2              = alpha * alpha;
  const gdouble Tm_pow_2            = Tm * Tm;
  const gdouble dalpha_dTm          = nc_recomb_seager_hummer_HeI_case_B_dTm (recomb_seager, cosmo, Tm);
  const gdouble dB_HeI_1s_1S0_dTm   = B_HeI_1s_1S0 * (ncm_c_HeI_ion_E_1s_1S0 () / (ncm_c_kb () * Tm_pow_2) + 1.5 / Tm);

  const gdouble dsaha_factor_dXHII  = XHeII * n_H / Hx;
  const gdouble dsaha_factor_dTm    = -XHeI * dB_HeI_1s_1S0_dTm / Hx;
  const gdouble dsaha_factor_dXHeII = ((XHeII + Xe) * n_H + B_HeI_1s_1S0) / Hx;

  const gdouble dR_sd_dKX              = -exp_E_HeI_2s_2p * B_HeI_2s_1S0 * n_H * alpha2 / R_sd_d_pow_2;
  const gdouble dR_sd_dexp_E_HeI_2s_2p = nKX * B_HeI_2s_1S0 * alpha2 / R_sd_d_pow_2;
  const gdouble dR_sd_dBalpha          = -R_sd * nKX / R_sd_d;

  const gdouble dexp_E_HeI_2s_2p_dTm   = exp_E_HeI_2s_2p * ncm_c_HeI_Balmer_E_kb_2p_1P1_2s_1S0 () / Tm_pow_2;
  const gdouble dBalpha_Tm             = B_HeI_2s_1S0 * alpha * (ncm_c_HeI_ion_E_2s_1S0 () / (ncm_c_kb () * Tm_pow_2) + 1.5 / Tm + dalpha_dTm / alpha);

  gdouble grad_temp[3];

  self->KX_HeI_2p_1P1_grad (recomb_seager, cosmo, x, XHI, Tm, XHeI, H, n_H, grad_temp);
  {
    const gdouble dR_sd_dXHII  = dR_sd_dKX * grad_temp[0];
    const gdouble dR_sd_dTm    = dR_sd_dKX * grad_temp[1] + dR_sd_dexp_E_HeI_2s_2p * dexp_E_HeI_2s_2p_dTm + dR_sd_dBalpha * dBalpha_Tm + R_sd * dalpha_dTm / alpha;
    const gdouble dR_sd_dXHeII = dR_sd_dKX * grad_temp[2];

    if (self->KX_HeI_2p_3Pmean == NULL)
    {
      grad[0] = R_sd * dsaha_factor_dXHII + saha_factor * dR_sd_dXHII;
      grad[1] = R_sd * dsaha_factor_dTm + saha_factor * dR_sd_dTm;
      grad[2] = R_sd * dsaha_factor_dXHeII + saha_factor * dR_sd_dXHeII;
    }
    else
    {
      const gdouble alpha_t         = nc_recomb_seager_hummer_HeI_case_B_trip (recomb_seager, cosmo, Tm);
      const gdouble KX_HeI_t        = self->KX_HeI_2p_3Pmean (recomb_seager, cosmo, x, XHI, Tm, XHeI, H, n_H);
      const gdouble nKX_t           = n_H * KX_HeI_t;
      const gdouble B_HeI_2p_3Pmean = 4.0 * ncm_c_boltzmann_factor_HeI_2p_3Pmean (Tm) * Tm3_2;

      const gdouble R_t_n           = alpha_t;
      const gdouble R_t_d           = nKX_t * B_HeI_2p_3Pmean * alpha_t + 1.0;
      const gdouble R_t             = R_t_n / R_t_d;
      const gdouble R               = R_sd + R_t;

      const gdouble R_t_d_pow_2     = R_t_d * R_t_d;
      const gdouble dalpha_t_dTm    = nc_recomb_seager_hummer_HeI_case_B_trip_dTm (recomb_seager, cosmo, Tm);

      const gdouble dR_t_dKX_t      = -R_t_n * B_HeI_2p_3Pmean * n_H * alpha_t / R_t_d_pow_2;
      const gdouble dR_t_dBalpha_t  = -R_t * nKX_t / R_t_d;
      const gdouble dBalpha_t_Tm    = B_HeI_2p_3Pmean * alpha_t * (ncm_c_HeI_ion_E_2p_3Pmean () / (ncm_c_kb () * Tm_pow_2) + 1.5 / Tm + dalpha_t_dTm / alpha_t);

      self->KX_HeI_2p_3Pmean_grad (recomb_seager, cosmo, x, XHI, Tm, XHeI, H, n_H, grad_temp);
      {
        const gdouble dR_t_dXHII  = dR_t_dKX_t * grad_temp[0];
        const gdouble dR_t_dTm    = dR_t_dKX_t * grad_temp[1] + dR_t_dBalpha_t * dBalpha_t_Tm + R_t * dalpha_t_dTm / alpha_t;
        const gdouble dR_t_dXHeII = dR_t_dKX_t * grad_temp[2];

        grad[0] = R * dsaha_factor_dXHII + saha_factor * (dR_sd_dXHII + dR_t_dXHII);
        grad[1] = R * dsaha_factor_dTm + saha_factor * (dR_sd_dTm + dR_t_dTm);
        grad[2] = R * dsaha_factor_dXHeII + saha_factor * (dR_sd_dXHeII + dR_t_dXHeII);
      }
    }
  }
}

static gdouble
nc_recomb_seager_Tm_dx (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHI, const gdouble XHII, const gdouble Tm, const gdouble XHeI, const gdouble XHeII, const gdouble x)
{
  const gdouble T0  = nc_hicosmo_T_gamma0 (cosmo);
  const gdouble T   = T0 * x;
  const gdouble T4  = gsl_pow_4 (T);
  const gdouble H   = nc_hicosmo_H (cosmo, x - 1.0) / ncm_c_kpc ();

  const gdouble f1  = (8.0 * ncm_c_thomson_cs () * ncm_c_AR () * T4 / (3.0 * ncm_c_c () * ncm_c_mass_e ())) / (H * x);

  const gdouble Xe  = XHII + XHeII;
  const gdouble XHe = nc_hicosmo_XHe (cosmo);
  const gdouble f2  = Xe * (Tm - T) / (1.0 + XHe + Xe);

  const gdouble f3  = 2.0 * Tm / x;

  return f1 * f2 + f3;
}

static void
nc_recomb_seager_Tm_dx_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHI, const gdouble XHII, const gdouble Tm, const gdouble XHeI, const gdouble XHeII, const gdouble x, gdouble* grad)
{
  const gdouble T0      = nc_hicosmo_T_gamma0 (cosmo);
  const gdouble T       = T0 * x;
  const gdouble H       = nc_hicosmo_H (cosmo, x - 1.0) / ncm_c_kpc ();
  const gdouble T4      = gsl_pow_4 (T);
  const gdouble f1      = (8.0 * ncm_c_thomson_cs () * ncm_c_AR () * T4 / (3.0 * ncm_c_c () * ncm_c_mass_e ())) / (H * x);

  const gdouble Xe      = XHII + XHeII;
  const gdouble XHe     = nc_hicosmo_XHe (cosmo);
  const gdouble f2      = Xe * (Tm - T) / (1.0 + XHe + Xe);

  const gdouble df2_dX  = (1.0 / Xe - 1.0 / (1.0 + XHe + Xe)) * f2;
  const gdouble ddXHII  = f1 * df2_dX;

  const gdouble ddTm    = f1 * Xe / (1.0 + XHe + Xe) + 2.0 / x;

  const gdouble ddXHeII = ddXHII;

  grad[0] = ddXHII;
  grad[1] = ddTm;
  grad[2] = ddXHeII;

  return;
}

/**
 * nc_recomb_seager_new:
 *
 * Creates a new #NcRecombSeager using default properties.
 *
 * Returns: (transfer full): a new #NcRecombSeager.
 */
NcRecombSeager*
nc_recomb_seager_new (void)
{
  return g_object_new (NC_TYPE_RECOMB_SEAGER,
                       NULL);
}

/**
 * nc_recomb_seager_new_full:
 * @init_frac: inital fraction of $X_{\HeIII}/X_{\He}$ where to start numerical integration
 * @zi: inital redshift
 * @prec: integration precision
 *
 * Creates a new #NcRecombSeager using @init_frac, @zi and @prec.
 *
 * Returns: (transfer full): a new #NcRecombSeager.
 */
NcRecombSeager*
nc_recomb_seager_new_full (gdouble init_frac, gdouble zi, gdouble prec)
{
  return g_object_new (NC_TYPE_RECOMB_SEAGER,
                       "init-frac", init_frac,
                       "zi", zi,
                       "prec", prec,
                       NULL);
}

/**
 * nc_recomb_seager_ref:
 * @recomb_seager: a #NcRecombSeager
 *
 * Increases the reference count of @recomb_seager.
 *
 * Returns: (transfer full): @recomb_seager.
 */
NcRecombSeager*
nc_recomb_seager_ref (NcRecombSeager *recomb_seager)
{
  return NC_RECOMB_SEAGER (g_object_ref (recomb_seager));
}

/**
 * nc_recomb_seager_free:
 * @recomb_seager: a #NcRecombSeager.
 *
 * Decreases the reference count of @recomb_seager.
 *
 */
void nc_recomb_seager_free (NcRecombSeager *recomb_seager)
{
  g_object_unref (recomb_seager);
}

/**
 * nc_recomb_seager_clear:
 * @recomb_seager: a #NcRecombSeager.
 *
 * Decreases the reference count of *@recomb_seager if
 * *@recomb_seager is not NULL, then sets *@recomb_seager to NULL.
 *
 */
void nc_recomb_seager_clear (NcRecombSeager** recomb_seager)
{
  g_clear_object (recomb_seager);
}

/**
 * nc_recomb_seager_set_options:
 * @recomb_seager: a #NcRecombSeager
 * @opts: a #NcRecombSeagerOpt
 *
 * Sets integration options #NcRecombSeagerOpt. To set the integration
 * options using the recfast compatible flags use nc_recomb_seager_set_switch().
 *
 */
void nc_recomb_seager_set_options (NcRecombSeager *recomb_seager, NcRecombSeagerOpt opts)
{
  NcRecombSeagerPrivate * const self = recomb_seager->priv;
  if (self->opts != opts)
  {
    NcRecomb *recomb = NC_RECOMB (recomb_seager);

    if (opts & NC_RECOM_SEAGER_OPT_HII_FUDGE_GAUSS_COR)
    {
      self->H_fudge = 1.125;
      self->K_HI_2p_2Pmean = _nc_recomb_seager_K_HI_2p_2Pmean_gcor;
    }
    else if (opts & NC_RECOM_SEAGER_OPT_HII_FUDGE)
    {
      self->H_fudge = 1.14;
      self->K_HI_2p_2Pmean = _nc_recomb_seager_K_HI_2p_2Pmean;
    }
    else
    {
      self->H_fudge = 1.0;
      self->K_HI_2p_2Pmean = _nc_recomb_seager_K_HI_2p_2Pmean;
    }

    if (opts & NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1)
    {
      if (opts & NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1_CO)
      {
        self->KX_HeI_2p_1P1 = &_nc_recomb_seager_KX_HeI_2p_1P1_sobolev_cont;
        self->KX_HeI_2p_1P1_grad = &_nc_recomb_seager_KX_HeI_2p_1P1_sobolev_cont_grad;
      }
      else
      {
        self->KX_HeI_2p_1P1 = &_nc_recomb_seager_KX_HeI_2p_1P1_sobolev;
        self->KX_HeI_2p_1P1_grad = &_nc_recomb_seager_KX_HeI_2p_1P1_sobolev_grad;
      }
    }
    else
    {
      self->KX_HeI_2p_1P1 = &_nc_recomb_seager_KX_HeI_2p_1P1;
      self->KX_HeI_2p_1P1_grad = &_nc_recomb_seager_KX_HeI_2p_1P1_grad;
    }

    if (opts & NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_3P012)
    {
      if (opts & NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_3P012_CO)
      {
        self->KX_HeI_2p_3Pmean = &_nc_recomb_seager_KX_HeI_2p_3Pmean_sobolev_cont;
        self->KX_HeI_2p_3Pmean_grad = &_nc_recomb_seager_KX_HeI_2p_3Pmean_sobolev_cont_grad;
      }
      else
      {
        self->KX_HeI_2p_3Pmean = &_nc_recomb_seager_KX_HeI_2p_3Pmean_sobolev;
        self->KX_HeI_2p_3Pmean_grad = &_nc_recomb_seager_KX_HeI_2p_3Pmean_sobolev_grad;
      }
    }
    else
    {
      self->KX_HeI_2p_3Pmean = NULL;
      self->KX_HeI_2p_3Pmean_grad = NULL;
    }

    ncm_model_ctrl_force_update (recomb->ctrl_cosmo);
    self->opts = opts;
  }
}

/**
 * nc_recomb_seager_set_switch:
 * @recomb_seager: a #NcRecombSeager
 * @H_switch: an integer between 0 and 1
 * @He_switch: an integer between 0 and 6
 *
 * Sets integration options #NcRecombSeagerOpt using the following map:
 *
 * - @H_switch == 0 => #NC_RECOM_SEAGER_OPT_HII_FUDGE;
 * - @H_switch == 1 => #NC_RECOM_SEAGER_OPT_HII_FUDGE_GAUSS_COR;
 *
 * - @He_switch == 0 => no additional flag;
 * - @He_switch == 1 => #NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1;
 * - @He_switch == 2 => #NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1 | #NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1_CO;
 * - @He_switch == 3 => #NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1 | #NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_3P012;
 * - @He_switch == 4 => #NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1 | #NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_3P012 | #NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_3P012_CO;
 * - @He_switch == 5 => #NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1 | #NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1_CO | #NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_3P012;
 * - @He_switch == 6 => #NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1 | #NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1_CO | #NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_3P012 | #NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_3P012_CO;
 *
 */
void nc_recomb_seager_set_switch (NcRecombSeager *recomb_seager, guint H_switch, guint He_switch)
{
  NcRecombSeagerOpt opts = 0;
  switch (H_switch)
  {
    case 0:
      opts |= NC_RECOM_SEAGER_OPT_HII_FUDGE;
      break;
    case 1:
      opts |= NC_RECOM_SEAGER_OPT_HII_FUDGE_GAUSS_COR;
      break;
    default:
      g_error ("nc_recomb_seager_set_options: unknown H_switch %u.\n", H_switch);
      break;
  }
  switch (He_switch)
  {
    case 0:
      break;
    case 1:
      opts |= NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1;
      break;
    case 2:
      opts |= NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1 | NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1_CO;
      break;
    case 3:
      opts |= NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1 | NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_3P012;
      break;
    case 4:
      opts |= NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1 | NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_3P012 | NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_3P012_CO;
      break;
    case 5:
      opts |= NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1 | NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1_CO | NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_3P012;
      break;
    case 6:
      opts |= NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1 | NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1_CO | NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_3P012 | NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_3P012_CO;
      break;
    default:
      g_error ("nc_recomb_seager_set_options: unknown He_switch %u.\n", He_switch);
      break;
  }
  nc_recomb_seager_set_options (recomb_seager, opts);
}

/**
 * nc_recomb_seager_get_options:
 * @recomb_seager: a #NcRecombSeager
 *
 * Gets integration options.
 *
 * Returns: currently used integration options.
 */
NcRecombSeagerOpt
nc_recomb_seager_get_options (NcRecombSeager *recomb_seager)
{
  NcRecombSeagerPrivate * const self = recomb_seager->priv;
  return self->opts;
}


/**
 * nc_recomb_seager_pequignot_HI_case_B:
 * @recomb_seager: a #NcRecombSeager
 * @cosmo: a #NcHICosmo
 * @Tm: the matter (baryons) temperature $T_m$
 *
 * The case B $\HyII$ recombination coefficient.
 *
 * The fitting formula of the case B recombination coefficient for $\HyII$ as
 * in [Pequignot (1991)][XPequignot1991].
 *
 * Returns: the value of the case B recombination coefficient for
 * $\HyII$, $\alpha_H$ .
 */
gdouble
nc_recomb_seager_pequignot_HI_case_B (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, gdouble Tm)
{
  NcRecombSeagerPrivate * const self = recomb_seager->priv;
  const gdouble F   = self->H_fudge; /* fudge factor */
  const gdouble G   = 1e-19;
  const gdouble a   = 4.309;
  const gdouble b   = -0.6166;
  const gdouble c   = 0.6703;
  const gdouble d   = 0.5300;
  const gdouble t   = Tm * 1.0e-4;
  const gdouble res = F * G * a * pow (t, b) / (1.0 + c * pow (t, d));
  NCM_UNUSED (cosmo);
  return res;
}

/**
 * nc_recomb_seager_pequignot_HI_case_B_dTm:
 * @recomb_seager: a #NcRecombSeager
 * @cosmo: a #NcHICosmo
 * @Tm: the matter (baryons) temperature $T_m$
 *
 * The case B $\HyII$ recombination coefficient derivative with respect to $T_m$.
 *
 * The derivative of the fitting formula of the case B recombination coefficient for $\HyII$
 * nc_recomb_seager_pequignot_HI_case_B ().
 *
 * Returns: the value of the case B recombination coefficient for $\HyII$, $d\alpha_H/dT_m$.
 */
gdouble
nc_recomb_seager_pequignot_HI_case_B_dTm (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, gdouble Tm)
{
  NcRecombSeagerPrivate * const self = recomb_seager->priv;
  const gdouble F   = self->H_fudge; /* fudge factor */
  const gdouble G   = 1e-19;
  const gdouble a   = 4.309;
  const gdouble b   = -0.6166;
  const gdouble c   = 0.6703;
  const gdouble d   = 0.5300;
  const gdouble t   = Tm * 1e-4;
  const gdouble t_b = pow (t, b);
  const gdouble t_d = pow (t, d);
  const gdouble res = a * F * G * (b + c * (b - d) * t_d) * t_b / (Tm * gsl_pow_2 (1.0 + c * t_d));
  NCM_UNUSED (cosmo);
  return res;
}

/**
 * nc_recomb_seager_hummer_HeI_case_B:
 * @recomb_seager: a #NcRecombSeager
 * @cosmo: a #NcHICosmo
 * @Tm: the matter (baryons) temperature $T_m$
 *
 * The case B $\HeII$ recombination coefficient.
 *
 * The fitting formula of the case B recombination coefficient for $\HeII$ as
 * in [Hummer (1998)][XHummer1998].
 *
 * Returns: the value of the case B recombination coefficient for $\HeII$, $\alpha_H$ .
 */
gdouble
nc_recomb_seager_hummer_HeI_case_B (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, gdouble Tm)
{
  const gdouble sqrt_Tm    = sqrt (Tm);
  const gdouble T1         = NC_RECOMB_SEAGER_HUMMER_HEI_CASE_B_T1;
  const gdouble T2         = NC_RECOMB_SEAGER_HUMMER_HEI_CASE_B_T2; /* original Seager paper: 3.0 */
  const gdouble p          = NC_RECOMB_SEAGER_HUMMER_HEI_CASE_B_P;
  const gdouble q          = NC_RECOMB_SEAGER_HUMMER_HEI_CASE_B_Q;
  const gdouble sqrt_T1    = sqrt (T1);
  const gdouble sqrt_T2    = sqrt (T2);
  const gdouble sqrt_Tm_T2 = sqrt_Tm / sqrt_T2;
  const gdouble sqrt_Tm_T1 = sqrt_Tm / sqrt_T1;
  const gdouble res        = q / (sqrt_Tm_T2 * pow (1.0 + sqrt_Tm_T2, 1.0 - p) * pow (1.0 + sqrt_Tm_T1, 1.0 + p));
  
  NCM_UNUSED (cosmo);

  return res;
}

/**
 * nc_recomb_seager_hummer_HeI_case_B_dTm:
 * @recomb_seager: a #NcRecombSeager
 * @cosmo: a #NcHICosmo
 * @Tm: the matter (baryons) temperature $T_m$
 *
 * The case B $\HeII$ recombination coefficient derivative with respect to Tm.
 *
 * The derivative of the fitting formula of the case B recombination coefficient for $\HeII$
 * nc_recomb_seager_hummer_HeI_case_B ().
 *
 * Returns: the value of the case B recombination coefficient for $\HeII$, $d\alpha_H/dT_m$.
 */
gdouble
nc_recomb_seager_hummer_HeI_case_B_dTm (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, gdouble Tm)
{
  const gdouble sqrt_Tm    = sqrt (Tm);
  const gdouble T1         = NC_RECOMB_SEAGER_HUMMER_HEI_CASE_B_T1;
  const gdouble T2         = NC_RECOMB_SEAGER_HUMMER_HEI_CASE_B_T2; /* original Seager paper: 3.0 */
  const gdouble p          = NC_RECOMB_SEAGER_HUMMER_HEI_CASE_B_P;
  const gdouble q          = NC_RECOMB_SEAGER_HUMMER_HEI_CASE_B_Q;
  const gdouble sqrt_T1    = sqrt (T1);
  const gdouble sqrt_T2    = sqrt (T2);
  const gdouble T1_2       = T1 * T1;
  const gdouble sqrt_Tm_T2 = sqrt_Tm / sqrt_T2;
  const gdouble sqrt_Tm_T1 = sqrt_Tm / sqrt_T1;
  const gdouble Tm_T1_3_2  = gsl_pow_3 (sqrt_Tm_T1);
  const gdouble res        = -q *
    (Tm * (2.0 + p + 3.0 * sqrt_Tm_T2) +
     T1 * sqrt_Tm_T1 * (1.0 + (2.0 - p) * sqrt_Tm_T2)) /
    (2.0 * T1_2 * Tm_T1_3_2 * sqrt_Tm_T2 *
     pow (1.0 + sqrt_Tm_T2, 2.0 - p) *
     pow (1.0 + sqrt_Tm_T1, 2.0 + p));
  
  NCM_UNUSED (cosmo);

  return res;
}

/**
 * nc_recomb_seager_hummer_HeI_case_B_trip:
 * @recomb_seager: a #NcRecombSeager
 * @cosmo: a #NcHICosmo
 * @Tm: the matter (baryons) temperature $T_m$
 *
 * The case B via triplets $\HeII$ recombination coefficient.
 *
 * The fitting formula of the case B via triplets recombination coefficient for $\HeII$ as
 * in [Hummer (1998)][XHummer1998].
 *
 * Returns: the value of the case B via triplets recombination coefficient for $\HeII$, $\alpha_H$ .
 */
gdouble
nc_recomb_seager_hummer_HeI_case_B_trip (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, gdouble Tm)
{
  const gdouble sqrt_Tm    = sqrt (Tm);
  const gdouble T1         = NC_RECOMB_SEAGER_HUMMER_HEI_CASE_B_T1;
  const gdouble T2         = NC_RECOMB_SEAGER_HUMMER_HEI_CASE_B_T2; /* original Seager paper: 3.0 */
  const gdouble p          = NC_RECOMB_SEAGER_HUMMER_HEI_CASE_B_P_TRIP;
  const gdouble q          = NC_RECOMB_SEAGER_HUMMER_HEI_CASE_B_Q_TRIP;
  const gdouble sqrt_T1    = sqrt (T1);
  const gdouble sqrt_T2    = sqrt (T2);
  const gdouble sqrt_Tm_T2 = sqrt_Tm / sqrt_T2;
  const gdouble sqrt_Tm_T1 = sqrt_Tm / sqrt_T1;
  const gdouble res        = q / (sqrt_Tm_T2 * pow (1.0 + sqrt_Tm_T2, 1.0 - p) * pow (1.0 + sqrt_Tm_T1, 1.0 + p));
  
  NCM_UNUSED (cosmo);

  return res;
}

/**
 * nc_recomb_seager_hummer_HeI_case_B_trip_dTm:
 * @recomb_seager: a #NcRecombSeager
 * @cosmo: a #NcHICosmo
 * @Tm: the matter (baryons) temperature $T_m$
 *
 * The case B via triplets $\HeII$ recombination coefficient derivative with respect to Tm.
 *
 * The derivative of the fitting formula of the case B via triplets recombination coefficient for $\HeII$
 * nc_recomb_seager_hummer_HeI_case_B_trip().
 *
 * Returns: the value of the case B via triplets recombination coefficient for $\HeII$, $d\alpha_H/dT_m$.
 */
gdouble
nc_recomb_seager_hummer_HeI_case_B_trip_dTm (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, gdouble Tm)
{
  const gdouble sqrt_Tm    = sqrt (Tm);
  const gdouble T1         = NC_RECOMB_SEAGER_HUMMER_HEI_CASE_B_T1;
  const gdouble T2         = NC_RECOMB_SEAGER_HUMMER_HEI_CASE_B_T2; /* original Seager paper: 3.0 */
  const gdouble p          = NC_RECOMB_SEAGER_HUMMER_HEI_CASE_B_P_TRIP;
  const gdouble q          = NC_RECOMB_SEAGER_HUMMER_HEI_CASE_B_Q_TRIP;
  const gdouble sqrt_T1    = sqrt (T1);
  const gdouble sqrt_T2    = sqrt (T2);
  const gdouble T1_2       = T1 * T1;
  const gdouble sqrt_Tm_T2 = sqrt_Tm / sqrt_T2;
  const gdouble sqrt_Tm_T1 = sqrt_Tm / sqrt_T1;
  const gdouble Tm_T1_3_2  = gsl_pow_3 (sqrt_Tm_T1);
  const gdouble res        = -q *
    (Tm * (2.0 + p + 3.0 * sqrt_Tm_T2) +
     T1 * sqrt_Tm_T1 * (1.0 + (2.0 - p) * sqrt_Tm_T2)) /
    (2.0 * T1_2 * Tm_T1_3_2 * sqrt_Tm_T2 *
     pow (1.0 + sqrt_Tm_T2, 2.0 - p) *
     pow (1.0 + sqrt_Tm_T1, 2.0 + p));
  
  NCM_UNUSED (cosmo);

  return res;
}

/**
 * nc_recomb_seager_weinberg_HII_ion_rate:
 * @recomb_seager: a #NcRecombSeager
 * @cosmo: a #NcHICosmo
 * @XHII: FIXME
 * @Tm: FIXME
 * @XHeII: FIXME
 * @x: normalized scale factor inverse $x = 1 + z = a_0/a$
 *
 * $dX_\e/dx$ implemented using Weinbergs book
 *
 * Returns: FIXME
 */
gdouble
nc_recomb_seager_weinberg_HII_ion_rate (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, gdouble XHII, gdouble Tm, gdouble XHeII, gdouble x)
{
  const gdouble x2    = x * x;
  const gdouble x3    = x2 * x;
  const gdouble Xe    = XHII + XHeII;
  const gdouble alpha = nc_recomb_seager_pequignot_HI_case_B (recomb_seager, cosmo, Tm);
  const gdouble n_H0  = nc_hicosmo_H_number_density (cosmo);
  const gdouble n_H   = n_H0 * x3;
  const gdouble H     = nc_hicosmo_H (cosmo, x - 1.0) / ncm_c_kpc ();
  const gdouble f1    = alpha * n_H / (H * x);

  const gdouble Tm3_2 = sqrt (gsl_pow_3 (Tm));
  const gdouble f2na  = ncm_c_decay_H_rate_2s_1s () * (1.0 - XHII);
  const gdouble f2nb  = H / (n_H * ncm_c_HI_Lyman_wl3_8pi_2p_2Pmean ());
  const gdouble f2n   = f2na + f2nb;
  const gdouble f2da  = ncm_c_boltzmann_factor_HI_2s_2S0_5 (Tm) * Tm3_2 * alpha * (1.0 - XHII);
  const gdouble f2d   = f2na + f2nb + f2da;
  const gdouble f2    = f2n / f2d;

  const gdouble S     = ncm_c_boltzmann_factor_HI_1s_2S0_5 (Tm) * Tm3_2 / n_H;
  const gdouble f3    = (XHII * Xe - (1.0 - XHII) * S);

  return f1 * f2 * f3;
}
