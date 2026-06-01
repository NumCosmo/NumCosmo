/***************************************************************************
 *            nc_transfer_func_eh_no_baryon.c
 *
 *  Mon Nov 03 17:46:27 2025
 *  Copyright  2025  Mariana Penna-Lima <pennalima@unb.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna-Lima 2025 <pennalima@unb.br>
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
 * NcTransferFuncEHNoBaryon:
 *
 * Eisenstein-Hu fitting function for the transfer function with no baryons.
 *
 * This objects implements the Eisenstein-Hu fitting function for the transfer function.
 * See [Eisenstein and Hu (1998)][XEisenstein1998]
 * [[arXiv](https://arxiv.org/abs/astro-ph/9709112)] for more details.
 *
 * FIXME (update this doc)
 * The transfer function is divided into a sum of two picies,
 * \begin{equation*}
 *  T(k) = \frac{\Omega_b}{\Omega_m} T_b(k) + \frac{\Omega_{c}}{\Omega_m} T_c(k) \, ,
 * \end{equation*}
 * where $\Omega_b$, $\Omega_m$ and $\Omega_{c}$ are the baryons (#nc_hicosmo_Omega_b0),
 * matter (#nc_hicosmo_Omega_m0) and cold dark matter (#nc_hicosmo_Omega_c0) density parameters today, respectevely.
 * $T_b(k)$ and $T_c(k)$ are the weights from baryons and cold dark matter to the transfer function $T(k)$.
 *
 * The cold dark matter term is defined as,
 * \begin{equation*}
 *  T_c(k) = f \, \widetilde{T}_0(k, 1, \beta_c) + (1-f) \, \widetilde{T}_0(k, \alpha_c, \beta_c)
 * \end{equation*}
 * and
 * $$ f = \frac{1}{1+ (ks/5.4)^4} \, ,$$
 * with
 * $$\widetilde{T}_0(k, \alpha_c, \beta_c) = \frac{\ln\left( e + 1.8 \beta_c \, q  \right)}{\ln\left( e + 1.8 \beta_c \, q  \right) + C q^2}  $$
 * and
 * $$ C = \frac{14.2}{\alpha_c} + \frac{386}{1+ 69.9 \, q^{1.08}} \, .$$
 * The paremeter $q$, is defined as,
 * $$ q = \frac{k}{13.41 k_{eq}} \, .$$
 * $\alpha_c$ and $\beta_c$ are fit by,
 * $$\alpha_c = a_1^{-\Omega_b/\Omega_m} \, a_2^{-(\Omega_b/\Omega_m)^3} \, ,$$
 * $$ a_1 = (46.9 \, \Omega_m h^2)^{0.670} \left[ 1 + (32.1 \, \Omega_m h^2)^{-0.532}  \right] \, ,$$
 * $$ a_2 = (12.0 \, \Omega_m h^2)^{0.424} \left[ 1 + (45.0 \, \Omega_m h^2)^{-0.582}  \right] \, ,$$
 * $$ \beta_c^{-1} = 1 + b_1 \left[ (\Omega_c/\Omega_m)^{b_{2}}  -1\right] \, ,$$
 * $$ b_1 = 0.944 \left[ 1+ (458\, \Omega_m h^2)^{-0.708}  \right]^{-1} \, , $$
 * $$ b_2 = (0.395 \, \Omega_m h^2)^{-0.0266}  \, .$$
 *
 * The baryon term is defined as,
 * \begin{equation*}
 *   T_b(k) = \left[ \frac{\widetilde{T}_0(k, 1, 1)}{1 + (ks/5.2)^2} + \frac{\alpha_b}{1+ (\beta_b/ks)^3}\mathrm{e}^{-(k/k_{\mathrm{Silk}})^{1.4}}  \right] \, \frac{\sin (k \tilde{s})}{k \tilde{s}} \, .
 * \end{equation*}
 * where, $s$ is the sound horizon scale, given by,
 * $$ s = \int_0^{z_d} c_s (1+z) \mathrm{d}t  = \frac{2}{3k_{eq}}\sqrt{\frac{6}{R(z_{eq})}} \ln \left( \frac{\sqrt{1+R(z_d)} + \sqrt{R(z_{d}) + R(z_{eq})}}{1 + \sqrt{R(z_{eq})}} \right) \, ,$$
 * with $z_d$ as the drag redshift (#nc_distance_drag_redshift) and $c_s$ is the baryon-photon plasma speed of sound (see #nc_hicosmo_bgp_cs2). The quantity $R$ is defined as,
 * $$ R(z) \equiv 31.5 \, \Omega_b h^2 \, \left( \frac{T_{\mathrm{CMB}}}{2.7} \right)^{-4} \left( \frac{z}{10^3}  \right)^{-1} \, .$$
 * The Silk damping scale is well fit by the approximation,
 * $$ k_{Silk} = 1.6 (\Omega_b h^2)^{0.52}(\Omega_m h^2)^{0.73} \left[ 1 + (10.4 \, \Omega_m h^2)^{-0.95}  \right] \, \mathrm{Mpc}^{-1} \, ,$$
 * and
 * $$ \alpha_b = 2.07 \,  k_{eq} \, s (1+ R_d )^{-3/4} G \left( \frac{1+ z_{eq}}{1+z_d}  \right) \, ,$$
 * $$ G(y) = y \left[-6 \sqrt{1+y} + (2+3y) \ln \left( \frac{\sqrt{1+y} + 1}{\sqrt{1+y} - 1} \right) \right]  \, ,$$
 * $$ k_{eq} = (2 \Omega_m H_0^2 z_{eq})^{1/2}  \, ,$$
 * $$ z_{eq} = 2.5 \times 10^4 \Omega_m h^2 \, (T_{\mathrm{CMB}} / 2.7)^{-4} \, ,$$
 * $$ R(z_d) = 31.5 \, \Omega_b h^2 (T_{cmb} / 2.7)^{-4} (z_d/10^3)^{-1} \, ,$$
 * $$ \tilde{s}(k) = \frac{s}{\left[ 1 + (\beta_{node}/ks)^3 \right]^{1/3}} \, ,$$
 * $$\beta_{node} = 8.41 \, (\Omega_m h^2)^{0.435} \, , $$
 * an finally the last parameter of the fitting function,
 * $$\beta_b = 0.5 + \frac{\Omega_b}{\Omega_m} + \left( 3  -2 \frac{\Omega_b}{\Omega_m} \right) \sqrt{(17.2 \, \Omega_m h^2)^2 +1} \, .$$
 *
 * With those relations in hand it is possible to evaluate the transfer function from [Eisenstein and Hu (1998)][XEisenstein1998] fitting formula.
 */


#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_transfer_func_eh_no_baryon.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_bessel.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcTransferFuncEHNoBaryonPrivate
{
  gdouble c2_2;
  gdouble h;
  gdouble wm;
  gdouble s;
  gdouble alphaGamma;
} NcTransferFuncEHNoBaryonPrivate;

struct _NcTransferFuncEHNoBaryon
{
  /*< private >*/
  NcTransferFunc parent_instance;
};

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcTransferFuncEHNoBaryon, nc_transfer_func_eh_no_baryon, NC_TYPE_TRANSFER_FUNC)

static void
nc_transfer_func_eh_no_baryon_init (NcTransferFuncEHNoBaryon *tf_eh)
{
  NcTransferFuncEHNoBaryonPrivate * const self = nc_transfer_func_eh_no_baryon_get_instance_private (tf_eh);

  self->c2_2       = 0.0;
  self->h          = 0.0;
  self->wm         = 0.0;
  self->s          = 0.0;
  self->alphaGamma = 0.0;
}

static void
_nc_transfer_func_eh_no_baryon_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_transfer_func_eh_no_baryon_parent_class)->finalize (object);
}

static void _nc_transfer_func_eh_no_baryon_prepare (NcTransferFunc *tf, NcHICosmo *cosmo);
static gdouble _nc_transfer_func_eh_no_baryon_calc (NcTransferFunc *tf, gdouble kh);

static void
nc_transfer_func_eh_no_baryon_class_init (NcTransferFuncEHNoBaryonClass *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcTransferFuncClass *parent_class = NC_TRANSFER_FUNC_CLASS (klass);

  object_class->finalize = &_nc_transfer_func_eh_no_baryon_finalize;

  parent_class->prepare = &_nc_transfer_func_eh_no_baryon_prepare;
  parent_class->calc    = &_nc_transfer_func_eh_no_baryon_calc;
}

static void
_nc_transfer_func_eh_no_baryon_prepare (NcTransferFunc *tf, NcHICosmo *cosmo)
{
  NcTransferFuncEHNoBaryon *tf_eh              = NC_TRANSFER_FUNC_EH_NO_BARYON (tf);
  NcTransferFuncEHNoBaryonPrivate * const self = nc_transfer_func_eh_no_baryon_get_instance_private (tf_eh);

  const gdouble T_0    = nc_hicosmo_T_gamma0 (cosmo);
  const gdouble c2     = T_0 / 2.7; /* \theta = 2.725/2.7 where 2.725 is the CMB temperature */
  const gdouble h      = nc_hicosmo_h (cosmo);
  const gdouble h2     = h * h;
  const gdouble wb     = nc_hicosmo_Omega_b0 (cosmo) * h2;
  const gdouble wm     = nc_hicosmo_Omega_m0 (cosmo) * h2;
  const gdouble wb_wm  = wb / wm;
  const gdouble wb_wm2 = wb_wm * wb_wm;
  const gdouble wb_3_4 = pow (wb, 0.75);

  const gdouble c2_2       = c2 * c2;
  const gdouble s          = 44.5 * log (9.83 / wm) / sqrt (1.0 + 10.0 * wb_3_4);
  const gdouble alphaGamma = 1.0 - 0.328 * log (431.0 * wm) * wb_wm + 0.38 * log (22.3 * wm) * wb_wm2;

  self->c2_2       = c2_2;
  self->h          = h;
  self->wm         = wm;
  self->s          = s;
  self->alphaGamma = alphaGamma;

/*  printf("s = %g\n, alphaGamma = %g\n", self->s, self->alphaGamma); */
}

static gdouble
_nc_transfer_func_eh_no_baryon_calc (NcTransferFunc *tf, gdouble kh)
{
  NcTransferFuncEHNoBaryon *tf_eh              = NC_TRANSFER_FUNC_EH_NO_BARYON (tf);
  NcTransferFuncEHNoBaryonPrivate * const self = nc_transfer_func_eh_no_baryon_get_instance_private (tf_eh);

  const gdouble k         = kh * self->h; /* [Mpc^-1] */
  const gdouble Gamma_eff = self->wm * (self->alphaGamma + (1.0 - self->alphaGamma) / (1.0 + pow (0.43 * k * self->s, 4)));
  const gdouble q         = k * self->c2_2 / Gamma_eff;
  const gdouble q2        = q * q;
  const gdouble Lo        = log (2.0 * M_E + 1.8 * q);
  const gdouble Co        = 14.2 + 731.0 / (1.0 + 62.5 * q);
  const gdouble To        = Lo / (Lo + Co * q2);

  return To;
}

/**
 * nc_transfer_func_eh_no_baryon_new:
 *
 * Creates a new #NcTransferFunc of the #NcTransferFuncEHNoBaryon type.
 *
 * Returns: A new #NcTransferFunc.
 */
NcTransferFunc *
nc_transfer_func_eh_no_baryon_new ()
{
  return g_object_new (NC_TYPE_TRANSFER_FUNC_EH_NO_BARYON, NULL);
}

