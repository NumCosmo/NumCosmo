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
 * This object implements the Eisenstein-Hu (1998) zero-baryon ("no-wiggle")
 * fitting function, which smooths over the baryon acoustic oscillations while
 * retaining the broadband shape:
 * \begin{equation*}
 *  T_0(k) = \frac{L_0}{L_0 + C_0\,q^2}, \qquad q = \frac{k\,\theta^2}{\Gamma_\mathrm{eff}},
 * \end{equation*}
 * where $\Gamma_\mathrm{eff}$ is an effective shape parameter built from the
 * matter (#nc_hicosmo_Omega_m0) and baryon (#nc_hicosmo_Omega_b0) densities and
 * $\theta = T_\mathrm{CMB}/2.7$.
 *
 * For the full definitions of $L_0$, $C_0$, and $\Gamma_\mathrm{eff}$, see the
 * theoretical background page (Zero-Baryon Variant section):
 * <a href="../../theory/transfer_func_eh.html">Eisenstein–Hu Transfer Function</a>.
 * Reference: [Eisenstein and Hu (1998)](https://arxiv.org/abs/astro-ph/9709112).
 */


#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/powspec/nc_transfer_func_eh_no_baryon.h"

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

