/***************************************************************************
 *            nc_transfer_func_eh.c
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
 * SECTION:nc_transfer_func_eh
 * @title: NcTransferFuncEH
 * @short_description: Eisenstein-Hu fitting function for the transfer function.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_transfer_func_eh.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_bessel.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcTransferFuncEHPrivate
{
  gdouble h;
  gdouble s;
  gdouble keq_1341;                /* keq * 13.41 */
  gdouble ksilk;
  gdouble b_node3;
  gdouble ab, bc, bb, bb3, ac_142; /* ac_142 = 14.2/ac */
  gdouble wb_wm, wc_wm;
  gboolean CCL_comp;
};

enum
{
  PROP_0,
  PROP_CCL_COMP,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcTransferFuncEH, nc_transfer_func_eh, NC_TYPE_TRANSFER_FUNC);

static void
nc_transfer_func_eh_init (NcTransferFuncEH *tf_eh)
{
  NcTransferFuncEHPrivate * const self = tf_eh->priv = nc_transfer_func_eh_get_instance_private (tf_eh);
  
  self->h        = 0.0;
  self->s        = 0.0;
  self->keq_1341 = 0.0;
  self->ksilk    = 0.0;
  self->b_node3  = 0.0;
  self->ab       = 0.0;
  self->bc       = 0.0;
  self->bb       = 0.0;
  self->bb3      = 0.0;
  self->ac_142   = 0.0;
  self->wb_wm    = 0.0;
  self->wc_wm    = 0.0;
  self->CCL_comp = FALSE;
}

static void
_nc_transfer_func_eh_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcTransferFuncEH *tf_eh = NC_TRANSFER_FUNC_EH (object);
  g_return_if_fail (NC_IS_TRANSFER_FUNC_EH (object));

  switch (prop_id)
  {
    case PROP_CCL_COMP:
      nc_transfer_func_eh_set_CCL_comp (tf_eh, g_value_get_boolean (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_transfer_func_eh_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcTransferFuncEH *tf_eh = NC_TRANSFER_FUNC_EH (object);
  NcTransferFuncEHPrivate * const self = tf_eh->priv;
  g_return_if_fail (NC_IS_TRANSFER_FUNC_EH (object));

  switch (prop_id)
  {
    case PROP_CCL_COMP:
      g_value_set_boolean (value, self->CCL_comp);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_transfer_func_eh_finalize (GObject *object)
{
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_transfer_func_eh_parent_class)->finalize (object);
}

static void _nc_transfer_func_eh_prepare (NcTransferFunc *tf, NcHICosmo *cosmo);
static gdouble _nc_transfer_func_eh_calc (NcTransferFunc *tf, gdouble kh);

static void
nc_transfer_func_eh_class_init (NcTransferFuncEHClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcTransferFuncClass* parent_class = NC_TRANSFER_FUNC_CLASS (klass);

  object_class->set_property = &_nc_transfer_func_eh_set_property;
  object_class->get_property = &_nc_transfer_func_eh_get_property;
  object_class->finalize     = &_nc_transfer_func_eh_finalize;

  g_object_class_install_property (object_class,
                                   PROP_CCL_COMP,
                                   g_param_spec_boolean ("CCL-comp",
                                                         NULL,
                                                         "Whether to use CCL compatible mode",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  parent_class->prepare       = &_nc_transfer_func_eh_prepare;
  parent_class->calc          = &_nc_transfer_func_eh_calc;
}

static void
_nc_transfer_func_eh_prepare (NcTransferFunc *tf, NcHICosmo *cosmo)
{
  NcTransferFuncEH *tf_eh = NC_TRANSFER_FUNC_EH (tf);
  NcTransferFuncEHPrivate * const self = tf_eh->priv;
  
  const gdouble T_0   = nc_hicosmo_T_gamma0 (cosmo);
  const gdouble c1    = sqrt (6.0);
  const gdouble c2    = T_0 / 2.7;    /* \theta = 2.725/2.7 where 2.725 is the CMB temperature */
  const gdouble h     = nc_hicosmo_h (cosmo);
  const gdouble h2    = h * h;
  const gdouble wm    = nc_hicosmo_Omega_m0 (cosmo) * h2;
  const gdouble wb    = nc_hicosmo_Omega_b0 (cosmo) * h2;
  const gdouble wc    = nc_hicosmo_Omega_c0 (cosmo) * h2;
  const gdouble wb_wm = wb/wm;     /* \frac{\Omega_{b0}}{\Omega_{m0}} */
  const gdouble wc_wm = wc/wm;    /* \frac{\Omega_{c0}}{\Omega_{m0}} */

  const gdouble a1      = pow (46.9 * wm, 0.670) * (1.0 + pow (32.1 * wm, -0.532));
  const gdouble a2      = pow (12.0 * wm, 0.424) * (1.0 + pow (45.0 * wm, -0.582));
  const gdouble b1      = 0.313 * pow (wm, -0.419) * (1.0 + 0.607 * pow (wm, 0.674));
  const gdouble b2      = 0.238 * pow (wm, 0.223);
  const gdouble b3      = 0.944 / (1.0 + pow (458.0 * wm, -0.708));
  const gdouble b4      = pow (0.395 * wm, -0.0266);
  const gdouble wm_pow  = pow (wm, 0.435);
  const gdouble wm_pow3 = wm_pow * wm_pow * wm_pow;
  const gdouble b_node3 = 594.823321  * wm_pow3;      /* b_node3 =  {\beta_node}^3  8.41^3 = 594.823321  */

  const gdouble c2_2     = c2 * c2;
  const gdouble keq      = 7.46e-2 * wm / c2_2; /* unit: [Mpc^{-1}] */
  const gdouble keq_1341 = keq * 13.41;
  const gdouble ksilk    = 1.6 * pow (wb, 0.52) * pow (wm, 0.73) * (1.0 + pow (10.4 * wm, -0.95)); /* unit: [Mpc^{-1}] */

  const gdouble c2_4 = c2_2 * c2_2;
  const gdouble zeq  = 2.5e4 * wm / c2_4;
  const gdouble zd   = 1291.0 * (pow (wm, 0.251) / (1.0 + 0.659 * pow (wm, 0.828))) * (1.0 + b1 * pow (wb, b2));

  const gdouble c3  = 31.5e3 * wb / c2_4;
  const gdouble Rd  = self->CCL_comp ? c3 / (1.0 + zd) : c3 / zd; /* Rd is the ratio of the baryon to photon momentum density at redshift zd */
  const gdouble Req = 1.26 * wb_wm;                                /* Req = c3/zeq = 1.26 * wb/wm */

  const gdouble sqrt_Req = sqrt (Req);
  const gdouble s        = (2.0 * c1 / (3.0 * keq * sqrt_Req)) * log ((sqrt (1.0 + Rd) + sqrt (Rd + Req)) / (1.0 + sqrt_Req));  /* Sound horizon */
  const gdouble y        = self->CCL_comp ? (0.0 + zeq)/(1.0 + zd) : (1.0 + zeq)/(1.0 + zd);
  const gdouble y1       = sqrt(1.0 + y);
  const gdouble ab       = 2.07 * keq * s * pow (1.0 + Rd, -0.75) * y * (-6.0 * y1 + (2.0 + 3.0 * y) * log((y1 + 1.0)/(y1 - 1.0)));  /* \alpha_b */
  const gdouble bb       = 0.5 + wb_wm + (3.0 - 2.0 * wb_wm) * sqrt(1.0 + gsl_pow_2 (17.2 * wm));    /* \beta_b */
  const gdouble bb3      = gsl_pow_3(bb);
  const gdouble wb_wm3   = wb_wm * wb_wm * wb_wm;    /* {\frac{\Omega_{b0}}{\Omega_{m0}}}^3 */
  const gdouble ac       = pow (a1, - wb_wm) * pow (a2, -wb_wm3);  /* \alpha_c */
  const gdouble ac_142   = 14.2 / ac;
  const gdouble bc       = 1.0 / (1.0 + b3 * (pow (wc_wm, b4) - 1.0));  /* \beta_c */
  
  self->h        = h;
  self->s        = s;
  self->keq_1341 = keq_1341;
  self->ksilk    = ksilk;
  self->b_node3  = b_node3;
  self->ab       = ab;
  self->bc       = bc;
  self->bb       = bb;
  self->bb3      = bb3;
  self->ac_142   = ac_142;
  self->wb_wm    = wb_wm;
  self->wc_wm    = wc_wm;

//  printf("(14.2/ac) = %g, s = %g, bc = %g, (keq * 13.41) = %g\n", EH->ac_142, EH->s, EH->bc, EH->keq_1341);

}

static gdouble
_nc_transfer_func_eh_calc (NcTransferFunc *tf, gdouble kh)
{
  NcTransferFuncEH *tf_eh = NC_TRANSFER_FUNC_EH (tf);
  NcTransferFuncEHPrivate * const self = tf_eh->priv;

  const gdouble k   = kh * self->h;
  const gdouble ks  = k * self->s;
  const gdouble ks2 = ks * ks;
  const gdouble ks3 = ks2 * ks;
  const gdouble ks4 = ks3 * ks;

  const gdouble q        = k / (self->keq_1341);
  const gdouble q2       = q * q;
  const gdouble c4       = log (M_E + 1.8 * q);
  const gdouble k_ksilk  = pow (k / self->ksilk, 1.4);
  const gdouble c5       = exp (- k_ksilk);                                                   /* pow(M_E, -k_ksilk); */
  const gdouble s_tilda  = self->s / cbrt (1.0 + self->b_node3 / ks3);
  const gdouble ks_tilda = k * s_tilda;
  const gdouble jo       = gsl_sf_bessel_j0 (ks_tilda);                                     /* j0 = \sin(ks_tilda)/ks_tilda */
  const gdouble q_1_08   = pow (q, 1.08);
  const gdouble C        = 14.2 + 386.0 / (1.0 + 69.9 * q_1_08);                               /* C com ac = 1 */ /* <============================================= */
  const gdouble To       = c4 / (c4 + C * q2);
  const gdouble Tb       = (To / (1.0 + ks2 / 27.04) + (self->ab * c5) / (1.0 + self->bb3/ks3)) * jo; /* 27.04 = 5.2^2*/
  const gdouble f        = 1.0 / (1.0 + ks4 / 850.3056);                                         /* 850.3056 = 5.4^4 */
  const gdouble c6       = log (M_E + 1.8 * self->bc * q);
  const gdouble To1      = c6 / (c6 + C * q2);
  const gdouble C_ac     = self->ac_142 + 386.0/(1.0 + 69.9 * q_1_08);
  const gdouble To2      = c6 / (c6 + C_ac * q2);
  const gdouble Tc       = f * To1 + (1.0 - f) * To2;

  return self->wb_wm * Tb + self->wc_wm * Tc;
}

/**
 * nc_transfer_func_eh_new:
 *
 * FIXME
 *
 * Returns: A new #NcTransferFunc.
 */
NcTransferFunc *
nc_transfer_func_eh_new ()
{
  return g_object_new (NC_TYPE_TRANSFER_FUNC_EH, NULL);
}

/**
 * nc_transfer_func_eh_set_CCL_comp:
 * @tf_eh: a #NcTransferFuncEH
 * @CCL_comp: a boolean
 * 
 * (Un)Sets CCL compatibility mode.
 * 
 */
void 
nc_transfer_func_eh_set_CCL_comp (NcTransferFuncEH *tf_eh, gboolean CCL_comp)
{
  NcTransferFuncEHPrivate * const self = tf_eh->priv;
  self->CCL_comp = CCL_comp;
}

