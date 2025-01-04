/***************************************************************************
 *            nc_hiprim_atan.c
 *
 *  Thu October 29 15:14:14 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hiprim_atan.c
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcHIPrimAtan:
 *
 * Arctangent modification of the power law primordial spectrum.
 *
 * This object implement the arctangent modification of the power law primordial
 * spectrum inspired in the quantum equilibrium models. See:
 * - [Valentini (2010)][XValentini2010]
 * - [Colin (2015)][XColin2015]
 * - [Underwood (2015)][XUnderwood2015]
 * - [Valentini (2015)][XValentini2015]
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_hiprim_atan.h"

G_DEFINE_TYPE (NcHIPrimAtan, nc_hiprim_atan, NC_TYPE_HIPRIM)

enum
{
  PROP_0,
  PROP_SIZE,
};

static void
nc_hiprim_atan_init (NcHIPrimAtan *prim_atan)
{
}

static void
nc_hiprim_atan_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hiprim_atan_parent_class)->finalize (object);
}

static gdouble _nc_hiprim_atan_lnSA_powespec_lnk (NcHIPrim *prim, const gdouble lnk);
static gdouble _nc_hiprim_atan_lnT_powespec_lnk (NcHIPrim *prim, const gdouble lnk);

static void
nc_hiprim_atan_class_init (NcHIPrimAtanClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcHIPrimClass *prim_class  = NC_HIPRIM_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->finalize = nc_hiprim_atan_finalize;

  ncm_model_class_set_name_nick (model_class, "Atan model for primordial spectra", "Atan");
  ncm_model_class_add_params (model_class, NC_HIPRIM_ATAN_SPARAM_LEN, 0, PROP_SIZE);

  /* Set ln10e10ASA param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_ATAN_LN10E10ASA, "\\log(10^{10}A_{\\mathrm{SA}})", "ln10e10ASA",
                              0.0, 5.0, 1.0,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_ATAN_DEFAULT_LN10E10ASA,
                              NCM_PARAM_TYPE_FIXED);

  /* Set N_SA param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_ATAN_N_SA, "n_{\\mathrm{SA}}", "n_SA",
                              0.5, 1.5, 1.0e-1,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_ATAN_DEFAULT_N_SA,
                              NCM_PARAM_TYPE_FIXED);

  /* Set lnkc param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_ATAN_LNKC, "\\ln(k_\\mathrm{c})", "lnkc",
                              -18.0, 1.0, 1.0e0,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_ATAN_DEFAULT_LNKC,
                              NCM_PARAM_TYPE_FIXED);
  /* Set c2 param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_ATAN_C2, "c_2", "c2",
                              0.0, 0.99, 2.0e-1,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_ATAN_DEFAULT_C2,
                              NCM_PARAM_TYPE_FIXED);
  /* Set c3 param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_ATAN_C3, "c_3", "c3",
                              0.5, 2.0, 3.0e-1,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_ATAN_DEFAULT_C3,
                              NCM_PARAM_TYPE_FIXED);
  /* Set lambda param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_ATAN_LAMBDA, "\\lambda", "lambda",
                              0.0, 60.0, 1.0,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_ATAN_DEFAULT_LAMBDA,
                              NCM_PARAM_TYPE_FIXED);

  /* Set T_SA_ratio param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_ATAN_T_SA_RATIO, "A_T/A_{\\mathrm{SA}}", "T_SA_ratio",
                              0.0, 10.0, 1.0e-1,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_ATAN_DEFAULT_T_SA_RATIO,
                              NCM_PARAM_TYPE_FIXED);

  /* Set N_T param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_ATAN_N_T, "n_{\\mathrm{T}}", "n_T",
                              -0.5, 0.5, 1.0e-2,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_ATAN_DEFAULT_N_T,
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  nc_hiprim_set_lnSA_powspec_lnk_impl (prim_class, &_nc_hiprim_atan_lnSA_powespec_lnk);
  nc_hiprim_set_lnT_powspec_lnk_impl  (prim_class, &_nc_hiprim_atan_lnT_powespec_lnk);
}

/**
 * nc_hiprim_atan_new: (constructor)
 *
 * This function instantiates a new object of type #NcHIPrimAtan.
 *
 * Returns: (transfer full): A new #NcHIPrimAtan
 */
NcHIPrimAtan *
nc_hiprim_atan_new (void)
{
  NcHIPrimAtan *prim_pl = g_object_new (NC_TYPE_HIPRIM_ATAN,
                                        NULL);

  return prim_pl;
}

#define VECTOR     (NCM_MODEL (prim))
#define LN10E10ASA (ncm_model_orig_param_get (VECTOR, NC_HIPRIM_ATAN_LN10E10ASA))
#define N_SA       (ncm_model_orig_param_get (VECTOR, NC_HIPRIM_ATAN_N_SA))
#define LNKC       (ncm_model_orig_param_get (VECTOR, NC_HIPRIM_ATAN_LNKC))
#define C2         (ncm_model_orig_param_get (VECTOR, NC_HIPRIM_ATAN_C2))
#define C3         (ncm_model_orig_param_get (VECTOR, NC_HIPRIM_ATAN_C3))
#define LAMBDA     (ncm_model_orig_param_get (VECTOR, NC_HIPRIM_ATAN_LAMBDA))
#define T_SA_RATIO (ncm_model_orig_param_get (VECTOR, NC_HIPRIM_ATAN_T_SA_RATIO))
#define N_T        (ncm_model_orig_param_get (VECTOR, NC_HIPRIM_ATAN_N_T))

/****************************************************************************
 * Power spectrum
 ****************************************************************************/

static gdouble
_Datan_xpd (const gdouble x, const gdouble d)
{
  const gdouble x2     = x * x;
  const gdouble onepx2 = 1.0 + x2;
  const gdouble Y      = d / onepx2;

  return Y * (1.0 + Y * (-x + Y * ((x2 - 1.0 / 3.0) + Y * x * (1.0 - x2))));
}

static gdouble
_nc_hiprim_atan_lnSA_powespec_lnk (NcHIPrim *prim, const gdouble lnk)
{
  const gdouble ln_ka     = lnk - prim->lnk_pivot;
  const gdouble ln_k_kc   = lnk - LNKC;
  const gdouble lambda    = LAMBDA;
  const gdouble k_kc_l    = exp (lambda * ln_k_kc);
  const gdouble c2        = C2;
  const gdouble c3        = C3;
  const gdouble pi_2      = ncm_c_pi () * 0.5;
  const gdouble pi_2_m_c3 = pi_2 - c3;
  const gdouble a1        = pi_2_m_c3 + c3 * c2;
  const gdouble tan_a1    = tan (a1);

  const gdouble atan_fac    = (k_kc_l > 1.0e-4) ? atan2 (k_kc_l + tan_a1, 1.0) - pi_2_m_c3 : _Datan_xpd (tan_a1, k_kc_l) + c3 * c2;
  const gdouble ln_atan_fac = log (atan_fac / c3);

  return (N_SA - 1.0) * ln_ka + LN10E10ASA - 10.0 * M_LN10 + ln_atan_fac;
}

static gdouble
_nc_hiprim_atan_lnT_powespec_lnk (NcHIPrim *prim, const gdouble lnk)
{
  const gdouble ln_ka = lnk - prim->lnk_pivot;

  return N_T * ln_ka + LN10E10ASA - 10.0 * M_LN10 + log (T_SA_RATIO);
}

