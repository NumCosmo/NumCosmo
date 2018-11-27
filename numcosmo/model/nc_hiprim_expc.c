/***************************************************************************
 *            nc_hiprim_expc.c
 *
 *  Sun May 08 18:13:19 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hiprim_expc.c
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_hiprim_expc
 * @title: NcHIPrimExpc
 * @short_description: Exponential cutoff modification of the power law primordial spectrum
 * 
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hiprim_expc.h"

G_DEFINE_TYPE (NcHIPrimExpc, nc_hiprim_expc, NC_TYPE_HIPRIM);

enum {
  PROP_0,
  PROP_SIZE,
};

static void
nc_hiprim_expc_init (NcHIPrimExpc *nc_hiprim_expc)
{
}

static void
nc_hiprim_expc_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hiprim_expc_parent_class)->finalize (object);
}

static gdouble _nc_hiprim_expc_lnSA_powespec_lnk (NcHIPrim *prim, const gdouble lnk);
static gdouble _nc_hiprim_expc_lnT_powespec_lnk (NcHIPrim *prim, const gdouble lnk);

static void
nc_hiprim_expc_class_init (NcHIPrimExpcClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcHIPrimClass *prim_class  = NC_HIPRIM_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->finalize = nc_hiprim_expc_finalize;

  ncm_model_class_set_name_nick (model_class, "Exponential cut model for primordial spectra", "Expc");
  ncm_model_class_add_params (model_class, NC_HIPRIM_EXPC_SPARAM_LEN, 0, PROP_SIZE);

  /* Set ln10e10ASA param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_EXPC_LN10E10ASA, "\\log(10^{10}A_{\\mathrm{SA}})", "ln10e10ASA",
                              0.0, 5.0, 1.0e0,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_EXPC_DEFAULT_LN10E10ASA,
                              NCM_PARAM_TYPE_FIXED);

  /* Set N_SA param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_EXPC_N_SA, "n_{\\mathrm{SA}}", "n_SA",
                              0.5, 1.5, 1.0e-1,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_EXPC_DEFAULT_N_SA,
                              NCM_PARAM_TYPE_FIXED);

  /* Set lambdac param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_EXPC_LAMBDAC, "\\lambda_\\mathrm{c}", "lambdac",
                              0.0, 30.0, 1.0,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_EXPC_DEFAULT_LAMBDAC,
                              NCM_PARAM_TYPE_FIXED);
  /* Set lnkc param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_EXPC_LNKC, "\\ln(k_\\mathrm{c})", "lnkc",
                              -12.0, -3.0, 1.0e0,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_EXPC_DEFAULT_LNKC,
                              NCM_PARAM_TYPE_FIXED);
  /* Set c param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_EXPC_C, "c", "c",
                              0.0, 0.99, 2.0e-1,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_EXPC_DEFAULT_C,
                              NCM_PARAM_TYPE_FIXED);

  /* Set T_SA_ratio param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_EXPC_T_SA_RATIO, "A_T/A_{\\mathrm{SA}}", "T_SA_ratio",
                              0.0, 10.0, 1.0e-1,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_EXPC_DEFAULT_T_SA_RATIO,
                              NCM_PARAM_TYPE_FIXED);
  
  /* Set N_T param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_EXPC_N_T, "n_{\\mathrm{T}}", "n_T",
                              -0.5, 0.5, 1.0e-2,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_EXPC_DEFAULT_N_T,
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  nc_hiprim_set_lnSA_powspec_lnk_impl (prim_class, &_nc_hiprim_expc_lnSA_powespec_lnk);
  nc_hiprim_set_lnT_powspec_lnk_impl (prim_class, &_nc_hiprim_expc_lnT_powespec_lnk);
}

/**
 * nc_hiprim_expc_new: (constructor)
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcHIPrimExpc *
nc_hiprim_expc_new (void)
{
  NcHIPrimExpc *prim_expc = g_object_new (NC_TYPE_HIPRIM_EXPC,
                                          NULL);
  return prim_expc;
}

#define VECTOR     (NCM_MODEL (prim)->params)
#define LN10E10ASA (ncm_vector_get (VECTOR, NC_HIPRIM_EXPC_LN10E10ASA))
#define N_SA       (ncm_vector_get (VECTOR, NC_HIPRIM_EXPC_N_SA))
#define LAMBDAC    (ncm_vector_get (VECTOR, NC_HIPRIM_EXPC_LAMBDAC))
#define LNKC       (ncm_vector_get (VECTOR, NC_HIPRIM_EXPC_LNKC))
#define C          (ncm_vector_get (VECTOR, NC_HIPRIM_EXPC_C))
#define T_SA_RATIO (ncm_vector_get (VECTOR, NC_HIPRIM_EXPC_T_SA_RATIO))
#define N_T        (ncm_vector_get (VECTOR, NC_HIPRIM_EXPC_N_T))

/****************************************************************************
 * Power spectrum
 ****************************************************************************/

static gdouble
_nc_hiprim_expc_lnSA_powespec_lnk (NcHIPrim *prim, const gdouble lnk)
{
  const gdouble ln_ka = lnk - prim->lnk_pivot;
  const gdouble lambda_c = LAMBDAC;
  const gdouble lnk_c    = LNKC;
  const gdouble c        = C;

  const gdouble k_kc_lambda_c = exp (lambda_c * (lnk - lnk_c)) - log1p (-c);

  gdouble ln_expc_fac;

  if (k_kc_lambda_c > 1.0)
    ln_expc_fac = log1p (- exp (- k_kc_lambda_c));
  else
    ln_expc_fac = log (- expm1 (- k_kc_lambda_c));

  return (N_SA - 1.0) * ln_ka + LN10E10ASA - 10.0 * M_LN10 + ln_expc_fac;
}

static gdouble
_nc_hiprim_expc_lnT_powespec_lnk (NcHIPrim *prim, const gdouble lnk)
{
  const gdouble ln_ka = lnk - prim->lnk_pivot;
  return N_T * ln_ka + LN10E10ASA - 10.0 * M_LN10 + log (T_SA_RATIO);
}
