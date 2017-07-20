/***************************************************************************
 *            nc_hiprim_sbpl.c
 *
 *  Thu July 20 13:06:51 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hiprim_sbpl.c
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
 * SECTION:nc_hiprim_sbpl
 * @title: NcHIPrimSBPL
 * @short_description: Smooth Broken power law modification of the power law primordial spectrum
 * 
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hiprim_sbpl.h"

G_DEFINE_TYPE (NcHIPrimSBPL, nc_hiprim_sbpl, NC_TYPE_HIPRIM);

enum {
  PROP_0,
  PROP_SIZE,
};

static void
nc_hiprim_sbpl_init (NcHIPrimSBPL *nc_hiprim_sbpl)
{
}

static void
nc_hiprim_sbpl_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hiprim_sbpl_parent_class)->finalize (object);
}

static gdouble _nc_hiprim_sbpl_lnSA_powespec_lnk (NcHIPrim *prim, const gdouble lnk);
static gdouble _nc_hiprim_sbpl_lnT_powespec_lnk (NcHIPrim *prim, const gdouble lnk);

static void
nc_hiprim_sbpl_class_init (NcHIPrimSBPLClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcHIPrimClass *prim_class  = NC_HIPRIM_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->finalize = nc_hiprim_sbpl_finalize;

  ncm_model_class_set_name_nick (model_class, "Broken power law model for primordial spectra", "BPL");
  ncm_model_class_add_params (model_class, NC_HIPRIM_SBPL_SPARAM_LEN, 0, PROP_SIZE);

  /* Set ln10e10ASA param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_SBPL_LN10E10ASA, "\\log(10^{10}A_{SA})", "ln10e10ASA",
                              0.0, 5.0, 1.0e0,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_SBPL_DEFAULT_LN10E10ASA,
                              NCM_PARAM_TYPE_FIXED);
  /* Set N_SA param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_SBPL_N_SA, "n_{SA}", "n_SA",
                              0.5, 1.5, 1.0e-1,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_SBPL_DEFAULT_N_SA,
                              NCM_PARAM_TYPE_FIXED);
  /* Set lambdac param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_SBPL_DELTA, "\\delta", "delta",
                              0.0, 8.0, 1.0,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_SBPL_DEFAULT_DELTA,
                              NCM_PARAM_TYPE_FIXED);
  /* Set lambdac param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_SBPL_RA, "R_A", "RA",
                              0.0, 4.0, 1.0,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_SBPL_DEFAULT_RA,
                              NCM_PARAM_TYPE_FIXED);
  /* Set lnkc param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_SBPL_LNKB, "\\ln(k_b)", "lnkb",
                              -12.0, -3.0, 1.0,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_SBPL_DEFAULT_LNKB,
                              NCM_PARAM_TYPE_FIXED);
  /* Set lnkc param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_SBPL_LAMBDA, "\\lambda", "lambda",
                              1.0e-3, 30.0, 1.0,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_SBPL_DEFAULT_LAMBDA,
                              NCM_PARAM_TYPE_FIXED);
  /* Set T_SA_ratio param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_SBPL_T_SA_RATIO, "A_T/A_{SA}", "T_SA_ratio",
                              0.0, 10.0, 1.0e-1,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_SBPL_DEFAULT_T_SA_RATIO,
                              NCM_PARAM_TYPE_FIXED);
  /* Set N_T param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_SBPL_N_T, "n_{T}", "n_T",
                              -0.5, 0.5, 1.0e-2,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_SBPL_DEFAULT_N_T,
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  nc_hiprim_set_lnSA_powspec_lnk_impl (prim_class, &_nc_hiprim_sbpl_lnSA_powespec_lnk);
  nc_hiprim_set_lnT_powspec_lnk_impl  (prim_class, &_nc_hiprim_sbpl_lnT_powespec_lnk);
}

/**
 * nc_hiprim_sbpl_new: (constructor)
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcHIPrimSBPL *
nc_hiprim_sbpl_new (void)
{
  NcHIPrimSBPL *prim_bpl = g_object_new (NC_TYPE_HIPRIM_SBPL,
                                          NULL);
  return prim_bpl;
}

#define VECTOR     (NCM_MODEL (prim)->params)
#define LN10E10ASA (ncm_vector_get (VECTOR, NC_HIPRIM_SBPL_LN10E10ASA))
#define N_SA       (ncm_vector_get (VECTOR, NC_HIPRIM_SBPL_N_SA))
#define DELTA      (ncm_vector_get (VECTOR, NC_HIPRIM_SBPL_DELTA))
#define RA         (ncm_vector_get (VECTOR, NC_HIPRIM_SBPL_RA))
#define LNKB       (ncm_vector_get (VECTOR, NC_HIPRIM_SBPL_LNKB))
#define LAMBDA     (ncm_vector_get (VECTOR, NC_HIPRIM_SBPL_LAMBDA))
#define T_SA_RATIO (ncm_vector_get (VECTOR, NC_HIPRIM_SBPL_T_SA_RATIO))
#define N_T        (ncm_vector_get (VECTOR, NC_HIPRIM_SBPL_N_T))

/****************************************************************************
 * Power spectrum
 ****************************************************************************/

static gdouble
_nc_hiprim_sbpl_lnSA_powespec_lnk (NcHIPrim *prim, const gdouble lnk)
{
  const gdouble ln_ka   = lnk - prim->lnk_pivot;
  const gdouble delta   = DELTA;
  const gdouble ra      = RA;
  const gdouble lnk_b   = LNKB;
  const gdouble lambda  = LAMBDA;

  const gdouble lnA_eff = LN10E10ASA - 10.0 * M_LN10;

  const gdouble lnPplaw = (N_SA - 1.0) * ln_ka + lnA_eff;

  const gdouble lnX     = (lnk - lnk_b);

  gdouble one_1px_plambda, one_1px_mlambda;
    
  if (lnX > 0.0)
  {
    const gdouble x_mlambda = exp (-lambda * lnX);

    one_1px_plambda = x_mlambda / (1.0 + x_mlambda);
    one_1px_mlambda = 1.0 / (1.0 + x_mlambda);
  }
  else
  {
    const gdouble x_plambda = exp (+lambda * lnX);

    one_1px_plambda = 1.0 / (1.0 + x_plambda);
    one_1px_mlambda = x_plambda / (1.0 + x_plambda);
  }

  {
    const gdouble lnxi    = log (ra * exp (delta * ln_ka) * one_1px_plambda + one_1px_mlambda);
    return lnxi + lnPplaw;
  }
}

static gdouble
_nc_hiprim_sbpl_lnT_powespec_lnk (NcHIPrim *prim, const gdouble lnk)
{
  const gdouble ln_ka = lnk - prim->lnk_pivot;
  return N_T * ln_ka + LN10E10ASA - 10.0 * M_LN10 + log (T_SA_RATIO);
}

