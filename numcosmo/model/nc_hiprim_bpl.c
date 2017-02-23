/***************************************************************************
 *            nc_hiprim_bpl.c
 *
 *  Sun May 08 18:56:36 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hiprim_bpl.c
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
 * SECTION:nc_hiprim_bpl
 * @title: NcHIPrimBPL
 * @short_description: Broken power law modification of the power law primordial spectrum
 * 
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hiprim_bpl.h"

G_DEFINE_TYPE (NcHIPrimBPL, nc_hiprim_bpl, NC_TYPE_HIPRIM);

enum {
  PROP_0,
  PROP_SIZE,
};

static void
nc_hiprim_bpl_init (NcHIPrimBPL *nc_hiprim_bpl)
{
}

static void
nc_hiprim_bpl_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hiprim_bpl_parent_class)->finalize (object);
}

static gdouble _nc_hiprim_bpl_lnSA_powespec_lnk (NcHIPrim *prim, const gdouble lnk);

static void
nc_hiprim_bpl_class_init (NcHIPrimBPLClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcHIPrimClass *prim_class  = NC_HIPRIM_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->finalize = nc_hiprim_bpl_finalize;

  ncm_model_class_set_name_nick (model_class, "Broken power law model for primordial spectra", "BPL");
  ncm_model_class_add_params (model_class, NC_HIPRIM_BPL_SPARAM_LEN, 0, PROP_SIZE);

  /* Set ln10e10ASA param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_BPL_LN10E10ASA, "\\log(10^{10}A_{SA})", "ln10e10ASA",
                              0.0, 5.0, 1.0e0,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_BPL_DEFAULT_LN10E10ASA,
                              NCM_PARAM_TYPE_FIXED);

  /* Set N_SA param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_BPL_N_SA, "n_{SA}", "n_SA",
                              0.5, 1.5, 1.0e-1,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_BPL_DEFAULT_N_SA,
                              NCM_PARAM_TYPE_FIXED);

  /* Set lambdac param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_BPL_DELTA, "\\delta", "delta",
                              0.0, 8.0, 1.0,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_BPL_DEFAULT_DELTA,
                              NCM_PARAM_TYPE_FIXED);
  /* Set lnkc param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_BPL_LNKB, "\\ln(k_b)", "lnkb",
                              -12.0, -3.0, 1.0,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_BPL_DEFAULT_LNKB,
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  nc_hiprim_set_lnSA_powspec_lnk_impl (prim_class, &_nc_hiprim_bpl_lnSA_powespec_lnk);
}

/**
 * nc_hiprim_bpl_new: (constructor)
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcHIPrimBPL *
nc_hiprim_bpl_new (void)
{
  NcHIPrimBPL *prim_bpl = g_object_new (NC_TYPE_HIPRIM_BPL,
                                          NULL);
  return prim_bpl;
}

#define VECTOR     (NCM_MODEL (prim)->params)
#define LN10E10ASA (ncm_vector_get (VECTOR, NC_HIPRIM_BPL_LN10E10ASA))
#define N_SA       (ncm_vector_get (VECTOR, NC_HIPRIM_BPL_N_SA))
#define DELTA      (ncm_vector_get (VECTOR, NC_HIPRIM_BPL_DELTA))
#define LNKB       (ncm_vector_get (VECTOR, NC_HIPRIM_BPL_LNKB))

/****************************************************************************
 * Power spectrum
 ****************************************************************************/

static gdouble
_nc_hiprim_bpl_lnSA_powespec_lnk (NcHIPrim *prim, const gdouble lnk)
{
  const gdouble ln_ka = lnk - prim->lnk_pivot;
  const gdouble delta = DELTA;
  const gdouble lnk_b = LNKB;
  gdouble n_tot   = N_SA - 1.0;
  gdouble lnA_eff = LN10E10ASA - 10.0 * M_LN10;

  if (lnk < lnk_b)
  {
    n_tot   += delta;
    lnA_eff += delta * (prim->lnk_pivot - lnk_b);
  }

  return n_tot * ln_ka + lnA_eff;
}
