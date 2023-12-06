/***************************************************************************
 *            nc_hiprim_power_law.c
 *
 *  Tue October 27 14:13:46 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hiprim_power_law.c
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
 * SECTION:nc_hiprim_power_law
 * @title: NcHIPrimPowerLaw
 * @short_description: Power law implementation for primordial spectra.
 *
 * Primordial adiabatic scalar power spectrum:
 * $$ \mathcal{P}_{SA}(k) = \mathcal{A}_\mathrm{s}\left(\frac{k}{k_\star}\right)^{n_s -1 }.$$
 *
 * Primordial tensor power spectrum:
 * $$ \mathcal{P}_T(k) = r \mathcal{A}_\mathrm{s} \left(\frac{k}{k_\star}\right)^{n_T -1 }.$$
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_hiprim_power_law.h"

G_DEFINE_TYPE (NcHIPrimPowerLaw, nc_hiprim_power_law, NC_TYPE_HIPRIM)

enum
{
  PROP_0,
  PROP_SIZE,
};

static void
nc_hiprim_power_law_init (NcHIPrimPowerLaw *nc_hiprim_power_law)
{
}

static void
_nc_hiprim_power_law_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hiprim_power_law_parent_class)->finalize (object);
}

static gdouble _nc_hiprim_power_law_lnSA_powespec_lnk (NcHIPrim *prim, const gdouble lnk);
static gdouble _nc_hiprim_power_law_lnT_powespec_lnk (NcHIPrim *prim, const gdouble lnk);

static void
nc_hiprim_power_law_class_init (NcHIPrimPowerLawClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcHIPrimClass *prim_class  = NC_HIPRIM_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->finalize = &_nc_hiprim_power_law_finalize;

  ncm_model_class_set_name_nick (model_class, "Power Law model for primordial spectra", "PowerLaw");
  ncm_model_class_add_params (model_class, NC_HIPRIM_POWER_LAW_SPARAM_LEN, 0, PROP_SIZE);

  /* Set ln10e10ASA param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_POWER_LAW_LN10E10ASA, "\\log(10^{10}A_{\\mathrm{SA}})", "ln10e10ASA",
                              2.0, 5.0, 1.0e0,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_POWER_LAW_DEFAULT_LN10E10ASA,
                              NCM_PARAM_TYPE_FIXED);

  /* Set T_SA_ratio param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_POWER_LAW_T_SA_RATIO, "A_T/A_{\\mathrm{SA}}", "T_SA_ratio",
                              0.0, 10.0, 1.0e-1,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_POWER_LAW_DEFAULT_T_SA_RATIO,
                              NCM_PARAM_TYPE_FIXED);

  /* Set N_SA param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_POWER_LAW_N_SA, "n_{\\mathrm{SA}}", "n_SA",
                              0.5, 1.5, 1.0e-1,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_POWER_LAW_DEFAULT_N_SA,
                              NCM_PARAM_TYPE_FIXED);

  /* Set N_T param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_POWER_LAW_N_T, "n_{\\mathrm{T}}", "n_T",
                              -0.5, 0.5, 1.0e-2,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_POWER_LAW_DEFAULT_N_T,
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  nc_hiprim_set_lnSA_powspec_lnk_impl (prim_class, &_nc_hiprim_power_law_lnSA_powespec_lnk);
  nc_hiprim_set_lnT_powspec_lnk_impl  (prim_class, &_nc_hiprim_power_law_lnT_powespec_lnk);
}

/**
 * nc_hiprim_power_law_new: (constructor)
 *
 * This function instantiates a new object of type #NcHIPrimPowerLaw.
 *
 * Returns: (transfer full): A new #NcHIPrimPowerLaw
 */
NcHIPrimPowerLaw *
nc_hiprim_power_law_new (void)
{
  NcHIPrimPowerLaw *prim_pl = g_object_new (NC_TYPE_HIPRIM_POWER_LAW,
                                            NULL);

  return prim_pl;
}

#define VECTOR (ncm_model_orig_params_peek_vector (NCM_MODEL (prim)))
#define LN10E10ASA (ncm_vector_get (VECTOR, NC_HIPRIM_POWER_LAW_LN10E10ASA))
#define T_SA_RATIO (ncm_vector_get (VECTOR, NC_HIPRIM_POWER_LAW_T_SA_RATIO))
#define N_SA       (ncm_vector_get (VECTOR, NC_HIPRIM_POWER_LAW_N_SA))
#define N_T        (ncm_vector_get (VECTOR, NC_HIPRIM_POWER_LAW_N_T))

/****************************************************************************
 * Power spectra
 ****************************************************************************/

static gdouble
_nc_hiprim_power_law_lnSA_powespec_lnk (NcHIPrim *prim, const gdouble lnk)
{
  const gdouble ln_ka = lnk - prim->lnk_pivot;

  return (N_SA - 1.0) * ln_ka + LN10E10ASA - 10.0 * M_LN10;
}

static gdouble
_nc_hiprim_power_law_lnT_powespec_lnk (NcHIPrim *prim, const gdouble lnk)
{
  const gdouble ln_ka = lnk - prim->lnk_pivot;

  return N_T * ln_ka + LN10E10ASA - 10.0 * M_LN10 + log (T_SA_RATIO);
}

