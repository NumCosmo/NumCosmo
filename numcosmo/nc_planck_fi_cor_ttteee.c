/***************************************************************************
 *            nc_planck_fi_cor_ttteee.c
 *
 *  Fri April 22 14:47:22 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_planck_fi_cor_ttteee.c
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
 * SECTION:nc_planck_fi_cor_ttteee
 * @title: NcPlanckFICorTTTEEE
 * @short_description: Planck Foreground and Instrument model for TT correlation maps
 *
 * FIXME
 *
 * If you use this object, cite [Planck 2015 results XI (2015)][XPlanckCollaboration2015a]
 * and related papers.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_planck_fi_cor_ttteee.h"

enum
{
  PROP_0,
  PROP_SIZE
};

G_DEFINE_TYPE (NcPlanckFICorTTTEEE, nc_planck_fi_cor_ttteee, NC_TYPE_PLANCK_FI_COR_TT);

static void
nc_planck_fi_cor_ttteee_init (NcPlanckFICorTTTEEE *cor_ttteee)
{
}

static void
nc_planck_fi_cor_ttteee_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_planck_fi_cor_ttteee_parent_class)->finalize (object);
}

static void
nc_planck_fi_cor_ttteee_class_init (NcPlanckFICorTTTEEEClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->finalize    = nc_planck_fi_cor_ttteee_finalize;

  ncm_model_class_set_name_nick (model_class, "Planck Foreground and Instument Model -- TT, TE, EE", "PlanckFICorTTTEEE");
  ncm_model_class_add_params (model_class, NC_PLANCK_FI_COR_TTTEEE_SPARAM_LEN - NC_PLANCK_FI_COR_TT_SPARAM_LEN, 0, PROP_SIZE);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_100, "A^{\\mathrm{dust}EE}_{100}", "galf_EE_A_100",
                              0.0, 10.0, 1.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_EE_A_100,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_100_143, "A^{\\mathrm{dust}EE}_{100 \\times 143}", "galf_EE_A_100_143",
                              0.0, 10.0, 1.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_EE_A_100_143,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_100_217, "A^{\\mathrm{dust}EE}_{100 \\times 217}", "galf_EE_A_100_217",
                              0.0, 10.0, 1.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_EE_A_100_217,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_143, "A^{\\mathrm{dust}EE}_{143}", "galf_EE_A_143",
                              0.0, 10.0, 1.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_EE_A_143,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_143_217, "A^{\\mathrm{dust}EE}_{143 \\times 217}", "galf_EE_A_143_217",
                              0.0, 10.0, 1.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_EE_A_143_217,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_217, "A^{\\mathrm{dust}EE}_{217}", "galf_EE_A_217",
                              0.0, 10.0, 1.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_EE_A_217,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_galf_EE_index, "n^{\\mathrm{dust}EE}", "galf_EE_index",
                              -10.0, 10.0, 0.1,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_EE_index,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_100, "A^{\\mathrm{dust}TE}_{100}", "galf_TE_A_100",
                              0.0, 10.0, 1.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_TE_A_100,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_100_143, "A^{\\mathrm{dust}TE}_{100 \\times 143}", "galf_TE_A_100_143",
                              0.0, 10.0, 1.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_TE_A_100_143,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_100_217, "A^{\\mathrm{dust}TE}_{100 \\times 217}", "galf_TE_A_100_217",
                              0.0, 10.0, 1.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_TE_A_100_217,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_143, "A^{\\mathrm{dust}TE}_{143}", "galf_TE_A_143",
                              0.0, 10.0, 1.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_TE_A_143,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_143_217, "A^{\\mathrm{dust}TE}_{143 \\times 217}", "galf_TE_A_143_217",
                              0.0, 10.0, 1.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_TE_A_143_217,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_217, "A^{\\mathrm{dust}TE}_{217}", "galf_TE_A_217",
                              0.0, 10.0, 1.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_TE_A_217,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_galf_TE_index, "n^{\\mathrm{dust}TE}", "galf_TE_index",
                              -10.0, 10.0, 0.1,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_TE_index,
                              NCM_PARAM_TYPE_FIXED);
  /*******************************************************************************************************************************************************/
  /* GAL DUST EPSILON TE 100 */
  /*******************************************************************************************************************************************************/
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_0T_0E, "\\epsilon^{\\mathrm{bleak}TE}_{0, 100}", "bleak_epsilon_0_0T_0E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_0T_0E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_0T_0E, "\\epsilon^{\\mathrm{bleak}TE}_{1, 100}", "bleak_epsilon_1_0T_0E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_0T_0E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_0T_0E, "\\epsilon^{\\mathrm{bleak}TE}_{2, 100}", "bleak_epsilon_2_0T_0E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_0T_0E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_0T_0E, "\\epsilon^{\\mathrm{bleak}TE}_{3, 100}", "bleak_epsilon_3_0T_0E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_0T_0E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_0T_0E, "\\epsilon^{\\mathrm{bleak}TE}_{4, 100}", "bleak_epsilon_4_0T_0E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_0T_0E,
                              NCM_PARAM_TYPE_FIXED);
  /*******************************************************************************************************************************************************/
  /* GAL DUST EPSILON TE 100 x 143 */    
  /*******************************************************************************************************************************************************/
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_0T_1E, "\\epsilon^{\\mathrm{bleak}TE}_{0, 100 \\times 143}", "bleak_epsilon_0_0T_1E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_0T_1E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_0T_1E, "\\epsilon^{\\mathrm{bleak}TE}_{1, 100 \\times 143}", "bleak_epsilon_1_0T_1E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_0T_1E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_0T_1E, "\\epsilon^{\\mathrm{bleak}TE}_{2, 100 \\times 143}", "bleak_epsilon_2_0T_1E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_0T_1E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_0T_1E, "\\epsilon^{\\mathrm{bleak}TE}_{3, 100 \\times 143}", "bleak_epsilon_3_0T_1E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_0T_1E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_0T_1E, "\\epsilon^{\\mathrm{bleak}TE}_{4, 100 \\times 143}", "bleak_epsilon_4_0T_1E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_0T_1E,
                              NCM_PARAM_TYPE_FIXED);
  /*******************************************************************************************************************************************************/
  /* GAL DUST EPSILON TE 100 x 217 */  
  /*******************************************************************************************************************************************************/
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_0T_2E, "\\epsilon^{\\mathrm{bleak}TE}_{0, 100 \\times 217}", "bleak_epsilon_0_0T_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_0T_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_0T_2E, "\\epsilon^{\\mathrm{bleak}TE}_{1, 100 \\times 217}", "bleak_epsilon_1_0T_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_0T_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_0T_2E, "\\epsilon^{\\mathrm{bleak}TE}_{2, 100 \\times 217}", "bleak_epsilon_2_0T_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_0T_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_0T_2E, "\\epsilon^{\\mathrm{bleak}TE}_{3, 100 \\times 217}", "bleak_epsilon_3_0T_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_0T_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_0T_2E, "\\epsilon^{\\mathrm{bleak}TE}_{4, 100 \\times 217}", "bleak_epsilon_4_0T_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_0T_2E,
                              NCM_PARAM_TYPE_FIXED);
  /*******************************************************************************************************************************************************/
  /* GAL DUST EPSILON TE 143 */  
  /*******************************************************************************************************************************************************/
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_1T_1E, "\\epsilon^{\\mathrm{bleak}TE}_{0, 143}", "bleak_epsilon_0_1T_1E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_1T_1E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_1T_1E, "\\epsilon^{\\mathrm{bleak}TE}_{1, 143}", "bleak_epsilon_1_1T_1E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_1T_1E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_1T_1E, "\\epsilon^{\\mathrm{bleak}TE}_{2, 143}", "bleak_epsilon_2_1T_1E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_1T_1E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_1T_1E, "\\epsilon^{\\mathrm{bleak}TE}_{3, 143}", "bleak_epsilon_3_1T_1E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_1T_1E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_1T_1E, "\\epsilon^{\\mathrm{bleak}TE}_{4, 143}", "bleak_epsilon_4_1T_1E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_1T_1E,
                              NCM_PARAM_TYPE_FIXED);
  /*******************************************************************************************************************************************************/
  /* GAL DUST EPSILON TE 143 x 217 */  
  /*******************************************************************************************************************************************************/
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_1T_2E, "\\epsilon^{\\mathrm{bleak}TE}_{0, 143 \\times 217}", "bleak_epsilon_0_1T_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_1T_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_1T_2E, "\\epsilon^{\\mathrm{bleak}TE}_{1, 143 \\times 217}", "bleak_epsilon_1_1T_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_1T_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_1T_2E, "\\epsilon^{\\mathrm{bleak}TE}_{2, 143 \\times 217}", "bleak_epsilon_2_1T_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_1T_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_1T_2E, "\\epsilon^{\\mathrm{bleak}TE}_{3, 143 \\times 217}", "bleak_epsilon_3_1T_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_1T_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_1T_2E, "\\epsilon^{\\mathrm{bleak}TE}_{4, 143 \\times 217}", "bleak_epsilon_4_1T_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_1T_2E,
                              NCM_PARAM_TYPE_FIXED);
  /*******************************************************************************************************************************************************/
  /* GAL DUST EPSILON TE 217 */
  /*******************************************************************************************************************************************************/
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_2T_2E, "\\epsilon^{\\mathrm{bleak}TE}_{0, 217}", "bleak_epsilon_0_2T_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_2T_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_2T_2E, "\\epsilon^{\\mathrm{bleak}TE}_{1, 217}", "bleak_epsilon_1_2T_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_2T_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_2T_2E, "\\epsilon^{\\mathrm{bleak}TE}_{2, 217}", "bleak_epsilon_2_2T_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_2T_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_2T_2E, "\\epsilon^{\\mathrm{bleak}TE}_{3, 217}", "bleak_epsilon_3_2T_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_2T_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_2T_2E, "\\epsilon^{\\mathrm{bleak}TE}_{4, 217}", "bleak_epsilon_4_2T_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_2T_2E,
                              NCM_PARAM_TYPE_FIXED);
  /*******************************************************************************************************************************************************/
  /* GAL DUST EPSILON EE 100 */    
  /*******************************************************************************************************************************************************/
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_0E_0E, "\\epsilon^{\\mathrm{bleak}EE}_{0, 100}", "bleak_epsilon_0_0E_0E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_0E_0E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_0E_0E, "\\epsilon^{\\mathrm{bleak}EE}_{1, 100}", "bleak_epsilon_1_0E_0E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_0E_0E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_0E_0E, "\\epsilon^{\\mathrm{bleak}EE}_{2, 100}", "bleak_epsilon_2_0E_0E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_0E_0E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_0E_0E, "\\epsilon^{\\mathrm{bleak}EE}_{3, 100}", "bleak_epsilon_3_0E_0E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_0E_0E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_0E_0E, "\\epsilon^{\\mathrm{bleak}EE}_{4, 100}", "bleak_epsilon_4_0E_0E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_0E_0E,
                              NCM_PARAM_TYPE_FIXED);
  /*******************************************************************************************************************************************************/
  /* GAL DUST EPSILON EE 100 x 143 */    
  /*******************************************************************************************************************************************************/
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_0E_1E, "\\epsilon^{\\mathrm{bleak}EE}_{0, 100 \\times 143}", "bleak_epsilon_0_0E_1E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_0E_1E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_0E_1E, "\\epsilon^{\\mathrm{bleak}EE}_{1, 100 \\times 143}", "bleak_epsilon_1_0E_1E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_0E_1E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_0E_1E, "\\epsilon^{\\mathrm{bleak}EE}_{2, 100 \\times 143}", "bleak_epsilon_2_0E_1E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_0E_1E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_0E_1E, "\\epsilon^{\\mathrm{bleak}EE}_{3, 100 \\times 143}", "bleak_epsilon_3_0E_1E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_0E_1E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_0E_1E, "\\epsilon^{\\mathrm{bleak}EE}_{4, 100 \\times 143}", "bleak_epsilon_4_0E_1E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_0E_1E,
                              NCM_PARAM_TYPE_FIXED);
  /*******************************************************************************************************************************************************/
  /* GAL DUST EPSILON EE 100 x 217 */  
  /*******************************************************************************************************************************************************/
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_0E_2E, "\\epsilon^{\\mathrm{bleak}EE}_{0, 100 \\times 217}", "bleak_epsilon_0_0E_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_0E_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_0E_2E, "\\epsilon^{\\mathrm{bleak}EE}_{1, 100 \\times 217}", "bleak_epsilon_1_0E_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_0E_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_0E_2E, "\\epsilon^{\\mathrm{bleak}EE}_{2, 100 \\times 217}", "bleak_epsilon_2_0E_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_0E_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_0E_2E, "\\epsilon^{\\mathrm{bleak}EE}_{3, 100 \\times 217}", "bleak_epsilon_3_0E_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_0E_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_0E_2E, "\\epsilon^{\\mathrm{bleak}EE}_{4, 100 \\times 217}", "bleak_epsilon_4_0E_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_0E_2E,
                              NCM_PARAM_TYPE_FIXED);
  /*******************************************************************************************************************************************************/
  /* GAL DUST EPSILON EE 143 */  
  /*******************************************************************************************************************************************************/
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_1E_1E, "\\epsilon^{\\mathrm{bleak}EE}_{0, 143}", "bleak_epsilon_0_1E_1E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_1E_1E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_1E_1E, "\\epsilon^{\\mathrm{bleak}EE}_{1, 143}", "bleak_epsilon_1_1E_1E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_1E_1E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_1E_1E, "\\epsilon^{\\mathrm{bleak}EE}_{2, 143}", "bleak_epsilon_2_1E_1E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_1E_1E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_1E_1E, "\\epsilon^{\\mathrm{bleak}EE}_{3, 143}", "bleak_epsilon_3_1E_1E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_1E_1E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_1E_1E, "\\epsilon^{\\mathrm{bleak}EE}_{4, 143}", "bleak_epsilon_4_1E_1E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_1E_1E,
                              NCM_PARAM_TYPE_FIXED);
  /*******************************************************************************************************************************************************/
  /* GAL DUST EPSILON EE 143 x 217 */  
  /*******************************************************************************************************************************************************/
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_1E_2E, "\\epsilon^{\\mathrm{bleak}EE}_{0, 143 \\times 217}", "bleak_epsilon_0_1E_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_1E_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_1E_2E, "\\epsilon^{\\mathrm{bleak}EE}_{1, 143 \\times 217}", "bleak_epsilon_1_1E_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_1E_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_1E_2E, "\\epsilon^{\\mathrm{bleak}EE}_{2, 143 \\times 217}", "bleak_epsilon_2_1E_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_1E_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_1E_2E, "\\epsilon^{\\mathrm{bleak}EE}_{3, 143 \\times 217}", "bleak_epsilon_3_1E_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_1E_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_1E_2E, "\\epsilon^{\\mathrm{bleak}EE}_{4, 143 \\times 217}", "bleak_epsilon_4_1E_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_1E_2E,
                              NCM_PARAM_TYPE_FIXED);
  /*******************************************************************************************************************************************************/
  /* GAL DUST EPSILON EE 217 */
  /*******************************************************************************************************************************************************/
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_2E_2E, "\\epsilon^{\\mathrm{bleak}EE}_{0, 217}", "bleak_epsilon_0_2E_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_2E_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_2E_2E, "\\epsilon^{\\mathrm{bleak}EE}_{1, 217}", "bleak_epsilon_1_2E_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_2E_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_2E_2E, "\\epsilon^{\\mathrm{bleak}EE}_{2, 217}", "bleak_epsilon_2_2E_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_2E_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_2E_2E, "\\epsilon^{\\mathrm{bleak}EE}_{3, 217}", "bleak_epsilon_3_2E_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_2E_2E,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_2E_2E, "\\epsilon^{\\mathrm{bleak}EE}_{4, 217}", "bleak_epsilon_4_2E_2E",
                              -1.0, 1.0, 1.0e-4,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_2E_2E,
                              NCM_PARAM_TYPE_FIXED);
  /*******************************************************************************************************************************************************/
  /* Amplitude and calibration */
  /*******************************************************************************************************************************************************/
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_calib_100P, "c_{100P}", "calib_100P",
                              0.5, 1.5, 1.0e-1,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_calib_100P,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_calib_143P, "c_{143P}", "calib_143P",
                              0.5, 1.5, 1.0e-1,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_calib_143P,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_calib_217P, "c_{217P}", "calib_217P",
                              0.5, 1.5, 1.0e-1,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_calib_217P,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TTTEEE_A_pol, "A_{\\mathrm{pol}}", "A_pol",
                              0.5, 1.5, 1.0e-1,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TTTEEE_DEFAULT_A_pol,
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);
}

/**
 * nc_planck_fi_cor_ttteee_add_galf_priors:
 * @lh: a #NcmLikelihood
 * @mean: a vector containing the means
 * @sigma: a vector containing the standard deviations
 *
 * Add the galaxy dust priors (on polarization data) as described in [Planck 2015 results XI (2015)][XPlanckCollaboration2015a].
 *
 */
void
nc_planck_fi_cor_ttteee_add_galf_priors (NcmLikelihood *lh, NcmVector *mean, NcmVector *sigma)
{
  g_assert_cmpuint (ncm_vector_len (mean), ==, 12);
  g_assert_cmpuint (ncm_vector_len (sigma), ==, 12);

  ncm_likelihood_priors_add_gauss_param (lh, nc_planck_fi_id (), NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_100,     ncm_vector_get (mean, 0),  ncm_vector_get (sigma, 0));
  ncm_likelihood_priors_add_gauss_param (lh, nc_planck_fi_id (), NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_100_143, ncm_vector_get (mean, 1),  ncm_vector_get (sigma, 1));
  ncm_likelihood_priors_add_gauss_param (lh, nc_planck_fi_id (), NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_100_217, ncm_vector_get (mean, 2),  ncm_vector_get (sigma, 2));
  ncm_likelihood_priors_add_gauss_param (lh, nc_planck_fi_id (), NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_143,     ncm_vector_get (mean, 3),  ncm_vector_get (sigma, 3));
  ncm_likelihood_priors_add_gauss_param (lh, nc_planck_fi_id (), NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_143_217, ncm_vector_get (mean, 4),  ncm_vector_get (sigma, 4));
  ncm_likelihood_priors_add_gauss_param (lh, nc_planck_fi_id (), NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_217,     ncm_vector_get (mean, 5),  ncm_vector_get (sigma, 5));
  ncm_likelihood_priors_add_gauss_param (lh, nc_planck_fi_id (), NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_100,     ncm_vector_get (mean, 6),  ncm_vector_get (sigma, 6));
  ncm_likelihood_priors_add_gauss_param (lh, nc_planck_fi_id (), NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_100_143, ncm_vector_get (mean, 7),  ncm_vector_get (sigma, 7));
  ncm_likelihood_priors_add_gauss_param (lh, nc_planck_fi_id (), NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_100_217, ncm_vector_get (mean, 8),  ncm_vector_get (sigma, 8));
  ncm_likelihood_priors_add_gauss_param (lh, nc_planck_fi_id (), NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_143,     ncm_vector_get (mean, 9),  ncm_vector_get (sigma, 9));
  ncm_likelihood_priors_add_gauss_param (lh, nc_planck_fi_id (), NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_143_217, ncm_vector_get (mean, 10), ncm_vector_get (sigma, 10));
  ncm_likelihood_priors_add_gauss_param (lh, nc_planck_fi_id (), NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_217,     ncm_vector_get (mean, 11), ncm_vector_get (sigma, 11));
}

/**
 * nc_planck_fi_cor_ttteee_add_default_galf_priors:
 * @lh: a #NcmLikelihood
 *
 * Add the galaxy dust priors (on polarization data) as described in [Planck 2015 results XI (2015)][XPlanckCollaboration2015a].
 * It uses the default values.
 *
 */
void
nc_planck_fi_cor_ttteee_add_default_galf_priors (NcmLikelihood *lh)
{
  gdouble mean[12]  = {0.060, 0.050, 0.110, 0.10, 0.240, 0.72, 0.140, 0.120, 0.30, 0.240, 0.600, 1.80};
  gdouble sigma[12] = {0.012, 0.015, 0.033, 0.02, 0.048, 0.14, 0.042, 0.036, 0.09, 0.072, 0.180, 0.54};
  NcmVector *mean_vec = ncm_vector_new_data_static (mean, 12, 1);
  NcmVector *sigma_vec = ncm_vector_new_data_static (sigma, 12, 1);

  nc_planck_fi_cor_ttteee_add_galf_priors (lh, mean_vec, sigma_vec);

  ncm_vector_free (mean_vec);
  ncm_vector_free (sigma_vec);
}

/**
 * nc_planck_fi_cor_ttteee_add_all_default_priors:
 * @lh: a #NcmLikelihood
 *
 * Adds all default priors above:
 * - nc_planck_fi_cor_ttteee_add_default_gal_priors()
 * - nc_planck_fi_cor_ttteee_add_default_calib_priors()
 * - nc_planck_fi_cor_ttteee_add_default_sz_prior()
 *
 */
void
nc_planck_fi_cor_ttteee_add_all_default_priors (NcmLikelihood *lh)
{
  nc_planck_fi_cor_tt_add_all_default_priors (lh);
  nc_planck_fi_cor_ttteee_add_default_galf_priors (lh);
}
