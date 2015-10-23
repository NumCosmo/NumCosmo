/***************************************************************************
 *            nc_planck_fi_tt.c
 *
 *  Thu October 22 16:21:36 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_planck_fi_tt.c
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_planck_fi_tt
 * @title: NcPlanclFI_TT
 * @short_description: Planck Foreground and Instrument model for TT maps
 *
 * FIXME
 * see~<link linkend="XPlanckCollaboration2015a">Planck 2015 results
 *           XI (2015)</link>
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_planck_fi_tt.h"

enum
{
  PROP_0,
  PROP_CMP,
  PROP_SIZE
};

G_DEFINE_TYPE (NcPlanckFI_TT, nc_planck_fi_tt, NC_TYPE_PLANCK_FI);

static void
nc_planck_fi_tt_init (NcPlanckFI_TT *nc_planck_fi_tt)
{
}

static void
nc_planck_fi_tt_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_PLANCK_FI_TT (object));

  switch (prop_id)
  {
    case PROP_CMP:
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_planck_fi_tt_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_PLANCK_FI_TT (object));

  switch (prop_id)
  {
    case PROP_CMP:
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_planck_fi_tt_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_planck_fi_tt_parent_class)->finalize (object);
}

static void
nc_planck_fi_tt_class_init (NcPlanckFI_TTClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = nc_planck_fi_tt_set_property;
  model_class->get_property = nc_planck_fi_tt_get_property;
  object_class->finalize    = nc_planck_fi_tt_finalize;

  ncm_model_class_set_name_nick (model_class, "Planck Foreground and Instument Model -- TT", "PlanckFI_TT");
  ncm_model_class_add_params (model_class, NC_PLANCK_FI_TT_SPARAM_LEN, 0, PROP_SIZE);

  g_object_class_install_property (object_class,
                                   PROP_CMP,
                                   g_param_spec_uint ("cmp",
                                                      NULL,
                                                      "mmm",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_TT_A_cib_217, "A^{\\mathrm{CIB}}_{217}", "A_cib_217",
                              0.0, 400.0, 5.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_TT_DEFAULT_A_cib_217,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_TT_cib_index, "n^{\\mathrm{CIB}}", "cib_index",
                              -2.0, 0.0, 0.01,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_TT_DEFAULT_cib_index,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_TT_xi_sz_cib, "\\xi^{\\mathrm{tSZ}\\times \\mathrm{CIB}}", "xi_sz_cib",
                              0.0, 1.0, 0.1,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_TT_DEFAULT_xi_sz_cib,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_TT_A_sz, "A^{\\mathrm{tSZ}}", "A_sz",
                              0.0, 10.0, 0.5,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_TT_DEFAULT_A_sz,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_TT_ps_A_100_100, "A^{\\mathrm{PS}}_{100}", "ps_A_100_100",
                              0.0, 400.0, 5.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_TT_DEFAULT_ps_A_100_100,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_TT_ps_A_143_143, "A^{\\mathrm{PS}}_{143}", "ps_A_143_143",
                              0.0, 400.0, 5.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_TT_DEFAULT_ps_A_143_143,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_TT_ps_A_143_217, "A^{\\mathrm{PS}}_{143\\times 217}", "ps_A_143_217",
                              0.0, 400.0, 5.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_TT_DEFAULT_ps_A_143_217,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_TT_ps_A_217_217, "A^{\\mathrm{PS}}_{217}", "ps_A_217_217",
                              0.0, 400.0, 5.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_TT_DEFAULT_ps_A_217_217,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_TT_ksz_norm, "A^{\\mathrm{kSZ}}", "ksz_norm",
                              0.0, 10.0, 0.5,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_TT_DEFAULT_ksz_norm,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_TT_gal545_A_100, "A^{\\mathrm{dust}TT}_{100}", "gal545_A_100",
                              0.0, 50.0, 2.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_TT_DEFAULT_gal545_A_100,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_TT_gal545_A_143, "A^{\\mathrm{dust}TT}_{143}", "gal545_A_143",
                              0.0, 50.0, 2.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_TT_DEFAULT_gal545_A_143,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_TT_gal545_A_143_217, "A^{\\mathrm{dust}TT}_{143 \\times 217}", "gal545_A_143_217",
                              0.0, 100.0, 8.5,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_TT_DEFAULT_gal545_A_143_217,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_TT_gal545_A_217, "A^{\\mathrm{dust}TT}_{217}", "gal545_A_217",
                              0.0, 400.0, 20.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_TT_DEFAULT_gal545_A_217,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_TT_calib_100T, "c_{100}", "calib_100T",
                              0.0, 3.0, 0.001,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_TT_DEFAULT_calib_100T,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_TT_calib_217T, "c_{217}", "calib_217T",
                              0.0, 3.0, 0.002,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_TT_DEFAULT_calib_217T,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_TT_A_planck, "y_{\\mathrm{cal}}", "A_planck",
                              0.9, 1.1, 0.0025,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_TT_DEFAULT_A_planck,
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);
}
