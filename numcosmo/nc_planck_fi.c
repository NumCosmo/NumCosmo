/***************************************************************************
 *            nc_planck_fi.c
 *
 *  Thu October 22 15:48:37 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_planck_fi.c
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
 * SECTION:nc_planck_fi
 * @title: NcPlanckFI
 * @short_description: Abstract class for Planck Foreground and Instrument models.
 *
 * NcPlanclFI is the abstract class designed describe a generic
 * model for Planck foreground and instruments.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_planck_fi.h"

enum
{
  PROP_0,
  PROP_VERSION,
  PROP_SIZE
};

G_DEFINE_TYPE (NcPlanckFI, nc_planck_fi, NCM_TYPE_MODEL);

static void
nc_planck_fi_init (NcPlanckFI *nc_planck_fi)
{
}

static void
nc_planck_fi_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcPlanckFI *pfi = NC_PLANCK_FI (object);
  g_return_if_fail (NC_IS_PLANCK_FI (object));

  switch (prop_id)
  {
    case PROP_VERSION:
      pfi->version = g_value_get_uint (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_planck_fi_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcPlanckFI *pfi = NC_PLANCK_FI (object);
  g_return_if_fail (NC_IS_PLANCK_FI (object));

  switch (prop_id)
  {
    case PROP_VERSION:
      g_value_set_uint (value, pfi->version);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_planck_fi_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_planck_fi_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_planck_fi, NC_TYPE_PLANCK_FI);

static void
nc_planck_fi_class_init (NcPlanckFIClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = nc_planck_fi_set_property;
  model_class->get_property = nc_planck_fi_get_property;
  object_class->finalize = nc_planck_fi_finalize;

  ncm_model_class_set_name_nick (model_class, "Planck Foreground and Instrument Abstract Class", "PlanckFI");
  ncm_model_class_add_params (NCM_MODEL_CLASS (klass), 0, 0, PROP_SIZE);

  ncm_mset_model_register_id (NCM_MODEL_CLASS (klass),
                              "NcPlanckFI",
                              "Planck Foreground and Instrument models.",
                              NULL,
                              FALSE,
                              NCM_MSET_MODEL_MAIN);

  g_object_class_install_property (object_class,
                                   PROP_VERSION,
                                   g_param_spec_uint ("version",
                                                      NULL,
                                                      "Planck compatible version",
                                                      0, G_MAXUINT, 2013,
                                                      G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  ncm_model_class_check_params_info (NCM_MODEL_CLASS (klass));
}

static void
_nc_planck_fi_log_all_models_go (GType model_type, guint n)
{
  guint nc, i, j;
  GType *models = g_type_children (model_type, &nc);
  for (i = 0; i < nc; i++)
  {
    guint ncc;
    GType *modelsc = g_type_children (models[i], &ncc);

    g_message ("#  ");
    for (j = 0; j < n; j++) g_message (" ");
    g_message ("%s\n", g_type_name (models[i]));
    if (ncc)
      _nc_planck_fi_log_all_models_go (models[i], n + 2);

    g_free (modelsc);
  }
  g_free (models);
}

/**
 * nc_planck_fi_log_all_models:
 *
 * This function lists all implemented models of cluster mass distributions.
 *
 */
void
nc_planck_fi_log_all_models (void)
{
  g_message ("# Registred NcPlanckFI:%s are:\n", g_type_name (NC_TYPE_PLANCK_FI));
  _nc_planck_fi_log_all_models_go (NC_TYPE_PLANCK_FI, 0);
}



/**
 * nc_planck_fi_new_from_name:
 * @pfi_name: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcPlanckFI *
nc_planck_fi_new_from_name (gchar *pfi_name)
{
  GType parent_type = NC_TYPE_PLANCK_FI;
  GObject *obj = ncm_serialize_global_from_string (pfi_name);
  GType model_type = G_OBJECT_TYPE (obj);

  if (!g_type_is_a (model_type, parent_type))
    g_error ("nc_planck_fi_new_from_name: NcPlanckFI %s do not descend from %s.", pfi_name, g_type_name (parent_type));
  return NC_PLANCK_FI (obj);
}

/**
 * nc_planck_fi_ref:
 * @pfi: a #NcPlanckFI
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcPlanckFI *
nc_planck_fi_ref (NcPlanckFI *pfi)
{
  return g_object_ref (pfi);
}

/**
 * nc_planck_fi_free:
 * @pfi: a #NcPlanckFI
 *
 * FIXME
 *
 */
void
nc_planck_fi_free (NcPlanckFI *pfi)
{
  g_object_unref (pfi);
}

/**
 * nc_planck_fi_clear:
 * @pfi: a #NcPlanckFI
 *
 * The reference count of @pfi is decreased and the pointer is set to NULL.
 *
 */
void
nc_planck_fi_clear (NcPlanckFI **pfi)
{
  g_clear_object (pfi);
}
