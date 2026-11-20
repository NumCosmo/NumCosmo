/***************************************************************************
 *            ncm_mset_preset.h
 *
 *  Thu Nov 19 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mset_preset.h
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_MSET_PRESET_H_
#define _NCM_MSET_PRESET_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/core/ncm_obj_array.h>
#include <numcosmo/ncm/model/ncm_mset.h>

G_BEGIN_DECLS

#define NCM_TYPE_MSET_PRESET (ncm_mset_preset_get_type ())

G_DECLARE_FINAL_TYPE (NcmMSetPreset, ncm_mset_preset, NCM, MSET_PRESET, GObject)

/**
 * NcmMSetPresetError:
 * @NCM_MSET_PRESET_ERROR_NOT_FOUND: no preset file has the requested name
 * @NCM_MSET_PRESET_ERROR_INVALID_FILE: the file does not hold an #NcmMSetPreset
 * @NCM_MSET_PRESET_ERROR_PARAM_NOT_FOUND: the #NcmMSet has no model or no parameter with
 *   a full name used by the preset
 * @NCM_MSET_PRESET_ERROR_REPARAM_MISMATCH: a model touched by the preset has a
 *   reparametrization other than the one the preset requires
 * @NCM_MSET_PRESET_ERROR_INVALID_ENTRY: an entry has the wrong type, or the resulting
 *   value, bounds, scale or absolute tolerance are inconsistent
 *
 * Error codes returned by the #NcmMSetPreset class.
 *
 */
typedef enum _NcmMSetPresetError /*< prefix=NCM_MSET_PRESET_ERROR >*/
{
  NCM_MSET_PRESET_ERROR_NOT_FOUND,
  NCM_MSET_PRESET_ERROR_INVALID_FILE,
  NCM_MSET_PRESET_ERROR_PARAM_NOT_FOUND,
  NCM_MSET_PRESET_ERROR_REPARAM_MISMATCH,
  NCM_MSET_PRESET_ERROR_INVALID_ENTRY,
} NcmMSetPresetError;

GQuark ncm_mset_preset_error_quark (void);

#define NCM_MSET_PRESET_ERROR (ncm_mset_preset_error_quark ())

/**
 * NCM_MSET_PRESET_DATA_DIR:
 *
 * The subdirectory of the data directory that holds the preset files.
 */
#define NCM_MSET_PRESET_DATA_DIR "mset_presets"

NcmMSetPreset *ncm_mset_preset_new (const gchar *desc);
NcmMSetPreset *ncm_mset_preset_new_from_file (const gchar *filename, GError **error);
NcmMSetPreset *ncm_mset_preset_new_from_name (const gchar *name, GError **error);
NcmMSetPreset *ncm_mset_preset_ref (NcmMSetPreset *preset);
void ncm_mset_preset_free (NcmMSetPreset *preset);
void ncm_mset_preset_clear (NcmMSetPreset **preset);

GStrv ncm_mset_preset_list_names (void);

const gchar *ncm_mset_preset_peek_desc (NcmMSetPreset *preset);
const gchar *ncm_mset_preset_peek_data_combination (NcmMSetPreset *preset);
const gchar *ncm_mset_preset_peek_statistic (NcmMSetPreset *preset);

void ncm_mset_preset_set_reparam (NcmMSetPreset *preset, const gchar *model_name, const gchar *reparam_type);
void ncm_mset_preset_set_value (NcmMSetPreset *preset, const gchar *fullname, gdouble value);
void ncm_mset_preset_set_scale (NcmMSetPreset *preset, const gchar *fullname, gdouble scale);
void ncm_mset_preset_set_lower_bound (NcmMSetPreset *preset, const gchar *fullname, gdouble lower_bound);
void ncm_mset_preset_set_upper_bound (NcmMSetPreset *preset, const gchar *fullname, gdouble upper_bound);
void ncm_mset_preset_set_abstol (NcmMSetPreset *preset, const gchar *fullname, gdouble abstol);
void ncm_mset_preset_set_fit (NcmMSetPreset *preset, const gchar *fullname, gboolean fit);

gboolean ncm_mset_preset_remove (NcmMSetPreset *preset, const gchar *fullname);

void ncm_mset_preset_apply (NcmMSetPreset *preset, NcmMSet *mset, GError **error);

G_END_DECLS

#endif /* _NCM_MSET_PRESET_H_ */

