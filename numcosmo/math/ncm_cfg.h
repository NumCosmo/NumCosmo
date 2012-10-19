/***************************************************************************
 *            ncm_cfg.h
 *
 *  Wed Aug 13 20:58:50 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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

#ifndef _NCM_CFG_H
#define _NCM_CFG_H

#include <glib.h>
#ifdef NUMCOSMO_HAVE_SQLITE3
#include <sqlite3.h>
#endif

G_BEGIN_DECLS

void ncm_cfg_init (void);
void ncm_cfg_enable_gsl_err_handler (void);
void ncm_cfg_register_model (GType model);
gchar *ncm_cfg_get_fullpath (gchar *filename, ...);
void ncm_cfg_keyfile_to_arg (GKeyFile *kfile, gchar *group_name, GOptionEntry *entries, gchar **argv, gint *argc);
void ncm_cfg_entries_to_keyfile (GKeyFile *kfile, gchar *group_name, GOptionEntry *entries);
const GEnumValue *ncm_cfg_get_enum_by_id_name_nick (GType enum_type, const gchar *id_name_nick);
void ncm_cfg_enum_print_all (GType enum_type, gchar *header);

gboolean ncm_cfg_load_fftw_wisdom (gchar *filename, ...);
gboolean ncm_cfg_save_fftw_wisdom (gchar *filename, ...);
gboolean ncm_cfg_exists (gchar *filename, ...);

void ncm_cfg_set_logfile (gchar *filename);
void ncm_cfg_logfile (gboolean on);
void ncm_cfg_logfile_flush (gboolean on);

void nc_message (gchar *msg, ...);

FILE *ncm_cfg_fopen (gchar *filename, gchar *mode, ...);
FILE *ncm_cfg_vfopen (gchar *filename, gchar *mode, va_list ap);

gboolean ncm_cfg_load_spline (gchar *filename, const gsl_interp_type *stype, NcmSpline **s, ...);
gboolean ncm_cfg_save_spline (gchar *filename, NcmSpline *s, ...);

gboolean ncm_cfg_load_vector (gchar *filename, gsl_vector *v, ...);
gboolean ncm_cfg_save_vector (gchar *filename, gsl_vector *v, ...);
gboolean ncm_cfg_load_matrix (gchar *filename, gsl_matrix *M, ...);
gboolean ncm_cfg_save_matrix (gchar *filename, gsl_matrix *M, ...);

#define LOAD_SAVE_VECTOR_MATRIX_DEF(typen) \
gboolean ncm_cfg_load_vector_##typen (gchar *filename, gsl_vector_##typen *v, ...); \
gboolean ncm_cfg_save_vector_##typen (gchar *filename, gsl_vector_##typen *v, ...); \
gboolean ncm_cfg_load_matrix_##typen (gchar *filename, gsl_matrix_##typen *M, ...); \
gboolean ncm_cfg_save_matrix_##typen (gchar *filename, gsl_matrix_##typen *M, ...);

LOAD_SAVE_VECTOR_MATRIX_DEF(int)
LOAD_SAVE_VECTOR_MATRIX_DEF(float)
LOAD_SAVE_VECTOR_MATRIX_DEF(complex)

#ifdef NUMCOSMO_HAVE_SQLITE3
sqlite3 *ncm_cfg_get_default_sqlite3 (void);
#endif

#define NC_CFG_DEFAULT_SQLITE3_FILENAME "data_observation.sqlite3"

typedef union _NcmDoubleInt64
{
  gint64 i;
  gdouble x;
} NcmDoubleInt64;

#define NCM_WRITE_INT32(_ff,_ii) do { gint32 _temp_i = GINT32_TO_BE ((_ii)); if (fwrite (&_temp_i, sizeof(gint32), (1), _ff) != 1) g_error ("NCM_WRITE_INT32: io error"); } while (FALSE)
#define NCM_WRITE_UINT32(_ff,_ii) do { guint32 _temp_i = GUINT32_TO_BE ((_ii)); if (fwrite (&_temp_i, sizeof(guint32), (1), _ff) != 1) g_error ("NCM_WRITE_UINT32: io error"); } while (FALSE)

#define NCM_WRITE_INT64(_ff,_ii) do { gint64 _temp_i = GINT64_TO_BE ((_ii)); if (fwrite (&_temp_i, sizeof(gint64), (1), _ff) != 1) g_error ("NCM_WRITE_INT64: io error"); } while (FALSE)
#define NCM_WRITE_UINT64(_ff,_ii) do { guint64 _temp_i = GUINT64_TO_BE ((_ii)); if (fwrite (&_temp_i, sizeof(guint64), (1), _ff) != 1) g_error ("NCM_WRITE_INT64: io error"); } while (FALSE)

#define NCM_WRITE_DOUBLE(_ff,_ii) do { NcmDoubleInt64 _iii; _iii.x = _ii; _iii.i = GINT64_TO_BE ((_iii.i)); if (fwrite (&_iii.i, sizeof(gint64), (1), _ff) != 1) g_error ("NCM_WRITE_DOUBLE: io error"); } while (FALSE)

#define NCM_READ_INT32(_ff,_ii) do { gint32 _temp_i; if (fread (&_temp_i, sizeof(gint32), (1), _ff) != 1) g_error ("NCM_READ_INT32: io error"); _ii = GINT32_FROM_BE (_temp_i); } while (FALSE)
#define NCM_READ_UINT32(_ff,_ii) do { guint32 _temp_i; if (fread (&_temp_i, sizeof(guint32), (1), _ff) != 1) g_error ("NCM_READ_UINT32: io error"); _ii = GUINT32_FROM_BE (_temp_i); } while (FALSE)

#define NCM_READ_INT64(_ff,_ii) do { gint64 _temp_i; if (fread (&_temp_i, sizeof(gint64), (1), _ff) != 1) g_error ("NCM_READ_INT64: io error"); _ii = GINT64_FROM_BE (_temp_i); } while (FALSE)
#define NCM_READ_UINT64(_ff,_ii) do { guint64 _temp_i; if (fread (&_temp_i, sizeof(guint64), (1), _ff) != 1) g_error ("NCM_READ_UINT64: io error"); _ii = GUINT64_FROM_BE (_temp_i); } while (FALSE)

#define NCM_READ_DOUBLE(_ff,_ii) do { NcmDoubleInt64 _iii; if (fread (&_iii.i, sizeof(gint64), (1), _ff) != 1) g_error ("NCM_READ_DOUBLE: io error"); _iii.i = GINT64_FROM_BE (_iii.i); _ii = _iii.x; } while (FALSE)

gchar *ncm_cfg_command_line (gchar *argv[], gint argc);

GObject *ncm_cfg_create_from_name_params (const gchar *obj_name, GVariant *params);
GObject *ncm_cfg_create_from_string (const gchar *obj_ser);
GVariant *ncm_cfg_gvalue_to_gvariant (GValue *val);
GVariant *ncm_cfg_serialize_to_variant (GObject *obj);
gchar *ncm_cfg_serialize_to_string (GObject *obj, gboolean valid_variant);


G_END_DECLS

#endif /* _NCM_CFG_H */
