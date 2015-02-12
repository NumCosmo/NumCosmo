/***************************************************************************
 *            ncm_cfg.h
 *
 *  Wed Aug 13 20:58:50 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#include <stdio.h>
#include <glib.h>
#include <glib/gstdio.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#ifndef NUMCOSMO_GIR_SCAN
#include <gmp.h>
#endif /* NUMCOSMO_GIR_SCAN */
#ifdef NUMCOSMO_HAVE_SQLITE3
#include <sqlite3.h>
#endif
#include <numcosmo/math/ncm_spline.h>

G_BEGIN_DECLS

void ncm_cfg_init (void);
void ncm_cfg_enable_gsl_err_handler (void);
void ncm_cfg_register_obj (GType obj);
gchar *ncm_cfg_get_fullpath (const gchar *filename, ...);
void ncm_cfg_keyfile_to_arg (GKeyFile *kfile, gchar *group_name, GOptionEntry *entries, gchar **argv, gint *argc);
void ncm_cfg_entries_to_keyfile (GKeyFile *kfile, gchar *group_name, GOptionEntry *entries);
gchar *ncm_cfg_string_to_comment (const gchar *str);
const GEnumValue *ncm_cfg_get_enum_by_id_name_nick (GType enum_type, const gchar *id_name_nick);
const GEnumValue *ncm_cfg_enum_get_value (GType enum_type, guint n);
void ncm_cfg_enum_print_all (GType enum_type, gchar *header);

gboolean ncm_cfg_load_fftw_wisdom (gchar *filename, ...);
gboolean ncm_cfg_save_fftw_wisdom (gchar *filename, ...);
gboolean ncm_cfg_exists (gchar *filename, ...);

void ncm_cfg_set_logfile (gchar *filename);
void ncm_cfg_logfile (gboolean on);
void ncm_cfg_logfile_flush (gboolean on);
void ncm_cfg_logfile_flush_now (void);

void ncm_message (gchar *msg, ...);
gchar *ncm_string_ww (const gchar *msg, const gchar *first, const gchar *rest, guint ncols);
void ncm_message_ww (const gchar *msg, const gchar *first, const gchar *rest, guint ncols);
void ncm_cfg_msg_sepa (void);

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

gchar *ncm_cfg_get_data_filename (const gchar *filename, gboolean must_exist);

#ifdef NUMCOSMO_HAVE_SQLITE3
sqlite3 *ncm_cfg_get_default_sqlite3 (void);
#endif

typedef union _NcmDoubleInt64
{
  gint64 i;
  gdouble x;
} NcmDoubleInt64;

gchar *ncm_cfg_command_line (gchar *argv[], gint argc);

GArray *ncm_cfg_variant_to_array (GVariant *var, gsize esize);
void ncm_cfg_array_set_variant (GArray *a, GVariant *var);
GVariant *ncm_cfg_array_to_variant (GArray *a, const GVariantType *etype);

void ncm_cfg_set_fftw_default_flag (guint flag);

extern guint fftw_default_flags;

/* Macros */

#ifdef SUNDIALS_USES_LONG_INT
#define _NCM_SUNDIALS_INT_TYPE glong
#else
#define _NCM_SUNDIALS_INT_TYPE gint
#endif

#define NCM_CFG_DATA_DIR_ENV "NUMCOSMO_DATA_DIR"

#if (GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 32)
#define _NCM_MUTEX_LOCK(l) g_static_mutex_lock (l)
#define _NCM_MUTEX_UNLOCK(l) g_static_mutex_unlock (l)
#define _NCM_MUTEX_TRYLOCK(l) g_static_mutex_trylock (l)
#define _NCM_MUTEX_TYPE GStaticMutex
#define _NCM_STATIC_MUTEX_DECL(l) static GStaticMutex l = G_STATIC_MUTEX_INIT
#define _NCM_MUTEX_INIT(l) g_static_mutex_init (l)
#define _NCM_MUTEX_CLEAR(l) g_static_mutex_free (l)
#else
#define _NCM_MUTEX_LOCK(l) g_mutex_lock (l)
#define _NCM_MUTEX_UNLOCK(l) g_mutex_unlock (l)
#define _NCM_MUTEX_TRYLOCK(l) g_mutex_trylock (l)
#define _NCM_MUTEX_TYPE GMutex
#define _NCM_STATIC_MUTEX_DECL(l) static GMutex l
#define _NCM_MUTEX_INIT(l) g_mutex_init (l)
#define _NCM_MUTEX_CLEAR(l) g_mutex_clear (l)
#endif

#define NCM_ZERO_LIMIT 1e-13
#define NCM_DEFAULT_PRECISION 1e-7

#ifndef NCM_THREAD_POOL_MAX
#define NCM_THREAD_POOL_MAX 5
#endif

#ifndef mpz_inits
#define mpz_inits ncm_mpz_inits
#endif /* mpz_inits */

#ifndef mpz_clears
#define mpz_clears ncm_mpz_clears
#endif /* mpz_inits */

#define NCM_CFG_DEFAULT_SQLITE3_FILENAME "data_observation.sqlite3"

/* Workaround on g_clear_pointer */
#if ((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 34))
#define g_clear_pointer(ptr,freefunc) \
G_STMT_START { \
  if (*(ptr) != NULL) \
  { \
    freefunc (*(ptr)); \
    *(ptr) = NULL; \
  } \
} G_STMT_END
#endif /* Glib version < 2.34 */

/* Workaround on g_clear_object. */
#if ((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 28))
#define g_clear_object(obj) g_clear_pointer(ptr,g_object_unref)
#endif /* Glib version < 2.28 */

#define NCM_FITS_ERROR(status) \
G_STMT_START { \
  if (status) \
  { \
    gchar errormsg[30]; \
    fits_get_errstatus (status, errormsg); \
    g_error ("FITS: %s", errormsg); \
  } \
} G_STMT_END

G_END_DECLS

#endif /* _NCM_CFG_H */
