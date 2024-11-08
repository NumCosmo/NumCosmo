/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_cfg.h
 *
 *  Wed Aug 13 20:58:50 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/math/ncm_mpi_job.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_rng.h>
#include <gmp.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

GQuark ncm_cfg_error_quark (void);

#define NCM_CFG_ERROR (ncm_cfg_error_quark ())


/**
 * NcmCfgError:
 * @NCM_CFG_ERROR_INVALID_FFTW_FLAG: Invalid FFTW flag.
 * @NCM_CFG_ERROR_INVALID_FFTW_FLAG_STRING: Invalid FFTW flag string.
 * @NCM_CFG_ERROR_INVALID_FFTW_TIMELIMIT: Invalid FFTW planner timelimit.
 *
 * Error codes for the ncm_cfg namespace.
 *
 */
typedef enum _NcmCfgError
{
  NCM_CFG_ERROR_INVALID_FFTW_FLAG,
  NCM_CFG_ERROR_INVALID_FFTW_FLAG_STRING,
  NCM_CFG_ERROR_INVALID_FFTW_TIMELIMIT,
} NcmCfgError;

typedef void (*NcmCfgLoggerFunc) (const gchar *msg);

void ncm_cfg_init (void);
void ncm_cfg_init_full_ptr (gint *argc, gchar ***argv);
gchar **ncm_cfg_init_full (gint argc, gchar **argv);

void ncm_cfg_enable_gsl_err_handler (void);
void ncm_cfg_register_obj (GType obj);
guint ncm_cfg_mpi_nslaves (void);
gchar *ncm_cfg_get_fullpath (const gchar *filename, ...);
const gchar *ncm_cfg_get_fullpath_base (void);

void ncm_cfg_keyfile_to_arg (GKeyFile *kfile, const gchar *group_name, GOptionEntry *entries, gchar **argv, gint *argc);
void ncm_cfg_entries_to_keyfile (GKeyFile *kfile, const gchar *group_name, GOptionEntry *entries);
gchar *ncm_cfg_string_to_comment (const gchar *str);
const GEnumValue *ncm_cfg_get_enum_by_id_name_nick (GType enum_type, const gchar *id_name_nick);
const GEnumValue *ncm_cfg_enum_get_value (GType enum_type, guint n);

void ncm_cfg_enum_print_all (GType enum_type, const gchar *header);

void ncm_cfg_lock_plan_fftw (void);
void ncm_cfg_unlock_plan_fftw (void);
gboolean ncm_cfg_load_fftw_wisdom (const gchar *filename, ...);
gboolean ncm_cfg_save_fftw_wisdom (const gchar *filename, ...);
gboolean ncm_cfg_exists (const gchar *filename, ...);

void ncm_cfg_set_logfile (gchar *filename);
void ncm_cfg_set_logstream (FILE *stream);
void ncm_cfg_set_log_handler (NcmCfgLoggerFunc logger);
void ncm_cfg_set_error_log_handler (NcmCfgLoggerFunc logger);

void ncm_cfg_set_openmp_nthreads (gint n);
void ncm_cfg_set_openblas_nthreads (gint n);
void ncm_cfg_set_blis_nthreads (gint n);
void ncm_cfg_set_mkl_nthreads (gint n);

void ncm_cfg_logfile (gboolean on);
void ncm_cfg_logfile_flush (gboolean on);
void ncm_cfg_logfile_flush_now (void);

void ncm_message_str (const gchar *msg);
void ncm_message (const gchar *msg, ...) G_GNUC_PRINTF (1, 2);
gchar *ncm_string_ww (const gchar *msg, const gchar *first, const gchar *rest, guint ncols);

void ncm_message_ww (const gchar *msg, const gchar *first, const gchar *rest, guint ncols);
void ncm_cfg_msg_sepa (void);

gchar *ncm_cfg_get_data_filename (const gchar *filename, gboolean must_exist);
gchar *ncm_cfg_get_data_directory (void);

gchar *ncm_cfg_command_line (gchar *argv[], gint argc);

void ncm_cfg_array_set_variant (GArray *a, GVariant *var);
GVariant *ncm_cfg_array_to_variant (GArray *a, const GVariantType *etype);

void ncm_cfg_set_fftw_default_flag (guint flag, const gdouble timeout, GError **error);
void ncm_cfg_set_fftw_default_flag_str (const gchar *flag_str, const gdouble timeout, GError **error);
void ncm_cfg_set_fftw_default_from_env (guint fallback_flag, const gdouble fallback_timeout, GError **error);
guint ncm_cfg_get_fftw_default_flag (void);
const gchar *ncm_cfg_get_fftw_default_flag_str (void);
gdouble ncm_cfg_get_fftw_timelimit (void);

/* Macros */

#define NCM_CFG_DATA_DIR_ENV "NUMCOSMO_DATA_DIR"

#ifdef NUMCOSMO_CHECK_PREPARE
#define NCM_CHECK_PREPARED(obj, name)                                      \
        G_STMT_START {                                                     \
          if (!obj->prepared)                                              \
          g_error ("calling method %s on an unprepared instance.", #name); \
        } G_STMT_END
#else /* NUMCOSMO_CHECK_PREPARE */
#define NCM_CHECK_PREPARED(obj, name)
#endif /* NUMCOSMO_CHECK_PREPARE */

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

#ifndef NUMCOSMO_GIR_SCAN
#define NCM_FITS_ERROR(status)                     \
        G_STMT_START {                             \
          if (status)                              \
          {                                        \
            gchar errormsg[30];                    \
            fits_get_errstatus (status, errormsg); \
            g_error ("FITS: %s", errormsg);        \
          }                                        \
        } G_STMT_END
#endif /* NUMCOSMO_GIR_SCAN */

G_END_DECLS

#endif /* _NCM_CFG_H */

