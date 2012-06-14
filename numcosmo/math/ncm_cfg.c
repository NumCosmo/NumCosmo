/***************************************************************************
 *            ncm_cfg.c
 *
 *  Wed Aug 13 20:59:22 2008
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

/**
 * SECTION:ncm_cfg
 * @title: Library Configuration
 * @short_description: FIXME
 *
 * FIXME
 */


#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <string.h>
#include <glib.h>
#include <glib/gstdio.h>
#ifdef NUMCOSMO_HAVE_SQLITE3
#include <sqlite3.h>
#endif

static gchar *numcosmo_path = NULL;
static gboolean numcosmo_init = FALSE;
static FILE *_log_stream = NULL;
static guint _log_msg_id = 0;
static gboolean _enable_msg = TRUE;
static gboolean _enable_msg_flush = TRUE;
static gsl_error_handler_t *gsl_err = NULL;

static void
_ncm_cfg_log_message (const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer user_data)
{
  if (_enable_msg && _log_stream)
  {
	fprintf (_log_stream, message);
	if (_enable_msg_flush)
	  fflush (_log_stream);
  }
}

/**
 * ncm_cfg_init:
 *
 * FIXME
 */
void
ncm_cfg_init (void)
{
  const gchar *home;
  if (numcosmo_init)
	return;
  home = g_get_home_dir ();
  numcosmo_path = g_build_filename (home, ".numcosmo", NULL);
  if (!g_file_test (numcosmo_path, G_FILE_TEST_EXISTS))
	g_mkdir_with_parents (numcosmo_path, 0755);

  gsl_err = gsl_set_error_handler_off ();
  g_thread_init (NULL);
  g_type_init();

  _log_msg_id = g_log_set_handler (G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE | G_LOG_LEVEL_DEBUG, _ncm_cfg_log_message, NULL);
  _log_stream = stdout;

  ncm_cfg_register_model (NC_TYPE_MODEL_LCDM);
  ncm_cfg_register_model (NC_TYPE_MODEL_DE_XCDM);
  ncm_cfg_register_model (NC_TYPE_MODEL_DE_LINDER);
  ncm_cfg_register_model (NC_TYPE_MODEL_DE_PAD);
  ncm_cfg_register_model (NC_TYPE_MODEL_DE_QE);
  ncm_cfg_register_model (NC_TYPE_MODEL_QCONST);
  ncm_cfg_register_model (NC_TYPE_MODEL_QLINEAR);
  ncm_cfg_register_model (NC_TYPE_MODEL_QPW);
  ncm_cfg_register_model (NC_TYPE_MODEL_QSPLINE);

  ncm_cfg_register_model (NC_TYPE_WINDOW_TOPHAT);
  ncm_cfg_register_model (NC_TYPE_WINDOW_GAUSSIAN);

  ncm_cfg_register_model (NC_TYPE_TRANSFER_FUNC_EH);

  ncm_cfg_register_model (NC_TYPE_MULTIPLICITY_FUNC_TINKER_MEAN);

  numcosmo_init = TRUE;
  return;
}

/**
 * ncm_cfg_enable_gsl_err_handler:
 *
 * FIXME
 */
void
ncm_cfg_enable_gsl_err_handler (void)
{
  g_assert (numcosmo_init);
  gsl_set_error_handler (gsl_err);
}

static uint nreg_model = 0;

/**
 * ncm_cfg_register_model:
 * @model: FIXME
 *
 * FIXME
 */
void
ncm_cfg_register_model (GType model)
{
#ifdef NCM_DEBUG_MSGS
  GType model_type = NCM_TYPE_GMODEL;
  GType model_type = g_type_next_base (model, model_type);
  g_message ("# %s: %s:%s registred\n",
             g_type_name (model_type),
             g_type_name (model_type),
             g_type_name (model));
#endif /* NCM_DEBUG_MSGS */
  nreg_model++;
}

/**
 * ncm_cfg_set_logfile:
 * @filename: FIXME
 *
 * FIXME
 */
void
ncm_cfg_set_logfile (gchar *filename)
{
  FILE *out = fopen (filename, "w");
  if (out != NULL)
	_log_stream = out;
  else
  {
	fprintf (stderr, "Can't open logfile (%s)\n", filename);
	exit (-1);
  }
}

/**
 * ncm_cfg_logfile:
 * @on: FIXME
 *
 * FIXME
 */
void
ncm_cfg_logfile (gboolean on) { _enable_msg = on; }

/**
 * ncm_cfg_logfile_flush:
 * @on: FIXME
 *
 * FIXME
 */
void
ncm_cfg_logfile_flush (gboolean on) { _enable_msg_flush = on; }


/**
 * nc_message:
 * @msg: FIXME
 * @...: FIXME
 *
 * FIXME
 */
void
nc_message (gchar *msg, ...)
{
  va_list ap;
  va_start(ap, msg);
  g_logv (G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE, msg, ap);
  va_end (ap);
}

/**
 * ncm_cfg_get_fullpath:
 * @filename: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gchar *
ncm_cfg_get_fullpath (gchar *filename, ...)
{
  gchar *file, *full_filename;

  g_assert (numcosmo_init);
  va_list ap;
  va_start (ap, filename);
  file = g_strdup_vprintf (filename, ap);
  va_end (ap);

  full_filename = g_build_filename (numcosmo_path, file, NULL);

  g_free (file);
  return full_filename;
}

#ifdef NUMCOSMO_HAVE_FFTW3
/**
 * ncm_cfg_load_fftw_wisdom:
 * @filename: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_cfg_load_fftw_wisdom (gchar *filename, ...)
{
  FILE *wis;
  gchar *file;
  gchar *full_filename;
  va_list ap;
  gboolean ret = FALSE;
  g_assert (numcosmo_init);

  va_start (ap, filename);
  file = g_strdup_vprintf (filename, ap);
  va_end (ap);

  full_filename = g_build_filename (numcosmo_path, file, NULL);

  if (g_file_test (full_filename, G_FILE_TEST_EXISTS))
  {
	wis = g_fopen (full_filename, "r");
	fftw_import_wisdom_from_file(wis);
	fclose (wis);
	ret = TRUE;
  }
  g_free (file);
  g_free (full_filename);
  return ret;
}

/**
 * ncm_cfg_save_fftw_wisdom:
 * @filename: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_cfg_save_fftw_wisdom (gchar *filename, ...)
{
  FILE *wis;
  gchar *file;
  gchar *full_filename;
  va_list ap;

  g_assert (numcosmo_init);

  va_start (ap, filename);
  file = g_strdup_vprintf (filename, ap);
  va_end (ap);
  full_filename = g_build_filename (numcosmo_path, file, NULL);

  wis = g_fopen (full_filename, "w");
  fftw_export_wisdom_to_file(wis);
  fclose (wis);
  g_free (file);
  g_free (full_filename);
  return TRUE;
}
#endif /* NUMCOSMO_HAVE_FFTW3 */

/**
 * ncm_cfg_fopen: (skip)
 * @filename: FIXME
 * @mode: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
FILE *
ncm_cfg_fopen (gchar *filename, gchar *mode, ...)
{
  FILE *F;
  gchar *file;
  gchar *full_filename;
  va_list ap;

  g_assert (numcosmo_init);

  va_start (ap, mode);
  file = g_strdup_vprintf (filename, ap);
  va_end (ap);
  full_filename = g_build_filename (numcosmo_path, file, NULL);

  F = g_fopen (full_filename, mode);
  g_free (file);
  g_free (full_filename);
  return F;
}

/**
 * ncm_cfg_vfopen: (skip)
 * @filename: FIXME
 * @mode: FIXME
 * @ap: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
FILE *
ncm_cfg_vfopen (gchar *filename, gchar *mode, va_list ap)
{
  FILE *F;
  gchar *file;
  gchar *full_filename;

  g_assert (numcosmo_init);
  file = g_strdup_vprintf (filename, ap);
  full_filename = g_build_filename (numcosmo_path, file, NULL);

  F = g_fopen (full_filename, mode);
  g_free (file);
  g_free (full_filename);
  return F;
}

/**
 * ncm_cfg_exists:
 * @filename: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_cfg_exists (gchar *filename, ...)
{
  gboolean exists;
  gchar *file;
  gchar *full_filename;
  va_list ap;

  g_assert (numcosmo_init);

  va_start (ap, filename);
  file = g_strdup_vprintf (filename, ap);
  va_end (ap);

  full_filename = g_build_filename (numcosmo_path, file, NULL);
  exists = g_file_test (full_filename, G_FILE_TEST_EXISTS);
  g_free (file);
  g_free (full_filename);
  return exists;
}

/**
 * ncm_cfg_load_spline: (skip)
 * @filename: FIXME
 * @stype: FIXME
 * @s: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_cfg_load_spline (gchar *filename, const gsl_interp_type *stype, NcmSpline **s, ...)
{
  guint64 size;
  NcmVector *xv, *yv;
  va_list ap;

  va_start (ap, s);
  FILE *F = ncm_cfg_vfopen (filename, "r", ap);
  va_end (ap);
  if (F == NULL)
	return FALSE;

  if (fread (&size, sizeof (guint64), 1, F) != 1)
	g_error ("ncm_cfg_save_spline: fwrite error");

  if (*s == NULL)
  {
	xv = ncm_vector_new (size);
	yv = ncm_vector_new (size);
  }
  else
  {
	xv = (*s)->xv;
	yv = (*s)->yv;
	g_assert (size == ncm_vector_len(xv));
  }

  if (fread (ncm_vector_ptr (xv, 0), sizeof (gdouble), size, F) != size)
	g_error ("ncm_cfg_save_spline: fwrite error");
  if (fread (ncm_vector_ptr (yv, 0), sizeof (gdouble), size, F) != size)
	g_error ("ncm_cfg_save_spline: fwrite error");

  if (*s == NULL)
  {
	*s = ncm_spline_gsl_new (stype);
	ncm_spline_set (*s, xv, yv, TRUE);
  }
  else
	ncm_spline_prepare (*s);

  fclose (F);
  return TRUE;
}

/**
 * ncm_cfg_save_spline:
 * @filename: FIXME
 * @s: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_cfg_save_spline (gchar *filename, NcmSpline *s, ...)
{
  guint64 size;
  va_list ap;

  va_start (ap, s);
  FILE *F = ncm_cfg_vfopen (filename, "w", ap);
  va_end (ap);

  if (F == NULL)
	return FALSE;

  size = s->len;

  if (fwrite (&size, sizeof (guint64), 1, F) != 1)
	g_error ("ncm_cfg_save_spline: fwrite error");

  if (fwrite (ncm_vector_ptr (s->xv, 0), sizeof (gdouble), size, F) != size)
	g_error ("ncm_cfg_save_spline: fwrite error");
  if (fwrite (ncm_vector_ptr (s->yv, 0), sizeof (gdouble), size, F) != size)
	g_error ("ncm_cfg_save_spline: fwrite error");

  fclose (F);
  return TRUE;
}

/**
 * ncm_cfg_load_vector: (skip)
 * @filename: FIXME
 * @v: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_cfg_load_vector (gchar *filename, gsl_vector *v, ...)
{
  va_list ap;

  va_start (ap, v);
  FILE *F = ncm_cfg_vfopen (filename, "r", ap);
  va_end (ap);

  if (F == NULL)
	return FALSE;
  gsl_vector_fread (F, v);
  fclose (F);
  return TRUE;
}

/**
 * ncm_cfg_save_vector: (skip)
 * @filename: FIXME
 * @v: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_cfg_save_vector (gchar *filename, gsl_vector *v, ...)
{
  va_list ap;

  va_start(ap, v);
  FILE *F = ncm_cfg_vfopen (filename, "w", ap);
  va_end (ap);

  if (F == NULL)
	return FALSE;
  gsl_vector_fwrite (F, v);
  fclose (F);
  return TRUE;
}

/**
 * ncm_cfg_load_matrix: (skip)
 * @filename: FIXME
 * @M: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_cfg_load_matrix (gchar *filename, gsl_matrix *M, ...)
{
  va_list ap;

  va_start (ap, M);
  FILE *F = ncm_cfg_vfopen (filename, "r", ap);
  va_end (ap);

  if (F == NULL)
	return FALSE;
  gsl_matrix_fread (F, M);
  fclose (F);
  return TRUE;
}

/**
 * ncm_cfg_save_matrix: (skip)
 * @filename: FIXME
 * @M: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_cfg_save_matrix (gchar *filename, gsl_matrix *M, ...)
{
  va_list ap;

  va_start (ap, M);
  FILE *F = ncm_cfg_vfopen (filename, "w", ap);
  va_end (ap);

  if (F == NULL)
	return FALSE;
  gsl_matrix_fwrite (F, M);
  fclose (F);
  return TRUE;
}

#define load_save_vector_matrix(typen) \
gboolean \
ncm_cfg_load_vector_##typen (gchar *filename, gsl_vector_##typen *v, ...) \
{ \
va_list ap; \
va_start (ap, v); \
FILE *F = ncm_cfg_vfopen (filename, "r", ap); \
va_end (ap); \
if (F == NULL) \
return FALSE; \
gsl_vector_##typen##_fread (F, v); \
fclose (F); \
return TRUE; \
} \
gboolean \
ncm_cfg_save_vector_##typen (gchar *filename, gsl_vector_##typen *v, ...) \
{ \
va_list ap; \
va_start(ap, v); \
FILE *F = ncm_cfg_vfopen (filename, "w", ap); \
va_end (ap); \
if (F == NULL) \
return FALSE; \
gsl_vector_##typen##_fwrite (F, v); \
fclose (F); \
return TRUE; \
} \
gboolean \
ncm_cfg_load_matrix_##typen (gchar *filename, gsl_matrix_##typen *M, ...) \
{ \
va_list ap; \
va_start(ap, M); \
FILE *F = ncm_cfg_vfopen (filename, "r", ap); \
va_end (ap); \
if (F == NULL) \
return FALSE; \
gsl_matrix_##typen##_fread (F, M); \
fclose (F); \
return TRUE; \
} \
gboolean \
ncm_cfg_save_matrix_##typen (gchar *filename, gsl_matrix_##typen *M, ...) \
{ \
va_list ap; \
va_start(ap, M); \
FILE *F = ncm_cfg_vfopen (filename, "w", ap); \
va_end (ap); \
if (F == NULL) \
return FALSE; \
gsl_matrix_##typen##_fwrite (F, M); \
fclose (F); \
return TRUE; \
} \

/**
 * ncm_cfg_load_vector_int: (skip)
 * @filename: FIXME
 * @v: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_cfg_save_vector_int: (skip)
 * @filename: FIXME
 * @v: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_cfg_load_matrix_int: (skip)
 * @filename: FIXME
 * @M: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_cfg_save_matrix_int: (skip)
 * @filename: FIXME
 * @M: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
load_save_vector_matrix(int)

/**
   * ncm_cfg_load_vector_float: (skip)
   * @filename: FIXME
   * @v: FIXME
   * @...: FIXME
   *
   * FIXME
   *
   * Returns: FIXME
   */
/**
   * ncm_cfg_save_vector_float: (skip)
   * @filename: FIXME
   * @v: FIXME
   * @...: FIXME
   *
   * FIXME
   *
   * Returns: FIXME
   */
/**
   * ncm_cfg_load_matrix_float: (skip)
   * @filename: FIXME
   * @M: FIXME
   * @...: FIXME
   *
   * FIXME
   *
   * Returns: FIXME
   */
/**
   * ncm_cfg_save_matrix_float: (skip)
   * @filename: FIXME
   * @M: FIXME
   * @...: FIXME
   *
   * FIXME
   *
   * Returns: FIXME
   */
load_save_vector_matrix(float)

/**
   * ncm_cfg_load_vector_complex: (skip)
   * @filename: FIXME
   * @v: FIXME
   * @...: FIXME
   *
   * FIXME
   *
   * Returns: FIXME
   */
/**
   * ncm_cfg_save_vector_complex: (skip)
   * @filename: FIXME
   * @v: FIXME
   * @...: FIXME
   *
   * FIXME
   *
   * Returns: FIXME
   */
/**
   * ncm_cfg_load_matrix_complex: (skip)
   * @filename: FIXME
   * @M: FIXME
   * @...: FIXME
   *
   * FIXME
   *
   * Returns: FIXME
   */
/**
   * ncm_cfg_save_matrix_complex: (skip)
   * @filename: FIXME
   * @M: FIXME
   * @...: FIXME
   *
   * FIXME
   *
   * Returns: FIXME
   */
load_save_vector_matrix(complex)

#ifdef NUMCOSMO_HAVE_SQLITE3
/**
   * ncm_cfg_get_default_sqlite3: (skip)
   *
   * FIXME
   *
   * Returns: FIXME
   */
sqlite3 *
ncm_cfg_get_default_sqlite3 (void)
{
  sqlite3 *db;
  gchar *filename = g_build_filename (PACKAGE_DATA_DIR, NC_CFG_DEFAULT_SQLITE3_FILENAME, NULL);
  gint ret;

  if (!g_file_test (filename, G_FILE_TEST_EXISTS))
	g_error ("Default database not found (%s)", filename);
  ret = sqlite3_open(filename, &db);
  if (ret  != SQLITE_OK)
	g_error ("Connection to database failed: %s", sqlite3_errmsg (db));

  return db;
}
#endif

/**
 * ncm_cfg_command_line:
 * @argv: FIXME
 * @argc: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gchar *
ncm_cfg_command_line (gchar *argv[], gint argc)
{
  gchar *full_cmd_line;
  gchar *full_cmd_line_ptr;
  guint tsize = (argc - 1) + 1;
  gint argv_size, i;

  for (i = 0; i < argc; i++)
  {
	tsize += strlen (argv[i]);
	if (g_strrstr (argv[i], " ") != NULL)
	  tsize += 2;
  }

  full_cmd_line_ptr = full_cmd_line = g_new (gchar, tsize);

  argv_size = strlen(argv[0]);
  memcpy (full_cmd_line_ptr, argv[0], argv_size);
  full_cmd_line_ptr = &full_cmd_line_ptr[argv_size];

  for (i = 1; i < argc; i++)
  {
	gboolean has_space = FALSE;
	full_cmd_line_ptr[0] = ' ';
	full_cmd_line_ptr++;
	argv_size = strlen(argv[i]);
	has_space = (g_strrstr (argv[i], " ") != NULL);
	if (has_space)
	  (full_cmd_line_ptr++)[0] = '\'';
	memcpy (full_cmd_line_ptr, argv[i], argv_size);
	full_cmd_line_ptr = &full_cmd_line_ptr[argv_size];
	if (has_space)
	  (full_cmd_line_ptr++)[0] = '\'';
  }

  full_cmd_line_ptr[0] = '\0';

  return full_cmd_line;
}
