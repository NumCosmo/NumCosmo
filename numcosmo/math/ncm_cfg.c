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
#include "build_cfg.h"

#include "math/ncm_cfg.h"
#include "math/ncm_spline_gsl.h"
#include "math/ncm_spline_cubic.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_spline2d_bicubic.h"
#include "math/ncm_spline2d_gsl.h"
#include "math/ncm_spline2d_spline.h"
#include "math/ncm_model.h"
#include "math/ncm_model_ctrl.h"
#include "math/ncm_reparam_linear.h"
#include "model/nc_hicosmo_qconst.h"
#include "model/nc_hicosmo_qlinear.h"
#include "model/nc_hicosmo_qpw.h"
#include "model/nc_hicosmo_qspline.h"
#include "model/nc_hicosmo_lcdm.h"
#include "model/nc_hicosmo_de_xcdm.h"
#include "model/nc_hicosmo_de_linder.h"
#include "model/nc_hicosmo_de_pad.h"
#include "model/nc_hicosmo_de_qe.h"
#include "lss/nc_window_tophat.h"
#include "lss/nc_window_gaussian.h"
#include "lss/nc_growth_func.h"
#include "lss/nc_transfer_func.h"
#include "lss/nc_transfer_func_bbks.h"
#include "lss/nc_transfer_func_eh.h"
#include "lss/nc_transfer_func_camb.h"
#include "lss/nc_transfer_func_pert.h"
#include "lss/nc_matter_var.h"
#include "lss/nc_multiplicity_func.h"
#include "lss/nc_multiplicity_func_st.h"
#include "lss/nc_multiplicity_func_ps.h"
#include "lss/nc_multiplicity_func_jenkins.h"
#include "lss/nc_multiplicity_func_warren.h"
#include "lss/nc_multiplicity_func_tinker.h"
#include "lss/nc_multiplicity_func_tinker_mean.h"
#include "lss/nc_multiplicity_func_tinker_crit.h"
#include "lss/nc_mass_function.h"
#include "lss/nc_galaxy_acf.h"
#include "lss/nc_cluster_mass.h"
#include "lss/nc_cluster_mass_nodist.h"
#include "lss/nc_cluster_mass_lnnormal.h"
#include "lss/nc_cluster_mass_vanderlinde.h"
#include "lss/nc_cluster_mass_benson.h"
#include "lss/nc_cluster_mass_benson_xray.h"
#include "lss/nc_cluster_redshift.h"
#include "lss/nc_cluster_redshift_nodist.h"
#include "lss/nc_cluster_photoz_gauss_global.h"
#include "lss/nc_cluster_photoz_gauss.h"
#include "lss/nc_halo_bias_func.h"
#include "lss/nc_halo_bias_type_ps.h"
#include "lss/nc_halo_bias_type_st_ellip.h"
#include "lss/nc_halo_bias_type_st_spher.h"
#include "lss/nc_halo_bias_type_tinker.h"
#include "lss/nc_cluster_abundance.h"

#include <gio/gio.h>
#ifdef NUMCOSMO_HAVE_FFTW3
#include <fftw3.h>
#endif /* NUMCOSMO_HAVE_FFTW3 */

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
	fprintf (_log_stream, "%s", message);
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

  g_setenv ("CUBACORES", "0", TRUE);

  _log_msg_id = g_log_set_handler (G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE | G_LOG_LEVEL_DEBUG, _ncm_cfg_log_message, NULL);
  _log_stream = stdout;

  ncm_cfg_register_model (NCM_TYPE_SPLINE);
  ncm_cfg_register_model (NCM_TYPE_SPLINE_CUBIC);
  ncm_cfg_register_model (NCM_TYPE_SPLINE_CUBIC_NOTAKNOT);
  ncm_cfg_register_model (NCM_TYPE_SPLINE_GSL);

  ncm_cfg_register_model (NCM_TYPE_SPLINE2D);
  ncm_cfg_register_model (NCM_TYPE_SPLINE2D_BICUBIC);
  ncm_cfg_register_model (NCM_TYPE_SPLINE2D_GSL);
  ncm_cfg_register_model (NCM_TYPE_SPLINE2D_SPLINE);

  ncm_cfg_register_model (NCM_TYPE_MODEL);
  ncm_cfg_register_model (NCM_TYPE_MODEL_CTRL);

  ncm_cfg_register_model (NCM_TYPE_REPARAM);
  ncm_cfg_register_model (NCM_TYPE_REPARAM_LINEAR);

  ncm_cfg_register_model (NC_TYPE_HICOSMO_QCONST);
  ncm_cfg_register_model (NC_TYPE_HICOSMO_QLINEAR);
  ncm_cfg_register_model (NC_TYPE_HICOSMO_QPW);
  ncm_cfg_register_model (NC_TYPE_HICOSMO_QSPLINE);
  ncm_cfg_register_model (NC_TYPE_HICOSMO_LCDM);
  ncm_cfg_register_model (NC_TYPE_HICOSMO_DE_XCDM);
  ncm_cfg_register_model (NC_TYPE_HICOSMO_DE_LINDER);
  ncm_cfg_register_model (NC_TYPE_HICOSMO_DE_PAD);
  ncm_cfg_register_model (NC_TYPE_HICOSMO_DE_QE);

  ncm_cfg_register_model (NC_TYPE_WINDOW);
  ncm_cfg_register_model (NC_TYPE_WINDOW_TOPHAT);
  ncm_cfg_register_model (NC_TYPE_WINDOW_GAUSSIAN);

  ncm_cfg_register_model (NC_TYPE_GROWTH_FUNC);

  ncm_cfg_register_model (NC_TYPE_TRANSFER_FUNC);
  ncm_cfg_register_model (NC_TYPE_TRANSFER_FUNC_BBKS);
  ncm_cfg_register_model (NC_TYPE_TRANSFER_FUNC_EH);
  ncm_cfg_register_model (NC_TYPE_TRANSFER_FUNC_CAMB);
  ncm_cfg_register_model (NC_TYPE_TRANSFER_FUNC_PERT);

  ncm_cfg_register_model (NC_TYPE_MATTER_VAR);

  ncm_cfg_register_model (NC_TYPE_MULTIPLICITY_FUNC);
  ncm_cfg_register_model (NC_TYPE_MULTIPLICITY_FUNC_PS);
  ncm_cfg_register_model (NC_TYPE_MULTIPLICITY_FUNC_ST);
  ncm_cfg_register_model (NC_TYPE_MULTIPLICITY_FUNC_JENKINS);
  ncm_cfg_register_model (NC_TYPE_MULTIPLICITY_FUNC_WARREN);
  ncm_cfg_register_model (NC_TYPE_MULTIPLICITY_FUNC_TINKER);
  ncm_cfg_register_model (NC_TYPE_MULTIPLICITY_FUNC_TINKER_MEAN);
  ncm_cfg_register_model (NC_TYPE_MULTIPLICITY_FUNC_TINKER_CRIT);

  ncm_cfg_register_model (NC_TYPE_MASS_FUNCTION);

  ncm_cfg_register_model (NC_TYPE_GALAXY_ACF);

  ncm_cfg_register_model (NC_TYPE_CLUSTER_MASS);
  ncm_cfg_register_model (NC_TYPE_CLUSTER_MASS_NODIST);
  ncm_cfg_register_model (NC_TYPE_CLUSTER_MASS_LNNORMAL);
  ncm_cfg_register_model (NC_TYPE_CLUSTER_MASS_VANDERLINDE);
  ncm_cfg_register_model (NC_TYPE_CLUSTER_MASS_BENSON);
  ncm_cfg_register_model (NC_TYPE_CLUSTER_MASS_BENSON_XRAY);

  ncm_cfg_register_model (NC_TYPE_CLUSTER_REDSHIFT);
  ncm_cfg_register_model (NC_TYPE_CLUSTER_REDSHIFT_NODIST);
  ncm_cfg_register_model (NC_TYPE_CLUSTER_PHOTOZ_GAUSS_GLOBAL);
  ncm_cfg_register_model (NC_TYPE_CLUSTER_PHOTOZ_GAUSS);

  ncm_cfg_register_model (NC_TYPE_HALO_BIAS_FUNC);

  ncm_cfg_register_model (NC_TYPE_HALO_BIAS_TYPE);
  ncm_cfg_register_model (NC_TYPE_HALO_BIAS_TYPE_PS);
  ncm_cfg_register_model (NC_TYPE_HALO_BIAS_TYPE_ST_ELLIP);
  ncm_cfg_register_model (NC_TYPE_HALO_BIAS_TYPE_ST_SPHER);
  ncm_cfg_register_model (NC_TYPE_HALO_BIAS_TYPE_TINKER);

  ncm_cfg_register_model (NC_TYPE_CLUSTER_ABUNDANCE);

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
  GType model_type = NCM_TYPE_MODEL;
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
    g_error ("ncm_cfg_set_logfile: Can't open logfile (%s)", filename);
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

/**
 * ncm_cfg_keyfile_to_arg:
 * @kfile: FIXME
 * @group_name: FIXME
 * @entries: FIXME
 * @argv: FIXME
 * @argc: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
ncm_cfg_keyfile_to_arg (GKeyFile *kfile, gchar *group_name, GOptionEntry *entries, gchar **argv, gint *argc)
{
  if (g_key_file_has_group (kfile, group_name))
  {
	GError *error = NULL;
	gint i;
	for (i = 0; entries[i].long_name != NULL; i++)
	{
	  if (g_key_file_has_key (kfile, group_name, entries[i].long_name, &error))
	  {
		gchar *val = g_key_file_get_value (kfile, group_name, entries[i].long_name, &error);
		if (entries[i].arg == G_OPTION_ARG_NONE)
		{
		  if ((g_ascii_strcasecmp (val, "1") == 0) ||
		      (g_ascii_strcasecmp (val, "true") == 0))
		  {
			argv[argc[0]++] = g_strdup_printf ("--%s", entries[i].long_name);
		  }
		  g_free (val);
		}
		else if (strlen (val) > 0)
		{
		  argv[argc[0]++] = g_strdup_printf ("--%s", entries[i].long_name);
		  argv[argc[0]++] = val;
		}
	  }
	}
  }
}

/**
 * ncm_cfg_keyfile_to_arg:
 * @kfile: FIXME
 * @group_name: FIXME
 * @entries: FIXME
 * @argv: FIXME
 * @argc: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
ncm_cfg_entries_to_keyfile (GKeyFile *kfile, gchar *group_name, GOptionEntry *entries)
{
  GError *error = NULL;
  gint i;
  for (i = 0; entries[i].long_name != NULL; i++)
  {
	gboolean skip_comment = FALSE;
	switch (entries[i].arg)
	{
	  case G_OPTION_ARG_NONE:
	  {
		gboolean arg_b = ((gboolean *)entries[i].arg_data)[0];
		g_key_file_set_boolean (kfile, group_name, entries[i].long_name, arg_b);
		break;
	  }
	  case G_OPTION_ARG_STRING:
	  case G_OPTION_ARG_FILENAME:
	  {
		gchar **arg_s = (gchar **)entries[i].arg_data;
		g_key_file_set_string (kfile, group_name, entries[i].long_name, *arg_s != NULL ? *arg_s : "");
		break;
	  }
	  case G_OPTION_ARG_STRING_ARRAY:
	  case G_OPTION_ARG_FILENAME_ARRAY:
	  {
		const gchar ***arg_as = (const gchar ***)entries[i].arg_data;
		gchar ***arg_cas = (gchar ***)entries[i].arg_data;
		if (*arg_cas != NULL)
		  g_key_file_set_string_list (kfile, group_name, entries[i].long_name, *arg_as, g_strv_length (*arg_cas));
		else
		  g_key_file_set_string_list (kfile, group_name, entries[i].long_name, NULL, 0);
		break;
	  }
	  case G_OPTION_ARG_INT:
	  {
		gint arg_i = ((gint *)entries[i].arg_data)[0];
		g_key_file_set_integer (kfile, group_name, entries[i].long_name, arg_i);
		break;
	  }
	  case G_OPTION_ARG_INT64:
	  {
		gint64 arg_l = ((gint64 *)entries[i].arg_data)[0];
		g_key_file_set_int64 (kfile, group_name, entries[i].long_name, arg_l);
		break;
	  }
	  case G_OPTION_ARG_DOUBLE:
	  {
		gdouble arg_d = ((double *)entries[i].arg_data)[0];
		g_key_file_set_double (kfile, group_name, entries[i].long_name, arg_d);
		break;
	  }
	  default:
		skip_comment = TRUE;
		//g_error ("ncm_cfg_entries_to_keyfile: cannot convert entry type %d to keyfile", entries[i].arg);
		break;
	}
	if (!skip_comment)
	{
	  if (!g_key_file_set_comment (kfile, group_name, entries[i].long_name, entries[i].description, &error))
		g_error ("ncm_cfg_entries_to_keyfile: %s", error->message);
	}
  }
}


/**
 * ncm_cfg_get_enum_by_id_name_nick:
 * @enum_type: FIXME
 * @id_name_nick: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
const GEnumValue *
ncm_cfg_get_enum_by_id_name_nick (GType enum_type, const gchar *id_name_nick)
{
  GEnumValue *res = NULL;
  gchar *endptr = NULL;
  gint64 id = g_ascii_strtoll (id_name_nick, &endptr, 10);
  GEnumClass *enum_class = NULL;

  g_assert (G_TYPE_IS_ENUM (enum_type));

  enum_class = g_type_class_ref (enum_type);
  if ((endptr == id_name_nick) || (strlen (endptr) > 0))
  {
	res = g_enum_get_value_by_name (enum_class, id_name_nick);
	if (res == NULL)
	  res = g_enum_get_value_by_nick (enum_class, id_name_nick);
  }
  else
	res = g_enum_get_value (enum_class, id);

  g_type_class_unref (enum_class);

  return res;
}

/**
 * ncm_cfg_enum_print_all:
 * @enum_type: FIXME
 * @header: FIXME
 *
 * FIXME
 *
 */
void
ncm_cfg_enum_print_all (GType enum_type, gchar *header)
{
  GEnumClass *enum_class;
  GEnumValue *snia;
  gint i = 0;
  gint name_max_len = 4;
  gint nick_max_len = 4;
  gint pad;

  g_assert (G_TYPE_IS_ENUM (enum_type));

  enum_class = g_type_class_ref (enum_type);

  while ((snia = g_enum_get_value (enum_class, i++)) != NULL)
  {
	name_max_len = GSL_MAX (name_max_len, strlen (snia->value_name));
	nick_max_len = GSL_MAX (nick_max_len, strlen (snia->value_nick));
  }

  printf ("# %s:\n", header);
  pad = 10 + name_max_len + nick_max_len;
  printf ("#");
  while (pad-- != 0)
	printf ("-");
  printf ("#\n");
  printf ("# Id | %-*s | %-*s |\n", name_max_len, "Name", nick_max_len, "Nick");
  i = 0;
  while ((snia = g_enum_get_value (enum_class, i++)) != NULL)
  {
	printf ("# %02d | %-*s | %-*s |\n", snia->value, name_max_len, snia->value_name, nick_max_len, snia->value_nick);
  }
  pad = 10 + name_max_len + nick_max_len;
  printf ("#");
  while (pad-- != 0)
	printf ("-");
  printf ("#\n");

  g_type_class_unref (enum_class);
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
  static sqlite3 *db = NULL;
  if (db == NULL)
  {
	gchar *filename = g_build_filename (PACKAGE_DATA_DIR, NC_CFG_DEFAULT_SQLITE3_FILENAME, NULL);
	gint ret;

	if (!g_file_test (filename, G_FILE_TEST_EXISTS))
	  g_error ("Default database not found (%s)", filename);
	ret = sqlite3_open (filename, &db);
	if (ret  != SQLITE_OK)
	  g_error ("Connection to database failed: %s", sqlite3_errmsg (db));

	g_free (filename);
  }

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

/**
 * ncm_cfg_create_from_string:
 * @obj_ser: String containing the serialized version of the object.
 *
 * Parses the serialized and returns the newly created object.
 *
 * Returns: (transfer full): A new #GObject.
 */
GObject *
ncm_cfg_create_from_string (const gchar *obj_ser)
{
  GError *error = NULL;
  static gsize regex_init = FALSE;
  static GRegex *regex = NULL;
  GMatchInfo *match_info = NULL;
  GObject *obj = NULL;
  GVariant *var_obj = NULL;

  if (g_once_init_enter (&regex_init))
  {
	regex = g_regex_new ("^([A-Z][A-Za-z]*)\\s*(.*?)\\s*$", 0, 0, &error);
	g_once_init_leave (&regex_init, TRUE);
  }

  var_obj = g_variant_parse (G_VARIANT_TYPE ("{sa{sv}}"), obj_ser, NULL, NULL, &error);
  g_clear_error (&error);
  if (var_obj != NULL)
  {
	GVariant *obj_name = g_variant_get_child_value (var_obj, 0);
	GVariant *params = g_variant_get_child_value (var_obj, 1);
	obj = ncm_cfg_create_from_name_params (g_variant_get_string (obj_name, NULL), params);
	g_variant_unref (obj_name);
	g_variant_unref (params);
	g_variant_unref (var_obj);
  }
  else if (g_regex_match (regex, obj_ser, 0, &match_info))
  {
	gchar *obj_name = g_match_info_fetch (match_info, 1);
	gchar *obj_prop = g_match_info_fetch (match_info, 2);
	GVariant *params = NULL;

	if (obj_prop != NULL && strlen (obj_prop) > 0)
	{
	  params = g_variant_parse (G_VARIANT_TYPE_VARDICT, obj_prop, NULL, NULL, &error);
	  if (error != NULL)
	  {
		g_error ("Cannot parse object '%s' parameters '%s', error %s\n", obj_name, obj_prop, error->message);
		g_error_free (error);
	  }
	}

	obj = ncm_cfg_create_from_name_params (obj_name, params);

	if (params != NULL)
	  g_variant_unref (params);
	g_free (obj_name);
	g_free (obj_prop);
	g_match_info_free (match_info);
  }
  else
  {
    g_match_info_free (match_info);
	g_error ("Cannot indentify object in string '%s'\n", obj_ser);
  }

  return obj;
}

/**
 * ncm_cfg_create_from_name_params:
 * @obj_name: string containing the object name.
 * @params: a #GVariant containing the object parameters.
 *
 * FIXME
 *
 * Returns: (transfer full): A new #GObject.
 */
GObject *
ncm_cfg_create_from_name_params (const gchar *obj_name, GVariant *params)
{
  GObject *obj = NULL;
  GType gtype = g_type_from_name (obj_name);

  if (gtype == 0)
	g_error ("Object '%s' is not registered\n", obj_name);

  g_assert (params == NULL || g_variant_is_of_type (params, G_VARIANT_TYPE ("a{sv}")));

  if (params != NULL)
  {
	GVariantIter *p_iter = g_variant_iter_new (params);
	gsize nprop = g_variant_iter_n_children (p_iter);
	GParameter *gprop = g_new (GParameter, nprop);
	GVariant *var = NULL;
	guint i = 0;

	while ((var = g_variant_iter_next_value (p_iter)))
	{
	  GVariant *var_key = g_variant_get_child_value (var, 0);
	  GVariant *var_val = g_variant_get_child_value (var, 1);
	  GVariant *val = g_variant_get_variant (var_val);
	  gprop[i].name = g_variant_get_string (var_key, NULL);

	  if (g_variant_is_of_type (val, G_VARIANT_TYPE ("{sa{sv}}")))
	  {
		GVariant *nest_obj_key = g_variant_get_child_value (val, 0);
		GVariant *nest_obj_params = g_variant_get_child_value (val, 1);
		GValue lval = G_VALUE_INIT;
		GObject *nest_obj =
		  ncm_cfg_create_from_name_params (g_variant_get_string (nest_obj_key, NULL),
		                                   nest_obj_params);
		g_value_init (&lval, G_TYPE_OBJECT);
		gprop[i].value = lval;
		g_value_take_object (&gprop[i].value, nest_obj);
		g_variant_unref (nest_obj_key);
		g_variant_unref (nest_obj_params);
	  }
	  else
		g_dbus_gvariant_to_gvalue (val, &gprop[i].value);

	  i++;
	  g_variant_unref (var_key);
	  g_variant_unref (val);
	  g_variant_unref (var_val);
	  g_variant_unref (var);
	}

	obj = g_object_newv (gtype, nprop, gprop);

	g_free (gprop);
	g_variant_iter_free (p_iter);
  }
  else
	obj = g_object_new (gtype, NULL);

  return obj;
}

static const GVariantType *
_ncm_cfg_gtype_to_gvariant_type (GType t)
{
  switch (t)
  {
	case G_TYPE_CHAR:
	case G_TYPE_UCHAR:
	  return G_VARIANT_TYPE_BYTE;
	  break;
	case G_TYPE_BOOLEAN:
	  return G_VARIANT_TYPE_BOOLEAN;
	  break;
	case G_TYPE_INT:
	{
	  switch (sizeof (gint))
	  {
		case 2:
		  return G_VARIANT_TYPE_INT16;
		  break;
		case 4:
		  return G_VARIANT_TYPE_INT32;
		  break;
		case 8:
		  return G_VARIANT_TYPE_INT64;
		  break;
		default:
		  g_error ("Unknown gint size %"G_GSIZE_FORMAT"\n", sizeof(gint));
		  break;
	  }
	  break;
	}
	case G_TYPE_UINT:
	{
	  switch (sizeof (guint))
	  {
		case 2:
		  return G_VARIANT_TYPE_UINT16;
		  break;
		case 4:
		  return G_VARIANT_TYPE_UINT32;
		  break;
		case 8:
		  return G_VARIANT_TYPE_UINT64;
		  break;
		default:
		  g_error ("Unknown gint size %"G_GSIZE_FORMAT"\n", sizeof(guint));
		  break;
	  }
	  break;
	}
	case G_TYPE_LONG:
	{
	  switch (sizeof (glong))
	  {
		case 2:
		  return G_VARIANT_TYPE_INT16;
		  break;
		case 4:
		  return G_VARIANT_TYPE_INT32;
		  break;
		case 8:
		  return G_VARIANT_TYPE_INT64;
		  break;
		default:
		  g_error ("Unknown gint size %"G_GSIZE_FORMAT"\n", sizeof(glong));
		  break;
	  }
	  break;
	}
	case G_TYPE_ULONG:
	{
	  switch (sizeof (gulong))
	  {
		case 2:
		  return G_VARIANT_TYPE_UINT16;
		  break;
		case 4:
		  return G_VARIANT_TYPE_UINT32;
		  break;
		case 8:
		  return G_VARIANT_TYPE_UINT64;
		  break;
		default:
		  g_error ("Unknown gint size %"G_GSIZE_FORMAT"\n", sizeof(gulong));
		  break;
	  }
	  break;
	}
	case G_TYPE_INT64:
	  return G_VARIANT_TYPE_INT64;
	  break;
	case G_TYPE_UINT64:
	  return G_VARIANT_TYPE_UINT64;
	  break;
	case G_TYPE_FLOAT:
	case G_TYPE_DOUBLE:
	  return G_VARIANT_TYPE_DOUBLE;
	  break;
	case G_TYPE_STRING:
	  return G_VARIANT_TYPE_STRING;
	  break;
	case G_TYPE_VARIANT:
	  return G_VARIANT_TYPE_VARIANT;
	  break;
	default:
	  return NULL;
	  break;
  }
}

/**
 * ncm_cfg_gvalue_to_gvariant:
 * @val: a #GValue.
 *
 * FIXME
 *
 * Returns: (transfer full): A #GVariant convertion of @val.
 */
GVariant *
ncm_cfg_gvalue_to_gvariant (GValue *val)
{
  GType t = G_VALUE_TYPE (val);
  GType fund_t = G_TYPE_FUNDAMENTAL (t);
  const GVariantType *var_type = _ncm_cfg_gtype_to_gvariant_type (fund_t);
  GVariant *var = NULL;

  if (var_type == NULL)
  {
	switch (fund_t)
	{
	  case G_TYPE_OBJECT:
	  {
		GObject *nest_obj = g_value_get_object (val);
		if (nest_obj != NULL)
		  var = ncm_cfg_serialize_to_variant (nest_obj);
		break;
	  }
	  case G_TYPE_ENUM:
		var = g_variant_ref_sink (g_variant_new_int32 (g_value_get_enum (val)));
		break;
	  case G_TYPE_FLAGS:
		var = g_variant_ref_sink (g_variant_new_uint32 (g_value_get_flags (val)));
		break;
	  default:
		g_error ("Cannot covert GValue '%s' to GVariant.", g_type_name (t));
		break;
	}
  }
  else
	var = g_dbus_gvalue_to_gvariant (val, var_type);

  return var;
}

/**
 * ncm_cfg_serialize_to_variant:
 * @obj: a #GObject.
 *
 * FIXME
 *
 * Returns: (transfer full): A #GVariant dictionary describing the @obj.
 */
GVariant *
ncm_cfg_serialize_to_variant (GObject *obj)
{
  const gchar *obj_name = G_OBJECT_TYPE_NAME (obj);
  GObjectClass *klass = G_OBJECT_GET_CLASS (obj);
  guint n_properties, i;
  GParamSpec **prop = g_object_class_list_properties (klass, &n_properties);

  if (n_properties == 0)
  {
	g_free (prop);
	return g_variant_ref_sink (g_variant_new ("{sa{sv}}", obj_name, NULL));
  }
  else
  {
	GVariantBuilder *b = g_variant_builder_new (G_VARIANT_TYPE ("a{sv}"));
	GVariant *params, *ser_var;

	for (i = 0; i < n_properties; i++)
	{
	  GVariant *var = NULL;
	  GValue val = G_VALUE_INIT;

	  if (!(prop[i]->flags & G_PARAM_WRITABLE))
		continue;

	  g_value_init (&val, prop[i]->value_type);
	  g_object_get_property (obj, prop[i]->name, &val);

	  var = ncm_cfg_gvalue_to_gvariant (&val);
	  if (var == NULL)
	  {
		g_value_unset (&val);
		continue;
	  }

	  g_variant_builder_add (b, "{sv}", prop[i]->name, var);

	  g_variant_unref (var);
	  g_value_unset (&val);
	}

    params = g_variant_builder_end (b);
	g_variant_builder_unref (b);
	g_free (prop);

	ser_var = g_variant_ref_sink (g_variant_new ("{s@a{sv}}", obj_name, params));
	return ser_var;
  }
}

/**
 * ncm_cfg_serialize_to_string:
 * @obj: a #GObject.
 * @valid_variant: FIXME.
 *
 * FIXME
 *
 *
 * Returns: (transfer full): A string containing the serialized version of @obj.
 */
gchar *
ncm_cfg_serialize_to_string (GObject *obj, gboolean valid_variant)
{
  GVariant *ser_var = ncm_cfg_serialize_to_variant (obj);
  gchar *ser = NULL;

  if (valid_variant)
	ser = g_variant_print (ser_var, TRUE);
  else
  {
	GVariant *params = NULL;
	gchar *obj_name = NULL;
	gchar *params_str;
	g_variant_get (ser_var, "{s@a{sv}}", &obj_name, &params);
    if (g_variant_n_children (params) == 0)
	{
      ser = obj_name;
	}
	else
	{
	  params_str = g_variant_print (params, TRUE);
	  ser = g_strdup_printf ("%s%s", obj_name, params_str);
	  g_free (params_str);
	  g_free (obj_name);
	}
	g_variant_unref (params);
  }
  g_variant_unref (ser_var);
  return ser;
}
