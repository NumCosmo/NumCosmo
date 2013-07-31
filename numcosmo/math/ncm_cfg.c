/***************************************************************************
 *            ncm_cfg.c
 *
 *  Wed Aug 13 20:59:22 2008
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

/**
 * SECTION:ncm_cfg
 * @title: Library Configuration
 * @short_description: Library configuration and helper functions
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_cfg.h"
#include "math/ncm_rng.h"
#include "math/ncm_spline_gsl.h"
#include "math/ncm_spline_cubic.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_spline2d_bicubic.h"
#include "math/ncm_spline2d_gsl.h"
#include "math/ncm_spline2d_spline.h"
#include "math/ncm_model.h"
#include "math/ncm_model_ctrl.h"
#include "math/ncm_reparam_linear.h"
#include "math/ncm_stats_vec.h"
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
#include "nc_distance.h"
#include "nc_recomb.h"
#include "nc_recomb_seager.h"
#include "nc_snia_dist_cov.h"

#include <gio/gio.h>
#ifdef NUMCOSMO_HAVE_FFTW3
#include <fftw3.h>
#endif /* NUMCOSMO_HAVE_FFTW3 */

#ifndef G_VALUE_INIT
#define G_VALUE_INIT {0}
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
#if (GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 32)
  g_thread_init (NULL);
#endif

#if (GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 36)
  g_type_init();
#endif

  g_setenv ("CUBACORES", "0", TRUE);

  _log_msg_id = g_log_set_handler (G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE | G_LOG_LEVEL_DEBUG, _ncm_cfg_log_message, NULL);
  _log_stream = stdout;

  ncm_cfg_register_obj (NCM_TYPE_RNG);
  
  ncm_cfg_register_obj (NCM_TYPE_VECTOR);
  ncm_cfg_register_obj (NCM_TYPE_MATRIX);
  
  ncm_cfg_register_obj (NCM_TYPE_SPLINE);
  ncm_cfg_register_obj (NCM_TYPE_SPLINE_CUBIC);
  ncm_cfg_register_obj (NCM_TYPE_SPLINE_CUBIC_NOTAKNOT);
  ncm_cfg_register_obj (NCM_TYPE_SPLINE_GSL);

  ncm_cfg_register_obj (NCM_TYPE_SPLINE2D);
  ncm_cfg_register_obj (NCM_TYPE_SPLINE2D_BICUBIC);
  ncm_cfg_register_obj (NCM_TYPE_SPLINE2D_GSL);
  ncm_cfg_register_obj (NCM_TYPE_SPLINE2D_SPLINE);

  ncm_cfg_register_obj (NCM_TYPE_MODEL);
  ncm_cfg_register_obj (NCM_TYPE_MODEL_CTRL);

  ncm_cfg_register_obj (NCM_TYPE_REPARAM);
  ncm_cfg_register_obj (NCM_TYPE_REPARAM_LINEAR);

  ncm_cfg_register_obj (NCM_TYPE_BOOTSTRAP);
  ncm_cfg_register_obj (NCM_TYPE_STATS_VEC);

  ncm_cfg_register_obj (NC_TYPE_HICOSMO_QCONST);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_QLINEAR);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_QPW);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_QSPLINE);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_LCDM);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_DE_XCDM);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_DE_LINDER);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_DE_PAD);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_DE_QE);

  ncm_cfg_register_obj (NC_TYPE_WINDOW);
  ncm_cfg_register_obj (NC_TYPE_WINDOW_TOPHAT);
  ncm_cfg_register_obj (NC_TYPE_WINDOW_GAUSSIAN);

  ncm_cfg_register_obj (NC_TYPE_GROWTH_FUNC);

  ncm_cfg_register_obj (NC_TYPE_TRANSFER_FUNC);
  ncm_cfg_register_obj (NC_TYPE_TRANSFER_FUNC_BBKS);
  ncm_cfg_register_obj (NC_TYPE_TRANSFER_FUNC_EH);
  ncm_cfg_register_obj (NC_TYPE_TRANSFER_FUNC_CAMB);
  ncm_cfg_register_obj (NC_TYPE_TRANSFER_FUNC_PERT);

  ncm_cfg_register_obj (NC_TYPE_MATTER_VAR);

  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_PS);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_ST);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_JENKINS);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_WARREN);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_TINKER);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_TINKER_MEAN);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_TINKER_CRIT);

  ncm_cfg_register_obj (NC_TYPE_MASS_FUNCTION);

  ncm_cfg_register_obj (NC_TYPE_GALAXY_ACF);

  ncm_cfg_register_obj (NC_TYPE_CLUSTER_MASS);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_MASS_NODIST);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_MASS_LNNORMAL);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_MASS_VANDERLINDE);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_MASS_BENSON);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_MASS_BENSON_XRAY);

  ncm_cfg_register_obj (NC_TYPE_CLUSTER_REDSHIFT);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_REDSHIFT_NODIST);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_PHOTOZ_GAUSS_GLOBAL);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_PHOTOZ_GAUSS);

  ncm_cfg_register_obj (NC_TYPE_HALO_BIAS_FUNC);

  ncm_cfg_register_obj (NC_TYPE_HALO_BIAS_TYPE);
  ncm_cfg_register_obj (NC_TYPE_HALO_BIAS_TYPE_PS);
  ncm_cfg_register_obj (NC_TYPE_HALO_BIAS_TYPE_ST_ELLIP);
  ncm_cfg_register_obj (NC_TYPE_HALO_BIAS_TYPE_ST_SPHER);
  ncm_cfg_register_obj (NC_TYPE_HALO_BIAS_TYPE_TINKER);

  ncm_cfg_register_obj (NC_TYPE_CLUSTER_ABUNDANCE);

  ncm_cfg_register_obj (NC_TYPE_DISTANCE);

  ncm_cfg_register_obj (NC_TYPE_RECOMB);
  ncm_cfg_register_obj (NC_TYPE_RECOMB_SEAGER);

  ncm_cfg_register_obj (NC_TYPE_SNIA_DIST_COV);
  
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
 * ncm_cfg_register_obj:
 * @obj: FIXME
 *
 * FIXME
 */
void
ncm_cfg_register_obj (GType obj)
{
  gpointer obj_class = g_type_class_ref (obj);
  g_type_class_unref (obj_class);
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
ncm_message (gchar *msg, ...)
{
  va_list ap;
  va_start(ap, msg);
  g_logv (G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE, msg, ap);
  va_end (ap);
}

/**
 * ncm_string_ww:
 * @msg: FIXME
 * @first: FIXME
 * @rest: FIXME
 * @ncols: FIXME
 *
 * FIXME
 * 
 * Returns: (transfer full): word wraped string @msg.
 */
gchar *
ncm_string_ww (const gchar *msg, const gchar *first, const gchar *rest, guint ncols)
{
  gchar **msg_split = g_strsplit (msg, " ", 0);
  guint size_print = strlen (msg_split[0]);
  guint first_size = strlen (first);
  guint rest_size = strlen (rest);
  guint msg_len = strlen (msg);
  guint msg_final_len = msg_len + first_size + msg_len / ncols * rest_size;
  GString *msg_ww = g_string_sized_new (msg_final_len);
  guint i = 1;

  g_string_append_printf (msg_ww, "%s%s", first, msg_split[0]);
  while (msg_split[i] != NULL)
  {
    guint lsize = strlen (msg_split[i]);
    if (size_print + lsize + 1 > (ncols - first_size))
      break;
    g_string_append_printf (msg_ww, " %s", msg_split[i]);
    size_print += lsize + 1;
    i++;
  }
  g_string_append_printf (msg_ww, "\n");

  while (msg_split[i] != NULL)
  {
    size_print = 0;
    g_string_append_printf (msg_ww, "%s", rest);
    while (msg_split[i] != NULL)
    {
      guint lsize = strlen (msg_split[i]);
      if (size_print + lsize + 1 > (ncols - rest_size))
        break;
      g_string_append_printf (msg_ww, " %s", msg_split[i]);
      size_print += lsize + 1;
      i++;
    }
    g_string_append_printf (msg_ww, "\n");
  }

  g_strfreev (msg_split);

  return g_string_free (msg_ww, FALSE);
}

/**
 * ncm_message_ww:
 * @msg: FIXME
 * @first: FIXME
 * @rest: FIXME
 * @ncols: FIXME
 *
 * FIXME
 */
void 
ncm_message_ww (const gchar *msg, const gchar *first, const gchar *rest, guint ncols)
{
  gchar *msg_ww = ncm_string_ww (msg, first, rest, ncols);
  g_message ("%s", msg_ww);
  g_free (msg_ww);
}

/**
 * ncm_cfg_msg_sepa:
 *
 * Log a message separator.
 * 
 */
void 
ncm_cfg_msg_sepa (void)
{
  g_message ("#----------------------------------------------------------------------------------\n");
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
ncm_cfg_get_fullpath (const gchar *filename, ...)
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
        if (entries[i].arg == G_OPTION_ARG_STRING_ARRAY || entries[i].arg == G_OPTION_ARG_FILENAME_ARRAY)
        {
          gint j;
          gsize length;
          gchar **vals = g_key_file_get_string_list (kfile, group_name, entries[i].long_name, &length, &error);
          if (error != NULL)
            g_error ("ncm_cfg_keyfile_to_arg: Cannot parse key file[%s]", error->message);
          for (j = 0; j < length; j++)
          {
            argv[argc[0]++] = g_strdup_printf ("--%s", entries[i].long_name);
            argv[argc[0]++] = vals[j];
          }
          g_free (vals);
        }
        else
        {
          gchar *val = g_key_file_get_value (kfile, group_name, entries[i].long_name, &error);
          if (error != NULL)
            g_error ("ncm_cfg_keyfile_to_arg: Cannot parse key file[%s]", error->message);
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
}

/**
 * ncm_cfg_string_to_comment:
 * @str: FIXME
 * 
 * FIXME
 * 
 * Returns: (transfer full): FIXME
 */
gchar *
ncm_cfg_string_to_comment (const gchar *str)
{
  g_assert (str != NULL);
  {
    gchar *desc_ww = ncm_string_ww (str, "  ", "  ", 80);
    gchar *desc = g_strdup_printf ("###############################################################################\n\n%s\n", desc_ww);
    g_free (desc_ww);
    return desc;
  }
}

/**
 * ncm_cfg_entries_to_keyfile:
 * @kfile: FIXME
 * @group_name: FIXME
 * @entries: FIXME
 *
 * FIXME
 *
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
#if !((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 26))
        gint64 arg_l = ((gint64 *)entries[i].arg_data)[0];
        g_key_file_set_int64 (kfile, group_name, entries[i].long_name, arg_l);
#else
        g_error ("ncm_cfg_entries_to_keyfile: gint64 not supported, recompile "PACKAGE_NAME" with glib >= 2.26");
#endif
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
      gchar *desc = ncm_cfg_string_to_comment (entries[i].description); 
      if (!g_key_file_set_comment (kfile, group_name, entries[i].long_name, desc, &error))
        g_error ("ncm_cfg_entries_to_keyfile: %s", error->message);
      g_free (desc);
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
  g_assert (id_name_nick != NULL);
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
    xv = ncm_spline_get_xv (*s);
    yv = ncm_spline_get_yv (*s);
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
    
  ncm_vector_free (xv);
  ncm_vector_free (yv);

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
    gchar *filename = g_build_filename (PACKAGE_DATA_DIR "/data", NCM_CFG_DEFAULT_SQLITE3_FILENAME, NULL);
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
