/***************************************************************************
 *            nc_data_planck_lkl.c
 *
 *  Tue October 20 15:25:14 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_data_planck_lkl.c
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
 * SECTION:nc_data_planck_lkl
 * @title: NcDataPlankLKL
 * @short_description: Planck Likelihood interface.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_planck_lkl.h"
#include "plc/clik.h"

enum
{
  PROP_0,
  PROP_DATA_FILE,
  PROP_IS_LENSING,
  PROP_NPARAMS,
  PROP_CHKSUM,
  PROP_SIZE
};

G_DEFINE_TYPE (NcDataPlanckLKL, nc_data_planck_lkl, NCM_TYPE_DATA);

void _nc_data_planck_set_filename (NcDataPlanckLKL *plik, const gchar *filename);
#define CLIK_OBJ(obj) ((clik_object *)(obj))
#define CLIK_LENS_OBJ(obj) ((clik_lensing_object *)(obj))

#define CLIK_CHECK_ERROR(str,err) \
G_STMT_START { \
    if (isError (err)) \
    { \
      gchar error_msg[4096]; \
      stringError (error_msg, (err)); \
      g_error ("%s: %s.", (str), error_msg); \
      g_free (error_msg); \
    } \
} G_STMT_END

static void
nc_data_planck_lkl_init (NcDataPlanckLKL *plik)
{
  plik->filename    = NULL;
  plik->obj         = NULL;
  plik->is_lensing  = FALSE;
  plik->nparams     = 0;
  plik->pnames      = NULL;
  plik->chksum      = NULL;
  plik->cmb_data    = 0;
  plik->data_params = NULL;
  plik->data_TT     = NULL;
  plik->data_EE     = NULL;
  plik->data_BB     = NULL;
  plik->data_TE     = NULL;
  plik->data_TB     = NULL;
  plik->data_EB     = NULL;
  plik->params      = NULL;
}

static void
nc_data_planck_lkl_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataPlanckLKL *plik = NC_DATA_PLANCK_LKL (object);
  g_return_if_fail (NC_IS_DATA_PLANCK_LKL (object));

  switch (prop_id)
  {
    case PROP_DATA_FILE:
      _nc_data_planck_set_filename (plik, g_value_get_string (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_planck_lkl_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataPlanckLKL *plik = NC_DATA_PLANCK_LKL (object);
  g_return_if_fail (NC_IS_DATA_PLANCK_LKL (object));

  switch (prop_id)
  {
    case PROP_DATA_FILE:
      g_value_set_string (value, plik->filename);
      break;
    case PROP_IS_LENSING:
      g_value_set_boolean (value, plik->is_lensing);
      break;
    case PROP_NPARAMS:
      g_value_set_uint (value, plik->nparams);
      break;
    case PROP_CHKSUM:
      g_value_set_string (value, plik->chksum);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_planck_lkl_dispose (GObject *object)
{
  NcDataPlanckLKL *plik = NC_DATA_PLANCK_LKL (object);

  ncm_vector_clear (&plik->data_params);
  ncm_vector_clear (&plik->data_TT);
  ncm_vector_clear (&plik->data_EE);
  ncm_vector_clear (&plik->data_BB);
  ncm_vector_clear (&plik->data_TE);
  ncm_vector_clear (&plik->data_TB);
  ncm_vector_clear (&plik->data_EB);
  ncm_vector_clear (&plik->params);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_planck_lkl_parent_class)->dispose (object);
}

static void
nc_data_planck_lkl_finalize (GObject *object)
{
  NcDataPlanckLKL *plik = NC_DATA_PLANCK_LKL (object);

  g_clear_pointer (&plik->filename, g_free);
  g_clear_pointer (&plik->chksum, g_free);

  if (plik->obj != NULL)
  {
    if (plik->is_lensing)
    {
      clik_lensing_object *obj = CLIK_LENS_OBJ (plik->obj);
      clik_lensing_cleanup (&obj);
    }
    else
    {
      clik_object *obj = CLIK_OBJ (plik->obj);
      clik_cleanup (&obj);
    }
    plik->obj = NULL;
  }

  g_strfreev (plik->pnames);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_planck_lkl_parent_class)->finalize (object);
}

static void
nc_data_planck_lkl_class_init (NcDataPlanckLKLClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = nc_data_planck_lkl_set_property;
  object_class->get_property = nc_data_planck_lkl_get_property;
  object_class->dispose      = nc_data_planck_lkl_dispose;
  object_class->finalize     = nc_data_planck_lkl_finalize;

  g_object_class_install_property (object_class,
                                   PROP_DATA_FILE,
                                   g_param_spec_string ("data-file",
                                                        NULL,
                                                        "Data file",
                                                        "no-file",
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_IS_LENSING,
                                   g_param_spec_boolean ("is-lensing",
                                                        NULL,
                                                        "Whether the likelihood has lensing",
                                                        FALSE,
                                                        G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_NPARAMS,
                                   g_param_spec_uint ("nparams",
                                                      NULL,
                                                      "Number of expected params",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_CHKSUM,
                                   g_param_spec_string ("checksum",
                                                        NULL,
                                                        "Params names checksum",
                                                        NULL,
                                                        G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}


void
_nc_data_planck_set_filename (NcDataPlanckLKL *plik, const gchar *filename)
{
  g_assert (plik->filename == NULL);
  g_assert (plik->obj == NULL);
  g_assert (filename != NULL);

  plik->filename = g_strdup (filename);

  g_assert (g_file_test (filename, G_FILE_TEST_EXISTS));
  {
    gchar *names_array[256];
    parname *names = NULL;
    error *err = initError ();
    gint data_params_len = 0;
    guint i;

    plik->is_lensing = clik_try_lensing (plik->filename, &err) == 1;
    CLIK_CHECK_ERROR ("_nc_data_planck_set_filename[clik_try_lensing]", err);

#define N_LENS_CMP 7
#define N_CMP 6
    if (plik->is_lensing)
    {
      gint lmax[N_LENS_CMP];
      guint vec_pos = 0;
      NcmVector **data_vec[N_LENS_CMP] =
      {
        &plik->data_PHIPHI, &plik->data_TT, &plik->data_EE,
        &plik->data_BB, &plik->data_TE, &plik->data_TB,
        &plik->data_EB
      };
      NcDataCMBDataType data_type_vec[N_LENS_CMP] =
      {
        NC_DATA_CMB_TYPE_PHIPHI,
        NC_DATA_CMB_TYPE_TT,
        NC_DATA_CMB_TYPE_EE,
        NC_DATA_CMB_TYPE_BB,
        NC_DATA_CMB_TYPE_TE,
        NC_DATA_CMB_TYPE_TB,
        NC_DATA_CMB_TYPE_EB
      };

      clik_lensing_object *obj = clik_lensing_init (plik->filename, &err);
      CLIK_CHECK_ERROR ("_nc_data_planck_set_filename[clik_lensing_init]", err);
      plik->obj = obj;

      plik->nparams = clik_lensing_get_extra_parameter_names (obj, &names, &err);
      CLIK_CHECK_ERROR ("_nc_data_planck_set_filename[clik_lensing_get_extra_parameter_names]", err);

      data_params_len += plik->nparams;

      clik_lensing_get_lmaxs (obj, lmax, &err);
      CLIK_CHECK_ERROR ("_nc_data_planck_set_filename[clik_lensing_get_lmaxs]", err);

      for (i = 0; i < N_LENS_CMP; i++)
        data_params_len += lmax[i] >= 0 ? (lmax[i] + 1) : 0;

      plik->data_params = ncm_vector_new (data_params_len);

      for (i = 0; i < N_LENS_CMP; i++)
      {
        if (lmax[i] >= 0)
        {
          plik->cmb_data |= data_type_vec[i];
          data_vec[i][0] = ncm_vector_get_subvector (plik->data_params, vec_pos, lmax[i] + 1);
          vec_pos += lmax[i] + 1;
        }
      }
      plik->params = ncm_vector_get_subvector (plik->data_params, vec_pos, plik->nparams);
    }
    else
    {
      gint lmax[N_CMP];
      guint vec_pos = 0;
      NcmVector **data_vec[N_CMP] =
      {
        &plik->data_TT, &plik->data_EE,
        &plik->data_BB, &plik->data_TE, &plik->data_TB,
        &plik->data_EB
      };
      NcDataCMBDataType data_type_vec[N_CMP] =
      {
        NC_DATA_CMB_TYPE_TT,
        NC_DATA_CMB_TYPE_EE,
        NC_DATA_CMB_TYPE_BB,
        NC_DATA_CMB_TYPE_TE,
        NC_DATA_CMB_TYPE_TB,
        NC_DATA_CMB_TYPE_EB
      };

      clik_object *obj = clik_init (plik->filename, &err);
      CLIK_CHECK_ERROR ("_nc_data_planck_set_filename[clik_init]", err);
      plik->obj = obj;

      plik->nparams = clik_get_extra_parameter_names (obj, &names, &err);
      CLIK_CHECK_ERROR ("_nc_data_planck_set_filename[clik_get_extra_parameter_names]", err);

      data_params_len += plik->nparams;

      clik_get_lmax (obj, lmax, &err);
      CLIK_CHECK_ERROR ("_nc_data_planck_set_filename[clik_get_lmaxs]", err);

      printf ("###################################\n");
      printf ("# data_params_len %u\n", data_params_len);
      for (i = 0; i < N_CMP; i++)
      {
        data_params_len += lmax[i] >= 0 ? (lmax[i] + 1) : 0;
        printf ("# data_params_len %d, %d\n", lmax[i], data_params_len);
      }

      plik->data_params = ncm_vector_new (data_params_len);

      for (i = 0; i < N_CMP; i++)
      {
        if (lmax[i] >= 0)
        {
          plik->cmb_data |= data_type_vec[i];
          data_vec[i][0] = ncm_vector_get_subvector (plik->data_params, vec_pos, lmax[i] + 1);
          vec_pos += lmax[i] + 1;
        }
      }
      plik->params = ncm_vector_get_subvector (plik->data_params, vec_pos, plik->nparams);
    }

    if (plik->nparams > 0)
    {
      GChecksum *names_checksum = g_checksum_new (G_CHECKSUM_MD5);
      for (i = 0; i < plik->nparams; i++)
      {
        names_array[i] = names[i];
        g_checksum_update (names_checksum, (guchar *)names[i], strlen (names[i]));
      }
      names_array[plik->nparams] = NULL;

      plik->chksum = g_strdup (g_checksum_get_string (names_checksum));

      g_checksum_free (names_checksum);
      plik->pnames = g_strdupv (names_array);
    }

    g_clear_pointer (&names, g_free);
    endError (&err);
  }
}

/**
 * nc_data_planck_lkl_new:
 * @filename: a Planck likelihood file
 *
 * FIXME
 *
 * Returns: a new #NcDataPlanckLKL
 */
NcDataPlanckLKL *
nc_data_planck_lkl_new (const gchar *filename)
{
  NcDataPlanckLKL *plik = g_object_new (NC_TYPE_DATA_PLANCK_LKL,
                                        "data-file", filename,
                                        NULL);
  return plik;
}

/**
 * nc_data_planck_lkl_get_param_name:
 * @plik: a #NcDataPlanckLKL
 * @i: param index
 *
 * FIXME
 *
 * Returns: (transfer none): a string constaining the param name
 */
const gchar *
nc_data_planck_lkl_get_param_name (NcDataPlanckLKL *plik, guint i)
{
  g_assert_cmpuint (plik->nparams, >, i);
  return plik->pnames[i];
}

/**
 * nc_data_planck_lkl_get_param_names:
 * @plik: a #NcDataPlanckLKL
 *
 * FIXME
 *
 * Returns: (array zero-terminated=1) (element-type utf8) (transfer full): an array of strings constaining the param names
 */
gchar **
nc_data_planck_lkl_get_param_names (NcDataPlanckLKL *plik)
{
  return g_strdupv (plik->pnames);
}
