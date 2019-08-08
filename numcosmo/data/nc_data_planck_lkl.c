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
 * @title: NcDataPlanckLKL
 * @short_description: Planck Likelihood interface.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_mset.h"
#include "math/ncm_model.h"
#include "nc_planck_fi.h"
#include "nc_planck_fi_cor_tt.h"
#include "data/nc_data_planck_lkl.h"

#ifndef NUMCOSMO_GIR_SCAN
#include "plc/clik.h"
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_DATA_FILE,
  PROP_PERT_BOLTZMANN,
  PROP_IS_LENSING,
  PROP_NPARAMS,
  PROP_CHKSUM,
  PROP_SIZE
};

G_DEFINE_TYPE (NcDataPlanckLKL, nc_data_planck_lkl, NCM_TYPE_DATA);

static void _nc_data_planck_lkl_set_filename (NcDataPlanckLKL *plik, const gchar *filename);

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
  plik->filename           = NULL;
  plik->pb                = NULL;
  plik->obj               = NULL;
  plik->is_lensing        = FALSE;
  plik->nparams           = 0;
  plik->ndata_entry       = 0;
  plik->pnames            = NULL;
  plik->chksum            = NULL;
	plik->check_m2lnL       = 0.0;
  plik->cmb_data          = 0;
  plik->data_params       = NULL;
  plik->check_data_params = NULL;
	plik->data_PHIPHI       = NULL;
  plik->data_TT           = NULL;
  plik->data_EE           = NULL;
  plik->data_BB           = NULL;
  plik->data_TE           = NULL;
  plik->data_TB           = NULL;
  plik->data_EB           = NULL;
  plik->params            = NULL;
  plik->pfi_ctrl           = ncm_model_ctrl_new (NULL);
  plik->cosmo_ctrl        = ncm_model_ctrl_new (NULL);
  plik->cm2lnL            = 0.0;
  plik->A_planck          = 0.0;
  plik->param_map         = g_array_new (TRUE, TRUE, sizeof (guint));
}

static void
nc_data_planck_lkl_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataPlanckLKL *plik = NC_DATA_PLANCK_LKL (object);
  g_return_if_fail (NC_IS_DATA_PLANCK_LKL (object));

  switch (prop_id)
  {
    case PROP_DATA_FILE:
      _nc_data_planck_lkl_set_filename (plik, g_value_get_string (value));
      break;
    case PROP_PERT_BOLTZMANN:
      nc_data_planck_lkl_set_hipert_boltzmann (plik, g_value_get_object (value));
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
    case PROP_PERT_BOLTZMANN:
      g_value_set_object (value, plik->pb);
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

  nc_hipert_boltzmann_clear (&plik->pb);

  ncm_vector_clear (&plik->data_params);
  ncm_vector_clear (&plik->check_data_params);
  ncm_vector_clear (&plik->data_PHIPHI);
  ncm_vector_clear (&plik->data_TT);
  ncm_vector_clear (&plik->data_EE);
  ncm_vector_clear (&plik->data_BB);
  ncm_vector_clear (&plik->data_TE);
  ncm_vector_clear (&plik->data_TB);
  ncm_vector_clear (&plik->data_EB);
  ncm_vector_clear (&plik->params);

  ncm_model_ctrl_clear (&plik->pfi_ctrl);
  ncm_model_ctrl_clear (&plik->cosmo_ctrl);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_planck_lkl_parent_class)->dispose (object);
}
  
static void
nc_data_planck_lkl_finalize (GObject *object)
{
  NcDataPlanckLKL *plik = NC_DATA_PLANCK_LKL (object);

  g_clear_pointer (&plik->filename, g_free);
  g_clear_pointer (&plik->chksum, g_free);
  g_clear_pointer (&plik->param_map, g_array_unref);

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

static guint _nc_data_planck_lkl_get_length (NcmData *data);
static void _nc_data_planck_lkl_begin (NcmData *data);
static void _nc_data_planck_lkl_prepare (NcmData *data, NcmMSet *mset);
/*static void _nc_data_planck_lkl_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng);*/
static void _nc_data_planck_lkl_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);

static void
nc_data_planck_lkl_class_init (NcDataPlanckLKLClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);

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
                                   PROP_PERT_BOLTZMANN,
                                   g_param_spec_object ("hipert-boltzmann",
                                                        NULL,
                                                        "NcHIPertBoltzmann object",
                                                        NC_TYPE_HIPERT_BOLTZMANN,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

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

  data_class->get_length = &_nc_data_planck_lkl_get_length;
  data_class->begin      = &_nc_data_planck_lkl_begin;
  data_class->prepare    = &_nc_data_planck_lkl_prepare;
  /*data_class->resample   = &_nc_data_planck_lkl_resample;*/
  data_class->m2lnL_val  = &_nc_data_planck_lkl_m2lnL_val;
}

static guint
_nc_data_planck_lkl_get_length (NcmData *data)
{
  return NC_DATA_PLANCK_LKL (data)->ndata_entry;
}

static void
_nc_data_planck_lkl_begin (NcmData *data)
{
  /* Nothing to do. */
  NCM_UNUSED (data);
}

static void
_nc_data_planck_lkl_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataPlanckLKL *clik = NC_DATA_PLANCK_LKL (data);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));

  if (clik->pb == NULL)
    g_error ("_nc_data_planck_lkl_prepare: cannot prepare without a #NcHIPertBoltzmann object. Use nc_data_planck_lkl_set_hipert_boltzmann to set the perturbations object.");
  if (cosmo == NULL)
    g_error ("_nc_data_planck_lkl_prepare: cannot prepare without a #NcHICosmo object. Add one to the #NcmMSet.");

  nc_hipert_boltzmann_prepare_if_needed (clik->pb, cosmo);
}

/*static void
_nc_data_planck_lkl_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng)
{
}
*/

/*
 * Some parameters have different upper/lower case combination in different
 * likelihood, so we make a case insensitive search, if it is found
 * update the name to mach the model's.
 */
static gboolean
_nc_data_planck_lkl_find_param (NcmModel *model, gchar *name, guint *pi)
{
  guint nparams = ncm_model_len (model);
  guint i;

  for (i = 0; i < nparams; i++)
  {
    const gchar *mname = ncm_model_param_name (model, i);
    const guint ns = strlen (name);
    if (strncasecmp (name, mname, ns) == 0)
    {
      memcpy (name, mname, ns);
      *pi = i;
      return TRUE;
    }
  }
  return FALSE;
}

static void
_nc_data_planck_lkl_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcDataPlanckLKL *clik = NC_DATA_PLANCK_LKL (data);
  gdouble *cl_and_pars  = ncm_vector_ptr (clik->data_params, 0);
  error *err            = initError ();
  NcPlanckFI *pfi  = NC_PLANCK_FI (ncm_mset_peek (mset, nc_planck_fi_id ()));
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  gboolean pfi_model_up = FALSE;
  gboolean pfi_up       = FALSE;
  gboolean cosmo_up     = FALSE;
  guint i;

  if (pfi == NULL)
    g_error ("_nc_data_planck_lkl_m2lnL_val: a NcPlanckFI* model is required in NcmMSet.");

  if (clik->nparams > 0)
  {
    if ((pfi_model_up = ncm_model_ctrl_model_update (clik->pfi_ctrl, NCM_MODEL (pfi))))
    {
      g_array_set_size (clik->param_map, clik->nparams);
      
      for (i = 0; i < clik->nparams; i++)
      {
        guint pi = 0;
        gboolean pfound = ncm_model_param_index_from_name (NCM_MODEL (pfi), clik->pnames[i], &pi);
        if (!pfound)
        {
          gboolean pfound2 = _nc_data_planck_lkl_find_param (NCM_MODEL (pfi), clik->pnames[i], &pi);
          if (!pfound2)
          {
            g_error ("_nc_data_planck_lkl_m2lnL_val: cannot find parameter `%s' in models `%s'.",
                     clik->pnames[i], ncm_model_name (NCM_MODEL (pfi)));
          }
        }
        g_array_index (clik->param_map, guint, i) = pi;
      }
    }

    for (i = 0; i < clik->nparams; i++)
    {
      const gdouble p_i = ncm_model_param_get (NCM_MODEL (pfi), g_array_index (clik->param_map, guint, i));
      ncm_vector_set (clik->params, i, p_i);
    }
  }

  pfi_up    = ncm_model_ctrl_update (clik->pfi_ctrl, NCM_MODEL (pfi));
  cosmo_up = ncm_model_ctrl_update (clik->cosmo_ctrl, NCM_MODEL (cosmo));

  if (pfi_up && clik->nparams == 1 && g_array_index (clik->param_map, guint, 0) == NC_PLANCK_FI_COR_TT_A_planck)
  {
    /*printf ("It's the case!\n");*/
    if (clik->A_planck == ncm_vector_get (clik->params, 0))
    {
      /*printf ("Not changed % 20.16g == % 20.16g\n", clik->A_planck, ncm_vector_get (clik->params, 0));*/
      pfi_up = FALSE;
    }
    else
    {
      /*printf ("Updating    % 20.16g to % 20.16g\n", clik->A_planck, ncm_vector_get (clik->params, 0));*/
      clik->A_planck = ncm_vector_get (clik->params, 0);
    }
  }
  
  if (cosmo_up)
  {
    if (clik->cmb_data & NC_DATA_CMB_TYPE_PHIPHI)
      nc_hipert_boltzmann_get_PHIPHI_Cls (clik->pb, clik->data_PHIPHI);

		if (clik->cmb_data & NC_DATA_CMB_TYPE_TT)
      nc_hipert_boltzmann_get_TT_Cls (clik->pb, clik->data_TT);

    if (clik->cmb_data & NC_DATA_CMB_TYPE_EE)
      nc_hipert_boltzmann_get_EE_Cls (clik->pb, clik->data_EE);

    if (clik->cmb_data & NC_DATA_CMB_TYPE_BB)
      nc_hipert_boltzmann_get_BB_Cls (clik->pb, clik->data_BB);

    if (clik->cmb_data & NC_DATA_CMB_TYPE_TE)
      nc_hipert_boltzmann_get_TE_Cls (clik->pb, clik->data_TE);

    if (clik->cmb_data & NC_DATA_CMB_TYPE_TB)
      nc_hipert_boltzmann_get_TB_Cls (clik->pb, clik->data_TB);

    if (clik->cmb_data & NC_DATA_CMB_TYPE_EB)
      nc_hipert_boltzmann_get_EB_Cls (clik->pb, clik->data_EB);
  }
/*
  ncm_mset_pretty_log (mset);
  printf ("# data[%p] cosmo_up:%d nparams:%u pfi_up:%d pfi_model_up:%d is_lensing:%d \n", data, cosmo_up, clik->nparams, pfi_up, pfi_model_up, clik->is_lensing);
*/  
  if (cosmo_up || (clik->nparams > 0 && (pfi_up || pfi_model_up)))
  {
    /*ncm_vector_log_vals (clik->data_params, "cl and vals: ", "% 22.15g");*/
    
    if (clik->is_lensing)
    {
      *m2lnL = -2.0 * clik_lensing_compute (clik->obj, cl_and_pars, &err);
      CLIK_CHECK_ERROR ("_nc_data_planck_lkl_m2lnL_val[clik_lensing_compute]", err);
    }
    else
    {
      *m2lnL = -2.0 * clik_compute (clik->obj, cl_and_pars, &err);
      CLIK_CHECK_ERROR ("_nc_data_planck_lkl_m2lnL_val[clik_compute]", err);

    }
    /*printf ("# m2lnL % 22.15g\n", *m2lnL);*/
    clik->cm2lnL = *m2lnL;
  }
  else
  {
    *m2lnL = clik->cm2lnL;
  }

  endError (&err);
  return;
}

static void
_nc_data_planck_lkl_set_filename (NcDataPlanckLKL *plik, const gchar *filename)
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

    plik->is_lensing = (clik_try_lensing (plik->filename, &err) == 1);
    CLIK_CHECK_ERROR ("_nc_data_planck_lkl_set_filename[clik_try_lensing]", err);

    {
      const gchar *type_name = g_type_name (NC_TYPE_DATA_PLANCK_LKL);
      gchar *bfile     = g_path_get_basename (filename);
      gchar *data_desc = g_strdup_printf ("%s[%s]", type_name, bfile);
      ncm_data_set_desc (NCM_DATA (plik), data_desc);
      g_free (bfile);
      g_free (data_desc);
    }

    plik->ndata_entry = 0;

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
      CLIK_CHECK_ERROR ("_nc_data_planck_lkl_set_filename[clik_lensing_init]", err);
      plik->obj = obj;

      plik->nparams = clik_lensing_get_extra_parameter_names (obj, &names, &err);
      CLIK_CHECK_ERROR ("_nc_data_planck_lkl_set_filename[clik_lensing_get_extra_parameter_names]", err);

      clik_lensing_get_lmaxs (obj, lmax, &err);
      CLIK_CHECK_ERROR ("_nc_data_planck_lkl_set_filename[clik_lensing_get_lmaxs]", err);

      for (i = 0; i < N_LENS_CMP; i++)
        data_params_len += lmax[i] >= 0 ? (lmax[i] + 1) : 0;

      plik->ndata_entry = data_params_len;

      data_params_len += plik->nparams;
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
      CLIK_CHECK_ERROR ("_nc_data_planck_lkl_set_filename[clik_init]", err);
      plik->obj = obj;

      plik->nparams = clik_get_extra_parameter_names (obj, &names, &err);
      CLIK_CHECK_ERROR ("_nc_data_planck_lkl_set_filename[clik_get_extra_parameter_names]", err);

      clik_get_lmax (obj, lmax, &err);
      CLIK_CHECK_ERROR ("_nc_data_planck_lkl_set_filename[clik_get_lmaxs]", err);

      /*printf ("###################################\n");*/
      /*printf ("# data_params_len %u\n", data_params_len);*/
      for (i = 0; i < N_CMP; i++)
      {
        data_params_len += lmax[i] >= 0 ? (lmax[i] + 1) : 0;
        /*printf ("# data_params_len %d, %d\n", lmax[i], data_params_len);*/
      }

      plik->ndata_entry = data_params_len;

      data_params_len += plik->nparams;
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

    {
      gint npar_out       = 0;
      gdouble check_value = 0.0;
      gdouble *chkp       = NULL;

			if (plik->is_lensing)
				clik_lensing_get_check_param (plik->obj, plik->filename, &chkp, &check_value, &npar_out, &err);
			else
				clik_get_check_param (plik->obj, plik->filename, &chkp, &check_value, &npar_out, &err);
			
      ncm_vector_set_data (plik->data_params, chkp, npar_out);
      ncm_vector_clear (&plik->check_data_params);

      plik->check_data_params = ncm_vector_dup (plik->data_params);
			plik->check_m2lnL       = -2.0 * check_value;

			{
				gdouble check_m2lnL = 0.0;
				if (plik->is_lensing)
				{
					check_m2lnL = -2.0 * clik_lensing_compute (plik->obj, chkp, &err);
					CLIK_CHECK_ERROR ("_nc_data_planck_lkl_m2lnL_val[clik_lensing_compute]", err);
				}
				else
				{
					check_m2lnL = -2.0 * clik_compute (plik->obj, chkp, &err);
					CLIK_CHECK_ERROR ("_nc_data_planck_lkl_m2lnL_val[clik_compute]", err);
				}
				ncm_assert_cmpdouble_e (check_m2lnL, ==, plik->check_m2lnL, 1.0e-4, 0.0);
			}
			
      g_free (chkp);
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

  ncm_data_set_init (NCM_DATA (plik), TRUE);
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
 * nc_data_planck_lkl_full_new:
 * @filename: a Planck likelihood file
 * @pb: a #NcHIPertBoltzmann
 *
 * FIXME
 *
 * Returns: a new #NcDataPlanckLKL
 */
NcDataPlanckLKL *
nc_data_planck_lkl_full_new (const gchar *filename, NcHIPertBoltzmann *pb)
{
  NcDataPlanckLKL *plik = g_object_new (NC_TYPE_DATA_PLANCK_LKL,
                                        "data-file", filename,
                                        "hipert-boltzmann", pb,
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

/**
 * nc_data_planck_lkl_set_hipert_boltzmann:
 * @plik: a #NcDataPlanckLKL
 * @pb: a #NcHIPertBoltzmann
 *
 * Sets the #NcHIPertBoltzmann to be used in the likelihood calculation.
 *
 */
void
nc_data_planck_lkl_set_hipert_boltzmann (NcDataPlanckLKL *plik, NcHIPertBoltzmann *pb)
{
  nc_hipert_boltzmann_clear (&plik->pb);
  plik->pb = nc_hipert_boltzmann_ref (pb);

  nc_hipert_boltzmann_append_target_Cls (plik->pb, plik->cmb_data);

  if (plik->data_PHIPHI != NULL)
  {
    guint PHIPHI_lmax = ncm_vector_len (plik->data_PHIPHI) - 1;
    if (PHIPHI_lmax > nc_hipert_boltzmann_get_PHIPHI_lmax (plik->pb))
      nc_hipert_boltzmann_set_PHIPHI_lmax (plik->pb, PHIPHI_lmax);
  }
  if (plik->data_TT != NULL)
  {
    guint TT_lmax = ncm_vector_len (plik->data_TT) - 1;
    if (TT_lmax > nc_hipert_boltzmann_get_TT_lmax (plik->pb))
      nc_hipert_boltzmann_set_TT_lmax (plik->pb, TT_lmax);
  }
  if (plik->data_EE != NULL)
  {
    guint EE_lmax = ncm_vector_len (plik->data_EE) - 1;
    if (EE_lmax > nc_hipert_boltzmann_get_EE_lmax (plik->pb))
      nc_hipert_boltzmann_set_EE_lmax (plik->pb, EE_lmax);
  }
  if (plik->data_BB != NULL)
  {
    guint BB_lmax = ncm_vector_len (plik->data_BB) - 1;
    if (BB_lmax > nc_hipert_boltzmann_get_BB_lmax (plik->pb))
      nc_hipert_boltzmann_set_BB_lmax (plik->pb, BB_lmax);
  }
  if (plik->data_TE != NULL)
  {
    guint TE_lmax = ncm_vector_len (plik->data_TE) - 1;
    if (TE_lmax > nc_hipert_boltzmann_get_TE_lmax (plik->pb))
      nc_hipert_boltzmann_set_TE_lmax (plik->pb, TE_lmax);
  }
  if (plik->data_TB != NULL)
  {
    guint TB_lmax = ncm_vector_len (plik->data_TB) - 1;
    if (TB_lmax > nc_hipert_boltzmann_get_TB_lmax (plik->pb))
      nc_hipert_boltzmann_set_TB_lmax (plik->pb, TB_lmax);
  }
  if (plik->data_EB != NULL)
  {
    guint EB_lmax = ncm_vector_len (plik->data_EB) - 1;
    if (EB_lmax > nc_hipert_boltzmann_get_EB_lmax (plik->pb))
      nc_hipert_boltzmann_set_EB_lmax (plik->pb, EB_lmax);
  }
}
