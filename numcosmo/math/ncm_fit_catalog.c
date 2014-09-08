/***************************************************************************
 *            ncm_fit_catalog.c
 *
 *  Tue February 18 10:49:26 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_fit_catalog.c
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_fit_catalog
 * @title: Fit Catalog
 * @short_description: Ordered catalog of different fits or likelihood values.
 *
 * This class defines a catalog type object. This object can automatically syncronise
 * with a fits file (thought cfitsio). 
 * 
 * For Motecarlo studies like resampling from a fiducial model or bootstrap it is used
 * to save the values of the best-fit values for each realization. Since the order of
 * resampling is important due to the fact that we use the same pseudo-random number 
 * generator for all resamplings, this object also guarantees the order of the samples
 * added.
 * 
 * For Markov Chain Montecarlo this object saves the value of the same likelihood in
 * different points of the parameter space.
 * 
 * For both applications this object keeps an interactive mean and variance for the
 * parameters added, this allows a sample by sample analyses of the convergence.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fit_catalog.h"
#include "math/ncm_cfg.h"
#include "math/ncm_func_eval.h"
#include "ncm_enum_types.h"

G_DEFINE_TYPE (NcmFitCatalog, ncm_fit_catalog, G_TYPE_OBJECT);

enum
{
  PROP_0,
  PROP_FIT,
  PROP_FILE,
  PROP_RUN_TYPE_STR,
  PROP_FLUSH_MODE,
  PROP_FLUSH_INTERVAL,
};

static void
ncm_fit_catalog_init (NcmFitCatalog *fcat)
{
  fcat->fit           = NULL;
  fcat->pstats        = NULL;
  fcat->fmode         = NCM_FIT_CATALOG_FLUSH_LEN;
  fcat->flush_timer   = g_timer_new ();
  fcat->cur_id        = -1; /* Represents that no elements in the catalog. */
  fcat->first_id      = 0;  /* The element to be in the catalog will be the one with index == 0 */
  fcat->file_cur_id   = -1; /* Represents that no elements in the catalog file. */
  fcat->file_first_id = 0;  /* The element to be in the catalog file will be the one with index == 0 */
  fcat->file          = NULL;
  fcat->rtype_str     = NULL;
  fcat->porder        = g_array_new (FALSE, FALSE, sizeof (gint));
#ifdef NUMCOSMO_HAVE_CFITSIO
  fcat->fptr          = NULL;
#endif /* NUMCOSMO_HAVE_CFITSIO */
  fcat->pdf_i         = -1;
  fcat->h             = NULL;
  fcat->h_pdf         = NULL;
  fcat->params_max    = NULL;
  fcat->params_min    = NULL;
}

static void _ncm_fit_catalog_set_fit_obj (NcmFitCatalog *fcat, NcmFit *fit);

static void
_ncm_fit_catalog_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFitCatalog *fcat = NCM_FIT_CATALOG (object);
  g_return_if_fail (NCM_IS_FIT_CATALOG (object));

  switch (prop_id)
  {
    case PROP_FIT:
      _ncm_fit_catalog_set_fit_obj (fcat, g_value_get_object (value));
      break;
    case PROP_FILE:
      ncm_fit_catalog_set_file (fcat, g_value_get_string (value));
      break;    
    case PROP_RUN_TYPE_STR:
      ncm_fit_catalog_set_run_type (fcat, g_value_get_string (value));
      break;    
    case PROP_FLUSH_MODE:
      ncm_fit_catalog_set_flush_mode (fcat, g_value_get_enum (value));
      break;
    case PROP_FLUSH_INTERVAL:
      ncm_fit_catalog_set_flush_interval (fcat, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_fit_catalog_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFitCatalog *fcat = NCM_FIT_CATALOG (object);
  g_return_if_fail (NCM_IS_FIT_CATALOG (object));

  switch (prop_id)
  {
    case PROP_FIT:
      g_value_set_object (value, fcat->fit);
      break;
    case PROP_FILE:
      g_value_set_string (value, fcat->file);
      break;
    case PROP_RUN_TYPE_STR:
      g_value_set_string (value, fcat->rtype_str);
      break;
    case PROP_FLUSH_MODE:
      g_value_set_enum (value, fcat->fmode);
      break;
    case PROP_FLUSH_INTERVAL:
      g_value_set_double (value, fcat->flush_interval);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}
static void
_ncm_fit_catalog_dispose (GObject *object)
{
  NcmFitCatalog *fcat = NCM_FIT_CATALOG (object);

  ncm_fit_clear (&fcat->fit);
  ncm_stats_vec_clear (&fcat->pstats);
  ncm_vector_free (fcat->params_max);
  ncm_vector_free (fcat->params_min);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_catalog_parent_class)->dispose (object);
}

static void _ncm_fit_catalog_close_file (NcmFitCatalog *fcat);

static void
_ncm_fit_catalog_finalize (GObject *object)
{
  NcmFitCatalog *fcat = NCM_FIT_CATALOG (object);

  if (fcat->h != NULL)
    gsl_histogram_free (fcat->h);
  if (fcat->h_pdf != NULL)
    gsl_histogram_pdf_free (fcat->h_pdf);

  _ncm_fit_catalog_close_file (fcat);

  g_clear_pointer (&fcat->rtype_str, g_free);
  
  g_array_unref (fcat->porder);
  g_timer_destroy (fcat->flush_timer);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_catalog_parent_class)->finalize (object);
}

static void
ncm_fit_catalog_class_init (NcmFitCatalogClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_fit_catalog_set_property;
  object_class->get_property = &_ncm_fit_catalog_get_property;
  object_class->dispose      = &_ncm_fit_catalog_dispose;
  object_class->finalize     = &_ncm_fit_catalog_finalize;

  g_object_class_install_property (object_class,
                                   PROP_FIT,
                                   g_param_spec_object ("fit",
                                                        NULL,
                                                        "Fit object",
                                                        NCM_TYPE_FIT,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_FLUSH_MODE,
                                   g_param_spec_enum ("fmode",
                                                      NULL,
                                                      "Catalog flush mode",
                                                      NCM_TYPE_FIT_CATALOG_FLUSH, NCM_FIT_CATALOG_FLUSH_AUTO,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_FILE,
                                   g_param_spec_string ("filename",
                                                        NULL,
                                                        "Catalog filename",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RUN_TYPE_STR,
                                   g_param_spec_string ("run-type-string",
                                                        NULL,
                                                        "Run type string",
                                                        NCM_FIT_CATALOG_RTYPE_UNDEFINED,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_FLUSH_INTERVAL,
                                   g_param_spec_double ("flush-interval",
                                                        NULL,
                                                        "Data flush interval",
                                                        0.0, 1.0e3, 10.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

}

static void 
_ncm_fit_catalog_set_fit_obj (NcmFitCatalog *fcat, NcmFit *fit)
{
  guint free_params_len = ncm_mset_fparams_len (fit->mset);
  
  fcat->fit        = ncm_fit_ref (fit);
  fcat->pstats     = ncm_stats_vec_new (free_params_len + 1, NCM_STATS_VEC_COV, TRUE);
  fcat->params_max = ncm_vector_new (free_params_len + 1);
  fcat->params_min = ncm_vector_new (free_params_len + 1);

  ncm_vector_set_all (fcat->params_max, GSL_NEGINF);
  ncm_vector_set_all (fcat->params_min, GSL_POSINF);

  g_array_set_size (fcat->porder, free_params_len + 1);
}

/**
 * ncm_fit_catalog_new:
 * @fit: a #NcmFit
 *
 * Creates a new #NcmFitCatalog based on the #NcmFit object @fit. The catalog assumes that
 * the @fit object will remain with the same set of free parameters during his whole lifetime.
 *
 * Returns: (transfer full): a new #NcmFitCatalog
 */
NcmFitCatalog *
ncm_fit_catalog_new (NcmFit *fit)
{
  return g_object_new (NCM_TYPE_FIT_CATALOG, 
                       "fit", fit,
                       NULL);
}

/**
 * ncm_fit_catalog_free:
 * @fcat: a #NcmFitCatalog
 *
 * Decrese the reference count of @fcat atomically.
 *
 */
void 
ncm_fit_catalog_free (NcmFitCatalog *fcat)
{
  g_object_unref (fcat);
}

/**
 * ncm_fit_catalog_clear:
 * @fcat: a #NcmFitCatalog
 *
 * Decrese the reference count of *@fcat atomically and sets the pointer *@fcat to null.
 *
 */
void 
ncm_fit_catalog_clear (NcmFitCatalog **fcat)
{
  ncm_fit_catalog_sync (*fcat, TRUE);  
  g_clear_object (fcat);
}

static void
_ncm_fit_catalog_open_create_file (NcmFitCatalog *fcat)
{
  gint status = 0;
  guint fparam_len = ncm_mset_fparam_len (fcat->fit->mset);
  guint i;
  
  g_assert (fcat->file != NULL);
  g_assert (fcat->fptr == NULL);

  if (g_file_test (fcat->file, G_FILE_TEST_EXISTS))
  {
    glong nrows;
    gint m2lnL_index;
    gchar key_text[FLEN_VALUE];

    fits_open_file (&fcat->fptr, fcat->file, READWRITE, &status);
    NCM_FITS_ERROR (status);

    fits_movnam_hdu (fcat->fptr, BINARY_TBL, NCM_FIT_CATALOG_EXTNAME, 0, &status);
    NCM_FITS_ERROR (status);

    fits_read_key (fcat->fptr, TINT, NCM_FIT_CATALOG_FIRST_ID_LABEL,
                   &fcat->file_first_id, NULL, &status);
    NCM_FITS_ERROR (status);

    fits_read_key (fcat->fptr, TSTRING, NCM_FIT_CATALOG_RTYPE_LABEL,
                   key_text, NULL, &status);
    NCM_FITS_ERROR (status);

    if (strcmp (fcat->rtype_str, key_text) != 0)
      g_error ("_ncm_fit_catalog_open_create_file: incompatible run type strings from catalog and file, catalog: `%s' file: `%s'.",
               fcat->rtype_str, key_text);

    fits_get_num_rows (fcat->fptr, &nrows, &status);
    NCM_FITS_ERROR (status);

    if (fcat->file_first_id != fcat->first_id)
    {
      if (nrows == 0)
      {
        if (fcat->file_first_id != 0)
          g_warning ("_ncm_fit_catalog_open_create_file: Empty data file with "NCM_FIT_CATALOG_FIRST_ID_LABEL" different from first_id: %d != %d. Setting to first_id.\n", 
                     fcat->file_first_id, fcat->first_id);
        fcat->file_first_id = fcat->first_id;
      } 
      else if (fcat->cur_id < fcat->first_id)
      {
        if (fcat->first_id != 0)
          g_warning ("_ncm_fit_catalog_open_create_file: Empty memory catalog with first_id different from "NCM_FIT_CATALOG_FIRST_ID_LABEL": %d != %d. Setting to "NCM_FIT_CATALOG_FIRST_ID_LABEL".\n", 
                     fcat->first_id, fcat->file_first_id);
        fcat->first_id = fcat->file_first_id;
        fcat->cur_id = fcat->file_first_id - 1;
      }
    }
    fcat->file_cur_id = fcat->file_first_id + nrows - 1; 

    if (fits_get_colnum (fcat->fptr, CASESEN, NCM_FIT_CATALOG_M2LNL_COLNAME, &m2lnL_index, &status))
      g_error ("_ncm_fit_catalog_open_create_file: Column "NCM_FIT_CATALOG_M2LNL_COLNAME" not found, invalid fits file.");
    if (m2lnL_index != 1)
      g_error ("_ncm_fit_catalog_open_create_file: Column "NCM_FIT_CATALOG_M2LNL_COLNAME" is not the first column [%d], invalid fits file.", 
               m2lnL_index);
    g_array_index (fcat->porder, gint, 0) = m2lnL_index;

    for (i = 0; i < fparam_len; i++)
    {
      const gchar *fparam_fullname = ncm_mset_fparam_full_name (fcat->fit->mset, i);
      if (fits_get_colnum (fcat->fptr, CASESEN, (gchar *)fparam_fullname, &g_array_index (fcat->porder, gint, i + 1), &status))
      {                /* I don't like this too ^^^^^^^^^  */
        g_error ("_ncm_fit_catalog_open_create_file: Column %s not found, invalid fits file.", fparam_fullname);
      }
    }
  }
  else
  {
    GPtrArray *ttype_array = g_ptr_array_sized_new (10);
    GPtrArray *tform_array = g_ptr_array_sized_new (10);
    fits_create_file (&fcat->fptr, fcat->file, &status);
    NCM_FITS_ERROR (status);

    g_ptr_array_add (ttype_array, NCM_FIT_CATALOG_M2LNL_COLNAME);
    g_ptr_array_add (tform_array, "1D");
    g_array_index (fcat->porder, gint, 0) = tform_array->len;

    for (i = 0; i < fparam_len; i++)
    {
      const gchar *fparam_fullname = ncm_mset_fparam_full_name (fcat->fit->mset, i);
      g_ptr_array_add (ttype_array, (gchar *)fparam_fullname);
      /* I don't like this too ^^^^^^^^^  */
      g_ptr_array_add (tform_array, "1D");
      g_array_index (fcat->porder, gint, i + 1) = tform_array->len;
    }

    /* append a new empty binary table onto the FITS file */
    fits_create_tbl (fcat->fptr, BINARY_TBL, 0, fparam_len + 1, (gchar **)ttype_array->pdata, (gchar **)tform_array->pdata,
                     NULL, NCM_FIT_CATALOG_EXTNAME, &status);
    NCM_FITS_ERROR (status);

    fits_update_key (fcat->fptr, TSTRING, NCM_FIT_CATALOG_RTYPE_LABEL, fcat->rtype_str, "Run type string.", &status);
    NCM_FITS_ERROR (status);  

    for (i = 0; i < fparam_len; i++)
    {
      gchar *fsymbi = g_strdup_printf ("%s%d", NCM_FIT_CATALOG_FSYMB_LABEL, i + 1);
      gchar *fsymb_desc = g_strdup_printf ("Symbol for parameter %s[%d]", 
                                           ncm_mset_fparam_name (fcat->fit->mset, i), 
                                           i + 1);
      const gchar *fsymb  = ncm_mset_fparam_symbol (fcat->fit->mset, i);
      
      fits_update_key (fcat->fptr, TSTRING, fsymbi, fsymb, fsymb_desc, &status);
      NCM_FITS_ERROR (status);  
      g_free (fsymbi);
      g_free (fsymb_desc);
    }
    
    fcat->file_first_id = 0;
    fcat->file_cur_id   = -1;
  }

  fits_update_key (fcat->fptr, TINT, NCM_FIT_CATALOG_FIRST_ID_LABEL, &fcat->file_first_id, "Id of the first element.", &status);
  NCM_FITS_ERROR (status);

  fits_flush_file (fcat->fptr, &status);
  NCM_FITS_ERROR (status);
}

static void _ncm_fit_catalog_flush_file (NcmFitCatalog *fcat);

/**
 * ncm_fit_catalog_set_file:
 * @fcat: a #NcmFitCatalog
 * @filename: a filename
 *
 * Sets the data filename to be used to sync/save data.
 *
 */
void 
ncm_fit_catalog_set_file (NcmFitCatalog *fcat, const gchar *filename)
{
  if (fcat->file != NULL && strcmp (fcat->file, filename) == 0)
    return;
  
  _ncm_fit_catalog_close_file (fcat);

  g_clear_pointer (&fcat->file, g_free);
  fcat->file = g_strdup (filename);  

  _ncm_fit_catalog_open_create_file (fcat);
  ncm_fit_catalog_sync (fcat, TRUE);

  _ncm_fit_catalog_flush_file (fcat);
}

/**
 * ncm_fit_catalog_set_flush_mode:
 * @fcat: a #NcmFitCatalog
 * @fmode: flush mode
 *
 * Sets the flush mode to @fmode.
 *
 */
void 
ncm_fit_catalog_set_flush_mode (NcmFitCatalog *fcat, NcmFitCatalogFlush fmode)
{
  fcat->fmode = fmode;
}

/**
 * ncm_fit_catalog_set_flush_interval:
 * @fcat: a #NcmFitCatalog
 * @interval: Minimum time interval between flushs
 *
 * Sets the minimum time interval between flushs.
 *
 */
void 
ncm_fit_catalog_set_flush_interval (NcmFitCatalog *fcat, gdouble interval)
{
  fcat->flush_interval = interval;
}

/**
 * ncm_fit_catalog_set_first_id:
 * @fcat: a #NcmFitCatalog
 * @first_id: the id of the first item in the sample.
 * 
 * Sets the first id of the catalog, mainly used to inform in which realization the catalog starts. 
 * 
 */
void 
ncm_fit_catalog_set_first_id (NcmFitCatalog *fcat, gint first_id)
{
  if (first_id == fcat->first_id)
    return;

  g_assert_cmpint (fcat->file_first_id, ==, fcat->first_id);
  g_assert_cmpint (fcat->file_cur_id, ==, fcat->cur_id);
  
  if (fcat->cur_id != fcat->first_id - 1)
    g_error ("ncm_fit_catalog_set_first_id: cannot modify first_id to %d in a non-empty catalog, catalog first id: %d, catalog current id: %d.", 
             first_id, fcat->first_id, fcat->cur_id);

  fcat->first_id = first_id;
  fcat->cur_id   = first_id - 1;

  fcat->file_first_id = first_id;
  fcat->file_cur_id   = first_id - 1;

  if (fcat->fptr != NULL)
  {
    gint status = 0;
    fits_update_key (fcat->fptr, TINT, NCM_FIT_CATALOG_FIRST_ID_LABEL, &fcat->file_first_id, "Id of the first element.", &status);
    NCM_FITS_ERROR (status);
    
    ncm_fit_catalog_sync (fcat, TRUE);
  }
}

/**
 * ncm_fit_catalog_set_run_type:
 * @fcat: a #NcmFitCatalog
 * @rtype_str: the run type string.
 * 
 * Sets the run type string. 
 * 
 */
void 
ncm_fit_catalog_set_run_type (NcmFitCatalog *fcat, const gchar *rtype_str)
{
  gint status = 0;
  g_assert (rtype_str != NULL);
  
  if (fcat->rtype_str != NULL)
  {
    if (strcmp (fcat->rtype_str, rtype_str) != 0)
    {
      if (fcat->cur_id + 1 != fcat->first_id)
        g_error ("ncm_fit_catalog_set_run_type: cannot change run type string in a non-empty catalog, actual: `%s' new: `%s'.",
                 fcat->rtype_str, rtype_str);
      else
      {
        g_clear_pointer (&fcat->rtype_str, g_free);
      }
    }
    else
      return;
  }

  fcat->rtype_str = g_strdup (rtype_str);
  if (fcat->fptr != NULL)
  {
    fits_update_key (fcat->fptr, TSTRING, NCM_FIT_CATALOG_RTYPE_LABEL, fcat->rtype_str, NULL, &status);
    NCM_FITS_ERROR (status);
  }
}

static void
_ncm_fit_catalog_flush_file (NcmFitCatalog *fcat)
{
  gint status = 0;
  gint64 nrows = fcat->file_cur_id - fcat->file_first_id + 1;

  fits_update_key (fcat->fptr, TLONGLONG, NCM_FIT_CATALOG_NROWS_LABEL, &nrows, NULL, &status);
  NCM_FITS_ERROR (status);

  fits_flush_buffer (fcat->fptr, 0, &status);
  NCM_FITS_ERROR (status);
}

static void
_ncm_fit_catalog_close_file (NcmFitCatalog *fcat)
{
  gint status = 0;
#ifdef NUMCOSMO_HAVE_CFITSIO
  if (fcat->fptr != NULL)
  {
    _ncm_fit_catalog_flush_file (fcat);
    
    fits_close_file (fcat->fptr, &status);
    NCM_FITS_ERROR (status);

    fcat->fptr = NULL;
  }
#endif /* NUMCOSMO_HAVE_CFITSIO */
}

static void 
_ncm_fit_catalog_write_row (NcmFitCatalog *fcat, NcmVector *row, guint row_index)
{
  guint i;
  gint status = 0;
  
  for (i = 0; i < ncm_vector_len (row); i++)
  {
    fits_write_col (fcat->fptr, TDOUBLE, g_array_index (fcat->porder, gint, i), row_index, 
                    1, 1, ncm_vector_ptr (row, i), &status);
    NCM_FITS_ERROR (status);
  }
}

static void 
_ncm_fit_catalog_read_row (NcmFitCatalog *fcat, NcmVector *row, guint row_index)
{
  guint i;
  gint status = 0;
    
  for (i = 0; i < ncm_vector_len (row); i++)
  {
    fits_read_col (fcat->fptr, TDOUBLE, g_array_index (fcat->porder, gint, i), row_index, 
                   1, 1, NULL, ncm_vector_ptr (row, i), NULL, &status);
    NCM_FITS_ERROR (status);
  }
}

/**
 * ncm_fit_catalog_sync:
 * @fcat: a #NcmFitCatalog
 * @check: whether to check consistence between file and memory data.
 * 
 * Syncronise memory and data file. If no file was defined simply returns.
 * 
 */
void
ncm_fit_catalog_sync (NcmFitCatalog *fcat, gboolean check)
{
  gint status = 0;
  guint i;

  if (fcat->file == NULL)
    return;

#ifdef NUMCOSMO_HAVE_CFITSIO
  g_assert (fcat->fptr != NULL);
  
  if (check)
  {
    gchar fptr_filename[FLEN_FILENAME];
  
    fits_file_name (fcat->fptr, fptr_filename, &status);
    NCM_FITS_ERROR (status);

    g_assert_cmpstr (fptr_filename, ==, fcat->file);

    if ((fcat->file_cur_id < fcat->first_id - 1) || (fcat->cur_id < fcat->file_first_id - 1))
      g_error ("ncm_fit_catalog_sync: file data & catalog mismatch, they do not intersect each other: file data [%d, %d] catalog [%d, %d]",
               fcat->file_first_id, fcat->file_cur_id,
               fcat->first_id, fcat->cur_id);
  }

  if (fcat->file_first_id != fcat->first_id)
  {
    if (fcat->file_first_id > fcat->first_id)
    {
      guint rows_to_add = fcat->file_first_id - fcat->first_id;
      fits_insert_rows (fcat->fptr, 0, rows_to_add, &status);
      NCM_FITS_ERROR (status);

      for (i = 0; i < rows_to_add; i++)
      {
        NcmVector *row = ncm_stats_vec_peek_row (fcat->pstats, i);
        _ncm_fit_catalog_write_row (fcat, row, i + 1);
      }
      fcat->file_first_id = fcat->first_id;
      
      fits_update_key (fcat->fptr, TINT, NCM_FIT_CATALOG_FIRST_ID_LABEL, &fcat->file_first_id, "Id of the first element.", &status);
      NCM_FITS_ERROR (status);
    }
    else if (fcat->file_first_id < fcat->first_id)
    {
      guint rows_to_add = fcat->first_id - fcat->file_first_id;
      GPtrArray *rows = g_ptr_array_new ();
      g_ptr_array_set_size (rows, rows_to_add);
      for (i = 0; i < rows_to_add; i++)
      {
        NcmVector *row = ncm_vector_dup (ncm_stats_vec_peek_x (fcat->pstats));
        _ncm_fit_catalog_read_row (fcat, row, i + 1);
        g_ptr_array_index (rows, i) = row;
      }
      ncm_stats_vec_prepend_data (fcat->pstats, rows, FALSE);
      g_ptr_array_unref (rows);
      fcat->first_id = fcat->file_first_id;
    }
    g_assert_cmpint (fcat->file_first_id, ==, fcat->first_id);
  }

  if (fcat->file_cur_id != fcat->cur_id)
  {
    if (fcat->file_cur_id < fcat->cur_id)
    {
      guint rows_to_add = fcat->cur_id - fcat->file_cur_id;
      guint offset = fcat->file_cur_id + 1 - fcat->file_first_id;

      for (i = 0; i < rows_to_add; i++)
      {
        NcmVector *row = ncm_stats_vec_peek_row (fcat->pstats, offset + i);
        _ncm_fit_catalog_write_row (fcat, row, offset + i + 1);
      }
      fcat->file_cur_id = fcat->cur_id;
    }
    else if (fcat->file_cur_id > fcat->cur_id)
    {
      guint rows_to_add = fcat->file_cur_id - fcat->cur_id;
      guint offset = fcat->cur_id + 1 - fcat->first_id;
      GPtrArray *rows = g_ptr_array_new ();
      g_ptr_array_set_size (rows, rows_to_add);
      for (i = 0; i < rows_to_add; i++)
      {
        NcmVector *row = ncm_vector_dup (ncm_stats_vec_peek_x (fcat->pstats));
        _ncm_fit_catalog_read_row (fcat, row, offset + i + 1);
        g_ptr_array_index (rows, i) = row;
      }
      ncm_stats_vec_append_data (fcat->pstats, rows, FALSE);
      g_ptr_array_unref (rows);
      fcat->cur_id = fcat->file_cur_id;
    } 
  }

  switch (fcat->fmode)
  {
    case NCM_FIT_CATALOG_FLUSH_DISABLE:
      break;
    case NCM_FIT_CATALOG_FLUSH_AUTO:
      _ncm_fit_catalog_flush_file (fcat);
      break;
    case NCM_FIT_CATALOG_FLUSH_TIMED:
    {
      if (g_timer_elapsed (fcat->flush_timer, NULL) > fcat->flush_interval)
      {
        g_timer_start (fcat->flush_timer);
        _ncm_fit_catalog_flush_file (fcat);        
      }
      break;
    }
    default:
      g_assert_not_reached ();
      break;
  }
#endif /* NUMCOSMO_HAVE_CFITSIO */
}

/**
 * ncm_fit_catalog_reset:
 * @fcat: a #NcmFitCatalog
 * 
 * Clean all catalog data, close any opened file. Otherwise it does 
 * not change any object's parameter.
 * 
 */
void
ncm_fit_catalog_reset (NcmFitCatalog *fcat)
{  
  ncm_fit_catalog_sync (fcat, TRUE);
  _ncm_fit_catalog_flush_file (fcat);
  _ncm_fit_catalog_close_file (fcat);
  
  ncm_stats_vec_reset (fcat->pstats);
  ncm_vector_set_all (fcat->params_max, GSL_NEGINF);
  ncm_vector_set_all (fcat->params_min, GSL_POSINF);
  fcat->first_id  =  0;
  fcat->cur_id    = -1;

  _ncm_fit_catalog_open_create_file (fcat);
}

/**
 * ncm_fit_catalog_erase_data:
 * @fcat: a #NcmFitCatalog
 * 
 * This function erases all data from the fits file associated with the
 * catalog.
 * 
 */
void
ncm_fit_catalog_erase_data (NcmFitCatalog *fcat)
{
  gint status = 0;
  gint nrows = fcat->file_cur_id - fcat->file_first_id + 1; 
  g_assert (fcat->fptr != NULL);

  if (nrows > 0)
  {    
    fits_delete_rows (fcat->fptr, 1, nrows, &status);
    NCM_FITS_ERROR (status);

    fcat->file_cur_id = fcat->file_first_id - 1;
    _ncm_fit_catalog_flush_file (fcat);
  }
}

/**
 * ncm_fit_catalog_peek_filename:
 * @fcat: a #NcmFitCatalog
 * 
 * Gets the filename associated with @fcat.
 * 
 * Returns: (transfer none): filename or NULL.
 */
const gchar *
ncm_fit_catalog_peek_filename (NcmFitCatalog *fcat)
{
  return fcat->file;
}

/**
 * ncm_fit_catalog_is_empty:
 * @fcat: a #NcmFitCatalog
 * 
 * Returns: TRUE when the catalog is empty;
 */
gboolean 
ncm_fit_catalog_is_empty (NcmFitCatalog *fcat)
{
  return (fcat->cur_id < fcat->first_id) && 
    (fcat->file_cur_id < fcat->file_first_id);
}

/**
 * ncm_fit_catalog_largest_error:
 * @fcat: a #NcmFitCatalog
 * 
 * Calculates the largest proportional error of the parameters included, i.e., $\text{lre} = \sigma_{\hat{p}}/(|\hat{p}|\sqrt{n})$ 
 * where $n$ represents the number of samples in the catalog, $\hat{p}$ is the estimated mean of the parameter $p$
 * and $\sigma_{\hat{p}}$ its standard deviation. 
 * 
 * It tries to guess when $p = 0$. In this case $\sigma_{\hat{p}} \approx |\hat{p}|\sqrt{n}$, therefore for $n > 10$ it tests
 * if $\text{lre} \approx 1$ and if it is the case it returns $\text{lre} = \sigma_{\hat{p}}/\sqrt{n}$ instead.
 * 
 * Returns: $\text{lre}$.
 */
gdouble
ncm_fit_catalog_largest_error (NcmFitCatalog *fcat)
{
  guint i;
  guint free_params_len = ncm_mset_fparams_len (fcat->fit->mset);
  const gdouble n = ncm_stats_vec_get_weight (fcat->pstats);
  const gdouble sqrt_n = sqrt (n);
  gdouble lerror = 0.0;
  
  if (n < 10)
  {
    for (i = 0; i < free_params_len; i++)
    {
      const gdouble mu = ncm_stats_vec_get_mean (fcat->pstats, i);
      const gdouble sd = ncm_stats_vec_get_sd (fcat->pstats, i);
      gdouble lerror_i = fabs (sd / (mu * sqrt_n));
      lerror = GSL_MAX (lerror, lerror_i);
    }
  }
  else
  {
    for (i = 0; i < free_params_len; i++)
    {
      const gdouble mu = ncm_stats_vec_get_mean (fcat->pstats, i);
      const gdouble sd = ncm_stats_vec_get_sd (fcat->pstats, i);
      gdouble lerror_i = fabs (sd / (mu * sqrt_n));
      guint lerror_i_truc = lerror_i; 
      if (lerror_i_truc == 1)
        lerror_i = fabs (sd / sqrt_n);
      lerror = GSL_MAX (lerror, lerror_i);
    }    
  }
  return lerror;
}

/**
 * ncm_fit_catalog_len:
 * @fcat: a #NcmFitCatalog
 * 
 * Number of itens in the catalog.
 * 
 * Returns: number of itens in the catalog.
 */
guint
ncm_fit_catalog_len (NcmFitCatalog *fcat)
{
  return fcat->pstats->nitens;
}

/**
 * ncm_fit_catalog_get_prng:
 * @fcat: a #NcmFitCatalog
 * @prng_algo: (out callee-allocates) (transfer full) (allow-none): the name of the PRNG algorithm.
 * @seed: (out caller-allocates) (allow-none): the seed used to initialize the PRNG.
 * 
 * This function checks if any pseudo random number generator (PRNG) is registred in the catalog
 * and returns its properties in @prng_algo and @seed if so.
 * 
 * Returns: TRUE if a PRNG is defined in the catalog.
 */
gboolean
ncm_fit_catalog_get_prng (NcmFitCatalog *fcat, gchar **prng_algo, gulong *seed)
{
  if (fcat->file != NULL)
  {
    gint status = 0;
    gchar algo[FLEN_VALUE];
    gboolean has_algo = FALSE;
    gboolean has_seed = FALSE;

    g_assert (*prng_algo == NULL);
    g_assert (seed != NULL);

    fits_read_key (fcat->fptr, TSTRING, NCM_FIT_CATALOG_PRNG_ALGO_LABEL,
                   &algo, NULL, &status);
    if (status == 0)
      has_algo = TRUE;
    else if (status == KEY_NO_EXIST)
      status = 0;
    else
      NCM_FITS_ERROR (status);

    fits_read_key (fcat->fptr, TULONG, NCM_FIT_CATALOG_PRNG_SEED_LABEL,
                   seed, NULL, &status);
    if (status == 0)
      has_seed = TRUE;
    else if (status == KEY_NO_EXIST)
      status = 0;
    else
      NCM_FITS_ERROR (status);

    if (has_seed != has_algo)
      g_error ("ncm_fit_catalog_get_prng: catalog contains algorithm name but no seed or vice-versa, has algorithm: %s, has seed: %s",
               has_algo ? "true" : "false", has_seed ? "true" : "false");

    if (has_algo)
      *prng_algo = g_strdup (algo);

    return has_algo;
  }
  else
    return FALSE;
}

/**
 * ncm_fit_catalog_set_prng:
 * @fcat: a #NcmFitCatalog
 * @rng: the #NcmRNG used.
 * 
 * This function sets the pseudo random number generator (PRNG) description in the catalog.
 * 
 */
void
ncm_fit_catalog_set_prng (NcmFitCatalog *fcat, NcmRNG *rng)
{
  if (fcat->file != NULL)
  {
    gint status = 0;
    gulong seed = ncm_rng_get_seed (rng);

    fits_update_key (fcat->fptr, TSTRING, NCM_FIT_CATALOG_PRNG_ALGO_LABEL, (gchar *)ncm_rng_get_algo (rng), "PRNG Algorithm name.", &status);
    NCM_FITS_ERROR (status);  

    fits_update_key (fcat->fptr, TULONG, NCM_FIT_CATALOG_PRNG_SEED_LABEL, &seed, "PRNG Algorithm seed.", &status);
    NCM_FITS_ERROR (status);  

    fits_update_key (fcat->fptr, TINT, NCM_FIT_CATALOG_FIRST_ID_LABEL, &fcat->file_first_id, "Id of the first element.", &status);
    NCM_FITS_ERROR (status);

    fits_flush_file (fcat->fptr, &status);
    NCM_FITS_ERROR (status);
  }
}

static void
_ncm_fit_catalog_post_update (NcmFitCatalog *fcat)
{
  guint i;
  guint len = fcat->pstats->len;
  NcmVector *x = ncm_stats_vec_peek_x (fcat->pstats);

  for (i = 0; i < len; i++)
  {
    gdouble p_i = ncm_vector_get (x, i);
    gdouble cur_max_p_i = ncm_vector_get (fcat->params_max, i);
    gdouble cur_min_p_i = ncm_vector_get (fcat->params_min, i);

    ncm_vector_set (fcat->params_max, i, GSL_MAX (p_i, cur_max_p_i));
    ncm_vector_set (fcat->params_min, i, GSL_MIN (p_i, cur_min_p_i));
  }
  
  ncm_stats_vec_update (fcat->pstats);

  fcat->cur_id++;

  ncm_fit_catalog_sync (fcat, FALSE);
}

/**
 * ncm_fit_catalog_add_from_fit:
 * @fcat: a #NcmFitCatalog
 * @fit: a #NcmFit
 * 
 * Adds a new element to the catalog using the parameters and likelihood value
 * from @fit. It assumes that @fit is compatible with the catalog.
 * 
 */
void
ncm_fit_catalog_add_from_fit (NcmFitCatalog *fcat, NcmFit *fit)
{
  NcmVector *row_i = ncm_stats_vec_peek_x (fcat->pstats);
  const gdouble m2lnL = ncm_fit_state_get_m2lnL_curval (fit->fstate);

  ncm_vector_set (row_i, 0, m2lnL);
  ncm_mset_fparams_get_vector_offset (fit->mset, row_i, 1);

  _ncm_fit_catalog_post_update (fcat);
}

/**
 * ncm_fit_catalog_add_from_vector:
 * @fcat: a #NcmFitCatalog
 * @vals: a #NcmVector
 * 
 * Adds a new element to the catalog using the values from the vector
 * @vals.
 * 
 */
void
ncm_fit_catalog_add_from_vector (NcmFitCatalog *fcat, NcmVector *vals)
{
  NcmVector *row_i = ncm_stats_vec_peek_x (fcat->pstats);
  ncm_vector_memcpy (row_i, vals);
  _ncm_fit_catalog_post_update (fcat);
}

static gdouble 
_fvar (gdouble v_i, gpointer user_data)
{
  NcmStatsVec *pstats = NCM_STATS_VEC (user_data);
  return sqrt (v_i * pstats->bias_wt);
}

static gdouble 
_fmeanvar (gdouble v_i, gpointer user_data)
{
  NcmStatsVec *pstats = NCM_STATS_VEC (user_data);
  return sqrt (v_i * pstats->bias_wt / pstats->nitens);
}

/**
 * ncm_fit_catalog_log_current_stats:
 * @fcat: a #NcmFitCatalog
 * 
 * Logs the current means and standard deviations of the catalog's parameters.
 * 
 */
void
ncm_fit_catalog_log_current_stats (NcmFitCatalog *fcat)
{
  ncm_vector_log_vals (fcat->pstats->mean,     "# NcmFitCatalog: Current mean:  ", "% -12.5g");
  ncm_vector_log_vals_func (fcat->pstats->var, "# NcmFitCatalog: Current msd :  ", "% -12.5g", &_fmeanvar, fcat->pstats);
  ncm_vector_log_vals_func (fcat->pstats->var, "# NcmFitCatalog: Current sd  :  ", "% -12.5g", &_fvar, fcat->pstats);
  ncm_vector_log_vals_avpb (fcat->pstats->var, "# NcmFitCatalog: Current var :  ", "% -12.5g", fcat->pstats->bias_wt, 0.0);
}

/**
 * ncm_fit_catalog_set_fit_mean_covar:
 * @fcat: a #NcmFitCatalog
 *
 * Sets the #NcmFit object parameters using the means from @fcat and
 * the covariance with the covariance of @fcat.
 * 
 */
void
ncm_fit_catalog_set_fit_mean_covar (NcmFitCatalog *fcat)
{
  ncm_stats_vec_get_cov_matrix (fcat->pstats, fcat->fit->fstate->covar, 1);
  ncm_stats_vec_get_mean_vector (fcat->pstats, fcat->fit->fstate->fparams, 1);
  ncm_mset_fparams_set_vector (fcat->fit->mset, fcat->fit->fstate->fparams);
  fcat->fit->fstate->has_covar = TRUE;
}

/**
 * ncm_fit_catalog_param_pdf:
 * @fcat: a #NcmFitCatalog
 * @i: parameter index.
 * 
 * Bins and calculates the pdf associated with the parameter @i.
 * (not ready yet FIXME)
 * 
 */
void
ncm_fit_catalog_param_pdf (NcmFitCatalog *fcat, guint i)
{
  const guint n = fcat->pstats->nitens;
  const guint nbins = n / 10 >= 10 ? n / 10 : 10;
  const gdouble p_max = ncm_vector_get (fcat->params_max, i);
  const gdouble p_min = ncm_vector_get (fcat->params_min, i);
  guint k;

  fcat->pdf_i = i;
	
  if (fcat->h != NULL && fcat->h->n != nbins)
  {
    gsl_histogram_free (fcat->h);
    fcat->h = NULL;
  }
  if (fcat->h == NULL)
    fcat->h = gsl_histogram_alloc (nbins);

  if (fcat->h_pdf != NULL && fcat->h_pdf->n != nbins)
  {
    gsl_histogram_pdf_free (fcat->h_pdf);
    fcat->h_pdf = NULL;
  }
  if (fcat->h_pdf == NULL)
    fcat->h_pdf = gsl_histogram_pdf_alloc (nbins);

  gsl_histogram_set_ranges_uniform (fcat->h, p_min, p_max);

  for (k = 0; k < fcat->pstats->nitens; k++)
  {
    NcmVector *row = ncm_stats_vec_peek_row (fcat->pstats, k);
    gsl_histogram_increment (fcat->h, ncm_vector_get (row, i));
  }

  gsl_histogram_pdf_init (fcat->h_pdf, fcat->h);
}

/**
 * ncm_fit_catalog_param_pdf_pvalue:
 * @fcat: a #NcmFitCatalog
 * @pval: parameter value
 * @both: one or both sides p-value.
 *
 * Calculates the p-value associated with the parameter value @pval.
 * 
 * Returns: the p-value.
 */
gdouble 
ncm_fit_catalog_param_pdf_pvalue (NcmFitCatalog *fcat, gdouble pval, gboolean both)
{
  g_assert_cmpint (fcat->pdf_i, >=, 0);
  g_assert (fcat->h_pdf != NULL);

  {
    const gdouble p_max = ncm_vector_get (fcat->params_max, fcat->pdf_i);
    const gdouble p_min = ncm_vector_get (fcat->params_min, fcat->pdf_i);
    gsize i = 0;

    NCM_UNUSED (both);
    if (pval < p_min || pval > p_max)
    {
      /*
       g_message ("# NcmFitCatalog: value % 20.15g outside mc obtained interval [% 20.15g % 20.15g]. Assuming 0 pvalue.",
       m2lnL, mc->fcat->m2lnL_min, mc->fcat->m2lnL_max);
       */
      return 0.0;
    }
    gsl_histogram_find (fcat->h, pval, &i);
    g_assert_cmpint (i, <=, fcat->h_pdf->n);
    if (i == 0)
      return 1.0;
    else
      return (1.0 - fcat->h_pdf->sum[i - 1]);
  }
}
