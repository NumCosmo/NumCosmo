/***************************************************************************
 *            ncm_mset_catalog.c
 *
 *  Tue February 18 10:49:26 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_mset_catalog.c
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
 * SECTION:ncm_mset_catalog
 * @title: Fit Catalog
 * @short_description: Ordered catalog of different #NcmMSet parameter values.
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

#include "math/ncm_mset_catalog.h"
#include "math/ncm_cfg.h"
#include "math/ncm_func_eval.h"
#include "ncm_enum_types.h"

G_DEFINE_TYPE (NcmMSetCatalog, ncm_mset_catalog, G_TYPE_OBJECT);

enum
{
  PROP_0,
  PROP_MSET,
  PROP_NADD_VALS,
  PROP_WEIGHTED,
  PROP_RNG,
  PROP_FILE,
  PROP_RUN_TYPE_STR,
  PROP_FLUSH_MODE,
  PROP_FLUSH_INTERVAL,
};

static void
ncm_mset_catalog_init (NcmMSetCatalog *mcat)
{
  mcat->mset           = NULL;
  mcat->nadd_vals      = 0;
  mcat->add_vals_names = g_ptr_array_new_with_free_func (g_free);
  mcat->pstats         = NULL;
  mcat->fmode          = NCM_MSET_CATALOG_FLUSH_LEN;
  mcat->rng            = NULL;
  mcat->weighted       = FALSE;
  mcat->first_flush    = FALSE;
  mcat->rng_inis       = NULL;
  mcat->rng_stat       = NULL;
  mcat->flush_timer    = g_timer_new ();
  mcat->cur_id         = -1; /* Represents that no elements in the catalog. */
  mcat->first_id       = 0;  /* The element to be in the catalog will be the one with index == 0 */
  mcat->file_cur_id    = -1; /* Represents that no elements in the catalog file. */
  mcat->file_first_id  = 0;  /* The element to be in the catalog file will be the one with index == 0 */
  mcat->file           = NULL;
  mcat->rtype_str      = NULL;
  mcat->porder         = g_array_new (FALSE, FALSE, sizeof (gint));
#ifdef NUMCOSMO_HAVE_CFITSIO
  mcat->fptr           = NULL;
#endif /* NUMCOSMO_HAVE_CFITSIO */
  mcat->pdf_i          = -1;
  mcat->h              = NULL;
  mcat->h_pdf          = NULL;
  mcat->params_max     = NULL;
  mcat->params_min     = NULL;
}

static void
_ncm_mset_catalog_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_mset_catalog_parent_class)->constructed (object);
  {
    NcmMSetCatalog *mcat = NCM_MSET_CATALOG (object);
    guint free_params_len = ncm_mset_fparams_len (mcat->mset);
    guint total = free_params_len + mcat->nadd_vals + (mcat->weighted ? 1 : 0);
    guint i;

    mcat->pstats     = ncm_stats_vec_new (total, NCM_STATS_VEC_COV, TRUE);
    mcat->params_max = ncm_vector_new (total);
    mcat->params_min = ncm_vector_new (total);

    ncm_vector_set_all (mcat->params_max, GSL_NEGINF);
    ncm_vector_set_all (mcat->params_min, GSL_POSINF);

    g_array_set_size (mcat->porder, total);
    g_ptr_array_set_size (mcat->add_vals_names, 0);
    
    for (i = 0; i < mcat->nadd_vals; i++)
    {
      g_ptr_array_add (mcat->add_vals_names, g_strdup_printf ("NcmMSetCatalog:additional-param-%u", i + 1));
    }

    if (mcat->weighted)
    {
      g_ptr_array_add (mcat->add_vals_names, g_strdup ("NcmMSetCatalog:Row-weights"));
      mcat->nadd_vals++;
    }
  }
}

static void
_ncm_mset_catalog_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmMSetCatalog *mcat = NCM_MSET_CATALOG (object);
  g_return_if_fail (NCM_IS_MSET_CATALOG (object));

  switch (prop_id)
  {
    case PROP_MSET:
      mcat->mset = g_value_dup_object (value);
      break;
    case PROP_NADD_VALS:
      mcat->nadd_vals = g_value_get_uint (value);
      break;
    case PROP_WEIGHTED:
      mcat->weighted = g_value_get_boolean (value);
      break;
    case PROP_RNG:
      ncm_mset_catalog_set_rng (mcat, g_value_get_object (value));
      break;
    case PROP_FILE:
      ncm_mset_catalog_set_file (mcat, g_value_get_string (value));
      break;    
    case PROP_RUN_TYPE_STR:
      ncm_mset_catalog_set_run_type (mcat, g_value_get_string (value));
      break;    
    case PROP_FLUSH_MODE:
      ncm_mset_catalog_set_flush_mode (mcat, g_value_get_enum (value));
      break;
    case PROP_FLUSH_INTERVAL:
      ncm_mset_catalog_set_flush_interval (mcat, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_mset_catalog_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmMSetCatalog *mcat = NCM_MSET_CATALOG (object);
  g_return_if_fail (NCM_IS_MSET_CATALOG (object));

  switch (prop_id)
  {
    case PROP_MSET:
      g_value_set_object (value, mcat->mset);
      break;
    case PROP_NADD_VALS:
      g_value_set_uint (value, mcat->nadd_vals);
      break;
    case PROP_WEIGHTED:
      g_value_set_boolean (value, mcat->weighted);
      break;
    case PROP_RNG:
      g_value_set_object (value, mcat->rng);
      break;
    case PROP_FILE:
      g_value_set_string (value, mcat->file);
      break;
    case PROP_RUN_TYPE_STR:
      g_value_set_string (value, mcat->rtype_str);
      break;
    case PROP_FLUSH_MODE:
      g_value_set_enum (value, mcat->fmode);
      break;
    case PROP_FLUSH_INTERVAL:
      g_value_set_double (value, mcat->flush_interval);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}
static void
_ncm_mset_catalog_dispose (GObject *object)
{
  NcmMSetCatalog *mcat = NCM_MSET_CATALOG (object);

  ncm_mset_clear (&mcat->mset);
  ncm_rng_clear (&mcat->rng);
  ncm_stats_vec_clear (&mcat->pstats);
  ncm_vector_free (mcat->params_max);
  ncm_vector_free (mcat->params_min);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_catalog_parent_class)->dispose (object);
}

static void _ncm_mset_catalog_close_file (NcmMSetCatalog *mcat);

static void
_ncm_mset_catalog_finalize (GObject *object)
{
  NcmMSetCatalog *mcat = NCM_MSET_CATALOG (object);

  if (mcat->h != NULL)
    gsl_histogram_free (mcat->h);
  if (mcat->h_pdf != NULL)
    gsl_histogram_pdf_free (mcat->h_pdf);

  _ncm_mset_catalog_close_file (mcat);

  g_clear_pointer (&mcat->rtype_str, g_free);
  
  g_array_unref (mcat->porder);
  g_timer_destroy (mcat->flush_timer);

  g_ptr_array_unref (mcat->add_vals_names);

  g_clear_pointer (&mcat->rng_inis, g_free);
  g_clear_pointer (&mcat->rng_stat, g_free);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_catalog_parent_class)->finalize (object);
}

static void
ncm_mset_catalog_class_init (NcmMSetCatalogClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->constructed  = &_ncm_mset_catalog_constructed;
  object_class->set_property = &_ncm_mset_catalog_set_property;
  object_class->get_property = &_ncm_mset_catalog_get_property;
  object_class->dispose      = &_ncm_mset_catalog_dispose;
  object_class->finalize     = &_ncm_mset_catalog_finalize;

  g_object_class_install_property (object_class,
                                   PROP_MSET,
                                   g_param_spec_object ("mset",
                                                        NULL,
                                                        "NcmMSet object",
                                                        NCM_TYPE_MSET,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_NADD_VALS,
                                   g_param_spec_uint ("nadd-vals",
                                                      NULL,
                                                      "Number of additional values",
                                                      1, G_MAXUINT32, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  g_object_class_install_property (object_class,
                                   PROP_WEIGHTED,
                                   g_param_spec_boolean ("weighted",
                                                         NULL,
                                                         "Catalog with weighted rows",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  g_object_class_install_property (object_class,
                                   PROP_RNG,
                                   g_param_spec_object ("rng",
                                                        NULL,
                                                        "Random number generator object",
                                                        NCM_TYPE_RNG,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_FLUSH_MODE,
                                   g_param_spec_enum ("fmode",
                                                      NULL,
                                                      "Catalog flush mode",
                                                      NCM_TYPE_MSET_CATALOG_FLUSH, NCM_MSET_CATALOG_FLUSH_AUTO,
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
                                                        NCM_MSET_CATALOG_RTYPE_UNDEFINED,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_FLUSH_INTERVAL,
                                   g_param_spec_double ("flush-interval",
                                                        NULL,
                                                        "Data flush interval",
                                                        0.0, 1.0e3, 10.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

}

/**
 * ncm_mset_catalog_new:
 * @mset: a #NcmMSet.
 * @nadd_vals: number of additional values.
 * @weighted: set to TRUE whenever the catalog is weighted.
 * @...: additional values names.
 *
 * Creates a new #NcmMSetCatalog based on the #NcmFit object @fit. The catalog assumes that
 * the @fit object will remain with the same set of free parameters during his whole lifetime.
 *
 * Returns: (transfer full): a new #NcmMSetCatalog
 */
NcmMSetCatalog *
ncm_mset_catalog_new (NcmMSet *mset, guint nadd_vals, gboolean weighted, ...)
{
  va_list ap;
  NcmMSetCatalog *mcat = g_object_new (NCM_TYPE_MSET_CATALOG, 
                                       "mset", mset,
                                       "nadd-vals", nadd_vals,
                                       "weighted", weighted,
                                       NULL);
  guint i;

  va_start (ap, weighted);
  
  for (i = 0; i < nadd_vals; i++)
  {
    gchar *name = va_arg (ap, gchar *);
    g_assert (name != NULL);
    ncm_mset_catalog_set_add_val_name (mcat, i, name);
  }
  
  va_end (ap);

  return mcat;
}

/**
 * ncm_mset_catalog_free:
 * @mcat: a #NcmMSetCatalog
 *
 * Decrese the reference count of @mcat atomically.
 *
 */
void 
ncm_mset_catalog_free (NcmMSetCatalog *mcat)
{
  ncm_mset_catalog_sync (mcat, TRUE);  
  g_object_unref (mcat);
}

/**
 * ncm_mset_catalog_clear:
 * @mcat: a #NcmMSetCatalog
 *
 * Decrese the reference count of *@mcat atomically and sets the pointer *@mcat to null.
 *
 */
void 
ncm_mset_catalog_clear (NcmMSetCatalog **mcat)
{
  ncm_mset_catalog_sync (*mcat, TRUE);  
  g_clear_object (mcat);
}

static void
_ncm_mset_catalog_open_create_file (NcmMSetCatalog *mcat)
{
  guint fparam_len = ncm_mset_fparam_len (mcat->mset);
  gchar key_text[FLEN_VALUE];
  gint status = 0;
  guint i;
  
  g_assert (mcat->file != NULL);
  g_assert (mcat->fptr == NULL);

  if (g_file_test (mcat->file, G_FILE_TEST_EXISTS))
  {
    glong nrows;

    fits_open_file (&mcat->fptr, mcat->file, READWRITE, &status);
    NCM_FITS_ERROR (status);

    fits_movnam_hdu (mcat->fptr, BINARY_TBL, NCM_MSET_CATALOG_EXTNAME, 0, &status);
    NCM_FITS_ERROR (status);

    fits_read_key (mcat->fptr, TINT, NCM_MSET_CATALOG_FIRST_ID_LABEL,
                   &mcat->file_first_id, NULL, &status);
    NCM_FITS_ERROR (status);

    fits_read_key (mcat->fptr, TSTRING, NCM_MSET_CATALOG_RTYPE_LABEL,
                   key_text, NULL, &status);
    NCM_FITS_ERROR (status);
    
    if (strcmp (mcat->rtype_str, key_text) != 0)
      g_error ("_ncm_mset_catalog_open_create_file: incompatible run type strings from catalog and file, catalog: `%s' file: `%s'.",
               mcat->rtype_str, key_text);

    fits_get_num_rows (mcat->fptr, &nrows, &status);
    NCM_FITS_ERROR (status);

    if (mcat->file_first_id != mcat->first_id)
    {
      if (nrows == 0)
      {
        if (mcat->file_first_id != 0)
          g_warning ("_ncm_mset_catalog_open_create_file: Empty data file with "NCM_MSET_CATALOG_FIRST_ID_LABEL" different from first_id: %d != %d. Setting to first_id.\n", 
                     mcat->file_first_id, mcat->first_id);
        mcat->file_first_id = mcat->first_id;
      } 
      else if (mcat->cur_id < mcat->first_id)
      {
        if (mcat->first_id != 0)
          g_warning ("_ncm_mset_catalog_open_create_file: Empty memory catalog with first_id different from "NCM_MSET_CATALOG_FIRST_ID_LABEL": %d != %d. Setting to "NCM_MSET_CATALOG_FIRST_ID_LABEL".\n", 
                     mcat->first_id, mcat->file_first_id);
        mcat->first_id = mcat->file_first_id;
        mcat->cur_id = mcat->file_first_id - 1;
      }
    }
    mcat->file_cur_id = mcat->file_first_id + nrows - 1; 

    for (i = 0; i < mcat->nadd_vals; i++)
    {
      const gchar *cname = g_ptr_array_index (mcat->add_vals_names, i);
      gint cindex = 0;
      
      if (fits_get_colnum (mcat->fptr, CASESEN, (gchar *)cname, &cindex, &status))
        g_error ("_ncm_mset_catalog_open_create_file: Column %s not found, invalid fits file.", cname);

      if (cindex != i + 1)
      {
        g_error ("_ncm_mset_catalog_open_create_file: Column %s is not the %d-th column [%d], invalid fits file.", 
                 cname, i + 1, cindex);
      }
      g_array_index (mcat->porder, gint, i) = cindex;  
    }

    for (i = 0; i < fparam_len; i++)
    {
      const gchar *fparam_fullname = ncm_mset_fparam_full_name (mcat->mset, i);
      if (fits_get_colnum (mcat->fptr, CASESEN, (gchar *)fparam_fullname, &g_array_index (mcat->porder, gint, i + mcat->nadd_vals), &status))
      {                /* I don't like this too ^^^^^^^^^  */
        g_error ("_ncm_mset_catalog_open_create_file: Column %s not found, invalid fits file.", fparam_fullname);
      }
    }    
  }
  else
  {
    GPtrArray *ttype_array = g_ptr_array_sized_new (10);
    GPtrArray *tform_array = g_ptr_array_sized_new (10);
    fits_create_file (&mcat->fptr, mcat->file, &status);
    NCM_FITS_ERROR (status);

    for (i = 0; i < mcat->nadd_vals; i++)
    {
      const gchar *cname = g_ptr_array_index (mcat->add_vals_names, i);
      g_ptr_array_add (ttype_array, (gchar *)cname);
      g_ptr_array_add (tform_array, "1D");
      g_array_index (mcat->porder, gint, i) = tform_array->len;
    }

    for (i = 0; i < fparam_len; i++)
    {
      const gchar *fparam_fullname = ncm_mset_fparam_full_name (mcat->mset, i);
      g_ptr_array_add (ttype_array, (gchar *)fparam_fullname);
      /* I don't like this too ^^^^^^^^^  */
      g_ptr_array_add (tform_array, "1D");
      g_array_index (mcat->porder, gint, i + mcat->nadd_vals) = tform_array->len;
    }

    /* append a new empty binary table onto the FITS file */
    fits_create_tbl (mcat->fptr, BINARY_TBL, 0, fparam_len + mcat->nadd_vals, (gchar **)ttype_array->pdata, (gchar **)tform_array->pdata,
                     NULL, NCM_MSET_CATALOG_EXTNAME, &status);
    NCM_FITS_ERROR (status);

    fits_update_key (mcat->fptr, TSTRING, NCM_MSET_CATALOG_RTYPE_LABEL, mcat->rtype_str, "Run type string.", &status);
    NCM_FITS_ERROR (status);  

    for (i = 0; i < fparam_len; i++)
    {
      gchar *fsymbi = g_strdup_printf ("%s%d", NCM_MSET_CATALOG_FSYMB_LABEL, i + 1);
      gchar *fsymb_desc = g_strdup_printf ("Symbol for parameter %s[%d]", 
                                           ncm_mset_fparam_name (mcat->mset, i), 
                                           i + 1);
      const gchar *fsymb  = ncm_mset_fparam_symbol (mcat->mset, i);
      
      fits_update_key (mcat->fptr, TSTRING, fsymbi, (gchar *)fsymb, fsymb_desc, &status);
      NCM_FITS_ERROR (status);  
      g_free (fsymbi);
      g_free (fsymb_desc);
    }
    
    mcat->file_first_id = 0;
    mcat->file_cur_id   = -1;
  }

  fits_read_key (mcat->fptr, TSTRING, NCM_MSET_CATALOG_RNG_ALGO_LABEL, 
                 key_text, NULL, &status);
  if (status == 0)
  {
    gulong seed = 0;
    gchar *inis = NULL;
    fits_read_key (mcat->fptr, TULONG, NCM_MSET_CATALOG_RNG_SEED_LABEL, 
                   &seed, NULL, &status);
    NCM_FITS_ERROR (status);

    fits_read_key_longstr (mcat->fptr, NCM_MSET_CATALOG_RNG_INIS_LABEL, &inis, NULL, &status);
    NCM_FITS_ERROR (status);
        
    if (mcat->rng != NULL)
    {
      const gchar *cat_algo = ncm_rng_get_algo (mcat->rng);
      g_assert_cmpstr (cat_algo, ==, key_text);
      g_assert_cmpstr (inis, ==, mcat->rng_inis);
    }
    else
    {
      NcmRNG *rng = ncm_rng_new (key_text);
      ncm_rng_set_state (rng, inis);
      ncm_mset_catalog_set_rng (mcat, rng);
      ncm_rng_free (rng);
    }

    fits_free_memory (inis, &status);
    NCM_FITS_ERROR (status);
  }
  else if (status == KEY_NO_EXIST)
  {
    status = 0;
    if (mcat->rng != NULL)
    {
      gulong seed = ncm_rng_get_seed (mcat->rng);

      fits_update_key (mcat->fptr, TSTRING, NCM_MSET_CATALOG_RNG_ALGO_LABEL, (gchar *)ncm_rng_get_algo (mcat->rng), "RNG Algorithm name.", &status);
      NCM_FITS_ERROR (status);  

      fits_update_key (mcat->fptr, TULONG, NCM_MSET_CATALOG_RNG_SEED_LABEL, &seed, "RNG Algorithm seed.", &status);
      NCM_FITS_ERROR (status);  

      fits_update_key_longstr (mcat->fptr, NCM_MSET_CATALOG_RNG_INIS_LABEL, mcat->rng_inis, NULL, &status);
      NCM_FITS_ERROR (status);
    }
  }
  else
    NCM_FITS_ERROR (status);

  fits_update_key (mcat->fptr, TINT, NCM_MSET_CATALOG_FIRST_ID_LABEL, &mcat->file_first_id, "Id of the first element.", &status);
  NCM_FITS_ERROR (status);

  fits_flush_file (mcat->fptr, &status);
  NCM_FITS_ERROR (status);
}

static void _ncm_mset_catalog_flush_file (NcmMSetCatalog *mcat);

/**
 * ncm_mset_catalog_set_add_val_name:
 * @mcat: a #NcmMSetCatalog.
 * @i: index of the additional value.
 * @name: additional value name.
 *
 * Sets the @i-th additional value name.
 *
 */
void 
ncm_mset_catalog_set_add_val_name (NcmMSetCatalog *mcat, guint i, const gchar *name)
{
  g_free (g_ptr_array_index (mcat->add_vals_names, i));
  g_ptr_array_index (mcat->add_vals_names, i) = g_strdup (name);
}

/**
 * ncm_mset_catalog_set_file:
 * @mcat: a #NcmMSetCatalog
 * @filename: a filename
 *
 * Sets the data filename to be used to sync/save data.
 *
 */
void 
ncm_mset_catalog_set_file (NcmMSetCatalog *mcat, const gchar *filename)
{
  gint status = 0;
  if (mcat->file != NULL && strcmp (mcat->file, filename) == 0)
    return;
  
  _ncm_mset_catalog_close_file (mcat);

  g_clear_pointer (&mcat->file, g_free);
  mcat->file = g_strdup (filename);  

  _ncm_mset_catalog_open_create_file (mcat);
  ncm_mset_catalog_sync (mcat, TRUE);

  _ncm_mset_catalog_flush_file (mcat);

  fits_flush_file (mcat->fptr, &status);
  NCM_FITS_ERROR (status);

  mcat->first_flush = TRUE;
}

/**
 * ncm_mset_catalog_set_flush_mode:
 * @mcat: a #NcmMSetCatalog
 * @fmode: flush mode
 *
 * Sets the flush mode to @fmode.
 *
 */
void 
ncm_mset_catalog_set_flush_mode (NcmMSetCatalog *mcat, NcmMSetCatalogFlush fmode)
{
  mcat->fmode = fmode;
}

/**
 * ncm_mset_catalog_set_flush_interval:
 * @mcat: a #NcmMSetCatalog
 * @interval: Minimum time interval between flushs
 *
 * Sets the minimum time interval between flushs.
 *
 */
void 
ncm_mset_catalog_set_flush_interval (NcmMSetCatalog *mcat, gdouble interval)
{
  mcat->flush_interval = interval;
}

/**
 * ncm_mset_catalog_set_first_id:
 * @mcat: a #NcmMSetCatalog
 * @first_id: the id of the first item in the sample.
 * 
 * Sets the first id of the catalog, mainly used to inform in which realization the catalog starts. 
 * 
 */
void 
ncm_mset_catalog_set_first_id (NcmMSetCatalog *mcat, gint first_id)
{
  if (first_id == mcat->first_id)
    return;

  g_assert_cmpint (mcat->file_first_id, ==, mcat->first_id);
  g_assert_cmpint (mcat->file_cur_id, ==, mcat->cur_id);

  if (mcat->cur_id != mcat->first_id - 1)
    g_error ("ncm_mset_catalog_set_first_id: cannot modify first_id to %d in a non-empty catalog, catalog first id: %d, catalog current id: %d.", 
             first_id, mcat->first_id, mcat->cur_id);

  mcat->first_id = first_id;
  mcat->cur_id   = first_id - 1;

  mcat->file_first_id = first_id;
  mcat->file_cur_id   = first_id - 1;

  if (mcat->fptr != NULL)
  {
    gint status = 0;
    fits_update_key (mcat->fptr, TINT, NCM_MSET_CATALOG_FIRST_ID_LABEL, &mcat->file_first_id, "Id of the first element.", &status);
    NCM_FITS_ERROR (status);
    
    ncm_mset_catalog_sync (mcat, TRUE);
  }
}

/**
 * ncm_mset_catalog_set_run_type:
 * @mcat: a #NcmMSetCatalog
 * @rtype_str: the run type string.
 * 
 * Sets the run type string. 
 * 
 */
void 
ncm_mset_catalog_set_run_type (NcmMSetCatalog *mcat, const gchar *rtype_str)
{
  gint status = 0;
  g_assert (rtype_str != NULL);
  
  if (mcat->rtype_str != NULL)
  {
    if (strcmp (mcat->rtype_str, rtype_str) != 0)
    {
      if (mcat->cur_id + 1 != mcat->first_id)
        g_error ("ncm_mset_catalog_set_run_type: cannot change run type string in a non-empty catalog, actual: `%s' new: `%s'.",
                 mcat->rtype_str, rtype_str);
      else
      {
        g_clear_pointer (&mcat->rtype_str, g_free);
      }
    }
    else
      return;
  }

  mcat->rtype_str = g_strdup (rtype_str);
  if (mcat->fptr != NULL)
  {
    fits_update_key (mcat->fptr, TSTRING, NCM_MSET_CATALOG_RTYPE_LABEL, mcat->rtype_str, NULL, &status);
    NCM_FITS_ERROR (status);
  }
}

/**
 * ncm_mset_catalog_set_rng:
 * @mcat: a #NcmMSetCatalog
 * @rng: a #NcmRNG.
 * 
 * Sets the random number generator. 
 * 
 */
void 
ncm_mset_catalog_set_rng (NcmMSetCatalog *mcat, NcmRNG *rng)
{
  g_assert (rng != NULL);
  
  if (mcat->rng != NULL)
  {
    g_error ("ncm_mset_catalog_set_rng: random number generator alread set.");
  }

  if (mcat->cur_id + 1 != mcat->first_id)
    g_warning ("ncm_mset_catalog_set_rng: setting RNG in a non-empty catalog, catalog first id: %d, catalog current id: %d.",
             mcat->first_id, mcat->cur_id);

  mcat->rng = ncm_rng_ref (rng);

  mcat->rng_inis = ncm_rng_get_state (rng);
  mcat->rng_stat = g_strdup (mcat->rng_inis);

  if (mcat->file != NULL)
  {
    gint status = 0;
    gulong seed = ncm_rng_get_seed (rng);

    fits_update_key (mcat->fptr, TSTRING, NCM_MSET_CATALOG_RNG_ALGO_LABEL, (gchar *)ncm_rng_get_algo (rng), "RNG Algorithm name.", &status);
    NCM_FITS_ERROR (status);  

    fits_update_key (mcat->fptr, TULONG, NCM_MSET_CATALOG_RNG_SEED_LABEL, &seed, "RNG Algorithm seed.", &status);
    NCM_FITS_ERROR (status);  

    fits_update_key (mcat->fptr, TINT, NCM_MSET_CATALOG_FIRST_ID_LABEL, &mcat->file_first_id, "Id of the first element.", &status);
    NCM_FITS_ERROR (status);

    fits_update_key_longstr (mcat->fptr, NCM_MSET_CATALOG_RNG_INIS_LABEL, mcat->rng_inis, NULL, &status);
    NCM_FITS_ERROR (status);

    fits_flush_file (mcat->fptr, &status);
    NCM_FITS_ERROR (status);
  }
}

static void
_ncm_mset_catalog_flush_file (NcmMSetCatalog *mcat)
{
  gint status = 0;
  gint64 nrows = mcat->file_cur_id - mcat->file_first_id + 1;

  fits_update_key (mcat->fptr, TLONGLONG, NCM_MSET_CATALOG_NROWS_LABEL, &nrows, NULL, &status);
  NCM_FITS_ERROR (status);

  if (mcat->rng != NULL)
  {
    g_clear_pointer (&mcat->rng_stat, g_free);
    mcat->rng_stat = ncm_rng_get_state (mcat->rng);

    fits_update_key_longstr (mcat->fptr, NCM_MSET_CATALOG_RNG_STAT_LABEL, mcat->rng_stat, NULL, &status);
    NCM_FITS_ERROR (status);
  }

  if (G_UNLIKELY (mcat->first_flush))
  {
    fits_flush_file (mcat->fptr, &status);
    NCM_FITS_ERROR (status);
    mcat->first_flush = FALSE;
  }
  else
  {
    fits_flush_buffer (mcat->fptr, 0, &status);
    NCM_FITS_ERROR (status);
  }
}

static void
_ncm_mset_catalog_close_file (NcmMSetCatalog *mcat)
{
  gint status = 0;
#ifdef NUMCOSMO_HAVE_CFITSIO
  if (mcat->fptr != NULL)
  {
    _ncm_mset_catalog_flush_file (mcat);
    
    fits_close_file (mcat->fptr, &status);
    NCM_FITS_ERROR (status);

    mcat->fptr = NULL;
  }
#endif /* NUMCOSMO_HAVE_CFITSIO */
}

static void 
_ncm_mset_catalog_write_row (NcmMSetCatalog *mcat, NcmVector *row, guint row_index)
{
  guint i;
  gint status = 0;
  
  for (i = 0; i < ncm_vector_len (row); i++)
  {
    fits_write_col (mcat->fptr, TDOUBLE, g_array_index (mcat->porder, gint, i), row_index, 
                    1, 1, ncm_vector_ptr (row, i), &status);
    NCM_FITS_ERROR (status);
  }
}

static void 
_ncm_mset_catalog_read_row (NcmMSetCatalog *mcat, NcmVector *row, guint row_index)
{
  guint i;
  gint status = 0;
    
  for (i = 0; i < ncm_vector_len (row); i++)
  {
    fits_read_col (mcat->fptr, TDOUBLE, g_array_index (mcat->porder, gint, i), row_index, 
                   1, 1, NULL, ncm_vector_ptr (row, i), NULL, &status);
    NCM_FITS_ERROR (status);
  }
}

/**
 * ncm_mset_catalog_sync:
 * @mcat: a #NcmMSetCatalog
 * @check: whether to check consistence between file and memory data.
 * 
 * Syncronise memory and data file. If no file was defined simply returns.
 * 
 */
void
ncm_mset_catalog_sync (NcmMSetCatalog *mcat, gboolean check)
{
  gint status = 0;
  guint i;

  if (mcat->file == NULL)
    return;

#ifdef NUMCOSMO_HAVE_CFITSIO
  g_assert (mcat->fptr != NULL);
  
  if (check)
  {
    gchar fptr_filename[FLEN_FILENAME];
  
    fits_file_name (mcat->fptr, fptr_filename, &status);
    NCM_FITS_ERROR (status);

    g_assert_cmpstr (fptr_filename, ==, mcat->file);

    if ((mcat->file_cur_id < mcat->first_id - 1) || (mcat->cur_id < mcat->file_first_id - 1))
      g_error ("ncm_mset_catalog_sync: file data & catalog mismatch, they do not intersect each other: file data [%d, %d] catalog [%d, %d]",
               mcat->file_first_id, mcat->file_cur_id,
               mcat->first_id, mcat->cur_id);
  }

  if (mcat->file_first_id != mcat->first_id)
  {
    if (mcat->file_first_id > mcat->first_id)
    {
      guint rows_to_add = mcat->file_first_id - mcat->first_id;
      fits_insert_rows (mcat->fptr, 0, rows_to_add, &status);
      NCM_FITS_ERROR (status);

      for (i = 0; i < rows_to_add; i++)
      {
        NcmVector *row = ncm_stats_vec_peek_row (mcat->pstats, i);
        _ncm_mset_catalog_write_row (mcat, row, i + 1);
      }
      mcat->file_first_id = mcat->first_id;

      if (mcat->rng != NULL)
      {
        fits_update_key_longstr (mcat->fptr, NCM_MSET_CATALOG_RNG_INIS_LABEL, mcat->rng_inis, NULL, &status);
        NCM_FITS_ERROR (status);
      }
      
      fits_update_key (mcat->fptr, TINT, NCM_MSET_CATALOG_FIRST_ID_LABEL, &mcat->file_first_id, "Id of the first element.", &status);
      NCM_FITS_ERROR (status);
    }
    else if (mcat->file_first_id < mcat->first_id)
    {
      guint rows_to_add = mcat->first_id - mcat->file_first_id;
      GPtrArray *rows = g_ptr_array_new ();
      gchar *inis = NULL;
      
      g_ptr_array_set_size (rows, rows_to_add);
      for (i = 0; i < rows_to_add; i++)
      {
        NcmVector *row = ncm_vector_dup (ncm_stats_vec_peek_x (mcat->pstats));
        _ncm_mset_catalog_read_row (mcat, row, i + 1);
        g_ptr_array_index (rows, i) = row;
      }
      ncm_stats_vec_prepend_data (mcat->pstats, rows, FALSE);
      g_ptr_array_unref (rows);
      mcat->first_id = mcat->file_first_id;

      if (mcat->rng != NULL)
      {
        fits_read_key_longstr (mcat->fptr, NCM_MSET_CATALOG_RNG_INIS_LABEL, &inis, NULL, &status);
        NCM_FITS_ERROR (status);

        g_clear_pointer (&mcat->rng_inis, g_free);
        mcat->rng_inis = g_strdup (inis);

        fits_free_memory (inis, &status);
        NCM_FITS_ERROR (status);
      }
    }
    g_assert_cmpint (mcat->file_first_id, ==, mcat->first_id);
  }

  if (mcat->file_cur_id != mcat->cur_id)
  {
    if (mcat->file_cur_id < mcat->cur_id)
    {
      guint rows_to_add = mcat->cur_id - mcat->file_cur_id;
      guint offset = mcat->file_cur_id + 1 - mcat->file_first_id;

      for (i = 0; i < rows_to_add; i++)
      {
        NcmVector *row = ncm_stats_vec_peek_row (mcat->pstats, offset + i);
        _ncm_mset_catalog_write_row (mcat, row, offset + i + 1);
      }
      mcat->file_cur_id = mcat->cur_id;

      if (mcat->rng != NULL)
      {
        g_clear_pointer (&mcat->rng_stat, g_free);
        mcat->rng_stat = ncm_rng_get_state (mcat->rng);

        fits_update_key_longstr (mcat->fptr, NCM_MSET_CATALOG_RNG_STAT_LABEL, mcat->rng_stat, NULL, &status);
        NCM_FITS_ERROR (status);
      }
    }
    else if (mcat->file_cur_id > mcat->cur_id)
    {
      guint rows_to_add = mcat->file_cur_id - mcat->cur_id;
      guint offset = mcat->cur_id + 1 - mcat->first_id;
      GPtrArray *rows = g_ptr_array_new ();
      gchar *stat = NULL;
      
      g_ptr_array_set_size (rows, rows_to_add);

      for (i = 0; i < rows_to_add; i++)
      {
        NcmVector *row = ncm_vector_dup (ncm_stats_vec_peek_x (mcat->pstats));
        _ncm_mset_catalog_read_row (mcat, row, offset + i + 1);
        g_ptr_array_index (rows, i) = row;
      }
      ncm_stats_vec_append_data (mcat->pstats, rows, FALSE);
      mcat->cur_id = mcat->file_cur_id;
      
      g_ptr_array_unref (rows);
      
      if (mcat->rng != NULL)
      {
        fits_read_key_longstr (mcat->fptr, NCM_MSET_CATALOG_RNG_STAT_LABEL, &stat, NULL, &status);
        NCM_FITS_ERROR (status);

        g_clear_pointer (&mcat->rng_stat, g_free);
        mcat->rng_stat = g_strdup (stat);

        fits_free_memory (stat, &status);
        NCM_FITS_ERROR (status);

        ncm_rng_set_state (mcat->rng, mcat->rng_stat);
      }
    } 
  }

  switch (mcat->fmode)
  {
    case NCM_MSET_CATALOG_FLUSH_DISABLE:
      break;
    case NCM_MSET_CATALOG_FLUSH_AUTO:
      _ncm_mset_catalog_flush_file (mcat);
      break;
    case NCM_MSET_CATALOG_FLUSH_TIMED:
    {
      if (g_timer_elapsed (mcat->flush_timer, NULL) > mcat->flush_interval)
      {
        g_timer_start (mcat->flush_timer);
        _ncm_mset_catalog_flush_file (mcat);        
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
 * ncm_mset_catalog_reset:
 * @mcat: a #NcmMSetCatalog
 * 
 * Clean all catalog data from memory and file. Otherwise it does 
 * not change any object's parameter.
 * 
 */
void
ncm_mset_catalog_reset (NcmMSetCatalog *mcat)
{  
  ncm_mset_catalog_erase_data (mcat);
  
  ncm_stats_vec_reset (mcat->pstats);
  ncm_vector_set_all (mcat->params_max, GSL_NEGINF);
  ncm_vector_set_all (mcat->params_min, GSL_POSINF);

  mcat->cur_id    = mcat->first_id - 1;
  _ncm_mset_catalog_close_file (mcat);
}

/**
 * ncm_mset_catalog_erase_data:
 * @mcat: a #NcmMSetCatalog
 * 
 * This function erases all data from the fits file associated with the
 * catalog.
 * 
 */
void
ncm_mset_catalog_erase_data (NcmMSetCatalog *mcat)
{
  if (mcat->fptr != NULL)
  {
    gint status = 0;
    gint nrows = mcat->file_cur_id - mcat->file_first_id + 1; 

    if (nrows > 0)
    {    
      fits_delete_rows (mcat->fptr, 1, nrows, &status);
      NCM_FITS_ERROR (status);

      mcat->file_cur_id = mcat->file_first_id - 1;
      _ncm_mset_catalog_flush_file (mcat);
    }
  }
}

/**
 * ncm_mset_catalog_peek_filename:
 * @mcat: a #NcmMSetCatalog
 * 
 * Gets the filename associated with @mcat.
 * 
 * Returns: (transfer none): filename or NULL.
 */
const gchar *
ncm_mset_catalog_peek_filename (NcmMSetCatalog *mcat)
{
  return mcat->file;
}

/**
 * ncm_mset_catalog_is_empty:
 * @mcat: a #NcmMSetCatalog
 * 
 * Returns: TRUE when the catalog is empty;
 */
gboolean 
ncm_mset_catalog_is_empty (NcmMSetCatalog *mcat)
{
  return (mcat->cur_id < mcat->first_id) && 
    (mcat->file_cur_id < mcat->file_first_id);
}

/**
 * ncm_mset_catalog_largest_error:
 * @mcat: a #NcmMSetCatalog
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
ncm_mset_catalog_largest_error (NcmMSetCatalog *mcat)
{
  guint i;
  guint free_params_len = ncm_mset_fparams_len (mcat->mset);
  const gdouble n = ncm_stats_vec_get_weight (mcat->pstats);
  const gdouble sqrt_n = sqrt (n);
  const gdouble fpi = mcat->nadd_vals;
  const gdouble fpf = free_params_len + mcat->nadd_vals;
  gdouble lerror = 0.0;
  
  if (n < 10)
  {
    for (i = fpi; i < fpf; i++)
    {
      const gdouble mu = ncm_stats_vec_get_mean (mcat->pstats, i);
      const gdouble sd = ncm_stats_vec_get_sd (mcat->pstats, i);
      gdouble lerror_i = fabs (sd / (mu * sqrt_n));
      lerror = GSL_MAX (lerror, lerror_i);
    }
  }
  else
  {
    for (i = fpi; i < fpf; i++)
    {
      const gdouble mu = ncm_stats_vec_get_mean (mcat->pstats, i);
      const gdouble sd = ncm_stats_vec_get_sd (mcat->pstats, i);
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
 * ncm_mset_catalog_len:
 * @mcat: a #NcmMSetCatalog
 * 
 * Number of itens in the catalog.
 * 
 * Returns: number of itens in the catalog.
 */
guint
ncm_mset_catalog_len (NcmMSetCatalog *mcat)
{
  return mcat->pstats->nitens;
}

/**
 * ncm_mset_catalog_get_rng:
 * @mcat: a #NcmMSetCatalog
 * @rng_algo: (out callee-allocates) (transfer full) (allow-none): the name of the RNG algorithm.
 * @seed: (out caller-allocates) (allow-none): the seed used to initialize the RNG.
 * 
 * This function checks if any pseudo random number generator (RNG) is registred in the catalog
 * and returns its properties in @rng_algo and @seed if so.
 * 
 * Returns: TRUE if a RNG is defined in the catalog.
 */
gboolean
ncm_mset_catalog_get_rng (NcmMSetCatalog *mcat, gchar **rng_algo, gulong *seed)
{
  if (mcat->file != NULL)
  {
    gint status = 0;
    gchar algo[FLEN_VALUE];
    gboolean has_algo = FALSE;
    gboolean has_seed = FALSE;

    g_assert (*rng_algo == NULL);
    g_assert (seed != NULL);

    fits_read_key (mcat->fptr, TSTRING, NCM_MSET_CATALOG_RNG_ALGO_LABEL,
                   &algo, NULL, &status);
    if (status == 0)
      has_algo = TRUE;
    else if (status == KEY_NO_EXIST)
      status = 0;
    else
      NCM_FITS_ERROR (status);

    fits_read_key (mcat->fptr, TULONG, NCM_MSET_CATALOG_RNG_SEED_LABEL,
                   seed, NULL, &status);
    if (status == 0)
      has_seed = TRUE;
    else if (status == KEY_NO_EXIST)
      status = 0;
    else
      NCM_FITS_ERROR (status);

    if (has_seed != has_algo)
      g_error ("ncm_mset_catalog_get_rng: catalog contains algorithm name but no seed or vice-versa, has algorithm: %s, has seed: %s",
               has_algo ? "true" : "false", has_seed ? "true" : "false");

    if (has_algo)
      *rng_algo = g_strdup (algo);

    return has_algo;
  }
  else
    return FALSE;
}

static void
_ncm_mset_catalog_post_update (NcmMSetCatalog *mcat)
{
  guint i;
  guint len = mcat->pstats->len;
  NcmVector *x = ncm_stats_vec_peek_x (mcat->pstats);

  for (i = 0; i < len; i++)
  {
    gdouble p_i = ncm_vector_get (x, i);
    gdouble cur_max_p_i = ncm_vector_get (mcat->params_max, i);
    gdouble cur_min_p_i = ncm_vector_get (mcat->params_min, i);

    ncm_vector_set (mcat->params_max, i, GSL_MAX (p_i, cur_max_p_i));
    ncm_vector_set (mcat->params_min, i, GSL_MIN (p_i, cur_min_p_i));
  }

  if (mcat->weighted)
    ncm_stats_vec_update_weight (mcat->pstats, ncm_vector_get (x, mcat->nadd_vals - 1));
  else  
    ncm_stats_vec_update (mcat->pstats);

  mcat->cur_id++;

  ncm_mset_catalog_sync (mcat, FALSE);
}

/**
 * ncm_mset_catalog_add_from_mset:
 * @mcat: a #NcmMSetCatalog
 * @mset: a #NcmMSet
 * @...: additional values.
 * 
 * Adds a new element to the catalog using the parameters from @mset. 
 * It assumes that @mset is compatible with the catalog and expect the
 * right number of additional values.
 * 
 */
void
ncm_mset_catalog_add_from_mset (NcmMSetCatalog *mcat, NcmMSet *mset, ...)
{
  NcmVector *row_i = ncm_stats_vec_peek_x (mcat->pstats);
  va_list ap;
  guint i;

  va_start (ap, mset);
  
  for (i = 0; i < mcat->nadd_vals; i++)
  {
    gdouble val = va_arg (ap, gdouble);
    ncm_vector_set (row_i, i, val);
  }
  
  va_end (ap);

  ncm_mset_fparams_get_vector_offset (mset, row_i, mcat->nadd_vals);

  _ncm_mset_catalog_post_update (mcat);
}

/**
 * ncm_mset_catalog_add_from_mset_array:
 * @mcat: a #NcmMSetCatalog
 * @mset: a #NcmMSet
 * @ax: (array) (element-type double): additional values array.
 * 
 * Adds a new element to the catalog using the parameters from @mset. 
 * It assumes that @mset is compatible with the catalog and expect the
 * right number of additional values in the array @ax.
 * 
 */
void
ncm_mset_catalog_add_from_mset_array (NcmMSetCatalog *mcat, NcmMSet *mset, gdouble *ax)
{
  NcmVector *row_i = ncm_stats_vec_peek_x (mcat->pstats);
  guint i;
  
  for (i = 0; i < mcat->nadd_vals; i++)
    ncm_vector_set (row_i, i, ax[i]);

  ncm_mset_fparams_get_vector_offset (mset, row_i, mcat->nadd_vals);

  _ncm_mset_catalog_post_update (mcat);
}

/**
 * ncm_mset_catalog_add_from_vector:
 * @mcat: a #NcmMSetCatalog
 * @vals: a #NcmVector
 * 
 * Adds a new element to the catalog using the values from the vector
 * @vals.
 * 
 */
void
ncm_mset_catalog_add_from_vector (NcmMSetCatalog *mcat, NcmVector *vals)
{
  NcmVector *row_i = ncm_stats_vec_peek_x (mcat->pstats);
  ncm_vector_memcpy (row_i, vals);
  _ncm_mset_catalog_post_update (mcat);
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
 * ncm_mset_catalog_log_current_stats:
 * @mcat: a #NcmMSetCatalog
 * 
 * Logs the current means and standard deviations of the catalog's parameters.
 * 
 */
void
ncm_mset_catalog_log_current_stats (NcmMSetCatalog *mcat)
{
  ncm_vector_log_vals (mcat->pstats->mean,     "# NcmMSetCatalog: Current mean:  ", "% -12.5g");
  ncm_vector_log_vals_func (mcat->pstats->var, "# NcmMSetCatalog: Current msd :  ", "% -12.5g", &_fmeanvar, mcat->pstats);
  ncm_vector_log_vals_func (mcat->pstats->var, "# NcmMSetCatalog: Current sd  :  ", "% -12.5g", &_fvar, mcat->pstats);
  ncm_vector_log_vals_avpb (mcat->pstats->var, "# NcmMSetCatalog: Current var :  ", "% -12.5g", mcat->pstats->bias_wt, 0.0);
}

/**
 * ncm_mset_catalog_peek_row:
 * @mcat: a #NcmMSetCatalog
 * @i: the row index.
 * 
 * Gets the @i-th row.
 * 
 * Returns: (transfer none): the row with index @i or NULL if not available. 
 */
NcmVector *
ncm_mset_catalog_peek_row (NcmMSetCatalog *mcat, guint i)
{
  if (mcat->pstats->nitens == 0 || i >= mcat->pstats->nitens)
    return NULL;
  else
    return ncm_stats_vec_peek_row (mcat->pstats, i);
}

/**
 * ncm_mset_catalog_peek_current_row:
 * @mcat: a #NcmMSetCatalog
 * 
 * Gets the last added row.
 * 
 * Returns: (transfer none): last added row or NULL if empty. 
 */
NcmVector *
ncm_mset_catalog_peek_current_row (NcmMSetCatalog *mcat)
{
  if (mcat->pstats->nitens == 0)
    return NULL;
  else
    return ncm_stats_vec_peek_row (mcat->pstats, mcat->pstats->nitens - 1);
}

/**
 * ncm_mset_catalog_get_mean:
 * @mcat: a #NcmMSetCatalog
 * @mean: (inout) (allow-none) (transfer full): a #NcmVector
 *
 * Gets the current mean vector.
 * 
 */
void
ncm_mset_catalog_get_mean (NcmMSetCatalog *mcat, NcmVector **mean)
{
  if (*mean == NULL)
    *mean = ncm_vector_new (mcat->pstats->len - mcat->nadd_vals);
  ncm_stats_vec_get_mean_vector (mcat->pstats, *mean, mcat->nadd_vals);
}

/**
 * ncm_mset_catalog_get_covar:
 * @mcat: a #NcmMSetCatalog
 * @cov: (inout) (allow-none) (transfer full): a #NcmMatrix
 *
 * Gets the current covariance matrix.
 * 
 */
void
ncm_mset_catalog_get_covar (NcmMSetCatalog *mcat, NcmMatrix **cov)
{
  if (*cov == NULL)
    *cov = ncm_matrix_new (mcat->pstats->len - mcat->nadd_vals, mcat->pstats->len - mcat->nadd_vals);
  ncm_stats_vec_get_cov_matrix (mcat->pstats, *cov, mcat->nadd_vals);
}

/**
 * ncm_mset_catalog_param_pdf:
 * @mcat: a #NcmMSetCatalog
 * @i: parameter index.
 * 
 * Bins and calculates the pdf associated with the parameter @i.
 * (not ready yet FIXME)
 * 
 */
void
ncm_mset_catalog_param_pdf (NcmMSetCatalog *mcat, guint i)
{
  const guint n = mcat->pstats->nitens;
  const guint nbins = n / 10 >= 10 ? n / 10 : 10;
  const gdouble p_max = ncm_vector_get (mcat->params_max, i);
  const gdouble p_min = ncm_vector_get (mcat->params_min, i);
  guint k;

  mcat->pdf_i = i;
	
  if (mcat->h != NULL && mcat->h->n != nbins)
  {
    gsl_histogram_free (mcat->h);
    mcat->h = NULL;
  }
  if (mcat->h == NULL)
    mcat->h = gsl_histogram_alloc (nbins);

  if (mcat->h_pdf != NULL && mcat->h_pdf->n != nbins)
  {
    gsl_histogram_pdf_free (mcat->h_pdf);
    mcat->h_pdf = NULL;
  }
  if (mcat->h_pdf == NULL)
    mcat->h_pdf = gsl_histogram_pdf_alloc (nbins);

  gsl_histogram_set_ranges_uniform (mcat->h, p_min, p_max);

  for (k = 0; k < mcat->pstats->nitens; k++)
  {
    NcmVector *row = ncm_stats_vec_peek_row (mcat->pstats, k);
    gsl_histogram_increment (mcat->h, ncm_vector_get (row, i));
  }

  gsl_histogram_pdf_init (mcat->h_pdf, mcat->h);
}

/**
 * ncm_mset_catalog_param_pdf_pvalue:
 * @mcat: a #NcmMSetCatalog
 * @pval: parameter value
 * @both: one or both sides p-value.
 *
 * Calculates the p-value associated with the parameter value @pval.
 * 
 * Returns: the p-value.
 */
gdouble 
ncm_mset_catalog_param_pdf_pvalue (NcmMSetCatalog *mcat, gdouble pval, gboolean both)
{
  g_assert_cmpint (mcat->pdf_i, >=, 0);
  g_assert (mcat->h_pdf != NULL);

  {
    const gdouble p_max = ncm_vector_get (mcat->params_max, mcat->pdf_i);
    const gdouble p_min = ncm_vector_get (mcat->params_min, mcat->pdf_i);
    gsize i = 0;

    NCM_UNUSED (both);
    if (pval < p_min || pval > p_max)
    {
      g_warning ("ncm_mset_catalog_param_pdf_pvalue: value % 20.15g outside mc obtained interval [% 20.15g % 20.15g]. Assuming 0 pvalue.",
                 pval, p_min, p_max);
      return 0.0;
    }
    gsl_histogram_find (mcat->h, pval, &i);
    g_assert_cmpint (i, <=, mcat->h_pdf->n);
    if (i == 0)
      return 1.0;
    else
      return (1.0 - mcat->h_pdf->sum[i - 1]);
  }
}
