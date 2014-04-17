/***************************************************************************
 *            ncm_bootstrap.c
 *
 *  Fri August 16 11:09:01 2013
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_bootstrap.c
 * Copyright (C) 2013 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_bootstrap
 * @title: Bootstrap object
 * @short_description: Generic index bootstrap.
 * 
 * This object generate random samples of indexes. These samples are used 
 * to calculate statistics using different combinations of the actual data set.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_bootstrap.h"
#include "math/ncm_cfg.h"
#include "math/ncm_c.h"

enum
{
  PROP_0,
  PROP_FSIZE,
  PROP_BSIZE,
  PROP_INIT,
  PROP_REAL,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcmBootstrap, ncm_bootstrap, G_TYPE_OBJECT);

static void
ncm_bootstrap_init (NcmBootstrap *bstrap)
{
  bstrap->fsize            = 0;
  bstrap->bsize            = 0;
  bstrap->bootstrap_index  = g_array_new (FALSE, FALSE, sizeof (guint));
  bstrap->increasing_index = g_array_new (FALSE, FALSE, sizeof (guint));
  bstrap->init             = FALSE;
}

static void
_ncm_bootstrap_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmBootstrap *bstrap = NCM_BOOTSTRAP (object);
  g_return_if_fail (NCM_IS_BOOTSTRAP (object));

  switch (prop_id)
  {
    case PROP_FSIZE:
      ncm_bootstrap_set_fsize (bstrap, g_value_get_uint (value));
      break;
    case PROP_BSIZE:
      ncm_bootstrap_set_bsize (bstrap, g_value_get_uint (value));
      break;
    case PROP_REAL:
    {
      GVariant *var = g_value_get_variant (value);
      guint bsize;
      g_assert (g_variant_is_of_type (var, G_VARIANT_TYPE ("au")));
      bsize = g_variant_n_children (var);

      if (bsize != 0)
      {
        guint i;
        g_assert_cmpuint (bsize, ==, bstrap->bsize);
        for (i = 0; i < bsize; i++)
        {
          guint j = 0;
          g_variant_get_child (var, i, "u", j);
          g_array_index (bstrap->bootstrap_index, guint, i) = j;
        }
        bstrap->init = TRUE;
      }
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_bootstrap_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmBootstrap *bstrap = NCM_BOOTSTRAP (object);
  g_return_if_fail (NCM_IS_BOOTSTRAP (object));

  switch (prop_id)
  {
    case PROP_FSIZE:
      g_value_set_uint (value, bstrap->fsize);
      break;
    case PROP_BSIZE:
      g_value_set_uint (value, bstrap->bsize);
      break;
    case PROP_INIT:
      g_value_set_boolean (value, bstrap->init);
      break;
    case PROP_REAL:
    {
      GVariant *var;
      if (bstrap->init)
      {
        gsize msize = sizeof (guint) * bstrap->bootstrap_index->len;
        gpointer mem = g_memdup (bstrap->bootstrap_index->data, msize);
        var = g_variant_new_from_data (G_VARIANT_TYPE ("au"),
                                       mem, msize, TRUE, &g_free, mem);
      }
      else
        var = g_variant_new ("au", NULL);
      
      g_value_take_variant (value, var);
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_bootstrap_finalize (GObject *object)
{
  NcmBootstrap *bstrap = NCM_BOOTSTRAP (object);
  
  if (bstrap->bootstrap_index != NULL)
  {
    g_array_unref (bstrap->bootstrap_index);
    bstrap->bootstrap_index = NULL;
  }

  if (bstrap->increasing_index != NULL)
  {
    g_array_unref (bstrap->increasing_index);
    bstrap->increasing_index = NULL;
  }

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_bootstrap_parent_class)->finalize (object);
}

static void
ncm_bootstrap_class_init (NcmBootstrapClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_bootstrap_set_property;
  object_class->get_property = &_ncm_bootstrap_get_property;
  object_class->finalize     = &_ncm_bootstrap_finalize;


  g_object_class_install_property (object_class,
                                   PROP_FSIZE,
                                   g_param_spec_uint ("full-size",
                                                      NULL,
                                                      "Data sample size",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_BSIZE,
                                   g_param_spec_uint ("bootstrap-size",
                                                      NULL,
                                                      "Bootstrap size",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_INIT,
                                   g_param_spec_boolean ("init",
                                                         NULL,
                                                         "Bootstrap initialization status",
                                                         FALSE,
                                                         G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_REAL,
                                   g_param_spec_variant ("realization",
                                                         NULL,
                                                         "Bootstrap current realization",
                                                         G_VARIANT_TYPE ("au"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_bootstrap_new:
 * 
 * Creates a new zero sized #NcmBootstrap object.
 * 
 * Returns: (transfer full): a #NcmBootstrap.
 */
NcmBootstrap *
ncm_bootstrap_new (void)
{
  NcmBootstrap *bstrap = g_object_new (NCM_TYPE_BOOTSTRAP, NULL);
  return bstrap;
}

/**
 * ncm_bootstrap_sized_new:
 * @fsize: sample size.
 * 
 * Creates a new #NcmBootstrap object for a sample of size @fsize. This object
 * will sample with replacement all indexes @fsize times.
 * 
 * Returns: (transfer full): a #NcmBootstrap.
 */
NcmBootstrap *
ncm_bootstrap_sized_new (guint fsize)
{
  NcmBootstrap *bstrap = g_object_new (NCM_TYPE_BOOTSTRAP, 
                                       "full-size", fsize, 
                                       "bootstrap-size", fsize,
                                       NULL);
  return bstrap;
}

/**
 * ncm_bootstrap_full_new:
 * @fsize: sample size.
 * @bsize: bootstrap size.
 * 
 * Creates a new #NcmBootstrap object for a sample of size @fsize. This object
 * will sample with replacement all indexes @bsize times.
 * 
 * Returns: (transfer full): a #NcmBootstrap.
 */
NcmBootstrap *
ncm_bootstrap_full_new (guint fsize, guint bsize)
{
  NcmBootstrap *bstrap = g_object_new (NCM_TYPE_BOOTSTRAP, 
                                       "full-size", fsize,
                                       "bootstrap-size", bsize, 
                                       NULL);
  return bstrap;
}

/**
 * ncm_bootstrap_ref:
 * @bstrap: a #NcmBootstrap.
 * 
 * Incresases the reference count of @bstrap by one.
 * 
 * Returns: (transfer full): a #NcmBootstrap.
 */
NcmBootstrap *
ncm_bootstrap_ref (NcmBootstrap *bstrap)
{
  return g_object_ref (bstrap);
}

/**
 * ncm_bootstrap_free:
 * @bstrap: a #NcmBootstrap.
 * 
 * Decreases the reference count of @bstrap by one.
 * 
 */
void 
ncm_bootstrap_free (NcmBootstrap *bstrap)
{
  g_object_unref (bstrap);
}

/**
 * ncm_bootstrap_clear:
 * @bstrap: a #NcmBootstrap.
 * 
 * Decreases the reference count of *@bstrap by one and sets *@bstrap tp NULL.
 * 
 */
void 
ncm_bootstrap_clear (NcmBootstrap **bstrap)
{
  g_clear_object (bstrap);
}

/**
 * ncm_bootstrap_set_fsize:
 * @bstrap: a #NcmBootstrap.
 * @fsize: full sample size. 
 * 
 * Sets the full sample size, it also sets the bsize to the same value @fsize.
 * 
 */
void 
ncm_bootstrap_set_fsize (NcmBootstrap *bstrap, guint fsize)
{
  g_array_set_size (bstrap->increasing_index, fsize);

  bstrap->fsize = fsize;
  
  if (fsize > 0)
  {
    guint i;
    for (i = 0; i < fsize; i++)
      g_array_index (bstrap->increasing_index, guint, i) = i;
  }
}

/**
 * ncm_bootstrap_get_fsize:
 * @bstrap: a #NcmBootstrap.
 * 
 * Gets the full sample size.
 * 
 * Returns: the full sample size.
 */
guint 
ncm_bootstrap_get_fsize (NcmBootstrap *bstrap)
{
  return bstrap->fsize;
}

/**
 * ncm_bootstrap_set_bsize:
 * @bstrap: a #NcmBootstrap.
 * @bsize: bootstrap size. 
 * 
 * Sets the bootstrap size.
 * 
 */
void 
ncm_bootstrap_set_bsize (NcmBootstrap *bstrap, guint bsize)
{
  bstrap->bsize = bsize;
  g_array_set_size (bstrap->bootstrap_index, bsize);
}

/**
 * ncm_bootstrap_get_bsize:
 * @bstrap: a #NcmBootstrap.
 * 
 * Gets the bootstrap size.
 * 
 * Returns: the bootstrap size.
 */
guint 
ncm_bootstrap_get_bsize (NcmBootstrap *bstrap)
{
  return bstrap->bsize;
}

/**
 * ncm_bootstrap_resample:
 * @bstrap: a #NcmBootstrap.
 * @rng: a #NcmRNG.
 * 
 * Sample with replacement #NcmBootstrap:bootstrap-size from the 
 * #NcmBootstrap:full-size indexes.
 * 
 */
/**
 * ncm_bootstrap_remix:
 * @bstrap: a #NcmBootstrap.
 * @rng: a #NcmRNG.
 * 
 * Sample without replacement #NcmBootstrap:bootstrap-size from the 
 * #NcmBootstrap:full-size indexes. Note that in this case
 * #NcmBootstrap:bootstrap-size must be equal or smaller than 
 * #NcmBootstrap:full-size.
 * 
 */
/**
 * ncm_bootstrap_get:
 * @bstrap: a #NcmBootstrap.
 * @i: index in [0, #NcmBootstrap:bootstrap-size - 1].
 * 
 * Gets the index associated with the @i-th resampled index.
 * 
 * Returns: the @i-th resampled index.
 */
/**
 * ncm_bootstrap_is_init:
 * @bstrap: a #NcmBootstrap.
 * 
 * Checks if the bootstrap object was initialized (remix or resample).
 * 
 * Returns: whether @bstrap is initialized.
 */
