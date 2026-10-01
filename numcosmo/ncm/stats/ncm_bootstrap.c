/***************************************************************************
 *            ncm_bootstrap.c
 *
 *  Fri August 16 11:09:01 2013
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_bootstrap.c
 * Copyright (C) 2013 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmBootstrap:
 *
 * Random samples of the indexes of a data set.
 *
 * A realization is an array of #NcmBootstrap:bootstrap-size indexes in [0,
 * #NcmBootstrap:full-size), drawn with replacement by ncm_bootstrap_resample() or without
 * replacement by ncm_bootstrap_remix(). Changing either size discards the realization.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/stats/ncm_bootstrap.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_c.h"

enum
{
  PROP_0,
  PROP_FSIZE,
  PROP_BSIZE,
  PROP_INIT,
  PROP_REAL,
};

struct _NcmBootstrap
{
  /*< private >*/
  GObject parent_instance;
  guint fsize;
  guint bsize;
  GArray *bootstrap_index;
  GArray *increasing_index;
  gboolean init;
};


G_DEFINE_TYPE (NcmBootstrap, ncm_bootstrap, G_TYPE_OBJECT)

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

        if (bsize != bstrap->bsize)
          g_error ("_ncm_bootstrap_set_property: realization has %u indexes, but the bootstrap size is %u.",
                   bsize, bstrap->bsize);

        for (i = 0; i < bsize; i++)
        {
          guint j = 0;

          g_variant_get_child (var, i, "u", &j);

          if (j >= bstrap->fsize)
            g_error ("_ncm_bootstrap_set_property: realization index %u is not smaller than the full size %u.",
                     j, bstrap->fsize);

          g_array_index (bstrap->bootstrap_index, guint, i) = j;
        }

        bstrap->init = TRUE;
      }

      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
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

#if GLIB_CHECK_VERSION (2, 68, 0)
        gpointer mem = g_memdup2 (bstrap->bootstrap_index->data, msize);

#else
        gpointer mem = g_memdup (bstrap->bootstrap_index->data, msize);
#endif /* GLIB_CHECK_VERSION(2,68,0) */
        var = g_variant_new_from_data (G_VARIANT_TYPE ("au"),
                                       mem, msize, TRUE, &g_free, mem);
      }
      else
      {
        var = g_variant_new ("au", NULL);
      }

      g_value_take_variant (value, var);
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
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
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

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
 * Creates a new #NcmBootstrap with both sizes zero.
 *
 * Returns: (transfer full): a new #NcmBootstrap.
 */
NcmBootstrap *
ncm_bootstrap_new (void)
{
  NcmBootstrap *bstrap = g_object_new (NCM_TYPE_BOOTSTRAP, NULL);

  return bstrap;
}

/**
 * ncm_bootstrap_sized_new:
 * @fsize: full sample size
 *
 * Creates a new #NcmBootstrap with #NcmBootstrap:full-size and
 * #NcmBootstrap:bootstrap-size both equal to @fsize.
 *
 * Returns: (transfer full): a new #NcmBootstrap.
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
 * @fsize: full sample size
 * @bsize: bootstrap size
 *
 * Creates a new #NcmBootstrap drawing @bsize indexes from [0, @fsize).
 *
 * Returns: (transfer full): a new #NcmBootstrap.
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
 * @bstrap: a #NcmBootstrap
 *
 * Increases the reference count of @bstrap by one.
 *
 * Returns: (transfer full): @bstrap.
 */
NcmBootstrap *
ncm_bootstrap_ref (NcmBootstrap *bstrap)
{
  return g_object_ref (bstrap);
}

/**
 * ncm_bootstrap_free:
 * @bstrap: a #NcmBootstrap
 *
 * Decreases the reference count of @bstrap by one.
 */
void
ncm_bootstrap_free (NcmBootstrap *bstrap)
{
  g_object_unref (bstrap);
}

/**
 * ncm_bootstrap_clear:
 * @bstrap: a #NcmBootstrap
 *
 * Decreases the reference count of *@bstrap by one and sets *@bstrap to %NULL.
 */
void
ncm_bootstrap_clear (NcmBootstrap **bstrap)
{
  g_clear_object (bstrap);
}

/**
 * ncm_bootstrap_set_fsize:
 * @bstrap: a #NcmBootstrap
 * @fsize: full sample size
 *
 * Sets #NcmBootstrap:full-size to @fsize. The bootstrap size is not changed. A new
 * size discards the realization.
 */
void
ncm_bootstrap_set_fsize (NcmBootstrap *bstrap, guint fsize)
{
  if (fsize != bstrap->fsize)
    bstrap->init = FALSE;

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
 * @bstrap: a #NcmBootstrap
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
 * @bstrap: a #NcmBootstrap
 * @bsize: bootstrap size
 *
 * Sets #NcmBootstrap:bootstrap-size to @bsize. A new size discards the realization.
 */
void
ncm_bootstrap_set_bsize (NcmBootstrap *bstrap, guint bsize)
{
  if (bsize != bstrap->bsize)
    bstrap->init = FALSE;

  bstrap->bsize = bsize;
  g_array_set_size (bstrap->bootstrap_index, bsize);
}

/**
 * ncm_bootstrap_get_bsize:
 * @bstrap: a #NcmBootstrap
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
 * @bstrap: a #NcmBootstrap
 * @rng: a #NcmRNG
 *
 * Draws a new realization of #NcmBootstrap:bootstrap-size indexes from [0,
 * #NcmBootstrap:full-size) with replacement. Locks @rng while drawing. Aborts if the
 * full size is zero and the bootstrap size is not.
 */
void
ncm_bootstrap_resample (NcmBootstrap *bstrap, NcmRNG *rng)
{
  gpointer bdata           = bstrap->bootstrap_index->data;
  gpointer idata           = bstrap->increasing_index->data;
  const gsize fsize        = bstrap->fsize;
  const gsize bsize        = bstrap->bsize;
  const gsize element_size = g_array_get_element_size (bstrap->bootstrap_index);

  if ((fsize == 0) && (bsize > 0))
    g_error ("ncm_bootstrap_resample: cannot draw %zu indexes from an empty sample.", bsize);

  ncm_rng_lock (rng);
  ncm_rng_sample (rng, bdata, bsize, idata, fsize, element_size);
  ncm_rng_unlock (rng);
  bstrap->init = TRUE;
}

/**
 * ncm_bootstrap_remix:
 * @bstrap: a #NcmBootstrap
 * @rng: a #NcmRNG
 *
 * Draws a new realization of #NcmBootstrap:bootstrap-size distinct indexes from [0,
 * #NcmBootstrap:full-size) without replacement, stored in increasing order. Locks @rng
 * while drawing. Aborts if the bootstrap size is larger than the full size.
 */
void
ncm_bootstrap_remix (NcmBootstrap *bstrap, NcmRNG *rng)
{
  gpointer bdata           = bstrap->bootstrap_index->data;
  gpointer idata           = bstrap->increasing_index->data;
  const gsize fsize        = bstrap->fsize;
  const gsize bsize        = bstrap->bsize;
  const gsize element_size = g_array_get_element_size (bstrap->bootstrap_index);

  if (bsize > fsize)
    g_error ("ncm_bootstrap_remix: cannot draw %zu distinct indexes from a sample of %zu.", bsize, fsize);

  ncm_rng_lock (rng);
  ncm_rng_choose (rng, bdata, bsize, idata, fsize, element_size);
  ncm_rng_unlock (rng);
  bstrap->init = TRUE;
}

/**
 * ncm_bootstrap_get:
 * @bstrap: a #NcmBootstrap
 * @i: position in the realization, in [0, #NcmBootstrap:bootstrap-size)
 *
 * Gets the index at position @i of the current realization. @i is not checked.
 *
 * Returns: the @i-th index of the realization.
 */
guint
ncm_bootstrap_get (NcmBootstrap *bstrap, guint i)
{
  return g_array_index (bstrap->bootstrap_index, guint, i);
}

static gint _ncm_bootstrap_get_sort (gconstpointer a, gconstpointer b);

/**
 * ncm_bootstrap_get_sortncomp:
 * @bstrap: a #NcmBootstrap
 *
 * Counts the distinct indexes of the current realization. The result holds the pairs
 * (index, number of occurrences), in increasing order of index; it is empty when the
 * bootstrap size is zero. The realization is not changed. Aborts if @bstrap has no
 * realization.
 *
 * Returns: (array) (element-type guint) (transfer full): the index and count pairs.
 */
GArray *
ncm_bootstrap_get_sortncomp (NcmBootstrap *bstrap)
{
  GArray *res     = g_array_sized_new (FALSE, TRUE, sizeof (guint), 2 * bstrap->bsize);
  GArray *sorted  = g_array_sized_new (FALSE, FALSE, sizeof (guint), bstrap->bsize);
  const guint one = 1;
  guint i, j, n_c;

  if (!bstrap->init)
    g_error ("ncm_bootstrap_get_sortncomp: the bootstrap has no realization, call ncm_bootstrap_resample() or ncm_bootstrap_remix() first.");

  if (bstrap->bsize == 0)
  {
    g_array_unref (sorted);

    return res;
  }

  g_array_append_vals (sorted, bstrap->bootstrap_index->data, bstrap->bsize);
  g_array_sort (sorted, &_ncm_bootstrap_get_sort);

  n_c = g_array_index (sorted, guint, 0);

  j = 0;
  g_array_append_val (res, n_c);
  g_array_append_val (res, one);

  for (i = 1; i < bstrap->bsize; i++)
  {
    const guint n_i = g_array_index (sorted, guint, i);

    if (n_i == n_c)
    {
      g_array_index (res, guint, 2 * j + 1)++;
    }
    else
    {
      g_array_append_val (res, n_i);
      g_array_append_val (res, one);
      n_c = n_i;
      j++;
    }
  }

  g_array_unref (sorted);

  return res;
}

/**
 * ncm_bootstrap_is_init:
 * @bstrap: a #NcmBootstrap
 *
 * Checks whether @bstrap holds a realization, drawn by ncm_bootstrap_resample() or
 * ncm_bootstrap_remix() or set through #NcmBootstrap:realization.
 *
 * Returns: %TRUE if @bstrap holds a realization.
 */
gboolean
ncm_bootstrap_is_init (NcmBootstrap *bstrap)
{
  return bstrap->init;
}

static gint
_ncm_bootstrap_get_sort (gconstpointer a, gconstpointer b)
{
  return (*((guint *) a) < *((guint *) b)) ? -1 : ((*((guint *) a) > *((guint *) b)) ? 1 : 0);
}

