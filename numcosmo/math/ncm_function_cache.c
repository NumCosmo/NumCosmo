/***************************************************************************
 *            ncm_function_cache.c
 *
 *  Wed Aug 13 21:18:48 2008
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
 * SECTION:ncm_function_cache
 * @title: NcmFunctionCache
 * @short_description: A generic cache for functions values
 * @stability: Stable
 * @include: numcosmo/math/ncm_function_cache.h
 *
 * A simple cache that saves function values at different argument values.
 * It can be used to find an already computed value or the value of the
 * function closest to an already computed point.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_function_cache.h"
#include "math/ncm_util.h"

struct _NcmFunctionCachePrivate
{
  GTree *tree;
  GMutex lock;
  gboolean clear;
  guint n;
  gdouble abstol;
  gdouble reltol;
};

enum
{
  PROP_0,
  PROP_DIM,
  PROP_RELTOL,
  PROP_ABSTOL,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmFunctionCache, ncm_function_cache, G_TYPE_OBJECT);

static gint gdouble_compare (gconstpointer a, gconstpointer b, gpointer user_data);
static void gdouble_free (gpointer data);
static gboolean cache_clean (NcmFunctionCache *cache);

static void
ncm_function_cache_init (NcmFunctionCache *cache)
{
  NcmFunctionCachePrivate * const self = cache->priv = ncm_function_cache_get_instance_private (cache);

  self->tree   = g_tree_new_full (&gdouble_compare, NULL, &gdouble_free, (GDestroyNotify)&ncm_vector_free);
  self->clear  = FALSE;
  self->n      = 0;
  self->abstol = 0.0;
  self->reltol = 0.0;

  g_mutex_init (&self->lock);
}

static void
_ncm_function_cache_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFunctionCache *cache = NCM_FUNCTION_CACHE (object);
  NcmFunctionCachePrivate * const self = cache->priv;
  g_return_if_fail (NCM_IS_FUNCTION_CACHE (object));

  switch (prop_id)
  {
    case PROP_DIM:
      self->n = g_value_get_uint (value);
      break;
    case PROP_RELTOL:
      self->reltol = g_value_get_double (value);
      break;
    case PROP_ABSTOL:
      self->abstol = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_function_cache_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFunctionCache *cache = NCM_FUNCTION_CACHE (object);
  NcmFunctionCachePrivate * const self = cache->priv;
  g_return_if_fail (NCM_IS_FUNCTION_CACHE (object));

  switch (prop_id)
  {
    case PROP_DIM:
      g_value_set_uint (value, self->n);
      break;
    case PROP_RELTOL:
      g_value_set_double (value, ncm_function_cache_get_reltol (cache));
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, ncm_function_cache_get_abstol (cache));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_function_cache_dispose (GObject *object)
{
  /*NcmFunctionCache *cache = NCM_FUNCTION_CACHE (object);*/
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_function_cache_parent_class)->dispose (object);
}

static void
_ncm_function_cache_finalize (GObject *object)
{
  /*NcmFunctionCache *cache = NCM_FUNCTION_CACHE (object);*/

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_function_cache_parent_class)->finalize (object);
}

static void
ncm_function_cache_class_init (NcmFunctionCacheClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_function_cache_set_property;
  object_class->get_property = &_ncm_function_cache_get_property;
  object_class->dispose      = &_ncm_function_cache_dispose;
  object_class->finalize     = &_ncm_function_cache_finalize;

  g_object_class_install_property (object_class,
                                   PROP_DIM,
                                   g_param_spec_uint ("dimension",
                                                      NULL,
                                                      "Function dimension",
                                                      1, G_MAXUINT, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                      NULL,
                                                      "Relative tolerance",
                                                      1.0e-15, 1.0, 1.0e-7,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ABSTOL,
                                   g_param_spec_double ("abstol",
                                                      NULL,
                                                      "Absolute tolerance",
                                                      0.0, G_MAXDOUBLE, 0.0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}


/**
 * ncm_function_cache_new:
 * @n: function dimension
 * @abstol: relative tolerance
 * @reltol: absolute tolerance
 *
 * Creates a new #NcmFunctionCache for a @n dimensional function.
 * The points are considered the same within the tolerance 
 * described by @abstol and @reltol.
 *
 * Returns: the newly created #NcmFunctionCache.
 */
NcmFunctionCache *
ncm_function_cache_new (guint n, gdouble abstol, gdouble reltol)
{
  NcmFunctionCache *cache = g_object_new (NCM_TYPE_FUNCTION_CACHE,
                                          "dimension", n,
                                          "reltol",    reltol,
                                          "abstol",    abstol,
                                          NULL);
  return cache;
}

/**
 * ncm_function_cache_ref:
 * @cache: a #NcmFunctionCache
 *
 * Increase the reference of @cache by one.
 *
 * Returns: (transfer full): @cache.
 */
NcmFunctionCache *
ncm_function_cache_ref (NcmFunctionCache *cache)
{
  return g_object_ref (cache);
}

/**
 * ncm_function_cache_free:
 * @cache: a #NcmFunctionCache
 *
 * Decrease the reference count of @cache by one.
 *
 */
void
ncm_function_cache_free (NcmFunctionCache *cache)
{
  g_object_unref (cache);
}

/**
 * ncm_function_cache_clear:
 * @cache: a #NcmFunctionCache
 *
 * Decrease the reference count of @cache by one, and sets the pointer *@cache to
 * NULL.
 *
 */
void
ncm_function_cache_clear (NcmFunctionCache **cache)
{
  g_clear_object (cache);
}

/**
 * ncm_function_cache_get_reltol:
 * @cache: a #NcmFunctionCache
 * 
 * Returns: the relative tolerance.
 */
gdouble 
ncm_function_cache_get_reltol (NcmFunctionCache *cache)
{
  NcmFunctionCachePrivate * const self = cache->priv;
  return self->reltol;
}

/**
 * ncm_function_cache_get_abstol:
 * @cache: a #NcmFunctionCache
 * 
 * Returns: the relative tolerance.
 */
gdouble 
ncm_function_cache_get_abstol (NcmFunctionCache *cache)
{
  NcmFunctionCachePrivate * const self = cache->priv;
  return self->abstol;
}

/**
 * ncm_function_cache_empty_cache:
 * @cache: a #NcmFunctionCache
 * 
 * Empties the content of @cache.
 * 
 */
void 
ncm_function_cache_empty_cache (NcmFunctionCache *cache)
{
  NcmFunctionCachePrivate * const self = cache->priv;
  self->clear = TRUE;
}

/**
 * ncm_function_cache_insert_vector:
 * @cache: a #NcmFunctionCache
 * @x: the argument $x$
 * @p: function value at $x$
 *
 * Insert a new point in the cache.
 *
 */
void
ncm_function_cache_insert_vector (NcmFunctionCache *cache, gdouble x, NcmVector *p)
{
  NcmFunctionCachePrivate * const self = cache->priv;
  gdouble *x_ptr = g_slice_new (gdouble);

  g_assert_cmpuint (self->n, ==, ncm_vector_len (p));

  *x_ptr = x;
  
  g_mutex_lock (&self->lock);
  cache_clean (cache);
  if (g_tree_lookup (self->tree, &x) != NULL)
  {
    g_mutex_unlock (&self->lock);
    return;
  }
  
  g_tree_insert (self->tree, x_ptr, ncm_vector_ref (p));
  
  g_mutex_unlock (&self->lock);
}

/**
 * ncm_function_cache_insert: (skip)
 * @cache: a #NcmFunctionCache
 * @x: the argument $x$
 * @...: function value at $x$
 *
 * Insert a new point in the cache. 
 *
 */
void
ncm_function_cache_insert (NcmFunctionCache *cache, gdouble x, ...)
{
  NcmFunctionCachePrivate * const self = cache->priv;
  gdouble *x_ptr;
  NcmVector *v;
  guint i;
  va_list ap;
  
  g_mutex_lock (&self->lock);
  
  cache_clean (cache);
  if (g_tree_lookup (self->tree, &x) != NULL)
  {
    g_mutex_unlock (&self->lock);
    return;
  }

  v     = ncm_vector_new (self->n);
  x_ptr = g_slice_new (gdouble);

  va_start (ap, x);

  for (i = 0; i < self->n; i++)
    ncm_vector_set (v, i, va_arg (ap, gdouble));

  *x_ptr = x;

  g_tree_insert (self->tree, x_ptr, v);

  g_mutex_unlock (&self->lock);
}

typedef struct _NcmFunctionCacheSearch
{
  gboolean found;
  gdouble near_x;
  gdouble x;
  gdouble diff;
  gint dir;
  gdouble reltol;
  NcmFunctionCacheSearchType type;
} NcmFunctionCacheSearch;

static gint gdouble_search_near (gconstpointer a, gconstpointer b);

/**
 * ncm_function_cache_get_near:
 * @cache: a #NcmFunctionCache
 * @x: the argument $x$
 * @x_found_ptr: Whether a point $x_c$ close to $x$ was found in the cache
 * @v: (out callee-allocates) (transfer full): the function at $x_c$ or NULL if no point was not found
 * @type: a #NcmFunctionCacheSearchType
 * 
 * Searches the @cache and returns the value of the function closest to $x$, $x_c$.
 * 
 * Returns: whether a point $x_c$ was found.
 */
gboolean
ncm_function_cache_get_near (NcmFunctionCache *cache, gdouble x, gdouble *x_found_ptr, NcmVector **v, NcmFunctionCacheSearchType type)
{
  NcmFunctionCachePrivate * const self = cache->priv;
  NcmFunctionCacheSearch search = {FALSE, x, 0.0, GSL_POSINF, 0, NCM_FUNCTION_CACHE_SEARCH_BOTH, self->reltol};
  NcmVector *res                = NULL;

  g_mutex_lock (&self->lock);
  search.type = type;

  if (cache_clean (cache))
  {
    g_mutex_unlock (&self->lock);
    return FALSE;
  }

  g_tree_search (self->tree, &gdouble_search_near, &search);
  if (search.found)
  {
    res          = g_tree_lookup (self->tree, &search.x);
    *x_found_ptr = search.x;
    *v           = ncm_vector_ref (res);
  }
  if (res == NULL)
  {
    g_mutex_unlock (&self->lock);
    return FALSE;
  }
  *v = res;

  g_mutex_unlock (&self->lock);
  return TRUE;
}

/**
 * ncm_function_cache_get: (skip)
 * @cache: a #NcmFunctionCache
 * @x_ptr: FIXME
 * @v: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_function_cache_get (NcmFunctionCache *cache, gdouble *x_ptr, NcmVector **v)
{
  NcmFunctionCachePrivate * const self = cache->priv;
  NcmVector *res;

  g_mutex_lock (&self->lock);

  if (cache_clean (cache))
  {
    g_mutex_unlock (&self->lock);
    return FALSE;
  }
  res = g_tree_lookup (self->tree, x_ptr);
  if (res == NULL)
  {
    g_mutex_unlock (&self->lock);
    return FALSE;
  }

  *v = res;
  g_mutex_unlock (&self->lock);
  return TRUE;
}

gboolean
cache_clean (NcmFunctionCache *cache)
{
  NcmFunctionCachePrivate * const self = cache->priv;
  if (self->clear)
  {
    g_tree_destroy (self->tree);
    self->tree = g_tree_new_full (&gdouble_compare, NULL, &gdouble_free, (GDestroyNotify)&ncm_vector_free);
    self->clear = FALSE;
    return TRUE;
  }
  return FALSE;
}

static gint
gdouble_compare (gconstpointer a, gconstpointer b, gpointer user_data)
{
  return gsl_fcmp ( *((gdouble *)a), *((gdouble *)b), 1e-15);
}

#define SX (*((gdouble *)a))

static gint
gdouble_search_near (gconstpointer a, gconstpointer b)
{
  NcmFunctionCacheSearch *search = (NcmFunctionCacheSearch *)b;
  gboolean found = FALSE;
  gdouble diff;
  
  search->dir = gsl_fcmp (search->near_x, SX, search->reltol);
  diff        = fabs (SX - search->near_x);

  switch (search->type)
  {
    case NCM_FUNCTION_CACHE_SEARCH_BOTH:
      found = TRUE;
      break;
    case NCM_FUNCTION_CACHE_SEARCH_GT:
      if (search->dir <= 0)
        found = TRUE;
      break;
    case NCM_FUNCTION_CACHE_SEARCH_LT:
      if (search->dir >= 0)
        found = TRUE;
      break;
  }

  if ((diff < search->diff) && found)
  {
    search->found = TRUE;
    search->x     = SX;
    search->diff  = diff;
  }

  return search->dir;
}

static void
gdouble_free (gpointer data)
{
  g_slice_free (gdouble, data);
}
