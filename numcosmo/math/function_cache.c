/***************************************************************************
 *            function_cache.c
 *
 *  Wed Aug 13 21:18:48 2008
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
 * SECTION:function_cache
 * @title: Function Cache
 * @short_description: A generic cache for functions values
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <glib.h>


/***************************************************************************
 *
 *
 ****************************************************************************/

static gint gdouble_compare (gconstpointer a, gconstpointer b, gpointer user_data);
static void gdouble_free (gpointer data);
static gboolean cache_clean (NcFunctionCache *cache);

/**
 * nc_function_cache_new: (skip)
 * @n: FIXME
 * @abstol: FIXME
 * @reltol: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcFunctionCache *
nc_function_cache_new (guint n, gdouble abstol, gdouble reltol)
{
  NcFunctionCache *cache = g_slice_new (NcFunctionCache);
  cache->tree = g_tree_new_full (&gdouble_compare, NULL, &gdouble_free, (GDestroyNotify)&gsl_vector_free);
  g_static_mutex_init(&cache->lock);
  cache->clear = FALSE;
  cache->n = n;
  cache->abstol = abstol;
  cache->reltol = reltol;

  return cache;
}

/**
 * nc_function_cache_free:
 * @cache: a #NcFunctionCache
 *
 * FIXME
 *
*/
void
nc_function_cache_free (NcFunctionCache *cache)
{
  g_static_mutex_free (&cache->lock);
  g_tree_destroy (cache->tree);
  g_slice_free (NcFunctionCache, cache);
  return;
}

/**
 * nc_function_cache_insert_vector: (skip)
 * @cache: a #NcFunctionCache
 * @x: FIXME
 * @p: FIXME
 *
 * FIXME
 *
*/
void
nc_function_cache_insert_vector (NcFunctionCache *cache, gdouble x, gsl_vector *p)
{
  gdouble *x_ptr = g_slice_new (gdouble);
  g_assert (cache->n == p->size);

  *x_ptr = x;
  g_static_mutex_lock (&cache->lock);
  cache_clean (cache);
  g_tree_insert (cache->tree, x_ptr, p);
  g_static_mutex_unlock (&cache->lock);
}

void
nc_function_cache_insert (NcFunctionCache *cache, gdouble x, ...)
{
  gdouble *x_ptr;
  gsl_vector *v;
  gint i;
  va_list ap;

  g_static_mutex_lock (&cache->lock);
  cache_clean (cache);

  if (g_tree_lookup (cache->tree, &x) != NULL)
  {
    g_static_mutex_unlock (&cache->lock);
    return;
  }

  v = gsl_vector_alloc (cache->n);
  x_ptr = g_slice_new (gdouble);

  va_start(ap, x);

  for (i = 0; i < cache->n; i++)
    gsl_vector_set (v, i, va_arg(ap, gdouble));
  *x_ptr = x;

  g_tree_insert (cache->tree, x_ptr, v);
  g_static_mutex_unlock (&cache->lock);
}

typedef struct _NcParamsCacheSearch
{
  gboolean found;
  gdouble near_x;
  gdouble x;
  gdouble diff;
  gint dir;
  NcFunctionCacheSearchType type;
} NcParamsCacheSearch;

static gint gdouble_search_near (gconstpointer a, gconstpointer b);

/**
 * nc_function_cache_get_near: (skip)
 * @cache: a #NcFunctionCache
 * @x: FIXME
 * @x_found_ptr: FIXME
 * @v: FIXME
 * @type: a #NcFunctionCacheSearchType
*/
gboolean
nc_function_cache_get_near (NcFunctionCache *cache, gdouble x, gdouble *x_found_ptr, gsl_vector **v, NcFunctionCacheSearchType type)
{
  NcParamsCacheSearch search = {FALSE, x, 0.0, GSL_POSINF, 0, NC_FUNCTION_CACHE_SEARCH_BOTH};
//  gdouble *res = NULL;
  gsl_vector *res = NULL;
  g_static_mutex_lock (&cache->lock);

  search.type = type;

//  printf ("--== %g %g %g %d\n", *x_ptr, search.near_x, search.x, search.dir);

  if(cache_clean (cache))
  {
    g_static_mutex_unlock (&cache->lock);
    return FALSE;
  }

  g_tree_search (cache->tree, &gdouble_search_near, &search);
  if (search.found)
  {
    res = g_tree_lookup (cache->tree, &search.x);
    *x_found_ptr = search.x;
    *v = res;
    //printf ("[%p]Found!! %g ... %g [%p] \n", g_thread_self (), search.x, gsl_vector_get (*v, 0), *v);
  }
  if (res == NULL)
  {
    g_static_mutex_unlock (&cache->lock);
    return FALSE;
  }
  //*f_ptr = *res;
  *v = res;

  g_static_mutex_unlock (&cache->lock);
  //printf ("[%p]Found2! %g ... %g [%p] \n", g_thread_self (), search.x, gsl_vector_get (*v, 0), *v);
  return TRUE;
}

/**
 * nc_function_cache_get: (skip)
 * @cache: a #NcFunctionCache
 * @x_ptr: FIXME
 * @v: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gboolean
nc_function_cache_get (NcFunctionCache *cache, gdouble *x_ptr, gsl_vector **v)
{
  gsl_vector *res;
  //gdouble *res;
  g_static_mutex_lock (&cache->lock);
  if(cache_clean (cache))
  {
    g_static_mutex_unlock (&cache->lock);
    return FALSE;
  }
  res = g_tree_lookup (cache->tree, x_ptr);
  if (res == NULL)
  {
    g_static_mutex_unlock (&cache->lock);
    return FALSE;
  }
  //*val_ptr = *res;
  *v = res;
  g_static_mutex_unlock (&cache->lock);
  return TRUE;
}

gboolean
cache_clean (NcFunctionCache *cache)
{
  if (cache->clear)
  {
    g_tree_destroy (cache->tree);
    cache->tree = g_tree_new_full (&gdouble_compare, NULL, &gdouble_free, (GDestroyNotify)&gsl_vector_free);
    cache->clear = FALSE;
    return TRUE;
  }
  return FALSE;
}

static gint
gdouble_compare (gconstpointer a, gconstpointer b, gpointer user_data)
{
//  printf ("BLOB %g %g\n", *((gdouble *)a), *((gdouble *)b));
  return gsl_fcmp ( *((gdouble *)a), *((gdouble *)b), 1e-15);
}

#define SX (*((gdouble *)a))

static gint
gdouble_search_near (gconstpointer a, gconstpointer b)
{
  NcParamsCacheSearch *search = (NcParamsCacheSearch *)b;
  gdouble diff;
  gboolean found = FALSE;
  search->dir = gsl_fcmp (search->near_x, SX, NC_ZERO_LIMIT);
  diff = fabs (SX - search->near_x);

  switch (search->type)
  {
    case NC_FUNCTION_CACHE_SEARCH_BOTH:
      found = TRUE;
      break;
    case NC_FUNCTION_CACHE_SEARCH_GT:
      if (search->dir <= 0)
        found = TRUE;
      break;
    case NC_FUNCTION_CACHE_SEARCH_LT:
      if (search->dir >= 0)
        found = TRUE;
      break;
  }

  if ((diff < search->diff) && found)
  {
    search->found = TRUE;
    search->x = SX;
    search->diff = diff;
  }

  return search->dir;
}

static void
gdouble_free (gpointer data)
{
  g_slice_free (gdouble, data);
}
