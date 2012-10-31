/***************************************************************************
 *            function_cache.h
 *
 *  Wed Aug 13 21:18:28 2008
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

#ifndef _NC_FUNCTION_CACHE_H
#define _NC_FUNCTION_CACHE_H

#include <glib.h>
#include <glib-object.h>

#include <numcosmo/nc_macros.h>
#include <gsl/gsl_vector.h>

G_BEGIN_DECLS

/**
 * NcFunctionCacheSearchType:
 * @NC_FUNCTION_CACHE_SEARCH_BOTH: FIXME
 * @NC_FUNCTION_CACHE_SEARCH_GT: FIXME
 * @NC_FUNCTION_CACHE_SEARCH_LT: FIXME
 *
 * FIXME
 */
typedef enum _NcFunctionCacheSearchType
{
  NC_FUNCTION_CACHE_SEARCH_BOTH = 0,
  NC_FUNCTION_CACHE_SEARCH_GT,
  NC_FUNCTION_CACHE_SEARCH_LT,
} NcFunctionCacheSearchType;

typedef struct _NcFunctionCache NcFunctionCache;

struct _NcFunctionCache
{
  /*< private >*/
  GTree *tree;
  _NCM_MUTEX_TYPE lock;
  gboolean clear;
  guint n;
  gdouble abstol;
  gdouble reltol;
};

NcFunctionCache *nc_function_cache_new (guint n, gdouble abstol, gdouble reltol);
void nc_function_cache_free (NcFunctionCache *cache);
void nc_function_cache_insert (NcFunctionCache *cache, gdouble x, ...);
void nc_function_cache_insert_vector (NcFunctionCache *cache, gdouble x, gsl_vector *p);
gboolean nc_function_cache_get (NcFunctionCache *cache, gdouble *x_ptr, gsl_vector **v);
gboolean nc_function_cache_get_near (NcFunctionCache *cache, gdouble x, gdouble *x_found_ptr, gsl_vector **v, NcFunctionCacheSearchType type);

#define NC_FUNCTION_CACHE(p) ((NcFunctionCache *)(p))

G_END_DECLS

#endif /* _NC_FUNCTION_CACHE_H */
