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
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_cfg.h>
#include <gsl/gsl_vector.h>

G_BEGIN_DECLS

/**
 * NcmFunctionCacheSearchType:
 * @NC_FUNCTION_CACHE_SEARCH_BOTH: FIXME
 * @NC_FUNCTION_CACHE_SEARCH_GT: FIXME
 * @NC_FUNCTION_CACHE_SEARCH_LT: FIXME
 *
 * FIXME
 */
typedef enum _NcmFunctionCacheSearchType
{
  NC_FUNCTION_CACHE_SEARCH_BOTH = 0,
  NC_FUNCTION_CACHE_SEARCH_GT,
  NC_FUNCTION_CACHE_SEARCH_LT,
} NcmFunctionCacheSearchType;

typedef struct _NcmFunctionCache NcmFunctionCache;

struct _NcmFunctionCache
{
  /*< private >*/
  GTree *tree;
  _NCM_MUTEX_TYPE lock;
  gboolean clear;
  guint n;
  gdouble abstol;
  gdouble reltol;
};

NcmFunctionCache *ncm_function_cache_new (guint n, gdouble abstol, gdouble reltol);
void ncm_function_cache_free (NcmFunctionCache *cache);
void ncm_function_cache_clear (NcmFunctionCache **cache);
void ncm_function_cache_insert (NcmFunctionCache *cache, gdouble x, ...);
void ncm_function_cache_insert_vector (NcmFunctionCache *cache, gdouble x, gsl_vector *p);
gboolean ncm_function_cache_get (NcmFunctionCache *cache, gdouble *x_ptr, gsl_vector **v);
gboolean ncm_function_cache_get_near (NcmFunctionCache *cache, gdouble x, gdouble *x_found_ptr, gsl_vector **v, NcmFunctionCacheSearchType type);

#define NC_FUNCTION_CACHE(p) ((NcmFunctionCache *)(p))

G_END_DECLS

#endif /* _NC_FUNCTION_CACHE_H */
