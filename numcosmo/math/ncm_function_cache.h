/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_function_cache.h
 *
 *  Wed Aug 13 21:18:28 2008
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

#ifndef _NCM_FUNCTION_CACHE_H_
#define _NCM_FUNCTION_CACHE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>

G_BEGIN_DECLS

#define NCM_TYPE_FUNCTION_CACHE             (ncm_function_cache_get_type ())
#define NCM_FUNCTION_CACHE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_FUNCTION_CACHE, NcmFunctionCache))
#define NCM_FUNCTION_CACHE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_FUNCTION_CACHE, NcmFunctionCacheClass))
#define NCM_IS_FUNCTION_CACHE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_FUNCTION_CACHE))
#define NCM_IS_FUNCTION_CACHE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_FUNCTION_CACHE))
#define NCM_FUNCTION_CACHE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_FUNCTION_CACHE, NcmFunctionCacheClass))

typedef struct _NcmFunctionCacheClass NcmFunctionCacheClass;
typedef struct _NcmFunctionCache NcmFunctionCache;
typedef struct _NcmFunctionCachePrivate NcmFunctionCachePrivate;

struct _NcmFunctionCacheClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcmFunctionCache
{
  /*< private >*/
  GObject parent_instance;
  NcmFunctionCachePrivate *priv;
};

/**
 * NcmFunctionCacheSearchType:
 * @NCM_FUNCTION_CACHE_SEARCH_BOTH: Searches both directions.
 * @NCM_FUNCTION_CACHE_SEARCH_GT: Searches cache upwards.
 * @NCM_FUNCTION_CACHE_SEARCH_LT: Searches cache downwards.
 *
 * Cache search direction.
 */
typedef enum _NcmFunctionCacheSearchType
{
  NCM_FUNCTION_CACHE_SEARCH_BOTH = 0,
  NCM_FUNCTION_CACHE_SEARCH_GT,
  NCM_FUNCTION_CACHE_SEARCH_LT,
} NcmFunctionCacheSearchType;

GType ncm_function_cache_get_type (void) G_GNUC_CONST;

NcmFunctionCache *ncm_function_cache_new (guint n, gdouble abstol, gdouble reltol);
NcmFunctionCache *ncm_function_cache_ref (NcmFunctionCache *cache);
void ncm_function_cache_free (NcmFunctionCache *cache);
void ncm_function_cache_clear (NcmFunctionCache **cache);

gdouble ncm_function_cache_get_reltol (NcmFunctionCache *cache);
gdouble ncm_function_cache_get_abstol (NcmFunctionCache *cache);

void ncm_function_cache_empty_cache (NcmFunctionCache *cache);

void ncm_function_cache_insert (NcmFunctionCache *cache, gdouble x, ...);
void ncm_function_cache_insert_vector (NcmFunctionCache *cache, gdouble x, NcmVector *p);
gboolean ncm_function_cache_get (NcmFunctionCache *cache, gdouble *x_ptr, NcmVector **v);
gboolean ncm_function_cache_get_near (NcmFunctionCache *cache, gdouble x, gdouble *x_found_ptr, NcmVector **v, NcmFunctionCacheSearchType type);

G_END_DECLS

#endif /* _NCM_FUNCTION_CACHE_H_ */
