/***************************************************************************
 *            memory_pool.c
 *
 *  Wed June 15 18:53:30 2011
 *  Copyright  2011 Sandro Dias Pinto Vitenti
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
 * SECTION:memory_pool
 * @title: Memory Pool
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/memory_pool.h"

/**
 * ncm_memory_pool_new: (skip)
 * @mp_alloc: a #NcmMemoryPoolAlloc, function used to alloc memory.
 * @userdata: userdata pointer for @mp_alloc function.
 * @mp_free: function used to free memory alloced by mp_alloc.
 *
 * This function prepare a memory pool which allocate memory
 * using mp_alloc and save it for future use, the memory must
 * be returned to the pool using #ncm_memory_pool_return.
 * These functions are thread safe.
 *
 * Returns: the memory pool #NcmMemoryPool
 */
NcmMemoryPool *
ncm_memory_pool_new (NcmMemoryPoolAlloc mp_alloc, gpointer userdata, GDestroyNotify mp_free)
{
  NcmMemoryPool *mp = g_slice_new (NcmMemoryPool);
  _NCM_MUTEX_INIT (&mp->append);
  mp->slices = g_ptr_array_new ();
  mp->alloc = mp_alloc;
  mp->free = mp_free;
  mp->userdata = userdata;
  return mp;
}

/**
 * ncm_memory_pool_free:
 * @mp: a #NcmMemoryPool, memory pool to be freed
 * @free_slices: if true and the pool was built with a free function, free the slices
 *
 * This function free the memory pool and also
 * the slices if free_slices == TRUE and the
 * pool was built with a free function
 *
 */
void
ncm_memory_pool_free (NcmMemoryPool *mp, gboolean free_slices)
{
  guint i;
  _NCM_MUTEX_LOCK (&mp->append);
  for (i = 0; i < mp->slices->len; i++)
  {
    NcmMemoryPoolSlice *slice = g_ptr_array_index (mp->slices, i);
    _NCM_MUTEX_LOCK (&slice->lock);
    if (free_slices && mp->free)
      mp->free (slice->p);
    _NCM_MUTEX_UNLOCK (&slice->lock);
    _NCM_MUTEX_CLEAR (&slice->lock);
    g_slice_free (NcmMemoryPoolSlice, slice);
  }
  g_ptr_array_free (mp->slices, FALSE);
  mp->slices = NULL;
  _NCM_MUTEX_UNLOCK (&mp->append);
  _NCM_MUTEX_CLEAR (&mp->append);
  g_slice_free (NcmMemoryPool, mp);
}

/**
 * ncm_memory_pool_set_min_size:
 * @mp: a #NcmMemoryPool
 * @n: minimun number of slices contained in mp
 *
 * if n grater than number of slices then allocate new slices until
 * n == slices.
 *
 */
void
ncm_memory_pool_set_min_size (NcmMemoryPool *mp, gsize n)
{
  _NCM_MUTEX_LOCK (&mp->append);
  while (mp->slices->len < n)
  {
    NcmMemoryPoolSlice *slice = g_slice_new (NcmMemoryPoolSlice);
    _NCM_MUTEX_INIT (&slice->lock);
    slice->p = mp->alloc (mp->userdata);
    slice->mp = mp;
    g_ptr_array_add (mp->slices, slice);
  }
  _NCM_MUTEX_UNLOCK (&mp->append);
}

/**
 * ncm_memory_pool_add:
 * @mp: a #NcmMemoryPool.
 * @p: a pointer to an object compatible with the pool.
 * 
 * Adds an already allocated pointer @p to the pool. It will be freed with
 * #NcmMemoryPool->free during the pool destruction.
 * 
 */
void
ncm_memory_pool_add (NcmMemoryPool *mp, gpointer p)
{
  _NCM_MUTEX_LOCK (&mp->append);
  {
    NcmMemoryPoolSlice *slice = g_slice_new (NcmMemoryPoolSlice);
    _NCM_MUTEX_INIT (&slice->lock);
    slice->p = p;
    slice->mp = mp;
    g_ptr_array_add (mp->slices, slice);
  }
  _NCM_MUTEX_UNLOCK (&mp->append);
}

/**
 * ncm_memory_pool_get:
 * @mp: a #NcmMemoryPool
 *
 * Search in the pool for a non used slice
 * and return the first finded. If none
 * allocate a new one add to the pool and
 * return it.
 *
 * Returns: (transfer full): a pointer to an unused #NcmMemoryPoolSlice
 */
gpointer
ncm_memory_pool_get (NcmMemoryPool *mp)
{
  guint i;
  NcmMemoryPoolSlice *slice = NULL;

  _NCM_MUTEX_LOCK (&mp->append);
  for (i = 0; i < mp->slices->len; i++)
  {
    if (_NCM_MUTEX_TRYLOCK (&((NcmMemoryPoolSlice *)g_ptr_array_index (mp->slices, i))->lock))
    {
      slice = (NcmMemoryPoolSlice *)g_ptr_array_index (mp->slices, i);
      break;
    }
  }
  if (slice == NULL)
  {
    slice = g_slice_new (NcmMemoryPoolSlice);
    _NCM_MUTEX_INIT (&slice->lock);
    _NCM_MUTEX_LOCK (&slice->lock);
    slice->p = mp->alloc (mp->userdata);
    slice->mp = mp;
    g_ptr_array_add (mp->slices, slice);
  }
  _NCM_MUTEX_UNLOCK (&mp->append);

  return slice;
}

/**
 * ncm_memory_pool_return:
 * @p: slice to be returned to the pool
 *
 * Returns: the slice pointed by slice to the pool
 */
void
ncm_memory_pool_return (gpointer p)
{
  NcmMemoryPoolSlice *slice = (NcmMemoryPoolSlice *)p;
  _NCM_MUTEX_LOCK (&slice->mp->append);
  _NCM_MUTEX_UNLOCK (&slice->lock);
  _NCM_MUTEX_UNLOCK (&slice->mp->append);
}
