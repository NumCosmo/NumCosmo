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
 * @title: NcmMemoryPool
 * @short_description: Generic memory pool.
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
  g_mutex_init (&mp->update);
  g_cond_init (&mp->finish);
  mp->slices_in_use = 0;
  mp->slices        = g_ptr_array_new ();
  mp->alloc         = mp_alloc;
  mp->free          = mp_free;
  mp->userdata      = userdata;
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

  g_mutex_lock (&mp->update);
  while (mp->slices_in_use != 0)
    g_cond_wait (&mp->finish, &mp->update);
  
  for (i = 0; i < mp->slices->len; i++)
  {
    NcmMemoryPoolSlice *slice = g_ptr_array_index (mp->slices, i);
    if (free_slices && mp->free)
      mp->free (slice->p);
    g_slice_free (NcmMemoryPoolSlice, slice);
  }
  g_ptr_array_free (mp->slices, TRUE);
  mp->slices = NULL;

  g_mutex_unlock (&mp->update);

  g_mutex_clear (&mp->update);
  g_cond_clear (&mp->finish);

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
  g_mutex_lock (&mp->update);
  while (mp->slices->len < n)
  {
    NcmMemoryPoolSlice *slice = g_slice_new (NcmMemoryPoolSlice);
    slice->in_use = FALSE;
    slice->p      = mp->alloc (mp->userdata);
    slice->mp     = mp;
    g_ptr_array_add (mp->slices, slice);
  }
  g_mutex_unlock (&mp->update);
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
  g_mutex_lock (&mp->update);
  {
    NcmMemoryPoolSlice *slice = g_slice_new (NcmMemoryPoolSlice);
    slice->in_use = FALSE;
    slice->p      = p;
    slice->mp     = mp;
    g_ptr_array_add (mp->slices, slice);
  }
  g_mutex_unlock (&mp->update);
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

  g_mutex_lock (&mp->update);
  for (i = 0; i < mp->slices->len; i++)
  {
    NcmMemoryPoolSlice *lslice = (NcmMemoryPoolSlice *)g_ptr_array_index (mp->slices, i);
    if (!lslice->in_use)
    {
      slice         = lslice;
      slice->in_use = TRUE;
      break;
    }
  }
  if (slice == NULL)
  {
    slice         = g_slice_new (NcmMemoryPoolSlice);
    slice->in_use = TRUE;
    slice->p      = mp->alloc (mp->userdata);
    slice->mp     = mp;
    g_ptr_array_add (mp->slices, slice);
  }
  mp->slices_in_use++;
  g_mutex_unlock (&mp->update);

  return slice;
}

/**
 * ncm_memory_pool_return:
 * @p: slice to be returned to the pool
 *
 * Put the slice pointed by slice back to the pool.
 */
void
ncm_memory_pool_return (gpointer p)
{
  NcmMemoryPoolSlice *slice = (NcmMemoryPoolSlice *)p;

  g_mutex_lock (&slice->mp->update);
  slice->in_use = FALSE;
  slice->mp->slices_in_use--;

  if (slice->mp->slices_in_use == 0)
    g_cond_signal (&slice->mp->finish);

  g_mutex_unlock (&slice->mp->update);
}
