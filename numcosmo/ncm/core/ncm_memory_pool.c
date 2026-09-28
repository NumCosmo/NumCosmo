/***************************************************************************
 *            ncm_memory_pool.c
 *
 *  Wed June 15 18:53:30 2011
 *  Copyright  2011 Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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
 * NcmMemoryPool:
 *
 * Thread-safe pool of reusable objects.
 *
 * Objects are allocated with a user-provided function, kept for reuse, and
 * released with the corresponding free function.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/core/ncm_memory_pool.h"

/**
 * ncm_memory_pool_new: (skip)
 * @mp_alloc: allocation function
 * @userdata: user data
 * @mp_free: (nullable): free function
 *
 * Creates a pool that allocates objects with @mp_alloc (called with
 * @userdata), reuses returned objects and frees them with @mp_free. Objects
 * must be returned with ncm_memory_pool_return().
 *
 * Returns: a new #NcmMemoryPool.
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
 * ncm_memory_pool_empty:
 * @mp: a #NcmMemoryPool
 * @free_slices: whether to free the pooled objects
 *
 * Removes every slice from @mp. The memory each slice points to is released
 * with the pool's free function only when @free_slices is %TRUE and such a
 * function was given to ncm_memory_pool_new(). Blocks until all slices
 * currently checked out have been returned.
 *
 * Returns: the number of slices removed
 */
guint
ncm_memory_pool_empty (NcmMemoryPool *mp, gboolean free_slices)
{
  guint i, n;

  g_mutex_lock (&mp->update);

  while (mp->slices_in_use != 0)
    g_cond_wait (&mp->finish, &mp->update);

  n = mp->slices->len;

  for (i = 0; i < mp->slices->len; i++)
  {
    NcmMemoryPoolSlice *slice = g_ptr_array_index (mp->slices, i);

    if (free_slices && mp->free)
      mp->free (slice->p);

    g_slice_free (NcmMemoryPoolSlice, slice);
  }

  g_ptr_array_free (mp->slices, TRUE);
  mp->slices = g_ptr_array_new ();

  g_mutex_unlock (&mp->update);

  return n;
}

/**
 * ncm_memory_pool_free:
 * @mp: a #NcmMemoryPool
 * @free_slices: whether to free the pooled objects
 *
 * Frees @mp. The memory each slice points to is released with the pool's free
 * function only when @free_slices is %TRUE and such a function was given to
 * ncm_memory_pool_new(). Blocks until all slices currently checked out have
 * been returned.
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
 * @n: minimum number of slices
 *
 * Allocates slices until @mp holds at least @n.
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
 * @mp: a #NcmMemoryPool
 * @p: an object
 *
 * Adds @p, allocated as by the pool's allocation function, to @mp. The pool
 * owns @p under the same rules as the objects it allocates.
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
 * Checks out the first unused slice, allocating a new one if all are in use.
 *
 * Returns: (transfer full): a #NcmMemoryPoolSlice, to be given back with ncm_memory_pool_return().
 */
gpointer
ncm_memory_pool_get (NcmMemoryPool *mp)
{
  guint i;
  NcmMemoryPoolSlice *slice = NULL;

  g_mutex_lock (&mp->update);

  for (i = 0; i < mp->slices->len; i++)
  {
    NcmMemoryPoolSlice *lslice = (NcmMemoryPoolSlice *) g_ptr_array_index (mp->slices, i);

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
 * @p: a #NcmMemoryPoolSlice
 *
 * Returns @p, obtained from ncm_memory_pool_get(), to its pool.
 */
void
ncm_memory_pool_return (gpointer p)
{
  NcmMemoryPoolSlice *slice = (NcmMemoryPoolSlice *) p;

  g_mutex_lock (&slice->mp->update);
  slice->in_use = FALSE;
  slice->mp->slices_in_use--;

  if (slice->mp->slices_in_use == 0)
    g_cond_signal (&slice->mp->finish);

  g_mutex_unlock (&slice->mp->update);
}

