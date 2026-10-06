/***************************************************************************
 *            ncm_prefetch_private.h
 *
 *  Wed September 23 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * ncm_prefetch_private.h
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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

#ifndef _NCM_PREFETCH_PRIVATE_H_
#define _NCM_PREFETCH_PRIVATE_H_

#include <glib.h>

G_BEGIN_DECLS

/*
 * Software prefetch hints for loops whose per-iteration data is scattered
 * over the heap (one small allocation per galaxy, say), so the hardware
 * prefetcher cannot follow it. A prefetch never faults, so a stale or NULL
 * address is harmless; it only costs the instruction. Without GCC/Clang
 * builtins, or with NCM_DISABLE_PREFETCH defined (for A/B timing), both
 * helpers compile to nothing.
 */

#if (defined (__GNUC__) || defined (__clang__)) && !defined (NCM_DISABLE_PREFETCH)
#define NCM_PREFETCH_ENABLED 1
#endif

#define NCM_PREFETCH_LINE 64

static inline void
ncm_prefetch (const void *p)
{
#ifdef NCM_PREFETCH_ENABLED
  __builtin_prefetch (p, 0, 3);
#endif
}

/* Every cache line overlapping [p, p + bytes). */
static inline void
ncm_prefetch_span (const void *p, const gsize bytes)
{
#ifdef NCM_PREFETCH_ENABLED

  if ((p != NULL) && (bytes > 0))
  {
    guintptr a         = ((guintptr) p) & ~((guintptr) NCM_PREFETCH_LINE - 1);
    const guintptr end = ((guintptr) p) + bytes;

    for ( ; a < end; a += NCM_PREFETCH_LINE)
      __builtin_prefetch ((const void *) a, 0, 3);
  }

#endif
}

G_END_DECLS

#endif /* _NCM_PREFETCH_PRIVATE_H_ */

