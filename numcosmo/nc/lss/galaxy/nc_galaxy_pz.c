/***************************************************************************
 *            nc_galaxy_pz.c
 *
 *  Tue Nov 17 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_pz.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcGalaxyPz:
 *
 * Abstract class for the redshift probability density $p(z)$ of one galaxy.
 *
 * An object is data, not a model: it has no parameters and does not change after
 * construction, so one object may be shared by several catalogs and read concurrently.
 * The density is not normalized: nc_galaxy_pz_get_norm() gives its integral over the
 * range of nc_galaxy_pz_get_range(), and nc_galaxy_pz_gen() draws from the normalized
 * density. The evaluation functions check nothing: the caller keeps $z$ inside the
 * range.
 *
 * A subclass implements the virtual functions of nc_galaxy_pz_get_range(),
 * nc_galaxy_pz_eval_p(), nc_galaxy_pz_eval_ln_p(), nc_galaxy_pz_get_norm() and
 * nc_galaxy_pz_gen(). The one of nc_galaxy_pz_cheb_moments() is optional: without it the
 * call aborts.
 *
 * The object keeps the moments of its last %NC_GALAXY_PZ_CHEB_CACHE_LEN calls to
 * nc_galaxy_pz_cheb_moments(), keyed on the interval and the order, and returns them
 * without computing them again. The cache is guarded by a lock, so concurrent callers
 * compute each entry once; it is not serialized.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_pz.h"

#include <string.h>

/* One computed set of moments: m_0..m_n on [a, b] */
typedef struct _NcGalaxyPzChebEntry
{
  gdouble a;
  gdouble b;
  guint n;
  GArray *moments;
} NcGalaxyPzChebEntry;

typedef struct _NcGalaxyPzPrivate
{
  GMutex cheb_lock;
  GArray *cheb_cache; /* element type: NcGalaxyPzChebEntry, oldest first */
} NcGalaxyPzPrivate;

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcGalaxyPz, nc_galaxy_pz, G_TYPE_OBJECT)

static void
_nc_galaxy_pz_cheb_entry_clear (gpointer entry)
{
  g_clear_pointer (&((NcGalaxyPzChebEntry *) entry)->moments, g_array_unref);
}

static void
nc_galaxy_pz_init (NcGalaxyPz *gpz)
{
  NcGalaxyPzPrivate * const self = nc_galaxy_pz_get_instance_private (gpz);

  g_mutex_init (&self->cheb_lock);
  self->cheb_cache = g_array_new (FALSE, FALSE, sizeof (NcGalaxyPzChebEntry));
  g_array_set_clear_func (self->cheb_cache, &_nc_galaxy_pz_cheb_entry_clear);
}

static void
_nc_galaxy_pz_finalize (GObject *object)
{
  NcGalaxyPzPrivate * const self = nc_galaxy_pz_get_instance_private (NC_GALAXY_PZ (object));

  g_clear_pointer (&self->cheb_cache, g_array_unref);
  g_mutex_clear (&self->cheb_lock);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_pz_parent_class)->finalize (object);
}

static void
_nc_galaxy_pz_cheb_moments_not_implemented (NcGalaxyPz *gpz, const gdouble a, const gdouble b, const guint n, GArray *moments)
{
  g_error ("nc_galaxy_pz_cheb_moments: %s does not implement the Chebyshev moments.", G_OBJECT_TYPE_NAME (gpz));
}

static void
nc_galaxy_pz_class_init (NcGalaxyPzClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->finalize = &_nc_galaxy_pz_finalize;
  klass->cheb_moments    = &_nc_galaxy_pz_cheb_moments_not_implemented;
}

/**
 * nc_galaxy_pz_ref:
 * @gpz: a #NcGalaxyPz
 *
 * Increases the reference count of @gpz by one.
 *
 * Returns: (transfer full): @gpz.
 */
NcGalaxyPz *
nc_galaxy_pz_ref (NcGalaxyPz *gpz)
{
  return g_object_ref (gpz);
}

/**
 * nc_galaxy_pz_free:
 * @gpz: a #NcGalaxyPz
 *
 * Decreases the reference count of @gpz by one.
 */
void
nc_galaxy_pz_free (NcGalaxyPz *gpz)
{
  g_object_unref (gpz);
}

/**
 * nc_galaxy_pz_clear:
 * @gpz: a #NcGalaxyPz
 *
 * Decreases the reference count of *@gpz by one and sets *@gpz to %NULL; does nothing
 * when *@gpz is %NULL.
 */
void
nc_galaxy_pz_clear (NcGalaxyPz **gpz)
{
  g_clear_object (gpz);
}

/**
 * nc_galaxy_pz_get_range:
 * @gpz: a #NcGalaxyPz
 * @z_min: (out): the lower end of the range
 * @z_max: (out): the upper end of the range
 *
 * Gets the redshift range where the density is defined.
 */
void
nc_galaxy_pz_get_range (NcGalaxyPz *gpz, gdouble *z_min, gdouble *z_max)
{
  NC_GALAXY_PZ_GET_CLASS (gpz)->get_range (gpz, z_min, z_max);
}

/**
 * nc_galaxy_pz_eval_p:
 * @gpz: a #NcGalaxyPz
 * @z: the redshift
 *
 * Evaluates the density at @z, not normalized. The redshift is not checked against the
 * range.
 *
 * Returns: $p(z)$.
 */
gdouble
nc_galaxy_pz_eval_p (NcGalaxyPz *gpz, const gdouble z)
{
  return NC_GALAXY_PZ_GET_CLASS (gpz)->eval_p (gpz, z);
}

/**
 * nc_galaxy_pz_eval_ln_p:
 * @gpz: a #NcGalaxyPz
 * @z: the redshift
 *
 * Evaluates the logarithm of the density at @z, not normalized. The redshift is not
 * checked against the range.
 *
 * Returns: $\ln p(z)$.
 */
gdouble
nc_galaxy_pz_eval_ln_p (NcGalaxyPz *gpz, const gdouble z)
{
  return NC_GALAXY_PZ_GET_CLASS (gpz)->eval_ln_p (gpz, z);
}

/**
 * nc_galaxy_pz_get_norm:
 * @gpz: a #NcGalaxyPz
 *
 * Gets the integral of the density over its range.
 *
 * Returns: the normalization.
 */
gdouble
nc_galaxy_pz_get_norm (NcGalaxyPz *gpz)
{
  return NC_GALAXY_PZ_GET_CLASS (gpz)->get_norm (gpz);
}

/**
 * nc_galaxy_pz_gen:
 * @gpz: a #NcGalaxyPz
 * @rng: a #NcmRNG
 *
 * Draws a redshift from the normalized density with @rng.
 *
 * Returns: the redshift drawn.
 */
gdouble
nc_galaxy_pz_gen (NcGalaxyPz *gpz, NcmRNG *rng)
{
  return NC_GALAXY_PZ_GET_CLASS (gpz)->gen (gpz, rng);
}

/**
 * nc_galaxy_pz_cheb_moments:
 * @gpz: a #NcGalaxyPz
 * @a: the lower end of the interval
 * @b: the upper end of the interval
 * @n: the highest order
 * @moments: (out callee-allocates) (transfer full) (element-type gdouble): the moments
 *
 * Computes the Chebyshev moments of the density on $[a, b]$,
 * $$
 * m_k = \int_a^b p(z)\, T_k\!\left(\frac{2z - a - b}{b - a}\right) \mathrm{d}z, \qquad k = 0, \dots, n,
 * $$
 * where $T_k$ is the Chebyshev polynomial of degree $k$, so $m_0$ is the mass of $p$ on
 * $[a, b]$. A C caller may pass in @moments an array to reuse, or a pointer to %NULL for a
 * new one; it receives $n + 1$ values. The moments of the same $a$, $b$ and $n$ as one of
 * the last %NC_GALAXY_PZ_CHEB_CACHE_LEN calls are copied from the cache. An interval that
 * is empty or not inside the range of nc_galaxy_pz_get_range() aborts, and so does a
 * subclass that does not implement the moments.
 */
void
nc_galaxy_pz_cheb_moments (NcGalaxyPz *gpz, const gdouble a, const gdouble b, const guint n, GArray **moments)
{
  NcGalaxyPzPrivate * const self = nc_galaxy_pz_get_instance_private (gpz);
  NcGalaxyPzChebEntry *found     = NULL;
  gdouble z_min, z_max;
  guint i;

  nc_galaxy_pz_get_range (gpz, &z_min, &z_max);

  if (!((z_min <= a) && (a < b) && (b <= z_max)))
    g_error ("nc_galaxy_pz_cheb_moments: the interval [%g, %g] is empty or not inside the range [%g, %g].",
             a, b, z_min, z_max);

  if (*moments == NULL)
    *moments = g_array_new (FALSE, FALSE, sizeof (gdouble));

  g_array_set_size (*moments, n + 1);

  /* The lock is held while an entry is computed, so a concurrent caller of the same entry
   * waits for it instead of computing it again */
  g_mutex_lock (&self->cheb_lock);

  for (i = 0; i < self->cheb_cache->len; i++)
  {
    NcGalaxyPzChebEntry *entry = &g_array_index (self->cheb_cache, NcGalaxyPzChebEntry, i);

    if ((entry->a == a) && (entry->b == b) && (entry->n == n))
    {
      found = entry;
      break;
    }
  }

  if (found == NULL)
  {
    NcGalaxyPzChebEntry entry = {a, b, n, g_array_sized_new (FALSE, FALSE, sizeof (gdouble), n + 1)};

    g_array_set_size (entry.moments, n + 1);
    NC_GALAXY_PZ_GET_CLASS (gpz)->cheb_moments (gpz, a, b, n, entry.moments);

    if (self->cheb_cache->len == NC_GALAXY_PZ_CHEB_CACHE_LEN)
      g_array_remove_index (self->cheb_cache, 0);

    g_array_append_val (self->cheb_cache, entry);
    found = &g_array_index (self->cheb_cache, NcGalaxyPzChebEntry, self->cheb_cache->len - 1);
  }

  memcpy ((*moments)->data, found->moments->data, sizeof (gdouble) * (n + 1));
  g_mutex_unlock (&self->cheb_lock);
}

