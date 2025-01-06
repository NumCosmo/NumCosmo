/***************************************************************************
 *            nc_powspec_mnl.c
 *
 *  Thu February 18 12:32:13 2016
 *  Copyright  2016  Cyrille Doux
 *  <cdoux@apc.in2p3.fr>
 ****************************************************************************/
/*
 * nc_powspec_mnl.c
 * Copyright (C) 2016 Cyrille Doux <cdoux@apc.in2p3.fr>
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
 * NcPowspecMNL:
 *
 * Abstrac class for non-linear matter power spectrum implementation.
 *
 * This module comprises the set of functions to compute a power spectrum and derived
 * quantities.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_powspec_mnl.h"

G_DEFINE_ABSTRACT_TYPE (NcPowspecMNL, nc_powspec_mnl, NCM_TYPE_POWSPEC)

static void
nc_powspec_mnl_init (NcPowspecMNL *nc_powspec_mnl)
{
}

static void
nc_powspec_mnl_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_powspec_mnl_parent_class)->finalize (object);
}

static void
nc_powspec_mnl_class_init (NcPowspecMNLClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->finalize = nc_powspec_mnl_finalize;
}

/**
 * nc_powspec_mnl_ref:
 * @ps_mnl: a #NcmMSetCatalog
 *
 * Increases the reference count of @ps_mnl by one atomically.
 *
 * Returns: (transfer full): @ps_mnl.
 */
NcPowspecMNL *
nc_powspec_mnl_ref (NcPowspecMNL *ps_mnl)
{
  return g_object_ref (ps_mnl);
}

/**
 * nc_powspec_mnl_free:
 * @ps_mnl: a #NcmMSetCatalog
 *
 * Atomically decrements the reference count of @ps_mnl by one.
 * If the reference count drops to 0, all memory allocated by @ps_mnl is released.
 *
 */
void
nc_powspec_mnl_free (NcPowspecMNL *ps_mnl)
{
  g_object_unref (ps_mnl);
}

/**
 * nc_powspec_mnl_clear:
 * @ps_mnl: a #NcmMSetCatalog
 *
 * If @ps_mnl is different from NULL,
 * atomically decrements the reference count of @ps_mnl by one.
 * If the reference count drops to 0,
 * all memory allocated by @ps_mnl is released and @ps_mnl is set to NULL.
 *
 */
void
nc_powspec_mnl_clear (NcPowspecMNL **ps_mnl)
{
  g_clear_object (ps_mnl);
}

