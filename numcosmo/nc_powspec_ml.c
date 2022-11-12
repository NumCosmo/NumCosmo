/***************************************************************************
 *            nc_powspec_ml.c
 *
 *  Thu February 18 12:32:13 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_powspec_ml.c
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * SECTION:nc_powspec_ml
 * @title: NcPowspecML
 * @short_description: Abstract class for linear matter power spectrum implementation.
 * @stability: Stable
 * @include: numcosmo/nc_powspec_ml.h
 *
 * This module comprises the set of functions to compute the linear matter power spectrum and
 * derived quantities.
 *
 * Following the description presented in #NcmPowspec, in this case we have that the field $\delta(\vec{x})$
 * represents the matter density fluctuations, i.e.,
 * $$\delta(\vec{x}) = \frac{\rho(\vec{x}) - \bar{\rho}}{\bar{\rho}},$$
 * where $\rho$ is the cold matter density field and $\bar{\rho}$ its mean.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_powspec_ml.h"

G_DEFINE_ABSTRACT_TYPE (NcPowspecML, nc_powspec_ml, NCM_TYPE_POWSPEC);

static void
nc_powspec_ml_init (NcPowspecML *nc_powspec_ml)
{
}

static void
_nc_powspec_ml_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_powspec_ml_parent_class)->finalize (object);
}

static void
nc_powspec_ml_class_init (NcPowspecMLClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  
  object_class->finalize = &_nc_powspec_ml_finalize;
}

/**
 * nc_powspec_ml_new_from_name:
 * @ps_ml_name: string which specifies the linear matter power spectrum object to be used
 *
 * This function returns a new #NcPowspecML whose type is defined by @ps_ml_name.
 *
 * Returns: A new #NcPowspecML.
 */
NcPowspecML *
nc_powspec_ml_new_from_name (const gchar *ps_ml_name)
{
  GObject *obj = ncm_serialize_global_from_string (ps_ml_name);
  
  if (!NC_IS_POWSPEC_ML (obj))
    g_error ("nc_powspec_ml_new_from_name: NcPowspecML %s do not descend from %s.", ps_ml_name, g_type_name (NC_TYPE_POWSPEC_ML));
  
  return NC_POWSPEC_ML (obj);
}

/**
 * nc_powspec_ml_ref:
 * @ps_ml: a #NcPowspecML
 *
 * Increases the reference count of @ps_ml by one atomically.
 *
 * Returns: (transfer full): @ps_ml.
 */
NcPowspecML *
nc_powspec_ml_ref (NcPowspecML *ps_ml)
{
  return g_object_ref (ps_ml);
}

/**
 * nc_powspec_ml_free:
 * @ps_ml: a #NcPowspecML
 *
 * Atomically decrements the reference count of @ps_ml by one.
 * If the reference count drops to 0, all memory allocated by @ps_ml is released.
 *
 */
void
nc_powspec_ml_free (NcPowspecML *ps_ml)
{
  g_object_unref (ps_ml);
}

/**
 * nc_powspec_ml_clear:
 * @ps_ml: a #NcPowspecML
 *
 *  If @ps_ml is different from NULL,
 *  atomically decrements the reference count of @ps_ml by one.
 *  If the reference count drops to 0,
 *  all memory allocated by @ps_ml is released and @ps_ml is set to NULL.
 *
 */
void
nc_powspec_ml_clear (NcPowspecML **ps_ml)
{
  g_clear_object (ps_ml);
}

