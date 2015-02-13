/***************************************************************************
 *            nc_density_profile.c
 *
 *  Sat June 07 19:46:31 2014
 *  Copyright  2014  
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_density_profile.c
 * Copyright (C) 2014 Mariana Penna Lima <pennalima@gmail.com>
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
 * SECTION:nc_density_profile
 * @title: NcDensityProfile
 * @short_description: Abstract class for density profile functions. 
 * 
 * This module comprises the set of functions to compute the matter density profile in both real 
 * and Fourier spaces. 
 * 
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_density_profile.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"

G_DEFINE_TYPE (NcDensityProfile, nc_density_profile, G_TYPE_OBJECT);

/**
 * nc_density_profile_new_from_name:
 * @density_profile_name: "#NcDensityProfileNFW".
 * 
 * This function returns a new #NcDensityProfile whose type is defined by @density_profile_name string.
 * 
 * Returns: A new #NcDensityProfile.
 */
NcDensityProfile *
nc_density_profile_new_from_name (gchar *density_profile_name)
{
  GObject *obj = ncm_serialize_global_from_string (density_profile_name);
  GType density_profile_type = G_OBJECT_TYPE (obj);
  if (!g_type_is_a (density_profile_type, NC_TYPE_DENSITY_PROFILE))
    g_error ("nc_density_profile_new_from_name: NcDensityProfile %s do not descend from %s.", density_profile_name, g_type_name (NC_TYPE_DENSITY_PROFILE));
  return NC_DENSITY_PROFILE (obj);
}

/**
 * nc_density_profile_eval_fourier:
 * @dp: a #NcDensityProfile
 * @model: a #NcHICosmo
 * @k: mode 
 * @M: mass
 * @z: redshift
 *  
 * This function computes the density profile in the Fourier space. 
 *
 * Returns: The value of the density profile in the Fourier space.
 */
gdouble 
nc_density_profile_eval_fourier (NcDensityProfile *dp, NcHICosmo *model, const gdouble k, const gdouble M, const gdouble z) 
{ 
  return NC_DENSITY_PROFILE_GET_CLASS (dp)->eval_fourier (dp, model, k, M, z);
}

/**
 * nc_density_profile_free: 
 * @dp: a #NcDensityProfile.
 *
 * Atomically decrements the reference count of @dp by one. If the reference count drops to 0, 
 * all memory allocated by @dp is released.
 *
 */
void 
nc_density_profile_free (NcDensityProfile *dp)
{
  g_object_unref (dp);
}

/**
 * nc_density_profile_clear: 
 * @dp: a #NcDensityProfile.
 *
 * Atomically decrements the reference count of @dp by one. If the reference count drops to 0, 
 * all memory allocated by @dp is released. Set the pointer to NULL;
 *
 */
void 
nc_density_profile_clear (NcDensityProfile **dp)
{
  g_clear_object (dp);
}

static void
nc_density_profile_init (NcDensityProfile *nc_density_profile)
{
  NCM_UNUSED (nc_density_profile);  
}

static void
nc_density_profile_finalize (GObject *object)
{
  /* TODO: Add deinitalization code here */

  G_OBJECT_CLASS (nc_density_profile_parent_class)->finalize (object);
}

static void
nc_density_profile_class_init (NcDensityProfileClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->finalize = nc_density_profile_finalize;

  klass->eval_fourier = NULL;
}
