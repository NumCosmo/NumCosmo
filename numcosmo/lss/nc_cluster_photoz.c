/***************************************************************************
 *            nc_cluster_photoz.c
 *
 *  Tue Apr 20 10:59:01 2010
 *  Copyright  2010  Mariana Penna Lima & Sandro Dias Pinto Vitenti
 *  <pennalima@gmail.com> & <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
 * SECTION:nc_cluster_photoz
 * @title: Photoz Cluster
 * @short_description: Photometric redshift
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <glib.h>
#include <gsl/gsl_sf_erf.h>

G_DEFINE_TYPE (NcClusterPhotoz, nc_cluster_photoz, G_TYPE_OBJECT);

/**
 * nc_cluster_photoz_new_from_name:
 * @photoz_name: string which specifies the type of the photometric redshift distribution.
 *
 * This function returns a new #NcClusterPhotoz whose type is defined by @photoz_name.
 *
 * Returns: A new #NcClusterPhotoz.
 */
NcClusterPhotoz *
nc_cluster_photoz_new_from_name (gchar *photoz_name)
{
  GType photoz_type = g_type_from_name (photoz_name);
  if (photoz_type == 0)
  {
	g_message ("# Invalid photometric redshift name %s\n", photoz_name);
	g_error ("Aborting...");
  }
  else if (!g_type_is_a (photoz_type, NC_TYPE_CLUSTER_PHOTOZ))
	g_error ("nc_cluster_photoz_new_from_name: NcClusterPhotoz %s do not descend from %s\n", photoz_name, g_type_name (NC_TYPE_CLUSTER_PHOTOZ));
  return g_object_new (photoz_type, NULL);
}

/**
 * nc_cluster_photoz_dist_eval:
 * @photo: a #NcClusterPhotoz.
 * @pz_data: FIXME
 * @z_photo: photometric redshift.
 * @z_real: real redshift.
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_cluster_photoz_dist_eval (NcClusterPhotoz *photo, gdouble *pz_data, gdouble z_photo, gdouble z_real)
{
  return NC_CLUSTER_PHOTOZ_GET_CLASS (photo)->dist_eval (photo, pz_data, z_photo, z_real);
}

/**
 * nc_cluster_photoz_resample:
 * @photo: a #NcClusterPhotoz.
 * @pz_data: FIXME.
 * @z_real: real redshift.
 * @result: FIXME
 *
 * FIXME
 * The function which will call this one is responsible to allocate enough memory for @z_lower and @z_upper.
 *
*/
void
nc_cluster_photoz_resample (NcClusterPhotoz *photo, gdouble *pz_data, gdouble z_real, gdouble *result)
{
  NC_CLUSTER_PHOTOZ_GET_CLASS (photo)->resample (photo, pz_data, z_real, result);
}

/**
 * nc_cluster_photoz_integ_limits:
 * @photo: a #NcClusterPhotoz.
 * @pz_data: FIXME.
 * @z_photo: photometric redshift.
 * @z_lower: pointer to the lower limit of the real redshift integration.
 * @z_upper: pointer to the upper limit of the real redshift integration.
 *
 * FIXME
 * The function which will call this one is responsible to allocate memory for @z_lower and @z_upper.
*/
void
nc_cluster_photoz_integ_limits (NcClusterPhotoz *photo, gdouble *pz_data, gdouble z_photo, gdouble *z_lower, gdouble *z_upper)
{
  NC_CLUSTER_PHOTOZ_GET_CLASS (photo)->integ_limits (photo, pz_data, z_photo, z_lower, z_upper);
}

/**
 * nc_cluster_photoz_free:
 * @photo: a #NcClusterPhotoz.
 *
 * Atomically decrements the reference count of @photo by one. If the reference count drops to 0,
 * all memory allocated by @photo is released.
 *
 */
void
nc_cluster_photoz_free (NcClusterPhotoz *photo)
{
  g_clear_object (&photo);
}

static void
nc_cluster_photoz_init (NcClusterPhotoz *nc_cluster_photoz)
{
  /* TODO: Add initialization code here */
}

static void
_nc_cluster_photoz_finalize (GObject *object)
{
  /* TODO: Add deinitalization code here */

  G_OBJECT_CLASS (nc_cluster_photoz_parent_class)->finalize (object);
}

static void
nc_cluster_photoz_class_init (NcClusterPhotozClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  //GObjectClass* parent_class = G_OBJECT_CLASS (klass);

  object_class->finalize = _nc_cluster_photoz_finalize;
}

