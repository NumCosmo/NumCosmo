/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_galaxy_redshift.h
 *
 *  Tue April 17 11:13:56 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <sandro@isoftware.com.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_galaxy_redshift.h
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
 * Copyright (C) 2018 Mariana Penna Lima <pennalima@gmail.com>
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

#ifndef _NC_GALAXY_REDSHIFT_H_
#define _NC_GALAXY_REDSHIFT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_REDSHIFT             (nc_galaxy_redshift_get_type ())
#define NC_GALAXY_REDSHIFT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_GALAXY_REDSHIFT, NcGalaxyRedshift))
#define NC_GALAXY_REDSHIFT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_GALAXY_REDSHIFT, NcGalaxyRedshiftClass))
#define NC_IS_GALAXY_REDSHIFT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_GALAXY_REDSHIFT))
#define NC_IS_GALAXY_REDSHIFT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_GALAXY_REDSHIFT))
#define NC_GALAXY_REDSHIFT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_GALAXY_REDSHIFT, NcGalaxyRedshiftClass))

typedef struct _NcGalaxyRedshiftClass NcGalaxyRedshiftClass;
typedef struct _NcGalaxyRedshift NcGalaxyRedshift;
typedef struct _NcGalaxyRedshiftPrivate NcGalaxyRedshiftPrivate;

struct _NcGalaxyRedshiftClass
{
  /*< private >*/
	GObjectClass parent_class;
	gboolean (*has_dist) (NcGalaxyRedshift *gz);
	gdouble (*mode) (NcGalaxyRedshift *gz);
	guint (*nintervals) (NcGalaxyRedshift *gz);
	gdouble (*interval_weight) (NcGalaxyRedshift *gz, const guint di);
	void (*pdf_limits) (NcGalaxyRedshift *gz, const guint di, gdouble *zmin, gdouble *zmax);
	gdouble (*pdf) (NcGalaxyRedshift *gz, const guint di, const gdouble z);
	gdouble (*gen) (NcGalaxyRedshift *gz, NcmRNG *rng);
};

struct _NcGalaxyRedshift
{
	/*< private >*/
	GObject parent_instance;
	NcGalaxyRedshiftPrivate *priv;
};

GType nc_galaxy_redshift_get_type (void) G_GNUC_CONST;

NcGalaxyRedshift *nc_galaxy_redshift_ref (NcGalaxyRedshift *gz);

void nc_galaxy_redshift_free (NcGalaxyRedshift *gz);
void nc_galaxy_redshift_clear (NcGalaxyRedshift **gz);

G_INLINE_FUNC gboolean nc_galaxy_redshift_has_dist (NcGalaxyRedshift *gz);
G_INLINE_FUNC gdouble nc_galaxy_redshift_mode (NcGalaxyRedshift *gz);
G_INLINE_FUNC guint nc_galaxy_redshift_nintervals (NcGalaxyRedshift *gz);
G_INLINE_FUNC gdouble nc_galaxy_redshift_interval_weight (NcGalaxyRedshift *gz, const guint di);
G_INLINE_FUNC void nc_galaxy_redshift_pdf_limits (NcGalaxyRedshift *gz, const guint di, gdouble *zmin, gdouble *zmax);
G_INLINE_FUNC gdouble nc_galaxy_redshift_pdf (NcGalaxyRedshift *gz, const guint di, const gdouble z);
G_INLINE_FUNC gdouble nc_galaxy_redshift_gen (NcGalaxyRedshift *gz, NcmRNG *rng);

G_END_DECLS

#endif /* _NC_GALAXY_REDSHIFT_H_ */

#ifndef _NC_GALAXY_REDSHIFT_INLINE_H_
#define _NC_GALAXY_REDSHIFT_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

G_INLINE_FUNC gboolean 
nc_galaxy_redshift_has_dist (NcGalaxyRedshift *gz)
{
	return NC_GALAXY_REDSHIFT_GET_CLASS (gz)->has_dist (gz);
}

G_INLINE_FUNC gdouble 
nc_galaxy_redshift_mode (NcGalaxyRedshift *gz)
{
	return NC_GALAXY_REDSHIFT_GET_CLASS (gz)->mode (gz);
}

G_INLINE_FUNC guint 
nc_galaxy_redshift_nintervals (NcGalaxyRedshift *gz)
{
	return NC_GALAXY_REDSHIFT_GET_CLASS (gz)->nintervals (gz);
}

G_INLINE_FUNC gdouble 
nc_galaxy_redshift_interval_weight (NcGalaxyRedshift *gz, const guint di)
{
	return NC_GALAXY_REDSHIFT_GET_CLASS (gz)->interval_weight (gz, di);
}

G_INLINE_FUNC void 
nc_galaxy_redshift_pdf_limits (NcGalaxyRedshift *gz, const guint di, gdouble *zmin, gdouble *zmax)
{
	return NC_GALAXY_REDSHIFT_GET_CLASS (gz)->pdf_limits (gz, di, zmin, zmax);
}

G_INLINE_FUNC gdouble 
nc_galaxy_redshift_pdf (NcGalaxyRedshift *gz, const guint di, const gdouble z)
{
	return NC_GALAXY_REDSHIFT_GET_CLASS (gz)->pdf (gz, di, z);
}

G_INLINE_FUNC gdouble 
nc_galaxy_redshift_gen (NcGalaxyRedshift *gz, NcmRNG *rng)
{
	return NC_GALAXY_REDSHIFT_GET_CLASS (gz)->gen (gz, rng);
}

G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_GALAXY_REDSHIFT_INLINE_H_ */
