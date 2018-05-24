/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_galaxy_redshift_spec.h
 *
 *  Tue April 17 14:26:11 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <sandro@isoftware.com.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_spec.h
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

#ifndef _NC_GALAXY_REDSHIFT_SPEC_H_
#define _NC_GALAXY_REDSHIFT_SPEC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_galaxy_redshift.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_REDSHIFT_SPEC             (nc_galaxy_redshift_spec_get_type ())
#define NC_GALAXY_REDSHIFT_SPEC(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_GALAXY_REDSHIFT_SPEC, NcGalaxyRedshiftSpec))
#define NC_GALAXY_REDSHIFT_SPEC_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_GALAXY_REDSHIFT_SPEC, NcGalaxyRedshiftSpecClass))
#define NC_IS_GALAXY_REDSHIFT_SPEC(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_GALAXY_REDSHIFT_SPEC))
#define NC_IS_GALAXY_REDSHIFT_SPEC_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_GALAXY_REDSHIFT_SPEC))
#define NC_GALAXY_REDSHIFT_SPEC_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_GALAXY_REDSHIFT_SPEC, NcGalaxyRedshiftSpecClass))

typedef struct _NcGalaxyRedshiftSpecClass NcGalaxyRedshiftSpecClass;
typedef struct _NcGalaxyRedshiftSpec NcGalaxyRedshiftSpec;
typedef struct _NcGalaxyRedshiftSpecPrivate NcGalaxyRedshiftSpecPrivate;

struct _NcGalaxyRedshiftSpecClass
{
	/*< private >*/
	NcGalaxyRedshiftClass parent_class;
};

struct _NcGalaxyRedshiftSpec
{
	/*< private >*/
	NcGalaxyRedshift parent_instance;
	NcGalaxyRedshiftSpecPrivate *priv;
};

GType nc_galaxy_redshift_spec_get_type (void) G_GNUC_CONST;

NcGalaxyRedshiftSpec *nc_galaxy_redshift_spec_new (const gdouble z_spec);
NcGalaxyRedshiftSpec *nc_galaxy_redshift_spec_ref (NcGalaxyRedshiftSpec *gzs);

void nc_galaxy_redshift_spec_free (NcGalaxyRedshiftSpec *gzs);
void nc_galaxy_redshift_spec_clear (NcGalaxyRedshiftSpec **gzs);

void nc_galaxy_redshift_spec_set_z (NcGalaxyRedshiftSpec *gzs, const gdouble z_spec);
gdouble nc_galaxy_redshift_spec_get_z (NcGalaxyRedshiftSpec *gzs);

G_END_DECLS

#endif /* _NC_GALAXY_REDSHIFT_SPEC_H_ */
