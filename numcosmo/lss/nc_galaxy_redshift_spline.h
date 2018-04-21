/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_galaxy_redshift_spline.h
 *
 *  Thu April 19 14:39:04 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <sandro@isoftware.com.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_spline.h
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

#ifndef _NC_GALAXY_REDSHIFT_SPLINE_H_
#define _NC_GALAXY_REDSHIFT_SPLINE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_galaxy_redshift.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_REDSHIFT_SPLINE             (nc_galaxy_redshift_spline_get_type ())
#define NC_GALAXY_REDSHIFT_SPLINE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_GALAXY_REDSHIFT_SPLINE, NcGalaxyRedshiftSpline))
#define NC_GALAXY_REDSHIFT_SPLINE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_GALAXY_REDSHIFT_SPLINE, NcGalaxyRedshiftSplineClass))
#define NC_IS_GALAXY_REDSHIFT_SPLINE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_GALAXY_REDSHIFT_SPLINE))
#define NC_IS_GALAXY_REDSHIFT_SPLINE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_GALAXY_REDSHIFT_SPLINE))
#define NC_GALAXY_REDSHIFT_SPLINE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_GALAXY_REDSHIFT_SPLINE, NcGalaxyRedshiftSplineClass))

typedef struct _NcGalaxyRedshiftSplineClass NcGalaxyRedshiftSplineClass;
typedef struct _NcGalaxyRedshiftSpline NcGalaxyRedshiftSpline;
typedef struct _NcGalaxyRedshiftSplinePrivate NcGalaxyRedshiftSplinePrivate;

struct _NcGalaxyRedshiftSplineClass
{
	/*< private >*/
	NcGalaxyRedshiftClass parent_class;
};

struct _NcGalaxyRedshiftSpline
{
	/*< private >*/	
	NcGalaxyRedshift parent_instance;
	NcGalaxyRedshiftSplinePrivate *priv;
};

GType nc_galaxy_redshift_spline_get_type (void) G_GNUC_CONST;

NcGalaxyRedshiftSpline *nc_galaxy_redshift_spline_new (void);
NcGalaxyRedshiftSpline *nc_galaxy_redshift_spline_ref (NcGalaxyRedshiftSpline *gzs);

void nc_galaxy_redshift_spline_free (NcGalaxyRedshiftSpline *gzs);
void nc_galaxy_redshift_spline_clear (NcGalaxyRedshiftSpline **gzs);

void nc_galaxy_redshift_spline_set_z_best (NcGalaxyRedshiftSpline *gzs, const gdouble z_best);
gdouble nc_galaxy_redshift_spline_get_z_best (NcGalaxyRedshiftSpline *gzs);

void nc_galaxy_redshift_spline_init_from_vectors (NcGalaxyRedshiftSpline *gzs, NcmVector *zv, NcmVector *Pzv);

#define NC_GALAXY_REDSHIFT_SPLINE_LKNOT_DROP (-2.0 * M_LN10)

G_END_DECLS

#endif /* _NC_GALAXY_REDSHIFT_SPLINE_H_ */
