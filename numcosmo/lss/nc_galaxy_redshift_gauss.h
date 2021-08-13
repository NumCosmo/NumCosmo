/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_redshift_gauss.h
 *
 *  Mon July 27 14:39:04 2020
 *  Copyright  2020  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <sandro@isoftware.com.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_gauss.h
 * Copyright (C) 2020 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
 * Copyright (C) 2020 Mariana Penna Lima <pennalima@gmail.com>
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

#ifndef _NC_GALAXY_REDSHIFT_GAUSS_H_
#define _NC_GALAXY_REDSHIFT_GAUSS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/lss/nc_galaxy_redshift.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_REDSHIFT_GAUSS             (nc_galaxy_redshift_gauss_get_type ())
#define NC_GALAXY_REDSHIFT_GAUSS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_GALAXY_REDSHIFT_GAUSS, NcGalaxyRedshiftGauss))
#define NC_GALAXY_REDSHIFT_GAUSS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_GALAXY_REDSHIFT_GAUSS, NcGalaxyRedshiftGaussClass))
#define NC_IS_GALAXY_REDSHIFT_GAUSS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_GALAXY_REDSHIFT_GAUSS))
#define NC_IS_GALAXY_REDSHIFT_GAUSS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_GALAXY_REDSHIFT_GAUSS))
#define NC_GALAXY_REDSHIFT_GAUSS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_GALAXY_REDSHIFT_GAUSS, NcGalaxyRedshiftGaussClass))

typedef struct _NcGalaxyRedshiftGaussClass NcGalaxyRedshiftGaussClass;
typedef struct _NcGalaxyRedshiftGauss NcGalaxyRedshiftGauss;
typedef struct _NcGalaxyRedshiftGaussPrivate NcGalaxyRedshiftGaussPrivate;

struct _NcGalaxyRedshiftGaussClass
{
  /*< private >*/
  NcGalaxyRedshiftClass parent_class;
};

struct _NcGalaxyRedshiftGauss
{
  /*< private >*/
  NcGalaxyRedshift parent_instance;
  NcGalaxyRedshiftGaussPrivate *priv;
};

GType nc_galaxy_redshift_gauss_get_type (void) G_GNUC_CONST;

NcGalaxyRedshiftGauss *nc_galaxy_redshift_gauss_new (void);
NcGalaxyRedshiftGauss *nc_galaxy_redshift_gauss_ref (NcGalaxyRedshiftGauss *gzg);

void nc_galaxy_redshift_gauss_free (NcGalaxyRedshiftGauss *gzg);
void nc_galaxy_redshift_gauss_clear (NcGalaxyRedshiftGauss **gzg);

void nc_galaxy_redshift_gauss_set_obs (NcGalaxyRedshiftGauss *gzg, NcmMatrix *obs);
NcmMatrix *nc_galaxy_redshift_gauss_peek_obs (NcGalaxyRedshiftGauss *gzg);

guint nc_galaxy_redshift_gauss_len (NcGalaxyRedshiftGauss *gzg);

G_END_DECLS

#endif /* _NC_GALAXY_REDSHIFT_GAUSS_H_ */

