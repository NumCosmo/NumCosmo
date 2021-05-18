/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_wl_reduced_shear_gauss.h
 *
 *  Mon July 27 11:13:56 2020
 *  Copyright  2020  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <sandro@isoftware.com.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_galaxy_wl_reduced_shear_gauss.h
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

#ifndef _NC_GALAXY_WL_REDUCED_SHEAR_GAUSS_H_
#define _NC_GALAXY_WL_REDUCED_SHEAR_GAUSS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/lss/nc_galaxy_wl_dist.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_WL_REDUCED_SHEAR_GAUSS             (nc_galaxy_wl_reduced_shear_gauss_get_type ())
#define NC_GALAXY_WL_REDUCED_SHEAR_GAUSS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_GALAXY_WL_REDUCED_SHEAR_GAUSS, NcGalaxyWLReducedShearGauss))
#define NC_GALAXY_WL_REDUCED_SHEAR_GAUSS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_GALAXY_WL_REDUCED_SHEAR_GAUSS, NcGalaxyWLReducedShearGaussClass))
#define NC_IS_GALAXY_WL_REDUCED_SHEAR_GAUSS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_GALAXY_WL_REDUCED_SHEAR_GAUSS))
#define NC_IS_GALAXY_WL_REDUCED_SHEAR_GAUSS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_GALAXY_WL_REDUCED_SHEAR_GAUSS))
#define NC_GALAXY_WL_REDUCED_SHEAR_GAUSS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_GALAXY_WL_REDUCED_SHEAR_GAUSS, NcGalaxyWLReducedShearGaussClass))

typedef struct _NcGalaxyWLReducedShearGaussClass NcGalaxyWLReducedShearGaussClass;
typedef struct _NcGalaxyWLReducedShearGauss NcGalaxyWLReducedShearGauss;
typedef struct _NcGalaxyWLReducedShearGaussPrivate NcGalaxyWLReducedShearGaussPrivate;

struct _NcGalaxyWLReducedShearGaussClass
{
  /*< private >*/
  NcGalaxyWLDistClass parent_class;
};

struct _NcGalaxyWLReducedShearGauss
{
  /*< private >*/
  NcGalaxyWLDist parent_instance;
  NcGalaxyWLReducedShearGaussPrivate *priv;
};

/**
 * NcGalaxyWLReducedShearGaussPos:
 * @NC_GALAXY_WL_REDUCED_SHEAR_GAUSS_POS_ANG: FIXME
 * @NC_GALAXY_WL_REDUCED_SHEAR_GAUSS_POS_R: FIXME
 *
 * FIXME
 *
 */
typedef enum _NcGalaxyWLReducedShearGaussPos
{
  NC_GALAXY_WL_REDUCED_SHEAR_GAUSS_POS_ANG,
  NC_GALAXY_WL_REDUCED_SHEAR_GAUSS_POS_R,
  /* < private > */
  NC_GALAXY_WL_REDUCED_SHEAR_GAUSS_POS_LEN, /*< skip >*/
} NcGalaxyWLReducedShearGaussPos;

GType nc_galaxy_wl_reduced_shear_gauss_get_type (void) G_GNUC_CONST;

NcGalaxyWLReducedShearGauss *nc_galaxy_wl_reduced_shear_gauss_new (NcGalaxyWLReducedShearGaussPos pos);
NcGalaxyWLReducedShearGauss *nc_galaxy_wl_reduced_shear_gauss_ref (NcGalaxyWLReducedShearGauss *grsg);

void nc_galaxy_wl_reduced_shear_gauss_free (NcGalaxyWLReducedShearGauss *grsg);
void nc_galaxy_wl_reduced_shear_gauss_clear (NcGalaxyWLReducedShearGauss **grsg);

void nc_galaxy_wl_reduced_shear_gauss_set_pos (NcGalaxyWLReducedShearGauss *grsg, NcGalaxyWLReducedShearGaussPos pos);
NcGalaxyWLReducedShearGaussPos nc_galaxy_wl_reduced_shear_gauss_get_pos (NcGalaxyWLReducedShearGauss *grsg);

void nc_galaxy_wl_reduced_shear_gauss_set_obs (NcGalaxyWLReducedShearGauss *grsg, NcmMatrix *obs);
NcmMatrix *nc_galaxy_wl_reduced_shear_gauss_peek_obs (NcGalaxyWLReducedShearGauss *grsg);

G_END_DECLS

#endif /* _NC_GALAXY_WL_REDUCED_SHEAR_GAUSS_H_ */

