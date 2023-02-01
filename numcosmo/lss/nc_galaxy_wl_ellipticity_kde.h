/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_wl_ellipticity_kde.h
 *
 *  Mon July 27 11:13:56 2020
 *  Copyright  2020  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <sandro@isoftware.com.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_galaxy_wl_ellipticity_kde.h
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

#ifndef _NC_GALAXY_WL_ELLIPTICITY_KDE_H_
#define _NC_GALAXY_WL_ELLIPTICITY_KDE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/lss/nc_galaxy_wl_dist.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_WL_ELLIPTICITY_KDE             (nc_galaxy_wl_ellipticity_kde_get_type ())
#define NC_GALAXY_WL_ELLIPTICITY_KDE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_GALAXY_WL_ELLIPTICITY_KDE, NcGalaxyWLEllipticityKDE))
#define NC_GALAXY_WL_ELLIPTICITY_KDE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_GALAXY_WL_ELLIPTICITY_KDE, NcGalaxyWLEllipticityKDEClass))
#define NC_IS_GALAXY_WL_ELLIPTICITY_KDE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_GALAXY_WL_ELLIPTICITY_KDE))
#define NC_IS_GALAXY_WL_ELLIPTICITY_KDE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_GALAXY_WL_ELLIPTICITY_KDE))
#define NC_GALAXY_WL_ELLIPTICITY_KDE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_GALAXY_WL_ELLIPTICITY_KDE, NcGalaxyWLEllipticityKDEClass))

typedef struct _NcGalaxyWLEllipticityKDEClass NcGalaxyWLEllipticityKDEClass;
typedef struct _NcGalaxyWLEllipticityKDE NcGalaxyWLEllipticityKDE;
typedef struct _NcGalaxyWLEllipticityKDEPrivate NcGalaxyWLEllipticityKDEPrivate;

struct _NcGalaxyWLEllipticityKDEClass
{
  /*< private >*/
  NcGalaxyWLDistClass parent_class;
};

struct _NcGalaxyWLEllipticityKDE
{
  /*< private >*/
  NcGalaxyWLDist parent_instance;
  NcGalaxyWLEllipticityKDEPrivate *priv;
};

GType nc_galaxy_wl_ellipticity_kde_get_type (void) G_GNUC_CONST;

NcGalaxyWLEllipticityKDE *nc_galaxy_wl_ellipticity_kde_new ();
NcGalaxyWLEllipticityKDE *nc_galaxy_wl_ellipticity_kde_ref (NcGalaxyWLEllipticityKDE *gekde);

void nc_galaxy_wl_ellipticity_kde_free (NcGalaxyWLEllipticityKDE *gekde);
void nc_galaxy_wl_ellipticity_kde_clear (NcGalaxyWLEllipticityKDE **gekde);

void nc_galaxy_wl_ellipticity_kde_set_obs (NcGalaxyWLEllipticityKDE *gekde, NcmMatrix *obs);
NcmMatrix *nc_galaxy_wl_ellipticity_kde_peek_obs (NcGalaxyWLEllipticityKDE *gekde);

G_END_DECLS

#endif /* _NC_GALAXY_WL_ELLIPTICITY_KDE_H_ */

