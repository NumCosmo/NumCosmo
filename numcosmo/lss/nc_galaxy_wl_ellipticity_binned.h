/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_wl_ellipticity_binned.h
 *
 *  Fri February 24 10:19:37 2023
 *  Copyright  2023  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <vitenti@uel.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_galaxy_wl_ellipticity_binned.h
 * Copyright (C) 2023 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
 * Copyright (C) 2023 Mariana Penna Lima <pennalima@gmail.com>
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

#ifndef _NC_GALAXY_WL_ELLIPTICITY_BINNED_H_
#define _NC_GALAXY_WL_ELLIPTICITY_BINNED_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/lss/nc_galaxy_wl_dist.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_WL_ELLIPTICITY_BINNED             (nc_galaxy_wl_ellipticity_binned_get_type ())
#define NC_GALAXY_WL_ELLIPTICITY_BINNED(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_GALAXY_WL_ELLIPTICITY_BINNED, NcGalaxyWLEllipticityBinned))
#define NC_GALAXY_WL_ELLIPTICITY_BINNED_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_GALAXY_WL_ELLIPTICITY_BINNED, NcGalaxyWLEllipticityBinnedClass))
#define NC_IS_GALAXY_WL_ELLIPTICITY_BINNED(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_GALAXY_WL_ELLIPTICITY_BINNED))
#define NC_IS_GALAXY_WL_ELLIPTICITY_BINNED_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_GALAXY_WL_ELLIPTICITY_BINNED))
#define NC_GALAXY_WL_ELLIPTICITY_BINNED_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_GALAXY_WL_ELLIPTICITY_BINNED, NcGalaxyWLEllipticityBinnedClass))

typedef struct _NcGalaxyWLEllipticityBinnedClass NcGalaxyWLEllipticityBinnedClass;
typedef struct _NcGalaxyWLEllipticityBinned NcGalaxyWLEllipticityBinned;
typedef struct _NcGalaxyWLEllipticityBinnedPrivate NcGalaxyWLEllipticityBinnedPrivate;

struct _NcGalaxyWLEllipticityBinnedClass
{
  /*< private >*/
  NcGalaxyWLDistClass parent_class;
};

struct _NcGalaxyWLEllipticityBinned
{
  /*< private >*/
  NcGalaxyWLDist parent_instance;
  NcGalaxyWLEllipticityBinnedPrivate *priv;
};

GType nc_galaxy_wl_ellipticity_binned_get_type (void) G_GNUC_CONST;

NcGalaxyWLEllipticityBinned *nc_galaxy_wl_ellipticity_binned_new ();
NcGalaxyWLEllipticityBinned *nc_galaxy_wl_ellipticity_binned_ref (NcGalaxyWLEllipticityBinned *gebin);

void nc_galaxy_wl_ellipticity_binned_free (NcGalaxyWLEllipticityBinned *gebin);
void nc_galaxy_wl_ellipticity_binned_clear (NcGalaxyWLEllipticityBinned **gebin);

void nc_galaxy_wl_ellipticity_binned_set_binobs (NcGalaxyWLEllipticityBinned *gebin, NcmMatrix *obs, NcmVector *bins);
NcmObjArray *nc_galaxy_wl_ellipticity_binned_peek_binobs (NcGalaxyWLEllipticityBinned *gebin);
NcmVector *nc_galaxy_wl_ellipticity_binned_peek_bins (NcGalaxyWLEllipticityBinned *gebin);

G_END_DECLS

#endif /* _NC_GALAXY_WL_ELLIPTICITY_BINNED_H_ */

