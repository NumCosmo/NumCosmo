/***************************************************************************
 *            nc_galaxy_acf.c
 *
 *  Fri May 11 21:18:21 2012
 *  Copyright  2012  Fernando de Simoni & Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Fernando de Simoni & Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef _NC_GALAXY_ACF_H_
#define _NC_GALAXY_ACF_H_

#include <glib-object.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_ACF             (nc_galaxy_acf_get_type ())
#define NC_GALAXY_ACF(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_GALAXY_ACF, NcGalaxyAcf))
#define NC_GALAXY_ACF_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_GALAXY_ACF, NcGalaxyAcfClass))
#define NC_IS_GALAXY_ACF(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_GALAXY_ACF))
#define NC_IS_GALAXY_ACF_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_GALAXY_ACF))
#define NC_GALAXY_ACF_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_GALAXY_ACF, NcGalaxyAcfClass))

typedef struct _NcGalaxyAcfClass NcGalaxyAcfClass;
typedef struct _NcGalaxyAcf NcGalaxyAcf;

struct _NcGalaxyAcfClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcGalaxyAcf
{
  /*< private >*/
  GObject parent_instance;
  NcGrowthFunc *gf;
  NcDistance *dist;
  NcTransferFunc *tf;
  NcmSpline *s;
  NcMatterVar *vp;
  gdouble b;
};

GType nc_galaxy_acf_get_type (void) G_GNUC_CONST;

NcGalaxyAcf *nc_galaxy_acf_new (NcGrowthFunc *gf, NcDistance *dist, NcTransferFunc *tf);
gdouble ncm_galaxy_acf_psi (NcGalaxyAcf *acf, NcHICosmo *model, gdouble k, guint l);
void ncm_galaxy_acf_prepare_psi (NcGalaxyAcf *acf, NcHICosmo *model, guint l);


G_END_DECLS

#endif /* _NC_GALAXY_ACF_H_ */
