/***************************************************************************
 *            nc_powspec_mnl.h
 *
 *  Thu February 18 12:32:25 2016
 *  Copyright  2016  Cyrille Doux
 *  <cdoux@apc.in2p3.fr>
 ****************************************************************************/
/*
 * nc_powspec_mnl.h
 * Copyright (C) 2016 Cyrille Doux <cdoux@apc.in2p3.fr>
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

#ifndef _NC_POWSPEC_MNL_H_
#define _NC_POWSPEC_MNL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_powspec.h>

G_BEGIN_DECLS

#define NC_TYPE_POWSPEC_MNL             (nc_powspec_mnl_get_type ())
#define NC_POWSPEC_MNL(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_POWSPEC_MNL, NcPowspecMNL))
#define NC_POWSPEC_MNL_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_POWSPEC_MNL, NcPowspecMNLClass))
#define NC_IS_POWSPEC_MNL(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_POWSPEC_MNL))
#define NC_IS_POWSPEC_MNL_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_POWSPEC_MNL))
#define NC_POWSPEC_MNL_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_POWSPEC_MNL, NcPowspecMNLClass))

typedef struct _NcPowspecMNLClass NcPowspecMNLClass;
typedef struct _NcPowspecMNL NcPowspecMNL;

struct _NcPowspecMNLClass
{
  /*< private > */
  NcmPowspecClass parent_class;
};

struct _NcPowspecMNL
{
  /*< private > */
  NcmPowspec parent_instance;
};

GType nc_powspec_mnl_get_type (void) G_GNUC_CONST;

NcPowspecMNL *nc_powspec_mnl_new_from_name (const gchar *ps_mnl_name);
NcPowspecMNL *nc_powspec_mnl_ref (NcPowspecMNL *ps_mnl);

void nc_powspec_mnl_free (NcPowspecMNL *ps_mnl);
void nc_powspec_mnl_clear (NcPowspecMNL **ps_mnl);

G_END_DECLS

#endif /* _NC_POWSPEC_MNL_H_ */

