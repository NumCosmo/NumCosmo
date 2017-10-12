/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_hipert_grav.h
 *
 *  Thu October 12 14:32:08 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_grav.h
 * Copyright (C) 2017 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_HIPERT_GRAV_H_
#define _NC_HIPERT_GRAV_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_GRAV             (nc_hipert_grav_get_type ())
#define NC_HIPERT_GRAV(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT_GRAV, NcHIPertGrav))
#define NC_HIPERT_GRAV_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPERT_GRAV, NcHIPertGravClass))
#define NC_IS_HIPERT_GRAV(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT_GRAV))
#define NC_IS_HIPERT_GRAV_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPERT_GRAV))
#define NC_HIPERT_GRAV_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPERT_GRAV, NcHIPertGravClass))

typedef struct _NcHIPertGravClass NcHIPertGravClass;
typedef struct _NcHIPertGrav NcHIPertGrav;
typedef struct _NcHIPertGravPrivate NcHIPertGravPrivate;

struct _NcHIPertGravClass
{
  GObjectClass parent_class;
};

struct _NcHIPertGrav
{
  GObject parent_instance;
  NcHIPertGravPrivate *priv;
};

/**
 * NcHIPertGravGauge:
 * @NC_HIPERT_GRAV_GAUGE_NEWTONIAN: Newtonian gauge
 * @NC_HIPERT_GRAV_GAUGE_SYNCHRONOUS: Synchronous gauge
 * 
 * FIXME
 * 
 */
typedef enum /*< enum,underscore_name=NC_HIPERT_GRAV_GAUGE >*/
{
  NC_HIPERT_GRAV_GAUGE_NEWTONIAN,
  NC_HIPERT_GRAV_GAUGE_SYNCHRONOUS,
  /* < private > */
  NC_HIPERT_GRAV_GAUGE_LEN, /*< skip >*/
} NcHIPertGravGauge;

GType nc_hipert_grav_get_type (void) G_GNUC_CONST;

G_END_DECLS

#endif /* _NC_HIPERT_GRAV_H_ */

