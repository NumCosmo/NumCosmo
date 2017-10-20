/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_hipert_grav_einstein.h
 *
 *  Fri October 13 10:36:58 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_grav_einstein.h
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

#ifndef _NC_HIPERT_GRAV_EINSTEIN_H_
#define _NC_HIPERT_GRAV_EINSTEIN_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/perturbations/nc_hipert_grav.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_GRAV_EINSTEIN             (nc_hipert_grav_einstein_get_type ())
#define NC_HIPERT_GRAV_EINSTEIN(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT_GRAV_EINSTEIN, NcHIPertGravEinstein))
#define NC_HIPERT_GRAV_EINSTEIN_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPERT_GRAV_EINSTEIN, NcHIPertGravEinsteinClass))
#define NC_IS_HIPERT_GRAV_EINSTEIN(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT_GRAV_EINSTEIN))
#define NC_IS_HIPERT_GRAV_EINSTEIN_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPERT_GRAV_EINSTEIN))
#define NC_HIPERT_GRAV_EINSTEIN_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPERT_GRAV_EINSTEIN, NcHIPertGravEinsteinClass))

typedef struct _NcHIPertGravEinsteinClass NcHIPertGravEinsteinClass;
typedef struct _NcHIPertGravEinstein NcHIPertGravEinstein;
typedef struct _NcHIPertGravEinsteinPrivate NcHIPertGravEinsteinPrivate;

struct _NcHIPertGravEinsteinClass
{
  /*< private >*/
  NcHIPertGravClass parent_class;
};

struct _NcHIPertGravEinstein
{
  /*< private >*/
  NcHIPertGrav parent_instance;
  NcHIPertGravEinsteinPrivate *priv;
};

GType nc_hipert_grav_einstein_get_type (void) G_GNUC_CONST;

NC_HIPERT_BG_VAR_ID_FUNC_DECL (nc_hipert_grav_einstein);

NcHIPertGravEinstein *nc_hipert_grav_einstein_new (void);
NcHIPertGravEinstein *nc_hipert_grav_einstein_ref (NcHIPertGravEinstein *gr);

void nc_hipert_grav_einstein_free (NcHIPertGravEinstein *gr);
void nc_hipert_grav_einstein_clear (NcHIPertGravEinstein **gr);

G_END_DECLS

#endif /* _NC_HIPERT_GRAV_EINSTEIN_H_ */
