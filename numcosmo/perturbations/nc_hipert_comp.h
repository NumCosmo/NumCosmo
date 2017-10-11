/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_hipert_comp.h
 *
 *  Wed October 11 15:54:27 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_comp.h
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

#ifndef _NC_HIPERT_COMP_H_
#define _NC_HIPERT_COMP_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_COMP             (nc_hipert_comp_get_type ())
#define NC_HIPERT_COMP(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT_COMP, NcHIPertComp))
#define NC_HIPERT_COMP_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPERT_COMP, NcHIPertCompClass))
#define NC_IS_HIPERT_COMP(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT_COMP))
#define NC_IS_HIPERT_COMP_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPERT_COMP))
#define NC_HIPERT_COMP_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPERT_COMP, NcHIPertCompClass))

typedef struct _NcHIPertCompClass NcHIPertCompClass;
typedef struct _NcHIPertComp NcHIPertComp;
typedef struct _NcHIPertCompPrivate NcHIPertCompPrivate;

struct _NcHIPertCompClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcHIPertComp
{
  /*< private >*/
  GObject parent_instance;
  NcHIPertCompPrivate *priv;
};

GType nc_hipert_comp_get_type (void) G_GNUC_CONST;

G_END_DECLS

#endif /* _NC_HIPERT_COMP_H_ */

