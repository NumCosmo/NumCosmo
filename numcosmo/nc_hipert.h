/***************************************************************************
 *            nc_hipert.h
 *
 *  Tue June 03 15:48:13 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert.h
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_HIPERT_H_
#define _NC_HIPERT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <nvector/nvector_serial.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT             (nc_hipert_get_type ())
#define NC_HIPERT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT, NcHIPert))
#define NC_HIPERT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPERT, NcHIPertClass))
#define NC_IS_HIPERT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT))
#define NC_IS_HIPERT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPERT))
#define NC_HIPERT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPERT, NcHIPertClass))

typedef struct _NcHIPertClass NcHIPertClass;
typedef struct _NcHIPert NcHIPert;

struct _NcHIPertClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcHIPert
{
  /*< private >*/
  GObject parent_instance;
  gdouble alpha0;
  N_Vector y;
  gpointer cvode;
  gboolean cvode_init;
  gdouble reltol;
  gdouble abstol;
  gdouble k;
  gboolean prepared;
};

GType nc_hipert_get_type (void) G_GNUC_CONST;

void nc_hipert_set_mode_k (NcHIPert *pert, gdouble k);

G_END_DECLS

#endif /* _NC_HIPERT_H_ */
