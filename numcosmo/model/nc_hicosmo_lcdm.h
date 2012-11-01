/***************************************************************************
 *            nc_hicosmo_lcdm.h
 *
 *  Mon Aug 11 19:56:24 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef _NC_HICOSMO_LCDM_H_
#define _NC_HICOSMO_LCDM_H_

#include <glib-object.h>
#include <numcosmo/nc_hicosmo.h>

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_LCDM             (nc_hicosmo_lcdm_get_type ())
#define NC_HICOSMO_LCDM(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_LCDM, NcHICosmoLCDM))
#define NC_HICOSMO_LCDM_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_LCDM, NcHICosmoLCDMClass))
#define NC_IS_HICOSMO_LCDM(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_LCDM))
#define NC_IS_HICOSMO_LCDM_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_LCDM))
#define NC_HICOSMO_LCDM_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_LCDM, NcHICosmoLCDMClass))

typedef struct _NcHICosmoLCDMClass NcHICosmoLCDMClass;
typedef struct _NcHICosmoLCDM NcHICosmoLCDM;

struct _NcHICosmoLCDMClass
{
  /*< private >*/
  NcHICosmoClass parent_class;
};

struct _NcHICosmoLCDM
{
  /*< private >*/
  NcHICosmo parent_instance;
};

GType nc_hicosmo_lcdm_get_type (void) G_GNUC_CONST;

NcHICosmoLCDM *nc_hicosmo_lcdm_new (void);

G_END_DECLS

#endif /* _NC_HICOSMO_LCDM_H_ */
