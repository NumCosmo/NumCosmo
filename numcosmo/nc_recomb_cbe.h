/***************************************************************************
 *            nc_recomb_cbe.h
 *
 *  Wed November 18 15:27:12 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>>
 ****************************************************************************/
/*
 * nc_recomb_cbe.h
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_RECOMB_CBE_H_
#define _NC_RECOMB_CBE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_recomb.h>

#include <numcosmo/nc_cbe.h>

G_BEGIN_DECLS

#define NC_TYPE_RECOMB_CBE             (nc_recomb_cbe_get_type ())
#define NC_RECOMB_CBE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_RECOMB_CBE, NcRecombCBE))
#define NC_RECOMB_CBE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_RECOMB_CBE, NcRecombCBEClass))
#define NC_IS_RECOMB_CBE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_RECOMB_CBE))
#define NC_IS_RECOMB_CBE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_RECOMB_CBE))
#define NC_RECOMB_CBE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_RECOMB_CBE, NcRecombCBEClass))

typedef struct _NcRecombCBEClass NcRecombCBEClass;
typedef struct _NcRecombCBE NcRecombCBE;

struct _NcRecombCBEClass
{
  /*< private >*/
  NcRecombClass parent_class;
};

struct _NcRecombCBE
{
  /*< private >*/
  NcRecomb parent_instance;
  NcCBE *cbe;
  NcmSpline *Xe_s;
};

GType nc_recomb_cbe_get_type (void) G_GNUC_CONST;

NcRecombCBE *nc_recomb_cbe_new (void);
NcRecombCBE *nc_recomb_cbe_full_new (NcCBE *cbe);
NcRecombCBE *nc_recomb_cbe_ref (NcRecombCBE *recomb_cbe);
void nc_recomb_cbe_free (NcRecombCBE *recomb_cbe);
void nc_recomb_cbe_clear (NcRecombCBE **recomb_cbe);

void nc_recomb_cbe_set_cbe (NcRecombCBE *recomb_cbe, NcCBE *cbe);
NcCBE *nc_recomb_cbe_peek_cbe (NcRecombCBE *recomb_cbe);

G_END_DECLS

#endif /* _NC_RECOMB_CBE_H_ */

