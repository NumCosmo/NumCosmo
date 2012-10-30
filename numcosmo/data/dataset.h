/***************************************************************************
 *            dataset.h
 *
 *  Mon Jul 16 18:04:45 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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

#ifndef _NC_DATASET_H
#define _NC_DATASET_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/data/data.h>

G_BEGIN_DECLS

GType nc_dataset_get_type (void);

typedef struct _NcDataSet NcDataSet;

/**
 * NcDataSet:
 *
 * FIXME
 */
struct _NcDataSet
{
  /*< private >*/
  GList *data_list;
  gboolean clone;
};

NcDataSet *nc_dataset_new ();
NcDataSet *nc_dataset_copy (NcDataSet *ds_orig);
void nc_dataset_free (NcDataSet *ds);
void nc_dataset_free0 (NcDataSet *ds, gboolean free_all);

guint nc_dataset_get_n (NcDataSet *ds);
gboolean nc_dataset_all_init (NcDataSet *ds);

void nc_dataset_append_data (NcDataSet *ds, NcData *data);
NcData *nc_dataset_get_data (NcDataSet *ds, guint n);
guint nc_dataset_get_ndata (NcDataSet *ds);

gboolean nc_dataset_resample (NcDataSet *ds, NcmMSet *mset, gboolean save);
gboolean nc_dataset_set_orig (NcDataSet *ds);

void nc_dataset_log_info (NcDataSet *ds);

G_END_DECLS

#endif /* _DATASET_H */
