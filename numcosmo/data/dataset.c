/***************************************************************************
 *            dataset.c
 *
 *  Tue May 29 19:28:48 2007
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

/**
 * SECTION:dataset
 * @title: Data Set module
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <glib.h>

G_DEFINE_BOXED_TYPE (NcDataSet, nc_dataset, nc_dataset_copy, nc_dataset_free);

/**
 * nc_dataset_new:
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcDataSet *
nc_dataset_new ()
{
  NcDataSet *ds = g_slice_new (NcDataSet);
  ds->data_list = NULL;
  ds->clone = FALSE;
  return ds;
}

/**
 * nc_dataset_copy:
 * @ds_orig: pointer to type defined by #NcDataSet
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcDataSet *
nc_dataset_copy (NcDataSet *ds_orig)
{
  NcDataSet *ds = g_slice_new (NcDataSet);
  ds->clone = TRUE;
  ds->data_list = g_list_copy (ds_orig->data_list);
  return ds;
}


/**
 * nc_dataset_append_data:
 * @ds: pointer to type defined by #NcDataSet
 * @data: #NcData object to be appended to #NcDataSet
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
nc_dataset_append_data (NcDataSet *ds, NcData *data)
{
  ds->data_list = g_list_append (ds->data_list, data);
}

/**
 * nc_dataset_get_n:
 * @ds: pointer to type defined by #NcDataSet
 *
 * Calculate the total number of data set points
 *
 * Returns: FIXME
 */
guint
nc_dataset_get_n (NcDataSet *ds)
{
  guint n = 0;
  GList *data_list;

  if (ds->data_list != NULL)
  {
    data_list = g_list_first (ds->data_list);
    while (data_list)
    {
      NcData *data = (NcData *)data_list->data;
      n += NC_DATA_LENGTH (data);
      data_list = g_list_next (data_list);
    }
  }
  return n;
}

/**
 * nc_dataset_all_init:
 * @ds: pointer to type defined by #NcDataSet
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
nc_dataset_all_init (NcDataSet *ds)
{
  GList *data_list;
  gboolean init = TRUE;

  if (ds->data_list != NULL)
  {
    data_list = g_list_first (ds->data_list);
    while (data_list)
    {
      NcData *data = (NcData *)data_list->data;
      init = (init && data->init);
      data_list = g_list_next (data_list);
    }
  }
  else
    init = FALSE;

  return init;
}

/**
 * nc_dataset_get_ndata:
 * @ds: pointer to type defined by #NcDataSet
 *
 * FIXME
 *
 * Returns: number of NcData objects in the set
 */
guint
nc_dataset_get_ndata (NcDataSet *ds)
{
  g_assert (ds->data_list != NULL);
  return g_list_length (ds->data_list);
}

/**
 * nc_dataset_get_data:
 * @ds: pointer to type defined by #NcDataSet
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcData *
nc_dataset_get_data (NcDataSet *ds, guint n)
{
  g_assert (ds->data_list != NULL);
  g_assert (n < g_list_length (ds->data_list));
  return (NcData *) g_list_nth_data (ds->data_list, n);
}

/***************************************************************************
 *
 *
 ****************************************************************************/

/**
 * nc_dataset_free0:
 * @ds: pointer to type defined by #NcDataSet
 * @free_all: FIXME
 *
 * FIXME
 */
void
nc_dataset_free0 (NcDataSet *ds, gboolean free_all)
{
  GList *data_list;
  if (ds->clone)
  {
    g_slice_free (NcDataSet, ds);
    return;
  }

  if (free_all)
  {
    data_list = g_list_first (ds->data_list);
    while (data_list)
    {
      NcData *data = (NcData *)data_list->data;
      nc_data_free0 (data, TRUE);
      data_list = g_list_next (data_list);
    }
  }

  g_list_free (ds->data_list);
  g_slice_free (NcDataSet, ds);

  return;
}

/**
 * nc_dataset_free:
 * @ds: pointer to type defined by #NcDataSet
 *
 * FIXME
 */
void
nc_dataset_free (NcDataSet *ds)
{
  nc_dataset_free0 (ds, FALSE);
}

/**
 * nc_dataset_resample:
 * @ds: a #NcDataSet.
 * @mset: a #NcmMSet.
 * @save: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
nc_dataset_resample (NcDataSet *ds, NcmMSet *mset, gboolean save)
{
  GList *data_list;
  data_list = g_list_first (ds->data_list);
  while (data_list)
  {
    NcData *data = (NcData *)data_list->data;
    nc_data_resample (data, mset, save);
    data_list = g_list_next (data_list);
  }
  return TRUE;
}

/**
 * nc_dataset_set_orig:
 * @ds: pointer to type defined by #NcDataSet
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
nc_dataset_set_orig (NcDataSet *ds)
{
  GList *data_list;
  data_list = g_list_first (ds->data_list);
  while (data_list)
  {
    NcData *data = (NcData *)data_list->data;
    nc_data_set_orig (data);
    data_list = g_list_next (data_list);
  }
  return TRUE;
}

/**
 * nc_dataset_log_info:
 * @ds: pointer to type defined by #NcDataSet
 *
 * FIXME
 */
void
nc_dataset_log_info (NcDataSet *ds)
{
  guint ndata = nc_dataset_get_ndata (ds);
  guint i;

  g_message ("#----------------------------------------------------------------------------------\n");
  g_message ("# Data used:\n");
  for (i = 0; i < ndata; i++)
	g_message ("#   - %s\n", nc_dataset_get_data (ds, i)->name);

  return;
}
