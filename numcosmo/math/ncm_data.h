/***************************************************************************
 *            ncm_data.h
 *
 *  Sat Mar 29 15:51:59 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_DATA_H_
#define _NCM_DATA_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/math/ncm_mset.h>

G_BEGIN_DECLS

#define NCM_TYPE_DATA             (ncm_data_get_type ())
#define NCM_DATA(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_DATA, NcmData))
#define NCM_DATA_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_DATA, NcmDataClass))
#define NCM_IS_DATA(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_DATA))
#define NCM_IS_DATA_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_DATA))
#define NCM_DATA_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_DATA, NcmDataClass))

typedef struct _NcmDataClass NcmDataClass;
typedef struct _NcmData NcmData;

struct _NcmDataClass
{
  /*< private >*/
  GObjectClass parent_class;
  gchar *name;
  guint (*get_length) (NcmData *data);
  NcmData *(*dup) (NcmData *data);
  void (*copyto) (NcmData *data, NcmData *data_dest);
  void (*begin) (NcmData *data);
  void (*prepare) (NcmData *data, NcmMSet *mset);
  void (*resample) (NcmData *data, NcmMSet *mset);
  void (*leastsquares_f) (NcmData *data, NcmMSet *mset, NcmVector *f);
  void (*leastsquares_J) (NcmData *data, NcmMSet *mset, NcmMatrix *J);
  void (*leastsquares_f_J) (NcmData *data, NcmMSet *mset, NcmVector *f, NcmMatrix *J);
  void (*m2lnL_val) (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
  void (*m2lnL_grad) (NcmData *data, NcmMSet *mset, NcmVector *grad);
  void (*m2lnL_val_grad) (NcmData *data, NcmMSet *mset, gdouble *m2lnL, NcmVector *grad);
};

struct _NcmData
{
  /*< private >*/
  GObject parent_instance;
  gchar *desc;
  gboolean init;
  gboolean begin;
};

GType ncm_data_get_type (void) G_GNUC_CONST;

NcmData *ncm_data_ref (NcmData *data);
void ncm_data_free (NcmData *data);
void ncm_data_clear (NcmData **data);

NcmData *ncm_data_dup (NcmData *data);

guint ncm_data_get_length (NcmData *data);
void ncm_data_set_init (NcmData *data);

void ncm_data_prepare (NcmData *data, NcmMSet *mset);
void ncm_data_resample (NcmData *data, NcmMSet *mset);
void ncm_data_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *f);
void ncm_data_leastsquares_J (NcmData *data, NcmMSet *mset, NcmMatrix *J);
void ncm_data_leastsquares_f_J (NcmData *data, NcmMSet *mset, NcmVector *f, NcmMatrix *J);
void ncm_data_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
void ncm_data_m2lnL_grad (NcmData *data, NcmMSet *mset, NcmVector *grad);
void ncm_data_m2lnL_val_grad (NcmData *data, NcmMSet *mset, gdouble *m2lnL, NcmVector *grad);

gchar *ncm_data_get_desc (NcmData *data);

G_END_DECLS

#endif /* _NCM_DATA_H_ */

