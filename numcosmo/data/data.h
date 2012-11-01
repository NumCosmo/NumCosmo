/***************************************************************************
 *            data.h
 *
 *  Sat Mar 29 15:51:59 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
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

#ifndef _NC_DATA_H
#define _NC_DATA_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/math/ncm_mset.h>

G_BEGIN_DECLS

GType nc_data_get_type (void);
GType nc_data_struct_get_type (void);

typedef struct _NcDataStruct NcDataStruct;

/**
 * NcDataStruct:
 *
 * FIXME
 */
struct _NcDataStruct
{
  /*< private >*/
  gpointer data;
  gpointer (*dup) (gpointer data);
  void (*copy) (gpointer dest, gpointer data);
  void (*free) (gpointer data);
  void (*begin) (gpointer data);
  void (*end) (gpointer data);
  guint (*get_length) (gpointer data);
};

NcDataStruct *nc_data_struct_new (void);
NcDataStruct *nc_data_struct_copy (NcDataStruct *dts);
void nc_data_struct_free (NcDataStruct *dts);

#define NC_DATA_STRUCT_DATA(dts) ((dts)->data)
#define NC_DATA_STRUCT_COPY(dest,src) ((src)->copy ((dest)->data, (src)->data))
#define NC_DATA_STRUCT_BEGIN(dts) ((dts)->begin ((dts)->data))
#define NC_DATA_STRUCT_END(dts) ((dts)->end ((dts)->data))
#define NC_DATA_STRUCT_LENGTH(dts) ((dts)->get_length ((dts)->data))

#define NC_DATA_STRUCT_HAS_DUP(dts) ((dts)->dup)
#define NC_DATA_STRUCT_HAS_COPY(dts) ((dts)->copy)
#define NC_DATA_STRUCT_HAS_FREE(dts) ((dts)->free)
#define NC_DATA_STRUCT_HAS_BEGIN(dts) ((dts)->begin)
#define NC_DATA_STRUCT_HAS_END(dts) ((dts)->end)

typedef GDestroyNotify NcDataFree;
typedef gpointer (*NcDataRef) (gpointer model);
typedef void (*NcDataInitModel) (gpointer model, gpointer data);
typedef void (*NcDataPrepare) (NcmMSet *mset, gpointer model, gpointer data);
typedef void (*NcDataResample) (NcmMSet *mset, gpointer model, gpointer data);
typedef void (*NcDataLeastSquaresF) (NcmMSet *mset, gpointer model, gpointer data, NcmVector *f);
typedef void (*NcDataLeastSquaresJ) (NcmMSet *mset, gpointer model, gpointer data, NcmMatrix *J);
typedef void (*NcDataLeastSquaresFJ) (NcmMSet *mset, gpointer model, gpointer data, NcmVector *f, NcmMatrix *J);
typedef void (*NcDataM2lnLVal) (NcmMSet *mset, gpointer model, gpointer data, gdouble *m2lnL);
typedef void (*NcDataM2lnLGrad) (NcmMSet *mset, gpointer model, gpointer data, NcmVector *grad);
typedef void (*NcDataM2lnLValGrad) (NcmMSet *mset, gpointer model, gpointer data, gdouble *m2lnL, NcmVector *grad);

typedef struct _NcData NcData;

/**
 * NcData:
 *
 * FIXME
 */
struct _NcData
{
  /*< private >*/
  gchar *name;
  glong type;
  gboolean init;
  gboolean clone;
  gpointer model;
  NcDataRef model_ref;
  NcDataFree model_free;
  NcDataStruct *dts;
  NcDataResample resample;
  NcDataInitModel model_init;
  NcDataPrepare prepare;
  NcDataLeastSquaresF calc_leastsquares_f;
  NcDataLeastSquaresJ calc_leastsquares_J;
  NcDataLeastSquaresFJ calc_leastsquares_f_J;
  NcDataM2lnLVal calc_m2lnL_val;
  NcDataM2lnLGrad calc_m2lnL_grad;
  NcDataM2lnLValGrad calc_m2lnL_val_grad;
  NcDataStruct *orig;
};

#define NC_DATA_MODEL(data) ((data)->model)
#define NC_DATA_DATA(data) (NC_DATA_STRUCT_DATA((data)->dts))
#define NC_DATA_LENGTH(dt) (NC_DATA_STRUCT_LENGTH((dt)->dts))
#define NC_DATA_HAS_PREPARE(data) ((data)->prepare != NULL)
#define NC_DATA_HAS_MODEL_INIT(data) ((data)->model_init != NULL)
#define NC_DATA_PREPARE(data,cp) ((data)->prepare ((cp),NC_DATA_MODEL ((data)), NC_DATA_DATA ((data))))
#define NC_DATA_MODEL_INIT(data) ((data)->model_init (NC_DATA_MODEL ((data)), NC_DATA_DATA ((data))))

NcData *nc_data_new (void);
NcData *nc_data_copy (NcData *data);
void nc_data_free (NcData *data);
gboolean nc_data_model_init (NcData *data);
gboolean nc_data_begin (NcData *data);
gboolean nc_data_init (NcData *data);
gboolean nc_data_copyto (NcData *dest, NcData *data);
gboolean nc_data_resample (NcData *data, NcmMSet *mset, gboolean save);
gboolean nc_data_bootstrap (NcData *data);
gboolean nc_data_set_orig (NcData *data);

G_END_DECLS

#endif /* DATA_H */
