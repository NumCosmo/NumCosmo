/***************************************************************************
 *            ncm_dataset.h
 *
 *  Mon Jul 16 18:04:45 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
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

#ifndef _NCM_DATASET_H_
#define _NCM_DATASET_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data.h>
#include <numcosmo/math/ncm_obj_array.h>

G_BEGIN_DECLS

#define NCM_TYPE_DATASET             (ncm_dataset_get_type ())
#define NCM_DATASET(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_DATASET, NcmDataset))
#define NCM_DATASET_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_DATASET, NcmDatasetClass))
#define NCM_IS_DATASET(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_DATASET))
#define NCM_IS_DATASET_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_DATASET))
#define NCM_DATASET_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_DATASET, NcmDatasetClass))

typedef struct _NcmDatasetClass NcmDatasetClass;
typedef struct _NcmDataset NcmDataset;

/**
 * NcmDatasetBStrapType:
 * @NCM_DATASET_BSTRAP_DISABLE: Bootstrap disabled.
 * @NCM_DATASET_BSTRAP_PARTIAL: Partial bootstrap, each #NcmData is bootstraped individually.
 * @NCM_DATASET_BSTRAP_TOTAL: Total bootstrap, all data is bootstraped simultaneously.
 *
 */
typedef enum _NcmDatasetBStrapType
{
  NCM_DATASET_BSTRAP_DISABLE = 0,
  NCM_DATASET_BSTRAP_PARTIAL,
  NCM_DATASET_BSTRAP_TOTAL, /*< private >*/
  NCM_DATASET_BSTRAP_LEN,   /*< skip >*/
} NcmDatasetBStrapType;

struct _NcmDatasetClass
{
  /*< private >*/
  GObjectClass parent_class;
};

/**
 * NcmDataSet:
 *
 * FIXME
 */
struct _NcmDataset
{
  /*< private >*/
  GObject parent_instance;
  NcmObjArray *oa;
  NcmDatasetBStrapType bstype;
  GArray *data_prob;
  GArray *bstrap;
};

GType ncm_dataset_get_type (void) G_GNUC_CONST;

NcmDataset *ncm_dataset_new (void);
NcmDataset *ncm_dataset_dup (NcmDataset *dset, NcmSerialize *ser);
NcmDataset *ncm_dataset_ref (NcmDataset *dset);
NcmDataset *ncm_dataset_copy (NcmDataset *dset);
void ncm_dataset_free (NcmDataset *dset);
void ncm_dataset_clear (NcmDataset **dset);

guint ncm_dataset_get_length (NcmDataset *dset);
guint ncm_dataset_get_n (NcmDataset *dset);
guint ncm_dataset_get_dof (NcmDataset *dset);
gboolean ncm_dataset_all_init (NcmDataset *dset);

void ncm_dataset_append_data (NcmDataset *dset, NcmData *data);
NcmData *ncm_dataset_get_data (NcmDataset *dset, guint n);
NcmData *ncm_dataset_peek_data (NcmDataset *dset, guint n);
guint ncm_dataset_get_ndata (NcmDataset *dset);
void ncm_dataset_set_data_array (NcmDataset *dset, NcmObjArray *oa);
NcmObjArray *ncm_dataset_get_data_array (NcmDataset *dset);
NcmObjArray *ncm_dataset_peek_data_array (NcmDataset *dset);

void ncm_dataset_resample (NcmDataset *dset, NcmMSet *mset, NcmRNG *rng);
void ncm_dataset_bootstrap_set (NcmDataset *dset, NcmDatasetBStrapType bstype);
void ncm_dataset_bootstrap_resample (NcmDataset *dset, NcmRNG *rng);
void ncm_dataset_log_info (NcmDataset *dset);
gchar *ncm_dataset_get_info (NcmDataset *dset);

gboolean ncm_dataset_has_leastsquares_f (NcmDataset *dset);
gboolean ncm_dataset_has_leastsquares_J (NcmDataset *dset);
gboolean ncm_dataset_has_leastsquares_f_J (NcmDataset *dset);

gboolean ncm_dataset_has_m2lnL_val (NcmDataset *dset);
gboolean ncm_dataset_has_m2lnL_grad (NcmDataset *dset);
gboolean ncm_dataset_has_m2lnL_val_grad (NcmDataset *dset);

void ncm_dataset_leastsquares_f (NcmDataset *dset, NcmMSet *mset, NcmVector *f);
void ncm_dataset_leastsquares_J (NcmDataset *dset, NcmMSet *mset, NcmMatrix *J);
void ncm_dataset_leastsquares_f_J (NcmDataset *dset, NcmMSet *mset, NcmVector *f, NcmMatrix *J);

void ncm_dataset_m2lnL_val (NcmDataset *dset, NcmMSet *mset, gdouble *m2lnL);
void ncm_dataset_m2lnL_vec (NcmDataset *dset, NcmMSet *mset, NcmVector *m2lnL_v);
void ncm_dataset_m2lnL_grad (NcmDataset *dset, NcmMSet *mset, NcmVector *grad);
void ncm_dataset_m2lnL_val_grad (NcmDataset *dset, NcmMSet *mset, gdouble *m2lnL, NcmVector *grad);

void ncm_dataset_m2lnL_i_val (NcmDataset *dset, NcmMSet *mset, guint i, gdouble *m2lnL_i);

void ncm_dataset_fisher_matrix (NcmDataset *dset, NcmMSet *mset, NcmMatrix **IM);

G_END_DECLS

#endif /* _NCM_DATASET_H_ */
