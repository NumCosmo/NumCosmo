/***************************************************************************
 *            ncm_data.h
 *
 *  Sat Mar 29 15:51:59 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/math/ncm_diff.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_bootstrap.h>

G_BEGIN_DECLS

#define NCM_TYPE_DATA (ncm_data_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmData, ncm_data, NCM, DATA, GObject)

/**
 * NcmDataFisherMatrix:
 * @data: a #NcmData
 * @mset: a #NcmMSet
 * @IM: (out): The fisher matrix
 *
 * Calculates the Fisher-information matrix @I.
 *
 */
typedef void (*NcmDataFisherMatrix) (NcmData *data, NcmMSet *mset, NcmMatrix **IM);

struct _NcmDataClass
{
  /*< private >*/
  GObjectClass parent_class;
  gchar *name;
  gboolean bootstrap;

  guint (*get_length) (NcmData *data);
  guint (*get_dof) (NcmData *data);
  void (*begin) (NcmData *data);
  void (*prepare) (NcmData *data, NcmMSet *mset);
  void (*resample) (NcmData *data, NcmMSet *mset, NcmRNG *rng);
  void (*leastsquares_f) (NcmData *data, NcmMSet *mset, NcmVector *f);
  void (*leastsquares_J) (NcmData *data, NcmMSet *mset, NcmMatrix *J);
  void (*leastsquares_f_J) (NcmData *data, NcmMSet *mset, NcmVector *f, NcmMatrix *J);
  void (*m2lnL_val) (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
  void (*m2lnL_grad) (NcmData *data, NcmMSet *mset, NcmVector *grad);
  void (*m2lnL_val_grad) (NcmData *data, NcmMSet *mset, gdouble *m2lnL, NcmVector *grad);
  void (*mean_vector) (NcmData *data, NcmMSet *mset, NcmVector *mu);
  void (*inv_cov_UH) (NcmData *data, NcmMSet *mset, NcmMatrix *H);

  NcmDataFisherMatrix fisher_matrix;
  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[2];
};

NcmData *ncm_data_ref (NcmData *data);
void ncm_data_free (NcmData *data);
void ncm_data_clear (NcmData **data);

NcmData *ncm_data_dup (NcmData *data, NcmSerialize *ser_obj);
NcmData *ncm_data_new_from_file (const gchar *filename);

guint ncm_data_get_length (NcmData *data);
guint ncm_data_get_dof (NcmData *data);
void ncm_data_set_init (NcmData *data, gboolean state);
gboolean ncm_data_is_init (NcmData *data);

void ncm_data_set_desc (NcmData *data, const gchar *desc);
void ncm_data_take_desc (NcmData *data, gchar *desc);
const gchar *ncm_data_peek_desc (NcmData *data);
gchar *ncm_data_get_desc (NcmData *data);

void ncm_data_prepare (NcmData *data, NcmMSet *mset);
void ncm_data_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng);

void ncm_data_bootstrap_create (NcmData *data);
void ncm_data_bootstrap_remove (NcmData *data);
void ncm_data_bootstrap_set (NcmData *data, NcmBootstrap *bstrap);
void ncm_data_bootstrap_resample (NcmData *data, NcmRNG *rng);
gboolean ncm_data_bootstrap_enabled (NcmData *data);
NcmBootstrap *ncm_data_peek_bootstrap (NcmData *data);

void ncm_data_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *f);
void ncm_data_leastsquares_J (NcmData *data, NcmMSet *mset, NcmMatrix *J);
void ncm_data_leastsquares_f_J (NcmData *data, NcmMSet *mset, NcmVector *f, NcmMatrix *J);
void ncm_data_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
void ncm_data_m2lnL_grad (NcmData *data, NcmMSet *mset, NcmVector *grad);
void ncm_data_m2lnL_val_grad (NcmData *data, NcmMSet *mset, gdouble *m2lnL, NcmVector *grad);

void ncm_data_mean_vector (NcmData *data, NcmMSet *mset, NcmVector *mu);
void ncm_data_sigma_vector (NcmData *data, NcmMSet *mset, NcmVector *sigma);
void ncm_data_inv_cov_UH (NcmData *data, NcmMSet *mset, NcmMatrix *H);

void ncm_data_fisher_matrix (NcmData *data, NcmMSet *mset, NcmMatrix **IM);

#define NCM_DATA_RESAMPLE_RNG_NAME "data_resample"

G_END_DECLS

#endif /* _NCM_DATA_H_ */

