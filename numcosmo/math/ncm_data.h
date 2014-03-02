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
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_bootstrap.h>

G_BEGIN_DECLS

#define NCM_TYPE_DATA             (ncm_data_get_type ())
#define NCM_DATA(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_DATA, NcmData))
#define NCM_DATA_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_DATA, NcmDataClass))
#define NCM_IS_DATA(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_DATA))
#define NCM_IS_DATA_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_DATA))
#define NCM_DATA_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_DATA, NcmDataClass))

typedef struct _NcmDataClass NcmDataClass;
typedef struct _NcmData NcmData;

struct _NcmData
{
  /*< private >*/
  GObject parent_instance;
  gchar *desc;
  gboolean init;
  gboolean begin;
  NcmBootstrap *bstrap;
};

/**
 * NcmDataClass:
 * @bootstrap: sets whenever the #NcmData implementations supports bootstrap.
 * @get_length: return the length associated to the #NcmData object.
 * @get_dof: return the effective degrees of freedom related to the #NcmData
 * statistics (likelihood or $\chi^2$) this number does not represent 
 * necessarely the number of data points.
 * @begin: perform any model independent precalculation.
 * @prepare: perform any model dependent precalculation.
 * @resample: resample data from the models in #NcmMSet.
 * @leastsquares_f: calculates the least squares $\vec{f}$ vector, i.e., 
 * $\chi^2 \equiv \vec{f}\cdot\vec{f}$.
 * @leastsquares_J: calculates the least squares $\vec{f}$ vector derivatives
 * with respect to the free parameter of @mset.
 * @leastsquares_f_J: calculates both least squares vector and its derivatives.
 * @m2lnL_val: evaluate the minus two times the natural logarithim of the 
 * likelihood, i.e., $-2\ln(L)$.
 * @m2lnL_grad: evaluate the gradient of $-2\ln(L)$ with respect to the free
 * parameters in @mset.
 * @m2lnL_val_grad: evaluate the value and the gradient of $-2\ln(L)$.
 * 
 * Virtual table for the #NcmData abstract class.
 * 
 */
struct _NcmDataClass
{
  /*< private >*/
  GObjectClass parent_class;
  gchar *name;
  /*< public >*/
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
};

GType ncm_data_get_type (void) G_GNUC_CONST;

NcmData *ncm_data_ref (NcmData *data);
void ncm_data_free (NcmData *data);
void ncm_data_clear (NcmData **data);

NcmData *ncm_data_dup (NcmData *data, NcmSerialize *ser_obj);

guint ncm_data_get_length (NcmData *data);
guint ncm_data_get_dof (NcmData *data);
void ncm_data_set_init (NcmData *data, gboolean state);

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

void ncm_data_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *f);
void ncm_data_leastsquares_J (NcmData *data, NcmMSet *mset, NcmMatrix *J);
void ncm_data_leastsquares_f_J (NcmData *data, NcmMSet *mset, NcmVector *f, NcmMatrix *J);
void ncm_data_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
void ncm_data_m2lnL_grad (NcmData *data, NcmMSet *mset, NcmVector *grad);
void ncm_data_m2lnL_val_grad (NcmData *data, NcmMSet *mset, gdouble *m2lnL, NcmVector *grad);

#define NCM_DATA_RESAMPLE_RNG_NAME "data_resample"

G_END_DECLS

#endif /* _NCM_DATA_H_ */

