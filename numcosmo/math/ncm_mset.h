/***************************************************************************
 *            ncm_mset.h
 *
 *  Fri May 25 09:38:14 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef _NCM_MSET_H_
#define _NCM_MSET_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/math/ncm_model.h>

#include <stdio.h>

G_BEGIN_DECLS

#define NCM_TYPE_MSET             (ncm_mset_get_type ())
#define NCM_MSET(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_MSET, NcmMSet))
#define NCM_MSET_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_MSET, NcmMSetClass))
#define NCM_IS_MSET(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_MSET))
#define NCM_IS_MSET_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_MSET))
#define NCM_MSET_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_MSET, NcmMSetClass))

typedef struct _NcmMSetClass NcmMSetClass;
typedef struct _NcmMSet NcmMSet;
typedef struct _NcmMSetPIndex NcmMSetPIndex;

#define NCM_MODEL_MAX_ID 5

struct _NcmMSetClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcmMSet
{
  /*< private >*/
  GObject parent_instance;
  NcmModel *model[NCM_MODEL_MAX_ID];
  GArray *pi_array;
  GArray *fpi_array[NCM_MODEL_MAX_ID];
  gboolean valid_map;
  guint total_len;
  guint fparam_len;
};

struct _NcmMSetPIndex
{
  /*< private >*/
  NcmModelID gmid;
  guint pid;
};

GType ncm_mset_get_type (void) G_GNUC_CONST;
GType ncm_mset_pindex_get_type (void) G_GNUC_CONST;

NcmMSetPIndex *ncm_mset_pindex_new (void);
NcmMSetPIndex *ncm_mset_pindex_copy (NcmMSetPIndex *pi);
void ncm_mset_pindex_free (NcmMSetPIndex *pi);

NcmMSet *ncm_mset_empty_new (void);
NcmMSet *ncm_mset_new (NcmModel *model0, ...);
NcmMSet *ncm_mset_newv (NcmModel *model0, va_list ap);
NcmMSet *ncm_mset_new_array (NcmModel **model);
NcmMSet *ncm_mset_copy_all (NcmMSet *mset);
NcmMSet *ncm_mset_copy_ref (NcmMSet *mset);

void ncm_mset_copyto (NcmMSet *mset_src, NcmMSet *mset_dest);
void ncm_mset_free (NcmMSet *mset);

NcmModel *ncm_mset_peek (NcmMSet *mset, NcmModelID gmid);
NcmModel *ncm_mset_get (NcmMSet *mset, NcmModelID gmid);
void ncm_mset_remove (NcmMSet *mset, NcmModelID gmid);
void ncm_mset_set (NcmMSet *mset, NcmModel *model);

void ncm_mset_prepare_fparam_map (NcmMSet *mset);

guint ncm_mset_total_len (NcmMSet *mset);
guint ncm_mset_fparam_len (NcmMSet *mset);
guint ncm_mset_max_param_name (NcmMSet *mset);
guint ncm_mset_max_fparam_name (NcmMSet *mset);
guint ncm_mset_max_model_nick (NcmMSet *mset);

void ncm_mset_pretty_log (NcmMSet *mset);
void ncm_mset_params_pretty_print (NcmMSet *mset, FILE *out, gchar *header);
void ncm_mset_params_log_vals (NcmMSet *mset);
void ncm_mset_params_print_vals (NcmMSet *mset, FILE *out);

void ncm_mset_param_set (NcmMSet *mset, NcmModelID gmid, guint pid, const gdouble x);
gdouble ncm_mset_param_get (NcmMSet *mset, NcmModelID gmid, guint pid);
gdouble ncm_mset_orig_param_get (NcmMSet *mset, NcmModelID gmid, guint pid);

const gchar *ncm_mset_param_name (NcmMSet *mset, NcmModelID gmid, guint pid);
void ncm_mset_param_set_ftype (NcmMSet *mset, NcmModelID gmid, guint pid, NcmParamType ftype);
void ncm_mset_param_set_all_ftype (NcmMSet *mset, NcmParamType ftype);
gdouble ncm_mset_param_get_scale (NcmMSet *mset, NcmModelID gmid, guint pid);
gdouble ncm_mset_param_get_lower_bound (NcmMSet *mset, NcmModelID gmid, guint pid);
gdouble ncm_mset_param_get_upper_bound (NcmMSet *mset, NcmModelID gmid, guint pid);
gdouble ncm_mset_param_get_abstol (NcmMSet *mset, NcmModelID gmid, guint pid);
NcmParamType ncm_mset_param_get_ftype (NcmMSet *mset, NcmModelID gmid, guint pid);

void ncm_mset_param_set_pi (NcmMSet *mset, NcmMSetPIndex *pi, const gdouble *x, guint n);
void ncm_mset_param_get_pi (NcmMSet *mset, NcmMSetPIndex *pi, gdouble *x, guint n);

void ncm_mset_fparams_get_vector (NcmMSet *mset, NcmVector *x);
void ncm_mset_fparams_set_vector (NcmMSet *mset, const NcmVector *x);
void ncm_mset_fparams_set_array (NcmMSet *mset, const gdouble *x);
void ncm_mset_fparams_set_gsl_vector (NcmMSet *mset, const gsl_vector *x);

guint ncm_mset_fparams_len (NcmMSet *mset);
const gchar *ncm_mset_fparam_name (NcmMSet *mset, guint n);
gdouble ncm_mset_fparam_get_scale (NcmMSet *mset, guint n);
gdouble ncm_mset_fparam_get_lower_bound (NcmMSet *mset, guint n);
gdouble ncm_mset_fparam_get_upper_bound (NcmMSet *mset, guint n);
gdouble ncm_mset_fparam_get_abstol (NcmMSet *mset, guint n);

gdouble ncm_mset_fparam_get (NcmMSet *mset, guint n);
void ncm_mset_fparam_set (NcmMSet *mset, guint n, const gdouble x);

NcmMSetPIndex *ncm_mset_fparam_get_pi (NcmMSet *mset, guint n);
gint ncm_mset_fparam_get_fpi (NcmMSet *mset, NcmModelID gmid, guint pid);

G_END_DECLS

#endif /* _NCM_MSET_H_ */
