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
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>

#include <stdio.h>

G_BEGIN_DECLS

#define NCM_TYPE_MSET             (ncm_mset_get_type ())
#define NCM_MSET(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_MSET, NcmMSet))
#define NCM_MSET_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_MSET, NcmMSetClass))
#define NCM_IS_MSET(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_MSET))
#define NCM_IS_MSET_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_MSET))
#define NCM_MSET_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_MSET, NcmMSetClass))

#define NCM_TYPE_MSET_PINDEX (ncm_mset_pindex_get_type ())

typedef struct _NcmMSetClass NcmMSetClass;
typedef struct _NcmMSet NcmMSet;
typedef struct _NcmMSetPIndex NcmMSetPIndex;
typedef struct _NcmMSetModelDesc NcmMSetModelDesc; 

#define NCM_MODEL_MAX_ID 5

struct _NcmMSetModelDesc
{
  /*< private >*/
  gboolean init;
  gchar *ns;
  gchar *desc;
  gchar *long_desc;
};

struct _NcmMSet
{
  /*< private >*/
  GObject parent_instance;
  NcmModel *model[NCM_MODEL_MAX_ID];
  GArray *pi_array;
  GArray *fpi_array[NCM_MODEL_MAX_ID];
  GPtrArray *fullname_parray;
  gboolean valid_map;
  guint total_len;
  guint fparam_len;
};

struct _NcmMSetClass
{
  /*< private >*/
  GObjectClass parent_class;
  GHashTable *ns_table;
  NcmMSetModelDesc model_desc[NCM_MODEL_MAX_ID];
};

struct _NcmMSetPIndex
{
  /*< private >*/
  NcmModelID mid;
  guint pid;
};

GType ncm_mset_get_type (void) G_GNUC_CONST;
GType ncm_mset_pindex_get_type (void) G_GNUC_CONST;

void ncm_mset_model_register_id (NcmModelClass *model_class, gchar *ns, gchar *desc, gchar *long_desc);

/**
 * NCM_MSET_MODEL_ID_FUNC: (skip)
 * @model_ns: FIXME
 * 
 * FIXME
 * 
 */
#define NCM_MSET_MODEL_ID_FUNC(model_ns) model_ns##_id

/**
 * NCM_MSET_MODEL_DECLARE_ID: (skip)
 * @model_ns: FIXME
 * 
 * FIXME
 * 
 */
#define NCM_MSET_MODEL_DECLARE_ID(model_ns) gint32 NCM_MSET_MODEL_ID_FUNC(model_ns) (void) G_GNUC_CONST

/**
 * NCM_MSET_MODEL_REGISTER_ID: (skip)
 * @model_ns: FIXME
 * @typemacro: FIXME
 * 
 * FIXME
 * 
 */
#define NCM_MSET_MODEL_REGISTER_ID(model_ns,typemacro) \
gint32 NCM_MSET_MODEL_ID_FUNC(model_ns) (void) \
{ \
  static gint32 id = -1; \
  if (id == -1) \
  { \
    NcmModelClass *model_class = g_type_class_ref (typemacro); \
    id = model_class->model_id; \
    g_type_class_unref (model_class); \
  } \
  return id; \
}

NcmMSetPIndex *ncm_mset_pindex_new (NcmModelID mid, guint pid);
NcmMSetPIndex *ncm_mset_pindex_dup (NcmMSetPIndex *pi);
void ncm_mset_pindex_free (NcmMSetPIndex *pi);

NcmMSet *ncm_mset_empty_new (void);
NcmMSet *ncm_mset_new (NcmModel *model0, ...);
NcmMSet *ncm_mset_newv (NcmModel *model0, va_list ap);
NcmMSet *ncm_mset_new_array (NcmModel **model);
NcmMSet *ncm_mset_ref (NcmMSet *mset);
NcmMSet *ncm_mset_dup (NcmMSet *mset, NcmSerialize *ser);

void ncm_mset_free (NcmMSet *mset);
void ncm_mset_clear (NcmMSet **mset);

NcmModel *ncm_mset_peek (NcmMSet *mset, NcmModelID mid);
NcmModel *ncm_mset_get (NcmMSet *mset, NcmModelID mid);
void ncm_mset_remove (NcmMSet *mset, NcmModelID mid);
void ncm_mset_set (NcmMSet *mset, NcmModel *model);
gboolean ncm_mset_exists (NcmMSet *mset, NcmModel *model);

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
gboolean ncm_mset_params_valid (NcmMSet *mset);
gboolean ncm_mset_cmp (NcmMSet *mset0, NcmMSet *mset1, gboolean cmp_model);

void ncm_mset_param_set (NcmMSet *mset, NcmModelID mid, guint pid, const gdouble x);
gdouble ncm_mset_param_get (NcmMSet *mset, NcmModelID mid, guint pid);
gdouble ncm_mset_orig_param_get (NcmMSet *mset, NcmModelID mid, guint pid);

guint ncm_mset_param_len (NcmMSet *mset);
const gchar *ncm_mset_param_name (NcmMSet *mset, NcmModelID mid, guint pid);
void ncm_mset_param_set_ftype (NcmMSet *mset, NcmModelID mid, guint pid, NcmParamType ftype);
void ncm_mset_param_set_all_ftype (NcmMSet *mset, NcmParamType ftype);
void ncm_mset_param_set_vector (NcmMSet *mset, NcmVector *params);
void ncm_mset_param_get_vector (NcmMSet *mset, NcmVector *params);
gdouble ncm_mset_param_get_scale (NcmMSet *mset, NcmModelID mid, guint pid);
gdouble ncm_mset_param_get_lower_bound (NcmMSet *mset, NcmModelID mid, guint pid);
gdouble ncm_mset_param_get_upper_bound (NcmMSet *mset, NcmModelID mid, guint pid);
gdouble ncm_mset_param_get_abstol (NcmMSet *mset, NcmModelID mid, guint pid);
NcmParamType ncm_mset_param_get_ftype (NcmMSet *mset, NcmModelID mid, guint pid);

void ncm_mset_param_set_pi (NcmMSet *mset, NcmMSetPIndex *pi, const gdouble *x, guint n);
void ncm_mset_param_get_pi (NcmMSet *mset, NcmMSetPIndex *pi, gdouble *x, guint n);

void ncm_mset_fparams_get_vector (NcmMSet *mset, NcmVector *x);
void ncm_mset_fparams_get_vector_offset (NcmMSet *mset, NcmVector *x, guint offset);
void ncm_mset_fparams_set_vector (NcmMSet *mset, const NcmVector *x);
void ncm_mset_fparams_set_array (NcmMSet *mset, const gdouble *x);
void ncm_mset_fparams_set_gsl_vector (NcmMSet *mset, const gsl_vector *x);

guint ncm_mset_fparams_len (NcmMSet *mset);
const gchar *ncm_mset_fparam_name (NcmMSet *mset, guint n);
const gchar *ncm_mset_fparam_full_name (NcmMSet *mset, guint n);
gdouble ncm_mset_fparam_get_scale (NcmMSet *mset, guint n);
gdouble ncm_mset_fparam_get_lower_bound (NcmMSet *mset, guint n);
gdouble ncm_mset_fparam_get_upper_bound (NcmMSet *mset, guint n);
gdouble ncm_mset_fparam_get_abstol (NcmMSet *mset, guint n);

gdouble ncm_mset_fparam_get (NcmMSet *mset, guint n);
void ncm_mset_fparam_set (NcmMSet *mset, guint n, const gdouble x);

NcmMSetPIndex *ncm_mset_fparam_get_pi (NcmMSet *mset, guint n);
gint ncm_mset_fparam_get_fpi (NcmMSet *mset, NcmModelID mid, guint pid);

void ncm_mset_save (NcmMSet *mset, gchar *filename, gboolean save_comment);
NcmMSet *ncm_mset_load (gchar *filename);

G_END_DECLS

#endif /* _NCM_MSET_H_ */
