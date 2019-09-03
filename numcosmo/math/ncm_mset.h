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
#include <numcosmo/math/ncm_obj_array.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <stdio.h>
#endif /* NUMCOSMO_GIR_SCAN */

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

#define NCM_MSET_MAX_STACKSIZE 1000
#define NCM_MSET_INIT_MARRAY 32
#define NCM_MSET_GET_BASE_MID(mid) (mid / NCM_MSET_MAX_STACKSIZE)
#define NCM_MSET_MID(id,pos) ((id) + pos)

struct _NcmMSetModelDesc
{
  /*< private >*/
  gboolean init;
  gchar *ns;
  gchar *desc;
  gchar *long_desc;
};

struct _NcmMSetClass
{
  /*< private >*/
  GObjectClass parent_class;
  GHashTable *ns_table;
  GArray *model_desc_array;
};

struct _NcmMSet
{
  /*< private >*/
  GObject parent_instance;
  NcmObjArray *model_array;
  GHashTable *mid_item_hash;
  GHashTable *model_item_hash;
  GHashTable *fpi_hash;
  GPtrArray *fullname_parray;
  GArray *pi_array;
  GRegex *fullname_regex;
  gboolean valid_map;
  guint total_len;
  guint fparam_len;
  NcmVector *temp_fparams;
};

struct _NcmMSetPIndex
{
  /*< private >*/
  NcmModelID mid;
  guint pid;
};

GType ncm_mset_get_type (void) G_GNUC_CONST;
GType ncm_mset_pindex_get_type (void) G_GNUC_CONST;

void ncm_mset_model_register_id (NcmModelClass *model_class, const gchar *ns, const gchar *desc, const gchar *long_desc, gboolean can_stack, NcmModelID main_model_id);

/**
 * NCM_MSET_MODEL_MAIN: (skip)
 *
 * FIXME
 *
 */
#define NCM_MSET_MODEL_MAIN (-1)

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
#define NCM_MSET_MODEL_DECLARE_ID(model_ns) NcmModelID NCM_MSET_MODEL_ID_FUNC(model_ns) (void) G_GNUC_CONST

/**
 * NCM_MSET_MODEL_REGISTER_ID: (skip)
 * @model_ns: FIXME
 * @typemacro: FIXME
 *
 * FIXME
 *
 */
#define NCM_MSET_MODEL_REGISTER_ID(model_ns,typemacro) \
NcmModelID NCM_MSET_MODEL_ID_FUNC(model_ns) (void) \
{ \
  static NcmModelID id = -1; \
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
NcmMSet *ncm_mset_new (gpointer model0, ...) G_GNUC_NULL_TERMINATED;
NcmMSet *ncm_mset_newv (gpointer model0, va_list ap);
NcmMSet *ncm_mset_new_array (GPtrArray *model_array);
NcmMSet *ncm_mset_ref (NcmMSet *mset);
NcmMSet *ncm_mset_dup (NcmMSet *mset, NcmSerialize *ser);
NcmMSet *ncm_mset_shallow_copy (NcmMSet *mset);

void ncm_mset_free (NcmMSet *mset);
void ncm_mset_clear (NcmMSet **mset);

NcmModel *ncm_mset_peek (NcmMSet *mset, NcmModelID mid);
NcmModel *ncm_mset_peek_pos (NcmMSet *mset, NcmModelID base_mid, guint stackpos_id);
NcmModel *ncm_mset_get (NcmMSet *mset, NcmModelID mid);
NcmModel *ncm_mset_peek_array_pos (NcmMSet *mset, guint i);
NcmModelID ncm_mset_get_mid_array_pos (NcmMSet *mset, guint i);
void ncm_mset_remove (NcmMSet *mset, NcmModelID mid);
void ncm_mset_set (NcmMSet *mset, NcmModel *model);
void ncm_mset_push (NcmMSet *mset, NcmModel *model);
void ncm_mset_set_pos (NcmMSet *mset, NcmModel *model, guint stackpos_id);
gboolean ncm_mset_exists (NcmMSet *mset, NcmModel *model);
gboolean ncm_mset_exists_pos (NcmMSet *mset, NcmModel *model, guint stackpos_id);
gboolean ncm_mset_is_subset (NcmMSet *mset, NcmMSet *sub_mset);
gint ncm_mset_cmp_all (NcmMSet *mset0, NcmMSet *mset1);

NcmModelID ncm_mset_get_id_by_type (GType model_type);
NcmModelID ncm_mset_get_id_by_ns (const gchar *ns);
const gchar *ncm_mset_get_ns_by_id (NcmModelID id);
GType ncm_mset_get_type_by_id (NcmModelID id);

void ncm_mset_set_fmap (NcmMSet *mset, const gchar *const *fmap, gboolean update_models);
gchar **ncm_mset_get_fmap (NcmMSet *mset);
void ncm_mset_prepare_fparam_map (NcmMSet *mset);

guint ncm_mset_total_len (NcmMSet *mset);
guint ncm_mset_fparam_len (NcmMSet *mset);
guint ncm_mset_max_param_name (NcmMSet *mset);
guint ncm_mset_max_fparam_name (NcmMSet *mset);
guint ncm_mset_max_model_nick (NcmMSet *mset);
guint ncm_mset_nmodels (NcmMSet *mset);

void ncm_mset_pretty_log (NcmMSet *mset);
void ncm_mset_params_pretty_print (NcmMSet *mset, FILE *out, const gchar *header);
void ncm_mset_params_log_vals (NcmMSet *mset);
void ncm_mset_params_print_vals (NcmMSet *mset, FILE *out);
void ncm_mset_fparams_log_covar (NcmMSet *mset, NcmMatrix *covar);

gboolean ncm_mset_params_valid (NcmMSet *mset);
gboolean ncm_mset_params_valid_bounds (NcmMSet *mset);
gboolean ncm_mset_cmp (NcmMSet *mset0, NcmMSet *mset1, gboolean cmp_model);

void ncm_mset_param_set (NcmMSet *mset, NcmModelID mid, guint pid, const gdouble x);
gdouble ncm_mset_param_get (NcmMSet *mset, NcmModelID mid, guint pid);
gdouble ncm_mset_orig_param_get (NcmMSet *mset, NcmModelID mid, guint pid);

const gchar *ncm_mset_param_name (NcmMSet *mset, NcmModelID mid, guint pid);
const gchar *ncm_mset_param_symbol (NcmMSet *mset, NcmModelID mid, guint pid);
void ncm_mset_param_set_ftype (NcmMSet *mset, NcmModelID mid, guint pid, NcmParamType ftype);
void ncm_mset_param_set_all_ftype (NcmMSet *mset, NcmParamType ftype);
void ncm_mset_param_set_mid_ftype (NcmMSet *mset, NcmModelID mid, NcmParamType ftype);
void ncm_mset_param_set_all_but_mid_ftype (NcmMSet *mset, NcmModelID mid, NcmParamType ftype);
void ncm_mset_param_set_ftype_from_fmap (NcmMSet *mset);
void ncm_mset_param_set_vector (NcmMSet *mset, NcmVector *params);
void ncm_mset_param_get_vector (NcmMSet *mset, NcmVector *params);
void ncm_mset_param_set_mset (NcmMSet *mset_dest, NcmMSet *mset_src);
gdouble ncm_mset_param_get_scale (NcmMSet *mset, NcmModelID mid, guint pid);
gdouble ncm_mset_param_get_lower_bound (NcmMSet *mset, NcmModelID mid, guint pid);
gdouble ncm_mset_param_get_upper_bound (NcmMSet *mset, NcmModelID mid, guint pid);
gdouble ncm_mset_param_get_abstol (NcmMSet *mset, NcmModelID mid, guint pid);
NcmParamType ncm_mset_param_get_ftype (NcmMSet *mset, NcmModelID mid, guint pid);
void ncm_mset_param_set_scale (NcmMSet *mset, NcmModelID mid, guint pid, gdouble scale);

void ncm_mset_param_set_pi (NcmMSet *mset, NcmMSetPIndex *pi, const gdouble *x, guint n);
void ncm_mset_param_get_pi (NcmMSet *mset, NcmMSetPIndex *pi, gdouble *x, guint n);

void ncm_mset_fparams_get_vector (NcmMSet *mset, NcmVector *x);
void ncm_mset_fparams_get_vector_offset (NcmMSet *mset, NcmVector *x, guint offset);
void ncm_mset_fparams_set_vector (NcmMSet *mset, const NcmVector *x);
void ncm_mset_fparams_set_vector_offset (NcmMSet *mset, const NcmVector *x, guint offset);
void ncm_mset_fparams_set_array (NcmMSet *mset, const gdouble *x);
void ncm_mset_fparams_set_gsl_vector (NcmMSet *mset, const gsl_vector *x);

guint ncm_mset_fparams_len (NcmMSet *mset);
const gchar *ncm_mset_fparam_name (NcmMSet *mset, guint n);
const gchar *ncm_mset_fparam_symbol (NcmMSet *mset, guint n);
const gchar *ncm_mset_fparam_full_name (NcmMSet *mset, guint n);
NcmMSetPIndex *ncm_mset_param_get_by_full_name (NcmMSet *mset, const gchar *fullname);
gdouble ncm_mset_fparam_get_scale (NcmMSet *mset, guint n);
gdouble ncm_mset_fparam_get_lower_bound (NcmMSet *mset, guint n);
gdouble ncm_mset_fparam_get_upper_bound (NcmMSet *mset, guint n);
gdouble ncm_mset_fparam_get_abstol (NcmMSet *mset, guint n);
void ncm_mset_fparam_set_scale (NcmMSet *mset, guint n, gdouble scale);
gboolean ncm_mset_fparam_valid_bounds (NcmMSet *mset, NcmVector *theta);
gboolean ncm_mset_fparam_valid_bounds_offset (NcmMSet *mset, NcmVector *theta, guint offset);
gboolean ncm_mset_fparam_validate_all (NcmMSet *mset, NcmVector *theta);

gdouble ncm_mset_fparam_get (NcmMSet *mset, guint n);
void ncm_mset_fparam_set (NcmMSet *mset, guint n, const gdouble x);

const NcmMSetPIndex *ncm_mset_fparam_get_pi (NcmMSet *mset, guint n);
gint ncm_mset_fparam_get_fpi (NcmMSet *mset, NcmModelID mid, guint pid);
const NcmMSetPIndex *ncm_mset_fparam_get_pi_by_name (NcmMSet *mset, const gchar *name);

void ncm_mset_save (NcmMSet *mset, NcmSerialize *ser, const gchar *filename, gboolean save_comment);
NcmMSet *ncm_mset_load (const gchar *filename, NcmSerialize *ser);

G_END_DECLS

#endif /* _NCM_MSET_H_ */
