/***************************************************************************
 *            ncm_model.h
 *
 *  Fri February 24 21:18:21 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef _NCM_MODEL_H_
#define _NCM_MODEL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_obj_array.h>
#include <numcosmo/math/ncm_sparam.h>
#include <numcosmo/math/ncm_vparam.h>
#include <numcosmo/math/ncm_reparam.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NCM_TYPE_MODEL             (ncm_model_get_type ())
#define NCM_MODEL(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_MODEL, NcmModel))
#define NCM_MODEL_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_MODEL, NcmModelClass))
#define NCM_IS_MODEL(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_MODEL))
#define NCM_IS_MODEL_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_MODEL))
#define NCM_MODEL_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_MODEL, NcmModelClass))

typedef struct _NcmModelClass NcmModelClass;
typedef struct _NcmModel NcmModel;
typedef gint32 NcmModelID;

typedef void (*NcmModelAddSubmodel) (NcmModel *model, NcmModel *submodel);

struct _NcmModelClass
{
  /*< private >*/
  GObjectClass parent_class;
  void (*get_property) (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec);
  void (*set_property) (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec);
  NcmModelAddSubmodel add_submodel;
  gboolean (*valid) (NcmModel *model);
  NcmModelID model_id;
  gboolean can_stack;
  NcmModelID main_model_id;
  gboolean is_submodel;
  gchar *name;
  gchar *nick;
  guint64 impl_flag;
  guint nonparam_prop_len;
  guint sparam_len;
  guint vparam_len;
  guint parent_sparam_len;
  guint parent_vparam_len;
  GPtrArray *sparam;
  GPtrArray *vparam;
};

/**
 * NcmModel:
 *
 * Base class for models.
 */
struct _NcmModel
{
  /*< private >*/
  GObject parent_instance;
  NcmReparam *reparam;
  NcmObjArray *sparams;
  NcmVector *params;
  NcmVector *p;
  GArray *vparam_pos;
  GArray *vparam_len;
  GArray *ptypes;
  GHashTable *sparams_name_id;
  GHashTable *submodel_mid_pos;
  GPtrArray *submodel_array;
  guint total_len;
  guint64 pkey;
  guint64 skey;
};

typedef gdouble (*NcmModelFunc0) (NcmModel *model);
typedef gdouble (*NcmModelFunc1) (NcmModel *model, const gdouble x);
typedef gdouble (*NcmModelFunc2) (NcmModel *model, const gdouble x, const gdouble y);

typedef gdouble (*NcmModelVFunc0) (NcmModel *model, const guint n);
typedef gdouble (*NcmModelVFunc1) (NcmModel *model, const guint n, const gdouble x);
typedef gdouble (*NcmModelVFunc2) (NcmModel *model, const guint n, const gdouble x, const gdouble y);

GType ncm_model_get_type (void) G_GNUC_CONST;

/*void ncm_model_class_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec);*/
/*void ncm_model_class_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec);*/

void ncm_model_class_add_params (NcmModelClass *model_class, guint sparam_len, guint vparam_len, guint nonparam_prop_len);
void ncm_model_class_set_name_nick (NcmModelClass *model_class, const gchar *name, const gchar *nick);

void ncm_model_class_set_sparam_obj (NcmModelClass *model_class, guint sparam_id, NcmSParam *sparam);
void ncm_model_class_set_vparam_obj (NcmModelClass *model_class, guint vparam_id, NcmVParam *vparam);

void ncm_model_class_set_sparam (NcmModelClass *model_class, guint sparam_id, const gchar *symbol, const gchar *name, gdouble lower_bound, gdouble upper_bound, gdouble scale, gdouble abstol, gdouble default_value, NcmParamType ppt);
void ncm_model_class_set_vparam (NcmModelClass *model_class, guint vparam_id, guint default_length, const gchar *symbol, const gchar *name, gdouble lower_bound, gdouble upper_bound, gdouble scale, gdouble abstol, gdouble default_value, NcmParamType ppt);

void ncm_model_class_check_params_info (NcmModelClass *model_class);

void ncm_model_class_add_impl_opts (NcmModelClass *model_class, gint opt1, ...);
void ncm_model_class_add_impl_flag (NcmModelClass *model_class, guint64 flag);

NcmModel *ncm_model_dup (NcmModel *model, NcmSerialize *ser);
NcmModel *ncm_model_ref (NcmModel *model);

void ncm_model_free (NcmModel *model);
void ncm_model_clear (NcmModel **model);
void ncm_model_set_reparam (NcmModel *model, NcmReparam *reparam);
gboolean ncm_model_is_equal (NcmModel *model1, NcmModel *model2);

NCM_INLINE NcmModelID ncm_model_id (NcmModel *model);
NCM_INLINE NcmModelID ncm_model_id_by_type (GType model_type);
NCM_INLINE gboolean ncm_model_check_impl_flag (NcmModel *model, guint64 impl);
NCM_INLINE gboolean ncm_model_check_impl_opt (NcmModel *model, gint opt);
gboolean ncm_model_check_impl_opts (NcmModel *model, gint opt1, ...);
NCM_INLINE guint ncm_model_len (NcmModel *model);
NCM_INLINE gboolean ncm_model_state_is_update (NcmModel *model);
NCM_INLINE void ncm_model_state_set_update (NcmModel *model);
NCM_INLINE void ncm_model_state_mark_outdated (NcmModel *model);

NCM_INLINE guint ncm_model_sparam_len (NcmModel *model);
NCM_INLINE guint ncm_model_vparam_array_len (NcmModel *model);
NCM_INLINE guint ncm_model_vparam_index (NcmModel *model, guint n, guint i);
NCM_INLINE guint ncm_model_vparam_len (NcmModel *model, guint n);
NCM_INLINE const gchar *ncm_model_name (NcmModel *model);
NCM_INLINE const gchar *ncm_model_nick (NcmModel *model);
NCM_INLINE NcmReparam *ncm_model_peek_reparam (NcmModel *model);
NCM_INLINE gboolean ncm_model_params_finite (NcmModel *model);
NCM_INLINE gboolean ncm_model_param_finite (NcmModel *model, guint i);
NCM_INLINE void ncm_model_params_update (NcmModel *model);
NCM_INLINE void ncm_model_orig_params_update (NcmModel *model);
NCM_INLINE NcmVector *ncm_model_orig_params_peek_vector (NcmModel *model);
void ncm_model_orig_params_log_all (NcmModel *model);

NCM_INLINE void ncm_model_param_set (NcmModel *model, guint n, gdouble val);
NCM_INLINE void ncm_model_param_set_default (NcmModel *model, guint n);
NCM_INLINE void ncm_model_orig_param_set (NcmModel *model, guint n, gdouble val);
NCM_INLINE void ncm_model_orig_vparam_set (NcmModel *model, guint n, guint i, gdouble val);
NCM_INLINE void ncm_model_orig_vparam_set_vector (NcmModel *model, guint n, NcmVector *val);
NCM_INLINE NcmSParam *ncm_model_param_peek_desc (NcmModel *model, guint n);
NCM_INLINE gdouble ncm_model_param_get (NcmModel *model, guint n);
NCM_INLINE gdouble ncm_model_orig_param_get (NcmModel *model, guint n);
NCM_INLINE gdouble ncm_model_orig_vparam_get (NcmModel *model, guint n, guint i);
NCM_INLINE NcmVector *ncm_model_orig_vparam_get_vector (NcmModel *model, guint n);

void ncm_model_params_copyto (NcmModel *model, NcmModel *model_dest);
void ncm_model_params_set_default (NcmModel *model);
void ncm_model_params_save_as_default (NcmModel *model);
void ncm_model_params_set_all (NcmModel *model, ...) G_GNUC_NULL_TERMINATED;
void ncm_model_params_set_all_data (NcmModel *model, gdouble *data);
void ncm_model_params_set_vector (NcmModel *model, NcmVector *v);
void ncm_model_params_set_model (NcmModel *model, NcmModel *model_src);
void ncm_model_params_print_all (NcmModel *model, FILE *out);
void ncm_model_params_log_all (NcmModel *model);
NcmVector *ncm_model_params_get_all (NcmModel *model);
gboolean ncm_model_params_valid (NcmModel *model);
gboolean ncm_model_params_valid_bounds (NcmModel *model);

gboolean ncm_model_orig_param_index_from_name (NcmModel *model, const gchar *param_name, guint *i);
gboolean ncm_model_param_index_from_name (NcmModel *model, const gchar *param_name, guint *i);
const gchar *ncm_model_orig_param_name (NcmModel *model, guint n);
const gchar *ncm_model_param_name (NcmModel *model, guint n);
const gchar *ncm_model_orig_param_symbol (NcmModel *model, guint n);
const gchar *ncm_model_param_symbol (NcmModel *model, guint n);

void ncm_model_param_set_by_name (NcmModel *model, const gchar *param_name, gdouble val);
void ncm_model_orig_param_set_by_name (NcmModel *model, const gchar *param_name, gdouble val);
gdouble ncm_model_param_get_by_name (NcmModel *model, const gchar *param_name);
gdouble ncm_model_orig_param_get_by_name (NcmModel *model, const gchar *param_name);

gdouble ncm_model_orig_param_get_scale (NcmModel *model, guint n);
gdouble ncm_model_orig_param_get_lower_bound (NcmModel *model, guint n);
gdouble ncm_model_orig_param_get_upper_bound (NcmModel *model, guint n);
gdouble ncm_model_orig_param_get_abstol (NcmModel *model, guint n);

gdouble ncm_model_param_get_scale (NcmModel *model, guint n);
gdouble ncm_model_param_get_lower_bound (NcmModel *model, guint n);
gdouble ncm_model_param_get_upper_bound (NcmModel *model, guint n);
gdouble ncm_model_param_get_abstol (NcmModel *model, guint n);

NcmParamType ncm_model_param_get_ftype (NcmModel *model, guint n);

void ncm_model_param_set_scale (NcmModel *model, guint n, const gdouble scale);
void ncm_model_param_set_lower_bound (NcmModel *model, guint n, const gdouble lb);
void ncm_model_param_set_upper_bound (NcmModel *model, guint n, const gdouble ub);
void ncm_model_param_set_abstol (NcmModel *model, guint n, const gdouble abstol);
void ncm_model_param_set_ftype (NcmModel *model, guint n, const NcmParamType ptype);

void ncm_model_reparam_df (NcmModel *model, NcmVector *fv, NcmVector *v);
void ncm_model_reparam_J (NcmModel *model, NcmMatrix *fJ, NcmMatrix *J);

gboolean ncm_model_is_submodel (NcmModel *model);
NcmModelID ncm_model_main_model (NcmModel *model);
void ncm_model_add_submodel (NcmModel *model, NcmModel *submodel);
guint ncm_model_get_submodel_len (NcmModel *model);
NcmModel *ncm_model_peek_submodel (NcmModel *model, guint i);
NcmModel *ncm_model_peek_submodel_by_mid (NcmModel *model, NcmModelID mid);
gint ncm_model_peek_submodel_pos_by_mid (NcmModel *model, NcmModelID mid);

gboolean ncm_model_type_is_submodel (GType model_type);
NcmModelID ncm_model_type_main_model (GType model_type);

#define NCM_MODEL_CLASS_IMPL_ALL ((guint64)(~((guint64)0)))
#define NCM_MODEL_OPT2IMPL(opt) (((guint64)1) << ((guint64)(opt)))
#define NCM_MODEL_2OPT2IMPL(opt1,opt2) (NCM_MODEL_OPT2IMPL (opt1) | NCM_MODEL_OPT2IMPL (opt2))
#define NCM_MODEL_3OPT2IMPL(opt1,opt2,opt3) (NCM_MODEL_2OPT2IMPL (opt1, opt2) | NCM_MODEL_OPT2IMPL (opt3))
#define NCM_MODEL_4OPT2IMPL(opt1,opt2,opt3,opt4) (NCM_MODEL_2OPT2IMPL (opt1, opt2) | NCM_MODEL_2OPT2IMPL (opt3, opt4))

/*
 * Model set functions
 */
#define NCM_MODEL_SET_IMPL_FUNC(NS_NAME,NsName,ns_name,type,name) \
void \
ns_name##_set_##name##_impl (NsName##Class *model_class, type f) \
{ \
  ncm_model_class_add_impl_opts (NCM_MODEL_CLASS (model_class), NS_NAME##_IMPL_##name, -1); \
  model_class->name = f; \
}

/*
 * Constant model functions call accessor
 */
#define NCM_MODEL_FUNC0_IMPL(NS_NAME,NsName,ns_name,name) \
NCM_INLINE gdouble ns_name##_##name (NsName *m) \
{ \
  return NS_NAME##_GET_CLASS (m)->name (NS_NAME (m)); \
}

/*
 * Model functions call
 */
#define NCM_MODEL_FUNC1_IMPL(NS_NAME,NsName,ns_name,name,var) \
NCM_INLINE gdouble ns_name##_##name (NsName *m, const gdouble var) \
{ \
  return NS_NAME##_GET_CLASS (m)->name (NS_NAME (m), var); \
}

/*
 * Model functions 2d call
 */
#define NCM_MODEL_FUNC2_IMPL(NS_NAME,NsName,ns_name,name) \
NCM_INLINE gdouble ns_name##_##name (NsName *m, const gdouble x, const gdouble y) \
{ \
  return NS_NAME##_GET_CLASS (m)->name (NS_NAME (m), x, y); \
}

/*
 * Constant model vector functions call accessor
 */
#define NCM_MODEL_VFUNC0_IMPL(NS_NAME,NsName,ns_name,name) \
NCM_INLINE gdouble ns_name##_##name (NsName *m, const guint n) \
{ \
  return NS_NAME##_GET_CLASS (m)->name (NS_NAME (m), n); \
}

/*
 * Model vector functions call
 */
#define NCM_MODEL_VFUNC1_IMPL(NS_NAME,NsName,ns_name,name,var) \
NCM_INLINE gdouble ns_name##_##name (NsName *m, const guint n, const gdouble var) \
{ \
  return NS_NAME##_GET_CLASS (m)->name (NS_NAME (m), n, var); \
}

/*
 * Model functions 2d call
 */
#define NCM_MODEL_VFUNC2_IMPL(NS_NAME,NsName,ns_name,name) \
NCM_INLINE gdouble ns_name##_##name (NsName *m, const guint n, const gdouble x, const gdouble y) \
{ \
  return NS_NAME##_GET_CLASS (m)->name (NS_NAME (m), n, x, y); \
}

G_END_DECLS

#endif /* _NCM_MODEL_H_ */

#ifndef _NCM_MODEL_INLINE_H_
#define _NCM_MODEL_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE NcmModelID
ncm_model_id (NcmModel *model)
{
  return NCM_MODEL_GET_CLASS (model)->model_id;
}

NCM_INLINE NcmModelID
ncm_model_id_by_type (GType model_type)
{
  if (!g_type_is_a (model_type, NCM_TYPE_MODEL))
  {
    g_error ("ncm_model_id_by_type: type (%s) is not a %s", g_type_name (model_type), g_type_name (NCM_TYPE_MODEL));
    return 0;
  }
  else
  {
    NcmModelClass *model_class = NCM_MODEL_CLASS (g_type_class_ref (model_type));
    NcmModelID id = model_class->model_id;
    g_type_class_unref (model_class);
    return id;
  }
}

NCM_INLINE gboolean
ncm_model_check_impl_flag (NcmModel *model, guint64 impl)
{
  if (impl == 0)
    return TRUE;
  else
    return ((NCM_MODEL_GET_CLASS (model)->impl_flag & impl) != 0);
}

NCM_INLINE gboolean 
ncm_model_check_impl_opt (NcmModel *model, gint opt)
{
  guint64 flag = NCM_MODEL_OPT2IMPL (opt);

  return ncm_model_check_impl_flag (model, flag);
}

NCM_INLINE guint
ncm_model_len (NcmModel *model)
{
  return model->total_len;
}

NCM_INLINE gboolean
ncm_model_state_is_update (NcmModel *model)
{
  return model->pkey == model->skey;
}

NCM_INLINE void
ncm_model_state_set_update (NcmModel *model)
{
  model->skey = model->pkey;
}

NCM_INLINE void 
ncm_model_state_mark_outdated (NcmModel *model)
{
  model->pkey++;
}

NCM_INLINE guint
ncm_model_sparam_len (NcmModel *model)
{
  return NCM_MODEL_GET_CLASS (model)->sparam_len;
}

NCM_INLINE guint
ncm_model_vparam_array_len (NcmModel *model)
{
  return NCM_MODEL_GET_CLASS (model)->vparam_len;
}

NCM_INLINE const gchar *
ncm_model_name (NcmModel *model)
{
  return NCM_MODEL_GET_CLASS (model)->name;
}

NCM_INLINE const gchar *
ncm_model_nick (NcmModel *model)
{
  return NCM_MODEL_GET_CLASS (model)->nick;
}

NCM_INLINE NcmReparam *
ncm_model_peek_reparam (NcmModel *model)
{
  return NCM_MODEL (model)->reparam;
}

NCM_INLINE gboolean
ncm_model_param_finite (NcmModel *model, guint i)
{
  NcmVector *params = model->reparam ? model->reparam->new_params : model->params;
  return gsl_finite (ncm_vector_get (params, i));
}

NCM_INLINE gboolean
ncm_model_params_finite (NcmModel *model)
{
  guint i;
  for (i = 0; i < ncm_model_len (model); i++)
  {
    if (!gsl_finite (ncm_vector_get (model->params, i)))
      return FALSE;
  }
  return TRUE;
}

NCM_INLINE void
ncm_model_params_update (NcmModel *model)
{
  model->pkey++;
  if (model->reparam)
    ncm_reparam_new2old (model->reparam, model);
}

NCM_INLINE void
ncm_model_orig_params_update (NcmModel *model)
{
  model->pkey++;
  if (model->reparam)
    ncm_reparam_old2new (model->reparam, model);
}

NCM_INLINE NcmVector *
ncm_model_orig_params_peek_vector (NcmModel *model)
{
  return model->params;
}

NCM_INLINE guint
ncm_model_vparam_index (NcmModel *model, guint n, guint i)
{
  return g_array_index (model->vparam_pos, guint, n) + i;
}

NCM_INLINE guint
ncm_model_vparam_len (NcmModel *model, guint n)
{
  return g_array_index (model->vparam_len, guint, n);
}

NCM_INLINE void
ncm_model_param_set (NcmModel *model, guint n, gdouble val)
{
  ncm_vector_set (model->p, n, val);
  ncm_model_params_update (model);
  return;
}

NCM_INLINE void
ncm_model_param_set_default (NcmModel *model, guint n)
{
  ncm_model_param_set (model, n, ncm_sparam_get_default_value (ncm_model_param_peek_desc (model, n)));
}

NCM_INLINE NcmSParam *
ncm_model_orig_param_peek_desc (NcmModel *model, guint n)
{
  g_assert_cmpuint (n, <, model->total_len);
  return g_ptr_array_index (model->sparams, n);
}

NCM_INLINE NcmSParam *
ncm_model_param_peek_desc (NcmModel *model, guint n)
{
  NcmReparam *reparam = ncm_model_peek_reparam (model);
  g_assert_cmpuint (n, <, model->total_len);
  if (reparam != NULL)
  {
    NcmSParam *sp = ncm_reparam_peek_param_desc (reparam, n);
    if (sp != NULL)
      return sp;
  }
  return ncm_model_orig_param_peek_desc (model, n);
}

NCM_INLINE gdouble
ncm_model_param_get (NcmModel *model, guint n)
{
  return ncm_vector_get (model->p, n);
}

NCM_INLINE void
ncm_model_orig_param_set (NcmModel *model, guint n, gdouble val)
{
  ncm_vector_set (model->params, n, val);
  ncm_model_orig_params_update (model);
  return;
}

NCM_INLINE gdouble
ncm_model_orig_param_get (NcmModel *model, guint n)
{
  return ncm_vector_get (model->params, n);
}

NCM_INLINE void
ncm_model_orig_vparam_set (NcmModel *model, guint n, guint i, gdouble val)
{
  ncm_vector_set (model->params, ncm_model_vparam_index (model, n, i), val);
  ncm_model_orig_params_update (model);
  return;
}

NCM_INLINE gdouble
ncm_model_orig_vparam_get (NcmModel *model, guint n, guint i)
{
  return ncm_vector_get (model->params, ncm_model_vparam_index (model, n, i));
}

NCM_INLINE void
ncm_model_orig_vparam_set_vector (NcmModel *model, guint n, NcmVector *val)
{
  ncm_vector_memcpy2 (model->params, val,
                      ncm_model_vparam_index (model, n, 0), 0,
                      ncm_model_vparam_len (model, n));
  ncm_model_orig_params_update (model);
}

NCM_INLINE NcmVector *
ncm_model_orig_vparam_get_vector (NcmModel *model, guint n)
{
  const guint vparam_len = ncm_model_vparam_len (model, n);
  if (vparam_len > 0)
  {
    NcmVector *val = ncm_vector_new (vparam_len);
    ncm_vector_memcpy2 (val, model->params,
                        0, ncm_model_vparam_index (model, n, 0),
                        ncm_model_vparam_len (model, n));
    return val;
  }
  else
    return NULL;
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_MODEL_INLINE_H_ */
