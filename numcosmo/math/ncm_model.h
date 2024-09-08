/***************************************************************************
 *            ncm_model.h
 *
 *  Fri February 24 21:18:21 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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

#define NCM_TYPE_MODEL (ncm_model_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmModel, ncm_model, NCM, MODEL, GObject)

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

#define NCM_MODEL_MAX_STATES (10)

typedef gdouble (*NcmModelFunc0) (NcmModel *model);
typedef gdouble (*NcmModelFunc1) (NcmModel *model, const gdouble x);
typedef gdouble (*NcmModelFunc2) (NcmModel *model, const gdouble x, const gdouble y);

typedef gdouble (*NcmModelVFunc0) (NcmModel *model, const guint n);
typedef gdouble (*NcmModelVFunc1) (NcmModel *model, const guint n, const gdouble x);
typedef gdouble (*NcmModelVFunc2) (NcmModel *model, const guint n, const gdouble x, const gdouble y);

/* Error domain for the NcmModel class */

/**
 * NcmModelError:
 * @NCM_MODEL_ERROR_PARAM_NAME_NOT_FOUND: The parameter name was not found.
 * @NCM_MODEL_ERROR_PARAM_ID_OUT_OF_RANGE: The parameter id is out of range.
 * @NCM_MODEL_ERROR_PARAM_INVALID_TYPE: The parameter type is invalid.
 * @NCM_MODEL_ERROR_ORIG_PARAM_NAME_NOT_FOUND: The original parameter name was not found.
 * @NCM_MODEL_ERROR_REPARAM_INCOMPATIBLE: The reparam is incompatible.
 * @NCM_MODEL_ERROR_INVALID_TYPE: The type is invalid.
 * @NCM_MODEL_ERROR_PARAM_CHANGED: The parameter was changed.
 *
 * Error codes returned by the #NcmModel class.
 *
 */
typedef enum _NcmModelError
{
  NCM_MODEL_ERROR_PARAM_NAME_NOT_FOUND,
  NCM_MODEL_ERROR_PARAM_ID_OUT_OF_RANGE,
  NCM_MODEL_ERROR_PARAM_INVALID_TYPE,
  NCM_MODEL_ERROR_ORIG_PARAM_NAME_NOT_FOUND,
  NCM_MODEL_ERROR_REPARAM_INCOMPATIBLE,
  NCM_MODEL_ERROR_INVALID_TYPE,
  NCM_MODEL_ERROR_PARAM_CHANGED
} NcmModelError;

GQuark ncm_model_error_quark (void);

#define NCM_MODEL_ERROR (ncm_model_error_quark ())

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
void ncm_model_set_reparam (NcmModel *model, NcmReparam *reparam, GError **error);
gboolean ncm_model_is_equal (NcmModel *model1, NcmModel *model2);

NcmModelID ncm_model_id (NcmModel *model);
NcmModelID ncm_model_id_by_type (GType model_type, GError **error);
gboolean ncm_model_check_impl_flag (NcmModel *model, guint64 impl);
gboolean ncm_model_check_impl_opt (NcmModel *model, gint opt);

gboolean ncm_model_check_impl_opts (NcmModel *model, gint opt1, ...);

guint ncm_model_len (NcmModel *model);
gboolean ncm_model_state_is_update (NcmModel *model);
void ncm_model_state_set_update (NcmModel *model);
gboolean ncm_model_lstate_is_update (NcmModel *model, guint i);
void ncm_model_lstate_set_update (NcmModel *model, guint i);
void ncm_model_state_mark_outdated (NcmModel *model);
guint64 ncm_model_state_get_pkey (NcmModel *model);

guint ncm_model_sparam_len (NcmModel *model);
guint ncm_model_vparam_array_len (NcmModel *model);
void ncm_model_set_vparam_len (NcmModel *model, guint n, guint len);
guint ncm_model_vparam_index (NcmModel *model, guint n, guint i);
guint ncm_model_vparam_len (NcmModel *model, guint n);
const gchar *ncm_model_name (NcmModel *model);
const gchar *ncm_model_nick (NcmModel *model);
NcmReparam *ncm_model_peek_reparam (NcmModel *model);
gboolean ncm_model_params_finite (NcmModel *model);
gboolean ncm_model_param_finite (NcmModel *model, guint i);
void ncm_model_params_update (NcmModel *model);
void ncm_model_orig_params_update (NcmModel *model);
NcmVector *ncm_model_orig_params_peek_vector (NcmModel *model);

void ncm_model_orig_params_log_all (NcmModel *model);

void ncm_model_param_set (NcmModel *model, guint n, gdouble val);
void ncm_model_param_set0 (NcmModel *model, guint n, gdouble val);
void ncm_model_param_set_default (NcmModel *model, guint n);
void ncm_model_orig_param_set (NcmModel *model, guint n, gdouble val);
void ncm_model_orig_vparam_set (NcmModel *model, guint n, guint i, gdouble val);
void ncm_model_orig_vparam_set_vector (NcmModel *model, guint n, NcmVector *val);
gdouble ncm_model_param_get (NcmModel *model, guint n);
gdouble ncm_model_orig_param_get (NcmModel *model, guint n);
gdouble ncm_model_orig_vparam_get (NcmModel *model, guint n, guint i);
NcmVector *ncm_model_orig_vparam_get_vector (NcmModel *model, guint n);

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
gboolean ncm_model_param_index_from_name (NcmModel *model, const gchar *param_name, guint *i, GError **error);
const gchar *ncm_model_orig_param_name (NcmModel *model, guint n);
const gchar *ncm_model_param_name (NcmModel *model, guint n);
const gchar *ncm_model_orig_param_symbol (NcmModel *model, guint n);
const gchar *ncm_model_param_symbol (NcmModel *model, guint n);
GPtrArray *ncm_model_param_names (NcmModel *model);

void ncm_model_param_set_by_name (NcmModel *model, const gchar *param_name, gdouble val, GError **error);
void ncm_model_orig_param_set_by_name (NcmModel *model, const gchar *param_name, gdouble val, GError **error);
gdouble ncm_model_param_get_by_name (NcmModel *model, const gchar *param_name, GError **error);
gdouble ncm_model_orig_param_get_by_name (NcmModel *model, const gchar *param_name, GError **error);

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
void ncm_model_params_set_default_ftype (NcmModel *model);

GHashTable *ncm_model_param_get_desc (NcmModel *model, gchar *param, GError **error);
void ncm_model_param_set_desc (NcmModel *model, gchar *param, GHashTable *desc, GError **error);

gboolean ncm_model_is_submodel (NcmModel *model);
NcmModelID ncm_model_main_model (NcmModel *model);
void ncm_model_add_submodel (NcmModel *model, NcmModel *submodel);
guint ncm_model_get_submodel_len (NcmModel *model);
NcmModel *ncm_model_peek_submodel (NcmModel *model, guint i);
NcmModel *ncm_model_peek_submodel_by_mid (NcmModel *model, NcmModelID mid);
gint ncm_model_peek_submodel_pos_by_mid (NcmModel *model, NcmModelID mid);

gboolean ncm_model_type_is_submodel (GType model_type);
NcmModelID ncm_model_type_main_model (GType model_type);

/* pygobject dict */
gdouble ncm_model___getitem__ (NcmModel *model, gchar *param, GError **error);
void ncm_model___setitem__ (NcmModel *model, gchar *param, gdouble val, GError **error);

#define NCM_MODEL_CLASS_IMPL_ALL ((guint64) (~((guint64) 0)))
#define NCM_MODEL_OPT2IMPL(opt) (((guint64) 1) << ((guint64) (opt)))
#define NCM_MODEL_2OPT2IMPL(opt1, opt2) (NCM_MODEL_OPT2IMPL (opt1) | NCM_MODEL_OPT2IMPL (opt2))
#define NCM_MODEL_3OPT2IMPL(opt1, opt2, opt3) (NCM_MODEL_2OPT2IMPL (opt1, opt2) | NCM_MODEL_OPT2IMPL (opt3))
#define NCM_MODEL_4OPT2IMPL(opt1, opt2, opt3, opt4) (NCM_MODEL_2OPT2IMPL (opt1, opt2) | NCM_MODEL_2OPT2IMPL (opt3, opt4))

/*
 * Model set functions
 */
#define NCM_MODEL_SET_IMPL_FUNC(NS_NAME, NsName, ns_name, type, name)                                   \
        void                                                                                            \
        ns_name ## _set_ ## name ## _impl (NsName ## Class * model_class, type f)                       \
        {                                                                                               \
          ncm_model_class_add_impl_opts (NCM_MODEL_CLASS (model_class), NS_NAME ## _IMPL_ ## name, -1); \
          model_class->name = f;                                                                        \
        }

/*
 * Constant model functions call accessor
 */
#define NCM_MODEL_FUNC0_IMPL(NS_NAME, NsName, ns_name, name)    \
        gdouble ns_name ## _ ## name (NsName * m)               \
        {                                                       \
          return NS_NAME ## _GET_CLASS (m)->name (NS_NAME (m)); \
        }

/*
 * Model functions call
 */
#define NCM_MODEL_FUNC1_IMPL(NS_NAME, NsName, ns_name, name, var)    \
        gdouble ns_name ## _ ## name (NsName * m, const gdouble var) \
        {                                                            \
          return NS_NAME ## _GET_CLASS (m)->name (NS_NAME (m), var); \
        }

/*
 * Model functions 2d call
 */
#define NCM_MODEL_FUNC2_IMPL(NS_NAME, NsName, ns_name, name)                        \
        gdouble ns_name ## _ ## name (NsName * m, const gdouble x, const gdouble y) \
        {                                                                           \
          return NS_NAME ## _GET_CLASS (m)->name (NS_NAME (m), x, y);               \
        }

/*
 * Constant model vector functions call accessor
 */
#define NCM_MODEL_VFUNC0_IMPL(NS_NAME, NsName, ns_name, name)      \
        gdouble ns_name ## _ ## name (NsName * m, const guint n)   \
        {                                                          \
          return NS_NAME ## _GET_CLASS (m)->name (NS_NAME (m), n); \
        }

/*
 * Model vector functions call
 */
#define NCM_MODEL_VFUNC1_IMPL(NS_NAME, NsName, ns_name, name, var)                  \
        gdouble ns_name ## _ ## name (NsName * m, const guint n, const gdouble var) \
        {                                                                           \
          return NS_NAME ## _GET_CLASS (m)->name (NS_NAME (m), n, var);             \
        }

/*
 * Model functions 2d call
 */
#define NCM_MODEL_VFUNC2_IMPL(NS_NAME, NsName, ns_name, name)                                      \
        gdouble ns_name ## _ ## name (NsName * m, const guint n, const gdouble x, const gdouble y) \
        {                                                                                          \
          return NS_NAME ## _GET_CLASS (m)->name (NS_NAME (m), n, x, y);                           \
        }

G_END_DECLS

#endif /* _NCM_MODEL_H_ */

