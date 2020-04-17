/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_ode_eval.h
 *
 *  Thu December 13 11:13:17 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_ode_eval.h
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_ODE_EVAL_H_
#define _NCM_ODE_EVAL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>

G_BEGIN_DECLS

#define NCM_TYPE_ODE_EVAL             (ncm_ode_eval_get_type ())
#define NCM_ODE_EVAL(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_ODE_EVAL, NcmODEEval))
#define NCM_ODE_EVAL_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_ODE_EVAL, NcmODEEvalClass))
#define NCM_IS_ODE_EVAL(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_ODE_EVAL))
#define NCM_IS_ODE_EVAL_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_ODE_EVAL))
#define NCM_ODE_EVAL_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_ODE_EVAL, NcmODEEvalClass))

typedef struct _NcmODEEvalClass NcmODEEvalClass;
typedef struct _NcmODEEval NcmODEEval;
typedef struct _NcmODEEvalPrivate NcmODEEvalPrivate;

typedef gint (*NcmODEEvalF) (NcmODEEval *ode_eval, const guint sys_size, const gdouble t, const gdouble * restrict f, gdouble * restrict df);
typedef gint (*NcmODEEvalJDense) (NcmODEEval *ode_eval, const guint sys_size, const gdouble t, const gdouble * restrict f, gdouble ** restrict J_col);

struct _NcmODEEvalClass
{
  /*< private >*/
  GObjectClass parent_class;
  NcmODEEvalF df;
  NcmODEEvalJDense J_dense;
  void (*clean_ls) (NcmODEEval *ode_eval);
};

struct _NcmODEEval
{
  /*< private >*/
  GObject parent_instance;
  NcmODEEvalPrivate *priv;
};

/**
 * NcmODEEvalReturn:
 * @NCM_ODE_EVAL_RETURN_SUCCESS: computation done with success
 * 
 * NcmODEEval return codes.
 * 
 */
typedef enum _NcmODEEvalReturn
{
  NCM_ODE_EVAL_RETURN_SUCCESS = 0,
} NcmODEEvalReturn;

GType ncm_ode_eval_get_type (void) G_GNUC_CONST;

NcmODEEval *ncm_ode_eval_ref (NcmODEEval *ode_eval);

void ncm_ode_eval_free (NcmODEEval *ode_eval);
void ncm_ode_eval_clear (NcmODEEval **ode_eval);

NCM_INLINE gint ncm_ode_eval_df (NcmODEEval *ode_eval, const guint sys_size, const gdouble t, const gdouble * restrict f, gdouble * restrict df);
NCM_INLINE gint ncm_ode_eval_J_dense (NcmODEEval *ode_eval, const guint sys_size, const gdouble t, const gdouble * restrict f, gdouble ** restrict J_col);

#define NCM_ODE_EVAL_DECLARE_IMPL(ModuleObjName, module_obj_name, MODULE, OBJ_NAME, LocalStruct)         \
  GType module_obj_name##_get_type (void);                                                               \
  G_GNUC_BEGIN_IGNORE_DEPRECATIONS                                                                       \
  typedef struct _##ModuleObjName ModuleObjName;                                                         \
  typedef struct _##ModuleObjName##Class ModuleObjName##Class;                                           \
  struct _##ModuleObjName##Class { NcmODEEvalClass parent_class; };                                      \
                                                                                                         \
  NCM_INLINE ModuleObjName * MODULE##_##OBJ_NAME (gpointer ptr) {                                     \
    return G_TYPE_CHECK_INSTANCE_CAST (ptr, module_obj_name##_get_type (), ModuleObjName); }             \
  NCM_INLINE ModuleObjName##Class * MODULE##_##OBJ_NAME##_CLASS (gpointer ptr) {                      \
    return G_TYPE_CHECK_CLASS_CAST (ptr, module_obj_name##_get_type (), ModuleObjName##Class); }         \
  NCM_INLINE gboolean MODULE##_IS_##OBJ_NAME (gpointer ptr) {                                         \
    return G_TYPE_CHECK_INSTANCE_TYPE (ptr, module_obj_name##_get_type ()); }                            \
  NCM_INLINE gboolean MODULE##_IS_##OBJ_NAME##_CLASS (gpointer ptr) {                                 \
    return G_TYPE_CHECK_CLASS_TYPE (ptr, module_obj_name##_get_type ()); }                               \
  NCM_INLINE ModuleObjName##Class * MODULE##_##OBJ_NAME##_GET_CLASS (gpointer ptr) {                  \
    return G_TYPE_INSTANCE_GET_CLASS (ptr, module_obj_name##_get_type (), ModuleObjName##Class); }       \
  G_GNUC_END_IGNORE_DEPRECATIONS                                                                         \
  ModuleObjName *module_obj_name##_new (void);                                                           \
  LocalStruct *_##module_obj_name##_peek_ls (ModuleObjName *ode_eval1);

#define NCM_ODE_EVAL_DEFINE_IMPL(ModuleObjName, module_obj_name, MODULE, OBJ_NAME, LocalStruct, df0, J_dense0, clean_ls0) \
struct _##ModuleObjName { NcmODEEval parent_instance; LocalStruct ls; };                                                  \
G_DEFINE_TYPE (ModuleObjName, module_obj_name, NCM_TYPE_ODE_EVAL);                                                        \
static void module_obj_name##_init (ModuleObjName *ode_eval1) {                                                           \
  memset (&ode_eval1->ls, 0, sizeof (LocalStruct));                                                                       \
}                                                                                                                         \
static void module_obj_name##_class_init (ModuleObjName##Class *klass) {                                                  \
  NcmODEEvalClass *ode_eval_class = NCM_ODE_EVAL_CLASS (klass);                                                           \
  if ((df0) != NULL)       ode_eval_class->df       = (df0);                                                              \
  if ((J_dense0) != NULL)  ode_eval_class->J_dense  = (J_dense0);                                                         \
  if ((clean_ls0) != NULL) ode_eval_class->clean_ls = (clean_ls0);                                                        \
}                                                                                                                         \
ModuleObjName *module_obj_name##_new (void) {                                                                             \
  return g_object_new (module_obj_name##_get_type (), NULL);                                                              \
}                                                                                                                         \
LocalStruct *_##module_obj_name##_peek_ls (ModuleObjName *ode_eval1) {                                                    \
  return &ode_eval1->ls;                                                                                                  \
}

G_END_DECLS

#endif /* _NCM_ODE_EVAL_H_ */

#ifndef _NCM_ODE_EVAL_INLINE_H_
#define _NCM_ODE_EVAL_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE gint 
ncm_ode_eval_df (NcmODEEval *ode_eval, const guint sys_size, const gdouble t, const gdouble * restrict f, gdouble * restrict df)
{
  NCM_ODE_EVAL_GET_CLASS (ode_eval)->df (ode_eval, sys_size, t, f, df);
  return 0;
}

NCM_INLINE gint 
ncm_ode_eval_J_dense (NcmODEEval *ode_eval, const guint sys_size, const gdouble t, const gdouble * restrict f, gdouble ** restrict J_col)
{
  return NCM_ODE_EVAL_GET_CLASS (ode_eval)->J_dense (ode_eval, sys_size, t, f, J_col);
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_ODE_EVAL_INLINE_H_ */
