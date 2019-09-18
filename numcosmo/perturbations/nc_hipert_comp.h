/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_hipert_comp.h
 *
 *  Wed October 11 15:54:27 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_comp.h
 * Copyright (C) 2017 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_HIPERT_COMP_H_
#define _NC_HIPERT_COMP_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/perturbations/nc_hipert_bg_var.h>
#include <numcosmo/perturbations/nc_hipert_grav.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_COMP             (nc_hipert_comp_get_type ())
#define NC_HIPERT_COMP(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT_COMP, NcHIPertComp))
#define NC_HIPERT_COMP_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPERT_COMP, NcHIPertCompClass))
#define NC_IS_HIPERT_COMP(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT_COMP))
#define NC_IS_HIPERT_COMP_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPERT_COMP))
#define NC_HIPERT_COMP_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPERT_COMP, NcHIPertCompClass))

typedef struct _NcHIPertCompClass NcHIPertCompClass;
typedef struct _NcHIPertComp NcHIPertComp;
typedef struct _NcHIPertCompPrivate NcHIPertCompPrivate;

typedef guint (*NcHIPertCompNDynVar) (NcHIPertComp *comp);
typedef GArray *(*NcHIPertCompDeps) (NcHIPertComp *comp, guint vindex);

typedef void (*NcHIPertCompSetGauge) (NcHIPertComp *comp, NcHIPertGravGauge gauge);
typedef NcHIPertGravGauge (*NcHIPertCompGetGauge) (NcHIPertComp *comp);

typedef NcHIPertGravTScalarInfo *(*NcHIPertCompGetTScalarInfo) (NcHIPertComp *comp);

typedef void (*NcHIPertCompGetTScalar) (NcHIPertComp *comp, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar);
typedef void (*NcHIPertCompGetTVector) (NcHIPertComp *comp, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTVector *T_vector);
typedef void (*NcHIPertCompGetTTensor) (NcHIPertComp *comp, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTTensor *T_tensor);

typedef void (*NcHIPertCompGetDYScalar) (NcHIPertComp *comp, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar, NcHIPertGravScalar *G_scalar);

struct _NcHIPertCompClass
{
  /*< private >*/
  GObjectClass parent_class;
  NcHIPertCompNDynVar ndyn_var;
  NcHIPertCompDeps get_deps;
  NcHIPertCompSetGauge set_gauge;
  NcHIPertCompGetGauge get_gauge;
  NcHIPertCompGetTScalarInfo get_T_scalar_info;
  NcHIPertCompGetTScalar get_T_scalar;
  NcHIPertCompGetTVector get_T_vector;
  NcHIPertCompGetTTensor get_T_tensor;
  NcHIPertCompGetDYScalar get_dy_scalar;
};

struct _NcHIPertComp
{
  /*< private >*/
  GObject parent_instance;
  NcHIPertCompPrivate *priv;
};

GType nc_hipert_comp_get_type (void) G_GNUC_CONST;

NcHIPertComp *nc_hipert_comp_ref (NcHIPertComp *comp);
void nc_hipert_comp_free (NcHIPertComp *comp);
void nc_hipert_comp_clear (NcHIPertComp **comp);

NCM_INLINE NcHIPertBGVarID nc_hipert_comp_get_id (NcHIPertComp *comp);

guint nc_hipert_comp_ndyn_var (NcHIPertComp *comp);
GArray *nc_hipert_comp_get_deps (NcHIPertComp *comp, guint vindex);

void nc_hipert_comp_set_gauge (NcHIPertComp *comp, NcHIPertGravGauge gauge);
NcHIPertGravGauge nc_hipert_comp_get_gauge (NcHIPertComp *comp);

NcHIPertGravTScalarInfo *nc_hipert_comp_get_T_scalar_info (NcHIPertComp *comp);

NCM_INLINE void nc_hipert_comp_get_T_scalar (NcHIPertComp *comp, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar);
NCM_INLINE void nc_hipert_comp_get_dy_scalar (NcHIPertComp *comp, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar, NcHIPertGravScalar *G_scalar);

G_END_DECLS

#endif /* _NC_HIPERT_COMP_H_ */

#ifndef _NC_HIPERT_COMP_INLINE_H_
#define _NC_HIPERT_COMP_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE NcHIPertBGVarID 
nc_hipert_comp_get_id (NcHIPertComp *comp)
{
  const NcHIPertBGVarID id = nc_hipert_bg_var_class_get_id_by_ns (G_OBJECT_TYPE_NAME (comp));
  g_assert_cmpint (id, >=, 0);
  return id;
}

NCM_INLINE void 
nc_hipert_comp_get_T_scalar (NcHIPertComp *comp, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar)
{
  return NC_HIPERT_COMP_GET_CLASS (comp)->get_T_scalar (comp, bg_var, ydy, T_scalar);
}

NCM_INLINE void 
nc_hipert_comp_get_dy_scalar (NcHIPertComp *comp, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar, NcHIPertGravScalar *G_scalar)
{
  return NC_HIPERT_COMP_GET_CLASS (comp)->get_dy_scalar (comp, bg_var, ydy, T_scalar, G_scalar);
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_HIPERT_COMP_INLINE_H_ */
