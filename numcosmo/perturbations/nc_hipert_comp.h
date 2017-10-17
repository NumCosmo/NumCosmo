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

typedef struct _NcHIPertCompTScalar
{
  gdouble drho_m_Aphi;
  gdouble A;
  gdouble rhopp_v;
  gdouble dp;
  gdouble Pi;
} NcHIPertCompTScalar;

typedef struct _NcHIPertCompTVector
{
  gdouble a;
} NcHIPertCompTVector;

typedef struct _NcHIPertCompTTensor
{
  gdouble a;
} NcHIPertCompTTensor;

typedef struct _NcHIPertCompCoupling
{
  gint c_a;
  gint c_b;
} NcHIPertCompCoupling;

typedef guint (*NcHIPertCompNDynVar) (NcHIPertComp *comp);
typedef GArray *(*NcHIPertCompCouplingGraph) (NcHIPertComp *comp);
typedef void (*NcHIPertCompCalcDY) (NcHIPertComp *comp, NcHIPertBGVar *bg_var, const NcmVector *y, NcmVector *dy);
typedef void (*NcHIPertCompCalcJ) (NcHIPertComp *comp, NcHIPertBGVar *bg_var, const NcmVector *y, NcmMatrix *J);
typedef void (*NcHIPertCompCalcDYJ) (NcHIPertComp *comp, NcHIPertBGVar *bg_var, const NcmVector *y, NcmVector *dy, NcmMatrix *J);
typedef void (*NcHIPertCompCalcTScalar) (NcHIPertComp *comp, NcHIPertBGVar *bg_var, const NcmVector *y, NcHIPertCompTScalar *TScalar);
typedef void (*NcHIPertCompCalcTVector) (NcHIPertComp *comp, NcHIPertBGVar *bg_var, const NcmVector *y, NcHIPertCompTScalar *TVector);
typedef void (*NcHIPertCompCalcTTensor) (NcHIPertComp *comp, NcHIPertBGVar *bg_var, const NcmVector *y, NcHIPertCompTScalar *TTensor);

struct _NcHIPertCompClass
{
  /*< private >*/
  GObjectClass parent_class;
  NcHIPertBGVarID bg_var_id;
  NcHIPertCompNDynVar ndyn_var;
  NcHIPertCompCouplingGraph cgraph;
  NcHIPertCompCalcDY dy;
  NcHIPertCompCalcJ J;
  NcHIPertCompCalcDYJ dy_J;
  NcHIPertCompCalcTScalar TScalar;
  NcHIPertCompCalcTVector TVector;
  NcHIPertCompCalcTTensor TTensor;
};

struct _NcHIPertComp
{
  /*< private >*/
  GObject parent_instance;
  NcHIPertCompPrivate *priv;
};

GType nc_hipert_comp_T_scalar_get_type (void) G_GNUC_CONST;
GType nc_hipert_comp_T_vector_get_type (void) G_GNUC_CONST;
GType nc_hipert_comp_T_tensor_get_type (void) G_GNUC_CONST;
GType nc_hipert_comp_coupling_get_type (void) G_GNUC_CONST;
GType nc_hipert_comp_get_type (void) G_GNUC_CONST;

void nc_hipert_comp_register_bg_var_id (NcHIPertCompClass *comp_class, guint cstruct_size, const gchar *ns, const gchar *desc, const gchar *long_desc);

/**
 * NC_HIPERT_COMP_BG_VAR_ID_FUNC: (skip)
 * @comp_ns: FIXME
 *
 * FIXME
 *
 */
#define NC_HIPERT_COMP_BG_VAR_ID_FUNC(comp_ns) comp_ns##_id

/**
 * NC_HIPERT_COMP_DECLARE_BG_VAR_ID: (skip)
 * @comp_ns: FIXME
 *
 * FIXME
 *
 */
#define NC_HIPERT_COMP_DECLARE_BG_VAR_ID(comp_ns) gint32 NC_HIPERT_COMP_BG_VAR_ID_FUNC(comp_ns) (void) G_GNUC_CONST

/**
 * NC_HIPERT_COMP_REGISTER_BG_VAR_ID: (skip)
 * @comp_ns: FIXME
 * @typemacro: FIXME
 *
 * FIXME
 *
 */
#define NC_HIPERT_COMP_REGISTER_BG_VAR_ID(comp_ns,typemacro) \
NcHIPertBGVarID NC_HIPERT_COMP_BG_VAR_ID_FUNC(comp_ns) (void) \
{ \
  static NcHIPertBGVarID id = -1; \
  if (id == -1) \
  { \
    NcHIPertCompClass *comp_class = g_type_class_ref (typemacro); \
    id = comp_class->bg_var_id; \
    g_type_class_unref (comp_class); \
  } \
  return id; \
}

G_INLINE_FUNC NcHIPertCompTScalar *nc_hipert_comp_T_scalar_new (void);
G_INLINE_FUNC NcHIPertCompTScalar *nc_hipert_comp_T_scalar_dup (NcHIPertCompTScalar *Ts);
G_INLINE_FUNC void nc_hipert_comp_T_scalar_free (NcHIPertCompTScalar *Ts);
G_INLINE_FUNC void nc_hipert_comp_T_scalar_add (NcHIPertCompTScalar *Ts, const NcHIPertCompTScalar *Ts1, const NcHIPertCompTScalar *Ts2);
G_INLINE_FUNC void nc_hipert_comp_T_scalar_set_zero (NcHIPertCompTScalar *Ts);

G_INLINE_FUNC NcHIPertCompTVector *nc_hipert_comp_T_vector_new (void);
G_INLINE_FUNC NcHIPertCompTVector *nc_hipert_comp_T_vector_dup (NcHIPertCompTVector *Tv);
G_INLINE_FUNC void nc_hipert_comp_T_vector_free (NcHIPertCompTVector *Tv);
G_INLINE_FUNC void nc_hipert_comp_T_vector_add (NcHIPertCompTVector *Tv, const NcHIPertCompTVector *Tv1, const NcHIPertCompTVector *Tv2);
G_INLINE_FUNC void nc_hipert_comp_T_vector_set_zero (NcHIPertCompTVector *Tv);

G_INLINE_FUNC NcHIPertCompTTensor *nc_hipert_comp_T_tensor_new (void);
G_INLINE_FUNC NcHIPertCompTTensor *nc_hipert_comp_T_tensor_dup (NcHIPertCompTTensor *Tt);
G_INLINE_FUNC void nc_hipert_comp_T_tensor_free (NcHIPertCompTTensor *Tt);
G_INLINE_FUNC void nc_hipert_comp_T_tensor_add (NcHIPertCompTTensor *Tt, const NcHIPertCompTTensor *Tt1, const NcHIPertCompTTensor *Tt2);
G_INLINE_FUNC void nc_hipert_comp_T_tensor_set_zero (NcHIPertCompTTensor *Tt);

G_INLINE_FUNC NcHIPertCompCoupling *nc_hipert_comp_coupling_new (void);
G_INLINE_FUNC NcHIPertCompCoupling *nc_hipert_comp_coupling_dup (NcHIPertCompCoupling *c);
G_INLINE_FUNC void nc_hipert_comp_coupling_free (NcHIPertCompCoupling *c);

NcHIPertComp *nc_hipert_comp_ref (NcHIPertComp *comp);
void nc_hipert_comp_free (NcHIPertComp *comp);
void nc_hipert_comp_clear (NcHIPertComp **comp);

G_INLINE_FUNC NcHIPertBGVarID nc_hipert_comp_get_id (NcHIPertComp *comp);

guint nc_hipert_comp_ndyn_var (NcHIPertComp *comp);
GArray *nc_hipert_comp_coupling_graph (NcHIPertComp *comp);

G_END_DECLS

#endif /* _NC_HIPERT_COMP_H_ */

#ifndef _NC_HIPERT_COMP_INLINE_H_
#define _NC_HIPERT_COMP_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

G_INLINE_FUNC NcHIPertCompTScalar *
nc_hipert_comp_T_scalar_new (void)
{
  return g_new0 (NcHIPertCompTScalar, 1);
}

G_INLINE_FUNC NcHIPertCompTScalar *
nc_hipert_comp_T_scalar_dup (NcHIPertCompTScalar *Ts)
{
  NcHIPertCompTScalar *Ts_dup = g_new0 (NcHIPertCompTScalar, 1);

  Ts_dup[0] = Ts[0];

  return Ts_dup;
}

G_INLINE_FUNC void
nc_hipert_comp_T_scalar_free (NcHIPertCompTScalar *Ts)
{
  g_free (Ts);
}

G_INLINE_FUNC void
nc_hipert_comp_T_scalar_add (NcHIPertCompTScalar *Ts, const NcHIPertCompTScalar *Ts1, const NcHIPertCompTScalar *Ts2)
{
  Ts->drho_m_Aphi = Ts1->drho_m_Aphi + Ts2->drho_m_Aphi;
  Ts->A           = Ts1->A           + Ts2->A;
  Ts->rhopp_v     = Ts1->rhopp_v     + Ts2->rhopp_v;
  Ts->dp          = Ts1->dp          + Ts2->dp;
  Ts->Pi          = Ts1->Pi          + Ts2->Pi;
}

G_INLINE_FUNC void
nc_hipert_comp_T_scalar_set_zero (NcHIPertCompTScalar *Ts)
{
  Ts->drho_m_Aphi = 0.0;
  Ts->A           = 0.0;
  Ts->rhopp_v     = 0.0;
  Ts->dp          = 0.0;
  Ts->Pi          = 0.0;
}

G_INLINE_FUNC NcHIPertCompTVector *
nc_hipert_comp_T_vector_new (void)
{
  return g_new0 (NcHIPertCompTVector, 1);
}

G_INLINE_FUNC NcHIPertCompTVector *
nc_hipert_comp_T_vector_dup (NcHIPertCompTVector *Tv)
{
  NcHIPertCompTVector *Tv_dup = g_new0 (NcHIPertCompTVector, 1);

  Tv_dup[0] = Tv[0];

  return Tv_dup;
}

G_INLINE_FUNC void
nc_hipert_comp_T_vector_free (NcHIPertCompTVector *Tv)
{
  g_free (Tv);
}

G_INLINE_FUNC void
nc_hipert_comp_T_vector_add (NcHIPertCompTVector *Tv, const NcHIPertCompTVector *Tv1, const NcHIPertCompTVector *Tv2)
{
  Tv->a = Tv1->a + Tv2->a;
}

G_INLINE_FUNC void
nc_hipert_comp_T_vector_set_zero (NcHIPertCompTVector *Tv)
{
  Tv->a = 0.0;
}

G_INLINE_FUNC NcHIPertCompTTensor *
nc_hipert_comp_T_tensor_new (void)
{
  return g_new0 (NcHIPertCompTTensor, 1);
}

G_INLINE_FUNC NcHIPertCompTTensor *
nc_hipert_comp_T_tensor_dup (NcHIPertCompTTensor *Tt)
{
  NcHIPertCompTTensor *Tt_dup = g_new0 (NcHIPertCompTTensor, 1);

  Tt_dup[0] = Tt[0];

  return Tt_dup;
}

G_INLINE_FUNC void
nc_hipert_comp_T_tensor_free (NcHIPertCompTTensor *Tt)
{
  g_free (Tt);
}

G_INLINE_FUNC void
nc_hipert_comp_T_tensor_add (NcHIPertCompTTensor *Tt, const NcHIPertCompTTensor *Tt1, const NcHIPertCompTTensor *Tt2)
{
  Tt->a = Tt1->a + Tt2->a;
}

G_INLINE_FUNC void
nc_hipert_comp_T_tensor_set_zero (NcHIPertCompTTensor *Tt)
{
  Tt->a = 0.0;
}

G_INLINE_FUNC NcHIPertCompCoupling *
nc_hipert_comp_coupling_new (void)
{
  return g_new0 (NcHIPertCompCoupling, 1);
}

G_INLINE_FUNC NcHIPertCompCoupling *
nc_hipert_comp_coupling_dup (NcHIPertCompCoupling *c)
{
  NcHIPertCompCoupling *c_dup = g_new0 (NcHIPertCompCoupling, 1);

  c_dup[0] = c[0];

  return c_dup;
}

G_INLINE_FUNC void
nc_hipert_comp_coupling_free (NcHIPertCompCoupling *c)
{
  g_free (c);
}

G_INLINE_FUNC NcHIPertBGVarID 
nc_hipert_comp_get_id (NcHIPertComp *comp)
{
  const NcHIPertBGVarID id = NC_HIPERT_COMP_GET_CLASS (comp)->bg_var_id;
  g_assert_cmpint (id, >=, 0);
  return id;
}

G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_HIPERT_COMP_INLINE_H_ */
