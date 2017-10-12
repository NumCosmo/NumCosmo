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
#include <numcosmo/nc_hicosmo.h>

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

typedef void (*NcHIPertCompCalcDY) (NcHIPertComp *comp, NcHICosmo *cosmo, gconstpointer *bg_vals, const NcmVector *y, NcmVector *dy);
typedef void (*NcHIPertCompCalcJ) (NcHIPertComp *comp, NcHICosmo *cosmo, gconstpointer *bg_vals, const NcmVector *y, NcmMatrix *J);
typedef void (*NcHIPertCompCalcDYJ) (NcHIPertComp *comp, NcHICosmo *cosmo, gconstpointer *bg_vals, const NcmVector *y, NcmVector *dy, NcmMatrix *J);
typedef void (*NcHIPertCompCalcTScalar) (NcHIPertComp *comp, NcHICosmo *cosmo, gconstpointer *bg_vals, const NcmVector *y, NcHIPertCompTScalar *TScalar);
typedef void (*NcHIPertCompCalcTVector) (NcHIPertComp *comp, NcHICosmo *cosmo, gconstpointer *bg_vals, const NcmVector *y, NcHIPertCompTScalar *TVector);
typedef void (*NcHIPertCompCalcTTensor) (NcHIPertComp *comp, NcHICosmo *cosmo, gconstpointer *bg_vals, const NcmVector *y, NcHIPertCompTScalar *TTensor);

struct _NcHIPertCompClass
{
  /*< private >*/
  GObjectClass parent_class;
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
GType nc_hipert_comp_get_type (void) G_GNUC_CONST;

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

G_END_DECLS

#endif /* _NC_HIPERT_COMP_H_ */

#ifndef _NC_HIPERT_COMP_INLINE_H_
#define _NC_HIPERT_COMP_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

NcHIPertCompTScalar *
nc_hipert_comp_T_scalar_new (void)
{
  return g_new0 (NcHIPertCompTScalar, 1);
}

NcHIPertCompTScalar *
nc_hipert_comp_T_scalar_dup (NcHIPertCompTScalar *Ts)
{
  NcHIPertCompTScalar *Ts_dup = g_new0 (NcHIPertCompTScalar, 1);

  Ts_dup[0] = Ts[0];

  return Ts_dup;
}

void
nc_hipert_comp_T_scalar_free (NcHIPertCompTScalar *Ts)
{
  g_free (Ts);
}

void
nc_hipert_comp_T_scalar_add (NcHIPertCompTScalar *Ts, const NcHIPertCompTScalar *Ts1, const NcHIPertCompTScalar *Ts2)
{
  Ts->drho_m_Aphi = Ts1->drho_m_Aphi + Ts2->drho_m_Aphi;
  Ts->A           = Ts1->A           + Ts2->A;
  Ts->rhopp_v     = Ts1->rhopp_v     + Ts2->rhopp_v;
  Ts->dp          = Ts1->dp          + Ts2->dp;
  Ts->Pi          = Ts1->Pi          + Ts2->Pi;
}

void
nc_hipert_comp_T_scalar_set_zero (NcHIPertCompTScalar *Ts)
{
  Ts->drho_m_Aphi = 0.0;
  Ts->A           = 0.0;
  Ts->rhopp_v     = 0.0;
  Ts->dp          = 0.0;
  Ts->Pi          = 0.0;
}

NcHIPertCompTVector *
nc_hipert_comp_T_vector_new (void)
{
  return g_new0 (NcHIPertCompTVector, 1);
}

NcHIPertCompTVector *
nc_hipert_comp_T_vector_dup (NcHIPertCompTVector *Tv)
{
  NcHIPertCompTVector *Tv_dup = g_new0 (NcHIPertCompTVector, 1);

  Tv_dup[0] = Tv[0];

  return Tv_dup;
}

void
nc_hipert_comp_T_vector_free (NcHIPertCompTVector *Tv)
{
  g_free (Tv);
}

void
nc_hipert_comp_T_vector_add (NcHIPertCompTVector *Tv, const NcHIPertCompTVector *Tv1, const NcHIPertCompTVector *Tv2)
{
  Tv->a = Tv1->a + Tv2->a;
}

void
nc_hipert_comp_T_vector_set_zero (NcHIPertCompTVector *Tv)
{
  Tv->a = 0.0;
}

NcHIPertCompTTensor *
nc_hipert_comp_T_tensor_new (void)
{
  return g_new0 (NcHIPertCompTTensor, 1);
}

NcHIPertCompTTensor *
nc_hipert_comp_T_tensor_dup (NcHIPertCompTTensor *Tt)
{
  NcHIPertCompTTensor *Tt_dup = g_new0 (NcHIPertCompTTensor, 1);

  Tt_dup[0] = Tt[0];

  return Tt_dup;
}

void
nc_hipert_comp_T_tensor_free (NcHIPertCompTTensor *Tt)
{
  g_free (Tt);
}

void
nc_hipert_comp_T_tensor_add (NcHIPertCompTTensor *Tt, const NcHIPertCompTTensor *Tt1, const NcHIPertCompTTensor *Tt2)
{
  Tt->a = Tt1->a + Tt2->a;
}

void
nc_hipert_comp_T_tensor_set_zero (NcHIPertCompTTensor *Tt)
{
  Tt->a = 0.0;
}

G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_HIPERT_COMP_INLINE_H_ */
