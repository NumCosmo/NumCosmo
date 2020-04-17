/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_hipert_grav.h
 *
 *  Thu October 12 14:32:08 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_grav.h
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

#ifndef _NC_HIPERT_GRAV_H_
#define _NC_HIPERT_GRAV_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/perturbations/nc_hipert_bg_var.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_GRAV             (nc_hipert_grav_get_type ())
#define NC_HIPERT_GRAV(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT_GRAV, NcHIPertGrav))
#define NC_HIPERT_GRAV_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPERT_GRAV, NcHIPertGravClass))
#define NC_IS_HIPERT_GRAV(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT_GRAV))
#define NC_IS_HIPERT_GRAV_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPERT_GRAV))
#define NC_HIPERT_GRAV_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPERT_GRAV, NcHIPertGravClass))

typedef struct _NcHIPertGravClass NcHIPertGravClass;
typedef struct _NcHIPertGrav NcHIPertGrav;
typedef struct _NcHIPertGravPrivate NcHIPertGravPrivate;
typedef struct _NcHIPertGravScalar NcHIPertGravScalar;
typedef struct _NcHIPertGravVector NcHIPertGravVector;
typedef struct _NcHIPertGravTensor NcHIPertGravTensor;
typedef struct _NcHIPertGravTScalar NcHIPertGravTScalar;
typedef struct _NcHIPertGravTVector NcHIPertGravTVector;
typedef struct _NcHIPertGravTTensor NcHIPertGravTTensor;
typedef struct _NcHIPertGravTScalarInfo NcHIPertGravTScalarInfo;
typedef struct _NcHIPertGravInfo NcHIPertGravInfo;

/**
 * NcHIPertGravScalar:
 * 
 * Boxed object describing scalar modes of the gravitational perturbations.
 */ 
struct _NcHIPertGravScalar
{
  /*< private >*/
  gdouble phi;
  gdouble dsigma;
  gdouble psi;
};

/**
 * NcHIPertGravVector:
 * 
 * Boxed object describing vector modes of the gravitational perturbations.
 */ 
struct _NcHIPertGravVector
{
  /*< private >*/
  gdouble dsigma[2];
};

/**
 * NcHIPertGravTensor:
 * 
 * Boxed object describing tensor modes of the gravitational perturbations.
 */ 
struct _NcHIPertGravTensor
{
  /*< private >*/
  gdouble h[2];
};

/**
 * NcHIPertGravTScalar:
 * 
 * Boxed object describing scalar modes of the energy-momentum tensor perturbations.
 */ 
struct _NcHIPertGravTScalar
{
  gdouble drho_m_Aphi;
  gdouble A;
  gdouble rhopp_v;
  gdouble dp;
  gdouble Pi;
};

/**
 * NcHIPertGravTVector:
 * 
 * Boxed object describing vector modes of the energy-momentum tensor perturbations.
 */ 
struct _NcHIPertGravTVector
{
  gdouble a;
};

/**
 * NcHIPertGravTTensor:
 * 
 * Boxed object describing tensor modes of the energy-momentum tensor perturbations.
 */ 
struct _NcHIPertGravTTensor
{
  gdouble a;
};

/**
 * NcHIPertGravTScalarInfo:
 * 
 * Boxed object describing the dependencies of the scalar modes of the gravitational perturbations.
 */ 
struct _NcHIPertGravTScalarInfo
{
  GArray *drho_deps;
  GArray *rhoppv_deps;
  GArray *dp_deps;
  GArray *dPi_deps;
};

/**
 * NcHIPertGravSElem:
 * @NC_HIPERT_GRAV_SELEM_PHI: gravitation potential $\phi$
 * @NC_HIPERT_GRAV_SELEM_DSIGMA: gravitation potential $\delta\sigma$
 * @NC_HIPERT_GRAV_SELEM_PSI: gravitation potential $\psi$
 * @NC_HIPERT_GRAV_SELEM_DOTPSI: time derivative of the gravitation potential $\dot{\psi}$
 * @NC_HIPERT_GRAV_SELEM_DRHO: energy momentum tensor component $\delta\rho$
 * @NC_HIPERT_GRAV_SELEM_RHOPPV: energy momentum tensor component $(\rho+p)\mathcal{V}$
 * @NC_HIPERT_GRAV_SELEM_DP: energy momentum tensor component $\delta{}p$
 * @NC_HIPERT_GRAV_SELEM_DPI: energy momentum tensor component $\delta\Pi$
 * 
 * Elements present in the scalar equations.
 * 
 */
typedef enum /*< enum,underscore_name=NC_HIPERT_GRAV_SELEM >*/
{
  NC_HIPERT_GRAV_SELEM_PHI = G_MININT,
  NC_HIPERT_GRAV_SELEM_DSIGMA,
  NC_HIPERT_GRAV_SELEM_PSI,
  NC_HIPERT_GRAV_SELEM_DOTPSI,
  NC_HIPERT_GRAV_SELEM_DRHO,
  NC_HIPERT_GRAV_SELEM_RHOPPV,
  NC_HIPERT_GRAV_SELEM_DP,
  NC_HIPERT_GRAV_SELEM_DPI,
  /* < private > */
  NC_HIPERT_GRAV_SELEM_LEN = NC_HIPERT_GRAV_SELEM_DPI - NC_HIPERT_GRAV_SELEM_PHI + 1, /*< skip >*/
} NcHIPertGravSElem;

/**
 * NcHIPertGravGauge:
 * @NC_HIPERT_GRAV_GAUGE_SYNCHRONOUS: Synchronous gauge
 * @NC_HIPERT_GRAV_GAUGE_NEWTONIAN: Newtonian gauge
 * @NC_HIPERT_GRAV_GAUGE_CONST_CURV: Constant curvature gauge
 * @NC_HIPERT_GRAV_GAUGE_CONST_EXP: Constant expansion gauge
 * 
 * Gravitation gauges.
 * 
 */
typedef enum /*< enum,underscore_name=NC_HIPERT_GRAV_GAUGE >*/
{
  NC_HIPERT_GRAV_GAUGE_SYNCHRONOUS,
  NC_HIPERT_GRAV_GAUGE_NEWTONIAN,
  NC_HIPERT_GRAV_GAUGE_CONST_CURV,
  NC_HIPERT_GRAV_GAUGE_CONST_EXP,
  /* < private > */
  NC_HIPERT_GRAV_GAUGE_LEN, /*< skip >*/
} NcHIPertGravGauge;

#define NC_HIPERT_GRAV_DYN_VAR(n) (n)

/**
 * NcHIPertGravInfo:
 * @phi_deps: (array) (element-type gint): $\phi$ dependencies
 * @dsigma_deps: (array) (element-type gint): $\delta\sigma$ dependencies
 * @psi_deps: (array) (element-type gint): $\psi$ dependencies
 * @dotpsi_deps: (array) (element-type gint): $\psi$ dependencies
 * 
 * Gravitation section info
 * 
 */
struct _NcHIPertGravInfo
{
  GArray *phi_deps;
  GArray *dsigma_deps;
  GArray *psi_deps;
  GArray *dotpsi_deps; 
};

typedef guint (*NcHIPertGravNDynVar) (NcHIPertGrav *grav);
typedef GArray *(*NcHIPertGravDeps) (NcHIPertGrav *grav, guint vindex);

typedef void (*NcHIPertGravSetGauge) (NcHIPertGrav *grav, NcHIPertGravGauge gauge);
typedef NcHIPertGravGauge (*NcHIPertGravGetGauge) (NcHIPertGrav *grav);

typedef NcHIPertGravInfo *(*NcHIPertGravGetScalarInfo) (NcHIPertGrav *grav);
typedef void (*NcHIPertGravGetScalar) (NcHIPertGrav *grav, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar, NcHIPertGravScalar *G_scalar);

typedef void (*NcHIPertGravGetDYScalar) (NcHIPertGrav *grav, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar, NcHIPertGravScalar *G_scalar);

struct _NcHIPertGravClass
{
  /*< private >*/
  GObjectClass parent_class;
  NcHIPertGravNDynVar ndyn_var;
  NcHIPertGravDeps get_deps;
  NcHIPertGravSetGauge set_gauge;
  NcHIPertGravGetGauge get_gauge;
  NcHIPertGravGetScalarInfo get_G_scalar_info;
  NcHIPertGravGetScalar get_G_scalar;
  NcHIPertGravGetDYScalar get_dy_scalar;
};

struct _NcHIPertGrav
{
  /*< private >*/
  GObject parent_instance;
  NcHIPertGravPrivate *priv;
};

GType nc_hipert_grav_scalar_get_type (void) G_GNUC_CONST;
GType nc_hipert_grav_vector_get_type (void) G_GNUC_CONST;
GType nc_hipert_grav_tensor_get_type (void) G_GNUC_CONST;

GType nc_hipert_grav_T_scalar_get_type (void) G_GNUC_CONST;
GType nc_hipert_grav_T_vector_get_type (void) G_GNUC_CONST;
GType nc_hipert_grav_T_tensor_get_type (void) G_GNUC_CONST;

GType nc_hipert_grav_info_get_type (void) G_GNUC_CONST;
GType nc_hipert_grav_T_scalar_info_get_type (void) G_GNUC_CONST;

GType nc_hipert_grav_get_type (void) G_GNUC_CONST;

NCM_INLINE NcHIPertGravScalar *nc_hipert_grav_scalar_new (void);
NCM_INLINE NcHIPertGravScalar *nc_hipert_grav_scalar_dup (NcHIPertGravScalar *gs);
NCM_INLINE void nc_hipert_grav_scalar_free (NcHIPertGravScalar *gs);
NCM_INLINE void nc_hipert_grav_scalar_set_zero (NcHIPertGravScalar *gs);

NCM_INLINE NcHIPertGravVector *nc_hipert_grav_vector_new (void);
NCM_INLINE NcHIPertGravVector *nc_hipert_grav_vector_dup (NcHIPertGravVector *gv);
NCM_INLINE void nc_hipert_grav_vector_free (NcHIPertGravVector *gv);
NCM_INLINE void nc_hipert_grav_vector_set_zero (NcHIPertGravVector *gv);

NCM_INLINE NcHIPertGravTensor *nc_hipert_grav_tensor_new (void);
NCM_INLINE NcHIPertGravTensor *nc_hipert_grav_tensor_dup (NcHIPertGravTensor *gt);
NCM_INLINE void nc_hipert_grav_tensor_free (NcHIPertGravTensor *gt);
NCM_INLINE void nc_hipert_grav_tensor_set_zero (NcHIPertGravTensor *gt);

NCM_INLINE NcHIPertGravTScalar *nc_hipert_grav_T_scalar_new (void);
NCM_INLINE NcHIPertGravTScalar *nc_hipert_grav_T_scalar_dup (NcHIPertGravTScalar *Ts);
NCM_INLINE void nc_hipert_grav_T_scalar_free (NcHIPertGravTScalar *Ts);
NCM_INLINE void nc_hipert_grav_T_scalar_add (NcHIPertGravTScalar *Ts, const NcHIPertGravTScalar *Ts1, const NcHIPertGravTScalar *Ts2);
NCM_INLINE void nc_hipert_grav_T_scalar_set_zero (NcHIPertGravTScalar *Ts);

NCM_INLINE NcHIPertGravTVector *nc_hipert_grav_T_vector_new (void);
NCM_INLINE NcHIPertGravTVector *nc_hipert_grav_T_vector_dup (NcHIPertGravTVector *Tv);
NCM_INLINE void nc_hipert_grav_T_vector_free (NcHIPertGravTVector *Tv);
NCM_INLINE void nc_hipert_grav_T_vector_add (NcHIPertGravTVector *Tv, const NcHIPertGravTVector *Tv1, const NcHIPertGravTVector *Tv2);
NCM_INLINE void nc_hipert_grav_T_vector_set_zero (NcHIPertGravTVector *Tv);

NCM_INLINE NcHIPertGravTTensor *nc_hipert_grav_T_tensor_new (void);
NCM_INLINE NcHIPertGravTTensor *nc_hipert_grav_T_tensor_dup (NcHIPertGravTTensor *Tt);
NCM_INLINE void nc_hipert_grav_T_tensor_free (NcHIPertGravTTensor *Tt);
NCM_INLINE void nc_hipert_grav_T_tensor_add (NcHIPertGravTTensor *Tt, const NcHIPertGravTTensor *Tt1, const NcHIPertGravTTensor *Tt2);
NCM_INLINE void nc_hipert_grav_T_tensor_set_zero (NcHIPertGravTTensor *Tt);

NCM_INLINE NcHIPertGravInfo *nc_hipert_grav_info_new (void);
NCM_INLINE NcHIPertGravInfo *nc_hipert_grav_info_dup (NcHIPertGravInfo *ginfo);

NCM_INLINE void nc_hipert_grav_info_free (NcHIPertGravInfo *ginfo);
NCM_INLINE void nc_hipert_grav_info_set_zero (NcHIPertGravInfo *ginfo);
NCM_INLINE void nc_hipert_grav_info_set_phi_deps (NcHIPertGravInfo *ginfo, GArray *phi_deps); 
NCM_INLINE void nc_hipert_grav_info_set_dsigma_deps (NcHIPertGravInfo *ginfo, GArray *dsigma_deps); 
NCM_INLINE void nc_hipert_grav_info_set_psi_deps (NcHIPertGravInfo *ginfo, GArray *psi_deps);
NCM_INLINE void nc_hipert_grav_info_set_dotpsi_deps (NcHIPertGravInfo *ginfo, GArray *dotpsi_deps);

NCM_INLINE GArray *nc_hipert_grav_info_get_phi_deps (NcHIPertGravInfo *ginfo);
NCM_INLINE GArray *nc_hipert_grav_info_get_dsigma_deps (NcHIPertGravInfo *ginfo); 
NCM_INLINE GArray *nc_hipert_grav_info_get_psi_deps (NcHIPertGravInfo *ginfo);
NCM_INLINE GArray *nc_hipert_grav_info_get_dotpsi_deps (NcHIPertGravInfo *ginfo); 

NCM_INLINE NcHIPertGravTScalarInfo *nc_hipert_grav_T_scalar_info_new (void);
NCM_INLINE NcHIPertGravTScalarInfo *nc_hipert_grav_T_scalar_info_dup (NcHIPertGravTScalarInfo *Tsinfo);
NCM_INLINE void nc_hipert_grav_T_scalar_info_free (NcHIPertGravTScalarInfo *Tsinfo);
NCM_INLINE void nc_hipert_grav_T_scalar_info_append (NcHIPertGravTScalarInfo *Tsinfo, NcHIPertGravTScalarInfo *Tsinfo1);
NCM_INLINE void nc_hipert_grav_T_scalar_info_set_zero (NcHIPertGravTScalarInfo *Tsinfo);

NcHIPertGrav *nc_hipert_grav_ref (NcHIPertGrav *grav);
void nc_hipert_grav_free (NcHIPertGrav *grav);
void nc_hipert_grav_clear (NcHIPertGrav **grav);

NCM_INLINE NcHIPertBGVarID nc_hipert_grav_get_id (NcHIPertGrav *grav);

NCM_INLINE guint nc_hipert_grav_ndyn_var (NcHIPertGrav *grav);
NCM_INLINE GArray *nc_hipert_grav_get_deps (NcHIPertGrav *grav, guint vindex);

NCM_INLINE void nc_hipert_grav_set_gauge (NcHIPertGrav *grav, NcHIPertGravGauge gauge);
NCM_INLINE NcHIPertGravGauge nc_hipert_grav_get_gauge (NcHIPertGrav *grav);

NCM_INLINE NcHIPertGravInfo *nc_hipert_grav_get_G_scalar_info (NcHIPertGrav *grav);
NCM_INLINE void nc_hipert_grav_get_G_scalar (NcHIPertGrav *grav, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar, NcHIPertGravScalar *G_scalar);

NCM_INLINE void nc_hipert_grav_get_dy_scalar (NcHIPertGrav *grav, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar, NcHIPertGravScalar *G_scalar);

G_END_DECLS

#endif /* _NC_HIPERT_GRAV_H_ */

#ifndef _NC_HIPERT_GRAV_INLINE_H_
#define _NC_HIPERT_GRAV_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE NcHIPertGravScalar *
nc_hipert_grav_scalar_new (void)
{
  return g_new0 (NcHIPertGravScalar, 1);
}

NCM_INLINE NcHIPertGravScalar *
nc_hipert_grav_scalar_dup (NcHIPertGravScalar *gs)
{
  NcHIPertGravScalar *gs_dup = g_new0 (NcHIPertGravScalar, 1);

  gs_dup[0] = gs[0];

  return gs_dup;
}

NCM_INLINE void
nc_hipert_grav_scalar_free (NcHIPertGravScalar *gs)
{
  g_free (gs);
}

NCM_INLINE void
nc_hipert_grav_scalar_set_zero (NcHIPertGravScalar *gs)
{
  gs->phi    = 0.0;
  gs->dsigma = 0.0;
  gs->psi    = 0.0;
}

NCM_INLINE NcHIPertGravVector *
nc_hipert_grav_vector_new (void)
{
  return g_new0 (NcHIPertGravVector, 1);
}

NCM_INLINE NcHIPertGravVector *
nc_hipert_grav_vector_dup (NcHIPertGravVector *gv)
{
  NcHIPertGravVector *gv_dup = g_new0 (NcHIPertGravVector, 1);

  gv_dup[0] = gv[0];

  return gv_dup;
}

NCM_INLINE void
nc_hipert_grav_vector_free (NcHIPertGravVector *gv)
{
  g_free (gv);
}

NCM_INLINE void
nc_hipert_grav_vector_set_zero (NcHIPertGravVector *gv)
{
  gv->dsigma[0] = 0.0;
  gv->dsigma[1] = 0.0;
}

NCM_INLINE NcHIPertGravTensor *
nc_hipert_grav_tensor_new (void)
{
  return g_new0 (NcHIPertGravTensor, 1);
}

NCM_INLINE NcHIPertGravTensor *
nc_hipert_grav_tensor_dup (NcHIPertGravTensor *gt)
{
  NcHIPertGravTensor *gt_dup = g_new0 (NcHIPertGravTensor, 1);

  gt_dup[0] = gt[0];

  return gt_dup;
}

NCM_INLINE void
nc_hipert_grav_tensor_free (NcHIPertGravTensor *gt)
{
  g_free (gt);
}

NCM_INLINE void
nc_hipert_grav_tensor_set_zero (NcHIPertGravTensor *gt)
{
  gt->h[0] = 0.0;
  gt->h[1] = 0.0;
}

NCM_INLINE NcHIPertGravTScalar *
nc_hipert_grav_T_scalar_new (void)
{
  return g_new0 (NcHIPertGravTScalar, 1);
}

NCM_INLINE NcHIPertGravTScalar *
nc_hipert_grav_T_scalar_dup (NcHIPertGravTScalar *Ts)
{
  NcHIPertGravTScalar *Ts_dup = g_new0 (NcHIPertGravTScalar, 1);

  Ts_dup[0] = Ts[0];

  return Ts_dup;
}

NCM_INLINE void
nc_hipert_grav_T_scalar_free (NcHIPertGravTScalar *Ts)
{
  g_free (Ts);
}

NCM_INLINE void
nc_hipert_grav_T_scalar_add (NcHIPertGravTScalar *Ts, const NcHIPertGravTScalar *Ts1, const NcHIPertGravTScalar *Ts2)
{
  Ts->drho_m_Aphi = Ts1->drho_m_Aphi + Ts2->drho_m_Aphi;
  Ts->A           = Ts1->A           + Ts2->A;
  Ts->rhopp_v     = Ts1->rhopp_v     + Ts2->rhopp_v;
  Ts->dp          = Ts1->dp          + Ts2->dp;
  Ts->Pi          = Ts1->Pi          + Ts2->Pi;
}

NCM_INLINE void
nc_hipert_grav_T_scalar_set_zero (NcHIPertGravTScalar *Ts)
{
  Ts->drho_m_Aphi = 0.0;
  Ts->A           = 0.0;
  Ts->rhopp_v     = 0.0;
  Ts->dp          = 0.0;
  Ts->Pi          = 0.0;
}

NCM_INLINE NcHIPertGravTVector *
nc_hipert_grav_T_vector_new (void)
{
  return g_new0 (NcHIPertGravTVector, 1);
}

NCM_INLINE NcHIPertGravTVector *
nc_hipert_grav_T_vector_dup (NcHIPertGravTVector *Tv)
{
  NcHIPertGravTVector *Tv_dup = g_new0 (NcHIPertGravTVector, 1);

  Tv_dup[0] = Tv[0];

  return Tv_dup;
}

NCM_INLINE void
nc_hipert_grav_T_vector_free (NcHIPertGravTVector *Tv)
{
  g_free (Tv);
}

NCM_INLINE void
nc_hipert_grav_T_vector_add (NcHIPertGravTVector *Tv, const NcHIPertGravTVector *Tv1, const NcHIPertGravTVector *Tv2)
{
  Tv->a = Tv1->a + Tv2->a;
}

NCM_INLINE void
nc_hipert_grav_T_vector_set_zero (NcHIPertGravTVector *Tv)
{
  Tv->a = 0.0;
}

NCM_INLINE NcHIPertGravTTensor *
nc_hipert_grav_T_tensor_new (void)
{
  return g_new0 (NcHIPertGravTTensor, 1);
}

NCM_INLINE NcHIPertGravTTensor *
nc_hipert_grav_T_tensor_dup (NcHIPertGravTTensor *Tt)
{
  NcHIPertGravTTensor *Tt_dup = g_new0 (NcHIPertGravTTensor, 1);

  Tt_dup[0] = Tt[0];

  return Tt_dup;
}

NCM_INLINE void
nc_hipert_grav_T_tensor_free (NcHIPertGravTTensor *Tt)
{
  g_free (Tt);
}

NCM_INLINE void
nc_hipert_grav_T_tensor_add (NcHIPertGravTTensor *Tt, const NcHIPertGravTTensor *Tt1, const NcHIPertGravTTensor *Tt2)
{
  Tt->a = Tt1->a + Tt2->a;
}

NCM_INLINE void
nc_hipert_grav_T_tensor_set_zero (NcHIPertGravTTensor *Tt)
{
  Tt->a = 0.0;
}

NCM_INLINE NcHIPertGravInfo *
nc_hipert_grav_info_new (void)
{
  NcHIPertGravInfo *ginfo = g_new0 (NcHIPertGravInfo, 1);
  
  ginfo->phi_deps     = g_array_new (TRUE, TRUE, sizeof (gint));
  ginfo->dsigma_deps  = g_array_new (TRUE, TRUE, sizeof (gint));
  ginfo->psi_deps     = g_array_new (TRUE, TRUE, sizeof (gint));
  ginfo->dotpsi_deps  = g_array_new (TRUE, TRUE, sizeof (gint));

  return ginfo;
}

NCM_INLINE NcHIPertGravInfo *
nc_hipert_grav_info_dup (NcHIPertGravInfo *ginfo)
{
  NcHIPertGravInfo *ginfo_dup = nc_hipert_grav_info_new ();

  g_array_append_vals (ginfo_dup->phi_deps,    ginfo->phi_deps->data,    ginfo->phi_deps->len);
  g_array_append_vals (ginfo_dup->dsigma_deps, ginfo->dsigma_deps->data, ginfo->dsigma_deps->len);
  g_array_append_vals (ginfo_dup->psi_deps,    ginfo->psi_deps->data,    ginfo->psi_deps->len);

  g_array_append_vals (ginfo_dup->dotpsi_deps, ginfo->dotpsi_deps->data, ginfo->dotpsi_deps->len);

  return ginfo_dup;
}

NCM_INLINE void
nc_hipert_grav_info_free (NcHIPertGravInfo *ginfo)
{
  g_array_unref (ginfo->phi_deps);
  g_array_unref (ginfo->dsigma_deps);
  g_array_unref (ginfo->psi_deps);
  g_array_unref (ginfo->dotpsi_deps);
  g_free (ginfo);
}

NCM_INLINE void
nc_hipert_grav_info_set_zero (NcHIPertGravInfo *ginfo)
{
  g_array_set_size (ginfo->phi_deps,    0);
  g_array_set_size (ginfo->dsigma_deps, 0);
  g_array_set_size (ginfo->psi_deps,    0);
  g_array_set_size (ginfo->dotpsi_deps, 0);
}

NCM_INLINE void
nc_hipert_grav_info_set_phi_deps (NcHIPertGravInfo *ginfo, GArray *phi_deps) 
{
  g_assert_cmpuint (g_array_get_element_size (ginfo->phi_deps), ==, g_array_get_element_size (phi_deps));
  g_array_set_size (ginfo->phi_deps,    0);
  g_array_append_vals (ginfo->phi_deps, phi_deps->data, phi_deps->len);
}

NCM_INLINE void
nc_hipert_grav_info_set_dsigma_deps (NcHIPertGravInfo *ginfo, GArray *dsigma_deps) 
{
  g_assert_cmpuint (g_array_get_element_size (ginfo->dsigma_deps), ==, g_array_get_element_size (dsigma_deps));
  g_array_set_size (ginfo->dsigma_deps,    0);
  g_array_append_vals (ginfo->dsigma_deps, dsigma_deps->data, dsigma_deps->len);
}

NCM_INLINE void
nc_hipert_grav_info_set_psi_deps (NcHIPertGravInfo *ginfo, GArray *psi_deps) 
{
  g_assert_cmpuint (g_array_get_element_size (ginfo->psi_deps), ==, g_array_get_element_size (psi_deps));
  g_array_set_size (ginfo->psi_deps,    0);
  g_array_append_vals (ginfo->psi_deps, psi_deps->data, psi_deps->len);
}

NCM_INLINE void
nc_hipert_grav_info_set_dotpsi_deps (NcHIPertGravInfo *ginfo, GArray *dotpsi_deps) 
{
  g_assert_cmpuint (g_array_get_element_size (ginfo->dotpsi_deps), ==, g_array_get_element_size (dotpsi_deps));
  g_array_set_size (ginfo->dotpsi_deps,    0);
  g_array_append_vals (ginfo->dotpsi_deps, dotpsi_deps->data, dotpsi_deps->len);
}

NCM_INLINE GArray *
nc_hipert_grav_info_get_phi_deps (NcHIPertGravInfo *ginfo) 
{
  return g_array_ref (ginfo->phi_deps);
}

NCM_INLINE GArray *
nc_hipert_grav_info_get_dsigma_deps (NcHIPertGravInfo *ginfo) 
{
  return g_array_ref (ginfo->dsigma_deps);
}

NCM_INLINE GArray *
nc_hipert_grav_info_get_psi_deps (NcHIPertGravInfo *ginfo) 
{
  return g_array_ref (ginfo->psi_deps);
}

NCM_INLINE GArray *
nc_hipert_grav_info_get_dotpsi_deps (NcHIPertGravInfo *ginfo) 
{
  return g_array_ref (ginfo->dotpsi_deps);
}

NCM_INLINE NcHIPertGravTScalarInfo *
nc_hipert_grav_T_scalar_info_new (void)
{
  NcHIPertGravTScalarInfo *Tsinfo = g_new0 (NcHIPertGravTScalarInfo, 1);

  Tsinfo->drho_deps   = g_array_new (TRUE, TRUE, sizeof (gint));
  Tsinfo->rhoppv_deps = g_array_new (TRUE, TRUE, sizeof (gint));
  Tsinfo->dp_deps     = g_array_new (TRUE, TRUE, sizeof (gint));
  Tsinfo->dPi_deps    = g_array_new (TRUE, TRUE, sizeof (gint));

  return Tsinfo;
}

NCM_INLINE NcHIPertGravTScalarInfo *
nc_hipert_grav_T_scalar_info_dup (NcHIPertGravTScalarInfo *Tsinfo)
{
  NcHIPertGravTScalarInfo *Tsinfo_dup = nc_hipert_grav_T_scalar_info_new ();
  nc_hipert_grav_T_scalar_info_append (Tsinfo_dup, Tsinfo);
  return Tsinfo_dup;
}

NCM_INLINE void 
nc_hipert_grav_T_scalar_info_free (NcHIPertGravTScalarInfo *Tsinfo)
{
  g_array_unref (Tsinfo->drho_deps);
  g_array_unref (Tsinfo->rhoppv_deps);
  g_array_unref (Tsinfo->dp_deps);
  g_array_unref (Tsinfo->dPi_deps);

  g_free (Tsinfo);
}

NCM_INLINE void 
nc_hipert_grav_T_scalar_info_append (NcHIPertGravTScalarInfo *Tsinfo, NcHIPertGravTScalarInfo *Tsinfo1)
{
  g_array_append_vals (Tsinfo->drho_deps,   Tsinfo1->drho_deps->data,   Tsinfo1->drho_deps->len);
  g_array_append_vals (Tsinfo->rhoppv_deps, Tsinfo1->rhoppv_deps->data, Tsinfo1->rhoppv_deps->len);
  g_array_append_vals (Tsinfo->dp_deps,     Tsinfo1->dp_deps->data,     Tsinfo1->dp_deps->len);
  g_array_append_vals (Tsinfo->dPi_deps,    Tsinfo1->dPi_deps->data,    Tsinfo1->dPi_deps->len);
}

#define __NC_HIPERT_COMP_ADD_PAD(a,pad) \
G_BEGIN_DECLS { \
  gint __i; \
  for (__i = 0; __i < (a)->len; __i++) \
  { \
    if (g_array_index ((a), gint, __i) >= 0) \
      { \
        g_array_index ((a), gint, __i) += pad; \
      } \
  } \
} G_END_DECLS

NCM_INLINE void 
nc_hipert_grav_T_scalar_info_add_pad (NcHIPertGravTScalarInfo *Tsinfo, gint pad)
{
  if (pad == 0)
    return;

  __NC_HIPERT_COMP_ADD_PAD (Tsinfo->drho_deps,   pad);
  __NC_HIPERT_COMP_ADD_PAD (Tsinfo->rhoppv_deps, pad);
  __NC_HIPERT_COMP_ADD_PAD (Tsinfo->dp_deps,     pad);
  __NC_HIPERT_COMP_ADD_PAD (Tsinfo->dPi_deps,    pad);
}

NCM_INLINE void 
nc_hipert_grav_T_scalar_info_set_zero (NcHIPertGravTScalarInfo *Tsinfo)
{
  g_array_set_size (Tsinfo->drho_deps,   0);
  g_array_set_size (Tsinfo->rhoppv_deps, 0);
  g_array_set_size (Tsinfo->dp_deps,     0);
  g_array_set_size (Tsinfo->dPi_deps,    0);
}


NCM_INLINE NcHIPertBGVarID 
nc_hipert_grav_get_id (NcHIPertGrav *grav)
{
  const NcHIPertBGVarID id = nc_hipert_bg_var_class_get_id_by_ns (G_OBJECT_TYPE_NAME (grav));
  g_assert_cmpint (id, >=, 0);
  return id;
}

NCM_INLINE guint 
nc_hipert_grav_ndyn_var (NcHIPertGrav *grav)
{
  return NC_HIPERT_GRAV_GET_CLASS (grav)->ndyn_var (grav);
}

NCM_INLINE GArray *
nc_hipert_grav_get_deps (NcHIPertGrav *grav, guint vindex)
{
  return NC_HIPERT_GRAV_GET_CLASS (grav)->get_deps (grav, vindex);
}

NCM_INLINE void 
nc_hipert_grav_set_gauge (NcHIPertGrav *grav, NcHIPertGravGauge gauge)
{
  NC_HIPERT_GRAV_GET_CLASS (grav)->set_gauge (grav, gauge);
}

NCM_INLINE NcHIPertGravGauge 
nc_hipert_grav_get_gauge (NcHIPertGrav *grav)
{
  return NC_HIPERT_GRAV_GET_CLASS (grav)->get_gauge (grav);
}

NCM_INLINE NcHIPertGravInfo *
nc_hipert_grav_get_G_scalar_info (NcHIPertGrav *grav)
{
  return NC_HIPERT_GRAV_GET_CLASS (grav)->get_G_scalar_info (grav);
}

NCM_INLINE void
nc_hipert_grav_get_G_scalar (NcHIPertGrav *grav, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar, NcHIPertGravScalar *G_scalar)
{
  return NC_HIPERT_GRAV_GET_CLASS (grav)->get_G_scalar (grav, bg_var, ydy, T_scalar, G_scalar);
}

NCM_INLINE void
nc_hipert_grav_get_dy_scalar (NcHIPertGrav *grav, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar, NcHIPertGravScalar *G_scalar)
{
  return NC_HIPERT_GRAV_GET_CLASS (grav)->get_dy_scalar (grav, bg_var, ydy, T_scalar, G_scalar);
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_HIPERT_GRAV_INLINE_H_ */
