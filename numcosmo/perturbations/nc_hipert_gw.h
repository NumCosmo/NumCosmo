/***************************************************************************
 *            nc_hipert_gw.h
 *
 *  Fri December 09 11:24:53 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_gw.h
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_HIPERT_GW_H_
#define _NC_HIPERT_GW_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_ode_spline.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/math/ncm_hoaa.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_IGW               (nc_hipert_igw_get_type ())
#define NC_HIPERT_IGW(obj)               (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT_IGW, NcHIPertIGW))
#define NC_IS_HIPERT_IGW(obj)            (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT_IGW))
#define NC_HIPERT_IGW_GET_INTERFACE(obj) (G_TYPE_INSTANCE_GET_INTERFACE ((obj), NC_TYPE_HIPERT_IGW, NcHIPertIGWInterface))

typedef struct _NcHIPertIGW NcHIPertIGW;
typedef struct _NcHIPertIGWInterface NcHIPertIGWInterface;

struct _NcHIPertIGWInterface
{
  /*< private >*/
  GTypeInterface parent;
  gdouble (*eval_mnu) (NcHIPertIGW *igw, const gdouble tau, const gdouble k);
  gdouble (*eval_nu) (NcHIPertIGW *igw, const gdouble tau, const gdouble k);
  gdouble (*eval_dlnmnu) (NcHIPertIGW *igw, const gdouble tau, const gdouble k);
  void (*eval_system) (NcHIPertIGW *igw, const gdouble tau, const gdouble k, gdouble *nu, gdouble *dlnmnu);
  guint (*nsing) (NcHIPertIGW *igw, const gdouble k);
  void (*get_sing_info) (NcHIPertIGW *igw, const gdouble k, const guint sing, gdouble *ts, gdouble *dts_i, gdouble *dts_f, NcmHOAASingType *st);
  gdouble (*eval_sing_mnu) (NcHIPertIGW *igw, const gdouble tau_m_taus, const gdouble k, const guint sing);
  gdouble (*eval_sing_dlnmnu) (NcHIPertIGW *igw, const gdouble tau_m_taus, const gdouble k, const guint sing);
  void (*eval_sing_system) (NcHIPertIGW *igw, const gdouble tau_m_taus, const gdouble k, const guint sing, gdouble *nu, gdouble *dlnmnu);
	gdouble (*eval_powspec_factor) (NcHIPertIGW *igw);
};

#define NC_TYPE_HIPERT_GW             (nc_hipert_gw_get_type ())
#define NC_HIPERT_GW(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT_GW, NcHIPertGW))
#define NC_HIPERT_GW_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPERT_GW, NcHIPertGWClass))
#define NC_IS_HIPERT_GW(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT_GW))
#define NC_IS_HIPERT_GW_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPERT_GW))
#define NC_HIPERT_GW_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPERT_GW, NcHIPertGWClass))

typedef struct _NcHIPertGWClass NcHIPertGWClass;
typedef struct _NcHIPertGW NcHIPertGW;

struct _NcHIPertGWClass
{
  /*< private >*/
  NcmHOAAClass parent_class;
};

/**
 * NcHIPertGWVars:
 * @NC_HIPERT_GW_RE_ZETA: $\text{Re}(\zeta)$
 * @NC_HIPERT_GW_IM_ZETA: $\text{Im}(\zeta)$
 * @NC_HIPERT_GW_RE_PZETA: $\text{Re}(P_\zeta)$
 * @NC_HIPERT_GW_IM_PZETA: $\text{Im}(P_\zeta)$
 * 
 * Perturbation variables enumerator.
 * 
 */
typedef enum /*< enum,underscore_name=NC_HIPERT_GW_VARS  >*/
{
  NC_HIPERT_GW_RE_ZETA = 0,
  NC_HIPERT_GW_IM_ZETA,
  NC_HIPERT_GW_RE_PZETA,
  NC_HIPERT_GW_IM_PZETA,
} NcHIPertGWVars;

struct _NcHIPertGW
{
  /*< private >*/
  NcmHOAA parent_instance;
};

GType nc_hipert_igw_get_type (void) G_GNUC_CONST;
GType nc_hipert_gw_get_type (void) G_GNUC_CONST;

NCM_INLINE gdouble nc_hipert_igw_eval_mnu (NcHIPertIGW *igw, const gdouble tau, const gdouble k);
NCM_INLINE gdouble nc_hipert_igw_eval_nu (NcHIPertIGW *igw, const gdouble tau, const gdouble k);
NCM_INLINE gdouble nc_hipert_igw_eval_dlnmnu (NcHIPertIGW *igw, const gdouble tau, const gdouble k);
NCM_INLINE void nc_hipert_igw_eval_system (NcHIPertIGW *igw, const gdouble tau, const gdouble k, gdouble *nu, gdouble *dlnmnu);

NCM_INLINE guint nc_hipert_igw_nsing (NcHIPertIGW *igw, const gdouble k);
NCM_INLINE void nc_hipert_igw_get_sing_info (NcHIPertIGW *igw, const gdouble k, const guint sing, gdouble *ts, gdouble *dts_i, gdouble *dts_f, NcmHOAASingType *st);
NCM_INLINE gdouble nc_hipert_igw_eval_sing_mnu (NcHIPertIGW *igw, const gdouble tau_m_taus, const gdouble k, const guint sing);
NCM_INLINE gdouble nc_hipert_igw_eval_sing_dlnmnu (NcHIPertIGW *igw, const gdouble tau_m_taus, const gdouble k, const guint sing);
NCM_INLINE void nc_hipert_igw_eval_sing_system (NcHIPertIGW *igw, const gdouble tau_m_taus, const gdouble k, const guint sing, gdouble *nu, gdouble *dlnmnu);

NCM_INLINE gdouble nc_hipert_igw_eval_powspec_factor (NcHIPertIGW *igw);

NcHIPertGW *nc_hipert_gw_new (void);
NcHIPertGW *nc_hipert_gw_ref (NcHIPertGW *pa);
void nc_hipert_gw_free (NcHIPertGW *pa);
void nc_hipert_gw_clear (NcHIPertGW **pa);

G_END_DECLS

#endif /* _NC_HIPERT_GW_H_ */

#ifndef _NC_HIPERT_GW_INLINE_H_
#define _NC_HIPERT_GW_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE gdouble
nc_hipert_igw_eval_mnu (NcHIPertIGW *igw, const gdouble tau, const gdouble k)
{
  return NC_HIPERT_IGW_GET_INTERFACE (igw)->eval_mnu (igw, tau, k);
}

NCM_INLINE gdouble
nc_hipert_igw_eval_nu (NcHIPertIGW *igw, const gdouble tau, const gdouble k)
{
  return NC_HIPERT_IGW_GET_INTERFACE (igw)->eval_nu (igw, tau, k);
}

NCM_INLINE gdouble
nc_hipert_igw_eval_dlnmnu (NcHIPertIGW *igw, const gdouble tau, const gdouble k)
{
  return NC_HIPERT_IGW_GET_INTERFACE (igw)->eval_dlnmnu (igw, tau, k);
}

NCM_INLINE void
nc_hipert_igw_eval_system (NcHIPertIGW *igw, const gdouble tau, const gdouble k, gdouble *nu, gdouble *dlnmnu)
{
  NC_HIPERT_IGW_GET_INTERFACE (igw)->eval_system (igw, tau, k, nu, dlnmnu);
}

NCM_INLINE guint 
nc_hipert_igw_nsing (NcHIPertIGW *igw, const gdouble k)
{
  return NC_HIPERT_IGW_GET_INTERFACE (igw)->nsing (igw, k);
}

NCM_INLINE void
nc_hipert_igw_get_sing_info (NcHIPertIGW *igw, const gdouble k, const guint sing, gdouble *ts, gdouble *dts_i, gdouble *dts_f, NcmHOAASingType *st)
{
  NC_HIPERT_IGW_GET_INTERFACE (igw)->get_sing_info (igw, k, sing, ts, dts_i, dts_f, st);
}

NCM_INLINE gdouble 
nc_hipert_igw_eval_sing_mnu (NcHIPertIGW *igw, const gdouble tau_m_taus, const gdouble k, const guint sing)
{
  return NC_HIPERT_IGW_GET_INTERFACE (igw)->eval_sing_mnu (igw, tau_m_taus, k, sing);
}

NCM_INLINE gdouble 
nc_hipert_igw_eval_sing_dlnmnu (NcHIPertIGW *igw, const gdouble tau_m_taus, const gdouble k, const guint sing)
{
  return NC_HIPERT_IGW_GET_INTERFACE (igw)->eval_sing_dlnmnu (igw, tau_m_taus, k, sing);
}

NCM_INLINE void 
nc_hipert_igw_eval_sing_system (NcHIPertIGW *igw, const gdouble tau_m_taus, const gdouble k, const guint sing, gdouble *nu, gdouble *dlnmnu)
{
  NC_HIPERT_IGW_GET_INTERFACE (igw)->eval_sing_system (igw, tau_m_taus, k, sing, nu, dlnmnu);
}

NCM_INLINE gdouble 
nc_hipert_igw_eval_powspec_factor (NcHIPertIGW *igw)
{
  return NC_HIPERT_IGW_GET_INTERFACE (igw)->eval_powspec_factor (igw);
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_HIPERT_GW_INLINE_H_ */
