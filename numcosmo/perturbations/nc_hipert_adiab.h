/***************************************************************************
 *            nc_hipert_adiab.h
 *
 *  Tue June 03 17:20:57 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_adiab.h
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

#ifndef _NC_HIPERT_ADIAB_H_
#define _NC_HIPERT_ADIAB_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_ode_spline.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/math/ncm_hoaa.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_IADIAB               (nc_hipert_iadiab_get_type ())
#define NC_HIPERT_IADIAB(obj)               (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT_IADIAB, NcHIPertIAdiab))
#define NC_IS_HIPERT_IADIAB(obj)            (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT_IADIAB))
#define NC_HIPERT_IADIAB_GET_INTERFACE(obj) (G_TYPE_INSTANCE_GET_INTERFACE ((obj), NC_TYPE_HIPERT_IADIAB, NcHIPertIAdiabInterface))

typedef struct _NcHIPertIAdiab NcHIPertIAdiab;
typedef struct _NcHIPertIAdiabInterface NcHIPertIAdiabInterface;

struct _NcHIPertIAdiabInterface
{
  /*< private >*/
  GTypeInterface parent;
  gdouble (*eval_mnu) (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
  gdouble (*eval_nu) (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
  gdouble (*eval_dlnmnu) (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
  void (*eval_system) (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k, gdouble *nu, gdouble *dlnmnu);
  guint (*nsing) (NcHIPertIAdiab *iad, const gdouble k);
  void (*get_sing_info) (NcHIPertIAdiab *iad, const gdouble k, const guint sing, gdouble *ts, gdouble *dts_i, gdouble *dts_f, NcmHOAASingType *st);
  gdouble (*eval_sing_mnu) (NcHIPertIAdiab *iad, const gdouble tau_m_taus, const gdouble k, const guint sing);
  gdouble (*eval_sing_dlnmnu) (NcHIPertIAdiab *iad, const gdouble tau_m_taus, const gdouble k, const guint sing);
  void (*eval_sing_system) (NcHIPertIAdiab *iad, const gdouble tau_m_taus, const gdouble k, const guint sing, gdouble *nu, gdouble *dlnmnu);
	gdouble (*eval_powspec_factor) (NcHIPertIAdiab *iad);
};

#define NC_TYPE_HIPERT_ADIAB             (nc_hipert_adiab_get_type ())
#define NC_HIPERT_ADIAB(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT_ADIAB, NcHIPertAdiab))
#define NC_HIPERT_ADIAB_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPERT_ADIAB, NcHIPertAdiabClass))
#define NC_IS_HIPERT_ADIAB(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT_ADIAB))
#define NC_IS_HIPERT_ADIAB_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPERT_ADIAB))
#define NC_HIPERT_ADIAB_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPERT_ADIAB, NcHIPertAdiabClass))

typedef struct _NcHIPertAdiabClass NcHIPertAdiabClass;
typedef struct _NcHIPertAdiab NcHIPertAdiab;

struct _NcHIPertAdiabClass
{
  /*< private >*/
  NcmHOAAClass parent_class;
};

/**
 * NcHIPertAdiabVars:
 * @NC_HIPERT_ADIAB_RE_ZETA: $\text{Re}(\zeta)$
 * @NC_HIPERT_ADIAB_IM_ZETA: $\text{Im}(\zeta)$
 * @NC_HIPERT_ADIAB_RE_PZETA: $\text{Re}(P_\zeta)$
 * @NC_HIPERT_ADIAB_IM_PZETA: $\text{Im}(P_\zeta)$
 * 
 * Perturbation variables enumerator.
 * 
 */
typedef enum /*< enum,underscore_name=NC_HIPERT_ADIAB_VARS  >*/
{
  NC_HIPERT_ADIAB_RE_ZETA = 0,
  NC_HIPERT_ADIAB_IM_ZETA,
  NC_HIPERT_ADIAB_RE_PZETA,
  NC_HIPERT_ADIAB_IM_PZETA,
} NcHIPertAdiabVars;

struct _NcHIPertAdiab
{
  /*< private >*/
  NcmHOAA parent_instance;
};

GType nc_hipert_iadiab_get_type (void) G_GNUC_CONST;
GType nc_hipert_adiab_get_type (void) G_GNUC_CONST;

NCM_INLINE gdouble nc_hipert_iadiab_eval_mnu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
NCM_INLINE gdouble nc_hipert_iadiab_eval_nu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
NCM_INLINE gdouble nc_hipert_iadiab_eval_dlnmnu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
NCM_INLINE void nc_hipert_iadiab_eval_system (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k, gdouble *nu, gdouble *dlnmnu);

NCM_INLINE guint nc_hipert_iadiab_nsing (NcHIPertIAdiab *iad, const gdouble k);
NCM_INLINE void nc_hipert_iadiab_get_sing_info (NcHIPertIAdiab *iad, const gdouble k, const guint sing, gdouble *ts, gdouble *dts_i, gdouble *dts_f, NcmHOAASingType *st);
NCM_INLINE gdouble nc_hipert_iadiab_eval_sing_mnu (NcHIPertIAdiab *iad, const gdouble tau_m_taus, const gdouble k, const guint sing);
NCM_INLINE gdouble nc_hipert_iadiab_eval_sing_dlnmnu (NcHIPertIAdiab *iad, const gdouble tau_m_taus, const gdouble k, const guint sing);
NCM_INLINE void nc_hipert_iadiab_eval_sing_system (NcHIPertIAdiab *iad, const gdouble tau_m_taus, const gdouble k, const guint sing, gdouble *nu, gdouble *dlnmnu);

NCM_INLINE gdouble nc_hipert_iadiab_eval_powspec_factor (NcHIPertIAdiab *iad);

NcHIPertAdiab *nc_hipert_adiab_new (void);
NcHIPertAdiab *nc_hipert_adiab_ref (NcHIPertAdiab *pa);
void nc_hipert_adiab_free (NcHIPertAdiab *pa);
void nc_hipert_adiab_clear (NcHIPertAdiab **pa);

G_END_DECLS

#endif /* _NC_HIPERT_ADIAB_H_ */

#ifndef _NC_HIPERT_ADIAB_INLINE_H_
#define _NC_HIPERT_ADIAB_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE gdouble
nc_hipert_iadiab_eval_mnu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  return NC_HIPERT_IADIAB_GET_INTERFACE (iad)->eval_mnu (iad, tau, k);
}

NCM_INLINE gdouble
nc_hipert_iadiab_eval_nu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  return NC_HIPERT_IADIAB_GET_INTERFACE (iad)->eval_nu (iad, tau, k);
}

NCM_INLINE gdouble
nc_hipert_iadiab_eval_dlnmnu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  return NC_HIPERT_IADIAB_GET_INTERFACE (iad)->eval_dlnmnu (iad, tau, k);
}

NCM_INLINE void
nc_hipert_iadiab_eval_system (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k, gdouble *nu, gdouble *dlnmnu)
{
  NC_HIPERT_IADIAB_GET_INTERFACE (iad)->eval_system (iad, tau, k, nu, dlnmnu);
}

NCM_INLINE guint 
nc_hipert_iadiab_nsing (NcHIPertIAdiab *iad, const gdouble k)
{
  return NC_HIPERT_IADIAB_GET_INTERFACE (iad)->nsing (iad, k);
}

NCM_INLINE void
nc_hipert_iadiab_get_sing_info (NcHIPertIAdiab *iad, const gdouble k, const guint sing, gdouble *ts, gdouble *dts_i, gdouble *dts_f, NcmHOAASingType *st)
{
  NC_HIPERT_IADIAB_GET_INTERFACE (iad)->get_sing_info (iad, k, sing, ts, dts_i, dts_f, st);
}

NCM_INLINE gdouble 
nc_hipert_iadiab_eval_sing_mnu (NcHIPertIAdiab *iad, const gdouble tau_m_taus, const gdouble k, const guint sing)
{
  return NC_HIPERT_IADIAB_GET_INTERFACE (iad)->eval_sing_mnu (iad, tau_m_taus, k, sing);
}

NCM_INLINE gdouble 
nc_hipert_iadiab_eval_sing_dlnmnu (NcHIPertIAdiab *iad, const gdouble tau_m_taus, const gdouble k, const guint sing)
{
  return NC_HIPERT_IADIAB_GET_INTERFACE (iad)->eval_sing_dlnmnu (iad, tau_m_taus, k, sing);
}

NCM_INLINE void 
nc_hipert_iadiab_eval_sing_system (NcHIPertIAdiab *iad, const gdouble tau_m_taus, const gdouble k, const guint sing, gdouble *nu, gdouble *dlnmnu)
{
  NC_HIPERT_IADIAB_GET_INTERFACE (iad)->eval_sing_system (iad, tau_m_taus, k, sing, nu, dlnmnu);
}

NCM_INLINE gdouble 
nc_hipert_iadiab_eval_powspec_factor (NcHIPertIAdiab *iad)
{
  return NC_HIPERT_IADIAB_GET_INTERFACE (iad)->eval_powspec_factor (iad);
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_HIPERT_ADIAB_INLINE_H_ */
