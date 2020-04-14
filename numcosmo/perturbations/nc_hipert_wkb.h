/***************************************************************************
 *            nc_hipert_wkb.h
 *
 *  Sun August 03 20:39:19 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_wkb.h
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

#ifndef _NC_HIPERT_WKB_H_
#define _NC_HIPERT_WKB_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_ode_spline.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/perturbations/nc_hipert.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_WKB             (nc_hipert_wkb_get_type ())
#define NC_HIPERT_WKB(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT_WKB, NcHIPertWKB))
#define NC_HIPERT_WKB_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPERT_WKB, NcHIPertWKBClass))
#define NC_IS_HIPERT_WKB(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT_WKB))
#define NC_IS_HIPERT_WKB_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPERT_WKB))
#define NC_HIPERT_WKB_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPERT_WKB, NcHIPertWKBClass))

typedef struct _NcHIPertWKBClass NcHIPertWKBClass;
typedef struct _NcHIPertWKB NcHIPertWKB;

typedef void (*NcHIPertWKBGetNuV) (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k, gdouble *nu, gdouble *V);
typedef void (*NcHIPertWKBGetMnuDmnu) (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k, gdouble *mnu, gdouble *dmnu);
typedef gdouble (*NcHIPertWKBGetVal) (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k);

struct _NcHIPertWKBClass
{
  /*< private >*/
  NcHIPertClass parent_class;
  NcHIPertWKBGetNuV get_nu_V;
  NcHIPertWKBGetMnuDmnu get_mnu_dmnu;
  NcHIPertWKBGetVal get_m;
  NcHIPertWKBGetVal get_nu2;
  NcHIPertWKBGetVal get_dVnu2;
};

/**
 * NcHIPertWKBVars:
 * @NC_HIPERT_WKB_RE_Q: $\text{Re}(q)$
 * @NC_HIPERT_WKB_IM_Q: $\text{Im}(q)$
 * @NC_HIPERT_WKB_RE_P: $\text{Re}(p)$
 * @NC_HIPERT_WKB_IM_P: $\text{Im}(p)$
 * 
 * Perturbation variables enumerator.
 * 
 */
typedef enum /*< enum,underscore_name=NC_HIPERT_WKB_VARS >*/
{
  NC_HIPERT_WKB_RE_Q = 0,
  NC_HIPERT_WKB_IM_Q,
  NC_HIPERT_WKB_RE_P,
  NC_HIPERT_WKB_IM_P,
} NcHIPertWKBVars;

/**
 * NcHIPertWKBCmp:
 * @NC_HIPERT_WKB_CMP_POTENTIAL: Compare with the potential.
 * @NC_HIPERT_WKB_CMP_ALPHA2: Compare with $\alpha^2$.
 * 
 * FIXME
 * 
 */
typedef enum _NcHIPertWKBCmp
{
  NC_HIPERT_WKB_CMP_POTENTIAL = 0,
  NC_HIPERT_WKB_CMP_ALPHA2,
} NcHIPertWKBCmp;

struct _NcHIPertWKB
{
  /*< private >*/
  NcHIPert parent_instance;
  NcmSpline *nuA;
  NcmSpline *lnF;
  NcmSpline *dlnF;
  gdouble alpha_phase;
  gdouble cur_phase;
  gdouble alpha_i;
  gdouble alpha_f;
  gdouble alpha_p;
};

GType nc_hipert_wkb_get_type (void) G_GNUC_CONST;

NcHIPertWKB *nc_hipert_wkb_new_by_name (const gchar *wkb_name);
NcHIPertWKB *nc_hipert_wkb_ref (NcHIPertWKB *wkb);
void nc_hipert_wkb_free (NcHIPertWKB *wkb);
void nc_hipert_wkb_clear (NcHIPertWKB **wkb);

void nc_hipert_wkb_set_interval (NcHIPertWKB *wkb, gdouble alpha_i, gdouble alpha_f);

NCM_INLINE void nc_hipert_wkb_get_nu_V (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k, gdouble *nu, gdouble *V);
NCM_INLINE void nc_hipert_wkb_get_mnu_dmnu (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k, gdouble *mnu, gdouble *dmnu);
NCM_INLINE gdouble nc_hipert_wkb_get_m (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k);
NCM_INLINE gdouble nc_hipert_wkb_get_nu2 (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k);
NCM_INLINE gdouble nc_hipert_wkb_get_dVnu2 (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k);

void nc_hipert_wkb_prepare (NcHIPertWKB *wkb, NcmModel *model);

void nc_hipert_wkb_q (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble *Re_q, gdouble *Im_q);
void nc_hipert_wkb_q_p (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble *Re_q, gdouble *Im_q, gdouble *Re_p, gdouble *Im_p);

gdouble nc_hipert_wkb_nuA (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha);
gdouble nc_hipert_wkb_phase (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha);

gdouble nc_hipert_wkb_maxtime (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha0, gdouble alpha1);
gdouble nc_hipert_wkb_maxtime_prec (NcHIPertWKB *wkb, NcmModel *model, NcHIPertWKBCmp cmp, gdouble alpha0, gdouble alpha1);

G_END_DECLS

#endif /* _NC_HIPERT_WKB_H_ */

#ifndef _NC_HIPERT_WKB_INLINE_H_
#define _NC_HIPERT_WKB_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE void 
nc_hipert_wkb_get_nu_V (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k, gdouble *nu, gdouble *V)
{
  NC_HIPERT_WKB_GET_CLASS (wkb)->get_nu_V (wkb, model, alpha, k, nu, V);
}

NCM_INLINE void 
nc_hipert_wkb_get_mnu_dmnu (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k, gdouble *mnu, gdouble *dmnu)
{
  NC_HIPERT_WKB_GET_CLASS (wkb)->get_mnu_dmnu (wkb, model, alpha, k, mnu, dmnu);
}

NCM_INLINE gdouble 
nc_hipert_wkb_get_m (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k)
{
  return NC_HIPERT_WKB_GET_CLASS (wkb)->get_m (wkb, model, alpha, k);
}

NCM_INLINE gdouble 
nc_hipert_wkb_get_nu2 (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k)
{
  return NC_HIPERT_WKB_GET_CLASS (wkb)->get_nu2 (wkb, model, alpha, k);
}

NCM_INLINE gdouble 
nc_hipert_wkb_get_dVnu2 (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k)
{
  return NC_HIPERT_WKB_GET_CLASS (wkb)->get_dVnu2 (wkb, model, alpha, k);
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_HIPERT_WKB_INLINE_H_ */
