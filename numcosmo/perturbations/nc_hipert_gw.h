/***************************************************************************
 *            nc_hipert_gw.h
 *
 *  Fri December 09 11:24:53 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hipert_gw.h
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
#include <numcosmo/math/ncm_csq1d.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_IGW (nc_hipert_igw_get_type ())
#define NC_TYPE_HIPERT_GW (nc_hipert_gw_get_type ())

G_DECLARE_INTERFACE (NcHIPertIGW, nc_hipert_igw, NC, HIPERT_IGW, GObject)
G_DECLARE_FINAL_TYPE (NcHIPertGW, nc_hipert_gw, NC, HIPERT_GW, NcmCSQ1D)

struct _NcHIPertIGWInterface
{
  /*< private >*/
  GTypeInterface parent;
  gdouble (*eval_xi) (NcHIPertIGW *igw, const gdouble tau, const gdouble k);
  gdouble (*eval_F1) (NcHIPertIGW *igw, const gdouble tau, const gdouble k);
  gdouble (*eval_nu) (NcHIPertIGW *igw, const gdouble tau, const gdouble k);
  gdouble (*eval_m) (NcHIPertIGW *igw, const gdouble tau, const gdouble k);
  gdouble (*eval_unit) (NcHIPertIGW *igw);
  gdouble (*eval_x) (NcHIPertIGW *igw, const gdouble tau);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[12];
};

/**
 * NcHIPertGWVars:
 * @NC_HIPERT_GW_RE_H: $\text{Re}(\zeta)$
 * @NC_HIPERT_GW_IM_H: $\text{Im}(\zeta)$
 * @NC_HIPERT_GW_RE_PH: $\text{Re}(P_\zeta)$
 * @NC_HIPERT_GW_IM_PH: $\text{Im}(P_\zeta)$
 *
 * Perturbation variables enumerator.
 *
 */
typedef enum /*< enum,underscore_name=NC_HIPERT_GW_VARS  >*/
{
  NC_HIPERT_GW_RE_H = 0,
  NC_HIPERT_GW_IM_H,
  NC_HIPERT_GW_RE_PH,
  NC_HIPERT_GW_IM_PH,
} NcHIPertGWVars;


gdouble nc_hipert_igw_eval_xi (NcHIPertIGW *igw, const gdouble tau, const gdouble k);
gdouble nc_hipert_igw_eval_F1 (NcHIPertIGW *igw, const gdouble tau, const gdouble k);
gdouble nc_hipert_igw_eval_nu (NcHIPertIGW *igw, const gdouble tau, const gdouble k);
gdouble nc_hipert_igw_eval_m (NcHIPertIGW *igw, const gdouble tau, const gdouble k);
gdouble nc_hipert_igw_eval_unit (NcHIPertIGW *igw);
gdouble nc_hipert_igw_eval_x (NcHIPertIGW *igw, const gdouble tau);

NcHIPertGW *nc_hipert_gw_new (void);
NcHIPertGW *nc_hipert_gw_ref (NcHIPertGW *pgw);
void nc_hipert_gw_free (NcHIPertGW *pgw);
void nc_hipert_gw_clear (NcHIPertGW **pgw);

void nc_hipert_gw_set_k (NcHIPertGW *pgw, const gdouble k);
gdouble nc_hipert_gw_get_k (NcHIPertGW *pgw);

G_END_DECLS

#endif /* _NC_HIPERT_GW_H_ */

