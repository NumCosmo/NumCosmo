/***************************************************************************
 *            nc_hipert_two_fluids.h
 *
 *  Tue June 10 19:15:11 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_two_fluids.h
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

#ifndef _NC_HIPERT_TWO_FLUIDS_H_
#define _NC_HIPERT_TWO_FLUIDS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_ode_spline.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/perturbations/nc_hipert.h>
#include <numcosmo/perturbations/nc_hipert_wkb.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_TWO_FLUIDS             (nc_hipert_two_fluids_get_type ())
#define NC_HIPERT_TWO_FLUIDS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT_TWO_FLUIDS, NcHIPertTwoFluids))
#define NC_HIPERT_TWO_FLUIDS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPERT_TWO_FLUIDS, NcHIPertTwoFluidsClass))
#define NC_IS_HIPERT_TWO_FLUIDS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT_TWO_FLUIDS))
#define NC_IS_HIPERT_TWO_FLUIDS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPERT_TWO_FLUIDS))
#define NC_HIPERT_TWO_FLUIDS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPERT_TWO_FLUIDS, NcHIPertTwoFluidsClass))

typedef struct _NcHIPertTwoFluidsClass NcHIPertTwoFluidsClass;
typedef struct _NcHIPertTwoFluids NcHIPertTwoFluids;

struct _NcHIPertTwoFluidsClass
{
  /*< private >*/
  NcHIPertClass parent_class;
};

struct _NcHIPertTwoFluids
{
  /*< private >*/
  NcHIPert parent_instance;
  NcHIPertWKB *wkb_zeta;
  NcHIPertWKB *wkb_S;
};

/**
 * NcHIPertTwoFluidsVars:
 * @NC_HIPERT_TWO_FLUIDS_RE_ZETA: $\text{Re}(\zeta)$
 * @NC_HIPERT_TWO_FLUIDS_IM_ZETA: $\text{Im}(\zeta)$
 * @NC_HIPERT_TWO_FLUIDS_RE_PZETA: $\text{Re}(P_\zeta)$
 * @NC_HIPERT_TWO_FLUIDS_IM_PZETA: $\text{Im}(P_\zeta)$
 * @NC_HIPERT_TWO_FLUIDS_RE_Q: $\text{Re}(Q)$
 * @NC_HIPERT_TWO_FLUIDS_IM_Q: $\text{Im}(Q)$
 * @NC_HIPERT_TWO_FLUIDS_RE_PQ: $\text{Re}(P_Q)$
 * @NC_HIPERT_TWO_FLUIDS_IM_PQ: $\text{Im}(P_Q)$
 * 
 * Perturbation variables enumerator.
 * 
 */
typedef enum _NcHIPertTwoFluidsVars
{
  NC_HIPERT_TWO_FLUIDS_RE_ZETA = 0,
  NC_HIPERT_TWO_FLUIDS_IM_ZETA,
  NC_HIPERT_TWO_FLUIDS_RE_PZETA,
  NC_HIPERT_TWO_FLUIDS_IM_PZETA,
  NC_HIPERT_TWO_FLUIDS_RE_Q,
  NC_HIPERT_TWO_FLUIDS_IM_Q,
  NC_HIPERT_TWO_FLUIDS_RE_PQ,
  NC_HIPERT_TWO_FLUIDS_IM_PQ, /*< private >*/
  NC_HIPERT_TWO_FLUIDS_LEN,   /*< skip >*/
} NcHIPertTwoFluidsVars;

GType nc_hipert_two_fluids_get_type (void) G_GNUC_CONST;

NcHIPertTwoFluids *nc_hipert_two_fluids_new (void);
NcHIPertTwoFluids *nc_hipert_two_fluids_ref (NcHIPertTwoFluids *ptf);
void nc_hipert_two_fluids_free (NcHIPertTwoFluids *ptf);
void nc_hipert_two_fluids_clear (NcHIPertTwoFluids **ptf);

void nc_hipert_two_fluids_prepare_wkb_zeta (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble prec, gdouble alpha_i, gdouble alpha_f);
void nc_hipert_two_fluids_prepare_wkb_S (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble prec, gdouble alpha_i, gdouble alpha_f);

gdouble nc_hipert_two_fluids_nuA (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha);
gdouble nc_hipert_two_fluids_nuB (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha);
gdouble nc_hipert_two_fluids_YAB (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha);

void nc_hipert_two_fluids_wkb_zeta (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, gdouble *Re_zeta, gdouble *Im_zeta);
void nc_hipert_two_fluids_wkb_zeta_Pzeta (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, gdouble *Re_zeta, gdouble *Im_zeta, gdouble *Re_Pzeta, gdouble *Im_Pzeta);

void nc_hipert_two_fluids_wkb_Q (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, gdouble *Re_Q, gdouble *Im_Q);
void nc_hipert_two_fluids_wkb_Q_PQ (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, gdouble *Re_Q, gdouble *Im_Q, gdouble *Re_PQ, gdouble *Im_PQ);

void nc_hipert_two_fluids_wkb_full_zeta (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, gdouble **vars);
void nc_hipert_two_fluids_wkb_full_Q (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, gdouble **vars);

gdouble nc_hipert_two_fluids_wkb_zeta_maxtime (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha0, gdouble alpha1);
gdouble nc_hipert_two_fluids_wkb_S_maxtime (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha0, gdouble alpha1);
gdouble nc_hipert_two_fluids_wkb_zeta_maxtime_prec (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, NcHIPertWKBCmp cmp, gdouble prec, gdouble alpha0, gdouble alpha1);
gdouble nc_hipert_two_fluids_wkb_S_maxtime_prec (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, NcHIPertWKBCmp cmp, gdouble prec, gdouble alpha0, gdouble alpha1);

void nc_hipert_two_fluids_set_init_cond (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alphai, gdouble *vars);

void nc_hipert_two_fluids_set_init_cond_wkb_zeta (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alphai);

void nc_hipert_two_fluids_set_init_cond_wkb_Q (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alphai);

void nc_hipert_two_fluids_evolve (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alphaf);
void nc_hipert_two_fluids_get_values (NcHIPertTwoFluids *ptf, gdouble *alphai, gdouble **vars);

G_END_DECLS

#endif /* _NC_HIPERT_TWO_FLUIDS_H_ */
