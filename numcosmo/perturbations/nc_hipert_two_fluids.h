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
#ifndef NUMCOSMO_GIR_SCAN
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_ode_spline.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/perturbations/nc_hipert.h>
#include <numcosmo/perturbations/nc_hipert_wkb.h>
#include <numcosmo/perturbations/nc_hipert_itwo_fluids.h>

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
  N_Vector abstol;
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

void nc_hipert_two_fluids_eom (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, NcHIPertITwoFluidsEOM **eom);

void nc_hipert_two_fluids_get_init_cond (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, guint main_mode, const gdouble beta_R, gdouble *state);

void nc_hipert_two_fluids_set_init_cond (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alphai, gdouble *vars);

void nc_hipert_two_fluids_set_init_cond_wkb_zeta (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alphai);

void nc_hipert_two_fluids_set_init_cond_wkb_Q (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alphai);

void nc_hipert_two_fluids_to_zeta_s (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, gdouble *state);

void nc_hipert_two_fluids_evolve (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alphaf);
void nc_hipert_two_fluids_get_values (NcHIPertTwoFluids *ptf, gdouble *alphai, gdouble **vars);

#define NC_HIPERT_TWO_FLUIDS_A2Q(Ai) (cimag (Ai)) 
#define NC_HIPERT_TWO_FLUIDS_A2P(Ai) (creal (Ai)) 
#define NC_HIPERT_TWO_FLUIDS_QP2A(Q,P) ((P) + I * (Q))

G_END_DECLS

#endif /* _NC_HIPERT_TWO_FLUIDS_H_ */
