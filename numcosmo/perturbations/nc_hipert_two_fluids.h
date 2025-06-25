/***************************************************************************
 *            nc_hipert_two_fluids.h
 *
 *  Tue June 10 19:15:11 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hipert_two_fluids.h
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

#define NC_TYPE_HIPERT_TWO_FLUIDS (nc_hipert_two_fluids_get_type ())

G_DECLARE_FINAL_TYPE (NcHIPertTwoFluids, nc_hipert_two_fluids, NC, HIPERT_TWO_FLUIDS, NcHIPert);

/* State interpolation */

typedef struct _NcHIPertTwoFluidsStateInterp NcHIPertTwoFluidsStateInterp;

struct _NcHIPertTwoFluidsStateInterp
{
  NcmSpline *mode1_splines[NC_HIPERT_ITWO_FLUIDS_VARS_LEN];
  NcmSpline *mode2_splines[NC_HIPERT_ITWO_FLUIDS_VARS_LEN];
  NcHIPertITwoFluidsState state;
  gint interp_mode; /* 1 == time, 2 == wavenumber */
};

GType nc_hipert_two_fluids_state_interp_get_type (void) G_GNUC_CONST;

NcHIPertTwoFluidsStateInterp *nc_hipert_two_fluids_state_interp_dup (NcHIPertTwoFluidsStateInterp *sinterp);
void nc_hipert_two_fluids_state_interp_free (NcHIPertTwoFluidsStateInterp *sinterp);

NcHIPertITwoFluidsState *nc_hipert_two_fluids_state_interp_eval (NcHIPertTwoFluidsStateInterp *sinterp, NcHICosmo *cosmo, gdouble x);

/* Two fluids */

NcHIPertTwoFluids *nc_hipert_two_fluids_new (void);
NcHIPertTwoFluids *nc_hipert_two_fluids_ref (NcHIPertTwoFluids *ptf);
void nc_hipert_two_fluids_free (NcHIPertTwoFluids *ptf);
void nc_hipert_two_fluids_clear (NcHIPertTwoFluids **ptf);

void nc_hipert_two_fluids_set_wkb_reltol (NcHIPertTwoFluids *ptf, gdouble reltol);
void nc_hipert_two_fluids_set_initial_time (NcHIPertTwoFluids *ptf, gdouble alpha_i);
void nc_hipert_two_fluids_set_final_time (NcHIPertTwoFluids *ptf, gdouble alpha_f);

gdouble nc_hipert_two_fluids_get_wkb_reltol (NcHIPertTwoFluids *ptf);
gdouble nc_hipert_two_fluids_get_initial_time (NcHIPertTwoFluids *ptf);
gdouble nc_hipert_two_fluids_get_final_time (NcHIPertTwoFluids *ptf);

void nc_hipert_two_fluids_eom (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, NcHIPertITwoFluidsEOM **eom);
void nc_hipert_two_fluids_wkb (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, NcHIPertITwoFluidsWKB **wkb);
gdouble nc_hipert_two_fluids_get_wkb_limit (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, guint main_mode, gdouble alpha_i, gdouble prec);

void nc_hipert_two_fluids_get_init_cond_QP (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, guint main_mode, const gdouble beta_R, NcmVector *init_cond);
void nc_hipert_two_fluids_get_init_cond_zetaS (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, guint main_mode, const gdouble beta_R, NcmVector *init_cond);
void nc_hipert_two_fluids_set_init_cond (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, guint main_mode, gboolean useQP, NcmVector *init_cond);
void nc_hipert_two_fluids_to_zeta_s (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, NcmVector *state);
void nc_hipert_two_fluids_evolve (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alphaf);

NcmVector *nc_hipert_two_fluids_peek_state (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble *alpha);
NcmMatrix *nc_hipert_two_fluids_evolve_array (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alphaf, gdouble step_reltol, gdouble step_abstol);
NcmSpline *nc_hipert_two_fluids_compute_zeta_spectrum (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, guint mode, gdouble alpha_i, gdouble alpha, gdouble ki, gdouble kf, guint nnodes);

NcHIPertTwoFluidsStateInterp *nc_hipert_two_fluids_evol_mode (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo);
NcHIPertTwoFluidsStateInterp *nc_hipert_two_fluids_state_compute_spectrum (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, const gdouble alpha, const gdouble k_i, const gdouble k_f, guint nnodes);

#define NC_HIPERT_TWO_FLUIDS_A2Q(Ai) (cimag (Ai))
#define NC_HIPERT_TWO_FLUIDS_A2P(Ai) (creal (Ai))
#define NC_HIPERT_TWO_FLUIDS_QP2A(Q, P) ((P) + I * (Q))

G_END_DECLS

#endif /* _NC_HIPERT_TWO_FLUIDS_H_ */

