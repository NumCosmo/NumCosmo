/***************************************************************************
 *            nc_hipert_itwo_fluids.h
 *
 *  Tue July 22 17:37:11 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hipert_itwo_fluids.h
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

#ifndef _NC_HIPERT_ITWO_FLUIDS_H_
#define _NC_HIPERT_ITWO_FLUIDS_H_

#include <glib-object.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/perturbations/nc_hipert_wkb.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_ITWO_FLUIDS (nc_hipert_itwo_fluids_get_type ())

G_DECLARE_INTERFACE (NcHIPertITwoFluids, nc_hipert_itwo_fluids, NC, HIPERT_ITWO_FLUIDS, GObject)

typedef struct _NcHIPertITwoFluidsTV NcHIPertITwoFluidsTV;
typedef struct _NcHIPertITwoFluidsState NcHIPertITwoFluidsState;
typedef struct _NcHIPertITwoFluidsEOM NcHIPertITwoFluidsEOM;
typedef struct _NcHIPertITwoFluidsWKB NcHIPertITwoFluidsWKB;

typedef NcHIPertITwoFluidsEOM *(*NcHIPertITwoFluidsFuncEOM) (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);
typedef NcHIPertITwoFluidsWKB *(*NcHIPertITwoFluidsFuncWKB) (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);
typedef NcHIPertITwoFluidsTV *(*NcHIPertITwoFluidsFuncTV) (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);

struct _NcHIPertITwoFluidsInterface
{
  /*< private >*/
  GTypeInterface parent;
  NcHIPertITwoFluidsFuncEOM eom;
  NcHIPertITwoFluidsFuncWKB wkb;
  NcHIPertITwoFluidsFuncTV tv;
  gdouble (*eval_unit) (NcHIPertITwoFluids *itf);
};

/**
 * NcHIPertITwoFluidsEOM:
 *
 * Structure used to store the perturbations' equations of motion coefficients for the
 * two-fluid model.
 *
 */
struct _NcHIPertITwoFluidsEOM
{
  /*< private >*/
  guint64 skey;
  gdouble alpha;
  gdouble k;
  gdouble nu1;
  gdouble nu2;
  gdouble gammabar11;
  gdouble gammabar22;
  gdouble gammabar12;
  gdouble taubar;
  gdouble m_zeta;
  gdouble m_s;
  gdouble mnu2_zeta;
  gdouble mnu2_s;
  gdouble y;
  gdouble sin2phi;
  gdouble cos2phi;
  gdouble cs2;
  gdouble cm2;
  gdouble gw1;
  gdouble gw2;
  gdouble Fnu;
};

/**
 * NcHIPertITwoFluidsVars:
 * @NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R: Real part of the adiabatic perturbation
 * @NC_HIPERT_ITWO_FLUIDS_VARS_S_R:  Real part of the entropy perturbation
 * @NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R: Real part of the adiabatic perturbation's momentum
 * @NC_HIPERT_ITWO_FLUIDS_VARS_PS_R:  Real part of the entropy perturbation's momentum
 * @NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I: Imaginary part of the adiabatic perturbation
 * @NC_HIPERT_ITWO_FLUIDS_VARS_S_I:   Imaginary part of the entropy perturbation
 * @NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I: Imaginary part of the adiabatic perturbation's momentum
 * @NC_HIPERT_ITWO_FLUIDS_VARS_PS_I:   Imaginary part of the entropy perturbation's momentum
 *
 * Enumeration of variables used in the two-fluid model.
 *
 */
typedef enum /*< enum,underscore_name=NC_HIPERT_ITWO_FLUIDS_VARS >*/
{
  NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R = 0,
  NC_HIPERT_ITWO_FLUIDS_VARS_S_R,
  NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R,
  NC_HIPERT_ITWO_FLUIDS_VARS_PS_R,
  NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I,
  NC_HIPERT_ITWO_FLUIDS_VARS_S_I,
  NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I,
  NC_HIPERT_ITWO_FLUIDS_VARS_PS_I,
  /* < private > */
  NC_HIPERT_ITWO_FLUIDS_VARS_LEN, /*< skip >*/
} NcHIPertITwoFluidsVars;

#define NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1 NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R
#define NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2 NC_HIPERT_ITWO_FLUIDS_VARS_S_R
#define NC_HIPERT_ITWO_FLUIDS_VARS_P_R1 NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R
#define NC_HIPERT_ITWO_FLUIDS_VARS_P_R2 NC_HIPERT_ITWO_FLUIDS_VARS_PS_R
#define NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1 NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I
#define NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2 NC_HIPERT_ITWO_FLUIDS_VARS_S_I
#define NC_HIPERT_ITWO_FLUIDS_VARS_P_I1 NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I
#define NC_HIPERT_ITWO_FLUIDS_VARS_P_I2 NC_HIPERT_ITWO_FLUIDS_VARS_PS_I


/**
 * NcHIPertITwoFluidsObs:
 * @NC_HIPERT_ITWO_FLUIDS_OBS_ZETA: Adiabatic perturbation (curvature perturbation), $\\zeta$.
 * @NC_HIPERT_ITWO_FLUIDS_OBS_FKU_TOT: Total velocity potential, $F_k\\mathcal{V}$.
 * @NC_HIPERT_ITWO_FLUIDS_OBS_FKU_DIFF: Velocity potential difference, $k_\\mathrm{phys}\\mathcal{V}_\\mathrm{diff}$.
 * @NC_HIPERT_ITWO_FLUIDS_OBS_DELTA_TOT: Total density contrast, $\\delta_\\mathrm{tot}$.
 * @NC_HIPERT_ITWO_FLUIDS_OBS_DELTA_DIFF: Density contrast difference (isocurvature perturbation), $\\delta_\\mathrm{diff}$.
 * @NC_HIPERT_ITWO_FLUIDS_OBS_FKU_R: Velocity potential of the radiation component, $k_\\mathrm{phys}\\mathcal{V}_r$.
 * @NC_HIPERT_ITWO_FLUIDS_OBS_FKU_W: Velocity potential of the matter component, $k_\\mathrm{phys}\\mathcal{V}_w$.
 * @NC_HIPERT_ITWO_FLUIDS_OBS_DELTA_R: Density contrast of the radiation component, $\\delta_r$.
 * @NC_HIPERT_ITWO_FLUIDS_OBS_DELTA_W: Density contrast of the matter component, $\\delta_w$.
 *
 * Enumeration of physical observables computed from the two-fluid perturbation state.
 *
 * Notes:
 * - The quantity $k_\\mathrm{phys}\\mathcal{V}_\\mathrm{diff} =
 *   k_\\mathrm{phys}(\\mathcal{V}_r - \\mathcal{V}_w)$ reflects differences in velocity
 *   potential.
 * - The isocurvature perturbation is $\delta_{\text{diff}} = \delta_r - \delta_w$.
 * - The total curvature perturbation $k_\\mathrm{phys}\\mathcal{V}_{\text{tot}}$ is
 *   computed as a weighted average of the components.
 * - The total density contrast $\delta_{\text{tot}}$ may be computed as a weighted
 *   average of the components.
 *
 * These observables cover adiabatic, isocurvature, and component-specific quantities,
 * and are suitable for correlation and power spectrum evaluations.
 */
typedef enum /*< enum,underscore_name=NC_HIPERT_ITWO_FLUIDS_OBS >*/
{
  NC_HIPERT_ITWO_FLUIDS_OBS_ZETA,
  NC_HIPERT_ITWO_FLUIDS_OBS_FKU_TOT,
  NC_HIPERT_ITWO_FLUIDS_OBS_FKU_DIFF,
  NC_HIPERT_ITWO_FLUIDS_OBS_DELTA_TOT,
  NC_HIPERT_ITWO_FLUIDS_OBS_DELTA_DIFF,
  NC_HIPERT_ITWO_FLUIDS_OBS_FKU_R,
  NC_HIPERT_ITWO_FLUIDS_OBS_FKU_W,
  NC_HIPERT_ITWO_FLUIDS_OBS_DELTA_R,
  NC_HIPERT_ITWO_FLUIDS_OBS_DELTA_W,
  /* < private > */
  NC_HIPERT_ITWO_FLUIDS_OBS_LEN, /*< skip >*/
} NcHIPertITwoFluidsObs;

/**
 * NcHIPertITwoFluidsObsMode:
 * @NC_HIPERT_ITWO_FLUIDS_OBS_MODE_ONE: Include only the first quantized mode (associated with the first complex solution).
 * @NC_HIPERT_ITWO_FLUIDS_OBS_MODE_TWO: Include only the second quantized mode (associated with the second complex solution).
 * @NC_HIPERT_ITWO_FLUIDS_OBS_MODE_BOTH: Include both quantized modes (sum of their separate contributions).
 *
 * Specifies which quantized mode(s) to include when computing two-point observables
 * such as power spectra and correlations.
 *
 * In the quantized theory, the system is expanded in two linearly independent complex
 * solutions, each associated with a distinct pair of creation and annihilation
 * operators. These operator pairs correspond to independent degrees of freedom and
 * commute with each other. As a result, selecting @NC_HIPERT_ITWO_FLUIDS_OBS_MODE_BOTH
 * computes the sum of the individual contributions of both modes.
 */
typedef enum /*< enum,underscore_name=NC_HIPERT_ITWO_FLUIDS_OBS_MODE >*/
{
  NC_HIPERT_ITWO_FLUIDS_OBS_MODE_ONE,
  NC_HIPERT_ITWO_FLUIDS_OBS_MODE_TWO,
  NC_HIPERT_ITWO_FLUIDS_OBS_MODE_BOTH,
  /* < private > */
  NC_HIPERT_ITWO_FLUIDS_OBS_MODE_LEN /*< skip >*/
} NcHIPertITwoFluidsObsMode;

/**
 * NcHIPertITwoFluidsTV:
 *
 * Structure used to store the perturbations' transfer functions for the two-fluid model.
 *
 */
struct _NcHIPertITwoFluidsTV
{
  /*< private >*/
  guint64 skey;
  gdouble alpha;
  gdouble k;
  gdouble zeta[NC_HIPERT_ITWO_FLUIDS_VARS_LEN / 2];
  gdouble s[NC_HIPERT_ITWO_FLUIDS_VARS_LEN / 2];
  gdouble Pzeta[NC_HIPERT_ITWO_FLUIDS_VARS_LEN / 2];
  gdouble Ps[NC_HIPERT_ITWO_FLUIDS_VARS_LEN / 2];
};

/**
 * NcHIPertITwoFluidsState:
 *
 * Structure used to store the perturbations' state for the two-fluid model.
 */
struct _NcHIPertITwoFluidsState
{
  gdouble alpha;
  gdouble k;
  gdouble gw1;
  gdouble gw2;
  gdouble Fnu;
  gdouble norma;
#ifndef NUMCOSMO_GIR_SCAN
  complex double zeta1;
  complex double Q1;
  complex double Pzeta1;
  complex double PQ1;

  complex double zeta2;
  complex double Q2;
  complex double Pzeta2;
  complex double PQ2;
#endif /* NUMCOSMO_GIR_SCAN */
};

/**
 * NcHIPertITwoFluidsWKB:
 *
 * Structure representing the WKB approximations of perturbations in the two-fluid
 * model. Stores two approximations, one for each eigenmode, along with their
 * corresponding estimated WKB scales.
 *
 */
struct _NcHIPertITwoFluidsWKB
{
  /*< private >*/
  gdouble mode1_scale;
  gdouble mode2_scale;
  gdouble mode1_zeta_scale;
  gdouble mode2_zeta_scale;
  gdouble mode1_Q_scale;
  gdouble mode2_Q_scale;
  gdouble mode1_Pzeta_scale;
  gdouble mode2_Pzeta_scale;
  gdouble mode1_PQ_scale;
  gdouble mode2_PQ_scale;
  NcHIPertITwoFluidsState state;
};

GType nc_hipert_itwo_fluids_tv_get_type (void) G_GNUC_CONST;
GType nc_hipert_itwo_fluids_state_get_type (void) G_GNUC_CONST;
GType nc_hipert_itwo_fluids_eom_get_type (void) G_GNUC_CONST;
GType nc_hipert_itwo_fluids_wkb_get_type (void) G_GNUC_CONST;

NcHIPertITwoFluidsTV *nc_hipert_itwo_fluids_tv_dup (NcHIPertITwoFluidsTV *tf_tv);
void nc_hipert_itwo_fluids_tv_free (NcHIPertITwoFluidsTV *tf_tv);

NcHIPertITwoFluidsState *nc_hipert_itwo_fluids_state_dup (NcHIPertITwoFluidsState *tf_state);
void nc_hipert_itwo_fluids_state_free (NcHIPertITwoFluidsState *tf_state);

NcHIPertITwoFluidsEOM *nc_hipert_itwo_fluids_eom_dup (NcHIPertITwoFluidsEOM *tf_eom);
void nc_hipert_itwo_fluids_eom_free (NcHIPertITwoFluidsEOM *tf_eom);

NcHIPertITwoFluidsWKB *nc_hipert_itwo_fluids_wkb_dup (NcHIPertITwoFluidsWKB *tf_wkb);
void nc_hipert_itwo_fluids_wkb_free (NcHIPertITwoFluidsWKB *tf_wkb);

NCM_INLINE NcHIPertITwoFluidsEOM *nc_hipert_itwo_fluids_eom_eval (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);
NCM_INLINE NcHIPertITwoFluidsWKB *nc_hipert_itwo_fluids_wkb_eval (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);
NCM_INLINE NcHIPertITwoFluidsTV *nc_hipert_itwo_fluids_tv_eval (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);
NCM_INLINE gdouble nc_hipert_itwo_fluids_eval_unit (NcHIPertITwoFluids *itf);

gdouble nc_hipert_itwo_fluids_state_eval_obs (NcHIPertITwoFluidsState *tf_state, NcHIPertITwoFluidsObsMode obs_mode, NcHIPertITwoFluidsObs obs_a, NcHIPertITwoFluidsObs obs_b);

G_END_DECLS

#endif /* _NC_HIPERT_ITWO_FLUIDS_H_ */

#ifndef _NC_HIPERT_ITWO_FLUIDS_INLINE_H_
#define _NC_HIPERT_ITWO_FLUIDS_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE NcHIPertITwoFluidsEOM *
nc_hipert_itwo_fluids_eom_eval (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k)
{
  return NC_HIPERT_ITWO_FLUIDS_GET_IFACE (itf)->eom (itf, alpha, k);
}

NCM_INLINE NcHIPertITwoFluidsWKB *
nc_hipert_itwo_fluids_wkb_eval (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k)
{
  return NC_HIPERT_ITWO_FLUIDS_GET_IFACE (itf)->wkb (itf, alpha, k);
}

NCM_INLINE NcHIPertITwoFluidsTV *
nc_hipert_itwo_fluids_tv_eval (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k)
{
  return NC_HIPERT_ITWO_FLUIDS_GET_IFACE (itf)->tv (itf, alpha, k);
}

NCM_INLINE gdouble
nc_hipert_itwo_fluids_eval_unit (NcHIPertITwoFluids *itf)
{
  return NC_HIPERT_ITWO_FLUIDS_GET_IFACE (itf)->eval_unit (itf);
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_HIPERT_ITWO_FLUIDS_INLINE_H_ */

