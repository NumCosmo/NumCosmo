/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_csq1d.h
 *
 *  Mon September 09 13:56:11 2019
 *  Copyright  2019  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_csq1d.h
 * Copyright (C) 2019 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_CSQ1D_H_
#define _NCM_CSQ1D_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>

G_BEGIN_DECLS

#define NCM_TYPE_CSQ1D (ncm_csq1d_get_type ())
#define NCM_TYPE_CSQ1D_STATE (ncm_csq1d_state_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmCSQ1D, ncm_csq1d, NCM, CSQ1D, GObject)

GType ncm_csq1d_state_get_type (void) G_GNUC_CONST;

struct _NcmCSQ1DClass
{
  /*< private >*/
  GObjectClass parent_class;

  gdouble (*eval_xi)         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
  gdouble (*eval_nu)         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
  gdouble (*eval_nu2)        (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
  gdouble (*eval_m)          (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
  gdouble (*eval_int_1_m)    (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
  gdouble (*eval_int_mnu2)   (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
  gdouble (*eval_int_qmnu2)  (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
  gdouble (*eval_int_q2mnu2) (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
  gdouble (*eval_F1)         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
  gdouble (*eval_F2)         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
  gdouble (*eval_FN)         (NcmCSQ1D *csq1d, NcmModel *model, const gint n, const gdouble t);
  gdouble (*eval_powspec_factor) (NcmCSQ1D *csq1d, NcmModel *model);
  void (*prepare) (NcmCSQ1D *csq1d, NcmModel *model);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[3];
};

typedef struct _NcmCSQ1DSingFitUp NcmCSQ1DSingFitUp;
typedef struct _NcmCSQ1DSingFitUm NcmCSQ1DSingFitUm;

/**
 * NcmCSQ1DEvolState:
 * @NCM_CSQ1D_EVOL_STATE_INVALID: Invalid state
 * @NCM_CSQ1D_EVOL_STATE_ADIABATIC: Adiabatic state variables $(\alpha,\, \delta\gamma)$
 * @NCM_CSQ1D_EVOL_STATE_UP: $(\chi,\, U_+)$ state variables
 * @NCM_CSQ1D_EVOL_STATE_UM: $(\chi,\, U_-)$ state variables
 *
 * Variables describing the system evolution state. The state @NCM_CSQ1D_EVOL_STATE_ADIABATIC
 * is used to describe the adiabatic evolution of the system, this state use the frame
 * #NCM_CSQ1D_FRAME_ADIAB1 with the variables $(\alpha,\, \delta\gamma)$ to compute the evolution.
 * The state @NCM_CSQ1D_EVOL_STATE_UP and @NCM_CSQ1D_EVOL_STATE_UM are used to describe the
 * non-adiabatic evolution of the system, these states use the frame #NCM_CSQ1D_FRAME_ORIG with
 * the variables $(\chi,\, U_+)$ and $(\chi,\, U_-)$ to compute the evolution, respectively.
 *
 */
typedef enum _NcmCSQ1DEvolState
{
  NCM_CSQ1D_EVOL_STATE_INVALID = 0,
  NCM_CSQ1D_EVOL_STATE_ADIABATIC,
  NCM_CSQ1D_EVOL_STATE_UP,
  NCM_CSQ1D_EVOL_STATE_UM,
} NcmCSQ1DEvolState;

/**
 * NcmCSQ1DFrame:
 * @NCM_CSQ1D_FRAME_ORIG: Original frame
 * @NCM_CSQ1D_FRAME_ADIAB1: Adiabatic frame 1
 * @NCM_CSQ1D_FRAME_ADIAB2: Adiabatic frame 2
 * @NCM_CSQ1D_FRAME_NONADIAB1: Non-adiabatic frame 1
 * @NCM_CSQ1D_FRAME_NONADIAB2: Non-adiabatic frame 2
 *
 * Frames for the system.
 *
 */
typedef enum _NcmCSQ1DFrame
{
  NCM_CSQ1D_FRAME_ORIG = 0,
  NCM_CSQ1D_FRAME_ADIAB1,
  NCM_CSQ1D_FRAME_ADIAB2,
  NCM_CSQ1D_FRAME_NONADIAB1,
  NCM_CSQ1D_FRAME_NONADIAB2,
} NcmCSQ1DFrame;

/**
 * NcmCSQ1DState:
 *
 * Represents the state of the system.
 */
typedef struct _NcmCSQ1DState
{
  /*< private >*/
  NcmCSQ1DFrame frame;
  gdouble alpha;
  gdouble gamma;
  gdouble t;
} NcmCSQ1DState;

/* State related functions */

NcmCSQ1DState *ncm_csq1d_state_new (void);
NcmCSQ1DState *ncm_csq1d_state_copy (NcmCSQ1DState *state);
void ncm_csq1d_state_free (NcmCSQ1DState *state);

void ncm_csq1d_state_set_ag (NcmCSQ1DState *state, const NcmCSQ1DFrame frame, const gdouble t, const gdouble alpha, const gdouble gamma);
void ncm_csq1d_state_set_up (NcmCSQ1DState *state, const NcmCSQ1DFrame frame, const gdouble t, const gdouble chi, const gdouble Up);
void ncm_csq1d_state_set_um (NcmCSQ1DState *state, const NcmCSQ1DFrame frame, const gdouble t, const gdouble chi, const gdouble Um);

gdouble ncm_csq1d_state_get_time (NcmCSQ1DState *state);
NcmCSQ1DFrame ncm_csq1d_state_get_frame (NcmCSQ1DState *state);

void ncm_csq1d_state_get_ag (NcmCSQ1DState *state, gdouble *alpha, gdouble *gamma);
void ncm_csq1d_state_get_up (NcmCSQ1DState *state, gdouble *chi, gdouble *Up);
void ncm_csq1d_state_get_um (NcmCSQ1DState *state, gdouble *chi, gdouble *Um);
void ncm_csq1d_state_get_J (NcmCSQ1DState *state, gdouble *J11, gdouble *J12, gdouble *J22);

void ncm_csq1d_state_get_phi_Pphi (NcmCSQ1DState *state, gdouble *phi, gdouble *Pphi);
void ncm_csq1d_state_get_poincare_half_plane (NcmCSQ1DState *state, gdouble *x, gdouble *lny);
void ncm_csq1d_state_get_poincare_disc (NcmCSQ1DState *state, gdouble *x, gdouble *y);
void ncm_csq1d_state_get_minkowski (NcmCSQ1DState *state, gdouble *x1, gdouble *x2);
void ncm_csq1d_state_get_circle (NcmCSQ1DState *state, const gdouble r, const gdouble theta, NcmCSQ1DState *cstate);

/* CSQ1D methods */

NcmCSQ1D *ncm_csq1d_ref (NcmCSQ1D *csq1d);
void ncm_csq1d_free (NcmCSQ1D *csq1d);
void ncm_csq1d_clear (NcmCSQ1D **csq1d);

void ncm_csq1d_set_reltol (NcmCSQ1D *csq1d, const gdouble reltol);
void ncm_csq1d_set_abstol (NcmCSQ1D *csq1d, const gdouble abstol);
void ncm_csq1d_set_ti (NcmCSQ1D *csq1d, const gdouble ti);
void ncm_csq1d_set_tf (NcmCSQ1D *csq1d, const gdouble tf);
void ncm_csq1d_set_adiab_threshold (NcmCSQ1D *csq1d, const gdouble adiab_threshold);
void ncm_csq1d_set_prop_threshold (NcmCSQ1D *csq1d, const gdouble prop_threshold);
void ncm_csq1d_set_save_evol (NcmCSQ1D *csq1d, const gboolean save);
void ncm_csq1d_set_max_order_2 (NcmCSQ1D *csq1d, const gboolean truncate);
void ncm_csq1d_set_init_cond (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DEvolState evol_state, NcmCSQ1DState *initial_state);
void ncm_csq1d_set_init_cond_adiab (NcmCSQ1D *csq1d, NcmModel *model, const gdouble ti);

gdouble ncm_csq1d_get_reltol (NcmCSQ1D *csq1d);
gdouble ncm_csq1d_get_abstol (NcmCSQ1D *csq1d);
gdouble ncm_csq1d_get_ti (NcmCSQ1D *csq1d);
gdouble ncm_csq1d_get_tf (NcmCSQ1D *csq1d);
gdouble ncm_csq1d_get_adiab_threshold (NcmCSQ1D *csq1d);
gdouble ncm_csq1d_get_prop_threshold (NcmCSQ1D *csq1d);
gboolean ncm_csq1d_get_save_evol (NcmCSQ1D *csq1d);
gboolean ncm_csq1d_get_max_order_2 (NcmCSQ1D *csq1d);

NCM_INLINE gdouble ncm_csq1d_eval_xi         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
NCM_INLINE gdouble ncm_csq1d_eval_nu         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
NCM_INLINE gdouble ncm_csq1d_eval_nu2        (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
NCM_INLINE gdouble ncm_csq1d_eval_m          (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
NCM_INLINE gdouble ncm_csq1d_eval_int_1_m    (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
NCM_INLINE gdouble ncm_csq1d_eval_int_mnu2   (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
NCM_INLINE gdouble ncm_csq1d_eval_int_qmnu2  (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
NCM_INLINE gdouble ncm_csq1d_eval_int_q2mnu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
NCM_INLINE gdouble ncm_csq1d_eval_dm         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
NCM_INLINE gdouble ncm_csq1d_eval_F1         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
NCM_INLINE gdouble ncm_csq1d_eval_F2         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
NCM_INLINE gdouble ncm_csq1d_eval_FN         (NcmCSQ1D *csq1d, NcmModel *model, const gint n, const gdouble t);
NCM_INLINE gdouble ncm_csq1d_eval_powspec_factor (NcmCSQ1D *csq1d, NcmModel *model);

void ncm_csq1d_prepare (NcmCSQ1D *csq1d, NcmModel *model);

GArray *ncm_csq1d_get_time_array (NcmCSQ1D *csq1d, gdouble *smallest_t);

gboolean ncm_csq1d_find_adiab_time_limit (NcmCSQ1D *csq1d, NcmModel *model, gdouble t0, gdouble t1, const gdouble reltol, gdouble *ti);
gdouble ncm_csq1d_find_adiab_max (NcmCSQ1D *csq1d, NcmModel *model, gdouble t0, gdouble t1, const gdouble border_eps, gdouble *F1_min, gdouble *t_Bl, gdouble *t_Bu);

NcmCSQ1DState *ncm_csq1d_compute_adiab (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, NcmCSQ1DState *state, gdouble *alpha_reltol, gdouble *dgamma_reltol);
NcmCSQ1DState *ncm_csq1d_compute_adiab_frame (NcmCSQ1D *csq1d, NcmModel *model, const NcmCSQ1DFrame frame, const gdouble t, NcmCSQ1DState *state, gdouble *alpha_reltol, gdouble *dgamma_reltol);
NcmCSQ1DState *ncm_csq1d_compute_nonadiab (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, NcmCSQ1DState *state);
NcmCSQ1DState *ncm_csq1d_compute_H (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, NcmCSQ1DState *state);
NcmCSQ1DState *ncm_csq1d_eval_at (NcmCSQ1D *csq1d, const gdouble t, NcmCSQ1DState *state);
NcmCSQ1DState *ncm_csq1d_eval_at_frame (NcmCSQ1D *csq1d, NcmModel *model, const NcmCSQ1DFrame frame, const gdouble t, NcmCSQ1DState *state);

NcmCSQ1DState *ncm_csq1d_change_frame (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DState *state, NcmCSQ1DFrame frame);

void ncm_csq1d_prepare_prop (NcmCSQ1D *csq1d, NcmModel *model, const gdouble ti, const gdouble tii, const gdouble tf);

gdouble ncm_csq1d_get_tf_prop (NcmCSQ1D *csq1d);
NcmCSQ1DState *ncm_csq1d_compute_prop_vector (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, NcmCSQ1DState *state);

NcmCSQ1DState *ncm_csq1d_evolve_prop_vector (NcmCSQ1D *csq1d, NcmModel *model, NcmCSQ1DState *initial_state, const NcmCSQ1DFrame frame, const gdouble t, NcmCSQ1DState *state);

G_END_DECLS

#endif /* _NCM_CSQ1D_H_ */
#ifndef _NCM_CSQ1D_INLINE_H_
#define _NCM_CSQ1D_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE gdouble
ncm_csq1d_eval_xi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_xi (csq1d, model, t);
}

NCM_INLINE gdouble
ncm_csq1d_eval_nu (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_nu (csq1d, model, t);
}

NCM_INLINE gdouble
ncm_csq1d_eval_nu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_nu2 (csq1d, model, t);
}

NCM_INLINE gdouble
ncm_csq1d_eval_m (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_m (csq1d, model, t);
}

NCM_INLINE gdouble
ncm_csq1d_eval_int_1_m (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_int_1_m (csq1d, model, t);
}

NCM_INLINE gdouble
ncm_csq1d_eval_int_mnu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_int_mnu2 (csq1d, model, t);
}

NCM_INLINE gdouble
ncm_csq1d_eval_int_qmnu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_int_qmnu2 (csq1d, model, t);
}

NCM_INLINE gdouble
ncm_csq1d_eval_int_q2mnu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_int_q2mnu2 (csq1d, model, t);
}

NCM_INLINE gdouble
ncm_csq1d_eval_F1 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_F1 (csq1d, model, t);
}

NCM_INLINE gdouble
ncm_csq1d_eval_F2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_F2 (csq1d, model, t);
}

NCM_INLINE gdouble
ncm_csq1d_eval_FN (NcmCSQ1D *csq1d, NcmModel *model, const gint n, const gdouble t)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_FN (csq1d, model, n, t);
}

NCM_INLINE gdouble
ncm_csq1d_eval_powspec_factor (NcmCSQ1D *csq1d, NcmModel *model)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_powspec_factor (csq1d, model);
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_CSQ1D_INLINE_H_ */

