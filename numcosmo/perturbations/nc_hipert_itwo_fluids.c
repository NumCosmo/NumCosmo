/***************************************************************************
 *            nc_hipert_itwo_fluids.c
 *
 *  Tue July 22 17:36:57 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hipert_itwo_fluids.c
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

/**
 * NcHIPertITwoFluids:
 *
 * Perturbation interface for two fluids system.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_hipert_itwo_fluids.h"

G_DEFINE_INTERFACE (NcHIPertITwoFluids, nc_hipert_itwo_fluids, G_TYPE_OBJECT)
G_DEFINE_BOXED_TYPE (NcHIPertITwoFluidsTV, nc_hipert_itwo_fluids_tv, nc_hipert_itwo_fluids_tv_dup, nc_hipert_itwo_fluids_tv_free)
G_DEFINE_BOXED_TYPE (NcHIPertITwoFluidsState, nc_hipert_itwo_fluids_state, nc_hipert_itwo_fluids_state_dup, nc_hipert_itwo_fluids_state_free)
G_DEFINE_BOXED_TYPE (NcHIPertITwoFluidsEOM, nc_hipert_itwo_fluids_eom, nc_hipert_itwo_fluids_eom_dup, nc_hipert_itwo_fluids_eom_free)
G_DEFINE_BOXED_TYPE (NcHIPertITwoFluidsWKB, nc_hipert_itwo_fluids_wkb, nc_hipert_itwo_fluids_wkb_dup, nc_hipert_itwo_fluids_wkb_free)

static void
nc_hipert_itwo_fluids_default_init (NcHIPertITwoFluidsInterface *iface)
{
  g_assert_cmpuint ((NC_HIPERT_ITWO_FLUIDS_VARS_LEN % 2), ==, 0);

  iface->eom       = NULL;
  iface->wkb       = NULL;
  iface->tv        = NULL;
  iface->eval_unit = NULL;
}

/**
 * nc_hipert_itwo_fluids_tv_dup:
 * @tf_tv: a #NcHIPertITwoFluidsTV
 *
 * Duplicates @tf_tv.
 *
 * Returns: (transfer full): a copy of @tf_tv.
 */
NcHIPertITwoFluidsTV *
nc_hipert_itwo_fluids_tv_dup (NcHIPertITwoFluidsTV *tf_tv)
{
  NcHIPertITwoFluidsTV *tf_tv_dup = g_new (NcHIPertITwoFluidsTV, 1);

  *tf_tv_dup = *tf_tv;

  return tf_tv_dup;
}

/**
 * nc_hipert_itwo_fluids_tv_free:
 * @tf_tv: a #NcHIPertITwoFluidsTV
 *
 * Frees @tf_tv.
 *
 */
void
nc_hipert_itwo_fluids_tv_free (NcHIPertITwoFluidsTV *tf_tv)
{
  g_free (tf_tv);
}

/**
 * nc_hipert_itwo_fluids_state_dup:
 * @tf_state: a #NcHIPertITwoFluidsState
 *
 * Duplicates @tf_state.
 *
 * Returns: (transfer full): a copy of @tf_state.
 */
NcHIPertITwoFluidsState *
nc_hipert_itwo_fluids_state_dup (NcHIPertITwoFluidsState *tf_state)
{
  NcHIPertITwoFluidsState *tf_state_dup = g_new (NcHIPertITwoFluidsState, 1);

  *tf_state_dup = *tf_state;

  return tf_state_dup;
}

/**
 * nc_hipert_itwo_fluids_state_free:
 * @tf_state: a #NcHIPertITwoFluidsState
 *
 * Frees @tf_state.
 *
 */
void
nc_hipert_itwo_fluids_state_free (NcHIPertITwoFluidsState *tf_state)
{
  g_free (tf_state);
}

/**
 * nc_hipert_itwo_fluids_eom_dup:
 * @tf_eom: a #NcHIPertITwoFluidsEOM
 *
 * Duplicates @tf_eom.
 *
 * Returns: (transfer full): a copy of @tf_eom.
 */
NcHIPertITwoFluidsEOM *
nc_hipert_itwo_fluids_eom_dup (NcHIPertITwoFluidsEOM *tf_eom)
{
  NcHIPertITwoFluidsEOM *tf_eom_dup = g_new (NcHIPertITwoFluidsEOM, 1);

  *tf_eom_dup = *tf_eom;

  return tf_eom_dup;
}

/**
 * nc_hipert_itwo_fluids_eom_free:
 * @tf_eom: a #NcHIPertITwoFluidsEOM
 *
 * Frees @tf_eom.
 *
 */
void
nc_hipert_itwo_fluids_eom_free (NcHIPertITwoFluidsEOM *tf_eom)
{
  g_free (tf_eom);
}

/**
 * nc_hipert_itwo_fluids_wkb_dup:
 * @tf_wkb: a #NcHIPertITwoFluidsWKB
 *
 * Duplicates @tf_wkb.
 *
 * Returns: (transfer full): a copy of @tf_wkb.
 */
NcHIPertITwoFluidsWKB *
nc_hipert_itwo_fluids_wkb_dup (NcHIPertITwoFluidsWKB *tf_wkb)
{
  NcHIPertITwoFluidsWKB *tf_wkb_dup = g_new (NcHIPertITwoFluidsWKB, 1);

  *tf_wkb_dup = *tf_wkb;

  return tf_wkb_dup;
}

/**
 * nc_hipert_itwo_fluids_wkb_free:
 * @tf_wkb: a #NcHIPertITwoFluidsWKB
 *
 * Frees @tf_wkb.
 *
 */
void
nc_hipert_itwo_fluids_wkb_free (NcHIPertITwoFluidsWKB *tf_wkb)
{
  g_free (tf_wkb);
}

/**
 * nc_hipert_itwo_fluids_eom_eval:
 * @itf: a #NcHIPertITwoFluids
 * @alpha: time in log of scale factor
 * @k: wave number
 *
 * Computes the coefficients of the differential equation for the
 * perturbations of the two fluids system.
 *
 *
 * Returns: (transfer none): a #NcHIPertITwoFluidsEOM.
 */

/**
 * nc_hipert_itwo_fluids_wkb_eval:
 * @itf: a #NcHIPertITwoFluids
 * @alpha: time in log of scale factor
 * @k: wave number
 *
 * Computes the WKB approximations of the perturbations of the
 * two fluids system.
 *
 * Returns: (transfer none): a #NcHIPertITwoFluidsWKB.
 */

typedef struct _NcHIPertITwoFluidsArgs
{
  const gdouble epsilon;
  const gdouble gw1;
  const gdouble gw2;
  const gdouble Fnu;
} NcHIPertITwoFluidsArgs;

/**
 * nc_hipert_itwo_fluids_tv_eval:
 * @itf: a #NcHIPertITwoFluids
 * @alpha: time in log of scale factor
 * @k: wave number
 *
 * Computes the transformation matrix between the perturbations of the
 * two fluids system and the variables used in the differential
 * equation.
 *
 * Returns: (transfer none): a #NcHIPertITwoFluidsTV.
 */

static complex double
_nc_hipert_itwo_fluids_state_eval_obs_helper (const complex double zeta, const complex double Q, const complex double Pzeta, const complex double PQ, NcHIPertITwoFluidsArgs *args, NcHIPertITwoFluidsObs obs)
{
  const gdouble epsilon = args->epsilon;
  const gdouble gw1     = args->gw1;
  const gdouble gw2     = args->gw2;
  const gdouble Fnu     = args->Fnu;

  switch (obs)
  {
    case NC_HIPERT_ITWO_FLUIDS_OBS_ZETA:
      return zeta;

      break;
    case NC_HIPERT_ITWO_FLUIDS_OBS_FKU_TOT:
      return Fnu * zeta;

      break;
    case NC_HIPERT_ITWO_FLUIDS_OBS_FKU_DIFF:
      return Fnu * (epsilon * (gw1 + gw2) / (gw1 * gw2) * Q);

      break;
    case NC_HIPERT_ITWO_FLUIDS_OBS_DELTA_TOT:
      return 3.0 * zeta - epsilon * Pzeta / (gw1 + gw2);

      break;
    case NC_HIPERT_ITWO_FLUIDS_OBS_DELTA_DIFF:
      return PQ;

      break;
    case NC_HIPERT_ITWO_FLUIDS_OBS_FKU_R:
      return Fnu * (zeta + epsilon * Q / gw1);

      break;
    case NC_HIPERT_ITWO_FLUIDS_OBS_FKU_W:
      return Fnu * (zeta - epsilon * Q / gw2);

      break;
    case NC_HIPERT_ITWO_FLUIDS_OBS_DELTA_R:
      return 3.0 * zeta + (gw2 * PQ - epsilon * Pzeta) / (gw1 + gw2);

      break;
    case NC_HIPERT_ITWO_FLUIDS_OBS_DELTA_W:
      return 3.0 * zeta - (gw1 * PQ + epsilon * Pzeta) / (gw1 + gw2);

      break;
    default:
      g_assert_not_reached ();
      break;
  }
}

/**
 * nc_hipert_itwo_fluids_state_eval_obs:
 * @tf_state: a #NcHIPertITwoFluidsState
 * @obs_a: a #NcHIPertITwoFluidsObs
 * @obs_b: a #NcHIPertITwoFluidsObs
 *
 * Computes the observable covariance between @obs_a and @obs_b.
 *
 * Returns: the value of the observable covariance.
 */
gdouble
nc_hipert_itwo_fluids_state_eval_obs (NcHIPertITwoFluidsState *tf_state, NcHIPertITwoFluidsObsMode obs_mode, NcHIPertITwoFluidsObs obs_a, NcHIPertITwoFluidsObs obs_b)
{
  const gdouble norma         = tf_state->norma * tf_state->norma * gsl_pow_3 (tf_state->k) / ncm_c_two_pi_2 ();
  const gdouble epsilon       = GSL_SIGN (tf_state->alpha);
  NcHIPertITwoFluidsArgs args = {epsilon, tf_state->gw1, tf_state->gw2, tf_state->Fnu};
  complex double obs1a        = 0.0;
  complex double obs1b        = 0.0;
  complex double obs2a        = 0.0;
  complex double obs2b        = 0.0;

  switch (obs_mode)
  {
    case NC_HIPERT_ITWO_FLUIDS_OBS_MODE_ONE:
      obs1a = _nc_hipert_itwo_fluids_state_eval_obs_helper (tf_state->zeta1, tf_state->Q1, tf_state->Pzeta1, tf_state->PQ1, &args, obs_a);
      obs1b = _nc_hipert_itwo_fluids_state_eval_obs_helper (tf_state->zeta1, tf_state->Q1, tf_state->Pzeta1, tf_state->PQ1, &args, obs_b);
      break;

    case NC_HIPERT_ITWO_FLUIDS_OBS_MODE_TWO:
      obs2a = _nc_hipert_itwo_fluids_state_eval_obs_helper (tf_state->zeta2, tf_state->Q2, tf_state->Pzeta2, tf_state->PQ2, &args, obs_a);
      obs2b = _nc_hipert_itwo_fluids_state_eval_obs_helper (tf_state->zeta2, tf_state->Q2, tf_state->Pzeta2, tf_state->PQ2, &args, obs_b);
      break;

    case NC_HIPERT_ITWO_FLUIDS_OBS_MODE_BOTH:
      obs1a = _nc_hipert_itwo_fluids_state_eval_obs_helper (tf_state->zeta1, tf_state->Q1, tf_state->Pzeta1, tf_state->PQ1, &args, obs_a);
      obs2a = _nc_hipert_itwo_fluids_state_eval_obs_helper (tf_state->zeta2, tf_state->Q2, tf_state->Pzeta2, tf_state->PQ2, &args, obs_a);
      obs1b = _nc_hipert_itwo_fluids_state_eval_obs_helper (tf_state->zeta1, tf_state->Q1, tf_state->Pzeta1, tf_state->PQ1, &args, obs_b);
      obs2b = _nc_hipert_itwo_fluids_state_eval_obs_helper (tf_state->zeta2, tf_state->Q2, tf_state->Pzeta2, tf_state->PQ2, &args, obs_b);
      break;

    default:
      g_assert_not_reached ();
      break;
  }

  return norma * creal (obs1a * conj (obs1b) + obs2a * conj (obs2b));
}

