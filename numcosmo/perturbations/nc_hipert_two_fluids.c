/***************************************************************************
 *            nc_hipert_two_fluids.c
 *
 *  Tue June 10 19:14:57 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hipert_two_fluids.c
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
 * NcHIPertTwoFluids:
 *
 * Perturbation object for a two fluids system.
 *
 * This object provides the computation of the two fluid system of cosmological
 * perturbations. This problem is described by two fluids with energy density and
 * pressure given respectively by $\bar{\rho}_i$ and $\bar{p}_i$ for $i = 1,2$.
 *
 * The system is written in terms of the gauge invariant variable
 * $$\zeta \equiv \Psi - \frac{2\bar{K}}{\kappa(\bar{\rho} + \bar{p})} + E\mathcal{V},$$
 * and the entropy mode $$S = \frac{\kappa\varpi}{x^3 H}(\mathcal{U}_1 - \mathcal{U}_2),$$
 * where $\mathcal{U}_i \equiv \psi + E\mathcal{V}_i$ and
 * $$\varpi \equiv \frac{(\bar{\rho}_1+\bar{p}_1)(\bar{\rho}_2+\bar{p}_2)}{\bar{\rho}+\bar{p}}.$$
 *
 * Their momentum are
 * \begin{split}
 * P_\zeta &= \frac{2\bar{D}^2_\bar{K}\Psi}{x^3E}, \\\\
 * P_S &= \frac{\delta\rho_2}{\bar{\rho}_2+\bar{p}_2} - \frac{\delta\rho_1}{\bar{\rho}_1 + \bar{p}_1}.
 * \end{split}
 *
 * The equations of motion in their first order form are
 * \begin{align}
 * \zeta^\prime &= \frac{P_\zeta}{m_\zeta} + Y S, \\\\
 * P_\zeta^\prime &= -m_\zeta\mu_\zeta^2\zeta, \\\\
 * S^\prime &= \frac{P_S}{m_S} + Y \zeta, \\\\
 * P_S^\prime &= -m_S\mu_S^2S.
 * \end{align}
 * The mass $m_\zeta$ and the frequency $\mu_\zeta$ are defined by
 * \begin{align}
 * m_\zeta     &= \frac{3\Delta_\bar{K}(\bar{\rho} + \bar{p})}{\rho_\text{crit0} N x^3 c_s^2 E^2}, \\\\
 * \mu_\zeta^2 &= x^2N^2c_s^2k^2, \\\\
 * m_S         &= \frac{x^3}{c_m^2\varpi N}, \\\\
 * \mu_S^2     &= x^2N^2c_m^2k^2, \\\
 * Y           &= \frac{c_n^2}{c_s^2c_m^2}\frac{1}{m_\zeta m_S \Delta_\bar{K} N E}.
 * \end{align}
 * where $\bar{\rho} + \bar{p}$ is the background total energy density plus pressure,
 * $E^2 = H^2/H_0^2$ is the dimensionless Hubble function squared (nc_hicosmo_E2()), $c_s^2$ the speed of sound,
 * $N$ is the lapse function that in this case (using $\alpha$ as time variable) is $N \equiv \vert{}E\vert^{-1}$, $\rho_\text{crit0}$
 * is the critical density today defined by $\rho_\text{crit0} \equiv 3H_0^2/\kappa$ and $$\Delta_\bar{K} \equiv \frac{k^2}{k^2 + \Omega_{k0}}.$$
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_spline_func.h"
#include "perturbations/nc_hipert_two_fluids.h"
#include "perturbations/nc_hipert_itwo_fluids.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <cvode/cvode.h>
#include <cvode/cvode_ls.h>

#include <nvector/nvector_serial.h>
#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>

#include <gsl/gsl_roots.h>
#include <gsl/gsl_odeiv2.h>

#include <arkode/arkode.h>
#include <arkode/arkode_ls.h>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>
#include <arkode/arkode_mristep.h>
#include <arkode/arkode_sprk.h>
#include <arkode/arkode_sprkstep.h>

#define SUN_DENSE_ACCESS SM_ELEMENT_D
#define SUN_BAND_ACCESS SM_ELEMENT_D

#endif /* NUMCOSMO_GIR_SCAN */
#undef HAVE_SUNDIALS_ARKODE
#include "perturbations/nc_hipert_private.h"

typedef struct _NcHIPertTwoFluidsPrivate
{
  NcHIPertWKB *wkb_zeta;
  NcHIPertWKB *wkb_S;
  N_Vector abstol;
  gboolean useQP;
  NcmVector *state;
  gpointer arg;
  gpointer arkode;
  gdouble wkb_reltol;
  gdouble alpha_i;
  gdouble alpha_f;
} NcHIPertTwoFluidsPrivate;

struct _NcHIPertTwoFluids
{
  /*< private >*/
  NcHIPert parent_instance;
  NcHIPertTwoFluidsPrivate *priv;
};

enum
{
  PROP_0,
  PROP_WKB_RELTOL,
  PROP_ALPHA_I,
  PROP_ALPHA_F,
  PROP_SIZE,
};

/* *INDENT-OFF* */
G_DEFINE_QUARK (nc-hipert-two-fluids-error, nc_hipert_two_fluids_error)
/* *INDENT-ON* */

G_DEFINE_TYPE_WITH_PRIVATE (NcHIPertTwoFluids, nc_hipert_two_fluids, NC_TYPE_HIPERT)
G_DEFINE_BOXED_TYPE (NcHIPertTwoFluidsStateInterp, nc_hipert_two_fluids_state_interp, nc_hipert_two_fluids_state_interp_dup, nc_hipert_two_fluids_state_interp_free)

typedef struct _NcHIPertTwoFluidsArg
{
  NcHICosmo *cosmo;
  NcHIPertTwoFluids *ptf;
  gdouble prec;
  guint mode;
  gdouble alpha_i;
  gdouble alpha;
} NcHIPertTwoFluidsArg;

static void
nc_hipert_two_fluids_init (NcHIPertTwoFluids *ptf)
{
  NcHIPertTwoFluidsPrivate * const self = ptf->priv = nc_hipert_two_fluids_get_instance_private (ptf);

  self->wkb_zeta = NULL;
  self->wkb_S    = NULL;
  self->abstol   = NULL;
  self->useQP    = FALSE;
  self->state    = ncm_vector_new (NC_HIPERT_ITWO_FLUIDS_VARS_LEN);

  self->arg = g_new0 (NcHIPertTwoFluidsArg, 1);

  self->alpha_i    = 0.0;
  self->alpha_f    = 0.0;
  self->wkb_reltol = 0.0;

#ifdef HAVE_SUNDIALS_ARKODE
  self->arkode = NULL;
#endif /* HAVE_SUNDIALS_ARKODE */
}

static void
_nc_hipert_two_fluids_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHIPertTwoFluids *ptf = NC_HIPERT_TWO_FLUIDS (object);

  g_return_if_fail (NC_IS_HIPERT_TWO_FLUIDS (object));

  switch (prop_id)
  {
    case PROP_WKB_RELTOL:
      nc_hipert_two_fluids_set_wkb_reltol (ptf, g_value_get_double (value));
      break;
    case PROP_ALPHA_I:
      nc_hipert_two_fluids_set_initial_time (ptf, g_value_get_double (value));
      break;
    case PROP_ALPHA_F:
      nc_hipert_two_fluids_set_final_time (ptf, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_hipert_two_fluids_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /* NcHIPertTwoFluids *ptf = NC_HIPERT_TWO_FLUIDS (object); */
  g_return_if_fail (NC_IS_HIPERT_TWO_FLUIDS (object));

  switch (prop_id)
  {
    case PROP_WKB_RELTOL:
      g_value_set_double (value, nc_hipert_two_fluids_get_wkb_reltol (NC_HIPERT_TWO_FLUIDS (object)));
      break;
    case PROP_ALPHA_I:
      g_value_set_double (value, nc_hipert_two_fluids_get_initial_time (NC_HIPERT_TWO_FLUIDS (object)));
      break;
    case PROP_ALPHA_F:
      g_value_set_double (value, nc_hipert_two_fluids_get_final_time (NC_HIPERT_TWO_FLUIDS (object)));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_hipert_two_fluids_dispose (GObject *object)
{
  NcHIPertTwoFluids *ptf                = NC_HIPERT_TWO_FLUIDS (object);
  NcHIPertTwoFluidsPrivate * const self = ptf->priv;

  nc_hipert_wkb_clear (&self->wkb_zeta);
  nc_hipert_wkb_clear (&self->wkb_S);

  ncm_vector_clear (&self->state);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_two_fluids_parent_class)->dispose (object);
}

static void
_nc_hipert_two_fluids_finalize (GObject *object)
{
  NcHIPertTwoFluids *ptf                = NC_HIPERT_TWO_FLUIDS (object);
  NcHIPertTwoFluidsPrivate * const self = ptf->priv;

#ifdef HAVE_SUNDIALS_ARKODE

  if (self->arkode != NULL)
  {
    ARKodeFree (&self->arkode);
    self->arkode = NULL;
  }

#endif /* HAVE_SUNDIALS_ARKODE */

  g_free (self->arg);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_two_fluids_parent_class)->finalize (object);
}

static void _nc_hipert_two_fluids_set_mode_k (NcHIPert *pert, gdouble k);
static void _nc_hipert_two_fluids_set_abstol (NcHIPert *pert, gdouble abstol);
static void _nc_hipert_two_fluids_set_reltol (NcHIPert *pert, gdouble reltol);

static void
nc_hipert_two_fluids_class_init (NcHIPertTwoFluidsClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_hipert_two_fluids_set_property;
  object_class->get_property = &_nc_hipert_two_fluids_get_property;
  object_class->dispose      = &_nc_hipert_two_fluids_dispose;
  object_class->finalize     = &_nc_hipert_two_fluids_finalize;

  g_object_class_install_property (object_class,
                                   PROP_WKB_RELTOL,
                                   g_param_spec_double ("wkb-reltol",
                                                        NULL,
                                                        "WKB reltol",
                                                        GSL_DBL_EPSILON, 1.0e-1, 1.0e-2,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_ALPHA_I,
                                   g_param_spec_double ("alpha-i",
                                                        NULL,
                                                        "alpha_i",
                                                        -G_MAXDOUBLE, +G_MAXDOUBLE, -100.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_ALPHA_F,
                                   g_param_spec_double ("alpha-f",
                                                        NULL,
                                                        "alpha_f",
                                                        -G_MAXDOUBLE, +G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  NC_HIPERT_CLASS (klass)->set_mode_k = &_nc_hipert_two_fluids_set_mode_k;
  NC_HIPERT_CLASS (klass)->set_abstol = &_nc_hipert_two_fluids_set_abstol;
  NC_HIPERT_CLASS (klass)->set_reltol = &_nc_hipert_two_fluids_set_reltol;
}

static void
_nc_hipert_two_fluids_set_mode_k (NcHIPert *pert, gdouble k)
{
  NC_HIPERT_CLASS (nc_hipert_two_fluids_parent_class)->set_mode_k (pert, k);

  /* Chain up : start */
  if (!nc_hipert_prepared (pert))
  {
/*
 *   NcHIPertTwoFluids *ptf = NC_HIPERT_TWO_FLUIDS (pert);
 *   nc_hipert_set_mode_k (NC_HIPERT (self->wkb_zeta), k);
 *   nc_hipert_set_mode_k (NC_HIPERT (self->wkb_S), k);
 */
  }
}

static void
_nc_hipert_two_fluids_set_abstol (NcHIPert *pert, gdouble abstol)
{
  NC_HIPERT_CLASS (nc_hipert_two_fluids_parent_class)->set_abstol (pert, abstol);

  /* Chain up : start */
  if (!nc_hipert_prepared (pert))
  {
/*
 *   NcHIPertTwoFluids *ptf = NC_HIPERT_TWO_FLUIDS (pert);
 *   nc_hipert_set_abstol (NC_HIPERT (self->wkb_zeta), abstol);
 *   nc_hipert_set_abstol (NC_HIPERT (self->wkb_S), abstol);
 */
  }
}

static void
_nc_hipert_two_fluids_set_reltol (NcHIPert *pert, gdouble reltol)
{
  NC_HIPERT_CLASS (nc_hipert_two_fluids_parent_class)->set_reltol (pert, reltol);

  /* Chain up : start */
  if (!nc_hipert_prepared (pert))
  {
/*
 *   NcHIPertTwoFluids *ptf = NC_HIPERT_TWO_FLUIDS (pert);
 *   nc_hipert_set_reltol (NC_HIPERT (self->wkb_zeta), reltol);
 *   nc_hipert_set_reltol (NC_HIPERT (self->wkb_S), reltol);
 */
  }
}

/**
 * nc_hipert_two_fluids_state_interp_dup:
 * @sinterp: a #NcHIPertTwoFluidsStateInterp
 *
 * Creates a shallow copy of @sinterp. The internal splines for both modes are not
 * duplicated but shared via reference counting (i.e., the references are increased
 * using ncm_spline_ref()).
 *
 * This function is useful when multiple components need read-only access to the same
 * interpolation structure without duplicating memory-heavy spline data.
 *
 * Returns: (transfer full): a newly allocated #NcHIPertTwoFluidsStateInterp with shared
 * spline references.
 */
NcHIPertTwoFluidsStateInterp *
nc_hipert_two_fluids_state_interp_dup (NcHIPertTwoFluidsStateInterp *sinterp)
{
  NcHIPertTwoFluidsStateInterp *sinterp_dup = g_new0 (NcHIPertTwoFluidsStateInterp, 1);
  gint i;

  for (i = 0; i < NC_HIPERT_ITWO_FLUIDS_VARS_LEN; i++)
  {
    sinterp_dup->mode1_splines[i] = ncm_spline_ref (sinterp->mode1_splines[i]);
    sinterp_dup->mode2_splines[i] = ncm_spline_ref (sinterp->mode2_splines[i]);
  }

  sinterp_dup->state = sinterp->state;

  return sinterp_dup;
}

/**
 * nc_hipert_two_fluids_state_interp_free:
 * @sinterp: a #NcHIPertTwoFluidsStateInterp
 *
 * Frees the memory allocated for @sinterp and releases the references to all internally
 * held splines for both modes.
 *
 * This function must be called when the interpolation structure is no longer needed.
 */
void
nc_hipert_two_fluids_state_interp_free (NcHIPertTwoFluidsStateInterp *sinterp)
{
  gint i;

  for (i = 0; i < NC_HIPERT_ITWO_FLUIDS_VARS_LEN; i++)
  {
    ncm_spline_clear (&sinterp->mode1_splines[i]);
    ncm_spline_clear (&sinterp->mode2_splines[i]);
  }

  g_free (sinterp);
}

/**
 * nc_hipert_two_fluids_state_interp_eval:
 * @sinterp: a #NcHIPertTwoFluidsStateInterp
 * @cosmo: a #NcHICosmo
 * @x: the interpolation point (e.g., time or wavenumber)
 *
 * Interpolates the perturbation state at the point @x using the splines stored in
 * @sinterp. The result contains the complex values of the variables $\zeta$, $Q$,
 * $P_\zeta$, and $P_Q$ for both quantized modes.
 *
 * The interpolated values are stored internally in @sinterp->state and are valid until
 * the next call to this function.
 *
 * Returns: (transfer none): a pointer to the internally stored interpolated
 * #NcHIPertITwoFluidsState. The caller must not free the returned pointer.
 */
NcHIPertITwoFluidsState *
nc_hipert_two_fluids_state_interp_eval (NcHIPertTwoFluidsStateInterp *sinterp, NcHICosmo *cosmo, gdouble x)
{
  sinterp->state.zeta1 = ncm_spline_eval (sinterp->mode1_splines[NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R], x)
                         + I * ncm_spline_eval (sinterp->mode1_splines[NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I], x);
  sinterp->state.zeta2 = ncm_spline_eval (sinterp->mode2_splines[NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R], x)
                         + I * ncm_spline_eval (sinterp->mode2_splines[NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I], x);

  sinterp->state.Q1 = ncm_spline_eval (sinterp->mode1_splines[NC_HIPERT_ITWO_FLUIDS_VARS_S_R], x)
                      + I * ncm_spline_eval (sinterp->mode1_splines[NC_HIPERT_ITWO_FLUIDS_VARS_S_I], x);
  sinterp->state.Q2 = ncm_spline_eval (sinterp->mode2_splines[NC_HIPERT_ITWO_FLUIDS_VARS_S_R], x)
                      + I * ncm_spline_eval (sinterp->mode2_splines[NC_HIPERT_ITWO_FLUIDS_VARS_S_I], x);

  sinterp->state.Pzeta1 = ncm_spline_eval (sinterp->mode1_splines[NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R], x)
                          + I * ncm_spline_eval (sinterp->mode1_splines[NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I], x);
  sinterp->state.Pzeta2 = ncm_spline_eval (sinterp->mode2_splines[NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R], x)
                          + I * ncm_spline_eval (sinterp->mode2_splines[NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I], x);

  sinterp->state.PQ1 = ncm_spline_eval (sinterp->mode1_splines[NC_HIPERT_ITWO_FLUIDS_VARS_PS_R], x)
                       + I * ncm_spline_eval (sinterp->mode1_splines[NC_HIPERT_ITWO_FLUIDS_VARS_PS_I], x);
  sinterp->state.PQ2 = ncm_spline_eval (sinterp->mode2_splines[NC_HIPERT_ITWO_FLUIDS_VARS_PS_R], x)
                       + I * ncm_spline_eval (sinterp->mode2_splines[NC_HIPERT_ITWO_FLUIDS_VARS_PS_I], x);


  switch (sinterp->interp_mode)
  {
    case 1:
      sinterp->state.alpha = x;
      break;
    case 2:
      sinterp->state.k = x;
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  {
    NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (cosmo), sinterp->state.alpha, sinterp->state.k);

    sinterp->state.gw1 = eom->gw1;
    sinterp->state.gw2 = eom->gw2;
    sinterp->state.Fnu = eom->Fnu;
  }

  return &sinterp->state;
}

/**
 * nc_hipert_two_fluids_new:
 *
 * Creates a new #NcHIPertTwoFluids object.
 *
 * Returns: (transfer full): a new #NcHIPertTwoFluids.
 */
NcHIPertTwoFluids *
nc_hipert_two_fluids_new (void)
{
  NcHIPertTwoFluids *ptf = g_object_new (NC_TYPE_HIPERT_TWO_FLUIDS,
                                         "sys-size", NC_HIPERT_ITWO_FLUIDS_VARS_LEN,
                                         NULL);

  return ptf;
}

/**
 * nc_hipert_two_fluids_ref:
 * @ptf: a #NcHIPertTwoFluids
 *
 * Increases the reference count of @ptf.
 *
 * Returns: (transfer full): @ptf.
 */
NcHIPertTwoFluids *
nc_hipert_two_fluids_ref (NcHIPertTwoFluids *ptf)
{
  return g_object_ref (ptf);
}

/**
 * nc_hipert_two_fluids_free:
 * @ptf: a #NcHIPertTwoFluids
 *
 * Decreases the reference count of @ptf.
 *
 */
void
nc_hipert_two_fluids_free (NcHIPertTwoFluids *ptf)
{
  g_object_unref (ptf);
}

/**
 * nc_hipert_two_fluids_clear:
 * @ptf: a #NcHIPertTwoFluids
 *
 * Decreases the reference count of *@ptf and sets *@ptf to NULL.
 *
 */
void
nc_hipert_two_fluids_clear (NcHIPertTwoFluids **ptf)
{
  g_clear_object (ptf);
}

/**
 * nc_hipert_two_fluids_set_wkb_reltol:
 * @ptf: a #NcHIPertTwoFluids
 * @reltol: Relative tolerance
 *
 * Sets the relative tolerance used to determine the validity range of the WKB
 * approximation. Beyond this threshold, the numerical solver is activated.
 *
 */
void
nc_hipert_two_fluids_set_wkb_reltol (NcHIPertTwoFluids *ptf, gdouble reltol)
{
  NcHIPertTwoFluidsPrivate * const self = ptf->priv;

  self->wkb_reltol = reltol;
}

/**
 * nc_hipert_two_fluids_set_initial_time:
 * @ptf: a #NcHIPertTwoFluids
 * @alpha_i: initial log-redshift time
 *
 * Sets the initial log-redshift time.
 *
 */
void
nc_hipert_two_fluids_set_initial_time (NcHIPertTwoFluids *ptf, gdouble alpha_i)
{
  NcHIPertTwoFluidsPrivate * const self = ptf->priv;

  self->alpha_i = alpha_i;
}

/**
 * nc_hipert_two_fluids_set_final_time:
 * @ptf: a #NcHIPertTwoFluids
 * @alpha_f: final log-redshift time
 *
 * Sets the final log-redshift time.
 *
 */
void
nc_hipert_two_fluids_set_final_time (NcHIPertTwoFluids *ptf, gdouble alpha_f)
{
  NcHIPertTwoFluidsPrivate * const self = ptf->priv;

  self->alpha_f = alpha_f;
}

/**
 * nc_hipert_two_fluids_get_wkb_reltol:
 * @ptf: a #NcHIPertTwoFluids
 *
 * Retrieves the relative tolerance that defines the transition point between the WKB
 * approximation and the numerical solver.
 *
 * Returns: the current relative tolerance.
 */
gdouble
nc_hipert_two_fluids_get_wkb_reltol (NcHIPertTwoFluids *ptf)
{
  NcHIPertTwoFluidsPrivate * const self = ptf->priv;

  return self->wkb_reltol;
}

/**
 * nc_hipert_two_fluids_get_initial_time:
 * @ptf: a #NcHIPertTwoFluids
 *
 * Retrieves the initial log-redshift time.
 *
 * Returns: the current initial log-redshift time.
 */
gdouble
nc_hipert_two_fluids_get_initial_time (NcHIPertTwoFluids *ptf)
{
  NcHIPertTwoFluidsPrivate * const self = ptf->priv;

  return self->alpha_i;
}

/**
 * nc_hipert_two_fluids_get_final_time:
 * @ptf: a #NcHIPertTwoFluids
 *
 * Retrieves the final log-redshift time.
 *
 * Returns: the current final log-redshift time.
 */
gdouble
nc_hipert_two_fluids_get_final_time (NcHIPertTwoFluids *ptf)
{
  NcHIPertTwoFluidsPrivate * const self = ptf->priv;

  return self->alpha_f;
}

/**
 * nc_hipert_two_fluids_eom:
 * @ptf: a #NcHIPertTwoFluids
 * @cosmo: a #NcHICosmo
 * @alpha: log-redshift time
 * @eom: (out callee-allocates) (transfer none): Equation of motion variables
 *
 * Calculates the equation of motion coefficients for the $(\zeta,\,S)$ system.
 *
 */
void
nc_hipert_two_fluids_eom (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, NcHIPertITwoFluidsEOM **eom)
{
  NcHIPert *pert = NC_HIPERT (ptf);

  *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, nc_hipert_get_mode_k (pert));

  return;
}

/**
 * nc_hipert_two_fluids_wkb:
 * @ptf: a #NcHIPertTwoFluids
 * @cosmo: a #NcHICosmo
 * @alpha: the log-redshift time
 * @wkb: (out callee-allocates) (transfer none): WKB variables
 *
 * Calculates the WKB approximation for the $(\zeta,\,S)$ system.
 *
 */
void
nc_hipert_two_fluids_wkb (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, NcHIPertITwoFluidsWKB **wkb)
{
  *wkb = nc_hipert_itwo_fluids_wkb_eval (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, nc_hipert_get_mode_k (NC_HIPERT (ptf)));
}

/**
 * nc_hipert_two_fluids_get_init_cond_QP:
 * @ptf: a #NcHIPertTwoFluids
 * @cosmo: a #NcHICosmo
 * @alpha: the log-redshift time
 * @main_mode: main mode
 * @beta_R: mode $R$ initial phase
 * @init_cond: a #NcmVector (size >= 8) where to put the initial conditions
 *
 * Calculates the initial condition for the $(Q,\,P)$ system with initial phase for the
 * R solution $\beta_R = $  @beta_R. The variable @main_mode chooses which mode is
 * excited (1 or 2).
 *
 */
void
nc_hipert_two_fluids_get_init_cond_QP (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, guint main_mode, const gdouble beta_R, NcmVector *init_cond)
{
  NcHIPert *pert             = NC_HIPERT (ptf);
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, nc_hipert_get_mode_k (pert));
  const gdouble beta_I       = beta_R + 0.5 * M_PI;

  switch (main_mode)
  {
    case 1:
    {
      const gdouble dsigma1     = eom->nu1 - 0.5 * eom->gammabar11;
      const gdouble dsigma2     = eom->nu2 - 0.5 * eom->gammabar22;
      const gdouble Deltadsigma = dsigma1 - dsigma2;

      const complex double a_R_1  = cexp (+I * beta_R) / M_SQRT2;
      const complex double ac_R_1 = cexp (-I * beta_R) / M_SQRT2;
      const complex double a_I_1  = cexp (+I * beta_I) / M_SQRT2;
      const complex double ac_I_1 = cexp (-I * beta_I) / M_SQRT2;

      const complex double A_R11 = a_R_1 - 0.25 * eom->gammabar11 * ac_R_1 / dsigma1;
      const complex double A_R12 = 0.5 * (a_R_1 * (-eom->gammabar12 + I * eom->taubar) - ac_R_1 * eom->gammabar12) / Deltadsigma;

      const complex double A_I11 = a_I_1 - 0.25 * eom->gammabar11 * ac_I_1 / dsigma1;
      const complex double A_I12 = 0.5 * (a_I_1 * (-eom->gammabar12 + I * eom->taubar) - ac_I_1 * eom->gammabar12) / Deltadsigma;

      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1, NC_HIPERT_TWO_FLUIDS_A2Q (A_R11));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1, NC_HIPERT_TWO_FLUIDS_A2P (A_R11));

      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2, NC_HIPERT_TWO_FLUIDS_A2Q (A_R12));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2, NC_HIPERT_TWO_FLUIDS_A2P (A_R12));

      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1, NC_HIPERT_TWO_FLUIDS_A2Q (A_I11));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1, NC_HIPERT_TWO_FLUIDS_A2P (A_I11));

      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2, NC_HIPERT_TWO_FLUIDS_A2Q (A_I12));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2, NC_HIPERT_TWO_FLUIDS_A2P (A_I12));

      break;
    }
    case 2:
    {
      const gdouble dsigma1     = eom->nu1 - 0.5 * eom->gammabar11;
      const gdouble dsigma2     = eom->nu2 - 0.5 * eom->gammabar22;
      const gdouble Deltadsigma = dsigma1 - dsigma2;

      const complex double a_R_2  = cexp (+I * beta_R) / M_SQRT2;
      const complex double ac_R_2 = cexp (-I * beta_R) / M_SQRT2;
      const complex double a_I_2  = cexp (+I * beta_I) / M_SQRT2;
      const complex double ac_I_2 = cexp (-I * beta_I) / M_SQRT2;

      const complex double A_R22 = a_R_2 - 0.25 * eom->gammabar22 * ac_R_2 / dsigma2;
      const complex double A_R21 = 0.5 * (a_R_2 * (I * eom->taubar + eom->gammabar12) - ac_R_2 * eom->gammabar12) / Deltadsigma;

      const complex double A_I22 = a_I_2 - 0.25 * eom->gammabar22 * ac_I_2 / dsigma2;
      const complex double A_I21 = 0.5 * (a_I_2 * (I * eom->taubar + eom->gammabar12) - ac_I_2 * eom->gammabar12) / Deltadsigma;

      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1, NC_HIPERT_TWO_FLUIDS_A2Q (A_R21));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1, NC_HIPERT_TWO_FLUIDS_A2P (A_R21));

      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2, NC_HIPERT_TWO_FLUIDS_A2Q (A_R22));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2, NC_HIPERT_TWO_FLUIDS_A2P (A_R22));

      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1, NC_HIPERT_TWO_FLUIDS_A2Q (A_I21));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1, NC_HIPERT_TWO_FLUIDS_A2P (A_I21));

      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2, NC_HIPERT_TWO_FLUIDS_A2Q (A_I22));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2, NC_HIPERT_TWO_FLUIDS_A2P (A_I22));

      break;
    }
    default:
      g_error ("nc_hipert_two_fluids_get_init_cond: Unknown main mode %u.", main_mode);
      break;
  }
}

/**
 * nc_hipert_two_fluids_get_init_cond_zetaS:
 * @ptf: a #NcHIPertTwoFluids
 * @cosmo: a #NcHICosmo
 * @alpha: the log-redshift time
 * @main_mode: main mode
 * @beta_R: mode $R$ initial phase
 * @init_cond: a #NcmVector (size >= 8) where to put the initial conditions
 *
 * Calculates the initial condition for the $\zeta{}S$ system with initial phase for the
 * R solution $\beta_R = $ @beta_R. The variable @main_mode chooses which mode is
 * excited (1 or 2).
 *
 */
void
nc_hipert_two_fluids_get_init_cond_zetaS (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, guint main_mode, const gdouble beta_R, NcmVector *init_cond)
{
  NcHIPert *pert                = NC_HIPERT (ptf);
  NcHIPertITwoFluidsWKB *tf_wkb = nc_hipert_itwo_fluids_wkb_eval (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, nc_hipert_get_mode_k (pert));

  switch (main_mode)
  {
    case 1:
    {
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R, creal (tf_wkb->state.zeta1));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I, cimag (tf_wkb->state.zeta1));

      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_S_R, creal (tf_wkb->state.Q1));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_S_I, cimag (tf_wkb->state.Q1));

      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R, creal (tf_wkb->state.Pzeta1));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I, cimag (tf_wkb->state.Pzeta1));

      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R, creal (tf_wkb->state.PQ1));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I, cimag (tf_wkb->state.PQ1));
      break;
    }
    case 2:
    {
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R, creal (tf_wkb->state.zeta2));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I, cimag (tf_wkb->state.zeta2));

      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_S_R, creal (tf_wkb->state.Q2));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_S_I, cimag (tf_wkb->state.Q2));

      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R, creal (tf_wkb->state.Pzeta2));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I, cimag (tf_wkb->state.Pzeta2));

      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R, creal (tf_wkb->state.PQ2));
      ncm_vector_set (init_cond, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I, cimag (tf_wkb->state.PQ2));
      break;
    }
    default:
      g_error ("nc_hipert_two_fluids_get_init_cond_zetaS: Unknown main mode %u.", main_mode);
      break;
  }
}

/**
 * nc_hipert_two_fluids_to_zeta_s:
 * @ptf: a #NcHIPertTwoFluids
 * @cosmo: a #NcHICosmo
 * @alpha: the log-redshift time
 * @state: a #NcmVector (size >= 8) current state in $(Q,\,P)$ variables
 *
 * Transform in-place the variables @init_cond from $(Q,\,P)$ to $(\zeta,\,S)$, assuming
 * they are calculated at $\alpha$ = @alpha.
 *
 */
void
nc_hipert_two_fluids_to_zeta_s (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, NcmVector *state)
{
  NcHIPert *pert                                 = NC_HIPERT (ptf);
  NcHIPertITwoFluidsTV *tv                       = nc_hipert_itwo_fluids_tv_eval (NC_HIPERT_ITWO_FLUIDS (cosmo), alpha, nc_hipert_get_mode_k (pert));
  const guint syssize                            = NC_HIPERT_ITWO_FLUIDS_VARS_LEN / 2;
  gdouble zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_LEN] = {
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0
  };
  guint i;

  for (i = 0; i < syssize; i++)
  {
    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R]  += tv->zeta[i]  * ncm_vector_get (state, i);
    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_S_R]     += tv->s[i]     * ncm_vector_get (state, i);
    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R] += tv->Pzeta[i] * ncm_vector_get (state, i);
    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_PS_R]    += tv->Ps[i]    * ncm_vector_get (state, i);

    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I]  += tv->zeta[i]  * ncm_vector_get (state, syssize + i);
    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_S_I]     += tv->s[i]     * ncm_vector_get (state, syssize + i);
    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I] += tv->Pzeta[i] * ncm_vector_get (state, syssize + i);
    zeta_s[NC_HIPERT_ITWO_FLUIDS_VARS_PS_I]    += tv->Ps[i]    * ncm_vector_get (state, syssize + i);
  }

  for (i = 0; i < NC_HIPERT_ITWO_FLUIDS_VARS_LEN; i++)
    ncm_vector_set (state, i, zeta_s[i]);
}

static gint
_nc_hipert_two_fluids_f_QP (sunrealtype alpha, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHIPertTwoFluidsArg *arg  = (NcHIPertTwoFluidsArg *) f_data;
  const gdouble k            = nc_hipert_get_mode_k (NC_HIPERT (arg->ptf));
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);

  const gdouble Q_R1 = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1);
  const gdouble P_R1 = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1);
  const gdouble Q_R2 = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2);
  const gdouble P_R2 = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2);

  const gdouble Q_I1 = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1);
  const gdouble P_I1 = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1);
  const gdouble Q_I2 = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2);
  const gdouble P_I2 = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2);

  const complex double A_R1 = NC_HIPERT_TWO_FLUIDS_QP2A (Q_R1, P_R1);
  const complex double A_R2 = NC_HIPERT_TWO_FLUIDS_QP2A (Q_R2, P_R2);

  const complex double A_I1 = NC_HIPERT_TWO_FLUIDS_QP2A (Q_I1, P_I1);
  const complex double A_I2 = NC_HIPERT_TWO_FLUIDS_QP2A (Q_I2, P_I2);

  const gdouble gammabar11 = eom->gammabar11;
  const gdouble gammabar22 = eom->gammabar22;
  const gdouble gammabar12 = eom->gammabar12;
  const gdouble taubar12   = eom->taubar;

  const complex double dA_R1 = I * eom->nu1 * A_R1 + gammabar11 * Q_R1 + gammabar12 * Q_R2 + 0.5 * taubar12 * A_R2;
  const complex double dA_R2 = I * eom->nu2 * A_R2 + gammabar22 * Q_R2 + gammabar12 * Q_R1 - 0.5 * taubar12 * A_R1;

  const complex double dA_I1 = I * eom->nu1 * A_I1 + gammabar11 * Q_I1 + gammabar12 * Q_I2 + 0.5 * taubar12 * A_I2;
  const complex double dA_I2 = I * eom->nu2 * A_I2 + gammabar22 * Q_I2 + gammabar12 * Q_I1 - 0.5 * taubar12 * A_I1;

  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1) = NC_HIPERT_TWO_FLUIDS_A2Q (dA_R1);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1) = NC_HIPERT_TWO_FLUIDS_A2P (dA_R1);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2) = NC_HIPERT_TWO_FLUIDS_A2Q (dA_R2);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2) = NC_HIPERT_TWO_FLUIDS_A2P (dA_R2);

  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1) = NC_HIPERT_TWO_FLUIDS_A2Q (dA_I1);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1) = NC_HIPERT_TWO_FLUIDS_A2P (dA_I1);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2) = NC_HIPERT_TWO_FLUIDS_A2Q (dA_I2);
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2) = NC_HIPERT_TWO_FLUIDS_A2P (dA_I2);

  return 0;
}

static gint
_nc_hipert_two_fluids_f_zetaS (sunrealtype alpha, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHIPertTwoFluidsArg *arg  = (NcHIPertTwoFluidsArg *) f_data;
  const gdouble k            = nc_hipert_get_mode_k (NC_HIPERT (arg->ptf));
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);

  const gdouble zeta_R  = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R);
  const gdouble S_R     = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_S_R);
  const gdouble Pzeta_R = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R);
  const gdouble PS_R    = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R);

  const gdouble zeta_I  = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I);
  const gdouble S_I     = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_S_I);
  const gdouble Pzeta_I = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I);
  const gdouble PS_I    = NV_Ith_S (y, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I);

  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = Pzeta_R / eom->m_zeta + eom->y * PS_R;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = PS_R    / eom->m_s    + eom->y * Pzeta_R;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = -eom->mnu2_zeta * zeta_R;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = -eom->mnu2_s    * S_R;

  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = Pzeta_I / eom->m_zeta + eom->y * PS_I;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = PS_I    / eom->m_s    + eom->y * Pzeta_I;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = -eom->mnu2_zeta * zeta_I;
  NV_Ith_S (ydot, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = -eom->mnu2_s    * S_I;

#ifdef TWO_FLUIDS_DEBUG_WKB
  {
    NcHIPertITwoFluidsWKB *wkb = nc_hipert_itwo_fluids_wkb_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);
    const gdouble Avzeta1      = cabs (wkb->zeta1);
    const gdouble AvQ1         = cabs (wkb->Q1);
    const gdouble AvPzeta1     = cabs (wkb->Pzeta1);
    const gdouble AvPQ1        = cabs (wkb->PQ1);

    const gdouble Avzeta2  = cabs (wkb->zeta2);
    const gdouble AvQ2     = cabs (wkb->Q2);
    const gdouble AvPzeta2 = cabs (wkb->Pzeta2);
    const gdouble AvPQ2    = cabs (wkb->PQ2);

    const gdouble Azeta  = hypot (zeta_R, zeta_I);
    const gdouble APzeta = hypot (Pzeta_R, Pzeta_I);
    const gdouble AS     = hypot (S_R, S_I);
    const gdouble APS    = hypot (PS_R, PS_I);

    printf ("% .4f % 18.11e % 18.11e % 18.11e % 18.11e, % 18.11e % 18.11e % 18.11e % 18.11e % 18.11e % 18.11e\n",
            alpha,
            wkb->mode1_err,
            wkb->mode2_err,
            Azeta / Avzeta1,
            AS / AvQ1,
            APzeta / AvPzeta1,
            APS / AvPQ1,
            Azeta / Avzeta2,
            AS / AvQ2,
            APzeta / AvPzeta2,
            APS / AvPQ2
           );
    fflush (stdout);
    fflush (stderr);
  }
#endif /* DEBUG_WKB */

  return 0;
}

static gint
_nc_hipert_two_fluids_J_zetaS (sunrealtype alpha, N_Vector y, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHIPertTwoFluidsArg *arg  = (NcHIPertTwoFluidsArg *) jac_data;
  const gdouble k            = nc_hipert_get_mode_k (NC_HIPERT (arg->ptf));
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);

  /************************************************************************************************************/
  /* Solution R */
  /************************************************************************************************************/

  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R,  NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R,  NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R,  NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = 1.0 / eom->m_zeta;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R,  NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = eom->y;

  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_R,     NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_R,     NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_R,     NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = eom->y;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_R,     NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = 1.0 / eom->m_s;

  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = -eom->mnu2_zeta;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R, NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = 0.0;

  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R,    NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R)  = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R,    NC_HIPERT_ITWO_FLUIDS_VARS_S_R)     = -eom->mnu2_s;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R,    NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R) = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_R,    NC_HIPERT_ITWO_FLUIDS_VARS_PS_R)    = 0.0;

  /************************************************************************************************************/
  /* Solution I */
  /************************************************************************************************************/

  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I,  NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I,  NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I,  NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = 1.0 / eom->m_zeta;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I,  NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = eom->y;

  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_I,     NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_I,     NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_I,     NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = eom->y;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_S_I,     NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = 1.0 / eom->m_s;

  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = -eom->mnu2_zeta;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I, NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = 0.0;

  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I,    NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I)  = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I,    NC_HIPERT_ITWO_FLUIDS_VARS_S_I)     = -eom->mnu2_s;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I,    NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I) = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_PS_I,    NC_HIPERT_ITWO_FLUIDS_VARS_PS_I)    = 0.0;

  return 0;
}

#ifdef HAVE_SUNDIALS_ARKODE

static gint
_nc_hipert_two_fluids_J_QP (sunrealtype alpha, N_Vector y, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHIPertTwoFluidsArg *arg  = (NcHIPertTwoFluidsArg *) jac_data;
  const gdouble k            = nc_hipert_get_mode_k (NC_HIPERT (arg->ptf));
  NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);

  const gdouble gammabar11 = eom->gammabar11;
  const gdouble gammabar22 = eom->gammabar22;
  const gdouble gammabar12 = eom->gammabar12;
  const gdouble taubar12   = eom->taubar;

  /************************************************************************************************************/
  /* Solution R */
  /************************************************************************************************************/

  /* eom->nu1 * P1 +  0.5 * taubar12 * Q2 */

  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1) = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R1) = eom->nu1;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2) = 0.5 * taubar12;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R2) = 0.0;

  /* - eom->nu1 * Q1 + gammabar11 * Q1 + gammabar12 * Q2 + 0.5 * taubar12 * P2 */

  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1) = -eom->nu1 + gammabar11;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R1) = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2) = gammabar12;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R2) = 0.5 * taubar12;

  /* eom->nu2 * P2 - 0.5 * taubar12 * Q1 */

  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1) = -0.5 * taubar12;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R1) = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2) = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R2) = eom->nu2;

  /* - eom->nu2 * Q2 + gammabar22 * Q2 + gammabar12 * Q1 - 0.5 * taubar12 * P1 */

  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1) = gammabar12;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R1) = -0.5 * taubar12;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2) = -eom->nu2 + gammabar22;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_R2) = 0.0;

  /************************************************************************************************************/
  /* Solution I */
  /************************************************************************************************************/

  /* eom->nu1 * P1 +  0.5 * taubar12 * Q2 */

  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1) = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I1) = eom->nu1;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2) = 0.5 * taubar12;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I2) = 0.0;

/* - eom->nu1 * Q1 + gammabar11 * Q1 + gammabar12 * Q2 + 0.5 * taubar12 * P2 */

  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1) = -eom->nu1 + gammabar11;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I1) = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2) = gammabar12;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I2) = 0.5 * taubar12;

/* eom->nu2 * P2 - 0.5 * taubar12 * Q1 */

  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1) = -0.5 * taubar12;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I1) = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2) = 0.0;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I2) = eom->nu2;

/* - eom->nu2 * Q2 + gammabar22 * Q2 + gammabar12 * Q1 - 0.5 * taubar12 * P1 */

  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1) = gammabar12;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I1) = -0.5 * taubar12;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2) = -eom->nu2 + gammabar22;
  SUN_DENSE_ACCESS (J, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2,  NC_HIPERT_ITWO_FLUIDS_VARS_P_I2) = 0.0;

  return 0;
}

#endif /* HAVE_SUNDIALS_ARKODE */

/**
 * nc_hipert_two_fluids_set_init_cond:
 * @ptf: a #NcHIPertTwoFluids
 * @cosmo: a #NcHICosmo
 * @alpha: the log-redshift time
 * @main_mode: main mode
 * @useQP: whether to use the $(Q,\,P)$ system
 * @init_cond: a #NcmVector (size >= 8) containing the initial conditions
 *
 * Sets the initial conditions for the two fluids system evolution.
 *
 */
void
nc_hipert_two_fluids_set_init_cond (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alpha, guint main_mode, gboolean useQP, NcmVector *init_cond)
{
  NcHIPertTwoFluidsPrivate * const self = ptf->priv;
  NcHIPert *pert                        = NC_HIPERT (ptf);
  NcHIPertPrivate * const pself         = nc_hipert_get_private (pert);
  gint vtype                            = useQP ? 1 : 0;
  gint c_vtype                          = self->useQP ? 1 : 0;
  gint flag;
  guint i;

  if (vtype != c_vtype)
    nc_hipert_reset_solver (pert);

  nc_hipert_set_sys_size (pert, NC_HIPERT_ITWO_FLUIDS_VARS_LEN);

  pself->alpha0 = alpha;

  for (i = 0; i < NC_HIPERT_ITWO_FLUIDS_VARS_LEN; i++)
  {
    NV_Ith_S (pself->y, i) = ncm_vector_get (init_cond, i);
  }

  if (!pself->cvode_init)
  {
#ifdef HAVE_SUNDIALS_ARKODE
    self->arkode = ERKStepCreate (_nc_hipert_two_fluids_f_zetaS, alpha, pself->y, pself->sunctx);
    NCM_CVODE_CHECK (self->arkode, "ERKStepCreate", 0, );
#endif /* HAVE_SUNDIALS_ARKODE */

    if (useQP)
    {
      flag = CVodeInit (pself->cvode, &_nc_hipert_two_fluids_f_QP, alpha, pself->y);
      NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

      pself->cvode_init = TRUE;
    }
    else
    {
      flag = CVodeInit (pself->cvode, &_nc_hipert_two_fluids_f_zetaS, alpha, pself->y);
      NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

      pself->cvode_init = TRUE;
    }
  }
  else
  {
    flag = CVodeReInit (pself->cvode, alpha, pself->y);
    NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );

#ifdef HAVE_SUNDIALS_ARKODE
    flag = ERKStepReInit (self->arkode, _nc_hipert_two_fluids_f_zetaS, alpha, pself->y);
    NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );
#endif /* HAVE_SUNDIALS_ARKODE */
  }

  flag = CVodeSStolerances (pself->cvode, pself->reltol, 0.0);
  NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );

  flag = CVodeSetMaxNumSteps (pself->cvode, G_MAXUINT32);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

  flag = CVodeSetLinearSolver (pself->cvode, pself->LS, pself->A);
  NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, );

  flag = CVodeSetJacFn (pself->cvode, &_nc_hipert_two_fluids_J_zetaS);
  NCM_CVODE_CHECK (&flag, "CVodeSetJacFn", 1, );

  flag = CVodeSetMaxHnilWarns (pself->cvode, -1);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxHnilWarns", 1, );

#ifdef HAVE_SUNDIALS_ARKODE

  flag = ARKodeSStolerances (self->arkode, pself->reltol, 0.0);
  NCM_CVODE_CHECK (&flag, "ARKodeSStolerances", 1, );

  flag = ARKodeSetInitStep (self->arkode, 1.0e-11);
  NCM_CVODE_CHECK (&flag, "ARKodeSetInitStep", 1, );

  flag = ARKodeSetMaxNumSteps (self->arkode, G_MAXUINT32);
  NCM_CVODE_CHECK (&flag, "ARKodeSetMaxNumSteps", 1, );

  flag = ARKodeSetAdaptController (self->arkode, NULL);
  NCM_CVODE_CHECK (&flag, "ARKodeSetAdaptController", 1, );

  /* flag = ARKodeSetLinear (self->arkode, 1); */
  /* NCM_CVODE_CHECK (&flag, "ARKodeSetLinear", 1, ); */

  if (useQP)
  {
    flag = ARKodeSetJacFn (self->arkode, &_nc_hipert_two_fluids_J_QP);
    NCM_CVODE_CHECK (&flag, "ARKodeSetJacFn", 1, );
  }

#endif /* HAVE_SUNDIALS_ARKODE */

  self->useQP = useQP;
}

/**
 * nc_hipert_two_fluids_evolve:
 * @ptf: a #NcHIPertTwoFluids
 * @cosmo: a #NcHICosmo
 * @alphaf: the final log-redshift time
 *
 * Evolve the system until @alphaf.
 *
 */
void
nc_hipert_two_fluids_evolve (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alphaf)
{
  NcHIPertTwoFluidsPrivate * const self = ptf->priv;
  NcHIPert *pert                        = NC_HIPERT (ptf);
  NcHIPertPrivate * const pself         = nc_hipert_get_private (pert);
  NcHIPertTwoFluidsArg *arg             = self->arg;
  gdouble alpha_i                       = 0.0;
  gint flag;

  arg->cosmo = cosmo;
  arg->ptf   = ptf;

  if (NV_LENGTH_S (pself->y) != NC_HIPERT_ITWO_FLUIDS_VARS_LEN)
    g_error ("nc_hipert_two_fluids_evolve: cannot evolve subsidiary approximated system, use the appropriated evolve function.");

#ifdef HAVE_SUNDIALS_ARKODE

  if (TRUE)
  {
    flag = ARKodeSetUserData (self->arkode, arg);
    NCM_CVODE_CHECK (&flag, "ARKodeSetUserData", 1, );

    flag = ARKodeSetStopTime (self->arkode, alphaf);
    NCM_CVODE_CHECK (&flag, "ARKodeSetStopTime", 1, );

    flag = ARKodeEvolve (self->arkode, alphaf, pself->y, &alpha_i, ARK_NORMAL);
    NCM_CVODE_CHECK (&flag, "ARKodeEvolve[nc_hipert_two_fluids_evolve]", 1, );

    pself->alpha0 = alpha_i;
  }
  else
#endif /* HAVE_SUNDIALS_ARKODE */
  {
    flag = CVodeSetUserData (pself->cvode, arg);
    NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

    flag = CVodeSetStopTime (pself->cvode, alphaf);
    NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

    flag = CVode (pself->cvode, alphaf, pself->y, &alpha_i, CV_NORMAL);
    NCM_CVODE_CHECK (&flag, "CVode[nc_hipert_two_fluids_evolve]", 1, );

    pself->alpha0 = alpha_i;
  }
}

/**
 * nc_hipert_two_fluids_evolve_array:
 * @ptf: a #NcHIPertTwoFluids
 * @cosmo: a #NcHICosmo
 * @alphaf: the final log-redshift time
 * @step_reltol: the step size relative tolerance
 * @step_abstol: the step size absolute tolerance
 *
 * Evolve the system until @alphaf and store the results in an array. Only steps with
 * difference satisfying the relative and absolute tolerances are stored.
 *
 * Returns: (transfer full): a #NcmMatrix
 */
NcmMatrix *
nc_hipert_two_fluids_evolve_array (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble alphaf, gdouble step_reltol, gdouble step_abstol)
{
  NcHIPertTwoFluidsPrivate * const self = ptf->priv;
  NcHIPert *pert                        = NC_HIPERT (ptf);
  NcHIPertPrivate * const pself         = nc_hipert_get_private (pert);
  NcHIPertTwoFluidsArg *arg             = self->arg;
  gdouble alpha_i                       = 0.0;
  GArray *array                         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  gdouble last_step                     = pself->alpha0;
  gint flag;

  g_assert_cmpfloat (step_reltol, >, GSL_DBL_EPSILON);
  g_assert_cmpfloat (step_abstol, >=, 0.0);

  arg->cosmo = cosmo;
  arg->ptf   = ptf;

  if (NV_LENGTH_S (pself->y) != NC_HIPERT_ITWO_FLUIDS_VARS_LEN)
    g_error ("nc_hipert_two_fluids_evolve: cannot evolve subsidiary approximated system, use the appropriated evolve function.");

#ifdef HAVE_SUNDIALS_ARKODE

  if (TRUE)
  {
    flag = ARKodeSetUserData (self->arkode, arg);
    NCM_CVODE_CHECK (&flag, "ARKodeSetUserData", 1, NULL);

    flag = ARKodeSetStopTime (self->arkode, alphaf);
    NCM_CVODE_CHECK (&flag, "ARKodeSetStopTime", 1, NULL);

    g_array_append_val (array, pself->alpha0);
    g_array_append_vals (array, NV_DATA_S (pself->y), NC_HIPERT_ITWO_FLUIDS_VARS_LEN);

    while (alphaf > pself->alpha0)
    {
      flag = ARKodeEvolve (self->arkode, alphaf, pself->y, &alpha_i, ARK_ONE_STEP);
      NCM_CVODE_CHECK (&flag, "ARKodeEvolve[nc_hipert_two_fluids_evolve]", 1, NULL);

      pself->alpha0 = alpha_i;

      if (fabs (alpha_i - last_step) > step_reltol * fabs (last_step) + step_abstol)
      {
        g_array_append_val (array, alpha_i);
        g_array_append_vals (array, NV_DATA_S (pself->y), NC_HIPERT_ITWO_FLUIDS_VARS_LEN);
        last_step = alpha_i;
      }
    }
  }
  else
#endif /* HAVE_SUNDIALS_ARKODE */
  {
    flag = CVodeSetUserData (pself->cvode, arg);
    NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, NULL);

    g_array_append_val (array, pself->alpha0);
    g_array_append_vals (array, NV_DATA_S (pself->y), NC_HIPERT_ITWO_FLUIDS_VARS_LEN);

    CVodeSetStopTime (pself->cvode, alphaf);

    while (alphaf > pself->alpha0)
    {
      flag = CVode (pself->cvode, alphaf, pself->y, &alpha_i, CV_ONE_STEP);
      NCM_CVODE_CHECK (&flag, "CVode[nc_hipert_two_fluids_evolve_array]", 1, NULL);

      pself->alpha0 = alpha_i;

      if (fabs (alpha_i - last_step) > step_reltol * fabs (last_step) + step_abstol)
      {
        g_array_append_val (array, alpha_i);
        g_array_append_vals (array, NV_DATA_S (pself->y), NC_HIPERT_ITWO_FLUIDS_VARS_LEN);
        last_step = alpha_i;
      }
    }
  }

  {
    const guint ncols = NC_HIPERT_ITWO_FLUIDS_VARS_LEN + 1;
    NcmMatrix *res    = ncm_matrix_new_array (array, ncols);

    g_array_unref (array);

    return res;
  }
}

/**
 * nc_hipert_two_fluids_peek_state:
 * @ptf: a #NcHIPertTwoFluids
 * @cosmo: a #NcHICosmo
 * @alpha: (out): current time
 *
 * Get the current time and values of the numerical solution.
 *
 * Returns: (transfer none): current solution state.
 */
NcmVector *
nc_hipert_two_fluids_peek_state (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, gdouble *alpha)
{
  NcHIPertTwoFluidsPrivate * const self = ptf->priv;
  NcHIPert *pert                        = NC_HIPERT (ptf);
  NcHIPertPrivate * const pself         = nc_hipert_get_private (pert);
  const guint len                       = NV_LENGTH_S (pself->y);
  guint i;

  if (len == NC_HIPERT_ITWO_FLUIDS_VARS_LEN)
  {
    for (i = 0; i < len; i++)
    {
      ncm_vector_set (self->state, i, NV_Ith_S (pself->y, i));
    }
  }
  else if (len == 4)
  {
    const gdouble k            = nc_hipert_get_mode_k (pert);
    NcHIPertITwoFluidsEOM *eom = nc_hipert_itwo_fluids_eom_eval (NC_HIPERT_ITWO_FLUIDS (cosmo), pself->alpha0, k);

    const gdouble Q_R2 = NV_Ith_S (pself->y, 0);
    const gdouble P_R2 = NV_Ith_S (pself->y, 1);

    const gdouble Q_I2 = NV_Ith_S (pself->y, 2);
    const gdouble P_I2 = NV_Ith_S (pself->y, 3);

    const complex double A_R2 = NC_HIPERT_TWO_FLUIDS_QP2A (Q_R2, P_R2);
    const complex double A_I2 = NC_HIPERT_TWO_FLUIDS_QP2A (Q_I2, P_I2);

    const gdouble gammabar11 = eom->gammabar11;
    const gdouble gammabar22 = eom->gammabar22;
    const gdouble gammabar12 = eom->gammabar12;
    const gdouble taubar12   = eom->taubar;
    const gdouble dsigma1    = eom->nu1 - 0.5 * gammabar11;
    const gdouble dsigma2    = eom->nu2 - 0.5 * gammabar22;

    const complex double Lambda = gammabar12 + I * taubar12;

    const complex double A_R1 = (
      +A_R2        * (Lambda     / (2.0 * (dsigma1 - dsigma2)) + gammabar11    * gammabar12 / (8.0 * dsigma1 * (dsigma1 + dsigma2)))
      - conj (A_R2) * (gammabar12 / (2.0 * (dsigma1 + dsigma2)) + conj (Lambda) * gammabar11 / (8.0 * dsigma1 * (dsigma1 - dsigma2)))
                                ) / (1.0 - gsl_pow_2 (gammabar11 / (4.0 * dsigma1)));
    const complex double A_I1 = (
      +A_I2        * (Lambda     / (2.0 * (dsigma1 - dsigma2)) + gammabar11    * gammabar12 / (8.0 * dsigma1 * (dsigma1 + dsigma2)))
      - conj (A_I2) * (gammabar12 / (2.0 * (dsigma1 + dsigma2)) + conj (Lambda) * gammabar11 / (8.0 * dsigma1 * (dsigma1 - dsigma2)))
                                ) / (1.0 - gsl_pow_2 (gammabar11 / (4.0 * dsigma1)));

    const gdouble Q_R1 = NC_HIPERT_TWO_FLUIDS_A2Q (A_R1);
    const gdouble P_R1 = NC_HIPERT_TWO_FLUIDS_A2P (A_R1);

    const gdouble Q_I1 = NC_HIPERT_TWO_FLUIDS_A2Q (A_I1);
    const gdouble P_I1 = NC_HIPERT_TWO_FLUIDS_A2P (A_I1);

    ncm_vector_set (self->state, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1, Q_R1);
    ncm_vector_set (self->state, NC_HIPERT_ITWO_FLUIDS_VARS_P_R1, P_R1);

    ncm_vector_set (self->state, NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2, Q_R2);
    ncm_vector_set (self->state, NC_HIPERT_ITWO_FLUIDS_VARS_P_R2, P_R2);

    ncm_vector_set (self->state, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1, Q_I1);
    ncm_vector_set (self->state, NC_HIPERT_ITWO_FLUIDS_VARS_P_I1, P_I1);

    ncm_vector_set (self->state, NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2, Q_I2);
    ncm_vector_set (self->state, NC_HIPERT_ITWO_FLUIDS_VARS_P_I2, P_I2);
  }
  else
  {
    g_assert_not_reached ();
  }


  alpha[0] = pself->alpha0;

  return self->state;
}

static gdouble
_nc_hipert_two_fluids_wkb_limit_mode1_root (gdouble alpha, gpointer userdata)
{
  NcHIPertTwoFluidsArg *arg = (NcHIPertTwoFluidsArg *) userdata;

  const gdouble k            = nc_hipert_get_mode_k (NC_HIPERT (arg->ptf));
  NcHIPertITwoFluidsWKB *wkb = nc_hipert_itwo_fluids_wkb_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);

  return log (wkb->mode1_scale / arg->prec);
}

static gdouble
_nc_hipert_two_fluids_wkb_limit_mode2_root (gdouble alpha, gpointer userdata)
{
  NcHIPertTwoFluidsArg *arg  = (NcHIPertTwoFluidsArg *) userdata;
  const gdouble k            = nc_hipert_get_mode_k (NC_HIPERT (arg->ptf));
  NcHIPertITwoFluidsWKB *wkb = nc_hipert_itwo_fluids_wkb_eval (NC_HIPERT_ITWO_FLUIDS (arg->cosmo), alpha, k);

  return log (wkb->mode2_scale / arg->prec);
}

/**
 * nc_hipert_two_fluids_get_wkb_limit:
 * @ptf: a #NcHIPertTwoFluids
 * @cosmo: a #NcHICosmo
 * @main_mode: main mode
 * @alpha_i: initial try
 * @prec: precision
 *
 * Get the initial time where the WKB approximation is valid up to @prec. The algorithm
 * uses @alpha_i as the initial guess.
 *
 * Returns: initial time $\alpha_i$.
 */
gdouble
nc_hipert_two_fluids_get_wkb_limit (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, guint main_mode, gdouble alpha_i, gdouble prec)
{
  NcHIPert *pert                        = NC_HIPERT (ptf);
  NcHIPertPrivate * const pself         = nc_hipert_get_private (pert);
  NcHIPertTwoFluidsPrivate * const self = ptf->priv;
  NcHIPertTwoFluidsArg *arg             = self->arg;
  gint status;
  gint iter = 0, max_iter = 1000000;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;
  gdouble alpha0, alpha1, alpha = alpha_i;

  switch (main_mode)
  {
    case 1:
      F.function = &_nc_hipert_two_fluids_wkb_limit_mode1_root;
      break;
    case 2:
      F.function = &_nc_hipert_two_fluids_wkb_limit_mode2_root;
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  F.params = arg;

  arg->cosmo = cosmo;
  arg->ptf   = ptf;
  arg->prec  = prec;

  alpha0 = alpha_i;

  while (GSL_FN_EVAL (&F, alpha0) > 0.0)
  {
    alpha0 -= 10.0;
  }

  alpha1 = alpha_i;

  while (GSL_FN_EVAL (&F, alpha1) < 0.0)
  {
    alpha1 += 10.0;
  }

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, alpha0, alpha1);

  do {
    iter++;
    status = gsl_root_fsolver_iterate (s);

    if (status)
    {
      g_warning ("%s", gsl_strerror (status));
      break;
    }

    alpha  = gsl_root_fsolver_root (s);
    alpha0 = gsl_root_fsolver_x_lower (s);
    alpha1 = gsl_root_fsolver_x_upper (s);

    status = gsl_root_test_residual (GSL_FN_EVAL (&F, alpha), pself->reltol);

    if (status == GSL_CONTINUE)
      status = gsl_root_test_interval (alpha0, alpha1, 0, prec);
  } while (status == GSL_CONTINUE && iter < max_iter);

  if (status)
    g_warning ("%s", gsl_strerror (status));

  if (iter >= max_iter)
    g_warning ("nc_hipert_two_fluids_get_wkb_limit: maximum number of iterations reached.");

  gsl_root_fsolver_free (s);

  return alpha;
}

static gdouble
_nc_hipert_two_fluids_zeta_spectrum (const gdouble lnk, gpointer userdata)
{
  NcHIPertTwoFluidsArg *arg = (NcHIPertTwoFluidsArg *) userdata;
  NcHIPert *pert            = NC_HIPERT (arg->ptf);
  NcmVector *init_cond      = ncm_vector_new (NC_HIPERT_ITWO_FLUIDS_VARS_LEN);
  const gdouble k           = exp (lnk);
  gdouble alpha_i0, alpha_i;

  nc_hipert_set_mode_k (pert, k);

  alpha_i = nc_hipert_two_fluids_get_wkb_limit (arg->ptf, arg->cosmo, arg->mode - 1, arg->alpha_i, arg->prec);
  nc_hipert_two_fluids_get_init_cond_zetaS (arg->ptf, arg->cosmo, alpha_i, arg->mode, M_PI * 0.25, init_cond);
  nc_hipert_two_fluids_set_init_cond (arg->ptf, arg->cosmo, alpha_i, arg->mode, FALSE, init_cond);

  nc_hipert_two_fluids_evolve (arg->ptf, arg->cosmo, arg->alpha);

  ncm_vector_clear (&init_cond);

  {
    NcmVector *res       = nc_hipert_two_fluids_peek_state (arg->ptf, arg->cosmo, &alpha_i0);
    const gdouble zeta_r = ncm_vector_get (res, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R);
    const gdouble zeta_i = ncm_vector_get (res, NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I);

    return log (gsl_pow_3 (k) * (gsl_pow_2 (zeta_r) + gsl_pow_2 (zeta_i)));
  }
}

/**
 * nc_hipert_two_fluids_compute_zeta_spectrum:
 * @ptf: a #NcHIPertTwoFluids
 * @cosmo: a #NcHICosmo
 * @mode: the mode
 * @alpha_i: the initial log-redshift time
 * @alpha: the log-redshift time
 * @ki: the initial k
 * @kf: the final k
 * @nnodes: number of knots
 *
 * Compute the spectra for the two fluids system.
 *
 * Returns: (transfer full): a #NcmSpline
 */
NcmSpline *
nc_hipert_two_fluids_compute_zeta_spectrum (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, guint mode, gdouble alpha_i, gdouble alpha, gdouble ki, gdouble kf, guint nnodes)
{
  NcHIPertTwoFluidsPrivate * const self = ptf->priv;
  NcHIPertTwoFluidsArg *arg             = self->arg;
  NcmSpline *spline                     = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  const gdouble lnki                    = log (ki);
  const gdouble lnkf                    = log (kf);
  NcmVector *xv                         = ncm_vector_new (nnodes);
  NcmVector *yv                         = ncm_vector_new (nnodes);
  gsl_function F;
  guint i;

  F.params = arg;

  arg->cosmo   = cosmo;
  arg->ptf     = ptf;
  arg->prec    = 1.0e-5;
  arg->alpha_i = alpha_i;
  arg->alpha   = alpha;
  arg->mode    = mode;

  F.function = &_nc_hipert_two_fluids_zeta_spectrum;

  for (i = 0; i < nnodes; i++)
  {
    const gdouble lnk = lnki + (lnkf - lnki) * i / (nnodes - 1);
    const gdouble y   = F.function (lnk, arg);

    ncm_vector_set (xv, i, lnk);
    ncm_vector_set (yv, i, y);
  }

  ncm_spline_set (spline, xv, yv, TRUE);

  ncm_vector_free (xv);
  ncm_vector_free (yv);

  return spline;
}

/**
 * nc_hipert_two_fluids_state_evol_mode:
 * @ptf: a #NcHIPertTwoFluids
 * @cosmo: a #NcHICosmo
 *
 * Computes the evolution of the two-fluid perturbation state across a predefined
 * interval in $\alpha$, as specified by the calls to
 * nc_hipert_two_fluids_set_initial_time() and nc_hipert_two_fluids_set_final_time().
 *
 * The evolution is performed per quantized mode. If the initial time $\alpha_i$ lies
 * before the WKB approximation limit, the initial conditions are obtained using the WKB
 * approximation. If the final time $\alpha_f$ extends beyond this limit, the evolution
 * is continued numerically.
 *
 * This function returns an interpolator for the complex-valued state variables $(\zeta,
 * Q, P_\zeta, P_Q)$ for each mode.
 *
 * Returns: (transfer full): a newly allocated #NcHIPertTwoFluidsStateInterp containing
 * spline interpolators for the perturbation evolution.
 */
NcHIPertTwoFluidsStateInterp *
nc_hipert_two_fluids_evol_mode (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo)
{
  NcHIPert *pert                        = NC_HIPERT (ptf);
  NcHIPertTwoFluidsPrivate * const self = ptf->priv;
  const gdouble alpha_try               = 0.5 * (self->alpha_i + self->alpha_f);
  const gdouble alpha_i1                = nc_hipert_two_fluids_get_wkb_limit (ptf, cosmo, 1, alpha_try, self->wkb_reltol);
  const gdouble alpha_i2                = nc_hipert_two_fluids_get_wkb_limit (ptf, cosmo, 2, alpha_try, self->wkb_reltol);
  const guint prealloc_n                = 1000;
  GArray *alpha1                        = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), prealloc_n);
  GArray *alpha2                        = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), prealloc_n);
  NcmVector *init_cond_vec              = ncm_vector_new (NC_HIPERT_ITWO_FLUIDS_VARS_LEN);
  GArray *state1[NC_HIPERT_ITWO_FLUIDS_VARS_LEN];
  GArray *state2[NC_HIPERT_ITWO_FLUIDS_VARS_LEN];
  guint i;

  for (i = 0; i < NC_HIPERT_ITWO_FLUIDS_VARS_LEN; i++)
  {
    state1[i] = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    state2[i] = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  }

#define G_ARRAY_APPEND_FUNCVAL(array, expr)    \
        G_STMT_START {                         \
          const gdouble __tmp = (expr);        \
          g_array_append_val ((array), __tmp); \
        } G_STMT_END

  if (self->alpha_i < alpha_i1)
  {
    NcHIPertITwoFluidsWKB *wkb_state;
    gdouble alpha_f    = GSL_MIN (self->alpha_f, alpha_i1);
    guint n_wkb_states = 20 * (fabs (alpha_f - self->alpha_i) + 10);

    /*
     *  Mode 1 starts before WKB limit. Computing states using WKB approximation.
     */
    for (i = 0; i < n_wkb_states; i++)
    {
      const gdouble alpha = self->alpha_i + i * (alpha_f - self->alpha_i) / (n_wkb_states * 1.0);

      nc_hipert_two_fluids_wkb (ptf, cosmo, alpha, &wkb_state);

      g_array_append_val (alpha1, alpha);
      G_ARRAY_APPEND_FUNCVAL (state1[NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R], creal (wkb_state->state.zeta1));
      G_ARRAY_APPEND_FUNCVAL (state1[NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I], cimag (wkb_state->state.zeta1));
      G_ARRAY_APPEND_FUNCVAL (state1[NC_HIPERT_ITWO_FLUIDS_VARS_S_R], creal (wkb_state->state.Q1));
      G_ARRAY_APPEND_FUNCVAL (state1[NC_HIPERT_ITWO_FLUIDS_VARS_S_I], cimag (wkb_state->state.Q1));
      G_ARRAY_APPEND_FUNCVAL (state1[NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I], creal (wkb_state->state.Pzeta1));
      G_ARRAY_APPEND_FUNCVAL (state1[NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R], cimag (wkb_state->state.Pzeta1));
      G_ARRAY_APPEND_FUNCVAL (state1[NC_HIPERT_ITWO_FLUIDS_VARS_PS_I], creal (wkb_state->state.PQ1));
      G_ARRAY_APPEND_FUNCVAL (state1[NC_HIPERT_ITWO_FLUIDS_VARS_PS_R], cimag (wkb_state->state.PQ1));
    }
  }

  if (self->alpha_f > alpha_i1)
  {
    NcmMatrix *evol_mat;

    /*
     *  Mode 1 extends beyond WKB limit. Computing states using numerical integration.
     */
    nc_hipert_two_fluids_get_init_cond_zetaS (ptf, cosmo, alpha_i1, 1, 0.25 * M_PI, init_cond_vec);
    nc_hipert_two_fluids_set_init_cond (ptf, cosmo, alpha_i1, 0, FALSE, init_cond_vec);
    evol_mat = nc_hipert_two_fluids_evolve_array (ptf, cosmo, self->alpha_f, 1.0e-3, 0.0);

    for (i = 0; i < ncm_matrix_nrows (evol_mat); i++)
    {
      gint j;

      G_ARRAY_APPEND_FUNCVAL (alpha1, ncm_matrix_get (evol_mat, i, 0));

      for (j = 0; j < NC_HIPERT_ITWO_FLUIDS_VARS_LEN; j++)
      {
        G_ARRAY_APPEND_FUNCVAL (state1[j], ncm_matrix_get (evol_mat, i, j + 1));
      }
    }

    ncm_matrix_free (evol_mat);
  }

  if (self->alpha_i < alpha_i2)
  {
    NcHIPertITwoFluidsWKB *wkb_state;
    gdouble alpha_f    = GSL_MIN (self->alpha_f, alpha_i2);
    guint n_wkb_states = 20 * (fabs (alpha_f - self->alpha_i) + 10);

    /*
     *  Mode 2 starts before WKB limit. Computing states using WKB approximation.
     */
    for (i = 0; i < n_wkb_states; i++)
    {
      const gdouble alpha = self->alpha_i + i * (alpha_f - self->alpha_i) / (n_wkb_states * 1.0);

      nc_hipert_two_fluids_wkb (ptf, cosmo, alpha, &wkb_state);

      g_array_append_val (alpha2, alpha);
      G_ARRAY_APPEND_FUNCVAL (state2[NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R], creal (wkb_state->state.zeta2));
      G_ARRAY_APPEND_FUNCVAL (state2[NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I], cimag (wkb_state->state.zeta2));
      G_ARRAY_APPEND_FUNCVAL (state2[NC_HIPERT_ITWO_FLUIDS_VARS_S_R], creal (wkb_state->state.Q2));
      G_ARRAY_APPEND_FUNCVAL (state2[NC_HIPERT_ITWO_FLUIDS_VARS_S_I], cimag (wkb_state->state.Q2));
      G_ARRAY_APPEND_FUNCVAL (state2[NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I], creal (wkb_state->state.Pzeta2));
      G_ARRAY_APPEND_FUNCVAL (state2[NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R], cimag (wkb_state->state.Pzeta2));
      G_ARRAY_APPEND_FUNCVAL (state2[NC_HIPERT_ITWO_FLUIDS_VARS_PS_I], creal (wkb_state->state.PQ2));
      G_ARRAY_APPEND_FUNCVAL (state2[NC_HIPERT_ITWO_FLUIDS_VARS_PS_R], cimag (wkb_state->state.PQ2));
    }
  }

  if (self->alpha_f > alpha_i2)
  {
    NcmMatrix *evol_mat;

    /*
     *  Mode 2 extends beyond WKB limit. Computing states using numerical integration.
     */
    nc_hipert_two_fluids_get_init_cond_zetaS (ptf, cosmo, alpha_i2, 2, 0.25 * M_PI, init_cond_vec);
    nc_hipert_two_fluids_set_init_cond (ptf, cosmo, alpha_i2, 0, FALSE, init_cond_vec);
    evol_mat = nc_hipert_two_fluids_evolve_array (ptf, cosmo, self->alpha_f, 1.0e-3, 0.0);

    for (i = 0; i < ncm_matrix_nrows (evol_mat); i++)
    {
      gint j;

      G_ARRAY_APPEND_FUNCVAL (alpha2, ncm_matrix_get (evol_mat, i, 0));

      for (j = 0; j < NC_HIPERT_ITWO_FLUIDS_VARS_LEN; j++)
      {
        G_ARRAY_APPEND_FUNCVAL (state2[j], ncm_matrix_get (evol_mat, i, j + 1));
      }
    }

    ncm_matrix_free (evol_mat);
  }

  {
    NcHIPertTwoFluidsStateInterp *sinterp = g_new0 (NcHIPertTwoFluidsStateInterp, 1);

    sinterp->interp_mode = 1;
    sinterp->state.k     = nc_hipert_get_mode_k (pert);
    sinterp->state.norma = nc_hipert_itwo_fluids_eval_unit (NC_HIPERT_ITWO_FLUIDS (cosmo));

    for (i = 0; i < NC_HIPERT_ITWO_FLUIDS_VARS_LEN; i++)
    {
      sinterp->mode1_splines[i] = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
      sinterp->mode2_splines[i] = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());

      ncm_spline_set_array (sinterp->mode1_splines[i], alpha1, state1[i], TRUE);
      ncm_spline_set_array (sinterp->mode2_splines[i], alpha2, state2[i], TRUE);
    }

    g_array_unref (alpha1);
    g_array_unref (alpha2);

    for (i = 0; i < NC_HIPERT_ITWO_FLUIDS_VARS_LEN; i++)
    {
      g_array_unref (state1[i]);
      g_array_unref (state2[i]);
    }

    ncm_vector_free (init_cond_vec);

    return sinterp;
  }
}

/**
 * nc_hipert_two_fluids_compute_spectrum:
 * @ptf: a #NcHIPertTwoFluids
 * @cosmo: a #NcHICosmo
 * @alpha: the scale factor
 * @k_a: (element-type gdouble): an array of k
 * @logger: (nullable) (scope call): a #NcHIPertTwoFluidsLogger
 *
 * Computes the spectrum at a given scale factor.
 *
 * Returns: a #NcHIPertTwoFluidsStateInterp
 */
NcHIPertTwoFluidsStateInterp *
nc_hipert_two_fluids_compute_spectrum (NcHIPertTwoFluids *ptf, NcHICosmo *cosmo, const gdouble alpha, GArray *k_a, NcHIPertTwoFluidsLogger logger)
{
  NcHIPert *pert                        = NC_HIPERT (ptf);
  NcHIPertTwoFluidsPrivate * const self = ptf->priv;
  const gdouble alpha_try               = 0.5 * (self->alpha_i + self->alpha_f);
  NcmVector *initial_condition          = ncm_vector_new (NC_HIPERT_ITWO_FLUIDS_VARS_LEN);
  GArray *k_array                       = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), k_a->len);
  GArray *state1[NC_HIPERT_ITWO_FLUIDS_VARS_LEN];
  GArray *state2[NC_HIPERT_ITWO_FLUIDS_VARS_LEN];
  NcmVector *state;
  gint mode;
  guint i;

  g_assert_cmpint (k_a->len, >, 0);
  g_assert_cmpint (g_array_get_element_size (k_a), ==, sizeof (gdouble));

  g_assert_cmpfloat (g_array_index (k_a, gdouble, 0), >, 0.0);

  for (i = 0; i < k_a->len - 1; i++)
  {
    g_assert_cmpfloat (g_array_index (k_a, gdouble, i + 1), >, g_array_index (k_a, gdouble, i));
  }

  for (i = 0; i < NC_HIPERT_ITWO_FLUIDS_VARS_LEN; i++)
  {
    state1[i] = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    state2[i] = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  }

  for (i = 0; i < k_a->len; i++)
  {
    const gdouble k = g_array_index (k_a, gdouble, i);

    nc_hipert_set_mode_k (pert, k);
    g_array_append_val (k_array, k);

    mode = 1;
    {
      const gdouble alpha_i = nc_hipert_two_fluids_get_wkb_limit (ptf, cosmo, mode, alpha_try, self->wkb_reltol);
      gdouble alpha_last;
      gint j;

      nc_hipert_two_fluids_get_init_cond_zetaS (ptf, cosmo, alpha_i, mode, M_PI * 0.25, initial_condition);
      nc_hipert_two_fluids_set_init_cond (ptf, cosmo, alpha_i, mode, FALSE, initial_condition);

      nc_hipert_two_fluids_evolve (ptf, cosmo, alpha);

      state = nc_hipert_two_fluids_peek_state (ptf, cosmo, &alpha_last);

      for (j = 0; j < NC_HIPERT_ITWO_FLUIDS_VARS_LEN; j++)
      {
        G_ARRAY_APPEND_FUNCVAL (state1[j], ncm_vector_get (state, j));
      }
    }

    mode = 2;
    {
      const gdouble alpha_i = nc_hipert_two_fluids_get_wkb_limit (ptf, cosmo, mode, alpha_try, self->wkb_reltol);
      gdouble alpha_last;
      gint j;

      nc_hipert_two_fluids_get_init_cond_zetaS (ptf, cosmo, alpha_i, mode, M_PI * 0.25, initial_condition);
      nc_hipert_two_fluids_set_init_cond (ptf, cosmo, alpha_i, mode, FALSE, initial_condition);

      nc_hipert_two_fluids_evolve (ptf, cosmo, alpha);

      state = nc_hipert_two_fluids_peek_state (ptf, cosmo, &alpha_last);

      for (j = 0; j < NC_HIPERT_ITWO_FLUIDS_VARS_LEN; j++)
      {
        G_ARRAY_APPEND_FUNCVAL (state2[j], ncm_vector_get (state, j));
      }
    }

    if (logger)
      logger (i, k_a->len);
  }

  ncm_vector_free (initial_condition);

  {
    NcHIPertTwoFluidsStateInterp *sinterp = g_new0 (NcHIPertTwoFluidsStateInterp, 1);

    sinterp->interp_mode = 2;
    sinterp->state.alpha = alpha;
    sinterp->state.norma = nc_hipert_itwo_fluids_eval_unit (NC_HIPERT_ITWO_FLUIDS (cosmo));

    for (i = 0; i < NC_HIPERT_ITWO_FLUIDS_VARS_LEN; i++)
    {
      sinterp->mode1_splines[i] = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
      sinterp->mode2_splines[i] = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());

      ncm_spline_set_array (sinterp->mode1_splines[i], k_array, state1[i], TRUE);
      ncm_spline_set_array (sinterp->mode2_splines[i], k_array, state2[i], TRUE);
    }

    g_array_unref (k_array);

    for (i = 0; i < NC_HIPERT_ITWO_FLUIDS_VARS_LEN; i++)
    {
      g_array_unref (state1[i]);
      g_array_unref (state2[i]);
    }

    return sinterp;
  }
}

