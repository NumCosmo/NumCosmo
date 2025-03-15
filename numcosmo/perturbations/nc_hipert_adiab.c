/***************************************************************************
 *            nc_hipert_adiab.c
 *
 *  Tue June 03 17:20:42 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hipert_adiab.c
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
 * NcHIPertAdiab:
 *
 * Perturbation object for adiabatic mode only.
 *
 * This object provides the computation of the adiabatic mode for the cosmological
 * perturbations. It solves the equation of motion for the gauge invariant variable
 * (see [Vitenti (2013)][XVitenti2013] for notation and details)
 * $$
 * \zeta \equiv \Psi - \frac{2\bar{K}}{\kappa(\bar{\rho} + \bar{p})} + H\mathcal{V}.
 * $$
 * Its conjugated
 * momentum is give by
 * \begin{split}
 * P_\zeta &= \frac{2\bar{D}^2_\bar{K}\Psi}{x^3H},
 * \end{split}
 *
 * The equations of motion in their first order form are
 * \begin{align}
 * \zeta^\prime &= \frac{P_\zeta}{m_\zeta}, \\\\
 * P_\zeta^\prime &= -m_\zeta\mu_\zeta^2\zeta. \end{align} The mass $m_\zeta$ and the
 * frequency $\mu_\zeta$ are defined by
 * \begin{align}
 * m_\zeta &= \frac{3\Delta_\bar{K}(\bar{\rho} + \bar{p})}{\rho_\text{crit0} N x^3
 * c_s^2 E^2},
 * \\\\
 * \mu_\zeta^2 &= x^2N^2c_s^2k^2,
 * \end{align}
 * where $\bar{\rho} + \bar{p}$ is the background total energy density plus pressure,
 * $E^2 = H^2/H_0^2$ is the dimensionless Hubble function squared (nc_hicosmo_E2()),
 * $c_s^2$ the speed of sound, $N$ is the lapse function that in this case (using
 * $\alpha$ as time variable) is $N \equiv \vert{}E\vert^{-1}$, $\rho_\text{crit0}$ is
 * the critical density today defined by $\rho_\text{crit0} \equiv 3H_0^2/\kappa$ and
 * $$
 * \Delta_\bar{K} \equiv \frac{k^2}{k^2 + \Omega_{k0}}.
 * $$
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "perturbations/nc_hipert_adiab.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_spline2d_bicubic.h"
#include "math/ncm_ode_spline.h"
#include "math/ncm_model_ctrl.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <gsl/gsl_roots.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcHIPertAdiab
{
  NcmCSQ1D parent_instance;
  gdouble k;
  NcmSpline2d *powspec_alpha;
  NcmSpline2d *powspec_gamma;
  NcmOdeSpline *ctime_forward;
  NcmOdeSpline *ctime_backward;
  NcmModelCtrl *model_ctrl;
};

G_DEFINE_INTERFACE (NcHIPertIAdiab, nc_hipert_iadiab, G_TYPE_OBJECT)
G_DEFINE_TYPE (NcHIPertAdiab, nc_hipert_adiab, NCM_TYPE_CSQ1D)

static void
nc_hipert_iadiab_default_init (NcHIPertIAdiabInterface *iface)
{
  iface->eval_xi         = NULL;
  iface->eval_F1         = NULL;
  iface->eval_nu         = NULL;
  iface->eval_m          = NULL;
  iface->eval_unit       = NULL;
  iface->eval_x          = NULL;
  iface->eval_p2Psi      = NULL;
  iface->eval_p2drho     = NULL;
  iface->eval_lapse      = NULL;
  iface->eval_tau_hubble = NULL;
  iface->eval_tau_jeans  = NULL;
  iface->eval_hubble     = NULL;
}

enum
{
  PROP_0,
  PROP_K,
  PROP_SIZE,
};

typedef struct _NcHIPertAdiabArg
{
  NcHICosmo *cosmo;
  NcHIPertAdiab *pa;
} NcHIPertAdiabArg;

static gdouble _nc_hipert_adiab_cosmic_time_integ_forward (const gdouble t, gdouble tau, void *params);
static gdouble _nc_hipert_adiab_cosmic_time_integ_backward (const gdouble t, gdouble tau, void *params);

static void
nc_hipert_adiab_init (NcHIPertAdiab *pa)
{
  NcmSpline *s1 = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  NcmSpline *s2 = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());

  pa->k              = 0.0;
  pa->powspec_alpha  = NULL;
  pa->powspec_gamma  = NULL;
  pa->ctime_forward  = ncm_ode_spline_new_full (s1, _nc_hipert_adiab_cosmic_time_integ_forward, 0.0, 0.0, 500.0);
  pa->ctime_backward = ncm_ode_spline_new_full (s2, _nc_hipert_adiab_cosmic_time_integ_backward, 0.0, 0.0, 500.0);
  pa->model_ctrl     = ncm_model_ctrl_new (NULL);

  ncm_ode_spline_auto_abstol (pa->ctime_forward, TRUE);
  ncm_ode_spline_auto_abstol (pa->ctime_backward, TRUE);

  ncm_spline_free (s1);
  ncm_spline_free (s2);
}

static void
_nc_hipert_adiab_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHIPertAdiab *pa = NC_HIPERT_ADIAB (object);

  g_return_if_fail (NC_IS_HIPERT_ADIAB (object));

  switch (prop_id)
  {
    case PROP_K:
      nc_hipert_adiab_set_k (pa, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_hipert_adiab_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHIPertAdiab *pa = NC_HIPERT_ADIAB (object);

  g_return_if_fail (NC_IS_HIPERT_ADIAB (object));

  switch (prop_id)
  {
    case PROP_K:
      g_value_set_double (value, nc_hipert_adiab_get_k (pa));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_hipert_adiab_dispose (GObject *object)
{
  NcHIPertAdiab *pa = NC_HIPERT_ADIAB (object);

  ncm_spline2d_clear (&pa->powspec_alpha);
  ncm_spline2d_clear (&pa->powspec_gamma);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_adiab_parent_class)->dispose (object);
}

static void
_nc_hipert_adiab_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_adiab_parent_class)->finalize (object);
}

static gdouble _nc_hipert_adiab_eval_xi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
static gdouble _nc_hipert_adiab_eval_F1 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
static gdouble _nc_hipert_adiab_eval_nu (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
static gdouble _nc_hipert_adiab_eval_m (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);

static void _nc_hipert_adiab_prepare (NcmCSQ1D *csq1d, NcmModel *model);

static void
nc_hipert_adiab_class_init (NcHIPertAdiabClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmCSQ1DClass *csq1d_class = NCM_CSQ1D_CLASS (klass);

  object_class->set_property = &_nc_hipert_adiab_set_property;
  object_class->get_property = &_nc_hipert_adiab_get_property;
  object_class->dispose      = &_nc_hipert_adiab_dispose;
  object_class->finalize     = &_nc_hipert_adiab_finalize;

  g_object_class_install_property (object_class,
                                   PROP_K,
                                   g_param_spec_double ("k",
                                                        NULL,
                                                        "Wave number",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  csq1d_class->eval_xi = &_nc_hipert_adiab_eval_xi;
  csq1d_class->eval_F1 = &_nc_hipert_adiab_eval_F1;
  csq1d_class->eval_nu = &_nc_hipert_adiab_eval_nu;
  csq1d_class->eval_m  = &_nc_hipert_adiab_eval_m;
  csq1d_class->prepare = &_nc_hipert_adiab_prepare;
}

static gdouble
_nc_hipert_adiab_eval_xi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  NcHIPertAdiab *pa = NC_HIPERT_ADIAB (csq1d);
  const gdouble k   = pa->k;

  return nc_hipert_iadiab_eval_xi (NC_HIPERT_IADIAB (model), t, k);
}

static gdouble
_nc_hipert_adiab_eval_F1 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  NcHIPertAdiab *pa = NC_HIPERT_ADIAB (csq1d);
  const gdouble k   = pa->k;

  return nc_hipert_iadiab_eval_F1 (NC_HIPERT_IADIAB (model), t, k);
}

static gdouble
_nc_hipert_adiab_eval_nu (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  NcHIPertAdiab *pa = NC_HIPERT_ADIAB (csq1d);
  const gdouble k   = pa->k;

  return nc_hipert_iadiab_eval_nu (NC_HIPERT_IADIAB (model), t, k);
}

static gdouble
_nc_hipert_adiab_eval_m (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  NcHIPertAdiab *pa = NC_HIPERT_ADIAB (csq1d);
  const gdouble k   = pa->k;

  return nc_hipert_iadiab_eval_m (NC_HIPERT_IADIAB (model), t, k);
}

static void
_nc_hipert_adiab_prepare (NcmCSQ1D *csq1d, NcmModel *model)
{
  g_assert (NC_IS_HIPERT_IADIAB (model));
}

/**
 * nc_hipert_iadiab_eval_xi:
 * @iad: a #NcHIPertIAdiab
 * @tau: $\tau$
 * @k: $k$
 *
 * Computes the value of $\xi = \ln(m\nu)$.
 *
 * Returns: $\xi$.
 */
gdouble
nc_hipert_iadiab_eval_xi (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  return NC_HIPERT_IADIAB_GET_IFACE (iad)->eval_xi (iad, tau, k);
}

/**
 * nc_hipert_iadiab_eval_F1:
 * @iad: a #NcHIPertIAdiab
 * @tau: $\tau$
 * @k: $k$
 *
 * Computes the value of $F_1 = \dot{\xi}/(2\nu)$.
 *
 * Returns: $F_1$.
 */
gdouble
nc_hipert_iadiab_eval_F1 (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  return NC_HIPERT_IADIAB_GET_IFACE (iad)->eval_F1 (iad, tau, k);
}

/**
 * nc_hipert_iadiab_eval_nu:
 * @iad: a #NcHIPertIAdiab
 * @tau: $\tau$
 * @k: $k$
 *
 * Computes the value of $\nu$.
 *
 * Returns: $\nu$.
 */
gdouble
nc_hipert_iadiab_eval_nu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  return NC_HIPERT_IADIAB_GET_IFACE (iad)->eval_nu (iad, tau, k);
}

/**
 * nc_hipert_iadiab_eval_m:
 * @iad: a #NcHIPertIAdiab
 * @tau: $\tau$
 * @k: $k$
 *
 * Computes the value of $m$.
 *
 * Returns: $m$.
 */
gdouble
nc_hipert_iadiab_eval_m (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  return NC_HIPERT_IADIAB_GET_IFACE (iad)->eval_m (iad, tau, k);
}

/**
 * nc_hipert_iadiab_eval_unit:
 * @iad: a #NcHIPertIAdiab
 *
 * Numerical factor for the power spectrum of the adiabatic mode.
 *
 * Returns: the numerical factor for the power spectrum of the adiabatic mode.
 */
gdouble
nc_hipert_iadiab_eval_unit (NcHIPertIAdiab *iad)
{
  return NC_HIPERT_IADIAB_GET_IFACE (iad)->eval_unit (iad);
}

/**
 * nc_hipert_iadiab_eval_x:
 * @iad: a #NcHIPertIAdiab
 * @tau: $\tau$
 *
 * Evaluates the value of $x = a_0/a$ at a given time $\tau$.
 *
 * Returns: $x$.
 */
gdouble
nc_hipert_iadiab_eval_x (NcHIPertIAdiab *iad, const gdouble tau)
{
  return NC_HIPERT_IADIAB_GET_IFACE (iad)->eval_x (iad, tau);
}

/**
 * nc_hipert_iadiab_eval_p2Psi:
 * @iad: a #NcHIPertIAdiab
 * @tau: $\tau$
 * @k: $k$
 *
 * Evaluates the conversion factor to convert the momentum of the adiabatic mode to the
 * gauge invariant variable $\Psi$.
 *
 * Returns: the conversion factor.
 */
gdouble
nc_hipert_iadiab_eval_p2Psi (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  return NC_HIPERT_IADIAB_GET_IFACE (iad)->eval_p2Psi (iad, tau, k);
}

/**
 * nc_hipert_iadiab_eval_p2drho:
 * @iad: a #NcHIPertIAdiab
 * @tau: $\tau$
 * @k: $k$
 *
 * Evaluates the conversion factor to convert the momentum of the adiabatic mode to the
 * gauge invariant variable $\delta\rho$.
 *
 * Returns: the conversion factor.
 */
gdouble
nc_hipert_iadiab_eval_p2drho (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  return NC_HIPERT_IADIAB_GET_IFACE (iad)->eval_p2drho (iad, tau, k);
}

/**
 * nc_hipert_iadiab_eval_lapse:
 * @iad: a #NcHIPertIAdiab
 * @tau: $\tau$
 *
 * Evaluates the lapse function at a given time $\tau$.
 *
 * Returns: the lapse function.
 */
gdouble
nc_hipert_iadiab_eval_lapse (NcHIPertIAdiab *iad, const gdouble tau)
{
  return NC_HIPERT_IADIAB_GET_IFACE (iad)->eval_lapse (iad, tau);
}

/**
 * nc_hipert_iadiab_eval_tau_hubble:
 * @iad: a #NcHIPertIAdiab
 * @k: $k$
 *
 * Evaluates the time at where the Hubble radius is equal to the wave number $k$.
 *
 * Returns: the time at where the Hubble radius is equal to the wave number $k$.
 */
gdouble
nc_hipert_iadiab_eval_tau_hubble (NcHIPertIAdiab *iad, const gdouble k)
{
  return NC_HIPERT_IADIAB_GET_IFACE (iad)->eval_tau_hubble (iad, k);
}

/**
 * nc_hipert_iadiab_eval_tau_jeans:
 * @iad: a #NcHIPertIAdiab
 * @k: $k$
 *
 * Evaluates the time at where the Jeans scale is equal to the wave number $k$.
 *
 * Returns: the time at where the Jeans scale is equal to the wave number $k$.
 */
gdouble
nc_hipert_iadiab_eval_tau_jeans (NcHIPertIAdiab *iad, const gdouble k)
{
  return NC_HIPERT_IADIAB_GET_IFACE (iad)->eval_tau_jeans (iad, k);
}

/**
 * nc_hipert_iadiab_eval_hubble:
 * @iad: a #NcHIPertIAdiab
 * @tau: $\tau$
 *
 * Evaluates the Hubble function at a given time $\tau$.
 *
 * Returns: the Hubble function.
 */
gdouble
nc_hipert_iadiab_eval_hubble (NcHIPertIAdiab *iad, const gdouble tau)
{
  return NC_HIPERT_IADIAB_GET_IFACE (iad)->eval_hubble (iad, tau);
}

/**
 * nc_hipert_adiab_new:
 *
 * Creates a new #NcHIPertAdiab object.
 *
 * Returns: (transfer full): a new #NcHIPertAdiab.
 */
NcHIPertAdiab *
nc_hipert_adiab_new (void)
{
  NcHIPertAdiab *pa = g_object_new (NC_TYPE_HIPERT_ADIAB,
                                    NULL);

  return pa;
}

/**
 * nc_hipert_adiab_ref:
 * @pa: a #NcHIPertAdiab.
 *
 * Increases the reference count of @pa.
 *
 * Returns: (transfer full): @pa.
 */
NcHIPertAdiab *
nc_hipert_adiab_ref (NcHIPertAdiab *pa)
{
  return g_object_ref (pa);
}

/**
 * nc_hipert_adiab_free:
 * @pa: a #NcHIPertAdiab.
 *
 * Decreases the reference count of @pa.
 *
 */
void
nc_hipert_adiab_free (NcHIPertAdiab *pa)
{
  g_object_unref (pa);
}

/**
 * nc_hipert_adiab_clear:
 * @pa: a #NcHIPertAdiab.
 *
 * Decreases the reference count of *@pa and sets *@pa to NULL.
 *
 */
void
nc_hipert_adiab_clear (NcHIPertAdiab **pa)
{
  g_clear_object (pa);
}

/**
 * nc_hipert_adiab_set_k:
 * @adiab: a #NcHIPertAdiab
 * @k: $k$
 *
 * Sets the wave number $k$.
 */
void
nc_hipert_adiab_set_k (NcHIPertAdiab *adiab, const gdouble k)
{
  adiab->k = k;
}

/**
 * nc_hipert_adiab_get_k:
 * @adiab: a #NcHIPertAdiab
 *
 * Gets the wave number $k$.
 *
 * Returns: $k$
 */
gdouble
nc_hipert_adiab_get_k (NcHIPertAdiab *adiab)
{
  return adiab->k;
}

typedef struct _NcHIPertAdiabCosmicTimeIntegArg
{
  NcHIPertAdiab *adiab;
  NcmModel *model;
} NcHIPertAdiabCosmicTimeIntegArg;

static gdouble
_nc_hipert_adiab_cosmic_time_integ_forward (gdouble t, gdouble tau, void *params)
{
  NcHIPertAdiabCosmicTimeIntegArg *arg = (NcHIPertAdiabCosmicTimeIntegArg *) params;

  return nc_hipert_iadiab_eval_lapse (NC_HIPERT_IADIAB (arg->model), tau);
}

static gdouble
_nc_hipert_adiab_cosmic_time_integ_backward (gdouble t, gdouble tau, void *params)
{
  NcHIPertAdiabCosmicTimeIntegArg *arg = (NcHIPertAdiabCosmicTimeIntegArg *) params;

  return -nc_hipert_iadiab_eval_lapse (NC_HIPERT_IADIAB (arg->model), -tau);
}

/**
 * nc_hipert_adiab_eval_cosmic_time:
 * @adiab: a #NcHIPertAdiab
 * @model: a #NcmModel
 * @tau: $\tau$
 *
 * Evaluates the cosmic time at a given conformal time $\tau$.
 *
 * Returns: the cosmic time.
 */
gdouble
nc_hipert_adiab_eval_cosmic_time (NcHIPertAdiab *adiab, NcmModel *model, const gdouble tau)
{
  if (ncm_model_ctrl_update (adiab->model_ctrl, model))
  {
    NcHIPertAdiabCosmicTimeIntegArg arg;

    arg.adiab = adiab;
    arg.model = model;

    ncm_ode_spline_prepare (adiab->ctime_forward, &arg);
    ncm_ode_spline_prepare (adiab->ctime_backward, &arg);
  }

  if (tau > 0)
  {
    NcmSpline *s = ncm_ode_spline_peek_spline (adiab->ctime_forward);

    return ncm_spline_eval (s, tau);
  }
  else
  {
    NcmSpline *s = ncm_ode_spline_peek_spline (adiab->ctime_backward);

    return ncm_spline_eval (s, -tau);
  }
}

/**
 * nc_hipert_adiab_eval_delta_critial:
 * @adiab: a #NcHIPertAdiab
 * @model: a #NcmModel
 * @tau: $\tau$
 *
 * Evaluates the critical density contrast at a given conformal time $\tau$.
 *
 * Returns: the critical density contrast.
 */
gdouble
nc_hipert_adiab_eval_delta_critial (NcHIPertAdiab *adiab, NcmModel *model, const gdouble tau)
{
  const gdouble tau_hubble = nc_hipert_iadiab_eval_tau_hubble (NC_HIPERT_IADIAB (model), adiab->k);
  const gdouble t_hubble   = nc_hipert_adiab_eval_cosmic_time (adiab, model, -tau_hubble);
  const gdouble t          = nc_hipert_adiab_eval_cosmic_time (adiab, model, tau);

  /* const gdouble E          = nc_hipert_iadiab_eval_hubble (NC_HIPERT_IADIAB (model), tau); */
  /* const gdouble delta_t    = t_hubble - t; */
  /* const gdouble F          = fabs (2.0 * E * delta_t / M_PI); */
  /* const gdouble cbrt_2     = cbrt (2.0); */
  /* const gdouble factor1    = cbrt (2.0 + 27.0 * gsl_pow_2 (F) - 3.0 * F * sqrt (12.0 + 81.0 * gsl_pow_2 (F))); */
  /* const gdouble factor2    = (1.0 + cbrt_2 / factor1 + factor1 / cbrt_2) / (3.0 * F); */
  /* printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", tau, tau_hubble, t, t_hubble, fabs (E * delta_t), 5.0 * t_hubble / t); */
  /*return sqrt (factor2); */

  return 5.0 * t_hubble / t;
}

static gdouble
_nc_hipert_eval_powspec_zeta_from_state (NcHIPertAdiab *adiab, NcmModel *model, NcmCSQ1DState *state, const gdouble k)
{
  const gdouble unit = nc_hipert_iadiab_eval_unit (NC_HIPERT_IADIAB (model));
  const gdouble n0   = 1.0 / ncm_c_two_pi_2 ();

  gdouble J11, J12, J22;

  ncm_csq1d_state_get_J (state, &J11, &J12, &J22);

  return gsl_pow_2 (unit) * n0 * (J11 / 2.0);
}

static gdouble
_nc_hipert_eval_powspec_Psi_from_state (NcHIPertAdiab *adiab, NcmModel *model, NcmCSQ1DState *state, const gdouble k)
{
  const gdouble unit  = nc_hipert_iadiab_eval_unit (NC_HIPERT_IADIAB (model));
  const gdouble tau   = ncm_csq1d_state_get_time (state);
  const gdouble p2Psi = nc_hipert_iadiab_eval_p2Psi (NC_HIPERT_IADIAB (model), tau, k);
  const gdouble n0    = 1.0 / ncm_c_two_pi_2 ();
  gdouble J11, J12, J22;

  ncm_csq1d_state_get_J (state, &J11, &J12, &J22);

  return gsl_pow_2 (unit * p2Psi) * n0 * (J22 / 2.0);
}

static gdouble
_nc_hipert_eval_powspec_drho_from_state (NcHIPertAdiab *adiab, NcmModel *model, NcmCSQ1DState *state, const gdouble k)
{
  const gdouble unit   = nc_hipert_iadiab_eval_unit (NC_HIPERT_IADIAB (model));
  const gdouble tau    = ncm_csq1d_state_get_time (state);
  const gdouble p2drho = nc_hipert_iadiab_eval_p2drho (NC_HIPERT_IADIAB (model), tau, k);
  const gdouble n0     = 1.0 / ncm_c_two_pi_2 ();
  gdouble J11, J12, J22;

  ncm_csq1d_state_get_J (state, &J11, &J12, &J22);

  return gsl_pow_2 (unit * p2drho) * n0 * (J22 / 2.0);
}

/**
 * nc_hipert_adiab_eval_powspec_zeta_at:
 * @adiab: a #NcHIPertAdiab
 * @model: a #NcmModel
 * @tau: $\tau$
 *
 * Evaluates the power spectrum of the gauge invariant variable $\zeta$ at a given time
 * $\tau$. The power spectrum is given by
 * $$
 * P_\zeta = u^2\frac{k^3}{2\pi^2} \frac{J_{11}}{2}.
 * $$
 * where $u$ is the numerical factor for the power spectrum of the adiabatic mode, $k$ is
 * the wave number.
 *
 * Returns: the power spectrum of the gauge invariant variable $\zeta$.
 */
gdouble
nc_hipert_adiab_eval_powspec_zeta_at (NcHIPertAdiab *adiab, NcmModel *model, const gdouble tau)
{
  const gdouble k = nc_hipert_adiab_get_k (adiab);
  NcmCSQ1DState state;

  ncm_csq1d_eval_at (NCM_CSQ1D (adiab), model, tau, &state);

  return _nc_hipert_eval_powspec_zeta_from_state (adiab, model, &state, k);
}

/**
 * nc_hipert_adiab_eval_powspec_Psi_at:
 * @adiab: a #NcHIPertAdiab
 * @model: a #NcmModel
 * @tau: $\tau$
 *
 * Evaluates the power spectrum of the gauge invariant variable $\Psi$ at a given time
 * $\tau$. The power spectrum is given by
 * $$
 * P_\Psi = u^2\frac{2\pi^2}{k^3} \frac{J_{22}}{2}.
 * $$
 *
 * Returns: the power spectrum of the gauge invariant variable $\Psi$.
 */
gdouble
nc_hipert_adiab_eval_powspec_Psi_at (NcHIPertAdiab *adiab, NcmModel *model, const gdouble tau)
{
  const gdouble k = nc_hipert_adiab_get_k (adiab);
  NcmCSQ1DState state;

  ncm_csq1d_eval_at (NCM_CSQ1D (adiab), model, tau, &state);

  return _nc_hipert_eval_powspec_Psi_from_state (adiab, model, &state, k);
}

/**
 * nc_hipert_adiab_eval_powspec_drho_at:
 * @adiab: a #NcHIPertAdiab
 * @model: a #NcmModel
 * @tau: $\tau$
 *
 * Evaluates the power spectrum of the gauge invariant variable $\delta\rho$ at a given time
 * $\tau$. The power spectrum is given by
 * $$
 * P_{\delta\rho} = u^2\frac{2\pi^2}{k^3} \frac{J_{22}}{2}.
 * $$
 *
 * Returns: the power spectrum of the gauge invariant variable $\delta\rho$.
 */
gdouble
nc_hipert_adiab_eval_powspec_drho_at (NcHIPertAdiab *adiab, NcmModel *model, const gdouble tau)
{
  const gdouble k = nc_hipert_adiab_get_k (adiab);
  NcmCSQ1DState state;

  ncm_csq1d_eval_at (NCM_CSQ1D (adiab), model, tau, &state);

  return _nc_hipert_eval_powspec_drho_from_state (adiab, model, &state, k);
}

/**
 * nc_hipert_adiab_prepare_spectrum:
 * @adiab: a #NcHIPertAdiab
 * @model: a #NcmModel
 * @k_array: (element-type gdouble): array of wave numbers
 * @tau_array: (element-type gdouble): array of times to evaluate the power spectrum
 *
 * Prepares the computation of the power spectrum of the adiabatic mode.
 *
 */
void
nc_hipert_adiab_prepare_spectrum (NcHIPertAdiab *adiab, NcmModel *model, GArray *k_array, GArray *tau_array)
{
  NcmCSQ1D *csq1d          = NCM_CSQ1D (adiab);
  NcmMatrix *powspec_alpha = ncm_matrix_new (k_array->len, tau_array->len);
  NcmMatrix *powspec_gamma = ncm_matrix_new (k_array->len, tau_array->len);
  guint i;

  for (i = 0; i < k_array->len; i++)
  {
    NcmCSQ1DState state;
    guint j;

    nc_hipert_adiab_set_k (adiab, g_array_index (k_array, gdouble, i));
    ncm_csq1d_prepare (csq1d, model);

    for (j = 0; j < tau_array->len; j++)
    {
      const gdouble tau = g_array_index (tau_array, gdouble, j);
      gdouble alpha, gamma;

      ncm_csq1d_eval_at (csq1d, model, tau, &state);
      ncm_csq1d_state_get_ag (&state, &alpha, &gamma);

      ncm_matrix_set (powspec_alpha, i, j, alpha);
      ncm_matrix_set (powspec_gamma, i, j, gamma);
    }
  }

  {
    NcmVector *tau_vec = ncm_vector_new_array (tau_array);
    NcmVector *k_vec   = ncm_vector_new_array (k_array);

    ncm_spline2d_clear (&adiab->powspec_alpha);
    ncm_spline2d_clear (&adiab->powspec_gamma);

    adiab->powspec_alpha = ncm_spline2d_bicubic_notaknot_new ();
    adiab->powspec_gamma = ncm_spline2d_bicubic_notaknot_new ();

    ncm_spline2d_set (adiab->powspec_alpha, tau_vec, k_vec, powspec_alpha, TRUE);
    ncm_spline2d_set (adiab->powspec_gamma, tau_vec, k_vec, powspec_gamma, TRUE);

    ncm_matrix_free (powspec_alpha);
    ncm_matrix_free (powspec_gamma);
    ncm_vector_free (tau_vec);
    ncm_vector_free (k_vec);
  }
}

static NcmPowspecSpline2d *
_nc_hipert_adiab_eval_powspec_func (NcHIPertAdiab *adiab, NcmModel *model,
                                    gdouble (*eval_from_state)(NcHIPertAdiab *adiab, NcmModel *model, NcmCSQ1DState *state, const gdouble k))
{
  if (!ncm_spline2d_is_init (adiab->powspec_alpha) || !ncm_spline2d_is_init (adiab->powspec_gamma))
  {
    g_error ("Power spectrum not prepared.");

    return NULL;
  }
  else
  {
    NcmVector *tau_vec     = ncm_spline2d_peek_xv (adiab->powspec_alpha);
    NcmVector *k_vec       = ncm_spline2d_peek_yv (adiab->powspec_alpha);
    const guint tau_len    = ncm_vector_len (tau_vec);
    const guint k_len      = ncm_vector_len (k_vec);
    NcmMatrix *powspec_mat = ncm_matrix_new (k_len, tau_len);
    guint i;

    for (i = 0; i < tau_len; i++)
    {
      const gdouble tau = ncm_vector_get (tau_vec, i);
      guint j;

      for (j = 0; j < k_len; j++)
      {
        const gdouble k = ncm_vector_get (k_vec, j);
        NcmCSQ1DState state;

        state.t     = tau;
        state.frame = NCM_CSQ1D_FRAME_ORIG;
        state.alpha = ncm_spline2d_eval (adiab->powspec_alpha, tau, k);
        state.gamma = ncm_spline2d_eval (adiab->powspec_gamma, tau, k);

        ncm_matrix_set (powspec_mat, j, i, log (eval_from_state (adiab, model, &state, k)));
      }
    }

    {
      NcmSpline2d *powspec_spline = ncm_spline2d_bicubic_notaknot_new ();
      NcmVector *lnk_vec          = ncm_vector_dup (k_vec);
      NcmPowspecSpline2d *powspec;
      guint i;

      for (i = 0; i < k_len; i++)
        ncm_vector_set (lnk_vec, i, log (ncm_vector_get (lnk_vec, i)));

      ncm_spline2d_set (powspec_spline, tau_vec, lnk_vec, powspec_mat, TRUE);

      powspec = ncm_powspec_spline2d_new (powspec_spline);

      ncm_matrix_free (powspec_mat);

      return powspec;
    }
  }
}

/**
 * nc_hipert_adiab_eval_powspec_zeta:
 * @adiab: a #NcHIPertAdiab
 * @model: a #NcmModel
 *
 * Evaluates the power spectrum for the gauge invariant variable $\zeta$.
 *
 * Returns: (transfer full): the power spectrum of $\zeta$.
 */
NcmPowspecSpline2d *
nc_hipert_adiab_eval_powspec_zeta (NcHIPertAdiab *adiab, NcmModel *model)
{
  return _nc_hipert_adiab_eval_powspec_func (adiab, model, _nc_hipert_eval_powspec_zeta_from_state);
}

/**
 * nc_hipert_adiab_eval_powspec_Psi:
 * @adiab: a #NcHIPertAdiab
 * @model: a #NcmModel
 *
 * Evaluates the power spectrum for the gauge invariant variable $\Psi$.
 *
 * Returns: (transfer full): the power spectrum of $\Psi$.
 */
NcmPowspecSpline2d *
nc_hipert_adiab_eval_powspec_Psi (NcHIPertAdiab *adiab, NcmModel *model)
{
  return _nc_hipert_adiab_eval_powspec_func (adiab, model, _nc_hipert_eval_powspec_Psi_from_state);
}

/**
 * nc_hipert_adiab_eval_powspec_drho:
 * @adiab: a #NcHIPertAdiab
 * @model: a #NcmModel
 *
 * Evaluates the power spectrum for the gauge invariant variable $\delta\rho$.
 *
 * Returns: (transfer full): the power spectrum of $\delta\rho$.
 */
NcmPowspecSpline2d *
nc_hipert_adiab_eval_powspec_drho (NcHIPertAdiab *adiab, NcmModel *model)
{
  return _nc_hipert_adiab_eval_powspec_func (adiab, model, _nc_hipert_eval_powspec_drho_from_state);
}

