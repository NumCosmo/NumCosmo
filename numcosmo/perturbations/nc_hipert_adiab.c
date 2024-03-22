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
 * SECTION:nc_hipert_adiab
 * @title: NcHIPertAdiab
 * @short_description: Perturbation object for adiabatic mode only.
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

#include "math/ncm_spline_cubic_notaknot.h"
#include "perturbations/nc_hipert_adiab.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <gsl/gsl_roots.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcHIPertAdiab
{
  NcmCSQ1D parent_instance;
  gdouble k;
};

G_DEFINE_INTERFACE (NcHIPertIAdiab, nc_hipert_iadiab, G_TYPE_OBJECT)
G_DEFINE_TYPE (NcHIPertAdiab, nc_hipert_adiab, NCM_TYPE_CSQ1D)

static void
nc_hipert_iadiab_default_init (NcHIPertIAdiabInterface *iface)
{
  iface->eval_xi     = NULL;
  iface->eval_F1     = NULL;
  iface->eval_nu     = NULL;
  iface->eval_m      = NULL;
  iface->eval_unit   = NULL;
  iface->eval_x      = NULL;
  iface->eval_p2Psi  = NULL;
  iface->eval_p2drho = NULL;
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

static void
nc_hipert_adiab_init (NcHIPertAdiab *pa)
{
  pa->k = 0.0;
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
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
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
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_adiab_dispose (GObject *object)
{
  /*NcHIPertAdiab *pa = NC_HIPERT_ADIAB (object);*/

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
static gdouble _nc_hipert_adiab_eval_unit (NcmCSQ1D *csq1d, NcmModel *model);

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


  csq1d_class->eval_xi   = &_nc_hipert_adiab_eval_xi;
  csq1d_class->eval_F1   = &_nc_hipert_adiab_eval_F1;
  csq1d_class->eval_nu   = &_nc_hipert_adiab_eval_nu;
  csq1d_class->eval_m    = &_nc_hipert_adiab_eval_m;
  csq1d_class->prepare   = &_nc_hipert_adiab_prepare;
  csq1d_class->eval_unit = &_nc_hipert_adiab_eval_unit;
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

static gdouble
_nc_hipert_adiab_eval_unit (NcmCSQ1D *csq1d, NcmModel *model)
{
  NcHIPertAdiab *pa = NC_HIPERT_ADIAB (csq1d);
  const gdouble k   = pa->k;

  return nc_hipert_iadiab_eval_unit (NC_HIPERT_IADIAB (model));
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

/**
 * nc_hipert_adiab_eval_powspec_zeta_at:
 * @adiab: a #NcHIPertAdiab
 * @model: a #NcmModel
 * @tau: $\tau$
 *
 * Evaluates the power spectrum of the gauge invariant variable $\zeta$ at a given time
 * $\tau$. The power spectrum is given by
 * $$
 * P_\zeta = u^2\frac{2\pi^2}{k^3} \frac{J_{11}}{2}.
 * $$
 * where $u$ is the numerical factor for the power spectrum of the adiabatic mode, $k$ is
 * the wave number.
 *
 * Returns: the power spectrum of the gauge invariant variable $\zeta$.
 */
gdouble
nc_hipert_adiab_eval_powspec_zeta_at (NcHIPertAdiab *adiab, NcmModel *model, const gdouble tau)
{
  const gdouble unit = nc_hipert_iadiab_eval_unit (NC_HIPERT_IADIAB (model));
  const gdouble fact = ncm_c_two_pi_2 ();
  const gdouble k    = nc_hipert_adiab_get_k (adiab);
  NcmCSQ1DState state;
  gdouble J11, J12, J22;

  ncm_csq1d_eval_at (NCM_CSQ1D (adiab), tau, &state);

  ncm_csq1d_state_get_J (&state, &J11, &J12, &J22);

  return gsl_pow_2 (unit) * gsl_pow_3 (k) * (J11 / 2.0) / fact;
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
  const gdouble unit  = nc_hipert_iadiab_eval_unit (NC_HIPERT_IADIAB (model));
  const gdouble fact  = ncm_c_two_pi_2 ();
  const gdouble k     = nc_hipert_adiab_get_k (adiab);
  const gdouble p2Psi = nc_hipert_iadiab_eval_p2Psi (NC_HIPERT_IADIAB (model), tau, k);
  NcmCSQ1DState state;
  gdouble J11, J12, J22;

  ncm_csq1d_eval_at (NCM_CSQ1D (adiab), tau, &state);

  ncm_csq1d_state_get_J (&state, &J11, &J12, &J22);

  return gsl_pow_2 (unit * p2Psi) * gsl_pow_3 (k) * (J22 / 2.0) / fact;
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
  const gdouble unit   = nc_hipert_iadiab_eval_unit (NC_HIPERT_IADIAB (model));
  const gdouble fact   = ncm_c_two_pi_2 ();
  const gdouble k      = nc_hipert_adiab_get_k (adiab);
  const gdouble p2drho = nc_hipert_iadiab_eval_p2drho (NC_HIPERT_IADIAB (model), tau, k);
  NcmCSQ1DState state;
  gdouble J11, J12, J22;

  ncm_csq1d_eval_at (NCM_CSQ1D (adiab), tau, &state);

  ncm_csq1d_state_get_J (&state, &J11, &J12, &J22);

  return gsl_pow_2 (unit * p2drho) * gsl_pow_3 (k) * (J22 / 2.0) / fact;
}

/**
 * nc_hipert_adiab_eval_powspec_zeta:
 * @adiab: a #NcHIPertAdiab
 * @model: a #NcmModel
 * @tAi: $\tau_i$ initial time to search for adiabatic limit
 * @tAf: $\tau_f$ final time to search for adiabatic limit
 * @adiab_reltol: relative tolerance for the adiabatic limit search
 * @k_array: (element-type gdouble): array of wave numbers
 * @tau: $\tau$ time to evaluate the power spectrum
 *
 * Evaluates the power spectrum of the adiabatic mode at a given time $\tau$.
 *
 * Returns: (transfer full) (element-type gdouble): the power spectrum of the adiabatic mode.
 */
GArray *
nc_hipert_adiab_eval_powspec_zeta (NcHIPertAdiab *adiab, NcmModel *model, const gdouble tAi, const gdouble tAf, const gdouble adiab_reltol, GArray *k_array, const gdouble tau)
{
  NcmCSQ1D *csq1d    = NCM_CSQ1D (adiab);
  const gdouble unit = nc_hipert_iadiab_eval_unit (NC_HIPERT_IADIAB (model));
  GArray *powspec    = g_array_new (FALSE, FALSE, sizeof (gdouble));
  gdouble tAd;
  gboolean found;
  guint i;

  for (i = 0; i < k_array->len; i++)
  {
    nc_hipert_adiab_set_k (adiab, g_array_index (k_array, gdouble, i));

    if (!ncm_csq1d_find_adiab_time_limit (csq1d, model, tAi, tAf, adiab_reltol, &tAd))
    {
      g_array_free (powspec, TRUE);

      return NULL;
    }
    else
    {
      NcmCSQ1DState state;
      gdouble zeta_k_tau, alpha_reltol, dgamma_reltol;

      ncm_csq1d_compute_adiab (csq1d, model, tAd, &state, &alpha_reltol, &dgamma_reltol);

      ncm_csq1d_set_init_cond_adiab (csq1d, model, tAd);
      ncm_csq1d_prepare (csq1d, model);

      zeta_k_tau = nc_hipert_adiab_eval_powspec_zeta_at (adiab, model, tau);

      g_array_append_val (powspec, zeta_k_tau);
    }
  }

  return powspec;
}

/**
 * nc_hipert_adiab_eval_powspec_Psi:
 * @adiab: a #NcHIPertAdiab
 * @model: a #NcmModel
 * @tAi: $\tau_i$ initial time to search for adiabatic limit
 * @tAf: $\tau_f$ final time to search for adiabatic limit
 * @adiab_reltol: relative tolerance for the adiabatic limit search
 * @k_array: (element-type gdouble): array of wave numbers
 * @tau: $\tau$ time to evaluate the power spectrum
 *
 * Evaluates the power spectrum of the adiabatic mode at a given time $\tau$.
 *
 * Returns: (transfer full) (element-type gdouble): the power spectrum of the adiabatic mode.
 */
GArray *
nc_hipert_adiab_eval_powspec_Psi (NcHIPertAdiab *adiab, NcmModel *model, const gdouble tAi, const gdouble tAf, const gdouble adiab_reltol, GArray *k_array, const gdouble tau)
{
  NcmCSQ1D *csq1d    = NCM_CSQ1D (adiab);
  const gdouble unit = nc_hipert_iadiab_eval_unit (NC_HIPERT_IADIAB (model));
  GArray *powspec    = g_array_new (FALSE, FALSE, sizeof (gdouble));
  gdouble tAd;
  gboolean found;
  guint i;

  for (i = 0; i < k_array->len; i++)
  {
    nc_hipert_adiab_set_k (adiab, g_array_index (k_array, gdouble, i));

    if (!ncm_csq1d_find_adiab_time_limit (csq1d, model, tAi, tAf, adiab_reltol, &tAd))
    {
      g_array_free (powspec, TRUE);

      return NULL;
    }
    else
    {
      NcmCSQ1DState state;
      gdouble Psi_k_tau, alpha_reltol, dgamma_reltol;

      ncm_csq1d_compute_adiab (csq1d, model, tAd, &state, &alpha_reltol, &dgamma_reltol);

      ncm_csq1d_set_init_cond_adiab (csq1d, model, tAd);
      ncm_csq1d_prepare (csq1d, model);

      Psi_k_tau = nc_hipert_adiab_eval_powspec_Psi_at (adiab, model, tau);

      g_array_append_val (powspec, Psi_k_tau);
    }
  }

  return powspec;
}

/**
 * nc_hipert_adiab_eval_powspec_drho:
 * @adiab: a #NcHIPertAdiab
 * @model: a #NcmModel
 * @tAi: $\tau_i$ initial time to search for adiabatic limit
 * @tAf: $\tau_f$ final time to search for adiabatic limit
 * @adiab_reltol: relative tolerance for the adiabatic limit search
 * @k_array: (element-type gdouble): array of wave numbers
 * @tau: $\tau$ time to evaluate the power spectrum
 *
 * Evaluates the power spectrum of the adiabatic mode at a given time $\tau$.
 *
 * Returns: (transfer full) (element-type gdouble): the power spectrum of the adiabatic mode.
 */
GArray *
nc_hipert_adiab_eval_powspec_drho (NcHIPertAdiab *adiab, NcmModel *model, const gdouble tAi, const gdouble tAf, const gdouble adiab_reltol, GArray *k_array, const gdouble tau)
{
  NcmCSQ1D *csq1d    = NCM_CSQ1D (adiab);
  const gdouble unit = nc_hipert_iadiab_eval_unit (NC_HIPERT_IADIAB (model));
  GArray *powspec    = g_array_new (FALSE, FALSE, sizeof (gdouble));
  gdouble tAd;
  gboolean found;
  guint i;

  for (i = 0; i < k_array->len; i++)
  {
    nc_hipert_adiab_set_k (adiab, g_array_index (k_array, gdouble, i));

    if (!ncm_csq1d_find_adiab_time_limit (csq1d, model, tAi, tAf, adiab_reltol, &tAd))
    {
      g_array_free (powspec, TRUE);

      return NULL;
    }
    else
    {
      NcmCSQ1DState state;
      gdouble drho_k_tau, alpha_reltol, dgamma_reltol;

      ncm_csq1d_compute_adiab (csq1d, model, tAd, &state, &alpha_reltol, &dgamma_reltol);

      ncm_csq1d_set_init_cond_adiab (csq1d, model, tAd);
      ncm_csq1d_prepare (csq1d, model);

      drho_k_tau = nc_hipert_adiab_eval_powspec_drho_at (adiab, model, tau);

      g_array_append_val (powspec, drho_k_tau);
    }
  }

  return powspec;
}

