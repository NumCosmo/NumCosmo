/***************************************************************************
 *            nc_hcosmo_bi.c
 *
 *  Sat Oct 17 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hcosmo_bi.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcHCosmoBI:
 *
 * Bianchi I cosmology with pressureless matter and a cosmological constant.
 *
 * The metric is $\mathrm{d}s^2 = -\mathrm{d}t^2 + \sum_i a_i^2(t)\,\mathrm{d}x_i^2$ with
 * $a_i(t_0) = 1$. The parameters are the Hubble rates of the three axes today, $H_i$,
 * $\Omega_{m0}$ and the orientation parameters $u_i$. The mean Hubble rate today, the shear today
 * and its density parameter are
 * $$
 * H_0 \equiv \frac{H_1 + H_2 + H_3}{3},
 * \qquad
 * s_i \equiv \frac{H_i}{H_0} - 1,
 * \qquad
 * \Omega_{s0} \equiv \frac{1}{6}\sum_i s_i^2,
 * $$
 * and the constraint equation gives $\Omega_{\Lambda0} = 1 - \Omega_{m0} - \Omega_{s0}$.
 * The model requires $\Omega_{\Lambda0} > 0$: ncm_model_params_valid() returns %FALSE
 * otherwise.
 *
 * The time variable is the lookback time in units of $1/H_0$, $t_{lb} \equiv H_0 (t_0 - t)$,
 * as in nc_distance_lookback_time(). With $\kappa \equiv 3\sqrt{\Omega_{\Lambda0}}$,
 * $w \equiv e^{-\kappa t_{lb}}$ and
 * $$
 * D \equiv \frac{\Omega_{m0}}{2\Omega_{\Lambda0}},
 * \qquad
 * A \equiv \frac{1}{2}\left(1 + D + \frac{1}{\sqrt{\Omega_{\Lambda0}}}\right),
 * \qquad
 * B \equiv \frac{1}{2}\left(1 + D - \frac{1}{\sqrt{\Omega_{\Lambda0}}}\right),
 * $$
 * the volume $V \equiv a_1 a_2 a_3$, the mean Hubble rate $E \equiv H / H_0$ and the scale
 * factors are
 * $$
 * V = \frac{A w^2 - D w + B}{w} = \frac{A (w - w_+)(w - w_-)}{w},
 * \qquad
 * E = \sqrt{\Omega_{\Lambda0}}\left(\frac{w}{w - w_+} + \frac{w_-}{w - w_-}\right),
 * \qquad
 * \ln a_i = \frac{1}{3}\ln V - s_i J,
 * $$
 * where $\Delta \equiv \sqrt{\Omega_{s0}/\Omega_{\Lambda0}}$, $w_\pm \equiv (D \pm \Delta)/(2A)$
 * and $J \equiv \int_0^{t_{lb}} \mathrm{d}t_{lb}^\prime / V$ is
 * $$
 * J = \frac{\ln(1 + \Delta g)}{\kappa\Delta},
 * \qquad
 * g \equiv \frac{2(1 - w)}{(Q + \Delta)(w - w_+)},
 * \qquad
 * Q \equiv 2A - D.
 * $$
 * The code evaluates $\ln(1 + \Delta g)/(\Delta g)$, so $J = g/\kappa$ at $\Omega_{s0} = 0$.
 * The Hubble rates of the axes are $H_i / H_0 = E + s_i / V$.
 *
 * The orientation parameters $u_i \in [0, 1]$ give the unit quaternion
 * $$
 * q = \left(\sqrt{u_1}\cos 2\pi u_3;\; \sqrt{1 - u_1}\sin 2\pi u_2,\;
 * \sqrt{1 - u_1}\cos 2\pi u_2,\; \sqrt{u_1}\sin 2\pi u_3\right),
 * $$
 * returned by nc_hcosmo_bi_get_quaternion(); a unit vector $\hat n$ in the equatorial frame
 * has components $q\,\hat n\,q^\dagger$ along the three axes. The uniform distribution
 * on the box $[0, 1]^3$ gives the uniform distribution on rotations, and the defaults
 * $(1, 0, 0)$ give the identity. The same quaternion in hyperspherical angles,
 * $q = (\cos\alpha;\; \sin\alpha\cos\beta,\; \sin\alpha\sin\beta\cos\gamma,\;
 * \sin\alpha\sin\beta\sin\gamma)$ with $\alpha \in [0, \pi/2]$, $\beta \in [0, \pi]$,
 * $\gamma \in [0, 2\pi)$, is read by nc_hcosmo_bi_get_angles() and set by
 * nc_hcosmo_bi_set_angles(); $u_1 = a^2 + d^2$, $2\pi u_2 = \arg(c + ib)$ and
 * $2\pi u_3 = \arg(a + id)$ for $q = (a; b, c, d)$, see nc_hcosmo_bi_set_quaternion().
 * The parametrization, the uniform measure on rotations and the symmetry of the axes are
 * described in <a href="../../theory/nc/background/hcosmo_bi.html">Bianchi I Orientation</a>.
 *
 * The volume vanishes at $w = w_+$, that is at $t_{lb} = H_0 t_0 = -\kappa^{-1}\ln w_+$,
 * returned by nc_hcosmo_bi_age(); nc_hcosmo_bi_ln_a() and nc_hcosmo_bi_E() require $t_{lb}$
 * below it. The code computes $w - w_+ = w_+\,\mathrm{expm1}[\kappa(H_0 t_0 - t_{lb})]$, so
 * the relative error of $V$ near $t_{lb} = H_0 t_0$ is the one set by rounding $t_{lb}$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/background/nc_hcosmo_bi.h"
#include "ncm/core/ncm_c.h"
#include "ncm/model/ncm_mset_func_list.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcHCosmoBI
{
  /*< private >*/
  NcmModel parent_instance;
};

enum
{
  PROP_0,
  PROP_SIZE,
};

typedef struct _NcHCosmoBIConst
{
  gdouble s[3];
  gdouble sqrt_Omega_L0;
  gdouble kappa;
  gdouble A;
  gdouble B;
  gdouble D;
  gdouble Delta;
  gdouble Q;
  gdouble w_p;
  gdouble w_m;
  gdouble ln_w_p;
} NcHCosmoBIConst;

G_DEFINE_TYPE (NcHCosmoBI, nc_hcosmo_bi, NCM_TYPE_MODEL)

#define VECTOR  (NCM_MODEL (cosmo))
#define H1     (ncm_model_orig_param_get (VECTOR, NC_HCOSMO_BI_H1))
#define H2     (ncm_model_orig_param_get (VECTOR, NC_HCOSMO_BI_H2))
#define H3     (ncm_model_orig_param_get (VECTOR, NC_HCOSMO_BI_H3))
#define OMEGA_M (ncm_model_orig_param_get (VECTOR, NC_HCOSMO_BI_OMEGA_M))
#define U1     (ncm_model_orig_param_get (VECTOR, NC_HCOSMO_BI_U1))
#define U2     (ncm_model_orig_param_get (VECTOR, NC_HCOSMO_BI_U2))
#define U3     (ncm_model_orig_param_get (VECTOR, NC_HCOSMO_BI_U3))

static gboolean _nc_hcosmo_bi_valid (NcmModel *model);
static void _nc_hcosmo_bi_const (NcHCosmoBI *cosmo, NcHCosmoBIConst *c);

static void
nc_hcosmo_bi_init (NcHCosmoBI *cosmo)
{
}

static void
_nc_hcosmo_bi_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_HCOSMO_BI (object));

  switch (prop_id)
  {
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_hcosmo_bi_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_HCOSMO_BI (object));

  switch (prop_id)
  {
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_hcosmo_bi_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hcosmo_bi_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_hcosmo_bi, NC_TYPE_HCOSMO_BI);

static void
nc_hcosmo_bi_class_init (NcHCosmoBIClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_hcosmo_bi_set_property;
  model_class->get_property = &_nc_hcosmo_bi_get_property;
  object_class->finalize    = &_nc_hcosmo_bi_finalize;

  ncm_model_class_set_name_nick (model_class, "Bianchi I, matter and cosmological constant", "BI");
  ncm_model_class_add_params (model_class, NC_HCOSMO_BI_SPARAM_LEN, 0, PROP_SIZE);

  ncm_mset_model_register_id (model_class,
                              "NcHCosmoBI",
                              "Bianchi I homogeneous cosmology.",
                              NULL,
                              FALSE,
                              NCM_MSET_MODEL_MAIN);

  ncm_model_class_set_sparam (model_class, NC_HCOSMO_BI_H1, "H_1", "H1",
                              10.0, 500.0, 1.0,
                              NC_HCOSMO_BI_DEFAULT_PARAMS_ABSTOL, NC_HCOSMO_BI_DEFAULT_H1,
                              NCM_PARAM_TYPE_FREE);
  ncm_model_class_set_sparam (model_class, NC_HCOSMO_BI_H2, "H_2", "H2",
                              10.0, 500.0, 1.0,
                              NC_HCOSMO_BI_DEFAULT_PARAMS_ABSTOL, NC_HCOSMO_BI_DEFAULT_H2,
                              NCM_PARAM_TYPE_FREE);
  ncm_model_class_set_sparam (model_class, NC_HCOSMO_BI_H3, "H_3", "H3",
                              10.0, 500.0, 1.0,
                              NC_HCOSMO_BI_DEFAULT_PARAMS_ABSTOL, NC_HCOSMO_BI_DEFAULT_H3,
                              NCM_PARAM_TYPE_FREE);
  ncm_model_class_set_sparam (model_class, NC_HCOSMO_BI_OMEGA_M, "\\Omega_{m0}", "Omegam",
                              1.0e-8, 1.0, 1.0e-2,
                              NC_HCOSMO_BI_DEFAULT_PARAMS_ABSTOL, NC_HCOSMO_BI_DEFAULT_OMEGA_M,
                              NCM_PARAM_TYPE_FREE);
  ncm_model_class_set_sparam (model_class, NC_HCOSMO_BI_U1, "u_1", "u1",
                              0.0, 1.0, 1.0e-1,
                              NC_HCOSMO_BI_DEFAULT_PARAMS_ABSTOL, NC_HCOSMO_BI_DEFAULT_U1,
                              NCM_PARAM_TYPE_FREE);
  ncm_model_class_set_sparam (model_class, NC_HCOSMO_BI_U2, "u_2", "u2",
                              0.0, 1.0, 1.0e-1,
                              NC_HCOSMO_BI_DEFAULT_PARAMS_ABSTOL, NC_HCOSMO_BI_DEFAULT_U2,
                              NCM_PARAM_TYPE_FREE);
  ncm_model_class_set_sparam (model_class, NC_HCOSMO_BI_U3, "u_3", "u3",
                              0.0, 1.0, 1.0e-1,
                              NC_HCOSMO_BI_DEFAULT_PARAMS_ABSTOL, NC_HCOSMO_BI_DEFAULT_U3,
                              NCM_PARAM_TYPE_FREE);

  ncm_model_class_check_params_info (model_class);

  model_class->valid = &_nc_hcosmo_bi_valid;
}

static gboolean
_nc_hcosmo_bi_valid (NcmModel *model)
{
  NcHCosmoBI *cosmo = NC_HCOSMO_BI (model);

  return nc_hcosmo_bi_Omega_L0 (cosmo) > 0.0;
}

/**
 * nc_hcosmo_bi_new:
 *
 * Creates a new #NcHCosmoBI with the default parameters.
 *
 * Returns: (transfer full): a new #NcHCosmoBI
 */
NcHCosmoBI *
nc_hcosmo_bi_new (void)
{
  NcHCosmoBI *cosmo = g_object_new (NC_TYPE_HCOSMO_BI,
                                    NULL);

  return cosmo;
}

/**
 * nc_hcosmo_bi_ref:
 * @cosmo: a #NcHCosmoBI
 *
 * Increases the reference count of @cosmo by one.
 *
 * Returns: (transfer full): @cosmo
 */
NcHCosmoBI *
nc_hcosmo_bi_ref (NcHCosmoBI *cosmo)
{
  return g_object_ref (cosmo);
}

/**
 * nc_hcosmo_bi_free:
 * @cosmo: a #NcHCosmoBI
 *
 * Decreases the reference count of @cosmo by one.
 *
 */
void
nc_hcosmo_bi_free (NcHCosmoBI *cosmo)
{
  g_object_unref (cosmo);
}

/**
 * nc_hcosmo_bi_clear:
 * @cosmo: a #NcHCosmoBI
 *
 * If *@cosmo is not %NULL, decreases its reference count and sets *@cosmo to %NULL.
 *
 */
void
nc_hcosmo_bi_clear (NcHCosmoBI **cosmo)
{
  g_clear_object (cosmo);
}

/**
 * nc_hcosmo_bi_H0:
 * @cosmo: a #NcHCosmoBI
 *
 * Returns: the mean Hubble rate today $H_0 = (H_1 + H_2 + H_3)/3$ [km/s/Mpc]
 */
gdouble
nc_hcosmo_bi_H0 (NcHCosmoBI *cosmo)
{
  return (H1 + H2 + H3) / 3.0;
}

/**
 * nc_hcosmo_bi_RH_Mpc:
 * @cosmo: a #NcHCosmoBI
 *
 * Returns: the Hubble radius of the mean Hubble rate today, $c / H_0$ [Mpc]
 */
gdouble
nc_hcosmo_bi_RH_Mpc (NcHCosmoBI *cosmo)
{
  return ncm_c_c () / (1.0e3 * nc_hcosmo_bi_H0 (cosmo));
}

/**
 * nc_hcosmo_bi_Omega_s0:
 * @cosmo: a #NcHCosmoBI
 *
 * Returns: the shear density parameter today $\Omega_{s0}$
 */
gdouble
nc_hcosmo_bi_Omega_s0 (NcHCosmoBI *cosmo)
{
  const gdouble H0 = nc_hcosmo_bi_H0 (cosmo);
  const gdouble s1 = (H1 - H0) / H0;
  const gdouble s2 = (H2 - H0) / H0;
  const gdouble s3 = (H3 - H0) / H0;

  return (s1 * s1 + s2 * s2 + s3 * s3) / 6.0;
}

/**
 * nc_hcosmo_bi_Omega_L0:
 * @cosmo: a #NcHCosmoBI
 *
 * Returns: the cosmological constant density parameter today
 *   $\Omega_{\Lambda0} = 1 - \Omega_{m0} - \Omega_{s0}$
 */
gdouble
nc_hcosmo_bi_Omega_L0 (NcHCosmoBI *cosmo)
{
  return 1.0 - OMEGA_M - nc_hcosmo_bi_Omega_s0 (cosmo);
}

/**
 * nc_hcosmo_bi_age:
 * @cosmo: a #NcHCosmoBI
 *
 * Computes the age of the universe in units of $1/H_0$, the lookback time at which the
 * volume vanishes.
 *
 * Returns: $H_0 t_0$
 */
gdouble
nc_hcosmo_bi_age (NcHCosmoBI *cosmo)
{
  NcHCosmoBIConst c;

  _nc_hcosmo_bi_const (cosmo, &c);

  return -c.ln_w_p / c.kappa;
}

/**
 * nc_hcosmo_bi_ln_a:
 * @cosmo: a #NcHCosmoBI
 * @t_lb: lookback time $t_{lb}$ in units of $1/H_0$
 * @ln_a1: (out): $\ln a_1$
 * @ln_a2: (out): $\ln a_2$
 * @ln_a3: (out): $\ln a_3$
 *
 * Computes the logarithms of the three scale factors at @t_lb, which must be smaller than
 * nc_hcosmo_bi_age().
 *
 */
void
nc_hcosmo_bi_ln_a (NcHCosmoBI *cosmo, const gdouble t_lb, gdouble *ln_a1, gdouble *ln_a2, gdouble *ln_a3)
{
  NcHCosmoBIConst c;

  _nc_hcosmo_bi_const (cosmo, &c);

  {
    const gdouble one_m_w = -expm1 (-c.kappa * t_lb);
    const gdouble w_m_w_p = c.w_p * expm1 (-c.ln_w_p - c.kappa * t_lb);
    const gdouble w_m_w_m = w_m_w_p + (c.w_p - c.w_m);
    const gdouble g       = 2.0 * one_m_w / ((c.Q + c.Delta) * w_m_w_p);
    const gdouble u       = c.Delta * g;
    const gdouble J       = g * ((u != 0.0) ? log1p (u) / u : 1.0) / c.kappa;
    const gdouble ln_V_3  = (log (c.A * w_m_w_p * w_m_w_m) + c.kappa * t_lb) / 3.0;

    ln_a1[0] = ln_V_3 - c.s[0] * J;
    ln_a2[0] = ln_V_3 - c.s[1] * J;
    ln_a3[0] = ln_V_3 - c.s[2] * J;
  }
}

/**
 * nc_hcosmo_bi_E:
 * @cosmo: a #NcHCosmoBI
 * @t_lb: lookback time $t_{lb}$ in units of $1/H_0$
 * @E1: (out): $H_1(t) / H_0$
 * @E2: (out): $H_2(t) / H_0$
 * @E3: (out): $H_3(t) / H_0$
 *
 * Computes the Hubble rates of the three axes at @t_lb in units of $H_0$; @t_lb must be
 * smaller than nc_hcosmo_bi_age().
 *
 */
void
nc_hcosmo_bi_E (NcHCosmoBI *cosmo, const gdouble t_lb, gdouble *E1, gdouble *E2, gdouble *E3)
{
  NcHCosmoBIConst c;

  _nc_hcosmo_bi_const (cosmo, &c);

  {
    const gdouble w       = exp (-c.kappa * t_lb);
    const gdouble w_m_w_p = c.w_p * expm1 (-c.ln_w_p - c.kappa * t_lb);
    const gdouble w_m_w_m = w_m_w_p + (c.w_p - c.w_m);
    const gdouble E       = c.sqrt_Omega_L0 * (w / w_m_w_p + c.w_m / w_m_w_m);
    const gdouble w_V     = w / (c.A * w_m_w_p * w_m_w_m);

    E1[0] = E + c.s[0] * w_V;
    E2[0] = E + c.s[1] * w_V;
    E3[0] = E + c.s[2] * w_V;
  }
}

/**
 * nc_hcosmo_bi_get_quaternion:
 * @cosmo: a #NcHCosmoBI
 * @q: (out caller-allocates): a #NcmQuaternion
 *
 * Sets @q to the unit quaternion of the orientation parameters $u_i$, see #NcHCosmoBI.
 *
 */
void
nc_hcosmo_bi_get_quaternion (NcHCosmoBI *cosmo, NcmQuaternion *q)
{
  const gdouble r1 = sqrt (U1);
  const gdouble r2 = sqrt (1.0 - U1);

  q->s      = r1 * cos (2.0 * M_PI * U3);
  q->v.c[0] = r2 * sin (2.0 * M_PI * U2);
  q->v.c[1] = r2 * cos (2.0 * M_PI * U2);
  q->v.c[2] = r1 * sin (2.0 * M_PI * U3);
}

/**
 * nc_hcosmo_bi_set_quaternion:
 * @cosmo: a #NcHCosmoBI
 * @q: a unit #NcmQuaternion
 *
 * Sets the orientation parameters $u_i$ to those of @q, see #NcHCosmoBI; $q$ and $-q$ give
 * different points of the box with the same rotation. The result is in $[0, 1]^3$, which may
 * lie outside the fundamental domain when the parameter bounds restrict it.
 *
 */
void
nc_hcosmo_bi_set_quaternion (NcHCosmoBI *cosmo, const NcmQuaternion *q)
{
  const gdouble a  = q->s;
  const gdouble b  = q->v.c[0];
  const gdouble c  = q->v.c[1];
  const gdouble d  = q->v.c[2];
  const gdouble u1 = a * a + d * d;
  gdouble u2       = atan2 (b, c) / (2.0 * M_PI);
  gdouble u3       = atan2 (d, a) / (2.0 * M_PI);

  u2 -= floor (u2);
  u3 -= floor (u3);

  ncm_model_orig_param_set (VECTOR, NC_HCOSMO_BI_U1, u1);
  ncm_model_orig_param_set (VECTOR, NC_HCOSMO_BI_U2, u2);
  ncm_model_orig_param_set (VECTOR, NC_HCOSMO_BI_U3, u3);
}

/**
 * nc_hcosmo_bi_get_angles:
 * @cosmo: a #NcHCosmoBI
 * @alpha: (out): $\alpha \in [0, \pi/2]$
 * @beta: (out): $\beta \in [0, \pi]$
 * @gamma: (out): $\gamma \in [0, 2\pi)$
 *
 * Computes the hyperspherical angles of the orientation quaternion, see #NcHCosmoBI, after
 * changing its sign if needed so that its scalar part is not negative.
 *
 */
void
nc_hcosmo_bi_get_angles (NcHCosmoBI *cosmo, gdouble *alpha, gdouble *beta, gdouble *gamma)
{
  NcmQuaternion q;
  gdouble sign;

  nc_hcosmo_bi_get_quaternion (cosmo, &q);
  sign = (q.s < 0.0) ? -1.0 : 1.0;

  {
    const gdouble a = sign * q.s;
    const gdouble b = sign * q.v.c[0];
    const gdouble c = sign * q.v.c[1];
    const gdouble d = sign * q.v.c[2];
    gdouble g       = atan2 (d, c);

    if (g < 0.0)
      g += 2.0 * M_PI;

    alpha[0] = acos (CLAMP (a, -1.0, 1.0));
    beta[0]  = atan2 (hypot (c, d), b);
    gamma[0] = g;
  }
}

/**
 * nc_hcosmo_bi_set_angles:
 * @cosmo: a #NcHCosmoBI
 * @alpha: $\alpha$
 * @beta: $\beta$
 * @gamma: $\gamma$
 *
 * Sets the orientation parameters from the hyperspherical angles of the quaternion, see
 * #NcHCosmoBI and nc_hcosmo_bi_set_quaternion().
 *
 */
void
nc_hcosmo_bi_set_angles (NcHCosmoBI *cosmo, const gdouble alpha, const gdouble beta, const gdouble gamma)
{
  NcmQuaternion q;

  q.s      = cos (alpha);
  q.v.c[0] = sin (alpha) * cos (beta);
  q.v.c[1] = sin (alpha) * sin (beta) * cos (gamma);
  q.v.c[2] = sin (alpha) * sin (beta) * sin (gamma);

  nc_hcosmo_bi_set_quaternion (cosmo, &q);
}

#define _NC_HCOSMO_BI_FUNC0_TO_FLIST(fname)                                                                                   \
        static void _nc_hcosmo_bi_flist_ ## fname (NcmMSetFuncList * flist, NcmMSet * mset, const gdouble * x, gdouble * res) \
        {                                                                                                                     \
          NcHCosmoBI *cosmo = NC_HCOSMO_BI (ncm_mset_peek (mset, nc_hcosmo_bi_id ()));                                        \
          res[0] = nc_hcosmo_bi_ ## fname (cosmo);                                                                            \
        }

_NC_HCOSMO_BI_FUNC0_TO_FLIST (H0)
_NC_HCOSMO_BI_FUNC0_TO_FLIST (RH_Mpc)
_NC_HCOSMO_BI_FUNC0_TO_FLIST (Omega_s0)
_NC_HCOSMO_BI_FUNC0_TO_FLIST (Omega_L0)
_NC_HCOSMO_BI_FUNC0_TO_FLIST (age)

static void
_nc_hcosmo_bi_flist_H2_minus_H1 (NcmMSetFuncList *flist, NcmMSet *mset, const gdouble *x, gdouble *res)
{
  NcHCosmoBI *cosmo = NC_HCOSMO_BI (ncm_mset_peek (mset, nc_hcosmo_bi_id ()));

  res[0] = H2 - H1;
}

static void
_nc_hcosmo_bi_flist_H3_minus_H2 (NcmMSetFuncList *flist, NcmMSet *mset, const gdouble *x, gdouble *res)
{
  NcHCosmoBI *cosmo = NC_HCOSMO_BI (ncm_mset_peek (mset, nc_hcosmo_bi_id ()));

  res[0] = H3 - H2;
}

void _nc_hcosmo_bi_register_functions (void);

void
_nc_hcosmo_bi_register_functions (void)
{
  ncm_mset_func_list_register ("H0",          "H_0",                   "NcHCosmoBI", "Mean Hubble rate today",           G_TYPE_NONE, _nc_hcosmo_bi_flist_H0,          0, 1);
  ncm_mset_func_list_register ("RH_Mpc",      "R_{H} [\\mathrm{Mpc}]", "NcHCosmoBI", "Hubble radius today in Mpc",       G_TYPE_NONE, _nc_hcosmo_bi_flist_RH_Mpc,      0, 1);
  ncm_mset_func_list_register ("Omega_s0",    "\\Omega_{s0}",          "NcHCosmoBI", "Shear density parameter today",    G_TYPE_NONE, _nc_hcosmo_bi_flist_Omega_s0,    0, 1);
  ncm_mset_func_list_register ("Omega_L0",    "\\Omega_{\\Lambda0}",   "NcHCosmoBI", "Cosmological constant today",      G_TYPE_NONE, _nc_hcosmo_bi_flist_Omega_L0,    0, 1);
  ncm_mset_func_list_register ("age",         "H_0 t_0",               "NcHCosmoBI", "Age in units of 1/H_0",            G_TYPE_NONE, _nc_hcosmo_bi_flist_age,         0, 1);
  ncm_mset_func_list_register ("H2_minus_H1", "H_2 - H_1",             "NcHCosmoBI", "Difference of the axis 2 and 1 Hubble rates", G_TYPE_NONE, _nc_hcosmo_bi_flist_H2_minus_H1, 0, 1);
  ncm_mset_func_list_register ("H3_minus_H2", "H_3 - H_2",             "NcHCosmoBI", "Difference of the axis 3 and 2 Hubble rates", G_TYPE_NONE, _nc_hcosmo_bi_flist_H3_minus_H2, 0, 1);
}

/*
 * Computes s_i, sqrt(Omega_L0), kappa, A, B, D, Delta, Q, w_+, w_- and ln w_+ of the
 * NcHCosmoBI description from the current parameters and writes them to @c.
 */
static void
_nc_hcosmo_bi_const (NcHCosmoBI *cosmo, NcHCosmoBIConst *c)
{
  const gdouble H0       = nc_hcosmo_bi_H0 (cosmo);
  const gdouble Omega_s0 = nc_hcosmo_bi_Omega_s0 (cosmo);
  const gdouble Omega_L0 = 1.0 - OMEGA_M - Omega_s0;

  c->s[0]          = (H1 - H0) / H0;
  c->s[1]          = (H2 - H0) / H0;
  c->s[2]          = (H3 - H0) / H0;
  c->sqrt_Omega_L0 = sqrt (Omega_L0);
  c->kappa         = 3.0 * c->sqrt_Omega_L0;
  c->D             = OMEGA_M / (2.0 * Omega_L0);
  c->A             = 0.5 * (1.0 + c->D + 1.0 / c->sqrt_Omega_L0);
  c->B             = 0.5 * (1.0 + c->D - 1.0 / c->sqrt_Omega_L0);
  c->Delta         = sqrt (Omega_s0 / Omega_L0);
  c->Q             = 2.0 * c->A - c->D;
  c->w_p           = (c->D + c->Delta) / (2.0 * c->A);
  c->w_m           = (c->D - c->Delta) / (2.0 * c->A);
  c->ln_w_p        = log (c->w_p);
}

