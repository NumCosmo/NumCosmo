/***************************************************************************
 *            nc_hcosmo_bi_dist.c
 *
 *  Sat Oct 17 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hcosmo_bi_dist.c
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
 * NcHCosmoBIDist:
 *
 * Redshift and luminosity distance along a light ray in a #NcHCosmoBI.
 *
 * A photon observed today from the unit direction $\hat n$ (components $n_i$ along the
 * axes of the model) has energy $E(t) = (\sum_i n_i^2 u_i)^{1/2}$ with $u_i \equiv a_i^{-2}$,
 * so its emission lookback time $t_e$ solves $E(t_e) = 1 + z$. The ray map from observed
 * directions to source positions at $t_e$ gives the angular diameter distance
 * $$
 * D_A^2 = a_1 a_2 a_3 (1 + z) \left(n_1^2 K_{12} K_{13} + n_2^2 K_{12} K_{23} + n_3^2 K_{13} K_{23}\right),
 * \qquad
 * K_{ij} \equiv \int_0^{t_e} \frac{u_i u_j}{E^3}\,\mathrm{d}t_{lb},
 * $$
 * with $a_i$ at $t_e$, and $D_L = (1 + z)^2 D_A$, both in units of $c / H_0$.
 *
 * nc_hcosmo_bi_dist_prepare() finds the lookback time $t_\mathrm{max}$ at which every
 * direction has $E = 1 + z_\mathrm{max}$ (#NcHCosmoBIDist:zmax) or more, chooses the number
 * of Chebyshev-Lobatto nodes on $[0, t_\mathrm{max}]$ at which the integrands of $K_{ij}$
 * along seven test directions (the axes, the face diagonals and the body diagonal) have
 * converged to #NcHCosmoBIDist:reltol, together with $u_i$ and $\ln V$, $V \equiv a_1 a_2 a_3$.
 * It stores $u_i$ at these nodes and the Chebyshev coefficients of $u_i$ and $\ln V$. For
 * each ray the coefficients of $E^2 = \sum_i n_i^2 u_i$ are a linear combination of the
 * stored ones; $t_e$ comes from Newton iterations on that series, $V(t_e)$ from the series
 * of $\ln V$, and each $K_{ij}$ is the integral of the Chebyshev interpolant of its
 * integrand on the stored nodes. The evaluation functions read only these tables, so they
 * take no model; they write to scratch space of the object and may not be called
 * concurrently.
 *
 * Every axis must expand over $[0, t_\mathrm{max}]$, which makes $E$ increasing along every
 * ray and $t_e$ unique. nc_hcosmo_bi_dist_prepare() stops with an error when an axis
 * contracts before $z_\mathrm{max}$ is reached (large shear), and the evaluation functions
 * stop with an error for a redshift above #NcHCosmoBIDist:zmax.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/background/nc_hcosmo_bi_dist.h"
#include "ncm/model/ncm_model_ctrl.h"
#include "ncm/algebra/ncm_spectral.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

#define NC_HCOSMO_BI_DIST_NTEST 7

struct _NcHCosmoBIDist
{
  /*< private >*/
  GObject parent_instance;
  gdouble zmax;
  gdouble reltol;
  NcmModelCtrl *ctrl;
  NcmSpectral *spectral;
  NcmQuaternion q;
  gdouble t_lb_max;
  gdouble E0[3];
  guint n;
  GArray *u;
  GArray *cu;
  GArray *C;
  GArray *A;
  GArray *W;
  GArray *e;
};

enum
{
  PROP_0,
  PROP_ZMAX,
  PROP_RELTOL,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcHCosmoBIDist, nc_hcosmo_bi_dist, G_TYPE_OBJECT)

static void
nc_hcosmo_bi_dist_init (NcHCosmoBIDist *dist)
{
  dist->zmax     = 0.0;
  dist->reltol   = 0.0;
  dist->ctrl     = ncm_model_ctrl_new (NULL);
  dist->spectral = ncm_spectral_new ();
  dist->t_lb_max = 0.0;
  dist->n        = 0;
  dist->u        = g_array_new (FALSE, FALSE, sizeof (gdouble));
  dist->cu       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  dist->C        = g_array_new (FALSE, FALSE, sizeof (gdouble));
  dist->A        = g_array_new (FALSE, FALSE, sizeof (gdouble));
  dist->W        = g_array_new (FALSE, FALSE, sizeof (gdouble));
  dist->e        = g_array_new (FALSE, FALSE, sizeof (gdouble));

  ncm_quaternion_set_I (&dist->q);
}

static void
_nc_hcosmo_bi_dist_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHCosmoBIDist *dist = NC_HCOSMO_BI_DIST (object);

  g_return_if_fail (NC_IS_HCOSMO_BI_DIST (object));

  switch (prop_id)
  {
    case PROP_ZMAX:
      nc_hcosmo_bi_dist_set_zmax (dist, g_value_get_double (value));
      break;
    case PROP_RELTOL:
      nc_hcosmo_bi_dist_set_reltol (dist, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_hcosmo_bi_dist_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHCosmoBIDist *dist = NC_HCOSMO_BI_DIST (object);

  g_return_if_fail (NC_IS_HCOSMO_BI_DIST (object));

  switch (prop_id)
  {
    case PROP_ZMAX:
      g_value_set_double (value, dist->zmax);
      break;
    case PROP_RELTOL:
      g_value_set_double (value, dist->reltol);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_hcosmo_bi_dist_dispose (GObject *object)
{
  NcHCosmoBIDist *dist = NC_HCOSMO_BI_DIST (object);

  ncm_model_ctrl_clear (&dist->ctrl);
  ncm_spectral_clear (&dist->spectral);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hcosmo_bi_dist_parent_class)->dispose (object);
}

static void
_nc_hcosmo_bi_dist_finalize (GObject *object)
{
  NcHCosmoBIDist *dist = NC_HCOSMO_BI_DIST (object);

  g_array_unref (dist->u);
  g_array_unref (dist->cu);
  g_array_unref (dist->C);
  g_array_unref (dist->A);
  g_array_unref (dist->W);
  g_array_unref (dist->e);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hcosmo_bi_dist_parent_class)->finalize (object);
}

static void
nc_hcosmo_bi_dist_class_init (NcHCosmoBIDistClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_hcosmo_bi_dist_set_property;
  object_class->get_property = &_nc_hcosmo_bi_dist_get_property;
  object_class->dispose      = &_nc_hcosmo_bi_dist_dispose;
  object_class->finalize     = &_nc_hcosmo_bi_dist_finalize;

  /**
   * NcHCosmoBIDist:zmax:
   *
   * Largest redshift the object computes.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_ZMAX,
                                   g_param_spec_double ("zmax",
                                                        NULL,
                                                        "Maximum redshift",
                                                        0.0, G_MAXDOUBLE, 0.2,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcHCosmoBIDist:reltol:
   *
   * Relative tolerance of the Chebyshev interpolation of the integrands.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance",
                                                        GSL_DBL_EPSILON, 1.0, NC_HCOSMO_BI_DIST_DEFAULT_RELTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_hcosmo_bi_dist_new:
 * @zmax: largest redshift
 *
 * Creates a new #NcHCosmoBIDist.
 *
 * Returns: (transfer full): a new #NcHCosmoBIDist
 */
NcHCosmoBIDist *
nc_hcosmo_bi_dist_new (const gdouble zmax)
{
  NcHCosmoBIDist *dist = g_object_new (NC_TYPE_HCOSMO_BI_DIST,
                                       "zmax", zmax,
                                       NULL);

  return dist;
}

/**
 * nc_hcosmo_bi_dist_ref:
 * @dist: a #NcHCosmoBIDist
 *
 * Increases the reference count of @dist by one.
 *
 * Returns: (transfer full): @dist
 */
NcHCosmoBIDist *
nc_hcosmo_bi_dist_ref (NcHCosmoBIDist *dist)
{
  return g_object_ref (dist);
}

/**
 * nc_hcosmo_bi_dist_free:
 * @dist: a #NcHCosmoBIDist
 *
 * Decreases the reference count of @dist by one.
 *
 */
void
nc_hcosmo_bi_dist_free (NcHCosmoBIDist *dist)
{
  g_object_unref (dist);
}

/**
 * nc_hcosmo_bi_dist_clear:
 * @dist: a #NcHCosmoBIDist
 *
 * If *@dist is not %NULL, decreases its reference count and sets *@dist to %NULL.
 *
 */
void
nc_hcosmo_bi_dist_clear (NcHCosmoBIDist **dist)
{
  g_clear_object (dist);
}

/**
 * nc_hcosmo_bi_dist_set_zmax:
 * @dist: a #NcHCosmoBIDist
 * @zmax: largest redshift
 *
 * Sets #NcHCosmoBIDist:zmax; the next prepare recomputes the tables.
 *
 */
void
nc_hcosmo_bi_dist_set_zmax (NcHCosmoBIDist *dist, const gdouble zmax)
{
  if (dist->zmax != zmax)
  {
    dist->zmax = zmax;
    ncm_model_ctrl_force_update (dist->ctrl);
  }
}

/**
 * nc_hcosmo_bi_dist_get_zmax:
 * @dist: a #NcHCosmoBIDist
 *
 * Returns: #NcHCosmoBIDist:zmax
 */
gdouble
nc_hcosmo_bi_dist_get_zmax (NcHCosmoBIDist *dist)
{
  return dist->zmax;
}

/**
 * nc_hcosmo_bi_dist_set_reltol:
 * @dist: a #NcHCosmoBIDist
 * @reltol: relative tolerance
 *
 * Sets #NcHCosmoBIDist:reltol; the next prepare recomputes the tables.
 *
 */
void
nc_hcosmo_bi_dist_set_reltol (NcHCosmoBIDist *dist, const gdouble reltol)
{
  if (dist->reltol != reltol)
  {
    dist->reltol = reltol;
    ncm_model_ctrl_force_update (dist->ctrl);
  }
}

/**
 * nc_hcosmo_bi_dist_get_reltol:
 * @dist: a #NcHCosmoBIDist
 *
 * Returns: #NcHCosmoBIDist:reltol
 */
gdouble
nc_hcosmo_bi_dist_get_reltol (NcHCosmoBIDist *dist)
{
  return dist->reltol;
}

typedef struct _NcHCosmoBIDistTest
{
  NcHCosmoBI *cosmo;
  gdouble n2[NC_HCOSMO_BI_DIST_NTEST][3];
} NcHCosmoBIDistTest;

static gdouble _nc_hcosmo_bi_dist_find_t_lb_max (NcHCosmoBI *cosmo, const gdouble zmax);
static void _nc_hcosmo_bi_dist_test_integrands (gpointer user_data, gdouble t_lb, NcmVector *fval);

/**
 * nc_hcosmo_bi_dist_prepare:
 * @dist: a #NcHCosmoBIDist
 * @cosmo: a #NcHCosmoBI
 *
 * Computes the node table of @cosmo, see #NcHCosmoBIDist, and stores its orientation.
 *
 */
void
nc_hcosmo_bi_dist_prepare (NcHCosmoBIDist *dist, NcHCosmoBI *cosmo)
{
  NcHCosmoBIDistTest test = {
    cosmo,
    {
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
      {0.5, 0.5, 0.0},
      {0.5, 0.0, 0.5},
      {0.0, 0.5, 0.5},
      {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
    }
  };
  NcmMatrix *coeffs = NULL;
  guint level, n, k, m;

  dist->t_lb_max = _nc_hcosmo_bi_dist_find_t_lb_max (cosmo, dist->zmax);

  level = ncm_spectral_compute_chebyshev_coeffs_batch_adaptive (dist->spectral,
                                                                &_nc_hcosmo_bi_dist_test_integrands,
                                                                3 * NC_HCOSMO_BI_DIST_NTEST + 4,
                                                                0.0, dist->t_lb_max,
                                                                3, dist->reltol, 0.0,
                                                                &coeffs, &test);
  ncm_matrix_free (coeffs);

  n       = 1 << level;
  dist->n = n;

  g_array_set_size (dist->u, 4 * (n + 1));
  g_array_set_size (dist->cu, 4 * (n + 1));
  g_array_set_size (dist->C, (n + 1) * (n + 1));
  g_array_set_size (dist->A, n + 1);
  g_array_set_size (dist->W, n + 1);
  g_array_set_size (dist->e, n + 1);

  for (m = 0; m <= n; m++)
  {
    const gdouble t_lb = 0.5 * dist->t_lb_max * (1.0 + cos (M_PI * m / n));
    gdouble ln_a1, ln_a2, ln_a3, E1, E2, E3;

    nc_hcosmo_bi_E (cosmo, t_lb, &E1, &E2, &E3);

    if (GSL_MIN (E1, GSL_MIN (E2, E3)) <= 0.0)
      g_error ("nc_hcosmo_bi_dist_prepare: an axis contracts at t_lb = %.15g, before "
               "zmax = %.15g is reached; lower NcHCosmoBIDist:zmax.", t_lb, dist->zmax);

    nc_hcosmo_bi_ln_a (cosmo, t_lb, &ln_a1, &ln_a2, &ln_a3);
    g_array_index (dist->u, gdouble, 0 * (n + 1) + m) = exp (-2.0 * ln_a1);
    g_array_index (dist->u, gdouble, 1 * (n + 1) + m) = exp (-2.0 * ln_a2);
    g_array_index (dist->u, gdouble, 2 * (n + 1) + m) = exp (-2.0 * ln_a3);
    g_array_index (dist->u, gdouble, 3 * (n + 1) + m) = ln_a1 + ln_a2 + ln_a3;
  }

  /* C maps values at s_m = cos(pi m / n) to the coefficients of sum_k c_k T_k(s) */
  for (k = 0; k <= n; k++)
  {
    const gdouble wk = ((k == 0) || (k == n)) ? 0.5 : 1.0;

    for (m = 0; m <= n; m++)
    {
      const gdouble wm = ((m == 0) || (m == n)) ? 0.5 : 1.0;

      g_array_index (dist->C, gdouble, k * (n + 1) + m) = 2.0 / n * wk * wm * cos (M_PI * ((k * m) % (2 * n)) / n);
    }
  }

  /* Chebyshev coefficients of u_1, u_2, u_3 and ln V */
  for (k = 0; k < 4 * (n + 1); k++)
  {
    const guint f  = k / (n + 1);
    const guint kk = k % (n + 1);
    gdouble c      = 0.0;

    for (m = 0; m <= n; m++)
      c += g_array_index (dist->C, gdouble, kk * (n + 1) + m) * g_array_index (dist->u, gdouble, f * (n + 1) + m);

    g_array_index (dist->cu, gdouble, k) = c;
  }

  nc_hcosmo_bi_E (cosmo, 0.0, &dist->E0[0], &dist->E0[1], &dist->E0[2]);

  nc_hcosmo_bi_get_quaternion (cosmo, &dist->q);
  ncm_model_ctrl_update (dist->ctrl, NCM_MODEL (cosmo));
}

/**
 * nc_hcosmo_bi_dist_prepare_if_needed:
 * @dist: a #NcHCosmoBIDist
 * @cosmo: a #NcHCosmoBI
 *
 * Calls nc_hcosmo_bi_dist_prepare() if @cosmo or the properties of @dist changed since the
 * last preparation.
 *
 */
void
nc_hcosmo_bi_dist_prepare_if_needed (NcHCosmoBIDist *dist, NcHCosmoBI *cosmo)
{
  if (ncm_model_ctrl_update (dist->ctrl, NCM_MODEL (cosmo)))
    nc_hcosmo_bi_dist_prepare (dist, cosmo);
}

/**
 * nc_hcosmo_bi_dist_get_t_lb_max:
 * @dist: a #NcHCosmoBIDist
 *
 * Returns: the upper end $t_\mathrm{max}$ of the node interval
 */
gdouble
nc_hcosmo_bi_dist_get_t_lb_max (NcHCosmoBIDist *dist)
{
  return dist->t_lb_max;
}

/**
 * nc_hcosmo_bi_dist_get_n_nodes:
 * @dist: a #NcHCosmoBIDist
 *
 * Returns: the number of Chebyshev-Lobatto nodes chosen by the last prepare
 */
guint
nc_hcosmo_bi_dist_get_n_nodes (NcHCosmoBIDist *dist)
{
  return dist->n + 1;
}

/**
 * nc_hcosmo_bi_dist_sky_to_axes:
 * @dist: a #NcHCosmoBIDist
 * @px: x component in the equatorial frame
 * @py: y component in the equatorial frame
 * @pz: z component in the equatorial frame
 * @n1: (out): component along axis 1
 * @n2: (out): component along axis 2
 * @n3: (out): component along axis 3
 *
 * Rotates a direction from the equatorial frame to the axes of the model, with the
 * orientation stored by the last prepare.
 *
 */
void
nc_hcosmo_bi_dist_sky_to_axes (NcHCosmoBIDist *dist, const gdouble px, const gdouble py, const gdouble pz, gdouble *n1, gdouble *n2, gdouble *n3)
{
  NcmTriVec v = {
    {
      px, py, pz
    }
  };

  ncm_quaternion_rotate (&dist->q, &v);

  n1[0] = v.c[0];
  n2[0] = v.c[1];
  n3[0] = v.c[2];
}

static gdouble _nc_hcosmo_bi_dist_t_lb_e (NcHCosmoBIDist *dist, const gdouble z, const gdouble nsq[3], gdouble *ln_V);

/**
 * nc_hcosmo_bi_dist_t_lb_emission:
 * @dist: a #NcHCosmoBIDist
 * @z: redshift
 * @n1: direction component along axis 1
 * @n2: direction component along axis 2
 * @n3: direction component along axis 3
 *
 * Computes the emission lookback time of a photon observed with redshift @z from the
 * unit direction $(n_1, n_2, n_3)$, in units of $1/H_0$.
 *
 * Returns: $t_e$
 */
gdouble
nc_hcosmo_bi_dist_t_lb_emission (NcHCosmoBIDist *dist, const gdouble z, const gdouble n1, const gdouble n2, const gdouble n3)
{
  const gdouble nsq[3] = {n1 * n1, n2 * n2, n3 * n3};
  gdouble ln_V;

  return _nc_hcosmo_bi_dist_t_lb_e (dist, z, nsq, &ln_V);
}

/**
 * nc_hcosmo_bi_dist_luminosity:
 * @dist: a #NcHCosmoBIDist
 * @z: redshift
 * @n1: direction component along axis 1
 * @n2: direction component along axis 2
 * @n3: direction component along axis 3
 *
 * Computes the luminosity distance $D_L$ of a source observed with redshift @z from the
 * unit direction $(n_1, n_2, n_3)$, see #NcHCosmoBIDist.
 *
 * Returns: $D_L$ in units of $c / H_0$
 */
gdouble
nc_hcosmo_bi_dist_luminosity (NcHCosmoBIDist *dist, const gdouble z, const gdouble n1, const gdouble n2, const gdouble n3)
{
  const guint n        = dist->n;
  const gdouble nsq[3] = {n1 * n1, n2 * n2, n3 * n3};
  const gdouble *u1    = &g_array_index (dist->u, gdouble, 0);
  const gdouble *u2    = &g_array_index (dist->u, gdouble, n + 1);
  const gdouble *u3    = &g_array_index (dist->u, gdouble, 2 * (n + 1));
  const gdouble *C     = &g_array_index (dist->C, gdouble, 0);
  gdouble *A           = &g_array_index (dist->A, gdouble, 0);
  gdouble *W           = &g_array_index (dist->W, gdouble, 0);
  gdouble ln_V;
  const gdouble t_e = _nc_hcosmo_bi_dist_t_lb_e (dist, z, nsq, &ln_V);
  gdouble K12 = 0.0, K13 = 0.0, K23 = 0.0;
  guint k, m;

  /*
   * Integrals of T_k(s) over s in [-1, s_e], with sin^2(phi) = (1 + s_e) / 2 = t_e / t_max,
   * using T_k(s_e) - (-1)^k = -2 (-1)^k sin^2(k phi), which keeps the relative precision
   * for t_e much smaller than t_max.
   */
  {
    const gdouble sp2 = t_e / dist->t_lb_max;
    const gdouble sp  = sqrt (sp2);
    const gdouble cp  = sqrt (1.0 - sp2);
    gdouble skm1      = sp;            /* sin((k - 1) phi) at k = 2 */
    gdouble sk        = 2.0 * sp * cp; /* sin(k phi) at k = 2 */

    A[0] = 2.0 * sp2;

    if (n >= 1)
      A[1] = -2.0 * sp2 * cp * cp;

    for (k = 2; k <= n; k++)
    {
      /* sin((k + 1) phi) = 2 cos(phi) sin(k phi) - sin((k - 1) phi) */
      const gdouble skp1 = 2.0 * cp * sk - skm1;
      const gdouble sgn  = (k % 2 == 0) ? 1.0 : -1.0;

      /* (-1)^(k+1) = (-1)^(k-1) = -sgn */
      A[k] = sgn * (skp1 * skp1 / (k + 1.0) - skm1 * skm1 / (k - 1.0));
      skm1 = sk;
      sk   = skp1;
    }
  }

  /* W = C^T A: integration weights of the node values */
  for (m = 0; m <= n; m++)
    W[m] = 0.0;

  for (k = 0; k <= n; k++)
  {
    const gdouble *Ck = &C[k * (n + 1)];
    const gdouble Ak  = A[k];

    for (m = 0; m <= n; m++)
      W[m] += Ck[m] * Ak;
  }

  for (m = 0; m <= n; m++)
  {
    const gdouble E2  = nsq[0] * u1[m] + nsq[1] * u2[m] + nsq[2] * u3[m];
    const gdouble wE3 = W[m] / (E2 * sqrt (E2));

    K12 += wE3 * u1[m] * u2[m];
    K13 += wE3 * u1[m] * u3[m];
    K23 += wE3 * u2[m] * u3[m];
  }

  {
    const gdouble h     = 0.5 * dist->t_lb_max;
    const gdouble Delta = h * h * (nsq[0] * K12 * K13 + nsq[1] * K12 * K23 + nsq[2] * K13 * K23);
    const gdouble DA2   = exp (ln_V) * (1.0 + z) * Delta;

    return (1.0 + z) * (1.0 + z) * sqrt (DA2);
  }
}

/**
 * nc_hcosmo_bi_dist_dmodulus_sky:
 * @dist: a #NcHCosmoBIDist
 * @z: redshift
 * @px: x component in the equatorial frame
 * @py: y component in the equatorial frame
 * @pz: z component in the equatorial frame
 *
 * Computes $5\log_{10}(D_L) + 25$ with $D_L$ from nc_hcosmo_bi_dist_luminosity() in units
 * of $c / H_0$, for the equatorial unit direction $(p_x, p_y, p_z)$. The distance modulus
 * is this value plus $5\log_{10}$ nc_hcosmo_bi_RH_Mpc().
 *
 * Returns: $5\log_{10}(D_L) + 25$
 */
gdouble
nc_hcosmo_bi_dist_dmodulus_sky (NcHCosmoBIDist *dist, const gdouble z, const gdouble px, const gdouble py, const gdouble pz)
{
  gdouble n1, n2, n3;

  nc_hcosmo_bi_dist_sky_to_axes (dist, px, py, pz, &n1, &n2, &n3);

  return 5.0 * log10 (nc_hcosmo_bi_dist_luminosity (dist, z, n1, n2, n3)) + 25.0;
}

/*
 * Solves E^2(t_e) = (1 + z)^2 by Newton iterations, safeguarded by bisection, on the
 * Chebyshev series of E^2 - 1 = sum_i n_i^2 u_i - 1 in delta = 1 + s = t_lb / h, h = t_max / 2;
 * the series is the combination of the stored coefficients of u_i, with c_0 lowered by one
 * so that the target is z (2 + z). Writes ln V(t_e) from its series to @ln_V and returns
 * t_e = h delta.
 */
static gdouble
_nc_hcosmo_bi_dist_t_lb_e (NcHCosmoBIDist *dist, const gdouble z, const gdouble nsq[3], gdouble *ln_V)
{
  const guint n        = dist->n;
  const gdouble h      = 0.5 * dist->t_lb_max;
  const gdouble *cu    = &g_array_index (dist->cu, gdouble, 0);
  const gdouble *cV    = &cu[3 * (n + 1)];
  const gdouble target = z * (2.0 + z);
  const guint max_iter = 200;
  gdouble *e           = &g_array_index (dist->e, gdouble, 0);
  gdouble lo           = 0.0;
  gdouble hi           = 2.0;
  gdouble delta;
  guint iter, k;

  if (z > dist->zmax)
    g_error ("nc_hcosmo_bi_dist: redshift %.15g above zmax = %.15g.", z, dist->zmax);

  for (k = 0; k <= n; k++)
    e[k] = nsq[0] * cu[k] + nsq[1] * cu[(n + 1) + k] + nsq[2] * cu[2 * (n + 1) + k];

  e[0] -= 1.0;

  /* E^2 - 1 ~ 2 t sum_i n_i^2 E_i(0) near t = 0 */
  delta = target / (2.0 * h * (nsq[0] * dist->E0[0] + nsq[1] * dist->E0[1] + nsq[2] * dist->E0[2]));

  if (!((delta > lo) && (delta < hi)))
    delta = 0.5 * (lo + hi);

  for (iter = 0; iter < max_iter; iter++)
  {
    const gdouble sv = delta - 1.0;
    gdouble Tkm1     = 1.0, Tk = sv;   /* T_0, T_1 */
    gdouble dTkm1    = 0.0, dTk = 1.0; /* T_0', T_1' */
    gdouble f  = e[0] + e[1] * sv;
    gdouble df = e[1];
    gdouble delta_new;

    for (k = 2; k <= n; k++)
    {
      const gdouble Tkp1  = 2.0 * sv * Tk - Tkm1;
      const gdouble dTkp1 = 2.0 * Tk + 2.0 * sv * dTk - dTkm1;

      f    += e[k] * Tkp1;
      df   += e[k] * dTkp1;
      Tkm1  = Tk;
      Tk    = Tkp1;
      dTkm1 = dTk;
      dTk   = dTkp1;
    }

    f -= target;

    if (f > 0.0)
      hi = delta;
    else
      lo = delta;

    delta_new = delta - f / df;

    if (!((delta_new > lo) && (delta_new < hi)))
      delta_new = 0.5 * (lo + hi);

    if (fabs (delta_new - delta) <= GSL_DBL_EPSILON * delta)
    {
      delta = delta_new;
      break;
    }

    delta = delta_new;
  }

  if (iter == max_iter)
    g_error ("nc_hcosmo_bi_dist: emission time for z = %.15g did not converge.", z);

  {
    const gdouble sv = delta - 1.0;
    gdouble Tkm1     = 1.0, Tk = sv;
    gdouble lV = cV[0] + cV[1] * sv;

    for (k = 2; k <= n; k++)
    {
      const gdouble Tkp1 = 2.0 * sv * Tk - Tkm1;

      lV  += cV[k] * Tkp1;
      Tkm1 = Tk;
      Tk   = Tkp1;
    }

    ln_V[0] = lV;
  }

  return h * delta;
}

/*
 * Finds the lookback time t_max at which max_i ln a_i = -ln(1 + zmax), so that every
 * direction has E >= 1 + zmax there, and returns the upper end of the final bisection
 * bracket. While every axis expands max_i ln a_i decreases, so the bracket is found by
 * probing t_lb = age (1 - 2^-k), k = 1, 2, ...; an axis contracting before the target is
 * reached (large shear) is an error.
 */
static gdouble
_nc_hcosmo_bi_dist_find_t_lb_max (NcHCosmoBI *cosmo, const gdouble zmax)
{
  const gdouble target = -log1p (zmax);
  const gdouble age    = nc_hcosmo_bi_age (cosmo);
  const guint max_iter = 200;
  gdouble lo           = 0.0;
  gdouble hi           = 0.5 * age;
  gdouble ln_a1, ln_a2, ln_a3, E1, E2, E3;
  guint iter;

  for (iter = 0; iter < max_iter; iter++)
  {
    nc_hcosmo_bi_E (cosmo, hi, &E1, &E2, &E3);

    if (GSL_MIN (E1, GSL_MIN (E2, E3)) <= 0.0)
      g_error ("nc_hcosmo_bi_dist_prepare: an axis contracts at t_lb = %.15g, before "
               "zmax = %.15g is reached; lower NcHCosmoBIDist:zmax.", hi, zmax);

    nc_hcosmo_bi_ln_a (cosmo, hi, &ln_a1, &ln_a2, &ln_a3);

    if (GSL_MAX (ln_a1, GSL_MAX (ln_a2, ln_a3)) <= target)
      break;

    lo = hi;
    hi = 0.5 * (hi + age);
  }

  if (iter == max_iter)
    g_error ("nc_hcosmo_bi_dist_prepare: zmax = %.15g is not reached before the "
             "singularity; lower NcHCosmoBIDist:zmax.", zmax);

  for (iter = 0; iter < max_iter; iter++)
  {
    const gdouble mid = 0.5 * (lo + hi);

    nc_hcosmo_bi_ln_a (cosmo, mid, &ln_a1, &ln_a2, &ln_a3);

    if (GSL_MAX (ln_a1, GSL_MAX (ln_a2, ln_a3)) > target)
      lo = mid;
    else
      hi = mid;

    if (hi - lo <= GSL_DBL_EPSILON * hi)
      break;
  }

  return hi;
}

/*
 * Integrands u_i u_j / E^3 of K_12, K_13 and K_23 along the test directions, then u_1, u_2,
 * u_3 and ln V, for the choice of the number of nodes in nc_hcosmo_bi_dist_prepare().
 */
static void
_nc_hcosmo_bi_dist_test_integrands (gpointer user_data, gdouble t_lb, NcmVector *fval)
{
  NcHCosmoBIDistTest *test = user_data;
  gdouble ln_a1, ln_a2, ln_a3;
  gdouble u1, u2, u3;
  guint d;

  nc_hcosmo_bi_ln_a (test->cosmo, t_lb, &ln_a1, &ln_a2, &ln_a3);
  u1 = exp (-2.0 * ln_a1);
  u2 = exp (-2.0 * ln_a2);
  u3 = exp (-2.0 * ln_a3);

  for (d = 0; d < NC_HCOSMO_BI_DIST_NTEST; d++)
  {
    const gdouble E2  = test->n2[d][0] * u1 + test->n2[d][1] * u2 + test->n2[d][2] * u3;
    const gdouble iE3 = 1.0 / (E2 * sqrt (E2));

    ncm_vector_set (fval, 3 * d + 0, u1 * u2 * iE3);
    ncm_vector_set (fval, 3 * d + 1, u1 * u3 * iE3);
    ncm_vector_set (fval, 3 * d + 2, u2 * u3 * iE3);
  }

  ncm_vector_set (fval, 3 * NC_HCOSMO_BI_DIST_NTEST + 0, u1);
  ncm_vector_set (fval, 3 * NC_HCOSMO_BI_DIST_NTEST + 1, u2);
  ncm_vector_set (fval, 3 * NC_HCOSMO_BI_DIST_NTEST + 2, u3);
  ncm_vector_set (fval, 3 * NC_HCOSMO_BI_DIST_NTEST + 3, ln_a1 + ln_a2 + ln_a3);
}

