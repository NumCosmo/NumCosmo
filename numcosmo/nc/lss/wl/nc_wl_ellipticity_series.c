/***************************************************************************
 *            nc_wl_ellipticity_series.c
 *
 *  Wed Jul 15 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_wl_ellipticity_series.c
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcWLEllipticitySeriesTrace:
 * NcWLEllipticitySeriesTraceDet:
 *
 * Truncated $g$-Taylor series of nc_wl_ellipticity.h's own closed-form
 * reduced-shear kernels: $\chi_I(\chi_L,g)$ (via apply_shear_inv) and the
 * *linear* (not log) Jacobian $\mathrm{Jac}(\chi_L,g)=|\det J|$ (via
 * det_jac), one object per ellipticity convention, mirroring that file's
 * own `_trace`/`_trace_det` split.
 *
 * Both conventions have closed forms in $g$ (no finite differences
 * anywhere) -- see docs/theory/wl_shape_marginalization_series.qmd and
 * `dev-notes/wl_shape_series_marginalization_derivation.py` sections 10-11
 * (verify_chi_closed_form_pieces, verify_eps_closed_form_pieces) for the
 * derivation and symbolic/numeric verification these two classes mirror:
 *
 * - TRACE_DET (epsilon): $\chi_I=(\chi_L-g)/(1-g\chi_L)$, holomorphic, so
 *   $\mathrm{Jac}=|f'(\chi_L,-g)|^2$; both closed-form power series in $g$
 *   with no recursion needed at all.
 * - TRACE (chi): $\chi_I(\chi_L,g)=f_\chi(\chi_L,-g)$ is a ratio of two
 *   QUADRATICS in $g$ ($\chi_L,\bar\chi_L$ fixed), so its Taylor series
 *   follows a 2-term power-series-division recursion; the chi map is NOT
 *   holomorphic in $\chi_L$ (involves $\bar\chi_L$ too), so the Jacobian
 *   needs the full real $2\times2$ determinant, itself a closed-form
 *   reciprocal-series-cubed construction. Both closed forms match the
 *   already-shipped
 *   nc_wl_ellipticity_apply_shear_inv_trace_c/nc_wl_ellipticity_det_jac_trace_c
 *   to machine precision (see dev-notes'
 *   verify_inverse_map_and_jacobian_match_production and this module's own
 *   test suite, tests/c/nc/lss/wl/test_nc_wl_ellipticity_series.c, which
 *   checks eval()+Horner-evaluate-at-small-g directly against those
 *   already-shipped closed forms).
 *
 * Every intermediate #NcmLaurentSeriesTPS/#NcmLaurentSeries this needs is
 * owned directly and allocated once at construction (order is
 * CONSTRUCT_ONLY): eval() only ever refills already-existing storage, never
 * allocates. TRACE_DET's chi and TRACE_DET's own $\chi_L$/$\bar\chi_L$/d1
 * quantities (shared between its chi and jac computation, since both run
 * sequentially against the same $\rho$ within one eval() call) are the only
 * scratch either class needs beyond the named #NcmLaurentSeriesTPS series
 * themselves.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/wl/nc_wl_ellipticity_series.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

/* ===========================================================================
 * NcWLEllipticitySeriesTraceDet (epsilon convention)
 * ===========================================================================
 */

typedef struct _NcWLEllipticitySeriesTraceDetPrivate
{
  guint order;
  NcmLaurentSeriesTPS *chi;
  NcmLaurentSeriesTPS *jac;
  NcmLaurentSeriesTPS *abs_sq;
  NcmLaurentSeriesTPS *chi_conj;
  NcmLaurentSeriesTPS *binom;
  NcmLaurentSeriesTPS *deriv;
  NcmLaurentSeriesTPS *deriv_conj;
} NcWLEllipticitySeriesTraceDetPrivate;

struct _NcWLEllipticitySeriesTraceDet
{
  GObject parent_instance;
};

enum
{
  PROP_TRACE_DET_0,
  PROP_TRACE_DET_ORDER,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcWLEllipticitySeriesTraceDet, nc_wl_ellipticity_series_trace_det, G_TYPE_OBJECT)

static void
nc_wl_ellipticity_series_trace_det_init (NcWLEllipticitySeriesTraceDet *ser)
{
  NcWLEllipticitySeriesTraceDetPrivate * const self = nc_wl_ellipticity_series_trace_det_get_instance_private (ser);

  self->order      = 0;
  self->chi        = NULL;
  self->jac        = NULL;
  self->abs_sq     = NULL;
  self->chi_conj   = NULL;
  self->binom      = NULL;
  self->deriv      = NULL;
  self->deriv_conj = NULL;
}

static void
_nc_wl_ellipticity_series_trace_det_set_order (NcWLEllipticitySeriesTraceDet *ser, guint order)
{
  NcWLEllipticitySeriesTraceDetPrivate * const self = nc_wl_ellipticity_series_trace_det_get_instance_private (ser);

  self->order      = order;
  self->chi        = ncm_laurent_series_tps_new (order);
  self->jac        = ncm_laurent_series_tps_new (order);
  self->abs_sq     = ncm_laurent_series_tps_new (order);
  self->chi_conj   = ncm_laurent_series_tps_new (order);
  self->binom      = ncm_laurent_series_tps_new (order + 2);
  self->deriv      = ncm_laurent_series_tps_new (order);
  self->deriv_conj = ncm_laurent_series_tps_new (order);
}

static void
_nc_wl_ellipticity_series_trace_det_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcWLEllipticitySeriesTraceDet *ser = NC_WL_ELLIPTICITY_SERIES_TRACE_DET (object);

  switch (prop_id)
  {
    case PROP_TRACE_DET_ORDER:
      _nc_wl_ellipticity_series_trace_det_set_order (ser, g_value_get_uint (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_wl_ellipticity_series_trace_det_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcWLEllipticitySeriesTraceDet *ser = NC_WL_ELLIPTICITY_SERIES_TRACE_DET (object);

  switch (prop_id)
  {
    case PROP_TRACE_DET_ORDER:
      g_value_set_uint (value, nc_wl_ellipticity_series_trace_det_get_order (ser));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_wl_ellipticity_series_trace_det_finalize (GObject *object)
{
  NcWLEllipticitySeriesTraceDetPrivate * const self = nc_wl_ellipticity_series_trace_det_get_instance_private (NC_WL_ELLIPTICITY_SERIES_TRACE_DET (object));

  ncm_laurent_series_tps_clear (&self->chi);
  ncm_laurent_series_tps_clear (&self->jac);
  ncm_laurent_series_tps_clear (&self->abs_sq);
  ncm_laurent_series_tps_clear (&self->chi_conj);
  ncm_laurent_series_tps_clear (&self->binom);
  ncm_laurent_series_tps_clear (&self->deriv);
  ncm_laurent_series_tps_clear (&self->deriv_conj);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_wl_ellipticity_series_trace_det_parent_class)->finalize (object);
}

static void
nc_wl_ellipticity_series_trace_det_class_init (NcWLEllipticitySeriesTraceDetClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_wl_ellipticity_series_trace_det_set_property;
  object_class->get_property = &_nc_wl_ellipticity_series_trace_det_get_property;
  object_class->finalize     = &_nc_wl_ellipticity_series_trace_det_finalize;

  /**
   * NcWLEllipticitySeriesTraceDet:order:
   *
   * Truncation order $N$ of the $g$-power series. Immutable for the
   * object's whole life -- every #NcmLaurentSeriesTPS this needs is sized
   * and allocated once, here.
   */
  g_object_class_install_property (object_class,
                                   PROP_TRACE_DET_ORDER,
                                   g_param_spec_uint ("order",
                                                      NULL,
                                                      "Truncation order",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_wl_ellipticity_series_trace_det_new:
 * @order: truncation order $N$ of the $g$-power series
 *
 * Returns: (transfer full): a new #NcWLEllipticitySeriesTraceDet
 */
NcWLEllipticitySeriesTraceDet *
nc_wl_ellipticity_series_trace_det_new (guint order)
{
  return g_object_new (NC_TYPE_WL_ELLIPTICITY_SERIES_TRACE_DET, "order", order, NULL);
}

/**
 * nc_wl_ellipticity_series_trace_det_ref:
 * @ser: a #NcWLEllipticitySeriesTraceDet
 *
 * Returns: (transfer full): @ser
 */
NcWLEllipticitySeriesTraceDet *
nc_wl_ellipticity_series_trace_det_ref (NcWLEllipticitySeriesTraceDet *ser)
{
  return g_object_ref (ser);
}

void
nc_wl_ellipticity_series_trace_det_free (NcWLEllipticitySeriesTraceDet *ser)
{
  g_object_unref (ser);
}

void
nc_wl_ellipticity_series_trace_det_clear (NcWLEllipticitySeriesTraceDet **ser)
{
  g_clear_object (ser);
}

/**
 * nc_wl_ellipticity_series_trace_det_get_order:
 * @ser: a #NcWLEllipticitySeriesTraceDet
 *
 * Returns: @ser's truncation order $N$
 */
guint
nc_wl_ellipticity_series_trace_det_get_order (NcWLEllipticitySeriesTraceDet *ser)
{
  NcWLEllipticitySeriesTraceDetPrivate * const self = nc_wl_ellipticity_series_trace_det_get_instance_private (ser);

  return self->order;
}

/**
 * nc_wl_ellipticity_series_trace_det_eval:
 * @ser: a #NcWLEllipticitySeriesTraceDet
 * @rho: $\rho=|\chi_L|$
 *
 * Refills @ser's chi and jac series (see nc_wl_ellipticity_series_trace_det_get_chi()/
 * _get_jac()) at $\chi_L=\rho\,w$.
 *
 * eps: $\chi_I(\chi_L,g)=(\chi_L-g)/(1-g\chi_L)$, closed form, no
 * recursion: $c_0=\chi_L=\{1{:}\rho\}$; $c_k=\chi_L^{k-1}(\chi_L^2-1)=
 * \rho^{k+1}w^{k+1}-\rho^{k-1}w^{k-1}$ for $k\ge1$ -- a 2-term Laurent
 * series at every order, written directly into @ser's own stable slots (no
 * scratch needed at all).
 *
 * eps: $f_\mathrm{eps}'(\chi_L,-g)=(1-g^2)/(1-g\chi_L)^2$ Taylor-in-$g$
 * coefficients -- closed form: $1/(1-g\chi_L)^2=\sum_n(n+1)\chi_L^n g^n$
 * (literal binomial series, single-term at harmonic $n$ each), times
 * $(1-g^2)$. $\mathrm{Jac}=|\mathrm{that}|^2$ is this series convolved with
 * its own conjugate. Also fills nc_wl_ellipticity_series_trace_det_get_abs_sq()
 * ($|\chi_I(\chi_L,g)|^2=\chi_I\overline{\chi_I}$, population-independent,
 * exactly the $\rho^2(g)$ input nc_galaxy_shape_pop_eval_p_rho2_g_series()
 * needs), since it's a pure function of chi already computed here.
 */
void
nc_wl_ellipticity_series_trace_det_eval (NcWLEllipticitySeriesTraceDet *ser, gdouble rho)
{
  NcWLEllipticitySeriesTraceDetPrivate * const self = nc_wl_ellipticity_series_trace_det_get_instance_private (ser);
  const guint order                                 = self->order;
  guint k, n;

  ncm_laurent_series_set_single_into (ncm_laurent_series_tps_get (self->chi, 0), 1, rho);

  for (k = 1; k <= order; k++)
  {
    gdouble rho_km1      = pow (rho, (gdouble) (k - 1));
    gdouble rho_kp1      = rho_km1 * rho * rho;
    NcmLaurentSeries *ck = ncm_laurent_series_tps_get (self->chi, k);

    ncm_laurent_series_reset (ck, (gint) k - 1, (gint) k + 1);
    ncm_laurent_series_set_c (ck, (gint) k - 1, -rho_km1);
    ncm_laurent_series_set_c (ck, (gint) k + 1, rho_kp1);
  }

  ncm_laurent_series_tps_conj (self->chi_conj, self->chi);
  ncm_laurent_series_tps_conv (self->abs_sq, self->chi, self->chi_conj);

  for (n = 0; n <= order + 2; n++)
    ncm_laurent_series_set_single_into (ncm_laurent_series_tps_get (self->binom, n), (gint) n, (n + 1.0) * pow (rho, (gdouble) n));

  for (n = 0; n <= order; n++)
  {
    if (n >= 2)
      ncm_laurent_series_add_into (ncm_laurent_series_tps_get (self->deriv, n), ncm_laurent_series_tps_get (self->binom, n), ncm_laurent_series_tps_get (self->binom, n - 2), -1.0);
    else
      ncm_laurent_series_scale_c_into (ncm_laurent_series_tps_get (self->deriv, n), ncm_laurent_series_tps_get (self->binom, n), 1.0);
  }

  ncm_laurent_series_tps_conj (self->deriv_conj, self->deriv);
  ncm_laurent_series_tps_conv (self->jac, self->deriv, self->deriv_conj);
}

/**
 * nc_wl_ellipticity_series_trace_det_get_chi:
 * @ser: a #NcWLEllipticitySeriesTraceDet
 *
 * Returns: (transfer none): $\chi_I(\chi_L,g)$'s $g$-Taylor coefficients as
 * of the most recent nc_wl_ellipticity_series_trace_det_eval() call
 */
NcmLaurentSeriesTPS *
nc_wl_ellipticity_series_trace_det_get_chi (NcWLEllipticitySeriesTraceDet *ser)
{
  NcWLEllipticitySeriesTraceDetPrivate * const self = nc_wl_ellipticity_series_trace_det_get_instance_private (ser);

  return self->chi;
}

/**
 * nc_wl_ellipticity_series_trace_det_get_jac:
 * @ser: a #NcWLEllipticitySeriesTraceDet
 *
 * Returns: (transfer none): $\mathrm{Jac}(\chi_L,g)$'s $g$-Taylor
 * coefficients as of the most recent
 * nc_wl_ellipticity_series_trace_det_eval() call
 */
NcmLaurentSeriesTPS *
nc_wl_ellipticity_series_trace_det_get_jac (NcWLEllipticitySeriesTraceDet *ser)
{
  NcWLEllipticitySeriesTraceDetPrivate * const self = nc_wl_ellipticity_series_trace_det_get_instance_private (ser);

  return self->jac;
}

/**
 * nc_wl_ellipticity_series_trace_det_get_abs_sq:
 * @ser: a #NcWLEllipticitySeriesTraceDet
 *
 * Returns: (transfer none): $|\chi_I(\chi_L,g)|^2$'s $g$-Taylor
 * coefficients (population-independent shear-map output) as of the most
 * recent nc_wl_ellipticity_series_trace_det_eval() call
 */
NcmLaurentSeriesTPS *
nc_wl_ellipticity_series_trace_det_get_abs_sq (NcWLEllipticitySeriesTraceDet *ser)
{
  NcWLEllipticitySeriesTraceDetPrivate * const self = nc_wl_ellipticity_series_trace_det_get_instance_private (ser);

  return self->abs_sq;
}

/* ===========================================================================
 * NcWLEllipticitySeriesTrace (chi/distortion convention)
 * ===========================================================================
 */

typedef struct _NcWLEllipticitySeriesTracePrivate
{
  guint order;
  NcmLaurentSeriesTPS *chi;
  NcmLaurentSeriesTPS *jac;
  NcmLaurentSeriesTPS *abs_sq;
  NcmLaurentSeriesTPS *chi_conj;
  NcmLaurentSeriesTPS *r;
  NcmLaurentSeriesTPS *r2;
  NcmLaurentSeriesTPS *r3;
  NcmLaurentSeriesTPS *num;   /* (1-g^2)^3: order-only-dependent, filled once at construction */
  NcmLaurentSeries *chibar_L; /* recomputed fresh every eval() (depends on rho) */
  NcmLaurentSeries *two_s;
  NcmLaurentSeries *d1;
  NcmLaurentSeries *scratch_term; /* transient scratch shared by the chi and jac
                                   * recursions below (never concurrent within
                                   * one eval() call) */
  NcmLaurentSeries *scratch_base;
  NcmLaurentSeries *scratch_acc[2];
} NcWLEllipticitySeriesTracePrivate;

struct _NcWLEllipticitySeriesTrace
{
  GObject parent_instance;
};

enum
{
  PROP_TRACE_0,
  PROP_TRACE_ORDER,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcWLEllipticitySeriesTrace, nc_wl_ellipticity_series_trace, G_TYPE_OBJECT)

static void
nc_wl_ellipticity_series_trace_init (NcWLEllipticitySeriesTrace *ser)
{
  NcWLEllipticitySeriesTracePrivate * const self = nc_wl_ellipticity_series_trace_get_instance_private (ser);

  self->order          = 0;
  self->chi            = NULL;
  self->jac            = NULL;
  self->abs_sq         = NULL;
  self->chi_conj       = NULL;
  self->r              = NULL;
  self->r2             = NULL;
  self->r3             = NULL;
  self->num            = NULL;
  self->chibar_L       = NULL;
  self->two_s          = NULL;
  self->d1             = NULL;
  self->scratch_term   = NULL;
  self->scratch_base   = NULL;
  self->scratch_acc[0] = NULL;
  self->scratch_acc[1] = NULL;
}

static void
_nc_wl_ellipticity_series_trace_set_order (NcWLEllipticitySeriesTrace *ser, guint order)
{
  NcWLEllipticitySeriesTracePrivate * const self = nc_wl_ellipticity_series_trace_get_instance_private (ser);
  guint k;

  self->order          = order;
  self->chi            = ncm_laurent_series_tps_new (order);
  self->jac            = ncm_laurent_series_tps_new (order);
  self->abs_sq         = ncm_laurent_series_tps_new (order);
  self->chi_conj       = ncm_laurent_series_tps_new (order);
  self->r              = ncm_laurent_series_tps_new (order);
  self->r2             = ncm_laurent_series_tps_new (order);
  self->r3             = ncm_laurent_series_tps_new (order);
  self->num            = ncm_laurent_series_tps_new (order);
  self->chibar_L       = ncm_laurent_series_new (0, 0);
  self->two_s          = ncm_laurent_series_new (0, 0);
  self->d1             = ncm_laurent_series_new (0, 0);
  self->scratch_term   = ncm_laurent_series_new (0, 0);
  self->scratch_base   = ncm_laurent_series_new (0, 0);
  self->scratch_acc[0] = ncm_laurent_series_new (0, 0);
  self->scratch_acc[1] = ncm_laurent_series_new (0, 0);

  /* (1-g^2)^3 = 1 - 3g^2 + 3g^4 - g^6: order-only-dependent (no rho), so
   * filled once here rather than on every eval(). */
  for (k = 0; k <= order; k++)
  {
    gdouble coeff;

    if (k == 0)
      coeff = 1.0;
    else if (k == 2)
      coeff = -3.0;
    else if (k == 4)
      coeff = 3.0;
    else if (k == 6)
      coeff = -1.0;
    else
      coeff = 0.0;

    ncm_laurent_series_set_single_into (ncm_laurent_series_tps_get (self->num, k), 0, coeff);
  }
}

static void
_nc_wl_ellipticity_series_trace_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcWLEllipticitySeriesTrace *ser = NC_WL_ELLIPTICITY_SERIES_TRACE (object);

  switch (prop_id)
  {
    case PROP_TRACE_ORDER:
      _nc_wl_ellipticity_series_trace_set_order (ser, g_value_get_uint (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_wl_ellipticity_series_trace_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcWLEllipticitySeriesTrace *ser = NC_WL_ELLIPTICITY_SERIES_TRACE (object);

  switch (prop_id)
  {
    case PROP_TRACE_ORDER:
      g_value_set_uint (value, nc_wl_ellipticity_series_trace_get_order (ser));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_wl_ellipticity_series_trace_finalize (GObject *object)
{
  NcWLEllipticitySeriesTracePrivate * const self = nc_wl_ellipticity_series_trace_get_instance_private (NC_WL_ELLIPTICITY_SERIES_TRACE (object));

  ncm_laurent_series_tps_clear (&self->chi);
  ncm_laurent_series_tps_clear (&self->jac);
  ncm_laurent_series_tps_clear (&self->abs_sq);
  ncm_laurent_series_tps_clear (&self->chi_conj);
  ncm_laurent_series_tps_clear (&self->r);
  ncm_laurent_series_tps_clear (&self->r2);
  ncm_laurent_series_tps_clear (&self->r3);
  ncm_laurent_series_tps_clear (&self->num);
  ncm_laurent_series_clear (&self->chibar_L);
  ncm_laurent_series_clear (&self->two_s);
  ncm_laurent_series_clear (&self->d1);
  ncm_laurent_series_clear (&self->scratch_term);
  ncm_laurent_series_clear (&self->scratch_base);
  ncm_laurent_series_clear (&self->scratch_acc[0]);
  ncm_laurent_series_clear (&self->scratch_acc[1]);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_wl_ellipticity_series_trace_parent_class)->finalize (object);
}

static void
nc_wl_ellipticity_series_trace_class_init (NcWLEllipticitySeriesTraceClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_wl_ellipticity_series_trace_set_property;
  object_class->get_property = &_nc_wl_ellipticity_series_trace_get_property;
  object_class->finalize     = &_nc_wl_ellipticity_series_trace_finalize;

  /**
   * NcWLEllipticitySeriesTrace:order:
   *
   * Truncation order $N$ of the $g$-power series. Immutable for the
   * object's whole life -- every #NcmLaurentSeriesTPS/#NcmLaurentSeries
   * this needs is sized and allocated once, here.
   */
  g_object_class_install_property (object_class,
                                   PROP_TRACE_ORDER,
                                   g_param_spec_uint ("order",
                                                      NULL,
                                                      "Truncation order",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_wl_ellipticity_series_trace_new:
 * @order: truncation order $N$ of the $g$-power series
 *
 * Returns: (transfer full): a new #NcWLEllipticitySeriesTrace
 */
NcWLEllipticitySeriesTrace *
nc_wl_ellipticity_series_trace_new (guint order)
{
  return g_object_new (NC_TYPE_WL_ELLIPTICITY_SERIES_TRACE, "order", order, NULL);
}

/**
 * nc_wl_ellipticity_series_trace_ref:
 * @ser: a #NcWLEllipticitySeriesTrace
 *
 * Returns: (transfer full): @ser
 */
NcWLEllipticitySeriesTrace *
nc_wl_ellipticity_series_trace_ref (NcWLEllipticitySeriesTrace *ser)
{
  return g_object_ref (ser);
}

void
nc_wl_ellipticity_series_trace_free (NcWLEllipticitySeriesTrace *ser)
{
  g_object_unref (ser);
}

void
nc_wl_ellipticity_series_trace_clear (NcWLEllipticitySeriesTrace **ser)
{
  g_clear_object (ser);
}

/**
 * nc_wl_ellipticity_series_trace_get_order:
 * @ser: a #NcWLEllipticitySeriesTrace
 *
 * Returns: @ser's truncation order $N$
 */
guint
nc_wl_ellipticity_series_trace_get_order (NcWLEllipticitySeriesTrace *ser)
{
  NcWLEllipticitySeriesTracePrivate * const self = nc_wl_ellipticity_series_trace_get_instance_private (ser);

  return self->order;
}

/* chi: chi_I(chi_L,g)=f_chi(chi_L,-g) is Num(g)/Den(g), ratio of two
 * quadratics in g (chi_L,chibar_L fixed): Num=chi_L-2g+chibar_L*g^2,
 * Den=1-(chi_L+chibar_L)*g+g^2 -- 2-term power-series-division recursion,
 * c_0=chi_L=n_0; c_k=n_k-d1*c_{k-1}-d2*c_{k-2} (k>=1), n_1=-2, n_2=chibar_L,
 * n_k=0 otherwise, d1=-(chi_L+chibar_L), d2=1 (so the d2 term is just
 * c_{k-2} itself, no convolution needed). c_0=chi_L is set by the caller
 * (nc_wl_ellipticity_series_trace_eval(), which also fills @self->d1);
 * every other c_k is assembled in @self's shared scratch and copied into
 * @self->chi's slot k at the end (that slot has stable identity, so a
 * plain pointer handoff is not an option). */
static void
_nc_wl_ellipticity_series_trace_fill_chi (NcWLEllipticitySeriesTracePrivate * const self)
{
  const guint order = self->order;
  guint k;

  for (k = 1; k <= order; k++)
  {
    NcmLaurentSeries *base;
    NcmLaurentSeries *acc;
    NcmLaurentSeries *ck = ncm_laurent_series_tps_get (self->chi, k);

    if (k == 1)
    {
      base = self->scratch_base;
      ncm_laurent_series_set_single_into (base, 0, -2.0);
    }
    else if (k == 2)
    {
      base = self->chibar_L;
    }
    else
    {
      base = self->scratch_base;
      ncm_laurent_series_reset (base, 0, 0);
    }

    ncm_laurent_series_conv_into (self->scratch_term, self->d1, ncm_laurent_series_tps_get (self->chi, k - 1));

    acc = self->scratch_acc[0];
    ncm_laurent_series_add_into (acc, base, self->scratch_term, -1.0);

    if (k >= 2)
    {
      NcmLaurentSeries *acc2 = self->scratch_acc[1];

      ncm_laurent_series_add_into (acc2, acc, ncm_laurent_series_tps_get (self->chi, k - 2), -1.0);
      acc = acc2;
    }

    ncm_laurent_series_scale_c_into (ck, acc, 1.0);
  }
}

/* chi: Jac(chi_L,g) = (1-g^2)^3 / D(g)^3, D(g) the same denominator as
 * _fill_chi above (the chi map is not holomorphic in chi_L, so the full
 * real 2x2 Jacobian is needed -- derived by hand, confirmed to machine
 * precision against the already-shipped nc_wl_ellipticity_det_jac_trace_c,
 * see the class docs). r = 1/D(g), D(g)=1+d1*g+g^2 (d0=1, d2=1); r2=r*r,
 * r3=r2*r via ncm_laurent_series_tps_conv(); num=(1-g^2)^3 was already
 * filled once at construction (order-only-dependent). */
static void
_nc_wl_ellipticity_series_trace_fill_jac (NcWLEllipticitySeriesTracePrivate * const self)
{
  const guint order = self->order;
  guint k;

  ncm_laurent_series_set_single_into (ncm_laurent_series_tps_get (self->r, 0), 0, 1.0);

  for (k = 1; k <= order; k++)
  {
    NcmLaurentSeries *acc;

    ncm_laurent_series_conv_into (self->scratch_term, self->d1, ncm_laurent_series_tps_get (self->r, k - 1));

    acc = self->scratch_acc[0];
    ncm_laurent_series_scale_c_into (acc, self->scratch_term, -1.0);

    if (k >= 2)
    {
      NcmLaurentSeries *acc2 = self->scratch_acc[1];

      ncm_laurent_series_add_into (acc2, acc, ncm_laurent_series_tps_get (self->r, k - 2), -1.0);
      acc = acc2;
    }

    ncm_laurent_series_scale_c_into (ncm_laurent_series_tps_get (self->r, k), acc, 1.0);
  }

  ncm_laurent_series_tps_conv (self->r2, self->r, self->r);
  ncm_laurent_series_tps_conv (self->r3, self->r2, self->r);
  ncm_laurent_series_tps_conv (self->jac, self->num, self->r3);
}

/**
 * nc_wl_ellipticity_series_trace_eval:
 * @ser: a #NcWLEllipticitySeriesTrace
 * @rho: $\rho=|\chi_L|$
 *
 * Refills @ser's chi and jac series (see nc_wl_ellipticity_series_trace_get_chi()/
 * _get_jac()) at $\chi_L=\rho\,w$. $\chi_L,\bar\chi_L,d_1=-(\chi_L+\bar\chi_L)$
 * are computed once here and shared by both the chi and jac recursions
 * (identical quantities, same $\rho$, never needed concurrently). Also
 * fills nc_wl_ellipticity_series_trace_get_abs_sq() ($|\chi_I(\chi_L,g)|^2=
 * \chi_I\overline{\chi_I}$, population-independent, exactly the $\rho^2(g)$
 * input nc_galaxy_shape_pop_eval_p_rho2_g_series() needs).
 */
void
nc_wl_ellipticity_series_trace_eval (NcWLEllipticitySeriesTrace *ser, gdouble rho)
{
  NcWLEllipticitySeriesTracePrivate * const self = nc_wl_ellipticity_series_trace_get_instance_private (ser);
  NcmLaurentSeries *chi_L                        = ncm_laurent_series_tps_get (self->chi, 0);

  ncm_laurent_series_set_single_into (chi_L, 1, rho);
  ncm_laurent_series_set_single_into (self->chibar_L, -1, rho);
  ncm_laurent_series_add_into (self->two_s, chi_L, self->chibar_L, 1.0);
  ncm_laurent_series_scale_c_into (self->d1, self->two_s, -1.0);

  _nc_wl_ellipticity_series_trace_fill_chi (self);
  _nc_wl_ellipticity_series_trace_fill_jac (self);

  ncm_laurent_series_tps_conj (self->chi_conj, self->chi);
  ncm_laurent_series_tps_conv (self->abs_sq, self->chi, self->chi_conj);
}

/**
 * nc_wl_ellipticity_series_trace_get_chi:
 * @ser: a #NcWLEllipticitySeriesTrace
 *
 * Returns: (transfer none): $\chi_I(\chi_L,g)$'s $g$-Taylor coefficients as
 * of the most recent nc_wl_ellipticity_series_trace_eval() call
 */
NcmLaurentSeriesTPS *
nc_wl_ellipticity_series_trace_get_chi (NcWLEllipticitySeriesTrace *ser)
{
  NcWLEllipticitySeriesTracePrivate * const self = nc_wl_ellipticity_series_trace_get_instance_private (ser);

  return self->chi;
}

/**
 * nc_wl_ellipticity_series_trace_get_jac:
 * @ser: a #NcWLEllipticitySeriesTrace
 *
 * Returns: (transfer none): $\mathrm{Jac}(\chi_L,g)$'s $g$-Taylor
 * coefficients as of the most recent nc_wl_ellipticity_series_trace_eval()
 * call
 */
NcmLaurentSeriesTPS *
nc_wl_ellipticity_series_trace_get_jac (NcWLEllipticitySeriesTrace *ser)
{
  NcWLEllipticitySeriesTracePrivate * const self = nc_wl_ellipticity_series_trace_get_instance_private (ser);

  return self->jac;
}

/**
 * nc_wl_ellipticity_series_trace_get_abs_sq:
 * @ser: a #NcWLEllipticitySeriesTrace
 *
 * Returns: (transfer none): $|\chi_I(\chi_L,g)|^2$'s $g$-Taylor
 * coefficients (population-independent shear-map output) as of the most
 * recent nc_wl_ellipticity_series_trace_eval() call
 */
NcmLaurentSeriesTPS *
nc_wl_ellipticity_series_trace_get_abs_sq (NcWLEllipticitySeriesTrace *ser)
{
  NcWLEllipticitySeriesTracePrivate * const self = nc_wl_ellipticity_series_trace_get_instance_private (ser);

  return self->abs_sq;
}

