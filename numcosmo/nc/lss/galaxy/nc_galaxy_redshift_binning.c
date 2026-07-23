/***************************************************************************
 *            nc_galaxy_redshift_binning.c
 *
 *  Tue Jul 1 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_binning.c
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
 * NcGalaxyRedshiftBinning:
 *
 * Binning calculator: the true-redshift distribution of a photometric bin.
 *
 * A calculator (a plain #GObject, NOT an #NcmModel and NOT held in an #NcmMSet)
 * that produces the true-redshift distribution $\mathrm{d}n/\mathrm{d}z$ of the
 * galaxies selected into a photometric window $[z_{p,\min}, z_{p,\max}]$, from a
 * population model $P(z\mid I)$ (a #NcGalaxyRedshiftPop) and a population photo-z
 * observable $P(z_p\mid z)$ (a #NcGalaxyRedshiftObsSel):
 * $$ \frac{\mathrm{d}n}{\mathrm{d}z}(z) \propto P(z\mid I)\,
 *    \frac{W(z; z_{p,\min}, z_{p,\max})}{N(z)}, \qquad
 *    N(z) = \int_0^\infty P(z_p\mid z)\,\mathrm{d}z_p, $$
 * where $W$ is the observable's selection mass in the window and $N$ the physical
 * (photo-z $\ge 0$) normalization. Both come from the observable's generic
 * nc_galaxy_redshift_obs_sel_window_mass(), so the calculator is scheme-free (no
 * per-kernel subclass).
 *
 * The photometric window is NOT object state: a single calculator produces
 * $\mathrm{d}n/\mathrm{d}z$ for arbitrarily many bins. The window is an argument
 * to nc_galaxy_redshift_binning_compute_dndz() (and the on-nodes variant), which
 * are pure producers: they return a freshly-built, normalized #NcmSpline on each
 * call and cache nothing, so they do NOT require nc_galaxy_redshift_binning_prepare().
 *
 * nc_galaxy_redshift_binning_prepare() builds only the window-free marginal
 * photo-z density $P(z_p)$ used by nc_galaxy_redshift_binning_eval_pzp() and the
 * equal-area edges nc_galaxy_redshift_binning_compute_equal_area_photoz_bins();
 * those two read the cached marginal and take no models. Re-call prepare after
 * changing the model parameters.
 *
 * Following the #NcDistance convention, the calculator does NOT hold the models:
 * they are consumed only inside the methods that need them and are passed as
 * arguments. The cached marginal is built with an adaptive-knot scheme that
 * guarantees a relative interpolation-error tolerance
 * (#NcGalaxyRedshiftBinning:reltol). Photo-z systematics (shift/stretch) are
 * intentionally NOT here: they belong to a future n(z) model layered on top of
 * this producer.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_redshift_binning.h"
#include "nc/lss/galaxy/nc_galaxy_redshift_obs_sel_gauss.h"
#include "ncm/algebra/ncm_vector.h"
#include "ncm/spline/ncm_spline_cubic_notaknot.h"
#include "ncm/spline/ncm_spline_func.h"
#include "ncm/integration/ncm_integral1d_ptr.h"
#include "ncm/stats/ncm_stats_dist1d_spline.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <float.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

#define NC_GALAXY_REDSHIFT_BINNING_DEFAULT_RELTOL (1.0e-7)
#define NC_GALAXY_REDSHIFT_BINNING_DEFAULT_ZP_SUPPORT_MAX (20.0)

typedef struct _NcGalaxyRedshiftBinningPrivate
{
  gdouble reltol;
  gdouble zp_support_max;
  NcmSpline *pzp;
  NcmStatsDist1d *pzp_stats;
} NcGalaxyRedshiftBinningPrivate;

struct _NcGalaxyRedshiftBinning
{
  GObject parent_instance;
};

enum
{
  PROP_0,
  PROP_RELTOL,
  PROP_ZP_SUPPORT_MAX,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyRedshiftBinning, nc_galaxy_redshift_binning, G_TYPE_OBJECT);

static void
nc_galaxy_redshift_binning_init (NcGalaxyRedshiftBinning *gsdrb)
{
  NcGalaxyRedshiftBinningPrivate * const self = nc_galaxy_redshift_binning_get_instance_private (gsdrb);

  self->reltol         = NC_GALAXY_REDSHIFT_BINNING_DEFAULT_RELTOL;
  self->zp_support_max = NC_GALAXY_REDSHIFT_BINNING_DEFAULT_ZP_SUPPORT_MAX;
  self->pzp            = NULL;
  self->pzp_stats      = NULL;
}

static void
_nc_galaxy_redshift_binning_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyRedshiftBinning *gsdrb = NC_GALAXY_REDSHIFT_BINNING (object);

  g_return_if_fail (NC_IS_GALAXY_REDSHIFT_BINNING (gsdrb));

  switch (prop_id)
  {
    case PROP_RELTOL:
      nc_galaxy_redshift_binning_set_reltol (gsdrb, g_value_get_double (value));
      break;
    case PROP_ZP_SUPPORT_MAX:
      nc_galaxy_redshift_binning_set_zp_support_max (gsdrb, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_redshift_binning_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyRedshiftBinning *gsdrb              = NC_GALAXY_REDSHIFT_BINNING (object);
  NcGalaxyRedshiftBinningPrivate * const self = nc_galaxy_redshift_binning_get_instance_private (gsdrb);

  g_return_if_fail (NC_IS_GALAXY_REDSHIFT_BINNING (gsdrb));

  switch (prop_id)
  {
    case PROP_RELTOL:
      g_value_set_double (value, self->reltol);
      break;
    case PROP_ZP_SUPPORT_MAX:
      g_value_set_double (value, self->zp_support_max);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_redshift_binning_dispose (GObject *object)
{
  NcGalaxyRedshiftBinning *gsdrb              = NC_GALAXY_REDSHIFT_BINNING (object);
  NcGalaxyRedshiftBinningPrivate * const self = nc_galaxy_redshift_binning_get_instance_private (gsdrb);

  ncm_spline_clear (&self->pzp);
  ncm_stats_dist1d_clear (&self->pzp_stats);

  G_OBJECT_CLASS (nc_galaxy_redshift_binning_parent_class)->dispose (object);
}

static void
_nc_galaxy_redshift_binning_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_redshift_binning_parent_class)->finalize (object);
}

static void
nc_galaxy_redshift_binning_class_init (NcGalaxyRedshiftBinningClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_redshift_binning_set_property;
  object_class->get_property = &_nc_galaxy_redshift_binning_get_property;
  object_class->dispose      = &_nc_galaxy_redshift_binning_dispose;
  object_class->finalize     = &_nc_galaxy_redshift_binning_finalize;

  /**
   * NcGalaxyRedshiftBinning:reltol:
   *
   * Relative interpolation-error tolerance for the cached adaptive splines.
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative interpolation tolerance",
                                                        DBL_EPSILON, 1.0e-2, NC_GALAXY_REDSHIFT_BINNING_DEFAULT_RELTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyRedshiftBinning:zp-support-max:
   *
   * Maximum photometric redshift for the marginal $P(z_p)$ support, used by
   * nc_galaxy_redshift_binning_eval_pzp() and the equal-area bin edges.
   */
  g_object_class_install_property (object_class,
                                   PROP_ZP_SUPPORT_MAX,
                                   g_param_spec_double ("zp-support-max",
                                                        NULL,
                                                        "Maximum photometric redshift for P(zp) support",
                                                        0.1, 100.0, NC_GALAXY_REDSHIFT_BINNING_DEFAULT_ZP_SUPPORT_MAX,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/* dn/dz machinery ------------------------------------------------------- */

/* Transient closure over the models and window for the adaptive-knot dn/dz
 * build; lives only for the duration of a compute call (nothing is held). */
typedef struct _NcGalaxyRedshiftBinningDndzCtx
{
  NcGalaxyRedshiftPop *population;
  NcGalaxyRedshiftObsSel *observable_population;
  gdouble zp_min;
  gdouble zp_max;
} NcGalaxyRedshiftBinningDndzCtx;

static gdouble
_nc_galaxy_redshift_binning_dndz_raw (NcGalaxyRedshiftBinningDndzCtx *ctx, const gdouble z)
{
  const gdouble Pz = nc_galaxy_redshift_pop_eval (ctx->population, z);
  const gdouble W  = nc_galaxy_redshift_obs_sel_window_mass (ctx->observable_population, z, ctx->zp_min, ctx->zp_max);
  /* Physical photo-z domain zp >= 0; N(z) = P(zp in [0, inf) | z). */
  const gdouble N = nc_galaxy_redshift_obs_sel_window_mass (ctx->observable_population, z, 0.0, INFINITY);

  return Pz * W / N;
}

static gdouble
_nc_galaxy_redshift_binning_dndz_raw_gsl (gdouble z, gpointer user_data)
{
  return _nc_galaxy_redshift_binning_dndz_raw ((NcGalaxyRedshiftBinningDndzCtx *) user_data, z);
}

/* Locate where the integrand drops below @threshold of its peak, so the
 * adaptive spline does not waste knots on the empty tails of the support. */
static void
_nc_galaxy_redshift_binning_effective_support (gsl_function *F, gdouble z_min, gdouble z_max, gdouble threshold, gdouble *z_min_eff, gdouble *z_max_eff)
{
  const guint n_search = 500;
  gdouble f_max        = -G_MAXDOUBLE;
  guint i_max          = 0;
  gdouble z_array[n_search];
  gdouble f_array[n_search];
  guint i;

  for (i = 0; i < n_search; i++)
  {
    z_array[i] = z_min + (z_max - z_min) * i / (n_search - 1.0);
    f_array[i] = GSL_FN_EVAL (F, z_array[i]);

    if (f_array[i] > f_max)
    {
      f_max = f_array[i];
      i_max = i;
    }
  }

  *z_min_eff = z_min;

  for (i = i_max; i > 0; i--)
  {
    if (f_array[i] < threshold * f_max)
    {
      *z_min_eff = z_array[i];
      break;
    }
  }

  *z_max_eff = z_max;

  for (i = i_max; i < n_search; i++)
  {
    if (f_array[i] < threshold * f_max)
    {
      *z_max_eff = z_array[i];
      break;
    }
  }

  if (*z_min_eff >= *z_max_eff)
  {
    *z_min_eff = z_min;
    *z_max_eff = z_max;
  }
}

/* Build a fresh, unit-normalized adaptive dn/dz spline for the window
 * [zp_min, zp_max]. No object state is touched; caller owns the result. */
static NcmSpline *
_nc_galaxy_redshift_binning_build_dndz (NcGalaxyRedshiftBinning *gsdrb, NcGalaxyRedshiftPop *population, NcGalaxyRedshiftObsSel *observable_population, const gdouble zp_min, const gdouble zp_max)
{
  NcGalaxyRedshiftBinningPrivate * const self = nc_galaxy_redshift_binning_get_instance_private (gsdrb);
  NcGalaxyRedshiftBinningDndzCtx ctx          = { population, observable_population, zp_min, zp_max };
  gdouble z_min, z_max, z_min_eff, z_max_eff, norm;
  NcmSpline *dndz;
  gsl_function F;

  g_assert_nonnull (population);
  g_assert_nonnull (observable_population);
  g_assert_cmpfloat (zp_min, <, zp_max);

  nc_galaxy_redshift_pop_get_lim (population, &z_min, &z_max);

  F.function = &_nc_galaxy_redshift_binning_dndz_raw_gsl;
  F.params   = &ctx;

  _nc_galaxy_redshift_binning_effective_support (&F, z_min, z_max, 1.0e-20, &z_min_eff, &z_max_eff);

  dndz = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  ncm_spline_set_func (dndz, NCM_SPLINE_FUNCTION_SPLINE, &F, z_min_eff, z_max_eff, 0, self->reltol);

  norm = ncm_spline_eval_integ (dndz, z_min_eff, z_max_eff);
  g_assert_cmpfloat (norm, >, 0.0);

  ncm_vector_scale (ncm_spline_peek_yv (dndz), 1.0 / norm);
  ncm_spline_prepare (dndz);

  return dndz;
}

/* Marginal photo-z density machinery ------------------------------------ */

typedef struct _NcGalaxyRedshiftBinningPzpIntegData
{
  NcGalaxyRedshiftPop *population;
  NcGalaxyRedshiftObsSel *observable_population;
  NcmIntegral1dPtr *integrator;
  gdouble zp;
  gdouble z_min;
  gdouble z_max;
} NcGalaxyRedshiftBinningPzpIntegData;

/* Integrand of P(zp) = int P(z) f(zp|z) / N(z) dz over the true redshift z,
 * with N(z) = window_mass(z, 0, inf) the physical zp >= 0 normalization. */
static gdouble
_nc_galaxy_redshift_binning_pzp_integrand (gpointer user_data, const gdouble z, const gdouble w)
{
  NcGalaxyRedshiftBinningPzpIntegData *data = (NcGalaxyRedshiftBinningPzpIntegData *) user_data;
  const gdouble Pz                          = nc_galaxy_redshift_pop_eval (data->population, z);
  const gdouble f                           = nc_galaxy_redshift_obs_sel_eval (data->observable_population, z, data->zp);
  const gdouble N                           = nc_galaxy_redshift_obs_sel_window_mass (data->observable_population, z, 0.0, INFINITY);

  return Pz * f / N;
}

static gdouble
_nc_galaxy_redshift_binning_m2lnpzp_gsl (gdouble zp, gpointer user_data)
{
  NcGalaxyRedshiftBinningPzpIntegData *data = (NcGalaxyRedshiftBinningPzpIntegData *) user_data;
  gdouble error                             = 0.0;
  gdouble result;

  data->zp = zp;
  ncm_integral1d_ptr_set_userdata (data->integrator, data);
  result = ncm_integral1d_eval (NCM_INTEGRAL1D (data->integrator), data->z_min, data->z_max, &error);

  g_assert_cmpfloat (result, >, 0.0);

  return -2.0 * log (fabs (result));
}

static void
_nc_galaxy_redshift_binning_assert_pzp_prepared (NcGalaxyRedshiftBinning *gsdrb)
{
  NcGalaxyRedshiftBinningPrivate * const self = nc_galaxy_redshift_binning_get_instance_private (gsdrb);

  if (self->pzp_stats == NULL)
    g_error ("nc_galaxy_redshift_binning: the marginal P(zp) has not been prepared; call nc_galaxy_redshift_binning_prepare () first.");
}

/* Public API ------------------------------------------------------------ */

/**
 * nc_galaxy_redshift_binning_new:
 *
 * Creates a new #NcGalaxyRedshiftBinning. The photometric window is not held by
 * the calculator: it is supplied to each dn/dz method. The population and
 * observable models are likewise not held; they are passed to every method that
 * needs them.
 *
 * Returns: (transfer full): a new #NcGalaxyRedshiftBinning.
 */
NcGalaxyRedshiftBinning *
nc_galaxy_redshift_binning_new (void)
{
  NcGalaxyRedshiftBinning *gsdrb = g_object_new (NC_TYPE_GALAXY_REDSHIFT_BINNING,
                                                 NULL);

  return gsdrb;
}

/**
 * nc_galaxy_redshift_binning_ref:
 * @gsdrb: a #NcGalaxyRedshiftBinning
 *
 * Increases the reference count of @gsdrb by one.
 *
 * Returns: (transfer full): @gsdrb.
 */
NcGalaxyRedshiftBinning *
nc_galaxy_redshift_binning_ref (NcGalaxyRedshiftBinning *gsdrb)
{
  return g_object_ref (gsdrb);
}

/**
 * nc_galaxy_redshift_binning_free:
 * @gsdrb: a #NcGalaxyRedshiftBinning
 *
 * Decreases the reference count of @gsdrb by one.
 *
 */
void
nc_galaxy_redshift_binning_free (NcGalaxyRedshiftBinning *gsdrb)
{
  g_object_unref (gsdrb);
}

/**
 * nc_galaxy_redshift_binning_clear:
 * @gsdrb: a #NcGalaxyRedshiftBinning
 *
 * Decreases the reference count of @gsdrb by one, and sets the pointer *@gsdrb
 * to NULL.
 *
 */
void
nc_galaxy_redshift_binning_clear (NcGalaxyRedshiftBinning **gsdrb)
{
  g_clear_object (gsdrb);
}

/**
 * nc_galaxy_redshift_binning_set_reltol:
 * @gsdrb: a #NcGalaxyRedshiftBinning
 * @reltol: the relative interpolation tolerance
 *
 * Sets the relative interpolation tolerance and invalidates the cached marginal
 * $P(z_p)$.
 *
 */
void
nc_galaxy_redshift_binning_set_reltol (NcGalaxyRedshiftBinning *gsdrb, const gdouble reltol)
{
  NcGalaxyRedshiftBinningPrivate * const self = nc_galaxy_redshift_binning_get_instance_private (gsdrb);

  self->reltol = reltol;

  ncm_spline_clear (&self->pzp);
  ncm_stats_dist1d_clear (&self->pzp_stats);
}

/**
 * nc_galaxy_redshift_binning_get_reltol:
 * @gsdrb: a #NcGalaxyRedshiftBinning
 *
 * Returns: the relative interpolation tolerance.
 */
gdouble
nc_galaxy_redshift_binning_get_reltol (NcGalaxyRedshiftBinning *gsdrb)
{
  NcGalaxyRedshiftBinningPrivate * const self = nc_galaxy_redshift_binning_get_instance_private (gsdrb);

  return self->reltol;
}

/**
 * nc_galaxy_redshift_binning_set_zp_support_max:
 * @gsdrb: a #NcGalaxyRedshiftBinning
 * @zp_support_max: the maximum photometric redshift for the P(zp) support
 *
 * Sets the maximum photometric redshift over which the marginal $P(z_p)$ is
 * tabulated and invalidates the cached marginal distribution. Does not affect the
 * bin $\mathrm{d}n/\mathrm{d}z$, which is independent of the P(zp) support.
 *
 */
void
nc_galaxy_redshift_binning_set_zp_support_max (NcGalaxyRedshiftBinning *gsdrb, const gdouble zp_support_max)
{
  NcGalaxyRedshiftBinningPrivate * const self = nc_galaxy_redshift_binning_get_instance_private (gsdrb);

  self->zp_support_max = zp_support_max;

  ncm_spline_clear (&self->pzp);
  ncm_stats_dist1d_clear (&self->pzp_stats);
}

/**
 * nc_galaxy_redshift_binning_get_zp_support_max:
 * @gsdrb: a #NcGalaxyRedshiftBinning
 *
 * Returns: the maximum photometric redshift for the P(zp) support.
 */
gdouble
nc_galaxy_redshift_binning_get_zp_support_max (NcGalaxyRedshiftBinning *gsdrb)
{
  NcGalaxyRedshiftBinningPrivate * const self = nc_galaxy_redshift_binning_get_instance_private (gsdrb);

  return self->zp_support_max;
}

/**
 * nc_galaxy_redshift_binning_compute_dndz:
 * @gsdrb: a #NcGalaxyRedshiftBinning
 * @population: a #NcGalaxyRedshiftPop
 * @observable_population: a #NcGalaxyRedshiftObsSel
 * @zp_min: the minimum photometric redshift of the bin
 * @zp_max: the maximum photometric redshift of the bin
 *
 * Builds the normalized bin $\mathrm{d}n/\mathrm{d}z$ for the photometric window
 * [@zp_min, @zp_max] as a fresh adaptive-knot #NcmSpline meeting
 * #NcGalaxyRedshiftBinning:reltol over its effective support, normalized to unit
 * integral. Pure producer: consumes the models on the call, holds nothing, and
 * does NOT require nc_galaxy_redshift_binning_prepare().
 *
 * Returns: (transfer full): a #NcmSpline of $\mathrm{d}n/\mathrm{d}z$.
 */
NcmSpline *
nc_galaxy_redshift_binning_compute_dndz (NcGalaxyRedshiftBinning *gsdrb, NcGalaxyRedshiftPop *population, NcGalaxyRedshiftObsSel *observable_population, const gdouble zp_min, const gdouble zp_max)
{
  return _nc_galaxy_redshift_binning_build_dndz (gsdrb, population, observable_population, zp_min, zp_max);
}

/**
 * nc_galaxy_redshift_binning_compute_dndz_on_nodes:
 * @gsdrb: a #NcGalaxyRedshiftBinning
 * @population: a #NcGalaxyRedshiftPop
 * @observable_population: a #NcGalaxyRedshiftObsSel
 * @zp_min: the minimum photometric redshift of the bin
 * @zp_max: the maximum photometric redshift of the bin
 * @z_nodes: the true-redshift nodes to tabulate on
 *
 * As nc_galaxy_redshift_binning_compute_dndz(), but returns the normalized bin
 * $\mathrm{d}n/\mathrm{d}z$ tabulated on the given @z_nodes (0 outside the
 * effective support). Same pure-producer semantics: holds nothing and does NOT
 * require nc_galaxy_redshift_binning_prepare().
 *
 * Returns: (transfer full): a #NcmSpline of $\mathrm{d}n/\mathrm{d}z$ on @z_nodes.
 */
NcmSpline *
nc_galaxy_redshift_binning_compute_dndz_on_nodes (NcGalaxyRedshiftBinning *gsdrb, NcGalaxyRedshiftPop *population, NcGalaxyRedshiftObsSel *observable_population, const gdouble zp_min, const gdouble zp_max, NcmVector *z_nodes)
{
  NcmSpline *dndz      = _nc_galaxy_redshift_binning_build_dndz (gsdrb, population, observable_population, zp_min, zp_max);
  NcmVector *dndz_xv   = ncm_spline_peek_xv (dndz);
  const guint dndz_len = ncm_vector_len (dndz_xv);
  const gdouble z_lo   = ncm_vector_get (dndz_xv, 0);
  const gdouble z_hi   = ncm_vector_get (dndz_xv, dndz_len - 1);
  const guint len      = ncm_vector_len (z_nodes);
  NcmVector *xv        = ncm_vector_dup (z_nodes);
  NcmVector *yv        = ncm_vector_new (len);
  NcmSpline *dndz_nodes;
  guint i;

  for (i = 0; i < len; i++)
  {
    const gdouble z = ncm_vector_get (z_nodes, i);
    const gdouble y = ((z < z_lo) || (z > z_hi)) ? 0.0 : ncm_spline_eval (dndz, z);

    ncm_vector_set (yv, i, y);
  }

  dndz_nodes = NCM_SPLINE (ncm_spline_cubic_notaknot_new_full (xv, yv, TRUE));

  ncm_vector_free (xv);
  ncm_vector_free (yv);
  ncm_spline_free (dndz);

  return dndz_nodes;
}

/**
 * nc_galaxy_redshift_binning_prepare:
 * @gsdrb: a #NcGalaxyRedshiftBinning
 * @population: a #NcGalaxyRedshiftPop
 * @observable_population: a #NcGalaxyRedshiftObsSel
 *
 * Builds the window-free marginal photo-z density $P(z_p)$ over
 * $[0, \mathtt{zp\_support\_max}]$ used by nc_galaxy_redshift_binning_eval_pzp()
 * and nc_galaxy_redshift_binning_compute_equal_area_photoz_bins(). The models are
 * consumed here and NOT retained; re-call this after changing their parameters.
 *
 * The dn/dz producers (nc_galaxy_redshift_binning_compute_dndz() and the on-nodes
 * variant) are stateless and do NOT depend on this preparation.
 *
 */
void
nc_galaxy_redshift_binning_prepare (NcGalaxyRedshiftBinning *gsdrb, NcGalaxyRedshiftPop *population, NcGalaxyRedshiftObsSel *observable_population)
{
  NcGalaxyRedshiftBinningPrivate * const self = nc_galaxy_redshift_binning_get_instance_private (gsdrb);
  NcGalaxyRedshiftBinningPzpIntegData integ_data;
  gdouble z_min, z_max;
  gsl_function F;

  g_assert_nonnull (population);
  g_assert_nonnull (observable_population);

  nc_galaxy_redshift_pop_get_lim (population, &z_min, &z_max);

  integ_data.population            = population;
  integ_data.observable_population = observable_population;
  integ_data.integrator            = ncm_integral1d_ptr_new (&_nc_galaxy_redshift_binning_pzp_integrand, NULL);
  integ_data.zp                    = 0.0;
  integ_data.z_min                 = z_min;
  integ_data.z_max                 = z_max;

  F.function = &_nc_galaxy_redshift_binning_m2lnpzp_gsl;
  F.params   = &integ_data;

  ncm_spline_clear (&self->pzp);
  self->pzp = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  ncm_spline_set_func (self->pzp, NCM_SPLINE_FUNCTION_SPLINE, &F, 0.0, self->zp_support_max, 0, self->reltol);

  ncm_integral1d_ptr_free (integ_data.integrator);

  ncm_stats_dist1d_clear (&self->pzp_stats);
  self->pzp_stats = NCM_STATS_DIST1D (ncm_stats_dist1d_spline_new (self->pzp));
  g_object_set (self->pzp_stats, "abstol", 1.0e-50, NULL);
  ncm_stats_dist1d_prepare (self->pzp_stats);
}

/**
 * nc_galaxy_redshift_binning_eval_pzp:
 * @gsdrb: a #NcGalaxyRedshiftBinning
 * @zp: the photometric redshift
 *
 * Evaluates the marginal photometric-redshift density
 * $$ P(z_p) = \int P(z\mid I)\,\frac{f(z_p\mid z)}{N(z)}\,\mathrm{d}z, \qquad
 *    N(z) = \int_0^\infty f(z_p'\mid z)\,\mathrm{d}z_p', $$
 * the distribution of photo-z over the whole population (independent of any bin
 * window). Requires a prior nc_galaxy_redshift_binning_prepare().
 *
 * Returns: the marginal density $P(z_p)$ at @zp.
 */
gdouble
nc_galaxy_redshift_binning_eval_pzp (NcGalaxyRedshiftBinning *gsdrb, const gdouble zp)
{
  NcGalaxyRedshiftBinningPrivate * const self = nc_galaxy_redshift_binning_get_instance_private (gsdrb);

  _nc_galaxy_redshift_binning_assert_pzp_prepared (gsdrb);

  return ncm_stats_dist1d_eval_p (self->pzp_stats, zp);
}

/**
 * nc_galaxy_redshift_binning_compute_equal_area_photoz_bins:
 * @gsdrb: a #NcGalaxyRedshiftBinning
 * @n_bins: the number of bins
 * @zp_max: the maximum photometric redshift to bin up to
 *
 * Computes @n_bins + 1 photometric-redshift edges that split the marginal
 * $P(z_p)$ (see nc_galaxy_redshift_binning_eval_pzp()) into @n_bins slices of
 * equal integrated probability, by inverting its CDF. Suitable for equal-area
 * source binning. Requires a prior nc_galaxy_redshift_binning_prepare(); @zp_max
 * must not exceed #NcGalaxyRedshiftBinning:zp-support-max.
 *
 * Returns: (transfer full): a #NcmVector with @n_bins + 1 photo-z edges.
 */
NcmVector *
nc_galaxy_redshift_binning_compute_equal_area_photoz_bins (NcGalaxyRedshiftBinning *gsdrb, const guint n_bins, const gdouble zp_max)
{
  NcGalaxyRedshiftBinningPrivate * const self = nc_galaxy_redshift_binning_get_instance_private (gsdrb);
  NcmVector *bin_edges;
  gdouble total_area;
  guint i;

  if (zp_max > self->zp_support_max)
    g_error ("nc_galaxy_redshift_binning_compute_equal_area_photoz_bins: requested zp_max=%.2f exceeds zp_support_max=%.2f.", zp_max, self->zp_support_max);

  _nc_galaxy_redshift_binning_assert_pzp_prepared (gsdrb);

  bin_edges  = ncm_vector_new (n_bins + 1);
  total_area = ncm_stats_dist1d_eval_pdf (self->pzp_stats, zp_max);

  for (i = 0; i <= n_bins; i++)
  {
    const gdouble target_area = (i * total_area) / n_bins;
    const gdouble zp_edge     = ncm_stats_dist1d_eval_inv_pdf (self->pzp_stats, target_area);

    ncm_vector_set (bin_edges, i, zp_edge);
  }

  return bin_edges;
}

/**
 * nc_galaxy_redshift_binning_lsst_srd_edges:
 * @type: a #NcGalaxyRedshiftPopLSSTSRDType
 * @population: (out) (transfer full): the LSST SRD population model for @type
 * @observable_population: (out) (transfer full): a Gaussian photo-z observable
 *   with the per-type scatter $\sigma_0$
 *
 * Computes the LSST SRD tomographic photo-z bin edges for @type and hands back
 * the matching population + Gaussian observable models (both configured for
 * @type) via the out-parameters. The edges follow the LSST DESC SRD recipe:
 *
 * - Lens (Y1, Y10): linearly spaced edges over $[0.2, 1.2]$
 *   (Y1: 5 bins, Y10: 10 bins), $\sigma_0 = 0.03$.
 * - Source (Y1, Y10): equal-area edges up to $z_p = 3.5$
 *   (5 bins), $\sigma_0 = 0.05$; computed by inverting the marginal $P(z_p)$.
 *
 * The returned edges feed nc_galaxy_redshift_binning_compute_dndz() (per
 * consecutive pair) with the handed-back models. This is a factory: it holds
 * nothing.
 *
 * Returns: (transfer full): a #NcmVector with n_bins + 1 photo-z edges.
 */
NcmVector *
nc_galaxy_redshift_binning_lsst_srd_edges (NcGalaxyRedshiftPopLSSTSRDType type, NcGalaxyRedshiftPop **population, NcGalaxyRedshiftObsSel **observable_population)
{
  NcGalaxyRedshiftPopLSSTSRD *pop        = nc_galaxy_redshift_pop_lsst_srd_new_from_type (type);
  NcGalaxyRedshiftObsSelGauss *obs_gauss = nc_galaxy_redshift_obs_sel_gauss_new ();
  NcGalaxyRedshiftObsSel *obs_sel        = NC_GALAXY_REDSHIFT_OBS_SEL (obs_gauss);
  NcmVector *bin_edges;
  gboolean is_lens;
  gdouble sigma0;
  guint n_bins;

  switch (type)
  {
    case NC_GALAXY_REDSHIFT_POP_LSST_SRD_Y1_SOURCE:
      n_bins  = 5;
      sigma0  = 0.05;
      is_lens = FALSE;
      break;
    case NC_GALAXY_REDSHIFT_POP_LSST_SRD_Y1_LENS:
      n_bins  = 5;
      sigma0  = 0.03;
      is_lens = TRUE;
      break;
    case NC_GALAXY_REDSHIFT_POP_LSST_SRD_Y10_SOURCE:
      n_bins  = 5;
      sigma0  = 0.05;
      is_lens = FALSE;
      break;
    case NC_GALAXY_REDSHIFT_POP_LSST_SRD_Y10_LENS:
      n_bins  = 10;
      sigma0  = 0.03;
      is_lens = TRUE;
      break;
    default:                                                                        /* LCOV_EXCL_LINE */
      g_error ("nc_galaxy_redshift_binning_lsst_srd_edges: invalid type %d", type); /* LCOV_EXCL_LINE */

      return NULL; /* LCOV_EXCL_LINE */
  }

  ncm_model_param_set (NCM_MODEL (obs_gauss), NC_GALAXY_REDSHIFT_OBS_SEL_GAUSS_SIGMA0, sigma0);

  if (is_lens)
  {
    const gdouble zp_min = 0.2;
    const gdouble zp_max = 1.2;
    guint i;

    bin_edges = ncm_vector_new (n_bins + 1);

    for (i = 0; i <= n_bins; i++)
      ncm_vector_set (bin_edges, i, zp_min + i * (zp_max - zp_min) / n_bins);
  }
  else
  {
    NcGalaxyRedshiftBinning *gsdrb = nc_galaxy_redshift_binning_new ();

    nc_galaxy_redshift_binning_prepare (gsdrb, NC_GALAXY_REDSHIFT_POP (pop), obs_sel);
    bin_edges = nc_galaxy_redshift_binning_compute_equal_area_photoz_bins (gsdrb, n_bins, 3.5);

    nc_galaxy_redshift_binning_free (gsdrb);
  }

  *population            = NC_GALAXY_REDSHIFT_POP (pop);
  *observable_population = obs_sel;

  return bin_edges;
}

