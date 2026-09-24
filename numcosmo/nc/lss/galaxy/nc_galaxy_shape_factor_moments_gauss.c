/***************************************************************************
 *            nc_galaxy_shape_factor_moments_gauss.c
 *
 *  Mon Sep 22 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_moments_gauss.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * NcGalaxyShapeFactorMomentsGauss:
 *
 * A Gaussian matched to the exact moments of the lensed, noisy
 * intrinsic-ellipticity marginal.
 *
 * In the frame rotated by $-\arg g$, with $(x,y)$ the observed ellipticity
 * in that frame,
 * $$
 *   \ln P(\chi_\mathrm{obs}\mid g) = -\frac{(x-\mu)^2}{2C_t} - \frac{y^2}{2C_x}
 *     - \frac{1}{2}\ln(C_t C_x) - \ln 2\pi,
 * $$
 * with $\mu = \mathrm{E}[x]$, $C_t = \mathrm{Var}(x)$ and
 * $C_x = \mathrm{E}[y^2]$ of the exact marginal. These are the same exact
 * moments #NcGalaxyShapeFactorMomentsTilt matches its tilt to, computed in
 * closed form for both ellipticity conventions; this class stops at the
 * Gaussian those moments define, instead of solving for the tilt.
 *
 * ## Tables
 *
 * The noise variance $\sigma_\nu^2$ enters the second moments only
 * additively, so the tables hold the moments of the lensed source alone,
 * $(\mu, \mathrm{Var}_\mathrm{src}(x), \mathrm{E}_\mathrm{src}[y^2])$, and
 * each galaxy's $\sigma_\nu^2$ is added at evaluation. A table therefore
 * depends only on the population, and every galaxy of a population with a
 * single dispersion shares one. The variance is tabulated rather than
 * $\mathrm{E}[x^2]$: near $\hat g = 1$ the difference
 * $\mathrm{E}[x^2] - \mu^2$ cancels down to a small number, and
 * interpolating the raw moment would carry its error into $C_t$.
 *
 * The moments are tabulated on the same dyadic Chebyshev mesh in
 * $\hat g = \min(|g|, 1/|g|)$ as #NcGalaxyShapeFactorMomentsTilt, all
 * panels, and each panel is resolved to the absolute tolerance
 * #NcGalaxyShapeFactorMomentsGauss:moment-tol. A node is one radial
 * quadrature, so the tolerance can be set far below anything that matters
 * rather than tuned: the error it leaves in $\ln P$ scales as
 * $1/\sigma_\nu^2$, and a loose tolerance would tie the accuracy to the
 * smallest $\sigma_\nu$ in the catalog.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_factor_moments_gauss.h"
#include "nc/lss/galaxy/nc_galaxy_shape_factor_moments_private.h"
#include "ncm/core/ncm_memory_pool.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

/* Three interleaved series: mean, source variance along the shear, source
 * second moment across it. */
#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS_N_COMP 3

/* Panels of a galaxy's table data_prefetch() fetches: the tables span the
 * full mesh, but weak-lensing shears keep nearly every evaluation in the
 * first two (ghat < 3/4). A later panel is still read correctly, only not
 * prefetched. */
#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS_PREFETCH_PANELS 2
#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS_PREFETCH_CAP 4096

struct _NcGalaxyShapeFactorMomentsGauss
{
  NcGalaxyShapeFactor parent_instance;
};

typedef struct _NcGalaxyShapeFactorMomentsGaussPrivate
{
  NcGalaxyWLObsEllipConv ellip_conv;
  gdouble moment_tol;
  guint max_degree;
  guint64 pop_hash;

  /* (e_rms, n_panels) -> table, valid for pop_hash == tab_cache_hash. The
   * key's sigma_nu slot is always zero: the tables do not depend on it. */
  GHashTable *tab_cache;
  guint64 tab_cache_hash;
  NcmMemoryPool *spectral_pool;
  GMutex cache_lock;
  gint table_build_count;
} NcGalaxyShapeFactorMomentsGaussPrivate;

typedef struct _NcGalaxyShapeFactorMomentsGaussLData
{
  NcGalaxyShapeFactorMomentsTable *table;
  guint64 pop_hash_seen;
  gboolean valid;

  /* Bytes of @table data_prefetch() fetches, set with @table so the
   * prefetch never reads the table header. */
  gsize table_span;
} NcGalaxyShapeFactorMomentsGaussLData;

enum
{
  PROP_0,
  PROP_MOMENT_TOL,
  PROP_MAX_DEGREE,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyShapeFactorMomentsGauss, nc_galaxy_shape_factor_moments_gauss, NC_TYPE_GALAXY_SHAPE_FACTOR)

/*
 * ---- The table build ----
 */
typedef struct _NcGalaxyShapeFactorMomentsGaussBuild
{
  NcGalaxyWLObsEllipConv ellip_conv;
  NcGalaxyShapePop *pop;
  NcGalaxyShapePopData *pop_data;
  GArray *r_arr;
  GArray *p_arr;
} NcGalaxyShapeFactorMomentsGaussBuild;

static gboolean
_moments_gauss_node_cb (gpointer user_data, const gdouble ghat, const NcGalaxyShapeFactorMomentsMemo *memo, gdouble *vals)
{
  NcGalaxyShapeFactorMomentsGaussBuild *b = (NcGalaxyShapeFactorMomentsGaussBuild *) user_data;
  gdouble t[3];

  /* Source moments: no noise. */
  _nc_galaxy_shape_factor_moments_exact_t (b->ellip_conv, b->pop, b->pop_data, ghat, 0.0, b->r_arr, &b->p_arr, t);

  vals[0] = t[0];
  vals[1] = t[1] - t[0] * t[0];
  vals[2] = t[2];

  return TRUE;
}

static void
_moments_gauss_panel_ctx_cb (gpointer user_data, const gdouble lo, const gdouble hi, gdouble *ctx)
{
}

/* Each series is resolved to the same absolute tolerance, so the size of a
 * coefficient vector is its largest component. */
static gdouble
_moments_gauss_probe_cb (gpointer user_data, const gdouble *v, const gdouble *ctx)
{
  return MAX (fabs (v[0]), MAX (fabs (v[1]), fabs (v[2])));
}

static NcGalaxyShapeFactorMomentsTable *
_moments_gauss_build_table (NcGalaxyShapeFactorMomentsGaussPrivate * const self, NcmSpectral *spectral,
                            NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data)
{
  NcGalaxyShapeFactorMomentsGaussBuild b;
  NcGalaxyShapeFactorMomentsTable *table;

  b.ellip_conv = self->ellip_conv;
  b.pop        = pop;
  b.pop_data   = pop_data;
  b.r_arr      = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), NC_GALAXY_SHAPE_FACTOR_MOMENTS_MOMENT_NNODES);
  b.p_arr      = NULL;

  g_array_set_size (b.r_arr, NC_GALAXY_SHAPE_FACTOR_MOMENTS_MOMENT_NNODES);

  table = _nc_galaxy_shape_factor_moments_build (spectral, NC_GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS_N_COMP,
                                                 NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_PANELS, 1.0,
                                                 self->max_degree, self->moment_tol, self->moment_tol,
                                                 &_moments_gauss_node_cb, &_moments_gauss_panel_ctx_cb,
                                                 &_moments_gauss_probe_cb, &b);

  g_array_unref (b.r_arr);

  if (b.p_arr != NULL)
    g_array_unref (b.p_arr);

  return table;
}

/* The pool keeps the returned pointer in a slot of its own and hands out
 * the slot: ncm_memory_pool_get() returns an NcmSpectral **. */
static gpointer
_moments_gauss_spectral_alloc (gpointer userdata)
{
  return ncm_spectral_new_with_max_order (6);
}

static void
_moments_gauss_spectral_free (gpointer p)
{
  ncm_spectral_free (NCM_SPECTRAL (p));
}

static void
nc_galaxy_shape_factor_moments_gauss_init (NcGalaxyShapeFactorMomentsGauss *gsfmg)
{
  NcGalaxyShapeFactorMomentsGaussPrivate * const self = nc_galaxy_shape_factor_moments_gauss_get_instance_private (gsfmg);

  self->ellip_conv = NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE;
  self->moment_tol = 1.0e-10;
  self->max_degree = 32;
  self->pop_hash   = 0;
  self->tab_cache  = g_hash_table_new_full (&_nc_galaxy_shape_factor_moments_key_hash,
                                            &_nc_galaxy_shape_factor_moments_key_equal,
                                            &g_free, (GDestroyNotify) & _nc_galaxy_shape_factor_moments_table_unref);
  self->tab_cache_hash = 0;
  self->spectral_pool  = ncm_memory_pool_new (&_moments_gauss_spectral_alloc, NULL, &_moments_gauss_spectral_free);
  g_mutex_init (&self->cache_lock);
  self->table_build_count = 0;
}

/* The convention is a construct-only property of the parent, so it is read
 * once here: the moments branch on it. */
static void
_nc_galaxy_shape_factor_moments_gauss_constructed (GObject *object)
{
  /* Chain up: start */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_moments_gauss_parent_class)->constructed (object);
  {
    NcGalaxyShapeFactorMomentsGaussPrivate * const self =
      nc_galaxy_shape_factor_moments_gauss_get_instance_private (NC_GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS (object));

    self->ellip_conv = nc_galaxy_shape_factor_get_ellip_conv (NC_GALAXY_SHAPE_FACTOR (object));
  }
}

static void
_nc_galaxy_shape_factor_moments_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactorMomentsGaussPrivate * const self =
    nc_galaxy_shape_factor_moments_gauss_get_instance_private (NC_GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS (object));

  switch (prop_id)
  {
    case PROP_MOMENT_TOL:
      self->moment_tol = g_value_get_double (value);
      break;
    case PROP_MAX_DEGREE:
      self->max_degree = g_value_get_uint (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_factor_moments_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactorMomentsGaussPrivate * const self =
    nc_galaxy_shape_factor_moments_gauss_get_instance_private (NC_GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS (object));

  switch (prop_id)
  {
    case PROP_MOMENT_TOL:
      g_value_set_double (value, self->moment_tol);
      break;
    case PROP_MAX_DEGREE:
      g_value_set_uint (value, self->max_degree);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_factor_moments_gauss_ldata_destroy (gpointer p)
{
  NcGalaxyShapeFactorMomentsGaussLData *ldata = (NcGalaxyShapeFactorMomentsGaussLData *) p;

  g_clear_pointer (&ldata->table, _nc_galaxy_shape_factor_moments_table_unref);
  g_free (ldata);
}

static void
_nc_galaxy_shape_factor_moments_gauss_ldata_noop (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
}

/* A per-galaxy population reads its dispersion from the catalog row, so a
 * new row invalidates the table without any model pkey moving. */
static void
_nc_galaxy_shape_factor_moments_gauss_ldata_read_row (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxyShapeFactorMomentsGaussLData *ldata = (NcGalaxyShapeFactorMomentsGaussLData *) data->ldata;

  ldata->valid = FALSE;
}

static void
_nc_galaxy_shape_factor_moments_gauss_ldata_required_columns (NcGalaxyShapeFactorData *data, GList **columns)
{
}

static void
_nc_galaxy_shape_factor_moments_gauss_data_init (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data)
{
  NcGalaxyShapeFactorMomentsGaussLData *ldata = g_new0 (NcGalaxyShapeFactorMomentsGaussLData, 1);

  data->ldata                  = ldata;
  data->ldata_destroy          = &_nc_galaxy_shape_factor_moments_gauss_ldata_destroy;
  data->ldata_read_row         = &_nc_galaxy_shape_factor_moments_gauss_ldata_read_row;
  data->ldata_write_row        = &_nc_galaxy_shape_factor_moments_gauss_ldata_noop;
  data->ldata_required_columns = &_nc_galaxy_shape_factor_moments_gauss_ldata_required_columns;
}

static void
_nc_galaxy_shape_factor_moments_gauss_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset)
{
  NcGalaxyShapeFactorMomentsGaussPrivate * const self =
    nc_galaxy_shape_factor_moments_gauss_get_instance_private (NC_GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS (gsf));

  /* No capability gate: only nc_galaxy_shape_pop_eval_p_array() and
   * nc_galaxy_shape_pop_moment_2k() are used, which every NcGalaxyShapePop
   * provides. */
  self->pop_hash = nc_galaxy_shape_factor_get_pop_hash (gsf);

  if (self->tab_cache_hash != self->pop_hash)
  {
    g_hash_table_remove_all (self->tab_cache);
    self->tab_cache_hash = self->pop_hash;
  }
}

/* The cold path of peek_table(), kept out of line so the warm check
 * inlines into every evaluation. */
G_GNUC_NO_INLINE static void
_nc_galaxy_shape_factor_moments_gauss_refresh_table (NcGalaxyShapeFactorMomentsGaussPrivate * const self,
                                                     NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data)
{
  NcGalaxyShapeFactorMomentsGaussLData *ldata = (NcGalaxyShapeFactorMomentsGaussLData *) data->ldata;

  const NcGalaxyShapeFactorMomentsKey key = { 0.0, data->pop_data->e_rms, NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_PANELS };
  NcGalaxyShapeFactorMomentsTable *table;

  /* Lookup and insertion under the lock, the build outside it; the
   * reference is taken while the lock is held, so a concurrent insert of
   * the same key cannot free the table first. */
  g_mutex_lock (&self->cache_lock);
  table = g_hash_table_lookup (self->tab_cache, &key);

  if (table != NULL)
    table = _nc_galaxy_shape_factor_moments_table_ref (table);

  g_mutex_unlock (&self->cache_lock);

  if (table == NULL)
  {
    NcmSpectral **sp                       = ncm_memory_pool_get (self->spectral_pool);
    NcGalaxyShapeFactorMomentsTable *built = _moments_gauss_build_table (self, *sp, pop, data->pop_data);

    ncm_memory_pool_return (sp);
    g_atomic_int_inc (&self->table_build_count);

    g_mutex_lock (&self->cache_lock);
    table = g_hash_table_lookup (self->tab_cache, &key);

    /* LCOV_EXCL_START: only when another thread built the same table meanwhile. */
    if (table != NULL)
    {
      table = _nc_galaxy_shape_factor_moments_table_ref (table);
      _nc_galaxy_shape_factor_moments_table_unref (built);
    }
    /* LCOV_EXCL_STOP */
    else
    {
      NcGalaxyShapeFactorMomentsKey *key_copy = g_new (NcGalaxyShapeFactorMomentsKey, 1);

      *key_copy = key;
      g_hash_table_insert (self->tab_cache, key_copy, built);
      table = _nc_galaxy_shape_factor_moments_table_ref (built);
    }

    g_mutex_unlock (&self->cache_lock);
  }

  g_clear_pointer (&ldata->table, _nc_galaxy_shape_factor_moments_table_unref);
  ldata->table      = table; /* already a reference */
  ldata->table_span = _nc_galaxy_shape_factor_moments_table_span (table,
                                                                  NC_GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS_PREFETCH_PANELS,
                                                                  NC_GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS_PREFETCH_CAP);
  ldata->pop_hash_seen = self->pop_hash;
  ldata->valid         = TRUE;
}

static inline const NcGalaxyShapeFactorMomentsTable *
_nc_galaxy_shape_factor_moments_gauss_peek_table (NcGalaxyShapeFactorMomentsGaussPrivate * const self,
                                                  NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data)
{
  NcGalaxyShapeFactorMomentsGaussLData *ldata = (NcGalaxyShapeFactorMomentsGaussLData *) data->ldata;

  if (G_UNLIKELY (!ldata->valid || (ldata->pop_hash_seen != self->pop_hash)))
    _nc_galaxy_shape_factor_moments_gauss_refresh_table (self, pop, data);

  return ldata->table;
}

/* Source moments at the folded shear, from the table: the three
 * interleaved series in one pass over the panel's block, which is the
 * evaluation hot path. */
static inline void
_nc_galaxy_shape_factor_moments_gauss_moments (const NcGalaxyShapeFactorMomentsTable *table, const gdouble ghat, gdouble *m)
{
  const guint j       = _nc_galaxy_shape_factor_moments_panel_index (table, ghat);
  const gdouble t     = _nc_galaxy_shape_factor_moments_panel_arg (table, j, ghat);
  const gdouble *coef = &table->coef[table->off[j]];
  const guint d       = table->deg[j];
  const gdouble t2    = 2.0 * t;
  gdouble a1 = 0.0, a2 = 0.0, b1 = 0.0, b2 = 0.0, c1 = 0.0, c2 = 0.0;
  guint k;

  /* x0 = t2 x1 + (ck - x2): ck and x2 are known before x1 is, so only
   * the multiply-add sits on each loop-carried chain. */
  for (k = d; k >= 1; k--)
  {
    const gdouble *ck = &coef[NC_GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS_N_COMP * k];
    const gdouble a0  = t2 * a1 + (ck[0] - a2);
    const gdouble b0  = t2 * b1 + (ck[1] - b2);
    const gdouble c0  = t2 * c1 + (ck[2] - c2);

    a2 = a1;
    a1 = a0;
    b2 = b1;
    b1 = b0;
    c2 = c1;
    c1 = c0;
  }

  m[0] = t * a1 - a2 + coef[0];
  m[1] = t * b1 - b2 + coef[1];
  m[2] = t * c1 - c2 + coef[2];
}

/*
 * Gauge-fixes (g, eps_obs) together by -arg(g), folds |g| to
 * ghat = min(|g|, 1/|g|) -- exact in both conventions for an isotropic
 * population -- and evaluates the Gaussian. The folded map conjugates the
 * observed ellipticity for |g| > 1, i.e. flips y, which this density is even
 * in.
 */
static gdouble
_nc_galaxy_shape_factor_moments_gauss_eval (const NcGalaxyShapeFactorMomentsTable *table, const gdouble sn2,
                                            const gdouble g_1, const gdouble g_2,
                                            const gdouble epsilon_obs_1, const gdouble epsilon_obs_2,
                                            const gboolean want_log)
{
  const gdouble g_mag  = sqrt (g_1 * g_1 + g_2 * g_2);
  const gdouble cos_pg = (g_mag > 0.0) ? g_1 / g_mag : 1.0;
  const gdouble sin_pg = (g_mag > 0.0) ? g_2 / g_mag : 0.0;
  const gdouble x      = epsilon_obs_1 * cos_pg + epsilon_obs_2 * sin_pg;
  const gdouble y      = -epsilon_obs_1 * sin_pg + epsilon_obs_2 * cos_pg;
  const gdouble ghat   = (g_mag <= 1.0) ? g_mag : 1.0 / g_mag;
  gdouble m[NC_GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS_N_COMP];
  gdouble Ct, Cx, dx, Q;

  _nc_galaxy_shape_factor_moments_gauss_moments (table, ghat, m);

  Ct = m[1] + sn2;
  Cx = m[2] + sn2;
  dx = x - m[0];
  Q  = dx * dx / Ct + y * y / Cx;

  if (want_log)
    return -0.5 * Q - 0.5 * log (Ct * Cx) - M_LN2 - M_LNPI;
  else
    return exp (-0.5 * Q) / (2.0 * M_PI * sqrt (Ct * Cx));
}

static void
_nc_galaxy_shape_factor_moments_gauss_data_prefetch (NcGalaxyShapeFactor *gsf, NcGalaxyShapeFactorData *data, const guint stage)
{
  const NcGalaxyShapeFactorMomentsGaussLData *ldata = (const NcGalaxyShapeFactorMomentsGaussLData *) data->ldata;

  if (stage == 1)
    ncm_prefetch_span (ldata, sizeof (NcGalaxyShapeFactorMomentsGaussLData));
  else if (stage == 2)
    ncm_prefetch_span (ldata->table, ldata->table_span);
}

static gdouble
_nc_galaxy_shape_factor_moments_gauss_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorMomentsGaussPrivate * const self =
    nc_galaxy_shape_factor_moments_gauss_get_instance_private (NC_GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS (gsf));
  const NcGalaxyShapeFactorMomentsTable *table = _nc_galaxy_shape_factor_moments_gauss_peek_table (self, pop, data);

  return _nc_galaxy_shape_factor_moments_gauss_eval (table, data->std_noise * data->std_noise, g_1, g_2, epsilon_obs_1, epsilon_obs_2, FALSE);
}

static gdouble
_nc_galaxy_shape_factor_moments_gauss_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorMomentsGaussPrivate * const self =
    nc_galaxy_shape_factor_moments_gauss_get_instance_private (NC_GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS (gsf));
  const NcGalaxyShapeFactorMomentsTable *table = _nc_galaxy_shape_factor_moments_gauss_peek_table (self, pop, data);

  return _nc_galaxy_shape_factor_moments_gauss_eval (table, data->std_noise * data->std_noise, g_1, g_2, epsilon_obs_1, epsilon_obs_2, TRUE);
}

static gchar *
_nc_galaxy_shape_factor_moments_gauss_get_desc (NcGalaxyShapeFactor *gsf)
{
  NcGalaxyShapeFactorMomentsGaussPrivate * const self =
    nc_galaxy_shape_factor_moments_gauss_get_instance_private (NC_GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS (gsf));
  gchar *parent_desc = NC_GALAXY_SHAPE_FACTOR_CLASS (nc_galaxy_shape_factor_moments_gauss_parent_class)->get_desc (gsf);
  gchar *desc        = g_strdup_printf ("%s, moment_tol=%g, max_degree=%u", parent_desc, self->moment_tol, self->max_degree);

  g_free (parent_desc);

  return desc;
}

static void
_nc_galaxy_shape_factor_moments_gauss_dispose (GObject *object)
{
  NcGalaxyShapeFactorMomentsGaussPrivate * const self =
    nc_galaxy_shape_factor_moments_gauss_get_instance_private (NC_GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS (object));

  g_clear_pointer (&self->tab_cache, g_hash_table_unref);

  if (self->spectral_pool != NULL)
  {
    ncm_memory_pool_free (self->spectral_pool, TRUE);
    self->spectral_pool = NULL;
  }

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_moments_gauss_parent_class)->dispose (object);
}

static void
_nc_galaxy_shape_factor_moments_gauss_finalize (GObject *object)
{
  NcGalaxyShapeFactorMomentsGaussPrivate * const self =
    nc_galaxy_shape_factor_moments_gauss_get_instance_private (NC_GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS (object));

  g_mutex_clear (&self->cache_lock);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_moments_gauss_parent_class)->finalize (object);
}

static void
nc_galaxy_shape_factor_moments_gauss_class_init (NcGalaxyShapeFactorMomentsGaussClass *klass)
{
  NcGalaxyShapeFactorClass *gsf_class = NC_GALAXY_SHAPE_FACTOR_CLASS (klass);
  GObjectClass *object_class          = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_shape_factor_moments_gauss_set_property;
  object_class->get_property = &_nc_galaxy_shape_factor_moments_gauss_get_property;
  object_class->constructed  = &_nc_galaxy_shape_factor_moments_gauss_constructed;
  object_class->dispose      = &_nc_galaxy_shape_factor_moments_gauss_dispose;
  object_class->finalize     = &_nc_galaxy_shape_factor_moments_gauss_finalize;

  /**
   * NcGalaxyShapeFactorMomentsGauss:moment-tol:
   *
   * Absolute tolerance to which each tabulated moment is resolved. The
   * error it leaves in $\ln P$ scales as $1/\sigma_\nu^2$, so the default of
   * $10^{-10}$ is set far below what any realistic $\sigma_\nu$ can
   * notice; a node costs one radial quadrature, so there is little to gain
   * from loosening it.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_MOMENT_TOL,
                                   g_param_spec_double ("moment-tol",
                                                        "Moment tolerance",
                                                        "Absolute tolerance of the tabulated moments",
                                                        1.0e-15, 1.0e-2, 1.0e-10,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyShapeFactorMomentsGauss:max-degree:
   *
   * Largest local degree a panel may reach. Refinement doubles from four,
   * so the default of 32 allows the ladder 4, 8, 16, 32.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_MAX_DEGREE,
                                   g_param_spec_uint ("max-degree",
                                                      "Maximum local degree",
                                                      "Largest Chebyshev degree a single panel may reach",
                                                      4, 64, 32,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  gsf_class->data_init = &_nc_galaxy_shape_factor_moments_gauss_data_init;

  /* The table cache is locked, and everything else a galaxy touches is its
   * own. */
  gsf_class->prepare          = &_nc_galaxy_shape_factor_moments_gauss_prepare;
  gsf_class->eval_marginal    = &_nc_galaxy_shape_factor_moments_gauss_eval_marginal;
  gsf_class->data_prefetch    = &_nc_galaxy_shape_factor_moments_gauss_data_prefetch;
  gsf_class->eval_ln_marginal = &_nc_galaxy_shape_factor_moments_gauss_eval_ln_marginal;
  gsf_class->get_desc         = &_nc_galaxy_shape_factor_moments_gauss_get_desc;
}

/**
 * nc_galaxy_shape_factor_moments_gauss_new:
 * @ellip_conv: the ellipticity convention #NcGalaxyWLObsEllipConv
 *
 * Creates a new #NcGalaxyShapeFactorMomentsGauss.
 *
 * Returns: (transfer full): a new #NcGalaxyShapeFactorMomentsGauss
 */
NcGalaxyShapeFactorMomentsGauss *
nc_galaxy_shape_factor_moments_gauss_new (NcGalaxyWLObsEllipConv ellip_conv)
{
  return g_object_new (NC_TYPE_GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS,
                       "ellip-conv", ellip_conv,
                       NULL);
}

/**
 * nc_galaxy_shape_factor_moments_gauss_ref:
 * @gsfmg: a #NcGalaxyShapeFactorMomentsGauss
 *
 * Increases the reference count of @gsfmg by one.
 *
 * Returns: (transfer full): @gsfmg
 */
NcGalaxyShapeFactorMomentsGauss *
nc_galaxy_shape_factor_moments_gauss_ref (NcGalaxyShapeFactorMomentsGauss *gsfmg)
{
  return g_object_ref (gsfmg);
}

/**
 * nc_galaxy_shape_factor_moments_gauss_free:
 * @gsfmg: a #NcGalaxyShapeFactorMomentsGauss
 *
 * Decreases the reference count of @gsfmg by one.
 *
 */
void
nc_galaxy_shape_factor_moments_gauss_free (NcGalaxyShapeFactorMomentsGauss *gsfmg)
{
  g_object_unref (gsfmg);
}

/**
 * nc_galaxy_shape_factor_moments_gauss_clear:
 * @gsfmg: a #NcGalaxyShapeFactorMomentsGauss
 *
 * If *@gsfmg is not %NULL, decreases its reference count by one and sets
 * *@gsfmg to %NULL.
 *
 */
void
nc_galaxy_shape_factor_moments_gauss_clear (NcGalaxyShapeFactorMomentsGauss **gsfmg)
{
  g_clear_object (gsfmg);
}

/**
 * nc_galaxy_shape_factor_moments_gauss_get_table_build_count:
 * @gsfmg: a #NcGalaxyShapeFactorMomentsGauss
 *
 * Number of tables actually built since the last reset, as opposed to
 * served from the instance cache. Under a parallel
 * nc_data_cluster_wl_factor_data_prepare(), two threads can build the same
 * table at once and one copy is discarded; both are counted.
 *
 * Returns: the running count of table builds
 */
guint
nc_galaxy_shape_factor_moments_gauss_get_table_build_count (NcGalaxyShapeFactorMomentsGauss *gsfmg)
{
  NcGalaxyShapeFactorMomentsGaussPrivate * const self = nc_galaxy_shape_factor_moments_gauss_get_instance_private (gsfmg);

  return (guint) g_atomic_int_get (&self->table_build_count);
}

/**
 * nc_galaxy_shape_factor_moments_gauss_reset_table_build_count:
 * @gsfmg: a #NcGalaxyShapeFactorMomentsGauss
 *
 * Zeroes the table-build counter.
 *
 */
void
nc_galaxy_shape_factor_moments_gauss_reset_table_build_count (NcGalaxyShapeFactorMomentsGauss *gsfmg)
{
  NcGalaxyShapeFactorMomentsGaussPrivate * const self = nc_galaxy_shape_factor_moments_gauss_get_instance_private (gsfmg);

  g_atomic_int_set (&self->table_build_count, 0);
}

/**
 * nc_galaxy_shape_factor_moments_gauss_exact_moments:
 * @gsfmg: a #NcGalaxyShapeFactorMomentsGauss
 * @pop: a #NcGalaxyShapePop
 * @data: a #NcGalaxyShapeFactorData
 * @ghat: the folded shear $\hat g=\min(g,1/g)\in[0,1]$
 * @mu: (out): $\mathrm{E}[x]$
 * @var_x: (out): $\mathrm{Var}(x)$, noise included
 * @e_y2: (out): $\mathrm{E}[y^2]$, noise included
 *
 * The exact moments of the observed marginal at @ghat, computed directly
 * rather than from the table. Exposed for validation against
 * nc_galaxy_shape_factor_moments_gauss_eval_moments().
 *
 */
void
nc_galaxy_shape_factor_moments_gauss_exact_moments (NcGalaxyShapeFactorMomentsGauss *gsfmg, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble ghat, gdouble *mu, gdouble *var_x, gdouble *e_y2)
{
  NcGalaxyShapeFactorMomentsGaussPrivate * const self = nc_galaxy_shape_factor_moments_gauss_get_instance_private (gsfmg);
  GArray *r_arr                                       = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), NC_GALAXY_SHAPE_FACTOR_MOMENTS_MOMENT_NNODES);
  GArray *p_arr                                       = NULL;
  const gdouble sn2                                   = data->std_noise * data->std_noise;
  gdouble t[3];

  g_array_set_size (r_arr, NC_GALAXY_SHAPE_FACTOR_MOMENTS_MOMENT_NNODES);
  _nc_galaxy_shape_factor_moments_exact_t (self->ellip_conv, pop, data->pop_data, ghat, 0.0, r_arr, &p_arr, t);

  *mu    = t[0];
  *var_x = t[1] - t[0] * t[0] + sn2;
  *e_y2  = t[2] + sn2;

  g_array_unref (r_arr);

  if (p_arr != NULL)
    g_array_unref (p_arr);
}

/**
 * nc_galaxy_shape_factor_moments_gauss_eval_moments:
 * @gsfmg: a #NcGalaxyShapeFactorMomentsGauss
 * @pop: a #NcGalaxyShapePop
 * @data: a #NcGalaxyShapeFactorData
 * @g_1: first reduced-shear component
 * @g_2: second reduced-shear component
 * @mu: (out): $\mathrm{E}[x]$
 * @var_x: (out): $\mathrm{Var}(x)$, noise included
 * @e_y2: (out): $\mathrm{E}[y^2]$, noise included
 *
 * The moments the evaluation uses at @g_1, @g_2, read from the table and
 * building it if needed: the same definitions as
 * nc_galaxy_shape_factor_moments_gauss_exact_moments(), so the two can be
 * compared directly.
 *
 */
void
nc_galaxy_shape_factor_moments_gauss_eval_moments (NcGalaxyShapeFactorMomentsGauss *gsfmg, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, gdouble *mu, gdouble *var_x, gdouble *e_y2)
{
  NcGalaxyShapeFactorMomentsGaussPrivate * const self = nc_galaxy_shape_factor_moments_gauss_get_instance_private (gsfmg);
  const NcGalaxyShapeFactorMomentsTable *table        = _nc_galaxy_shape_factor_moments_gauss_peek_table (self, pop, data);
  const gdouble g_mag                                 = hypot (g_1, g_2);
  const gdouble ghat                                  = (g_mag <= 1.0) ? g_mag : 1.0 / g_mag;
  const gdouble sn2                                   = data->std_noise * data->std_noise;
  gdouble m[NC_GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS_N_COMP];

  _nc_galaxy_shape_factor_moments_gauss_moments (table, ghat, m);

  *mu    = m[0];
  *var_x = m[1] + sn2;
  *e_y2  = m[2] + sn2;
}

/**
 * nc_galaxy_shape_factor_moments_gauss_peek_layout:
 * @gsfmg: a #NcGalaxyShapeFactorMomentsGauss
 * @pop: a #NcGalaxyShapePop
 * @data: a #NcGalaxyShapeFactorData
 * @n_panels: (out): number of panels
 * @top: (out): upper end of the mesh
 * @degrees: (out) (element-type gdouble) (transfer full): per-panel degree
 *
 * The mesh this galaxy's table uses, building it if needed.
 *
 */
void
nc_galaxy_shape_factor_moments_gauss_peek_layout (NcGalaxyShapeFactorMomentsGauss *gsfmg, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, guint *n_panels, gdouble *top, GArray **degrees)
{
  NcGalaxyShapeFactorMomentsGaussPrivate * const self = nc_galaxy_shape_factor_moments_gauss_get_instance_private (gsfmg);
  const NcGalaxyShapeFactorMomentsTable *table        = _nc_galaxy_shape_factor_moments_gauss_peek_table (self, pop, data);
  guint j;

  *n_panels = table->n_panels;
  *top      = table->top;
  *degrees  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), table->n_panels);

  for (j = 0; j < table->n_panels; j++)
  {
    const gdouble d = (gdouble) table->deg[j];

    g_array_append_val (*degrees, d);
  }
}

