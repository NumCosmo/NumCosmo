/***************************************************************************
 *            nc_galaxy_redshift_factor_spline.c
 *
 *  Sun Jul 13 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_factor_spline.c
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
 * NcGalaxyRedshiftFactorSpline:
 *
 * Spline redshift calculator scheme: a per-galaxy pre-tabulated $p(z)$.
 *
 * Each galaxy carries its own $p(z)$ as an #NcmSpline (set via
 * nc_galaxy_redshift_factor_spline_data_set(), stored on the #NcGalaxyWLObs
 * row via nc_galaxy_wl_obs_set_pz()/nc_galaxy_wl_obs_peek_pz()); the joint
 * density this scheme hands out is exactly that spline, evaluated (or its
 * natural logarithm taken) at the true redshift $z$. Unlike
 * #NcGalaxyRedshiftFactorComposed, no #NcmMSet model is involved: the
 * spline is per-galaxy data, not a fitted population/observable pair, so
 * this scheme's density never changes except when the underlying data
 * (i.e. the spline itself) changes.
 *
 * Sampling draws from an #NcmStatsDist1dSpline built from @pz (lazily, on
 * the first nc_galaxy_redshift_factor_update_data() call after a fresh
 * @pz is set).
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_redshift_factor_spline.h"
#include "ncm/spline/ncm_spline_cubic_notaknot.h"
#include "ncm/stats/ncm_stats_dist1d_spline.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcGalaxyRedshiftFactorSplinePrivate
{
  gint placeholder;
} NcGalaxyRedshiftFactorSplinePrivate;

struct _NcGalaxyRedshiftFactorSpline
{
  NcGalaxyRedshiftFactor parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyRedshiftFactorSpline, nc_galaxy_redshift_factor_spline, NC_TYPE_GALAXY_REDSHIFT_FACTOR);

/* Per-galaxy fragment: the spline itself, plus the sampler built from it
 * (NULL until the first update_data() call after @pz is set). */
typedef struct _SplineLData
{
  NcmSpline *pz;
  NcmStatsDist1d *dist;
} SplineLData;

static void
nc_galaxy_redshift_factor_spline_init (NcGalaxyRedshiftFactorSpline *gsdrs)
{
}

static void
_nc_galaxy_redshift_factor_spline_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_redshift_factor_spline_parent_class)->finalize (object);
}

static void _nc_galaxy_redshift_factor_spline_data_init (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data);
static void _nc_galaxy_redshift_factor_spline_gen (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, NcmRNG *rng);
static gboolean _nc_galaxy_redshift_factor_spline_gen1 (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, NcmRNG *rng);
static NcGalaxyRedshiftFactorIntegrand *_nc_galaxy_redshift_factor_spline_integ (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, gboolean use_lnp);
static void _nc_galaxy_redshift_factor_spline_get_integ_lim (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, gdouble *z_min, gdouble *z_max);
static gdouble _nc_galaxy_redshift_factor_spline_norm (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data);
static NcmIntegralFixed *_nc_galaxy_redshift_factor_spline_make_fixed_nodes (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, gdouble z_lo, gdouble z_hi, guint n_nodes, guint rule_n);

static void
nc_galaxy_redshift_factor_spline_class_init (NcGalaxyRedshiftFactorSplineClass *klass)
{
  GObjectClass *object_class              = G_OBJECT_CLASS (klass);
  NcGalaxyRedshiftFactorClass *gsdr_class = NC_GALAXY_REDSHIFT_FACTOR_CLASS (klass);

  object_class->finalize = &_nc_galaxy_redshift_factor_spline_finalize;

  gsdr_class->data_init        = &_nc_galaxy_redshift_factor_spline_data_init;
  gsdr_class->gen              = &_nc_galaxy_redshift_factor_spline_gen;
  gsdr_class->gen1             = &_nc_galaxy_redshift_factor_spline_gen1;
  gsdr_class->integ            = &_nc_galaxy_redshift_factor_spline_integ;
  gsdr_class->get_integ_lim    = &_nc_galaxy_redshift_factor_spline_get_integ_lim;
  gsdr_class->norm             = &_nc_galaxy_redshift_factor_spline_norm;
  gsdr_class->make_fixed_nodes = &_nc_galaxy_redshift_factor_spline_make_fixed_nodes;

  /* prepare/get_hash: base defaults (no-op / constant 0) are exactly right
   * here -- this scheme has no NcmMSet model dependency, so its density
   * only ever changes when the per-galaxy data (the spline itself) does,
   * which the orchestrator's own obs_changed already covers.
   *
   * update_data: base default (no-op) is right too -- this scheme's only
   * per-galaxy cache derived from update_data(), the resample-only
   * inverse-CDF sampler, is now built lazily inside _gen() itself (see
   * _nc_galaxy_redshift_factor_spline_ensure_dist()), not eagerly on every
   * prepare() cycle. */
}

/* Fragment vtable --------------------------------------------------------- */

static void
_spline_ldata_free (gpointer ldata)
{
  SplineLData *sldata = (SplineLData *) ldata;

  ncm_spline_clear (&sldata->pz);
  ncm_stats_dist1d_clear (&sldata->dist);
  g_free (sldata);
}

static void
_spline_ldata_read_row (NcGalaxyRedshiftFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  SplineLData *sldata = (SplineLData *) data->ldata;

  ncm_spline_clear (&sldata->pz);
  ncm_stats_dist1d_clear (&sldata->dist);

  sldata->pz = ncm_spline_ref (nc_galaxy_wl_obs_peek_pz (obs, i));
}

static void
_spline_ldata_write_row (NcGalaxyRedshiftFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  SplineLData *sldata = (SplineLData *) data->ldata;

  nc_galaxy_wl_obs_set_pz (obs, i, sldata->pz);
}

static void
_spline_ldata_required_columns (NcGalaxyRedshiftFactorData *data, GList **columns)
{
  /* The spline is carried on the NcGalaxyWLObs row's own dedicated pz slot,
   * not a regular scalar column. */
}

static void
_nc_galaxy_redshift_factor_spline_data_init (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data)
{
  SplineLData *sldata = g_new0 (SplineLData, 1);

  data->ldata                  = sldata;
  data->ldata_destroy          = &_spline_ldata_free;
  data->ldata_read_row         = &_spline_ldata_read_row;
  data->ldata_write_row        = &_spline_ldata_write_row;
  data->ldata_required_columns = &_spline_ldata_required_columns;
}

/* Sampling ----------------------------------------------------------------- */

/* Builds sldata->dist (an inverse-CDF sampler over @pz) the first time a
 * galaxy is actually sampled -- resample()-only state, never read by
 * integ/norm/get_integ_lim/make_fixed_nodes (all of which use sldata->pz
 * directly), so building it eagerly for every galaxy on every prepare()
 * cycle (as update_data() used to do) paid its full cost -- dominated by
 * ncm_stats_dist1d_prepare()'s own inverse-CDF construction -- even for
 * fits/MCMC runs that never call resample() at all. */
static void
_nc_galaxy_redshift_factor_spline_ensure_dist (SplineLData *sldata)
{
  if (sldata->dist != NULL)
    return;

  {
    NcmVector *xv     = ncm_spline_peek_xv (sldata->pz);
    NcmVector *yv     = ncm_spline_peek_yv (sldata->pz);
    NcmVector *m2lnyv = ncm_vector_new (ncm_vector_len (yv));
    NcmSpline *m2lnp;
    gdouble z_min, z_max;
    guint j;

    ncm_spline_get_bounds (sldata->pz, &z_min, &z_max);

    for (j = 0; j < ncm_vector_len (yv); j++)
    {
      const gdouble y = -2.0 * log (ncm_vector_fast_get (yv, j) + 1.0e-5);

      ncm_vector_set (m2lnyv, j, y);
    }

    m2lnp = NCM_SPLINE (ncm_spline_cubic_notaknot_new_full (xv, m2lnyv, TRUE));

    sldata->dist = NCM_STATS_DIST1D (ncm_stats_dist1d_spline_new (m2lnp));
    g_object_set (G_OBJECT (sldata->dist), "reltol", 1.0e-5, NULL);
    ncm_stats_dist1d_set_xi (sldata->dist, z_min);
    ncm_stats_dist1d_set_xf (sldata->dist, z_max);
    ncm_stats_dist1d_prepare (sldata->dist);

    ncm_spline_free (m2lnp);
    ncm_vector_free (m2lnyv);
  }
}

static void
_nc_galaxy_redshift_factor_spline_gen (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, NcmRNG *rng)
{
  SplineLData *sldata = (SplineLData *) data->ldata;
  gdouble z_min, z_max, z;

  g_assert_nonnull (sldata->pz);
  _nc_galaxy_redshift_factor_spline_ensure_dist (sldata);
  ncm_spline_get_bounds (sldata->pz, &z_min, &z_max);

  do {
    z = ncm_stats_dist1d_gen (sldata->dist, rng);
  } while ((z < z_min) || (z > z_max));

  data->z = z;
}

static gboolean
_nc_galaxy_redshift_factor_spline_gen1 (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, NcmRNG *rng)
{
  _nc_galaxy_redshift_factor_spline_gen (gsdr, mset, data, rng);

  return TRUE;
}

/* Integrand ---------------------------------------------------------------- */

static gdouble
_spline_integ_f (gpointer callback_data, const gdouble z, NcGalaxyRedshiftFactorData *data)
{
  SplineLData *sldata = (SplineLData *) data->ldata;

  return ncm_spline_eval (sldata->pz, z);
}

static gdouble
_spline_ln_integ_f (gpointer callback_data, const gdouble z, NcGalaxyRedshiftFactorData *data)
{
  SplineLData *sldata = (SplineLData *) data->ldata;

  return log (ncm_spline_eval (sldata->pz, z));
}

/* No closure payload needed: the integrand reads the density straight from
 * @data->ldata at eval time, unlike Composed which must capture the
 * mset-resolved population/observable. callback_data stays NULL throughout;
 * these two just satisfy the callback machinery's non-NULL requirement on
 * the free/copy function pointers themselves. */
static gpointer
_spline_integ_data_copy (gpointer callback_data)
{
  return NULL;
}

static void
_spline_integ_data_free (gpointer callback_data)
{
}

static NcGalaxyRedshiftFactorIntegrand *
_nc_galaxy_redshift_factor_spline_integ (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, gboolean use_lnp)
{
  return nc_galaxy_redshift_factor_integrand_new (use_lnp ? _spline_ln_integ_f : _spline_integ_f,
                                                  &_spline_integ_data_free,
                                                  &_spline_integ_data_copy,
                                                  NULL,
                                                  NULL);
}

/* Integration limits and normalization -------------------------------------- */

static void
_nc_galaxy_redshift_factor_spline_get_integ_lim (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, gdouble *z_min, gdouble *z_max)
{
  SplineLData *sldata = (SplineLData *) data->ldata;

  g_assert_nonnull (sldata->pz);
  ncm_spline_get_bounds (sldata->pz, z_min, z_max);
}

static gdouble
_nc_galaxy_redshift_factor_spline_norm (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data)
{
  SplineLData *sldata = (SplineLData *) data->ldata;
  gdouble z_min, z_max;

  g_assert_nonnull (sldata->pz);
  ncm_spline_get_bounds (sldata->pz, &z_min, &z_max);

  return ncm_spline_eval_integ (sldata->pz, z_min, z_max);
}

static gdouble
_spline_fixed_nodes_gsl_f (gdouble z, gpointer user_data)
{
  NcGalaxyRedshiftFactorData *data = (NcGalaxyRedshiftFactorData *) user_data;

  return _spline_integ_f (NULL, z, data);
}

static NcmIntegralFixed *
_nc_galaxy_redshift_factor_spline_make_fixed_nodes (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, gdouble z_lo, gdouble z_hi, guint n_nodes, guint rule_n)
{
  NcmIntegralFixed *intf;
  gsl_function F;

  intf = ncm_integral_fixed_new (n_nodes, rule_n, z_lo, z_hi);

  F.function = &_spline_fixed_nodes_gsl_f;
  F.params   = data;

  ncm_integral_fixed_calc_nodes (intf, &F);

  return intf;
}

/* Public API ----------------------------------------------------------------- */

/**
 * nc_galaxy_redshift_factor_spline_new:
 *
 * Creates a new #NcGalaxyRedshiftFactorSpline calculator scheme.
 *
 * Returns: (transfer full): a new #NcGalaxyRedshiftFactorSpline.
 */
NcGalaxyRedshiftFactorSpline *
nc_galaxy_redshift_factor_spline_new (void)
{
  return g_object_new (NC_TYPE_GALAXY_REDSHIFT_FACTOR_SPLINE, NULL);
}

/**
 * nc_galaxy_redshift_factor_spline_ref:
 * @gsdrs: a #NcGalaxyRedshiftFactorSpline
 *
 * Increases the reference count of @gsdrs by one.
 *
 * Returns: (transfer full): @gsdrs.
 */
NcGalaxyRedshiftFactorSpline *
nc_galaxy_redshift_factor_spline_ref (NcGalaxyRedshiftFactorSpline *gsdrs)
{
  return g_object_ref (gsdrs);
}

/**
 * nc_galaxy_redshift_factor_spline_free:
 * @gsdrs: a #NcGalaxyRedshiftFactorSpline
 *
 * Decreases the reference count of @gsdrs by one.
 *
 */
void
nc_galaxy_redshift_factor_spline_free (NcGalaxyRedshiftFactorSpline *gsdrs)
{
  g_object_unref (gsdrs);
}

/**
 * nc_galaxy_redshift_factor_spline_clear:
 * @gsdrs: a #NcGalaxyRedshiftFactorSpline
 *
 * Decreases the reference count of @gsdrs by one, and sets the pointer
 * *@gsdrs to NULL.
 *
 */
void
nc_galaxy_redshift_factor_spline_clear (NcGalaxyRedshiftFactorSpline **gsdrs)
{
  g_clear_object (gsdrs);
}

/**
 * nc_galaxy_redshift_factor_spline_data_set:
 * @gsdrs: a #NcGalaxyRedshiftFactorSpline
 * @data: a #NcGalaxyRedshiftFactorData
 * @pz: an #NcmSpline
 *
 * Sets @data's per-galaxy $p(z)$ to @pz.
 *
 */
void
nc_galaxy_redshift_factor_spline_data_set (NcGalaxyRedshiftFactorSpline *gsdrs, NcGalaxyRedshiftFactorData *data, NcmSpline *pz)
{
  SplineLData *sldata = (SplineLData *) data->ldata;

  ncm_spline_clear (&sldata->pz);
  ncm_stats_dist1d_clear (&sldata->dist);

  sldata->pz = ncm_spline_ref (pz);
}

/**
 * nc_galaxy_redshift_factor_spline_data_peek:
 * @gsdrs: a #NcGalaxyRedshiftFactorSpline
 * @data: a #NcGalaxyRedshiftFactorData
 *
 * Returns: (transfer none): @data's per-galaxy $p(z)$.
 */
NcmSpline *
nc_galaxy_redshift_factor_spline_data_peek (NcGalaxyRedshiftFactorSpline *gsdrs, NcGalaxyRedshiftFactorData *data)
{
  SplineLData *sldata = (SplineLData *) data->ldata;

  return sldata->pz;
}

