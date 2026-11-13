/***************************************************************************
 *            nc_galaxy_redshift_factor_spline.c
 *
 *  Sun Jul 13 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_factor_spline.c
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
 * NcGalaxyRedshiftFactorSpline:
 *
 * Spline redshift calculator scheme: a per-galaxy tabulated $p(z)$.
 *
 * Each galaxy carries its own table of $p(z)$ as an #NcmSpline (set via
 * nc_galaxy_redshift_factor_spline_data_set(), stored on the #NcGalaxyWLObs row via
 * nc_galaxy_wl_obs_set_pz()/nc_galaxy_wl_obs_peek_pz()). The scheme uses only the knots
 * $z_j$ and values $p_j \geq 0$ of that spline, not its interpolant: the density is
 * $$
 * p(z) = s(z)^2,
 * $$
 * where $s$ is the order-6 (degree 5) #NcmSplineBSpline through $(z_j, \sqrt{p_j})$.
 * The density is non-negative and $C^4$ everywhere on the table; its logarithm is
 * $2\ln\vert s(z)\vert$. The integration limits are the first and last knots, the
 * normalization is the exact integral of $s^2$ (6-point Gauss-Legendre on each polynomial
 * piece), and sampling inverts the same integral.
 *
 * Reading a row with a negative value, or with fewer than 6 knots, aborts.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_redshift_factor_spline.h"
#include "ncm/spline/ncm_spline_bspline.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <gsl/gsl_integration.h>
#endif /* NUMCOSMO_GIR_SCAN */

#define _NC_GALAXY_REDSHIFT_FACTOR_SPLINE_ORDER (6)
#define _NC_GALAXY_REDSHIFT_FACTOR_SPLINE_GL_N (6)

typedef struct _NcGalaxyRedshiftFactorSplinePrivate
{
  gint placeholder;
} NcGalaxyRedshiftFactorSplinePrivate;

struct _NcGalaxyRedshiftFactorSpline
{
  NcGalaxyRedshiftFactor parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyRedshiftFactorSpline, nc_galaxy_redshift_factor_spline, NC_TYPE_GALAXY_REDSHIFT_FACTOR);

/* Per-galaxy fragment: the table @pz as given, the B-spline @s of sqrt(p) through its
 * knots, the breakpoints @brk of @s, and @cum, the integral of s^2 from brk[0] to each
 * breakpoint (cum[0] = 0, the last entry is the normalization). */
typedef struct _SplineLData
{
  NcmSpline *pz;
  NcmSpline *s;
  NcmVector *brk;
  GArray *cum;
} SplineLData;

static void _spline_ldata_set_pz (SplineLData *sldata, NcmSpline *pz);
static gdouble _spline_s2_integ (NcmSpline *s, const gdouble a, const gdouble b);
static const gsl_integration_glfixed_table *_spline_gl (void);

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

  /* prepare, get_hash and update_data keep the base no-ops: the density depends on no
   * model, and everything derived from the table is built when the row is read. */
}

/* Fragment vtable --------------------------------------------------------- */

static void
_spline_ldata_free (gpointer ldata)
{
  SplineLData *sldata = (SplineLData *) ldata;

  ncm_spline_clear (&sldata->pz);
  ncm_spline_clear (&sldata->s);
  ncm_vector_clear (&sldata->brk);
  g_clear_pointer (&sldata->cum, g_array_unref);
  g_free (sldata);
}

static void
_spline_ldata_read_row (NcGalaxyRedshiftFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  SplineLData *sldata = (SplineLData *) data->ldata;

  _spline_ldata_set_pz (sldata, nc_galaxy_wl_obs_peek_pz (obs, i));
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
  /* The spline is carried on the NcGalaxyWLObs row's own dedicated pz slot, not a
   * regular scalar column. */
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

/* Density s^2 ---------------------------------------------------------------- */

/* Sets @pz, builds @s through (z_j, sqrt(p_j)), its breakpoints and the cumulative
 * integrals of s^2. Aborts on a negative value or fewer than ORDER knots. */
static void
_spline_ldata_set_pz (SplineLData *sldata, NcmSpline *pz)
{
  NcmVector *xv = ncm_spline_peek_xv (pz);
  NcmVector *yv = ncm_spline_peek_yv (pz);
  const guint n = ncm_vector_len (yv);
  NcmVector *sv;
  gdouble total = 0.0;
  guint j;

  if (n < _NC_GALAXY_REDSHIFT_FACTOR_SPLINE_ORDER)
    g_error ("nc_galaxy_redshift_factor_spline: the P(z) table has %u knots, at least %d are required.",
             n, _NC_GALAXY_REDSHIFT_FACTOR_SPLINE_ORDER);

  sv = ncm_vector_new (n);

  for (j = 0; j < n; j++)
  {
    const gdouble p_j = ncm_vector_get (yv, j);

    if (!(p_j >= 0.0))
      g_error ("nc_galaxy_redshift_factor_spline: the P(z) table value %u is %g at z = %g; values must be non-negative.",
               j, p_j, ncm_vector_get (xv, j));

    ncm_vector_fast_set (sv, j, sqrt (p_j));
  }

  ncm_spline_clear (&sldata->s);
  ncm_vector_clear (&sldata->brk);
  g_clear_pointer (&sldata->cum, g_array_unref);
  ncm_spline_ref (pz);
  ncm_spline_clear (&sldata->pz);

  sldata->pz  = pz;
  sldata->s   = NCM_SPLINE (ncm_spline_bspline_new_full (_NC_GALAXY_REDSHIFT_FACTOR_SPLINE_ORDER, xv, sv, TRUE));
  sldata->brk = ncm_spline_bspline_get_breakpoints (NCM_SPLINE_BSPLINE (sldata->s));
  ncm_spline_bspline_release_workspace (NCM_SPLINE_BSPLINE (sldata->s));
  sldata->cum = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), ncm_vector_len (sldata->brk));

  g_array_append_val (sldata->cum, total);

  for (j = 0; j + 1 < ncm_vector_len (sldata->brk); j++)
  {
    total += _spline_s2_integ (sldata->s, ncm_vector_get (sldata->brk, j), ncm_vector_get (sldata->brk, j + 1));
    g_array_append_val (sldata->cum, total);
  }

  ncm_vector_free (sv);
}

/* Integral of s^2 over [a, b] inside one polynomial piece of s, where s^2 has degree 10
 * and the 6-point Gauss-Legendre rule is exact. */
static gdouble
_spline_s2_integ (NcmSpline *s, const gdouble a, const gdouble b)
{
  const gsl_integration_glfixed_table *gl = _spline_gl ();
  const gdouble c                         = 0.5 * (a + b);
  const gdouble h                         = 0.5 * (b - a);
  gdouble res                             = 0.0;
  guint j;

  for (j = 0; j < _NC_GALAXY_REDSHIFT_FACTOR_SPLINE_GL_N / 2; j++)
  {
    const gdouble s_m = ncm_spline_eval (s, c - h * gl->x[j]);
    const gdouble s_p = ncm_spline_eval (s, c + h * gl->x[j]);

    res += gl->w[j] * (s_m * s_m + s_p * s_p);
  }

  return h * res;
}

static const gsl_integration_glfixed_table *
_spline_gl (void)
{
  static gsize init                        = 0;
  static gsl_integration_glfixed_table *gl = NULL;

  if (g_once_init_enter (&init))
  {
    gl = gsl_integration_glfixed_table_alloc (_NC_GALAXY_REDSHIFT_FACTOR_SPLINE_GL_N);
    g_once_init_leave (&init, 1);
  }

  return gl;
}

/* Sampling ----------------------------------------------------------------- */

/* Draws z with density s^2 / norm: picks the polynomial piece by its share of the
 * normalization, then solves int_{brk_j}^{z} s^2 = target in it by Newton steps kept
 * inside a bisection bracket. */
static void
_nc_galaxy_redshift_factor_spline_gen (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, NcmRNG *rng)
{
  SplineLData *sldata  = (SplineLData *) data->ldata;
  const guint n_brk    = ncm_vector_len (sldata->brk);
  const gdouble total  = g_array_index (sldata->cum, gdouble, n_brk - 1);
  const gdouble target = ncm_rng_uniform_gen (rng, 0.0, total);
  guint lo             = 0;
  guint hi             = n_brk - 1;
  gdouble a, b, z, rem;
  guint iter;

  g_assert_nonnull (sldata->s);

  while (hi - lo > 1)
  {
    const guint mid = (lo + hi) / 2;

    if (target < g_array_index (sldata->cum, gdouble, mid))
      hi = mid;
    else
      lo = mid;
  }

  a   = ncm_vector_get (sldata->brk, lo);
  b   = ncm_vector_get (sldata->brk, lo + 1);
  rem = target - g_array_index (sldata->cum, gdouble, lo);

  {
    const gdouble x0 = a;
    gdouble z_lo     = a;
    gdouble z_hi     = b;

    z = 0.5 * (a + b);

    for (iter = 0; iter < 200; iter++)
    {
      const gdouble F  = _spline_s2_integ (sldata->s, x0, z) - rem;
      const gdouble sz = ncm_spline_eval (sldata->s, z);
      const gdouble dF = sz * sz;
      gdouble z_new;

      if (F > 0.0)
        z_hi = z;
      else
        z_lo = z;

      z_new = (dF > 0.0) ? z - F / dF : 0.5 * (z_lo + z_hi);

      if (!((z_new > z_lo) && (z_new < z_hi)))
        z_new = 0.5 * (z_lo + z_hi);

      if (fabs (z_new - z) <= 4.0 * GSL_DBL_EPSILON * fabs (z))
      {
        z = z_new;
        break;
      }

      z = z_new;
    }
  }

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
  const gdouble s_z   = ncm_spline_eval (sldata->s, z);

  return s_z * s_z;
}

static gdouble
_spline_ln_integ_f (gpointer callback_data, const gdouble z, NcGalaxyRedshiftFactorData *data)
{
  SplineLData *sldata = (SplineLData *) data->ldata;

  return 2.0 * log (fabs (ncm_spline_eval (sldata->s, z)));
}

/* No closure payload needed. */
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

  g_assert_nonnull (sldata->cum);

  return g_array_index (sldata->cum, gdouble, sldata->cum->len - 1);
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
 * Sets @data's per-galaxy $p(z)$ table to the knots and values of @pz, see
 * #NcGalaxyRedshiftFactorSpline.
 *
 */
void
nc_galaxy_redshift_factor_spline_data_set (NcGalaxyRedshiftFactorSpline *gsdrs, NcGalaxyRedshiftFactorData *data, NcmSpline *pz)
{
  SplineLData *sldata = (SplineLData *) data->ldata;

  _spline_ldata_set_pz (sldata, pz);
}

/**
 * nc_galaxy_redshift_factor_spline_data_peek:
 * @gsdrs: a #NcGalaxyRedshiftFactorSpline
 * @data: a #NcGalaxyRedshiftFactorData
 *
 * Returns: (transfer none): @data's per-galaxy $p(z)$ table, as set.
 */
NcmSpline *
nc_galaxy_redshift_factor_spline_data_peek (NcGalaxyRedshiftFactorSpline *gsdrs, NcGalaxyRedshiftFactorData *data)
{
  SplineLData *sldata = (SplineLData *) data->ldata;

  return sldata->pz;
}

