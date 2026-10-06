/***************************************************************************
 *            nc_galaxy_shape_factor_moments_private.c
 *
 *  Mon Sep 22 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_moments_private.c
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_factor_moments_private.h"
#include "nc/lss/halo/nc_halo_mass_summary.h"
#include "nc/lss/wl/nc_wl_surface_mass_density.h"
#include "ncm/algebra/ncm_matrix.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#endif /* NUMCOSMO_GIR_SCAN */

/*
 * ---- Gauss-Legendre node table on [-1,1], built once ----
 *
 * The node count is a compile-time constant, so the table depends on nothing
 * per galaxy: gsl_integration_glfixed_table_alloc() must not be called per
 * evaluation (it allocates).
 */
static gdouble _nc_galaxy_shape_factor_moments_gl64_x[NC_GALAXY_SHAPE_FACTOR_MOMENTS_MOMENT_NNODES];
static gdouble _nc_galaxy_shape_factor_moments_gl64_w[NC_GALAXY_SHAPE_FACTOR_MOMENTS_MOMENT_NNODES];

static void _nc_galaxy_shape_factor_moments_gl64_init (void);

void
_nc_galaxy_shape_factor_moments_gl_fill (gdouble *x, gdouble *w, const guint n)
{
  gsl_integration_glfixed_table *table = gsl_integration_glfixed_table_alloc (n);
  guint i;

  for (i = 0; i < n; i++)
    gsl_integration_glfixed_point (-1.0, 1.0, i, &x[i], &w[i], table);

  gsl_integration_glfixed_table_free (table);
}

/*
 * ---- Stage 1: the exact target moments ----
 *
 * (E[x], E[x^2], E[y^2]) of the exact lensed-and-noisy marginal at the
 * folded shear ghat = min(|g|, 1/|g|), in the gauge frame where g is real.
 * The noise variance @sn2 is included in the second moments; pass zero for
 * the moments of the lensed source alone. The fold is exact in both
 * conventions for an isotropic population (see the TRACE_DET branch).
 *
 * @r_arr and @p_arr are caller-owned scratch reused across the nodes of one
 * table build (r_arr sized to MOMENT_NNODES; *p_arr may start NULL).
 */

/*
 * TRACE (chi) convention. The chi map depends on g only through
 * delta = 2 ghat / (1 + ghat^2), and the circle averages are closed form by
 * the Moebius lemma once r = sin(psi)/delta is substituted; the psi
 * substitution and the grouping of `num` and `wc2` are what keep the
 * integrand finite and cancellation-free as delta -> 1.
 */
static void
_nc_galaxy_shape_factor_moments_exact_t_trace (NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data,
                                               const gdouble ghat, const gdouble sn2,
                                               GArray *r_arr, GArray **p_arr, gdouble *t)
{
  const guint n_nodes = NC_GALAXY_SHAPE_FACTOR_MOMENTS_MOMENT_NNODES;
  const gdouble gh2   = ghat * ghat;
  const gdouble d     = 2.0 * ghat / (1.0 + gh2);
  const gdouble c2    = gsl_pow_2 ((1.0 - gh2) / (1.0 + gh2));
  gdouble *r_data     = (gdouble *) r_arr->data;
  const gdouble *p_data;
  gdouble w_data[NC_GALAXY_SHAPE_FACTOR_MOMENTS_MOMENT_NNODES];
  gdouble psi_max, sum_mu, sum_u2a, sum_u2b, sum_v2;
  guint i;

  _nc_galaxy_shape_factor_moments_gl64_init ();

  if (d < 1.0e-12)
  {
    /* g = 0: the marginal is P_0 itself and the moments are the
     * population's own, with the noise added. */
    const gdouble M2 = nc_galaxy_shape_pop_moment_2k (pop, pop_data, 1);

    t[0] = 0.0;
    t[1] = 0.5 * M2 + sn2;
    t[2] = 0.5 * M2 + sn2;

    return;
  }

  /* r = sin(psi)/delta, w = cos(psi): absorbs the sqrt singularity L and B
   * develop at r -> 1 as delta -> 1, which is the whole reason g = 1 is
   * reachable at all here. */
  psi_max = asin (MIN (d, 1.0));

  /* sin and cos of the same psi side by side, so the compiler can fuse
   * them into one sincos call. */
  for (i = 0; i < n_nodes; i++)
  {
    const gdouble psi = 0.5 * psi_max * (_nc_galaxy_shape_factor_moments_gl64_x[i] + 1.0);

    r_data[i] = sin (psi) / d;
    w_data[i] = cos (psi);
  }

  nc_galaxy_shape_pop_eval_p_array (pop, pop_data, r_arr, p_arr);
  p_data = (const gdouble *) (*p_arr)->data;

  sum_mu  = 0.0;
  sum_u2a = 0.0;
  sum_u2b = 0.0;
  sum_v2  = 0.0;

  for (i = 0; i < n_nodes; i++)
  {
    const gdouble w  = w_data[i];
    const gdouble r  = r_data[i];
    const gdouble r2 = r * r;

    /* measure in psi: int P F dr = (1/delta) int P F w dpsi, the 1/delta
     * folded into the per-term factors below. */
    const gdouble P   = p_data[i] * 0.5 * psi_max * _nc_galaxy_shape_factor_moments_gl64_w[i];
    const gdouble num = 1.0 + c2 - r2;
    const gdouble wc2 = w + c2;

    sum_mu  += P * (num / wc2); /* delta * int P a dr */
    sum_u2a += P * num * num / (w * wc2 * wc2);
    sum_u2b += P * r2 / (w * w * (1.0 + w));
    sum_v2  += P * r2 / (1.0 + w);
  }

  t[0] = sum_mu;
  t[1] = d * sum_u2a + (c2 * c2 / d) * sum_u2b + sn2;
  t[2] = (c2 / d) * sum_v2 + sn2;
}

/*
 * TRACE_DET (epsilon) convention. For |ghat| <= 1 the map is the disc
 * automorphism w = (z + ghat)/(1 + ghat z), with z = r e^{i phi} the source
 * ellipticity and its pole at -1/ghat outside the disc. Over each source
 * circle:
 *
 * - w is analytic in z, so by the mean-value property <w> = ghat and
 *   <w^2> = ghat^2, whatever r;
 * - the automorphism identity 1 - |w|^2 = (1 - ghat^2)(1 - r^2)/|1 + ghat z|^2
 *   and the Poisson integral <|1 + ghat z|^-2> = 1/(1 - ghat^2 r^2) give
 *   <|w|^2> = A(r) = 1 - (1 - ghat^2)(1 - r^2)/(1 - ghat^2 r^2).
 *
 * With x = Re w and y = Im w, x^2 = (w^2 + conj(w)^2 + 2|w|^2)/4 and
 * y^2 = (2|w|^2 - w^2 - conj(w)^2)/4, so
 *
 *   E[x]   = ghat                              (exactly: epsilon is unbiased)
 *   E[x^2] = (ghat^2 + <A>)/2 + sn2
 *   E[y^2] = (<A> - ghat^2)/2 + sn2
 *
 * with <A> the average of A(r) over the population's radial marginal: one
 * one-dimensional quadrature. The denominator is written as
 * (1 - ghat^2) + ghat^2 (1 - r^2), with each factor as a product of
 * (1 - .)(1 + .), so nothing cancels as ghat or r approaches one.
 *
 * The fold to ghat <= 1 is exact here too. With g real, the |g| > 1 map
 * (1 + g conj(eps))/(conj(eps) + g) is the |g| <= 1 map at 1/g applied to
 * conj(eps), and conjugation preserves an isotropic population's measure.
 * It conjugates the observed ellipticity, i.e. flips y, which no moment used
 * here can see.
 */
static void
_nc_galaxy_shape_factor_moments_exact_t_trace_det (NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data,
                                                   const gdouble ghat, const gdouble sn2,
                                                   GArray *r_arr, GArray **p_arr, gdouble *t)
{
  const guint n_nodes    = NC_GALAXY_SHAPE_FACTOR_MOMENTS_MOMENT_NNODES;
  const gdouble gh2      = ghat * ghat;
  const gdouble one_m_g2 = (1.0 - ghat) * (1.0 + ghat);
  gdouble *r_data        = (gdouble *) r_arr->data;
  const gdouble *p_data;
  gdouble mean_A = 0.0;
  guint i;

  _nc_galaxy_shape_factor_moments_gl64_init ();

  for (i = 0; i < n_nodes; i++)
    r_data[i] = 0.5 * (_nc_galaxy_shape_factor_moments_gl64_x[i] + 1.0);

  nc_galaxy_shape_pop_eval_p_array (pop, pop_data, r_arr, p_arr);
  p_data = (const gdouble *) (*p_arr)->data;

  for (i = 0; i < n_nodes; i++)
  {
    const gdouble r        = r_data[i];
    const gdouble one_m_r2 = (1.0 - r) * (1.0 + r);
    const gdouble den      = one_m_g2 + gh2 * one_m_r2;
    const gdouble A        = (den > 0.0) ? 1.0 - one_m_g2 * one_m_r2 / den : 1.0;

    mean_A += 0.5 * _nc_galaxy_shape_factor_moments_gl64_w[i] * p_data[i] * A;
  }

  t[0] = ghat;
  t[1] = 0.5 * (gh2 + mean_A) + sn2;
  t[2] = 0.5 * (mean_A - gh2) + sn2;
}

void
_nc_galaxy_shape_factor_moments_exact_t (const NcGalaxyWLObsEllipConv ellip_conv,
                                         NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data,
                                         const gdouble ghat, const gdouble sn2,
                                         GArray *r_arr, GArray **p_arr, gdouble *t)
{
  switch (ellip_conv)
  {
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
      _nc_galaxy_shape_factor_moments_exact_t_trace (pop, pop_data, ghat, sn2, r_arr, p_arr, t);
      break;
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
      _nc_galaxy_shape_factor_moments_exact_t_trace_det (pop, pop_data, ghat, sn2, r_arr, p_arr, t);
      break;
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
  }
}

NcGalaxyShapeFactorMomentsTable *
_nc_galaxy_shape_factor_moments_table_ref (NcGalaxyShapeFactorMomentsTable *table)
{
  g_atomic_ref_count_inc (&table->ref_count);

  return table;
}

void
_nc_galaxy_shape_factor_moments_table_unref (NcGalaxyShapeFactorMomentsTable *table)
{
  if (g_atomic_ref_count_dec (&table->ref_count))
    g_free (table);  /* coef lives in the same block */
}

void
_nc_galaxy_shape_factor_moments_panel_bounds (const NcGalaxyShapeFactorMomentsTable *table, const guint j, gdouble *lo, gdouble *hi)
{
  *lo = 1.0 - ldexp (1.0, -(gint) j);
  *hi = (j + 1 == table->n_panels) ? table->top : 1.0 - ldexp (1.0, -(gint) (j + 1));
}

/* Ascending Chebyshev-Lobatto node k of degree d, mapped to [lo, hi]:
 * k = 0 is lo and k = d is hi. Under d -> 2d node k becomes node 2k, which
 * is what makes refinement reuse every solve already done. */
gdouble
_nc_galaxy_shape_factor_moments_node (const gdouble lo, const gdouble hi, const guint d, const guint k)
{
  const gdouble t = cos (M_PI * (gdouble) (d - k) / (gdouble) d);

  return lo + 0.5 * (hi - lo) * (t + 1.0);
}

/* The @n_comp series separately, for the validation entry points and for
 * the build's own check that a truncation did not move the nodes. */
void
_nc_galaxy_shape_factor_moments_clenshaw_raw (const guint n_comp, const gdouble *coef, const guint d, const gdouble t, gdouble *out)
{
  const gdouble t2 = 2.0 * t;
  guint a;

  for (a = 0; a < n_comp; a++)
  {
    gdouble b1 = 0.0;
    gdouble b2 = 0.0;
    guint k;

    for (k = d; k >= 1; k--)
    {
      const gdouble b0 = t2 * b1 - b2 + coef[n_comp * k + a];

      b2 = b1;
      b1 = b0;
    }

    out[a] = t * b1 - b2 + coef[a];
  }
}

/* Grows *@table to hold @n_coef coefficients right after the header and
 * points coef there: one allocation per table, so a table can be prefetched
 * as a single span and its coefficients are never a pointer hop away from
 * the header. Moves the block; nothing may hold the old pointer. */
static void
_nc_galaxy_shape_factor_moments_table_alloc_coef (NcGalaxyShapeFactorMomentsTable **table, const guint n_coef)
{
  NcGalaxyShapeFactorMomentsTable *t = g_realloc (*table, sizeof (NcGalaxyShapeFactorMomentsTable) + n_coef * sizeof (gdouble));

  G_STATIC_ASSERT (sizeof (NcGalaxyShapeFactorMomentsTable) % sizeof (gdouble) == 0);

  t->n_coef = n_coef;
  t->coef   = (gdouble *) (t + 1);
  *table    = t;
}

/* Caches panel_bounds() in @table; called once n_panels and top are set. */
void
_nc_galaxy_shape_factor_moments_table_set_bounds (NcGalaxyShapeFactorMomentsTable *table)
{
  guint j;

  for (j = 0; j < table->n_panels; j++)
  {
    gdouble lo, hi;

    _nc_galaxy_shape_factor_moments_panel_bounds (table, j, &lo, &hi);
    table->lo[j]    = lo;
    table->width[j] = hi - lo;
  }
}

void
_nc_galaxy_shape_factor_moments_clenshaw (const NcGalaxyShapeFactorMomentsTable *table, const guint j, const gdouble t, gdouble *out)
{
  _nc_galaxy_shape_factor_moments_clenshaw_raw (table->n_comp, &table->coef[table->off[j]], table->deg[j], t, out);
}

guint
_nc_galaxy_shape_factor_moments_key_hash (gconstpointer p)
{
  const NcGalaxyShapeFactorMomentsKey *k = (const NcGalaxyShapeFactorMomentsKey *) p;

  return g_double_hash (&k->sn) ^ (g_double_hash (&k->e_rms) * 31U) ^ (k->n_panels * 41U);
}

gboolean
_nc_galaxy_shape_factor_moments_key_equal (gconstpointer a, gconstpointer b)
{
  const NcGalaxyShapeFactorMomentsKey *ka = (const NcGalaxyShapeFactorMomentsKey *) a;
  const NcGalaxyShapeFactorMomentsKey *kb = (const NcGalaxyShapeFactorMomentsKey *) b;

  /* Bit-exact: these are cache keys, not physical comparisons. */
  return (ka->sn == kb->sn) && (ka->e_rms == kb->e_rms) && (ka->n_panels == kb->n_panels);
}

NcmVector *
_nc_galaxy_shape_factor_moments_table_to_vector (const NcGalaxyShapeFactorMomentsKey *key, const NcGalaxyShapeFactorMomentsTable *table)
{
  const guint len = NC_GALAXY_SHAPE_FACTOR_MOMENTS_TABLES_HEADER + table->n_panels + table->n_coef;
  NcmVector *v    = ncm_vector_new (len);
  guint i;

  ncm_vector_set (v, 0, (gdouble) table->n_comp);
  ncm_vector_set (v, 1, key->sn);
  ncm_vector_set (v, 2, key->e_rms);
  ncm_vector_set (v, 3, (gdouble) table->n_panels);
  ncm_vector_set (v, 4, table->top);
  ncm_vector_set (v, 5, (gdouble) table->n_failed);

  for (i = 0; i < table->n_panels; i++)
    ncm_vector_set (v, NC_GALAXY_SHAPE_FACTOR_MOMENTS_TABLES_HEADER + i, (gdouble) table->deg[i]);

  for (i = 0; i < table->n_coef; i++)
    ncm_vector_set (v, NC_GALAXY_SHAPE_FACTOR_MOMENTS_TABLES_HEADER + table->n_panels + i, table->coef[i]);

  return v;
}

/* Rebuilds one table from its vector, rejecting anything whose shape does
 * not add up rather than trusting the length. Returns FALSE and leaves the
 * caller to rebuild. */
gboolean
_nc_galaxy_shape_factor_moments_table_from_vector (NcmVector *v, NcGalaxyShapeFactorMomentsKey *key, NcGalaxyShapeFactorMomentsTable **table_out)
{
  NcGalaxyShapeFactorMomentsTable *table;
  guint n_comp, n_panels, n_coef, off, i;

  if (ncm_vector_len (v) < NC_GALAXY_SHAPE_FACTOR_MOMENTS_TABLES_HEADER + 1)
    return FALSE;

  n_comp   = (guint) ncm_vector_get (v, 0);
  n_panels = (guint) ncm_vector_get (v, 3);

  if ((n_comp < 1) || (n_comp > NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_COMP))
    return FALSE;

  if ((n_panels < 1) || (n_panels > NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_PANELS))
    return FALSE;

  /* The last panel starts at 1 - 2^-(n_panels - 1) and must have width. */
  if (!((ncm_vector_get (v, 4) > 1.0 - ldexp (1.0, -(gint) (n_panels - 1))) && (ncm_vector_get (v, 4) <= 1.0)))
    return FALSE;

  table           = g_new0 (NcGalaxyShapeFactorMomentsTable, 1);
  table->n_comp   = n_comp;
  table->n_panels = n_panels;
  table->top      = ncm_vector_get (v, 4);
  table->n_failed = (guint) ncm_vector_get (v, 5);
  g_atomic_ref_count_init (&table->ref_count);
  _nc_galaxy_shape_factor_moments_table_set_bounds (table);

  off = 0;

  for (i = 0; i < n_panels; i++)
  {
    const gdouble d = ncm_vector_get (v, NC_GALAXY_SHAPE_FACTOR_MOMENTS_TABLES_HEADER + i);

    if ((d < 1.0) || (d > 255.0))
    {
      _nc_galaxy_shape_factor_moments_table_unref (table);

      return FALSE;
    }

    table->deg[i] = (guint8) d;
    table->off[i] = (guint16) off;
    off          += n_comp * (table->deg[i] + 1);
  }

  n_coef = off;

  if (ncm_vector_len (v) != NC_GALAXY_SHAPE_FACTOR_MOMENTS_TABLES_HEADER + n_panels + n_coef)
  {
    _nc_galaxy_shape_factor_moments_table_unref (table);

    return FALSE;
  }

  _nc_galaxy_shape_factor_moments_table_alloc_coef (&table, n_coef);

  for (i = 0; i < n_coef; i++)
    table->coef[i] = ncm_vector_get (v, NC_GALAXY_SHAPE_FACTOR_MOMENTS_TABLES_HEADER + n_panels + i);

  key->sn       = ncm_vector_get (v, 1);
  key->e_rms    = ncm_vector_get (v, 2);
  key->n_panels = n_panels;

  *table_out = table;

  return TRUE;
}

/*
 * Largest displacement the halo centre's prior permits, as a physical
 * transverse length at the cluster redshift. Zero when the centre is fixed.
 * The corners are fed through nc_halo_position_projected_radius_from_ra_dec(),
 * which is the same mapping the pipeline uses for a galaxy's own radius, so
 * the cos(dec) factor and the angular-diameter conversion come with it.
 *
 * This is a bound over the PRIOR BOX, not over the posterior: walkers are
 * initialised across the whole box, so a bound taken from the current point
 * would not hold for the chain. Its being a property of the prior is also
 * what lets a table be built once and stay valid throughout.
 */
gdouble
_nc_galaxy_shape_factor_moments_centre_delta (NcHaloPosition *hp, NcHICosmo *cosmo)
{
  NcmModel *model         = NCM_MODEL (hp);
  const gboolean ra_free  = ncm_model_param_get_ftype (model, NC_HALO_POSITION_RA) == NCM_PARAM_TYPE_FREE;
  const gboolean dec_free = ncm_model_param_get_ftype (model, NC_HALO_POSITION_DEC) == NCM_PARAM_TYPE_FREE;
  gdouble ra_c[2], dec_c[2];
  gdouble delta = 0.0;
  guint i, j;

  if (!ra_free && !dec_free)
    return 0.0;

  ra_c[0] = ra_free ? ncm_model_param_get_lower_bound (model, NC_HALO_POSITION_RA) : ncm_model_param_get (model, NC_HALO_POSITION_RA);
  ra_c[1] = ra_free ? ncm_model_param_get_upper_bound (model, NC_HALO_POSITION_RA) : ra_c[0];

  dec_c[0] = dec_free ? ncm_model_param_get_lower_bound (model, NC_HALO_POSITION_DEC) : ncm_model_param_get (model, NC_HALO_POSITION_DEC);
  dec_c[1] = dec_free ? ncm_model_param_get_upper_bound (model, NC_HALO_POSITION_DEC) : dec_c[0];

  for (i = 0; i < 2; i++)
  {
    for (j = 0; j < 2; j++)
      delta = MAX (delta, nc_halo_position_projected_radius_from_ra_dec (hp, cosmo, ra_c[i], dec_c[j]));
  }

  return delta;
}

/* Sets @dp_copy's mass to the top of its prior box and its concentration to
 * @c_bound, leaving the profile in the mset untouched. Returns FALSE when
 * the parameters cannot be resolved by name, which is not an error: the
 * caller then skips the restriction and builds the full mesh. */
gboolean
_nc_galaxy_shape_factor_moments_set_box_corner (NcHaloDensityProfile *dp_copy, NcHaloDensityProfile *dp_src, const gboolean c_upper)
{
  NcHaloMassSummary *hms_copy = nc_halo_density_profile_peek_mass_summary (dp_copy);
  NcHaloMassSummary *hms_src  = nc_halo_density_profile_peek_mass_summary (dp_src);
  NcmModel *m_copy            = NCM_MODEL (hms_copy);
  NcmModel *m_src             = NCM_MODEL (hms_src);
  GError *error               = NULL;
  guint idx;

  /* Resolve by name: the parameter set differs per NcHaloMassSummary
   * implementation, and several have no free concentration at all. */
  /* LCOV_EXCL_START: every NcHaloMassSummary in the library has log10MDelta. */
  if (!ncm_model_param_index_from_name (m_src, "log10MDelta", &idx, &error))
  {
    g_clear_error (&error);

    return FALSE;
  }

  /* LCOV_EXCL_STOP */

  ncm_model_param_set (m_copy, idx, ncm_model_param_get_upper_bound (m_src, idx));

  if (ncm_model_param_index_from_name (m_src, "cDelta", &idx, &error))
    ncm_model_param_set (m_copy, idx, c_upper ? ncm_model_param_get_upper_bound (m_src, idx)
                                              : ncm_model_param_get_lower_bound (m_src, idx));
  else
    g_clear_error (&error);


  return TRUE;
}

static void
_nc_galaxy_shape_factor_moments_gl64_init (void)
{
  static gsize init = 0;

  if (g_once_init_enter (&init))
  {
    _nc_galaxy_shape_factor_moments_gl_fill (_nc_galaxy_shape_factor_moments_gl64_x,
                                             _nc_galaxy_shape_factor_moments_gl64_w,
                                             NC_GALAXY_SHAPE_FACTOR_MOMENTS_MOMENT_NNODES);
    g_once_init_leave (&init, 1);
  }
}

/*
 * ---- Memo of solved nodes ----
 *
 * Every node the build evaluates, sorted by folded shear. It serves three
 * purposes. NcmSpectral asks for nodes in DESCENDING order within a panel,
 * while a continuation wants to climb from ghat = 0, so the builder solves
 * each panel's first level itself, ascending, and NcmSpectral's requests
 * then land in the memo. A refinement node always has solved neighbours on
 * both sides to start from. And the first node of panel j is the last node
 * of panel j-1, so it is found here rather than solved twice.
 *
 * Lookups use a tolerance, not bitwise equality: NcmSpectral computes its
 * cosines from its own table, and the nodes are far apart compared to it.
 */
typedef struct _NcGalaxyShapeFactorMomentsMemoEntry
{
  gdouble x;
  gdouble v[NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_COMP];
} NcGalaxyShapeFactorMomentsMemoEntry;

struct _NcGalaxyShapeFactorMomentsMemo
{
  GArray *e;
};

guint
_nc_galaxy_shape_factor_moments_memo_len (const NcGalaxyShapeFactorMomentsMemo *memo)
{
  return memo->e->len;
}

gdouble
_nc_galaxy_shape_factor_moments_memo_x (const NcGalaxyShapeFactorMomentsMemo *memo, const guint i)
{
  return g_array_index (memo->e, NcGalaxyShapeFactorMomentsMemoEntry, i).x;
}

const gdouble *
_nc_galaxy_shape_factor_moments_memo_vals (const NcGalaxyShapeFactorMomentsMemo *memo, const guint i)
{
  return g_array_index (memo->e, NcGalaxyShapeFactorMomentsMemoEntry, i).v;
}

/* Index of the largest entry strictly below @x, or MEMO_NONE. */
guint
_nc_galaxy_shape_factor_moments_memo_below (const NcGalaxyShapeFactorMomentsMemo *memo, const gdouble x)
{
  guint lo = 0, hi = memo->e->len;

  while (lo < hi)
  {
    const guint mid = (lo + hi) / 2;

    if (g_array_index (memo->e, NcGalaxyShapeFactorMomentsMemoEntry, mid).x < x)
      lo = mid + 1;
    else
      hi = mid;
  }

  return (lo == 0) ? NC_GALAXY_SHAPE_FACTOR_MOMENTS_MEMO_NONE : lo - 1;
}

typedef struct _NcGalaxyShapeFactorMomentsBuildCtx
{
  NcGalaxyShapeFactorMomentsMemo memo;
  guint n_comp;
  gdouble tol_x;
  NcGalaxyShapeFactorMomentsNodeF node;
  gpointer user_data;
  guint n_failed;
  guint comp;
} NcGalaxyShapeFactorMomentsBuildCtx;

/* The values at @x, from the memo when already there, otherwise evaluated
 * and inserted in order. Copied out to @out, because an insertion can move
 * the array. */
static void
_nc_galaxy_shape_factor_moments_memo_get (NcGalaxyShapeFactorMomentsBuildCtx *ctx, const gdouble x, gdouble *out)
{
  const guint below = _nc_galaxy_shape_factor_moments_memo_below (&ctx->memo, x);
  const guint at    = (below == NC_GALAXY_SHAPE_FACTOR_MOMENTS_MEMO_NONE) ? 0 : below + 1;
  guint c;

  if ((at < ctx->memo.e->len) &&
      (fabs (g_array_index (ctx->memo.e, NcGalaxyShapeFactorMomentsMemoEntry, at).x - x) <= ctx->tol_x))
  {
    memcpy (out, g_array_index (ctx->memo.e, NcGalaxyShapeFactorMomentsMemoEntry, at).v, sizeof (gdouble) * ctx->n_comp);

    return;
  }

  if ((below != NC_GALAXY_SHAPE_FACTOR_MOMENTS_MEMO_NONE) &&
      (fabs (g_array_index (ctx->memo.e, NcGalaxyShapeFactorMomentsMemoEntry, below).x - x) <= ctx->tol_x))
  {
    memcpy (out, g_array_index (ctx->memo.e, NcGalaxyShapeFactorMomentsMemoEntry, below).v, sizeof (gdouble) * ctx->n_comp);

    return;
  }

  {
    NcGalaxyShapeFactorMomentsMemoEntry entry;

    entry.x = x;

    for (c = 0; c < NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_COMP; c++)
      entry.v[c] = 0.0;

    if (!ctx->node (ctx->user_data, x, &ctx->memo, entry.v))
      ctx->n_failed++;

    g_array_insert_val (ctx->memo.e, at, entry);
    memcpy (out, entry.v, sizeof (gdouble) * ctx->n_comp);
  }
}

static void
_nc_galaxy_shape_factor_moments_batch_f (gpointer user_data, gdouble x, NcmVector *y)
{
  NcGalaxyShapeFactorMomentsBuildCtx *ctx = (NcGalaxyShapeFactorMomentsBuildCtx *) user_data;
  gdouble v[NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_COMP];
  guint c;

  _nc_galaxy_shape_factor_moments_memo_get (ctx, x, v);

  for (c = 0; c < ctx->n_comp; c++)
    ncm_vector_set (y, c, v[c]);
}

static gdouble
_nc_galaxy_shape_factor_moments_scalar_f (gpointer user_data, gdouble x)
{
  NcGalaxyShapeFactorMomentsBuildCtx *ctx = (NcGalaxyShapeFactorMomentsBuildCtx *) user_data;
  gdouble v[NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_COMP];

  _nc_galaxy_shape_factor_moments_memo_get (ctx, x, v);

  return v[ctx->comp];
}

/*
 * Builds a table of @n_panels panels ending at @top.
 *
 * Per panel: solve the first level (degree 4) ascending, so the continuation
 * runs up the interval; let NcmSpectral refine until every component has
 * converged to @abstol, capped at @max_degree; then truncate to the smallest
 * degree whose tail, measured by @probe, is inside @trunc_tol, and confirm
 * that truncation at the panel's own nodes before accepting it.
 *
 * NcmSpectral's convergence test is per component and compares successive
 * levels, so it measures the error BETWEEN nodes, which is where an
 * unresolved series goes wrong. The probe alone cannot be trusted with that:
 * the exact tilt satisfies dW/dghat = t . dlambda/dghat, so over a panel
 * where t barely moves, a probe weighted by t sees the coefficients of W
 * cancel those of t . lambda, and reports a small error for a series that is
 * not resolved. The probe is used only to trim what NcmSpectral has already
 * resolved.
 */
NcGalaxyShapeFactorMomentsTable *
_nc_galaxy_shape_factor_moments_build (NcmSpectral *spectral, const guint n_comp, const guint n_panels, const gdouble top,
                                       const guint max_degree, const gdouble abstol, const gdouble trunc_tol,
                                       NcGalaxyShapeFactorMomentsNodeF node, NcGalaxyShapeFactorMomentsPanelCtxF panel_ctx,
                                       NcGalaxyShapeFactorMomentsProbeF probe, gpointer user_data)
{
  const guint k_min                      = 2; /* degree 4 */
  guint k_cap                            = k_min;
  NcGalaxyShapeFactorMomentsTable *table = g_new0 (NcGalaxyShapeFactorMomentsTable, 1);
  GArray *all_coef                       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcGalaxyShapeFactorMomentsBuildCtx ctx;
  guint j;

  g_assert_cmpuint (n_comp, >=, 1);
  g_assert_cmpuint (n_comp, <=, NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_COMP);
  g_assert_cmpuint (n_panels, >=, 1);
  g_assert_cmpuint (n_panels, <=, NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_PANELS);

  while ((1U << (k_cap + 1)) <= max_degree)
    k_cap++;

  table->n_comp   = n_comp;
  table->n_panels = n_panels;
  table->top      = top;
  _nc_galaxy_shape_factor_moments_table_set_bounds (table);
  g_atomic_ref_count_init (&table->ref_count);

  ctx.memo.e    = g_array_new (FALSE, FALSE, sizeof (NcGalaxyShapeFactorMomentsMemoEntry));
  ctx.n_comp    = n_comp;
  ctx.node      = node;
  ctx.user_data = user_data;
  ctx.n_failed  = 0;
  ctx.comp      = 0;

  for (j = 0; j < n_panels; j++)
  {
    gdouble lo, hi;
    gdouble pctx[NC_GALAXY_SHAPE_FACTOR_MOMENTS_CTX_LEN];
    gdouble *pcoef;
    NcmMatrix *cm = NULL;
    guint d, m, k, c;

    _nc_galaxy_shape_factor_moments_panel_bounds (table, j, &lo, &hi);
    ctx.tol_x = 1.0e-13 * (hi - lo);

    panel_ctx (user_data, lo, hi, pctx);

    /* First level, ascending, so each node starts from the one below it. */
    {
      const guint d0 = 1U << k_min;
      gdouble v[NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_COMP];

      for (k = 0; k <= d0; k++)
        _nc_galaxy_shape_factor_moments_memo_get (&ctx, _nc_galaxy_shape_factor_moments_node (lo, hi, d0, k), v);
    }

    if (ncm_spectral_compute_chebyshev_coeffs_batch_adaptive_cap (spectral, &_nc_galaxy_shape_factor_moments_batch_f,
                                                                  n_comp, lo, hi, k_min, k_cap, 0.0, abstol, FALSE,
                                                                  &cm, &ctx) > 0)
    {
      d = ncm_matrix_ncols (cm) - 1;
    }
    else
    {
      /* Not converged at the cap: take the cap level itself. Every node of
       * it is already in the memo, so this solves nothing. */
      d = 1U << k_cap;
      g_clear_pointer (&cm, ncm_matrix_free);
      cm = ncm_matrix_new (n_comp, d + 1);

      for (c = 0; c < n_comp; c++)
      {
        GArray *cc = NULL;

        ctx.comp = c;
        ncm_spectral_compute_chebyshev_coeffs (spectral, &_nc_galaxy_shape_factor_moments_scalar_f, lo, hi, d + 1, &cc, &ctx);

        for (k = 0; k <= d; k++)
          ncm_matrix_set (cm, c, k, g_array_index (cc, gdouble, k));

        g_array_unref (cc);
      }
    }

    pcoef = g_new (gdouble, n_comp * (d + 1));

    for (k = 0; k <= d; k++)
      for (c = 0; c < n_comp; c++)
        pcoef[n_comp * k + c] = ncm_matrix_get (cm, c, k);

    ncm_matrix_free (cm);

    /* Truncate to the smallest degree whose remaining tail is inside the
     * margin... */
    {
      gdouble tail = 0.0;

      m = d;

      while (m > 1)
      {
        const gdouble next = tail + probe (user_data, &pcoef[n_comp * m], pctx);

        if (next >= trunc_tol)
          break;

        tail = next;
        m--;
      }
    }

    /* ...then confirm it at the panel's own nodes rather than trusting the
     * tail estimate. The full series reproduces every node exactly, so any
     * disagreement there is the truncation's own, and the check solves
     * nothing: the nodes are all in the memo. */
    while (m < d)
    {
      gdouble worst = 0.0;

      for (k = 0; k <= d; k++)
      {
        const gdouble t_node = cos (M_PI * (gdouble) (d - k) / (gdouble) d);
        gdouble lw[NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_COMP];
        gdouble vn[NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_COMP];
        gdouble delta[NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_COMP];

        _nc_galaxy_shape_factor_moments_clenshaw_raw (n_comp, pcoef, m, t_node, lw);
        _nc_galaxy_shape_factor_moments_memo_get (&ctx, _nc_galaxy_shape_factor_moments_node (lo, hi, d, k), vn);

        for (c = 0; c < n_comp; c++)
          delta[c] = lw[c] - vn[c];

        worst = MAX (worst, probe (user_data, delta, pctx));
      }

      if (worst < trunc_tol)
        break;

      /* LCOV_EXCL_START: safety net; no scanned population, noise or gate
       * has had the tail estimate under-call the truncation error. */
      m++;
      /* LCOV_EXCL_STOP */
    }

    table->deg[j] = (guint8) m;
    table->off[j] = (guint16) all_coef->len;
    g_array_append_vals (all_coef, pcoef, n_comp * (m + 1));
    g_free (pcoef);
  }

  table->n_failed = ctx.n_failed;
  _nc_galaxy_shape_factor_moments_table_alloc_coef (&table, all_coef->len);
  memcpy (table->coef, all_coef->data, all_coef->len * sizeof (gdouble));
  g_array_unref (all_coef);

  g_array_unref (ctx.memo.e);

  return table;
}

