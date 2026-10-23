/***************************************************************************
 *            nc_galaxy_shape_factor_moments_private.h
 *
 *  Mon Sep 22 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_moments_private.h
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

/*
 * Machinery shared by the shape factors built on the exact target moments
 * of the lensed intrinsic-ellipticity distribution
 * (#NcGalaxyShapeFactorMomentsTilt and #NcGalaxyShapeFactorMomentsGauss):
 * the closed-form target moments, the dyadic Chebyshev mesh in
 * ghat = min(|g|, 1/|g|) and its interleaved coefficient tables, the
 * instance-cache key, the tables' serialized form, and the reachable-range
 * bound. Not public API.
 */
#ifndef _NC_GALAXY_SHAPE_FACTOR_MOMENTS_PRIVATE_H_
#define _NC_GALAXY_SHAPE_FACTOR_MOMENTS_PRIVATE_H_

#ifndef NUMCOSMO_GIR_SCAN

#include <glib.h>
#include <string.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/background/nc_hicosmo.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_pop.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_wl_obs.h>
#include <numcosmo/nc/lss/halo/nc_halo_position.h>
#include <numcosmo/nc/lss/halo/nc_halo_density_profile.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/ncm/algebra/ncm_spectral.h>
#include "ncm/core/ncm_prefetch_private.h"

G_BEGIN_DECLS

/* Node count of the radial quadrature of the target moments: in the psi
 * variable of r = sin(psi)/delta for the TRACE convention, in r for
 * TRACE_DET. */
#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_MOMENT_NNODES 64

/* Largest number of panels. K = round(log2(1/sigma_nu)) is 7 at
 * sigma_nu = 0.008, the smallest per-galaxy shape dispersion in production
 * catalogs; the cap is one above that. */
#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_PANELS 8

/* Largest number of interleaved series a table holds. */
#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_COMP 4

/* Tolerance of the upper-end guard. */
#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_TOP_TOL 1.0e-12

/* Leading slots of a serialized table: n_comp, sigma_nu, e_rms, n_panels,
 * top, n_failed. */
#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_TABLES_HEADER 6

/*
 * A table: n_panels dyadic panels, panel j on [1 - 2^-j, 1 - 2^-(j+1)] except
 * the last, which ends at @top, each with its own Chebyshev series of degree
 * deg[j] for @n_comp interleaved series.
 *
 * off[j] indexes @coef in DOUBLES and already includes the factor n_comp:
 * off[0] = 0 and off[j+1] = off[j] + n_comp * (deg[j] + 1). Build and
 * evaluation have to agree on that, since a mismatch reads a valid but wrong
 * part of the block rather than running off its end.
 *
 * @coef points just past the struct, into the same allocation (see
 * _nc_galaxy_shape_factor_moments_table_alloc_coef()): a table is one block,
 * header then coefficients, and table_unref() frees it as one.
 *
 * lo[j] and width[j] = hi - lo are panel_bounds() cached at construction
 * (_nc_galaxy_shape_factor_moments_table_set_bounds()), so the evaluation
 * forms the panel argument without calling into libm.
 */
typedef struct _NcGalaxyShapeFactorMomentsTable
{
  guint n_comp;
  guint n_panels;
  gdouble top;
  guint8 deg[NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_PANELS];
  guint16 off[NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_PANELS];
  gdouble lo[NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_PANELS];
  gdouble width[NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_PANELS];
  gdouble *coef;
  guint n_coef;
  guint n_failed;
  gatomicrefcount ref_count;
} NcGalaxyShapeFactorMomentsTable;

/* Instance-cache key. Bit-exact equality: these are cache keys, not
 * physical comparisons. */
typedef struct _NcGalaxyShapeFactorMomentsKey
{
  gdouble sn;
  gdouble e_rms;
  guint n_panels;
} NcGalaxyShapeFactorMomentsKey;

void _nc_galaxy_shape_factor_moments_gl_fill (gdouble *x, gdouble *w, const guint n);
void _nc_galaxy_shape_factor_moments_exact_t (const NcGalaxyWLObsEllipConv ellip_conv, NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data, const gdouble ghat, const gdouble sn2, GArray *r_arr, GArray **p_arr, gdouble *t);
NcGalaxyShapeFactorMomentsTable *_nc_galaxy_shape_factor_moments_table_ref (NcGalaxyShapeFactorMomentsTable *table);
void _nc_galaxy_shape_factor_moments_table_unref (NcGalaxyShapeFactorMomentsTable *table);
void _nc_galaxy_shape_factor_moments_panel_bounds (const NcGalaxyShapeFactorMomentsTable *table, const guint j, gdouble *lo, gdouble *hi);
gdouble _nc_galaxy_shape_factor_moments_node (const gdouble lo, const gdouble hi, const guint d, const guint k);
void _nc_galaxy_shape_factor_moments_clenshaw_raw (const guint n_comp, const gdouble *coef, const guint d, const gdouble t, gdouble *out);
void _nc_galaxy_shape_factor_moments_table_set_bounds (NcGalaxyShapeFactorMomentsTable *table);
void _nc_galaxy_shape_factor_moments_clenshaw (const NcGalaxyShapeFactorMomentsTable *table, const guint j, const gdouble t, gdouble *out);
guint _nc_galaxy_shape_factor_moments_key_hash (gconstpointer p);
gboolean _nc_galaxy_shape_factor_moments_key_equal (gconstpointer a, gconstpointer b);
NcmVector *_nc_galaxy_shape_factor_moments_table_to_vector (const NcGalaxyShapeFactorMomentsKey *key, const NcGalaxyShapeFactorMomentsTable *table);
gboolean _nc_galaxy_shape_factor_moments_table_from_vector (NcmVector *v, NcGalaxyShapeFactorMomentsKey *key, NcGalaxyShapeFactorMomentsTable **table_out);
gdouble _nc_galaxy_shape_factor_moments_centre_delta (NcHaloPosition *hp, NcHICosmo *cosmo);
gboolean _nc_galaxy_shape_factor_moments_set_box_corner (NcHaloDensityProfile *dp_copy, NcHaloDensityProfile *dp_src);

/*
 * Panel index of @ghat, in O(1) and without a search: panel j covers
 * 1 - ghat in (2^-(j+1), 2^-j], so j is minus the frexp() exponent of
 * u = 1 - ghat. That exponent is read from the bits directly: u is a
 * difference of doubles in [0, 1], hence a multiple of 2^-53 and never
 * subnormal, so frexp's e is the biased exponent minus 1022.
 *
 * A ghat exactly on a breakpoint selects the lower panel; a breakpoint is a
 * Lobatto endpoint of that panel, where its series is exact. At ghat = 1,
 * u = 0 is not a normal double (frexp(0) reports exponent zero, which the
 * formula would map to panel 0), so that case is branched explicitly to the
 * last panel.
 */
static inline guint
_nc_galaxy_shape_factor_moments_panel_index (const NcGalaxyShapeFactorMomentsTable *table, const gdouble ghat)
{
  const gdouble u = 1.0 - ghat;

  if (u <= 0.0)
  {
    return table->n_panels - 1;
  }
  else
  {
    guint64 bits;
    gint jj;

    memcpy (&bits, &u, sizeof (bits));
    jj = 1022 - (gint) ((bits >> 52) & 0x7ff);

    if (jj < 0)
      jj = 0;

    return MIN ((guint) jj, table->n_panels - 1);
  }
}

/* Bytes from the start of @table through the end of its first @n_used
 * panels' coefficients (clamped to the table), capped at @cap: the span a
 * per-galaxy prefetch fetches. Reads the header, so it is computed when the
 * table is assigned to a galaxy, never at prefetch time. */
static inline gsize
_nc_galaxy_shape_factor_moments_table_span (const NcGalaxyShapeFactorMomentsTable *table, const guint n_used, const gsize cap)
{
  const guint j     = MIN (MAX (n_used, 1), table->n_panels) - 1;
  const gsize bytes = sizeof (NcGalaxyShapeFactorMomentsTable) + (table->off[j] + table->n_comp * (table->deg[j] + 1)) * sizeof (gdouble);

  return MIN (bytes, cap);
}

/* Argument of panel @j's series at @ghat, mapped to [-1, 1]. Every panel
 * has positive width: built tables by construction, loaded ones by the check
 * in _nc_galaxy_shape_factor_moments_table_from_vector(). */
static inline gdouble
_nc_galaxy_shape_factor_moments_panel_arg (const NcGalaxyShapeFactorMomentsTable *table, const guint j, const gdouble ghat)
{
  return 2.0 * (ghat - table->lo[j]) / table->width[j] - 1.0;
}

/*
 * ---- The generic table builder ----
 *
 * The build is driven by NcmSpectral, which places the Chebyshev-Lobatto
 * nodes, refines 2^k -> 2^(k+1) by evaluating only the new odd nodes, tests
 * convergence per component, and runs the DCT. The cost of a node and the
 * error a coefficient contributes differ between schemes, so those are
 * callbacks:
 *
 * - @node fills the n_comp values at a folded shear. Every node it has
 *   produced so far is in @memo, sorted, so it can take its starting point
 *   from the solved neighbors. Returns FALSE when it could not produce an
 *   accurate value; the value it wrote is kept and the failure counted.
 * - @panel_ctx fills a small per-panel context (e.g. the target moments at
 *   the panel midpoint) passed back to @probe.
 * - @probe measures the error a coefficient vector (or a nodal difference)
 *   of n_comp entries contributes to ln P. It is used for the truncation and
 *   its verification after NcmSpectral has converged.
 */
typedef struct _NcGalaxyShapeFactorMomentsMemo NcGalaxyShapeFactorMomentsMemo;

#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_MEMO_NONE G_MAXUINT

typedef gboolean (*NcGalaxyShapeFactorMomentsNodeF) (gpointer user_data, const gdouble ghat, const NcGalaxyShapeFactorMomentsMemo *memo, gdouble *vals);
typedef void (*NcGalaxyShapeFactorMomentsPanelCtxF) (gpointer user_data, const gdouble lo, const gdouble hi, gdouble *ctx);
typedef gdouble (*NcGalaxyShapeFactorMomentsProbeF) (gpointer user_data, const gdouble *v, const gdouble *ctx);

#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_CTX_LEN 4

guint _nc_galaxy_shape_factor_moments_memo_len (const NcGalaxyShapeFactorMomentsMemo *memo);
gdouble _nc_galaxy_shape_factor_moments_memo_x (const NcGalaxyShapeFactorMomentsMemo *memo, const guint i);
const gdouble *_nc_galaxy_shape_factor_moments_memo_vals (const NcGalaxyShapeFactorMomentsMemo *memo, const guint i);
guint _nc_galaxy_shape_factor_moments_memo_below (const NcGalaxyShapeFactorMomentsMemo *memo, const gdouble x);

NcGalaxyShapeFactorMomentsTable *_nc_galaxy_shape_factor_moments_build (NcmSpectral *spectral, const guint n_comp, const guint n_panels, const gdouble top, const guint max_degree, const gdouble abstol, const gdouble trunc_tol, NcGalaxyShapeFactorMomentsNodeF node, NcGalaxyShapeFactorMomentsPanelCtxF panel_ctx, NcGalaxyShapeFactorMomentsProbeF probe, gpointer user_data);

G_END_DECLS

#endif /* NUMCOSMO_GIR_SCAN */

#endif /* _NC_GALAXY_SHAPE_FACTOR_MOMENTS_PRIVATE_H_ */

