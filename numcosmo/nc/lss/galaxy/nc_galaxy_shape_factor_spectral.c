/***************************************************************************
 *            nc_galaxy_shape_factor_spectral.c
 *
 *  Wed Nov 11 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_spectral.c
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
 * NcGalaxyShapeFactorSpectral:
 *
 * Intrinsic-ellipticity marginal by adaptive Chebyshev panels in polar coordinates.
 *
 * Computes
 * $$P(\epsilon_\mathrm{obs} \mid g) = \int_{|\chi_I|<1} \mathrm{d}^2\chi_I\,
 *   P_\mathrm{pop}(\chi_I)\, N_2\big(\epsilon_\mathrm{obs} - f_g(\chi_I);
 *   \sigma_\mathrm{noise}^2\big)$$
 * with #NcmSpectralPolar in the polar coordinates $(r, \theta)$ of the variable of
 * integration, for any isotropic population that implements
 * nc_galaxy_shape_pop_exponent_at_origin(), nc_galaxy_shape_pop_exponent_at_one() and
 * nc_galaxy_shape_pop_eval_p_reduced(). The two end power laws of the density of $r$ are
 * the weight of the integral and the reduced density times the noise kernel is its
 * smooth factor.
 *
 * The variable of integration is $\chi_I$ for %NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET. For
 * %NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE it is $\epsilon_I = T^{-1}(\chi_I)$,
 * $T(\epsilon) = 2\epsilon / (1 + |\epsilon|^2)$, the trace--determinant ellipticity of
 * the same shape, in which the shear map is the composition of the trace--determinant
 * map with $T$. In it the integrand varies on the scale $1 - |g|$ near $|g| = 1$, and in
 * $\chi_I$ on the scale $(1 - |g|)^2$.
 *
 * The shear enters through its modulus. At $|g| = 1$ every point of the open disk maps
 * to $g$ and $P = N_2(\epsilon_\mathrm{obs} - g; \sigma_\mathrm{noise}^2)$. For
 * $|g| > 1$ the integral equals the one at $g / |g|^2$, which is used. The frame is then
 * rotated by the phase of $g$, so that the shear is the real number $t = |g|$.
 *
 * For $1 - t < 10^{-5}$ the integral is not computed at $t$, where the integrand varies
 * on the scale $1 - t$ next to the unit circle. The value is the sum of the marginal to
 * second order in the distance of the sheared ellipticity to its mean and of a term
 * proportional to a power of $(1 - t) / (1 + t)$, whose coefficients are set by the
 * integral at $1 - t = 10^{-5}$, and for %NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET also at
 * $1 - t = 10^{-4}$. It assumes $10^{-5} |\epsilon_\mathrm{obs} - g| \ll
 * \sigma_\mathrm{noise}^2$.
 *
 * The panels are those of ncm_spectral_polar_peak_edges() for a peak located on the ray
 * through the preimage of $\epsilon_\mathrm{obs}$, at the radius that maximizes the
 * logarithm of the density of $r$ minus
 * $|\epsilon_\mathrm{obs} - f_t(\chi_I)|^2 / 2\sigma_\mathrm{noise}^2$ on a scan of that ray, with half widths $\sigma_\mathrm{noise}$ times the smaller
 * stretch factor of the inverse map there. The scan places the panels only: the integral is
 * accepted by the test of #NcmSpectralPolar with the relative tolerance
 * #NcGalaxyShapeFactorSpectral:reltol on panels that cover the whole disk.
 *
 * The class computes $\ln P$, with the smallest $|\epsilon_\mathrm{obs} - f_t(\chi_I)|^2$
 * of the scan taken out of the exponent, and nc_galaxy_shape_factor_eval_marginal()
 * returns its exponential. A panel that reaches #NcmSpectralPolar:max-depth, or a
 * non-finite $\ln P$, is an error that aborts with the point in the message.
 *
 * With #NcGalaxyShapeFactorSpectral:use-interp, nc_galaxy_shape_factor_eval_marginal() and
 * nc_galaxy_shape_factor_eval_ln_marginal() take $\ln P$ for
 * $|g| \le$ #NcGalaxyShapeFactorSpectral:interp-g-max from an approximation of $\ln P$ in
 * $g$ built for each galaxy by ncm_spectral_polar_build() from the integral, to the
 * tolerance #NcGalaxyShapeFactorSpectral:interp-abstol, and kept in the
 * #NcGalaxyShapeFactorData of the galaxy. nc_galaxy_shape_factor_data_prepare() builds it,
 * and so does the first evaluation of a galaxy that has none valid.
 *
 * The interpolants can be stored in the catalog: nc_galaxy_shape_factor_data_write_row()
 * writes them (#NcGalaxyWLObs:shape-interp) and
 * nc_galaxy_shape_factor_spectral_set_obs_source() records the shape factor and the
 * population that built them. Reading the rows of such a catalog loads the interpolants;
 * nc_galaxy_shape_factor_data_prepare() then compares the record with the fit (the
 * convention, reltol, the interp-* settings and the population parameters within a relative
 * tolerance of $10^{-15}$, none of them free; for an expansion over the population
 * parameters, see #NcGalaxyShapeFactorSpectral:pop-order) and aborts on any difference. A loaded
 * interpolant is never rebuilt: an evaluation that finds it invalid, or a catalog with
 * interpolants used without #NcGalaxyShapeFactorSpectral:use-interp, aborts.
 *
 * A call takes its scratch space from a pool, so one object may be used for several
 * galaxies in several threads. Two threads must not evaluate the same galaxy when
 * #NcGalaxyShapeFactorSpectral:use-interp is set, since the first evaluation builds its
 * interpolant.
 *
 * See <a href="../../theory/nc/lss/galaxy/wl_shape_marginalization_spectral.html">Shape
 * Marginalization by Spectral Panels</a> for the formulas.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_factor_spectral.h"
#include "nc/lss/galaxy/nc_galaxy_shape_factor_reach_private.h"
#include "nc/lss/wl/nc_wl_ellipticity.h"
#include "ncm/algebra/ncm_spectral.h"
#include "ncm/algebra/ncm_spectral_polar.h"
#include "ncm/algebra/ncm_spectral_polar_panels.h"
#include "ncm/algebra/ncm_spectral_polar_private.h"
#include "ncm/algebra/ncm_spectral_polar_panels_array.h"
#include "ncm/core/ncm_memory_pool.h"
#include "ncm/core/ncm_prefetch_private.h"
#include "ncm/core/ncm_serialize.h"
#include "ncm/core/ncm_vmath_private.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <string.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

/* The radii of the scan of the ray: a uniform grid, the two ends, and the peak of the noise */
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_UNIFORM (401)
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_END_DECADES (11)
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_OFFSETS (10)
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_EXTRA (2 * NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_END_DECADES + 4 * NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_OFFSETS)

/*
 * Below this distance of min (|g|, 1 / |g|) to 1 the integral is not computed at the shear;
 * the value is formed from integrals at this distance. The errors measured against the Arb
 * reference hold for this value.
 */
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NEAR_UNIT_DELTA (1.0e-5)

/*
 * The table along the line of the shears of a galaxy: the first grid has 2^LEVEL_MIN + 1
 * points, the finest 2^LEVEL_MAX + 1 (at most the largest level with a generated kernel,
 * NCM_SPECTRAL_CHEB_MAX_LEVEL, so that no FFTW plan is made inside the threaded
 * likelihood), and the largest |g| of the line it covers.
 */
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_LEVEL_MIN (1)
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_LEVEL_MAX (7)
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_G_MAX (0.9)

/* The smallest radius in |g| an interpolant is built for with restrict-range */
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_REACH_MIN (0.02)

/* The largest number of populations of the expansion of pop-order */
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_POP_NODES_MAX (4096)


/* Steps of the golden-section search of the maximum of the scan */
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_GOLDEN_STEPS (40)

struct _NcGalaxyShapeFactorSpectral
{
  NcGalaxyShapeFactor parent_instance;
};

typedef struct _NcGalaxyShapeFactorSpectralKernel NcGalaxyShapeFactorSpectralKernel;
typedef struct _NcGalaxyShapeFactorSpectralNodes NcGalaxyShapeFactorSpectralNodes;

/*
 * The box of an expansion over the population parameters: its order K, the indices of the
 * n_free parameters it spans and their bounds, n_nodes = K^n_free. Shared, read-only once
 * made.
 */
typedef struct _NcGalaxyShapeFactorSpectralBox
{
  guint K;
  guint n_free;
  guint n_nodes;
  guint *idx;
  gdouble *lb;
  gdouble *ub;
} NcGalaxyShapeFactorSpectralBox;

static void _nc_galaxy_shape_factor_spectral_nodes_free (NcGalaxyShapeFactorSpectralNodes *nodes);

/*
 * The functions of the convention act on the variable of integration z: chi_I for
 * TRACE_DET and T^{-1} (chi_I) for TRACE, T (e) = 2 e / (1 + |e|^2). Those of the integrals
 * take the convention from NcGalaxyShapeFactorSpectralKernel:trace, once per call.
 */
typedef struct _NcGalaxyShapeFactorSpectralPrivate
{
  /* the sheared ellipticity e_obs to z, at the shear g */
  complex double (*map_inv) (complex double g, complex double e_obs);

  gboolean axis_ratio_var; /* TRUE for TRACE: z = T^{-1} (chi_I) */

  NcmSpectralPolar *polar;        /* the integrals */
  NcmSpectralPolar *polar_interp; /* the build of the interpolant in g */
  NcmMemoryPool *work_pool;       /* NcGalaxyShapeFactorSpectralWork, one per concurrent evaluation */
  NcmMemoryPool *line_pool;       /* NcmSpectral of the tables along a line, one per concurrent eval_line */
  gdouble reltol;
  gdouble line_abstol;
  gboolean use_interp;
  gdouble interp_g_max;
  gdouble interp_abstol;
  guint interp_min_radial_panels;

  /*
   * restrict-range: the bound of the reachable range over the box of the fit, taken by
   * prepare() (_nc_galaxy_shape_factor_reach_prepare(), with a free cDelta sampled). It
   * applies when range_on: restrict-range and use-interp are set, the models are in the fit
   * and the box covers every free parameter. n_short counts the galaxies whose loaded
   * interpolant does not cover their reach in the box of the fit, warned about once per
   * catalog (_nc_galaxy_shape_factor_spectral_warn_once()).
   */
  gboolean restrict_range;
  gboolean range_on;
  NcGalaxyShapeFactorReach reach;
  gint n_short;

  /*
   * pop-order: the galaxies built by this object expand ln P over the box of the free
   * parameters of the population, pop_box (pop_order^n_free nodes, pop_nodes, made by
   * prepare() when the box changes, pop_box_key). pop_fit_key is the key of the free
   * parameters, their bounds and the values of the fixed ones, the state a loaded expansion
   * was compared with. pop_n_unconverged counts the galaxies whose expansion did not
   * converge since the box was set; they are integrated. A loaded expansion compared with a
   * fit whose population reaches outside its box is warned about once per catalog
   * (_nc_galaxy_shape_factor_spectral_warn_once()).
   */
  guint pop_order;
  gboolean pop_expand;
  NcGalaxyShapeFactorSpectralBox *pop_box;
  guint64 pop_box_key;
  guint64 pop_fit_key;
  NcGalaxyShapeFactorSpectralNodes *pop_nodes;
  gdouble *pop_T; /* the modes from the values at the nodes, n_nodes x n_nodes */
  gint pop_n_unconverged;
  gint pop_warned;
} NcGalaxyShapeFactorSpectralPrivate;

/* The scratch space of one evaluation */
typedef struct _NcGalaxyShapeFactorSpectralWork
{
  /*
   * pop-order: the panels combined from the modes of the galaxy whose modes_id is comb_id
   * with the weights of generation comb_w_gen
   */
  NcmSpectralPolarPanelsComb *comb;
  guint64 comb_id;
  guint64 comb_w_gen;
  GArray *r_edges;
  GArray *theta_edges;
  GArray *r_edges_all; /* edges of a peak before those outside a region are left out */
  GArray *theta_edges_all;
} NcGalaxyShapeFactorSpectralWork;

/*
 * Relative tolerance of the comparison of a loaded interpolant's population parameters and
 * settings with those of the fit: a few ulps, which a decimal round trip of the record can
 * change; a parameter difference this small changes ln P far below interp-abstol.
 */
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SOURCE_RELTOL (1.0e-15)

/*
 * The interpolant of ln P in g of one galaxy and what it was built for. It is rebuilt when
 * one of these values differs from the one of the evaluation, unless it was loaded from the
 * catalog (loaded): then src_gsfs and src_pop are the record that built it, which
 * data_prepare() compares with the fit, and a loaded interpolant is never rebuilt.
 */
typedef struct _NcGalaxyShapeFactorSpectralData
{
  gboolean valid;
  gboolean loaded;
  NcGalaxyShapeFactorSpectral *src_gsfs;
  NcGalaxyShapePop *src_pop;
  gdouble epsilon_obs_1;
  gdouble epsilon_obs_2;
  gdouble std_noise;
  guint64 pop_pkey;
  gdouble g_max;
  gdouble abstol;
  gdouble reltol;
  guint min_radial_panels;
  gdouble cos_e; /* direction of eps_obs */
  gdouble sin_e;
  gdouble g_c; /* center of the coordinates of the interpolant */

  /*
   * The radius in |g| the interpolant is built for: the reachable range of the galaxy that
   * data_prepare() found with restrict-range (memo holds the bound and what it was computed
   * for), or a negative value for interp-g-max. g_max above is the radius of the interpolant
   * held; shears beyond it are integrated. g_max_setting is the interp-g-max a loaded
   * interpolant was compared with.
   */
  gdouble g_reach;

  /*
   * The interpolant: the Chebyshev modes of its expansion over the population parameters,
   * panels of one layout packed in mode_spp and evaluated by combining them on the panel of
   * each run of points, and the box they span: that of the fit for a built expansion, that
   * of the record for a loaded one. A fixed population is the expansion with no free
   * parameter: one mode, box NULL, evaluated as it is (weight one). box_key is the
   * pop_box_key (built) or the pop_fit_key (loaded) they were made or compared for. w are
   * the weights of the modes at the parameters of w_pop at its pkey w_pkey, w_inside
   * whether those are in the box, w_gen a value unique to those weights among those of
   * modes_id (w_count of them so far).
   * expand_failed: the expansion of the galaxy did not converge, and it is integrated.
   */
  NcmSpectralPolarPanelsArray *mode_spp;
  const NcmSpectralPolarPanelsPack *mode_pack; /* the pack of mode_spp, which owns it */
  NcGalaxyShapeFactorSpectralBox *box;
  guint64 box_key;
  guint64 modes_id; /* unique to the modes in mode_spp, never 0 once they are set */
  gdouble *w;       /* n_nodes weights followed by the n_free K values of the polynomials */
  gconstpointer w_pop;
  guint64 w_pkey;
  guint64 w_gen;
  guint64 w_count; /* the weights made so far for the modes of modes_id */
  gboolean w_inside;
  gboolean expand_failed;
  gboolean box_is_fit; /* a loaded expansion over the box and order the fit expands over */
  NcGalaxyShapeFactorReachMemo reach_memo;
  gdouble g_max_setting;
} NcGalaxyShapeFactorSpectralData;

/*
 * The functions that a loop of the integrals calls with the convention as a constant: inlined
 * always, so that each copy of the loop holds a single map and vectorizes
 */
#ifdef __GNUC__
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_INLINE static inline __attribute__ ((always_inline))
#else
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_INLINE static inline
#endif

/* A new value of NcGalaxyShapeFactorSpectralData:modes_id, never 0 nor repeated */
static guint64
_nc_galaxy_shape_factor_spectral_id_new (void)
{
  static gsize last = 0;

  return (guint64) g_atomic_pointer_add (&last, 1) + 1;
}

static void
_nc_galaxy_shape_factor_spectral_box_clear_func (gpointer p)
{
  NcGalaxyShapeFactorSpectralBox *box = p;

  g_free (box->idx);
  g_free (box->lb);
  g_free (box->ub);
}

static void
_nc_galaxy_shape_factor_spectral_box_unref (NcGalaxyShapeFactorSpectralBox *box)
{
  g_atomic_rc_box_release_full (box, &_nc_galaxy_shape_factor_spectral_box_clear_func);
}

/* Whether a parameter of pop is free */
static gboolean
_nc_galaxy_shape_factor_spectral_pop_has_free (NcGalaxyShapePop *pop)
{
  NcmModel *model = NCM_MODEL (pop);
  guint pid;

  for (pid = 0; pid < ncm_model_len (model); pid++)
    if (ncm_model_param_get_ftype (model, pid) == NCM_PARAM_TYPE_FREE)
      return TRUE;

  return FALSE;
}

/*
 * The box of the expansion of order K over the free parameters of pop and their bounds,
 * NULL when no parameter is free. Aborts when a bound is not finite or the box has more
 * than NC_GALAXY_SHAPE_FACTOR_SPECTRAL_POP_NODES_MAX nodes.
 */
static NcGalaxyShapeFactorSpectralBox *
_nc_galaxy_shape_factor_spectral_box_new (NcGalaxyShapePop *pop, const guint K)
{
  NcmModel *model = NCM_MODEL (pop);
  const guint len = ncm_model_len (model);
  NcGalaxyShapeFactorSpectralBox *box;
  guint pid, i, n_free = 0;

  for (pid = 0; pid < len; pid++)
    if (ncm_model_param_get_ftype (model, pid) == NCM_PARAM_TYPE_FREE)
      n_free++;

  if (n_free == 0)
    return NULL;

  box          = g_atomic_rc_box_new0 (NcGalaxyShapeFactorSpectralBox);
  box->K       = K;
  box->n_free  = n_free;
  box->n_nodes = 1;
  box->idx     = g_new (guint, n_free);
  box->lb      = g_new (gdouble, n_free);
  box->ub      = g_new (gdouble, n_free);

  for (i = 0; i < n_free; i++)
  {
    box->n_nodes *= K;

    if (box->n_nodes > NC_GALAXY_SHAPE_FACTOR_SPECTRAL_POP_NODES_MAX)
      g_error ("nc_galaxy_shape_factor_spectral: pop-order %u with %u free population parameters makes more than %u nodes.",
               K, n_free, NC_GALAXY_SHAPE_FACTOR_SPECTRAL_POP_NODES_MAX);
  }

  for (pid = 0, i = 0; pid < len; pid++)
  {
    if (ncm_model_param_get_ftype (model, pid) == NCM_PARAM_TYPE_FREE)
    {
      box->idx[i] = pid;
      box->lb[i]  = ncm_model_param_get_lower_bound (model, pid);
      box->ub[i]  = ncm_model_param_get_upper_bound (model, pid);

      if (!(box->lb[i] < box->ub[i]) || !gsl_finite (box->ub[i] - box->lb[i]))
        g_error ("nc_galaxy_shape_factor_spectral: pop-order needs finite bounds for the free parameter %s of %s.",
                 ncm_model_param_name (model, pid), G_OBJECT_TYPE_NAME (pop));

      i++;
    }
  }

  return box;
}

/*
 * Tolerance of the test of a point on the box of an expansion: a parameter set to a bound
 * of the box can differ from it by a few ulps after a decimal round trip.
 */
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_BOX_TOL (1.0e-12)

static void _nc_galaxy_shape_factor_spectral_cheb_values (const guint K, const gdouble x, gdouble *T);

/*
 * The weights of the modes of the galaxy of ldata at the current parameters of pop, made
 * again only when pop or its parameters change. Returns whether the parameters are in the
 * box of the modes; outside of it the galaxy is integrated.
 */
static gboolean
_nc_galaxy_shape_factor_spectral_weights (NcGalaxyShapeFactorSpectralData *ldata, NcGalaxyShapePop *pop)
{
  const guint64 pkey = ncm_model_state_get_pkey (NCM_MODEL (pop));

  if ((ldata->w_gen == 0) || (ldata->w_pop != (gconstpointer) pop) || (ldata->w_pkey != pkey))
  {
    const NcGalaxyShapeFactorSpectralBox *box = ldata->box;
    const guint K                             = box->K;
    gdouble *L                                = ldata->w + box->n_nodes;
    guint i, q;

    ldata->w_inside = TRUE;

    for (i = 0; i < box->n_free; i++)
    {
      const gdouble theta = ncm_model_param_get (NCM_MODEL (pop), box->idx[i]);
      const gdouble x     = (2.0 * theta - box->lb[i] - box->ub[i]) / (box->ub[i] - box->lb[i]);

      if (!(fabs (x) <= 1.0 + NC_GALAXY_SHAPE_FACTOR_SPECTRAL_BOX_TOL))
        ldata->w_inside = FALSE;

      _nc_galaxy_shape_factor_spectral_cheb_values (K, CLAMP (x, -1.0, 1.0), L + (gsize) i * K);
    }

    for (q = 0; q < box->n_nodes; q++)
    {
      guint r   = q;
      gdouble w = 1.0;

      for (i = 0; i < box->n_free; i++)
      {
        w *= L[(gsize) i * K + r % K];
        r /= K;
      }

      ldata->w[q] = w;
    }

    /*
     * Unique with modes_id, which the caches of combined panels are keyed by too: a count
     * of the galaxy instead of a global one, which every thread would write for every
     * galaxy. It restarts only with new modes, and so with a new modes_id.
     */
    ldata->w_pop  = pop;
    ldata->w_pkey = pkey;
    ldata->w_gen  = ++ldata->w_count;
  }

  return ldata->w_inside;
}

/* The weight of the single mode of a fixed population */
static const gdouble _nc_galaxy_shape_factor_spectral_one = 1.0;

/*
 * Whether the interpolant of the galaxy of ldata can be used at the current parameters of
 * pop: it converged, and an expansion has pop inside its box
 * (_nc_galaxy_shape_factor_spectral_weights()). Sets *w to the weights of its modes.
 */
static gboolean
_nc_galaxy_shape_factor_spectral_usable (NcGalaxyShapeFactorSpectralData *ldata, NcGalaxyShapePop *pop, const gdouble **w)
{
  if (ldata->expand_failed)
    return FALSE;

  if (ldata->box == NULL)
  {
    *w = &_nc_galaxy_shape_factor_spectral_one;

    return TRUE;
  }

  *w = ldata->w;

  return _nc_galaxy_shape_factor_spectral_weights (ldata, pop);
}

/* Sets the modes of ldata, of the expansion over box (NULL for a fixed population), and the space of their weights */
static void
_nc_galaxy_shape_factor_spectral_ldata_set_modes (NcGalaxyShapeFactorSpectralData *ldata, NcmSpectralPolarPanelsArray *modes,
                                                  NcGalaxyShapeFactorSpectralBox *box)
{
  if (modes != ldata->mode_spp)
  {
    ncm_spectral_polar_panels_array_clear (&ldata->mode_spp);
    ldata->mode_spp  = (modes != NULL) ? ncm_spectral_polar_panels_array_ref (modes) : NULL;
    ldata->mode_pack = (modes != NULL) ? ncm_spectral_polar_panels_array_peek_pack (modes) : NULL;
    ldata->modes_id  = (modes != NULL) ? _nc_galaxy_shape_factor_spectral_id_new () : 0;
    ldata->w_count   = 0;
  }

  if (box != ldata->box)
  {
    g_clear_pointer (&ldata->box, _nc_galaxy_shape_factor_spectral_box_unref);
    ldata->box = (box != NULL) ? g_atomic_rc_box_acquire (box) : NULL;
    ldata->w   = (box != NULL) ? g_renew (gdouble, ldata->w, box->n_nodes + (gsize) box->n_free * box->K) : ldata->w;
  }

  ldata->w_gen = 0;
}

/* The work's cache of combined panels for the modes of ldata and the current weights */
static NcmSpectralPolarPanelsComb *
_nc_galaxy_shape_factor_spectral_work_comb (NcGalaxyShapeFactorSpectralWork *work, const guint64 modes_id, const guint64 w_gen)
{
  if ((work->comb_id != modes_id) || (work->comb_w_gen != w_gen))
  {
    ncm_spectral_polar_panels_comb_reset (work->comb);
    work->comb_id    = modes_id;
    work->comb_w_gen = w_gen;
  }

  return work->comb;
}

/*
 * A set of populations evaluated together (_nc_galaxy_shape_factor_spectral_eval_ln_nodes()):
 * the layout population, whose adaptive integral places the panels, and the n_w populations
 * integrated on them, with their prepared data, their exponents and the quadrature rules of
 * those exponents. Read-only once made, so that threads share it.
 */
struct _NcGalaxyShapeFactorSpectralNodes
{
  NcGalaxyShapePop *layout_pop;
  NcGalaxyShapePopData *layout_data;
  gdouble layout_a;
  gdouble layout_b;
  gdouble layout_dens_norm;
  gdouble layout_dens_pow;
  guint n_w;
  NcGalaxyShapePop **pops;
  NcGalaxyShapePopData **pop_data;
  gdouble *a;
  gdouble *b;
  gdouble *dens_norm;
  gdouble *dens_pow;
  NcmSpectralPolarMulti *multi;
};


/* What the integrand needs, in the frame where the shear is the real number t */
struct _NcGalaxyShapeFactorSpectralKernel
{
  NcGalaxyShapePop *pop;
  NcGalaxyShapePopData *pop_data;
  gdouble t;
  gdouble e_1;
  gdouble e_2;
  gdouble c2;        /* 1 / (2 std_noise^2) */
  gdouble norm;      /* 1 / (2 pi std_noise^2) / (2 pi) */
  gdouble shift_c2;  /* c2 times the offset taken out of the exponent */
  gdouble dens_norm; /* TRACE: 2^(alpha_o + 1), alpha_o and alpha_1 the exponents of the population */
  gdouble dens_pow;  /* TRACE: alpha_o + alpha_1 + 2 */
  gboolean trace;    /* the convention: TRUE for TRACE (NcGalaxyShapeFactorSpectralPrivate:axis_ratio_var) */
};

/* The peak of the integrand on the ray through the preimage of eps_obs */
typedef struct _NcGalaxyShapeFactorSpectralSeed
{
  gdouble r_star;  /* radius of the maximum of h on the scan */
  gdouble theta_0; /* angle of the ray */
  gdouble J_1;     /* radial stretch factor of the inverse map at eps_obs */
  gdouble J_2;     /* tangential stretch factor */
  gdouble shift;   /* smallest squared distance of eps_obs to the sheared ellipticity on the scan */
} NcGalaxyShapeFactorSpectralSeed;

enum
{
  PROP_0,
  PROP_RELTOL,
  PROP_USE_INTERP,
  PROP_INTERP_G_MAX,
  PROP_INTERP_ABSTOL,
  PROP_INTERP_MIN_RADIAL_PANELS,
  PROP_LINE_ABSTOL,
  PROP_RESTRICT_RANGE,
  PROP_POP_ORDER,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyShapeFactorSpectral, nc_galaxy_shape_factor_spectral, NC_TYPE_GALAXY_SHAPE_FACTOR)

static gpointer _nc_galaxy_shape_factor_spectral_work_new (gpointer userdata);
static void _nc_galaxy_shape_factor_spectral_work_free (gpointer p);
static gpointer _nc_galaxy_shape_factor_spectral_line_new (gpointer userdata);
static void _nc_galaxy_shape_factor_spectral_line_free (gpointer p);

static void
nc_galaxy_shape_factor_spectral_init (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  self->map_inv        = NULL;
  self->axis_ratio_var = FALSE;
  self->polar          = ncm_spectral_polar_new ();
  self->polar_interp   = ncm_spectral_polar_new ();
  self->work_pool      = ncm_memory_pool_new (_nc_galaxy_shape_factor_spectral_work_new, NULL,
                                              _nc_galaxy_shape_factor_spectral_work_free);
  self->line_pool = ncm_memory_pool_new (_nc_galaxy_shape_factor_spectral_line_new, NULL,
                                         _nc_galaxy_shape_factor_spectral_line_free);
  self->reltol        = 0.0;
  self->line_abstol   = 0.0;
  self->use_interp    = FALSE;
  self->interp_g_max  = 0.0;
  self->interp_abstol = 0.0;

  self->interp_min_radial_panels = 0;

  self->restrict_range = FALSE;
  self->range_on       = FALSE;
  self->n_short        = 0;
  _nc_galaxy_shape_factor_reach_init (&self->reach, TRUE);

  self->pop_order         = 0;
  self->pop_expand        = FALSE;
  self->pop_box           = NULL;
  self->pop_box_key       = 0;
  self->pop_fit_key       = 0;
  self->pop_nodes         = NULL;
  self->pop_T             = NULL;
  self->pop_n_unconverged = 0;
  self->pop_warned        = 0;
}

static void
_nc_galaxy_shape_factor_spectral_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactorSpectral *gsfs = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (object);

  switch (prop_id)
  {
    case PROP_RELTOL:
      nc_galaxy_shape_factor_spectral_set_reltol (gsfs, g_value_get_double (value));
      break;
    case PROP_USE_INTERP:
      nc_galaxy_shape_factor_spectral_set_use_interp (gsfs, g_value_get_boolean (value));
      break;
    case PROP_INTERP_G_MAX:
      nc_galaxy_shape_factor_spectral_set_interp_g_max (gsfs, g_value_get_double (value));
      break;
    case PROP_INTERP_ABSTOL:
      nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs, g_value_get_double (value));
      break;
    case PROP_INTERP_MIN_RADIAL_PANELS:
      nc_galaxy_shape_factor_spectral_set_interp_min_radial_panels (gsfs, g_value_get_uint (value));
      break;
    case PROP_LINE_ABSTOL:
      nc_galaxy_shape_factor_spectral_set_line_abstol (gsfs, g_value_get_double (value));
      break;
    case PROP_RESTRICT_RANGE:
      nc_galaxy_shape_factor_spectral_set_restrict_range (gsfs, g_value_get_boolean (value));
      break;
    case PROP_POP_ORDER:
      nc_galaxy_shape_factor_spectral_set_pop_order (gsfs, g_value_get_uint (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_factor_spectral_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactorSpectral *gsfs = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (object);

  switch (prop_id)
  {
    case PROP_RELTOL:
      g_value_set_double (value, nc_galaxy_shape_factor_spectral_get_reltol (gsfs));
      break;
    case PROP_USE_INTERP:
      g_value_set_boolean (value, nc_galaxy_shape_factor_spectral_get_use_interp (gsfs));
      break;
    case PROP_INTERP_G_MAX:
      g_value_set_double (value, nc_galaxy_shape_factor_spectral_get_interp_g_max (gsfs));
      break;
    case PROP_INTERP_ABSTOL:
      g_value_set_double (value, nc_galaxy_shape_factor_spectral_get_interp_abstol (gsfs));
      break;
    case PROP_INTERP_MIN_RADIAL_PANELS:
      g_value_set_uint (value, nc_galaxy_shape_factor_spectral_get_interp_min_radial_panels (gsfs));
      break;
    case PROP_LINE_ABSTOL:
      g_value_set_double (value, nc_galaxy_shape_factor_spectral_get_line_abstol (gsfs));
      break;
    case PROP_RESTRICT_RANGE:
      g_value_set_boolean (value, nc_galaxy_shape_factor_spectral_get_restrict_range (gsfs));
      break;
    case PROP_POP_ORDER:
      g_value_set_uint (value, nc_galaxy_shape_factor_spectral_get_pop_order (gsfs));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static complex double _nc_galaxy_shape_factor_spectral_map_inv_trace (complex double g, complex double e_obs);

static void _nc_galaxy_shape_factor_spectral_map_trace (gdouble t, gdouble x, gdouble y, gdouble * restrict f_x, gdouble * restrict f_y);
static void _nc_galaxy_shape_factor_spectral_kernel (gpointer user_data, gdouble r, guint n, const gdouble *theta, const gdouble *cos_theta, const gdouble *sin_theta, gdouble *fval);
static void _nc_galaxy_shape_factor_spectral_radial (gpointer user_data, guint n, const gdouble *r, gdouble *rho);
static void _nc_galaxy_shape_factor_spectral_density (const NcGalaxyShapeFactorSpectralKernel *k, guint n, const gdouble *r, gdouble *out);
static void _nc_galaxy_shape_factor_spectral_dist2 (const NcGalaxyShapeFactorSpectralKernel *k, guint n, const gdouble *r, gdouble cos_0, gdouble sin_0, gdouble *d);

static void
_nc_galaxy_shape_factor_spectral_constructed (GObject *object)
{
  /* Chain up: start */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_spectral_parent_class)->constructed (object);
  {
    NcGalaxyShapeFactorSpectral *gsfs               = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (object);
    NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);
    const NcGalaxyWLObsEllipConv ellip_conv         = nc_galaxy_shape_factor_get_ellip_conv (NC_GALAXY_SHAPE_FACTOR (gsfs));

    switch (ellip_conv)
    {
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
        self->map_inv        = &_nc_galaxy_shape_factor_spectral_map_inv_trace;
        self->axis_ratio_var = TRUE;
        break;
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
        self->map_inv        = &nc_wl_ellipticity_apply_shear_inv_trace_det;
        self->axis_ratio_var = FALSE;
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
    }
  }
}

static void
_nc_galaxy_shape_factor_spectral_dispose (GObject *object)
{
  NcGalaxyShapeFactorSpectral *gsfs               = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (object);
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  ncm_spectral_polar_clear (&self->polar);
  ncm_spectral_polar_clear (&self->polar_interp);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_spectral_parent_class)->dispose (object);
}

static void
_nc_galaxy_shape_factor_spectral_finalize (GObject *object)
{
  NcGalaxyShapeFactorSpectral *gsfs               = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (object);
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  ncm_memory_pool_free (self->work_pool, TRUE);
  ncm_memory_pool_free (self->line_pool, TRUE);
  _nc_galaxy_shape_factor_reach_clear (&self->reach);
  g_clear_pointer (&self->pop_nodes, _nc_galaxy_shape_factor_spectral_nodes_free);
  g_clear_pointer (&self->pop_box, _nc_galaxy_shape_factor_spectral_box_unref);
  g_free (self->pop_T);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_spectral_parent_class)->finalize (object);
}

static void _nc_galaxy_shape_factor_spectral_data_init (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data);
static void _nc_galaxy_shape_factor_spectral_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset);
static void _nc_galaxy_shape_factor_spectral_data_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data, const gdouble z_max);
static gdouble _nc_galaxy_shape_factor_spectral_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2);
static gdouble _nc_galaxy_shape_factor_spectral_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2);
static void _nc_galaxy_shape_factor_spectral_data_prefetch (NcGalaxyShapeFactor *gsf, NcGalaxyShapeFactorData *data, const guint stage);
static void _nc_galaxy_shape_factor_spectral_eval_nodes (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const guint n, const gdouble *g_1, const gdouble *g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gboolean want_log, gdouble *out);
static gboolean _nc_galaxy_shape_factor_spectral_line_table (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble t_min, const gdouble t_max, const gdouble u_1, const gdouble u_2, const gdouble b_1, const gdouble b_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const guint max_points, GArray **coeffs);

static void
nc_galaxy_shape_factor_spectral_class_init (NcGalaxyShapeFactorSpectralClass *klass)
{
  NcGalaxyShapeFactorClass *gsf_class = NC_GALAXY_SHAPE_FACTOR_CLASS (klass);
  GObjectClass *object_class          = G_OBJECT_CLASS (klass);

  object_class->constructed  = &_nc_galaxy_shape_factor_spectral_constructed;
  object_class->set_property = &_nc_galaxy_shape_factor_spectral_set_property;
  object_class->get_property = &_nc_galaxy_shape_factor_spectral_get_property;
  object_class->dispose      = &_nc_galaxy_shape_factor_spectral_dispose;
  object_class->finalize     = &_nc_galaxy_shape_factor_spectral_finalize;

  /**
   * NcGalaxyShapeFactorSpectral:reltol:
   *
   * The relative tolerance of the integral, passed to ncm_spectral_polar_integrate().
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance of the integral",
                                                        GSL_DBL_EPSILON, 1.0, 1.0e-5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyShapeFactorSpectral:use-interp:
   *
   * Whether nc_galaxy_shape_factor_eval_marginal() and
   * nc_galaxy_shape_factor_eval_ln_marginal() take $\ln P$ from the interpolant in $g$ of
   * the galaxy for $|g| \le$ #NcGalaxyShapeFactorSpectral:interp-g-max. The interpolant is
   * built at the first such evaluation of a galaxy and again when its observed
   * ellipticity, its noise, the parameters of the population or a tolerance changes. A
   * build computes of the order of $10^3$ integrals, so it pays when the population
   * parameters are fixed.
   */
  g_object_class_install_property (object_class,
                                   PROP_USE_INTERP,
                                   g_param_spec_boolean ("use-interp",
                                                         NULL,
                                                         "Whether ln P is interpolated in g for each galaxy",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyShapeFactorSpectral:interp-g-max:
   *
   * The largest $|g|$ computed from the interpolant. It is below $1/3$: the region of the
   * interpolant reaches $|g| < 3$ times this value and stays inside $|g| < 1$.
   */
  g_object_class_install_property (object_class,
                                   PROP_INTERP_G_MAX,
                                   g_param_spec_double ("interp-g-max",
                                                        NULL,
                                                        "Largest |g| computed from the interpolant",
                                                        1.0e-3, 0.33, 0.3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyShapeFactorSpectral:interp-abstol:
   *
   * The absolute tolerance on $\ln P$ of the interpolant, passed to
   * ncm_spectral_polar_build(): a panel is accepted when $\ln P$ changes by at most this
   * value between its two grids. Each panel of the interpolant is then trimmed
   * (ncm_spectral_polar_panels_trim()) to the coefficients that change $\ln P$ by more than
   * $10^{-4}$ times this value, in memory and in a catalog prepared with them.
   */
  g_object_class_install_property (object_class,
                                   PROP_INTERP_ABSTOL,
                                   g_param_spec_double ("interp-abstol",
                                                        NULL,
                                                        "Absolute tolerance on ln P of the interpolant",
                                                        GSL_DBL_EPSILON, 1.0, 1.0e-5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyShapeFactorSpectral:interp-min-radial-panels:
   *
   * The #NcmSpectralPolar:min-radial-panels of the build of the interpolant: the radial
   * range of the interpolant is split into at least this many equal intervals. More
   * intervals give panels of lower degree, so a faster evaluation of the interpolant, at the
   * cost of more integrals in the build.
   */
  g_object_class_install_property (object_class,
                                   PROP_INTERP_MIN_RADIAL_PANELS,
                                   g_param_spec_uint ("interp-min-radial-panels",
                                                      NULL,
                                                      "Minimum number of radial intervals of the interpolant",
                                                      1, G_MAXUINT, 2,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyShapeFactorSpectral:line-abstol:
   *
   * The absolute tolerance on $\ln P$ of the table along the line of the shears of the
   * nodes of a galaxy (see nc_galaxy_shape_factor_eval_line()). The nodes of a galaxy whose
   * shears are integrated (without #NcGalaxyShapeFactorSpectral:use-interp; with it the
   * interpolant is used, except for a galaxy that cannot use it anywhere on the line, see
   * nc_galaxy_shape_factor_spectral_get_fallback_stats()) take $\ln P$ from a Chebyshev series
   * in the position on that line, built from integrals on nested Chebyshev-Lobatto grids
   * and accepted when its coefficients change by at most this value in $\ell_1$ norm
   * between two grids. The table is built only when it needs at most half as many integrals
   * as there are nodes and every shear of the line has $|g| \le 0.9$; otherwise each node
   * is integrated, as when the table does not reach the tolerance. Each integral is
   * accurate to about
   * #NcGalaxyShapeFactorSpectral:reltol relative, which is its accuracy in $\ln P$: a
   * negative value, the default, takes reltol / 10, so that the table adds an error below
   * that of the integrals. Zero integrates every node.
   */
  g_object_class_install_property (object_class,
                                   PROP_LINE_ABSTOL,
                                   g_param_spec_double ("line-abstol",
                                                        NULL,
                                                        "Absolute tolerance on ln P of the table along the line of the shears (negative: reltol / 10)",
                                                        -1.0, 1.0, -1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyShapeFactorSpectral:restrict-range:
   *
   * With #NcGalaxyShapeFactorSpectral:use-interp, whether nc_galaxy_shape_factor_data_prepare()
   * builds the interpolant of each galaxy only over the shears it can reach in the box of
   * the fit, instead of the whole disk $|g| \le$ #NcGalaxyShapeFactorSpectral:interp-g-max.
   *
   * The box is that of the free parameters: the halo center over the bounds of ra and dec,
   * the mass up to the upper bound of log10MDelta and a free cDelta over its bounds. The
   * bound of a galaxy is the calibrated reduced shear $|1 + m| |g_t| + |c|$ at the upper
   * bound of the mass, the top of the redshift support of the galaxy and the smallest
   * projected radius the box of the center allows; a free cDelta is sampled at 9 log-spaced
   * points of its box and the largest value taken times 1.05, since the shear is not
   * monotonic in the concentration. The interpolant of a galaxy covers
   * $|g| \le$ min (bound, interp-g-max), at least 0.02, or interp-g-max where the convergence
   * reaches one at the corner of the box. The bound depends on the box only, not on the
   * current point, so it is computed once per galaxy and the interpolants are not rebuilt as
   * the fit moves.
   *
   * The restriction applies only when the box covers every free parameter: log10MDelta, ra,
   * dec, cDelta and the parameters of the population may be free. With any other free
   * parameter (the cosmology, the redshift of the halo) every galaxy gets interp-g-max.
   *
   * A shear beyond the radius of the interpolant of its galaxy is integrated, as one beyond
   * interp-g-max always is, so the values are those of the integrals there: a box too small
   * costs time, not accuracy. At that radius $\ln P$ changes from the interpolant to the
   * integral, which differ by up to #NcGalaxyShapeFactorSpectral:interp-abstol. The
   * interpolants of a catalog cover the box it was prepared for; a fit with a larger box
   * integrates the shears beyond it, with a warning
   * (nc_galaxy_shape_factor_spectral_get_n_loaded_short()).
   *
   * #NcGalaxyShapeFactorMomentsTilt:restrict-range takes the same bound
   * at a fixed concentration only, since its tables have no integral to fall back on.
   */
  g_object_class_install_property (object_class,
                                   PROP_RESTRICT_RANGE,
                                   g_param_spec_boolean ("restrict-range",
                                                         NULL,
                                                         "Whether to build each interpolant only over the reachable shears",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyShapeFactorSpectral:pop-order:
   *
   * With #NcGalaxyShapeFactorSpectral:use-interp, the number $K$ of nodes per free
   * parameter of the population over which the interpolant of each galaxy is expanded; 0 or
   * 1 for none. With $d$ free parameters, each between its lower and upper bound,
   * nc_galaxy_shape_factor_prepare() makes the $K^d$ populations at the tensor
   * Chebyshev-Lobatto nodes of that box, and nc_galaxy_shape_factor_data_prepare() builds
   * for each galaxy the interpolants in $g$ of their $\ln P$ on one panel layout, that of
   * the population at the center of the box. They are kept as the Chebyshev modes of the
   * interpolation over the box, $\ln P = \sum_m T_m (\theta) L_m (g)$ with $T_m$ the
   * products of Chebyshev polynomials of the parameters mapped to $[-1, 1]$, each mode cut
   * to the coefficients it needs: $\ln P$ changes by at most $10^{-4}$
   * #NcGalaxyShapeFactorSpectral:interp-abstol (ncm_spectral_polar_panels_trim_modes()), the
   * tolerance of the cut of the interpolant of one population. The error of the expansion
   * over the box is estimated from the modes of the highest degree in some parameter and,
   * for $K \ge 4$, from their prediction by the two degrees below (which does not vanish
   * for a function even or odd in a parameter); a galaxy where it exceeds
   * #NcGalaxyShapeFactorSpectral:interp-abstol is integrated instead, and counted
   * (nc_galaxy_shape_factor_spectral_get_pop_n_unconverged()), with a warning at the first,
   * so that a too low order or a too wide box costs time, not accuracy.
   * The modes are combined with the weights $T_m (\theta)$ of the current parameters, once
   * per change of the parameters, and the galaxy is evaluated as with a fixed population. A
   * population outside its box is integrated. Any #NcGalaxyShapePop works, with any number of
   * free parameters, for both ellipticity conventions. The populations of the nodes are
   * prepared once, for every galaxy, so the population must not depend on the galaxy (its
   * #NcGalaxyShapePopData holds no data of the galaxy).
   *
   * The modes are written to prepared catalogs (#NcGalaxyWLObs:shape-interp), with the
   * box in the record of the catalog. A fit that reads them evaluates the expansion over that
   * box, whatever its own pop-order: it may free any of the parameters of the box, over any
   * bounds, or fix them at any value, and must fix the other parameters at the values of the
   * record. The bounds of the fit are its prior and are never changed: a population outside
   * the box is integrated, at the cost of the integrals, and a warning is given once.
   */
  g_object_class_install_property (object_class,
                                   PROP_POP_ORDER,
                                   g_param_spec_uint ("pop-order",
                                                      NULL,
                                                      "Nodes per free population parameter of the expansion of the interpolants",
                                                      0, 64, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  gsf_class->data_init        = &_nc_galaxy_shape_factor_spectral_data_init;
  gsf_class->prepare          = &_nc_galaxy_shape_factor_spectral_prepare;
  gsf_class->data_prepare     = &_nc_galaxy_shape_factor_spectral_data_prepare;
  gsf_class->eval_marginal    = &_nc_galaxy_shape_factor_spectral_eval_marginal;
  gsf_class->eval_ln_marginal = &_nc_galaxy_shape_factor_spectral_eval_ln_marginal;
  gsf_class->eval_nodes       = &_nc_galaxy_shape_factor_spectral_eval_nodes;
  gsf_class->data_prefetch    = &_nc_galaxy_shape_factor_spectral_data_prefetch;
  gsf_class->line_table       = &_nc_galaxy_shape_factor_spectral_line_table;
}

static gpointer
_nc_galaxy_shape_factor_spectral_work_new (gpointer userdata)
{
  NcGalaxyShapeFactorSpectralWork *work = g_new0 (NcGalaxyShapeFactorSpectralWork, 1);

  work->comb            = ncm_spectral_polar_panels_comb_new ();
  work->comb_id         = 0;
  work->comb_w_gen      = 0;
  work->r_edges         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->theta_edges     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->r_edges_all     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->theta_edges_all = g_array_new (FALSE, FALSE, sizeof (gdouble));

  return work;
}

static gpointer
_nc_galaxy_shape_factor_spectral_line_new (gpointer userdata)
{
  return ncm_spectral_new_with_max_level (NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_LEVEL_MAX);
}

static void
_nc_galaxy_shape_factor_spectral_line_free (gpointer p)
{
  ncm_spectral_free (p);
}

static void
_nc_galaxy_shape_factor_spectral_work_free (gpointer p)
{
  NcGalaxyShapeFactorSpectralWork *work = p;

  ncm_spectral_polar_panels_comb_free (work->comb);
  g_array_unref (work->r_edges);
  g_array_unref (work->theta_edges);
  g_array_unref (work->r_edges_all);
  g_array_unref (work->theta_edges_all);

  g_free (work);
}

/*
 * Takes the interpolant of row i of obs, when the catalog has one for its current rows, as a
 * loaded interpolant (shared with obs, never rebuilt in place). Otherwise the galaxy has no
 * valid interpolant.
 */
static void
_nc_galaxy_shape_factor_spectral_ldata_read_row (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxyShapeFactorSpectralData *ldata = data->ldata;
  NcmSpectralPolarPanelsArray *modes     = nc_galaxy_wl_obs_peek_shape_interp (obs, i);

  if (ldata->loaded)
  {
    _nc_galaxy_shape_factor_spectral_ldata_set_modes (ldata, NULL, NULL);
    g_clear_object (&ldata->src_gsfs);
    nc_galaxy_shape_pop_clear (&ldata->src_pop);
    ldata->loaded = FALSE;
  }

  ldata->valid         = FALSE;
  ldata->expand_failed = FALSE;

  if (modes != NULL)
  {
    NcmObjDictStr *source = nc_galaxy_wl_obs_peek_shape_interp_source (obs);
    GObject *src_gsfs     = ncm_obj_dict_str_peek (source, "shape-factor");
    GObject *src_pop      = ncm_obj_dict_str_peek (source, "shape-pop");

    if (!NC_IS_GALAXY_SHAPE_FACTOR_SPECTRAL (src_gsfs) || !NC_IS_GALAXY_SHAPE_POP (src_pop))
      g_error ("nc_galaxy_shape_factor_spectral: the catalog has interpolants, but its record "
               "lacks the NcGalaxyShapeFactorSpectral `shape-factor' or the NcGalaxyShapePop `shape-pop'.");

    if (ncm_spectral_polar_panels_array_len (modes) == 0)
      g_error ("nc_galaxy_shape_factor_spectral: row %u of the catalog has an interpolant with no modes.", i);

    /* The modes, whose box data_prepare() reads from the record */
    _nc_galaxy_shape_factor_spectral_ldata_set_modes (ldata, modes, NULL);

    ldata->src_gsfs = g_object_ref (NC_GALAXY_SHAPE_FACTOR_SPECTRAL (src_gsfs));
    ldata->src_pop  = nc_galaxy_shape_pop_ref (NC_GALAXY_SHAPE_POP (src_pop));
    ldata->loaded   = TRUE;
  }
}

/*
 * Writes the interpolant of the galaxy, the modes of its expansion over the population
 * parameters, to row i of obs when it is valid for the observed ellipticity and noise of the
 * galaxy. A galaxy whose expansion did not converge has none, and its row none.
 * nc_galaxy_shape_factor_spectral_set_obs_source() makes the written interpolants current.
 */
static void
_nc_galaxy_shape_factor_spectral_ldata_write_row (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxyShapeFactorSpectralData *ldata = data->ldata;

  if (ldata->valid && (ldata->epsilon_obs_1 == data->epsilon_obs_1) && (ldata->epsilon_obs_2 == data->epsilon_obs_2) &&
      (ldata->std_noise == data->std_noise))
  {
    /*
     * The record written with the rows is that of the fit
     * (nc_galaxy_shape_factor_spectral_set_obs_source()), so a loaded expansion over another
     * box or order would be read over the wrong box.
     */
    if ((ldata->box != NULL) && ldata->loaded && !ldata->box_is_fit)
      g_error ("nc_galaxy_shape_factor_spectral: row %u holds an expansion loaded over a box other than the one the fit "
               "expands over (the free population parameters, their bounds, the fixed values and pop-order); "
               "prepare the catalog again from rows without interpolants.", i);
    else if (ldata->mode_spp != NULL)
      nc_galaxy_wl_obs_set_shape_interp (obs, i, ldata->mode_spp);
  }
}

static void
_nc_galaxy_shape_factor_spectral_ldata_required_columns (NcGalaxyShapeFactorData *data, GList **columns)
{
}

static void _nc_galaxy_shape_factor_spectral_ldata_free (gpointer p);

/*
 * The state per galaxy is the interpolant of ln P in g, empty until
 * NcGalaxyShapeFactorSpectral:use-interp builds it or a row of a catalog that carries one
 * is read. Nothing of it is a column of the catalog.
 */
static void
_nc_galaxy_shape_factor_spectral_data_init (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data)
{
  NcGalaxyShapeFactorSpectralData *ldata = g_new0 (NcGalaxyShapeFactorSpectralData, 1);

  ldata->valid    = FALSE;
  ldata->loaded   = FALSE;
  ldata->src_gsfs = NULL;
  ldata->g_reach  = -1.0;

  ldata->mode_spp       = NULL;
  ldata->mode_pack      = NULL;
  ldata->box            = NULL;
  ldata->box_key        = 0;
  ldata->modes_id       = 0;
  ldata->w              = NULL;
  ldata->w_pop          = NULL;
  ldata->w_pkey         = 0;
  ldata->w_gen          = 0;
  ldata->w_count        = 0;
  ldata->w_inside       = FALSE;
  ldata->expand_failed  = FALSE;
  ldata->box_is_fit     = FALSE;
  ldata->reach_memo.gen = 0;
  ldata->g_max_setting  = 0.0;
  ldata->src_pop        = NULL;

  data->ldata                  = ldata;
  data->ldata_destroy          = &_nc_galaxy_shape_factor_spectral_ldata_free;
  data->ldata_read_row         = &_nc_galaxy_shape_factor_spectral_ldata_read_row;
  data->ldata_write_row        = &_nc_galaxy_shape_factor_spectral_ldata_write_row;
  data->ldata_required_columns = &_nc_galaxy_shape_factor_spectral_ldata_required_columns;
}

static void
_nc_galaxy_shape_factor_spectral_ldata_free (gpointer p)
{
  NcGalaxyShapeFactorSpectralData *ldata = p;

  _nc_galaxy_shape_factor_spectral_ldata_set_modes (ldata, NULL, NULL);
  g_clear_object (&ldata->src_gsfs);
  nc_galaxy_shape_pop_clear (&ldata->src_pop);
  g_free (ldata->w);
  g_free (ldata);
}

static NcGalaxyShapeFactorSpectralNodes *_nc_galaxy_shape_factor_spectral_nodes_new (NcGalaxyShapeFactorSpectralPrivate * const self, NcGalaxyShapePop *layout_pop, const guint n_w, NcGalaxyShapePop **pops);

/* T_0 (x), ..., T_{K-1} (x) */
static void
_nc_galaxy_shape_factor_spectral_cheb_values (const guint K, const gdouble x, gdouble *T)
{
  guint k;

  T[0] = 1.0;

  if (K > 1)
    T[1] = x;

  for (k = 2; k < K; k++)
    T[k] = 2.0 * x * T[k - 1] - T[k - 2];
}

/*
 * The Chebyshev coefficients of the interpolant of degree K - 1 from its values at the K
 * Chebyshev-Lobatto nodes cos (pi j / (K - 1)): c_a = sum_j V[a K + j] f_j, with
 * V[a K + j] = 2 / (K - 1) h_a h_j T_a (x_j), h = 1/2 at the ends.
 */
static void
_nc_galaxy_shape_factor_spectral_cheb_transform (const guint K, gdouble *V)
{
  guint a, j;

  for (a = 0; a < K; a++)
  {
    const gdouble h_a = ((a == 0) || (a == K - 1)) ? 0.5 : 1.0;

    for (j = 0; j < K; j++)
    {
      const gdouble h_j = ((j == 0) || (j == K - 1)) ? 0.5 : 1.0;

      V[a * K + j] = 2.0 / (K - 1.0) * h_a * h_j * cos (M_PI * (gdouble) ((a * j) % (2 * (K - 1))) / (K - 1.0));
    }
  }
}

/* A 64-bit FNV-1a step over n bytes */
static guint64
_nc_galaxy_shape_factor_spectral_fnv (guint64 h, gconstpointer p, gsize n)
{
  const guint8 *c = p;
  gsize i;

  for (i = 0; i < n; i++)
  {
    h ^= c[i];
    h *= 1099511628211ULL;
  }

  return h;
}

/*
 * The key of the box of the expansion of order K over the free parameters of pop: its type,
 * K, which parameters are free, the bounds of the free ones and the values of the others.
 */
static guint64
_nc_galaxy_shape_factor_spectral_pop_box_key (NcGalaxyShapePop *pop, const guint K)
{
  NcmModel *model = NCM_MODEL (pop);
  const guint len = ncm_model_len (model);
  guint64 key     = 14695981039346656037ULL;
  guint pid;

  key = _nc_galaxy_shape_factor_spectral_fnv (key, G_OBJECT_TYPE_NAME (pop), strlen (G_OBJECT_TYPE_NAME (pop)));
  key = _nc_galaxy_shape_factor_spectral_fnv (key, &K, sizeof (K));

  for (pid = 0; pid < len; pid++)
  {
    const gboolean is_free = ncm_model_param_get_ftype (model, pid) == NCM_PARAM_TYPE_FREE;

    key = _nc_galaxy_shape_factor_spectral_fnv (key, &is_free, sizeof (is_free));

    if (is_free)
    {
      const gdouble lb = ncm_model_param_get_lower_bound (model, pid);
      const gdouble ub = ncm_model_param_get_upper_bound (model, pid);

      key = _nc_galaxy_shape_factor_spectral_fnv (key, &lb, sizeof (lb));
      key = _nc_galaxy_shape_factor_spectral_fnv (key, &ub, sizeof (ub));
    }
    else
    {
      const gdouble v = ncm_model_param_get (model, pid);

      key = _nc_galaxy_shape_factor_spectral_fnv (key, &v, sizeof (v));
    }
  }

  return key;
}

/*
 * pop-order: the key of the state of the population a loaded expansion is compared with,
 * and, with use-interp and free parameters in the population, the populations of the nodes
 * of their box when the box (the type of the population, its fixed parameters, the bounds of
 * the free ones, the order) changes.
 */
static void
_nc_galaxy_shape_factor_spectral_prepare_expansion (NcGalaxyShapeFactorSpectralPrivate * const self, NcmMSet *mset)
{
  NcGalaxyShapePop *pop = NC_GALAXY_SHAPE_POP (ncm_mset_peek (mset, nc_galaxy_shape_pop_id ()));
  const guint K         = self->pop_order;
  NcGalaxyShapeFactorSpectralBox *box;
  NcmModel *model;
  guint n_free, n_nodes, i, q;
  guint64 key;

  self->pop_expand = FALSE;

  if (pop == NULL)
    return;

  self->pop_fit_key = _nc_galaxy_shape_factor_spectral_pop_box_key (pop, 0);

  if (!self->use_interp || (K < 2))
    return;

  model = NCM_MODEL (pop);
  key   = _nc_galaxy_shape_factor_spectral_pop_box_key (pop, K);

  if ((self->pop_box != NULL) && (key == self->pop_box_key))
  {
    self->pop_expand = TRUE;

    return;
  }

  box = _nc_galaxy_shape_factor_spectral_box_new (pop, K);

  if (box == NULL)
    return;

  n_free  = box->n_free;
  n_nodes = box->n_nodes;

  {
    NcmSerialize *ser        = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
    NcGalaxyShapePop **nodes = g_new (NcGalaxyShapePop *, n_nodes);
    NcGalaxyShapePop *center = NC_GALAXY_SHAPE_POP (ncm_model_dup (model, ser));

    g_clear_pointer (&self->pop_nodes, _nc_galaxy_shape_factor_spectral_nodes_free);
    g_clear_pointer (&self->pop_box, _nc_galaxy_shape_factor_spectral_box_unref);
    self->pop_box = box;

    for (i = 0; i < n_free; i++)
      ncm_model_param_set (NCM_MODEL (center), box->idx[i], 0.5 * (box->lb[i] + box->ub[i]));

    for (q = 0; q < n_nodes; q++)
    {
      guint r = q;

      nodes[q] = NC_GALAXY_SHAPE_POP (ncm_model_dup (model, ser));

      for (i = 0; i < n_free; i++)
      {
        const guint j     = r % K;
        const gdouble mid = 0.5 * (box->lb[i] + box->ub[i]);
        const gdouble hw  = 0.5 * (box->ub[i] - box->lb[i]);

        ncm_model_param_set (NCM_MODEL (nodes[q]), box->idx[i], mid + hw * cos (M_PI * j / (K - 1.0)));
        r /= K;
      }
    }

    self->pop_nodes         = _nc_galaxy_shape_factor_spectral_nodes_new (self, center, n_nodes, nodes);
    self->pop_box_key       = key;
    self->pop_n_unconverged = 0;
    self->pop_warned        = 0;

    /* The modes from the values at the nodes: the tensor product of the transform of each parameter */
    {
      gdouble *V = g_new (gdouble, (gsize) K * K);
      guint m;

      _nc_galaxy_shape_factor_spectral_cheb_transform (K, V);
      self->pop_T = g_renew (gdouble, self->pop_T, (gsize) n_nodes * n_nodes);

      for (m = 0; m < n_nodes; m++)
      {
        for (q = 0; q < n_nodes; q++)
        {
          guint rm  = m;
          guint rq  = q;
          gdouble t = 1.0;

          for (i = 0; i < n_free; i++)
          {
            t  *= V[(rm % K) * K + rq % K];
            rm /= K;
            rq /= K;
          }

          self->pop_T[(gsize) m * n_nodes + q] = t;
        }
      }

      g_free (V);
    }

    for (q = 0; q < n_nodes; q++)
      nc_galaxy_shape_pop_free (nodes[q]);

    nc_galaxy_shape_pop_free (center);
    ncm_serialize_free (ser);
    g_free (nodes);
  }

  self->pop_expand = TRUE;
}

/*
 * With restrict-range and use-interp, takes the bound of the reachable range over the box of
 * the fit (_nc_galaxy_shape_factor_reach_prepare()); it applies when the models are in the
 * fit and the box covers every free parameter, and otherwise every galaxy gets interp-g-max.
 */
static void
_nc_galaxy_shape_factor_spectral_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (NC_GALAXY_SHAPE_FACTOR_SPECTRAL (gsf));

  self->range_on = FALSE;

  if (self->restrict_range && self->use_interp)
  {
    _nc_galaxy_shape_factor_reach_prepare (&self->reach, mset);
    self->range_on = self->reach.models_ready && self->reach.holds;
  }

  _nc_galaxy_shape_factor_spectral_prepare_expansion (self, mset);
}

/* The radius in |g| an interpolant of the galaxy of ldata is built for */
static inline gdouble
_nc_galaxy_shape_factor_spectral_target (NcGalaxyShapeFactorSpectralPrivate * const self, const NcGalaxyShapeFactorSpectralData *ldata)
{
  return (ldata->g_reach >= 0.0) ? ldata->g_reach : self->interp_g_max;
}

static gdouble _nc_galaxy_shape_factor_spectral_eval_ln (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data, const gdouble std_noise, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, NcmMatrix **panels);
static gboolean _nc_galaxy_shape_factor_spectral_interp_valid (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2);
static void _nc_galaxy_shape_factor_spectral_interp_build (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2);
static void _nc_galaxy_shape_factor_spectral_interp_record (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2);
static void _nc_galaxy_shape_factor_spectral_check_loaded (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data);
static gdouble _nc_galaxy_shape_factor_spectral_interp_center (NcGalaxyShapeFactorSpectralPrivate *self, const gdouble e);

/*
 * With NcGalaxyShapeFactorSpectral:use-interp, builds the interpolant of the galaxy of data
 * for its observed ellipticity, unless it already has a valid one. The likelihood calls this
 * for different galaxies in parallel (nc_data_cluster_wl_factor_data_prepare()), so the
 * builds of a catalog run there instead of at the first evaluation of each galaxy. An
 * interpolant loaded from the catalog is compared with the fit instead
 * (_nc_galaxy_shape_factor_spectral_check_loaded()); a loaded interpolant without
 * use-interp aborts.
 */
static void
_nc_galaxy_shape_factor_spectral_data_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data, const gdouble z_max)
{
  NcGalaxyShapeFactorSpectral *gsfs               = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (gsf);
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);
  NcGalaxyShapeFactorSpectralData *ldata          = data->ldata;

  if (self->use_interp)
  {
    NcGalaxyShapePop *pop = NC_GALAXY_SHAPE_POP (ncm_mset_peek (mset, nc_galaxy_shape_pop_id ()));

    ldata->g_reach = -1.0;

    if (self->range_on)
    {
      const gdouble reach = _nc_galaxy_shape_factor_reach_eval (&self->reach, &ldata->reach_memo, data, z_max);

      if (reach < self->interp_g_max)
        ldata->g_reach = MAX (reach, MIN (NC_GALAXY_SHAPE_FACTOR_SPECTRAL_REACH_MIN, self->interp_g_max));
    }

    if (!_nc_galaxy_shape_factor_spectral_interp_valid (gsfs, pop, data, data->epsilon_obs_1, data->epsilon_obs_2))
    {
      if (ldata->loaded)
        _nc_galaxy_shape_factor_spectral_check_loaded (gsfs, pop, data);
      else
        _nc_galaxy_shape_factor_spectral_interp_build (gsfs, pop, data, data->epsilon_obs_1, data->epsilon_obs_2);
    }
  }
  else if (ldata->loaded)
  {
    g_error ("nc_galaxy_shape_factor_spectral: the catalog carries interpolants, but use-interp is not set.");
  }
}

static gdouble
_nc_galaxy_shape_factor_spectral_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  /* The exponential of eval_nodes, so that both take the exponential of the same function */
  return ncm_vmath_exp (_nc_galaxy_shape_factor_spectral_eval_ln_marginal (gsf, pop, data, g_1, g_2, epsilon_obs_1, epsilon_obs_2));
}

/*
 * The coordinates of the interpolant of the galaxy of ldata at the shear (g_1, g_2): the
 * shear in the frame where eps_obs is on the positive real axis, relative to g_c, in polar
 * coordinates, with the angle in [0, pi] (the interpolant covers psi < 0 by reflection).
 */
static inline void
_nc_galaxy_shape_factor_spectral_interp_coords (const NcGalaxyShapeFactorSpectralData *ldata, const gdouble g_1, const gdouble g_2,
                                                gdouble *rho, gdouble *psi)
{
  const gdouble gf_1 = g_1 * ldata->cos_e + g_2 * ldata->sin_e - ldata->g_c;
  const gdouble gf_2 = g_2 * ldata->cos_e - g_1 * ldata->sin_e;

  *rho = sqrt (gf_1 * gf_1 + gf_2 * gf_2);
  *psi = ncm_vmath_atan2 (fabs (gf_2), gf_1);
}

/*
 * Makes the interpolant of the galaxy of data valid for these eps_obs: builds it when it is
 * not; a loaded one, which a fit does not rebuild, is compared with the fit again
 * (_nc_galaxy_shape_factor_spectral_check_loaded(), as nc_galaxy_shape_factor_data_prepare()
 * does, so that a fixed population parameter may change inside the box of a loaded
 * expansion), and the call aborts when it is still not valid.
 */
static void
_nc_galaxy_shape_factor_spectral_interp_ensure (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                                const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorSpectralData *ldata = data->ldata;

  if (_nc_galaxy_shape_factor_spectral_interp_valid (gsfs, pop, data, epsilon_obs_1, epsilon_obs_2))
    return;

  if (ldata->loaded)
    _nc_galaxy_shape_factor_spectral_check_loaded (gsfs, pop, data);

  if (ldata->loaded && !_nc_galaxy_shape_factor_spectral_interp_valid (gsfs, pop, data, epsilon_obs_1, epsilon_obs_2))
    g_error ("nc_galaxy_shape_factor_spectral: the interpolant loaded from the catalog is not valid for this evaluation "
             "(eps_obs = (% 22.15g, % 22.15g)); a loaded interpolant is not rebuilt.", epsilon_obs_1, epsilon_obs_2);

  _nc_galaxy_shape_factor_spectral_interp_build (gsfs, pop, data, epsilon_obs_1, epsilon_obs_2);
}

/*
 * The evaluations made without the interpolant of a galaxy while use-interp is set, by
 * reason, for nc_galaxy_shape_factor_spectral_get_fallback_stats(): lines tabulated
 * (nc_galaxy_shape_factor_eval_line()), shears integrated one by one, and the time they took,
 * in microseconds. They are counted for the process, since each copy of a likelihood a
 * sampler makes has its own shape factor; atomics, since the copies run in threads. The
 * first fallback outside the coverage of the interpolants (POPULATION or RANGE) is reported
 * once with g_message(): CAP, the shears above interp-g-max, happens at ordinary points.
 */
typedef struct _NcGalaxyShapeFactorSpectralFallbackStats
{
  gsize lines;
  gsize shears;
  gsize usec;
} NcGalaxyShapeFactorSpectralFallbackStats;

G_STATIC_ASSERT (sizeof (gsize) == sizeof (gpointer));

static NcGalaxyShapeFactorSpectralFallbackStats _nc_galaxy_shape_factor_spectral_fallback_stats[NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_LEN];
static gint _nc_galaxy_shape_factor_spectral_fallback_noted = 0;

static void
_nc_galaxy_shape_factor_spectral_fallback_record (const NcGalaxyShapeFactorSpectralFallback reason, const gsize lines, const gsize shears,
                                                  const gint64 usec)
{
  NcGalaxyShapeFactorSpectralFallbackStats *st = &_nc_galaxy_shape_factor_spectral_fallback_stats[reason];

  g_atomic_pointer_add (&st->lines, (gssize) lines);
  g_atomic_pointer_add (&st->shears, (gssize) shears);
  g_atomic_pointer_add (&st->usec, (gssize) MAX (usec, 0));

  if ((reason != NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_CAP) &&
      g_atomic_int_compare_and_exchange (&_nc_galaxy_shape_factor_spectral_fallback_noted, 0, 1))
    g_message ("# nc_galaxy_shape_factor_spectral: ln P of some galaxies is computed without their interpolants, %s: "
               "from a table along their line of shears where possible, otherwise by the integrals, at a higher cost. "
               "nc_galaxy_shape_factor_spectral_get_fallback_stats() counts them.\n",
               (reason == NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_POPULATION) ?
               "with the population outside the box of their expansion (or an expansion that did not converge)" :
               "at shears beyond the range of their interpolants");
}

/* ln P by the integral, at a shear the interpolant does not cover, counted */
static gdouble
_nc_galaxy_shape_factor_spectral_fallback_eval_ln (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                                   const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2,
                                                   const NcGalaxyShapeFactorSpectralFallback reason)
{
  const gint64 t0    = g_get_monotonic_time ();
  const gdouble ln_P = _nc_galaxy_shape_factor_spectral_eval_ln (gsfs, pop, data->pop_data, data->std_noise, g_1, g_2, epsilon_obs_1, epsilon_obs_2, NULL);

  _nc_galaxy_shape_factor_spectral_fallback_record (reason, 0, 1, g_get_monotonic_time () - t0);

  return ln_P;
}

/* Flags of _nc_galaxy_shape_factor_spectral_warn_once() */
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_WARNED_SHORT (1u << 0)
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_WARNED_OUTSIDE (1u << 1)

G_LOCK_DEFINE_STATIC (nc_galaxy_shape_factor_spectral_warned);

/*
 * Whether the warning flag has not been given yet for the catalog whose record is record, and
 * marks it given. The flags are kept on the record, which the copies of a likelihood share
 * with the catalog, so each catalog is warned about once in a process.
 */
static gboolean
_nc_galaxy_shape_factor_spectral_warn_once (gpointer record, const guint flag)
{
  static GQuark quark = 0;
  gboolean first;
  guint flags;

  G_LOCK (nc_galaxy_shape_factor_spectral_warned);

  if (G_UNLIKELY (quark == 0))
    quark = g_quark_from_static_string ("nc-galaxy-shape-factor-spectral-warned");

  flags = GPOINTER_TO_UINT (g_object_get_qdata (G_OBJECT (record), quark));
  first = !(flags & flag);

  if (first)
    g_object_set_qdata (G_OBJECT (record), quark, GUINT_TO_POINTER (flags | flag));

  G_UNLOCK (nc_galaxy_shape_factor_spectral_warned);

  return first;
}

/*
 * ln P (eps_obs | g) for the galaxy of data. With NcGalaxyShapeFactorSpectral:use-interp
 * and |g| <= NcGalaxyShapeFactorSpectral:interp-g-max it is the value of the interpolant
 * of the galaxy in g, which is built first when the galaxy has none for these eps_obs,
 * std_noise, population parameters and tolerances. Otherwise it is the value of
 * _nc_galaxy_shape_factor_spectral_eval_ln(). The interpolant is evaluated at
 * (rho, psi), g = g_c + rho exp (i psi) in the frame where eps_obs is on the positive
 * real axis, with the sign of the imaginary part of g dropped: ln P is even under that
 * reflection.
 */
static gdouble
_nc_galaxy_shape_factor_spectral_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorSpectral *gsfs               = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (gsf);
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  if (!self->use_interp)
    return _nc_galaxy_shape_factor_spectral_eval_ln (gsfs, pop, data->pop_data, data->std_noise, g_1, g_2, epsilon_obs_1, epsilon_obs_2, NULL);

  if (g_1 * g_1 + g_2 * g_2 > gsl_pow_2 (self->interp_g_max))
    return _nc_galaxy_shape_factor_spectral_fallback_eval_ln (gsfs, pop, data, g_1, g_2, epsilon_obs_1, epsilon_obs_2,
                                                              NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_CAP);

  {
    NcGalaxyShapeFactorSpectralData *ldata = data->ldata;
    NcGalaxyShapeFactorSpectralWork **ptr;
    const gdouble *w;
    gdouble rho, psi, val;

    _nc_galaxy_shape_factor_spectral_interp_ensure (gsfs, pop, data, epsilon_obs_1, epsilon_obs_2);

    /* A galaxy whose expansion did not converge, or a population outside the box of the expansion */
    if (!_nc_galaxy_shape_factor_spectral_usable (ldata, pop, &w))
      return _nc_galaxy_shape_factor_spectral_fallback_eval_ln (gsfs, pop, data, g_1, g_2, epsilon_obs_1, epsilon_obs_2,
                                                                NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_POPULATION);

    if (g_1 * g_1 + g_2 * g_2 > gsl_pow_2 (ldata->g_max))
      return _nc_galaxy_shape_factor_spectral_fallback_eval_ln (gsfs, pop, data, g_1, g_2, epsilon_obs_1, epsilon_obs_2,
                                                                NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_RANGE);

    /* The modes combined on the panel of the point only */
    ptr = ncm_memory_pool_get (self->work_pool);
    _nc_galaxy_shape_factor_spectral_interp_coords (ldata, g_1, g_2, &rho, &psi);
    ncm_spectral_polar_panels_pack_eval_on_array_combine (ldata->mode_pack, w,
                                                          _nc_galaxy_shape_factor_spectral_work_comb (*ptr, ldata->modes_id, ldata->w_gen),
                                                          ncm_spectral_polar_panels_pack_find_panel_near (ldata->mode_pack, G_MAXUINT, rho, psi),
                                                          1, &rho, &psi, &val);
    ncm_memory_pool_return (ptr);

    return val;
  }
}

/*
 * Most bytes of the pack of a galaxy data_prefetch() fetches: the structure, the panels, the
 * edges and the first terms. On HWL16a-002 (about 20 kB a galaxy) the whole pack is no
 * faster, and one line of it 8% slower; the hardware prefetcher follows the rest.
 */
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_PREFETCH_CAP 8192

/*
 * data_prefetch: stage 1 the state of the galaxy, stage 2 the pack of its modes and their
 * weights, what _nc_galaxy_shape_factor_spectral_eval_nodes() reads for a pop-order
 * expansion. The pack is one allocation, so that it is reached from the state in one step.
 */
static void
_nc_galaxy_shape_factor_spectral_data_prefetch (NcGalaxyShapeFactor *gsf, NcGalaxyShapeFactorData *data, const guint stage)
{
  const NcGalaxyShapeFactorSpectralData *ldata = (const NcGalaxyShapeFactorSpectralData *) data->ldata;

  if (stage == 1)
  {
    ncm_prefetch_span (ldata, sizeof (NcGalaxyShapeFactorSpectralData));
  }
  else if ((stage == 2) && (ldata->mode_pack != NULL))
  {
    ncm_prefetch_span (ldata->mode_pack, MIN (NC_GALAXY_SHAPE_FACTOR_SPECTRAL_PREFETCH_CAP,
                                              ncm_spectral_polar_panels_pack_get_size (ldata->mode_pack)));

    if (ldata->box != NULL)
      ncm_prefetch_span (ldata->w, sizeof (gdouble) * ldata->box->n_nodes);
  }
}

/* Points per pass of _nc_galaxy_shape_factor_spectral_eval_nodes() */
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NODE_CHUNK 64

/*
 * eval_nodes: ln P (or P) at n shears of one galaxy, the values of
 * _nc_galaxy_shape_factor_spectral_eval_ln_marginal() at each shear within rounding (see
 * nc_galaxy_shape_factor_eval_nodes()). The interpolant is checked (and built when it is
 * not valid) once, at the first shear inside interp-g-max, instead of once per shear. The
 * shears are taken in passes: the ones the interpolant does not cover are integrated
 * (_nc_galaxy_shape_factor_spectral_fallback_eval_ln(): above interp-g-max, beyond the
 * range of the interpolant, or every shear of a galaxy that cannot use it), the coordinates
 * of the others are computed in one loop that vectorizes, and then the interpolant is
 * evaluated at them.
 */
static void
_nc_galaxy_shape_factor_spectral_eval_nodes (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const guint n,
                                             const gdouble *g_1, const gdouble *g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2,
                                             const gboolean want_log, gdouble *out)
{
  NcGalaxyShapeFactorSpectral *gsfs               = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (gsf);
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);
  NcGalaxyShapeFactorSpectralData *ldata          = data->ldata;
  const gdouble cap2                              = gsl_pow_2 (self->interp_g_max);
  gboolean checked                                = FALSE;
  gboolean usable                                 = TRUE;
  const gdouble *w                                = NULL;
  guint panel                                     = G_MAXUINT;
  guint j0;

  if (!self->use_interp)
  {
    guint j;

    for (j = 0; j < n; j++)
    {
      const gdouble ln_P = _nc_galaxy_shape_factor_spectral_eval_ln (gsfs, pop, data->pop_data, data->std_noise, g_1[j], g_2[j],
                                                                     epsilon_obs_1, epsilon_obs_2, NULL);

      out[j] = want_log ? ln_P : ncm_vmath_exp (ln_P);
    }

    return;
  }

  for (j0 = 0; j0 < n; j0 += NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NODE_CHUNK)
  {
    const guint j1 = MIN (n, j0 + NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NODE_CHUNK);
    gdouble ig_1[NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NODE_CHUNK];
    gdouble ig_2[NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NODE_CHUNK];
    gdouble rho[NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NODE_CHUNK];
    gdouble psi[NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NODE_CHUNK];
    guint idx[NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NODE_CHUNK];
    guint c = 0;
    guint j, k, pass;

    for (j = j0; j < j1; j++)
    {
      ig_1[c] = g_1[j];
      ig_2[c] = g_2[j];
      idx[c]  = j;
      c++;
    }

    /*
     * Pass 0 leaves out the shears above interp-g-max; pass 1, once the interpolant is
     * checked, those beyond its range. Every shear of a galaxy that cannot use its
     * interpolant is left out (lim2 < 0).
     */
    for (pass = 0; pass < 2; pass++)
    {
      NcGalaxyShapeFactorSpectralFallback reason;
      gdouble lim2;
      guint c_in = 0;

      if (pass == 1)
      {
        if (c == 0)
          break;

        if (!checked)
        {
          _nc_galaxy_shape_factor_spectral_interp_ensure (gsfs, pop, data, epsilon_obs_1, epsilon_obs_2);
          usable  = _nc_galaxy_shape_factor_spectral_usable (ldata, pop, &w);
          checked = TRUE;
        }
      }

      if (!usable)
      {
        lim2   = -1.0;
        reason = NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_POPULATION;
      }
      else if (pass == 0)
      {
        lim2   = cap2;
        reason = NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_CAP;
      }
      else
      {
        lim2   = gsl_pow_2 (ldata->g_max);
        reason = NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_RANGE;
      }

      for (k = 0; k < c; k++)
      {
        if (ig_1[k] * ig_1[k] + ig_2[k] * ig_2[k] > lim2)
        {
          const gdouble ln_P = _nc_galaxy_shape_factor_spectral_fallback_eval_ln (gsfs, pop, data, ig_1[k], ig_2[k], epsilon_obs_1, epsilon_obs_2, reason);

          out[idx[k]] = want_log ? ln_P : ncm_vmath_exp (ln_P);
        }
        else
        {
          ig_1[c_in] = ig_1[k];
          ig_2[c_in] = ig_2[k];
          idx[c_in]  = idx[k];
          c_in++;
        }
      }

      c = c_in;
    }

    if (c == 0)
      continue;

    #pragma omp simd

    for (k = 0; k < c; k++)
      _nc_galaxy_shape_factor_spectral_interp_coords (ldata, ig_1[k], ig_2[k], &rho[k], &psi[k]);

    /*
     * ln P, kept in rho, which is no longer needed. The shears of the nodes of a galaxy
     * lie on a line, so consecutive nodes are usually in the same panel (a run of about 30
     * on average): each run is evaluated in one call, which loads the coefficients of the
     * panel once for several nodes and combines the modes on that panel only.
     */
    {
      const NcmSpectralPolarPanelsPack *pack = ldata->mode_pack;
      NcGalaxyShapeFactorSpectralWork **ptr  = ncm_memory_pool_get (self->work_pool);
      NcmSpectralPolarPanelsComb *comb       = _nc_galaxy_shape_factor_spectral_work_comb (*ptr, ldata->modes_id, ldata->w_gen);

      k     = 0;
      panel = ncm_spectral_polar_panels_pack_find_panel_near (pack, panel, rho[0], psi[0]);

      while (k < c)
      {
        guint k1   = k + 1;
        guint next = panel;

        while (k1 < c)
        {
          next = ncm_spectral_polar_panels_pack_find_panel_near (pack, panel, rho[k1], psi[k1]);

          if (next != panel)
            break;

          k1++;
        }

        ncm_spectral_polar_panels_pack_eval_on_array_combine (pack, w, comb, panel, k1 - k, &rho[k], &psi[k], &rho[k]);

        k     = k1;
        panel = next;
      }

      ncm_memory_pool_return (ptr);
    }

    if (!want_log)
    {
      #pragma omp simd

      for (k = 0; k < c; k++)
        rho[k] = ncm_vmath_exp (rho[k]);
    }

    for (k = 0; k < c; k++)
      out[idx[k]] = rho[k];
  }
}

/* What the table along a line needs */
typedef struct _NcGalaxyShapeFactorSpectralLine
{
  NcGalaxyShapeFactorSpectral *gsfs;
  NcGalaxyShapePop *pop;
  NcGalaxyShapeFactorData *data;
  gdouble u_1;
  gdouble u_2;
  gdouble b_1;
  gdouble b_2;
  gdouble epsilon_obs_1;
  gdouble epsilon_obs_2;
} NcGalaxyShapeFactorSpectralLine;

/* ln P at the position t of the line, by the integral */
static gdouble
_nc_galaxy_shape_factor_spectral_line_ln_P (gpointer user_data, gdouble t)
{
  const NcGalaxyShapeFactorSpectralLine *line = user_data;

  return _nc_galaxy_shape_factor_spectral_eval_ln (line->gsfs, line->pop, line->data->pop_data, line->data->std_noise,
                                                   line->b_1 + t * line->u_1, line->b_2 + t * line->u_2,
                                                   line->epsilon_obs_1, line->epsilon_obs_2, NULL);
}

/*
 * The Chebyshev series of ln P at b + t u in t on [t_min, t_max], in *coeffs, expanded from
 * integrals on nested Chebyshev-Lobatto grids from level NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_LEVEL_MIN
 * to level_cap and accepted at NcGalaxyShapeFactorSpectral:line-abstol (reltol / 10 when
 * negative). FALSE, with *coeffs NULL, with the interpolant, without line-abstol, for a
 * level_cap too low or a single point, when the line leaves
 * |g| <= NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_G_MAX, or when the series does not reach the
 * tolerance. Nothing is kept from one call to the next, so the values depend only on the
 * line given. With fallback the table is built with the interpolant too, for a galaxy that
 * cannot use it on this line.
 */
static gboolean
_nc_galaxy_shape_factor_spectral_line_expand (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                              const gdouble t_min, const gdouble t_max, const gdouble u_1, const gdouble u_2,
                                              const gdouble b_1, const gdouble b_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2,
                                              const guint level_cap, const gboolean fallback, GArray **coeffs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);
  const gdouble line_abstol                       = (self->line_abstol < 0.0) ? 0.1 * self->reltol : self->line_abstol;

  *coeffs = NULL;

  /* |b + t u| is largest at an end of the segment */
  if ((fallback || !self->use_interp) && (line_abstol > 0.0) &&
      (level_cap > NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_LEVEL_MIN) && (t_max > t_min) &&
      (gsl_hypot (b_1 + t_min * u_1, b_2 + t_min * u_2) <= NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_G_MAX) &&
      (gsl_hypot (b_1 + t_max * u_1, b_2 + t_max * u_2) <= NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_G_MAX))
  {
    NcGalaxyShapeFactorSpectralLine line = {gsfs, pop, data, u_1, u_2, b_1, b_2, epsilon_obs_1, epsilon_obs_2};
    NcmSpectral **spectral               = ncm_memory_pool_get (self->line_pool);
    gboolean converged                   = FALSE;

    ncm_spectral_compute_chebyshev_coeffs_adaptive_try (*spectral, &_nc_galaxy_shape_factor_spectral_line_ln_P, t_min, t_max,
                                                        NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_LEVEL_MIN, level_cap,
                                                        0.0, line_abstol, coeffs, &line, &converged);
    ncm_memory_pool_return (spectral);

    if (converged)
      return TRUE;

    g_clear_pointer (coeffs, g_array_unref);
  }

  return FALSE;
}

/*
 * line_table: the series of _nc_galaxy_shape_factor_spectral_line_expand() up to the
 * largest level whose grid has at most max_points points, and at most
 * NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_LEVEL_MAX. With use-interp the interpolant is used,
 * and there is no table, except for a galaxy that cannot use it anywhere on the line: a
 * population outside the box of its expansion, or an expansion that did not converge. Its
 * line is then tabulated as without use-interp, and counted
 * (nc_galaxy_shape_factor_spectral_get_fallback_stats()).
 */
static gboolean
_nc_galaxy_shape_factor_spectral_line_table (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                             const gdouble t_min, const gdouble t_max, const gdouble u_1, const gdouble u_2,
                                             const gdouble b_1, const gdouble b_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2,
                                             const guint max_points, GArray **coeffs)
{
  NcGalaxyShapeFactorSpectral *gsfs               = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (gsf);
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);
  guint level_cap                                 = 0;

  while ((level_cap < NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_LEVEL_MAX) && ((1u << (level_cap + 1)) + 1 <= max_points))
    level_cap++;

  if (self->use_interp)
  {
    const gdouble *w;
    gboolean done;
    gint64 t0;

    *coeffs = NULL;
    _nc_galaxy_shape_factor_spectral_interp_ensure (gsfs, pop, data, epsilon_obs_1, epsilon_obs_2);

    if (_nc_galaxy_shape_factor_spectral_usable (data->ldata, pop, &w))
      return FALSE;

    t0   = g_get_monotonic_time ();
    done = _nc_galaxy_shape_factor_spectral_line_expand (gsfs, pop, data, t_min, t_max, u_1, u_2, b_1, b_2,
                                                         epsilon_obs_1, epsilon_obs_2, level_cap, TRUE, coeffs);
    _nc_galaxy_shape_factor_spectral_fallback_record (NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_POPULATION, done ? 1 : 0, 0,
                                                      g_get_monotonic_time () - t0);

    return done;
  }

  return _nc_galaxy_shape_factor_spectral_line_expand (gsfs, pop, data, t_min, t_max, u_1, u_2, b_1, b_2,
                                                       epsilon_obs_1, epsilon_obs_2, level_cap, FALSE, coeffs);
}

/*
 * Whether the galaxy of data has an interpolant built for these eps_obs, its std_noise,
 * the current state of pop and the current settings of the class.
 */
static gboolean
_nc_galaxy_shape_factor_spectral_interp_valid (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                               const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);
  const NcGalaxyShapeFactorSpectralData *ldata    = data->ldata;

  gboolean pop_ok;

  if (!ldata->valid || (ldata->epsilon_obs_1 != epsilon_obs_1) || (ldata->epsilon_obs_2 != epsilon_obs_2) ||
      (ldata->std_noise != data->std_noise))
    return FALSE;

  /*
   * The population: a loaded expansion holds while the free parameters, their bounds and the
   * fixed values it was compared with do not change; a built expansion (or one that did not
   * converge) while the box of the fit does not; the interpolant of a fixed population (no
   * box) while the parameters do not. A galaxy built one way is not valid for the other.
   */
  if (ldata->box != NULL)
    pop_ok = (ldata->box_key == (ldata->loaded ? self->pop_fit_key : self->pop_box_key)) && (ldata->loaded || self->pop_expand);
  else if (ldata->expand_failed)
    pop_ok = self->pop_expand && (ldata->box_key == self->pop_box_key);
  else
    pop_ok = (ldata->mode_spp != NULL) && (ldata->loaded || !self->pop_expand) &&
             (ldata->pop_pkey == ncm_model_state_get_pkey (NCM_MODEL (pop)));

  return pop_ok &&
         (ldata->loaded ? (ldata->g_max_setting == self->interp_g_max) : (ldata->g_max == _nc_galaxy_shape_factor_spectral_target (self, ldata))) &&
         (ldata->abstol == self->interp_abstol) && (ldata->reltol == self->reltol) &&
         (ldata->min_radial_panels == self->interp_min_radial_panels);
}

static void _nc_galaxy_shape_factor_spectral_interp_kernel (gpointer user_data, gdouble rho, guint n, const gdouble *psi, const gdouble *cos_psi, const gdouble *sin_psi, gdouble *fval);

/* What the function approximated in g needs */
typedef struct _NcGalaxyShapeFactorSpectralInterp
{
  NcGalaxyShapeFactorSpectral *gsfs;
  NcGalaxyShapePop *pop;
  NcGalaxyShapePopData *pop_data;
  gdouble std_noise;
  gdouble e;   /* |eps_obs| */
  gdouble g_c; /* center of the coordinates */
  GHashTable *memo;
} NcGalaxyShapeFactorSpectralInterp;

static void _nc_galaxy_shape_factor_spectral_eval_ln_nodes (NcGalaxyShapeFactorSpectral *gsfs, const NcGalaxyShapeFactorSpectralNodes *nodes, const gdouble std_noise, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, gdouble *out, gdouble *layout);

/*
 * The function the interpolants of the nodes of pop-order expand: ln P of every population
 * of the nodes at g = g_c + rho exp (i psi), in the frame where eps_obs is real.
 */
static void
_nc_galaxy_shape_factor_spectral_interp_kernel_nodes (gpointer user_data, gdouble rho, guint n, const gdouble *psi,
                                                      const gdouble *cos_psi, const gdouble *sin_psi, guint n_out, gdouble *fval)
{
  const NcGalaxyShapeFactorSpectralInterp *interp = user_data;
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (interp->gsfs);
  gdouble *entry                                  = g_new (gdouble, n_out + 2);
  guint j, q;

  for (j = 0; j < n; j++)
  {
    const gdouble *val;

    entry[0] = rho;
    entry[1] = psi[j];
    val      = (interp->memo != NULL) ? g_hash_table_lookup (interp->memo, entry) : NULL;

    if (val == NULL)
    {
      _nc_galaxy_shape_factor_spectral_eval_ln_nodes (interp->gsfs, self->pop_nodes, interp->std_noise,
                                                      interp->g_c + rho * cos_psi[j], rho * sin_psi[j], interp->e, 0.0, entry + 2, NULL);
      val = entry;
    }

    for (q = 0; q < n_out; q++)
      fval[(gsize) q * n + j] = val[q + 2];
  }

  g_free (entry);
}

/*
 * The values of the populations of the nodes at the points the build samples, so that the
 * interpolants of the nodes on its layout need no integral: each entry is
 * (rho, psi, ln P of the n_w populations), keyed by the bits of (rho, psi).
 */
static guint
_nc_galaxy_shape_factor_spectral_memo_hash (gconstpointer key)
{
  const gdouble *x = key;
  guint64 a, b;

  memcpy (&a, &x[0], sizeof (guint64));
  memcpy (&b, &x[1], sizeof (guint64));
  a = (a * G_GUINT64_CONSTANT (0x9E3779B97F4A7C15)) ^ b;
  a = (a ^ (a >> 31)) * G_GUINT64_CONSTANT (0xBF58476D1CE4E5B9);

  return (guint) (a ^ (a >> 32));
}

static gboolean
_nc_galaxy_shape_factor_spectral_memo_equal (gconstpointer a, gconstpointer b)
{
  return memcmp (a, b, 2 * sizeof (gdouble)) == 0;
}

/*
 * Largest change of ln P, relative to interp-abstol, of the coefficients trimmed from each
 * panel of an interpolant: with the default interp-abstol (1e-5) the trimmed coefficients
 * change ln P by at most 1e-9, against a measured interpolation error of 1e-8.
 */
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_INTERP_TRIM (1.0e-4)

/*
 * The edges of the region of the interpolant of a galaxy with |eps_obs| giving the center
 * g_c, covering the shears |g| <= g_max (see _nc_galaxy_shape_factor_spectral_interp_build()),
 * in work->r_edges and work->theta_edges.
 */
static void
_nc_galaxy_shape_factor_spectral_interp_region (NcGalaxyShapeFactorSpectralPrivate * const self, NcGalaxyShapeFactorSpectralWork *work,
                                                const gdouble std_noise, const gdouble g_c, const gdouble g_max)
{
  const gdouble L       = (self->axis_ratio_var ? 0.5 : 1.0) * std_noise;
  const gdouble rho_min = (g_c < g_max) ? 0.0 : g_c - g_max;
  const gdouble rho_max = g_c + g_max;
  const gdouble psi_min = (g_c < g_max) ? 0.0 : M_PI - asin (g_max / g_c);
  const gdouble psi_max = M_PI;
  guint i;

  ncm_spectral_polar_peak_edges (self->polar_interp, rho_max, 0.0, 0.5 * M_PI, L, M_PI, work->r_edges_all, work->theta_edges_all);

  g_array_set_size (work->r_edges, 0);
  g_array_append_val (work->r_edges, rho_min);

  for (i = 0; i < work->r_edges_all->len; i++)
  {
    const gdouble rho = g_array_index (work->r_edges_all, gdouble, i);

    if ((rho > rho_min + 0.125 * std_noise) && (rho < rho_max))
      g_array_append_val (work->r_edges, rho);
  }

  /* The region ends at rho_max even when it is narrower than std_noise / 8 */
  g_array_append_val (work->r_edges, rho_max);

  g_array_set_size (work->theta_edges, 0);
  g_array_append_val (work->theta_edges, psi_min);

  if (psi_max - psi_min > 0.5 * M_PI)
  {
    const gdouble psi_mid = 0.5 * (psi_min + psi_max);

    g_array_append_val (work->theta_edges, psi_mid);
  }

  g_array_append_val (work->theta_edges, psi_max);
}

/*
 * The truncation of an expansion of order K over n_free parameters, the n_nodes modes in
 * modes. On each panel D_d is the sum of the absolute values of the coefficients of the
 * modes whose highest degree in a parameter is d, and the error of the expansion over the
 * box is estimated by the larger of D_{K-1}, the last terms, as the last coefficients of a
 * Chebyshev series are, and of D_{K-2}^2 / D_{K-3}, the terms of degree K - 1 predicted
 * from the two degrees before: a function even or odd in a parameter has every other D_d
 * zero, and the prediction then does not vanish with D_{K-1}. Below K = 4 D_{K-1} alone is
 * taken. Returns the largest estimate over the panels.
 */
static gdouble
_nc_galaxy_shape_factor_spectral_expansion_trunc (const guint K, const guint n_free, const guint n_nodes, NcmSpectralPolarPanels **modes)
{
  const guint n_panel = ncm_spectral_polar_panels_get_n_panels (modes[0]);
  gdouble *D          = g_new (gdouble, K);
  gdouble trunc       = 0.0;
  guint i, m;

  for (i = 0; i < n_panel; i++)
  {
    gdouble tail;
    guint d;

    for (d = 0; d < K; d++)
      D[d] = 0.0;

    for (m = 0; m < n_nodes; m++)
    {
      gdouble ra, rb, ta, tb;
      guint n_r, n_theta, k;
      const gdouble *c = ncm_spectral_polar_panels_peek_panel (modes[m], i, &ra, &rb, &ta, &tb, &n_r, &n_theta);
      guint r          = m;
      guint deg        = 0;

      for (d = 0; d < n_free; d++)
      {
        deg = MAX (deg, r % K);
        r  /= K;
      }

      for (k = 0; k < n_r * n_theta; k++)
        D[deg] += fabs (c[k]);
    }

    tail = D[K - 1];

    if ((K >= 4) && (D[K - 2] > 0.0))
      tail = MAX (tail, (D[K - 3] > 0.0) ? D[K - 2] * (D[K - 2] / D[K - 3]) : GSL_POSINF);

    trunc = MAX (trunc, tail);
  }

  g_free (D);

  return trunc;
}

/*
 * pop-order: the expansion over the box of the population of the interpolant whose layout
 * (the population at the center of the box) is layout. The values at the nodes on that
 * layout, then the Chebyshev modes in the population parameters, cut together so that ln P
 * changes by at most the tolerance, packed. Returns NULL, with the estimated error in
 * *trunc, when the expansion does not converge within interp-abstol, the accuracy of the
 * interpolant in g itself.
 */
static NcmSpectralPolarPanelsArray *
_nc_galaxy_shape_factor_spectral_interp_expand (NcGalaxyShapeFactorSpectralPrivate * const self, NcGalaxyShapeFactorSpectralInterp *interp,
                                                NcmSpectralPolarPanels *layout, gdouble *trunc)
{
  const guint n_nodes                   = self->pop_nodes->n_w;
  NcmSpectralPolarPanels **nodes_spp    = g_new (NcmSpectralPolarPanels *, n_nodes);
  NcmSpectralPolarPanels **modes        = g_new (NcmSpectralPolarPanels *, n_nodes);
  NcmSpectralPolarPanelsArray *mode_spp = NULL;
  guint q;

  for (q = 0; q < n_nodes; q++)
  {
    nodes_spp[q] = ncm_spectral_polar_panels_new ();
    modes[q]     = ncm_spectral_polar_panels_new ();
  }

  _ncm_spectral_polar_refit_multi (self->polar_interp, layout, &_nc_galaxy_shape_factor_spectral_interp_kernel_nodes, interp,
                                   n_nodes, nodes_spp, NULL);

  for (q = 0; q < n_nodes; q++)
    ncm_spectral_polar_panels_combine (modes[q], n_nodes, nodes_spp, self->pop_T + (gsize) q * n_nodes);

  *trunc = _nc_galaxy_shape_factor_spectral_expansion_trunc (self->pop_order, self->pop_box->n_free, n_nodes, modes);

  if (*trunc <= self->interp_abstol)
  {
    mode_spp = ncm_spectral_polar_panels_array_new ();

    /* The modes are cut as the interpolant of one population is (ncm_spectral_polar_panels_trim()) */
    ncm_spectral_polar_panels_trim_modes (n_nodes, modes, NC_GALAXY_SHAPE_FACTOR_SPECTRAL_INTERP_TRIM * self->interp_abstol);

    for (q = 0; q < n_nodes; q++)
      ncm_spectral_polar_panels_array_add (mode_spp, modes[q]);

    /* The coefficients are kept once, packed; the array is not yet shared */
    ncm_spectral_polar_panels_array_pack (mode_spp);
  }

  for (q = 0; q < n_nodes; q++)
  {
    ncm_spectral_polar_panels_free (nodes_spp[q]);
    ncm_spectral_polar_panels_free (modes[q]);
  }

  g_free (nodes_spp);
  g_free (modes);

  return mode_spp;
}

/*
 * Builds the interpolant of ln P in g of the galaxy of data, for the observed ellipticity
 * (epsilon_obs_1, epsilon_obs_2), on a region that contains the disk
 * |g| <= NcGalaxyShapeFactorSpectral:interp-g-max, and records in the data of the galaxy
 * what it was built for. In the frame where eps_obs is on the positive real axis,
 * g = g_c + rho exp (i psi), with g_c the shear that maps chi_I = 0 to eps_obs:
 * e / (1 + sqrt (1 - e^2)) for TRACE (1 for e >= 1) and e for TRACE_DET, e = |eps_obs|
 * (_nc_galaxy_shape_factor_spectral_interp_center()). ln P has its narrowest structure at
 * g_c, of width std_noise / 2 (TRACE) or std_noise (TRACE_DET). The region is
 * [0, g_max + g_c] x [0, pi] for g_c < g_max and otherwise
 * [g_c - g_max, g_c + g_max] x [pi - asin (g_max / g_c), pi]; psi < 0 is covered by the
 * reflection. The radial edges are those of ncm_spectral_polar_peak_edges() for a peak at
 * rho = 0 of that width that are inside the region, and the angular panels span pi / 2 at
 * most. ncm_spectral_polar_build() accepts a panel when the change of ln P between its two
 * grids is at most NcGalaxyShapeFactorSpectral:interp-abstol. A build that reaches the
 * maximum depth aborts with the galaxy in the message.
 *
 * With pop-order the interpolant is expanded over the box of the population
 * (_nc_galaxy_shape_factor_spectral_interp_expand()), and a galaxy whose expansion does not
 * converge keeps none: it is integrated. Otherwise it is the expansion of one mode, the
 * interpolant of the current population.
 */
static void
_nc_galaxy_shape_factor_spectral_interp_build (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                               const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);
  NcGalaxyShapeFactorSpectralData *ldata          = data->ldata;
  NcGalaxyShapeFactorSpectralWork **ptr           = ncm_memory_pool_get (self->work_pool);
  const gdouble g_max                             = _nc_galaxy_shape_factor_spectral_target (self, ldata);
  const gdouble e                                 = hypot (epsilon_obs_1, epsilon_obs_2);
  const gdouble g_c                               = _nc_galaxy_shape_factor_spectral_interp_center (self, e);
  NcGalaxyShapeFactorSpectralInterp interp        = {
    gsfs,
    self->pop_expand ? self->pop_nodes->layout_pop : pop,
    self->pop_expand ? self->pop_nodes->layout_data : data->pop_data,
    data->std_noise, e, g_c,
    self->pop_expand ? g_hash_table_new_full (&_nc_galaxy_shape_factor_spectral_memo_hash, &_nc_galaxy_shape_factor_spectral_memo_equal, g_free, NULL) : NULL
  };
  NcmSpectralPolarPanels *layout = ncm_spectral_polar_panels_new ();
  GError *error                  = NULL;

  _nc_galaxy_shape_factor_spectral_interp_region (self, *ptr, data->std_noise, g_c, g_max);

  ldata->valid = FALSE;
  ncm_spectral_polar_build (self->polar_interp, &_nc_galaxy_shape_factor_spectral_interp_kernel, &interp, (*ptr)->r_edges, (*ptr)->theta_edges,
                            0.0, self->interp_abstol, layout, NULL, &error);
  ncm_memory_pool_return (ptr);

  if (error != NULL)
    g_error ("nc_galaxy_shape_factor_spectral: the interpolant in g: %s for the population %s, std_noise = % 22.15g, "
             "eps_obs = (% 22.15g, % 22.15g), interp-g-max = % 22.15g.",
             error->message, G_OBJECT_TYPE_NAME (pop), data->std_noise, epsilon_obs_1, epsilon_obs_2, g_max);

  if (self->pop_expand)
  {
    gdouble trunc;
    NcmSpectralPolarPanelsArray *mode_spp = _nc_galaxy_shape_factor_spectral_interp_expand (self, &interp, layout, &trunc);

    ldata->box_key       = self->pop_box_key;
    ldata->expand_failed = (mode_spp == NULL);
    _nc_galaxy_shape_factor_spectral_ldata_set_modes (ldata, mode_spp, (mode_spp != NULL) ? self->pop_box : NULL);

    if (mode_spp != NULL)
    {
      ncm_spectral_polar_panels_array_free (mode_spp);
    }
    else
    {
      g_atomic_int_inc (&self->pop_n_unconverged);

      if (g_atomic_int_compare_and_exchange (&self->pop_warned, 0, 1))
        g_warning ("nc_galaxy_shape_factor_spectral: the expansion over the population parameters does not converge "
                   "for some galaxies, which are integrated instead (see "
                   "nc_galaxy_shape_factor_spectral_get_pop_n_unconverged()): for the first, eps_obs = (% 22.15g, % 22.15g), "
                   "std_noise = % 22.15g, the estimated error of its expansion reaches % 22.15g in ln P, above interp-abstol = "
                   "% 22.15g. Raise pop-order or narrow the box of the free population parameters.",
                   epsilon_obs_1, epsilon_obs_2, data->std_noise, trunc, self->interp_abstol);
    }
  }
  else
  {
    /* A fixed population: the expansion of one mode, the interpolant itself */
    NcmSpectralPolarPanelsArray *mode_spp = ncm_spectral_polar_panels_array_new ();

    ncm_spectral_polar_panels_trim (layout, NC_GALAXY_SHAPE_FACTOR_SPECTRAL_INTERP_TRIM * self->interp_abstol);
    ncm_spectral_polar_panels_array_add (mode_spp, layout);
    ncm_spectral_polar_panels_array_pack (mode_spp);

    ldata->expand_failed = FALSE;
    _nc_galaxy_shape_factor_spectral_ldata_set_modes (ldata, mode_spp, NULL);
    ncm_spectral_polar_panels_array_free (mode_spp);
  }

  ncm_spectral_polar_panels_free (layout);
  g_clear_pointer (&interp.memo, g_hash_table_unref);

  _nc_galaxy_shape_factor_spectral_interp_record (gsfs, pop, data, epsilon_obs_1, epsilon_obs_2);
}

/*
 * Marks the interpolant of the galaxy of data as valid and records what it is valid for:
 * these eps_obs, the noise of the galaxy, the current state of pop and the current settings
 * of the class, and the frame of the interpolant (the direction of eps_obs and the center).
 */
static void
_nc_galaxy_shape_factor_spectral_interp_record (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                                const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);
  NcGalaxyShapeFactorSpectralData *ldata          = data->ldata;
  const gdouble e                                 = hypot (epsilon_obs_1, epsilon_obs_2);

  ldata->valid             = TRUE;
  ldata->epsilon_obs_1     = epsilon_obs_1;
  ldata->epsilon_obs_2     = epsilon_obs_2;
  ldata->std_noise         = data->std_noise;
  ldata->pop_pkey          = ncm_model_state_get_pkey (NCM_MODEL (pop));
  ldata->g_max             = _nc_galaxy_shape_factor_spectral_target (self, ldata);
  ldata->abstol            = self->interp_abstol;
  ldata->reltol            = self->reltol;
  ldata->min_radial_panels = self->interp_min_radial_panels;
  ldata->cos_e             = (e > 0.0) ? epsilon_obs_1 / e : 1.0;
  ldata->sin_e             = (e > 0.0) ? epsilon_obs_2 / e : 0.0;
  ldata->g_c               = _nc_galaxy_shape_factor_spectral_interp_center (self, e);
  ldata->g_max_setting     = self->interp_g_max;

  /* A loaded interpolant covers |g| <= rho_max - g_c, whatever box built it */
  if (ldata->loaded)
  {
    gdouble ra, rb, ta, tb;

    ncm_spectral_polar_panels_pack_get_range (ldata->mode_pack, &ra, &rb, &ta, &tb);

    ldata->g_max = MIN (rb - ldata->g_c, self->interp_g_max);
  }
}

/*
 * Compares a loaded expansion over the population parameters with the population of the fit,
 * appending the differences to diff. The box of the expansion is that of the record: the
 * free parameters of src_pop and their bounds, at the pop-order of src. The fit may free
 * any of those parameters, over any bounds, or fix them at any value; outside the box the
 * galaxy is integrated (_nc_galaxy_shape_factor_spectral_weights()), and outside is set
 * when the fit reaches there. The parameters fixed in the record must be fixed at the same
 * values. Returns the box, or NULL with a difference.
 */
static NcGalaxyShapeFactorSpectralBox *
_nc_galaxy_shape_factor_spectral_check_loaded_box (NcGalaxyShapePop *pop, NcGalaxyShapePop *src_pop, const guint src_K,
                                                   NcmSpectralPolarPanelsArray *modes, GString *diff, gboolean *outside)
{
  const gdouble tol = NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SOURCE_RELTOL;
  NcmModel *fit     = NCM_MODEL (pop);
  NcmModel *src     = NCM_MODEL (src_pop);
  const gsize len0  = diff->len;
  NcGalaxyShapeFactorSpectralBox *box;
  guint pid;

  if (G_OBJECT_TYPE (pop) != G_OBJECT_TYPE (src_pop))
  {
    g_string_append_printf (diff, " the population (%s, built with %s);", G_OBJECT_TYPE_NAME (pop), G_OBJECT_TYPE_NAME (src_pop));

    return NULL;
  }

  box = (src_K >= 2) ? _nc_galaxy_shape_factor_spectral_box_new (src_pop, src_K) : NULL;

  if ((box == NULL) || (box->n_nodes != ncm_spectral_polar_panels_array_len (modes)))
  {
    g_string_append_printf (diff, " the record of the expansion over the population parameters (pop-order %u, %u modes);",
                            src_K, ncm_spectral_polar_panels_array_len (modes));
    g_clear_pointer (&box, _nc_galaxy_shape_factor_spectral_box_unref);

    return NULL;
  }

  for (pid = 0; pid < ncm_model_len (fit); pid++)
  {
    const gchar *name     = ncm_model_param_name (fit, pid);
    const gboolean f_free = ncm_model_param_get_ftype (fit, pid) == NCM_PARAM_TYPE_FREE;
    const gdouble v       = ncm_model_param_get (fit, pid);

    if (ncm_model_param_get_ftype (src, pid) == NCM_PARAM_TYPE_FREE)
    {
      const gdouble lb   = ncm_model_param_get_lower_bound (src, pid);
      const gdouble ub   = ncm_model_param_get_upper_bound (src, pid);
      const gdouble f_lb = f_free ? ncm_model_param_get_lower_bound (fit, pid) : v;
      const gdouble f_ub = f_free ? ncm_model_param_get_upper_bound (fit, pid) : v;

      if ((ncm_cmp (f_lb, lb, tol, 0.0) < 0) || (ncm_cmp (f_ub, ub, tol, 0.0) > 0))
        *outside = TRUE;
    }
    else if (f_free)
    {
      g_string_append_printf (diff, " the population parameter %s is free, and the expansion holds for it fixed at % 22.15g;",
                              name, ncm_model_param_get (src, pid));
    }
    else if (ncm_cmp (v, ncm_model_param_get (src, pid), tol, 0.0) != 0)
    {
      g_string_append_printf (diff, " the population parameter %s (% 22.15g, built with it fixed at % 22.15g);",
                              name, v, ncm_model_param_get (src, pid));
    }
  }

  if (diff->len > len0)
    g_clear_pointer (&box, _nc_galaxy_shape_factor_spectral_box_unref);

  return box;
}

/*
 * Compares the record of the loaded interpolant of the galaxy of data with the fit: the
 * convention, reltol and the three interp-* settings of the class, and the population: for
 * the interpolant of one population its parameters (ncm_model_params_equal()), none of them
 * free; for an expansion over the population parameters, see
 * _nc_galaxy_shape_factor_spectral_check_loaded_box(). When all agree the interpolant is
 * recorded as valid for the current state of pop; otherwise the call aborts with every
 * difference in the message.
 */
static void
_nc_galaxy_shape_factor_spectral_check_loaded (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);
  NcGalaxyShapeFactorSpectralData *ldata          = data->ldata;
  NcGalaxyShapeFactorSpectralPrivate * const src  = nc_galaxy_shape_factor_spectral_get_instance_private (ldata->src_gsfs);
  const gdouble tol                               = NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SOURCE_RELTOL;
  GString *diff                                   = g_string_new (NULL);
  NcGalaxyShapeFactorSpectralBox *box             = NULL;
  gboolean outside                                = FALSE;
  guint k;

  if (nc_galaxy_shape_factor_get_ellip_conv (NC_GALAXY_SHAPE_FACTOR (gsfs)) != nc_galaxy_shape_factor_get_ellip_conv (NC_GALAXY_SHAPE_FACTOR (ldata->src_gsfs)))
    g_string_append (diff, " the ellipticity convention;");

  if (ncm_cmp (self->reltol, src->reltol, tol, 0.0) != 0)
    g_string_append_printf (diff, " reltol (% 22.15g, built with % 22.15g);", self->reltol, src->reltol);

  if (ncm_cmp (self->interp_g_max, src->interp_g_max, tol, 0.0) != 0)
    g_string_append_printf (diff, " interp-g-max (% 22.15g, built with % 22.15g);", self->interp_g_max, src->interp_g_max);

  if (ncm_cmp (self->interp_abstol, src->interp_abstol, tol, 0.0) != 0)
    g_string_append_printf (diff, " interp-abstol (% 22.15g, built with % 22.15g);", self->interp_abstol, src->interp_abstol);

  if (self->interp_min_radial_panels != src->interp_min_radial_panels)
    g_string_append_printf (diff, " interp-min-radial-panels (%u, built with %u);", self->interp_min_radial_panels, src->interp_min_radial_panels);

  /* An expansion when the record expands over free population parameters, otherwise the one mode of a fixed population */
  if ((src->pop_order >= 2) && _nc_galaxy_shape_factor_spectral_pop_has_free (ldata->src_pop))
  {
    box = _nc_galaxy_shape_factor_spectral_check_loaded_box (pop, ldata->src_pop, src->pop_order, ldata->mode_spp, diff, &outside);
  }
  else
  {
    if (ncm_spectral_polar_panels_array_len (ldata->mode_spp) != 1)
      g_string_append_printf (diff, " the interpolant of a fixed population has %u modes;", ncm_spectral_polar_panels_array_len (ldata->mode_spp));

    if (!ncm_model_params_equal (NCM_MODEL (pop), NCM_MODEL (ldata->src_pop), tol, 0.0))
      g_string_append_printf (diff, " the population (%s, built with %s or other parameter values);",
                              G_OBJECT_TYPE_NAME (pop), G_OBJECT_TYPE_NAME (ldata->src_pop));

    for (k = 0; k < ncm_model_len (NCM_MODEL (pop)); k++)
    {
      if (ncm_model_param_get_ftype (NCM_MODEL (pop), k) == NCM_PARAM_TYPE_FREE)
        g_string_append_printf (diff, " the population parameter %s is free;", ncm_model_param_name (NCM_MODEL (pop), k));
    }
  }

  if (diff->len > 0)
    g_error ("nc_galaxy_shape_factor_spectral: the interpolants of the catalog were built for other settings:%s "
             "prepare the catalog again for this fit.", diff->str);

  g_string_free (diff, TRUE);

  _nc_galaxy_shape_factor_spectral_interp_record (gsfs, pop, data, data->epsilon_obs_1, data->epsilon_obs_2);

  /* The bounds of the fit are its prior: a population outside the box is integrated */
  if (outside && _nc_galaxy_shape_factor_spectral_warn_once (ldata->src_gsfs, NC_GALAXY_SHAPE_FACTOR_SPECTRAL_WARNED_OUTSIDE))
    g_warning ("nc_galaxy_shape_factor_spectral: the population of this fit reaches outside the box of the expansion of the "
               "catalog (its free parameters over their bounds, or a fixed value outside it); there the galaxies are "
               "computed without the expansion (from a table along their line of shears where possible, otherwise by the "
               "integrals), at a higher cost. Prepare the catalog over a box that covers the fit.");

  if (box != NULL)
  {
    _nc_galaxy_shape_factor_spectral_ldata_set_modes (ldata, ldata->mode_spp, box);
    _nc_galaxy_shape_factor_spectral_box_unref (box);
    ldata->box_key    = self->pop_fit_key;
    ldata->box_is_fit = self->pop_expand && (src->pop_order == self->pop_order) &&
                        (_nc_galaxy_shape_factor_spectral_pop_box_key (ldata->src_pop, src->pop_order) == self->pop_box_key);
  }

  /*
   * A catalog built for a smaller box than that of the fit: the shears beyond it are
   * integrated. The range of a loaded interpolant is (g_c + g_max) - g_c, rounded at the
   * scale of g_c + g_max: a few of its ulps are allowed besides the relative tolerance, or a
   * galaxy with a small reach and a large g_c would be short of its own reach.
   */
  if (ldata->g_max < (1.0 - NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SOURCE_RELTOL) * _nc_galaxy_shape_factor_spectral_target (self, ldata) -
      4.0 * GSL_DBL_EPSILON * (ldata->g_c + _nc_galaxy_shape_factor_spectral_target (self, ldata)))
  {
    g_atomic_int_inc (&self->n_short);

    if (_nc_galaxy_shape_factor_spectral_warn_once (ldata->src_gsfs, NC_GALAXY_SHAPE_FACTOR_SPECTRAL_WARNED_SHORT))
      g_warning ("nc_galaxy_shape_factor_spectral: the interpolants of the catalog do not cover the shears some galaxies "
                 "reach in the box of this fit (for the first, |g| <= % 22.15g against % 22.15g), which are integrated "
                 "instead: the values are those of the integrals, at their cost. Prepare the catalog for the box of the "
                 "fit (see nc_galaxy_shape_factor_spectral_get_n_loaded_short()).",
                 ldata->g_max, _nc_galaxy_shape_factor_spectral_target (self, ldata));
  }
}

/*
 * The center g_c of the coordinates of the interpolant for |eps_obs| = e: the shear that maps
 * chi_I = 0 to eps_obs, e / (1 + sqrt (1 - e^2)) for TRACE (1 for e >= 1) and e for
 * TRACE_DET.
 */
static gdouble
_nc_galaxy_shape_factor_spectral_interp_center (NcGalaxyShapeFactorSpectralPrivate *self, const gdouble e)
{
  if (!self->axis_ratio_var)
    return e;

  return (e < 1.0) ? e / (1.0 + sqrt (1.0 - e * e)) : 1.0;
}

/*
 * The function approximated by ncm_spectral_polar_build(): ln P at the shear
 * g_c + rho exp (i psi) for the observed ellipticity (e, 0), by the integral.
 */
static void
_nc_galaxy_shape_factor_spectral_interp_kernel (gpointer user_data, gdouble rho, guint n, const gdouble *psi,
                                                const gdouble *cos_psi, const gdouble *sin_psi, gdouble *fval)
{
  const NcGalaxyShapeFactorSpectralInterp *interp = user_data;
  guint j;

  if (interp->memo != NULL)
  {
    /* pop-order: every population of the nodes from the samples of the layout population */
    NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (interp->gsfs);
    const guint n_w                                 = self->pop_nodes->n_w;

    for (j = 0; j < n; j++)
    {
      gdouble *entry = g_new (gdouble, n_w + 2);

      entry[0] = rho;
      entry[1] = psi[j];
      _nc_galaxy_shape_factor_spectral_eval_ln_nodes (interp->gsfs, self->pop_nodes, interp->std_noise,
                                                      interp->g_c + rho * cos_psi[j], rho * sin_psi[j], interp->e, 0.0, entry + 2, &fval[j]);
      g_hash_table_add (interp->memo, entry);
    }
  }
  else
  {
    for (j = 0; j < n; j++)
      fval[j] = _nc_galaxy_shape_factor_spectral_eval_ln (interp->gsfs, interp->pop, interp->pop_data, interp->std_noise,
                                                          interp->g_c + rho * cos_psi[j], rho * sin_psi[j], interp->e, 0.0, NULL);
  }
}

static gdouble _nc_galaxy_shape_factor_spectral_one_minus_abs2 (const gdouble g_1, const gdouble g_2);
static gdouble _nc_galaxy_shape_factor_spectral_integrate_ln (NcGalaxyShapeFactorSpectralPrivate * const self, NcGalaxyShapeFactorSpectralKernel *k, const gdouble std_noise, const gdouble a, const gdouble b, NcGalaxyShapeFactorSpectralWork *work, NcmMatrix **panels, GError **error);
static gdouble _nc_galaxy_shape_factor_spectral_near_unit_ln (NcGalaxyShapeFactorSpectralPrivate * const self, NcGalaxyShapeFactorSpectralKernel *k, const gdouble std_noise, const gdouble a, const gdouble b, const gdouble delta, const gdouble lambda_inv, NcGalaxyShapeFactorSpectralWork *work, NcmMatrix **panels, GError **error);

/*
 * ln P (eps_obs | g) for the population pop with resolved data pop_data: the closed form
 * at |g| = 1, the integral at t = min (|g|, 1 / |g|) for
 * 1 - t >= NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NEAR_UNIT_DELTA, and below it the second-order
 * term plus the term scaled from the integral at that distance from 1. With panels not
 * NULL, also the accepted panels of the integral computed, with their angles in the frame
 * of g and eps_obs, or NULL at |g| = 1, where nothing is integrated. A failed integral or
 * a non-finite ln P aborts with the point in the message.
 */
static gdouble
_nc_galaxy_shape_factor_spectral_eval_ln (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data,
                                          const gdouble std_noise, const gdouble g_1, const gdouble g_2,
                                          const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, NcmMatrix **panels)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);
  const gdouble noise_var                         = gsl_pow_2 (std_noise);
  const gdouble abs_g2                            = g_1 * g_1 + g_2 * g_2;

  if (abs_g2 == 1.0)
  {
    if (panels != NULL)
      *panels = NULL;

    return -(gsl_pow_2 (epsilon_obs_1 - g_1) + gsl_pow_2 (epsilon_obs_2 - g_2)) / (2.0 * noise_var) - log (2.0 * M_PI * noise_var);
  }

  {
    /* For |g| > 1 the integral is the one at g / |g|^2, of the same phase */
    const gdouble abs_g                   = sqrt (abs_g2);
    const gdouble u                       = fabs (_nc_galaxy_shape_factor_spectral_one_minus_abs2 (g_1, g_2));
    const gdouble delta                   = (abs_g2 > 1.0) ? u / (abs_g * (1.0 + abs_g)) : u / (1.0 + abs_g);
    const gdouble lambda_inv              = u / gsl_pow_2 (1.0 + abs_g);
    const gdouble phase                   = (abs_g2 > 0.0) ? atan2 (g_2, g_1) : 0.0;
    const gdouble cos_p                   = cos (phase);
    const gdouble sin_p                   = sin (phase);
    const gdouble a                       = nc_galaxy_shape_pop_exponent_at_origin (pop);
    const gdouble b_pop                   = nc_galaxy_shape_pop_exponent_at_one (pop);
    const gdouble b                       = self->axis_ratio_var ? 2.0 * b_pop + 1.0 : b_pop;
    NcGalaxyShapeFactorSpectralWork **ptr = ncm_memory_pool_get (self->work_pool);
    NcGalaxyShapeFactorSpectralWork *work = *ptr;
    NcGalaxyShapeFactorSpectralKernel k   = {
      .pop       = pop,
      .pop_data  = pop_data,
      .t         = (abs_g2 > 1.0) ? 1.0 / abs_g : abs_g,
      .e_1       = epsilon_obs_1 * cos_p + epsilon_obs_2 * sin_p,
      .e_2       = epsilon_obs_2 * cos_p - epsilon_obs_1 * sin_p,
      .c2        = 0.5 / noise_var,
      .norm      = 1.0 / (2.0 * M_PI * noise_var) / (2.0 * M_PI),
      .shift_c2  = 0.0,
      .dens_norm = exp2 (a + 1.0),
      .dens_pow  = a + b_pop + 2.0,
      .trace     = self->axis_ratio_var
    };
    GError *error = NULL;
    gdouble ln_P;

    if (delta < NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NEAR_UNIT_DELTA)
      ln_P = _nc_galaxy_shape_factor_spectral_near_unit_ln (self, &k, std_noise, a, b, delta, lambda_inv, work, panels, &error);
    else
      ln_P = _nc_galaxy_shape_factor_spectral_integrate_ln (self, &k, std_noise, a, b, work, panels, &error);

    ncm_memory_pool_return (ptr);

    if ((error != NULL) || !gsl_finite (ln_P))
      g_error ("nc_galaxy_shape_factor_spectral: %s for the population %s, std_noise = % 22.15g, "
               "g = (% 22.15g, % 22.15g), eps_obs = (% 22.15g, % 22.15g).",
               (error != NULL) ? error->message : "non-finite ln P",
               G_OBJECT_TYPE_NAME (pop), std_noise, g_1, g_2, epsilon_obs_1, epsilon_obs_2);

    if ((panels != NULL) && (phase != 0.0))
    {
      guint i;

      for (i = 0; i < ncm_matrix_nrows (*panels); i++)
      {
        ncm_matrix_addto (*panels, i, 2, phase);
        ncm_matrix_addto (*panels, i, 3, phase);
      }
    }

    return ln_P;
  }
}

/*
 * 1 - g_1^2 - g_2^2 with the rounding errors of the two squares and of their sum removed:
 * the squares are split into their rounded value and its error by fma(), the larger
 * rounded square is subtracted from 1 first, and the errors last. Next to |g| = 1 both
 * subtractions are exact, so the result keeps the digits of the inputs.
 */
static gdouble
_nc_galaxy_shape_factor_spectral_one_minus_abs2 (const gdouble g_1, const gdouble g_2)
{
  const gdouble p   = g_1 * g_1;
  const gdouble q   = g_2 * g_2;
  const gdouble e_p = fma (g_1, g_1, -p);
  const gdouble e_q = fma (g_2, g_2, -q);
  const gdouble hi  = GSL_MAX (p, q);
  const gdouble lo  = GSL_MIN (p, q);

  return ((1.0 - hi) - lo) - (e_p + e_q);
}

static void _nc_galaxy_shape_factor_spectral_seed (NcGalaxyShapeFactorSpectralPrivate * const self, const NcGalaxyShapeFactorSpectralKernel *k, const gdouble std_noise, const gdouble a, const gdouble b, NcGalaxyShapeFactorSpectralSeed *seed);

/*
 * ln of the integral at the real shear k->t, by _ncm_spectral_polar_integrate_radial() on the
 * panels placed by the seed, with the offset of the seed set in k->shift_c2. NaN with
 * error set when a panel reaches the maximum depth.
 */
static gdouble
_nc_galaxy_shape_factor_spectral_integrate_ln (NcGalaxyShapeFactorSpectralPrivate * const self, NcGalaxyShapeFactorSpectralKernel *k,
                                               const gdouble std_noise, const gdouble a, const gdouble b,
                                               NcGalaxyShapeFactorSpectralWork *work, NcmMatrix **panels, GError **error)
{
  NcGalaxyShapeFactorSpectralSeed seed;
  gdouble L, integral;

  k->shift_c2 = 0.0;
  _nc_galaxy_shape_factor_spectral_seed (self, k, std_noise, a, b, &seed);
  k->shift_c2 = seed.shift * k->c2;

  L = std_noise * GSL_MIN (seed.J_1, seed.J_2);
  ncm_spectral_polar_peak_edges (self->polar, 1.0, seed.r_star, seed.theta_0, L, L / GSL_MAX (seed.r_star, GSL_DBL_MIN),
                                 work->r_edges, work->theta_edges);

  integral = _ncm_spectral_polar_integrate_radial (self->polar, &_nc_galaxy_shape_factor_spectral_kernel, &_nc_galaxy_shape_factor_spectral_radial,
                                                   k, 1.0, a, b, work->r_edges, work->theta_edges, self->reltol, 0.0, NULL, panels, error);

  return log (integral) - k->shift_c2;
}

static gdouble _nc_galaxy_shape_factor_spectral_second_order_ln (NcGalaxyShapeFactorSpectralPrivate * const self, const NcGalaxyShapeFactorSpectralKernel *k, const gdouble std_noise, const gdouble a, const gdouble b, const gdouble delta, NcGalaxyShapeFactorSpectralWork *work, GError **error);

/*
 * ln P at the real shear t = 1 - delta, delta < delta_s =
 * NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NEAR_UNIT_DELTA, without an integral at t:
 * P (t) = A (t) + B (t), with A the second-order term of
 * _nc_galaxy_shape_factor_spectral_second_order_ln() and
 * B (t) = lambda^-(b + 2) (k_0 + k_1 / lambda), lambda = (1 + t) / (1 - t) = 1 / lambda_inv.
 * For TRACE k_1 = 0 and k_0 is set by B = P - A at 1 - delta_s, with P the integral there.
 * For TRACE_DET k_0 and k_1 are set by B at 1 - delta_s and at 1 - 10 delta_s. Every
 * quantity is divided by P (1 - delta_s), so that none leaves the double range. panels,
 * when not NULL, receives the panels of the integral at 1 - delta_s. NaN when an integral
 * fails or when A + B is not positive.
 */
static gdouble
_nc_galaxy_shape_factor_spectral_near_unit_ln (NcGalaxyShapeFactorSpectralPrivate * const self, NcGalaxyShapeFactorSpectralKernel *k,
                                               const gdouble std_noise, const gdouble a, const gdouble b,
                                               const gdouble delta, const gdouble lambda_inv,
                                               NcGalaxyShapeFactorSpectralWork *work, NcmMatrix **panels, GError **error)
{
  const gdouble t        = k->t;
  const gdouble t_1      = 1.0 - NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NEAR_UNIT_DELTA;
  const gdouble delta_1  = 1.0 - t_1;
  const gdouble linv_1   = delta_1 / (1.0 + t_1);
  const gdouble exponent = b + 2.0;
  gdouble ln_P_1, ln_A_1, ln_A, b_1, k_0, k_1;

  k->t   = t_1;
  ln_P_1 = _nc_galaxy_shape_factor_spectral_integrate_ln (self, k, std_noise, a, b, work, panels, error);

  if (*error != NULL)
    return GSL_NAN;

  ln_A_1 = _nc_galaxy_shape_factor_spectral_second_order_ln (self, k, std_noise, a, b, delta_1, work, error);

  if (*error != NULL)
    return GSL_NAN;

  b_1 = -expm1 (ln_A_1 - ln_P_1) * pow (linv_1, -exponent);

  if (self->axis_ratio_var)
  {
    k_0 = b_1;
    k_1 = 0.0;
  }
  else
  {
    const gdouble t_2     = 1.0 - 10.0 * NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NEAR_UNIT_DELTA;
    const gdouble delta_2 = 1.0 - t_2;
    const gdouble linv_2  = delta_2 / (1.0 + t_2);
    gdouble ln_P_2, ln_A_2, b_2;

    k->t   = t_2;
    ln_P_2 = _nc_galaxy_shape_factor_spectral_integrate_ln (self, k, std_noise, a, b, work, NULL, error);

    if (*error != NULL)
      return GSL_NAN;

    ln_A_2 = _nc_galaxy_shape_factor_spectral_second_order_ln (self, k, std_noise, a, b, delta_2, work, error);

    if (*error != NULL)
      return GSL_NAN;

    b_2 = (exp (ln_P_2 - ln_P_1) - exp (ln_A_2 - ln_P_1)) * pow (linv_2, -exponent);
    k_1 = (b_1 - b_2) / (linv_1 - linv_2);
    k_0 = b_1 - k_1 * linv_1;
  }

  k->t = t;
  ln_A = _nc_galaxy_shape_factor_spectral_second_order_ln (self, k, std_noise, a, b, delta, work, error);

  if (*error != NULL)
    return GSL_NAN;

  return ln_P_1 + log (exp (ln_A - ln_P_1) + pow (lambda_inv, exponent) * (k_0 + k_1 * lambda_inv));
}

static void _nc_galaxy_shape_factor_spectral_var_kernel (gpointer user_data, gdouble v, guint n, const gdouble *theta, const gdouble *cos_theta, const gdouble *sin_theta, gdouble *fval);
static void _nc_galaxy_shape_factor_spectral_var_radial (gpointer user_data, guint n, const gdouble *v, gdouble *rho);

/* What the integrand of the variance needs */
typedef struct _NcGalaxyShapeFactorSpectralVar
{
  const NcGalaxyShapeFactorSpectralKernel *k;
  gdouble delta; /* 1 - t */
  gdouble *r;    /* the radii r = 1 - v of a call of the radial factor, 2^level-cap + 1 long */
} NcGalaxyShapeFactorSpectralVar;

/*
 * ln of A (t) = Phi (t) + s (t) Lap Phi (t) / 2 at t = k->t = 1 - delta, the marginal to
 * second order in the distance of the sheared ellipticity to its mean. w = m_t (z) is the
 * TRACE_DET map of the variable of integration z, which has the mean t and the variance
 * s (t) per component,
 *   s (t) = (1 - t^2)^2 / 2 INT_0^1 (density of r) r^2 / (1 - t^2 r^2) dr,
 * and Phi (w) = N_2 (eps_obs - C (w)), C the identity for TRACE_DET and T for TRACE. The
 * integral is computed in v = 1 - r by _ncm_spectral_polar_integrate_radial(), where the factor
 * 1 / (delta + t v) of 1 / (1 - t^2 r^2) varies on the scale delta, with radial edges
 * placed for that scale at v = 0. The integrand depends only on r, so all of it is the
 * radial factor, computed once per radius. Differences to 1 are formed from delta.
 */
static gdouble
_nc_galaxy_shape_factor_spectral_second_order_ln (NcGalaxyShapeFactorSpectralPrivate * const self, const NcGalaxyShapeFactorSpectralKernel *k,
                                                  const gdouble std_noise, const gdouble a, const gdouble b, const gdouble delta,
                                                  NcGalaxyShapeFactorSpectralWork *work, GError **error)
{
  const gdouble t                    = k->t;
  const gdouble noise_var            = gsl_pow_2 (std_noise);
  const gdouble one_minus_t2         = delta * (1.0 + t);
  NcGalaxyShapeFactorSpectralVar var = {k, delta, g_new (gdouble, (1u << ncm_spectral_polar_get_level_cap (self->polar)) + 1)};
  gdouble s, d, lap_ratio;

  ncm_spectral_polar_peak_edges (self->polar, 1.0, 0.0, 0.0, delta, M_PI, work->r_edges, work->theta_edges);
  s = 0.5 * gsl_pow_2 (one_minus_t2) *
      _ncm_spectral_polar_integrate_radial (self->polar, &_nc_galaxy_shape_factor_spectral_var_kernel, &_nc_galaxy_shape_factor_spectral_var_radial,
                                            &var, 1.0, b, a, work->r_edges, work->theta_edges, self->reltol, 0.0, NULL, NULL, error);
  g_free (var.r);

  if (self->axis_ratio_var)
  {
    /* C (t) = 2 t / (1 + t^2) = 1 - delta^2 / (1 + t^2) */
    const gdouble q     = 1.0 + t * t;
    const gdouble d_1   = (k->e_1 - 1.0) + delta * delta / q;
    const gdouble d_2   = k->e_2;
    const gdouble gr_1  = -4.0 * d_1 * one_minus_t2 / (q * q);
    const gdouble gr_2  = -4.0 * d_2 / q;
    const gdouble lap_d = 8.0 * gsl_pow_2 (one_minus_t2 / (q * q)) + 8.0 / (q * q) + 32.0 * t * d_1 / gsl_pow_3 (q);

    d         = d_1 * d_1 + d_2 * d_2;
    lap_ratio = -0.5 * lap_d / noise_var + 0.25 * (gr_1 * gr_1 + gr_2 * gr_2) / gsl_pow_2 (noise_var);
  }
  else
  {
    const gdouble d_1 = (k->e_1 - 1.0) + delta;
    const gdouble d_2 = k->e_2;

    d         = d_1 * d_1 + d_2 * d_2;
    lap_ratio = d / gsl_pow_2 (noise_var) - 2.0 / noise_var;
  }

  return -0.5 * d / noise_var - log (2.0 * M_PI * noise_var) + log1p (0.5 * s * lap_ratio);
}

/*
 * The integrand of the variance for _ncm_spectral_polar_integrate_radial() in v = 1 - r, with
 * the weight v^b (1 - v)^a, is (reduced density of r) r^2 / ((delta + t v) (1 + t r)) / (2 pi):
 * it does not depend on the angle, and all of it is the radial factor
 * (_nc_galaxy_shape_factor_spectral_var_radial()). The rest is 1.
 */
static void
_nc_galaxy_shape_factor_spectral_var_kernel (gpointer user_data, gdouble v, guint n, const gdouble *theta,
                                             const gdouble *cos_theta, const gdouble *sin_theta, gdouble *fval)
{
  guint j;

  for (j = 0; j < n; j++)
    fval[j] = 1.0;
}

/* The integrand of the variance at the n values v[i] of v = 1 - r, in rho[i] */
static void
_nc_galaxy_shape_factor_spectral_var_radial (gpointer user_data, guint n, const gdouble *v, gdouble *rho)
{
  const NcGalaxyShapeFactorSpectralVar *var = user_data;
  const gdouble t                           = var->k->t;
  const gdouble delta                       = var->delta;
  gdouble *r                                = var->r;
  guint i;

  for (i = 0; i < n; i++)
    r[i] = 1.0 - v[i];

  _nc_galaxy_shape_factor_spectral_density (var->k, n, r, rho);

  for (i = 0; i < n; i++)
    rho[i] *= r[i] * r[i] / ((delta + t * v[i]) * (1.0 + t * r[i])) / (2.0 * M_PI);
}

static void _nc_galaxy_shape_factor_spectral_h (const NcGalaxyShapeFactorSpectralKernel *k, const gdouble a, const gdouble b, const guint n, const gdouble *r, const gdouble cos_0, const gdouble sin_0, gdouble *h, gdouble *d, gdouble *x);
static void _nc_galaxy_shape_factor_spectral_sort (gdouble *x, const guint n);

/*
 * Locates the peak of the integrand on the ray through c, the preimage of eps_obs by the
 * map of the variable of integration, or of the point of the unit circle nearest to
 * eps_obs when |eps_obs| >= 1: the angle of the ray, the stretch factors of the inverse
 * map at that point by finite differences, the smallest d found and the radius r* that
 * maximizes h (r) = ln (density of r) - d / 2 std_noise^2, d the squared distance of
 * eps_obs to the sheared ellipticity. h is scanned
 * on a uniform grid of radii merged with radii next to both ends and at |c| plus and minus
 * multiples of the widths std_noise J_1 and std_noise J_2, and r* is refined by golden
 * section between the neighbors of the best radius.
 */
static void
_nc_galaxy_shape_factor_spectral_seed (NcGalaxyShapeFactorSpectralPrivate * const self, const NcGalaxyShapeFactorSpectralKernel *k,
                                       const gdouble std_noise, const gdouble a, const gdouble b,
                                       NcGalaxyShapeFactorSpectralSeed *seed)
{
  static const gdouble offsets[NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_OFFSETS] = {0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 25.0};
  const complex double eps                                                   = k->e_1 + I * k->e_2;
  const gdouble abs_eps                                                      = cabs (eps);
  const complex double eps_0                                                 = (abs_eps < 1.0) ? eps : eps / abs_eps * (1.0 - 1.0e-12);
  const complex double c                                                     = self->map_inv (k->t, eps_0);
  const gdouble abs_c                                                        = cabs (c);
  const gdouble cos_0                                                        = (abs_c > 0.0) ? creal (c) / abs_c : 1.0;
  const gdouble sin_0                                                        = (abs_c > 0.0) ? cimag (c) / abs_c : 0.0;
  const gdouble dh                                                           = 1.0e-6;
  const guint n_uniform                                                      = NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_UNIFORM;
  gdouble extra[NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_EXTRA];
  gdouble scan_r[NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_UNIFORM + NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_EXTRA];
  gdouble scan_h[NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_UNIFORM + NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_EXTRA];
  gdouble scan_d[NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_UNIFORM + NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_EXTRA];
  gdouble scan_x[NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_UNIFORM + NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_EXTRA];
  guint n_scan     = 0;
  gdouble best_h   = GSL_NEGINF;
  gdouble r_before = 0.0;
  gdouble r_after  = 0.0;
  gdouble r_prev   = 0.0;
  gdouble shift    = GSL_POSINF;
  gboolean pending = FALSE;
  guint n_extra    = 0;
  guint iu         = 0;
  guint ie         = 0;
  guint i;

  seed->theta_0 = carg (c);

  if (cabs (eps_0) > 0.0)
  {
    seed->J_1 = cabs (self->map_inv (k->t, eps_0 * (1.0 - dh)) - c) / (dh * cabs (eps_0));
    seed->J_2 = cabs (self->map_inv (k->t, eps_0 * cexp (I * dh)) - c) / (dh * cabs (eps_0));
  }
  else
  {
    seed->J_1 = seed->J_2 = cabs (self->map_inv (k->t, dh) - c) / dh;
  }

  for (i = 2; i < 2 + NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_END_DECADES; i++)
  {
    extra[n_extra++] = pow (10.0, -(gdouble) i);
    extra[n_extra++] = 1.0 - pow (10.0, -(gdouble) i);
  }

  for (i = 0; i < NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_OFFSETS; i++)
  {
    extra[n_extra++] = abs_c + std_noise * seed->J_1 * offsets[i];
    extra[n_extra++] = abs_c - std_noise * seed->J_1 * offsets[i];
    extra[n_extra++] = abs_c + std_noise * seed->J_2 * offsets[i];
    extra[n_extra++] = abs_c - std_noise * seed->J_2 * offsets[i];
  }

  _nc_galaxy_shape_factor_spectral_sort (extra, n_extra);

  /* The uniform radii and the extra ones, in ascending order, each radius once */
  while ((iu < n_uniform) || (ie < n_extra))
  {
    const gdouble r_u = (iu < n_uniform) ? 1.0e-4 + (1.0 - 2.0e-4) * iu / (n_uniform - 1.0) : GSL_POSINF;
    const gdouble r_e = (ie < n_extra) ? extra[ie] : GSL_POSINF;
    const gdouble r   = GSL_MIN (r_u, r_e);

    if (r_u <= r_e)
      iu++;
    else
      ie++;

    if (!(r > 0.0) || !(r < 1.0) || ((n_scan > 0) && (r == scan_r[n_scan - 1])))
      continue;

    scan_r[n_scan++] = r;
  }

  _nc_galaxy_shape_factor_spectral_h (k, a, b, n_scan, scan_r, cos_0, sin_0, scan_h, scan_d, scan_x);

  for (i = 0; i < n_scan; i++)
  {
    const gdouble r = scan_r[i];
    const gdouble h = scan_h[i];

    shift = GSL_MIN (shift, scan_d[i]);

    if (pending)
    {
      r_after = r;
      pending = FALSE;
    }

    if (h > best_h)
    {
      best_h   = h;
      r_before = (r_prev > 0.0) ? r_prev : r;
      r_after  = r;
      pending  = TRUE;
    }

    r_prev = r;
  }

  /* The maximum of h on [r_before, r_after] by golden section */
  {
    const gdouble gr = 0.5 * (sqrt (5.0) - 1.0);
    gdouble lo       = r_before;
    gdouble hi       = r_after;
    gdouble x1       = hi - gr * (hi - lo);
    gdouble x2       = lo + gr * (hi - lo);
    gdouble d, xs, h1, h2;

    _nc_galaxy_shape_factor_spectral_h (k, a, b, 1, &x1, cos_0, sin_0, &h1, &d, &xs);
    _nc_galaxy_shape_factor_spectral_h (k, a, b, 1, &x2, cos_0, sin_0, &h2, &d, &xs);

    for (i = 0; i < NC_GALAXY_SHAPE_FACTOR_SPECTRAL_GOLDEN_STEPS; i++)
    {
      if (h1 > h2)
      {
        hi = x2;
        x2 = x1;
        h2 = h1;
        x1 = hi - gr * (hi - lo);
        _nc_galaxy_shape_factor_spectral_h (k, a, b, 1, &x1, cos_0, sin_0, &h1, &d, &xs);
      }
      else
      {
        lo = x1;
        x1 = x2;
        h1 = h2;
        x2 = lo + gr * (hi - lo);
        _nc_galaxy_shape_factor_spectral_h (k, a, b, 1, &x2, cos_0, sin_0, &h2, &d, &xs);
      }
    }

    seed->r_star = 0.5 * (lo + hi);
  }

  seed->shift = shift;
}

/* log1p (x) for x > -1, in a loop that vectorizes: Kahan's log (u) x / (u - 1), u = 1 + x */
static inline gdouble
_nc_galaxy_shape_factor_spectral_log1p (const gdouble x)
{
  const gdouble u = 1.0 + x;

  return ncm_vmath_select (u == 1.0, x, ncm_vmath_log (u) * (x / (u - 1.0)));
}

/*
 * h (r) = ln (density of r) - d / (2 std_noise^2) at the n radii r[i], 0 < r < 1, of the
 * ray of direction (cos_0, sin_0), in h[i], with d the squared distance of eps_obs to the
 * sheared ellipticity, in d[i]. The logarithm of the density is the two end power laws,
 * of exponents a and b, plus the logarithm of the density of the convention divided by
 * them (in x, scratch). Both pieces of the convention are computed at the n points at once.
 */
static void
_nc_galaxy_shape_factor_spectral_h (const NcGalaxyShapeFactorSpectralKernel *k, const gdouble a, const gdouble b, const guint n,
                                    const gdouble *r, const gdouble cos_0, const gdouble sin_0, gdouble *h, gdouble *d, gdouble *x)
{
  const gdouble c2 = k->c2;
  guint i;

  _nc_galaxy_shape_factor_spectral_dist2 (k, n, r, cos_0, sin_0, d);
  _nc_galaxy_shape_factor_spectral_density (k, n, r, x);

  #pragma omp simd

  for (i = 0; i < n; i++)
    h[i] = a * ncm_vmath_log (r[i]) + b * _nc_galaxy_shape_factor_spectral_log1p (-r[i]) + ncm_vmath_log (x[i]) - d[i] * c2;
}

/* Sorts x[0..n-1] in ascending order, by insertion */
static void
_nc_galaxy_shape_factor_spectral_sort (gdouble *x, const guint n)
{
  guint i;

  for (i = 1; i < n; i++)
  {
    const gdouble v = x[i];
    guint j         = i;

    while ((j > 0) && (x[j - 1] > v))
    {
      x[j] = x[j - 1];
      j--;
    }

    x[j] = v;
  }
}

/*
 * The map of the convention at the real shear t, in the variable of integration z = x + i y:
 * T (m_t (z)) for TRACE (_nc_galaxy_shape_factor_spectral_map_trace()) and m_t (z) for
 * TRACE_DET, m_t the TRACE_DET map. Its callers pass trace as a constant, so that each loop
 * that calls it holds a single map once inlined, and vectorizes.
 */
NC_GALAXY_SHAPE_FACTOR_SPECTRAL_INLINE void
_nc_galaxy_shape_factor_spectral_map (const gboolean trace, const gdouble t, const gdouble x, const gdouble y,
                                      gdouble * restrict f_x, gdouble * restrict f_y)
{
  if (trace)
    _nc_galaxy_shape_factor_spectral_map_trace (t, x, y, f_x, f_y);
  else
    nc_wl_ellipticity_apply_shear_trace_det_real (t, x, y, f_x, f_y);
}

/* The loop of _nc_galaxy_shape_factor_spectral_kernel() for one convention */
NC_GALAXY_SHAPE_FACTOR_SPECTRAL_INLINE void
_nc_galaxy_shape_factor_spectral_kernel_loop (const NcGalaxyShapeFactorSpectralKernel *k, const gboolean trace, const gdouble r, const guint n,
                                              const gdouble *cos_theta, const gdouble *sin_theta, gdouble *fval)
{
  const gdouble t        = k->t;
  const gdouble e_1      = k->e_1;
  const gdouble e_2      = k->e_2;
  const gdouble shift_c2 = k->shift_c2;
  const gdouble c2       = k->c2;
  guint j;

  /* The fields of k in locals: fval might alias them, which keeps the loop scalar */
  #pragma omp simd

  for (j = 0; j < n; j++)
  {
    gdouble f_x, f_y;

    _nc_galaxy_shape_factor_spectral_map (trace, t, r * cos_theta[j], r * sin_theta[j], &f_x, &f_y);
    fval[j] = ncm_vmath_exp (shift_c2 - (gsl_pow_2 (e_1 - f_x) + gsl_pow_2 (e_2 - f_y)) * c2);
  }
}

/*
 * The integrand of _ncm_spectral_polar_integrate_radial(), in the polar coordinates of the
 * variable of integration z: G (r, theta) = (reduced density of r = |z|)
 * N_2 (eps_obs - F_t (z)) / (2 pi), F_t the map of the convention
 * (_nc_galaxy_shape_factor_spectral_map()), with the offset of the kernel taken out of the
 * exponent of N_2. The factor that depends only on r, the reduced density times norm, is
 * _nc_galaxy_shape_factor_spectral_radial(): here only the exponential.
 */
static void
_nc_galaxy_shape_factor_spectral_kernel (gpointer user_data, gdouble r, guint n, const gdouble *theta,
                                         const gdouble *cos_theta, const gdouble *sin_theta, gdouble *fval)
{
  const NcGalaxyShapeFactorSpectralKernel *k = user_data;

  if (k->trace)
    _nc_galaxy_shape_factor_spectral_kernel_loop (k, TRUE, r, n, cos_theta, sin_theta, fval);
  else
    _nc_galaxy_shape_factor_spectral_kernel_loop (k, FALSE, r, n, cos_theta, sin_theta, fval);
}

/*
 * The factor of the integrand that depends only on r, given to
 * _ncm_spectral_polar_integrate_radial() so that it is computed once per radius: norm times
 * the density of the convention (_nc_galaxy_shape_factor_spectral_density()).
 */
static void
_nc_galaxy_shape_factor_spectral_radial (gpointer user_data, guint n, const gdouble *r, gdouble *rho)
{
  const NcGalaxyShapeFactorSpectralKernel *k = user_data;
  guint i;

  _nc_galaxy_shape_factor_spectral_density (k, n, r, rho);

  for (i = 0; i < n; i++)
    rho[i] *= k->norm;
}

/*
 * The TRACE map at the real shear t in the variable e_I = x + i y: T (m_t (e_I)), with
 * m_t the TRACE_DET map. It equals f_t (T (e_I)), f_t the TRACE map. With
 * m_t (e_I) = N / D, N = (x + t) + i y and D = (1 + t x) + i t y, and T (m) = 2 m / (1 + |m|^2),
 * it is 2 N conj (D) / (|D|^2 + |N|^2): one division instead of the two of m_t followed by
 * T, the slowest operations of the loop of the kernel.
 */
static void
_nc_galaxy_shape_factor_spectral_map_trace (gdouble t, gdouble x, gdouble y, gdouble * restrict f_x, gdouble * restrict f_y)
{
  const gdouble num_x = x + t;
  const gdouble den_x = 1.0 + t * x;
  const gdouble den_y = t * y;
  const gdouble fac   = 2.0 / ((num_x * num_x + y * y) + (den_x * den_x + den_y * den_y));

  *f_x = (num_x * den_x + y * den_y) * fac;
  *f_y = (y * den_x - num_x * den_y) * fac;
}

/* The inverse of the map above at the shear g: m_g^{-1} (T^{-1} (e_obs)), for |e_obs| <= 1 */
static complex double
_nc_galaxy_shape_factor_spectral_map_inv_trace (complex double g, complex double e_obs)
{
  const gdouble abs_e2 = gsl_pow_2 (creal (e_obs)) + gsl_pow_2 (cimag (e_obs));

  return nc_wl_ellipticity_apply_shear_inv_trace_det (g, e_obs / (1.0 + sqrt (1.0 - abs_e2)));
}

/*
 * The density of r = |z| divided by its two end power laws, at the n radii r[i] in out[i].
 * TRACE_DET: the reduced density of the population at r = |chi_I|. TRACE: the density of
 * r = |e_I| divided by r^alpha_o (1 - r)^(2 alpha_1 + 1), for a population that is a
 * density of |chi_I| = 2 r / (1 + r^2) with the exponents alpha_o and alpha_1,
 * 2^(alpha_o + 1) (1 + r) (1 + r^2)^-(alpha_o + alpha_1 + 2) times the reduced density of
 * the population at |chi_I|, the factor of the convention applied in a loop that
 * vectorizes.
 */
static void
_nc_galaxy_shape_factor_spectral_density (const NcGalaxyShapeFactorSpectralKernel *k, guint n, const gdouble *r, gdouble *out)
{
  guint i;

  if (k->trace)
  {
    const gdouble c  = k->dens_norm;
    const gdouble mp = -k->dens_pow;

    for (i = 0; i < n; i++)
      out[i] = nc_galaxy_shape_pop_eval_p_reduced (k->pop, k->pop_data, 2.0 * r[i] / (1.0 + r[i] * r[i]));

    #pragma omp simd

    for (i = 0; i < n; i++)
      out[i] *= c * (1.0 + r[i]) * ncm_vmath_pow_pos (1.0 + r[i] * r[i], mp);
  }
  else
  {
    for (i = 0; i < n; i++)
      out[i] = nc_galaxy_shape_pop_eval_p_reduced (k->pop, k->pop_data, r[i]);
  }
}

/* The loop of _nc_galaxy_shape_factor_spectral_dist2() for one convention */
NC_GALAXY_SHAPE_FACTOR_SPECTRAL_INLINE void
_nc_galaxy_shape_factor_spectral_dist2_loop (const NcGalaxyShapeFactorSpectralKernel *k, const gboolean trace, const guint n, const gdouble *r,
                                             const gdouble cos_0, const gdouble sin_0, gdouble *d)
{
  const gdouble t   = k->t;
  const gdouble e_1 = k->e_1;
  const gdouble e_2 = k->e_2;
  guint i;

  #pragma omp simd

  for (i = 0; i < n; i++)
  {
    gdouble f_x, f_y;

    _nc_galaxy_shape_factor_spectral_map (trace, t, r[i] * cos_0, r[i] * sin_0, &f_x, &f_y);
    d[i] = gsl_pow_2 (e_1 - f_x) + gsl_pow_2 (e_2 - f_y);
  }
}

/*
 * The squared distance of eps_obs to the map of the convention of z = r (cos_0, sin_0)
 * (_nc_galaxy_shape_factor_spectral_map()), at the n radii r[i] in d[i]
 */
static void
_nc_galaxy_shape_factor_spectral_dist2 (const NcGalaxyShapeFactorSpectralKernel *k, guint n, const gdouble *r,
                                        gdouble cos_0, gdouble sin_0, gdouble *d)
{
  if (k->trace)
    _nc_galaxy_shape_factor_spectral_dist2_loop (k, TRUE, n, r, cos_0, sin_0, d);
  else
    _nc_galaxy_shape_factor_spectral_dist2_loop (k, FALSE, n, r, cos_0, sin_0, d);
}

/**
 * nc_galaxy_shape_factor_spectral_new:
 * @ellip_conv: a #NcGalaxyWLObsEllipConv
 *
 * Creates a new #NcGalaxyShapeFactorSpectral with the ellipticity convention @ellip_conv.
 *
 * Returns: (transfer full): a new #NcGalaxyShapeFactorSpectral.
 */
NcGalaxyShapeFactorSpectral *
nc_galaxy_shape_factor_spectral_new (NcGalaxyWLObsEllipConv ellip_conv)
{
  return g_object_new (NC_TYPE_GALAXY_SHAPE_FACTOR_SPECTRAL,
                       "ellip-conv", ellip_conv,
                       NULL);
}

/**
 * nc_galaxy_shape_factor_spectral_ref:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Increases the reference count of @gsfs by one.
 *
 * Returns: (transfer full): @gsfs.
 */
NcGalaxyShapeFactorSpectral *
nc_galaxy_shape_factor_spectral_ref (NcGalaxyShapeFactorSpectral *gsfs)
{
  return g_object_ref (gsfs);
}

/**
 * nc_galaxy_shape_factor_spectral_free:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Decreases the reference count of @gsfs by one.
 */
void
nc_galaxy_shape_factor_spectral_free (NcGalaxyShapeFactorSpectral *gsfs)
{
  g_object_unref (gsfs);
}

/**
 * nc_galaxy_shape_factor_spectral_clear:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Decreases the reference count of *@gsfs by one and sets *@gsfs to %NULL.
 */
void
nc_galaxy_shape_factor_spectral_clear (NcGalaxyShapeFactorSpectral **gsfs)
{
  g_clear_object (gsfs);
}

/**
 * nc_galaxy_shape_factor_spectral_set_reltol:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @reltol: relative tolerance of the integral
 *
 * Sets #NcGalaxyShapeFactorSpectral:reltol.
 */
void
nc_galaxy_shape_factor_spectral_set_reltol (NcGalaxyShapeFactorSpectral *gsfs, const gdouble reltol)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  self->reltol = reltol;
}

/**
 * nc_galaxy_shape_factor_spectral_get_reltol:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Returns: #NcGalaxyShapeFactorSpectral:reltol.
 */
gdouble
nc_galaxy_shape_factor_spectral_get_reltol (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  return self->reltol;
}

/**
 * nc_galaxy_shape_factor_spectral_set_obs_source:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @pop: the #NcGalaxyShapePop that built the interpolants
 * @obs: a #NcGalaxyWLObs
 *
 * Sets the record of the interpolants written into @obs by
 * nc_galaxy_shape_factor_data_write_row() (nc_galaxy_wl_obs_set_shape_interp_source()):
 * copies of @gsfs under `shape-factor` and of @pop with its parameter values under
 * `shape-pop`. Call it after the last row is written. A fit that reads the rows compares its
 * shape factor and population with the record in nc_galaxy_shape_factor_data_prepare().
 */
void
nc_galaxy_shape_factor_spectral_set_obs_source (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyWLObs *obs)
{
  NcmSerialize *ser     = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcmObjDictStr *source = ncm_obj_dict_str_new ();
  GObject *gsfs_dup     = ncm_serialize_dup_obj (ser, G_OBJECT (gsfs));
  GObject *pop_dup;

  ncm_serialize_reset (ser, FALSE);
  pop_dup = ncm_serialize_dup_obj (ser, G_OBJECT (pop));

  ncm_obj_dict_str_add (source, "shape-factor", gsfs_dup);
  ncm_obj_dict_str_add (source, "shape-pop", pop_dup);
  nc_galaxy_wl_obs_set_shape_interp_source (obs, source);

  g_object_unref (gsfs_dup);
  g_object_unref (pop_dup);
  ncm_obj_dict_str_unref (source);
  ncm_serialize_free (ser);
}

/**
 * nc_galaxy_shape_factor_spectral_set_use_interp:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @use_interp: whether $\ln P$ is interpolated in $g$ for each galaxy
 *
 * Sets #NcGalaxyShapeFactorSpectral:use-interp.
 */
void
nc_galaxy_shape_factor_spectral_set_use_interp (NcGalaxyShapeFactorSpectral *gsfs, const gboolean use_interp)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  self->use_interp = use_interp;
}

/**
 * nc_galaxy_shape_factor_spectral_get_use_interp:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Returns: #NcGalaxyShapeFactorSpectral:use-interp.
 */
gboolean
nc_galaxy_shape_factor_spectral_get_use_interp (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  return self->use_interp;
}

/**
 * nc_galaxy_shape_factor_spectral_set_interp_g_max:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @g_max: largest $|g|$ computed from the interpolant, in $(0, 1/3)$
 *
 * Sets #NcGalaxyShapeFactorSpectral:interp-g-max.
 */
void
nc_galaxy_shape_factor_spectral_set_interp_g_max (NcGalaxyShapeFactorSpectral *gsfs, const gdouble g_max)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  g_assert_cmpfloat (g_max, >, 0.0);
  g_assert_cmpfloat (3.0 * g_max, <, 1.0);

  self->interp_g_max = g_max;
}

/**
 * nc_galaxy_shape_factor_spectral_get_interp_g_max:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Returns: #NcGalaxyShapeFactorSpectral:interp-g-max.
 */
gdouble
nc_galaxy_shape_factor_spectral_get_interp_g_max (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  return self->interp_g_max;
}

/**
 * nc_galaxy_shape_factor_spectral_set_interp_abstol:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @abstol: absolute tolerance on $\ln P$ of the interpolant
 *
 * Sets #NcGalaxyShapeFactorSpectral:interp-abstol.
 */
void
nc_galaxy_shape_factor_spectral_set_interp_abstol (NcGalaxyShapeFactorSpectral *gsfs, const gdouble abstol)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  g_assert_cmpfloat (abstol, >, 0.0);

  self->interp_abstol = abstol;
}

/**
 * nc_galaxy_shape_factor_spectral_get_interp_abstol:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Returns: #NcGalaxyShapeFactorSpectral:interp-abstol.
 */
gdouble
nc_galaxy_shape_factor_spectral_get_interp_abstol (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  return self->interp_abstol;
}

/**
 * nc_galaxy_shape_factor_spectral_set_interp_min_radial_panels:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @n: minimum number of radial intervals of the interpolant
 *
 * Sets #NcGalaxyShapeFactorSpectral:interp-min-radial-panels, which must be at least 1.
 */
void
nc_galaxy_shape_factor_spectral_set_interp_min_radial_panels (NcGalaxyShapeFactorSpectral *gsfs, const guint n)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  ncm_spectral_polar_set_min_radial_panels (self->polar_interp, n);

  self->interp_min_radial_panels = n;
}

/**
 * nc_galaxy_shape_factor_spectral_get_interp_min_radial_panels:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Returns: #NcGalaxyShapeFactorSpectral:interp-min-radial-panels.
 */
guint
nc_galaxy_shape_factor_spectral_get_interp_min_radial_panels (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  return self->interp_min_radial_panels;
}

/**
 * nc_galaxy_shape_factor_spectral_set_line_abstol:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @abstol: the absolute tolerance on $\ln P$ in [-1, 1], zero to integrate every node,
 *   negative for #NcGalaxyShapeFactorSpectral:reltol / 10
 *
 * Sets #NcGalaxyShapeFactorSpectral:line-abstol.
 */
void
nc_galaxy_shape_factor_spectral_set_line_abstol (NcGalaxyShapeFactorSpectral *gsfs, const gdouble abstol)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  g_assert_cmpfloat (abstol, >=, -1.0);
  g_assert_cmpfloat (abstol, <=, 1.0);

  self->line_abstol = abstol;
}

/**
 * nc_galaxy_shape_factor_spectral_get_line_abstol:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Returns: #NcGalaxyShapeFactorSpectral:line-abstol.
 */
gdouble
nc_galaxy_shape_factor_spectral_get_line_abstol (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  return self->line_abstol;
}

/**
 * nc_galaxy_shape_factor_spectral_set_restrict_range:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @restrict_range: whether to build each interpolant only over the reachable shears
 *
 * Sets #NcGalaxyShapeFactorSpectral:restrict-range. It applies from the next
 * nc_galaxy_shape_factor_prepare() and nc_galaxy_shape_factor_data_prepare().
 */
void
nc_galaxy_shape_factor_spectral_set_restrict_range (NcGalaxyShapeFactorSpectral *gsfs, const gboolean restrict_range)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  self->restrict_range = restrict_range;
}

/**
 * nc_galaxy_shape_factor_spectral_get_restrict_range:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Returns: #NcGalaxyShapeFactorSpectral:restrict-range.
 */
gboolean
nc_galaxy_shape_factor_spectral_get_restrict_range (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  return self->restrict_range;
}

/**
 * nc_galaxy_shape_factor_spectral_set_pop_order:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @pop_order: the number of nodes per free population parameter, 0 or 1 for none
 *
 * Sets #NcGalaxyShapeFactorSpectral:pop-order. It applies from the next
 * nc_galaxy_shape_factor_prepare(); the interpolants are built again for the new box.
 */
void
nc_galaxy_shape_factor_spectral_set_pop_order (NcGalaxyShapeFactorSpectral *gsfs, const guint pop_order)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  self->pop_order = pop_order;
}

/**
 * nc_galaxy_shape_factor_spectral_get_pop_order:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Returns: #NcGalaxyShapeFactorSpectral:pop-order.
 */
guint
nc_galaxy_shape_factor_spectral_get_pop_order (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  return self->pop_order;
}

/**
 * nc_galaxy_shape_factor_spectral_get_pop_n_unconverged:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * The number of expansions over the population parameters built by @gsfs that did not
 * converge (see #NcGalaxyShapeFactorSpectral:pop-order) since its box was last set: those
 * galaxies are integrated. A galaxy is counted at each build, so a galaxy built again for
 * another observed ellipticity or noise is counted again.
 *
 * Returns: the number of expansions that did not converge.
 */
guint
nc_galaxy_shape_factor_spectral_get_pop_n_unconverged (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  return (guint) g_atomic_int_get (&self->pop_n_unconverged);
}

/**
 * nc_galaxy_shape_factor_spectral_get_n_loaded_short:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * The number of galaxies whose interpolant, loaded from a catalog, does not cover the shears
 * they reach in the box of the fit (see #NcGalaxyShapeFactorSpectral:restrict-range): the
 * catalog was prepared for a smaller box, and the shears beyond it are integrated. A galaxy
 * is counted each time its interpolant is compared with the fit.
 *
 * Returns: the number of loaded interpolants short of the box of the fit.
 */
guint
nc_galaxy_shape_factor_spectral_get_n_loaded_short (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  return (guint) g_atomic_int_get (&self->n_short);
}

/**
 * nc_galaxy_shape_factor_spectral_get_fallback_stats:
 * @reason: a #NcGalaxyShapeFactorSpectralFallback
 * @lines: (out) (optional): the lines of shears tabulated
 * @shears: (out) (optional): the shears integrated one by one
 * @seconds: (out) (optional): the time spent on them, in seconds
 *
 * The evaluations of $\ln P$ that could not use the interpolant of their galaxy, with
 * #NcGalaxyShapeFactorSpectral:use-interp set, for @reason, since the process started or
 * nc_galaxy_shape_factor_spectral_reset_fallback_stats() was called. A galaxy that cannot use
 * its interpolant anywhere on the line of its shears (nc_galaxy_shape_factor_eval_line()) is
 * tabulated along it, from integrals on at most half as many points as it has shears, as
 * without use-interp (#NcGalaxyShapeFactorSpectral:line-abstol); any other shear the
 * interpolant does not cover is integrated. The counts are summed over every
 * #NcGalaxyShapeFactorSpectral of the process, so that they include each copy of a likelihood
 * a sampler makes.
 */
void
nc_galaxy_shape_factor_spectral_get_fallback_stats (NcGalaxyShapeFactorSpectralFallback reason, guint64 *lines, guint64 *shears, gdouble *seconds)
{
  NcGalaxyShapeFactorSpectralFallbackStats *st;

  g_assert_cmpint (reason, <, NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_LEN);
  st = &_nc_galaxy_shape_factor_spectral_fallback_stats[reason];

  if (lines != NULL)
    *lines = GPOINTER_TO_SIZE (g_atomic_pointer_get (&st->lines));

  if (shears != NULL)
    *shears = GPOINTER_TO_SIZE (g_atomic_pointer_get (&st->shears));

  if (seconds != NULL)
    *seconds = 1.0e-6 * GPOINTER_TO_SIZE (g_atomic_pointer_get (&st->usec));
}

/**
 * nc_galaxy_shape_factor_spectral_reset_fallback_stats:
 *
 * Sets the counts of nc_galaxy_shape_factor_spectral_get_fallback_stats() to zero, for
 * every reason, and lets the next fallback outside the coverage of the interpolants be
 * reported again.
 */
void
nc_galaxy_shape_factor_spectral_reset_fallback_stats (void)
{
  guint i;

  for (i = 0; i < NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_LEN; i++)
  {
    g_atomic_pointer_and (&_nc_galaxy_shape_factor_spectral_fallback_stats[i].lines, 0);
    g_atomic_pointer_and (&_nc_galaxy_shape_factor_spectral_fallback_stats[i].shears, 0);
    g_atomic_pointer_and (&_nc_galaxy_shape_factor_spectral_fallback_stats[i].usec, 0);
  }

  g_atomic_int_set (&_nc_galaxy_shape_factor_spectral_fallback_noted, 0);
}

/**
 * nc_galaxy_shape_factor_spectral_eval_ln_direct:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @pop: a #NcGalaxyShapePop
 * @pop_data: a resolved #NcGalaxyShapePopData of @pop
 * @g_1: reduced shear, real component
 * @g_2: reduced shear, imaginary component
 * @epsilon_obs_1: observed ellipticity, real component
 * @epsilon_obs_2: observed ellipticity, imaginary component
 * @std_noise: measurement noise standard deviation
 *
 * Computes $\ln P(\epsilon_\mathrm{obs} \mid g)$ from the population and the values of
 * one galaxy, without a #NcGalaxyShapeFactorData. @pop_data must have been prepared with
 * nc_galaxy_shape_pop_prepare() since the parameters of @pop last changed.
 * nc_galaxy_shape_factor_eval_ln_marginal() computes the same value for the galaxy of a
 * #NcGalaxyShapeFactorData.
 *
 * Returns: $\ln P(\epsilon_\mathrm{obs} \mid g)$.
 */
gdouble
nc_galaxy_shape_factor_spectral_eval_ln_direct (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data,
                                                const gdouble g_1, const gdouble g_2,
                                                const gdouble epsilon_obs_1, const gdouble epsilon_obs_2,
                                                const gdouble std_noise)
{
  return _nc_galaxy_shape_factor_spectral_eval_ln (gsfs, pop, pop_data, std_noise, g_1, g_2, epsilon_obs_1, epsilon_obs_2, NULL);
}

/**
 * nc_galaxy_shape_factor_spectral_eval_panels:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @pop: a #NcGalaxyShapePop
 * @pop_data: a resolved #NcGalaxyShapePopData of @pop
 * @g_1: reduced shear, real component
 * @g_2: reduced shear, imaginary component
 * @epsilon_obs_1: observed ellipticity, real component
 * @epsilon_obs_2: observed ellipticity, imaginary component
 * @std_noise: measurement noise standard deviation
 * @ln_marginal: (out) (optional): $\ln P(\epsilon_\mathrm{obs} \mid g)$
 *
 * Computes $\ln P(\epsilon_\mathrm{obs} \mid g)$ as
 * nc_galaxy_shape_factor_spectral_eval_ln_direct() and returns the panels the integral
 * was accepted on, for inspection: one row per panel with the columns of
 * ncm_spectral_polar_integrate(), $r_a$, $r_b$, $\theta_a$, $\theta_b$, the two levels,
 * the value and the change of the panel. $(r, \theta)$ are the polar coordinates of the
 * variable of integration in the frame of @g_1, @g_2 and of the observed ellipticity:
 * $\chi_I$ for %NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, and for
 * %NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE the variable with $|\chi_I| = 2 r / (1 + r^2)$ and
 * the angle of $\chi_I$. For $1 - \min(|g|, 1 / |g|) < 10^{-5}$ the panels are those of
 * the integral at $1 - 10^{-5}$. For $|g| > 1$ the
 * panels are those of the integral at $g / |g|^2$, and at $|g| = 1$ there is no integral
 * and the function returns %NULL.
 *
 * Returns: (transfer full) (nullable): the accepted panels.
 */
NcmMatrix *
nc_galaxy_shape_factor_spectral_eval_panels (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data,
                                             const gdouble g_1, const gdouble g_2,
                                             const gdouble epsilon_obs_1, const gdouble epsilon_obs_2,
                                             const gdouble std_noise, gdouble *ln_marginal)
{
  NcmMatrix *panels  = NULL;
  const gdouble ln_P = _nc_galaxy_shape_factor_spectral_eval_ln (gsfs, pop, pop_data, std_noise, g_1, g_2, epsilon_obs_1, epsilon_obs_2, &panels);

  if (ln_marginal != NULL)
    *ln_marginal = ln_P;

  return panels;
}

/* The argument of the kernels of one multi-population integral: k first, so that the kernel reads it */
typedef struct _NcGalaxyShapeFactorSpectralMultiArg
{
  NcGalaxyShapeFactorSpectralKernel k;
  guint n_w;
  NcGalaxyShapePop **pops;
  NcGalaxyShapePopData **pop_data;
  gdouble *dens_norm;
  gdouble *dens_pow;
} NcGalaxyShapeFactorSpectralMultiArg;

/* The exponents of pop in the variable of integration of self, and the factors of the TRACE density */
static void
_nc_galaxy_shape_factor_spectral_pop_exponents (NcGalaxyShapeFactorSpectralPrivate * const self, NcGalaxyShapePop *pop,
                                                gdouble *a, gdouble *b, gdouble *dens_norm, gdouble *dens_pow)
{
  const gdouble b_pop = nc_galaxy_shape_pop_exponent_at_one (pop);

  *a         = nc_galaxy_shape_pop_exponent_at_origin (pop);
  *b         = self->axis_ratio_var ? 2.0 * b_pop + 1.0 : b_pop;
  *dens_norm = exp2 (*a + 1.0);
  *dens_pow  = *a + b_pop + 2.0;
}

/* The node set of layout_pop and the n_w populations pops, each prepared here */
static NcGalaxyShapeFactorSpectralNodes *
_nc_galaxy_shape_factor_spectral_nodes_new (NcGalaxyShapeFactorSpectralPrivate * const self, NcGalaxyShapePop *layout_pop,
                                            const guint n_w, NcGalaxyShapePop **pops)
{
  NcGalaxyShapeFactorSpectralNodes *nodes = g_new0 (NcGalaxyShapeFactorSpectralNodes, 1);
  guint q;

  nodes->layout_pop  = nc_galaxy_shape_pop_ref (layout_pop);
  nodes->layout_data = nc_galaxy_shape_pop_data_new (layout_pop);
  nc_galaxy_shape_pop_prepare (layout_pop, nodes->layout_data);
  _nc_galaxy_shape_factor_spectral_pop_exponents (self, layout_pop, &nodes->layout_a, &nodes->layout_b,
                                                  &nodes->layout_dens_norm, &nodes->layout_dens_pow);

  nodes->n_w       = n_w;
  nodes->pops      = g_new (NcGalaxyShapePop *, n_w);
  nodes->pop_data  = g_new (NcGalaxyShapePopData *, n_w);
  nodes->a         = g_new (gdouble, n_w);
  nodes->b         = g_new (gdouble, n_w);
  nodes->dens_norm = g_new (gdouble, n_w);
  nodes->dens_pow  = g_new (gdouble, n_w);

  for (q = 0; q < n_w; q++)
  {
    nodes->pops[q]     = nc_galaxy_shape_pop_ref (pops[q]);
    nodes->pop_data[q] = nc_galaxy_shape_pop_data_new (pops[q]);
    nc_galaxy_shape_pop_prepare (pops[q], nodes->pop_data[q]);
    _nc_galaxy_shape_factor_spectral_pop_exponents (self, pops[q], &nodes->a[q], &nodes->b[q], &nodes->dens_norm[q], &nodes->dens_pow[q]);
  }

  nodes->multi = _ncm_spectral_polar_multi_new (ncm_spectral_polar_get_level_cap (self->polar), n_w, nodes->a, nodes->b);

  return nodes;
}

static void
_nc_galaxy_shape_factor_spectral_nodes_free (NcGalaxyShapeFactorSpectralNodes *nodes)
{
  guint q;

  for (q = 0; q < nodes->n_w; q++)
  {
    nc_galaxy_shape_pop_free (nodes->pops[q]);
    nc_galaxy_shape_pop_data_unref (nodes->pop_data[q]);
  }

  nc_galaxy_shape_pop_free (nodes->layout_pop);
  nc_galaxy_shape_pop_data_unref (nodes->layout_data);
  _ncm_spectral_polar_multi_free (nodes->multi);
  g_free (nodes->pops);
  g_free (nodes->pop_data);
  g_free (nodes->a);
  g_free (nodes->b);
  g_free (nodes->dens_norm);
  g_free (nodes->dens_pow);
  g_free (nodes);
}

/*
 * The radial factor of the integrand for each population: norm times the density of the
 * convention (_nc_galaxy_shape_factor_spectral_density()) with the population q and its
 * factors, in w[q n + i].
 */
static void
_nc_galaxy_shape_factor_spectral_radial_multi (gpointer user_data, guint n, const gdouble *r, guint n_w, gdouble *w)
{
  const NcGalaxyShapeFactorSpectralMultiArg *arg = user_data;
  guint i, q;

  for (q = 0; q < n_w; q++)
  {
    NcGalaxyShapeFactorSpectralKernel kq = arg->k;
    gdouble *w_q                         = w + (gsize) q * n;

    kq.pop       = arg->pops[q];
    kq.pop_data  = arg->pop_data[q];
    kq.dens_norm = arg->dens_norm[q];
    kq.dens_pow  = arg->dens_pow[q];

    _nc_galaxy_shape_factor_spectral_density (&kq, n, r, w_q);

    for (i = 0; i < n; i++)
      w_q[i] *= kq.norm;
  }
}

/*
 * ln P at the shear (g_1, g_2) for each population of nodes, in out[0..n_w - 1]: the panels
 * of the adaptive integral of the layout population, and every population integrated on
 * them from the samples of the integrand that integral takes
 * (_ncm_spectral_polar_integrate_radial_multi()). Within 1e-5 of |g| = 1 each value is
 * _nc_galaxy_shape_factor_spectral_eval_ln(). With layout not NULL, also ln P of the layout
 * population, as _nc_galaxy_shape_factor_spectral_eval_ln() computes it.
 */
static void
_nc_galaxy_shape_factor_spectral_eval_ln_nodes (NcGalaxyShapeFactorSpectral *gsfs, const NcGalaxyShapeFactorSpectralNodes *nodes,
                                                const gdouble std_noise, const gdouble g_1, const gdouble g_2,
                                                const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, gdouble *out, gdouble *layout)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);
  const guint n_w                                 = nodes->n_w;
  const gdouble noise_var                         = gsl_pow_2 (std_noise);
  const gdouble abs_g2                            = g_1 * g_1 + g_2 * g_2;
  const gdouble abs_g                             = sqrt (abs_g2);
  const gdouble u                                 = fabs (_nc_galaxy_shape_factor_spectral_one_minus_abs2 (g_1, g_2));
  const gdouble delta                             = (abs_g2 > 1.0) ? u / (abs_g * (1.0 + abs_g)) : u / (1.0 + abs_g);
  guint q;

  if ((abs_g2 == 1.0) || (delta < NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NEAR_UNIT_DELTA))
  {
    for (q = 0; q < n_w; q++)
      out[q] = _nc_galaxy_shape_factor_spectral_eval_ln (gsfs, nodes->pops[q], nodes->pop_data[q], std_noise, g_1, g_2, epsilon_obs_1, epsilon_obs_2, NULL);

    if (layout != NULL)
      *layout = _nc_galaxy_shape_factor_spectral_eval_ln (gsfs, nodes->layout_pop, nodes->layout_data, std_noise, g_1, g_2, epsilon_obs_1, epsilon_obs_2, NULL);
  }
  else
  {
    const gdouble phase                   = (abs_g2 > 0.0) ? atan2 (g_2, g_1) : 0.0;
    const gdouble cos_p                   = cos (phase);
    const gdouble sin_p                   = sin (phase);
    NcGalaxyShapeFactorSpectralWork **ptr = ncm_memory_pool_get (self->work_pool);
    NcGalaxyShapeFactorSpectralWork *work = *ptr;
    NcGalaxyShapeFactorSpectralMultiArg arg;
    NcGalaxyShapeFactorSpectralSeed seed;
    GError *error = NULL;
    gdouble L, integral;

    arg.k = (NcGalaxyShapeFactorSpectralKernel) {
      .pop       = nodes->layout_pop,
      .pop_data  = nodes->layout_data,
      .t         = (abs_g2 > 1.0) ? 1.0 / abs_g : abs_g,
      .e_1       = epsilon_obs_1 * cos_p + epsilon_obs_2 * sin_p,
      .e_2       = epsilon_obs_2 * cos_p - epsilon_obs_1 * sin_p,
      .c2        = 0.5 / noise_var,
      .norm      = 1.0 / (2.0 * M_PI * noise_var) / (2.0 * M_PI),
      .shift_c2  = 0.0,
      .dens_norm = nodes->layout_dens_norm,
      .dens_pow  = nodes->layout_dens_pow,
      .trace     = self->axis_ratio_var
    };
    arg.n_w       = n_w;
    arg.pops      = nodes->pops;
    arg.pop_data  = nodes->pop_data;
    arg.dens_norm = nodes->dens_norm;
    arg.dens_pow  = nodes->dens_pow;

    /*
     * The adaptive integral of the layout population, as _nc_galaxy_shape_factor_spectral_integrate_ln(),
     * and on each panel it accepts the integrals of every population from the same samples
     */
    _nc_galaxy_shape_factor_spectral_seed (self, &arg.k, std_noise, nodes->layout_a, nodes->layout_b, &seed);
    arg.k.shift_c2 = seed.shift * arg.k.c2;
    L              = std_noise * GSL_MIN (seed.J_1, seed.J_2);
    ncm_spectral_polar_peak_edges (self->polar, 1.0, seed.r_star, seed.theta_0, L, L / GSL_MAX (seed.r_star, GSL_DBL_MIN),
                                   work->r_edges, work->theta_edges);
    integral = _ncm_spectral_polar_integrate_radial_multi (self->polar, &_nc_galaxy_shape_factor_spectral_kernel, &_nc_galaxy_shape_factor_spectral_radial, &arg.k, 1.0, nodes->layout_a, nodes->layout_b,
                                                           work->r_edges, work->theta_edges, self->reltol, 0.0, nodes->multi,
                                                           &_nc_galaxy_shape_factor_spectral_radial_multi, out, NULL, &error);
    ncm_memory_pool_return (ptr);

    if (error != NULL)
      g_error ("nc_galaxy_shape_factor_spectral: the integral of the layout population: %s, std_noise = % 22.15g, "
               "g = (% 22.15g, % 22.15g), eps_obs = (% 22.15g, % 22.15g).",
               error->message, std_noise, g_1, g_2, epsilon_obs_1, epsilon_obs_2);

    for (q = 0; q < n_w; q++)
      out[q] = log (out[q]) - arg.k.shift_c2;

    if (layout != NULL)
      *layout = log (integral) - arg.k.shift_c2;
  }
}

