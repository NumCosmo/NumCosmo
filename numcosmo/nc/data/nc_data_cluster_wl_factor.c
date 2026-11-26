/***************************************************************************
 *            nc_data_cluster_wl_factor.c
 *
 *  Sun Jul 5 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_data_cluster_wl_factor.c
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
 * NcDataClusterWLFactor:
 *
 * Cluster weak-lensing likelihood built on the per-galaxy Factor
 * calculators (#NcGalaxyPositionFactor / #NcGalaxyRedshiftFactor /
 * #NcGalaxyShapeFactor).
 *
 * The per-galaxy shape factor is
 * $$P(\epsilon_\mathrm{obs} \mid g) = \int_{|\chi_I|<1} \mathrm{d}^2\chi_I\,
 *   P_\mathrm{pop}(\chi_I)\, N_2\big(\epsilon_\mathrm{obs} - f_g(\chi_I);
 *   \sigma_\mathrm{noise}^2\big),$$
 * combined with the position and redshift factors and integrated over source
 * redshift.
 *
 * The redshift integral is done once, by this class, against every other
 * $z$-dependent factor. This follows #NcGalaxyRedshiftFactor's contract: the
 * calculator never integrates $z$ itself.
 *
 * Three methods are available: `CHEB_WEIGHTS`, the default, nested Chebyshev-Lobatto
 * levels whose weights integrate $p(z)$ exactly (nc_galaxy_redshift_factor_cheb_moments()),
 * with the level of each galaxy chosen at each evaluation, see
 * #NcDataClusterWLFactor:cheb-tol and #NcDataClusterWLFactor:cheb-level-min; `LNINT`,
 * an adaptive log-domain 1D integral; and `CUBATURE`, an adaptive #NcmIntegralND over the
 * linear-domain product. All three converge to the same integral and differ only in
 * numerical strategy and cost, and all three support bootstrap resampling.
 *
 * Under `CHEB_WEIGHTS` the nodes and weights of every galaxy form an
 * #NcGalaxyWLObsChebGrid. The dataset reads the one stored with the catalog
 * (#NcGalaxyWLObs:cheb-grid, see nc_data_cluster_wl_factor_write_cheb_grid()) when it
 * describes the current rows, so the copies a sampler makes of a dataset on a shared catalog
 * read one grid; otherwise it builds its own, in parallel, at prepare.
 *
 * nc_data_cluster_wl_factor_data_prepare() precomputes the galaxies in
 * parallel with OpenMP (thread count from `OMP_NUM_THREADS`), including the
 * redshift grid, for every position, redshift and shape factor. The
 * per-cycle prepare() of a fit runs inside the sampler's own threads and is
 * serial.
 *
 * No `r_min`/`r_max` weighting is applied at fit time, because applying it
 * biases the mass estimate.
 *
 * All preparation happens in prepare(): resolving models, refreshing each
 * Factor's caches, and re-preparing the $z$-integrands. m2lnL_val() and
 * resample() only use already-prepared state. The three Factor objects are
 * held as construct-only references, since calculators are designed to be
 * shared across likelihoods. This class never computes a Factor's
 * change-detection hash, never calls ncm_model_state_get_pkey(), and never
 * knows which #NcmModel a Factor depends on; it only compares the opaque
 * #guint64 values each Factor's prepare() returns. See the Factor classes for
 * that get_hash() and update_data() contract.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/data/nc_data_cluster_wl_factor.h"
#include "nc/lss/galaxy/nc_galaxy_shape_factor_moments_tilt.h"
#include "nc/lss/galaxy/nc_galaxy_shape_factor_spectral.h"
#include "nc_enum_types.h"
#include "nc/background/nc_hicosmo.h"
#include "nc/lss/halo/nc_halo_position.h"
#include "nc/lss/wl/nc_wl_surface_mass_density.h"
#include "nc/background/nc_distance.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_serialize.h"
#include "ncm/integration/ncm_integral1d_ptr.h"
#include "ncm/integration/ncm_integral_nd.h"
#include "ncm/integration/ncm_integrate.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#ifdef _OPENMP
#include <omp.h>
#endif /* _OPENMP */
#endif /* NUMCOSMO_GIR_SCAN */

/* The dataset that registered an #NcGalaxyWLObs as shared, stored as qdata
 * on that @obs by register_shared(). It is compared by address only and
 * never dereferenced. A clone holding the same @obs is not its owner: before
 * its first resample() writes rows it switches to a private copy of the
 * columns (see _obs_make_private()), so the realization one #NcmFitMC
 * worker draws does not reach the @obs the other workers read. */
/* *INDENT-OFF* */
G_DEFINE_QUARK (nc-data-cluster-wl-factor-shared-obs, _nc_data_cluster_wl_factor_shared_obs)
/* *INDENT-ON* */

/* Fallback -2lnP for a galaxy whose marginal is not positive. */
#define NC_GALAXY_LOW_PROB 1.0e6

typedef struct _NcDataClusterWLFactorIntArg
{
  NcGalaxyPositionFactorIntegrand *integ_pos;
  NcGalaxyRedshiftFactorIntegrand *integ_z;
  NcGalaxyShapeFactorIntegrand *integ_shape;
  NcGalaxyShapeFactorData *s_data;
} NcDataClusterWLFactorIntArg;

/* Linear-domain argument for the cubature product integrand. */
typedef struct _NcDataClusterWLFactorCubatureIntArg
{
  NcGalaxyPositionFactorIntegrand *integ_pos_lin;
  NcGalaxyRedshiftFactorIntegrand *integ_z_lin;
  NcGalaxyShapeFactorIntegrand *integ_shape_lin;
  NcGalaxyShapeFactorData *s_data;
} NcDataClusterWLFactorCubatureIntArg;

static void _nc_data_cluster_wl_factor_cubature_integrand (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);
static void _nc_data_cluster_wl_factor_cubature_dim (NcmIntegralND *intnd, guint *dim, guint *fdim);

NCM_INTEGRAL_ND_DEFINE_TYPE (NC, DATA_CLUSTER_WL_FACTOR_CUBATURE_INTEGRAND, NcDataClusterWLFactorCubatureIntegrand, nc_data_cluster_wl_factor_cubature_integrand, _nc_data_cluster_wl_factor_cubature_dim, _nc_data_cluster_wl_factor_cubature_integrand, NcDataClusterWLFactorCubatureIntArg);

/* Linear-domain product int_pos*int_z*int_shape at each of the @npoints
 * values of z in @x. */
static void
_nc_data_cluster_wl_factor_cubature_integrand (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval)
{
  NcDataClusterWLFactorCubatureIntegrand *lh_int = NC_DATA_CLUSTER_WL_FACTOR_CUBATURE_INTEGRAND (intnd);
  NcDataClusterWLFactorCubatureIntArg *arg       = &lh_int->data;
  guint i;

  for (i = 0; i < npoints; i++)
  {
    const gdouble z         = ncm_vector_fast_get (x, i);
    const gdouble int_pos   = nc_galaxy_position_factor_integrand_eval (arg->integ_pos_lin, arg->s_data->pos_data);
    const gdouble int_z     = nc_galaxy_redshift_factor_integrand_eval (arg->integ_z_lin, z, arg->s_data->z_data);
    const gdouble int_shape = nc_galaxy_shape_factor_integrand_eval (arg->integ_shape_lin, z, arg->s_data);

    ncm_vector_set (fval, i, int_pos * int_z * int_shape);
  }
}

static void
_nc_data_cluster_wl_factor_cubature_dim (NcmIntegralND *intnd, guint *dim, guint *fdim)
{
  *dim  = 1;
  *fdim = 1;
}

typedef struct _NcDataClusterWLFactorPrivate
{
  NcGalaxyWLObs *obs;

  /* keep: one gboolean per row of obs, whether the likelihood uses it, or NULL for every
   * row (NcDataClusterWLFactor:mask); mask: the rows it keeps, increasing (guint), derived
   * from keep. set_mask() replaces both arrays and never changes them. */
  GArray *keep;
  GArray *mask;
  NcGalaxyPositionFactor *position_factor;
  NcGalaxyRedshiftFactor *redshift_factor;
  NcGalaxyShapeFactor *shape_factor;

  gdouble r_min;
  gdouble r_max;
  gdouble prec;
  guint len;
  NcDataClusterWLResampleFlag resample_flag;

  NcDataClusterWLIntegMethod integ_method;

  /* CHEB_WEIGHTS: the tolerance of the nested levels, relative to the numerator of the
   * marginal of a galaxy, see _nc_data_cluster_wl_factor_cheb_levels(). */
  gdouble cheb_tol;

  /* CHEB_WEIGHTS: the level whose nodes are evaluated in the first call, 2 to 7 */
  guint cheb_level_min;

  /* Set by nc_data_cluster_wl_factor_data_prepare(): the next prepare() runs
   * the per-galaxy steps even when the staleness checks would skip them, in
   * parallel, and clears the flag. */
  gboolean force_data_prepare;

  /* TRUE while ncm_data_resample() drives prepare(), read from
   * ncm_data_is_resampling() at the top of every prepare(). That prepare()
   * only feeds gen(), which never evaluates the marginal, so the shape
   * factor's data_prepare() and the redshift grid are skipped. */
  gboolean resampling;

  /* How many NC_GALAXY_LOW_PROB substitutions happened during the LAST eval_m2lnP. Zeroed at the top of every evaluation. */
  guint low_prob_count;

  /* CHEB_WEIGHTS, for the LAST eval_m2lnP: the galaxies that reached the last level without
   * meeting cheb_tol, and the worst |B_7 - B_6| / |numerator| among them. */
  guint cheb_unconverged;
  gdouble cheb_worst_relerr;

  /* CHEB_WEIGHTS: the redshift grid of the rows of obs, the one stored with obs
   * (#NcGalaxyWLObs:cheb-grid) when it is current, otherwise one this dataset built;
   * NULL until prepare() needs one. z_nodes holds, for each galaxy, a view of its row of
   * the grid: the anchor z_lo followed by its nodes. grid_z_cl is the cluster redshift of
   * the grid. */
  NcGalaxyWLObsChebGrid *cheb_grid;
  GPtrArray *z_nodes; /* element type: NcmVector* */
  gdouble grid_z_cl;

  /* Linear-domain twins of integ_pos/integ_z/integ_shape below, used by the
   * CHEB_WEIGHTS evaluation (integ_pos_lin) and by CUBATURE's product integrand
   * (integ_z_lin/integ_shape_lin) -- all inherently linear-domain, unlike
   * integ_pos/integ_z/integ_shape (built with use_lnp=TRUE for LNINT's
   * log-domain sum). Built lazily alongside integ_pos/integ_z/integ_shape
   * and re-prepared unconditionally every cycle exactly like them. */
  NcGalaxyPositionFactorIntegrand *integ_pos_lin;
  NcGalaxyRedshiftFactorIntegrand *integ_z_lin;
  NcGalaxyShapeFactorIntegrand *integ_shape_lin;

  /* Per-galaxy data cache; element type NcGalaxyShapeFactorData (a boxed
   * type, ref/unref via nc_galaxy_shape_factor_data_{ref,unref}, NOT a
   * GObject -- hence a plain GPtrArray here rather than an NcmObjArray). */
  GPtrArray *shape_data;

  /* Rebuild trigger: set whenever obs is replaced; this is the *only*
   * structural-rebuild trigger (orthogonal to any Factor's hash -- a
   * Factor's model-derived state changing never requires reallocating
   * shape_data, only refreshing what is cached in it). */
  gboolean obs_changed;

  /*
   * The largest upper limit of the redshift integration over the galaxies of obs, set
   * when obs is loaded and required of the distance of the surface mass density at every
   * prepare (_nc_data_cluster_wl_factor_require_zf()).
   */
  gdouble z_top;

  /* Last-seen hashes, compared against each Factor's own get_hash() every
   * prepare() cycle. */
  guint64 pos_hash;
  guint64 z_hash;
  guint64 shape_radius_hash;
  guint64 shape_optzs_hash;
  guint64 shape_crit_hash;
  guint64 shape_pop_hash;

  /* Cached by prepare(), used by m2lnL_val()/resample() -- neither ever
   * peeks mset itself. */
  NcHICosmo *cosmo;
  NcHaloPosition *halo_position;
  gdouble z_cl;

  /* z-integrands, built once (lazily) and re-prepared every prepare()
   * cycle: they only capture mset-level model references, so their
   * freshness is unrelated to obs_changed or any Factor hash -- see the
   * class documentation. */
  NcGalaxyPositionFactorIntegrand *integ_pos;
  NcGalaxyRedshiftFactorIntegrand *integ_z;
  NcGalaxyShapeFactorIntegrand *integ_shape;
  NcmIntegral1dPtr *int1d;
  NcDataClusterWLFactorIntArg int_arg;

  /* CUBATURE-only state: built lazily alongside int1d, reusing
   * integ_pos_lin/integ_z_lin/integ_shape_lin (no re-preparation needed
   * beyond what those three already get every cycle). zpi/zpf/res/err are
   * length-1 vectors reused across every ncm_integral_nd_eval() call. */
  NcmIntegralND *cub_intnd;
  NcDataClusterWLFactorCubatureIntArg *cub_arg; /* owned by cub_intnd, not freed separately */
  NcmVector *cub_zpi;
  NcmVector *cub_zpf;
  NcmVector *cub_res;
  NcmVector *cub_err;

  gboolean constructed;
} NcDataClusterWLFactorPrivate;

enum
{
  PROP_0,
  PROP_OBS,
  PROP_POSITION_FACTOR,
  PROP_REDSHIFT_FACTOR,
  PROP_SHAPE_FACTOR,
  PROP_R_MIN,
  PROP_R_MAX,
  PROP_PREC,
  PROP_LEN,
  PROP_RESAMPLE_FLAG,
  PROP_INTEG_METHOD,
  PROP_CHEB_TOL,
  PROP_CHEB_LEVEL_MIN,
  PROP_MASK,
};

struct _NcDataClusterWLFactor
{
  /*< private >*/
  NcmData parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcDataClusterWLFactor, nc_data_cluster_wl_factor, NCM_TYPE_DATA)

/* Update-step dispatch: "specialize and choose once" -- the flags below are
 * decided once per prepare() cycle, then a fixed, pre-built list of steps is
 * applied uniformly per galaxy, with no per-galaxy branch on which of them
 * needed updating. Deliberately not a generic "N update levels" mechanism:
 * this is not a generic factor orchestrator, it has exactly three fixed
 * participants (position/redshift/shape), and only shape currently has more
 * than one independent update level -- see the companion Factor classes'
 * documentation for the full rationale. `mset`/`gal_i` are only used by the
 * CHEB_WEIGHTS steps below (the grid build needs mset; the grid and node steps index
 * the galaxy by gal_i) -- the five other steps ignore both, kept uniform across all
 * steps for one shared loop. */
typedef void (*NcDataClusterWLFactorStep) (NcDataClusterWLFactorPrivate *self, NcmMSet *mset, NcGalaxyShapeFactorData *s_data, guint gal_i);

static inline guint _nc_data_cluster_wl_factor_row (NcDataClusterWLFactorPrivate *self, const guint gal_i);
static NcmVector *_cheb_grid_view (NcDataClusterWLFactorPrivate *self, const guint gal_i);

static void
_step_pos (NcDataClusterWLFactorPrivate *self, NcmMSet *mset, NcGalaxyShapeFactorData *s_data, guint gal_i)
{
  nc_galaxy_position_factor_update_data (self->position_factor, s_data->pos_data);
}

static void
_step_z (NcDataClusterWLFactorPrivate *self, NcmMSet *mset, NcGalaxyShapeFactorData *s_data, guint gal_i)
{
  nc_galaxy_redshift_factor_update_data (self->redshift_factor, s_data->z_data);
}

static void
_step_shape_radius (NcDataClusterWLFactorPrivate *self, NcmMSet *mset, NcGalaxyShapeFactorData *s_data, guint gal_i)
{
  nc_galaxy_shape_factor_update_data_radius (self->shape_factor, s_data);
}

static void
_step_shape_optzs (NcDataClusterWLFactorPrivate *self, NcmMSet *mset, NcGalaxyShapeFactorData *s_data, guint gal_i)
{
  nc_galaxy_shape_factor_update_data_optzs (self->shape_factor, s_data);
}

static void
_step_shape_pop (NcDataClusterWLFactorPrivate *self, NcmMSet *mset, NcGalaxyShapeFactorData *s_data, guint gal_i)
{
  nc_galaxy_shape_factor_update_data_pop (self->shape_factor, s_data);
}

/*
 * CHEB_WEIGHTS layout. The background [a, z_hi], a = max (z_lo, z_cl), of a galaxy has the
 * nested Chebyshev-Lobatto nodes of levels 2 to 7: level L has the 2^L + 1 nodes
 * x_j = cos (pi j / 2^L), and its nodes are those of level 7 at j 2^(7 - L). The nodes are
 * stored by level, those of level 2 first and then the new ones of each level in ascending j,
 * so the 2^L + 1 nodes of level L are the first ones. _cheb_j_of_pos maps a stored position to
 * the index of the node at level 7, _cheb_pos_of_j the reverse, and _cheb_cos holds
 * cos (pi q / 128), q < 256.
 */
#define NC_DATA_CLUSTER_WL_CHEB_LEVEL_MIN (2)
#define NC_DATA_CLUSTER_WL_CHEB_LEVEL_MAX (7)
#define NC_DATA_CLUSTER_WL_CHEB_N_MAX (1 << NC_DATA_CLUSTER_WL_CHEB_LEVEL_MAX)
#define NC_DATA_CLUSTER_WL_CHEB_N_NODES (NC_DATA_CLUSTER_WL_CHEB_N_MAX + 1)

/* The weights of levels 2 to 7, 5 + 9 + 17 + 33 + 65 + 129 values, and the offset of each level */
#define NC_DATA_CLUSTER_WL_CHEB_N_WEIGHTS (258)
static const guint _cheb_w_offset[NC_DATA_CLUSTER_WL_CHEB_LEVEL_MAX + 1] = {0, 0, 0, 5, 14, 31, 64, 129};

static guint _cheb_j_of_pos[NC_DATA_CLUSTER_WL_CHEB_N_NODES];
static guint _cheb_pos_of_j[NC_DATA_CLUSTER_WL_CHEB_N_NODES];
static gdouble _cheb_cos[2 * NC_DATA_CLUSTER_WL_CHEB_N_MAX];

static void
_cheb_layout_init (void)
{
  static gsize init = 0;

  if (g_once_init_enter (&init))
  {
    guint pos = 0;
    guint L, j, q;

    for (L = NC_DATA_CLUSTER_WL_CHEB_LEVEL_MIN; L <= NC_DATA_CLUSTER_WL_CHEB_LEVEL_MAX; L++)
    {
      const guint step = NC_DATA_CLUSTER_WL_CHEB_N_MAX >> L;

      for (j = 0; j <= NC_DATA_CLUSTER_WL_CHEB_N_MAX; j += step)
      {
        /* A node of a lower level has j a multiple of 2 step */
        if ((L == NC_DATA_CLUSTER_WL_CHEB_LEVEL_MIN) || ((j / step) % 2 == 1))
        {
          _cheb_j_of_pos[pos] = j;
          _cheb_pos_of_j[j]   = pos;
          pos++;
        }
      }
    }

    g_assert_cmpuint (pos, ==, NC_DATA_CLUSTER_WL_CHEB_N_NODES);

    for (q = 0; q < 2 * NC_DATA_CLUSTER_WL_CHEB_N_MAX; q++)
      _cheb_cos[q] = cos (M_PI * q / NC_DATA_CLUSTER_WL_CHEB_N_MAX);

    g_once_init_leave (&init, 1);
  }
}

/*
 * Builds galaxy gal_i's row of the CHEB_WEIGHTS grid, at its row of obs: the anchor z_lo at
 * index 0 and the nodes of levels 2 to 7 after it (see _cheb_layout_init), the weights of
 * every level, and the exact masses, then takes the view of its nodes. With m_k = int_a^z_hi p T_k dz, k <= 128 (nc_galaxy_redshift_factor_cheb_moments()),
 * the weight of node j of level L, N = 2^L, is
 *   w_j = (2 / N) g_j sum_k c_k m_k cos (pi k j / N),  c_k, g_j = 1/2 at 0 and N, 1 otherwise,
 * so sum_j w_j F (z_j) = int_a^z_hi p F_N dz for F_N the interpolant of F on those nodes. The
 * background mass is m_0 and the foreground one the moment of order 0 on [z_lo, z_cl].
 */
static void
_step_cheb_grid (NcDataClusterWLFactorPrivate *self, NcmMSet *mset, NcGalaxyShapeFactorData *s_data, guint gal_i)
{
  GArray *w       = NULL;
  gdouble bg_norm = 0.0;
  gdouble fg_norm = 0.0;
  guint n_total   = 0;
  gdouble z_lo, z_hi, norm;
  NcmVector *z_nodes;

  _cheb_layout_init ();

  nc_galaxy_redshift_factor_get_integ_lim (self->redshift_factor, mset, s_data->z_data, &z_lo, &z_hi);
  norm = nc_galaxy_redshift_factor_norm (self->redshift_factor, mset, s_data->z_data);

  if (z_hi > self->z_cl)
  {
    const gdouble a = MAX (z_lo, self->z_cl);
    const gdouble c = 0.5 * (a + z_hi);
    const gdouble h = 0.5 * (z_hi - a);
    GArray *m_arr   = NULL;
    const gdouble *m;
    gdouble *wd;
    guint L, p;

    nc_galaxy_redshift_factor_cheb_moments (self->redshift_factor, mset, s_data->z_data, a, z_hi, NC_DATA_CLUSTER_WL_CHEB_N_MAX, &m_arr);
    m       = &g_array_index (m_arr, gdouble, 0);
    bg_norm = m[0];

    w = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), NC_DATA_CLUSTER_WL_CHEB_N_WEIGHTS);
    g_array_set_size (w, NC_DATA_CLUSTER_WL_CHEB_N_WEIGHTS);
    wd = &g_array_index (w, gdouble, 0);

    for (L = NC_DATA_CLUSTER_WL_CHEB_LEVEL_MIN; L <= NC_DATA_CLUSTER_WL_CHEB_LEVEL_MAX; L++)
    {
      const guint N    = 1u << L;
      const guint step = NC_DATA_CLUSTER_WL_CHEB_N_MAX / N;
      guint jL, k;

      for (jL = 0; jL <= N; jL++)
      {
        const gdouble g_j = ((jL == 0) || (jL == N)) ? 0.5 : 1.0;
        gdouble acc       = 0.5 * m[0];

        /* cos (pi k jL / N) = _cheb_cos[(k jL step) mod 256] */
        for (k = 1; k < N; k++)
          acc += m[k] * _cheb_cos[(k * jL * step) % (2 * NC_DATA_CLUSTER_WL_CHEB_N_MAX)];

        acc += 0.5 * m[N] * _cheb_cos[(N * jL * step) % (2 * NC_DATA_CLUSTER_WL_CHEB_N_MAX)];

        wd[_cheb_w_offset[L] + _cheb_pos_of_j[jL * step]] = (2.0 / N) * g_j * acc;
      }
    }

    n_total = NC_DATA_CLUSTER_WL_CHEB_N_NODES;
    z_nodes = ncm_vector_new (n_total + 1);

    for (p = 0; p < n_total; p++)
      ncm_vector_set (z_nodes, p + 1, c + h * _cheb_cos[_cheb_j_of_pos[p]]);

    g_array_unref (m_arr);

    if (z_lo < self->z_cl)
    {
      GArray *m_fg = NULL;

      nc_galaxy_redshift_factor_cheb_moments (self->redshift_factor, mset, s_data->z_data, z_lo, self->z_cl, 0, &m_fg);
      fg_norm = g_array_index (m_fg, gdouble, 0);
      g_array_unref (m_fg);
    }
  }
  else
  {
    z_nodes = ncm_vector_new (1);
  }

  ncm_vector_set (z_nodes, 0, z_lo);

  nc_galaxy_wl_obs_cheb_grid_set_galaxy (self->cheb_grid, _nc_data_cluster_wl_factor_row (self, gal_i), n_total, w, z_nodes, norm, bg_norm, fg_norm);
  g_ptr_array_index (self->z_nodes, gal_i) = _cheb_grid_view (self, gal_i);

  ncm_vector_free (z_nodes);

  if (w != NULL)
    g_array_unref (w);
}

static void
_step_nodes_crit (NcDataClusterWLFactorPrivate *self, NcmMSet *mset, NcGalaxyShapeFactorData *s_data, guint gal_i)
{
  nc_galaxy_shape_factor_update_data_at_nodes_crit (self->shape_factor, s_data, g_ptr_array_index (self->z_nodes, gal_i));
}

static void
_step_nodes_sigma (NcDataClusterWLFactorPrivate *self, NcmMSet *mset, NcGalaxyShapeFactorData *s_data, guint gal_i)
{
  nc_galaxy_shape_factor_update_data_at_nodes_sigma (self->shape_factor, s_data);
}

/* Calls nc_galaxy_shape_factor_data_prepare() for one galaxy with the mset
 * and the galaxy's data, so a shape factor that caches something depending
 * on both builds it here instead of lazily in the evaluation. Runs after
 * the radius/optzs/pop steps, whose results it may read. z_hi bounds the
 * galaxy's redshift support under every integration method. */
static void
_step_shape_data_prepare (NcDataClusterWLFactorPrivate *self, NcmMSet *mset, NcGalaxyShapeFactorData *s_data, guint gal_i)
{
  gdouble z_lo, z_hi;

  nc_galaxy_redshift_factor_get_integ_lim (self->redshift_factor, mset, s_data->z_data, &z_lo, &z_hi);

  nc_galaxy_shape_factor_data_prepare (self->shape_factor, mset, s_data, z_hi);
}

/* z_nodes entries are NULL until the grid step sets them */
static void
_z_nodes_free_entry (gpointer z_nodes)
{
  if (z_nodes != NULL)
    ncm_vector_free ((NcmVector *) z_nodes);
}

/* The number of galaxies of the likelihood: the rows of the mask, or of obs */
static guint
_nc_data_cluster_wl_factor_n_rows (NcDataClusterWLFactorPrivate *self)
{
  if (self->mask != NULL)
    return self->mask->len;

  return (self->obs != NULL) ? ncm_catalog_len (NCM_CATALOG (self->obs)) : 0;
}

/* The row of obs of galaxy gal_i */
static inline guint
_nc_data_cluster_wl_factor_row (NcDataClusterWLFactorPrivate *self, const guint gal_i)
{
  return (self->mask != NULL) ? g_array_index (self->mask, guint, gal_i) : gal_i;
}

/* Aborts unless the mask has one entry per row of obs */
static void
_nc_data_cluster_wl_factor_check_mask (NcDataClusterWLFactorPrivate *self)
{
  if ((self->keep != NULL) && (self->obs != NULL) && (self->keep->len != ncm_catalog_len (NCM_CATALOG (self->obs))))
    g_error ("nc_data_cluster_wl_factor: the mask has %u entries, but the catalog has %u rows.",
             self->keep->len, ncm_catalog_len (NCM_CATALOG (self->obs)));
}

/* The view of galaxy gal_i's nodes in its row of the grid: the anchor z_lo and its nodes */
static NcmVector *
_cheb_grid_view (NcDataClusterWLFactorPrivate *self, const guint gal_i)
{
  const gdouble *w, *z_nodes, *norms;
  const guint n_total = nc_galaxy_wl_obs_cheb_grid_peek_galaxy (self->cheb_grid, _nc_data_cluster_wl_factor_row (self, gal_i), &w, &z_nodes, &norms);

  return (NcmVector *) ncm_vector_const_new_data (z_nodes, n_total + 1, 1);
}

static void
_cheb_grid_clear (NcDataClusterWLFactorPrivate *self)
{
  g_ptr_array_set_size (self->z_nodes, 0);
  nc_galaxy_wl_obs_cheb_grid_clear (&self->cheb_grid);
  self->grid_z_cl = GSL_NAN;
}

/*
 * Takes the CHEB_WEIGHTS grid of the rows of obs: the one stored with obs when it is
 * current and the redshift factor does not depend on the model, with a view per galaxy;
 * otherwise a new grid for this dataset to build (_step_cheb_grid). Returns whether the
 * grid must be built. A stored grid of another redshift factor or cluster redshift aborts.
 */
static gboolean
_cheb_grid_take (NcDataClusterWLFactorPrivate *self)
{
  const guint n_gal           = self->shape_data->len;
  NcGalaxyWLObsChebGrid *grid = (nc_galaxy_redshift_factor_get_hash (self->redshift_factor) == 0) ?
                                nc_galaxy_wl_obs_peek_cheb_grid (self->obs) : NULL;
  guint i;

  g_ptr_array_set_size (self->z_nodes, n_gal);
  self->grid_z_cl = self->z_cl;

  if (grid == NULL)
  {
    gchar *data_hash = ncm_catalog_get_data_hash (NCM_CATALOG (self->obs));

    self->cheb_grid = nc_galaxy_wl_obs_cheb_grid_new (ncm_catalog_len (NCM_CATALOG (self->obs)),
                                                      NC_DATA_CLUSTER_WL_CHEB_N_WEIGHTS,
                                                      NC_DATA_CLUSTER_WL_CHEB_N_NODES + 1,
                                                      self->z_cl,
                                                      G_OBJECT_TYPE_NAME (self->redshift_factor),
                                                      data_hash);
    g_free (data_hash);

    return TRUE;
  }

  if (g_strcmp0 (nc_galaxy_wl_obs_cheb_grid_peek_z_type (grid), G_OBJECT_TYPE_NAME (self->redshift_factor)) != 0)
    g_error ("nc_data_cluster_wl_factor: the Chebyshev grid of the catalog was built by a %s, the dataset has a %s; "
             "rebuild the tables of the catalog.", nc_galaxy_wl_obs_cheb_grid_peek_z_type (grid), G_OBJECT_TYPE_NAME (self->redshift_factor));

  if (nc_galaxy_wl_obs_cheb_grid_get_z_cl (grid) != self->z_cl)
    g_error ("nc_data_cluster_wl_factor: the Chebyshev grid of the catalog was built for the cluster redshift %.15g, "
             "the dataset has %.15g; rebuild the tables of the catalog for this redshift.",
             nc_galaxy_wl_obs_cheb_grid_get_z_cl (grid), self->z_cl);

  self->cheb_grid = nc_galaxy_wl_obs_cheb_grid_ref (grid);

  for (i = 0; i < n_gal; i++)
    g_ptr_array_index (self->z_nodes, i) = _cheb_grid_view (self, i);

  return FALSE;
}

/*
 * Requires of the distance of the surface mass density a final redshift of at least
 * self->z_top, so that NcDistance evaluates the distance to every source redshift from
 * its spline. It is called at every prepare, before the prepare of the factors: a raised
 * final redshift takes effect when the shape factor prepares the surface mass density,
 * and a distance that already reaches self->z_top is not changed. Calling it at every
 * prepare covers a surface mass density replaced in the model set under the same obs.
 */
static void
_nc_data_cluster_wl_factor_require_zf (NcDataClusterWLFactorPrivate *self, NcmMSet *mset)
{
  NcWLSurfaceMassDensity *smd = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));

  g_assert_nonnull (smd);
  g_assert_nonnull (smd->dist);

  nc_distance_require_zf (smd->dist, self->z_top);
}

/*
 * Returns the largest upper limit of the redshift integration over the galaxies of obs,
 * from the redshift data that _nc_data_cluster_wl_factor_load_obs() read. The limits come
 * from the rows of obs and from the limits of the population, so prepare computes it once,
 * when obs is loaded. A redshift above it is still evaluated, by the integral of
 * NcDistance above its final redshift.
 */
static gdouble
_nc_data_cluster_wl_factor_z_top (NcDataClusterWLFactorPrivate *self, NcmMSet *mset)
{
  gdouble z_top = 0.0;
  guint i;

  for (i = 0; i < self->shape_data->len; i++)
  {
    NcGalaxyShapeFactorData *s_data = g_ptr_array_index (self->shape_data, i);
    gdouble z_lo, z_hi;

    nc_galaxy_redshift_factor_get_integ_lim (self->redshift_factor, mset, s_data->z_data, &z_lo, &z_hi);
    z_top = MAX (z_top, z_hi);
  }

  return z_top;
}

static void
nc_data_cluster_wl_factor_init (NcDataClusterWLFactor *dcwlf)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  self->obs             = NULL;
  self->keep            = NULL;
  self->mask            = NULL;
  self->position_factor = NULL;
  self->redshift_factor = NULL;
  self->shape_factor    = NULL;

  self->r_min         = 0.0;
  self->r_max         = 5.0;
  self->prec          = 1.0e-6;
  self->len           = 0;
  self->resample_flag = NC_DATA_CLUSTER_WL_RESAMPLE_FLAG_ALL;

  self->shape_data = g_ptr_array_new_with_free_func ((GDestroyNotify) nc_galaxy_shape_factor_data_unref);

  self->obs_changed       = TRUE; /* the very first prepare() always rebuilds */
  self->z_top             = 0.0;
  self->pos_hash          = 0;
  self->z_hash            = 0;
  self->shape_radius_hash = 0;
  self->shape_optzs_hash  = 0;
  self->shape_crit_hash   = 0;
  self->shape_pop_hash    = 0;

  self->cosmo         = NULL;
  self->halo_position = NULL;
  self->z_cl          = 0.0;

  self->integ_pos      = NULL;
  self->integ_z        = NULL;
  self->integ_shape    = NULL;
  self->int1d          = NULL;
  self->int_arg.s_data = NULL;

  self->integ_method   = NC_DATA_CLUSTER_WL_INTEG_METHOD_LNINT;
  self->low_prob_count = 0;
  self->cheb_tol       = 0.0;
  self->cheb_level_min = 0;

  self->cheb_unconverged  = 0;
  self->cheb_worst_relerr = 0.0;

  self->force_data_prepare = FALSE;

  self->cheb_grid = NULL;
  self->z_nodes   = g_ptr_array_new_with_free_func (&_z_nodes_free_entry);
  self->grid_z_cl = GSL_NAN;

  self->integ_pos_lin   = NULL;
  self->integ_z_lin     = NULL;
  self->integ_shape_lin = NULL;

  self->cub_intnd = NULL;
  self->cub_arg   = NULL;
  self->cub_zpi   = NULL;
  self->cub_zpf   = NULL;
  self->cub_res   = NULL;
  self->cub_err   = NULL;

  self->constructed = FALSE;
}

static void
_nc_data_cluster_wl_factor_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataClusterWLFactor *dcwlf              = NC_DATA_CLUSTER_WL_FACTOR (object);
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  g_return_if_fail (NC_IS_DATA_CLUSTER_WL_FACTOR (object));

  switch (prop_id)
  {
    case PROP_OBS:
      nc_data_cluster_wl_factor_set_obs (dcwlf, g_value_get_object (value));
      break;
    case PROP_POSITION_FACTOR:
      g_assert_null (self->position_factor); /* construct-only */
      self->position_factor = g_value_dup_object (value);
      break;
    case PROP_REDSHIFT_FACTOR:
      g_assert_null (self->redshift_factor); /* construct-only */
      self->redshift_factor = g_value_dup_object (value);
      break;
    case PROP_SHAPE_FACTOR:
      g_assert_null (self->shape_factor); /* construct-only */
      self->shape_factor = g_value_dup_object (value);
      break;
    case PROP_R_MIN:
      nc_data_cluster_wl_factor_set_cut (dcwlf, g_value_get_double (value), self->r_max);
      break;
    case PROP_R_MAX:
      nc_data_cluster_wl_factor_set_cut (dcwlf, self->r_min, g_value_get_double (value));
      break;
    case PROP_PREC:
      nc_data_cluster_wl_factor_set_prec (dcwlf, g_value_get_double (value));
      break;
    case PROP_RESAMPLE_FLAG:
      nc_data_cluster_wl_factor_set_resample_flag (dcwlf, g_value_get_flags (value));
      break;
    case PROP_INTEG_METHOD:
      nc_data_cluster_wl_factor_set_integ_method (dcwlf, g_value_get_enum (value));
      break;
    case PROP_CHEB_TOL:
      nc_data_cluster_wl_factor_set_cheb_tol (dcwlf, g_value_get_double (value));
      break;
    case PROP_CHEB_LEVEL_MIN:
      nc_data_cluster_wl_factor_set_cheb_level_min (dcwlf, g_value_get_uint (value));
      break;
    case PROP_MASK:
    {
      GVariant *var = g_value_get_variant (value);

      if (var == NULL)
      {
        nc_data_cluster_wl_factor_set_mask (dcwlf, NULL);
      }
      else
      {
        /* The booleans of a GVariant array are one byte each */
        gsize n_rows        = 0;
        const guchar *bytes = g_variant_get_fixed_array (var, &n_rows, sizeof (guchar));
        GArray *keep        = g_array_sized_new (FALSE, FALSE, sizeof (gboolean), n_rows);
        gsize i;

        for (i = 0; i < n_rows; i++)
        {
          const gboolean keep_i = (bytes[i] != 0);

          g_array_append_val (keep, keep_i);
        }

        nc_data_cluster_wl_factor_set_mask (dcwlf, keep);
        g_array_unref (keep);
      }

      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_data_cluster_wl_factor_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataClusterWLFactor *dcwlf              = NC_DATA_CLUSTER_WL_FACTOR (object);
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  g_return_if_fail (NC_IS_DATA_CLUSTER_WL_FACTOR (object));

  switch (prop_id)
  {
    case PROP_OBS:
      g_value_set_object (value, nc_data_cluster_wl_factor_peek_obs (dcwlf));
      break;
    case PROP_POSITION_FACTOR:
      g_value_set_object (value, self->position_factor);
      break;
    case PROP_REDSHIFT_FACTOR:
      g_value_set_object (value, self->redshift_factor);
      break;
    case PROP_SHAPE_FACTOR:
      g_value_set_object (value, self->shape_factor);
      break;
    case PROP_R_MIN:
      g_value_set_double (value, self->r_min);
      break;
    case PROP_R_MAX:
      g_value_set_double (value, self->r_max);
      break;
    case PROP_PREC:
      g_value_set_double (value, self->prec);
      break;
    case PROP_LEN:
      g_value_set_uint (value, self->len);
      break;
    case PROP_RESAMPLE_FLAG:
      g_value_set_flags (value, nc_data_cluster_wl_factor_get_resample_flag (dcwlf));
      break;
    case PROP_INTEG_METHOD:
      g_value_set_enum (value, self->integ_method);
      break;
    case PROP_CHEB_TOL:
      g_value_set_double (value, self->cheb_tol);
      break;
    case PROP_CHEB_LEVEL_MIN:
      g_value_set_uint (value, self->cheb_level_min);
      break;
    case PROP_MASK:

      if (self->keep == NULL)
      {
        g_value_set_variant (value, NULL);
      }
      else
      {
        guchar *bytes = g_new (guchar, MAX (self->keep->len, 1));
        guint i;

        for (i = 0; i < self->keep->len; i++)
          bytes[i] = g_array_index (self->keep, gboolean, i) ? 1 : 0;

        g_value_set_variant (value, g_variant_new_fixed_array (G_VARIANT_TYPE_BOOLEAN, bytes, self->keep->len, sizeof (guchar)));
        g_free (bytes);
      }

      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

/* Removes the owner mark register_shared() set on @obs, if this dataset set
 * it, so that a dataset later holding @obs alone resamples it in place. */
static void
_shared_obs_owner_clear (NcDataClusterWLFactor *dcwlf, NcGalaxyWLObs *obs)
{
  if ((obs != NULL) && (g_object_get_qdata (G_OBJECT (obs), _nc_data_cluster_wl_factor_shared_obs_quark ()) == dcwlf))
    g_object_set_qdata (G_OBJECT (obs), _nc_data_cluster_wl_factor_shared_obs_quark (), NULL);
}

/*
 * Replaces @obs by a copy of its columns that references its p(z) objects.
 * The columns, of the same names and types, are the rows resample() writes
 * (about 100 bytes per galaxy); the p(z) objects are most of the catalog and
 * are only read, so they are added with nc_galaxy_wl_obs_set_pz(), which
 * references them. The NcmCatalog calls get a NULL error: the columns come
 * from @obs itself, so a failure aborts.
 */
static void
_obs_make_private (NcDataClusterWLFactorPrivate *self)
{
  NcGalaxyWLObs *shared  = self->obs;
  NcmCatalog *shared_cat = NCM_CATALOG (shared);
  const guint len        = ncm_catalog_len (shared_cat);
  GStrv columns          = ncm_catalog_peek_columns (shared_cat);
  const guint ncols      = g_strv_length (columns);
  GArray *types          = g_array_sized_new (FALSE, FALSE, sizeof (NcmCatalogColType), ncols);
  NcGalaxyWLObs *copy;
  NcmCatalog *copy_cat;
  guint i;

  for (i = 0; i < ncols; i++)
  {
    const NcmCatalogColType t = ncm_catalog_get_col_type (shared_cat, columns[i], NULL);

    g_array_append_val (types, t);
  }

  copy = nc_galaxy_wl_obs_new (nc_galaxy_wl_obs_get_ellip_conv (shared),
                               nc_galaxy_wl_obs_get_coord (shared),
                               len, columns,
                               (const NcmCatalogColType *) types->data, ncols);
  copy_cat = NCM_CATALOG (copy);

  ncm_matrix_memcpy (ncm_catalog_peek_data (copy_cat), ncm_catalog_peek_data (shared_cat));

  for (i = 0; i < ncols; i++)
  {
    GArray *vals = NULL;

    switch (g_array_index (types, NcmCatalogColType, i))
    {
      case NCM_CATALOG_COL_TYPE_INT:
        ncm_catalog_get_int_column (shared_cat, columns[i], &vals, NULL);
        ncm_catalog_set_int_column (copy_cat, columns[i], vals, NULL);
        break;
      case NCM_CATALOG_COL_TYPE_BOOL:
        ncm_catalog_get_bool_column (shared_cat, columns[i], &vals, NULL);
        ncm_catalog_set_bool_column (copy_cat, columns[i], vals, NULL);
        break;
      default:
        break;
    }

    g_clear_pointer (&vals, g_array_unref);
  }

  ncm_catalog_mark_outdated (copy_cat);
  g_array_unref (types);

  for (i = 0; i < len; i++)
  {
    NcGalaxyPz *pz = nc_galaxy_wl_obs_peek_pz (shared, i);

    if (pz != NULL)
      nc_galaxy_wl_obs_set_pz (copy, i, pz);
  }

  nc_galaxy_wl_obs_free (shared);
  self->obs = copy;
}

static void
_nc_data_cluster_wl_factor_dispose (GObject *object)
{
  NcDataClusterWLFactor *dcwlf              = NC_DATA_CLUSTER_WL_FACTOR (object);
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  _shared_obs_owner_clear (dcwlf, self->obs);
  nc_galaxy_wl_obs_clear (&self->obs);
  nc_galaxy_position_factor_clear (&self->position_factor);
  nc_galaxy_redshift_factor_clear (&self->redshift_factor);
  nc_galaxy_shape_factor_clear (&self->shape_factor);

  g_clear_pointer (&self->shape_data, g_ptr_array_unref);

  nc_hicosmo_clear (&self->cosmo);
  nc_halo_position_clear (&self->halo_position);

  g_clear_pointer (&self->integ_pos, nc_galaxy_position_factor_integrand_free);
  g_clear_pointer (&self->integ_z, nc_galaxy_redshift_factor_integrand_free);
  g_clear_pointer (&self->integ_shape, nc_galaxy_shape_factor_integrand_free);
  g_clear_pointer (&self->int1d, ncm_integral1d_ptr_free);
  g_clear_pointer (&self->integ_pos_lin, nc_galaxy_position_factor_integrand_free);
  g_clear_pointer (&self->integ_z_lin, nc_galaxy_redshift_factor_integrand_free);
  g_clear_pointer (&self->integ_shape_lin, nc_galaxy_shape_factor_integrand_free);

  /* cub_arg is owned by cub_intnd (embedded, not a separate allocation);
   * freeing cub_intnd is enough. */
  self->cub_arg = NULL;
  g_clear_pointer (&self->cub_intnd, ncm_integral_nd_free);
  g_clear_pointer (&self->cub_zpi, ncm_vector_free);
  g_clear_pointer (&self->cub_zpf, ncm_vector_free);
  g_clear_pointer (&self->cub_res, ncm_vector_free);
  g_clear_pointer (&self->cub_err, ncm_vector_free);

  g_clear_pointer (&self->z_nodes, g_ptr_array_unref);
  nc_galaxy_wl_obs_cheb_grid_clear (&self->cheb_grid);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_data_cluster_wl_factor_parent_class)->dispose (object);
}

static void
_nc_data_cluster_wl_factor_finalize (GObject *object)
{
  NcDataClusterWLFactor *dcwlf              = NC_DATA_CLUSTER_WL_FACTOR (object);
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  g_clear_pointer (&self->keep, g_array_unref);
  g_clear_pointer (&self->mask, g_array_unref);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_data_cluster_wl_factor_parent_class)->finalize (object);
}

static void
_nc_data_cluster_wl_factor_constructed (GObject *object)
{
  /* Chain up: start */
  G_OBJECT_CLASS (nc_data_cluster_wl_factor_parent_class)->constructed (object);
  {
    NcDataClusterWLFactor *dcwlf              = NC_DATA_CLUSTER_WL_FACTOR (object);
    NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

    g_assert_nonnull (self->position_factor);
    g_assert_nonnull (self->redshift_factor);
    g_assert_nonnull (self->shape_factor);
    g_assert_cmpfloat (self->r_min, <, self->r_max);

    self->constructed = TRUE;
  }
}

/* Rebuild cascade, per galaxy (near-direct translation of
 * _nc_data_cluster_wl_load_obs): builds the linked pos_data/z_data/s_data
 * triple and reads the catalog row through the shape level's own cascading
 * read_row (which already reaches down through pos_data/z_data/pop_data/
 * ldata in one call). */
static void
_nc_data_cluster_wl_factor_load_obs (NcDataClusterWLFactorPrivate * const self, NcmMSet *mset, NcGalaxyWLObs *obs)
{
  guint i;

  g_ptr_array_set_size (self->shape_data, 0);
  _nc_data_cluster_wl_factor_check_mask (self);

  for (i = 0; i < _nc_data_cluster_wl_factor_n_rows (self); i++)
  {
    NcGalaxyPositionFactorData *pos_data = nc_galaxy_position_factor_data_new (self->position_factor, mset);
    NcGalaxyRedshiftFactorData *z_data   = nc_galaxy_redshift_factor_data_new (self->redshift_factor, mset);
    NcGalaxyShapeFactorData *s_data      = nc_galaxy_shape_factor_data_new (self->shape_factor, mset, pos_data, z_data);

    nc_galaxy_shape_factor_data_read_row (s_data, obs, _nc_data_cluster_wl_factor_row (self, i));
    g_ptr_array_add (self->shape_data, s_data);

    nc_galaxy_position_factor_data_unref (pos_data);
    nc_galaxy_redshift_factor_data_unref (z_data);
  }
}

static gdouble
_nc_data_cluster_wl_factor_lnint_integrand (gpointer user_data, const gdouble z, const gdouble w)
{
  NcDataClusterWLFactorIntArg *arg = (NcDataClusterWLFactorIntArg *) user_data;
  const gdouble int_pos            = nc_galaxy_position_factor_integrand_eval (arg->integ_pos, arg->s_data->pos_data);
  const gdouble int_z              = nc_galaxy_redshift_factor_integrand_eval (arg->integ_z, z, arg->s_data->z_data);
  const gdouble int_shape          = nc_galaxy_shape_factor_integrand_eval (arg->integ_shape, z, arg->s_data);

  return int_pos + int_z + int_shape;
}

static void
_nc_data_cluster_wl_factor_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataClusterWLFactor *dcwlf              = NC_DATA_CLUSTER_WL_FACTOR (data);
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);
  GError *error                             = NULL;
  gboolean pos_changed, z_changed, radius_changed, optzs_changed, crit_changed, pop_changed, grid_changed;
  NcDataClusterWLFactorStep steps[9];
  guint n_steps       = 0;
  gboolean cheb_build = FALSE;

  self->resampling = ncm_data_is_resampling (data);

  /* Orchestrator's own direct needs only -- z_cl for the split-integration
   * logic below, halo_position+cosmo for resample()'s radius rejection. */
  nc_hicosmo_clear (&self->cosmo);
  nc_halo_position_clear (&self->halo_position);
  self->cosmo         = nc_hicosmo_ref (NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ())));
  self->halo_position = nc_halo_position_ref (NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ())));
  g_assert_nonnull (self->cosmo);
  g_assert_nonnull (self->halo_position);

  nc_halo_position_prepare_if_needed (self->halo_position, self->cosmo);
  self->z_cl = nc_halo_position_get_redshift (self->halo_position);

  if (!nc_galaxy_shape_factor_check_obs (self->shape_factor, self->obs, &error))
    g_error ("nc_data_cluster_wl_factor_prepare: %s", error->message);

  /* Captured before self->obs_changed is reset to FALSE just below: the change flags
   * computed further down all need to know whether obs was replaced *this* cycle, and
   * must read that fact through this local -- reading self->obs_changed itself at that
   * point would always see FALSE. A pure obs swap changes no NcmModel, so obs_was_changed
   * is the only thing that signals that the grid must be taken again. */
  const gboolean obs_was_changed = self->obs_changed;

  if (obs_was_changed)
  {
    _nc_data_cluster_wl_factor_load_obs (self, mset, self->obs);
    self->obs_changed = FALSE;
    self->z_top       = _nc_data_cluster_wl_factor_z_top (self, mset);
  }

  _nc_data_cluster_wl_factor_require_zf (self, mset);

  nc_galaxy_position_factor_prepare (self->position_factor, mset);
  nc_galaxy_redshift_factor_prepare (self->redshift_factor, mset);
  nc_galaxy_shape_factor_prepare (self->shape_factor, mset);

  /* z-integrands: built once, re-prepared unconditionally every cycle (see
   * the class documentation for why this is decoupled from obs_changed and
   * from every Factor's hash). */
  if (self->integ_pos == NULL)
  {
    self->integ_pos   = nc_galaxy_position_factor_integ (self->position_factor, mset, TRUE);
    self->integ_z     = nc_galaxy_redshift_factor_integ (self->redshift_factor, mset, TRUE);
    self->integ_shape = nc_galaxy_shape_factor_integ (self->shape_factor, mset, TRUE);

    /* Linear-domain twins of integ_pos/integ_z/integ_shape, needed by the
     * CHEB_WEIGHTS evaluation (integ_pos_lin) and CUBATURE's product integrand
     * (integ_z_lin/integ_shape_lin) -- all inherently linear-domain, see
     * the class documentation. */
    self->integ_pos_lin   = nc_galaxy_position_factor_integ (self->position_factor, mset, FALSE);
    self->integ_z_lin     = nc_galaxy_redshift_factor_integ (self->redshift_factor, mset, FALSE);
    self->integ_shape_lin = nc_galaxy_shape_factor_integ (self->shape_factor, mset, FALSE);

    self->int_arg.integ_pos   = self->integ_pos;
    self->int_arg.integ_z     = self->integ_z;
    self->int_arg.integ_shape = self->integ_shape;

    self->int1d = ncm_integral1d_ptr_new (&_nc_data_cluster_wl_factor_lnint_integrand, NULL);
    ncm_integral1d_ptr_set_userdata (self->int1d, &self->int_arg);

    /* CUBATURE-only state: an NcmIntegralND wrapping the linear-domain
     * product int_pos_lin*int_z_lin*int_shape_lin, built once here alongside
     * int1d -- no re-preparation needed beyond what the three _lin twins
     * already get every cycle, since the ND integrator object itself is
     * stateless across calls. */
    self->cub_intnd = g_object_new (nc_data_cluster_wl_factor_cubature_integrand_get_type (), NULL);
    self->cub_arg   = &NC_DATA_CLUSTER_WL_FACTOR_CUBATURE_INTEGRAND (self->cub_intnd)->data;

    self->cub_arg->integ_pos_lin   = self->integ_pos_lin;
    self->cub_arg->integ_z_lin     = self->integ_z_lin;
    self->cub_arg->integ_shape_lin = self->integ_shape_lin;

    ncm_integral_nd_set_reltol (self->cub_intnd, self->prec);
    ncm_integral_nd_set_abstol (self->cub_intnd, 0.0);
    ncm_integral_nd_set_method (self->cub_intnd, NCM_INTEGRAL_ND_METHOD_CUBATURE_H_V);

    self->cub_zpi = ncm_vector_new (1);
    self->cub_zpf = ncm_vector_new (1);
    self->cub_res = ncm_vector_new (1);
    self->cub_err = ncm_vector_new (1);
  }
  else
  {
    nc_galaxy_position_factor_integrand_prepare (self->integ_pos, mset);
    nc_galaxy_redshift_factor_integrand_prepare (self->integ_z, mset);
    nc_galaxy_shape_factor_integrand_prepare (self->integ_shape, mset);
    nc_galaxy_position_factor_integrand_prepare (self->integ_pos_lin, mset);
    nc_galaxy_redshift_factor_integrand_prepare (self->integ_z_lin, mset);
    nc_galaxy_shape_factor_integrand_prepare (self->integ_shape_lin, mset);
  }

  ncm_integral1d_set_reltol (NCM_INTEGRAL1D (self->int1d), self->prec);
  ncm_integral1d_set_abstol (NCM_INTEGRAL1D (self->int1d), 0.0);
  ncm_integral_nd_set_reltol (self->cub_intnd, self->prec);

  /* radius_changed: galaxy geometry changed (halo_position/cosmo).
   * optzs_changed:  lens model changed (cosmo/density_profile/mass_summary/
   *                 surface_mass_density/halo_position).
   * grid_changed:   CHEB_WEIGHTS needs its redshift grid taken again. The grid depends on
   *                 the rows of obs, the redshift factor and z_cl only, so it is dropped
   *                 when one of those changes; radius_hash would also bump on RA/Dec of a
   *                 miscentering model, which leave the grid unchanged. The resampling
   *                 prepare() feeds gen() only and takes no grid. */
  pos_changed    = obs_was_changed || (nc_galaxy_position_factor_get_hash (self->position_factor) != self->pos_hash);
  z_changed      = obs_was_changed || (nc_galaxy_redshift_factor_get_hash (self->redshift_factor) != self->z_hash);
  radius_changed = obs_was_changed || (nc_galaxy_shape_factor_get_radius_hash (self->shape_factor) != self->shape_radius_hash);
  optzs_changed  = obs_was_changed || (nc_galaxy_shape_factor_get_optzs_hash (self->shape_factor) != self->shape_optzs_hash);
  crit_changed   = obs_was_changed || (nc_galaxy_shape_factor_get_crit_hash (self->shape_factor) != self->shape_crit_hash);
  pop_changed    = obs_was_changed || (nc_galaxy_shape_factor_get_pop_hash (self->shape_factor) != self->shape_pop_hash);

  if (z_changed || (self->z_cl != self->grid_z_cl))
    _cheb_grid_clear (self);

  grid_changed = (self->integ_method == NC_DATA_CLUSTER_WL_INTEG_METHOD_CHEB_WEIGHTS) &&
                 (self->cheb_grid == NULL) && !self->resampling;

  self->pos_hash          = nc_galaxy_position_factor_get_hash (self->position_factor);
  self->z_hash            = nc_galaxy_redshift_factor_get_hash (self->redshift_factor);
  self->shape_radius_hash = nc_galaxy_shape_factor_get_radius_hash (self->shape_factor);
  self->shape_optzs_hash  = nc_galaxy_shape_factor_get_optzs_hash (self->shape_factor);
  self->shape_crit_hash   = nc_galaxy_shape_factor_get_crit_hash (self->shape_factor);
  self->shape_pop_hash    = nc_galaxy_shape_factor_get_pop_hash (self->shape_factor);

  if (!pos_changed && !z_changed && !radius_changed && !optzs_changed && !crit_changed && !pop_changed && !grid_changed &&
      !self->force_data_prepare)
    return;  /* nothing to do this cycle */

  if (grid_changed)
    cheb_build = _cheb_grid_take (self);

  /* Specialize and choose once: no per-galaxy branch on any of these. Ordering
   * matters: _step_shape_radius before _step_nodes_sigma (sigma reads the projected
   * radius) and _step_cheb_grid before _step_nodes_crit (crit reads the nodes). */
  if (pos_changed)
    steps[n_steps++] = &_step_pos;

  if (z_changed)
    steps[n_steps++] = &_step_z;

  if (radius_changed)
    steps[n_steps++] = &_step_shape_radius;

  if (optzs_changed)
    steps[n_steps++] = &_step_shape_optzs;

  if (pop_changed)
    steps[n_steps++] = &_step_shape_pop;

  /* _step_shape_data_prepare runs after radius/optzs/pop, whose results it reads. It is
   * not gated on radius_changed or optzs_changed: what a shape factor precomputes here
   * holds over the whole parameter box, so running it on every halo move would repeat
   * per-galaxy work for a result that does not change. A shape factor whose cache
   * depends on the current point must key it by value and refresh lazily, as it must for
   * resample() anyway.
   *
   * Skipped while resampling: gen() never evaluates the marginal, and resample() sets
   * obs_changed, so the next prepare() runs this step on the new data. */
  if ((z_changed || pop_changed || grid_changed || self->force_data_prepare) && !self->resampling)
    steps[n_steps++] = &_step_shape_data_prepare;

  if (cheb_build)
    steps[n_steps++] = &_step_cheb_grid;

  if (self->cheb_grid != NULL)
  {
    if (grid_changed || crit_changed)
      steps[n_steps++] = &_step_nodes_crit;

    if (grid_changed || radius_changed || optzs_changed)
      steps[n_steps++] = &_step_nodes_sigma;
  }

  {
    const guint n_gal       = self->shape_data->len;
    const gboolean parallel = self->force_data_prepare || cheb_build;
    guint gal_i;

    /* Serial loop: runs every step for one galaxy at a time until a galaxy has run every
     * step (in a grid build, a galaxy with background support, which reaches the
     * redshift factor's moments), or for every galaxy when not parallel. Several models
     * the factors reach update shared state lazily on the first call after their
     * parameters change: the halo position's rotation, the redshift population's derived
     * constants, the splines of the numerically integrated density profiles and Cuba's
     * settings. Those first calls must not come from several threads at once: the writes
     * would race, and in the profile's case corrupt a shared spline. The serial galaxy
     * performs all of them. */
    for (gal_i = 0; gal_i < n_gal; gal_i++)
    {
      NcGalaxyShapeFactorData *s_data = g_ptr_array_index (self->shape_data, gal_i);
      guint j;

      for (j = 0; j < n_steps; j++)
        steps[j](self, mset, s_data, gal_i);

      if (parallel && (!cheb_build || (ncm_vector_len (g_ptr_array_index (self->z_nodes, gal_i)) > 1)))
      {
        gal_i++;
        break;
      }
    }

    /* Parallel loop over the remaining galaxies, from
     * nc_data_cluster_wl_factor_data_prepare() or a CHEB_WEIGHTS grid build. Every
     * position, redshift and shape factor supports this: their per-galaxy methods write
     * only the galaxy's own data (see the thread safety note of #NcGalaxyShapeFactor,
     * #NcGalaxyPositionFactor and #NcGalaxyRedshiftFactor), and each galaxy writes only
     * its own row of the grid and its own view. Any other per-cycle prepare() of a fit is
     * serial: it runs inside the sampler's own threads, where a nested region adds only
     * overhead. */
    if (parallel && (gal_i < n_gal))
    {
      const gint first = (gint) gal_i;
      const gint last  = (gint) n_gal;
      gint k;

#ifdef _OPENMP
      #pragma omp parallel for schedule(dynamic, 1) if (!omp_in_parallel ())
#endif /* _OPENMP */

      for (k = first; k < last; k++)
      {
        NcGalaxyShapeFactorData *s_data = g_ptr_array_index (self->shape_data, k);
        guint j;

        for (j = 0; j < n_steps; j++)
          steps[j](self, mset, s_data, (guint) k);
      }
    }

    self->force_data_prepare = FALSE;
  }
}

/* Combine two adaptive sub-integrals split at z_cl. Each panel returns
 * -2 ln I_panel; recombine in the linear domain I = I_lo + I_hi. */
static gdouble
_nc_data_cluster_wl_factor_integ_combine (gdouble m2_lo, gdouble m2_hi)
{
  if (!gsl_finite (m2_lo))
    return m2_hi;

  if (!gsl_finite (m2_hi))
    return m2_lo;

  {
    const gdouble a = -0.5 * m2_lo;
    const gdouble b = -0.5 * m2_hi;
    const gdouble m = MAX (a, b);

    return -2.0 * (m + log (exp (a - m) + exp (b - m)));
  }
}

static gdouble
_nc_data_cluster_wl_factor_integrate (NcDataClusterWLFactorPrivate * const self, gdouble zmin, gdouble zmax)
{
  gdouble err;

  return -2.0 * ncm_integral1d_eval_lnint (NCM_INTEGRAL1D (self->int1d), zmin, zmax, &err);
}

static gdouble
_nc_data_cluster_wl_factor_eval_m2lnP_lnint (NcDataClusterWLFactor *dcwlf, NcmMSet *mset, NcmVector *m2lnP_gal)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);
  NcmData *data                             = NCM_DATA (dcwlf);
  const gboolean use_bootstrap              = ncm_data_bootstrap_enabled (data);
  NcmBootstrap *bstrap                      = use_bootstrap ? ncm_data_peek_bootstrap (data) : NULL;
  const guint n_iter                        = use_bootstrap ? ncm_bootstrap_get_bsize (bstrap) : self->shape_data->len;
  gdouble result                            = 0.0;
  guint i;

  /* LNINT never substitutes NC_GALAXY_LOW_PROB -- it skips a non-finite galaxy and
   * stores NaN -- so its count is zero. Reset it anyway, or get_low_prob_count() would
   * report a stale count from an earlier CHEB_WEIGHTS/CUBATURE evaluation on the same
   * object. */
  self->low_prob_count = 0;

  if (m2lnP_gal != NULL)
    g_assert_cmpuint (ncm_vector_len (m2lnP_gal), ==, self->len);

  for (i = 0; i < n_iter; i++)
  {
    const guint gal_i               = use_bootstrap ? ncm_bootstrap_get (bstrap, i) : i;
    NcGalaxyShapeFactorData *s_data = g_ptr_array_index (self->shape_data, gal_i);
    gdouble zpi, zpf, m2lnP_gal_i;

    self->int_arg.s_data = s_data;

    nc_galaxy_redshift_factor_get_integ_lim (self->redshift_factor, mset, s_data->z_data, &zpi, &zpf);

    /* Split at z_cl when it falls strictly inside the support, since the
     * reduced shear has a kink there. */
    if ((self->z_cl > zpi) && (self->z_cl < zpf))
    {
      const gdouble m2_lo = _nc_data_cluster_wl_factor_integrate (self, zpi, self->z_cl);
      const gdouble m2_hi = _nc_data_cluster_wl_factor_integrate (self, self->z_cl, zpf);

      m2lnP_gal_i = _nc_data_cluster_wl_factor_integ_combine (m2_lo, m2_hi);
    }
    else
    {
      m2lnP_gal_i = _nc_data_cluster_wl_factor_integrate (self, zpi, zpf);
    }

    if (!gsl_finite (m2lnP_gal_i))
    {
      g_warning ("_nc_data_cluster_wl_factor_eval_m2lnP_lnint: galaxy %d has undefined likelihood [%g]. Skipping it.", gal_i, m2lnP_gal_i);

      if (m2lnP_gal != NULL)
        ncm_vector_set (m2lnP_gal, gal_i, GSL_NAN);

      continue;
    }

    if (m2lnP_gal != NULL)
      ncm_vector_set (m2lnP_gal, gal_i, m2lnP_gal_i);

    result += m2lnP_gal_i;
  }

  return result;
}

/*
 * The numerator F_0 m_fg + B of a galaxy's marginal under CHEB_WEIGHTS, where B is the
 * background integral sum_j w_j F (z_j) of a level with the weights @w, F_0 the value at the
 * anchor and m_fg = @fg the foreground mass. The
 * first call evaluates the nodes of level L_0 = cheb_level_min, and the levels from 3 to L_0
 * are summed on those values; the result is B_{L_0} when |B_{L_0} - B_{L_0 - 1}| <=
 * cheb_tol |F_0 m_fg + B_{L_0}|, and otherwise the levels above L_0 are evaluated one call
 * each, through the range evaluation of the shape factor, until a level meets the same test
 * or level 7 is reached. With @want_log the shape factor gives ln F, and the values are taken
 * as exp (ln F - s), s the largest ln F of the first call, given in @s; the quadrature is
 * linear in them. Writes in @converged whether a level met the tolerance and in @relerr the
 * last relative difference.
 */
static gdouble
_nc_data_cluster_wl_factor_cheb_levels (NcDataClusterWLFactorPrivate * const self, NcmMSet *mset, NcGalaxyShapeFactorData *s_data, const gdouble *w,
                                        const gdouble fg, NcmVector *z_nodes, NcmVector *vals, const gboolean want_log, gdouble *s, gboolean *converged,
                                        gdouble *relerr)
{
  const guint L_0 = self->cheb_level_min;
  gdouble *f      = ncm_vector_data (vals);
  guint n_prev    = (1u << L_0) + 1;
  gdouble B_prev  = 0.0;
  gdouble num     = 0.0;
  guint L, k;

  if (want_log)
    nc_galaxy_shape_factor_eval_ln_at_nodes_range (self->shape_factor, mset, s_data, z_nodes, 0, 1 + n_prev, vals);
  else
    nc_galaxy_shape_factor_eval_at_nodes_range (self->shape_factor, mset, s_data, z_nodes, 0, 1 + n_prev, vals);

  *s         = 0.0;
  *converged = FALSE;
  *relerr    = 0.0;

  if (want_log)
  {
    *s = GSL_NEGINF;

    for (k = 0; k < 1 + n_prev; k++)
      *s = MAX (*s, f[k]);

    for (k = 0; k < 1 + n_prev; k++)
      f[k] = exp (f[k] - *s);
  }

  for (L = NC_DATA_CLUSTER_WL_CHEB_LEVEL_MIN; L <= NC_DATA_CLUSTER_WL_CHEB_LEVEL_MAX; L++)
  {
    const guint n = (1u << L) + 1;
    gdouble B     = 0.0;

    /* Above L_0 the level's new nodes are evaluated now */
    if (L > L_0)
    {
      if (want_log)
      {
        nc_galaxy_shape_factor_eval_ln_at_nodes_range (self->shape_factor, mset, s_data, z_nodes, 1 + n_prev, 1 + n, vals);

        for (k = 1 + n_prev; k < 1 + n; k++)
          f[k] = exp (f[k] - *s);
      }
      else
      {
        nc_galaxy_shape_factor_eval_at_nodes_range (self->shape_factor, mset, s_data, z_nodes, 1 + n_prev, 1 + n, vals);
      }

      n_prev = n;
    }

    for (k = 0; k < n; k++)
      B += w[_cheb_w_offset[L] + k] * f[1 + k];

    num = f[0] * fg + B;

    /* Levels below L_0 only lead to the comparison at L_0 */
    if ((L >= L_0) && (L > NC_DATA_CLUSTER_WL_CHEB_LEVEL_MIN))
    {
      const gdouble diff = fabs (B - B_prev);

      *relerr = (num != 0.0) ? diff / fabs (num) : ((diff == 0.0) ? 0.0 : GSL_POSINF);

      if (diff <= self->cheb_tol * fabs (num))
      {
        *converged = TRUE;

        return num;
      }
    }

    B_prev = B;
  }

  return num;
}

static gdouble
_nc_data_cluster_wl_factor_eval_m2lnP_cheb (NcDataClusterWLFactor *dcwlf, NcmMSet *mset, NcmVector *m2lnP_gal)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);
  NcmData *data                             = NCM_DATA (dcwlf);
  const gboolean use_bootstrap              = ncm_data_bootstrap_enabled (data);
  NcmBootstrap *bstrap                      = use_bootstrap ? ncm_data_peek_bootstrap (data) : NULL;
  const guint n_iter                        = use_bootstrap ? ncm_bootstrap_get_bsize (bstrap) : self->shape_data->len;
  NcmVector *vals                           = ncm_vector_new (NC_DATA_CLUSTER_WL_CHEB_N_NODES + 1);
  gdouble result                            = 0.0;
  guint i;

  self->low_prob_count    = 0;
  self->cheb_unconverged  = 0;
  self->cheb_worst_relerr = 0.0;

  if (m2lnP_gal != NULL)
    g_assert_cmpuint (ncm_vector_len (m2lnP_gal), ==, self->len);

  for (i = 0; i < n_iter; i++)
  {
    const guint gal_i               = use_bootstrap ? ncm_bootstrap_get (bstrap, i) : i;
    NcGalaxyShapeFactorData *s_data = g_ptr_array_index (self->shape_data, gal_i);
    NcmVector *z_nodes              = g_ptr_array_index (self->z_nodes, gal_i);
    const gdouble int_pos           = nc_galaxy_position_factor_integrand_eval (self->integ_pos_lin, s_data->pos_data);
    const gdouble *w, *z, *norms;
    const guint n_total = nc_galaxy_wl_obs_cheb_grid_peek_galaxy (self->cheb_grid, _nc_data_cluster_wl_factor_row (self, gal_i), &w, &z, &norms);
    const gdouble norm  = norms[0];
    const gdouble mass  = norms[1] + norms[2];
    gdouble m2lnP_gal_i = NC_GALAXY_LOW_PROB;

    if (n_total == 0)
    {
      /* Fully in front of the lens: the shape likelihood is the constant at the anchor */
      nc_galaxy_shape_factor_eval_at_nodes_range (self->shape_factor, mset, s_data, z_nodes, 0, 1, vals);

      if (ncm_vector_get (vals, 0) * norm * int_pos > 0.0)
      {
        m2lnP_gal_i = -2.0 * log (ncm_vector_get (vals, 0) * norm * int_pos);
      }
      else
      {
        nc_galaxy_shape_factor_eval_ln_at_nodes_range (self->shape_factor, mset, s_data, z_nodes, 0, 1, vals);

        if (gsl_finite (ncm_vector_get (vals, 0)) && (norm * int_pos > 0.0))
          m2lnP_gal_i = -2.0 * (log (norm * int_pos) + ncm_vector_get (vals, 0));
        else
          self->low_prob_count++;
      }
    }
    else
    {
      gboolean converged;
      gdouble s, relerr, num;

      num = _nc_data_cluster_wl_factor_cheb_levels (self, mset, s_data, w, norms[2], z_nodes, vals, FALSE, &s, &converged, &relerr);

      /* The shape likelihood underflows at every node: the same levels on ln F, rescaled */
      if (!((mass > 0.0) && (norm * num * int_pos / mass > 0.0)))
        num = _nc_data_cluster_wl_factor_cheb_levels (self, mset, s_data, w, norms[2], z_nodes, vals, TRUE, &s, &converged, &relerr);

      if (!converged)
      {
        self->cheb_unconverged++;
        self->cheb_worst_relerr = MAX (self->cheb_worst_relerr, relerr);
      }

      if ((mass > 0.0) && gsl_finite (s) && (norm * num * int_pos / mass > 0.0))
        m2lnP_gal_i = -2.0 * (log (norm * num * int_pos / mass) + s);
      else
        self->low_prob_count++;
    }

    if (!gsl_finite (m2lnP_gal_i))
    {
      g_warning ("_nc_data_cluster_wl_factor_eval_m2lnP_cheb: galaxy %u has undefined likelihood [%g]. Skipping it.", gal_i, m2lnP_gal_i);

      if (m2lnP_gal != NULL)
        ncm_vector_set (m2lnP_gal, gal_i, GSL_NAN);

      continue;
    }

    if (m2lnP_gal != NULL)
      ncm_vector_set (m2lnP_gal, gal_i, m2lnP_gal_i);

    result += m2lnP_gal_i;
  }

  ncm_vector_free (vals);

  return result;
}

/* One panel of the CUBATURE integral over [zmin, zmax], as -2lnP. */
static gdouble
_nc_data_cluster_wl_factor_cubature_integrate (NcDataClusterWLFactorPrivate * const self, gdouble zmin, gdouble zmax)
{
  gdouble P_gal;

  ncm_vector_fast_set (self->cub_zpi, 0, zmin);
  ncm_vector_fast_set (self->cub_zpf, 0, zmax);

  ncm_integral_nd_eval (self->cub_intnd, self->cub_zpi, self->cub_zpf, self->cub_res, self->cub_err);

  P_gal = ncm_vector_fast_get (self->cub_res, 0);

  if (P_gal > 0.0)
    return -2.0 * log (P_gal);

  self->low_prob_count++;

  return NC_GALAXY_LOW_PROB;
}

static gdouble
_nc_data_cluster_wl_factor_eval_m2lnP_cubature (NcDataClusterWLFactor *dcwlf, NcmMSet *mset, NcmVector *m2lnP_gal)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);
  NcmData *data                             = NCM_DATA (dcwlf);
  const gboolean use_bootstrap              = ncm_data_bootstrap_enabled (data);
  NcmBootstrap *bstrap                      = use_bootstrap ? ncm_data_peek_bootstrap (data) : NULL;
  const guint n_iter                        = use_bootstrap ? ncm_bootstrap_get_bsize (bstrap) : self->shape_data->len;
  gdouble result                            = 0.0;
  guint gal_i, i;

  self->low_prob_count = 0;

  if (m2lnP_gal != NULL)
    g_assert_cmpuint (ncm_vector_len (m2lnP_gal), ==, self->len);

  for (i = 0; i < n_iter; i++)
  {
    gal_i = use_bootstrap ? ncm_bootstrap_get (bstrap, i) : i;

    NcGalaxyShapeFactorData *s_data = g_ptr_array_index (self->shape_data, gal_i);
    gdouble zpi, zpf, m2lnP_gal_i;

    self->cub_arg->s_data = s_data;

    nc_galaxy_redshift_factor_get_integ_lim (self->redshift_factor, mset, s_data->z_data, &zpi, &zpf);

    /* Split at z_cl when it falls strictly inside the support, exactly like
     * the LNINT path above. */
    if ((self->z_cl > zpi) && (self->z_cl < zpf))
    {
      const gdouble m2_lo = _nc_data_cluster_wl_factor_cubature_integrate (self, zpi, self->z_cl);
      const gdouble m2_hi = _nc_data_cluster_wl_factor_cubature_integrate (self, self->z_cl, zpf);

      m2lnP_gal_i = _nc_data_cluster_wl_factor_integ_combine (m2_lo, m2_hi);
    }
    else
    {
      m2lnP_gal_i = _nc_data_cluster_wl_factor_cubature_integrate (self, zpi, zpf);
    }

    if (!gsl_finite (m2lnP_gal_i))
    {
      g_warning ("_nc_data_cluster_wl_factor_eval_m2lnP_cubature: galaxy %d has undefined likelihood [%g]. Skipping it.", gal_i, m2lnP_gal_i);

      if (m2lnP_gal != NULL)
        ncm_vector_set (m2lnP_gal, gal_i, GSL_NAN);

      continue;
    }

    if (m2lnP_gal != NULL)
      ncm_vector_set (m2lnP_gal, gal_i, m2lnP_gal_i);

    result += m2lnP_gal_i;
  }

  return result;
}

/* Dispatch: decided once per call, not per galaxy -- mirrors the steps[]
 * "specialize and choose once" convention above. */
static gdouble
_nc_data_cluster_wl_factor_eval_m2lnP (NcDataClusterWLFactor *dcwlf, NcmMSet *mset, NcmVector *m2lnP_gal)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  switch (self->integ_method)
  {
    case NC_DATA_CLUSTER_WL_INTEG_METHOD_LNINT:
      return _nc_data_cluster_wl_factor_eval_m2lnP_lnint (dcwlf, mset, m2lnP_gal);

    case NC_DATA_CLUSTER_WL_INTEG_METHOD_CUBATURE:
      return _nc_data_cluster_wl_factor_eval_m2lnP_cubature (dcwlf, mset, m2lnP_gal);

    case NC_DATA_CLUSTER_WL_INTEG_METHOD_CHEB_WEIGHTS:
      return _nc_data_cluster_wl_factor_eval_m2lnP_cheb (dcwlf, mset, m2lnP_gal);

    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */

      return 0.0; /* LCOV_EXCL_LINE */
  }
}

static void
_nc_data_cluster_wl_factor_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcDataClusterWLFactor *dcwlf = NC_DATA_CLUSTER_WL_FACTOR (data);

  m2lnL[0] = _nc_data_cluster_wl_factor_eval_m2lnP (dcwlf, mset, NULL);
}

static void
_nc_data_cluster_wl_factor_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng)
{
  NcDataClusterWLFactor *dcwlf              = NC_DATA_CLUSTER_WL_FACTOR (data);
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);
  guint gal_i;
  gpointer owner = g_object_get_qdata (G_OBJECT (self->obs), _nc_data_cluster_wl_factor_shared_obs_quark ());

  /* @obs registered as shared by another dataset is read by every clone of
   * it: this realization goes to a private copy of the columns instead. */
  if ((owner != NULL) && (owner != dcwlf))
    _obs_make_private (self);

  for (gal_i = 0; gal_i < self->shape_data->len; gal_i++)
  {
    NcGalaxyShapeFactorData *s_data      = g_ptr_array_index (self->shape_data, gal_i);
    NcGalaxyPositionFactorData *pos_data = s_data->pos_data;
    NcGalaxyRedshiftFactorData *z_data   = s_data->z_data;

    if (self->resample_flag & NC_DATA_CLUSTER_WL_RESAMPLE_FLAG_REDSHIFT)
      nc_galaxy_redshift_factor_gen (self->redshift_factor, mset, z_data, rng);

    if (self->resample_flag & NC_DATA_CLUSTER_WL_RESAMPLE_FLAG_POSITION)
    {
      gdouble radius;

      do {
        nc_galaxy_position_factor_gen (self->position_factor, mset, pos_data, rng);
        radius = nc_halo_position_projected_radius_from_ra_dec (self->halo_position, self->cosmo, pos_data->ra, pos_data->dec);
      } while ((radius < self->r_min) || (radius > self->r_max));
    }

    nc_galaxy_shape_factor_gen (self->shape_factor, mset, s_data, rng);
    nc_galaxy_shape_factor_data_write_row (s_data, self->obs, _nc_data_cluster_wl_factor_row (self, gal_i));
  }

  /* ncm_data_resample() calls prepare() *before* this function runs, using whatever
   * per-galaxy data was current at that point. This loop overwrote every galaxy's raw
   * data in place but touches none of the NcmModel-side hashes the change flags of
   * prepare() are keyed on, so every cache built from the pre-resample data, the
   * redshift grid included, would stay stale. Marking obs as changed makes the next
   * prepare() redo everything from the new rows, as after set_obs(). */
  self->obs_changed = TRUE;
}

static guint
_nc_data_cluster_wl_factor_get_len (NcmData *data)
{
  NcDataClusterWLFactor *dcwlf              = NC_DATA_CLUSTER_WL_FACTOR (data);
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  return _nc_data_cluster_wl_factor_n_rows (self);
}

static void
_nc_data_cluster_wl_factor_register_shared (NcmData *data, NcmSerialize *ser)
{
  NcDataClusterWLFactor *dcwlf              = NC_DATA_CLUSTER_WL_FACTOR (data);
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  /* Registers @obs so that clones share it read-only, and marks this dataset
   * as its owner: a clone's resample() then writes to a private copy of the
   * columns (see _obs_make_private()), while this dataset resamples @obs in
   * place. */
  if (self->obs != NULL)
  {
    gchar *name = g_strdup_printf ("dcwlfObs%p", (void *) dcwlf);

    ncm_serialize_set (ser, self->obs, name, FALSE);
    g_object_set_qdata (G_OBJECT (self->obs), _nc_data_cluster_wl_factor_shared_obs_quark (), dcwlf);
    g_free (name);
  }

  /* Registers the tables of an #NcGalaxyShapeFactorMomentsTilt shape factor
   * as shared, as @obs is above, so that every worker's clone references
   * the built tables instead of serializing and rebuilding them (see
   * nc_galaxy_shape_factor_moments_tilt_register_shared()). The tables are
   * read-only, one Chebyshev table per distinct (population, shape
   * dispersion, panel count), and do not depend on the model parameters an
   * ESMCMC/MC worker samples. The other #NcGalaxyShapeFactor implementations
   * have no such cache and are deep-copied per clone. */
  if (NC_IS_GALAXY_SHAPE_FACTOR_MOMENTS_TILT (self->shape_factor))
    nc_galaxy_shape_factor_moments_tilt_register_shared (NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT (self->shape_factor), ser);
}

static void
nc_data_cluster_wl_factor_class_init (NcDataClusterWLFactorClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);

  object_class->set_property = &_nc_data_cluster_wl_factor_set_property;
  object_class->get_property = &_nc_data_cluster_wl_factor_get_property;
  object_class->dispose      = &_nc_data_cluster_wl_factor_dispose;
  object_class->finalize     = &_nc_data_cluster_wl_factor_finalize;
  object_class->constructed  = &_nc_data_cluster_wl_factor_constructed;

  /**
   * NcDataClusterWLFactor:obs:
   *
   * Galaxy weak lensing observables.
   */
  g_object_class_install_property (object_class,
                                   PROP_OBS,
                                   g_param_spec_object ("obs",
                                                        NULL,
                                                        "Galaxy weak lensing observables",
                                                        NC_TYPE_GALAXY_WL_OBS,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWLFactor:position-factor:
   *
   * The #NcGalaxyPositionFactor calculator, held and shared across
   * likelihoods -- construct-only, since this class tracks its identity as
   * fixed for the lifetime of the instance (see the class documentation).
   */
  g_object_class_install_property (object_class,
                                   PROP_POSITION_FACTOR,
                                   g_param_spec_object ("position-factor",
                                                        NULL,
                                                        "Galaxy position factor calculator",
                                                        NC_TYPE_GALAXY_POSITION_FACTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWLFactor:redshift-factor:
   *
   * The #NcGalaxyRedshiftFactor calculator, held and shared across
   * likelihoods -- construct-only (see :position-factor).
   */
  g_object_class_install_property (object_class,
                                   PROP_REDSHIFT_FACTOR,
                                   g_param_spec_object ("redshift-factor",
                                                        NULL,
                                                        "Galaxy redshift factor calculator",
                                                        NC_TYPE_GALAXY_REDSHIFT_FACTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWLFactor:shape-factor:
   *
   * The #NcGalaxyShapeFactor calculator, held and shared across
   * likelihoods -- construct-only (see :position-factor).
   */
  g_object_class_install_property (object_class,
                                   PROP_SHAPE_FACTOR,
                                   g_param_spec_object ("shape-factor",
                                                        NULL,
                                                        "Galaxy shape factor calculator",
                                                        NC_TYPE_GALAXY_SHAPE_FACTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWLFactor:r-min:
   *
   * Minimum radius of the weak lensing observables, enforced only in
   * resample()'s position rejection-sampling loop -- there is no fit-time
   * weighting.
   */
  g_object_class_install_property (object_class,
                                   PROP_R_MIN,
                                   g_param_spec_double ("r-min",
                                                        NULL,
                                                        "Minimum radius of the weak lensing observables",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWLFactor:r-max:
   *
   * Maximum radius of the weak lensing observables (see :r-min).
   */
  g_object_class_install_property (object_class,
                                   PROP_R_MAX,
                                   g_param_spec_double ("r-max",
                                                        NULL,
                                                        "Maximum radius of the weak lensing observables",
                                                        0.0, G_MAXDOUBLE, 5.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWLFactor:prec:
   *
   * Relative tolerance for the z-integral.
   */
  g_object_class_install_property (object_class,
                                   PROP_PREC,
                                   g_param_spec_double ("prec",
                                                        NULL,
                                                        "Precision for the z-integral",
                                                        0.0, 1.0, 1.0e-6,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWLFactor:len:
   *
   * Number of galaxies in the catalog.
   */
  g_object_class_install_property (object_class,
                                   PROP_LEN,
                                   g_param_spec_uint ("len",
                                                      NULL,
                                                      "Number of galaxies",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWLFactor:resample-flag:
   *
   * Flags selecting which per-galaxy quantities resample() regenerates.
   */
  g_object_class_install_property (object_class,
                                   PROP_RESAMPLE_FLAG,
                                   g_param_spec_flags ("resample-flag",
                                                       NULL,
                                                       "Resample flag",
                                                       NC_TYPE_DATA_CLUSTER_WL_RESAMPLE_FLAG,
                                                       NC_DATA_CLUSTER_WL_RESAMPLE_FLAG_ALL,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWLFactor:integ-method:
   *
   * Integration method for the redshift integral: `CHEB_WEIGHTS` (the default), `LNINT` or
   * `CUBATURE`.
   */
  g_object_class_install_property (object_class,
                                   PROP_INTEG_METHOD,
                                   g_param_spec_enum ("integ-method",
                                                      NULL,
                                                      "Integration method for the redshift integral",
                                                      NC_TYPE_DATA_CLUSTER_WL_INTEG_METHOD,
                                                      NC_DATA_CLUSTER_WL_INTEG_METHOD_CHEB_WEIGHTS,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWLFactor:cheb-tol:
   *
   * The tolerance of the nested Chebyshev-Lobatto levels under `CHEB_WEIGHTS` (ignored
   * otherwise): a galaxy stops at the first level whose background integral differs from the
   * one of the level before by at most this value times the numerator of its marginal.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHEB_TOL,
                                   g_param_spec_double ("cheb-tol",
                                                        NULL,
                                                        "Tolerance of the nested Chebyshev levels (CHEB_WEIGHTS only)",
                                                        GSL_DBL_EPSILON, 1.0, 1.0e-8,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWLFactor:cheb-level-min:
   *
   * The level whose 2^L + 1 nodes `CHEB_WEIGHTS` evaluates in its first call to the shape
   * factor (ignored otherwise), from 2 to 7. The levels above it are evaluated one call each,
   * adaptively. A shape factor whose nodes each cost the same, interpolants or one integral
   * per node, runs fastest from 2, the default. With the line table of
   * #NcGalaxyShapeFactorSpectral:line-abstol a call costs about the same whatever its number
   * of nodes, as long as it has enough of them to build the table: 5 (33 nodes) is then about
   * three times faster than 2, and 7 slightly slower than 5.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHEB_LEVEL_MIN,
                                   g_param_spec_uint ("cheb-level-min",
                                                      NULL,
                                                      "First level of the nested Chebyshev levels (CHEB_WEIGHTS only)",
                                                      2, 7, 2,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWLFactor:mask:
   *
   * Whether the likelihood uses each row of #NcDataClusterWLFactor:obs, as an
   * array of booleans ("ab") with one entry per row, or %NULL for every row
   * (see nc_data_cluster_wl_factor_set_mask()).
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_MASK,
                                   g_param_spec_variant ("mask",
                                                         NULL,
                                                         "Whether the likelihood uses each row of the catalog",
                                                         G_VARIANT_TYPE ("ab"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->bootstrap       = TRUE;
  data_class->resample        = &_nc_data_cluster_wl_factor_resample;
  data_class->m2lnL_val       = &_nc_data_cluster_wl_factor_m2lnL_val;
  data_class->get_length      = &_nc_data_cluster_wl_factor_get_len;
  data_class->prepare         = &_nc_data_cluster_wl_factor_prepare;
  data_class->register_shared = &_nc_data_cluster_wl_factor_register_shared;
}

/**
 * nc_data_cluster_wl_factor_new:
 * @position_factor: a #NcGalaxyPositionFactor
 * @redshift_factor: a #NcGalaxyRedshiftFactor
 * @shape_factor: a #NcGalaxyShapeFactor
 *
 * Creates a new #NcDataClusterWLFactor holding @position_factor,
 * @redshift_factor and @shape_factor (construct-only: these three
 * cannot be replaced after construction).
 *
 * Returns: (transfer full): a new #NcDataClusterWLFactor.
 */
NcDataClusterWLFactor *
nc_data_cluster_wl_factor_new (NcGalaxyPositionFactor *position_factor, NcGalaxyRedshiftFactor *redshift_factor, NcGalaxyShapeFactor *shape_factor)
{
  return g_object_new (NC_TYPE_DATA_CLUSTER_WL_FACTOR,
                       "position-factor", position_factor,
                       "redshift-factor", redshift_factor,
                       "shape-factor", shape_factor,
                       NULL);
}

/**
 * nc_data_cluster_wl_factor_ref:
 * @dcwlf: a #NcDataClusterWLFactor
 *
 * Increases the reference count of @dcwlf by one.
 *
 * Returns: (transfer full): @dcwlf.
 */
NcDataClusterWLFactor *
nc_data_cluster_wl_factor_ref (NcDataClusterWLFactor *dcwlf)
{
  return g_object_ref (dcwlf);
}

/**
 * nc_data_cluster_wl_factor_free:
 * @dcwlf: a #NcDataClusterWLFactor
 *
 * Decreases the reference count of @dcwlf by one.
 */
void
nc_data_cluster_wl_factor_free (NcDataClusterWLFactor *dcwlf)
{
  g_object_unref (dcwlf);
}

/**
 * nc_data_cluster_wl_factor_clear:
 * @dcwlf: a #NcDataClusterWLFactor
 *
 * Decreases the reference count of *@dcwlf by one, and sets the pointer
 * *@dcwlf to NULL.
 */
void
nc_data_cluster_wl_factor_clear (NcDataClusterWLFactor **dcwlf)
{
  g_clear_object (dcwlf);
}

/* Builds NcmData's :desc from current config + each Factor's own get_desc();
 * refreshed by every setter below. */
static gchar *
_nc_data_cluster_wl_factor_build_desc (NcDataClusterWLFactor *dcwlf)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);
  const GEnumValue *integ_val               = ncm_cfg_enum_get_value (NC_TYPE_DATA_CLUSTER_WL_INTEG_METHOD, self->integ_method);
  gchar *pos_desc                           = nc_galaxy_position_factor_get_desc (self->position_factor);
  gchar *z_desc                             = nc_galaxy_redshift_factor_get_desc (self->redshift_factor);
  gchar *shape_desc                         = nc_galaxy_shape_factor_get_desc (self->shape_factor);
  GString *str                              = g_string_new (NULL);

  g_string_append_printf (str,
                          "NcDataClusterWLFactor: %u galaxies, r in [%.4g, %.4g] Mpc, prec=%.2g. "
                          "position: %s. redshift: %s. shape: %s. integ_method=%s",
                          self->len, self->r_min, self->r_max, self->prec,
                          pos_desc, z_desc, shape_desc, integ_val->value_nick);

  if (self->integ_method == NC_DATA_CLUSTER_WL_INTEG_METHOD_CHEB_WEIGHTS)
    g_string_append_printf (str, " (cheb_tol=%.2g, cheb_level_min=%u)", self->cheb_tol, self->cheb_level_min);

  g_free (pos_desc);
  g_free (z_desc);
  g_free (shape_desc);

  return g_string_free (str, FALSE);
}

/**
 * nc_data_cluster_wl_factor_set_prec:
 * @dcwlf: a #NcDataClusterWLFactor
 * @prec: relative tolerance for the z-integral
 *
 * Sets the relative tolerance for the z-integral.
 */
void
nc_data_cluster_wl_factor_set_prec (NcDataClusterWLFactor *dcwlf, gdouble prec)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  self->prec = prec;
  ncm_data_take_desc (NCM_DATA (dcwlf), _nc_data_cluster_wl_factor_build_desc (dcwlf));
}

/**
 * nc_data_cluster_wl_factor_set_obs:
 * @dcwlf: a #NcDataClusterWLFactor
 * @obs: a #NcGalaxyWLObs
 *
 * Sets the observables catalog @obs, triggering a structural rebuild of the
 * per-galaxy data cache on the next prepare() call.
 *
 * ncm_data_resample() writes the realization into @obs, except when @obs was
 * registered as shared (ncm_dataset_register_shared()) by another dataset,
 * as for the clones #NcmFitMC makes of a dataset: then this dataset first
 * switches to a copy of the columns of @obs that references its P(z)
 * splines, and nc_data_cluster_wl_factor_peek_obs() returns that copy.
 */
void
nc_data_cluster_wl_factor_set_obs (NcDataClusterWLFactor *dcwlf, NcGalaxyWLObs *obs)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  if (obs != self->obs)
    _shared_obs_owner_clear (dcwlf, self->obs);

  nc_galaxy_wl_obs_clear (&self->obs);
  self->obs         = nc_galaxy_wl_obs_ref (obs);
  self->len         = _nc_data_cluster_wl_factor_n_rows (self);
  self->obs_changed = TRUE;

  ncm_data_set_init (NCM_DATA (dcwlf), TRUE);
  ncm_data_take_desc (NCM_DATA (dcwlf), _nc_data_cluster_wl_factor_build_desc (dcwlf));
}

/**
 * nc_data_cluster_wl_factor_set_cut:
 * @dcwlf: a #NcDataClusterWLFactor
 * @r_min: minimum projected radius
 * @r_max: maximum projected radius
 *
 * Sets the radial cut enforced by resample()'s position rejection-sampling
 * loop (see the class documentation: there is no fit-time weighting).
 */
void
nc_data_cluster_wl_factor_set_cut (NcDataClusterWLFactor *dcwlf, const gdouble r_min, const gdouble r_max)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  g_assert_cmpfloat (r_min, <, r_max);

  self->r_min = r_min;
  self->r_max = r_max;
  ncm_data_take_desc (NCM_DATA (dcwlf), _nc_data_cluster_wl_factor_build_desc (dcwlf));
}

/**
 * nc_data_cluster_wl_factor_peek_obs:
 * @dcwlf: a #NcDataClusterWLFactor
 *
 * Returns: (transfer none): the observables catalog.
 */
NcGalaxyWLObs *
nc_data_cluster_wl_factor_peek_obs (NcDataClusterWLFactor *dcwlf)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  return self->obs;
}

/**
 * nc_data_cluster_wl_factor_set_mask:
 * @dcwlf: a #NcDataClusterWLFactor
 * @keep: (element-type gboolean) (nullable): whether the likelihood uses each row of the catalog
 *
 * Restricts the likelihood to the rows of the catalog (#NcDataClusterWLFactor:obs)
 * whose entry of @keep is %TRUE: only those galaxies are loaded, evaluated,
 * bootstrapped and resampled, in the order of their rows. The other rows are
 * kept in the catalog, untouched, so that one catalog and the caches of its rows
 * serve several selections. %NULL uses every row.
 *
 * @keep is copied. It has one entry per row of the catalog, checked when the
 * catalog is loaded, at the next prepare; it must keep at least one row, so an
 * empty array aborts. Setting the
 * mask, like setting the catalog, makes the next prepare rebuild every galaxy.
 *
 * The likelihood of a mask equals, bit for bit, that of the catalog made of
 * its rows alone.
 */
void
nc_data_cluster_wl_factor_set_mask (NcDataClusterWLFactor *dcwlf, GArray *keep)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);
  guint i;

  g_clear_pointer (&self->keep, g_array_unref);
  g_clear_pointer (&self->mask, g_array_unref);

  if (keep != NULL)
  {
    const guint n_rows = keep->len;

    self->keep = g_array_sized_new (FALSE, FALSE, sizeof (gboolean), n_rows);
    self->mask = g_array_new (FALSE, FALSE, sizeof (guint));

    for (i = 0; i < n_rows; i++)
    {
      const gboolean keep_i = (g_array_index (keep, gboolean, i) != FALSE);

      g_array_append_val (self->keep, keep_i);

      if (keep_i)
        g_array_append_val (self->mask, i);
    }

    if (self->mask->len == 0)
      g_error ("nc_data_cluster_wl_factor_set_mask: the mask keeps none of the %u rows.", n_rows);
  }

  self->len = _nc_data_cluster_wl_factor_n_rows (self);

  /* As set_obs(): every galaxy is rebuilt */
  self->obs_changed = TRUE;

  if (self->obs != NULL)
  {
    ncm_data_set_init (NCM_DATA (dcwlf), TRUE);
    ncm_data_take_desc (NCM_DATA (dcwlf), _nc_data_cluster_wl_factor_build_desc (dcwlf));
  }
}

/**
 * nc_data_cluster_wl_factor_get_mask:
 * @dcwlf: a #NcDataClusterWLFactor
 *
 * Gets whether the likelihood uses each row of the catalog
 * (nc_data_cluster_wl_factor_set_mask()).
 *
 * Without a mask the C function returns %NULL; the language bindings return an
 * empty list.
 *
 * Returns: (transfer full) (element-type gboolean) (nullable): a copy of the mask, one entry per row, or %NULL for every row
 */
GArray *
nc_data_cluster_wl_factor_get_mask (NcDataClusterWLFactor *dcwlf)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);
  GArray *keep;

  if (self->keep == NULL)
    return NULL;

  keep = g_array_sized_new (FALSE, FALSE, sizeof (gboolean), self->keep->len);
  g_array_append_vals (keep, self->keep->data, self->keep->len);

  return keep;
}

/**
 * nc_data_cluster_wl_factor_set_resample_flag:
 * @dcwlf: a #NcDataClusterWLFactor
 * @resample_flag: a #NcDataClusterWLResampleFlag
 *
 * Sets which per-galaxy quantities resample() regenerates.
 */
void
nc_data_cluster_wl_factor_set_resample_flag (NcDataClusterWLFactor *dcwlf, NcDataClusterWLResampleFlag resample_flag)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  if (resample_flag & ~NC_DATA_CLUSTER_WL_RESAMPLE_FLAG_ALL)
    g_error ("nc_data_cluster_wl_factor_set_resample_flag: unknown resample flag %d.", resample_flag);

  self->resample_flag = resample_flag;
}

/**
 * nc_data_cluster_wl_factor_get_resample_flag:
 * @dcwlf: a #NcDataClusterWLFactor
 *
 * Returns: the #NcDataClusterWLResampleFlag.
 */
NcDataClusterWLResampleFlag
nc_data_cluster_wl_factor_get_resample_flag (NcDataClusterWLFactor *dcwlf)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  return self->resample_flag;
}

/**
 * nc_data_cluster_wl_factor_set_integ_method:
 * @dcwlf: a #NcDataClusterWLFactor
 * @integ_method: a #NcDataClusterWLIntegMethod
 *
 * Sets the redshift-integral method.
 */
void
nc_data_cluster_wl_factor_set_integ_method (NcDataClusterWLFactor *dcwlf, NcDataClusterWLIntegMethod integ_method)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  /* The node caches of the galaxies are kept only under CHEB_WEIGHTS: a change of method
   * takes the grid again, with them */
  if (integ_method != self->integ_method)
    _cheb_grid_clear (self);

  self->integ_method = integ_method;
  ncm_data_take_desc (NCM_DATA (dcwlf), _nc_data_cluster_wl_factor_build_desc (dcwlf));
}

/**
 * nc_data_cluster_wl_factor_get_integ_method:
 * @dcwlf: a #NcDataClusterWLFactor
 *
 * Returns: the #NcDataClusterWLIntegMethod.
 */
NcDataClusterWLIntegMethod
nc_data_cluster_wl_factor_get_integ_method (NcDataClusterWLFactor *dcwlf)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  return self->integ_method;
}

/**
 * nc_data_cluster_wl_factor_get_low_prob_count:
 * @dcwlf: a #NcDataClusterWLFactor
 *
 * Number of galaxies whose marginal probability evaluated non-positive during
 * the most recent likelihood evaluation, and for which the flat
 * NC_GALAXY_LOW_PROB fallback was therefore substituted in place of
 * -2ln(P_gal).
 *
 * Under CHEB_WEIGHTS a galaxy whose shape likelihood underflows to zero is
 * evaluated again in log space and not counted; only a genuine zero is
 * counted, such as a grid with no P(z) mass. CUBATURE integrates in linear
 * space and substitutes the fallback on underflow.
 *
 * Reset to zero at the start of every evaluation, so it always describes the
 * last one.
 *
 * Returns: the number of low-probability substitutions in the last evaluation.
 */
guint
nc_data_cluster_wl_factor_get_low_prob_count (NcDataClusterWLFactor *dcwlf)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  return self->low_prob_count;
}

/**
 * nc_data_cluster_wl_factor_data_prepare:
 * @dcwlf: a #NcDataClusterWLFactor
 * @mset: a #NcmMSet
 *
 * Runs every per-galaxy precomputation of this dataset, whether or not the
 * staleness checks of prepare() would run it: the CHEB_WEIGHTS redshift grid, when the
 * catalog stores none, and what the shape factor builds per galaxy through
 * nc_galaxy_shape_factor_data_prepare(). A shape factor may store its results in
 * properties, as #NcGalaxyShapeFactorMomentsTilt:tables does, so serializing it afterwards
 * carries them.
 *
 * Galaxies are prepared in parallel with OpenMP for every position,
 * redshift and shape factor, with the thread count taken from
 * `OMP_NUM_THREADS`. The result does not depend on the thread count. Called
 * from inside an existing parallel region, it runs serially.
 *
 * Calling it is not required: everything it builds is also built on demand
 * by prepare(), with the same results.
 *
 */
void
nc_data_cluster_wl_factor_data_prepare (NcDataClusterWLFactor *dcwlf, NcmMSet *mset)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  self->force_data_prepare = TRUE;

  ncm_data_prepare (NCM_DATA (dcwlf), mset);
}

/**
 * nc_data_cluster_wl_factor_write_shape_interp:
 * @dcwlf: a #NcDataClusterWLFactor
 * @mset: a #NcmMSet
 * @obs: a #NcGalaxyWLObs with the rows of the catalog of @dcwlf
 *
 * Prepares the per-galaxy data with nc_data_cluster_wl_factor_data_prepare(), builds the
 * interpolant of every galaxy with nc_galaxy_shape_factor_spectral_build_interp(), in
 * parallel, writes them into @obs (#NcGalaxyWLObs:shape-interp) and records the shape factor
 * and the population of @mset that built them
 * (nc_galaxy_shape_factor_spectral_set_obs_source()). The columns of @obs are not written;
 * @obs may be the catalog of @dcwlf (nc_data_cluster_wl_factor_peek_obs()). A fit that reads
 * the rows of @obs uses the interpolants; it never builds them. The free parameters of the
 * population are variables of the interpolants, over their bounds.
 *
 * The shape factor must be a #NcGalaxyShapeFactorSpectral and @obs must have as many rows as
 * the catalog of @dcwlf; otherwise the call aborts.
 */
void
nc_data_cluster_wl_factor_write_shape_interp (NcDataClusterWLFactor *dcwlf, NcmMSet *mset, NcGalaxyWLObs *obs)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);
  NcmModel *pop                             = ncm_mset_peek (mset, nc_galaxy_shape_pop_id ());
  NcGalaxyShapeFactorSpectral *gsfs;
  gint gal_s;
  guint gal_i;

  if (!NC_IS_GALAXY_SHAPE_FACTOR_SPECTRAL (self->shape_factor))
    g_error ("nc_data_cluster_wl_factor_write_shape_interp: the shape factor %s is not a NcGalaxyShapeFactorSpectral.",
             G_OBJECT_TYPE_NAME (self->shape_factor));

  gsfs = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (self->shape_factor);
  nc_data_cluster_wl_factor_data_prepare (dcwlf, mset);

  if (ncm_catalog_len (NCM_CATALOG (obs)) != ncm_catalog_len (NCM_CATALOG (self->obs)))
    g_error ("nc_data_cluster_wl_factor_write_shape_interp: the catalog has %u rows, that of the dataset %u.",
             ncm_catalog_len (NCM_CATALOG (obs)), ncm_catalog_len (NCM_CATALOG (self->obs)));

  #pragma omp parallel for schedule(dynamic, 1) if (!omp_in_parallel ())

  for (gal_s = 0; gal_s < (gint) self->shape_data->len; gal_s++)
  {
    NcGalaxyShapeFactorData *s_data = g_ptr_array_index (self->shape_data, gal_s);
    gdouble z_lo, z_hi;

    nc_galaxy_redshift_factor_get_integ_lim (self->redshift_factor, mset, s_data->z_data, &z_lo, &z_hi);
    nc_galaxy_shape_factor_spectral_build_interp (gsfs, mset, s_data, z_hi);
  }

  /* Each galaxy goes to its own row: with a mask the other rows are left as they are */
  for (gal_i = 0; gal_i < self->shape_data->len; gal_i++)
  {
    NcGalaxyShapeFactorData *s_data = g_ptr_array_index (self->shape_data, gal_i);

    s_data->ldata_write_row (s_data, obs, _nc_data_cluster_wl_factor_row (self, gal_i));
  }

  nc_galaxy_shape_factor_spectral_set_obs_source (gsfs, NC_GALAXY_SHAPE_POP (pop), obs);
}

/**
 * nc_data_cluster_wl_factor_write_cheb_grid:
 * @dcwlf: a #NcDataClusterWLFactor
 * @mset: a #NcmMSet
 * @obs: a #NcGalaxyWLObs with the rows of the catalog of @dcwlf
 *
 * Builds the CHEB_WEIGHTS redshift grid of every galaxy with
 * nc_data_cluster_wl_factor_data_prepare(), in parallel, and stores it in @obs
 * (#NcGalaxyWLObs:cheb-grid) with the cluster redshift of @mset, the type of the redshift
 * factor and the data hash of the rows of @obs. A dataset on @obs under CHEB_WEIGHTS then
 * takes the grid from the catalog instead of building it. A grid already stored in @obs is
 * removed first; @obs may be the catalog of @dcwlf (nc_data_cluster_wl_factor_peek_obs()).
 *
 * The integration method of @dcwlf must be %NC_DATA_CLUSTER_WL_INTEG_METHOD_CHEB_WEIGHTS, its
 * redshift factor must not depend on the model, it must have no mask
 * (#NcDataClusterWLFactor:mask), and @obs must have as many rows as the catalog of @dcwlf;
 * otherwise the call aborts. A dataset with a mask takes the rows of the stored grid it keeps.
 */
void
nc_data_cluster_wl_factor_write_cheb_grid (NcDataClusterWLFactor *dcwlf, NcmMSet *mset, NcGalaxyWLObs *obs)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  if (self->integ_method != NC_DATA_CLUSTER_WL_INTEG_METHOD_CHEB_WEIGHTS)
    g_error ("nc_data_cluster_wl_factor_write_cheb_grid: the integration method is not CHEB_WEIGHTS.");

  if (nc_galaxy_redshift_factor_get_hash (self->redshift_factor) != 0)
    g_error ("nc_data_cluster_wl_factor_write_cheb_grid: the redshift factor %s depends on the model.",
             G_OBJECT_TYPE_NAME (self->redshift_factor));

  if (self->keep != NULL)
    g_error ("nc_data_cluster_wl_factor_write_cheb_grid: the dataset has a mask; the grid of a catalog "
             "covers every row, so write it from a dataset without one.");

  nc_galaxy_wl_obs_set_cheb_grid (obs, NULL);
  _cheb_grid_clear (self);
  nc_data_cluster_wl_factor_data_prepare (dcwlf, mset);

  if (ncm_catalog_len (NCM_CATALOG (obs)) != self->shape_data->len)
    g_error ("nc_data_cluster_wl_factor_write_cheb_grid: the catalog has %u rows, the dataset %u galaxies.",
             ncm_catalog_len (NCM_CATALOG (obs)), self->shape_data->len);

  nc_galaxy_wl_obs_set_cheb_grid (obs, self->cheb_grid);
}

/**
 * nc_data_cluster_wl_factor_write_shape_interp_file:
 * @dcwlf: a #NcDataClusterWLFactor
 * @mset: a #NcmMSet
 * @obs: a #NcGalaxyWLObs with the rows of the catalog of @dcwlf
 * @filename: the file of the interpolants
 *
 * For a catalog whose interpolants should not be held in memory together, writes to
 * @filename the interpolant of every galaxy, the record of the shape factor and of the
 * population of @mset that built them, and the data hash of the rows of @obs, through a
 * #NcGalaxyWLObsShapeInterp created by nc_galaxy_wl_obs_shape_interp_create_empty(), and sets
 * that object, then backed by @filename, as #NcGalaxyWLObs:shape-interp of @obs. Serializing
 * @obs then writes the name of @filename and the record, not the interpolants. The galaxies
 * are built in parallel and each interpolant is freed once it is appended to the records
 * file, so the call keeps in memory only the interpolants being built, one per thread.
 *
 * The shape factor must be a #NcGalaxyShapeFactorSpectral, @dcwlf must have no mask
 * (nc_data_cluster_wl_factor_set_mask()) and @obs must have as many rows as the catalog of
 * @dcwlf; otherwise, or when a file cannot be written, the call aborts. A dataset with a mask
 * uses the interpolants of the rows it keeps from a file written for the whole catalog.
 */
void
nc_data_cluster_wl_factor_write_shape_interp_file (NcDataClusterWLFactor *dcwlf, NcmMSet *mset, NcGalaxyWLObs *obs, const gchar *filename)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);
  NcmModel *pop                             = ncm_mset_peek (mset, nc_galaxy_shape_pop_id ());
  NcGalaxyShapeFactorSpectral *gsfs;
  NcGalaxyWLObsShapeInterp *si;
  gint gal_s;

  /* The file holds one interpolant per row of the catalog */
  if (self->mask != NULL)
    g_error ("nc_data_cluster_wl_factor_write_shape_interp_file: the dataset has a mask; the interpolants file "
             "is written from a dataset that uses every row of the catalog.");

  if (!NC_IS_GALAXY_SHAPE_FACTOR_SPECTRAL (self->shape_factor))
    g_error ("nc_data_cluster_wl_factor_write_shape_interp_file: the shape factor %s is not a NcGalaxyShapeFactorSpectral.",
             G_OBJECT_TYPE_NAME (self->shape_factor));

  gsfs = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (self->shape_factor);
  nc_data_cluster_wl_factor_data_prepare (dcwlf, mset);

  if (ncm_catalog_len (NCM_CATALOG (obs)) != self->shape_data->len)
    g_error ("nc_data_cluster_wl_factor_write_shape_interp_file: the catalog has %u rows, the dataset %u galaxies.",
             ncm_catalog_len (NCM_CATALOG (obs)), self->shape_data->len);

  si = nc_galaxy_wl_obs_shape_interp_create_empty (filename, NULL);

  #pragma omp parallel for schedule(dynamic, 1) if (!omp_in_parallel ())

  for (gal_s = 0; gal_s < (gint) self->shape_data->len; gal_s++)
  {
    NcGalaxyShapeFactorData *s_data = g_ptr_array_index (self->shape_data, gal_s);
    NcmSpectralTensor *st;
    gdouble z_lo, z_hi;

    nc_galaxy_redshift_factor_get_integ_lim (self->redshift_factor, mset, s_data->z_data, &z_lo, &z_hi);
    nc_galaxy_shape_factor_spectral_build_interp (gsfs, mset, s_data, z_hi);

    st = nc_galaxy_shape_factor_spectral_data_take_interp (s_data);

    if (st != NULL)
    {
      nc_galaxy_wl_obs_shape_interp_set_interp (si, gal_s, st, NULL);
      ncm_spectral_tensor_free (st);
    }
  }

  /* The record and the data hash are set on si by the shape factor, through obs */
  g_object_set (obs, "shape-interp", si, NULL);
  nc_galaxy_shape_factor_spectral_set_obs_source (gsfs, NC_GALAXY_SHAPE_POP (pop), obs);
  nc_galaxy_wl_obs_shape_interp_close (si, NULL);
  nc_galaxy_wl_obs_shape_interp_free (si);
}

/**
 * nc_data_cluster_wl_factor_is_data_prepared:
 * @dcwlf: a #NcDataClusterWLFactor
 *
 * Whether the per-galaxy state is present: a prepare() has built the
 * per-galaxy data.
 *
 * Returns: %TRUE when the per-galaxy state is present
 */
gboolean
nc_data_cluster_wl_factor_is_data_prepared (NcDataClusterWLFactor *dcwlf)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  return self->shape_data->len > 0;
}

/**
 * nc_data_cluster_wl_factor_peek_data_array:
 * @dcwlf: a #NcDataClusterWLFactor
 *
 * Returns: (transfer none) (element-type NcGalaxyShapeFactorData): the per-galaxy data array.
 */
GPtrArray *
nc_data_cluster_wl_factor_peek_data_array (NcDataClusterWLFactor *dcwlf)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  return self->shape_data;
}

/**
 * nc_data_cluster_wl_factor_eval_m2lnP_gal:
 * @dcwlf: a #NcDataClusterWLFactor
 * @mset: a #NcmMSet
 * @m2lnP_gal: a #NcmVector of length equal to the number of galaxies, filled in place
 *
 * Computes the per-galaxy $-2\ln P_i$ and stores them in @m2lnP_gal -- the
 * per-galaxy breakdown of the total $-2\ln L$ returned by
 * ncm_data_m2lnL_val(), primarily a parity-debugging diagnostic (see the
 * class documentation). Galaxies whose likelihood is not finite (and are
 * therefore skipped in the total) are left as NAN.
 */
void
nc_data_cluster_wl_factor_eval_m2lnP_gal (NcDataClusterWLFactor *dcwlf, NcmMSet *mset, NcmVector *m2lnP_gal)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  g_assert_nonnull (m2lnP_gal);

  ncm_data_prepare (NCM_DATA (dcwlf), mset);

  g_assert_cmpuint (ncm_vector_len (m2lnP_gal), ==, self->len);

  ncm_vector_set_all (m2lnP_gal, GSL_NAN);

  _nc_data_cluster_wl_factor_eval_m2lnP (dcwlf, mset, m2lnP_gal);
}

/**
 * nc_data_cluster_wl_factor_set_cheb_tol:
 * @dcwlf: a #NcDataClusterWLFactor
 * @cheb_tol: the tolerance
 *
 * Sets #NcDataClusterWLFactor:cheb-tol.
 */
void
nc_data_cluster_wl_factor_set_cheb_tol (NcDataClusterWLFactor *dcwlf, const gdouble cheb_tol)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  self->cheb_tol = cheb_tol;
  ncm_data_take_desc (NCM_DATA (dcwlf), _nc_data_cluster_wl_factor_build_desc (dcwlf));
}

/**
 * nc_data_cluster_wl_factor_get_cheb_tol:
 * @dcwlf: a #NcDataClusterWLFactor
 *
 * Returns: #NcDataClusterWLFactor:cheb-tol.
 */
gdouble
nc_data_cluster_wl_factor_get_cheb_tol (NcDataClusterWLFactor *dcwlf)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  return self->cheb_tol;
}

/**
 * nc_data_cluster_wl_factor_get_cheb_unconverged_count:
 * @dcwlf: a #NcDataClusterWLFactor
 *
 * Under `CHEB_WEIGHTS`, the number of galaxies of the last evaluation that reached level 7
 * (129 nodes) without meeting #NcDataClusterWLFactor:cheb-tol; their marginal is the one of
 * level 7. Reset to zero at the start of every evaluation.
 *
 * Returns: the number of galaxies not converged in the last evaluation.
 */
guint
nc_data_cluster_wl_factor_get_cheb_unconverged_count (NcDataClusterWLFactor *dcwlf)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  return self->cheb_unconverged;
}

/**
 * nc_data_cluster_wl_factor_get_cheb_worst_relerr:
 * @dcwlf: a #NcDataClusterWLFactor
 *
 * Under `CHEB_WEIGHTS`, the largest relative difference between levels 6 and 7 among the
 * galaxies of nc_data_cluster_wl_factor_get_cheb_unconverged_count(), relative to the
 * numerator of their marginal; zero when every galaxy converged.
 *
 * Returns: the worst relative difference of the last evaluation.
 */
gdouble
nc_data_cluster_wl_factor_get_cheb_worst_relerr (NcDataClusterWLFactor *dcwlf)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  return self->cheb_worst_relerr;
}

/**
 * nc_data_cluster_wl_factor_set_cheb_level_min:
 * @dcwlf: a #NcDataClusterWLFactor
 * @cheb_level_min: the first level, 2 to 7
 *
 * Sets #NcDataClusterWLFactor:cheb-level-min. A level outside 2 to 7 aborts.
 */
void
nc_data_cluster_wl_factor_set_cheb_level_min (NcDataClusterWLFactor *dcwlf, const guint cheb_level_min)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  if ((cheb_level_min < 2) || (cheb_level_min > 7))
    g_error ("nc_data_cluster_wl_factor_set_cheb_level_min: the level %u is outside 2 to 7.", cheb_level_min);

  self->cheb_level_min = cheb_level_min;
  ncm_data_take_desc (NCM_DATA (dcwlf), _nc_data_cluster_wl_factor_build_desc (dcwlf));
}

/**
 * nc_data_cluster_wl_factor_get_cheb_level_min:
 * @dcwlf: a #NcDataClusterWLFactor
 *
 * Returns: #NcDataClusterWLFactor:cheb-level-min.
 */
guint
nc_data_cluster_wl_factor_get_cheb_level_min (NcDataClusterWLFactor *dcwlf)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  return self->cheb_level_min;
}

