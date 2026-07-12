/***************************************************************************
 *            nc_data_cluster_wl_factor.c
 *
 *  Sun Jul 5 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_cluster_wl_factor.c
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
 * NcDataClusterWLFactor:
 *
 * Cluster weak-lensing likelihood built on the new per-galaxy Factor
 * calculators (#NcGalaxyPositionFactor / #NcGalaxyRedshiftFactor /
 * #NcGalaxyShapeFactor), instead of the legacy #NcGalaxySDPosition /
 * #NcGalaxySDObsRedshift / #NcGalaxySDShape classes #NcDataClusterWL is
 * built on. #NcDataClusterWL is left completely untouched: it remains the
 * parity oracle this class is validated against (see the test suite), and
 * this class is the additive sibling proving the new Factor framework can
 * actually be plugged into a fit before any legacy code is reconsidered.
 *
 * The likelihood factor per galaxy is
 * $$P(\epsilon_\mathrm{obs} \mid g) = \int_{|\chi_I|<1} \mathrm{d}^2\chi_I\,
 *   P_\mathrm{pop}(\chi_I)\, N_2\big(\epsilon_\mathrm{obs} - f_g(\chi_I);
 *   \sigma_\mathrm{noise}^2\big),$$
 * combined multiplicatively with the position and joint-redshift factors and
 * integrated once over the source redshift $z$ -- exactly
 * #NcGalaxyRedshiftFactor's documented contract ("the calculator never
 * integrates $z$ itself; the orchestrator integrates this against every
 * other $z$-dependent factor").
 *
 * Implements all three of legacy's redshift-integral methods: `LNINT`
 * (default, adaptive log-domain 1D), `FIXED_NODES` (fixed Gauss-Legendre,
 * with an optional per-galaxy `auto-nodes` calibration -- see
 * nc_data_cluster_wl_factor_set_auto_nodes()) and `CUBATURE` (adaptive
 * `NcmIntegralND` over the linear-domain product). All three are
 * mathematically exact, differing only in numerical strategy/cost. No
 * bootstrap, no Fisher-matrix support, no OpenMP parallelism (including for
 * `auto-nodes` calibration, which stays serial for the same reason), no
 * fit-time `r_min`/`r_max` weighting (matching #NcDataClusterWL's own
 * deliberate no-op there -- a previous attempt "introduced bias in the
 * estimates").
 *
 * Design principle followed throughout: all preparation work happens in
 * prepare() -- resolving models, refreshing each Factor's own caches,
 * re-preparing the z-integrands -- so m2lnL_val() and resample() only *use*
 * already-prepared state. The three Factor objects are held as direct,
 * construct-only references (confirmed the intended pattern: calculators
 * are designed to be held and shared across likelihoods); this orchestrator
 * never computes a Factor's change-detection hash itself, never calls
 * `ncm_model_state_get_pkey()`, and never knows which `NcmModel`s a given
 * Factor depends on -- it only compares the opaque `guint64` values each
 * Factor's own prepare() call produces. See the three Factor classes' own
 * documentation for the get_hash()/update_data() contract this relies on.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/data/nc_data_cluster_wl_factor.h"
#include "nc_enum_types.h"
#include "nc/background/nc_hicosmo.h"
#include "nc/lss/halo/nc_halo_position.h"
#include "ncm/integration/ncm_integral1d_ptr.h"
#include "ncm/integration/ncm_integral_nd.h"
#include "ncm/integration/ncm_integrate.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

/* Fallback -2lnP for a non-positive fixed-nodes P_gal; matches legacy's own
 * local #define at nc_data_cluster_wl.c:67 (not currently shared). */
#define NC_GALAXY_LOW_PROB 1.0e6

typedef struct _NcDataClusterWLFactorIntArg
{
  NcGalaxyPositionFactorIntegrand *integ_pos;
  NcGalaxyRedshiftFactorIntegrand *integ_z;
  NcGalaxyShapeFactorIntegrand *integ_shape;
  NcGalaxyShapeFactorData *s_data;
} NcDataClusterWLFactorIntArg;

/* Linear-domain arg for the auto-nodes calibration probes F(z)=P(z),
 * G(z)=P(eps_obs|z): the position factor is deliberately not included, since
 * it scales F*G and INT F equally and so cancels in every tolerance
 * ncm_integral_fixed_calibrate() checks. */
typedef struct _NcDataClusterWLFactorCalibArg
{
  NcGalaxyRedshiftFactorIntegrand *integ_z_lin;
  NcGalaxyShapeFactorIntegrand *integ_shape_lin;
  NcGalaxyRedshiftFactorData *z_data;
  NcGalaxyShapeFactorData *s_data;
} NcDataClusterWLFactorCalibArg;

/* Linear-domain arg for CUBATURE's product integrand int_pos*int_z*int_shape
 * -- reuses the same integ_pos_lin/integ_z_lin/integ_shape_lin twins the
 * auto-nodes calibration above already needs. */
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
 * values of z in @x -- direct translation of legacy's
 * nc_data_cluster_wl_cubature_integrand. */
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

/* gsl_function F(z) = P(z): the photo-z weight baked into the calibrated
 * quadrature nodes. Direct translation of legacy's
 * _nc_data_cluster_wl_calib_pz. */
static gdouble
_nc_data_cluster_wl_factor_calib_pz (gdouble z, gpointer user_data)
{
  NcDataClusterWLFactorCalibArg *a = (NcDataClusterWLFactorCalibArg *) user_data;

  return nc_galaxy_redshift_factor_integrand_eval (a->integ_z_lin, z, a->z_data);
}

/* gsl_function G(z) = P(eps_obs|z): the slowly-varying shape factor probed
 * at the candidate nodes. Direct translation of legacy's
 * _nc_data_cluster_wl_calib_shape. */
static gdouble
_nc_data_cluster_wl_factor_calib_shape (gdouble z, gpointer user_data)
{
  NcDataClusterWLFactorCalibArg *a = (NcDataClusterWLFactorCalibArg *) user_data;

  return nc_galaxy_shape_factor_integrand_eval (a->integ_shape_lin, z, a->s_data);
}

typedef struct _NcDataClusterWLFactorPrivate
{
  NcGalaxyWLObs *obs;
  NcGalaxyPositionFactor *position_factor;
  NcGalaxyRedshiftFactor *redshift_factor;
  NcGalaxyShapeFactor *shape_factor;

  gdouble r_min;
  gdouble r_max;
  gdouble prec;
  guint len;
  NcDataClusterWLResampleFlag resample_flag;

  /* FIXED_NODES/CUBATURE state. */
  NcDataClusterWLIntegMethod integ_method;
  guint n_nodes;
  guint rule_n;

  /* Per-galaxy auto-nodes calibration (FIXED_NODES only): instead of every
   * galaxy using the global (n_nodes, rule_n), calibrate each galaxy's own
   * minimal fixed Gauss-Legendre configuration reaching node_reltol, via
   * ncm_integral_fixed_calibrate(). Opt-in, default off (matches every
   * other orchestrator knob's own opt-in convention). */
  gboolean auto_nodes;
  gdouble node_reltol;
  guint max_total_nodes;

  /* Per-galaxy fixed-node grid: index i pairs with shape_data's element i.
   * fixed_bg_nodes[i] is NULL when galaxy i's redshift support lies entirely
   * in front of the lens (no background quadrature needed, see
   * _step_fixed_nodes_grid). */
  GPtrArray *fixed_bg_nodes;     /* element type: NcmIntegralFixed* */
  GPtrArray *z_nodes_per_galaxy; /* element type: NcmVector* */
  GArray *fixed_norm;            /* element type: gdouble */

  /* Per-galaxy background node count (n_nodes-1)*rule_n selected by the
   * auto-nodes calibration (or the global value when auto_nodes is off,
   * i.e. every galaxy the same); 0 for fully-foreground galaxies. Mirrors
   * fixed_norm's exact lifecycle (same reset points). */
  GArray *n_total_per_galaxy; /* element type: guint */

  /* Last z_cl the grid was actually built for (GSL_NAN forces the first
   * rebuild). The grid only depends on z_cl, not on halo_position as a
   * whole -- radius_hash/optzs_hash bump on *any* halo_position parameter
   * (e.g. RA/Dec in a miscentering model), which would needlessly rebuild
   * the whole per-galaxy quadrature grid on every fit iteration if used as
   * the rebuild trigger instead. Compared by plain !=, matching legacy's
   * own fixed_nodes_zcl. */
  gdouble fixed_nodes_zcl;

  /* Last n_nodes/rule_n the grid was actually built for. n_nodes/rule_n are
   * orchestrator properties, not NcmModel parameters -- changing them via
   * set_n_nodes()/set_rule_n() bumps no pkey at all, so without this
   * dedicated tracking fixed_grid_changed would never notice a resize,
   * leaving z_nodes_per_galaxy/fixed_bg_nodes/shape_at_nodes sized to a
   * stale panel count (found to corrupt the heap via
   * ncm_integral_fixed_integ_vec_mult on a resized subvector). 0 is not a
   * valid n_nodes/rule_n (properties are bounded to >= 2/>= 1), so it
   * forces the first rebuild same as fixed_nodes_zcl's GSL_NAN. */
  guint fixed_nodes_n_nodes_seen;
  guint fixed_nodes_rule_n_seen;

  /* Same mid-run-change-detection mechanism as fixed_nodes_n_nodes_seen/
   * fixed_nodes_rule_n_seen above, for the three auto-nodes properties
   * (also orchestrator properties bumping no pkey). */
  gboolean fixed_nodes_auto_nodes_seen;
  gdouble fixed_nodes_node_reltol_seen;
  guint fixed_nodes_max_total_nodes_seen;

  /* Linear-domain twins of integ_pos/integ_z/integ_shape below, used by the
   * fixed-nodes path's control-variate combination (integ_pos_lin) and by
   * the auto-nodes calibration probes / CUBATURE's product integrand
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

  /* Last-seen hashes, compared against each Factor's own get_hash() every
   * prepare() cycle. */
  guint64 pos_hash;
  guint64 z_hash;
  guint64 shape_radius_hash;
  guint64 shape_optzs_hash;
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
   * length-1 vectors reused across every ncm_integral_nd_eval() call,
   * matching legacy's own CubatureIntegrator. */
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
  PROP_N_NODES,
  PROP_RULE_N,
  PROP_AUTO_NODES,
  PROP_NODE_RELTOL,
  PROP_MAX_TOTAL_NODES,
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
 * FIXED_NODES-only steps below (grid rebuild needs mset; grid/crit index
 * into the per-galaxy fixed-node arrays by gal_i) -- the five LNINT-path
 * steps ignore both, kept uniform across all steps for one shared loop. */
typedef void (*NcDataClusterWLFactorStep) (NcDataClusterWLFactorPrivate *self, NcmMSet *mset, NcGalaxyShapeFactorData *s_data, guint gal_i);

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

/* Rebuilds galaxy gal_i's fixed-node grid (anchor at index 0, background
 * Gauss-Legendre nodes after) -- direct translation of legacy's per-galaxy
 * construction body (nc_data_cluster_wl.c:1296-1337). The foreground bulk is
 * carried analytically through fixed_norm, so no quadrature interval ever
 * crosses the non-smooth point at z_cl.
 *
 * When auto_nodes is on, (n_nodes_i, rule_n_i) are calibrated per galaxy via
 * ncm_integral_fixed_calibrate() instead of using the global (n_nodes,
 * rule_n) for every galaxy -- F(z)=P(z), G(z)=P(eps_obs|z), both linear
 * domain (see integ_z_lin/integ_shape_lin's own docs). exact_bg_norm is the
 * exact background P(z) mass, usable as the calibration's missed-mass guard
 * baseline only when the galaxy's support doesn't straddle z_cl (only then
 * is the background mass analytically the full norm); otherwise
 * calibrate() falls back to its own high-resolution reference mass.
 * calibrate()'s own returned NcmIntegralFixed is discarded and the grid is
 * rebuilt via make_fixed_nodes() at the chosen (n_nodes_i, rule_n_i): that
 * virtual method may bake in a scheme-specific construction (e.g. bypassing
 * the generic integrand-eval path) that calibrate()'s generic F closure
 * cannot reproduce, so reusing its result would risk a different answer
 * than the non-auto-nodes path gives at the same (n_nodes_i, rule_n_i). */
static void
_step_fixed_nodes_grid (NcDataClusterWLFactorPrivate *self, NcmMSet *mset, NcGalaxyShapeFactorData *s_data, guint gal_i)
{
  guint n_nodes_i = self->n_nodes;
  guint rule_n_i  = self->rule_n;
  guint n_total_i = 0;
  gdouble z_lo, z_hi, norm;
  gboolean has_bg;
  NcmIntegralFixed *bg_intf;
  NcmVector *z_nodes;

  nc_galaxy_redshift_factor_get_integ_lim (self->redshift_factor, mset, s_data->z_data, &z_lo, &z_hi);
  norm   = nc_galaxy_redshift_factor_norm (self->redshift_factor, mset, s_data->z_data);
  has_bg = z_hi > self->z_cl;

  if (has_bg)
  {
    const gdouble bg_lo = MAX (z_lo, self->z_cl);
    NcmVector *z_nodes_bg;

    if (self->auto_nodes)
    {
      NcDataClusterWLFactorCalibArg calib_arg = { self->integ_z_lin, self->integ_shape_lin, s_data->z_data, s_data };
      gsl_function F                          = { &_nc_data_cluster_wl_factor_calib_pz, &calib_arg };
      gsl_function G                          = { &_nc_data_cluster_wl_factor_calib_shape, &calib_arg };
      const gdouble exact_bg_norm             = (z_lo >= self->z_cl) ? norm : GSL_NAN;
      NcmIntegralFixed *cal                   = ncm_integral_fixed_calibrate (&F, &G, bg_lo, z_hi, self->node_reltol, exact_bg_norm, self->max_total_nodes, &n_nodes_i, &rule_n_i);

      ncm_integral_fixed_free (cal);
    }

    n_total_i  = (n_nodes_i - 1) * rule_n_i;
    bg_intf    = nc_galaxy_redshift_factor_make_fixed_nodes (self->redshift_factor, mset, s_data->z_data, bg_lo, z_hi, n_nodes_i, rule_n_i);
    z_nodes    = ncm_vector_new (n_total_i + 1);
    z_nodes_bg = ncm_vector_get_subvector (z_nodes, 1, n_total_i);

    ncm_integral_fixed_get_nodes (bg_intf, z_nodes_bg);
    ncm_vector_free (z_nodes_bg);
  }
  else
  {
    bg_intf = NULL;
    z_nodes = ncm_vector_new (1);
  }

  ncm_vector_set (z_nodes, 0, z_lo);

  g_ptr_array_index (self->fixed_bg_nodes, gal_i)        = bg_intf;
  g_ptr_array_index (self->z_nodes_per_galaxy, gal_i)    = z_nodes;
  g_array_index (self->fixed_norm, gdouble, gal_i)       = norm;
  g_array_index (self->n_total_per_galaxy, guint, gal_i) = n_total_i;
}

static void
_step_fixed_nodes_crit (NcDataClusterWLFactorPrivate *self, NcmMSet *mset, NcGalaxyShapeFactorData *s_data, guint gal_i)
{
  NcmVector *z_nodes = (NcmVector *) g_ptr_array_index (self->z_nodes_per_galaxy, gal_i);

  nc_galaxy_shape_factor_update_data_at_nodes_crit (self->shape_factor, s_data, z_nodes);
}

static void
_step_fixed_nodes_sigma (NcDataClusterWLFactorPrivate *self, NcmMSet *mset, NcGalaxyShapeFactorData *s_data, guint gal_i)
{
  nc_galaxy_shape_factor_update_data_at_nodes_sigma (self->shape_factor, s_data);
}

/* fixed_bg_nodes has no GDestroyNotify (NcmIntegralFixed is not NULL-safe to
 * free, and entries are legitimately NULL for foreground-only galaxies), so
 * both disposal and grid-rebuild resizing must free non-NULL entries by hand
 * through this one shared helper. */
static void
_fixed_bg_nodes_free_entries (GPtrArray *fixed_bg_nodes)
{
  guint i;

  for (i = 0; i < fixed_bg_nodes->len; i++)
  {
    NcmIntegralFixed *bg_intf = (NcmIntegralFixed *) g_ptr_array_index (fixed_bg_nodes, i);

    if (bg_intf != NULL)
      ncm_integral_fixed_free (bg_intf);
  }
}

/* Resizes the per-galaxy fixed-node arrays to the current galaxy count,
 * freeing whatever they held before -- called once per prepare() cycle
 * (not per galaxy) right before the steps[] loop, only when
 * fixed_grid_changed. After this, every slot is NULL/0 and
 * _step_fixed_nodes_grid's indexed writes are leak-free. */
static void
_fixed_nodes_grid_reset (NcDataClusterWLFactorPrivate *self)
{
  _fixed_bg_nodes_free_entries (self->fixed_bg_nodes);
  g_ptr_array_set_size (self->fixed_bg_nodes, 0);
  g_ptr_array_set_size (self->fixed_bg_nodes, self->shape_data->len);

  g_ptr_array_set_size (self->z_nodes_per_galaxy, 0);
  g_ptr_array_set_size (self->z_nodes_per_galaxy, self->shape_data->len);

  g_array_set_size (self->fixed_norm, self->shape_data->len);

  g_array_set_size (self->n_total_per_galaxy, 0);
  g_array_set_size (self->n_total_per_galaxy, self->shape_data->len);
}

static void
nc_data_cluster_wl_factor_init (NcDataClusterWLFactor *dcwlf)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  self->obs             = NULL;
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
  self->pos_hash          = 0;
  self->z_hash            = 0;
  self->shape_radius_hash = 0;
  self->shape_optzs_hash  = 0;
  self->shape_pop_hash    = 0;

  self->cosmo         = NULL;
  self->halo_position = NULL;
  self->z_cl          = 0.0;

  self->integ_pos      = NULL;
  self->integ_z        = NULL;
  self->integ_shape    = NULL;
  self->int1d          = NULL;
  self->int_arg.s_data = NULL;

  self->integ_method    = NC_DATA_CLUSTER_WL_INTEG_METHOD_LNINT;
  self->n_nodes         = 10;
  self->rule_n          = 5;
  self->auto_nodes      = FALSE;
  self->node_reltol     = 1.0e-4;
  self->max_total_nodes = 2000;

  /* No destroy-func on fixed_bg_nodes: entries are nullable (NcmIntegralFixed
   * is not NULL-safe to free), so clearing/resizing this array is handled
   * manually with an explicit NULL check -- see _fixed_nodes_grid_reset().
   * z_nodes_per_galaxy's entries are never NULL, so its destroy-func is
   * safe to rely on for both resizing and disposal. */
  self->fixed_bg_nodes                   = g_ptr_array_new ();
  self->z_nodes_per_galaxy               = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_vector_free);
  self->fixed_norm                       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->n_total_per_galaxy               = g_array_new (FALSE, FALSE, sizeof (guint));
  self->fixed_nodes_zcl                  = GSL_NAN;
  self->fixed_nodes_n_nodes_seen         = 0;
  self->fixed_nodes_rule_n_seen          = 0;
  self->fixed_nodes_auto_nodes_seen      = FALSE;
  self->fixed_nodes_node_reltol_seen     = 0.0;
  self->fixed_nodes_max_total_nodes_seen = 0;

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
    case PROP_N_NODES:
      nc_data_cluster_wl_factor_set_n_nodes (dcwlf, g_value_get_uint (value));
      break;
    case PROP_RULE_N:
      nc_data_cluster_wl_factor_set_rule_n (dcwlf, g_value_get_uint (value));
      break;
    case PROP_AUTO_NODES:
      nc_data_cluster_wl_factor_set_auto_nodes (dcwlf, g_value_get_boolean (value));
      break;
    case PROP_NODE_RELTOL:
      nc_data_cluster_wl_factor_set_node_reltol (dcwlf, g_value_get_double (value));
      break;
    case PROP_MAX_TOTAL_NODES:
      nc_data_cluster_wl_factor_set_max_total_nodes (dcwlf, g_value_get_uint (value));
      break;
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
    case PROP_N_NODES:
      g_value_set_uint (value, self->n_nodes);
      break;
    case PROP_RULE_N:
      g_value_set_uint (value, self->rule_n);
      break;
    case PROP_AUTO_NODES:
      g_value_set_boolean (value, self->auto_nodes);
      break;
    case PROP_NODE_RELTOL:
      g_value_set_double (value, self->node_reltol);
      break;
    case PROP_MAX_TOTAL_NODES:
      g_value_set_uint (value, self->max_total_nodes);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_data_cluster_wl_factor_dispose (GObject *object)
{
  NcDataClusterWLFactor *dcwlf              = NC_DATA_CLUSTER_WL_FACTOR (object);
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

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

  _fixed_bg_nodes_free_entries (self->fixed_bg_nodes);
  g_clear_pointer (&self->fixed_bg_nodes, g_ptr_array_unref);
  g_clear_pointer (&self->z_nodes_per_galaxy, g_ptr_array_unref);
  g_clear_pointer (&self->fixed_norm, g_array_unref);
  g_clear_pointer (&self->n_total_per_galaxy, g_array_unref);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_data_cluster_wl_factor_parent_class)->dispose (object);
}

static void
_nc_data_cluster_wl_factor_finalize (GObject *object)
{
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

  for (i = 0; i < nc_galaxy_wl_obs_len (obs); i++)
  {
    NcGalaxyPositionFactorData *pos_data = nc_galaxy_position_factor_data_new (self->position_factor, mset);
    NcGalaxyRedshiftFactorData *z_data   = nc_galaxy_redshift_factor_data_new (self->redshift_factor, mset);
    NcGalaxyShapeFactorData *s_data      = nc_galaxy_shape_factor_data_new (self->shape_factor, mset, pos_data, z_data);

    nc_galaxy_shape_factor_data_read_row (s_data, obs, i);
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
  gboolean pos_changed, z_changed, radius_changed, optzs_changed, pop_changed, fixed_grid_changed;
  NcDataClusterWLFactorStep steps[8];
  guint n_steps = 0;

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

  /* Captured before self->obs_changed is reset to FALSE just below: the
   * pos_changed/z_changed/radius_changed/optzs_changed/pop_changed/
   * fixed_grid_changed flags computed further down all need to know whether
   * obs was replaced *this* cycle, and must read that fact through this
   * local -- reading self->obs_changed itself at that point would always
   * see FALSE, silently turning every "self->obs_changed ||" term there
   * into dead code. Found as a real, pre-existing Phase-1 bug (this
   * ordering predates the FIXED_NODES work): it was masked for the other
   * five flags because their own hash fields start at 0, which any real
   * hash is virtually certain to differ from, so "first ever" was still
   * caught by the hash comparison alone -- but not masked for
   * fixed_grid_changed, whose own comparisons (z_cl, n_nodes, rule_n) all
   * still match their last-seen values after a pure obs swap that changes
   * no NcmModel at all, so the dead obs_changed term was the only thing
   * that could have signaled "rebuild anyway." */
  const gboolean obs_was_changed = self->obs_changed;

  if (obs_was_changed)
  {
    _nc_data_cluster_wl_factor_load_obs (self, mset, self->obs);
    self->obs_changed = FALSE;
  }

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
     * fixed-nodes control-variate combination (integ_pos_lin), the
     * auto-nodes calibration probes, and CUBATURE's product integrand
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

  /* radius_changed:      galaxy geometry changed (halo_position/cosmo).
   * optzs_changed:       lens model changed (cosmo/density_profile/
   *                      mass_summary/surface_mass_density/halo_position).
   * fixed_grid_changed:  integration discretization changed -- new z-nodes
   *                      needed, only meaningful under FIXED_NODES. Deliberately
   *                      *not* driven by radius_changed: the grid only
   *                      depends on z_cl, but radius_hash bumps on *any*
   *                      halo_position parameter (e.g. RA/Dec in a
   *                      miscentering model) -- using it here would rebuild
   *                      the whole per-galaxy quadrature grid every fit
   *                      iteration even when z_cl itself never moved. Tracked
   *                      instead via the dedicated fixed_nodes_zcl scalar
   *                      compare, matching legacy's own fixed_nodes_zcl --
   *                      plus n_nodes/rule_n, orchestrator properties (not
   *                      NcmModel parameters) that resize the grid without
   *                      bumping any pkey at all; missing this crashed via
   *                      a stale-sized subvector passed to
   *                      ncm_integral_fixed_integ_vec_mult(). */
  pos_changed    = obs_was_changed || (nc_galaxy_position_factor_get_hash (self->position_factor) != self->pos_hash);
  z_changed      = obs_was_changed || (nc_galaxy_redshift_factor_get_hash (self->redshift_factor) != self->z_hash);
  radius_changed = obs_was_changed || (nc_galaxy_shape_factor_get_radius_hash (self->shape_factor) != self->shape_radius_hash);
  optzs_changed  = obs_was_changed || (nc_galaxy_shape_factor_get_optzs_hash (self->shape_factor) != self->shape_optzs_hash);
  pop_changed    = obs_was_changed || (nc_galaxy_shape_factor_get_pop_hash (self->shape_factor) != self->shape_pop_hash);

  fixed_grid_changed = obs_was_changed || z_changed || (self->z_cl != self->fixed_nodes_zcl) ||
                       (self->n_nodes != self->fixed_nodes_n_nodes_seen) || (self->rule_n != self->fixed_nodes_rule_n_seen) ||
                       (self->auto_nodes != self->fixed_nodes_auto_nodes_seen) ||
                       (self->node_reltol != self->fixed_nodes_node_reltol_seen) ||
                       (self->max_total_nodes != self->fixed_nodes_max_total_nodes_seen);

  self->pos_hash          = nc_galaxy_position_factor_get_hash (self->position_factor);
  self->z_hash            = nc_galaxy_redshift_factor_get_hash (self->redshift_factor);
  self->shape_radius_hash = nc_galaxy_shape_factor_get_radius_hash (self->shape_factor);
  self->shape_optzs_hash  = nc_galaxy_shape_factor_get_optzs_hash (self->shape_factor);
  self->shape_pop_hash    = nc_galaxy_shape_factor_get_pop_hash (self->shape_factor);

  if (!pos_changed && !z_changed && !radius_changed && !optzs_changed && !pop_changed && !fixed_grid_changed)
    return;  /* nothing to do this cycle */

  if ((self->integ_method == NC_DATA_CLUSTER_WL_INTEG_METHOD_FIXED_NODES) && fixed_grid_changed)
  {
    _fixed_nodes_grid_reset (self);
    self->fixed_nodes_zcl                  = self->z_cl;
    self->fixed_nodes_n_nodes_seen         = self->n_nodes;
    self->fixed_nodes_rule_n_seen          = self->rule_n;
    self->fixed_nodes_auto_nodes_seen      = self->auto_nodes;
    self->fixed_nodes_node_reltol_seen     = self->node_reltol;
    self->fixed_nodes_max_total_nodes_seen = self->max_total_nodes;
  }

  /* Specialize and choose once: no per-galaxy branch on any of these. The
   * three fixed-nodes steps are only ever appended under FIXED_NODES;
   * ordering matters -- _step_shape_radius before _step_fixed_nodes_sigma
   * (sigma reads the projected radius) and _step_fixed_nodes_grid before
   * _step_fixed_nodes_crit (crit reads the current z_nodes).
   *
   * optzs/pop are moved ahead of the FIXED_NODES block (rather than after
   * it, as their own hash-gating alone would otherwise suggest) because
   * auto-nodes calibration inside _step_fixed_nodes_grid evaluates the
   * continuous-z shape integrand directly, which reads cdata->optzs (from
   * update_data_optzs) and data->pop_data (from update_data_pop) -- both
   * uninitialized on a galaxy's very first cycle if left in their original
   * post-FIXED_NODES position. Moving them earlier costs nothing when
   * auto_nodes is off (grid doesn't touch either), and optzs still only
   * needs radius (already placed before this block) to be current. */
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

  if (self->integ_method == NC_DATA_CLUSTER_WL_INTEG_METHOD_FIXED_NODES)
  {
    if (fixed_grid_changed)
      steps[n_steps++] = &_step_fixed_nodes_grid;

    if (fixed_grid_changed || optzs_changed)
      steps[n_steps++] = &_step_fixed_nodes_crit;

    /* fixed_grid_changed (not just radius_changed/optzs_changed) must gate
     * sigma too: its GSL_NAN-seeded fixed_nodes_zcl sentinel is what
     * detects "FIXED_NODES has never run before on this object" (e.g. a
     * user switching integ_method mid-run after cycles spent in LNINT,
     * during which radius_hash/optzs_hash were already kept current for
     * the LNINT path's own needs and so would not, on their own, signal
     * that sigma_cache -- a FIXED_NODES-only cache -- has never actually
     * been computed). Found via exactly this switch-mid-run scenario. */
    if (fixed_grid_changed || radius_changed || optzs_changed)
      steps[n_steps++] = &_step_fixed_nodes_sigma;
  }

  {
    guint gal_i;

    for (gal_i = 0; gal_i < self->shape_data->len; gal_i++)
    {
      NcGalaxyShapeFactorData *s_data = g_ptr_array_index (self->shape_data, gal_i);
      guint j;

      for (j = 0; j < n_steps; j++)
        steps[j](self, mset, s_data, gal_i);
    }
  }
}

/* Combine two adaptive sub-integrals split at z_cl. Each panel returns
 * -2 ln I_panel; recombine in the linear domain I = I_lo + I_hi -- identical
 * to NcDataClusterWL's own _nc_data_cluster_wl_integ_combine(). */
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
  gdouble result                            = 0.0;
  guint gal_i;

  if (m2lnP_gal != NULL)
    g_assert_cmpuint (ncm_vector_len (m2lnP_gal), ==, self->len);

  for (gal_i = 0; gal_i < self->shape_data->len; gal_i++)
  {
    NcGalaxyShapeFactorData *s_data = g_ptr_array_index (self->shape_data, gal_i);
    gdouble zpi, zpf, m2lnP_gal_i;

    self->int_arg.s_data = s_data;

    nc_galaxy_redshift_factor_get_integ_lim (self->redshift_factor, mset, s_data->z_data, &zpi, &zpf);

    /* Split at z_cl when it falls strictly inside the support (the reduced
     * shear has a kink there): same reasoning as NcDataClusterWL's own
     * adaptive-quadrature path. */
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

/* Integrates the shape values at @shape_at_nodes against galaxy @gal_i's
 * precomputed P(z)-weighted background quadrature using a control-variate
 * form: P_gal = p_a * norm + Int_bg P(z) [P(e_o,z) - p_a] dz, where @p_a is
 * the anchor value (index 0) and @norm is the exact full-support P(z)
 * normalization -- algebraically exact, and never places a quadrature
 * interval across the non-smooth point at z_cl (the foreground is handled
 * analytically through @norm). Direct translation of NcDataClusterWL's own
 * _nc_data_cluster_wl_fixed_panels_integ. */
static gdouble
_nc_data_cluster_wl_factor_fixed_panels_integ (NcDataClusterWLFactorPrivate * const self, guint gal_i, NcmVector *shape_at_nodes, NcmVector *sub)
{
  NcmIntegralFixed *bg_intf = (NcmIntegralFixed *) g_ptr_array_index (self->fixed_bg_nodes, gal_i);
  const gdouble norm        = g_array_index (self->fixed_norm, gdouble, gal_i);
  const gdouble p_a         = ncm_vector_get (shape_at_nodes, 0);
  gdouble P                 = p_a * norm;

  if (bg_intf != NULL)
  {
    ncm_vector_add_constant (sub, -p_a);
    P += ncm_integral_fixed_integ_vec_mult (bg_intf, sub);
  }

  return P;
}

static gdouble
_nc_data_cluster_wl_factor_eval_m2lnP_fixed (NcDataClusterWLFactor *dcwlf, NcmMSet *mset, NcmVector *m2lnP_gal)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);
  guint max_n_total                         = 0;
  NcmVector *shape_at_nodes;
  gdouble result = 0.0;
  guint gal_i;

  /* Per-galaxy node counts can differ under auto-nodes (and are all equal
   * to the global (n_nodes-1)*rule_n otherwise): size the reusable shape
   * buffer to the largest, then take a per-galaxy leading subvector below.
   * Built once per call, not once per galaxy, to avoid a GObject
   * alloc/dispose per galaxy (ncm_vector_get_subvector is a full
   * g_object_new). */
  for (gal_i = 0; gal_i < self->n_total_per_galaxy->len; gal_i++)
    max_n_total = MAX (max_n_total, g_array_index (self->n_total_per_galaxy, guint, gal_i));

  shape_at_nodes = ncm_vector_new (max_n_total + 1);

  if (m2lnP_gal != NULL)
    g_assert_cmpuint (ncm_vector_len (m2lnP_gal), ==, self->len);

  for (gal_i = 0; gal_i < self->shape_data->len; gal_i++)
  {
    NcGalaxyShapeFactorData *s_data = g_ptr_array_index (self->shape_data, gal_i);
    NcmVector *z_nodes              = (NcmVector *) g_ptr_array_index (self->z_nodes_per_galaxy, gal_i);
    const guint n_total_i           = g_array_index (self->n_total_per_galaxy, guint, gal_i);
    const gdouble int_pos           = nc_galaxy_position_factor_integrand_eval (self->integ_pos_lin, s_data->pos_data);

    /* eval_at_nodes needs an @out exactly as long as @z_nodes (which is
     * this galaxy's own n_total_i+1, possibly shorter than shape_at_nodes'
     * own max_n_total+1 buffer under auto-nodes) -- a leading [0,
     * n_total_i+1) view, distinct from @sub's [1, n_total_i) view used only
     * by _fixed_panels_integ below. Both are non-owning views of the same
     * underlying buffer, safe to hold simultaneously. */
    NcmVector *nodes_view = (n_total_i < max_n_total) ? ncm_vector_get_subvector (shape_at_nodes, 0, n_total_i + 1) : shape_at_nodes;
    NcmVector *sub        = (n_total_i > 0) ? ncm_vector_get_subvector (shape_at_nodes, 1, n_total_i) : NULL;
    gdouble P_gal, m2lnP_gal_i;

    nc_galaxy_shape_factor_eval_at_nodes (self->shape_factor, mset, s_data, z_nodes, nodes_view);

    P_gal       = _nc_data_cluster_wl_factor_fixed_panels_integ (self, gal_i, shape_at_nodes, sub) * int_pos;
    m2lnP_gal_i = P_gal > 0.0 ? -2.0 * log (P_gal) : NC_GALAXY_LOW_PROB;

    if (nodes_view != shape_at_nodes)
      ncm_vector_free (nodes_view);

    if (sub != NULL)
      ncm_vector_free (sub);

    if (!gsl_finite (m2lnP_gal_i))
    {
      g_warning ("_nc_data_cluster_wl_factor_eval_m2lnP_fixed: galaxy %d has undefined likelihood [%g]. Skipping it.", gal_i, m2lnP_gal_i);

      if (m2lnP_gal != NULL)
        ncm_vector_set (m2lnP_gal, gal_i, GSL_NAN);

      continue;
    }

    if (m2lnP_gal != NULL)
      ncm_vector_set (m2lnP_gal, gal_i, m2lnP_gal_i);

    result += m2lnP_gal_i;
  }

  ncm_vector_free (shape_at_nodes);

  return result;
}

/* One panel of the CUBATURE integral over [zmin, zmax], as -2lnP -- direct
 * translation of legacy's cubature_integrator_integrate. */
static gdouble
_nc_data_cluster_wl_factor_cubature_integrate (NcDataClusterWLFactorPrivate * const self, gdouble zmin, gdouble zmax)
{
  gdouble P_gal;

  ncm_vector_fast_set (self->cub_zpi, 0, zmin);
  ncm_vector_fast_set (self->cub_zpf, 0, zmax);

  ncm_integral_nd_eval (self->cub_intnd, self->cub_zpi, self->cub_zpf, self->cub_res, self->cub_err);

  P_gal = ncm_vector_fast_get (self->cub_res, 0);

  return P_gal > 0.0 ? -2.0 * log (P_gal) : NC_GALAXY_LOW_PROB;
}

static gdouble
_nc_data_cluster_wl_factor_eval_m2lnP_cubature (NcDataClusterWLFactor *dcwlf, NcmMSet *mset, NcmVector *m2lnP_gal)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);
  gdouble result                            = 0.0;
  guint gal_i;

  if (m2lnP_gal != NULL)
    g_assert_cmpuint (ncm_vector_len (m2lnP_gal), ==, self->len);

  for (gal_i = 0; gal_i < self->shape_data->len; gal_i++)
  {
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
    case NC_DATA_CLUSTER_WL_INTEG_METHOD_FIXED_NODES:
      return _nc_data_cluster_wl_factor_eval_m2lnP_fixed (dcwlf, mset, m2lnP_gal);

    case NC_DATA_CLUSTER_WL_INTEG_METHOD_LNINT:
      return _nc_data_cluster_wl_factor_eval_m2lnP_lnint (dcwlf, mset, m2lnP_gal);

    case NC_DATA_CLUSTER_WL_INTEG_METHOD_CUBATURE:
      return _nc_data_cluster_wl_factor_eval_m2lnP_cubature (dcwlf, mset, m2lnP_gal);

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
    nc_galaxy_shape_factor_data_write_row (s_data, self->obs, gal_i);
  }

  /* ncm_data_resample() calls prepare() *before* this function runs, using
   * whatever per-galaxy data was current at that point (e.g. set_obs()'s
   * placeholder zeros, on a fresh object). This loop just overwrote every
   * galaxy's raw data in place, but touches none of the NcmModel-side
   * hashes that pos_changed/z_changed/radius_changed/optzs_changed/
   * pop_changed/fixed_grid_changed are keyed on -- so without this, every
   * cache built from the pre-resample data (most dangerously the
   * FIXED_NODES z-node grid, which is gated only by z_cl/n_nodes/rule_n and
   * so is *not* naturally invalidated by a subsequent mass-only fit/scan)
   * would silently stay stale for the rest of the object's life. Marking
   * obs as changed forces the next prepare() cycle to redo everything from
   * the fresh data, exactly as if set_obs() had just been called. */
  self->obs_changed = TRUE;
}

static guint
_nc_data_cluster_wl_factor_get_len (NcmData *data)
{
  NcDataClusterWLFactor *dcwlf              = NC_DATA_CLUSTER_WL_FACTOR (data);
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  return (self->obs != NULL) ? nc_galaxy_wl_obs_len (self->obs) : 0;
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
   * resample()'s position rejection-sampling loop (matching
   * #NcDataClusterWL's own semantics -- there is no fit-time weighting).
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
   * Integration method for the redshift integral: `LNINT`, `FIXED_NODES` or
   * `CUBATURE`.
   */
  g_object_class_install_property (object_class,
                                   PROP_INTEG_METHOD,
                                   g_param_spec_enum ("integ-method",
                                                      NULL,
                                                      "Integration method for the redshift integral",
                                                      NC_TYPE_DATA_CLUSTER_WL_INTEG_METHOD,
                                                      NC_DATA_CLUSTER_WL_INTEG_METHOD_LNINT,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWLFactor:n-nodes:
   *
   * Number of fixed-quadrature panels per galaxy under `FIXED_NODES`
   * (ignored otherwise).
   */
  g_object_class_install_property (object_class,
                                   PROP_N_NODES,
                                   g_param_spec_uint ("n-nodes",
                                                      NULL,
                                                      "Number of fixed-quadrature panels (FIXED_NODES only)",
                                                      2, G_MAXUINT, 10,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWLFactor:rule-n:
   *
   * Gauss-Legendre rule order per fixed-quadrature panel under
   * `FIXED_NODES` (ignored otherwise).
   */
  g_object_class_install_property (object_class,
                                   PROP_RULE_N,
                                   g_param_spec_uint ("rule-n",
                                                      NULL,
                                                      "Gauss-Legendre rule order per panel (FIXED_NODES only)",
                                                      1, G_MAXUINT, 5,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWLFactor:auto-nodes:
   *
   * Whether to automatically select, per galaxy, the minimal fixed
   * Gauss-Legendre configuration reaching #NcDataClusterWLFactor:node-reltol
   * (`FIXED_NODES` only, ignored otherwise). When disabled (the default),
   * every galaxy uses the global #NcDataClusterWLFactor:n-nodes /
   * #NcDataClusterWLFactor:rule-n.
   */
  g_object_class_install_property (object_class,
                                   PROP_AUTO_NODES,
                                   g_param_spec_boolean ("auto-nodes",
                                                         NULL,
                                                         "Automatically select the per-galaxy fixed-node configuration (FIXED_NODES only)",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWLFactor:node-reltol:
   *
   * Target relative tolerance for the per-galaxy fixed-node selection (see
   * #NcDataClusterWLFactor:auto-nodes).
   */
  g_object_class_install_property (object_class,
                                   PROP_NODE_RELTOL,
                                   g_param_spec_double ("node-reltol",
                                                        NULL,
                                                        "Target relative tolerance for the per-galaxy fixed-node selection",
                                                        0.0, G_MAXDOUBLE, 1.0e-4,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWLFactor:max-total-nodes:
   *
   * Safety ceiling on the total background node count $(n_\mathrm{nodes}-1)\,
   * \mathrm{rule}_n$ explored by the #NcDataClusterWLFactor:auto-nodes
   * selection. If the tolerance is not met within this budget the best
   * configuration found is used and a warning is emitted.
   */
  g_object_class_install_property (object_class,
                                   PROP_MAX_TOTAL_NODES,
                                   g_param_spec_uint ("max-total-nodes",
                                                      NULL,
                                                      "Safety ceiling on the total background node count for auto-nodes selection",
                                                      2, G_MAXUINT, 2000,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->bootstrap  = FALSE; /* deferred to a follow-up, see the class doc */
  data_class->resample   = &_nc_data_cluster_wl_factor_resample;
  data_class->m2lnL_val  = &_nc_data_cluster_wl_factor_m2lnL_val;
  data_class->get_length = &_nc_data_cluster_wl_factor_get_len;
  data_class->prepare    = &_nc_data_cluster_wl_factor_prepare;
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
}

/**
 * nc_data_cluster_wl_factor_set_obs:
 * @dcwlf: a #NcDataClusterWLFactor
 * @obs: a #NcGalaxyWLObs
 *
 * Sets the observables catalog @obs, triggering a structural rebuild of the
 * per-galaxy data cache on the next prepare() call.
 */
void
nc_data_cluster_wl_factor_set_obs (NcDataClusterWLFactor *dcwlf, NcGalaxyWLObs *obs)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  nc_galaxy_wl_obs_clear (&self->obs);
  self->len         = nc_galaxy_wl_obs_len (obs);
  self->obs         = nc_galaxy_wl_obs_ref (obs);
  self->obs_changed = TRUE;

  ncm_data_set_init (NCM_DATA (dcwlf), TRUE);
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

  self->integ_method = integ_method;
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
 * nc_data_cluster_wl_factor_set_n_nodes:
 * @dcwlf: a #NcDataClusterWLFactor
 * @n_nodes: number of fixed-quadrature panels
 *
 * Sets the number of fixed-quadrature panels per galaxy (`FIXED_NODES`
 * only).
 */
void
nc_data_cluster_wl_factor_set_n_nodes (NcDataClusterWLFactor *dcwlf, guint n_nodes)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  self->n_nodes = n_nodes;
}

/**
 * nc_data_cluster_wl_factor_get_n_nodes:
 * @dcwlf: a #NcDataClusterWLFactor
 *
 * Returns: the number of fixed-quadrature panels.
 */
guint
nc_data_cluster_wl_factor_get_n_nodes (NcDataClusterWLFactor *dcwlf)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  return self->n_nodes;
}

/**
 * nc_data_cluster_wl_factor_set_rule_n:
 * @dcwlf: a #NcDataClusterWLFactor
 * @rule_n: Gauss-Legendre rule order per panel
 *
 * Sets the Gauss-Legendre rule order per fixed-quadrature panel
 * (`FIXED_NODES` only).
 */
void
nc_data_cluster_wl_factor_set_rule_n (NcDataClusterWLFactor *dcwlf, guint rule_n)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  self->rule_n = rule_n;
}

/**
 * nc_data_cluster_wl_factor_get_rule_n:
 * @dcwlf: a #NcDataClusterWLFactor
 *
 * Returns: the Gauss-Legendre rule order per panel.
 */
guint
nc_data_cluster_wl_factor_get_rule_n (NcDataClusterWLFactor *dcwlf)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  return self->rule_n;
}

/**
 * nc_data_cluster_wl_factor_set_auto_nodes:
 * @dcwlf: a #NcDataClusterWLFactor
 * @auto_nodes: whether to auto-calibrate the per-galaxy fixed-node configuration
 *
 * Sets whether to automatically select, per galaxy, the minimal fixed
 * Gauss-Legendre configuration reaching #NcDataClusterWLFactor:node-reltol
 * (`FIXED_NODES` only).
 */
void
nc_data_cluster_wl_factor_set_auto_nodes (NcDataClusterWLFactor *dcwlf, gboolean auto_nodes)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  self->auto_nodes = auto_nodes;
}

/**
 * nc_data_cluster_wl_factor_get_auto_nodes:
 * @dcwlf: a #NcDataClusterWLFactor
 *
 * Returns: whether per-galaxy auto-nodes calibration is enabled.
 */
gboolean
nc_data_cluster_wl_factor_get_auto_nodes (NcDataClusterWLFactor *dcwlf)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  return self->auto_nodes;
}

/**
 * nc_data_cluster_wl_factor_set_node_reltol:
 * @dcwlf: a #NcDataClusterWLFactor
 * @node_reltol: target relative tolerance
 *
 * Sets the target relative tolerance for the per-galaxy fixed-node
 * selection (see #NcDataClusterWLFactor:auto-nodes).
 */
void
nc_data_cluster_wl_factor_set_node_reltol (NcDataClusterWLFactor *dcwlf, gdouble node_reltol)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  self->node_reltol = node_reltol;
}

/**
 * nc_data_cluster_wl_factor_get_node_reltol:
 * @dcwlf: a #NcDataClusterWLFactor
 *
 * Returns: the target relative tolerance for the per-galaxy fixed-node
 * selection.
 */
gdouble
nc_data_cluster_wl_factor_get_node_reltol (NcDataClusterWLFactor *dcwlf)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  return self->node_reltol;
}

/**
 * nc_data_cluster_wl_factor_set_max_total_nodes:
 * @dcwlf: a #NcDataClusterWLFactor
 * @max_total_nodes: safety ceiling on the total background node count
 *
 * Sets the safety ceiling on the total background node count explored by
 * the auto-nodes selection (see #NcDataClusterWLFactor:auto-nodes).
 */
void
nc_data_cluster_wl_factor_set_max_total_nodes (NcDataClusterWLFactor *dcwlf, guint max_total_nodes)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  self->max_total_nodes = max_total_nodes;
}

/**
 * nc_data_cluster_wl_factor_get_max_total_nodes:
 * @dcwlf: a #NcDataClusterWLFactor
 *
 * Returns: the safety ceiling on the total background node count for
 * auto-nodes selection.
 */
guint
nc_data_cluster_wl_factor_get_max_total_nodes (NcDataClusterWLFactor *dcwlf)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  return self->max_total_nodes;
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

