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
 * v1 scope, explicit rather than silent: only the `LNINT`-equivalent
 * adaptive log-domain 1D z-integral is implemented (legacy's `CUBATURE` and
 * `FIXED_NODES` schemes are deferred, not silently omitted -- all three are
 * mathematically exact, differing only in numerical strategy/cost). No
 * bootstrap, no Fisher-matrix support, no OpenMP parallelism, no fit-time
 * `r_min`/`r_max` weighting (matching #NcDataClusterWL's own deliberate
 * no-op there -- a previous attempt "introduced bias in the estimates").
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

  /* FIXED_NODES-only state: CUBATURE is not implemented by this class in
   * either phase (see the class documentation). */
  NcDataClusterWLIntegMethod integ_method;
  guint n_nodes;
  guint rule_n;

  /* Per-galaxy fixed-node grid: index i pairs with shape_data's element i.
   * fixed_bg_nodes[i] is NULL when galaxy i's redshift support lies entirely
   * in front of the lens (no background quadrature needed, see
   * _step_fixed_nodes_grid). */
  GPtrArray *fixed_bg_nodes;     /* element type: NcmIntegralFixed* */
  GPtrArray *z_nodes_per_galaxy; /* element type: NcmVector* */
  GArray *fixed_norm;            /* element type: gdouble */

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

  /* Linear-domain position integrand, used only by the fixed-nodes path's
   * control-variate combination (which is inherently linear-domain) --
   * built lazily alongside integ_pos/integ_z/integ_shape and re-prepared
   * unconditionally every cycle exactly like them. */
  NcGalaxyPositionFactorIntegrand *integ_pos_lin;

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
 * crosses the non-smooth point at z_cl. */
static void
_step_fixed_nodes_grid (NcDataClusterWLFactorPrivate *self, NcmMSet *mset, NcGalaxyShapeFactorData *s_data, guint gal_i)
{
  const guint n_total_nodes = (self->n_nodes - 1) * self->rule_n;
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

    bg_intf    = nc_galaxy_redshift_factor_make_fixed_nodes (self->redshift_factor, mset, s_data->z_data, bg_lo, z_hi, self->n_nodes, self->rule_n);
    z_nodes    = ncm_vector_new (n_total_nodes + 1);
    z_nodes_bg = ncm_vector_get_subvector (z_nodes, 1, n_total_nodes);

    ncm_integral_fixed_get_nodes (bg_intf, z_nodes_bg);
    ncm_vector_free (z_nodes_bg);
  }
  else
  {
    bg_intf = NULL;
    z_nodes = ncm_vector_new (1);
  }

  ncm_vector_set (z_nodes, 0, z_lo);

  g_ptr_array_index (self->fixed_bg_nodes, gal_i)     = bg_intf;
  g_ptr_array_index (self->z_nodes_per_galaxy, gal_i) = z_nodes;
  g_array_index (self->fixed_norm, gdouble, gal_i)    = norm;
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

  self->integ_method = NC_DATA_CLUSTER_WL_INTEG_METHOD_LNINT;
  self->n_nodes      = 10;
  self->rule_n       = 5;

  /* No destroy-func on fixed_bg_nodes: entries are nullable (NcmIntegralFixed
   * is not NULL-safe to free), so clearing/resizing this array is handled
   * manually with an explicit NULL check -- see _fixed_nodes_grid_reset().
   * z_nodes_per_galaxy's entries are never NULL, so its destroy-func is
   * safe to rely on for both resizing and disposal. */
  self->fixed_bg_nodes           = g_ptr_array_new ();
  self->z_nodes_per_galaxy       = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_vector_free);
  self->fixed_norm               = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->fixed_nodes_zcl          = GSL_NAN;
  self->fixed_nodes_n_nodes_seen = 0;
  self->fixed_nodes_rule_n_seen  = 0;

  self->integ_pos_lin = NULL;

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

  _fixed_bg_nodes_free_entries (self->fixed_bg_nodes);
  g_clear_pointer (&self->fixed_bg_nodes, g_ptr_array_unref);
  g_clear_pointer (&self->z_nodes_per_galaxy, g_ptr_array_unref);
  g_clear_pointer (&self->fixed_norm, g_array_unref);

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

    /* Linear-domain twin of integ_pos, needed only by the fixed-nodes
     * control-variate combination (inherently linear-domain -- see the
     * class documentation). */
    self->integ_pos_lin = nc_galaxy_position_factor_integ (self->position_factor, mset, FALSE);

    self->int_arg.integ_pos   = self->integ_pos;
    self->int_arg.integ_z     = self->integ_z;
    self->int_arg.integ_shape = self->integ_shape;

    self->int1d = ncm_integral1d_ptr_new (&_nc_data_cluster_wl_factor_lnint_integrand, NULL);
    ncm_integral1d_ptr_set_userdata (self->int1d, &self->int_arg);
  }
  else
  {
    nc_galaxy_position_factor_integrand_prepare (self->integ_pos, mset);
    nc_galaxy_redshift_factor_integrand_prepare (self->integ_z, mset);
    nc_galaxy_shape_factor_integrand_prepare (self->integ_shape, mset);
    nc_galaxy_position_factor_integrand_prepare (self->integ_pos_lin, mset);
  }

  ncm_integral1d_set_reltol (NCM_INTEGRAL1D (self->int1d), self->prec);
  ncm_integral1d_set_abstol (NCM_INTEGRAL1D (self->int1d), 0.0);

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
                       (self->n_nodes != self->fixed_nodes_n_nodes_seen) || (self->rule_n != self->fixed_nodes_rule_n_seen);

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
    self->fixed_nodes_zcl          = self->z_cl;
    self->fixed_nodes_n_nodes_seen = self->n_nodes;
    self->fixed_nodes_rule_n_seen  = self->rule_n;
  }

  /* Specialize and choose once: no per-galaxy branch on any of these. The
   * three fixed-nodes steps are only ever appended under FIXED_NODES;
   * ordering matters -- _step_shape_radius before _step_fixed_nodes_sigma
   * (sigma reads the projected radius) and _step_fixed_nodes_grid before
   * _step_fixed_nodes_crit (crit reads the current z_nodes). */
  if (pos_changed)
    steps[n_steps++] = &_step_pos;

  if (z_changed)
    steps[n_steps++] = &_step_z;

  if (radius_changed)
    steps[n_steps++] = &_step_shape_radius;

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

  if (optzs_changed)
    steps[n_steps++] = &_step_shape_optzs;

  if (pop_changed)
    steps[n_steps++] = &_step_shape_pop;

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
  const guint n_total_nodes                 = (self->n_nodes - 1) * self->rule_n;
  NcmVector *shape_at_nodes                 = ncm_vector_new (n_total_nodes + 1);

  /* Always the same [1, n_total_nodes) view into shape_at_nodes (itself
   * overwritten, not reallocated, every galaxy) -- built once per call
   * instead of once per galaxy to avoid a GObject alloc/dispose per galaxy
   * (ncm_vector_get_subvector is a full g_object_new). */
  NcmVector *sub = (n_total_nodes > 0) ? ncm_vector_get_subvector (shape_at_nodes, 1, n_total_nodes) : NULL;
  gdouble result = 0.0;
  guint gal_i;

  if (m2lnP_gal != NULL)
    g_assert_cmpuint (ncm_vector_len (m2lnP_gal), ==, self->len);

  for (gal_i = 0; gal_i < self->shape_data->len; gal_i++)
  {
    NcGalaxyShapeFactorData *s_data = g_ptr_array_index (self->shape_data, gal_i);
    NcmVector *z_nodes              = (NcmVector *) g_ptr_array_index (self->z_nodes_per_galaxy, gal_i);
    const gdouble int_pos           = nc_galaxy_position_factor_integrand_eval (self->integ_pos_lin, s_data->pos_data);
    gdouble P_gal, m2lnP_gal_i;

    nc_galaxy_shape_factor_eval_at_nodes (self->shape_factor, mset, s_data, z_nodes, shape_at_nodes);

    P_gal       = _nc_data_cluster_wl_factor_fixed_panels_integ (self, gal_i, shape_at_nodes, sub) * int_pos;
    m2lnP_gal_i = P_gal > 0.0 ? -2.0 * log (P_gal) : NC_GALAXY_LOW_PROB;

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

  if (sub != NULL)
    ncm_vector_free (sub);

  ncm_vector_free (shape_at_nodes);

  return result;
}

/* Dispatch: decided once per call, not per galaxy -- mirrors the steps[]
 * "specialize and choose once" convention above. CUBATURE is not
 * implemented by this class in either phase (see the class doc). */
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
    default:
      g_error ("_nc_data_cluster_wl_factor_eval_m2lnP: CUBATURE is not implemented by NcDataClusterWLFactor.");

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
   * Integration method for the redshift integral. `CUBATURE` is not
   * implemented by this class (see the class documentation).
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
 * Sets the redshift-integral method. `CUBATURE` is not implemented by this
 * class in either phase.
 */
void
nc_data_cluster_wl_factor_set_integ_method (NcDataClusterWLFactor *dcwlf, NcDataClusterWLIntegMethod integ_method)
{
  NcDataClusterWLFactorPrivate * const self = nc_data_cluster_wl_factor_get_instance_private (dcwlf);

  if (integ_method == NC_DATA_CLUSTER_WL_INTEG_METHOD_CUBATURE)
    g_error ("nc_data_cluster_wl_factor_set_integ_method: CUBATURE is not implemented by NcDataClusterWLFactor.");

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

