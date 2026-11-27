/***************************************************************************
 *            nc_galaxy_shape_factor.c
 *
 *  Thu Jul 2 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor.c
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
 * NcGalaxyShapeFactor:
 *
 * Shape likelihood-factor calculator for the weak-lensing pipeline.
 *
 * The calculator provides the HSM measurement engine: the intrinsic
 * ellipticity $\chi_I$ is drawn from the #NcGalaxyShapePop resolved from the
 * #NcmMSet, mapped deterministically by the reduced-shear transformation with
 * calibration bias $\tilde g = (1+m)\,g + c$, and observed with additive
 * per-galaxy Gaussian pixel noise,
 * $\epsilon_\mathrm{obs} = f_{\tilde g}(\chi_I) + n$, $n \sim N_2(0, \sigma_\mathrm{noise}^2)$.
 *
 * The likelihood factor is the marginal over the latent intrinsic shape,
 * $$P(\epsilon_\mathrm{obs} \mid g) = \int_{|\chi_I| < 1} d^2\chi_I\,
 *   P_\mathrm{pop}(\chi_I)\, N_2\!\big(\epsilon_\mathrm{obs} - f_{\tilde g}(\chi_I);
 *   \sigma_\mathrm{noise}^2\big).$$
 * The calculator evaluates the marginal of each galaxy at a given reduced shear
 * $g = g_1 + i g_2$. The caller passes $g$, with the calibration bias applied, in the
 * frame of the galaxy's ellipticity. Subclasses implement the eval_marginal() and
 * eval_ln_marginal() hooks.
 *
 * This is a calculator, not a model: it holds no fitted parameters and does
 * not live in the #NcmMSet. The only model it uses is the shape population,
 * taken from the mset passed to nc_galaxy_shape_factor_prepare(). The reduced
 * shear comes from the caller, for example #NcDataClusterWLFactor through
 * #NcGalaxyHaloWL and #NcWLClusterLens.
 *
 * Thread safety: every subclass must support its per-galaxy methods being
 * called for different galaxies from different threads at the same time,
 * because nc_data_cluster_wl_factor_data_prepare() prepares galaxies in
 * parallel. A per-galaxy method may write only that galaxy's own data
 * fragment; shared instance state may change only in prepare(), which always
 * runs alone. Models reached through the #NcmMSet that update themselves
 * lazily on their first call after a parameter change are handled by
 * #NcDataClusterWLFactor, which prepares one galaxy alone before starting
 * the rest.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_enum_types.h"
#include "nc/lss/galaxy/nc_galaxy_wl_obs.h"
#include "nc/lss/galaxy/nc_galaxy_shape_factor.h"
#include "ncm/core/ncm_prefetch_private.h"
#include "nc/lss/galaxy/nc_galaxy_shape_pop.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/stats/ncm_stats_vec.h"
#include "ncm/model/ncm_model.h"
#include "ncm/model/ncm_model_ctrl.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <complex.h>
#include <string.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcGalaxyShapeFactorPrivate
{
  NcGalaxyWLObsEllipConv ellip_conv;
  NcmStatsVec *obs_stats;

  NcmModelCtrl *ctrl_pop;
  guint64 pop_hash; /* the number of prepares that changed the population (#NcmModelCtrl) */
} NcGalaxyShapeFactorPrivate;

enum
{
  PROP_0,
  PROP_ELLIP_CONV,
  PROP_LEN,
};

/* *INDENT-OFF* */
G_DEFINE_QUARK (nc-galaxy-shape-factor-error, nc_galaxy_shape_factor_error)
/* *INDENT-ON* */
G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcGalaxyShapeFactor, nc_galaxy_shape_factor, G_TYPE_OBJECT)
G_DEFINE_BOXED_TYPE (NcGalaxyShapeFactorData, nc_galaxy_shape_factor_data, nc_galaxy_shape_factor_data_ref, nc_galaxy_shape_factor_data_unref); /* LCOV_EXCL_LINE */

static void
nc_galaxy_shape_factor_init (NcGalaxyShapeFactor *gsf)
{
  NcGalaxyShapeFactorPrivate * const self = nc_galaxy_shape_factor_get_instance_private (gsf);

  self->ellip_conv = 0;
  self->obs_stats  = ncm_stats_vec_new (6, NCM_STATS_VEC_COV, FALSE);
  self->ctrl_pop   = ncm_model_ctrl_new (NULL);
  self->pop_hash   = 0;
}

static void _nc_galaxy_shape_factor_set_ellip_conv (NcGalaxyShapeFactor *gsf, NcGalaxyWLObsEllipConv ellip_conv);

static void
_nc_galaxy_shape_factor_set_property (GObject *object, guint property_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactor *gsf = NC_GALAXY_SHAPE_FACTOR (object);

  g_return_if_fail (NC_IS_GALAXY_SHAPE_FACTOR (gsf));

  switch (property_id)
  {
    case PROP_ELLIP_CONV:
      _nc_galaxy_shape_factor_set_ellip_conv (gsf, g_value_get_enum (value));
      break;
    default:                                                          /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                          /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_factor_get_property (GObject *object, guint property_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactor *gsf = NC_GALAXY_SHAPE_FACTOR (object);

  g_return_if_fail (NC_IS_GALAXY_SHAPE_FACTOR (gsf));

  switch (property_id)
  {
    case PROP_ELLIP_CONV:
      g_value_set_enum (value, nc_galaxy_shape_factor_get_ellip_conv (gsf));
      break;
    default:                                                          /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                          /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_factor_dispose (GObject *object)
{
  NcGalaxyShapeFactor *gsf                = NC_GALAXY_SHAPE_FACTOR (object);
  NcGalaxyShapeFactorPrivate * const self = nc_galaxy_shape_factor_get_instance_private (gsf);

  ncm_stats_vec_clear (&self->obs_stats);
  ncm_model_ctrl_clear (&self->ctrl_pop);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_parent_class)->dispose (object);
}

static void
_nc_galaxy_shape_factor_finalize (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_parent_class)->finalize (object);
}

/* LCOV_EXCL_START */
static void
_nc_galaxy_shape_factor_data_init (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data)
{
  g_error ("_nc_galaxy_shape_factor_data_init: method not implemented.");
}

static gdouble
_nc_galaxy_shape_factor_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  g_error ("_nc_galaxy_shape_factor_eval_marginal: method not implemented.");

  return 0.0;
}

static gdouble
_nc_galaxy_shape_factor_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  g_error ("_nc_galaxy_shape_factor_eval_ln_marginal: method not implemented.");

  return 0.0;
}

/* LCOV_EXCL_STOP */

static void
_nc_galaxy_shape_factor_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset)
{
  /* Default: the evaluation strategy needs no additional setup. */
}

static void
_nc_galaxy_shape_factor_data_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data, const gdouble z_max)
{
  /* Default: the evaluation strategy caches nothing per galaxy. */
}

/* Default eval_nodes: eval_marginal or eval_ln_marginal at each shear */
static void
_nc_galaxy_shape_factor_eval_nodes (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const guint n,
                                    const gdouble *g_1, const gdouble *g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2,
                                    const gboolean want_log, gdouble *out)
{
  NcGalaxyShapeFactorClass *klass = NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsf);
  guint j;

  if (want_log)
    for (j = 0; j < n; j++)
      out[j] = klass->eval_ln_marginal (gsf, pop, data, g_1[j], g_2[j], epsilon_obs_1, epsilon_obs_2);

  else
    for (j = 0; j < n; j++)
      out[j] = klass->eval_marginal (gsf, pop, data, g_1[j], g_2[j], epsilon_obs_1, epsilon_obs_2);
}

/* Default eval_line: no evaluation along the line */
static gboolean
_nc_galaxy_shape_factor_eval_line (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const guint n,
                                   const gdouble *t, const gdouble u_1, const gdouble u_2, const gdouble b_1, const gdouble b_2,
                                   const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gboolean want_log, gdouble *out)
{
  return FALSE;
}

static gchar *
_nc_galaxy_shape_factor_get_desc (NcGalaxyShapeFactor *gsf)
{
  const NcGalaxyWLObsEllipConv ellip_conv = nc_galaxy_shape_factor_get_ellip_conv (gsf);
  const GEnumValue *ellip_val             = ncm_cfg_enum_get_value (NC_TYPE_GALAXY_WL_OBS_ELLIP_CONV, ellip_conv);

  /* Default: only the property shared by every subclass. */
  return g_strdup_printf ("%s (ellip_conv=%s)", G_OBJECT_TYPE_NAME (gsf), ellip_val->value_nick);
}

static void
nc_galaxy_shape_factor_class_init (NcGalaxyShapeFactorClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_shape_factor_set_property;
  object_class->get_property = &_nc_galaxy_shape_factor_get_property;
  object_class->dispose      = &_nc_galaxy_shape_factor_dispose;
  object_class->finalize     = &_nc_galaxy_shape_factor_finalize;

  /**
   * NcGalaxyShapeFactor:ellip-conv:
   *
   * Weak lensing observables ellipticity convention #NcGalaxyWLObsEllipConv.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_ELLIP_CONV,
                                   g_param_spec_enum ("ellip-conv",
                                                      "Ellipticity convention",
                                                      "Weak lensing observables ellipticity convention",
                                                      NC_TYPE_GALAXY_WL_OBS_ELLIP_CONV,
                                                      NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  klass->data_init        = &_nc_galaxy_shape_factor_data_init;
  klass->prepare          = &_nc_galaxy_shape_factor_prepare;
  klass->data_prepare     = &_nc_galaxy_shape_factor_data_prepare;
  klass->eval_marginal    = &_nc_galaxy_shape_factor_eval_marginal;
  klass->eval_ln_marginal = &_nc_galaxy_shape_factor_eval_ln_marginal;
  klass->get_desc         = &_nc_galaxy_shape_factor_get_desc;
  klass->eval_nodes       = &_nc_galaxy_shape_factor_eval_nodes;
  klass->eval_line        = &_nc_galaxy_shape_factor_eval_line;
}

static void
_nc_galaxy_shape_factor_set_ellip_conv (NcGalaxyShapeFactor *gsf, NcGalaxyWLObsEllipConv ellip_conv)
{
  NcGalaxyShapeFactorPrivate * const self = nc_galaxy_shape_factor_get_instance_private (gsf);

  switch (ellip_conv)
  {
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
      break;
    default:                                                                                           /* LCOV_EXCL_LINE */
      g_error ("nc_galaxy_shape_factor_set_ellip_conv: ellipse type %d not implemented.", ellip_conv); /* LCOV_EXCL_LINE */
      break;                                                                                           /* LCOV_EXCL_LINE */
  }

  self->ellip_conv = ellip_conv;
}

/**
 * nc_galaxy_shape_factor_get_ellip_conv:
 * @gsf: a #NcGalaxyShapeFactor
 *
 * Gets the ellipticity convention of @gsf.
 *
 * Returns: a #NcGalaxyWLObsEllipConv
 */
NcGalaxyWLObsEllipConv
nc_galaxy_shape_factor_get_ellip_conv (NcGalaxyShapeFactor *gsf)
{
  NcGalaxyShapeFactorPrivate * const self = nc_galaxy_shape_factor_get_instance_private (gsf);

  return self->ellip_conv;
}

/**
 * nc_galaxy_shape_factor_check_obs:
 * @gsf: a #NcGalaxyShapeFactor
 * @obs: a #NcGalaxyWLObs
 * @error: return location for a #GError
 *
 * Checks that @obs's ellipticity convention matches @gsf's own
 * #NcGalaxyWLObsEllipConv. Validation the shape factor owns (it is the one
 * that knows what convention it expects), not its callers.
 *
 * Returns: %TRUE if the conventions match, %FALSE otherwise (with @error set).
 */
gboolean
nc_galaxy_shape_factor_check_obs (NcGalaxyShapeFactor *gsf, NcGalaxyWLObs *obs, GError **error)
{
  const NcGalaxyWLObsEllipConv obs_conv = nc_galaxy_wl_obs_get_ellip_conv (obs);
  const NcGalaxyWLObsEllipConv gsf_conv = nc_galaxy_shape_factor_get_ellip_conv (gsf);

  if (obs_conv != gsf_conv)
  {
    ncm_util_set_or_call_error (error, NC_GALAXY_SHAPE_FACTOR_ERROR, NC_GALAXY_SHAPE_FACTOR_ERROR_ELLIP_CONV_MISMATCH,
                                "nc_galaxy_shape_factor_check_obs: observation ellipticity convention (%d) "
                                "does not match the shape factor's (%d).", obs_conv, gsf_conv);

    return FALSE;
  }

  return TRUE;
}

/**
 * nc_galaxy_shape_factor_get_desc:
 * @gsf: a #NcGalaxyShapeFactor
 *
 * Returns: (transfer full): a human-readable one-line description of this
 * scheme's own configuration (default: the concrete type name plus
 * :ellip-conv).
 */
gchar *
nc_galaxy_shape_factor_get_desc (NcGalaxyShapeFactor *gsf)
{
  return NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsf)->get_desc (gsf);
}

/**
 * nc_galaxy_shape_factor_prepare:
 * @gsf: a #NcGalaxyShapeFactor
 * @mset: a #NcmMSet with a #NcGalaxyShapePop
 *
 * Increments the population counter (nc_galaxy_shape_factor_get_pop_hash()) when the
 * #NcGalaxyShapePop of @mset is another object or its state changed (#NcmModelCtrl), then
 * calls the subclass's `prepare` hook. No model is kept: the evaluations and the
 * per-galaxy updates take the population as an argument. It does nothing else when
 * nothing changed, so it can be called every prepare cycle.
 */
void
nc_galaxy_shape_factor_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset)
{
  NcGalaxyShapeFactorClass *klass         = NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsf);
  NcGalaxyShapeFactorPrivate * const self = nc_galaxy_shape_factor_get_instance_private (gsf);
  NcGalaxyShapePop *pop                   = NC_GALAXY_SHAPE_POP (ncm_mset_peek (mset, nc_galaxy_shape_pop_id ()));

  g_assert_nonnull (pop);

  if (ncm_model_ctrl_update (self->ctrl_pop, NCM_MODEL (pop)))
    self->pop_hash++;

  klass->prepare (gsf, mset);
}

/**
 * nc_galaxy_shape_factor_data_prepare:
 * @gsf: a #NcGalaxyShapeFactor
 * @mset: a #NcmMSet
 * @data: a #NcGalaxyShapeFactorData
 * @z_max: an upper bound on the galaxy's redshift support
 *
 * Per-galaxy prepare: calls the subclass's data_prepare with @mset and
 * @data, so that it precomputes what depends on the galaxy and on the
 * models together. Call it after nc_galaxy_shape_factor_update_data_pop() and
 * after setting the bound fields of @data, whose values it may read, and
 * before any redshift quadrature evaluates the marginal of the galaxy.
 *
 * It is an optimization, not a precondition: a subclass that caches here
 * must build the same state lazily when this function did not run for a
 * galaxy.
 */
void
nc_galaxy_shape_factor_data_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data, const gdouble z_max)
{
  NcGalaxyShapeFactorClass *klass = NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsf);

  klass->data_prepare (gsf, mset, data, z_max);
}

/**
 * nc_galaxy_shape_factor_data_prefetch:
 * @gsf: a #NcGalaxyShapeFactor
 * @data: a #NcGalaxyShapeFactorData
 * @stage: prefetch stage, 0, 1 or 2
 *
 * Issues software prefetches for the memory nc_galaxy_shape_factor_eval_shear_nodes()
 * reads for @data, so that a loop over galaxies fetches the next galaxies'
 * state while it evaluates the current one. Each galaxy's state is a chain
 * of small heap allocations, which the hardware prefetcher cannot follow;
 * with many galaxies it does not fit in cache, and without prefetching
 * every evaluation waits on memory at each link of the chain.
 *
 * The chain is walked one link per @stage, so that no stage reads a pointer
 * that is not already being fetched:
 *
 * - stage 0: @data itself;
 * - stage 1: the structures @data points to, read through @data;
 * - stage 2: the arrays those structures point to.
 *
 * The caller issues the stages for a galaxy on consecutive iterations, the
 * highest stage for the galaxy evaluated next: stage 2 for galaxy i+1,
 * stage 1 for i+2 and stage 0 for i+3 while evaluating galaxy i. A stage
 * called out of order changes nothing and stalls on the pointer the skipped
 * stage would have fetched.
 *
 * This only issues hints: it changes no state and has no effect on results.
 */
void
nc_galaxy_shape_factor_data_prefetch (NcGalaxyShapeFactor *gsf, NcGalaxyShapeFactorData *data, const guint stage)
{
  NcGalaxyShapeFactorClass *klass = NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsf);

  switch (stage)
  {
    case 0:
      ncm_prefetch_span (data, sizeof (NcGalaxyShapeFactorData));
      break;
    case 1:
    case 2:
      break;
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
      break;                   /* LCOV_EXCL_LINE */
  }

  if (klass->data_prefetch != NULL)
    klass->data_prefetch (gsf, data, stage);
}

/**
 * nc_galaxy_shape_factor_get_pop_hash:
 * @gsf: a #NcGalaxyShapeFactor
 *
 * Returns: the number of nc_galaxy_shape_factor_prepare() calls in which the population
 *   changed; a different value means a change.
 */
guint64
nc_galaxy_shape_factor_get_pop_hash (NcGalaxyShapeFactor *gsf)
{
  NcGalaxyShapeFactorPrivate * const self = nc_galaxy_shape_factor_get_instance_private (gsf);

  return self->pop_hash;
}

/**
 * nc_galaxy_shape_factor_update_data_pop:
 * @gsf: a #NcGalaxyShapeFactor
 * @pop: the #NcGalaxyShapePop of the last nc_galaxy_shape_factor_prepare()
 * @data: a #NcGalaxyShapeFactorData
 *
 * Resolves the current parameters of @pop into @data's population fragment.
 * Unconditional: call it when nc_galaxy_shape_factor_get_pop_hash() changed.
 */
void
nc_galaxy_shape_factor_update_data_pop (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data)
{
  nc_galaxy_shape_pop_prepare (pop, data->pop_data);
}

/**
 * nc_galaxy_shape_factor_apply_shear:
 * @gsf: a #NcGalaxyShapeFactor
 * @g: input reduced shear as a #NcmComplex
 * @E: input intrinsic ellipticity as a #NcmComplex
 * @E_obs: output observed ellipticity as a #NcmComplex
 *
 * Applies the reduced shear @g to the intrinsic ellipticity @E, storing the
 * resulting observed ellipticity in @E_obs. The transformation depends on the
 * #NcGalaxyWLObsEllipConv configured in @gsf.
 */
void
nc_galaxy_shape_factor_apply_shear (NcGalaxyShapeFactor *gsf, const NcmComplex *g, const NcmComplex *E, NcmComplex *E_obs)
{
  NcGalaxyShapeFactorPrivate * const self = nc_galaxy_shape_factor_get_instance_private (gsf);
  const complex double gn                 = ncm_complex_c (g);
  const complex double En                 = ncm_complex_c (E);

  switch (self->ellip_conv)
  {
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
      ncm_complex_set_c (E_obs, nc_wl_ellipticity_apply_shear_trace (gn, En));
      break;
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
      ncm_complex_set_c (E_obs, nc_wl_ellipticity_apply_shear_trace_det (gn, En));
      break;
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
  }
}

/**
 * nc_galaxy_shape_factor_apply_shear_inv:
 * @gsf: a #NcGalaxyShapeFactor
 * @g: input reduced shear as a #NcmComplex
 * @E_obs: input observed ellipticity as a #NcmComplex
 * @E: output intrinsic ellipticity as a #NcmComplex
 *
 * Applies the inverse shear transformation using @g to recover the intrinsic
 * ellipticity @E from the observed ellipticity @E_obs. The transformation
 * depends on the #NcGalaxyWLObsEllipConv configured in @gsf.
 */
void
nc_galaxy_shape_factor_apply_shear_inv (NcGalaxyShapeFactor *gsf, const NcmComplex *g, const NcmComplex *E_obs, NcmComplex *E)
{
  NcGalaxyShapeFactorPrivate * const self = nc_galaxy_shape_factor_get_instance_private (gsf);
  const complex double gn                 = ncm_complex_c (g);
  const complex double En_obs             = ncm_complex_c (E_obs);

  switch (self->ellip_conv)
  {
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
      ncm_complex_set_c (E, nc_wl_ellipticity_apply_shear_inv_trace (gn, En_obs));
      break;
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
      ncm_complex_set_c (E, nc_wl_ellipticity_apply_shear_inv_trace_det (gn, En_obs));
      break;
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
  }
}

/**
 * nc_galaxy_shape_factor_lndet_jac:
 * @gsf: a #NcGalaxyShapeFactor
 * @g: input reduced shear as a #NcmComplex
 * @E_obs: input observed ellipticity as a #NcmComplex
 *
 * Computes the natural logarithm of the absolute value of the Jacobian
 * determinant of the transformation from intrinsic to observed ellipticity.
 *
 * Returns: the log-determinant of the shear Jacobian.
 */
gdouble
nc_galaxy_shape_factor_lndet_jac (NcGalaxyShapeFactor *gsf, const NcmComplex *g, const NcmComplex *E_obs)
{
  NcGalaxyShapeFactorPrivate * const self = nc_galaxy_shape_factor_get_instance_private (gsf);
  const complex double gn                 = ncm_complex_c (g);
  const complex double En_obs             = ncm_complex_c (E_obs);

  switch (self->ellip_conv)
  {
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
      return nc_wl_ellipticity_lndet_jac_trace (gn, En_obs);

    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
      return nc_wl_ellipticity_lndet_jac_trace_det (gn, En_obs);

    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */

      return 0.0; /* LCOV_EXCL_LINE */
  }
}

/* LCOV_EXCL_START */

/**
 * nc_galaxy_shape_factor_data_ref:
 * @data: a #NcGalaxyShapeFactorData
 *
 * Increases the reference count of @data by one.
 *
 * Returns: (transfer full): @data.
 */
NcGalaxyShapeFactorData *
nc_galaxy_shape_factor_data_ref (NcGalaxyShapeFactorData *data)
{
  g_atomic_ref_count_inc (&data->ref_count);

  return data;
}

/* LCOV_EXCL_STOP */

/**
 * nc_galaxy_shape_factor_data_unref:
 * @data: a #NcGalaxyShapeFactorData
 *
 * Decreases the reference count of @data by one. If the reference count
 * reaches 0, the data is freed.
 *
 */
void
nc_galaxy_shape_factor_data_unref (NcGalaxyShapeFactorData *data)
{
  if (g_atomic_ref_count_dec (&data->ref_count))
  {
    g_assert_nonnull (data->ldata_destroy);
    data->ldata_destroy (data->ldata);
    nc_galaxy_shape_pop_data_unref (data->pop_data);
    g_free (data);
  }
}

/**
 * nc_galaxy_shape_factor_data_read_row:
 * @data: a #NcGalaxyShapeFactorData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Reads the row @i from @obs into @data, cascading to the population fragment and the
 * subclass fragment.
 */
void
nc_galaxy_shape_factor_data_read_row (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  data->coord         = nc_galaxy_wl_obs_get_coord (obs);
  data->epsilon_int_1 = ncm_catalog_get (NCM_CATALOG (obs), NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_INT_1, i, NULL);
  data->epsilon_int_2 = ncm_catalog_get (NCM_CATALOG (obs), NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_INT_2, i, NULL);
  data->epsilon_obs_1 = ncm_catalog_get (NCM_CATALOG (obs), NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_OBS_1, i, NULL);
  data->epsilon_obs_2 = ncm_catalog_get (NCM_CATALOG (obs), NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_OBS_2, i, NULL);
  data->std_noise     = ncm_catalog_get (NCM_CATALOG (obs), NC_GALAXY_SHAPE_FACTOR_COL_STD_NOISE, i, NULL);
  data->c1            = ncm_catalog_get (NCM_CATALOG (obs), NC_GALAXY_SHAPE_FACTOR_COL_C1, i, NULL);
  data->c2            = ncm_catalog_get (NCM_CATALOG (obs), NC_GALAXY_SHAPE_FACTOR_COL_C2, i, NULL);
  data->m             = ncm_catalog_get (NCM_CATALOG (obs), NC_GALAXY_SHAPE_FACTOR_COL_M, i, NULL);

  nc_galaxy_shape_pop_data_read_row (data->pop_data, obs, i);
  data->ldata_read_row (data, obs, i);
}

/**
 * nc_galaxy_shape_factor_data_write_row:
 * @data: a #NcGalaxyShapeFactorData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Writes the row @i of @data into @obs, cascading to the population fragment and the
 * subclass fragment.
 */
void
nc_galaxy_shape_factor_data_write_row (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  ncm_catalog_set (NCM_CATALOG (obs), NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_INT_1, i, data->epsilon_int_1, NULL);
  ncm_catalog_set (NCM_CATALOG (obs), NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_INT_2, i, data->epsilon_int_2, NULL);
  ncm_catalog_set (NCM_CATALOG (obs), NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_OBS_1, i, data->epsilon_obs_1, NULL);
  ncm_catalog_set (NCM_CATALOG (obs), NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_OBS_2, i, data->epsilon_obs_2, NULL);
  ncm_catalog_set (NCM_CATALOG (obs), NC_GALAXY_SHAPE_FACTOR_COL_STD_NOISE, i, data->std_noise, NULL);
  ncm_catalog_set (NCM_CATALOG (obs), NC_GALAXY_SHAPE_FACTOR_COL_C1, i, data->c1, NULL);
  ncm_catalog_set (NCM_CATALOG (obs), NC_GALAXY_SHAPE_FACTOR_COL_C2, i, data->c2, NULL);
  ncm_catalog_set (NCM_CATALOG (obs), NC_GALAXY_SHAPE_FACTOR_COL_M, i, data->m, NULL);

  nc_galaxy_shape_pop_data_write_row (data->pop_data, obs, i);
  data->ldata_write_row (data, obs, i);
}

/**
 * nc_galaxy_shape_factor_data_required_columns:
 * @data: a #NcGalaxyShapeFactorData
 *
 * Returns: (element-type utf8) (transfer full): the columns of the shape data: its own,
 *   those of the population fragment and those of the subclass fragment.
 */
GList *
nc_galaxy_shape_factor_data_required_columns (NcGalaxyShapeFactorData *data)
{
  GList *columns = NULL;

  columns = g_list_append (columns, g_strdup (NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_INT_1));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_INT_2));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_OBS_1));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_OBS_2));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SHAPE_FACTOR_COL_STD_NOISE));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SHAPE_FACTOR_COL_C1));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SHAPE_FACTOR_COL_C2));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SHAPE_FACTOR_COL_M));

  columns = g_list_concat (columns, nc_galaxy_shape_pop_data_required_columns (data->pop_data));
  data->ldata_required_columns (data, &columns);

  return columns;
}

/**
 * nc_galaxy_shape_factor_ref:
 * @gsf: a #NcGalaxyShapeFactor
 *
 * Increases the reference count of @gsf by one.
 *
 * Returns: (transfer full): @gsf.
 */
NcGalaxyShapeFactor *
nc_galaxy_shape_factor_ref (NcGalaxyShapeFactor *gsf)
{
  return g_object_ref (gsf);
}

/**
 * nc_galaxy_shape_factor_free:
 * @gsf: a #NcGalaxyShapeFactor
 *
 * Decreases the reference count of @gsf by one.
 *
 */
void
nc_galaxy_shape_factor_free (NcGalaxyShapeFactor *gsf)
{
  g_object_unref (gsf);
}

/**
 * nc_galaxy_shape_factor_clear:
 * @gsf: a #NcGalaxyShapeFactor
 *
 * Decreases the reference count of @gsf by one, and sets the pointer *@gsf to
 * NULL.
 *
 */
void
nc_galaxy_shape_factor_clear (NcGalaxyShapeFactor **gsf)
{
  g_clear_object (gsf);
}

/**
 * nc_galaxy_shape_factor_data_new:
 * @gsf: a #NcGalaxyShapeFactor
 * @mset: a #NcmMSet
 *
 * Creates a new per-galaxy shape data. The intrinsic-ellipticity fragment is
 * allocated by the #NcGalaxyShapePop resolved from @mset.
 *
 * Returns: (transfer full): a new #NcGalaxyShapeFactorData object.
 */
NcGalaxyShapeFactorData *
nc_galaxy_shape_factor_data_new (NcGalaxyShapeFactor *gsf, NcmMSet *mset)
{
  NcGalaxyShapeFactorData *data = g_new0 (NcGalaxyShapeFactorData, 1);
  NcGalaxyShapePop *pop         = NC_GALAXY_SHAPE_POP (ncm_mset_peek (mset, nc_galaxy_shape_pop_id ()));

  g_assert_nonnull (pop);

  data->pop_data = nc_galaxy_shape_pop_data_new (pop);
  data->coord    = NC_WL_ELLIPTICITY_FRAME_CELESTIAL;

  data->bound_ready = FALSE;
  data->bound_holds = FALSE;
  data->bound_z_cl  = GSL_NAN;
  data->bound_R_min = GSL_NAN;
  data->bound_kappa = GSL_NAN;
  data->bound_g_t   = GSL_NAN;

  data->ldata                  = NULL;
  data->ldata_destroy          = NULL;
  data->ldata_read_row         = NULL;
  data->ldata_write_row        = NULL;
  data->ldata_required_columns = NULL;

  g_atomic_ref_count_init (&data->ref_count);
  NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsf)->data_init (gsf, mset, data);

  g_assert_nonnull (data->ldata_destroy);
  g_assert_nonnull (data->ldata_read_row);
  g_assert_nonnull (data->ldata_write_row);
  g_assert_nonnull (data->ldata_required_columns);

  return data;
}

/**
 * nc_galaxy_shape_factor_data_set:
 * @gsf: a #NcGalaxyShapeFactor
 * @data: a #NcGalaxyShapeFactorData
 * @epsilon_obs_1: the observed ellipticity component 1
 * @epsilon_obs_2: the observed ellipticity component 2
 * @std_noise: the observational shape dispersion
 * @c1: the first additive bias parameter (in the @coord frame)
 * @c2: the second additive bias parameter (in the @coord frame)
 * @m: the multiplicative bias parameter
 * @coord: the ellipticity handedness frame #NcWLEllipticityFrame
 *
 * Sets the observed ellipticity and the per-galaxy measurement quantities.
 *
 */
void
nc_galaxy_shape_factor_data_set (NcGalaxyShapeFactor *gsf, NcGalaxyShapeFactorData *data, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gdouble std_noise, const gdouble c1, const gdouble c2, const gdouble m, NcWLEllipticityFrame coord)
{
  data->epsilon_obs_1 = epsilon_obs_1;
  data->epsilon_obs_2 = epsilon_obs_2;
  data->std_noise     = std_noise;
  data->c1            = c1;
  data->c2            = c2;
  data->m             = m;
  data->coord         = coord;
}

/**
 * nc_galaxy_shape_factor_data_get:
 * @gsf: a #NcGalaxyShapeFactor
 * @data: a #NcGalaxyShapeFactorData
 * @epsilon_obs_1: (out): the observed ellipticity component 1
 * @epsilon_obs_2: (out): the observed ellipticity component 2
 * @std_noise: (out): the observational shape dispersion
 * @c1: (out): the first additive bias parameter
 * @c2: (out): the second additive bias parameter
 * @m: (out): the multiplicative bias parameter
 *
 * Gets the observed ellipticity and the per-galaxy measurement quantities.
 *
 */
void
nc_galaxy_shape_factor_data_get (NcGalaxyShapeFactor *gsf, NcGalaxyShapeFactorData *data, gdouble *epsilon_obs_1, gdouble *epsilon_obs_2, gdouble *std_noise, gdouble *c1, gdouble *c2, gdouble *m)
{
  *epsilon_obs_1 = data->epsilon_obs_1;
  *epsilon_obs_2 = data->epsilon_obs_2;
  *std_noise     = data->std_noise;
  *c1            = data->c1;
  *c2            = data->c2;
  *m             = data->m;
}

/**
 * nc_galaxy_shape_factor_prepare_data_array:
 * @gsf: a #NcGalaxyShapeFactor
 * @mset: a #NcmMSet
 * @data_array: (element-type NcGalaxyShapeFactorData): a #GPtrArray of #NcGalaxyShapeFactorData
 *
 * Prepares @gsf (nc_galaxy_shape_factor_prepare()) and the population fragment
 * of every galaxy of @data_array (nc_galaxy_shape_factor_update_data_pop()).
 *
 * Returns: TRUE.
 */
gboolean
nc_galaxy_shape_factor_prepare_data_array (NcGalaxyShapeFactor *gsf, NcmMSet *mset, GPtrArray *data_array)
{
  NcGalaxyShapePop *pop = NC_GALAXY_SHAPE_POP (ncm_mset_peek (mset, nc_galaxy_shape_pop_id ()));
  guint i;

  nc_galaxy_shape_factor_prepare (gsf, mset);

  for (i = 0; i < data_array->len; i++)
    nc_galaxy_shape_factor_update_data_pop (gsf, pop, g_ptr_array_index (data_array, i));

  return TRUE;
}

/**
 * nc_galaxy_shape_factor_eval_shear:
 * @gsf: a #NcGalaxyShapeFactor
 * @pop: the #NcGalaxyShapePop of the last nc_galaxy_shape_factor_prepare()
 * @data: a #NcGalaxyShapeFactorData
 * @g_1: real part of the reduced shear, bias applied, in the frame of @data
 * @g_2: imaginary part of the reduced shear
 * @want_log: whether to return the logarithm
 *
 * Evaluates the marginal $P(\epsilon_\mathrm{obs} | g)$ of the galaxy of @data, or its
 * logarithm, with @pop.
 *
 * Returns: the marginal or its logarithm.
 */
gdouble
nc_galaxy_shape_factor_eval_shear (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gboolean want_log)
{
  NcGalaxyShapeFactorClass *klass = NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsf);

  return want_log ?
         klass->eval_ln_marginal (gsf, pop, data, g_1, g_2, data->epsilon_obs_1, data->epsilon_obs_2) :
         klass->eval_marginal (gsf, pop, data, g_1, g_2, data->epsilon_obs_1, data->epsilon_obs_2);
}

/**
 * nc_galaxy_shape_factor_eval_shear_nodes:
 * @gsf: a #NcGalaxyShapeFactor
 * @pop: the #NcGalaxyShapePop of the last nc_galaxy_shape_factor_prepare()
 * @data: a #NcGalaxyShapeFactorData
 * @n: the number of shears
 * @g_1: (array length=n): real parts of the reduced shears, bias applied
 * @g_2: (array length=n): imaginary parts of the reduced shears
 * @want_log: whether to compute the logarithms
 * @out: (array length=n) (out caller-allocates): the values
 *
 * Evaluates nc_galaxy_shape_factor_eval_shear() of the galaxy of @data at the @n shears,
 * together when the subclass supports it.
 */
void
nc_galaxy_shape_factor_eval_shear_nodes (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const guint n, const gdouble *g_1, const gdouble *g_2, const gboolean want_log, gdouble *out)
{
  NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsf)->eval_nodes (gsf, pop, data, n, g_1, g_2, data->epsilon_obs_1, data->epsilon_obs_2, want_log, out);
}

/**
 * nc_galaxy_shape_factor_eval_shear_line:
 * @gsf: a #NcGalaxyShapeFactor
 * @pop: the #NcGalaxyShapePop of the last nc_galaxy_shape_factor_prepare()
 * @data: a #NcGalaxyShapeFactorData
 * @n: the number of shears
 * @t: (array length=n): the positions $t_j$ along the line
 * @u_1: real part of the direction $u$
 * @u_2: imaginary part of the direction
 * @b_1: real part of the offset $b$
 * @b_2: imaginary part of the offset
 * @want_log: whether to compute the logarithms
 * @out: (array length=n): the values
 *
 * Evaluates nc_galaxy_shape_factor_eval_shear() of the galaxy of @data at the @n shears
 * $b + t_j u$ of one line, when the subclass computes them along the line. @out may be
 * @t.
 *
 * Returns: whether the subclass computed them; when FALSE, @out is unchanged unless it is
 *   @t.
 */
gboolean
nc_galaxy_shape_factor_eval_shear_line (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const guint n, const gdouble *t, const gdouble u_1, const gdouble u_2,
                                        const gdouble b_1, const gdouble b_2, const gboolean want_log, gdouble *out)
{
  return NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsf)->eval_line (gsf, pop, data, n, t, u_1, u_2, b_1, b_2, data->epsilon_obs_1, data->epsilon_obs_2, want_log, out);
}

/**
 * nc_galaxy_shape_factor_gen_shear:
 * @gsf: a #NcGalaxyShapeFactor
 * @pop: the #NcGalaxyShapePop of the last nc_galaxy_shape_factor_prepare()
 * @data: a #NcGalaxyShapeFactorData
 * @lensed: whether the galaxy is behind the lens
 * @g_1: real part of the reduced shear, bias applied, in the frame of @data
 * @g_2: imaginary part of the reduced shear
 * @rng: a #NcmRNG
 *
 * Draws the intrinsic and observed ellipticities of the galaxy of @data with @pop: the
 * noise, then the intrinsic ellipticity, both expressed in the frame of @data; the
 * observed one is the intrinsic one sheared by $g$ when @lensed and translated by $g$
 * (the bias alone) otherwise, plus the noise.
 */
void
nc_galaxy_shape_factor_gen_shear (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gboolean lensed, const gdouble g_1, const gdouble g_2, NcmRNG *rng)
{
  const gdouble noise1   = ncm_rng_gaussian_gen (rng, 0.0, data->std_noise);
  const gdouble noise2   = ncm_rng_gaussian_gen (rng, 0.0, data->std_noise);
  complex double noise   = noise1 + I * noise2;
  const complex double g = g_1 + I * g_2;
  gdouble e_int_1, e_int_2;
  complex double e_s, e_o;

  nc_galaxy_shape_pop_prepare (pop, data->pop_data);
  nc_galaxy_shape_pop_gen (pop, data->pop_data, rng, &e_int_1, &e_int_2);
  e_s = e_int_1 + I * e_int_2;

  /*
   * e_s and the noise are isotropic: expressing them in data->coord changes nothing
   * statistically and keeps them in the frame of g
   */
  e_s   = nc_wl_ellipticity_celestial_to_frame (data->coord, e_s);
  noise = nc_wl_ellipticity_celestial_to_frame (data->coord, noise);

  if (lensed)
  {
    NcmComplex cplx_E_s = NCM_COMPLEX_INIT (e_s);
    NcmComplex cplx_E_o, cplx_g;

    ncm_complex_set_c (&cplx_g, g);
    nc_galaxy_shape_factor_apply_shear (gsf, &cplx_g, &cplx_E_s, &cplx_E_o);

    e_o = ncm_complex_c (&cplx_E_o);
  }
  else
  {
    e_o = e_s + g;
  }

  e_o += noise;

  data->epsilon_int_1 = creal (e_s);
  data->epsilon_int_2 = cimag (e_s);
  data->epsilon_obs_1 = creal (e_o);
  data->epsilon_obs_2 = cimag (e_o);
}

/**
 * nc_galaxy_shape_factor_direct_estimate:
 * @gsf: a #NcGalaxyShapeFactor
 * @mset: a #NcmMSet
 * @data_array: (element-type NcGalaxyShapeFactorData): a #GPtrArray of #NcGalaxyShapeFactorData
 * @rotations: (element-type gdouble): the pairs $\cos 2\phi_i$, $\sin 2\phi_i$ of every galaxy
 * @gt: (out): the reduced tangential shear
 * @gx: (out): the cross shear
 * @sigma_t: (out): the tangential scatter
 * @sigma_x: (out): the cross scatter
 * @rho: (out): the correlation
 *
 * Computes the inverse-variance-weighted direct estimate of the reduced shear
 * from the observed ellipticities. The per-galaxy intrinsic scatter enters
 * through the population model's rms (see nc_galaxy_shape_pop_e_rms()). The
 * ellipticity and the additive bias c1, c2 of each galaxy are rotated to its
 * tangential/cross frame before they are averaged.
 *
 * The angle $\phi_i$ is the celestial position angle of galaxy $i$ around the
 * halo center (nc_galaxy_halo_wl_data_get_rotation()). @rotations holds
 * $2 n$ values for the $n$ galaxies of @data_array.
 */
void
nc_galaxy_shape_factor_direct_estimate (NcGalaxyShapeFactor *gsf, NcmMSet *mset, GPtrArray *data_array, GArray *rotations, gdouble *gt, gdouble *gx, gdouble *sigma_t, gdouble *sigma_x, gdouble *rho)
{
  NcGalaxyShapeFactorPrivate * const self = nc_galaxy_shape_factor_get_instance_private (gsf);
  NcGalaxyShapePop *pop                   = NC_GALAXY_SHAPE_POP (ncm_mset_peek (mset, nc_galaxy_shape_pop_id ()));
  NcGalaxyWLObsEllipConv ellip_conv       = nc_galaxy_shape_factor_get_ellip_conv (gsf);
  guint i;

  g_assert_nonnull (pop);
  g_assert_cmpuint (rotations->len, ==, 2 * data_array->len);

  ncm_stats_vec_reset (self->obs_stats, TRUE);

  for (i = 0; i < data_array->len; i++)
  {
    NcGalaxyShapeFactorData *data_i = g_ptr_array_index (data_array, i);
    const gdouble rot_re            = g_array_index (rotations, gdouble, 2 * i);
    const gdouble rot_im            = g_array_index (rotations, gdouble, 2 * i + 1);
    const gdouble e1                = data_i->epsilon_obs_1;
    const gdouble e2                = data_i->epsilon_obs_2;
    const gdouble std_noise         = data_i->std_noise;
    const gdouble m                 = data_i->m;
    const gdouble c1                = data_i->c1;
    const gdouble c2                = data_i->c2;
    complex double e_o              = e1 + I * e2;
    complex double c_o              = c1 + I * c2;
    gdouble std_shape;
    complex double hat_g, hat_c;
    gdouble var_tot, weight;
    NcmComplex rot;

    nc_galaxy_shape_pop_prepare (pop, data_i->pop_data);
    std_shape = nc_galaxy_shape_pop_e_rms (pop, data_i->pop_data);
    var_tot   = std_shape * std_shape + std_noise * std_noise;
    weight    = 1.0 / var_tot;

    /* Rotate the data-frame ellipticity and additive bias into the tangential/cross frame */
    rot   = nc_wl_ellipticity_celestial_to_frame (data_i->coord, rot_re + I * rot_im);
    hat_g = e_o * conj (rot);
    hat_c = c_o * conj (rot);

    ncm_stats_vec_set (self->obs_stats, 0, creal (hat_g));
    ncm_stats_vec_set (self->obs_stats, 1, cimag (hat_g));
    ncm_stats_vec_set (self->obs_stats, 2, std_shape * std_shape);
    ncm_stats_vec_set (self->obs_stats, 3, m);
    ncm_stats_vec_set (self->obs_stats, 4, creal (hat_c));
    ncm_stats_vec_set (self->obs_stats, 5, cimag (hat_c));

    ncm_stats_vec_update_weight (self->obs_stats, weight);
  }

  {
    const gdouble mean_gt          = ncm_stats_vec_get_mean (self->obs_stats, 0);
    const gdouble mean_gx          = ncm_stats_vec_get_mean (self->obs_stats, 1);
    const gdouble mean_sigma_true2 = ncm_stats_vec_get_mean (self->obs_stats, 2);
    const gdouble mean_m           = ncm_stats_vec_get_mean (self->obs_stats, 3);
    const gdouble mean_ct          = ncm_stats_vec_get_mean (self->obs_stats, 4);
    const gdouble mean_cx          = ncm_stats_vec_get_mean (self->obs_stats, 5);
    const gdouble R                = 1.0 - mean_sigma_true2;

    switch (ellip_conv)
    {
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
        /* epsilon convention: <epsilon> = g, no responsivity factor R */
        *gt      = (mean_gt - mean_ct) / (1.0 + mean_m);
        *gx      = (mean_gx - mean_cx) / (1.0 + mean_m);
        *sigma_t = ncm_stats_vec_get_sd (self->obs_stats, 0) / (1.0 + mean_m);
        *sigma_x = ncm_stats_vec_get_sd (self->obs_stats, 1) / (1.0 + mean_m);
        *rho     = ncm_stats_vec_get_cor (self->obs_stats, 0, 1);
        break;
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
        /* distortion convention: <chi> = 2 R g, responsivity applies to both components */
        *gt      = (0.5 * mean_gt / R - mean_ct) / (1.0 + mean_m);
        *gx      = (0.5 * mean_gx / R - mean_cx) / (1.0 + mean_m);
        *sigma_t = 0.5 * ncm_stats_vec_get_sd (self->obs_stats, 0) / R / (1.0 + mean_m);
        *sigma_x = 0.5 * ncm_stats_vec_get_sd (self->obs_stats, 1) / R / (1.0 + mean_m);
        *rho     = ncm_stats_vec_get_cor (self->obs_stats, 0, 1);
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
        break;                   /* LCOV_EXCL_LINE */
    }
  }
}

/**
 * nc_galaxy_shape_factor_eval_marginal: (virtual eval_marginal)
 * @gsf: a #NcGalaxyShapeFactor
 * @pop: the #NcGalaxyShapePop intrinsic-ellipticity model
 * @data: a #NcGalaxyShapeFactorData
 * @g_1: real part of the reduced shear (bias applied)
 * @g_2: imaginary part of the reduced shear (bias applied)
 * @epsilon_obs_1: observed ellipticity component 1 (same frame as @g_1, @g_2)
 * @epsilon_obs_2: observed ellipticity component 2 (same frame as @g_1, @g_2)
 *
 * Evaluates the intrinsic-ellipticity marginal
 * P(epsilon_obs | g) = int d^2chi_I P_pop(chi_I) N_2(epsilon_obs - f_g(chi_I); std_noise^2)
 * using the subclass evaluation strategy. The shear and the observed
 * ellipticity must be expressed in the same frame; @data supplies the
 * per-galaxy noise (std_noise) and the population fragment (pop_data).
 *
 * Returns: the marginal probability density.
 */
gdouble
nc_galaxy_shape_factor_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  return NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsf)->eval_marginal (gsf, pop, data, g_1, g_2, epsilon_obs_1, epsilon_obs_2);
}

/**
 * nc_galaxy_shape_factor_eval_ln_marginal: (virtual eval_ln_marginal)
 * @gsf: a #NcGalaxyShapeFactor
 * @pop: the #NcGalaxyShapePop intrinsic-ellipticity model
 * @data: a #NcGalaxyShapeFactorData
 * @g_1: real part of the reduced shear (bias applied)
 * @g_2: imaginary part of the reduced shear (bias applied)
 * @epsilon_obs_1: observed ellipticity component 1 (same frame as @g_1, @g_2)
 * @epsilon_obs_2: observed ellipticity component 2 (same frame as @g_1, @g_2)
 *
 * Logarithmic counterpart of nc_galaxy_shape_factor_eval_marginal().
 *
 * Returns: the natural logarithm of the marginal probability density.
 */
gdouble
nc_galaxy_shape_factor_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  return NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsf)->eval_ln_marginal (gsf, pop, data, g_1, g_2, epsilon_obs_1, epsilon_obs_2);
}

/**
 * nc_galaxy_shape_factor_eval_ln_marginal_array:
 * @gsf: a #NcGalaxyShapeFactor
 * @pop: the #NcGalaxyShapePop intrinsic-ellipticity model
 * @data: a #NcGalaxyShapeFactorData
 * @g_1: real parts of the reduced shears
 * @g_2: imaginary parts of the reduced shears
 * @epsilon_obs_1: observed ellipticities, component 1
 * @epsilon_obs_2: observed ellipticities, component 2
 * @out: a #NcmVector receiving the results
 *
 * Sets element $i$ of @out to nc_galaxy_shape_factor_eval_ln_marginal() at element $i$
 * of @g_1, @g_2, @epsilon_obs_1 and @epsilon_obs_2, in order and in one call, for the
 * galaxy of @data. The five vectors must have the same length. It evaluates many points
 * without a call per point from the caller, for example to time the evaluation of a
 * class from a language binding without the cost of the binding.
 */
void
nc_galaxy_shape_factor_eval_ln_marginal_array (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                               const NcmVector *g_1, const NcmVector *g_2,
                                               const NcmVector *epsilon_obs_1, const NcmVector *epsilon_obs_2,
                                               NcmVector *out)
{
  NcGalaxyShapeFactorClass *klass = NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsf);
  const guint n                   = ncm_vector_len (out);
  guint i;

  g_assert_cmpuint (ncm_vector_len (g_1), ==, n);
  g_assert_cmpuint (ncm_vector_len (g_2), ==, n);
  g_assert_cmpuint (ncm_vector_len (epsilon_obs_1), ==, n);
  g_assert_cmpuint (ncm_vector_len (epsilon_obs_2), ==, n);

  for (i = 0; i < n; i++)
    ncm_vector_set (out, i,
                    klass->eval_ln_marginal (gsf, pop, data,
                                             ncm_vector_get (g_1, i), ncm_vector_get (g_2, i),
                                             ncm_vector_get (epsilon_obs_1, i), ncm_vector_get (epsilon_obs_2, i)));
}

