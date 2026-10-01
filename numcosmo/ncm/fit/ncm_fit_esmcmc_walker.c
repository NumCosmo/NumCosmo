/***************************************************************************
 *            ncm_fit_esmcmc_walker.c
 *
 *  Wed March 16 13:07:31 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_fit_esmcmc_walker.c
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcmFitESMCMCWalker:
 *
 * Abstract proposal move of #NcmFitESMCMC.
 *
 * An ensemble sampler moves $L$ walkers $X_1, \dots, X_L$, each a point of the
 * free-parameter space, whose joint target is the product of the posterior at each
 * walker. The move of walker $k$ uses the positions of the others (Goodman and Weare,
 * [Ensemble samplers with affine invariance](https://msp.org/camcos/2010/5-1/camcos-v5-n1-p04-s.pdf)).
 * A walker class draws the proposals (ncm_fit_esmcmc_walker_setup() and
 * ncm_fit_esmcmc_walker_step()) and gives the proposal factor $q$ of the acceptance
 * probability $\min(1, q\,L^\star/L)$ (ncm_fit_esmcmc_walker_prob_norm()).
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/fit/ncm_fit_esmcmc_walker.h"
#include "ncm/core/ncm_serialize.h"

enum
{
  PROP_0,
  PROP_SIZE,
  PROP_NPARAMS,
};

G_DEFINE_ABSTRACT_TYPE (NcmFitESMCMCWalker, ncm_fit_esmcmc_walker, G_TYPE_OBJECT)

static void
ncm_fit_esmcmc_walker_init (NcmFitESMCMCWalker *walker)
{
}

static void
ncm_fit_esmcmc_walker_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_esmcmc_walker_parent_class)->finalize (object);
}

static void
ncm_fit_esmcmc_walker_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFitESMCMCWalker *walker = NCM_FIT_ESMCMC_WALKER (object);

  g_return_if_fail (NCM_IS_FIT_ESMCMC_WALKER (object));

  switch (prop_id)
  {
    case PROP_SIZE:
      ncm_fit_esmcmc_walker_set_size (walker, g_value_get_uint (value));
      break;
    case PROP_NPARAMS:
      ncm_fit_esmcmc_walker_set_nparams (walker, g_value_get_uint (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_fit_esmcmc_walker_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFitESMCMCWalker *walker = NCM_FIT_ESMCMC_WALKER (object);

  g_return_if_fail (NCM_IS_FIT_ESMCMC_WALKER (object));

  switch (prop_id)
  {
    case PROP_SIZE:
      g_value_set_uint (value, ncm_fit_esmcmc_walker_get_size (walker));
      break;
    case PROP_NPARAMS:
      g_value_set_uint (value, ncm_fit_esmcmc_walker_get_nparams (walker));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

/* LCOV_EXCL_START */

static void
_ncm_fit_esmcmc_walker_set_size (NcmFitESMCMCWalker *walker, guint size)
{
  g_error ("method set_size not implemented by %s.", G_OBJECT_TYPE_NAME (walker));
}

static guint
_ncm_fit_esmcmc_walker_get_size (NcmFitESMCMCWalker *walker)
{
  g_error ("method get_size not implemented by %s.", G_OBJECT_TYPE_NAME (walker));

  return 0;
}

static void
_ncm_fit_esmcmc_walker_set_nparams (NcmFitESMCMCWalker *walker, guint nparams)
{
  g_error ("method set_nparams not implemented by %s.", G_OBJECT_TYPE_NAME (walker));
}

static guint
_ncm_fit_esmcmc_walker_get_nparams (NcmFitESMCMCWalker *walker)
{
  g_error ("method get_nparams not implemented by %s.", G_OBJECT_TYPE_NAME (walker));

  return 0;
}

static void
_ncm_fit_esmcmc_walker_setup (NcmFitESMCMCWalker *walker, NcmMSet *mset, GPtrArray *theta, GPtrArray *m2lnL, guint ki, guint kf, NcmRNG *rng)
{
  g_error ("method setup not implemented by %s.", G_OBJECT_TYPE_NAME (walker));
}

static void
_ncm_fit_esmcmc_walker_step (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k)
{
  g_error ("method step not implemented by %s.", G_OBJECT_TYPE_NAME (walker));
}

static gdouble
_ncm_fit_esmcmc_walker_prob (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star)
{
  g_error ("method prob not implemented by %s.", G_OBJECT_TYPE_NAME (walker));

  return 0.0;
}

static gdouble
_ncm_fit_esmcmc_walker_prob_norm (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k)
{
  g_error ("method prob_norm not implemented by %s.", G_OBJECT_TYPE_NAME (walker));

  return 0.0;
}

static void
_ncm_fit_esmcmc_walker_clean (NcmFitESMCMCWalker *walker, guint ki, guint kf)
{
  g_error ("method clean not implemented by %s.", G_OBJECT_TYPE_NAME (walker));
}

static const gchar *
_ncm_fit_esmcmc_walker_desc (NcmFitESMCMCWalker *walker)
{
  g_error ("method desc not implemented by %s.", G_OBJECT_TYPE_NAME (walker));

  return NULL;
}

/* A walker with nothing to tune reports no options, which is not a missing method. */
static const gchar *
_ncm_fit_esmcmc_walker_opts (NcmFitESMCMCWalker *walker)
{
  return NULL;
}

static void
_ncm_fit_esmcmc_walker_start_run (NcmFitESMCMCWalker *walker, gboolean initial, guint exploration_done)
{
  /* Nothing to prepare by default. */
}

static void
_ncm_fit_esmcmc_walker_end_run (NcmFitESMCMCWalker *walker)
{
  /* Nothing to release by default. */
}

static gboolean
_ncm_fit_esmcmc_walker_is_markovian (NcmFitESMCMCWalker *walker)
{
  /* A walker that does not modify the acceptance is Markovian at every iteration. */
  return TRUE;
}

/* LCOV_EXCL_STOP */

static void
ncm_fit_esmcmc_walker_class_init (NcmFitESMCMCWalkerClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = ncm_fit_esmcmc_walker_set_property;
  object_class->get_property = ncm_fit_esmcmc_walker_get_property;

  object_class->finalize = ncm_fit_esmcmc_walker_finalize;

  g_object_class_install_property (object_class,
                                   PROP_SIZE,
                                   g_param_spec_uint ("size",
                                                      NULL,
                                                      "Number of walkers",
                                                      1, G_MAXUINT, 100,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NPARAMS,
                                   g_param_spec_uint ("nparams",
                                                      NULL,
                                                      "Number of parameters",
                                                      1, G_MAXUINT, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->set_size     = _ncm_fit_esmcmc_walker_set_size;
  klass->get_size     = _ncm_fit_esmcmc_walker_get_size;
  klass->set_nparams  = _ncm_fit_esmcmc_walker_set_nparams;
  klass->get_nparams  = _ncm_fit_esmcmc_walker_get_nparams;
  klass->setup        = _ncm_fit_esmcmc_walker_setup;
  klass->step         = _ncm_fit_esmcmc_walker_step;
  klass->prob         = _ncm_fit_esmcmc_walker_prob;
  klass->prob_norm    = _ncm_fit_esmcmc_walker_prob_norm;
  klass->clean        = _ncm_fit_esmcmc_walker_clean;
  klass->desc         = _ncm_fit_esmcmc_walker_desc;
  klass->opts         = _ncm_fit_esmcmc_walker_opts;
  klass->start_run    = &_ncm_fit_esmcmc_walker_start_run;
  klass->end_run      = &_ncm_fit_esmcmc_walker_end_run;
  klass->is_markovian = &_ncm_fit_esmcmc_walker_is_markovian;
}

/**
 * ncm_fit_esmcmc_walker_ref:
 * @walker: a #NcmFitESMCMCWalker
 *
 * Increases the reference count of @walker atomically.
 *
 * Returns: (transfer full): @walker.
 */
NcmFitESMCMCWalker *
ncm_fit_esmcmc_walker_ref (NcmFitESMCMCWalker *walker)
{
  return g_object_ref (walker);
}

/**
 * ncm_fit_esmcmc_walker_free:
 * @walker: a #NcmFitESMCMCWalker
 *
 * Decreases the reference count of @walker atomically.
 *
 */
void
ncm_fit_esmcmc_walker_free (NcmFitESMCMCWalker *walker)
{
  g_object_unref (walker);
}

/**
 * ncm_fit_esmcmc_walker_clear:
 * @walker: a #NcmFitESMCMCWalker
 *
 * Decreases the reference count of *@walker and sets it to %NULL.
 *
 */
void
ncm_fit_esmcmc_walker_clear (NcmFitESMCMCWalker **walker)
{
  g_clear_object (walker);
}

/**
 * ncm_fit_esmcmc_walker_set_size: (virtual set_size)
 * @walker: a #NcmFitESMCMCWalker
 * @size: number of walkers
 *
 * Sets the number of walkers.
 *
 */
void
ncm_fit_esmcmc_walker_set_size (NcmFitESMCMCWalker *walker, guint size)
{
  NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->set_size (walker, size);
}

/**
 * ncm_fit_esmcmc_walker_get_size: (virtual get_size)
 * @walker: a #NcmFitESMCMCWalker
 *
 * Returns: the number of walkers
 *
 */
guint
ncm_fit_esmcmc_walker_get_size (NcmFitESMCMCWalker *walker)
{
  return NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->get_size (walker);
}

/**
 * ncm_fit_esmcmc_walker_set_nparams: (virtual set_nparams)
 * @walker: a #NcmFitESMCMCWalker
 * @nparams: number of parameters
 *
 * Sets the number of free parameters; #NcmFitESMCMC sets it to the number of free
 * parameters of its fit.
 *
 */
void
ncm_fit_esmcmc_walker_set_nparams (NcmFitESMCMCWalker *walker, guint nparams)
{
  NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->set_nparams (walker, nparams);
}

/**
 * ncm_fit_esmcmc_walker_get_nparams: (virtual get_nparams)
 * @walker: a #NcmFitESMCMCWalker
 *
 * Returns: the number of free parameters
 *
 */
guint
ncm_fit_esmcmc_walker_get_nparams (NcmFitESMCMCWalker *walker)
{
  return NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->get_nparams (walker);
}

/**
 * ncm_fit_esmcmc_walker_setup: (virtual setup)
 * @walker: a #NcmFitESMCMCWalker
 * @mset: a #NcmMSet
 * @theta: (element-type NcmVector): array of walkers positions
 * @m2lnL: (element-type NcmVector): array of walkers $-2\ln(L)$
 * @ki: first walker index
 * @kf: last walker index
 * @rng: a #NcmRNG
 *
 * Draws the random numbers of the moves of the walkers @ki to @kf - 1.
 *
 */
void
ncm_fit_esmcmc_walker_setup (NcmFitESMCMCWalker *walker, NcmMSet *mset, GPtrArray *theta, GPtrArray *m2lnL, guint ki, guint kf, NcmRNG *rng)
{
  NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->setup (walker, mset, theta, m2lnL, ki, kf, rng);
}

/**
 * ncm_fit_esmcmc_walker_step: (virtual step)
 * @walker: a #NcmFitESMCMCWalker
 * @theta: (element-type NcmVector): array of walkers positions
 * @m2lnL: (element-type NcmVector): array of walkers $-2\ln(L)$
 * @thetastar: a #NcmVector
 * @k: index of the walker to move
 *
 * Computes the proposal of walker @k in @thetastar, from the numbers drawn by
 * ncm_fit_esmcmc_walker_setup().
 *
 */
void
ncm_fit_esmcmc_walker_step (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k)
{
  NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->step (walker, theta, m2lnL, thetastar, k);
}

/**
 * ncm_fit_esmcmc_walker_prob: (virtual prob)
 * @walker: a #NcmFitESMCMCWalker
 * @theta: (element-type NcmVector): array of walkers positions
 * @m2lnL: (element-type NcmVector): array of walkers $-2\ln(L)$
 * @thetastar: a #NcmVector
 * @k: index of the walker to move
 * @m2lnL_cur: current value of $-2\ln(L)$
 * @m2lnL_star: proposed value for $-2\ln(L^\star)$
 *
 * Computes the acceptance ratio $q\,L^\star/L$ of moving walker @k from $-2\ln L$ =
 * @m2lnL_cur to @thetastar with $-2\ln L^\star$ = @m2lnL_star.
 *
 * Returns: the acceptance ratio, before its minimum with one
 */
gdouble
ncm_fit_esmcmc_walker_prob (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star)
{
  return NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->prob (walker, theta, m2lnL, thetastar, k, m2lnL_cur, m2lnL_star);
}

/**
 * ncm_fit_esmcmc_walker_prob_norm: (virtual prob_norm)
 * @walker: a #NcmFitESMCMCWalker
 * @theta: (element-type NcmVector): array of walkers positions
 * @m2lnL: (element-type NcmVector): array of walkers $-2\ln(L)$
 * @thetastar: a #NcmVector
 * @k: index of the walker to move
 *
 * Computes $\ln q$, the proposal factor of the acceptance ratio of moving walker @k to
 * @thetastar (e.g. $(d - 1)\ln z$ for the stretch move).
 *
 * Returns: $\ln q$
 */
gdouble
ncm_fit_esmcmc_walker_prob_norm (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k)
{
  return NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->prob_norm (walker, theta, m2lnL, thetastar, k);
}

/**
 * ncm_fit_esmcmc_walker_clean: (virtual clean)
 * @walker: a #NcmFitESMCMCWalker
 * @ki: first walker index
 * @kf: last walker index
 *
 * Releases what the moves of the walkers @ki to @kf - 1 needed.
 *
 */
void
ncm_fit_esmcmc_walker_clean (NcmFitESMCMCWalker *walker, guint ki, guint kf)
{
  NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->clean (walker, ki, kf);
}

/**
 * ncm_fit_esmcmc_walker_desc: (virtual desc)
 * @walker: a #NcmFitESMCMCWalker
 *
 * Returns: (transfer none): walker description.
 */
const gchar *
ncm_fit_esmcmc_walker_desc (NcmFitESMCMCWalker *walker)
{
  return NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->desc (walker);
}

/**
 * ncm_fit_esmcmc_walker_opts: (virtual opts)
 * @walker: a #NcmFitESMCMCWalker
 *
 * The walker's tunable settings as a colon separated list of `name=value' pairs, the
 * companion of ncm_fit_esmcmc_walker_desc(), which names the structure alone. Two runs of
 * the same walker differ here and nowhere else, so a catalog records both.
 *
 * Returns: (transfer none) (nullable): the options, or %NULL when the walker has none.
 */
const gchar *
ncm_fit_esmcmc_walker_opts (NcmFitESMCMCWalker *walker)
{
  return NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->opts (walker);
}

/**
 * ncm_fit_esmcmc_walker_start_run: (virtual start_run)
 * @walker: a #NcmFitESMCMCWalker
 * @initial: whether the chain has no Markovian rows yet
 * @exploration_done: iterations of a non-Markovian phase already recorded in the catalog
 *
 * Called by #NcmFitESMCMC once when a run starts, before the first iteration. @initial is
 * TRUE when the run starts from the initial ensemble, or resumes a catalog whose
 * #NcmMSetCatalog:markovian-id lies beyond its last row, i.e. an interrupted exploration
 * phase; the walker then resumes the phase with @exploration_done iterations already
 * counted. @initial is FALSE when the catalog already holds Markovian rows, and the walker
 * must not start any non-Markovian phase.
 *
 */
void
ncm_fit_esmcmc_walker_start_run (NcmFitESMCMCWalker *walker, gboolean initial, guint exploration_done)
{
  NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->start_run (walker, initial, exploration_done);
}

/**
 * ncm_fit_esmcmc_walker_end_run: (virtual end_run)
 * @walker: a #NcmFitESMCMCWalker
 *
 * Called by #NcmFitESMCMC once when a run ends.
 *
 */
void
ncm_fit_esmcmc_walker_end_run (NcmFitESMCMCWalker *walker)
{
  NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->end_run (walker);
}

/**
 * ncm_fit_esmcmc_walker_is_markovian: (virtual is_markovian)
 * @walker: a #NcmFitESMCMCWalker
 *
 * Whether the acceptance of the last completed iteration was the exact Metropolis-Hastings
 * one for every walker. FALSE means some walker was moved by a rule that does not satisfy
 * detailed balance (an exploration phase), and #NcmFitESMCMC then moves the catalog's
 * #NcmMSetCatalog:markovian-id past that iteration.
 *
 * Returns: TRUE if the last iteration was Markovian.
 */
gboolean
ncm_fit_esmcmc_walker_is_markovian (NcmFitESMCMCWalker *walker)
{
  return NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->is_markovian (walker);
}

