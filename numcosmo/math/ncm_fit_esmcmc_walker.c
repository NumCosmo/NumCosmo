/***************************************************************************
 *            ncm_fit_esmcmc_walker.c
 *
 *  Wed March 16 13:07:31 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_fit_esmcmc_walker.c
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_fit_esmcmc_walker
 * @title: NcmFitESMCMCWalker
 * @short_description: Ensemble sampler Markov Chain Monte Carlo walker class.
 *
 * Abstract class for implementing walkers for #NcmFitESMCMC.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fit_esmcmc_walker.h"
#include "math/ncm_serialize.h"

enum
{
  PROP_0,
  PROP_SIZE,
  PROP_NPARAMS,
  
};

G_DEFINE_ABSTRACT_TYPE (NcmFitESMCMCWalker, ncm_fit_esmcmc_walker, G_TYPE_OBJECT);

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
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
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
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void _ncm_fit_esmcmc_walker_set_size (NcmFitESMCMCWalker *walker, guint size) { g_error ("_ncm_fit_esmcmc_walker_set_size: method not implemented."); }
static guint _ncm_fit_esmcmc_walker_get_size (NcmFitESMCMCWalker *walker) { g_error ("_ncm_fit_esmcmc_walker_get_size: method not implemented."); return 0; }
static void _ncm_fit_esmcmc_walker_set_nparams (NcmFitESMCMCWalker *walker, guint nparams) { g_error ("_ncm_fit_esmcmc_walker_set_nparams: method not implemented."); }
static guint _ncm_fit_esmcmc_walker_get_nparams (NcmFitESMCMCWalker *walker) { g_error ("_ncm_fit_esmcmc_walker_get_nparams: method not implemented."); return 0; }
static void _ncm_fit_esmcmc_walker_setup (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, guint ki, guint kf, NcmRNG *rng) { g_error ("_ncm_fit_esmcmc_walker_setup: method not implemented."); }
static void _ncm_fit_esmcmc_walker_step (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k) { g_error ("_ncm_fit_esmcmc_walker_step: method not implemented."); }
static gdouble _ncm_fit_esmcmc_walker_prob (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star) { g_error ("_ncm_fit_esmcmc_walker_prob: method not implemented."); return 0.0; }
static gdouble _ncm_fit_esmcmc_walker_prob_norm (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k) { g_error ("_ncm_fit_esmcmc_walker_prob_norm: method not implemented."); return 0.0; }
static void _ncm_fit_esmcmc_walker_clean (NcmFitESMCMCWalker *walker, guint ki, guint kf) { g_error ("_ncm_fit_esmcmc_walker_clean: method not implemented."); }
static const gchar *_ncm_fit_esmcmc_walker_desc (NcmFitESMCMCWalker *walker) { g_error ("_ncm_fit_esmcmc_walker_desc: method not implemented."); return NULL; }

static void
ncm_fit_esmcmc_walker_class_init (NcmFitESMCMCWalkerClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

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

  klass->set_size    = _ncm_fit_esmcmc_walker_set_size;
  klass->get_size    = _ncm_fit_esmcmc_walker_get_size;
  klass->set_nparams = _ncm_fit_esmcmc_walker_set_nparams;
  klass->get_nparams = _ncm_fit_esmcmc_walker_get_nparams;
  klass->setup       = _ncm_fit_esmcmc_walker_setup;
  klass->step        = _ncm_fit_esmcmc_walker_step;
  klass->prob        = _ncm_fit_esmcmc_walker_prob;
  klass->prob_norm   = _ncm_fit_esmcmc_walker_prob_norm;
  klass->clean       = _ncm_fit_esmcmc_walker_clean;
  klass->desc        = _ncm_fit_esmcmc_walker_desc;
}

/**
 * ncm_fit_esmcmc_walker_new_from_name:
 * @walker_name: string which specifies the walker object to be used
 *
 * This function returns a new #NcmFitESMCMCWalker whose type is defined by @walker_name.
 *
 * Returns: A new #NcmFitESMCMCWalker.
 */
NcmFitESMCMCWalker *
ncm_fit_esmcmc_walker_new_from_name (const gchar *walker_name)
{
  GObject *obj = ncm_serialize_global_from_string (walker_name);

  if (!NCM_IS_FIT_ESMCMC_WALKER (obj))
    g_error ("ncm_fit_esmcmc_walker_new_from_name: NcmFitESMCMCWalker %s do not descend from %s.", walker_name, g_type_name (NCM_TYPE_FIT_ESMCMC_WALKER));

  return NCM_FIT_ESMCMC_WALKER (obj);
}

/**
 * ncm_fit_esmcmc_walker_ref:
 * @walker: a #NcmMSetCatalog
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
 * @walker: a #NcmMSetCatalog
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
 * @walker: a #NcmMSetCatalog
 *
 * Decrese the reference count of *@walker atomically and sets the pointer *@walker to null.
 *
 */
void 
ncm_fit_esmcmc_walker_clear (NcmFitESMCMCWalker **walker)
{
  g_clear_object (walker);
}

/**
 * ncm_fit_esmcmc_walker_set_size: (virtual set_size)
 * @walker: a #NcmMSetCatalog
 * @size: new walker's size
 * 
 * Sets the walker's size.
 *
 */
void
ncm_fit_esmcmc_walker_set_size (NcmFitESMCMCWalker *walker, guint size)
{
  NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->set_size (walker, size);
}

/**
 * ncm_fit_esmcmc_walker_get_size: (virtual get_size)
 * @walker: a #NcmMSetCatalog
 *
 * Returns: the size of the @walker.
 *
 */
guint 
ncm_fit_esmcmc_walker_get_size (NcmFitESMCMCWalker *walker)
{
  return NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->get_size (walker);
}

/**
 * ncm_fit_esmcmc_walker_set_nparams: (virtual set_nparams)
 * @walker: a #NcmMSetCatalog
 * @nparams: number of parameters
 * 
 * Sets the number parameters of the walker.
 *
 */
void
ncm_fit_esmcmc_walker_set_nparams (NcmFitESMCMCWalker *walker, guint nparams)
{
  NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->set_nparams (walker, nparams);
}

/**
 * ncm_fit_esmcmc_walker_get_nparams: (virtual get_nparams)
 * @walker: a #NcmMSetCatalog
 *
 * Returns: the nparams of the @walker.
 *
 */
guint 
ncm_fit_esmcmc_walker_get_nparams (NcmFitESMCMCWalker *walker)
{
  return NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->get_nparams (walker);
}

/**
 * ncm_fit_esmcmc_walker_setup: (virtual setup)
 * @walker: a #NcmMSetCatalog
 * @theta: (element-type NcmVector): array of walkers positions
 * @m2lnL: (element-type NcmVector): array of walkers $-2\ln(L)$
 * @ki: first walker index
 * @kf: last walker index
 * @rng: a #NcmRNG
 *
 * Setup the walkers @ki to @kf (@kf not included). 
 *
 */
void 
ncm_fit_esmcmc_walker_setup (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, guint ki, guint kf, NcmRNG *rng)
{
  NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->setup (walker, theta, m2lnL, ki, kf, rng);
}

/**
 * ncm_fit_esmcmc_walker_step: (virtual step)
 * @walker: a #NcmMSetCatalog
 * @theta: (element-type NcmVector): array of walkers positions
 * @m2lnL: (element-type NcmVector): array of walkers $-2\ln(L)$
 * @thetastar: a #NcmVector 
 * @k: index of the walker to move
 * 
 * Move the @k-th walker and assign the new position in @thetastar.
 * 
 */
void 
ncm_fit_esmcmc_walker_step (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k)
{
  NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->step (walker, theta, m2lnL, thetastar, k);
}

/**
 * ncm_fit_esmcmc_walker_prob: (virtual prob)
 * @walker: a #NcmMSetCatalog
 * @theta: (element-type NcmVector): array of walkers positions
 * @m2lnL: (element-type NcmVector): array of walkers $-2\ln(L)$
 * @thetastar: a #NcmVector 
 * @k: index of the walker to move
 * @m2lnL_cur: current value of $-2\ln(L)$
 * @m2lnL_star: proposed value for $-2\ln(L^\star)$
 * 
 * Calculates the transition probability
 * 
 * Returns: the transition probability.
 */
gdouble
ncm_fit_esmcmc_walker_prob (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star)
{
  return NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->prob (walker, theta, m2lnL, thetastar, k, m2lnL_cur, m2lnL_star);
}

/**
 * ncm_fit_esmcmc_walker_prob_norm: (virtual prob_norm)
 * @walker: a #NcmMSetCatalog
 * @theta: (element-type NcmVector): array of walkers positions
 * @m2lnL: (element-type NcmVector): array of walkers $-2\ln(L)$
 * @thetastar: a #NcmVector 
 * @k: index of the walker to move
 * 
 * Calculates the transition probability norm, this method is used in the MPI implementation.
 * 
 * Returns: the transition probability log-norm.
 */
gdouble
ncm_fit_esmcmc_walker_prob_norm (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k)
{
  return NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->prob_norm (walker, theta, m2lnL, thetastar, k);
}

/**
 * ncm_fit_esmcmc_walker_clean: (virtual clean)
 * @walker: a #NcmMSetCatalog
 * @ki: first walker index
 * @kf: last walker index
 * 
 * Cleanup after moving walkers from @ki to @kf (@kf not included).
 * 
 */
void 
ncm_fit_esmcmc_walker_clean (NcmFitESMCMCWalker *walker, guint ki, guint kf)
{
  NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->clean (walker, ki, kf);
}

/**
 * ncm_fit_esmcmc_walker_desc: (virtual desc)
 * @walker: a #NcmMSetCatalog
 * 
 * Returns: (transfer none): walker description.
 */
const gchar * 
ncm_fit_esmcmc_walker_desc (NcmFitESMCMCWalker *walker)
{
  return NCM_FIT_ESMCMC_WALKER_GET_CLASS (walker)->desc (walker);
}
