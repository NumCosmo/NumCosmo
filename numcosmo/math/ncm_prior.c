/***************************************************************************
 *            ncm_prior.c
 *
 *  Wed August 03 10:08:51 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_prior.c
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
 * SECTION:ncm_prior
 * @title: NcmPrior
 * @short_description: Base class for prior distributions.
 *
 * This object defines a base class for priors used by NcmLikelihood. These objects
 * describe various prior distributions applicable to parameters or any derived
 * quantity. Two types of priors are supported:
 *
 * 1. Priors returning $-2\ln(L_\mathrm{prior})$
 * 2. Priors returning $f$ such that $-2\ln(P_\mathrm{prior}) = f^2$
 *
 * The second type is essential when conducting least-squares based analysis but
 * can also be used in any other type of analysis.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_prior.h"

G_DEFINE_TYPE (NcmPrior, ncm_prior, NCM_TYPE_MSET_FUNC)

static void
ncm_prior_init (NcmPrior *prior)
{
}

static void
_ncm_prior_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_prior_parent_class)->constructed (object);

  NcmMSetFunc *func = NCM_MSET_FUNC (object);

  ncm_mset_func_set_meta (func, NULL, NULL, NULL, NULL, 0, 1);
}

static void
_ncm_prior_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_prior_parent_class)->finalize (object);
}

static void
ncm_prior_class_init (NcmPriorClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  /* NcmMSetFuncClass *func_class = NCM_MSET_FUNC_CLASS (klass); */

  object_class->constructed = &_ncm_prior_constructed;
  object_class->finalize    = &_ncm_prior_finalize;

  klass->is_m2lnL = FALSE;
}

/**
 * ncm_prior_ref:
 * @prior: a #NcmPrior
 *
 * Increases the reference count of @prior atomically.
 *
 * Returns: (transfer full): @prior.
 */
NcmPrior *
ncm_prior_ref (NcmPrior *prior)
{
  return g_object_ref (prior);
}

/**
 * ncm_prior_free:
 * @prior: a #NcmPrior
 *
 * Decreases the reference count of @prior atomically.
 *
 */
void
ncm_prior_free (NcmPrior *prior)
{
  g_object_unref (prior);
}

/**
 * ncm_prior_clear:
 * @prior: a #NcmPrior
 *
 * Decreases the reference count of *@prior and sets *@prior to NULL.
 *
 */
void
ncm_prior_clear (NcmPrior **prior)
{
  g_clear_object (prior);
}

/**
 * ncm_prior_is_m2lnL:
 * @prior: a #NcmPrior
 *
 * Returns: TRUE if the prior calculates $-2\ln(L_\mathrm{prior})$ and FALSE
 * if it returns $f$ such that $-2\ln(L_\mathrm{prior}) = f^2$.
 */
gboolean
ncm_prior_is_m2lnL (NcmPrior *prior)
{
  return NCM_PRIOR_GET_CLASS (prior)->is_m2lnL;
}

