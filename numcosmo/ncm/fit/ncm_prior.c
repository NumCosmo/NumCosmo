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
 * NcmPrior:
 *
 * Base class for prior distributions.
 *
 * This object defines a base class for priors used by #NcmLikelihood. These objects
 * describe prior distributions applicable to parameters or any derived quantity. A
 * prior returns one of two quantities, see ncm_prior_is_m2lnL():
 *
 * 1. $-2\ln P$, added to $-2\ln L$ as it is;
 * 2. $f$ such that $-2\ln P = f^2$, added to $-2\ln L$ as $f^2$.
 *
 * The second form is also a least-squares residual, so those priors can be used by the
 * least-squares fits as well as by any other analysis.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/fit/ncm_prior.h"

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
ncm_prior_class_init (NcmPriorClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->constructed = &_ncm_prior_constructed;

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
 * Returns: TRUE if the prior returns $-2\ln P$ and FALSE if it returns $f$ such that
 * $-2\ln P = f^2$.
 */
gboolean
ncm_prior_is_m2lnL (NcmPrior *prior)
{
  return NCM_PRIOR_GET_CLASS (prior)->is_m2lnL;
}

