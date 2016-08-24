/***************************************************************************
 *            ncm_prior.c
 *
 *  Wed August 03 10:08:51 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_prior.c
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
 * SECTION:ncm_prior
 * @title: NcmPrior
 * @short_description: A prior for NcmLikelihood
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_prior.h"

G_DEFINE_TYPE (NcmPrior, ncm_prior, NCM_TYPE_MSET_FUNC);

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

  func->dim  = 1;
  func->nvar = 0;
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
  GObjectClass* object_class   = G_OBJECT_CLASS (klass);
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
 * Returns true if the prior calculates $-2\ln(L_\mathrm{prior})$ and false
 * if it returns $f$ such that $-2\ln(L_\mathrm{prior}) = f^2$.
 * 
 */
gboolean 
ncm_prior_is_m2lnL (NcmPrior *prior)
{
  return NCM_PRIOR_GET_CLASS (prior)->is_m2lnL;
}
