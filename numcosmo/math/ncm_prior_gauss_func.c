/***************************************************************************
 *            ncm_prior_gauss_func.c
 *
 *  Wed August 03 16:27:00 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>>
 ****************************************************************************/
/*
 * ncm_prior_gauss_func.c
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
 * SECTION:ncm_prior_gauss_func
 * @title: NcmPriorGaussFunc
 * @short_description: a gaussian prior on a parameter
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_prior_gauss_func.h"

enum
{
  PROP_0,
  PROP_MEAN_FUNC,
};

G_DEFINE_TYPE (NcmPriorGaussFunc, ncm_prior_gauss_func, NCM_TYPE_PRIOR_GAUSS)

static void
ncm_prior_gauss_func_init (NcmPriorGaussFunc *pgf)
{
  pgf->mean_func = NULL;
}

static void
_ncm_prior_gauss_func_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmPriorGaussFunc *pgf = NCM_PRIOR_GAUSS_FUNC (object);
  g_return_if_fail (NCM_IS_PRIOR_GAUSS_FUNC (object));

  switch (prop_id)
  {
    case PROP_MEAN_FUNC:
      ncm_mset_func_clear (&pgf->mean_func);
      pgf->mean_func = g_value_dup_object (value);

      g_assert (ncm_mset_func_is_scalar (pgf->mean_func));
      g_assert_cmpint (ncm_mset_func_get_nvar (pgf->mean_func), <=, 1);
      
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_prior_gauss_func_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmPriorGaussFunc *pgf = NCM_PRIOR_GAUSS_FUNC (object);
  g_return_if_fail (NCM_IS_PRIOR_GAUSS_FUNC (object));

  switch (prop_id)
  {
    case PROP_MEAN_FUNC:
      g_value_set_object (value, pgf->mean_func);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_prior_gauss_func_dispose (GObject *object)
{
  NcmPriorGaussFunc *pgf = NCM_PRIOR_GAUSS_FUNC (object);

  ncm_mset_func_clear (&pgf->mean_func);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_prior_gauss_func_parent_class)->dispose (object);
}

static void
_ncm_prior_gauss_func_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_prior_gauss_func_parent_class)->finalize (object);
}

static gdouble _ncm_prior_gauss_func_mean (NcmPriorGauss *pg, NcmMSet *mset);

static void
ncm_prior_gauss_func_class_init (NcmPriorGaussFuncClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmPriorGaussClass *pg_class = NCM_PRIOR_GAUSS_CLASS (klass);

  object_class->set_property = &_ncm_prior_gauss_func_set_property;
  object_class->get_property = &_ncm_prior_gauss_func_get_property;
  object_class->dispose      = &_ncm_prior_gauss_func_dispose;
  object_class->finalize     = &_ncm_prior_gauss_func_finalize;

  g_object_class_install_property (object_class,
                                   PROP_MEAN_FUNC,
                                   g_param_spec_object ("mean-func",
                                                        NULL,
                                                        "mean function",
                                                        NCM_TYPE_MSET_FUNC,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  pg_class->mean = &_ncm_prior_gauss_func_mean;
}

static gdouble 
_ncm_prior_gauss_func_mean (NcmPriorGauss *pg, NcmMSet *mset)
{
  NcmPriorGaussFunc *pgf = NCM_PRIOR_GAUSS_FUNC (pg);
  return ncm_mset_func_eval1 (pgf->mean_func, mset, pg->var);
}

/**
 * ncm_prior_gauss_func_new:
 * @mean_func: a #NcmMSetFunc
 * @mu: mean
 * @sigma: standard deviation
 * @var: variable
 * 
 * Creates a new Gaussiam prior for parameter @pid of model @mid.
 * 
 * Returns: (transfer full): @pgf.
 */
NcmPriorGaussFunc *
ncm_prior_gauss_func_new (NcmMSetFunc *mean_func, gdouble mu, gdouble sigma, gdouble var)
{
  NcmPriorGaussFunc *pgf = g_object_new (NCM_TYPE_PRIOR_GAUSS_FUNC, 
                                         "mean-func", mean_func,
                                         "mu", mu,
                                         "sigma", sigma,
                                         "variable", var,
                                         NULL);
  
  return pgf;
}

/**
 * ncm_prior_gauss_func_ref:
 * @pgf: a #NcmPriorGaussFunc
 * 
 * Increases the reference count of @pgf atomically.
 * 
 * Returns: (transfer full): @pgf.
 */
NcmPriorGaussFunc *
ncm_prior_gauss_func_ref (NcmPriorGaussFunc *pgf)
{
  return g_object_ref (pgf);
}

/**
 * ncm_prior_gauss_func_free:
 * @pgf: a #NcmPriorGaussFunc
 * 
 * Decreases the reference count of @pgf atomically.
 * 
 */
void 
ncm_prior_gauss_func_free (NcmPriorGaussFunc *pgf)
{
  g_object_unref (pgf);
}

/**
 * ncm_prior_gauss_func_clear:
 * @pgf: a #NcmPriorGaussFunc
 * 
 * Decreases the reference count of *@pgf and sets *@pgf to NULL.
 * 
 */
void 
ncm_prior_gauss_func_clear (NcmPriorGaussFunc **pgf)
{
  g_clear_object (pgf);
}
