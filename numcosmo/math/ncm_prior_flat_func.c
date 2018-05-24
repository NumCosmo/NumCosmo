/***************************************************************************
 *            ncm_prior_flat_func.c
 *
 *  Wed August 03 16:27:00 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>>
 ****************************************************************************/
/*
 * ncm_prior_flat_func.c
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
 * SECTION:ncm_prior_flat_func
 * @title: NcmPriorFlatFunc
 * @short_description: a flat prior on a parameter
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_prior_flat_func.h"

enum
{
  PROP_0,
  PROP_MEAN_FUNC,
};

G_DEFINE_TYPE (NcmPriorFlatFunc, ncm_prior_flat_func, NCM_TYPE_PRIOR_FLAT);

static void
ncm_prior_flat_func_init (NcmPriorFlatFunc *pff)
{
  pff->mean_func = NULL;
}

static void
_ncm_prior_flat_func_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmPriorFlatFunc *pff = NCM_PRIOR_FLAT_FUNC (object);
  g_return_if_fail (NCM_IS_PRIOR_FLAT_FUNC (object));

  switch (prop_id)
  {
    case PROP_MEAN_FUNC:
      ncm_mset_func_clear (&pff->mean_func);
      pff->mean_func = g_value_dup_object (value);

      g_assert (ncm_mset_func_is_scalar (pff->mean_func));
      g_assert_cmpint (ncm_mset_func_get_nvar (pff->mean_func), <=, 1);
      
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_prior_flat_func_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmPriorFlatFunc *pff = NCM_PRIOR_FLAT_FUNC (object);
  g_return_if_fail (NCM_IS_PRIOR_FLAT_FUNC (object));

  switch (prop_id)
  {
    case PROP_MEAN_FUNC:
      g_value_set_object (value, pff->mean_func);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_prior_flat_func_dispose (GObject *object)
{
  NcmPriorFlatFunc *pff = NCM_PRIOR_FLAT_FUNC (object);

  ncm_mset_func_clear (&pff->mean_func);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_prior_flat_func_parent_class)->dispose (object);
}

static void
_ncm_prior_flat_func_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_prior_flat_func_parent_class)->finalize (object);
}

static gdouble _ncm_prior_flat_func_mean (NcmPriorFlat *pf, NcmMSet *mset);

static void
ncm_prior_flat_func_class_init (NcmPriorFlatFuncClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmPriorFlatClass *pf_class = NCM_PRIOR_FLAT_CLASS (klass);

  object_class->set_property = &_ncm_prior_flat_func_set_property;
  object_class->get_property = &_ncm_prior_flat_func_get_property;
  object_class->dispose      = &_ncm_prior_flat_func_dispose;
  object_class->finalize     = &_ncm_prior_flat_func_finalize;

  g_object_class_install_property (object_class,
                                   PROP_MEAN_FUNC,
                                   g_param_spec_object ("mean-func",
                                                        NULL,
                                                        "mean function",
                                                        NCM_TYPE_MSET_FUNC,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  pf_class->mean = &_ncm_prior_flat_func_mean;
}

static gdouble 
_ncm_prior_flat_func_mean (NcmPriorFlat *pf, NcmMSet *mset)
{
  NcmPriorFlatFunc *pff = NCM_PRIOR_FLAT_FUNC (pf);
  return ncm_mset_func_eval1 (pff->mean_func, mset, pf->var);
}

/**
 * ncm_prior_flat_func_new:
 * @mean_func: a #NcmMSetFunc
 * @x_low: FIXME
 * @x_upp: FIXME
 * @scale: FIXME
 * @variable: FIXME
 * 
 * Creates a new Flat prior for parameter @pid of model @mid.
 * 
 * Returns: (transfer full): @pff.
 */
NcmPriorFlatFunc *
ncm_prior_flat_func_new (NcmMSetFunc *mean_func, gdouble x_low, gdouble x_upp, gdouble scale, gdouble variable)
{
  NcmPriorFlatFunc *pff = g_object_new (NCM_TYPE_PRIOR_FLAT_FUNC, 
                                        "mean-func", mean_func,
                                        "x-low", x_low,
                                        "x-upp", x_upp,
                                        "scale", scale,
                                        "variable", variable,
                                        NULL);
  
  return pff;
}

/**
 * ncm_prior_flat_func_ref:
 * @pff: a #NcmPriorFlatFunc
 * 
 * Increases the reference count of @pff atomically.
 * 
 * Returns: (transfer full): @pff.
 */
NcmPriorFlatFunc *
ncm_prior_flat_func_ref (NcmPriorFlatFunc *pff)
{
  return g_object_ref (pff);
}

/**
 * ncm_prior_flat_func_free:
 * @pff: a #NcmPriorFlatFunc
 * 
 * Decreases the reference count of @pff atomically.
 * 
 */
void 
ncm_prior_flat_func_free (NcmPriorFlatFunc *pff)
{
  g_object_unref (pff);
}

/**
 * ncm_prior_flat_func_clear:
 * @pff: a #NcmPriorFlatFunc
 * 
 * Decreases the reference count of *@pff and sets *@pff to NULL.
 * 
 */
void 
ncm_prior_flat_func_clear (NcmPriorFlatFunc **pff)
{
  g_clear_object (pff);
}
