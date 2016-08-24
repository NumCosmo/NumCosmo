/***************************************************************************
 *            ncm_prior_flat_param.c
 *
 *  Wed August 03 10:55:21 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_prior_flat_param.c
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
 * SECTION:ncm_prior_flat_param
 * @title: NcmPriorFlatParam
 * @short_description: a flat prior on a parameter
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_prior_flat_param.h"

enum
{
  PROP_0,
  PROP_MID,
  PROP_PID
};

G_DEFINE_TYPE (NcmPriorFlatParam, ncm_prior_flat_param, NCM_TYPE_PRIOR_FLAT);

static void
ncm_prior_flat_param_init (NcmPriorFlatParam *pfp)
{
}

static void
_ncm_prior_flat_param_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmPriorFlatParam *pfp = NCM_PRIOR_FLAT_PARAM (object);
  g_return_if_fail (NCM_IS_PRIOR_FLAT_PARAM (object));

  switch (prop_id)
  {
    case PROP_MID:
      pfp->mid = g_value_get_int (value);
      break;
    case PROP_PID:
      pfp->pid = g_value_get_uint (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_prior_flat_param_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmPriorFlatParam *pfp = NCM_PRIOR_FLAT_PARAM (object);
  g_return_if_fail (NCM_IS_PRIOR_FLAT_PARAM (object));

  switch (prop_id)
  {
    case PROP_MID:
      g_value_set_int (value, pfp->mid);
      break;
    case PROP_PID:
      g_value_set_uint (value, pfp->pid);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_prior_flat_param_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_prior_flat_param_parent_class)->finalize (object);
}

static gdouble _ncm_prior_flat_param_mean (NcmPriorFlat *pf, NcmMSet *mset);

static void
ncm_prior_flat_param_class_init (NcmPriorFlatParamClass *klass)
{
  GObjectClass* object_class   = G_OBJECT_CLASS (klass);
  NcmPriorFlatClass *pf_class = NCM_PRIOR_FLAT_CLASS (klass);

  object_class->set_property = &_ncm_prior_flat_param_set_property;
  object_class->get_property = &_ncm_prior_flat_param_get_property;
  object_class->finalize     = &_ncm_prior_flat_param_finalize;

  g_object_class_install_property (object_class,
                                   PROP_MID,
                                   g_param_spec_int ("mid",
                                                     NULL,
                                                     "model id",
                                                     0, G_MAXINT, 0,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_PID,
                                   g_param_spec_uint ("pid",
                                                     NULL,
                                                     "parameter id",
                                                     0, G_MAXUINT32, 0,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  pf_class->mean = &_ncm_prior_flat_param_mean;
}

static gdouble 
_ncm_prior_flat_param_mean (NcmPriorFlat *pf, NcmMSet *mset)
{
  NcmPriorFlatParam *pfp = NCM_PRIOR_FLAT_PARAM (pf);
  return ncm_mset_param_get (mset, pfp->mid, pfp->pid);
}

/**
 * ncm_prior_flat_param_new:
 * @mid: model id
 * @pid: parameter id
 * @x_low: FIXME
 * @x_upp: FIXME
 * @scale: FIXME
 * 
 * Creates a new Flat prior for parameter @pid of model @mid.
 * 
 * Returns: (transfer full): @pfp.
 */
NcmPriorFlatParam *
ncm_prior_flat_param_new (NcmModelID mid, guint pid, gdouble x_low, gdouble x_upp, gdouble scale)
{
  NcmPriorFlatParam *pfp = g_object_new (NCM_TYPE_PRIOR_FLAT_PARAM, 
                                         "x-low", x_low,
                                         "x-upp", x_upp,
                                         "scale", scale,
                                         "mid", mid,
                                         "pid", pid,
                                         NULL);
  
  return pfp;
}

/**
 * ncm_prior_flat_param_new_pindex:
 * @pi: a #NcmMSetPIndex
 * @x_low: FIXME
 * @x_upp: FIXME
 * @scale: FIXME
 * 
 * Creates a new Flat prior for parameter @pid of model @mid.
 * 
 * Returns: (transfer full): @pfp.
 */
NcmPriorFlatParam *
ncm_prior_flat_param_new_pindex (const NcmMSetPIndex *pi, gdouble x_low, gdouble x_upp, gdouble scale)
{
  return ncm_prior_flat_param_new (pi->mid, pi->pid, x_low, x_upp, scale);
}

/**
 * ncm_prior_flat_param_new_name:
 * @mset: a #NcmMSet
 * @name: parameter name
 * @x_low: FIXME
 * @x_upp: FIXME
 * @scale: FIXME
 * 
 * Creates a new Flat prior for parameter named @name in @mset.
 * 
 * Returns: (transfer full): @pfp.
 */
NcmPriorFlatParam *
ncm_prior_flat_param_new_name (NcmMSet *mset, const gchar *name, gdouble x_low, gdouble x_upp, gdouble scale)
{
  const NcmMSetPIndex *pi = NULL;
  if ((pi = ncm_mset_param_get_by_full_name (mset, name)) != NULL)
  {
    NcmPriorFlatParam *pfp = ncm_prior_flat_param_new_pindex (pi, x_low, x_upp, scale);

    ncm_mset_pindex_free ((NcmMSetPIndex *)pi);

    return pfp;
  }
  else if ((pi = ncm_mset_fparam_get_pi_by_name (mset, name)) != NULL)
  {
    return ncm_prior_flat_param_new_pindex (pi, x_low, x_upp, scale);
  }
  else
  {
    g_error ("ncm_prior_flat_param_new_name: cannot find parameter named `%s' in the mset.", name);
    return NULL;
  }
}

/**
 * ncm_prior_flat_param_ref:
 * @pfp: a #NcmPriorFlatParam
 * 
 * Increases the reference count of @pfp atomically.
 * 
 * Returns: (transfer full): @pfp.
 */
NcmPriorFlatParam *
ncm_prior_flat_param_ref (NcmPriorFlatParam *pfp)
{
  return g_object_ref (pfp);
}

/**
 * ncm_prior_flat_param_free:
 * @pfp: a #NcmPriorFlatParam
 * 
 * Decreases the reference count of @pfp atomically.
 * 
 */
void 
ncm_prior_flat_param_free (NcmPriorFlatParam *pfp)
{
  g_object_unref (pfp);
}

/**
 * ncm_prior_flat_param_clear:
 * @pfp: a #NcmPriorFlatParam
 * 
 * Decreases the reference count of *@pfp and sets *@pfp to NULL.
 * 
 */
void 
ncm_prior_flat_param_clear (NcmPriorFlatParam **pfp)
{
  g_clear_object (pfp);
}
