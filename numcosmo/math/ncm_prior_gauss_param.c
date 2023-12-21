/***************************************************************************
 *            ncm_prior_gauss_param.c
 *
 *  Wed August 03 10:55:21 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_prior_gauss_param.c
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
 * SECTION:ncm_prior_gauss_param
 * @title: NcmPriorGaussParam
 * @short_description: a gaussian prior on a parameter
 *
 * This object is a subclass of #NcmPriorGauss, precisely designed for a Gaussian
 * prior on a parameter. The prior relies on a parameter ID and
 * user-specified mean and standard deviation parameters.
 *
 * Users have flexibility in specifying the parameter in various ways:
 * - Using the pair NcmModelID and the parameter pid.
 * - Providing a single NcmMSetPIndex.
 * - Supplying a string consisting of a parameter full name "model:parameter".
 *
 * These options provide versatile and convenient methods for parameter identification.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_prior_gauss_param.h"

enum
{
  PROP_0,
  PROP_MID,
  PROP_PID
};

struct _NcmPriorGaussParam
{
  /*< private >*/
  NcmPriorGauss parent_instance;
  NcmModelID mid;
  guint pid;
};

G_DEFINE_TYPE (NcmPriorGaussParam, ncm_prior_gauss_param, NCM_TYPE_PRIOR_GAUSS)

static void
ncm_prior_gauss_param_init (NcmPriorGaussParam *pgp)
{
}

static void
_ncm_prior_gauss_param_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmPriorGaussParam *pgp = NCM_PRIOR_GAUSS_PARAM (object);

  g_return_if_fail (NCM_IS_PRIOR_GAUSS_PARAM (object));

  switch (prop_id)
  {
    case PROP_MID:
      pgp->mid = g_value_get_int (value);
      break;
    case PROP_PID:
      pgp->pid = g_value_get_uint (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_prior_gauss_param_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmPriorGaussParam *pgp = NCM_PRIOR_GAUSS_PARAM (object);

  g_return_if_fail (NCM_IS_PRIOR_GAUSS_PARAM (object));

  switch (prop_id)
  {
    case PROP_MID:
      g_value_set_int (value, pgp->mid);
      break;
    case PROP_PID:
      g_value_set_uint (value, pgp->pid);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_prior_gauss_param_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_prior_gauss_param_parent_class)->finalize (object);
}

static gdouble _ncm_prior_gauss_param_mean (NcmPriorGauss *pg, NcmMSet *mset);

static void
ncm_prior_gauss_param_class_init (NcmPriorGaussParamClass *klass)
{
  GObjectClass *object_class   = G_OBJECT_CLASS (klass);
  NcmPriorGaussClass *pg_class = NCM_PRIOR_GAUSS_CLASS (klass);

  object_class->set_property = &_ncm_prior_gauss_param_set_property;
  object_class->get_property = &_ncm_prior_gauss_param_get_property;
  object_class->finalize     = &_ncm_prior_gauss_param_finalize;

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

  pg_class->mean = &_ncm_prior_gauss_param_mean;
}

static gdouble
_ncm_prior_gauss_param_mean (NcmPriorGauss *pg, NcmMSet *mset)
{
  NcmPriorGaussParam *pgp = NCM_PRIOR_GAUSS_PARAM (pg);

  return ncm_mset_param_get (mset, pgp->mid, pgp->pid);
}

/**
 * ncm_prior_gauss_param_new:
 * @mid: model id
 * @pid: parameter id
 * @mu: mean
 * @sigma: standard deviation
 *
 * Creates a new Gaussian prior for parameter @pid of model @mid.
 *
 * Returns: (transfer full): @pgp.
 */
NcmPriorGaussParam *
ncm_prior_gauss_param_new (NcmModelID mid, guint pid, gdouble mu, gdouble sigma)
{
  NcmPriorGaussParam *pgp = g_object_new (NCM_TYPE_PRIOR_GAUSS_PARAM,
                                          "mu", mu,
                                          "sigma", sigma,
                                          "mid", mid,
                                          "pid", pid,
                                          NULL);

  return pgp;
}

/**
 * ncm_prior_gauss_param_new_pindex:
 * @pi: a #NcmMSetPIndex
 * @mu: mean
 * @sigma: standard deviation
 *
 * Creates a new Gaussian prior for parameter @pid of model @mid.
 *
 * Returns: (transfer full): @pgp.
 */
NcmPriorGaussParam *
ncm_prior_gauss_param_new_pindex (const NcmMSetPIndex *pi, gdouble mu, gdouble sigma)
{
  return ncm_prior_gauss_param_new (pi->mid, pi->pid, mu, sigma);
}

/**
 * ncm_prior_gauss_param_new_name:
 * @mset: a #NcmMSet
 * @name: parameter name
 * @mu: mean
 * @sigma: standard deviation
 *
 * Creates a new Gaussian prior for parameter named @name in @mset.
 *
 * Returns: (transfer full): @pgp.
 */
NcmPriorGaussParam *
ncm_prior_gauss_param_new_name (NcmMSet *mset, const gchar *name, gdouble mu, gdouble sigma)
{
  const NcmMSetPIndex *pi = NULL;

  if ((pi = ncm_mset_param_get_by_full_name (mset, name)) != NULL)
  {
    NcmPriorGaussParam *pgp = ncm_prior_gauss_param_new_pindex (pi, mu, sigma);

    ncm_mset_pindex_free ((NcmMSetPIndex *) pi);

    return pgp;
  }
  else if ((pi = ncm_mset_fparam_get_pi_by_name (mset, name)) != NULL)
  {
    return ncm_prior_gauss_param_new_pindex (pi, mu, sigma);
  }
  else
  {
    g_error ("ncm_prior_gauss_param_new_name: cannot find parameter named `%s' in the mset.", name);

    return NULL;
  }
}

/**
 * ncm_prior_gauss_param_ref:
 * @pgp: a #NcmPriorGaussParam
 *
 * Increases the reference count of @pgp atomically.
 *
 * Returns: (transfer full): @pgp.
 */
NcmPriorGaussParam *
ncm_prior_gauss_param_ref (NcmPriorGaussParam *pgp)
{
  return g_object_ref (pgp);
}

/**
 * ncm_prior_gauss_param_free:
 * @pgp: a #NcmPriorGaussParam
 *
 * Decreases the reference count of @pgp atomically.
 *
 */
void
ncm_prior_gauss_param_free (NcmPriorGaussParam *pgp)
{
  g_object_unref (pgp);
}

/**
 * ncm_prior_gauss_param_clear:
 * @pgp: a #NcmPriorGaussParam
 *
 * Decreases the reference count of *@pgp and sets *@pgp to NULL.
 *
 */
void
ncm_prior_gauss_param_clear (NcmPriorGaussParam **pgp)
{
  g_clear_object (pgp);
}

