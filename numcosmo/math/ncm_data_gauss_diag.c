/***************************************************************************
 *            ncm_data_gauss_diag.c
 *
 *  Fri Mar 19 14:57:35 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_data_gauss_diag
 * @title: Gaussian Data - DiagCov
 * @short_description: Gaussian data object, diagonal covariance
 *
 * Gaussian distribution which uses a diagonal covariance matrix as input.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_data_gauss_diag.h"
#include "math/ncm_cfg.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>

enum
{
  PROP_0,
  PROP_NPOINTS,
  PROP_WMEAN,
  PROP_MEAN,
  PROP_SIGMA,
  PROP_SIZE,
};

G_DEFINE_ABSTRACT_TYPE (NcmDataGaussDiag, ncm_data_gauss_diag, NCM_TYPE_DATA);

static void
ncm_data_gauss_diag_init (NcmDataGaussDiag *diag)
{
  diag->np         = 0;
  diag->y          = NULL;
  diag->v          = NULL;
  diag->sigma      = NULL;
  diag->weight     = NULL;
  diag->wt         = 0.0;
  diag->prepared_w = FALSE;
  diag->wmean      = FALSE;
}

static void
_ncm_data_gauss_diag_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_data_gauss_diag_parent_class)->constructed (object);

}

static void
ncm_data_gauss_diag_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (object);
  g_return_if_fail (NCM_IS_DATA_GAUSS_DIAG (object));

  switch (prop_id)
  {
    case PROP_NPOINTS:
      ncm_data_gauss_diag_set_size (diag, g_value_get_uint (value));
      break;
    case PROP_WMEAN:
      diag->wmean = g_value_get_boolean (value);
      break;
    case PROP_MEAN:
      ncm_vector_set_from_variant (diag->y, g_value_get_variant (value));
      break;
    case PROP_SIGMA:
      ncm_vector_set_from_variant (diag->sigma, g_value_get_variant (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_data_gauss_diag_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (object);
  g_return_if_fail (NCM_IS_DATA_GAUSS_DIAG (object));

  switch (prop_id)
  {
    case PROP_NPOINTS:
      g_value_set_uint (value, ncm_data_gauss_diag_get_size (diag));
      break;
    case PROP_WMEAN:
      g_value_set_boolean (value, diag->wmean);
      break;
    case PROP_MEAN:
      g_value_take_variant (value, ncm_vector_get_variant (diag->y));
      break;
    case PROP_SIGMA:
      g_value_take_variant (value, ncm_vector_get_variant (diag->sigma));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_data_gauss_diag_dispose (GObject *object)
{
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (object);

  ncm_vector_clear (&diag->y);
  ncm_vector_clear (&diag->v);
  ncm_vector_clear (&diag->sigma);
  ncm_vector_clear (&diag->weight);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_gauss_diag_parent_class)->dispose (object);
}

static void
ncm_data_gauss_diag_finalize (GObject *object)
{
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_gauss_diag_parent_class)->finalize (object);
}

static guint _ncm_data_gauss_diag_get_length (NcmData *data);
static guint _ncm_data_gauss_diag_get_dof (NcmData *data);
/* static void _ncm_data_gauss_diag_begin (NcmData *data); */
static void _ncm_data_gauss_diag_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng);
static void _ncm_data_gauss_diag_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
static void _ncm_data_gauss_diag_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *v);
static void _ncm_data_gauss_diag_set_size (NcmDataGaussDiag *diag, guint np);
static guint _ncm_data_gauss_diag_get_size (NcmDataGaussDiag *diag);

static void
ncm_data_gauss_diag_class_init (NcmDataGaussDiagClass *klass)
{
  GObjectClass* object_class     = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class       = NCM_DATA_CLASS (klass);
  NcmDataGaussDiagClass *gauss_diag_class = NCM_DATA_GAUSS_DIAG_CLASS (klass);

  object_class->constructed  = &_ncm_data_gauss_diag_constructed;
  object_class->set_property = &ncm_data_gauss_diag_set_property;
  object_class->get_property = &ncm_data_gauss_diag_get_property;
  object_class->dispose      = &ncm_data_gauss_diag_dispose;
  object_class->finalize     = &ncm_data_gauss_diag_finalize;

  g_object_class_install_property (object_class,
                                   PROP_NPOINTS,
                                   g_param_spec_uint ("n-points",
                                                      NULL,
                                                      "Data sample size",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_WMEAN,
                                   g_param_spec_boolean ("w-mean",
                                                         NULL,
                                                         "Whether to minimize analytically over the weighted mean",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_MEAN,
                                   g_param_spec_variant ("mean",
                                                         NULL,
                                                         "Data mean",
                                                         G_VARIANT_TYPE ("ad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  g_object_class_install_property (object_class,
                                   PROP_SIGMA,
                                   g_param_spec_variant ("sigma",
                                                         NULL,
                                                         "Data standard deviation",
                                                         G_VARIANT_TYPE ("ad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->bootstrap        = TRUE;
  data_class->get_length       = &_ncm_data_gauss_diag_get_length;
  data_class->get_dof          = &_ncm_data_gauss_diag_get_dof;
  data_class->resample         = &_ncm_data_gauss_diag_resample;
  data_class->m2lnL_val        = &_ncm_data_gauss_diag_m2lnL_val;
  data_class->leastsquares_f   = &_ncm_data_gauss_diag_leastsquares_f;

  gauss_diag_class->mean_func  = NULL;
  gauss_diag_class->sigma_func = NULL;
  gauss_diag_class->set_size   = &_ncm_data_gauss_diag_set_size;
  gauss_diag_class->get_size   = &_ncm_data_gauss_diag_get_size;
}

static guint 
_ncm_data_gauss_diag_get_length (NcmData *data) 
{ 
  return NCM_DATA_GAUSS_DIAG (data)->np; 
}

static guint 
_ncm_data_gauss_diag_get_dof (NcmData *data) 
{ 
  guint dof = NCM_DATA_GAUSS_DIAG (data)->np;
  if (NCM_DATA_GAUSS_DIAG (data)->wmean && dof > 0)
    dof--;
  return dof;
}

static void
_ncm_data_gauss_prepare_weight (NcmData *data)
{
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (data);
  guint i;

  if (diag->weight == NULL)
    diag->weight = ncm_vector_new (diag->np);

  diag->wt = 0.0;
  for (i = 0; i < diag->np; i++)
  {
    const gdouble sigma_i = ncm_vector_get (diag->sigma, i);
    const gdouble w_i = 1.0 / (sigma_i * sigma_i);
    ncm_vector_set (diag->weight, i, w_i);
    diag->wt += w_i;
  }

  diag->prepared_w = TRUE;
}

#include "numcosmo/nc_hicosmo.h"
#include "numcosmo/nc_data_dist_mu.h"

static void
_ncm_data_gauss_diag_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng)
{
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (data);
  NcmDataGaussDiagClass *gauss_diag_class = NCM_DATA_GAUSS_DIAG_GET_CLASS (diag);
  guint i;

  if (gauss_diag_class->sigma_func != NULL)
    gauss_diag_class->sigma_func (diag, mset, diag->sigma);

  gauss_diag_class->mean_func (diag, mset, diag->y); 
  
  ncm_rng_lock (rng);
  for (i = 0; i < diag->np; i++)
  {
    const gdouble sigma_i = ncm_vector_get (diag->sigma, i);
    const gdouble u_i = gsl_ran_gaussian_ziggurat (rng->r, sigma_i);
    ncm_vector_subfrom (diag->y, i, u_i);
  }
  ncm_rng_unlock (rng);
}

static void
_ncm_data_gauss_diag_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (data);
  NcmDataGaussDiagClass *gauss_diag_class = NCM_DATA_GAUSS_DIAG_GET_CLASS (diag);
  gboolean sigma_update = FALSE;
  guint i;

  *m2lnL = 0.0;

  if (gauss_diag_class->sigma_func != NULL)
    sigma_update = gauss_diag_class->sigma_func (diag, mset, diag->sigma);

  gauss_diag_class->mean_func (diag, mset, diag->v);

  if (diag->wmean)
  {
    if (sigma_update || !diag->prepared_w)
      _ncm_data_gauss_prepare_weight (data);

    if (!ncm_data_bootstrap_enabled (data))
    {
      gdouble tmp = 0.0;

      for (i = 0; i < diag->np; i++)
      {
        const gdouble y_i = ncm_vector_get (diag->y, i);
        const gdouble yt_i = ncm_vector_get (diag->v, i);
        const gdouble r_i = (yt_i - y_i);
        const gdouble w_i = ncm_vector_get (diag->weight, i);
        const gdouble r_i2 = r_i * r_i;
        *m2lnL += r_i2 * w_i;
        tmp += r_i * w_i;
      }
      *m2lnL -= tmp * tmp / diag->wt;
    }
    else
    {
      const guint bsize = ncm_bootstrap_get_bsize (data->bstrap);
      gdouble tmp = 0.0;
      gdouble wt = 0.0;
      
      for (i = 0; i < bsize; i++)
      {
        guint k = ncm_bootstrap_get (data->bstrap, i);
        const gdouble y_k = ncm_vector_get (diag->y, k);
        const gdouble yt_k = ncm_vector_get (diag->v, k);
        const gdouble r_k = (yt_k - y_k);
        const gdouble w_k = ncm_vector_get (diag->weight, k);
        const gdouble r_k2 = r_k * r_k;
        *m2lnL += r_k2 * w_k;
        tmp += r_k * w_k;
        wt += w_k;
      }
      *m2lnL -= tmp * tmp / wt;
    }
  }
  else
  {
    if (!ncm_data_bootstrap_enabled (data))
    {
      for (i = 0; i < diag->np; i++)
      {
        const gdouble y_i = ncm_vector_get (diag->y, i);
        const gdouble yt_i = ncm_vector_get (diag->v, i);
        const gdouble sigma_i = ncm_vector_get (diag->sigma, i);
        const gdouble r_i = (yt_i - y_i) / sigma_i;
        *m2lnL += r_i * r_i;
      }
    }
    else
    {
      const guint bsize = ncm_bootstrap_get_bsize (data->bstrap);
      for (i = 0; i < bsize; i++)
      {
        guint k = ncm_bootstrap_get (data->bstrap, i);
        const gdouble y_k = ncm_vector_get (diag->y, k);
        const gdouble yt_k = ncm_vector_get (diag->v, k);
        const gdouble sigma_k = ncm_vector_get (diag->sigma, k);
        const gdouble r_k = (yt_k - y_k) / sigma_k;
        *m2lnL += r_k * r_k;
      }
    }
  }
}

static void
_ncm_data_gauss_diag_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *v)
{
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (data);
  NcmDataGaussDiagClass *gauss_diag_class = NCM_DATA_GAUSS_DIAG_GET_CLASS (diag);
  gboolean sigma_update = FALSE;
  guint i;

  if (ncm_data_bootstrap_enabled (data))
    g_error ("NcmDataGaussDiag: does not support bootstrap with least squares");
  
  if (gauss_diag_class->sigma_func != NULL)
    sigma_update = gauss_diag_class->sigma_func (diag, mset, diag->sigma);

  gauss_diag_class->mean_func (diag, mset, v);

  if (diag->wmean)
  {
    gdouble wmean;
    if (sigma_update || !diag->prepared_w)
      _ncm_data_gauss_prepare_weight (data);

    ncm_vector_sub (v, diag->y);

    wmean = gsl_stats_wmean (ncm_vector_gsl (diag->weight)->data,
                             ncm_vector_gsl (diag->weight)->stride,
                             ncm_vector_gsl (v)->data,
                             ncm_vector_gsl (v)->stride,
                             ncm_vector_gsl (v)->size);

    for (i = 0; i < diag->np; i++)
    {
      const gdouble sigma_i = ncm_vector_get (diag->sigma, i);
      const gdouble v_i = ncm_vector_get (v, i);
      ncm_vector_set (v, i, (v_i - wmean) / sigma_i);
    }
  }
  else
  {
    for (i = 0; i < diag->np; i++)
    {
      const gdouble y_i = ncm_vector_get (diag->y, i);
      const gdouble yt_i = ncm_vector_get (v, i);
      const gdouble sigma_i = ncm_vector_get (diag->sigma, i);
      const gdouble r_i = (yt_i - y_i) / sigma_i;
      ncm_vector_set (v, i, r_i);
    }
  }
}

static void 
_ncm_data_gauss_diag_set_size (NcmDataGaussDiag *diag, guint np)
{
  NcmData *data = NCM_DATA (diag);
  if ((np == 0) || (np != diag->np))
  {
    diag->np = 0;
    ncm_vector_clear (&diag->y);
    ncm_vector_clear (&diag->v);
    ncm_vector_clear (&diag->sigma);
    ncm_vector_clear (&diag->weight);
    data->init = FALSE;
  }
  if ((np != 0) && (np != diag->np))
  {
    diag->np    = np;
    diag->y     = ncm_vector_new (diag->np);
    diag->v     = ncm_vector_new (diag->np);
    diag->sigma = ncm_vector_new (diag->np);
    if (ncm_data_bootstrap_enabled (data))
    {
      ncm_bootstrap_set_fsize (data->bstrap, np);
      ncm_bootstrap_set_bsize (data->bstrap, np);
    }
    data->init = FALSE;
  }
}

static guint 
_ncm_data_gauss_diag_get_size (NcmDataGaussDiag *diag)
{
  return diag->np;
}

/**
 * ncm_data_gauss_diag_set_size: (virtual set_size)
 * @diag: a #NcmDataGaussDiag
 * @np: data size.
 *
 * Sets the data size to @np.
 * 
 */
void 
ncm_data_gauss_diag_set_size (NcmDataGaussDiag *diag, guint np)
{
  NCM_DATA_GAUSS_DIAG_GET_CLASS (diag)->set_size (diag, np);
}

/**
 * ncm_data_gauss_diag_get_size: (virtual get_size)
 * @diag: a #NcmDataGaussDiag
 *
 * Gets the data size.
 * 
 * Returns: Data size.
 * 
 */
guint 
ncm_data_gauss_diag_get_size (NcmDataGaussDiag *diag)
{
  return NCM_DATA_GAUSS_DIAG_GET_CLASS (diag)->get_size (diag);
}
