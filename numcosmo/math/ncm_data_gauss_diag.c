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
 * @title: Gaussian Data
 * @short_description: Gaussian abstract class - Diagonal covariance
 *
 * Generic gaussian distribution which uses the covariance matrix
 * to perform the calculations.
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
  PROP_SIZE,
};

G_DEFINE_ABSTRACT_TYPE (NcmDataGaussDiag, ncm_data_gauss_diag, NCM_TYPE_DATA);

static void
ncm_data_gauss_diag_init (NcmDataGaussDiag *gauss)
{
  gauss->np         = 0;
  gauss->y          = NULL;
  gauss->v          = NULL;
  gauss->sigma      = NULL;
  gauss->weight     = NULL;
  gauss->wt         = 0.0;
  gauss->prepared_w = FALSE;
  gauss->wmean      = FALSE;
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
  NcmDataGaussDiag *gauss = NCM_DATA_GAUSS_DIAG (object);
  g_return_if_fail (NCM_IS_DATA_GAUSS_DIAG (object));

  switch (prop_id)
  {
    case PROP_NPOINTS:
      ncm_data_gauss_diag_set_size (gauss, g_value_get_uint (value));
      break;
    case PROP_WMEAN:
      gauss->wmean = g_value_get_boolean (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_data_gauss_diag_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmDataGaussDiag *gauss = NCM_DATA_GAUSS_DIAG (object);
  g_return_if_fail (NCM_IS_DATA_GAUSS_DIAG (object));

  switch (prop_id)
  {
    case PROP_NPOINTS:
      g_value_set_uint (value, ncm_data_gauss_diag_get_size (gauss));
      break;
    case PROP_WMEAN:
      g_value_set_boolean (value, gauss->wmean);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_data_gauss_diag_dispose (GObject *object)
{
  NcmDataGaussDiag *gauss = NCM_DATA_GAUSS_DIAG (object);

  ncm_vector_clear (&gauss->y);
  ncm_vector_clear (&gauss->v);
  ncm_vector_clear (&gauss->sigma);
  ncm_vector_clear (&gauss->weight);

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
static void _ncm_data_gauss_diag_copyto (NcmData *data, NcmData *data_dest);
/* static void _ncm_data_gauss_diag_begin (NcmData *data); */
static void _ncm_data_gauss_diag_resample (NcmData *data, NcmMSet *mset);
static void _ncm_data_gauss_diag_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
static void _ncm_data_gauss_diag_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *v);

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
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_WMEAN,
                                   g_param_spec_boolean ("w-mean",
                                                         NULL,
                                                         "Whether to minimize analytically over the weighted mean",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->get_length       = &_ncm_data_gauss_diag_get_length;
  data_class->get_dof          = &_ncm_data_gauss_diag_get_dof;
  data_class->copyto           = &_ncm_data_gauss_diag_copyto;

  data_class->resample         = &_ncm_data_gauss_diag_resample;
  data_class->m2lnL_val        = &_ncm_data_gauss_diag_m2lnL_val;
  data_class->leastsquares_f   = &_ncm_data_gauss_diag_leastsquares_f;

  gauss_diag_class->mean_func  = NULL;
  gauss_diag_class->sigma_func = NULL;
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
_ncm_data_gauss_diag_copyto (NcmData *data, NcmData *data_dest)
{
  NcmDataGaussDiag *src  = NCM_DATA_GAUSS_DIAG (data);
  NcmDataGaussDiag *dest = NCM_DATA_GAUSS_DIAG (data_dest);
  
  ncm_vector_memcpy (dest->y, src->y);
  ncm_vector_memcpy (dest->sigma, src->sigma);

  g_assert (dest->weight == NULL);
  
  if (src->weight != NULL)
    dest->weight = ncm_vector_dup (src->weight);

  ncm_data_set_init (data_dest);
}

static void
_ncm_data_gauss_prepare_weight (NcmData *data)
{
  NcmDataGaussDiag *gauss = NCM_DATA_GAUSS_DIAG (data);
  gint i;

  if (gauss->weight == NULL)
    gauss->weight = ncm_vector_new_sunk (gauss->np);

  gauss->wt = 0.0;
  for (i = 0; i < gauss->np; i++)
  {
    const gdouble sigma_i = ncm_vector_get (gauss->sigma, i);
    const gdouble w_i = 1.0 / (sigma_i * sigma_i);
    ncm_vector_set (gauss->weight, i, w_i);
    gauss->wt += w_i;
  }

  gauss->prepared_w = TRUE;
}

static void
_ncm_data_gauss_diag_resample (NcmData *data, NcmMSet *mset)
{
  NcmDataGaussDiag *gauss = NCM_DATA_GAUSS_DIAG (data);
  NcmDataGaussDiagClass *gauss_diag_class = NCM_DATA_GAUSS_DIAG_GET_CLASS (gauss);
  gsl_rng *rng = ncm_cfg_rng_get ();
  gint i;

  if (gauss_diag_class->sigma_func != NULL)
    gauss_diag_class->sigma_func (gauss, mset, gauss->sigma);

  gauss_diag_class->mean_func (gauss, mset, gauss->y);
  
  for (i = 0; i < gauss->np; i++)
  {
    const gdouble sigma_i = ncm_vector_get (gauss->sigma, i);
    const gdouble u_i = gsl_ran_gaussian_ziggurat (rng, sigma_i);
    ncm_vector_subfrom (gauss->y, i, u_i);
  }
}

static void
_ncm_data_gauss_diag_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcmDataGaussDiag *gauss = NCM_DATA_GAUSS_DIAG (data);
  NcmDataGaussDiagClass *gauss_diag_class = NCM_DATA_GAUSS_DIAG_GET_CLASS (gauss);
  gboolean sigma_update = FALSE;
  gint i;

  *m2lnL = 0.0;

  if (gauss_diag_class->sigma_func != NULL)
    sigma_update = gauss_diag_class->sigma_func (gauss, mset, gauss->sigma);

  gauss_diag_class->mean_func (gauss, mset, gauss->v);

  if (gauss->wmean)
  {
    gdouble tmp = 0.0;
    if (sigma_update || !gauss->prepared_w)
      _ncm_data_gauss_prepare_weight (data);

    for (i = 0; i < gauss->np; i++)
    {
      const gdouble y_i = ncm_vector_get (gauss->y, i);
      const gdouble yt_i = ncm_vector_get (gauss->v, i);
      const gdouble r_i = (yt_i - y_i);
      const gdouble w_i = ncm_vector_get (gauss->weight, i);
      const gdouble r_i2 = r_i * r_i;
      *m2lnL += r_i2 * w_i;
      tmp += r_i * w_i;
    }

    *m2lnL -= tmp * tmp / gauss->wt;
  }
  else
  {
    for (i = 0; i < gauss->np; i++)
    {
      const gdouble y_i = ncm_vector_get (gauss->y, i);
      const gdouble yt_i = ncm_vector_get (gauss->v, i);
      const gdouble sigma_i = ncm_vector_get (gauss->sigma, i);
      const gdouble r_i = (yt_i - y_i) / sigma_i;
      *m2lnL += r_i * r_i;
    }
  }
}

static void
_ncm_data_gauss_diag_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *v)
{
  NcmDataGaussDiag *gauss = NCM_DATA_GAUSS_DIAG (data);
  NcmDataGaussDiagClass *gauss_diag_class = NCM_DATA_GAUSS_DIAG_GET_CLASS (gauss);
  gboolean sigma_update = FALSE;
  gint i;

  if (gauss_diag_class->sigma_func != NULL)
    sigma_update = gauss_diag_class->sigma_func (gauss, mset, gauss->sigma);

  gauss_diag_class->mean_func (gauss, mset, v);

  if (gauss->wmean)
  {
    gdouble wmean;
    if (sigma_update || !gauss->prepared_w)
      _ncm_data_gauss_prepare_weight (data);

    ncm_vector_sub (v, gauss->y);

    wmean = gsl_stats_wmean (ncm_vector_gsl (gauss->weight)->data,
                             ncm_vector_gsl (gauss->weight)->stride,
                             ncm_vector_gsl (v)->data,
                             ncm_vector_gsl (v)->stride,
                             ncm_vector_gsl (v)->size);

    for (i = 0; i < gauss->np; i++)
    {
      const gdouble sigma_i = ncm_vector_get (gauss->sigma, i);
      const gdouble v_i = ncm_vector_get (v, i);
      ncm_vector_set (v, i, (v_i - wmean) / sigma_i);
    }
  }
  else
  {
    for (i = 0; i < gauss->np; i++)
    {
      const gdouble y_i = ncm_vector_get (gauss->y, i);
      const gdouble yt_i = ncm_vector_get (v, i);
      const gdouble sigma_i = ncm_vector_get (gauss->sigma, i);
      const gdouble r_i = (yt_i - y_i) / sigma_i;
      ncm_vector_set (v, i, r_i);
    }
  }
}

/**
 * ncm_data_gauss_diag_set_size:
 * @diag: a #NcmDataGauss
 * @np: FIXME
 *
 * FIXME
 * 
 */
void 
ncm_data_gauss_diag_set_size (NcmDataGaussDiag *diag, guint np)
{
  if ((np == 0) || (np != diag->np))
  {
    diag->np = 0;
    ncm_vector_clear (&diag->y);
    ncm_vector_clear (&diag->v);
    ncm_vector_clear (&diag->sigma);
    ncm_vector_clear (&diag->weight);
  }
  if ((np != 0) && (np != diag->np))
  {
    diag->np    = np;
    diag->y     = ncm_vector_new_sunk (diag->np);
    diag->v     = ncm_vector_new_sunk (diag->np);
    diag->sigma = ncm_vector_new_sunk (diag->np);    
  }
}

/**
 * ncm_data_gauss_diag_get_size:
 * @diag: a #NcmDataGauss
 *
 * FIXME
 * 
 * Returns: FIXME
 */
guint 
ncm_data_gauss_diag_get_size (NcmDataGaussDiag *diag)
{
  return diag->np;
}
