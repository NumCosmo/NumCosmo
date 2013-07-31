/***************************************************************************
 *            ncm_data_dist1d.c
 *
 *  Thu Apr 15 11:16:11 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
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
 * SECTION:ncm_data_dist1d
 * @title: One Variable Distribution Data
 * @short_description: Object representing a one variable distribution data
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_data_dist1d.h"
#include "math/ncm_rng.h"

enum
{
  PROP_0,
  PROP_NPOINTS,
  PROP_SIZE,
};

G_DEFINE_ABSTRACT_TYPE (NcmDataDist1d, ncm_data_dist1d, NCM_TYPE_DATA);

static void
ncm_data_dist1d_init (NcmDataDist1d *dist1d)
{
  dist1d->np = 0;
  dist1d->x  = NULL;
}

static void
_ncm_data_dist1d_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_data_dist1d_parent_class)->constructed (object);
}

static void
_ncm_data_dist1d_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmDataDist1d *dist1d = NCM_DATA_DIST1D (object);
  g_return_if_fail (NCM_IS_DATA_DIST1D (object));

  switch (prop_id)
  {
    case PROP_NPOINTS:
      ncm_data_dist1d_set_size (dist1d, g_value_get_uint (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_data_dist1d_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmDataDist1d *dist1d = NCM_DATA_DIST1D (object);

  g_return_if_fail (NCM_IS_DATA_DIST1D (object));

  switch (prop_id)
  {
    case PROP_NPOINTS:
      g_value_set_uint (value, ncm_data_dist1d_get_size (dist1d));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_data_dist1d_dispose (GObject *object)
{
  NcmDataDist1d *dist1d = NCM_DATA_DIST1D (object);

  ncm_vector_clear (&dist1d->x);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_dist1d_parent_class)->dispose (object);
}

static void
ncm_data_dist1d_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_dist1d_parent_class)->finalize (object);
}

static guint _ncm_data_dist1d_get_length (NcmData *data);
static void _ncm_data_dist1d_copyto (NcmData *data, NcmData *data_dest);
static void _ncm_data_dist1d_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
static void _ncm_data_dist1d_resample (NcmData *data, NcmMSet *mset);

static void
ncm_data_dist1d_class_init (NcmDataDist1dClass *klass)
{
  GObjectClass *object_class       = G_OBJECT_CLASS (klass);
  NcmDataDist1dClass *dist1d_class = NCM_DATA_DIST1D_CLASS (klass);
  NcmDataClass *data_class         = NCM_DATA_CLASS (klass);

  object_class->constructed  = &_ncm_data_dist1d_constructed;
  object_class->set_property = &_ncm_data_dist1d_set_property;
  object_class->get_property = &_ncm_data_dist1d_get_property;

  object_class->dispose      = &ncm_data_dist1d_dispose;
  object_class->finalize     = &ncm_data_dist1d_finalize;

  g_object_class_install_property (object_class,
                                   PROP_NPOINTS,
                                   g_param_spec_uint ("n-points",
                                                      NULL,
                                                      "Data sample size",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->get_length = &_ncm_data_dist1d_get_length;
  data_class->copyto     = &_ncm_data_dist1d_copyto;
  data_class->begin      = NULL;

  data_class->resample   = &_ncm_data_dist1d_resample;
  data_class->m2lnL_val  = &_ncm_data_dist1d_m2lnL_val;

  dist1d_class->dist    = NULL;
  dist1d_class->inv_pdf = NULL;

}

static guint 
_ncm_data_dist1d_get_length (NcmData *data)
{ 
  NcmDataDist1d *dist1d = NCM_DATA_DIST1D (data);
  return dist1d->np; 
}

static void
_ncm_data_dist1d_copyto (NcmData *data, NcmData *data_dest)
{
  /* Chain up : start */
  NCM_DATA_CLASS (ncm_data_dist1d_parent_class)->copyto (data, data_dest);
  {
    NcmDataDist1d *dist1d = NCM_DATA_DIST1D (data);
    NcmDataDist1d *dist1_dest = NCM_DATA_DIST1D (data_dest);
    ncm_vector_memcpy (dist1_dest->x, dist1d->x);
  }
}

static void
_ncm_data_dist1d_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcmDataDist1d *dist1d = NCM_DATA_DIST1D (data);
  NcmDataDist1dClass *dist1d_class = NCM_DATA_DIST1D_GET_CLASS (data);
  gint i;
  
  *m2lnL = 0.0;
  if (!data->bootstrap)
  {
    for (i = 0; i < dist1d->np; i++)
    {
      const gdouble x_i = ncm_vector_get (dist1d->x, i);
      *m2lnL += ncm_mset_func_eval1 (dist1d_class->dist, mset, x_i);
    }
  }
  else
  {
    const guint bsize = ncm_bootstrap_get_bsize (data->bstrap);
    for (i = 0; i < bsize; i++)
    {
      guint k = ncm_bootstrap_get (data->bstrap, i);
      const gdouble x_i = ncm_vector_get (dist1d->x, k);
      *m2lnL += ncm_mset_func_eval1 (dist1d_class->dist, mset, x_i);
    }    
  }
  return;
}

static void
_ncm_data_dist1d_resample (NcmData *data, NcmMSet *mset)
{
  NcmDataDist1d *dist1d = NCM_DATA_DIST1D (data);
  NcmDataDist1dClass *dist1d_class = NCM_DATA_DIST1D_GET_CLASS (data);
  NcmRNG *rng = ncm_rng_pool_get (NCM_DATA_RESAMPLE_RNG_NAME);
  gint i;

  if (dist1d_class->inv_pdf == NULL)
    g_error ("_ncm_data_dist1d_resample: This object do not implement the inverse of the pdf.");

  ncm_rng_lock (rng);
  for (i = 0; i < dist1d->np; i++)
  {
    const gdouble u_i = gsl_rng_uniform (rng->r);
    const gdouble x_i = ncm_mset_func_eval1 (dist1d_class->inv_pdf, mset, u_i);
    ncm_vector_set (dist1d->x, i, x_i);
  }
  ncm_rng_unlock (rng);
  ncm_rng_free (rng);
}

/**
 * ncm_data_dist1d_set_size:
 * @dist1d: a #NcmDataDist1d
 * @np: FIXME
 *
 * FIXME
 * 
 */
void 
ncm_data_dist1d_set_size (NcmDataDist1d *dist1d, guint np)
{
  NcmData *data = NCM_DATA (dist1d);
  if ((np == 0) || (np != dist1d->np))
  {
    dist1d->np = 0;
    ncm_vector_clear (&dist1d->x);
    data->init = FALSE;
  }
  if ((np != 0) && (np != dist1d->np))
  {
    dist1d->np = np;
    dist1d->x  = ncm_vector_new (dist1d->np);
    if (data->bootstrap)
      ncm_bootstrap_set_fsize (data->bstrap, np);
    data->init = FALSE;
  }
}

/**
 * ncm_data_dist1d_get_size:
 * @dist1d: a #NcmDataDist1d
 *
 * FIXME
 * 
 * Returns: FIXME
 */
guint 
ncm_data_dist1d_get_size (NcmDataDist1d *dist1d)
{
  return dist1d->np;
}

