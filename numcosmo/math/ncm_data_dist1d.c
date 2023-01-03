/***************************************************************************
 *            ncm_data_dist1d.c
 *
 *  Thu Apr 15 11:16:11 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * @title: NcmDataDist1d
 * @short_description: Abstract class for one variable distribution data.
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
  PROP_VECTOR,
  PROP_SIZE,
};

typedef struct _NcmDataDist1dPrivate
{
  guint np;
  NcmVector *x;
} NcmDataDist1dPrivate;

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmDataDist1d, ncm_data_dist1d, NCM_TYPE_DATA);

static void
ncm_data_dist1d_init (NcmDataDist1d *dist1d)
{
  NcmDataDist1dPrivate * const self = ncm_data_dist1d_get_instance_private (dist1d);

  self->np = 0;
  self->x  = NULL;
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
  NcmDataDist1d *dist1d             = NCM_DATA_DIST1D (object);
  NcmDataDist1dPrivate * const self = ncm_data_dist1d_get_instance_private (dist1d);

  g_return_if_fail (NCM_IS_DATA_DIST1D (object));

  switch (prop_id)
  {
    case PROP_NPOINTS:
      ncm_data_dist1d_set_size (dist1d, g_value_get_uint (value));
      break;
    case PROP_VECTOR:
      ncm_vector_substitute (&self->x, g_value_get_object (value), TRUE);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_data_dist1d_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmDataDist1d *dist1d             = NCM_DATA_DIST1D (object);
  NcmDataDist1dPrivate * const self = ncm_data_dist1d_get_instance_private (dist1d);

  g_return_if_fail (NCM_IS_DATA_DIST1D (object));

  switch (prop_id)
  {
    case PROP_NPOINTS:
      g_value_set_uint (value, ncm_data_dist1d_get_size (dist1d));
      break;
    case PROP_VECTOR:
      g_value_set_object (value, self->x);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_data_dist1d_dispose (GObject *object)
{
  NcmDataDist1d *dist1d             = NCM_DATA_DIST1D (object);
  NcmDataDist1dPrivate * const self = ncm_data_dist1d_get_instance_private (dist1d);

  ncm_vector_clear (&self->x);

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
static void _ncm_data_dist1d_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
static void _ncm_data_dist1d_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng);
static void _ncm_data_dist1d_set_size (NcmDataDist1d *dist1d, guint np);
static guint _ncm_data_dist1d_get_size (NcmDataDist1d *dist1d);

static void
ncm_data_dist1d_class_init (NcmDataDist1dClass *klass)
{
  GObjectClass *object_class       = G_OBJECT_CLASS (klass);
  NcmDataDist1dClass *dist1d_class = NCM_DATA_DIST1D_CLASS (klass);
  NcmDataClass *data_class         = NCM_DATA_CLASS (klass);

  object_class->constructed  = &_ncm_data_dist1d_constructed;
  object_class->set_property = &_ncm_data_dist1d_set_property;
  object_class->get_property = &_ncm_data_dist1d_get_property;

  object_class->dispose  = &ncm_data_dist1d_dispose;
  object_class->finalize = &ncm_data_dist1d_finalize;

  g_object_class_install_property (object_class,
                                   PROP_NPOINTS,
                                   g_param_spec_uint ("n-points",
                                                      NULL,
                                                      "Data sample size",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_VECTOR,
                                   g_param_spec_object ("vector",
                                                        NULL,
                                                        "Data vector",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->bootstrap  = TRUE;
  data_class->get_length = &_ncm_data_dist1d_get_length;
  data_class->begin      = NULL;

  data_class->resample  = &_ncm_data_dist1d_resample;
  data_class->m2lnL_val = &_ncm_data_dist1d_m2lnL_val;

  dist1d_class->m2lnL_val = NULL;
  dist1d_class->inv_pdf   = NULL;
  dist1d_class->set_size  = &_ncm_data_dist1d_set_size;
  dist1d_class->get_size  = &_ncm_data_dist1d_get_size;
}

static guint
_ncm_data_dist1d_get_length (NcmData *data)
{
  NcmDataDist1d *dist1d             = NCM_DATA_DIST1D (data);
  NcmDataDist1dPrivate * const self = ncm_data_dist1d_get_instance_private (dist1d);

  return self->np;
}

static void
_ncm_data_dist1d_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcmDataDist1d *dist1d             = NCM_DATA_DIST1D (data);
  NcmDataDist1dPrivate * const self = ncm_data_dist1d_get_instance_private (dist1d);
  NcmDataDist1dClass *dist1d_class  = NCM_DATA_DIST1D_GET_CLASS (data);
  guint i;

  *m2lnL = 0.0;

  if (!ncm_data_bootstrap_enabled (data))
  {
    for (i = 0; i < self->np; i++)
    {
      const gdouble x_i = ncm_vector_get (self->x, i);

      *m2lnL += dist1d_class->m2lnL_val (dist1d, mset, x_i);
    }
  }
  else
  {
    NcmBootstrap *bstrap = ncm_data_peek_bootstrap (data);
    const guint bsize    = ncm_bootstrap_get_bsize (bstrap);

    for (i = 0; i < bsize; i++)
    {
      guint k           = ncm_bootstrap_get (bstrap, i);
      const gdouble x_i = ncm_vector_get (self->x, k);

      *m2lnL += dist1d_class->m2lnL_val (dist1d, mset, x_i);
    }
  }

  return;
}

static void
_ncm_data_dist1d_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng)
{
  NcmDataDist1d *dist1d             = NCM_DATA_DIST1D (data);
  NcmDataDist1dPrivate * const self = ncm_data_dist1d_get_instance_private (dist1d);
  NcmDataDist1dClass *dist1d_class  = NCM_DATA_DIST1D_GET_CLASS (data);
  guint i;

  if (dist1d_class->inv_pdf == NULL)
    g_error ("_ncm_data_dist1d_resample: This object do not implement the inverse of the pdf.");

  ncm_rng_lock (rng);

  for (i = 0; i < self->np; i++)
  {
    const gdouble u_i = gsl_rng_uniform (rng->r);
    const gdouble x_i = dist1d_class->inv_pdf (dist1d, mset, u_i);

    ncm_vector_set (self->x, i, x_i);
  }

  ncm_rng_unlock (rng);
}

static void
_ncm_data_dist1d_set_size (NcmDataDist1d *dist1d, guint np)
{
  NcmData *data                     = NCM_DATA (dist1d);
  NcmDataDist1dPrivate * const self = ncm_data_dist1d_get_instance_private (dist1d);

  if ((np == 0) || (np != self->np))
  {
    self->np = 0;
    ncm_vector_clear (&self->x);
    ncm_data_set_init (data, FALSE);
  }

  if ((np != 0) && (np != self->np))
  {
    NcmBootstrap *bstrap = ncm_data_peek_bootstrap (data);

    self->np = np;
    self->x  = ncm_vector_new (self->np);

    if (ncm_data_bootstrap_enabled (data))
    {
      ncm_bootstrap_set_fsize (bstrap, np);
      ncm_bootstrap_set_bsize (bstrap, np);
    }

    ncm_data_set_init (data, FALSE);
  }
}

static guint
_ncm_data_dist1d_get_size (NcmDataDist1d *dist1d)
{
  NcmDataDist1dPrivate * const self = ncm_data_dist1d_get_instance_private (dist1d);

  return self->np;
}

/**
 * ncm_data_dist1d_set_size: (virtual set_size)
 * @dist1d: a #NcmDataDist1d
 * @np: data size.
 *
 * Sets the data size to @np.
 *
 */
void
ncm_data_dist1d_set_size (NcmDataDist1d *dist1d, guint np)
{
  NCM_DATA_DIST1D_GET_CLASS (dist1d)->set_size (dist1d, np);
}

/**
 * ncm_data_dist1d_get_size: (virtual get_size)
 * @dist1d: a #NcmDataDist1d
 *
 * Gets the data size.
 *
 * Returns: Data size.
 *
 */
guint
ncm_data_dist1d_get_size (NcmDataDist1d *dist1d)
{
  return NCM_DATA_DIST1D_GET_CLASS (dist1d)->get_size (dist1d);
}

/**
 * ncm_data_dist1d_get_data:
 * @dist1d: a #NcmDataDist1d
 *
 * Gets the data #NcmVector.
 *
 * Returns: (transfer full): Data vector.
 */
NcmVector *
ncm_data_dist1d_get_data (NcmDataDist1d *dist1d)
{
  NcmDataDist1dPrivate * const self = ncm_data_dist1d_get_instance_private (dist1d);

  return ncm_vector_ref (self->x);
}

