/***************************************************************************
 *            ncm_data_dist2d.c
 *
 *  Fri Sep 1 15:19:32 2017
 *  Copyright  2017  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2017 Mariana Penna Lima <pennalima@gmail.com>
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
 * SECTION:ncm_data_dist2d
 * @title: NcmDataDist2d
 * @short_description: Abstract class for two-variables distribution data.
 *
 * This object is designate to data that is described by a bivariate and arbitrary distribution.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_data_dist2d.h"
#include "math/ncm_rng.h"

enum
{
  PROP_0,
  PROP_NPOINTS,
  PROP_MATRIX,
  PROP_SIZE,
};

G_DEFINE_ABSTRACT_TYPE (NcmDataDist2d, ncm_data_dist2d, NCM_TYPE_DATA);

static void
ncm_data_dist2d_init (NcmDataDist2d *dist2d)
{
  dist2d->np = 0;
  dist2d->m  = NULL;
}

static void
_ncm_data_dist2d_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_data_dist2d_parent_class)->constructed (object);
}

static void
_ncm_data_dist2d_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmDataDist2d *dist2d = NCM_DATA_DIST2D (object);

  g_return_if_fail (NCM_IS_DATA_DIST2D (object));

  switch (prop_id)
  {
    case PROP_NPOINTS:
      ncm_data_dist2d_set_size (dist2d, g_value_get_uint (value));
      break;
    case PROP_MATRIX:
      ncm_matrix_substitute (&dist2d->m, g_value_get_object (value), TRUE);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_data_dist2d_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmDataDist2d *dist2d = NCM_DATA_DIST2D (object);

  g_return_if_fail (NCM_IS_DATA_DIST2D (object));

  switch (prop_id)
  {
    case PROP_NPOINTS:
      g_value_set_uint (value, ncm_data_dist2d_get_size (dist2d));
      break;
    case PROP_MATRIX:
      g_value_set_object (value, dist2d->m);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_data_dist2d_dispose (GObject *object)
{
  NcmDataDist2d *dist2d = NCM_DATA_DIST2D (object);

  ncm_matrix_clear (&dist2d->m);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_dist2d_parent_class)->dispose (object);
}

static void
ncm_data_dist2d_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_dist2d_parent_class)->finalize (object);
}

static guint _ncm_data_dist2d_get_length (NcmData *data);
static void _ncm_data_dist2d_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
static void _ncm_data_dist2d_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng);
static void _ncm_data_dist2d_set_size (NcmDataDist2d *dist2d, guint np);
static guint _ncm_data_dist2d_get_size (NcmDataDist2d *dist2d);

static void
ncm_data_dist2d_class_init (NcmDataDist2dClass *klass)
{
  GObjectClass *object_class       = G_OBJECT_CLASS (klass);
  NcmDataDist2dClass *dist2d_class = NCM_DATA_DIST2D_CLASS (klass);
  NcmDataClass *data_class         = NCM_DATA_CLASS (klass);

  object_class->constructed  = &_ncm_data_dist2d_constructed;
  object_class->set_property = &_ncm_data_dist2d_set_property;
  object_class->get_property = &_ncm_data_dist2d_get_property;

  object_class->dispose  = &ncm_data_dist2d_dispose;
  object_class->finalize = &ncm_data_dist2d_finalize;

  g_object_class_install_property (object_class,
                                   PROP_NPOINTS,
                                   g_param_spec_uint ("n-points",
                                                      NULL,
                                                      "Data sample size",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_MATRIX,
                                   g_param_spec_object ("matrix",
                                                        NULL,
                                                        "Data matrix",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->bootstrap  = TRUE;
  data_class->get_length = &_ncm_data_dist2d_get_length;
  data_class->begin      = NULL;

  data_class->resample  = &_ncm_data_dist2d_resample;
  data_class->m2lnL_val = &_ncm_data_dist2d_m2lnL_val;

  dist2d_class->m2lnL_val = NULL;
  dist2d_class->inv_pdf   = NULL;
  dist2d_class->set_size  = &_ncm_data_dist2d_set_size;
  dist2d_class->get_size  = &_ncm_data_dist2d_get_size;
}

static guint
_ncm_data_dist2d_get_length (NcmData *data)
{
  NcmDataDist2d *dist2d = NCM_DATA_DIST2D (data);

  return dist2d->np;
}

static void
_ncm_data_dist2d_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcmDataDist2d *dist2d            = NCM_DATA_DIST2D (data);
  NcmDataDist2dClass *dist2d_class = NCM_DATA_DIST2D_GET_CLASS (data);
  guint i;

  *m2lnL = 0.0;

  if (!ncm_data_bootstrap_enabled (data))
  {
    for (i = 0; i < dist2d->np; i++)
    {
      const gdouble x_i = ncm_matrix_get (dist2d->m, i, 0);
      const gdouble y_i = ncm_matrix_get (dist2d->m, i, 1);

      *m2lnL += dist2d_class->m2lnL_val (dist2d, mset, x_i, y_i);
    }
  }
  else
  {
    const guint bsize = ncm_bootstrap_get_bsize (data->bstrap);

    for (i = 0; i < bsize; i++)
    {
      guint k           = ncm_bootstrap_get (data->bstrap, i);
      const gdouble x_i = ncm_matrix_get (dist2d->m, k, 0);
      const gdouble y_i = ncm_matrix_get (dist2d->m, k, 1);

      *m2lnL += dist2d_class->m2lnL_val (dist2d, mset, x_i, y_i);
    }
  }

  return;
}

static void
_ncm_data_dist2d_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng)
{
  NcmDataDist2d *dist2d            = NCM_DATA_DIST2D (data);
  NcmDataDist2dClass *dist2d_class = NCM_DATA_DIST2D_GET_CLASS (data);
  guint i;

  if (dist2d_class->inv_pdf == NULL)
    g_error ("_ncm_data_dist2d_resample: This object do not implement the inverse of the pdf.");

  ncm_rng_lock (rng);

  for (i = 0; i < dist2d->np; i++)
  {
    gdouble x_i;
    gdouble y_i;
    const gdouble u_i = gsl_rng_uniform (rng->r);
    const gdouble v_i = gsl_rng_uniform (rng->r);

    dist2d_class->inv_pdf (dist2d, mset, u_i, v_i, &x_i, &y_i);
    ncm_matrix_set (dist2d->m, i, 0, x_i);
    ncm_matrix_set (dist2d->m, i, 1, y_i);
  }

  ncm_rng_unlock (rng);
}

static void
_ncm_data_dist2d_set_size (NcmDataDist2d *dist2d, guint np)
{
  NcmData *data = NCM_DATA (dist2d);

  if ((np == 0) || (np != dist2d->np))
  {
    dist2d->np = 0;
    ncm_matrix_clear (&dist2d->m);
    data->init = FALSE;
  }

  if ((np != 0) && (np != dist2d->np))
  {
    dist2d->np = np;
    dist2d->m  = ncm_matrix_new (dist2d->np, 2);

    if (ncm_data_bootstrap_enabled (data))
    {
      ncm_bootstrap_set_fsize (data->bstrap, np);
      ncm_bootstrap_set_bsize (data->bstrap, np);
    }

    data->init = FALSE;
  }
}

static guint
_ncm_data_dist2d_get_size (NcmDataDist2d *dist2d)
{
  return dist2d->np;
}

/**
 * ncm_data_dist2d_set_size: (virtual set_size)
 * @dist2d: a #NcmDataDist2d
 * @np: data size.
 *
 * Sets the data size to @np.
 *
 */
void
ncm_data_dist2d_set_size (NcmDataDist2d *dist2d, guint np)
{
  NCM_DATA_DIST2D_GET_CLASS (dist2d)->set_size (dist2d, np);
}

/**
 * ncm_data_dist2d_get_size: (virtual get_size)
 * @dist2d: a #NcmDataDist2d
 *
 * Gets the data size.
 *
 * Returns: Data size.
 *
 */
guint
ncm_data_dist2d_get_size (NcmDataDist2d *dist2d)
{
  return NCM_DATA_DIST2D_GET_CLASS (dist2d)->get_size (dist2d);
}

/**
 * ncm_data_dist2d_get_data:
 * @dist2d: a #NcmDataDist2d
 *
 * Gets the data #NcmMatrix.
 *
 * Returns: (transfer full): Data matrix.
 */
NcmMatrix *
ncm_data_dist2d_get_data (NcmDataDist2d *dist2d)
{
  return ncm_matrix_ref (dist2d->m);
}

