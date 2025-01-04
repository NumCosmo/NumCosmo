/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_data_gaussmix2d.c
 *
 *  Sat April 17 11:11:28 2021
 *  Copyright  2021  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_data_gaussmix2d.c
 * Copyright (C) 2021 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmDataGaussMix2D:
 *
 * Gaussian Mixture 2d distribution.
 *
 * Data object describing a Gaussian Mixture 2d distribution.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_data_gaussmix2d.h"
#include "math/ncm_model_rosenbrock.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcmDataGaussMix2DPrivate
{
  gint unused;
} NcmDataGaussMix2DPrivate;

struct _NcmDataGaussMix2D
{
  NcmData parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmDataGaussMix2D, ncm_data_gaussmix2d, NCM_TYPE_DATA)

static void
ncm_data_gaussmix2d_init (NcmDataGaussMix2D *gm2d)
{
}

static void
ncm_data_gaussmix2d_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_data_gaussmix2d_parent_class)->constructed (object);

  ncm_data_set_init (NCM_DATA (object), TRUE);
}

static void
ncm_data_gaussmix2d_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_gaussmix2d_parent_class)->finalize (object);
}

static guint _ncm_data_gaussmix2d_get_length (NcmData *data);
static guint _ncm_data_gaussmix2d_get_dof (NcmData *data);
static void _ncm_data_gaussmix2d_prepare (NcmData *data, NcmMSet *mset);
static void _ncm_data_gaussmix2d_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);

static void
ncm_data_gaussmix2d_class_init (NcmDataGaussMix2DClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);

  object_class->constructed = ncm_data_gaussmix2d_constructed;
  object_class->finalize    = ncm_data_gaussmix2d_finalize;

  data_class->get_length = &_ncm_data_gaussmix2d_get_length;
  data_class->get_dof    = &_ncm_data_gaussmix2d_get_dof;
  data_class->prepare    = &_ncm_data_gaussmix2d_prepare;
  data_class->m2lnL_val  = &_ncm_data_gaussmix2d_m2lnL_val;
}

static guint
_ncm_data_gaussmix2d_get_length (NcmData *data)
{
  return 10;
}

static guint
_ncm_data_gaussmix2d_get_dof (NcmData *data)
{
  return 10;
}

static void
_ncm_data_gaussmix2d_prepare (NcmData *data, NcmMSet *mset)
{
  NcmModelRosenbrock *mrb = NCM_MODEL_ROSENBROCK (ncm_mset_peek (mset, ncm_model_rosenbrock_id ()));

  g_assert (mrb != NULL);
  g_assert (NCM_IS_MODEL_ROSENBROCK (mrb));
}

static void
_ncm_data_gaussmix2d_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcmModelRosenbrock *mrb = NCM_MODEL_ROSENBROCK (ncm_mset_peek (mset, ncm_model_rosenbrock_id ()));
  const gdouble x1        = ncm_model_param_get (NCM_MODEL (mrb), NCM_MODEL_ROSENBROCK_X1);
  const gdouble x2        = ncm_model_param_get (NCM_MODEL (mrb), NCM_MODEL_ROSENBROCK_X2);
  const gdouble rho1      = +0.6;
  const gdouble rho2      = -0.6;
  const gdouble sigma1    = 0.4;
  const gdouble sigma2    = 0.2;
  const gdouble mu11      = -1.5;
  const gdouble mu21      = +0.0;
  const gdouble mu12      = +1.5;
  const gdouble mu22      = +0.0;
  const gdouble w1        = 0.5;
  const gdouble w2        = 0.5;
  const gdouble dx11      = (x1 - mu11) / sigma1;
  const gdouble dx21      = (x2 - mu21) / sigma1;
  const gdouble dx12      = (x1 - mu12) / sigma2;
  const gdouble dx22      = (x2 - mu22) / sigma2;
  const gdouble z1        = dx11 * dx11 - 2.0 * rho1 * dx11 * dx21 + dx21 * dx21;
  const gdouble z2        = dx12 * dx12 - 2.0 * rho2 * dx12 * dx22 + dx22 * dx22;
  const gdouble chi1      = z1 / (2.0 * (1.0 - rho1 * rho1));
  const gdouble chi2      = z2 / (2.0 * (1.0 - rho2 * rho2));
  const gdouble lnN1      = -log (2.0 * M_PI * sigma1 * sigma1 * sqrt (1.0 - rho1 * rho1));
  const gdouble lnN2      = -log (2.0 * M_PI * sigma2 * sigma2 * sqrt (1.0 - rho2 * rho2));
  const gdouble lnp1      = log (w1) + lnN1 - chi1;
  const gdouble lnp2      = log (w2) + lnN2 - chi2;

  if (lnp1 > lnp2)
    m2lnL[0] = -2.0 * (lnp1 + log1p (exp (lnp2 - lnp1)));
  else
    m2lnL[0] = -2.0 * (lnp2 + log1p (exp (lnp1 - lnp2)));
}

/**
 * ncm_data_gaussmix2d_new:
 *
 * Creates a new #NcmDataGaussMix2D.
 *
 * Returns: the newly created object.
 */
NcmDataGaussMix2D *
ncm_data_gaussmix2d_new (void)
{
  NcmDataGaussMix2D *gm2d = g_object_new (NCM_TYPE_DATA_GAUSSMIX2D,
                                          NULL);

  return gm2d;
}

/**
 * ncm_data_gaussmix2d_ref:
 * @gm2d: a #NcmDataGaussMix2D
 *
 * Increases the reference count of @gm2d by onG
 * Returns: (transfer full): @gm2d
 */
NcmDataGaussMix2D *
ncm_data_gaussmix2d_ref (NcmDataGaussMix2D *gm2d)
{
  return g_object_ref (gm2d);
}

/**
 * ncm_data_gaussmix2d_free:
 * @gm2d: a #NcmDataGaussMix2D
 *
 * Decreases the reference count of @gm2d by onG
 */
void
ncm_data_gaussmix2d_free (NcmDataGaussMix2D *gm2d)
{
  g_object_unref (gm2d);
}

/**
 * ncm_data_gaussmix2d_clear:
 * @gm2d: a #NcmDataGaussMix2D
 *
 * If @gm2d is different from NULL, decreases the reference count of
 * @gm2d by one and sets Gto NULL.
 *
 */
void
ncm_data_gaussmix2d_clear (NcmDataGaussMix2D **gm2d)
{
  g_clear_object (gm2d);
}

