/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_data_funnel.c
 *
 *  Wed May 12 21:32:23 2021
 *  Copyright  2021  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_data_funnel.c
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
 * SECTION:ncm_data_funnel
 * @title: NcmDataFunnel
 * @short_description: Funnel distribution.
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_data_funnel.h"
#include "math/ncm_model_funnel.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmDataFunnelPrivate
{
  gint unused;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmDataFunnel, ncm_data_funnel, NCM_TYPE_DATA);

static void
ncm_data_funnel_init (NcmDataFunnel *dfu)
{
  dfu->priv = ncm_data_funnel_get_instance_private (dfu);
}

static void
ncm_data_funnel_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_data_funnel_parent_class)->constructed (object);

  ncm_data_set_init (NCM_DATA (object), TRUE);
}

static void
ncm_data_funnel_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_funnel_parent_class)->finalize (object);
}

static guint _ncm_data_funnel_get_length (NcmData *data);
static guint _ncm_data_funnel_get_dof (NcmData *data);
static void _ncm_data_funnel_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);

static void
ncm_data_funnel_class_init (NcmDataFunnelClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);

  object_class->constructed = ncm_data_funnel_constructed;
  object_class->finalize    = ncm_data_funnel_finalize;

  data_class->get_length = &_ncm_data_funnel_get_length;
  data_class->get_dof    = &_ncm_data_funnel_get_dof;
  data_class->m2lnL_val  = &_ncm_data_funnel_m2lnL_val;
}

static guint _ncm_data_funnel_get_length (NcmData *data) { return 10; }
static guint _ncm_data_funnel_get_dof (NcmData *data) { return 10; }

static void
_ncm_data_funnel_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcmModelFunnel *mrb    = NCM_MODEL_FUNNEL (ncm_mset_peek (mset, ncm_model_funnel_id ()));
  const gdouble nu       = ncm_model_param_get (NCM_MODEL (mrb), NCM_MODEL_FUNNEL_NU);
  const gdouble sigma_nu = exp (0.5 * nu);
  const guint x_len      = ncm_model_vparam_len (NCM_MODEL (mrb), NCM_MODEL_FUNNEL_X);
  const guint x_0_i      = ncm_model_vparam_index (NCM_MODEL (mrb), NCM_MODEL_FUNNEL_X, 0);
  gint i;
  
  /*printf ("nu: % 22.15g\n", nu);*/
  m2lnL[0] = x_len * nu + gsl_pow_2 (nu / 3.0);
  for (i = 0; i < x_len; i++)
  {
    const gdouble x_i = ncm_model_param_get (NCM_MODEL (mrb), x_0_i + i);
    /*printf ("x[%d] = % 22.15g\n", i, x_i);*/
    m2lnL[0] += gsl_pow_2 (x_i / sigma_nu);
  }
}

/**
 * ncm_data_funnel_new:
 * 
 * Creates a new Funnel data.
 * 
 * Returns: the newly created object.
 */ 
NcmDataFunnel *
ncm_data_funnel_new (void)
{
  NcmDataFunnel *dfu = g_object_new (NCM_TYPE_DATA_FUNNEL,
                                         NULL);
  return dfu;
}

/**
 * ncm_data_funnel_ref:
 * @dfu: a #NcmDataFunnel
 * 
 * Increases the reference count of @dfu by one.
 * 
 * Returns: (transfer full): @dfu
 */
NcmDataFunnel *
ncm_data_funnel_ref (NcmDataFunnel *dfu)
{
  return g_object_ref (dfu);
}

/**
 * ncm_data_funnel_free:
 * @dfu: a #NcmDataFunnel
 * 
 * Decreases the reference count of @dfu by one.
 * 
 */
void 
ncm_data_funnel_free (NcmDataFunnel *dfu)
{
  g_object_unref (dfu);
}

/**
 * ncm_data_funnel_clear:
 * @dfu: a #NcmDataFunnel
 * 
 * If @dfu is different from NULL, decreases the reference count of
 * @dfu by one and sets @dfu to NULL.
 * 
 */
void 
ncm_data_funnel_clear (NcmDataFunnel **dfu)
{
  g_clear_object (dfu);
}
