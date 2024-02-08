/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_model_funnel.h
 *
 *  Wed May 12 21:25:45 2021
 *  Copyright  2021  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_model_funnel.h
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

#ifndef _NCM_MODEL_FUNNEL_H_
#define _NCM_MODEL_FUNNEL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/math/ncm_model.h>

G_BEGIN_DECLS

#define NCM_TYPE_MODEL_FUNNEL (ncm_model_funnel_get_type ())

G_DECLARE_FINAL_TYPE (NcmModelFunnel, ncm_model_funnel, NCM, MODEL_FUNNEL, NcmModel)

/**
 * NcmModelFunnelSParams:
 * @NCM_MODEL_FUNNEL_NU: $\nu$
 *
 * Funnel model parameters.
 *
 */
typedef enum _NcmModelFunnelSParams
{
  NCM_MODEL_FUNNEL_NU,
  /* < private > */
  NNCM_MODEL_FUNNEL_SPARAM_LEN, /*< skip >*/
} NcmModelFunnelSParams;

/**
 * NcmModelFunnelVParams:
 * @NCM_MODEL_FUNNEL_X: $x$
 *
 * Funnel model parameters.
 *
 */
typedef enum _NcmModelFunnelVParams
{
  NCM_MODEL_FUNNEL_X,
  /* < private > */
  NNCM_MODEL_FUNNEL_VPARAM_LEN, /*< skip >*/
} NcmModelFunnelVParams;

NCM_MSET_MODEL_DECLARE_ID (ncm_model_funnel);

NcmModelFunnel *ncm_model_funnel_new (guint n);
NcmModelFunnel *ncm_model_funnel_ref (NcmModelFunnel *mfu);
void ncm_model_funnel_free (NcmModelFunnel *mfu);
void ncm_model_funnel_clear (NcmModelFunnel **mfu);

G_END_DECLS

#endif /* _NCM_MODEL_FUNNEL_H_ */
