/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_model_mvnd.h
 *
 *  Sun February 04 15:31:36 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_model_mvnd.h
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_MODEL_MVND_H_
#define _NCM_MODEL_MVND_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/math/ncm_model.h>

G_BEGIN_DECLS

#define NCM_TYPE_MODEL_MVND (ncm_model_mvnd_get_type ())

G_DECLARE_FINAL_TYPE (NcmModelMVND, ncm_model_mvnd, NCM, MODEL_MVND, NcmModel)

/**
 * NcmModelMVNDVParams:
 * @NCM_MODEL_MVND_MEAN: Mean vector
 *
 * MVND model parameters
 *
 */
typedef enum _NcmModelMVNDVParams
{
  NCM_MODEL_MVND_MEAN,
  /* < private > */
  NNCM_MODEL_MVND_VPARAM_LEN, /*< skip >*/
} NcmModelMVNDVParams;

NCM_MSET_MODEL_DECLARE_ID (ncm_model_mvnd);

NcmModelMVND *ncm_model_mvnd_new (const guint dim);
NcmModelMVND *ncm_model_mvnd_ref (NcmModelMVND *model_mvnd);
void ncm_model_mvnd_free (NcmModelMVND *model_mvnd);
void ncm_model_mvnd_clear (NcmModelMVND **model_mvnd);

void ncm_model_mvnd_mean (NcmModelMVND *model_mvnd, NcmVector *y);

G_END_DECLS

#endif /* _NCM_MODEL_MVND_H_ */

