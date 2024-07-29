/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_model_rosenbrock.h
 *
 *  Sat April 17 10:57:36 2021
 *  Copyright  2021  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_model_rosenbrock.h
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

#ifndef _NCM_MODEL_ROSENBROCK_H_
#define _NCM_MODEL_ROSENBROCK_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/math/ncm_model.h>

G_BEGIN_DECLS

#define NCM_TYPE_MODEL_ROSENBROCK (ncm_model_rosenbrock_get_type ())

G_DECLARE_FINAL_TYPE (NcmModelRosenbrock, ncm_model_rosenbrock, NCM, MODEL_ROSENBROCK, NcmModel)

/**
 * NcmModelRosenbrockSParams:
 * @NCM_MODEL_ROSENBROCK_X1: $x_1$
 * @NCM_MODEL_ROSENBROCK_X2: $x_2$
 *
 * Rosenbrock model parameters
 *
 */
typedef enum _NcmModelRosenbrockSParams
{
  NCM_MODEL_ROSENBROCK_X1,
  NCM_MODEL_ROSENBROCK_X2,
  /* < private > */
  NNCM_MODEL_ROSENBROCK_SPARAM_LEN, /*< skip >*/
} NcmModelRosenbrockSParams;

NCM_MSET_MODEL_DECLARE_ID (ncm_model_rosenbrock);

NcmModelRosenbrock *ncm_model_rosenbrock_new (void);
NcmModelRosenbrock *ncm_model_rosenbrock_ref (NcmModelRosenbrock *mrb);
void ncm_model_rosenbrock_free (NcmModelRosenbrock *mrb);
void ncm_model_rosenbrock_clear (NcmModelRosenbrock **mrb);

G_END_DECLS

#endif /* _NCM_MODEL_ROSENBROCK_H_ */

