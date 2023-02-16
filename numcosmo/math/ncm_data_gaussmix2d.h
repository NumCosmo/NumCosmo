/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_data_gaussmix2d.h
 *
 *  Sat April 17 11:11:28 2021
 *  Copyright  2021  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_data_gaussmix2d.h
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

#ifndef _NCM_DATA_GAUSSMIX2D_H_
#define _NCM_DATA_GAUSSMIX2D_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data_gauss_cov.h>

G_BEGIN_DECLS

#define NCM_TYPE_DATA_GAUSSMIX2D             (ncm_data_gaussmix2d_get_type ())
#define NCM_DATA_GAUSSMIX2D(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_DATA_GAUSSMIX2D, NcmDataGaussMix2D))
#define NCM_DATA_GAUSSMIX2D_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_DATA_GAUSSMIX2D, NcmDataGaussMix2DClass))
#define NCM_IS_DATA_GAUSSMIX2D(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_DATA_GAUSSMIX2D))
#define NCM_IS_DATA_GAUSSMIX2D_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_DATA_GAUSSMIX2D))
#define NCM_DATA_GAUSSMIX2D_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_DATA_GAUSSMIX2D, NcmDataGaussMix2DClass))

typedef struct _NcmDataGaussMix2DClass NcmDataGaussMix2DClass;
typedef struct _NcmDataGaussMix2D NcmDataGaussMix2D;
typedef struct _NcmDataGaussMix2DPrivate NcmDataGaussMix2DPrivate;

struct _NcmDataGaussMix2DClass
{
  /*< private >*/
  NcmDataClass parent_class;
};

struct _NcmDataGaussMix2D
{
  /*< private >*/
  NcmData parent_instance;
  NcmDataGaussMix2DPrivate *priv;
};

GType ncm_data_gaussmix2d_get_type (void) G_GNUC_CONST;

NcmDataGaussMix2D *ncm_data_gaussmix2d_new (void);
NcmDataGaussMix2D *ncm_data_gaussmix2d_ref (NcmDataGaussMix2D *gm2d);
void ncm_data_gaussmix2d_free (NcmDataGaussMix2D *gm2d);
void ncm_data_gaussmix2d_clear (NcmDataGaussMix2D **gm2d);

G_END_DECLS

#endif /* _NCM_DATA_GAUSSMIX2D_H_ */

