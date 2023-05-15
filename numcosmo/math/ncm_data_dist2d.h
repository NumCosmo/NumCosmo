/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_data_dist2d.h
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

#ifndef _NCM_DATA_DIST2D_H_
#define _NCM_DATA_DIST2D_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_data.h>

G_BEGIN_DECLS

#define NCM_TYPE_DATA_DIST2D (ncm_data_dist2d_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmDataDist2d, ncm_data_dist2d, NCM, DATA_DIST2D, NcmData)

struct _NcmDataDist2dClass
{
  /*< private >*/
  NcmDataClass parent_class;

  gdouble (*m2lnL_val) (NcmDataDist2d *dist2d, NcmMSet *mset, gdouble x, gdouble y);
  void (*inv_pdf) (NcmDataDist2d *dist2d, NcmMSet *mset, gdouble u, gdouble v, gdouble *x, gdouble *y);
  void (*set_size) (NcmDataDist2d *dist2d, guint np);
  guint (*get_size) (NcmDataDist2d *dist2d);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[14];
};

void ncm_data_dist2d_set_size (NcmDataDist2d *dist2d, guint np);
guint ncm_data_dist2d_get_size (NcmDataDist2d *dist2d);
NcmMatrix *ncm_data_dist2d_get_data (NcmDataDist2d *dist2d);

G_END_DECLS

#endif /* _NCM_DATA_DIST2D_H_ */

