/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_data_dist1d.h
 *
 *  Thu Apr 15 11:17:25 2010
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

#ifndef _NCM_DATA_DIST1D_H_
#define _NCM_DATA_DIST1D_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_data.h>

G_BEGIN_DECLS

#define NCM_TYPE_DATA_DIST1D (ncm_data_dist1d_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmDataDist1d, ncm_data_dist1d, NCM, DATA_DIST1D, NcmData)

struct _NcmDataDist1dClass
{
  /*< private >*/
  NcmDataClass parent_class;

  gdouble (*dist1d_m2lnL_val) (NcmDataDist1d *dist1d, NcmMSet *mset, gdouble x);
  gdouble (*inv_pdf) (NcmDataDist1d *dist1d, NcmMSet *mset, gdouble u);
  void (*set_size) (NcmDataDist1d *dist1d, guint np);
  guint (*get_size) (NcmDataDist1d *dist1d);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[14];
};

void ncm_data_dist1d_set_size (NcmDataDist1d *dist1d, guint np);
guint ncm_data_dist1d_get_size (NcmDataDist1d *dist1d);
NcmVector *ncm_data_dist1d_get_data (NcmDataDist1d *dist1d);

G_END_DECLS

#endif /* _NCM_DATA_DIST1D_H_ */

