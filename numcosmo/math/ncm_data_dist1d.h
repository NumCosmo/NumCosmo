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

#define NCM_TYPE_DATA_DIST1D             (ncm_data_dist1d_get_type ())
#define NCM_DATA_DIST1D(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_DATA_DIST1D, NcmDataDist1d))
#define NCM_DATA_DIST1D_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_DATA_DIST1D, NcmDataDist1dClass))
#define NCM_IS_DATA_DIST1D(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_DATA_DIST1D))
#define NCM_IS_DATA_DIST1D_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_DATA_DIST1D))
#define NCM_DATA_DIST1D_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_DATA_DIST1D, NcmDataDist1dClass))

typedef struct _NcmDataDist1dClass NcmDataDist1dClass;
typedef struct _NcmDataDist1d NcmDataDist1d;

struct _NcmDataDist1dClass
{
  /*< private >*/
  NcmDataClass parent_class;
  gdouble (*m2lnL_val) (NcmDataDist1d *dist1d, NcmMSet *mset, gdouble x);
  gdouble (*inv_pdf) (NcmDataDist1d *dist1d, NcmMSet *mset, gdouble u);  
  void (*set_size) (NcmDataDist1d *dist1d, guint np);
  guint (*get_size) (NcmDataDist1d *dist1d);
};

struct _NcmDataDist1d
{
  /*< private >*/
  NcmData parent_instance;
  guint np;
  NcmVector *x;
};

GType ncm_data_dist1d_get_type (void) G_GNUC_CONST;

void ncm_data_dist1d_set_size (NcmDataDist1d *dist1d, guint np);
guint ncm_data_dist1d_get_size (NcmDataDist1d *dist1d);
NcmVector *ncm_data_dist1d_get_data (NcmDataDist1d *dist1d);

G_END_DECLS

#endif /* _NCM_DATA_DIST1D_H_ */

