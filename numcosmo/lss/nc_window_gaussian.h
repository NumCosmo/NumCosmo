/***************************************************************************
 *            nc_window_gaussian.h
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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

#ifndef _NC_WINDOW_GAUSSIAN_H_
#define _NC_WINDOW_GAUSSIAN_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_window.h>

G_BEGIN_DECLS

#define NC_TYPE_WINDOW_GAUSSIAN             (nc_window_gaussian_get_type ())
#define NC_WINDOW_GAUSSIAN(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_WINDOW_GAUSSIAN, NcWindowGaussian))
#define NC_WINDOW_GAUSSIAN_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_WINDOW_GAUSSIAN, NcWindowGaussianClass))
#define NC_IS_WINDOW_GAUSSIAN(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_WINDOW_GAUSSIAN))
#define NC_IS_WINDOW_GAUSSIAN_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_WINDOW_GAUSSIAN))
#define NC_WINDOW_GAUSSIAN_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_WINDOW_GAUSSIAN, NcWindowGaussianClass))

typedef struct _NcWindowGaussianClass NcWindowGaussianClass;
typedef struct _NcWindowGaussian NcWindowGaussian;

struct _NcWindowGaussianClass
{
  /*< private > */
  NcWindowClass parent_class;
};

struct _NcWindowGaussian
{
  /*< private > */
  NcWindow parent_instance;
};

GType nc_window_gaussian_get_type (void) G_GNUC_CONST;

NcWindow *nc_window_gaussian_new (void);

#define NC_WINDOW_VOLUME_GAUSSIAN (sqrt(2.0 * M_PI)*sqrt(2.0 * M_PI)*sqrt(2.0 * M_PI))  /* (2.0 \Pi)^(3/2) */

G_END_DECLS

#endif /* _NC_WINDOW_GAUSSIAN_H_ */
