/***************************************************************************
 *            nc_cbe_precision.h
 *
 *  Sun October 25 20:45:31 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_cbe_precision.h
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_CBE_PRECISION_H_
#define _NC_CBE_PRECISION_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

G_BEGIN_DECLS

#define NC_TYPE_CBE_PRECISION             (nc_cbe_precision_get_type ())
#define NC_CBE_PRECISION(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_CBE_PRECISION, NcCBEPrecision))
#define NC_CBE_PRECISION_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_CBE_PRECISION, NcCBEPrecisionClass))
#define NC_IS_CBE_PRECISION(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_CBE_PRECISION))
#define NC_IS_CBE_PRECISION_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_CBE_PRECISION))
#define NC_CBE_PRECISION_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_CBE_PRECISION, NcCBEPrecisionClass))

typedef struct _NcCBEPrecisionClass NcCBEPrecisionClass;
typedef struct _NcCBEPrecision NcCBEPrecision;

typedef struct _NcCBEPrecisionPrivate NcCBEPrecisionPrivate;

struct _NcCBEPrecisionClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcCBEPrecision
{
  /*< private >*/
  GObject parent_instance;
  NcCBEPrecisionPrivate *priv;
};

GType nc_cbe_precision_get_type (void) G_GNUC_CONST;

NcCBEPrecision *nc_cbe_precision_ref (NcCBEPrecision *cbe_prec);
NcCBEPrecision *nc_cbe_precision_new (void);
void nc_cbe_precision_free (NcCBEPrecision *cbe_prec);
void nc_cbe_precision_clear (NcCBEPrecision **cbe_prec);
void nc_cbe_precision_assert_default (NcCBEPrecision *cbe_prec);

G_END_DECLS

#endif /* _NC_CBE_PRECISION_H_ */

