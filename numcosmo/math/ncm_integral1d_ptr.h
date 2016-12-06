/***************************************************************************
 *            ncm_integral1d_ptr.h
 *
 *  Mon December 05 17:23:53 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_integral1d_ptr.h
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_INTEGRAL1D_PTR_H_
#define _NCM_INTEGRAL1D_PTR_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_integral1d.h>

G_BEGIN_DECLS

#define NCM_TYPE_INTEGRAL1D_PTR             (ncm_integral1d_ptr_get_type ())
#define NCM_INTEGRAL1D_PTR(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_INTEGRAL1D_PTR, NcmIntegral1dPtr))
#define NCM_INTEGRAL1D_PTR_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_INTEGRAL1D_PTR, NcmIntegral1dPtrClass))
#define NCM_IS_INTEGRAL1D_PTR(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_INTEGRAL1D_PTR))
#define NCM_IS_INTEGRAL1D_PTR_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_INTEGRAL1D_PTR))
#define NCM_INTEGRAL1D_PTR_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_INTEGRAL1D_PTR, NcmIntegral1dPtrClass))

typedef struct _NcmIntegral1dPtrClass NcmIntegral1dPtrClass;
typedef struct _NcmIntegral1dPtr NcmIntegral1dPtr;
typedef struct _NcmIntegral1dPtrPrivate NcmIntegral1dPtrPrivate;
typedef gdouble (*NcmIntegral1dPtrF) (gpointer userdata, const gdouble x, const gdouble w);

struct _NcmIntegral1dPtrClass
{
  /*< private >*/
  NcmIntegral1dClass parent_class;
};

struct _NcmIntegral1dPtr
{
  /*< private >*/
  NcmIntegral1d parent_instance;
  NcmIntegral1dPtrPrivate *priv;
};

GType ncm_integral1d_ptr_get_type (void) G_GNUC_CONST;

NcmIntegral1dPtr *ncm_integral1d_ptr_new (NcmIntegral1dPtrF F, GDestroyNotify userfree);
NcmIntegral1dPtr *ncm_integral1d_ptr_new_full (NcmIntegral1dPtrF F, GDestroyNotify userfree, const gdouble reltol, const gdouble abstol, const guint partition, const guint rule);

NcmIntegral1dPtr *ncm_integral1d_ptr_ref (NcmIntegral1dPtr *int1d_ptr);
void ncm_integral1d_ptr_free (NcmIntegral1dPtr *int1d_ptr);
void ncm_integral1d_ptr_clear (NcmIntegral1dPtr **int1d_ptr);

void ncm_integral1d_ptr_set_userdata (NcmIntegral1dPtr *int1d_ptr, gpointer userdata);

G_END_DECLS

#endif /* _NCM_INTEGRAL1D_PTR_H_ */
