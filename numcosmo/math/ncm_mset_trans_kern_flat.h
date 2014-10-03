/***************************************************************************
 *            ncm_mset_trans_kern_flat.h
 *
 *  Thu October 02 13:36:55 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_mset_trans_kern_flat.h
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_MSET_TRANS_KERN_FLAT_H_
#define _NCM_MSET_TRANS_KERN_FLAT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mset_trans_kern.h>

G_BEGIN_DECLS

#define NCM_TYPE_MSET_TRANS_KERN_FLAT             (ncm_mset_trans_kern_flat_get_type ())
#define NCM_MSET_TRANS_KERN_FLAT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_MSET_TRANS_KERN_FLAT, NcmMSetTransKernFlat))
#define NCM_MSET_TRANS_KERN_FLAT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_MSET_TRANS_KERN_FLAT, NcmMSetTransKernFlatClass))
#define NCM_IS_MSET_TRANS_KERN_FLAT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_MSET_TRANS_KERN_FLAT))
#define NCM_IS_MSET_TRANS_KERN_FLAT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_MSET_TRANS_KERN_FLAT))
#define NCM_MSET_TRANS_KERN_FLAT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_MSET_TRANS_KERN_FLAT, NcmMSetTransKernFlatClass))

typedef struct _NcmMSetTransKernFlatClass NcmMSetTransKernFlatClass;
typedef struct _NcmMSetTransKernFlat NcmMSetTransKernFlat;

struct _NcmMSetTransKernFlatClass
{
  /*< private >*/
  NcmMSetTransKernClass parent_class;
};

struct _NcmMSetTransKernFlat
{
  /*< private >*/
  NcmMSetTransKern parent_instance;
  gdouble parea;
};

GType ncm_mset_trans_kern_flat_get_type (void) G_GNUC_CONST;

NcmMSetTransKernFlat *ncm_mset_trans_kern_flat_new (void);

G_END_DECLS

#endif /* _NCM_MSET_TRANS_KERN_FLAT_H_ */
