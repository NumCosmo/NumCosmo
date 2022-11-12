/***************************************************************************
 *            ncm_mset_trans_kern_cat.h
 *
 *  Fri September 23 12:36:57 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mset_trans_kern_cat.h
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_MSET_TRANS_KERN_CAT_H_
#define _NCM_MSET_TRANS_KERN_CAT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mset_trans_kern.h>
#include <numcosmo/math/ncm_mset_catalog.h>
#include <numcosmo/math/ncm_stats_dist.h>

G_BEGIN_DECLS

#define NCM_TYPE_MSET_TRANS_KERN_CAT             (ncm_mset_trans_kern_cat_get_type ())
#define NCM_MSET_TRANS_KERN_CAT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_MSET_TRANS_KERN_CAT, NcmMSetTransKernCat))
#define NCM_MSET_TRANS_KERN_CAT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_MSET_TRANS_KERN_CAT, NcmMSetTransKernCatClass))
#define NCM_IS_MSET_TRANS_KERN_CAT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_MSET_TRANS_KERN_CAT))
#define NCM_IS_MSET_TRANS_KERN_CAT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_MSET_TRANS_KERN_CAT))
#define NCM_MSET_TRANS_KERN_CAT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_MSET_TRANS_KERN_CAT, NcmMSetTransKernCatClass))

typedef struct _NcmMSetTransKernCatClass NcmMSetTransKernCatClass;
typedef struct _NcmMSetTransKernCat NcmMSetTransKernCat;
typedef struct _NcmMSetTransKernCatPrivate NcmMSetTransKernCatPrivate;

struct _NcmMSetTransKernCatClass
{
  /*< private >*/
  NcmMSetTransKernClass parent_class;
};

struct _NcmMSetTransKernCat
{
  /*< private >*/
  NcmMSetTransKern parent_instance;
  NcmMSetTransKernCatPrivate *priv;
};

/**
 * NcmMSetTransKernCatSampling:
 * @NCM_MSET_TRANS_KERN_CAT_SAMPLING_CHOOSE: sampling by choosing radomly points from the catalog
 * @NCM_MSET_TRANS_KERN_CAT_SAMPLING_RBF_INTERP: sampling from the RBF interpolation built using the last ensemble from the catalog 
 * @NCM_MSET_TRANS_KERN_CAT_SAMPLING_KDE: sampling from the KDE estimate built using the last ensemble from the catalog 
 * 
 * Sampling methods.
 * 
 */ 
typedef enum _NcmMSetTransKernCatSampling
{
  NCM_MSET_TRANS_KERN_CAT_SAMPLING_CHOOSE = 0,
  NCM_MSET_TRANS_KERN_CAT_SAMPLING_RBF_INTERP,
  NCM_MSET_TRANS_KERN_CAT_SAMPLING_KDE,
  /* < private > */
  NCM_MSET_TRANS_KERN_CAT_SAMPLING_LEN, /*< skip >*/
} NcmMSetTransKernCatSampling;

GType ncm_mset_trans_kern_cat_get_type (void) G_GNUC_CONST;

NcmMSetTransKernCat *ncm_mset_trans_kern_cat_new (NcmMSetCatalog *mcat, NcmStatsDist *sd);

void ncm_mset_trans_kern_cat_set_sampling (NcmMSetTransKernCat *tcat, NcmMSetTransKernCatSampling sampling);
NcmMSetTransKernCatSampling ncm_mset_trans_kern_cat_get_sampling (NcmMSetTransKernCat *tcat);

G_END_DECLS

#endif /* _NCM_MSET_TRANS_KERN_CAT_H_ */
