/***************************************************************************
 *            ncm_mset_trans_kern.h
 *
 *  Fri August 29 18:56:39 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mset_trans_kern.h
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_MSET_TRANS_KERN_H_
#define _NCM_MSET_TRANS_KERN_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/math/ncm_rng.h>

G_BEGIN_DECLS

#define NCM_TYPE_MSET_TRANS_KERN             (ncm_mset_trans_kern_get_type ())
#define NCM_MSET_TRANS_KERN(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_MSET_TRANS_KERN, NcmMSetTransKern))
#define NCM_MSET_TRANS_KERN_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_MSET_TRANS_KERN, NcmMSetTransKernClass))
#define NCM_IS_MSET_TRANS_KERN(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_MSET_TRANS_KERN))
#define NCM_IS_MSET_TRANS_KERN_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_MSET_TRANS_KERN))
#define NCM_MSET_TRANS_KERN_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_MSET_TRANS_KERN, NcmMSetTransKernClass))

typedef struct _NcmMSetTransKernClass NcmMSetTransKernClass;
typedef struct _NcmMSetTransKern NcmMSetTransKern;

struct _NcmMSetTransKernClass
{
  /*< private >*/
  GObjectClass parent_class;
  gboolean bernoulli_scheme;

  void (*set_mset) (NcmMSetTransKern *tkern, NcmMSet *mset);
  void (*generate) (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng);
  gdouble (*pdf) (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar);
  void (*reset) (NcmMSetTransKern *tkern);

  const gchar *(*get_name) (NcmMSetTransKern *tkern);
};

struct _NcmMSetTransKern
{
  /*< private >*/
  GObject parent_instance;
  NcmMSet *mset;
  NcmVector *theta;
};

GType ncm_mset_trans_kern_get_type (void) G_GNUC_CONST;

NcmMSetTransKern *ncm_mset_trans_kern_ref (NcmMSetTransKern *tkern);
void ncm_mset_trans_kern_free (NcmMSetTransKern *tkern);
void ncm_mset_trans_kern_clear (NcmMSetTransKern **tkern);

void ncm_mset_trans_kern_set_mset (NcmMSetTransKern *tkern, NcmMSet *mset);
void ncm_mset_trans_kern_set_prior (NcmMSetTransKern *tkern, NcmVector *theta);
void ncm_mset_trans_kern_set_prior_from_mset (NcmMSetTransKern *tkern);
void ncm_mset_trans_kern_generate (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng);
gdouble ncm_mset_trans_kern_pdf (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar);
void ncm_mset_trans_kern_prior_sample (NcmMSetTransKern *tkern, NcmVector *thetastar, NcmRNG *rng);
gdouble ncm_mset_trans_kern_prior_pdf (NcmMSetTransKern *tkern, NcmVector *thetastar);
void ncm_mset_trans_kern_reset (NcmMSetTransKern *tkern);

const gchar *ncm_mset_trans_kern_get_name (NcmMSetTransKern *tkern);

G_END_DECLS

#endif /* _NCM_MSET_TRANS_KERN_H_ */

