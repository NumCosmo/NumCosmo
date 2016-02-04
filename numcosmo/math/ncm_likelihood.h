/***************************************************************************
 *            ncm_likelihood.h
 *
 *  Mon Jul 16 18:05:56 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_LIKELIHOOD_H_
#define _NCM_LIKELIHOOD_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_dataset.h>
#include <numcosmo/math/ncm_mset_func.h>

G_BEGIN_DECLS

#define NCM_TYPE_LIKELIHOOD             (ncm_likelihood_get_type ())
#define NCM_LIKELIHOOD(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_LIKELIHOOD, NcmLikelihood))
#define NCM_LIKELIHOOD_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_LIKELIHOOD, NcmLikelihoodClass))
#define NCM_IS_LIKELIHOOD(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_LIKELIHOOD))
#define NCM_IS_LIKELIHOOD_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_LIKELIHOOD))
#define NCM_LIKELIHOOD_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_LIKELIHOOD, NcmLikelihoodClass))

typedef struct _NcmLikelihoodClass NcmLikelihoodClass;
typedef struct _NcmLikelihood NcmLikelihood;

struct _NcmLikelihood
{
  /*< private >*/
  GObject parent_instance;
  NcmDataset *dset;
  GPtrArray *priors_f;
  GPtrArray *priors_m2lnL;
  NcmVector *m2lnL_v;
};

struct _NcmLikelihoodClass
{
  /*< private >*/
  GObjectClass parent_class;
};

GType ncm_likelihood_get_type (void) G_GNUC_CONST;

NcmLikelihood *ncm_likelihood_new (NcmDataset *dset);
NcmLikelihood *ncm_likelihood_ref (NcmLikelihood *lh);
NcmLikelihood *ncm_likelihood_dup (NcmLikelihood *lh, NcmSerialize *ser);
NcmLikelihood *ncm_likelihood_copy (NcmLikelihood *lh);
void ncm_likelihood_free (NcmLikelihood *lh);
void ncm_likelihood_clear (NcmLikelihood **lh);

void ncm_likelihood_priors_add (NcmLikelihood *lh, NcmMSetFunc *prior, gboolean is_m2lnL);
NcmMSetFunc *ncm_likelihood_priors_peek_f (NcmLikelihood *lh, guint i);
guint ncm_likelihood_priors_length_f (NcmLikelihood *lh);
NcmMSetFunc *ncm_likelihood_priors_peek_m2lnL (NcmLikelihood *lh, guint i);
guint ncm_likelihood_priors_length_m2lnL (NcmLikelihood *lh);

gboolean ncm_likelihood_has_leastsquares_J (NcmLikelihood *lh);
gboolean ncm_likelihood_has_m2lnL_grad (NcmLikelihood *lh);

void ncm_likelihood_priors_leastsquares_f (NcmLikelihood *lh, NcmMSet *mset, NcmVector *priors_f);
void ncm_likelihood_leastsquares_f (NcmLikelihood *lh, NcmMSet *mset, NcmVector *f);
void ncm_likelihood_leastsquares_J (NcmLikelihood *lh, NcmMSet *mset, NcmMatrix *J);
void ncm_likelihood_leastsquares_f_J (NcmLikelihood *lh, NcmMSet *mset, NcmVector *f, NcmMatrix *J);

void ncm_likelihood_priors_m2lnL_val (NcmLikelihood *lh, NcmMSet *mset, gdouble *priors_m2lnL);
void ncm_likelihood_priors_m2lnL_vec (NcmLikelihood *lh, NcmMSet *mset, NcmVector *priors_m2lnL_v);
void ncm_likelihood_m2lnL_val (NcmLikelihood *lh, NcmMSet *mset, gdouble *m2lnL);
void ncm_likelihood_m2lnL_grad (NcmLikelihood *lh, NcmMSet *mset, NcmVector *grad);
void ncm_likelihood_m2lnL_val_grad (NcmLikelihood *lh, NcmMSet *mset, gdouble *m2lnL, NcmVector *grad);

G_END_DECLS

#endif /* _NCM_LIKELIHOOD_H_ */

