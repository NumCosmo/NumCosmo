/***************************************************************************
 *            likelihood.h
 *
 *  Mon Jul 16 18:05:56 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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

#ifndef _NC_LIKELIHOOD_H
#define _NC_LIKELIHOOD_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/data/dataset.h>
#include <numcosmo/math/ncm_mset_func.h>

G_BEGIN_DECLS

GType nc_likelihood_get_type (void);

typedef struct _NcLikelihood NcLikelihood;

/**
 * NcLikelihood:
 *
 * FIXME
 */
struct _NcLikelihood
{
  /*< private >*/
  NcDataSet *ds;
  GList *priors;
  gboolean clone;
};

NcLikelihood *nc_likelihood_new (NcDataSet *ds);
NcLikelihood *nc_likelihood_copy (NcLikelihood *lh_orig);
void nc_likelihood_free (NcLikelihood *lh);

void nc_likelihood_priors_add (NcLikelihood *lh, NcmMSetFunc *prior);
gint nc_likelihood_priors_length (NcLikelihood *lh);

gboolean nc_likelihood_has_leastsquares_J (NcLikelihood *lh);
gboolean nc_likelihood_has_m2lnL_grad (NcLikelihood *lh);

void nc_likelihood_leastsquares_f (NcLikelihood *lh, NcmMSet *mset, NcmVector *f);
void nc_likelihood_data_leastsquares_f (NcLikelihood *lh, NcmMSet *mset, NcmVector *data_f);
void nc_likelihood_priors_leastsquares_f (NcLikelihood *lh, NcmMSet *mset, NcmVector *priors_f);
void nc_likelihood_leastsquares_J (NcLikelihood *lh, NcmMSet *mset, NcmMatrix *J);
void nc_likelihood_leastsquares_f_J (NcLikelihood *lh, NcmMSet *mset, NcmVector *f, NcmMatrix *J);

void nc_likelihood_m2lnL_val (NcLikelihood *lh, NcmMSet *mset, gdouble *m2lnL);
void nc_likelihood_data_m2lnL_val (NcLikelihood *lh, NcmMSet *mset, gdouble *data_m2lnL);
void nc_likelihood_priors_m2lnL_val (NcLikelihood *lh, NcmMSet *mset, gdouble *priors_m2lnL);
void nc_likelihood_m2lnL_grad (NcLikelihood *lh, NcmMSet *mset, NcmVector *grad);
void nc_likelihood_m2lnL_val_grad (NcLikelihood *lh, NcmMSet *mset, gdouble *m2lnL, NcmVector *grad);

G_END_DECLS

#endif /* _LIKELIHOOD_H */
