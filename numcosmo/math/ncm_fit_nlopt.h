/***************************************************************************
 *            ncm_fit_nlopt.h
 *
 *  Sat Apr  3 16:07:17 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
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

#ifndef _NCM_FIT_NLOPT_H_
#define _NCM_FIT_NLOPT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_fit.h>
#ifdef NUMCOSMO_HAVE_NLOPT
#include <nlopt.h>

G_BEGIN_DECLS

#define NCM_TYPE_FIT_NLOPT             (ncm_fit_nlopt_get_type ())
#define NCM_FIT_NLOPT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_FIT_NLOPT, NcmFitNLOpt))
#define NCM_FIT_NLOPT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_FIT_NLOPT, NcmFitNLOptClass))
#define NCM_IS_FIT_NLOPT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_FIT_NLOPT))
#define NCM_IS_FIT_NLOPT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_FIT_NLOPT))
#define NCM_FIT_NLOPT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_FIT_NLOPT, NcmFitNLOptClass))

typedef struct _NcmFitNLOptClass NcmFitNLOptClass;
typedef struct _NcmFitNLOpt NcmFitNLOpt;

struct _NcmFitNLOpt
{
  /*< private >*/
  NcmFit parent_instance;
#ifdef NUMCOSMO_HAVE_NLOPT
#ifdef HAVE_NLOPT_2_2
  nlopt_opt nlopt;
  nlopt_opt local_nlopt;
#endif /* HAVE_NLOPT_2_2 */  
  nlopt_algorithm nlopt_algo;
  nlopt_algorithm local_nlopt_algo;
#endif /* NUMCOSMO_HAVE_NLOPT */
  NcmVector *lb;
  NcmVector *ub;
  NcmVector *pabs;
  NcmVector *pscale;
  gchar *desc;
  guint fparam_len;
};

struct _NcmFitNLOptClass
{
  /*< private >*/
  NcmFitClass parent_class;
};

GType ncm_fit_nlopt_get_type (void) G_GNUC_CONST;

#ifndef NUMCOSMO_GIR_SCAN
NcmFit *ncm_fit_nlopt_new (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, nlopt_algorithm algo);
NcmFit *ncm_fit_nlopt_local_new (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, nlopt_algorithm algo, nlopt_algorithm local_algo);
void ncm_fit_nlopt_set_algo (NcmFitNLOpt *fit_nlopt, nlopt_algorithm algo);
void ncm_fit_nlopt_set_local_algo (NcmFitNLOpt *fit_nlopt, nlopt_algorithm algo);
#endif /* NUMCOSMO_GIR_SCAN */
NcmFit *ncm_fit_nlopt_new_by_name (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, gchar *algo_name);
NcmFit *ncm_fit_nlopt_new_default (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype);

G_END_DECLS

#endif /* NUMCOSMO_HAVE_NLOPT */
#endif /* _NCM_FIT_NLOPT_H_ */
