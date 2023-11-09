/***************************************************************************
 *            ncm_fit_nlopt.h
 *
 *  Sat Apr  3 16:07:17 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
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

#ifndef _NCM_FIT_NLOPT_H_
#define _NCM_FIT_NLOPT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_fit.h>

#ifdef NUMCOSMO_HAVE_NLOPT
#include <numcosmo/ncm_fit_nlopt_enum.h>

G_BEGIN_DECLS

#define NCM_TYPE_FIT_NLOPT (ncm_fit_nlopt_get_type ())

G_DECLARE_FINAL_TYPE (NcmFitNLOpt, ncm_fit_nlopt, NCM, FIT_NLOPT, NcmFit)

NcmFit *ncm_fit_nlopt_new (NcmLikelihood * lh, NcmMSet * mset, NcmFitGradType gtype, NcmFitNloptAlgorithm algo);
NcmFit *ncm_fit_nlopt_local_new (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, NcmFitNloptAlgorithm algo, NcmFitNloptAlgorithm local_algo);
void ncm_fit_nlopt_set_algo (NcmFitNLOpt *fit_nlopt, NcmFitNloptAlgorithm algo);
void ncm_fit_nlopt_set_local_algo (NcmFitNLOpt *fit_nlopt, NcmFitNloptAlgorithm algo);

NcmFit *ncm_fit_nlopt_new_by_name (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, gchar *algo_name);
NcmFit *ncm_fit_nlopt_new_default (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype);

G_END_DECLS

#endif /* NUMCOSMO_HAVE_NLOPT */
#endif /* _NCM_FIT_NLOPT_H_ */

