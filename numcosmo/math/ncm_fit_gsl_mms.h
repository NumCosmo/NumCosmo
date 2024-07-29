/***************************************************************************
 *            ncm_fit_gsl_mms.h
 *
 *  Sat Aug 16 19:57:28 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
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

#ifndef _NCM_FIT_GSL_MMS_H_
#define _NCM_FIT_GSL_MMS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_fit.h>

G_BEGIN_DECLS

#define NCM_TYPE_FIT_GSL_MMS (ncm_fit_gsl_mms_get_type ())

G_DECLARE_FINAL_TYPE (NcmFitGSLMMS, ncm_fit_gsl_mms, NCM, FIT_GSL_MMS, NcmFit)

/**
 * NcmFitGSLMMSAlgos:
 * @NCM_FIT_GSL_MMS_NMSIMPLEX2: Nelder-Mead Simplex algorithm - faster implementation
 * @NCM_FIT_GSL_MMS_NMSIMPLEX: Nelson-Mead Simplex algorithm - original implementation
 * @NCM_FIT_GSL_MMS_NMSIMPLEX2RAND: Nelder-Mead Simplex algorithm implementation with randomization
 *
 * GSL Multidimensional minimization algorithms without derivatives
 *
 */
typedef enum _NcmFitGSLMMSAlgos
{
  NCM_FIT_GSL_MMS_NMSIMPLEX2 = 0,
  NCM_FIT_GSL_MMS_NMSIMPLEX,
  NCM_FIT_GSL_MMS_NMSIMPLEX2RAND,
  /* < private > */
  NCM_FIT_GSL_MMS_NUM_ALGOS, /*< skip >*/
} NcmFitGSLMMSAlgos;

NcmFit *ncm_fit_gsl_mms_new (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, NcmFitGSLMMSAlgos algo);
NcmFit *ncm_fit_gsl_mms_new_default (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype);
NcmFit *ncm_fit_gsl_mms_new_by_name (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, gchar *algo_name);
void ncm_fit_gsl_mms_set_algo (NcmFitGSLMMS *fit_gsl_mms, NcmFitGSLMMSAlgos algo);

G_END_DECLS

#endif /* _NCM_FIT_GSL_MMS_H_ */

