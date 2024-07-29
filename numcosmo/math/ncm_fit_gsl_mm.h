/***************************************************************************
 *            ncm_fit_gsl_mm.h
 *
 *  Sat Aug 16 19:56:29 2008
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

#ifndef _NCM_FIT_GSL_MM_H_
#define _NCM_FIT_GSL_MM_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_fit.h>

G_BEGIN_DECLS

#define NCM_TYPE_FIT_GSL_MM (ncm_fit_gsl_mm_get_type ())

G_DECLARE_FINAL_TYPE (NcmFitGSLMM, ncm_fit_gsl_mm, NCM, FIT_GSL_MM, NcmFit)

/**
 * NcmFitGSLMMAlgos:
 * @NCM_FIT_GSL_MM_CONJUGATE_FR: Fletcher-Reeves conjugate gradient algorithm
 * @NCM_FIT_GSL_MM_CONJUGATE_PR: Polak-Ribiere conjugate gradient algorithm
 * @NCM_FIT_GSL_MM_VECTOR_BFGS: Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm
 * @NCM_FIT_GSL_MM_VECTOR_BFGS2: More efficient variation of BFGS algorithm
 * @NCM_FIT_GSL_MM_STEEPEST_DESCENT: Steepest descent algorithm
 *
 * GSL Multidimensional minimization algorithms
 *
 */
typedef enum _NcmFitGSLMMAlgos
{
  NCM_FIT_GSL_MM_CONJUGATE_FR = 0,
  NCM_FIT_GSL_MM_CONJUGATE_PR,
  NCM_FIT_GSL_MM_VECTOR_BFGS,
  NCM_FIT_GSL_MM_VECTOR_BFGS2,
  NCM_FIT_GSL_MM_STEEPEST_DESCENT,
  /* < private > */
  NCM_FIT_GSL_MM_NUM_ALGOS, /*< skip >*/
} NcmFitGSLMMAlgos;


NcmFit *ncm_fit_gsl_mm_new (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, NcmFitGSLMMAlgos algo);
NcmFit *ncm_fit_gsl_mm_new_default (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype);
NcmFit *ncm_fit_gsl_mm_new_by_name (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, gchar *algo_name);
void ncm_fit_gsl_mm_set_algo (NcmFitGSLMM *fit_gsl_mm, NcmFitGSLMMAlgos algo);


G_END_DECLS

#endif /* _NCM_FIT_GSL_MM_H_ */

