/***************************************************************************
 *            ncm_fit_gsl_mm.h
 *
 *  Sat Aug 16 19:56:29 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
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

#ifndef _NCM_FIT_GSL_MM_H_
#define _NCM_FIT_GSL_MM_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_fit.h>

G_BEGIN_DECLS

#define NCM_TYPE_FIT_GSL_MM             (ncm_fit_gsl_mm_get_type ())
#define NCM_FIT_GSL_MM(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_FIT_GSL_MM, NcmFitGSLMM))
#define NCM_FIT_GSL_MM_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_FIT_GSL_MM, NcmFitGSLMMClass))
#define NCM_IS_FIT_GSL_MM(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_FIT_GSL_MM))
#define NCM_IS_FIT_GSL_MM_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_FIT_GSL_MM))
#define NCM_FIT_GSL_MM_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_FIT_GSL_MM, NcmFitGSLMMClass))

typedef struct _NcmFitGSLMMClass NcmFitGSLMMClass;
typedef struct _NcmFitGSLMM NcmFitGSLMM;

/**
 * NcmFitGSLMMAlgos:
 * @NCM_FIT_GSL_MM_CONJUGATE_FR: FIXME
 * @NCM_FIT_GSL_MM_CONJUGATE_PR: FIXME
 * @NCM_FIT_GSL_MM_VECTOR_BFGS: FIXME
 * @NCM_FIT_GSL_MM_VECTOR_BFGS2: FIXME
 * @NCM_FIT_GSL_MM_STEEPEST_DESCENT: FIXME
 *
 * FIXME
 */
typedef enum _NcmFitGSLMMAlgos
{
  NCM_FIT_GSL_MM_CONJUGATE_FR = 0,
  NCM_FIT_GSL_MM_CONJUGATE_PR,
  NCM_FIT_GSL_MM_VECTOR_BFGS,
  NCM_FIT_GSL_MM_VECTOR_BFGS2,
  NCM_FIT_GSL_MM_STEEPEST_DESCENT,  /*< private >*/
  NCM_FIT_GSL_MM_NUM_ALGOS,         /*< skip >*/
} NcmFitGSLMMAlgos;

struct _NcmFitGSLMM
{
  /*< private >*/
  NcmFit parent_instance;
  gsl_multimin_fdfminimizer *mm;
  gsl_multimin_function_fdf f;
  NcmFitGSLMMAlgos algo;
  gchar *desc;
  gdouble err_a;
  gdouble err_b;
};

struct _NcmFitGSLMMClass
{
  /*< private >*/
  NcmFitClass parent_class;
};

GType ncm_fit_gsl_mm_get_type (void) G_GNUC_CONST;

NcmFit *ncm_fit_gsl_mm_new (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, NcmFitGSLMMAlgos algo);
NcmFit *ncm_fit_gsl_mm_new_default (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype);
NcmFit *ncm_fit_gsl_mm_new_by_name (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, gchar *algo_name);
void ncm_fit_gsl_mm_set_algo (NcmFitGSLMM *fit_gsl_mm, NcmFitGSLMMAlgos algo);


G_END_DECLS

#endif /* _NCM_FIT_GSL_MM_H_ */

