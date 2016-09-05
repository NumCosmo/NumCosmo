/***************************************************************************
 *            ncm_fit_gsl_ls.h
 *
 *  Sat Aug 16 19:58:39 2008
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

#ifndef _NCM_FIT_GSL_LS_H_
#define _NCM_FIT_GSL_LS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_fit.h>

G_BEGIN_DECLS

#define NCM_TYPE_FIT_GSL_LS             (ncm_fit_gsl_ls_get_type ())
#define NCM_FIT_GSL_LS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_FIT_GSL_LS, NcmFitGSLLS))
#define NCM_FIT_GSL_LS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_FIT_GSL_LS, NcmFitGSLLSClass))
#define NCM_IS_FIT_GSL_LS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_FIT_GSL_LS))
#define NCM_IS_FIT_GSL_LS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_FIT_GSL_LS))
#define NCM_FIT_GSL_LS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_FIT_GSL_LS, NcmFitGSLLSClass))

typedef struct _NcmFitGSLLSClass NcmFitGSLLSClass;
typedef struct _NcmFitGSLLS NcmFitGSLLS;

struct _NcmFitGSLLSClass
{
  /*< private >*/
  NcmFitClass parent_class;
};

struct _NcmFitGSLLS
{
  /*< private >*/
  NcmFit parent_instance;
  gsl_multifit_fdfsolver *ls;
  gsl_multifit_function_fdf f;
  const gsl_multifit_fdfsolver_type *T;
};

GType ncm_fit_gsl_ls_get_type (void) G_GNUC_CONST;

NcmFit *ncm_fit_gsl_ls_new (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype);

G_END_DECLS

#endif /* _NCM_FIT_GSL_LS_H_ */

