/***************************************************************************
 *            ncm_prior_gauss.h
 *
 *  Wed August 03 10:19:32 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_prior_gauss.h
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

#ifndef _NCM_PRIOR_GAUSS_H_
#define _NCM_PRIOR_GAUSS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_prior.h>

G_BEGIN_DECLS

#define NCM_TYPE_PRIOR_GAUSS             (ncm_prior_gauss_get_type ())
#define NCM_PRIOR_GAUSS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_PRIOR_GAUSS, NcmPriorGauss))
#define NCM_PRIOR_GAUSS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_PRIOR_GAUSS, NcmPriorGaussClass))
#define NCM_IS_PRIOR_GAUSS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_PRIOR_GAUSS))
#define NCM_IS_PRIOR_GAUSS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_PRIOR_GAUSS))
#define NCM_PRIOR_GAUSS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_PRIOR_GAUSS, NcmPriorGaussClass))

typedef struct _NcmPriorGaussClass NcmPriorGaussClass;
typedef struct _NcmPriorGauss NcmPriorGauss;

typedef gdouble (*NcmPriorGaussMean) (NcmPriorGauss *pg, NcmMSet *mset);

struct _NcmPriorGaussClass
{
  /*< private >*/
  NcmPriorClass parent_class;
  NcmPriorGaussMean mean;
};

struct _NcmPriorGauss
{
  /*< private >*/
  NcmPrior parent_instance;
  gdouble mu;
  gdouble sigma;
  gdouble var;
};

GType ncm_prior_gauss_get_type (void) G_GNUC_CONST;

NcmPriorGauss *ncm_prior_gauss_ref (NcmPriorGauss *pg);
void ncm_prior_gauss_free (NcmPriorGauss *pg);
void ncm_prior_gauss_clear (NcmPriorGauss **pg);

G_END_DECLS

#endif /* _NCM_PRIOR_GAUSS_H_ */
