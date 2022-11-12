/***************************************************************************
 *            ncm_prior_gauss_param.h
 *
 *  Wed August 03 10:55:33 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_prior_gauss_param.h
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

#ifndef _NCM_PRIOR_GAUSS_PARAM_H_
#define _NCM_PRIOR_GAUSS_PARAM_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_prior_gauss.h>

G_BEGIN_DECLS

#define NCM_TYPE_PRIOR_GAUSS_PARAM             (ncm_prior_gauss_param_get_type ())
#define NCM_PRIOR_GAUSS_PARAM(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_PRIOR_GAUSS_PARAM, NcmPriorGaussParam))
#define NCM_PRIOR_GAUSS_PARAM_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_PRIOR_GAUSS_PARAM, NcmPriorGaussParamClass))
#define NCM_IS_PRIOR_GAUSS_PARAM(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_PRIOR_GAUSS_PARAM))
#define NCM_IS_PRIOR_GAUSS_PARAM_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_PRIOR_GAUSS_PARAM))
#define NCM_PRIOR_GAUSS_PARAM_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_PRIOR_GAUSS_PARAM, NcmPriorGaussParamClass))

typedef struct _NcmPriorGaussParamClass NcmPriorGaussParamClass;
typedef struct _NcmPriorGaussParam NcmPriorGaussParam;

struct _NcmPriorGaussParamClass
{
  /*< private >*/
  NcmPriorGaussClass parent_class;
};

struct _NcmPriorGaussParam
{
  /*< private >*/
  NcmPriorGauss parent_instance;
  NcmModelID mid;
  guint pid;
};

GType ncm_prior_gauss_param_get_type (void) G_GNUC_CONST;

NcmPriorGaussParam *ncm_prior_gauss_param_new (NcmModelID mid, guint pid, gdouble mu, gdouble sigma);
NcmPriorGaussParam *ncm_prior_gauss_param_new_pindex (const NcmMSetPIndex *pi, gdouble mu, gdouble sigma);
NcmPriorGaussParam *ncm_prior_gauss_param_new_name (NcmMSet *mset, const gchar *name, gdouble mu, gdouble sigma);

NcmPriorGaussParam *ncm_prior_gauss_param_ref (NcmPriorGaussParam *pgp);

void ncm_prior_gauss_param_free (NcmPriorGaussParam *pgp);
void ncm_prior_gauss_param_clear (NcmPriorGaussParam **pgp);

G_END_DECLS

#endif /* _NCM_PRIOR_GAUSS_PARAM_H_ */
