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
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_prior_gauss.h>

G_BEGIN_DECLS

#define NCM_TYPE_PRIOR_GAUSS_PARAM (ncm_prior_gauss_param_get_type ())

G_DECLARE_FINAL_TYPE (NcmPriorGaussParam, ncm_prior_gauss_param, NCM, PRIOR_GAUSS_PARAM, NcmPriorGauss)

NcmPriorGaussParam *ncm_prior_gauss_param_new (NcmModel * model, guint pid, gdouble mu, gdouble sigma);
NcmPriorGaussParam *ncm_prior_gauss_param_new_name (const gchar *name, gdouble mu, gdouble sigma, GError **error);

NcmPriorGaussParam *ncm_prior_gauss_param_ref (NcmPriorGaussParam *pgp);

void ncm_prior_gauss_param_free (NcmPriorGaussParam *pgp);
void ncm_prior_gauss_param_clear (NcmPriorGaussParam **pgp);

void ncm_prior_gauss_param_set_model_ns (NcmPriorGaussParam *pgp, const gchar *model_ns, GError **error);
void ncm_prior_gauss_param_set_stack_pos (NcmPriorGaussParam *pgp, guint stack_pos);
void ncm_prior_gauss_param_set_param_name (NcmPriorGaussParam *pgp, const gchar *param_name);

const gchar *ncm_prior_gauss_param_peek_model_ns (NcmPriorGaussParam *pgp);
const gchar *ncm_prior_gauss_param_peek_param_name (NcmPriorGaussParam *pgp);
guint ncm_prior_gauss_param_get_stack_pos (NcmPriorGaussParam *pgp);

G_END_DECLS

#endif /* _NCM_PRIOR_GAUSS_PARAM_H_ */

