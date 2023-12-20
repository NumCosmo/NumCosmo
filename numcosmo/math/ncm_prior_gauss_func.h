/***************************************************************************
 *            ncm_prior_gauss_func.h
 *
 *  Wed August 03 16:26:48 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_prior_gauss_func.h
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

#ifndef _NCM_PRIOR_GAUSS_FUNC_H_
#define _NCM_PRIOR_GAUSS_FUNC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_prior_gauss.h>

G_BEGIN_DECLS

#define NCM_TYPE_PRIOR_GAUSS_FUNC (ncm_prior_gauss_func_get_type ())

G_DECLARE_FINAL_TYPE(NcmPriorGaussFunc, ncm_prior_gauss_func, NCM, PRIOR_GAUSS_FUNC, NcmPriorGauss)

NcmPriorGaussFunc *ncm_prior_gauss_func_new (NcmMSetFunc *mean_func, gdouble mu, gdouble sigma, gdouble var);
NcmPriorGaussFunc *ncm_prior_gauss_func_ref (NcmPriorGaussFunc *pgf);

void ncm_prior_gauss_func_free (NcmPriorGaussFunc *pgf);
void ncm_prior_gauss_func_clear (NcmPriorGaussFunc **pgf);

G_END_DECLS

#endif /* _NCM_PRIOR_GAUSS_FUNC_H_ */
