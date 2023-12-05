/***************************************************************************
 *            ncm_fit_gsl_ls.h
 *
 *  Sat Aug 16 19:58:39 2008
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

#ifndef _NCM_FIT_GSL_LS_H_
#define _NCM_FIT_GSL_LS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_fit.h>

G_BEGIN_DECLS

#define NCM_TYPE_FIT_GSL_LS (ncm_fit_gsl_ls_get_type ())

G_DECLARE_FINAL_TYPE (NcmFitGSLLS, ncm_fit_gsl_ls, NCM, FIT_GSL_LS, NcmFit)

NcmFit *ncm_fit_gsl_ls_new (NcmLikelihood * lh, NcmMSet * mset, NcmFitGradType gtype);

G_END_DECLS

#endif /* _NCM_FIT_GSL_LS_H_ */

