/***************************************************************************
 *            ncm_fit_levmar.h
 *
 *  Wed Feb 24 21:19:48 2010
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

#ifndef _NCM_FIT_LEVMAR_H_
#define _NCM_FIT_LEVMAR_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_fit.h>

G_BEGIN_DECLS

#define NCM_TYPE_FIT_LEVMAR (ncm_fit_levmar_get_type ())

G_DECLARE_FINAL_TYPE (NcmFitLevmar, ncm_fit_levmar, NCM, FIT_LEVMAR, NcmFit)

/**
 * NcmFitLevmarAlgos:
 * @NCM_FIT_LEVMAR_DER: with external derivatives.
 * @NCM_FIT_LEVMAR_DIF: with internal derivatives (inside levmar).
 * @NCM_FIT_LEVMAR_BC_DER: with box constraints and external derivatives.
 * @NCM_FIT_LEVMAR_BC_DIF: with box constraints and internal derivatives (inside levmar).
 *
 * Levmar algorithms.
 *
 */
typedef enum _NcmFitLevmarAlgos
{
  NCM_FIT_LEVMAR_DER = 0,
  NCM_FIT_LEVMAR_DIF,
  NCM_FIT_LEVMAR_BC_DER,
  NCM_FIT_LEVMAR_BC_DIF,
  /* < private > */
  NCM_FIT_LEVMAR_NUM_ALGOS, /*< skip >*/
} NcmFitLevmarAlgos;


NcmFit *ncm_fit_levmar_new (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, NcmFitLevmarAlgos algo);
NcmFit *ncm_fit_levmar_new_default (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype);
NcmFit *ncm_fit_levmar_new_by_name (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, gchar *algo_name);
void ncm_fit_levmar_set_algo (NcmFitLevmar *fit_levmar, NcmFitLevmarAlgos algo);

G_END_DECLS

#endif /* _NCM_FIT_LEVMAR_H_ */

