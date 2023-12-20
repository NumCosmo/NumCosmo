/***************************************************************************
 *            ncm_mset_trans_kern_gauss.h
 *
 *  Wed September 03 14:55:28 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mset_trans_kern_gauss.h
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_MSET_TRANS_KERN_GAUSS_H_
#define _NCM_MSET_TRANS_KERN_GAUSS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mset_trans_kern.h>

G_BEGIN_DECLS

#define NCM_TYPE_MSET_TRANS_KERN_GAUSS (ncm_mset_trans_kern_gauss_get_type ())

G_DECLARE_FINAL_TYPE (NcmMSetTransKernGauss, ncm_mset_trans_kern_gauss, NCM, MSET_TRANS_KERN_GAUSS, NcmMSetTransKern)

NcmMSetTransKernGauss *ncm_mset_trans_kern_gauss_new (guint len);

void ncm_mset_trans_kern_gauss_set_size (NcmMSetTransKernGauss *tkerng, guint len);
guint ncm_mset_trans_kern_gauss_get_size (NcmMSetTransKernGauss *tkerng);

void ncm_mset_trans_kern_gauss_set_cov (NcmMSetTransKernGauss *tkerng, const NcmMatrix *cov);
void ncm_mset_trans_kern_gauss_set_cov_variant (NcmMSetTransKernGauss *tkerng, GVariant *cov);
NcmMatrix *ncm_mset_trans_kern_gauss_get_cov (NcmMSetTransKernGauss *tkerng);

void ncm_mset_trans_kern_gauss_set_cov_from_scale (NcmMSetTransKernGauss *tkerng);
void ncm_mset_trans_kern_gauss_set_cov_from_rescale (NcmMSetTransKernGauss *tkerng, const gdouble epsilon);

G_END_DECLS

#endif /* _NCM_MSET_TRANS_KERN_GAUSS_H_ */

