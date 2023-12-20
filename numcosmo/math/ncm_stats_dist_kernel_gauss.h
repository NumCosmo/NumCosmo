/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_stats_dist_kernel_gauss.h
 *
 *  Wed November 07 17:41:38 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_kernel_gauss.h
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_STATS_DIST_KERNEL_GAUSS_H_
#define _NCM_STATS_DIST_KERNEL_GAUSS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_stats_dist_kernel.h>
#include <numcosmo/math/ncm_vector.h>

G_BEGIN_DECLS

#define NCM_TYPE_STATS_DIST_KERNEL_GAUSS (ncm_stats_dist_kernel_gauss_get_type ())

G_DECLARE_FINAL_TYPE (NcmStatsDistKernelGauss, ncm_stats_dist_kernel_gauss, NCM, STATS_DIST_KERNEL_GAUSS, NcmStatsDistKernel)

NcmStatsDistKernelGauss *ncm_stats_dist_kernel_gauss_new (const guint dim);

NcmStatsDistKernelGauss *ncm_stats_dist_kernel_gauss_ref (NcmStatsDistKernelGauss *sdkg);
void ncm_stats_dist_kernel_gauss_free (NcmStatsDistKernelGauss *sdkg);
void ncm_stats_dist_kernel_gauss_clear (NcmStatsDistKernelGauss **sdkg);

G_END_DECLS

#endif /* _NCM_STATS_DIST_KERNEL_GAUSS_H_ */

