/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_stats_dist_kernel_st.h
 *
 *  Wed November 07 17:41:38 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_kernel_st.h
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

#ifndef _NCM_STATS_DIST_KERNEL_ST_H_
#define _NCM_STATS_DIST_KERNEL_ST_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_stats_dist_kernel.h>
#include <numcosmo/math/ncm_vector.h>

G_BEGIN_DECLS

#define NCM_TYPE_STATS_DIST_KERNEL_ST             (ncm_stats_dist_kernel_st_get_type ())
#define NCM_STATS_DIST_KERNEL_ST(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_STATS_DIST_KERNEL_ST,  NcmStatsDistKernelST))
#define NCM_STATS_DIST_KERNEL_ST_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_STATS_DIST_KERNEL_ST,  NcmStatsDistKernelSTClass))
#define NCM_IS_STATS_DIST_KERNEL_ST(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_STATS_DIST_KERNEL_ST))
#define NCM_IS_STATS_DIST_KERNEL_ST_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_STATS_DIST_KERNEL_ST))
#define NCM_STATS_DIST_KERNEL_ST_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_STATS_DIST_KERNEL_ST,  NcmStatsDistKernelSTClass))

typedef struct _NcmStatsDistKernelSTClass NcmStatsDistKernelSTClass;
typedef struct _NcmStatsDistKernelST NcmStatsDistKernelST;
typedef struct _NcmStatsDistKernelSTPrivate NcmStatsDistKernelSTPrivate;


struct _NcmStatsDistKernelSTClass
{
  NcmStatsDistKernelClass parent_class;
};

struct _NcmStatsDistKernelST
{
  NcmStatsDistKernel parent_instance;
  
  NcmStatsDistKernelSTPrivate *priv;
};

GType ncm_stats_dist_kernel_st_get_type (void) G_GNUC_CONST;

NcmStatsDistKernelST *ncm_stats_dist_kernel_st_new (const guint dim, const gdouble nu);

NcmStatsDistKernelST *ncm_stats_dist_kernel_st_ref (NcmStatsDistKernelST *sdkst);
void ncm_stats_dist_kernel_st_free (NcmStatsDistKernelST *sdkst);
void ncm_stats_dist_kernel_st_clear (NcmStatsDistKernelST **sdkst);

void ncm_stats_dist_kernel_st_set_nu (NcmStatsDistKernelST *sdkst, const gdouble nu);
gdouble ncm_stats_dist_kernel_st_get_nu (NcmStatsDistKernelST *sdkst);

G_END_DECLS

#endif /* _NCM_STATS_DIST_KERNEL_ST_H_ */

