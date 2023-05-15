/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_stats_dist_kde_private.h
 *
 *  Thu July 22 15:12:38 2021
 *  Copyright  2021  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_kde_private.h
 * Copyright (C) 2021 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_STATS_DIST_KDE_PRIVATE_H_
#define _NCM_STATS_DIST_KDE_PRIVATE_H_

#include <glib.h>
#include "math/ncm_stats_dist_kde.h"
#include "math/ncm_memory_pool.h"

G_BEGIN_DECLS

struct _NcmStatsDistKDEPrivate
{
  /*< private >*/
  NcmStatsVec *sample;
  NcmStatsDistKDECovType cov_type;
  NcmMatrix *cov_fixed;
  NcmMatrix *cov_decomp;
  NcmMatrix *sample_matrix;
  NcmMatrix *invUsample_matrix;
  GPtrArray *invUsample_array;
  gdouble kernel_lnnorm;
  guint nearPD_maxiter;
  NcmMemoryPool *mp_eval_vars;
  guint mp_eval_vars_len;
};

G_END_DECLS

#endif /* _NCM_STATS_DIST_KDE_PRIVATE_H_ */

