/***************************************************************************
 *            ncm_stats_rows.h
 *
 *  Fri Oct 24 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_rows.h
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_STATS_ROWS_H_
#define _NCM_STATS_ROWS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/ncm/algebra/ncm_matrix.h>

G_BEGIN_DECLS

#define NCM_TYPE_STATS_ROWS (ncm_stats_rows_get_type ())

G_DECLARE_FINAL_TYPE (NcmStatsRows, ncm_stats_rows, NCM, STATS_ROWS, GObject)

NcmStatsRows *ncm_stats_rows_new (guint ncols);
NcmStatsRows *ncm_stats_rows_ref (NcmStatsRows *rows);
void ncm_stats_rows_free (NcmStatsRows *rows);
void ncm_stats_rows_clear (NcmStatsRows **rows);

void ncm_stats_rows_set_rows (NcmStatsRows *rows, GPtrArray *sample);
void ncm_stats_rows_set_nrows (NcmStatsRows *rows, guint nrows);
void ncm_stats_rows_reset (NcmStatsRows *rows);
guint ncm_stats_rows_get_nrows (NcmStatsRows *rows);
guint ncm_stats_rows_get_ncols (NcmStatsRows *rows);
NcmVector *ncm_stats_rows_peek_row (NcmStatsRows *rows, guint i);

NcmVector *ncm_stats_rows_max_ess_time (NcmStatsRows *rows, const guint ntests, gint *bindex, guint *wp, guint *wp_order, gdouble *wp_ess);
NcmVector *ncm_stats_rows_heidel_diag (NcmStatsRows *rows, const guint ntests, const gdouble pvalue, gint *bindex, guint *wp, guint *wp_order, gdouble *wp_pvalue);
NcmVector *ncm_stats_rows_visual_heidel_diag (NcmStatsRows *rows, const guint p, const guint fi, gdouble *mean, gdouble *var);
gdouble ncm_stats_rows_estimate_const_break (NcmStatsRows *rows, guint p);
void ncm_stats_rows_cov_robust_diag (NcmStatsRows *rows, NcmMatrix *cov);
void ncm_stats_rows_cov_robust_ogk (NcmStatsRows *rows, NcmMatrix *cov);

G_END_DECLS

#endif /* _NCM_STATS_ROWS_H_ */

