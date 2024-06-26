/***************************************************************************
 *            ncm_stats_dist1d_epdf.h
 *
 *  Sat March 14 19:31:53 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_dist1d_epdf.h
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_STATS_DIST1D_EPDF_H_
#define _NCM_STATS_DIST1D_EPDF_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_stats_dist1d.h>
#include <numcosmo/math/ncm_stats_vec.h>

G_BEGIN_DECLS

#define NCM_TYPE_STATS_DIST1D_EPDF (ncm_stats_dist1d_epdf_get_type ())

G_DECLARE_FINAL_TYPE (NcmStatsDist1dEPDF, ncm_stats_dist1d_epdf, NCM, STATS_DIST1D_EPDF, NcmStatsDist1d)

/**
 * NcmStatsDist1dEPDFBw:
 * @NCM_STATS_DIST1D_EPDF_BW_FIXED: Uses the given value of bandwidth.
 * @NCM_STATS_DIST1D_EPDF_BW_RoT: Uses the Silverman's rule of thumb to determine the bandwidth.
 * @NCM_STATS_DIST1D_EPDF_BW_AUTO: Uses Botev's et al method to automatically determine the best bandwidth.
 *
 * Gaussian kernel bandwidth type.
 *
 */
typedef enum _NcmStatsDist1dEPDFBw
{
  NCM_STATS_DIST1D_EPDF_BW_FIXED = 0,
  NCM_STATS_DIST1D_EPDF_BW_RoT,
  NCM_STATS_DIST1D_EPDF_BW_AUTO,
  /* < private > */
  NCM_STATS_DIST1D_EPDF_BW_LEN, /*< skip >*/
} NcmStatsDist1dEPDFBw;

NcmStatsDist1dEPDF *ncm_stats_dist1d_epdf_new_full (guint max_obs, NcmStatsDist1dEPDFBw bw, gdouble h_fixed, gdouble sd_min_scale);
NcmStatsDist1dEPDF *ncm_stats_dist1d_epdf_new (gdouble sd_min_scale);

NcmStatsDist1dEPDF *ncm_stats_dist1d_epdf_ref (NcmStatsDist1dEPDF *epdf1d);
void ncm_stats_dist1d_epdf_free (NcmStatsDist1dEPDF *epdf1d);
void ncm_stats_dist1d_epdf_clear (NcmStatsDist1dEPDF **epdf1d);

void ncm_stats_dist1d_epdf_set_bw_type (NcmStatsDist1dEPDF *epdf1d, NcmStatsDist1dEPDFBw bw);
NcmStatsDist1dEPDFBw ncm_stats_dist1d_epdf_get_bw_type (NcmStatsDist1dEPDF *epdf1d);

void ncm_stats_dist1d_epdf_set_h_fixed (NcmStatsDist1dEPDF *epdf1d, gdouble h_fixed);
gdouble ncm_stats_dist1d_epdf_get_h_fixed (NcmStatsDist1dEPDF *epdf1d);

void ncm_stats_dist1d_epdf_add_obs_weight (NcmStatsDist1dEPDF *epdf1d, const gdouble x, const gdouble w);
void ncm_stats_dist1d_epdf_add_obs (NcmStatsDist1dEPDF *epdf1d, const gdouble x);
void ncm_stats_dist1d_epdf_reset (NcmStatsDist1dEPDF *epdf1d);

void ncm_stats_dist1d_epdf_set_min (NcmStatsDist1dEPDF *epdf1d, const gdouble min);
void ncm_stats_dist1d_epdf_set_max (NcmStatsDist1dEPDF *epdf1d, const gdouble max);

gdouble ncm_stats_dist1d_epdf_get_obs_mean (NcmStatsDist1dEPDF *epdf1d);

G_END_DECLS

#endif /* _NCM_STATS_DIST1D_EPDF_H_ */

