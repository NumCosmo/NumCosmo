/***************************************************************************
 *            ncm_stats_dist1d_epdf.h
 *
 *  Sat March 14 19:31:53 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_stats_dist1d_epdf.h
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#define NCM_TYPE_STATS_DIST1D_EPDF             (ncm_stats_dist1d_epdf_get_type ())
#define NCM_STATS_DIST1D_EPDF(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_STATS_DIST1D_EPDF, NcmStatsDist1dEPDF))
#define NCM_STATS_DIST1D_EPDF_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_STATS_DIST1D_EPDF, NcmStatsDist1dEPDFClass))
#define NCM_IS_STATS_DIST1D_EPDF(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_STATS_DIST1D_EPDF))
#define NCM_IS_STATS_DIST1D_EPDF_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_STATS_DIST1D_EPDF))
#define NCM_STATS_DIST1D_EPDF_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_STATS_DIST1D_EPDF, NcmStatsDist1dEPDFClass))

typedef struct _NcmStatsDist1dEPDFClass NcmStatsDist1dEPDFClass;
typedef struct _NcmStatsDist1dEPDF NcmStatsDist1dEPDF;

struct _NcmStatsDist1dEPDFClass
{
  /*< private >*/
  NcmStatsDist1dClass parent_class;
};

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
  NCM_STATS_DIST1D_EPDF_BW_AUTO, /*< private >*/
  NCM_STATS_DIST1D_EPDF_BW_LEN,  /*< skip >*/
} NcmStatsDist1dEPDFBw;

struct _NcmStatsDist1dEPDF
{
  /*< private >*/
  NcmStatsDist1d parent_instance;
  NcmStatsVec *obs_stats;
  guint max_obs;
  NcmStatsDist1dEPDFBw bw;
  gdouble h_fixed;
  gdouble sd_min_scale;
  gdouble outliers_threshold;
  gdouble h;
  guint n_obs;
  guint np_obs;
  GArray *obs;
  gdouble min;
  gdouble max;
  gboolean list_sorted;
  guint fftsize;
  NcmVector *Iv;
  NcmVector *p_data;
  NcmVector *p_tilde;
  NcmVector *p_tilde2;
  NcmVector *p_est;
  NcmVector *xv;
  NcmVector *pv;
  gpointer fft_data_to_tilde; 
  gpointer fft_tilde_to_est;
  NcmSpline *p_spline;
  gboolean bw_set;
};

GType ncm_stats_dist1d_epdf_get_type (void) G_GNUC_CONST;

NcmStatsDist1dEPDF *ncm_stats_dist1d_epdf_new_full (guint max_obs, NcmStatsDist1dEPDFBw bw, gdouble h_fixed, gdouble sd_min_scale);
NcmStatsDist1dEPDF *ncm_stats_dist1d_epdf_new (gdouble sd_min_scale);

void ncm_stats_dist1d_epdf_add_obs (NcmStatsDist1dEPDF *epdf1d, gdouble x);
void ncm_stats_dist1d_epdf_reset (NcmStatsDist1dEPDF *epdf1d);

gdouble ncm_stats_dist1d_epdf_get_obs_mean (NcmStatsDist1dEPDF *epdf1d);

G_END_DECLS

#endif /* _NCM_STATS_DIST1D_EPDF_H_ */
