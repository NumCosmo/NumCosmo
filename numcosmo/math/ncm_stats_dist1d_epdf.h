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

struct _NcmStatsDist1dEPDF
{
  /*< private >*/
  NcmStatsDist1d parent_instance;
  NcmStatsVec *obs_stats;
  guint max_obs;
  guint min_knots;
  guint max_knots;
  gdouble sd_scale;
  gdouble sd_min_scale;
  gdouble outliers_threshold;
  gdouble sd;
  guint n_obs;
  guint np_obs;
  GArray *obs;
  gdouble min;
  gdouble max;
  gboolean list_sorted;
};

GType ncm_stats_dist1d_epdf_get_type (void) G_GNUC_CONST;

NcmStatsDist1dEPDF *ncm_stats_dist1d_epdf_new (guint max_obs, gdouble sd_scale, gdouble sd_min_scale);

void ncm_stats_dist1d_epdf_add_obs (NcmStatsDist1dEPDF *epdf1d, gdouble x);
void ncm_stats_dist1d_epdf_reset (NcmStatsDist1dEPDF *epdf1d);

gdouble ncm_stats_dist1d_epdf_get_obs_mean (NcmStatsDist1dEPDF *epdf1d);

G_END_DECLS

#endif /* _NCM_STATS_DIST1D_EPDF_H_ */
