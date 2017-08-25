/***************************************************************************
 *            ncm_stats_dist1d.h
 *
 *  Thu February 12 15:37:22 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_stats_dist1d.h
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

#ifndef _NCM_STATS_DIST1D_H_
#define _NCM_STATS_DIST1D_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_ode_spline.h>
#include <numcosmo/math/ncm_rng.h>
#include <gsl/gsl_min.h>

G_BEGIN_DECLS

#define NCM_TYPE_STATS_DIST1D             (ncm_stats_dist1d_get_type ())
#define NCM_STATS_DIST1D(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_STATS_DIST1D, NcmStatsDist1d))
#define NCM_STATS_DIST1D_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_STATS_DIST1D, NcmStatsDist1dClass))
#define NCM_IS_STATS_DIST1D(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_STATS_DIST1D))
#define NCM_IS_STATS_DIST1D_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_STATS_DIST1D))
#define NCM_STATS_DIST1D_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_STATS_DIST1D, NcmStatsDist1dClass))

typedef struct _NcmStatsDist1dClass NcmStatsDist1dClass;
typedef struct _NcmStatsDist1d NcmStatsDist1d;

struct _NcmStatsDist1dClass
{
  /*< private >*/
  GObjectClass parent_class;
  gdouble (*p) (NcmStatsDist1d *sd1, gdouble x);
  gdouble (*m2lnp) (NcmStatsDist1d *sd1, gdouble x);
  void (*prepare) (NcmStatsDist1d *sd1);
};

struct _NcmStatsDist1d
{
  /*< private >*/
  GObject parent_instance;
  gdouble xi;
  gdouble xf;
  gdouble norma;
  gdouble reltol;
  gdouble abstol;
  gdouble max_prob;
  NcmOdeSpline *inv_cdf;
  NcmOdeSpline *pdf;
  gsl_min_fminimizer *fmin;
};

GType ncm_stats_dist1d_get_type (void) G_GNUC_CONST;

NcmStatsDist1d *ncm_stats_dist1d_ref (NcmStatsDist1d *sd1);
void ncm_stats_dist1d_free (NcmStatsDist1d *sd1);
void ncm_stats_dist1d_clear (NcmStatsDist1d **sd1);

void ncm_stats_dist1d_prepare (NcmStatsDist1d *sd1);

gdouble ncm_stats_dist1d_get_xi (NcmStatsDist1d *sd1);
gdouble ncm_stats_dist1d_get_xf (NcmStatsDist1d *sd1);

gdouble ncm_stats_dist1d_eval_p (NcmStatsDist1d *sd1, gdouble x);
gdouble ncm_stats_dist1d_eval_m2lnp (NcmStatsDist1d *sd1, gdouble x);
gdouble ncm_stats_dist1d_eval_pdf (NcmStatsDist1d *sd1, gdouble x);
gdouble ncm_stats_dist1d_eval_norma (NcmStatsDist1d *sd1);
gdouble ncm_stats_dist1d_eval_inv_pdf (NcmStatsDist1d *sd1, gdouble u);
gdouble ncm_stats_dist1d_eval_inv_pdf_tail (NcmStatsDist1d *sd1, gdouble v);
gdouble ncm_stats_dist1d_gen (NcmStatsDist1d *sd1, NcmRNG *rng);
gdouble ncm_stats_dist1d_eval_mode (NcmStatsDist1d *sd1);

G_END_DECLS

#endif /* _NCM_STATS_DIST1D_H_ */
