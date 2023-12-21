/***************************************************************************
 *            ncm_stats_dist2d.h
 *
 *  Sat July 22 16:21:25 2017
 *  Copyright  2017  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * ncm_stats_dist2d.h
 * Copyright (C) 2017 Mariana Penna Lima <pennalima@gmail.com>
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

#ifndef _NCM_STATS_DIST2D_H_
#define _NCM_STATS_DIST2D_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_ode_spline.h>
#include <numcosmo/math/ncm_rng.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_min.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NCM_TYPE_STATS_DIST2D (ncm_stats_dist2d_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmStatsDist2d, ncm_stats_dist2d, NCM, STATS_DIST2D, GObject)

struct _NcmStatsDist2dClass
{
  /*< private >*/
  GObjectClass parent_class;

  void (*xbounds) (NcmStatsDist2d *sd2, gdouble *xi, gdouble *xf);
  void (*ybounds) (NcmStatsDist2d *sd2, gdouble *yi, gdouble *yf);
  gdouble (*pdf) (NcmStatsDist2d *sd2, const gdouble x, const gdouble y);
  gdouble (*m2lnp) (NcmStatsDist2d *sd2, const gdouble x, const gdouble y);
  gdouble (*cdf) (NcmStatsDist2d *sd2, const gdouble x, const gdouble y);
  gdouble (*marginal_pdf) (NcmStatsDist2d *sd2, const gdouble xy);
  gdouble (*marginal_cdf) (NcmStatsDist2d *sd2, const gdouble xy);
  gdouble (*marginal_inv_cdf) (NcmStatsDist2d *sd2, const gdouble u);
  gdouble (*inv_cond) (NcmStatsDist2d *sd2, const gdouble u, const gdouble xy);
  void (*prepare) (NcmStatsDist2d *sd2);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[8];
};

NcmStatsDist2d *ncm_stats_dist2d_ref (NcmStatsDist2d *sd2);
void ncm_stats_dist2d_free (NcmStatsDist2d *sd2);
void ncm_stats_dist2d_clear (NcmStatsDist2d **sd2);

void ncm_stats_dist2d_prepare (NcmStatsDist2d *sd2);

void ncm_stats_dist2d_xbounds (NcmStatsDist2d *sd2, gdouble *xi, gdouble *xf);
void ncm_stats_dist2d_ybounds (NcmStatsDist2d *sd2, gdouble *yi, gdouble *yf);

gdouble ncm_stats_dist2d_eval_pdf (NcmStatsDist2d *sd2, const gdouble x, const gdouble y);
gdouble ncm_stats_dist2d_eval_m2lnp (NcmStatsDist2d *sd2, const gdouble x, const gdouble y);
gdouble ncm_stats_dist2d_eval_cdf (NcmStatsDist2d *sd2, const gdouble x, const gdouble y);

gdouble ncm_stats_dist2d_eval_marginal_pdf (NcmStatsDist2d *sd2, const gdouble xy);
gdouble ncm_stats_dist2d_eval_marginal_cdf (NcmStatsDist2d *sd2, const gdouble xy);
gdouble ncm_stats_dist2d_eval_marginal_inv_cdf (NcmStatsDist2d *sd2, const gdouble xy);
gdouble ncm_stats_dist2d_eval_inv_cond (NcmStatsDist2d *sd2, const gdouble u, const gdouble xy);

G_END_DECLS

#endif /* _NCM_STATS_DIST2D_H_ */

