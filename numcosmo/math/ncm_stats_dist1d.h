/***************************************************************************
 *            ncm_stats_dist1d.h
 *
 *  Thu February 12 15:37:22 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_dist1d.h
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

#ifndef _NCM_STATS_DIST1D_H_
#define _NCM_STATS_DIST1D_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_ode_spline.h>
#include <numcosmo/math/ncm_rng.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_min.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NCM_TYPE_STATS_DIST1D (ncm_stats_dist1d_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmStatsDist1d, ncm_stats_dist1d, NCM, STATS_DIST1D, GObject)

struct _NcmStatsDist1dClass
{
  /*< private >*/
  GObjectClass parent_class;

  gdouble (*p) (NcmStatsDist1d *sd1, gdouble x);
  gdouble (*m2lnp) (NcmStatsDist1d *sd1, gdouble x);
  void (*prepare) (NcmStatsDist1d *sd1);
  gdouble (*get_current_h) (NcmStatsDist1d *sd1);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[14];
};

NcmStatsDist1d *ncm_stats_dist1d_ref (NcmStatsDist1d *sd1);
void ncm_stats_dist1d_free (NcmStatsDist1d *sd1);
void ncm_stats_dist1d_clear (NcmStatsDist1d **sd1);

void ncm_stats_dist1d_prepare (NcmStatsDist1d *sd1);
gdouble ncm_stats_dist1d_get_current_h (NcmStatsDist1d *sd1);

void ncm_stats_dist1d_set_compute_cdf (NcmStatsDist1d *sd1, gboolean compute_cdf);
gboolean ncm_stats_dist1d_get_compute_cdf (NcmStatsDist1d *sd1);

void ncm_stats_dist1d_set_xi (NcmStatsDist1d *sd1, gdouble xi);
void ncm_stats_dist1d_set_xf (NcmStatsDist1d *sd1, gdouble xf);

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

