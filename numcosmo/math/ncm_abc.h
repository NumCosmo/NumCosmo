/***************************************************************************
 *            ncm_abc.h
 *
 *  Tue September 30 15:46:33 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_abc.h
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_ABC_H_
#define _NCM_ABC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_timer.h>
#include <numcosmo/math/ncm_dataset.h>
#include <numcosmo/math/ncm_mset_trans_kern.h>
#include <numcosmo/math/ncm_mset_catalog.h>
#include <numcosmo/math/ncm_memory_pool.h>

G_BEGIN_DECLS

#define NCM_TYPE_ABC (ncm_abc_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmABC, ncm_abc, NCM, ABC, GObject)

struct _NcmABCClass
{
  /*< private >*/
  GObjectClass parent_class;

  gboolean (*data_summary) (NcmABC *abc);
  gdouble (*mock_distance) (NcmABC *abc, NcmDataset *dset, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng);
  gdouble (*distance_prob) (NcmABC *abc, gdouble distance);
  void (*update_tkern) (NcmABC *abc);

  const gchar *(*get_desc) (NcmABC *abc);
  const gchar *(*log_info) (NcmABC *abc);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[12];
};

void ncm_abc_free (NcmABC *abc);
void ncm_abc_clear (NcmABC **abc);

gboolean ncm_abc_data_summary (NcmABC *abc);
gdouble ncm_abc_mock_distance (NcmABC *abc, NcmDataset *dset, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng);
gdouble ncm_abc_distance_prob (NcmABC *abc, gdouble distance);
void ncm_abc_update_tkern (NcmABC *abc);
const gchar *ncm_abc_get_desc (NcmABC *abc);
const gchar *ncm_abc_log_info (NcmABC *abc);

void ncm_abc_set_mtype (NcmABC *abc, NcmFitRunMsgs mtype);
void ncm_abc_set_data_file (NcmABC *abc, const gchar *filename);
void ncm_abc_set_nthreads (NcmABC *abc, guint nthreads);
void ncm_abc_set_rng (NcmABC *abc, NcmRNG *rng);
void ncm_abc_set_trans_kern (NcmABC *abc, NcmMSetTransKern *tkern);

gdouble ncm_abc_get_dist_quantile (NcmABC *abc, gdouble p);
gdouble ncm_abc_get_accept_ratio (NcmABC *abc);
void ncm_abc_update_epsilon (NcmABC *abc, gdouble epsilon);
gdouble ncm_abc_get_epsilon (NcmABC *abc);
gdouble ncm_abc_get_depsilon (NcmABC *abc);
NcmFitRunMsgs ncm_abc_get_mtype (NcmABC *abc);

NcmMSetCatalog *ncm_abc_peek_catalog (NcmABC *abc);
NcmDataset *ncm_abc_peek_dataset (NcmABC *abc);
NcmMatrix *ncm_abc_peek_covar (NcmABC *abc);
NcmMSetTransKern *ncm_abc_peek_trans_kern (NcmABC *abc);

void ncm_abc_start_run (NcmABC *abc);
void ncm_abc_end_run (NcmABC *abc);
void ncm_abc_reset (NcmABC *abc);
void ncm_abc_run (NcmABC *abc, guint nparticles);
void ncm_abc_mean_covar (NcmABC *abc, NcmFit *fit);

void ncm_abc_start_update (NcmABC *abc);
void ncm_abc_end_update (NcmABC *abc);
void ncm_abc_update (NcmABC *abc);

#define NCM_ABC_MIN_SYNC_INTERVAL (10.0)

G_END_DECLS

#endif /* _NCM_ABC_H_ */

