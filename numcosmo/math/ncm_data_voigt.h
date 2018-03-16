/***************************************************************************
 *            ncm_data_voigt.h
 *
 *  Thu Mar  15 16:02:50 2018
 *  Copyright  2018  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Mariana Penna Lima <pennalima@gmail.com>
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

#ifndef _NCM_DATA_VOIGT_H_
#define _NCM_DATA_VOIGT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mset_func.h>
#include <numcosmo/math/ncm_data.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_histogram.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NCM_TYPE_DATA_VOIGT             (ncm_data_voigt_get_type ())
#define NCM_DATA_VOIGT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_DATA_VOIGT, NcmDataVoigt))
#define NCM_DATA_VOIGT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_DATA_VOIGT, NcmDataVoigtClass))
#define NCM_IS_DATA_VOIGT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_DATA_VOIGT))
#define NCM_IS_DATA_VOIGT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_DATA_VOIGT))
#define NCM_DATA_VOIGT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_DATA_VOIGT, NcmDataVoigtClass))

typedef struct _NcmDataVoigtClass NcmDataVoigtClass;
typedef struct _NcmDataVoigt NcmDataVoigt;

/**
 * NcmDataVoigtType:
 * @NCM_DATA_VOIGT_INT: FIXME
 *
 * FIXME
 */
typedef enum _NcmDataVoigtType
{
  NCM_DATA_VOIGT_INT,
} NcmDataVoigtType;

struct _NcmDataVoigtClass
{
  /* < private > */
  NcmDataClass parent_class;
  gdouble (*mean_func) (NcmDataVoigt *voigt, NcmMSet *mset, guint n);
  void (*set_size) (NcmDataVoigt *voigt, guint nbins);
  guint (*get_size) (NcmDataVoigt *voigt);
};

/**
 * NcmDataVoigt:
 *
 * FIXME
 */
struct _NcmDataVoigt
{
  /* < private > */
  NcmData parent_instance;
  gsl_histogram *h;
	NcmVector *means;
  NcmVector *log_Nfac;
  guint nbins;
};

GType ncm_data_voigt_get_type (void) G_GNUC_CONST;

void ncm_data_voigt_init_from_vector (NcmDataVoigt *voigt, NcmVector *nodes, NcmVector *N);
void ncm_data_voigt_init_from_histogram (NcmDataVoigt *voigt, gsl_histogram *h);
void ncm_data_voigt_init_zero (NcmDataVoigt *voigt, NcmVector *nodes);
void ncm_data_voigt_init_from_binning (NcmDataVoigt *voigt, NcmVector *nodes, NcmVector *x);

void ncm_data_voigt_set_size (NcmDataVoigt *voigt, guint nbins);
guint ncm_data_voigt_get_size (NcmDataVoigt *voigt);
gdouble ncm_data_voigt_get_sum (NcmDataVoigt *voigt);

NcmVector *ncm_data_voigt_get_hist_vals (NcmDataVoigt *voigt);
NcmVector *ncm_data_voigt_get_hist_means (NcmDataVoigt *voigt, NcmMSet *mset);

G_END_DECLS

#endif /* _NCM_DATA_VOIGT_H_ */

