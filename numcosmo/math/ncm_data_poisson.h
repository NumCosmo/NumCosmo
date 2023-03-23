/***************************************************************************
 *            ncm_data_poisson.h
 *
 *  Sun Apr  4 21:57:50 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_DATA_POISSON_H_
#define _NCM_DATA_POISSON_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mset_func.h>
#include <numcosmo/math/ncm_data.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_histogram.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NCM_TYPE_DATA_POISSON             (ncm_data_poisson_get_type ())
#define NCM_DATA_POISSON(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_DATA_POISSON, NcmDataPoisson))
#define NCM_DATA_POISSON_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_DATA_POISSON, NcmDataPoissonClass))
#define NCM_IS_DATA_POISSON(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_DATA_POISSON))
#define NCM_IS_DATA_POISSON_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_DATA_POISSON))
#define NCM_DATA_POISSON_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_DATA_POISSON, NcmDataPoissonClass))

typedef struct _NcmDataPoissonClass NcmDataPoissonClass;
typedef struct _NcmDataPoisson NcmDataPoisson;

/**
 * NcmDataPoissonType:
 * @NCM_DATA_POISSON_INT: FIXME
 *
 * FIXME
 */
typedef enum _NcmDataPoissonType
{
  NCM_DATA_POISSON_INT,
} NcmDataPoissonType;

struct _NcmDataPoissonClass
{
  /* < private > */
  NcmDataClass parent_class;
  gdouble (*mean_func) (NcmDataPoisson *poisson, NcmMSet *mset, guint n);
  void (*set_size) (NcmDataPoisson *poisson, guint nbins);
  guint (*get_size) (NcmDataPoisson *poisson);
};

/**
 * NcmDataPoisson:
 *
 * FIXME
 */
struct _NcmDataPoisson
{
  /* < private > */
  NcmData parent_instance;
  gsl_histogram *h;
	NcmVector *means;
  NcmVector *log_Nfac;
  guint nbins;
};

GType ncm_data_poisson_get_type (void) G_GNUC_CONST;

void ncm_data_poisson_init_from_vector (NcmDataPoisson *poisson, NcmVector *nodes, NcmVector *N);
void ncm_data_poisson_init_from_histogram (NcmDataPoisson *poisson, gsl_histogram *h);
void ncm_data_poisson_init_zero (NcmDataPoisson *poisson, NcmVector *nodes);
void ncm_data_poisson_init_from_binning (NcmDataPoisson *poisson, NcmVector *nodes, NcmVector *x);

void ncm_data_poisson_set_size (NcmDataPoisson *poisson, guint nbins);
guint ncm_data_poisson_get_size (NcmDataPoisson *poisson);
gdouble ncm_data_poisson_get_sum (NcmDataPoisson *poisson);

NcmVector *ncm_data_poisson_get_hist_vals (NcmDataPoisson *poisson);
NcmVector *ncm_data_poisson_get_hist_means (NcmDataPoisson *poisson, NcmMSet *mset);

G_END_DECLS

#endif /* _NCM_DATA_POISSON_H_ */

