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

G_BEGIN_DECLS

#define NCM_TYPE_DATA_POISSON (ncm_data_poisson_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmDataPoisson, ncm_data_poisson, NCM, DATA_POISSON, NcmData)

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

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[15];
};

void ncm_data_poisson_init_from_vector (NcmDataPoisson *poisson, NcmVector *nodes, NcmVector *N);
void ncm_data_poisson_init_zero (NcmDataPoisson *poisson, NcmVector *nodes);
void ncm_data_poisson_init_from_binning (NcmDataPoisson *poisson, NcmVector *nodes, NcmVector *x);

void ncm_data_poisson_set_size (NcmDataPoisson *poisson, guint nbins);
guint ncm_data_poisson_get_size (NcmDataPoisson *poisson);
gdouble ncm_data_poisson_get_sum (NcmDataPoisson *poisson);

NcmVector *ncm_data_poisson_get_hist_vals (NcmDataPoisson *poisson);
NcmVector *ncm_data_poisson_get_hist_means (NcmDataPoisson *poisson, NcmMSet *mset);

G_END_DECLS

#endif /* _NCM_DATA_POISSON_H_ */

