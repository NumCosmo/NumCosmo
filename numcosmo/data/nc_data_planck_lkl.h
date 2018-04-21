/***************************************************************************
 *            nc_data_planck_lkl.h
 *
 *  Tue October 20 15:24:54 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_data_planck_lkl.h
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

#ifndef _NC_DATA_PLANCK_LKL_H_
#define _NC_DATA_PLANCK_LKL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/data/nc_data_cmb.h>
#include <numcosmo/perturbations/nc_hipert_boltzmann.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_PLANCK_LKL             (nc_data_planck_lkl_get_type ())
#define NC_DATA_PLANCK_LKL(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DATA_PLANCK_LKL, NcDataPlanckLKL))
#define NC_DATA_PLANCK_LKL_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DATA_PLANCK_LKL, NcDataPlanckLKLClass))
#define NC_IS_DATA_PLANCK_LKL(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DATA_PLANCK_LKL))
#define NC_IS_DATA_PLANCK_LKL_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DATA_PLANCK_LKL))
#define NC_DATA_PLANCK_LKL_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DATA_PLANCK_LKL, NcDataPlanckLKLClass))

typedef struct _NcDataPlanckLKLClass NcDataPlanckLKLClass;
typedef struct _NcDataPlanckLKL NcDataPlanckLKL;

struct _NcDataPlanckLKLClass
{
  /*< private >*/
  NcmDataClass parent_class;
};

struct _NcDataPlanckLKL
{
  /*< private >*/
  NcmData parent_instance;
  NcHIPertBoltzmann *pb;
  gchar *filename;
  gpointer obj;
  gboolean is_lensing;
  guint nparams;
  guint ndata_entry;
  gchar **pnames;
  gchar *chksum;
	gdouble check_m2lnL;
  NcDataCMBDataType cmb_data;
  NcmVector *data_params;
  NcmVector *check_data_params;
  NcmVector *data_TT;
  NcmVector *data_EE;
  NcmVector *data_BB;
  NcmVector *data_TE;
  NcmVector *data_TB;
  NcmVector *data_EB;
  NcmVector *data_PHIPHI;
  NcmVector *params;
  NcmModelCtrl *pfi_ctrl;
  NcmModelCtrl *cosmo_ctrl;
  gdouble cm2lnL;
  gdouble A_planck;
  GArray *param_map;
};

GType nc_data_planck_lkl_get_type (void) G_GNUC_CONST;

NcDataPlanckLKL *nc_data_planck_lkl_new (const gchar *filename);
NcDataPlanckLKL *nc_data_planck_lkl_full_new (const gchar *filename, NcHIPertBoltzmann *pb);
const gchar *nc_data_planck_lkl_get_param_name (NcDataPlanckLKL *plik, guint i);
gchar **nc_data_planck_lkl_get_param_names (NcDataPlanckLKL *plik);

void nc_data_planck_lkl_set_hipert_boltzmann (NcDataPlanckLKL *plik, NcHIPertBoltzmann *pb);

G_END_DECLS

#endif /* _NC_DATA_PLANCK_LKL_H_ */

