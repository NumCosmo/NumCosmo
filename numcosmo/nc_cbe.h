/***************************************************************************
 *            nc_cbe.h
 *
 *  Sat October 24 11:57:37 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_cbe.h
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

#ifndef _NC_CBE_H_
#define _NC_CBE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model_ctrl.h>
#include <numcosmo/math/ncm_spline2d.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/nc_hireion.h>
#include <numcosmo/nc_hiprim.h>
#include <numcosmo/nc_scalefactor.h>
#include <numcosmo/nc_cbe_precision.h>
#include <numcosmo/data/nc_data_cmb.h>

G_BEGIN_DECLS

#define NC_TYPE_CBE             (nc_cbe_get_type ())
#define NC_CBE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_CBE, NcCBE))
#define NC_CBE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_CBE, NcCBEClass))
#define NC_IS_CBE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_CBE))
#define NC_IS_CBE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_CBE))
#define NC_CBE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_CBE, NcCBEClass))

typedef struct _NcCBEClass NcCBEClass;
typedef struct _NcCBE NcCBE;

typedef struct _NcCBEPrivate NcCBEPrivate;

struct _NcCBEClass
{
  /*< private >*/
  GObjectClass parent_class;
};

typedef void (*NcCBECall) (NcCBE *cbe, NcHICosmo *cosmo);
typedef void (*NcCBEFree) (NcCBE *cbe);

struct _NcCBE
{
  /*< private >*/
  GObject parent_instance;
  NcCBEPrecision *prec;
  NcCBEPrivate *priv;
  NcScalefactor *a;
  guint bg_verbose;
  guint thermo_verbose;
  guint pert_verbose;
  guint transfer_verbose;
  guint prim_verbose;
  guint spectra_verbose;
  guint nonlin_verbose;
  guint lensing_verbose;
  NcDataCMBDataType target_Cls;
  gboolean calc_transfer;
  gboolean use_halofit;
  gboolean use_lensed_Cls;
  gboolean use_tensor;
  gboolean use_thermodyn;
  guint scalar_lmax;
  guint vector_lmax;
  guint tensor_lmax;
  NcmModelCtrl *ctrl_cosmo;
  NcmModelCtrl *ctrl_prim;
  NcCBECall call;
  NcCBEFree free;
  gboolean allocated;
  gboolean thermodyn_prepared;
};

GType nc_cbe_get_type (void) G_GNUC_CONST;

NcCBE *nc_cbe_new (void);
NcCBE *nc_cbe_prec_new (NcCBEPrecision *cbe_prec);
NcCBE *nc_cbe_prec_file_new (gchar *prec_filename);
NcCBE *nc_cbe_ref (NcCBE *cbe);
void nc_cbe_free (NcCBE *cbe);
void nc_cbe_clear (NcCBE **cbe);

void nc_cbe_set_precision (NcCBE *cbe, NcCBEPrecision *cbe_prec);
void nc_cbe_set_target_Cls (NcCBE *cbe, NcDataCMBDataType target_Cls);
void nc_cbe_set_calc_transfer (NcCBE *cbe, gboolean calc_transfer);
void nc_cbe_set_halofit (NcCBE *cbe, gboolean use_halofit);
void nc_cbe_set_lensed_Cls (NcCBE *cbe, gboolean use_lensed_Cls);
void nc_cbe_set_tensor (NcCBE *cbe, gboolean use_tensor);
void nc_cbe_set_thermodyn (NcCBE *cbe, gboolean use_thermodyn);
void nc_cbe_set_scalar_lmax (NcCBE *cbe, guint scalar_lmax);
void nc_cbe_set_vector_lmax (NcCBE *cbe, guint vector_lmax);
void nc_cbe_set_tensor_lmax (NcCBE *cbe, guint tensor_lmax);

void nc_cbe_set_max_matter_pk_z (NcCBE *cbe, gdouble zmax);
void nc_cbe_set_max_matter_pk_k (NcCBE *cbe, gdouble kmax);
gdouble nc_cbe_get_max_matter_pk_z (NcCBE *cbe);
gdouble nc_cbe_get_max_matter_pk_k (NcCBE *cbe);

NcCBEPrecision *nc_cbe_peek_precision (NcCBE *cbe);
NcDataCMBDataType nc_cbe_get_target_Cls (NcCBE *cbe);
gboolean nc_cbe_calc_transfer (NcCBE *cbe);
gboolean nc_cbe_halofit (NcCBE *cbe);
gboolean nc_cbe_lensed_Cls (NcCBE *cbe);
gboolean nc_cbe_tensor (NcCBE *cbe);
gboolean nc_cbe_thermodyn (NcCBE *cbe);
guint nc_cbe_get_scalar_lmax (NcCBE *cbe);
guint nc_cbe_get_vector_lmax (NcCBE *cbe);
guint nc_cbe_get_tensor_lmax (NcCBE *cbe);

void nc_cbe_use_ppf (NcCBE *cbe, gboolean use_ppf);

void nc_cbe_thermodyn_prepare (NcCBE *cbe, NcHICosmo *cosmo);
void nc_cbe_thermodyn_prepare_if_needed (NcCBE *cbe, NcHICosmo *cosmo);
void nc_cbe_prepare (NcCBE *cbe, NcHICosmo *cosmo);
void nc_cbe_prepare_if_needed (NcCBE *cbe, NcHICosmo *cosmo);

gdouble nc_cbe_compare_bg (NcCBE *cbe, NcHICosmo *cosmo, gboolean log_cmp);

NcmSpline *nc_cbe_thermodyn_get_Xe (NcCBE *cbe);
gdouble nc_cbe_thermodyn_v_tau_max_z (NcCBE *cbe);
gdouble nc_cbe_thermodyn_z_d (NcCBE *cbe);

NcmSpline2d *nc_cbe_get_matter_ps (NcCBE *cbe);

gdouble nc_cbe_get_sigma8 (NcCBE *cbe);

void nc_cbe_get_all_Cls (NcCBE *cbe, NcmVector *PHIPHI_Cls, NcmVector *TT_Cls, NcmVector *EE_Cls, NcmVector *BB_Cls, NcmVector *TE_Cls);

G_END_DECLS

#endif /* _NC_CBE_H_ */

