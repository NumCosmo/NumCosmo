/***************************************************************************
 *            nc_data_xcor.h
 *
 *  Tue July 14 12:00:00 2015
 *  Copyright  2015  Cyrille Doux
 *  <cdoux@apc.in2p3.fr>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2015 Cyrille Doux <cdoux@apc.in2p3.fr>
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

#ifndef _NC_DATA_XCOR_H_
#define _NC_DATA_XCOR_H_

#include <glib-object.h>
#include <glib.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data_gauss_cov.h>
#include <numcosmo/math/ncm_model_ctrl.h>
#include <numcosmo/math/ncm_obj_array.h>
#include <numcosmo/xcor/nc_xcor.h>
#include <numcosmo/xcor/nc_xcor_AB.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_XCOR (nc_data_xcor_get_type ())
#define NC_DATA_XCOR(obj) (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DATA_XCOR, NcDataXcor))
#define NC_DATA_XCOR_CLASS(klass) (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DATA_XCOR, NcDataXcorClass))
#define NC_IS_DATA_XCOR(obj) (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DATA_XCOR))
#define NC_IS_DATA_XCOR_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DATA_XCOR))
#define NC_DATA_XCOR_GET_CLASS(obj) (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DATA_XCOR, NcDataXcorClass))

#define NC_DATA_XCOR_MAX (5)

// typedef struct _NcDataXcorAB
// {
// 	gint a;
// 	gint b;
//
// 	gint ell_th_cut_off;
// 	gint ell_lik_min;
// 	gint ell_lik_max;
// 	gint nell_lik;
//
// 	NcmMatrix* mixing;
// 	NcmMatrix* cl_th; //column 0 : C_l^th, 1 : C_l^th+N_l, 2 : mixed C_l
// 	NcmVector* cl_obs;
//
// } NcDataXcorAB;

typedef struct _NcDataXcorClass NcDataXcorClass;
typedef struct _NcDataXcor NcDataXcor;

struct _NcDataXcor
{
	/*< private >*/
	NcmDataGaussCov parent_instance;

	guint nobs;

	NcXcorAB* xcab[NC_DATA_XCOR_MAX][NC_DATA_XCOR_MAX];

	NcmObjArray* xcab_oa;
	// guint xcab_oa_idx[NC_DATA_XCOR_MAX * NC_DATA_XCOR_MAX][2];
	// guint xcab_oa_ctr;

	gint xcidx[NC_DATA_XCOR_MAX][NC_DATA_XCOR_MAX];
	guint xcidx_ctr;

	// gboolean X_init;
	NcmMatrix* X1; //[NC_DATA_XCOR_MAX_NOBS][NC_DATA_XCOR_MAX_NOBS][NC_DATA_XCOR_MAX_NOBS][NC_DATA_XCOR_MAX_NOBS];
	NcmMatrix* X2; //[NC_DATA_XCOR_MAX_NOBS][NC_DATA_XCOR_MAX_NOBS][NC_DATA_XCOR_MAX_NOBS][NC_DATA_XCOR_MAX_NOBS]; /* X matrices (=mask dependent, cosmology independent part of the covariances <C_l^{a,b}C_l'^{c,d}>) */

	NcmVector* pcl;
	NcmMatrix* pcov;

	NcXcor* xc;

	NcmModelCtrl* cosmo_ctrl;
	GPtrArray* xclk_ctrl;
};

struct _NcDataXcorClass
{
	/*< private >*/
	NcmDataGaussCovClass parent_class;
};

GType nc_data_xcor_get_type (void) G_GNUC_CONST;

// NcmData* nc_data_xcor_new (gboolean use_norma);
NcDataXcor* nc_data_xcor_new_full (const guint nobs, NcXcor* xc, gboolean use_norma); //, const gchar* xcname[]);
void nc_data_xcor_set_2 (NcDataXcor* dxc, guint a, guint b, guint ell_th_cut_off, guint ell_lik_min, guint ell_lik_max, const gchar* clobs_filename, const gchar* mixing_filename, const guint mixing_filelength);
void nc_data_xcor_set_AB (NcDataXcor* dxc, NcXcorAB* xcab);
void nc_data_xcor_set_3 (NcDataXcor* dxc);
void nc_data_xcor_set_4 (NcDataXcor* dxc, guint a, guint b, guint c, guint d, const gchar* X1_filename, const gchar* X2_filename, guint X_filelength);
void nc_data_xcor_set_5 (NcDataXcor* dxc);
void nc_data_xcor_mean_func_ab (NcDataXcor* dxc, NcmVector* vp, guint a, guint b);
void nc_data_xcor_get_cl_obs (NcDataXcor* dxc, NcmVector* vp, guint a, guint b);
void nc_data_xcor_cov_func_abcd (NcDataXcor* dxc, NcmMatrix* cov, guint a, guint b, guint c, guint d);

#define NC_DATA_XCOR_DL 10

G_END_DECLS

#endif /* _NC_DATA_XCOR_H_ */
