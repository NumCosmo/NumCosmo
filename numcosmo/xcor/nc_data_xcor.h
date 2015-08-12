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

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data_gauss_cov.h>
#include <numcosmo/math/ncm_model_ctrl.h>
#include <numcosmo/xcor/nc_xcor.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_XCOR (nc_data_xcor_get_type ())
#define NC_DATA_XCOR(obj) (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DATA_XCOR, NcDataXcor))
#define NC_DATA_XCOR_CLASS(klass) (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DATA_XCOR, NcDataXcorClass))
#define NC_IS_DATA_XCOR(obj) (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DATA_XCOR))
#define NC_IS_DATA_XCOR_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DATA_XCOR))
#define NC_DATA_XCOR_GET_CLASS(obj) (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DATA_XCOR, NcDataXcorClass))

typedef struct _NcDataXcorClass NcDataXcorClass;
typedef struct _NcDataXcor NcDataXcor;

struct _NcDataXcor
{
	/*< private >*/
	NcmDataGaussCov parent_instance;
	guint mu_len;

	guint nobs; /* number of observables */
	guint nell; /* number of multipoles */

	NcmVector* ell; /* mutipole vector */
	NcmVector* cldata; /* contains the pseudo C_l^{a,b} from data (size = ncl * nell) */
	NcmVector* clth; /* contains the theoretical pseudo C_l^{a,b} from model after mask-mixing (size = ncl * nell) */

	guint ncl; /* number of auto and cross spectra = nobs*(nobs+1)/2 */

	NcmMatrix* clorder; /* Healpix ordering of cross-spectra : AA BB CC AB BC AC */

	NcmMatrix* X_matrix_1;
	NcmMatrix* X_matrix_2; /* X matrices (=mask dependent, cosmology independent part of the covariances <C_l^{a,b}C_l'^{c,d}>) */
	NcmMatrix* mixing; /* mixing^{a,b}_{l,l'} (size=(ncl*nell, nell)) */

	NcXcor* xc;

	NcmModelCtrl* cosmo_ctrl;
	GPtrArray* xcl_ctrl_array;
};

struct _NcDataXcorClass
{
	/*< private >*/
	NcmDataGaussCovClass parent_class;
};

GType nc_data_xcor_get_type (void) G_GNUC_CONST;

NcmData* nc_data_xcor_new (gboolean use_norma);
NcDataXcor* nc_data_xcor_new_full (NcmVector* ell, const guint nobs, NcmMatrix* X_matrix_1, NcmMatrix* X_matrix_2, NcmMatrix* mixing, NcXcor* xc, NcmVector* Clobs, gboolean use_norma);

G_END_DECLS

#endif /* _NC_DATA_XCOR_H_ */
