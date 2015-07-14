/***************************************************************************
 *            nc_data_xcor.c
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "xcor/nc_data_xcor.h"
#include "xcor/nc_xcor.h"
#include "nc_snia_dist_cov.h"
#include "math/ncm_model_ctrl.h"
#include "math/ncm_lapack.h"
#include "math/ncm_cfg.h"

#include <glib/gstdio.h>
#ifdef NUMCOSMO_HAVE_CFITSIO
#include <fitsio.h>
#endif /* NUMCOSMO_HAVE_CFITSIO */

#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>

enum
{
	PROP_0,
	PROP_NOBS,
	PROP_ELL,
	PROP_CLDATA,
	PROP_OBS_POS,
	PROP_X_MATRIX_1,
	PROP_X_MATRIX_2,
	PROP_MIXING,
	PROP_XC,
	PROP_SIZE,
};

G_DEFINE_TYPE (NcDataXcor, nc_data_xcor, NCM_TYPE_DATA_GAUSS_COV);

static void
nc_data_xcor_init (NcDataXcor* xcdata)
{
	xcdata->mu_len = 0;
	xcdata->nell = 0;
	xcdata->ell = NULL;
	xcdata->nobs = 0; /* number of observables*/
	xcdata->ncl = 0; /* number of auto and cross spectra = nobs(nobs+1)/2 */

	xcdata->X_matrix_1 = NULL;
	xcdata->X_matrix_2 = NULL; /* X matrices (=mask dependent, cosmology independent part of the covariances <C_l^{a,b}C_l'^{c,d}>) */
	xcdata->mixing = NULL; /* mixing^{a,b}_{l,l'} (size=(ncl*nell, nell)) */
	xcdata->xc = NULL;

	xcdata->clorder = NULL;

	xcdata->cosmo_ctrl = ncm_model_ctrl_new (NULL);
	xcdata->xcl_ctrl_array = NULL;
}

static void
nc_data_xcor_set_property (GObject* object, guint prop_id, const GValue* value, GParamSpec* pspec)
{
	NcDataXcor* xcdata = NC_DATA_XCOR (object);

	switch (prop_id)
	{
	case PROP_NOBS:
	{
		xcdata->nobs = g_value_get_int (value);
		break;
	}
	case PROP_ELL:
	{
		fprintf (stderr, "courgette\n");
		xcdata->ell = g_value_get_object (value);
		break;
	}
	case PROP_CLDATA:
	{
		fprintf (stderr, "tomate\n");

		xcdata->cldata = g_value_get_object (value);
		break;
	}
	case PROP_X_MATRIX_1:
	{
		fprintf (stderr, "aubergine\n");

		xcdata->X_matrix_1 = g_value_get_object (value);
		break;
	}
	case PROP_X_MATRIX_2:
	{
		fprintf (stderr, "sel\n");

		xcdata->X_matrix_2 = g_value_get_object (value);
		break;
	}
	case PROP_MIXING:
	{
		fprintf (stderr, "poivre\n");

		xcdata->mixing = g_value_get_object (value);
		break;
	}
	case PROP_XC:
	{
		fprintf (stderr, "lardon\n");

		xcdata->xc = g_value_get_object (value);
		break;
	}
	default:
		G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
		break;
	}
}

static void
nc_data_xcor_get_property (GObject* object, guint prop_id, GValue* value, GParamSpec* pspec)
{
	NcDataXcor* xcdata = NC_DATA_XCOR (object);

	switch (prop_id)
	{
	case PROP_NOBS:
	{
		g_value_set_int (value, xcdata->nobs);
		break;
	}
	case PROP_ELL:
	{
		g_value_set_object (value, xcdata->ell);
		break;
	}
	case PROP_CLDATA:
	{
		g_value_set_object (value, xcdata->cldata);
		break;
	}
	case PROP_X_MATRIX_1:
	{
		g_value_set_object (value, xcdata->X_matrix_1);
		break;
	}
	case PROP_X_MATRIX_2:
	{
		g_value_set_object (value, xcdata->X_matrix_2);
		break;
	}
	case PROP_MIXING:
	{
		g_value_set_object (value, xcdata->mixing);
		break;
	}
	case PROP_XC:
	{
		g_value_set_object (value, xcdata->xc);
		break;
	}
	default:
		G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
		break;
	}
}

static void
nc_data_xcor_constructed (GObject* object)
{
	/* Chain up : start */
	G_OBJECT_CLASS (nc_data_xcor_parent_class)
	->constructed (object);
	{
	}
}

static void
nc_data_xcor_dispose (GObject* object)
{
	NcDataXcor* xcdata = NC_DATA_XCOR (object);
	NcmDataGaussCov* gauss = NCM_DATA_GAUSS_COV (object);
	ncm_data_gauss_cov_set_size (gauss, 0);

	ncm_model_ctrl_clear (&xcdata->cosmo_ctrl);
	g_ptr_array_unref (xcdata->xcl_ctrl_array);

	ncm_matrix_clear (&xcdata->X_matrix_1);
	ncm_matrix_clear (&xcdata->X_matrix_2);
	ncm_matrix_clear (&xcdata->mixing);
	ncm_matrix_clear (&xcdata->clorder);

	ncm_vector_clear (&xcdata->ell);

	nc_xcor_free (xcdata->xc);

	/* Chain up : end */
	G_OBJECT_CLASS (nc_data_xcor_parent_class)
	->dispose (object);
}

static void
nc_data_xcor_finalize (GObject* object)
{
	/* Chain up : end */
	G_OBJECT_CLASS (nc_data_xcor_parent_class)
	->finalize (object);
}

static void _nc_data_xcor_prepare (NcmData* data, NcmMSet* mset);
static void _nc_data_xcor_compute_cl (NcDataXcor* xcdata, NcmMSet* mset);
static void _nc_data_xcor_mean_func (NcmDataGaussCov* gauss, NcmMSet* mset, NcmVector* vp);
static gboolean _nc_data_xcor_cov_func (NcmDataGaussCov* gauss, NcmMSet* mset, NcmMatrix* cov);

static void
nc_data_xcor_class_init (NcDataXcorClass* klass)
{
	GObjectClass* object_class = G_OBJECT_CLASS (klass);
	NcmDataClass* data_class = NCM_DATA_CLASS (klass);
	NcmDataGaussCovClass* gauss_class = NCM_DATA_GAUSS_COV_CLASS (klass);

	object_class->set_property = &nc_data_xcor_set_property;
	object_class->get_property = &nc_data_xcor_get_property;
	object_class->constructed = &nc_data_xcor_constructed;
	object_class->dispose = &nc_data_xcor_dispose;
	object_class->finalize = &nc_data_xcor_finalize;

	g_object_class_install_property (object_class,
	                                 PROP_NOBS,
	                                 g_param_spec_int ("nobs",
	                                                   NULL,
	                                                   "Number of observables",
	                                                   0, 100, 0,
	                                                   G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_ELL,
	                                 g_param_spec_object ("ell",
	                                                      NULL,
	                                                      "Multipole values",
	                                                      NCM_TYPE_VECTOR,
	                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	g_object_class_install_property (object_class,
	                                 PROP_CLDATA,
	                                 g_param_spec_object ("cldata",
	                                                      NULL,
	                                                      "Cross and auto spectra from data",
	                                                      NCM_TYPE_VECTOR,
	                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


	g_object_class_install_property (object_class,
	                                 PROP_OBS_POS,
	                                 g_param_spec_pointer ("obs-pos",
	                                                       NULL,
	                                                       "Indices of the xcl objects in mset",
	                                                       // G_TYPE_POINTER,
	                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	g_object_class_install_property (object_class,
	                                 PROP_X_MATRIX_1,
	                                 g_param_spec_object ("X-matrix-1",
	                                                      NULL,
	                                                      "X matrix 1 for the covariance of partial sky pseudo spectra estimator",
	                                                      NCM_TYPE_MATRIX,
	                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	g_object_class_install_property (object_class,
	                                 PROP_X_MATRIX_2,
	                                 g_param_spec_object ("X-matrix-2",
	                                                      NULL,
	                                                      "X matrix 2 for the covariance of partial sky pseudo spectra estimator",
	                                                      NCM_TYPE_MATRIX,
	                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	g_object_class_install_property (object_class,
	                                 PROP_MIXING,
	                                 g_param_spec_object ("mixing",
	                                                      NULL,
	                                                      "Mixing matrix for partial sky pseudo spectra estimator",
	                                                      NCM_TYPE_MATRIX,
	                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	g_object_class_install_property (object_class,
	                                 PROP_XC,
	                                 g_param_spec_object ("xc",
	                                                      NULL,
	                                                      "Xcor object to compute theoretical spectra",
	                                                      NC_TYPE_XCOR,
	                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	data_class->prepare = &_nc_data_xcor_prepare;

	gauss_class->mean_func = &_nc_data_xcor_mean_func;
	gauss_class->cov_func = &_nc_data_xcor_cov_func;
}

static void
_nc_data_xcor_prepare (NcmData* data, NcmMSet* mset)
{
	fprintf (stderr, "pomme\n");

	NcDataXcor* xcdata = NC_DATA_XCOR (data);
	NcHICosmo* cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));

	fprintf (stderr, "poire\n");

	gboolean doprep = ncm_model_ctrl_update (xcdata->cosmo_ctrl, NCM_MODEL (cosmo)); /*if any xcl or cosmo has been updated */

	fprintf (stderr, "pruneau\n");

	if (doprep)
	{
		nc_xcor_prepare (xcdata->xc, cosmo);
	}
	fprintf (stderr, "peche\n");


	guint i;
	gboolean doprepxcl;
	for (i = 0; i < xcdata->nobs; i++)
	{
		NcmModelCtrl* ctrl = g_ptr_array_index (xcdata->xcl_ctrl_array, i);
		NcXcorLimber* xcl = NC_XCOR_LIMBER (ncm_mset_peek_pos (mset, nc_xcor_limber_id (), i));

		doprepxcl = ncm_model_ctrl_update (ctrl, NCM_MODEL (xcl));

		if (doprepxcl)
		{
			nc_xcor_limber_prepare (xcl, cosmo);
		}

		doprep = doprepxcl | doprep;
	}
	fprintf (stderr, "prune");


	if (doprep) _nc_data_xcor_compute_cl (xcdata, mset);
}

static void
_nc_data_xcor_compute_cl (NcDataXcor* xcdata, NcmMSet* mset)
{
	NcHICosmo* cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
	NcXcor* xc = xcdata->xc;
	guint nell = xcdata->nell;

	guint obs, diag, k;
	guint pos_ell = 0;
	NcmVector* clpure = ncm_vector_new (nell);
	NcmVector* clmixed = ncm_vector_new (nell);
	gdouble cltemp;

	guint i, vp_pos;
	gint ret;
	/* First get auto-spectra with noise */
	for (obs = 0; obs < xcdata->nobs; obs++)
	{
		NcXcorLimber* xcl = NC_XCOR_LIMBER (ncm_mset_peek_pos (mset, nc_xcor_limber_id (), obs));

		nc_xcor_limber_auto_cl (xc, xcl, cosmo, xcdata->ell, clpure, 0, TRUE); //NcDistance* dist, NcTransferFunc* tf, NcGrowthFunc* gf,
		NcmMatrix* sub_mixing = ncm_matrix_get_submatrix (xcdata->mixing, 0, pos_ell, nell, nell);
		ret = gsl_blas_dgemv (CblasNoTrans, 1.0, ncm_matrix_gsl (sub_mixing), ncm_vector_gsl (clpure), 0.0, ncm_vector_gsl (clmixed));
		if (ret != 0)
			g_error ("_nc_data_xcor_compute_cl : %i", ret);
		for (i = 0; i < nell; i++)
		{
			cltemp = ncm_vector_get (clmixed, i);
			vp_pos = pos_ell + i;
			ncm_vector_set (xcdata->clth, vp_pos, cltemp);
		}
		pos_ell += nell;
	}

	/* Then get cross-spectra */
	if (xcdata->nobs > 1)
	{
		for (diag = 1; diag < xcdata->nobs; diag++)
		{
			for (k = 0; k < xcdata->nobs - diag; k++)
			{
				NcXcorLimber* xcl1 = NC_XCOR_LIMBER (ncm_mset_peek_pos (mset, nc_xcor_limber_id (), k));
				NcXcorLimber* xcl2 = NC_XCOR_LIMBER (ncm_mset_peek_pos (mset, nc_xcor_limber_id (), diag + k));


				nc_xcor_limber_cross_cl (xc, xcl1, xcl2, cosmo, xcdata->ell, clpure, 0);
				NcmMatrix* sub_mixing = ncm_matrix_get_submatrix (xcdata->mixing, 0, pos_ell, nell, nell);
				ret = gsl_blas_dgemv (CblasNoTrans, 1.0, ncm_matrix_gsl (sub_mixing), ncm_vector_gsl (clpure), 0.0, ncm_vector_gsl (clmixed));
				if (ret != 0)
					g_error ("_nc_data_xcor_compute_cl : %i", ret);

				for (i = 0; i < xcdata->nell; i++)
				{
					cltemp = ncm_vector_get (clmixed, i);
					vp_pos = pos_ell + i;
					ncm_vector_set (xcdata->clth, vp_pos, cltemp);
				}
				pos_ell += nell;
			}
		}
	}
}

static void
_nc_data_xcor_mean_func (NcmDataGaussCov* gauss, NcmMSet* mset, NcmVector* vp)
{
	NcDataXcor* xcdata = NC_DATA_XCOR (gauss);
	guint i;
	gdouble cltemp;

	for (i = 0; i < xcdata->mu_len; i++) /* is it mu_len ?? */
	{
		cltemp = ncm_vector_get (xcdata->clth, i);
		ncm_vector_set (vp, i, cltemp);
	}
}

static gboolean
_nc_data_xcor_cov_func (NcmDataGaussCov* gauss, NcmMSet* mset, NcmMatrix* cov)
{
	NcDataXcor* xcdata = NC_DATA_XCOR (gauss);

	guint diag1, k1, diag2, k2;
	guint a, b, c, d, l, ll;
	guint ladidx, lladidx, lbcidx, llbcidx, lacidx, llacidx, lbdidx, llbdidx, Xlidx, Xllidx;
	gdouble res;

	/* For loop over AB */
	for (diag1 = 0; diag1 < xcdata->nobs; diag1++)
	{
		for (k1 = 0; k1 < xcdata->nobs - diag1; k1++)
		{
			a = k1;
			b = diag1 + k1;

			/* For loop over CD */
			for (diag2 = 0; diag2 < xcdata->nobs; diag2++)
			{
				for (k2 = 0; k2 < xcdata->nobs - diag2; k2++)
				{
					c = k2;
					d = diag2 + k2;

					/* For loop over multipoles */
					for (l = 0; l < xcdata->nell; l++)
					{
						for (ll = 0; ll < xcdata->nell; ll++)
						{
							ladidx = ncm_matrix_get (xcdata->clorder, a, d) * xcdata->nell + l;
							lladidx = ncm_matrix_get (xcdata->clorder, a, d) * xcdata->nell + ll;
							lbcidx = ncm_matrix_get (xcdata->clorder, b, c) * xcdata->nell + l;
							llbcidx = ncm_matrix_get (xcdata->clorder, b, c) * xcdata->nell + ll;

							lacidx = ncm_matrix_get (xcdata->clorder, a, c) * xcdata->nell + l;
							llacidx = ncm_matrix_get (xcdata->clorder, a, c) * xcdata->nell + ll;
							lbdidx = ncm_matrix_get (xcdata->clorder, b, d) * xcdata->nell + l;
							llbdidx = ncm_matrix_get (xcdata->clorder, b, d) * xcdata->nell + ll;

							Xlidx = ncm_matrix_get (xcdata->clorder, a, b) * xcdata->nell + l;
							Xllidx = ncm_matrix_get (xcdata->clorder, c, d) * xcdata->nell + ll;

							res = 0.0;

							res += sqrt (fabs (ncm_vector_get (xcdata->clth, ladidx) *
							                   ncm_vector_get (xcdata->clth, lladidx) *
							                   ncm_vector_get (xcdata->clth, lbcidx) *
							                   ncm_vector_get (xcdata->clth, llbcidx))) *
							       ncm_matrix_get (xcdata->X_matrix_1, Xlidx, Xllidx);

							res += sqrt (fabs (ncm_vector_get (xcdata->clth, lacidx) *
							                   ncm_vector_get (xcdata->clth, llacidx) *
							                   ncm_vector_get (xcdata->clth, lbdidx) *
							                   ncm_vector_get (xcdata->clth, llbdidx))) *
							       ncm_matrix_get (xcdata->X_matrix_2, Xlidx, Xllidx);

							ncm_matrix_set (cov, Xlidx, Xllidx, res);
						}
					}
				}
			}
		}
	}

	return TRUE;
}

/**
 * nc_data_xcor_new_full:
 * @ell: a #NcmVector
 * @nobs: a #guint
 * @X_matrix_1: a #NcmMatrix
 * @X_matrix_2: a #NcmMatrix
 * @mixing: a #NcmMatrix
 * @xc: a #NcXcor
 * @Clobs: a #NcmVector
 * @use_norma: a #gboolean
 *
 * FIXME
 *
 * Returns: FIXME
 *
*/
NcDataXcor* nc_data_xcor_new_full (NcmVector* ell, const guint nobs, NcmMatrix* X_matrix_1, NcmMatrix* X_matrix_2, NcmMatrix* mixing, NcXcor* xc, NcmVector* Clobs, gboolean use_norma)
{
	guint ncl = nobs * (nobs + 1) / 2;
	guint nell = ncm_vector_len (ell);
	guint mu_len = nell * ncl;

	if (ncm_vector_len (Clobs) != mu_len)
	{
		g_error ("\nThe size of Clobs doesn't match nobs and nell.\n");
	}

	if ((ncm_matrix_nrows (X_matrix_1) != mu_len) | (ncm_matrix_ncols (X_matrix_1) != mu_len) | (ncm_matrix_nrows (X_matrix_2) != mu_len) | (ncm_matrix_ncols (X_matrix_2) != mu_len))
	{
		g_error ("\nThe size of X_matrix_1 or 2 doesn't match nobs and nell.\n");
	}

	if ((ncm_matrix_nrows (mixing) != nell) | (ncm_matrix_ncols (mixing) != mu_len))
	{
		g_error ("\nThe size of mixing doesn't match nobs and nell.\n");
	}

	fprintf (stderr, "pingouin\n");

	NcDataXcor* xcdata = g_object_new (NC_TYPE_DATA_XCOR,
	                                   "use-norma", use_norma,
	                                   "nobs", nobs,
	                                   "ell", ell,
	                                   "cldata", Clobs,
	                                   "X_matrix_1", X_matrix_1,
	                                   "X_matrix_2", X_matrix_2,
	                                   "mixing", mixing,
	                                   "xc", xc, /* ADD EVRYTHING HERE ! */
	                                   NULL);

	fprintf (stderr, "requin\n");

	xcdata->mu_len = mu_len;
	xcdata->nell = nell;
	xcdata->clth = ncm_vector_new (mu_len);
	xcdata->ncl = ncl;
	xcdata->cosmo_ctrl = ncm_model_ctrl_new (NULL);

	/* xcl_ctrl_array initialization */
	GPtrArray* ctrl_array = g_ptr_array_new ();
	guint i;

	fprintf (stderr, "lion\n");


	g_ptr_array_set_free_func (ctrl_array, &g_object_unref);

	fprintf (stderr, "poule\n");


	for (i = 0; i < nobs; i++)
	{
		NcmModelCtrl* ctrl = ncm_model_ctrl_new (NULL);
		g_ptr_array_add (ctrl_array, ctrl);
	}

	fprintf (stderr, "kangourou\n");


	xcdata->xcl_ctrl_array = ctrl_array;

	/* clorder initialization */
	xcdata->clorder = ncm_matrix_new (nobs, nobs);
	i = 0;
	guint diag, k;
	for (diag = 0; diag < xcdata->nobs; diag++)
	{
		for (k = 0; k < xcdata->nobs - diag; k++)
		{
			ncm_matrix_set (xcdata->clorder, k, diag + k, i);
			ncm_matrix_set (xcdata->clorder, diag + k, k, i);
			i++;
		}
	}

	fprintf (stderr, "narwhal\n");

	NcmDataGaussCov* gauss = NCM_DATA_GAUSS_COV (xcdata);
	NcmData* data = NCM_DATA (gauss);
	gauss->np = mu_len;
	gauss->y = ncm_vector_ref (xcdata->cldata);
	gauss->v = ncm_vector_new (gauss->np);
	gauss->cov = ncm_matrix_new (gauss->np, gauss->np);
	data->init = TRUE;

	return xcdata;
}
