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

#include "data/nc_data_xcor.h"
#include "math/ncm_cfg.h"
#include "math/ncm_lapack.h"
#include "math/ncm_model_ctrl.h"
#include "nc_hireion.h"
#include "nc_snia_dist_cov.h"
#include "xcor/nc_xcor.h"
#include "xcor/nc_xcor_limber_kernel_gal.h"

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
	PROP_XCAB_OA,
	// PROP_ELL,
	// PROP_ELLTH,
	// PROP_CLDATA,
	// PROP_OBS_POS,
	PROP_X1,
	PROP_X2,
	// PROP_MIXING,
	PROP_XC,
	PROP_SIZE,
};

G_DEFINE_TYPE (NcDataXcor, nc_data_xcor, NCM_TYPE_DATA_GAUSS_COV);

static void
_nc_data_xcor_sort (const gint a, const gint b, gint* aa, gint* bb)
{
	*aa = (a <= b) ? a : b;
	*bb = (a <= b) ? b : a;
}

static void _nc_data_xcor_set_by_oa (NcDataXcor* dxc, NcmObjArray* oa);

static void
nc_data_xcor_init (NcDataXcor* dxc)
{
	dxc->nobs = 0;

	dxc->xcab_oa = NULL; //ncm_obj_array_new ();

	// dxc->xcab_oa_ctr = 0;

	dxc->xcidx_ctr = 0;

	dxc->X1 = NULL;
	dxc->X2 = NULL;

	dxc->pcl = NULL;
	dxc->pcov = NULL;

	dxc->xc = NULL;

	dxc->cosmo_ctrl = ncm_model_ctrl_new (NULL);
	dxc->xclk_ctrl = g_ptr_array_new ();

	guint a, b; //, c, d;
	for (a = 0; a < NC_DATA_XCOR_MAX; a++)
	{
		for (b = 0; b < NC_DATA_XCOR_MAX; b++)
		{
			dxc->xcidx[a][b] = -1;

			// NcXcorAB xcab = { -1, -1, -1, -1, -1, -1, NULL, NULL, NULL };
			dxc->xcab[a][b] = NULL;
		}
	}

	// guint i, j;
	// for (i = 0; i < NC_DATA_XCOR_MAX * NC_DATA_XCOR_MAX; i++)
	// {
	// 	for (j = 0; j < 2; j++)
	// 	{
	// 		dxc->xcab_oa_idx[i][j] = 99;
	// 	}
	// }
	// dxc->xcidx_ctr = 0;
	//
	// dxc->X1 = NULL;
	// dxc->X2 = NULL;
}

static void
nc_data_xcor_set_property (GObject* object, guint prop_id, const GValue* value, GParamSpec* pspec)
{
	NcDataXcor* dxc = NC_DATA_XCOR (object);

	switch (prop_id)
	{
	case PROP_NOBS:
	{
		dxc->nobs = g_value_get_uint (value);
		break;
	}
	case PROP_XCAB_OA:
	{
		// dxc->xcab_oa = g_value_dup_object (value);
		NcmObjArray* oa = (NcmObjArray*)g_value_get_boxed (value);
		if (oa != NULL)
		{
			_nc_data_xcor_set_by_oa (dxc, oa);
		}
		break;
	}
	case PROP_X1:
	{
		dxc->X1 = g_value_dup_object (value);
		break;
	}
	case PROP_X2:
	{
		dxc->X2 = g_value_dup_object (value);
		break;
	}
	case PROP_XC:
	{
		dxc->xc = g_value_dup_object (value);
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
	NcDataXcor* dxc = NC_DATA_XCOR (object);

	switch (prop_id)
	{
	case PROP_NOBS:
	{
		g_value_set_uint (value, dxc->nobs);
		break;
	}
	case PROP_XCAB_OA:
	{
		// g_value_set_object (value, dxc->xcab_oa);
		NcmObjArray* oa = ncm_obj_array_new ();
		// guint a, b;
		guint k, diag;
		// for (a = 0; a < dxc->nobs; a++)
		// {
		// 	for (b = a; b < dxc->nobs; b++)
		for (diag = 0; diag < dxc->nobs; diag++)
		{
			for (k = 0; k < dxc->nobs - diag; k++)
			{
				NcXcorAB* xcab = dxc->xcab[k][diag + k];
				ncm_obj_array_add (oa, G_OBJECT (xcab));
			}
		}
		g_value_take_boxed (value, oa);
		break;
	}
	case PROP_X1:
	{
		g_value_set_object (value, dxc->X1);
		break;
	}
	case PROP_X2:
	{
		g_value_set_object (value, dxc->X2);
		break;
	}
	case PROP_XC:
	{
		g_value_set_object (value, dxc->xc);
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
	G_OBJECT_CLASS (nc_data_xcor_parent_class)->constructed (object);
	{
		NcDataXcor* dxc = NC_DATA_XCOR (object);

		dxc->cosmo_ctrl = ncm_model_ctrl_new (NULL);

		/* xcl_ctrl_array initialization */
		g_ptr_array_set_free_func (dxc->xclk_ctrl, &g_object_unref);

		guint i;

		for (i = 0; i < dxc->nobs; i++)
		{
			NcmModelCtrl* ctrl = ncm_model_ctrl_new (NULL);
			g_ptr_array_add (dxc->xclk_ctrl, ctrl);
		}
	}
}

static void
nc_data_xcor_dispose (GObject* object)
{
	NcDataXcor* dxc = NC_DATA_XCOR (object);
	NcmDataGaussCov* gauss = NCM_DATA_GAUSS_COV (object);
	ncm_data_gauss_cov_set_size (gauss, 0);

	ncm_model_ctrl_clear (&dxc->cosmo_ctrl);
	g_ptr_array_unref (dxc->xclk_ctrl);

	ncm_matrix_clear (&dxc->X1);
	ncm_matrix_clear (&dxc->X2);

	ncm_matrix_clear (&dxc->pcov);
	ncm_vector_clear (&dxc->pcl);

	// ncm_matrix_clear (&dxc->mixing);
	// ncm_matrix_clear (&dxc->clorder);
	//
	// ncm_vector_clear (&dxc->ell);
	// ncm_vector_clear (&dxc->ellth);

	nc_xcor_clear (&dxc->xc);

	ncm_obj_array_clear (&dxc->xcab_oa);
	//
	guint a, b;
	for (a = 0; a < dxc->nobs; a++)
	{
		for (b = a; b < dxc->nobs; b++)
		{
			nc_xcor_AB_clear (&dxc->xcab[a][b]);
		}
	}


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
// static void _nc_data_xcor_compute_cl (NcDataXcor* dxc, NcmMSet* mset);
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
	                                 g_param_spec_uint ("nobs",
	                                                    NULL,
	                                                    "Number of observables",
	                                                    0, NC_DATA_XCOR_MAX, 0,
	                                                    G_PARAM_CONSTRUCT | G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_XCAB_OA,
	                                 g_param_spec_boxed ("xcab-oa",
	                                                     NULL,
	                                                     "NcXcorAB array",
	                                                     NCM_TYPE_OBJ_ARRAY,
	                                                     G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_X1,
	                                 g_param_spec_object ("X1",
	                                                      NULL,
	                                                      "X matrix",
	                                                      NCM_TYPE_MATRIX,
	                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_X2,
	                                 g_param_spec_object ("X2",
	                                                      NULL,
	                                                      "X matrix",
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
	gauss_class->cov_func  = &_nc_data_xcor_cov_func;
}

static gboolean
_nc_data_xcor_fast_update (NcDataXcor* dxc, NcXcorLimberKernel* xcl, guint a, NcmMSet* mset)
{
	if (NC_IS_XCOR_LIMBER_KERNEL_GAL (xcl))
	{
		NcXcorLimberKernelGal* xclkg = NC_XCOR_LIMBER_KERNEL_GAL (xcl);
		if (xclkg->fast_update)
		{
			const gdouble biasratio = *(xclkg->bias) / xclkg->bias_old;
			// printf ("fast bias thingy\n");
			// printf ("%g %g %g\n", *(xclkg->bias), xclkg->bias_old, biasratio);

			NcmVector* cl_th_0_aa = ncm_matrix_get_col (dxc->xcab[a][a]->cl_th, 0);
			NcmVector* cl_th_1_aa = ncm_matrix_get_col (dxc->xcab[a][a]->cl_th, 1);

			ncm_vector_add_constant (cl_th_0_aa, -1.0 * xclkg->noise_bias_old);
			ncm_vector_scale (cl_th_0_aa, gsl_pow_2 (biasratio));
			// ncm_vector_add_constant (ncm_matrix_get_col (dxc->xcab[a][a]->cl_th, 0), ncm_vector_get (NCM_MODEL (xclkg)->params, NC_XCOR_LIMBER_KERNEL_GAL_NOISE_BIAS));
			nc_xcor_limber_kernel_add_noise (xcl, cl_th_0_aa, cl_th_1_aa, 0);

			ncm_vector_free(cl_th_0_aa);
			ncm_vector_free(cl_th_1_aa);

			const guint nobs = dxc->nobs;
			guint b;

			if (nobs > 1)
			{
				for (b = 0; b < nobs; b++)
				{
					if (a < b)
					{
						NcmVector* cl_th_0_ab = ncm_matrix_get_col (dxc->xcab[a][b]->cl_th, 0);
						NcmVector* cl_th_1_ab = ncm_matrix_get_col (dxc->xcab[a][b]->cl_th, 1);

						ncm_vector_scale (cl_th_0_ab, biasratio);
						ncm_vector_memcpy (cl_th_1_ab, cl_th_0_ab);

						ncm_vector_free(cl_th_0_ab);
						ncm_vector_free(cl_th_1_ab);
					}
					else if (b < a)
					{
						NcmVector* cl_th_0_ba = ncm_matrix_get_col (dxc->xcab[b][a]->cl_th, 0);
						NcmVector* cl_th_1_ba = ncm_matrix_get_col (dxc->xcab[b][a]->cl_th, 1);

						ncm_vector_scale (cl_th_0_ba, biasratio);
						ncm_vector_memcpy (cl_th_1_ba, cl_th_0_ba);

						ncm_vector_free(cl_th_0_ba);
						ncm_vector_free(cl_th_1_ba);
					}
				}
			}

			xclkg->bias_old = *(xclkg->bias);
			xclkg->noise_bias_old = ncm_vector_get (NCM_MODEL (xclkg)->params, NC_XCOR_LIMBER_KERNEL_GAL_NOISE_BIAS);

			return TRUE;
		}
		else
			return FALSE;
	}
	else
		return FALSE;
}

static void
_nc_data_xcor_prepare (NcmData* data, NcmMSet* mset)
{
	NcDataXcor* dxc = NC_DATA_XCOR (data);
	NcHICosmo* cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
	// NcHIReion* reion = NC_HIREION (ncm_mset_peek (mset, nc_hireion_id ()));
	const guint nobs = dxc->nobs;

	/* Array of booleans where prep[a][b] = TRUE means that the theoretical spectrum C_l^{a,b} needs to
	 * be recalculated
	 * Only the upper diagonal should be used !
	 */
	gboolean prep[NC_DATA_XCOR_MAX][NC_DATA_XCOR_MAX] = { FALSE };

	guint a, b;

	/* If cosmo has changed, prepare everything */
	if (ncm_model_ctrl_update (dxc->cosmo_ctrl, NCM_MODEL (cosmo)))
	{
		/* Prepare xcor */
		nc_xcor_prepare (dxc->xc, cosmo);

		/* Prepare the kernels */
		for (a = 0; a < nobs; a++)
		{
			NcXcorLimberKernel* xcl = NC_XCOR_LIMBER_KERNEL (ncm_mset_peek_pos (mset, nc_xcor_limber_kernel_id (), a));
			nc_xcor_limber_kernel_prepare (xcl, cosmo);
			for (b = a; b < nobs; b++)
			{
				// if (dxc->xcidx[a][b] > -1) All of the theoretical spectra must be computed since they can appear in the cross-terms of the covariance ! :)
				prep[a][b] = TRUE;
			}
		}
	}

	/* Otherwise check if one of the xclk's has been changed, reprepare it, and compute the auto and cross-spectra involving this variable only */
	else
	{
		// guint aa, bb;
		for (a = 0; a < nobs; a++)
		{
			NcmModelCtrl* ctrl = g_ptr_array_index (dxc->xclk_ctrl, a);
			NcXcorLimberKernel* xcl = NC_XCOR_LIMBER_KERNEL (ncm_mset_peek_pos (mset, nc_xcor_limber_kernel_id (), a));

			if (ncm_model_ctrl_update (ctrl, NCM_MODEL (xcl)))
			{
				/* If the observable is gal, then simply scale with the bias and add the noise_bias term*/
				if (!(_nc_data_xcor_fast_update (dxc, xcl, a, mset)))
				{
					nc_xcor_limber_kernel_prepare (xcl, cosmo);
					for (b = 0; b < nobs; b++)
					{
						// if (dxc->xcidx[a][b] > -1) All of the theoretical spectra must be computed since they can appear in the cross-terms of the covariance ! :)
						prep[a][b] = TRUE;
						prep[b][a] = TRUE;
					}
				}
			}
		}
	}

	/* Compute all the Cl's that need to be updated */
	for (a = 0; a < nobs; a++)
	{
		NcXcorLimberKernel* xcl1 = NC_XCOR_LIMBER_KERNEL (ncm_mset_peek_pos (mset, nc_xcor_limber_kernel_id (), a));

		if (prep[a][a])
		{
			NcmVector* cl_th_0_aa = ncm_matrix_get_col (dxc->xcab[a][a]->cl_th, 0);
			NcmVector* cl_th_1_aa = ncm_matrix_get_col (dxc->xcab[a][a]->cl_th, 1);

			// nc_xcor_limber (dxc->xc, xcl1, NULL, cosmo, 0, dxc->xcab[a][a]->ell_th_cut_off, cl_th_0_aa, NC_XCOR_LIMBER_METHOD_CVODE);
			nc_xcor_limber (dxc->xc, xcl1, NULL, cosmo, 0, dxc->xcab[a][a]->ell_th_cut_off, cl_th_0_aa, NC_XCOR_LIMBER_METHOD_GSL);

			nc_xcor_limber_kernel_add_noise (xcl1, cl_th_0_aa, cl_th_1_aa, 0);

			ncm_vector_free(cl_th_0_aa);
			ncm_vector_free(cl_th_1_aa);
		}

		if (nobs > 1)
		{
			for (b = a + 1; b < nobs; b++)
			{
				// printf("%i %i \n", a, b);
				if (prep[a][b])
				{
					NcXcorLimberKernel* xcl2 = NC_XCOR_LIMBER_KERNEL (ncm_mset_peek_pos (mset, nc_xcor_limber_kernel_id (), b));

					NcmVector* cl_th_0_ab = ncm_matrix_get_col (dxc->xcab[a][b]->cl_th, 0);
					NcmVector* cl_th_1_ab = ncm_matrix_get_col (dxc->xcab[a][b]->cl_th, 1);

					// nc_xcor_limber (dxc->xc, xcl1, xcl2, cosmo, 0, dxc->xcab[a][b]->ell_th_cut_off, cl_th_0_ab, NC_XCOR_LIMBER_METHOD_CVODE);
					nc_xcor_limber (dxc->xc, xcl1, xcl2, cosmo, 0, dxc->xcab[a][b]->ell_th_cut_off, cl_th_0_ab, NC_XCOR_LIMBER_METHOD_GSL);

					ncm_vector_memcpy (cl_th_1_ab, cl_th_0_ab);

					ncm_vector_free(cl_th_0_ab);
					ncm_vector_free(cl_th_1_ab);
				}
			}
		}
	}
}

/**
 * nc_data_xcor_mean_func_ab :
 * @dxc: a #NcDataXcor
 * @vp: a #NcmVector
 * @a: a #guint
 * @b: a #guint
 *
 * FIXME
 *
 */
void
nc_data_xcor_mean_func_ab (NcDataXcor* dxc, NcmVector* vp, guint a, guint b)
{
	guint ret;
	gint ell_idx;

	ell_idx = dxc->xcidx[a][b];
	if (ell_idx > -1)
	{
		NcXcorAB* xcab = dxc->xcab[a][b];

		// NcmVector* vp_sub = ncm_vector_get_subvector (vp, ell_idx, xcab->nell_lik);
		NcmVector* cl_th = ncm_matrix_get_col (xcab->cl_th, 0);
		NcmMatrix* mixing_sub = ncm_matrix_get_submatrix (xcab->mixing, xcab->ell_lik_min, 0, xcab->nell_lik, xcab->ell_th_cut_off + 1);

		ret = gsl_blas_dgemv (CblasNoTrans, 1.0, ncm_matrix_gsl (mixing_sub), ncm_vector_gsl (cl_th), 0.0, ncm_vector_gsl (vp));

		if (ret != 0)
			g_error ("_nc_data_xcor_compute_cl : %i", ret);

		// ncm_vector_free (vp_sub);
		ncm_vector_free (cl_th);
		ncm_matrix_free (mixing_sub);
	}
}

static void
_nc_data_xcor_mean (NcDataXcor* dxc)
{
	guint nobs = dxc->nobs;

	guint a, b, ret;
	gint ell_idx;

	for (a = 0; a < nobs; a++)
	{
		for (b = a; b < nobs; b++)
		{
			ell_idx = dxc->xcidx[a][b];
			if (ell_idx > -1)
			{
				NcXcorAB* xcab = dxc->xcab[a][b];

				NcmVector* pcl_sub = ncm_vector_get_subvector (dxc->pcl, ell_idx, xcab->nell_lik);
				NcmVector* cl_th = ncm_matrix_get_col (xcab->cl_th, 0);
				NcmMatrix* mixing_sub = ncm_matrix_get_submatrix (xcab->mixing, xcab->ell_lik_min, 0, xcab->nell_lik, xcab->ell_th_cut_off + 1);

				ret = gsl_blas_dgemv (CblasNoTrans, 1.0, ncm_matrix_gsl (mixing_sub), ncm_vector_gsl (cl_th), 0.0, ncm_vector_gsl (pcl_sub));

				if (ret != 0)
					g_error ("_nc_data_xcor_compute_cl : %i", ret);

				ncm_vector_free (pcl_sub);
				ncm_vector_free (cl_th);
				ncm_matrix_free (mixing_sub);
			}
		}
	}
}

static void
_nc_data_xcor_bin_vector (NcmVector* vf, NcmVector* vb, const guint dl)
{
	const guint lvf = ncm_vector_len (vf);
	const guint lvb = ncm_vector_len (vb);
	if (lvf != lvb * dl)
	{
		g_error ("_nc_data_xcor_bin_vector : length issue");
	}

	guint L, l;

	for (L = 0; L < lvb; L++)
	{
		gdouble vbL = 0.0;
		for (l = dl * L; l < dl * (L + 1); l++)
		{
			vbL += ncm_vector_get (vf, l);
		}
		ncm_vector_set (vb, L, vbL);
	}
}

static void
_nc_data_xcor_mean_func (NcmDataGaussCov* gauss, NcmMSet* mset, NcmVector* vp)
{
	NcDataXcor* dxc = NC_DATA_XCOR (gauss);

	_nc_data_xcor_mean (dxc); // compute pseudo-C_l's and store them in dxc->pcl

	_nc_data_xcor_bin_vector (dxc->pcl, vp, NC_DATA_XCOR_DL);
}

/**
 * nc_data_xcor_cov_func_abcd :
 * @dxc: a #NcDataXcor
 * @cov: a #NcmMatrix
 * @a: a #guint
 * @b: a #guint
 * @c: a #guint
 * @d: a #guint
 *
 * FIXME
 *
 */
void
nc_data_xcor_cov_func_abcd (NcDataXcor* dxc, NcmMatrix* cov, guint a, guint b, guint c, guint d)
{
	gint aa = -1, bb = -1, cc = -1, dd = -1;
	gint ell_idx_ab, ell_idx_cd;
	guint l, ll, i, j, ix, jx, lmin, llmin;
	gdouble res;

	ell_idx_ab = dxc->xcidx[a][b];
	ell_idx_cd = dxc->xcidx[c][d];

	lmin = dxc->xcab[a][b]->ell_lik_min;
	llmin = dxc->xcab[c][d]->ell_lik_min;

	for (l = lmin; l < dxc->xcab[a][b]->ell_lik_max + 1; l++)
	{
		i = l - lmin;
		ix = ell_idx_ab + i;
		for (ll = llmin; ll < dxc->xcab[c][d]->ell_lik_max + 1; ll++)
		{
			j = ll - llmin;
			jx = ell_idx_cd + j;
			if (i <= j)
			{
				res = 0.0;

				_nc_data_xcor_sort (a, d, &aa, &dd);
				_nc_data_xcor_sort (b, c, &bb, &cc);

				// printf ("%i %i %i %i %i %i %i %i \n", a, b, c, d, aa, bb, cc, dd);

				res = sqrt (fabs (ncm_matrix_get (dxc->xcab[aa][dd]->cl_th, l, 1) *
				                  ncm_matrix_get (dxc->xcab[aa][dd]->cl_th, ll, 1) *
				                  ncm_matrix_get (dxc->xcab[bb][cc]->cl_th, l, 1) *
				                  ncm_matrix_get (dxc->xcab[bb][cc]->cl_th, ll, 1))) *
				      ncm_matrix_get (dxc->X1, ix, jx);

				// printf ("%g %g %g %g %g %g", ncm_matrix_get (dxc->xcab[aa][dd]->cl_th, l, 1),
				//         ncm_matrix_get (dxc->xcab[aa][dd]->cl_th, ll, 1),
				//         ncm_matrix_get (dxc->xcab[bb][cc]->cl_th, l, 1),
				//         ncm_matrix_get (dxc->xcab[bb][cc]->cl_th, ll, 1),
				//         ncm_matrix_get (dxc->X1, ix, jx),
				//         res);

				_nc_data_xcor_sort (a, c, &aa, &cc);
				_nc_data_xcor_sort (b, d, &bb, &dd);

				// if ((a != b) && (c != d))
				// {
				res += sqrt (fabs (ncm_matrix_get (dxc->xcab[aa][cc]->cl_th, l, 1) *
				                   ncm_matrix_get (dxc->xcab[aa][cc]->cl_th, ll, 1) *
				                   ncm_matrix_get (dxc->xcab[bb][dd]->cl_th, l, 1) *
				                   ncm_matrix_get (dxc->xcab[bb][dd]->cl_th, ll, 1))) *
				       ncm_matrix_get (dxc->X2, ix, jx);
				// }

				// printf ("%g\n", res);
				ncm_matrix_set (cov, i, j, res);
				ncm_matrix_set (cov, j, i, res);
			}
		}
	}
}

// static gboolean
// _nc_data_xcor_cov (NcDataXcor* dxc)
static gboolean
_nc_data_xcor_cov_func (NcmDataGaussCov* gauss, NcmMSet* mset, NcmMatrix* cov)
{
	NcDataXcor* dxc = NC_DATA_XCOR (gauss);

	ncm_matrix_set_all (cov, 0.0);

	const guint nobs = dxc->nobs;

	guint a, b, c, d;
	gint aa = -1, bb = -1, cc = -1, dd = -1;
	gint ell_idx_ab, ell_idx_cd;
	guint l, ll, i, j, lmin, llmin, L, LL;
	gdouble res;

	for (a = 0; a < nobs; a++)
	{
		for (b = a; b < nobs; b++)
		{
			ell_idx_ab = dxc->xcidx[a][b];

			if (ell_idx_ab > -1)
			{
				for (c = 0; c < nobs; c++)
				{
					for (d = c; d < nobs; d++)
					{
						ell_idx_cd = dxc->xcidx[c][d];

						if (ell_idx_cd > -1)
						{
							lmin = dxc->xcab[a][b]->ell_lik_min;
							llmin = dxc->xcab[c][d]->ell_lik_min;

							for (l = lmin; l < dxc->xcab[a][b]->ell_lik_max + 1; l++)
							{
								i = ell_idx_ab + l - lmin;

								for (ll = llmin; ll < dxc->xcab[c][d]->ell_lik_max + 1; ll++)
								{
									j = ell_idx_cd + ll - llmin;

									L = i / NC_DATA_XCOR_DL;
									LL = j / NC_DATA_XCOR_DL;

									if (L <= LL)
									{
										res = 0.0;

										_nc_data_xcor_sort (a, d, &aa, &dd);
										_nc_data_xcor_sort (b, c, &bb, &cc);

										res = sqrt (fabs (ncm_matrix_get (dxc->xcab[aa][dd]->cl_th, l, 1) *
										                  ncm_matrix_get (dxc->xcab[aa][dd]->cl_th, ll, 1) *
										                  ncm_matrix_get (dxc->xcab[bb][cc]->cl_th, l, 1) *
										                  ncm_matrix_get (dxc->xcab[bb][cc]->cl_th, ll, 1))) *
										      ncm_matrix_get (dxc->X1, i, j);

										_nc_data_xcor_sort (a, c, &aa, &cc);
										_nc_data_xcor_sort (b, d, &bb, &dd);

										res += sqrt (fabs (ncm_matrix_get (dxc->xcab[aa][cc]->cl_th, l, 1) *
										                   ncm_matrix_get (dxc->xcab[aa][cc]->cl_th, ll, 1) *
										                   ncm_matrix_get (dxc->xcab[bb][dd]->cl_th, l, 1) *
										                   ncm_matrix_get (dxc->xcab[bb][dd]->cl_th, ll, 1))) *
										       ncm_matrix_get (dxc->X2, i, j);

										ncm_matrix_addto (cov, L, LL, res);
									}
								}
							}
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
 * @nobs: a #guint
 * @xc: a #NcXcor
 * @use_norma: a #gboolean
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcDataXcor *
nc_data_xcor_new_full (const guint nobs, NcXcor* xc, const gboolean use_norma) //, const gchar* xcname[])
{
	return g_object_new (NC_TYPE_DATA_XCOR,
	                     "use-norma", use_norma,
	                     "nobs", nobs,
	                     "xc", xc, /* ADD EVRYTHING HERE ! */
	                     NULL);
}

/**
 * nc_data_xcor_set_AB:
 * @dxc: a #NcDataXcor
 * @xcab: a #NcXcorAB
 *
 * FIXME
 *
 */
void
nc_data_xcor_set_AB (NcDataXcor* dxc, NcXcorAB* xcab)
{
	const guint a = xcab->a;
	const guint b = xcab->b;

	/* check if it's not an empty object (used only in covariance matrix)... */
	if (xcab->cl_obs != NULL)
	{

		// Set C_l counters
		dxc->xcidx[a][b] = dxc->xcidx_ctr;
		dxc->xcidx_ctr += xcab->nell_lik;

		// Set OA indices
		// dxc->xcab_oa_idx[dxc->xcab_oa_ctr][0] = a;
		// dxc->xcab_oa_idx[dxc->xcab_oa_ctr][1] = b;
		// dxc->xcab_oa_ctr += 1;


		// ncm_obj_array_add (dxc->xcab_oa, xcab);

		// ncm_matrix_free (mixing_full);

		//
		// // dxc->ncl += 1;
	}

	// Ref the xcor_AB obejct
	dxc->xcab[a][b] = nc_xcor_AB_ref (xcab);
}

static void
_nc_data_xcor_set_by_oa (NcDataXcor* dxc, NcmObjArray* oa)
{
	guint a, b;
	for (a = 0; a < NC_DATA_XCOR_MAX; a++)
	{
		for (b = 0; b < NC_DATA_XCOR_MAX; b++)
		{
			dxc->xcidx[a][b] = -1;
			nc_xcor_AB_clear(&dxc->xcab[a][b]);
		}
	}

	// Reset index for ObjectArray serialization
	// guint i, j;
	// for (i = 0; i < NC_DATA_XCOR_MAX * NC_DATA_XCOR_MAX; i++)
	// {
	// 	for (j = 0; j < 2; j++)
	// 	{
	// 		dxc->xcab_oa_idx[i][j] = 99;
	// 	}
	// }

	// Reset the counters
	dxc->xcidx_ctr = 0;
	// dxc->xcab_oa_ctr = 0;

	// Ref all the xcor_AB objects
	guint i;
	for (i = 0; i < oa->len; i++)
	{
		nc_data_xcor_set_AB (dxc, NC_XCOR_AB (ncm_obj_array_peek (oa, i)));
	}

	if (dxc->pcl == NULL)
	{
		const guint ctr = dxc->xcidx_ctr;

		dxc->pcl = ncm_vector_new (ctr);
		ncm_vector_set_zero (dxc->pcl);
		dxc->pcov = ncm_matrix_new (ctr, ctr);
		ncm_matrix_set_all (dxc->pcov, 0.0);
	}
}

/**
 * nc_data_xcor_set_3:
 * @dxc: a #NcDataXcor
 *
 * FIXME
 *
 */
void
nc_data_xcor_set_3 (NcDataXcor* dxc)
{
	const guint ctr = dxc->xcidx_ctr;

	dxc->X1 = ncm_matrix_new (ctr, ctr);
	ncm_matrix_set_all (dxc->X1, 0.);
	dxc->X2 = ncm_matrix_new (ctr, ctr);
	ncm_matrix_set_all (dxc->X2, 0.);

	dxc->pcl = ncm_vector_new (ctr);
	ncm_vector_set_zero (dxc->pcl);
	dxc->pcov = ncm_matrix_new (ctr, ctr);
	ncm_matrix_set_all (dxc->pcov, 0.0);

	/* Set the parent instance : size, model vector and data vector*/
	NcmDataGaussCov* gauss = NCM_DATA_GAUSS_COV (dxc);
	const guint np = ctr / NC_DATA_XCOR_DL;
	gauss->np = np;
	gauss->v = ncm_vector_new (np);
	gauss->cov = ncm_matrix_new (np, np);
	gauss->y = ncm_vector_new (np);
	ncm_vector_set_zero (gauss->v);
	ncm_vector_set_zero (gauss->y);
	ncm_matrix_set_zero (gauss->cov);

	// ncm_data_gauss_cov_set_size (gauss, np);

	gint ell_idx;
	guint a, b;
	guint nobs = dxc->nobs;

	// Put the binned mean value (temporarilly borrowing pcl), find max of ell_lik_max
	guint max_ell_lik_max = 0;
	for (a = 0; a < nobs; a++)
	{
		for (b = a; b < nobs; b++)
		{
			ell_idx = dxc->xcidx[a][b];
			if (ell_idx > -1)
			{
				NcXcorAB* xcab = dxc->xcab[a][b];

				ncm_vector_memcpy2 (dxc->pcl,
				                    xcab->cl_obs,
				                    ell_idx,
				                    xcab->ell_lik_min,
				                    xcab->nell_lik);

				max_ell_lik_max = GSL_MAX (max_ell_lik_max, xcab->ell_lik_max);
			}
		}
	}

	// Create missing xcab that may be used in the covariance, and check that for the other xcab, ell_th_cut_off >= max_ell_lik_max
	for (a = 0; a < nobs; a++)
	{
		for (b = a; b < nobs; b++)
		{
			ell_idx = dxc->xcidx[a][b];
			if (ell_idx > -1)
			{
				g_assert_cmpuint (dxc->xcab[a][b]->ell_th_cut_off, >=, max_ell_lik_max);
			}
			else // C_l^{a,b} not used
			{
				dxc->xcab[a][b] = nc_xcor_AB_new (a, b, max_ell_lik_max, 0, 0, NULL, NULL, 0);
			}
		}
	}

	_nc_data_xcor_bin_vector (dxc->pcl, gauss->y, NC_DATA_XCOR_DL);

	ncm_vector_set_zero (dxc->pcl);
}

/**
 * nc_data_xcor_set_4:
 * @dxc: a #NcDataXcor
 * @a: a #guint
 * @b: a #guint
 * @c: a #guint
 * @d: a #guint
 * @X1_filename: a #gchar
 * @X2_filename: a #gchar
 * @X_filelength: a #guint
 *
 * FIXME
 *
 */
void
nc_data_xcor_set_4 (NcDataXcor* dxc, guint a, guint b, guint c, guint d, const gchar* X1_filename, const gchar* X2_filename, guint X_filelength)
{
	if ((b < a) | (d < c))
		g_error ("b,d must be greater or equal to a,c");

	// NcXcorAB* xcab = ncm_obj_array_peek (dxc->xcab_oa, dxc->xcab_oa_idx[a][b])
	// NcXcorAB* xccd = ncm_obj_array_peek (dxc->xcab_oa, dxc->xcab_oa_idx[c][d])
	NcXcorAB* xcab = dxc->xcab[a][b];
	NcXcorAB* xccd = dxc->xcab[c][d];

	NcmMatrix* X1f = ncm_matrix_new (X_filelength, X_filelength);
	FILE* f1 = fopen (X1_filename, "r");
	gsl_matrix_fscanf (f1, ncm_matrix_gsl (X1f));
	fclose (f1);

	NcmMatrix* X1f_sub = ncm_matrix_get_submatrix (X1f, xcab->ell_lik_min, xccd->ell_lik_min, xcab->nell_lik, xccd->nell_lik);
	// printf ("%i %i %i %i \n", dxc->xcab[a][b]->ell_lik_min, dxc->xcab[c][d]->ell_lik_min, dxc->xcab[a][b]->nell_lik, dxc->xcab[c][d]->nell_lik);
	NcmMatrix* X1_sub = ncm_matrix_get_submatrix (dxc->X1, dxc->xcidx[a][b], dxc->xcidx[c][d], xcab->nell_lik, xccd->nell_lik);
	// printf ("%i %i %i %i \n", dxc->xcidx[a][b], dxc->xcidx[c][d], dxc->xcab[a][b]->nell_lik, dxc->xcab[c][d]->nell_lik);

	ncm_matrix_memcpy (X1_sub, X1f_sub);

	ncm_matrix_free (X1f);
	ncm_matrix_free (X1f_sub);
	ncm_matrix_free (X1_sub);


	NcmMatrix* X2f = ncm_matrix_new (X_filelength, X_filelength);
	FILE* f2 = fopen (X2_filename, "r");
	gsl_matrix_fscanf (f2, ncm_matrix_gsl (X2f));
	fclose (f2);

	NcmMatrix* X2f_sub = ncm_matrix_get_submatrix (X2f, xcab->ell_lik_min, xccd->ell_lik_min, xcab->nell_lik, xccd->nell_lik);
	// printf ("%i %i %i %i \n", dxc->xcab[a][b]->ell_lik_min, dxc->xcab[c][d]->ell_lik_min, dxc->xcab[a][b]->nell_lik, dxc->xcab[c][d]->nell_lik);
	NcmMatrix* X2_sub = ncm_matrix_get_submatrix (dxc->X2, dxc->xcidx[a][b], dxc->xcidx[c][d], xcab->nell_lik, xccd->nell_lik);
	// printf ("%i %i %i %i \n", dxc->xcidx[a][b], dxc->xcidx[c][d], dxc->xcab[a][b]->nell_lik, dxc->xcab[c][d]->nell_lik);

	ncm_matrix_memcpy (X2_sub, X2f_sub);

	ncm_matrix_free (X2f);
	ncm_matrix_free (X2f_sub);
	ncm_matrix_free (X2_sub);
	// dxc->X1[a][b][c][d] = ncm_matrix_ref (X1);
	// dxc->X1[c][d][a][b] = ncm_matrix_ref (X1);
	//
	// dxc->X2[a][b][c][d] = ncm_matrix_ref (X2);
	// dxc->X2[c][d][a][b] = ncm_matrix_ref (X2);
}

/**
 * nc_data_xcor_set_5:
 * @dxc: a #NcDataXcor
 *
 * FIXME
 *
 */
void
nc_data_xcor_set_5 (NcDataXcor* dxc)
{
	NcmData *data = NCM_DATA (dxc);
	data->init = TRUE;
}

/**
 * nc_data_xcor_get_cl_obs:
 * @dxc: a #NcDataXcor
 * @vp: a #NcmVector
 * @a: a #guint
 * @b: a #guint
 *
 * FIXME
 *
 */
void
nc_data_xcor_get_cl_obs (NcDataXcor* dxc, NcmVector* vp, guint a, guint b)
{
	ncm_vector_memcpy2 (vp, dxc->xcab[a][b]->cl_obs, 0, dxc->xcab[a][b]->ell_lik_min, dxc->xcab[a][b]->nell_lik);
}
