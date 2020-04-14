/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_spline_rbf.c
 *
 *  Fri April 06 20:44:42 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_spline_rbf.c
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

/**
 * SECTION:ncm_spline_rbf
 * @title: NcmSplineRBF
 * @short_description: Radial Basis Function implementation of spline class.
 *
 * This object implements #NcmSpline, using Radial Basis Function (RBF) methods.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline_rbf.h"
#include "math/ncm_matrix.h"
#include "ncm_enum_types.h"

#include <gsl/gsl_min.h>
#include <nvector/nvector_serial.h>

struct _NcmSplineRBFPrivate
{
	NcmSplineRBFType type_id;
	gchar *name;
	NcmMatrix *interp_matrix;
	NcmVector *coeff;
	NcmVector *shape_params;
	void (*default_shape_params) (NcmSplineRBFPrivate * const self, NcmVector *xv);
	void (*prepare_matrix) (NcmSplineRBFPrivate * const self, NcmVector *xv, NcmMatrix *interp_matrix, const guint len);
	gdouble (*eval) (NcmSplineRBFPrivate * const self, NcmVector *xv, const gdouble x);
	gdouble (*deriv) (NcmSplineRBFPrivate * const self, NcmVector *xv, const gdouble x);
	gdouble (*deriv2) (NcmSplineRBFPrivate * const self, NcmVector *xv, const gdouble x);
	gdouble (*integ) (NcmSplineRBFPrivate * const self, NcmVector *xv, const gdouble x0, const gdouble x1);
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmSplineRBF, ncm_spline_rbf, NCM_TYPE_SPLINE);

enum
{
  PROP_0,
  PROP_TYPE_ID,
	PROP_SHAPE_PARAMS,
  PROP_SIZE,
};

static void _ncm_spline_rbf_default_shape_params (NcmSplineRBFPrivate * const self, NcmVector *xv) { g_error ("_ncm_spline_rbf_default_shape_params: method not set."); }
static void _ncm_spline_rbf_type_prepare_matrix (NcmSplineRBFPrivate * const self, NcmVector *xv, NcmMatrix *interp_matrix, const guint len) { g_error ("_ncm_spline_rbf_type_prepare_matrix: method not set."); }
static gdouble _ncm_spline_rbf_type_eval (NcmSplineRBFPrivate * const self, NcmVector *xv, const gdouble x) { g_error ("_ncm_spline_rbf_type_eval: method not set."); return 0.0; }
static gdouble _ncm_spline_rbf_type_deriv (NcmSplineRBFPrivate * const self, NcmVector *xv, const gdouble x) { g_error ("_ncm_spline_rbf_type_deriv: method not set."); return 0.0; }
static gdouble _ncm_spline_rbf_type_deriv2 (NcmSplineRBFPrivate * const self, NcmVector *xv, const gdouble x) { g_error ("_ncm_spline_rbf_type_deriv2: method not set."); return 0.0; }

static void
ncm_spline_rbf_init (NcmSplineRBF *rbf)
{
	NcmSplineRBFPrivate * const self = rbf->priv = ncm_spline_rbf_get_instance_private (rbf);

	self->type_id        = NCM_SPLINE_RBF_TYPE_LEN;
	self->name           = NULL;
	self->interp_matrix  = NULL;
	self->coeff          = NULL;
	self->shape_params   = NULL;
	
	self->default_shape_params = &_ncm_spline_rbf_default_shape_params;
	self->prepare_matrix       = &_ncm_spline_rbf_type_prepare_matrix;
	self->eval                 = &_ncm_spline_rbf_type_eval;
	self->deriv                = &_ncm_spline_rbf_type_deriv;
	self->deriv2               = &_ncm_spline_rbf_type_deriv2;
	
}

static void
_ncm_spline_rbf_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSplineRBF *rbf = NCM_SPLINE_RBF (object);
  g_return_if_fail (NCM_IS_SPLINE_RBF (object));

  switch (prop_id)
  {
    case PROP_TYPE_ID:
      ncm_spline_rbf_set_type (rbf, g_value_get_enum (value));
      break;
    case PROP_SHAPE_PARAMS:
      ncm_spline_rbf_set_shape_params (rbf, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_spline_rbf_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSplineRBF *rbf = NCM_SPLINE_RBF (object);
	NcmSplineRBFPrivate * const self = rbf->priv;
  g_return_if_fail (NCM_IS_SPLINE_RBF (object));

  switch (prop_id)
  {
    case PROP_TYPE_ID:
    {
      g_value_set_enum (value, self->type_id);
      break;
    }
		case PROP_SHAPE_PARAMS:
			g_value_set_object (value, self->shape_params);
			break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_spline_rbf_dispose (GObject *object)
{
	NcmSplineRBF *rbf = NCM_SPLINE_RBF (object);
	NcmSplineRBFPrivate * const self = rbf->priv;
	
	ncm_matrix_clear (&self->interp_matrix);
	ncm_vector_clear (&self->coeff);
	ncm_vector_clear (&self->shape_params);
	
	/* Chain up : end */
	G_OBJECT_CLASS (ncm_spline_rbf_parent_class)->dispose (object);
}

static void
_ncm_spline_rbf_finalize (GObject *object)
{
	NcmSplineRBF *rbf = NCM_SPLINE_RBF (object);
	NcmSplineRBFPrivate * const self = rbf->priv;
	
	g_clear_pointer (&self->name, g_free);

	/* Chain up : end */
	G_OBJECT_CLASS (ncm_spline_rbf_parent_class)->finalize (object);
}

static const gchar *_ncm_spline_rbf_name (NcmSpline *s);
static void _ncm_spline_rbf_reset (NcmSpline *s);
static void _ncm_spline_rbf_prepare (NcmSpline *s);
static gsize _ncm_spline_rbf_min_size (const NcmSpline *s);
static gdouble _ncm_spline_rbf_eval (const NcmSpline *s, const gdouble x);
static gdouble _ncm_spline_rbf_deriv (const NcmSpline *s, const gdouble x);
static gdouble _ncm_spline_rbf_deriv2 (const NcmSpline *s, const gdouble x);
/*static gdouble _ncm_spline_rbf_deriv_nmax (const NcmSpline *s, const gdouble x);*/
static gdouble _ncm_spline_rbf_integ (const NcmSpline *s, const gdouble x0, const gdouble x1);
static NcmSpline *_ncm_spline_rbf_copy_empty (const NcmSpline *s);

static void
ncm_spline_rbf_class_init (NcmSplineRBFClass *klass)
{
	GObjectClass* object_class = G_OBJECT_CLASS (klass);
	NcmSplineClass* s_class = NCM_SPLINE_CLASS (klass);

	object_class->set_property = &_ncm_spline_rbf_set_property;
  object_class->get_property = &_ncm_spline_rbf_get_property;
	object_class->dispose      = &_ncm_spline_rbf_dispose;
	object_class->finalize     = &_ncm_spline_rbf_finalize;

	g_object_class_install_property (object_class,
	                                 PROP_TYPE_ID,
	                                 g_param_spec_enum ("type-id",
	                                                    NULL,
	                                                    "Type ID",
	                                                    NCM_TYPE_SPLINE_RBF_TYPE, NCM_SPLINE_RBF_TYPE_POSDEF_GAUSS,
	                                                    G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	
	
	s_class->name         = &_ncm_spline_rbf_name;
	s_class->prepare      = &_ncm_spline_rbf_prepare;
	s_class->prepare_base = NULL;
	s_class->min_size     = &_ncm_spline_rbf_min_size;
	s_class->reset        = &_ncm_spline_rbf_reset;
	s_class->eval         = &_ncm_spline_rbf_eval;
	s_class->deriv        = &_ncm_spline_rbf_deriv;
	s_class->deriv2       = &_ncm_spline_rbf_deriv2;
  s_class->deriv_nmax   = &_ncm_spline_rbf_deriv2;
	s_class->integ        = &_ncm_spline_rbf_integ;
	s_class->copy_empty   = &_ncm_spline_rbf_copy_empty;	
}

static const gchar *
_ncm_spline_rbf_name (NcmSpline *s)
{
	NcmSplineRBF *rbf = NCM_SPLINE_RBF (s);
	NcmSplineRBFPrivate * const self = rbf->priv;

	return self->name;
}

static void 
_ncm_spline_rbf_reset (NcmSpline *s)
{
	NcmSplineRBF *rbf = NCM_SPLINE_RBF (s);
	NcmSplineRBFPrivate * const self = rbf->priv;

	if ((self->coeff == NULL) || (ncm_vector_len (self->coeff) != s->len))
	{
		ncm_vector_clear (&self->coeff);
		self->coeff = ncm_vector_new (s->len);
	}	

	ncm_vector_set_all (self->shape_params, GSL_NAN);
}

gint 
ATimes (gpointer user_data, N_Vector v_vec, N_Vector z_vec)
{
	NcmSplineRBFPrivate * const self = user_data;
	const gint len = ncm_matrix_nrows (self->interp_matrix);
	gdouble *v, *z;

	v = N_VGetArrayPointer (v_vec);
	z = N_VGetArrayPointer (z_vec);
	
	cblas_dsymv (CblasRowMajor, CblasLower, len, 
	             1.0, 
	             ncm_matrix_data (self->interp_matrix), ncm_matrix_tda (self->interp_matrix),
	             v, 1, 
	             0.0,
	             z, 1);

	return 0;
}

gint 
PSolve (gpointer user_data, N_Vector r_vec, N_Vector z_vec, realtype tol, int lr)
{
	NcmSplineRBFPrivate * const self = user_data;
	const gint len = ncm_matrix_nrows (self->interp_matrix);
  gdouble *r, *z;
  gint i;
  
  r = N_VGetArrayPointer (r_vec);
  z = N_VGetArrayPointer (z_vec);

  for (i = 0; i < len; i++) 
	{
		const gdouble Aiim1 = /*sqrt*/ (1.0 / ncm_matrix_get (self->interp_matrix, i, i));
    z[i] = r[i] * Aiim1;
	}

  return 0;
}

gint PSetup (gpointer Data) { return 0; }

static void
_ncm_spline_rbf_calc_coeff (NcmSplineRBFPrivate * const self, NcmVector *xv, NcmVector *yv, const guint len)
{
	if (TRUE)
	{
		gint ret;

		self->prepare_matrix (self, xv, self->interp_matrix, len);
		ncm_vector_memcpy (self->coeff, yv);
		ret = ncm_matrix_cholesky_solve (self->interp_matrix, self->coeff, 'L');
		if (ret != 0)
			g_error ("_ncm_spline_rbf_prepare[ncm_matrix_cholesky_solve]: %d.", ret);
	}
	else
	{
		NcmVector *Q = ncm_vector_new (len);
		self->prepare_matrix (self, xv, self->interp_matrix, len);

		ncm_vector_set_zero (Q);
		ncm_vector_memcpy (self->coeff, yv);
		
		cblas_dsymv (CblasRowMajor, CblasLower, len, 
		             1.0, 
		             ncm_matrix_data (self->interp_matrix), ncm_matrix_tda (self->interp_matrix),
		             ncm_vector_data (self->coeff), 1, 
		             0.0,
		             ncm_vector_data (Q), 1);

		
		
		
	}
}

static gdouble 
_ncm_spline_rbf_LOOCV_err2 (gdouble h, gpointer user_data)
{
	NcmSpline *s = user_data;
	NcmSplineRBF *rbf = NCM_SPLINE_RBF (s);
	NcmSplineRBFPrivate * const self = rbf->priv;
	
	const guint len = ncm_spline_get_len (s);
	gdouble err     = 0.0;
	guint i;
	gint ret;

	ncm_vector_set (self->shape_params, 0, h);
	self->prepare_matrix (self, s->xv, self->interp_matrix, len);

	ret = ncm_matrix_cholesky_decomp (self->interp_matrix, 'L');
	
	if (ret != 0)
	{
		/*printf ("h % 22.15g err % 22.15e\n", h, 1.0e20);*/
		return 1.0e200;
	}

	ncm_vector_memcpy (self->coeff, s->yv);
	ret = ncm_matrix_cholesky_solve2 (self->interp_matrix, self->coeff, 'L');
	if (ret != 0)
		g_error ("_ncm_spline_rbf_LOOCV_err2[ncm_matrix_cholesky_solve]: %d.", ret);

	ncm_matrix_cholesky_inverse (self->interp_matrix, 'L');
	if (ret != 0)
		g_error ("_ncm_spline_rbf_LOOCV_err2[ncm_matrix_cholesky_decomp]: %d.", ret);

	for (i = 0; i < len; i++)
	{
		/*err += gsl_pow_2 (ncm_vector_get (self->coeff, i) / (ncm_matrix_get (self->interp_matrix, i, i) * ncm_vector_get (s->yv, i)));*/
		err += gsl_pow_2 (ncm_vector_get (self->coeff, i) / (ncm_matrix_get (self->interp_matrix, i, i)));
	}

	/*printf ("h % 22.15g err % 22.15e\n", h, err);*/
	return err;	
}

static void 
_ncm_spline_rbf_prepare (NcmSpline *s)
{
	NcmSplineRBF *rbf = NCM_SPLINE_RBF (s);
	NcmSplineRBFPrivate * const self = rbf->priv;
	const guint len = ncm_spline_get_len (s);

	if (len == 0)
		g_error ("_ncm_spline_rbf_prepare: empty spline.");

	if (!ncm_vector_is_finite (self->shape_params))
		self->default_shape_params (self, s->xv);

	self->interp_matrix = ncm_matrix_new (len, len);

	if (FALSE)
	{
		gsl_min_fminimizer *fmin = gsl_min_fminimizer_alloc (gsl_min_fminimizer_brent);
		gdouble ht = ncm_vector_get (self->shape_params, 0);
		gdouble hi = ht / 1000.0;
		gdouble hf = ht * 1000.0;
		gdouble last_hi = hi;
		gdouble last_hf = hf;
		guint iter = 0;
		gint status;
		gsl_function F;

		F.params   = s;
		F.function = &_ncm_spline_rbf_LOOCV_err2;

		gsl_min_fminimizer_set (fmin, &F, ht, hi, hf);
		do {
			iter++;
			status = gsl_min_fminimizer_iterate (fmin);
			if (status)
				g_error ("_ncm_spline_rbf_prepare: Cannot find minimum (%s)", gsl_strerror (status));

			ht = gsl_min_fminimizer_x_minimum (fmin);
			hi = gsl_min_fminimizer_x_lower (fmin);
			hf = gsl_min_fminimizer_x_upper (fmin);

			status = gsl_min_test_interval (hi, hf, 0.0, 1.0e-2);

			if ((status == GSL_CONTINUE) && (hi == last_hi) && (hf == last_hf))
			{
				g_warning ("_ncm_spline_rbf_prepare: minimization not improving, giving up...");
				break;
			}

			printf ("[%d] % 22.15e % 22.15e % 22.15e\n", status, ht, hi, hf);

			last_hi = hi;
			last_hf = hf;

		} while (status == GSL_CONTINUE && iter < 1000);

		ncm_vector_set (self->shape_params, 0, gsl_min_fminimizer_x_minimum (fmin));
		gsl_min_fminimizer_free (fmin);
	}


	_ncm_spline_rbf_calc_coeff (self, s->xv, s->yv, len);
	
	ncm_matrix_clear (&self->interp_matrix);
}

static gsize 
_ncm_spline_rbf_min_size (const NcmSpline *s)
{
	return 1;
}

static gdouble 
_ncm_spline_rbf_eval (const NcmSpline *s, const gdouble x)
{
	NcmSplineRBF *rbf = NCM_SPLINE_RBF (s);
	NcmSplineRBFPrivate * const self = rbf->priv;

	return self->eval (self, s->xv, x);
}

static gdouble 
_ncm_spline_rbf_deriv (const NcmSpline *s, const gdouble x)
{
	NcmSplineRBF *rbf = NCM_SPLINE_RBF (s);
	NcmSplineRBFPrivate * const self = rbf->priv;

	return self->deriv (self, s->xv, x);	
}

static gdouble 
_ncm_spline_rbf_deriv2 (const NcmSpline *s, const gdouble x)
{
	NcmSplineRBF *rbf = NCM_SPLINE_RBF (s);
	NcmSplineRBFPrivate * const self = rbf->priv;

	return self->deriv2 (self, s->xv, x);	
}

static gdouble _ncm_spline_rbf_integ (const NcmSpline *s, const gdouble x0, const gdouble x1)
{
	NcmSplineRBF *rbf = NCM_SPLINE_RBF (s);
	NcmSplineRBFPrivate * const self = rbf->priv;

	return self->integ (self, s->xv, x0, x1);		
}

static NcmSpline *
_ncm_spline_rbf_copy_empty (const NcmSpline *s)
{
	NcmSplineRBF *rbf = NCM_SPLINE_RBF (s);
	NcmSplineRBFPrivate * const self = rbf->priv;
	return NCM_SPLINE (ncm_spline_rbf_new (self->type_id));
}

static void 
_ncm_spline_rbf_posdef_gauss_default_shape_params (NcmSplineRBFPrivate * const self, NcmVector *xv)
{
	const guint len  = ncm_vector_len (xv);
	const gdouble xi = ncm_vector_get (xv, 0);
	const gdouble xf = ncm_vector_get (xv, len - 1);
	
	ncm_vector_set (self->shape_params, 0, 1.0 / (2.5 * (xf - xi) / (1.0 * len)));
}

static gdouble
_ncm_spline_rbf_posdef_gauss_kern (const gdouble x1, const gdouble x2, const gdouble h)
{
	return x1 * x2 * exp (- 0.5 * gsl_pow_2 ((x1 - x2) * h));
}

static gdouble
_ncm_spline_rbf_posdef_gauss_kern_deriv (const gdouble x1, const gdouble x2, const gdouble h)
{
	return x2 * (1.0 - x1 * (x1 - x2) * h * h) * exp (- 0.5 * gsl_pow_2 ((x1 - x2) * h));
}

static gdouble
_ncm_spline_rbf_posdef_gauss_kern_deriv2 (const gdouble x1, const gdouble x2, const gdouble h)
{
	const gdouble h2 = h * h;
	return h2 * x2 * ( x1 * (h2 * gsl_pow_2 (x1 - x2) - 3.0) + 2.0 * x2) * exp (- 0.5 * gsl_pow_2 ((x1 - x2) * h));
}

static void 
_ncm_spline_rbf_posdef_gauss_prepare_matrix (NcmSplineRBFPrivate * const self, NcmVector *xv, NcmMatrix *interp_matrix, const guint len)
{ 
	const gdouble h = ncm_vector_get (self->shape_params, 0);
	guint i;

	for (i = 0; i < len; i++)
	{
		const gdouble x1 = ncm_vector_get (xv, i); 
		guint j;

		for (j = i; j < len; j++)
		{
			const gdouble x2    = ncm_vector_get (xv, j); 
			const gdouble Kx1x2 = _ncm_spline_rbf_posdef_gauss_kern (x1, x2, h);
			ncm_matrix_set (interp_matrix, j, i, Kx1x2);
		}
	}
}

static gdouble 
_ncm_spline_rbf_posdef_gauss_eval (NcmSplineRBFPrivate * const self, NcmVector *xv, const gdouble x)
{ 
	const gdouble h = ncm_vector_get (self->shape_params, 0);
	const guint len = ncm_vector_len (self->coeff);
	gdouble res = 0.0;
	guint i;

	for (i = 0; i < len; i++)
	{
		const gdouble xi = ncm_vector_get (xv, i); 
		const gdouble ci = ncm_vector_get (self->coeff, i); 
		res += ci * _ncm_spline_rbf_posdef_gauss_kern (x, xi, h);
	}

	return res;
}

static gdouble 
_ncm_spline_rbf_posdef_gauss_deriv (NcmSplineRBFPrivate * const self, NcmVector *xv, const gdouble x)
{ 
	const gdouble h = ncm_vector_get (self->shape_params, 0);
	const guint len = ncm_vector_len (self->coeff);
	gdouble res = 0.0;
	guint i;

	for (i = 0; i < len; i++)
	{
		const gdouble xi = ncm_vector_get (xv, i); 
		const gdouble ci = ncm_vector_get (self->coeff, i); 
		res += ci * _ncm_spline_rbf_posdef_gauss_kern_deriv (x, xi, h);
	}

	return res;
}

static gdouble 
_ncm_spline_rbf_posdef_gauss_deriv2 (NcmSplineRBFPrivate * const self, NcmVector *xv, const gdouble x)
{ 
	const gdouble h = ncm_vector_get (self->shape_params, 0);
	const guint len = ncm_vector_len (self->coeff);
	gdouble res = 0.0;
	guint i;

	for (i = 0; i < len; i++)
	{
		const gdouble xi = ncm_vector_get (xv, i); 
		const gdouble ci = ncm_vector_get (self->coeff, i); 
		res += ci * _ncm_spline_rbf_posdef_gauss_kern_deriv2 (x, xi, h);
	}

	return res;
}

static gdouble 
_ncm_spline_rbf_integ_placeholder (NcmSplineRBFPrivate * const self, NcmVector *xv, const gdouble x0, const gdouble x1)
{
	g_error ("_ncm_spline_rbf_integ_placeholder: method not implemented!");
	return 0.0;
}

static void 
_ncm_spline_rbf_gauss_default_shape_params (NcmSplineRBFPrivate * const self, NcmVector *xv)
{
	const guint len  = ncm_vector_len (xv);
	const gdouble xi = ncm_vector_get (xv, 0);
	const gdouble xf = ncm_vector_get (xv, len - 1);
	
	ncm_vector_set (self->shape_params, 0, 1.0 / (2.5 * (xf - xi) / (1.0 * len)));
}

static gdouble
_ncm_spline_rbf_gauss_kern (const gdouble x1, const gdouble x2, const gdouble h)
{
	return exp (- 0.5 * gsl_pow_2 ((x1 - x2) * h));
}

static gdouble
_ncm_spline_rbf_gauss_kern_deriv (const gdouble x1, const gdouble x2, const gdouble h)
{
	return - (x1 - x2) * h * h * exp (- 0.5 * gsl_pow_2 ((x1 - x2) * h));
}

static gdouble
_ncm_spline_rbf_gauss_kern_deriv2 (const gdouble x1, const gdouble x2, const gdouble h)
{
	const gdouble h2 = h * h;
	return h2 * (-1.0 + gsl_pow_2 ((x1 - x2) * h)) * exp (- 0.5 * gsl_pow_2 ((x1 - x2) * h));
}

static void 
_ncm_spline_rbf_gauss_prepare_matrix (NcmSplineRBFPrivate * const self, NcmVector *xv, NcmMatrix *interp_matrix, const guint len)
{ 
	const gdouble h = ncm_vector_get (self->shape_params, 0);
	guint i;

	for (i = 0; i < len; i++)
	{
		const gdouble x1 = ncm_vector_get (xv, i); 
		guint j;

		for (j = i; j < len; j++)
		{
			const gdouble x2    = ncm_vector_get (xv, j); 
			const gdouble Kx1x2 = _ncm_spline_rbf_gauss_kern (x1, x2, h);
			ncm_matrix_set (interp_matrix, j, i, Kx1x2);
		}
	}
}

static gdouble 
_ncm_spline_rbf_gauss_eval (NcmSplineRBFPrivate * const self, NcmVector *xv, const gdouble x)
{ 
	const gdouble h = ncm_vector_get (self->shape_params, 0);
	const guint len = ncm_vector_len (self->coeff);
	gdouble res = 0.0;
	guint i;

	for (i = 0; i < len; i++)
	{
		const gdouble xi = ncm_vector_get (xv, i); 
		const gdouble ci = ncm_vector_get (self->coeff, i); 
		res += ci * _ncm_spline_rbf_gauss_kern (x, xi, h);
	}

	return res;
}

static gdouble 
_ncm_spline_rbf_gauss_deriv (NcmSplineRBFPrivate * const self, NcmVector *xv, const gdouble x)
{ 
	const gdouble h = ncm_vector_get (self->shape_params, 0);
	const guint len = ncm_vector_len (self->coeff);
	gdouble res = 0.0;
	guint i;

	for (i = 0; i < len; i++)
	{
		const gdouble xi = ncm_vector_get (xv, i); 
		const gdouble ci = ncm_vector_get (self->coeff, i); 
		res += ci * _ncm_spline_rbf_gauss_kern_deriv (x, xi, h);
	}

	return res;
}

static gdouble 
_ncm_spline_rbf_gauss_deriv2 (NcmSplineRBFPrivate * const self, NcmVector *xv, const gdouble x)
{ 
	const gdouble h = ncm_vector_get (self->shape_params, 0);
	const guint len = ncm_vector_len (self->coeff);
	gdouble res = 0.0;
	guint i;

	for (i = 0; i < len; i++)
	{
		const gdouble xi = ncm_vector_get (xv, i); 
		const gdouble ci = ncm_vector_get (self->coeff, i); 
		res += ci * _ncm_spline_rbf_gauss_kern_deriv2 (x, xi, h);
	}

	return res;
}

/**
 * ncm_spline_rbf_new:
 * @type_id: a #NcmSplineRBFType
 * 
 * Creates a new RBF using @type_id. 
 * 
 * Returns: (transfer full): a newly created #NcmSplineRBF.
 */ 
NcmSplineRBF *
ncm_spline_rbf_new (NcmSplineRBFType type_id)
{
	NcmSplineRBF *rbf = g_object_new (NCM_TYPE_SPLINE_RBF, 
	                                  "type-id", type_id,
	                                  NULL);
	return rbf;
}

/**
 * ncm_spline_rbf_ref:
 * @rbf: a #NcmSplineRBF
 *
 * Increase the reference of @rbf by one.
 *
 * Returns: (transfer full): @rbf.
 */
NcmSplineRBF *
ncm_spline_rbf_ref (NcmSplineRBF *rbf)
{
  return g_object_ref (rbf);
}

/**
 * ncm_spline_rbf_free:
 * @rbf: a #NcmSplineRBF
 *
 * Decrease the reference count of @rbf by one.
 *
 */
void
ncm_spline_rbf_free (NcmSplineRBF *rbf)
{
  g_object_unref (rbf);
}

/**
 * ncm_spline_rbf_clear:
 * @rbf: a #NcmSplineRBF
 *
 * Decrease the reference count of @rbf by one, and sets the pointer *@rbf to
 * NULL.
 *
 */
void
ncm_spline_rbf_clear (NcmSplineRBF **rbf)
{
  g_clear_object (rbf);
}

/**
 * ncm_spline_rbf_set_type:
 * @rbf: a #NcmSplineRBF
 * @type_id: a #NcmSplineRBFType
 * 
 * Sets the RBF type function to @type_id. 
 * 
 */ 
void
ncm_spline_rbf_set_type (NcmSplineRBF *rbf, NcmSplineRBFType type_id)
{
	NcmSplineRBFPrivate * const self = rbf->priv;
	
	if (type_id != self->type_id)
	{
		g_clear_pointer (&self->name, g_free);
		ncm_vector_clear (&self->shape_params);
		
		switch (type_id)
		{
			case NCM_SPLINE_RBF_TYPE_POSDEF_GAUSS:
				self->default_shape_params = &_ncm_spline_rbf_posdef_gauss_default_shape_params;
				self->prepare_matrix       = &_ncm_spline_rbf_posdef_gauss_prepare_matrix;
				self->eval                 = &_ncm_spline_rbf_posdef_gauss_eval;
				self->deriv                = &_ncm_spline_rbf_posdef_gauss_deriv;
				self->deriv2               = &_ncm_spline_rbf_posdef_gauss_deriv2;
				self->integ                = &_ncm_spline_rbf_integ_placeholder;
				
				self->name                 = g_strdup ("NcmSplineRBF:Posdef-Gaussian");
				self->shape_params         = ncm_vector_new (1);
				
				ncm_vector_set (self->shape_params, 0, GSL_NAN);
				break;
			case NCM_SPLINE_RBF_TYPE_GAUSS:
				self->default_shape_params = &_ncm_spline_rbf_gauss_default_shape_params;
				self->prepare_matrix       = &_ncm_spline_rbf_gauss_prepare_matrix;
				self->eval                 = &_ncm_spline_rbf_gauss_eval;
				self->deriv                = &_ncm_spline_rbf_gauss_deriv;
				self->deriv2               = &_ncm_spline_rbf_gauss_deriv2;
				
				self->name                 = g_strdup ("NcmSplineRBF:Gaussian");
				self->shape_params         = ncm_vector_new (1);
				
				ncm_vector_set (self->shape_params, 0, GSL_NAN);
				break;
			default:
				g_assert_not_reached ();
				break;
		}
	}
}

/**
 * ncm_spline_rbf_set_shape_params:
 * @rbf: a #NcmSplineRBF
 * @shape_params: a #NcmSplineRBFType
 * 
 * Sets the RBF shape parameters to @shape_params. 
 * 
 */ 
void 
ncm_spline_rbf_set_shape_params (NcmSplineRBF *rbf, NcmVector *shape_params)
{
	NcmSplineRBFPrivate * const self = rbf->priv;
	ncm_vector_memcpy (self->shape_params, shape_params);
}


/*
		SUNLinearSolver LS;
		N_Vector        xhat, x, b;
		realtype        *vecdata;
		gint i;
		gint maxl = len;
		gint ret;
		
		x          = N_VNew_Serial (len);
		xhat       = N_VNew_Serial (len);
		b          = N_VNew_Serial (len);

		vecdata = N_VGetArrayPointer (x);
		for (i = 0; i < len; i++) 
			vecdata[i] = ncm_vector_get (yv, i);

		vecdata = N_VGetArrayPointer (b);
		for (i = 0; i < len; i++) 
			vecdata[i] = ncm_vector_get (yv, i);

		self->prepare_matrix (self, xv, self->interp_matrix, len);

		//LS = SUNPCG (x, PREC_BOTH, maxl);
		LS = SUNSPGMR (x, PREC_LEFT, maxl);
		//LS = SUNSPTFQMR (x, PREC_RIGHT, maxl);
		//LS = SUNSPBCGS (x, PREC_RIGHT, maxl);

		ret = SUNLinSolSetPreconditioner (LS, self, PSetup, PSolve);
		printf ("ret 0 = %d\n", ret);
		
		ret = SUNLinSolSetATimes (LS, self, ATimes);
		printf ("ret 1 = %d\n", ret);

		ret = SUNLinSolInitialize (LS);
		printf ("ret 2 = %d\n", ret);

		ret = SUNLinSolSolve (LS, NULL, x, b, 1.0e-2);
		printf ("ret 3 = %d\n", ret);

		vecdata = N_VGetArrayPointer (x);
		for (i = 0; i < len; i++) 
			ncm_vector_set (self->coeff, i, vecdata[i]);
 */ 
