/***************************************************************************
 *            quadrature.h
 *
 *  Wed May 13 16:35:10 2009
 *  Copyright  2009  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef _NC_QUADRATURE_H
#define _NC_QUADRATURE_H

#include <glib.h>
#include <glib-object.h>

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

G_BEGIN_DECLS

typedef struct _NcmQuadFilonError NcmQuadFilonError;

/**
 * NcmQuadFilonError:
 * 
 * FIXME
 */
struct _NcmQuadFilonError
{
  /*< private >*/
  gsl_matrix *vandermonde;
  gsl_permutation *lu_p;
  gint lu_s;
  gsl_vector *inter_point;
  gsl_matrix *dPn_1;
  gsl_matrix *dPn_m1;
};

typedef struct _NcmQuadFilon NcmQuadFilon;

/**
 * NcmQuadFilon:
 * 
 * FIXME
 */
struct _NcmQuadFilon
{
  /*< private >*/
  gdouble omega;
  gdouble psi;
  gdouble range;
  gint order;
  gint inter_n;
  gint n;
  gsl_matrix *vandermonde;
  gsl_matrix *vandermonde_lu;
  gsl_permutation *lu_p;
  gint lu_s;
  gsl_vector *inter_point;
  gsl_matrix *dxnm_1;
  gsl_matrix *dxnm_m1;
  gsl_matrix *dPn_1;
  gsl_matrix *dPn_m1;
  gsl_vector *Re_mu;
  gsl_vector *Im_mu;
  gsl_vector *Re_b;
  gsl_vector *Im_b;
  gsl_vector *residual;
  NcmQuadFilonError err;
};

NcmQuadFilon *ncm_quadrature_filon_new (gdouble omega, gint inter_n, gint order, gdouble range);
gboolean ncm_quadrature_filon_calc_mu_dxnm (NcmQuadFilon *quadf);
gboolean ncm_quadrature_filon_calc_inter_point (NcmQuadFilon *quadf, gdouble g);
gboolean ncm_quadrature_filon_calc_vandermonde (NcmQuadFilon *quadf);
gboolean ncm_quadrature_filon_solve_vandermonde (NcmQuadFilon *quadf);
gboolean ncm_quadrature_filon_eval (NcmQuadFilon *quadf, gsl_function *F, gdouble xi, gsl_complex *res, gdouble *err);

G_END_DECLS

#endif /* _NC_QUADRATURE_H */
