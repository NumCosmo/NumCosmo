/***************************************************************************
 *            quadrature.c
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

/**
 * SECTION:quadrature
 * @title: Quadrature Algorithims
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/quadrature.h"
#include "math/ncm_cfg.h"

#include <string.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#define LOCAL_MAX_ORDER 100

/**
 * ncm_quadrature_filon_new: (skip)
 * @omega: FIXME
 * @inter_n: FIXME
 * @order: FIXME
 * @range: FIXME
 * 
 * FIXME
 *
 * Returns: FIXME
*/ 
NcmQuadFilon *
ncm_quadrature_filon_new (gdouble omega, gint inter_n, gint order, gdouble range)
{
  NcmQuadFilon *quadf;
  gint n = inter_n + 2 * (++order);

  g_assert (n + 2 <= LOCAL_MAX_ORDER);
  g_assert (order > 1);
  g_assert (inter_n >= 0);

  quadf = (NcmQuadFilon *)g_slice_new (NcmQuadFilon);

  quadf->omega = omega;
  quadf->psi = omega * range / 2.0;
  quadf->range = range;
  quadf->order = order;
  quadf->inter_n = inter_n;
  quadf->n = n;
  /*********************************************************** 
   * Linear system matrix and interpolation points for order *
   *                                                         *
   ***********************************************************/
  quadf->lu_p = gsl_permutation_alloc (n);
  quadf->vandermonde = gsl_matrix_alloc (n, n);
  quadf->vandermonde_lu = gsl_matrix_alloc (n, n);
  quadf->inter_point = gsl_vector_alloc (n);

  quadf->err.lu_p = gsl_permutation_alloc (n - 2);
  quadf->err.vandermonde = gsl_matrix_alloc (n - 2, n - 2);
  quadf->err.inter_point = gsl_vector_alloc (n - 2);
  /***********************************************************/

  quadf->dxnm_1 = gsl_matrix_alloc (order, n);
  quadf->dxnm_m1 = gsl_matrix_alloc (order, n);
  quadf->dPn_1 = gsl_matrix_alloc (order, n);
  quadf->dPn_m1 = gsl_matrix_alloc (order, n);

  quadf->err.dPn_1 = gsl_matrix_alloc (order - 1, n - 2);
  quadf->err.dPn_m1 = gsl_matrix_alloc (order - 1, n - 2);
  
  quadf->Re_mu = gsl_vector_alloc (n);
  quadf->Im_mu = gsl_vector_alloc (n);
  quadf->Re_b = gsl_vector_alloc (n);
  quadf->Im_b = gsl_vector_alloc (n);

  quadf->residual = gsl_vector_alloc (n);

  ncm_quadrature_filon_calc_mu_dxnm (quadf);
  ncm_quadrature_filon_calc_inter_point (quadf, 1.0);
  ncm_quadrature_filon_calc_vandermonde (quadf);
  ncm_quadrature_filon_solve_vandermonde (quadf);
  
  return quadf;
}

/*********************************************************** 
 * Calculate the moments $\int_{-1}^{1}x^n\exp{iwx}dx$ and *
 * derivatives $\frac{d x^n}{dx^l}|_{1 \text{and} -1}$     *
***********************************************************/

/**
 * ncm_quadrature_filon_calc_mu_dxnm: 
 * @quadf: a #NcmQuadFilon
 *
 * FIXME
 *
 * Returns: FIXME 
*/ 
gboolean
ncm_quadrature_filon_calc_mu_dxnm (NcmQuadFilon *quadf)
{
  gdouble jl[LOCAL_MAX_ORDER];
  gint k;

  gsl_sf_bessel_jl_array (quadf->n, quadf->psi, jl);
  
  for (k = 0; k < quadf->n; k++)
  {
    if (GSL_IS_EVEN(k))
    {
      gint sign_jl = (k/2) % 2 == 1 ? -1 : 1;
      gsl_vector_set (quadf->Re_mu, k, sign_jl * 2.0 * jl[k]);
      gsl_vector_set (quadf->Im_mu, k, 0.0);
    }
    else
    {
      gint sign_jl = ((k-1)/2) % 2 == 1 ? -1 : 1;
      gsl_vector_set (quadf->Re_mu, k, 0.0);
      gsl_vector_set (quadf->Im_mu, k, sign_jl * 2.0 * jl[k]);
    }
  }

  return TRUE;
}

/***********************************************************
 * Initializing the interpolation points                   *
 ***********************************************************/

/**
 * ncm_quadrature_filon_calc_inter_point: 
 * @quadf: a #NcmQuadFilon
 * @g: FIXME 
 *
 * FIXME
 *
 * Returns: FIXME 
*/
gboolean
ncm_quadrature_filon_calc_inter_point (NcmQuadFilon *quadf, gdouble g)
{
  gint k;
  for (k = 0; k < quadf->order; k++)
  {
    gsl_vector_set (quadf->inter_point,            k, -1.0 + k * g / quadf->psi);
    gsl_vector_set (quadf->inter_point, quadf->n-k-1,  1.0 - k * g / quadf->psi);
    if (k < quadf->order - 1)
    {
      gsl_vector_set (quadf->err.inter_point,            k, -1.0 + k * g / quadf->psi);
      gsl_vector_set (quadf->err.inter_point, quadf->n-k-3,  1.0 - k * g / quadf->psi);
    }
  }
  for (k = 0; k < quadf->inter_n; k++)
  {
    gdouble chebchev_node_k = -cos ((2.0 * k + 1.0) * M_PI / (2.0 * quadf->inter_n));
    gsl_vector_set (quadf->inter_point, quadf->order + k, chebchev_node_k);
    gsl_vector_set (quadf->err.inter_point, quadf->order - 1 + k, chebchev_node_k);
  }
  
  return TRUE;
}

/***********************************************************
 * Initializing the Vandermonde matrix $c_l^n$ and         *
 * and aplying the LU decomposition to solve $c v = u$.    *
 ***********************************************************/

/**
 * ncm_quadrature_filon_calc_vandermonde: 
 * @quadf: a #NcmQuadFilon
 *
 * FIXME
 *
 * Returns: FIXME 
*/
gboolean
ncm_quadrature_filon_calc_vandermonde (NcmQuadFilon *quadf)
{
  gint k;

  for (k = 0; k < quadf->n; k++)
  {
    gsl_vector_view vrowv = gsl_matrix_row (quadf->vandermonde, k);
    gsl_sf_legendre_Pl_array (quadf->n-1, gsl_vector_get(quadf->inter_point, k), vrowv.vector.data);
  }

  gsl_matrix_transpose (quadf->vandermonde);
  
  gsl_matrix_memcpy (quadf->vandermonde_lu, quadf->vandermonde);
  gsl_linalg_LU_decomp (quadf->vandermonde_lu, quadf->lu_p, &quadf->lu_s);
  gsl_linalg_LU_decomp (quadf->err.vandermonde, quadf->err.lu_p, &quadf->err.lu_s);
  
  return TRUE;
}

/**
 * solve_vandermonde: 
 * @x: FIXME
 * @w: FIXME
 * @q: FIXME
 * @n: FIXME
 *
 * FIXME
 *
*/ 
void
solve_vandermonde (double x[], double w[], double q[], int n)
{
  int i, j, k;
  double b, s, t, xx;
  double c[LOCAL_MAX_ORDER];

  if (n == 1) 
    w[1] = q[1];
  else 
  {
    for (i=1;i<=n;i++) 
      c[i] = 0.0;

    c[n] = -x[1]; 
    for (i = 2; i <= n; i++)
    {
      xx = -x[i];
      for (j = (n+1-i); j <= (n-1); j++)
        c[j] += xx*c[j+1];
      c[n] += xx;
    }
    for (i = 1; i <=n ; i++)
    {
      xx = x[i];
      t = b = 1.0;
      s = q[n];
      for (k=n;k>=2;k--)
      {
        b = c[k] + xx * b;
        s += q[k-1] * b;
        t = xx * t + b;
      }
      w[i] = s / t;
    }
  }
}

/**
 * ncm_quadrature_filon_solve_vandermonde: 
 * @quadf: a #NcmQuadFilon
 *
 * FIXME
 *
 * Returns: FIXME 
*/
gboolean
ncm_quadrature_filon_solve_vandermonde (NcmQuadFilon *quadf)
{
  gint k, ret;

  ret = gsl_linalg_LU_solve (quadf->vandermonde_lu, quadf->lu_p, quadf->Re_mu, quadf->Re_b);
  NCM_TEST_GSL_RESULT("ncm_quadrature_filon_solve_vandermonde", ret);
  
  ret = gsl_linalg_LU_solve (quadf->vandermonde_lu, quadf->lu_p, quadf->Im_mu, quadf->Im_b);
  NCM_TEST_GSL_RESULT("ncm_quadrature_filon_solve_vandermonde", ret);
  
  for (k = 0; k < quadf->order; k++)
  {
    gsl_vector_view dxnm_1v  = gsl_matrix_row (quadf->dxnm_1,  k);
    gsl_vector_view dxnm_m1v  = gsl_matrix_row (quadf->dxnm_m1,  k);
    gsl_vector_view dPn_1v  = gsl_matrix_row (quadf->dPn_1,  k);
    gsl_vector_view dPn_m1v  = gsl_matrix_row (quadf->dPn_m1,  k);

    ret = gsl_linalg_LU_solve (quadf->vandermonde_lu, quadf->lu_p, &dxnm_1v.vector, &dPn_1v.vector);
    NCM_TEST_GSL_RESULT("ncm_quadrature_filon_solve_vandermonde", ret);
    ret = gsl_linalg_LU_solve (quadf->vandermonde_lu, quadf->lu_p, &dxnm_m1v.vector, &dPn_m1v.vector);
    NCM_TEST_GSL_RESULT("ncm_quadrature_filon_solve_vandermonde", ret);
    if (k < quadf->order - 1)
    {
      gsl_vector_view err_dxnm_1v = gsl_vector_subvector (&dxnm_1v.vector, 0, quadf->n - 2);
      gsl_vector_view err_dxnm_m1v = gsl_vector_subvector (&dxnm_m1v.vector, 0, quadf->n - 2);
      gsl_vector_view err_dPn_1v  = gsl_matrix_row (quadf->err.dPn_1, k);
      gsl_vector_view err_dPn_m1v  = gsl_matrix_row (quadf->err.dPn_m1, k);

      ret = gsl_linalg_LU_solve (quadf->err.vandermonde, quadf->err.lu_p, &err_dxnm_1v.vector, &err_dPn_1v.vector);
      NCM_TEST_GSL_RESULT("ncm_quadrature_filon_solve_vandermonde", ret);
      ret = gsl_linalg_LU_solve (quadf->err.vandermonde, quadf->err.lu_p, &err_dxnm_m1v.vector, &err_dPn_m1v.vector);
      NCM_TEST_GSL_RESULT("ncm_quadrature_filon_solve_vandermonde", ret);
    }
  }
  
  return TRUE;
}

/**
 * ncm_quadrature_filon_eval: (skip)
 * @quadf: a #NcmQuadFilon
 * @F: FIXME
 * @xi: FIXME
 * @res: FIXME
 * @err: FIXME 
 *
 * FIXME
 *
 * Returns: FIXME 
*/
gboolean
ncm_quadrature_filon_eval (NcmQuadFilon *quadf, gsl_function *F, gdouble xi, gsl_complex *res, gdouble *err)
{
  const gdouble u = xi + quadf->range / 2.0;
  const gdouble v = quadf->range / 2.0;
  gdouble dfa_1[LOCAL_MAX_ORDER];
  gdouble dfa_m1[LOCAL_MAX_ORDER];
  gdouble dfb_1[LOCAL_MAX_ORDER];
  gdouble dfb_m1[LOCAL_MAX_ORDER];
  gdouble psi_n;
  gsl_complex err_1, err_m1;
  gint k, kk;
  gint ip[4] = {1, -1, -1, 1};

  GSL_SET_COMPLEX(res, 0.0, 0.0);
  GSL_SET_COMPLEX(&err_1, 0.0, 0.0);
  GSL_SET_COMPLEX(&err_m1, 0.0, 0.0);

  memset (dfa_1,  0, sizeof(gdouble) * LOCAL_MAX_ORDER);
  memset (dfa_m1, 0, sizeof(gdouble) * LOCAL_MAX_ORDER);
  memset (dfb_1,  0, sizeof(gdouble) * LOCAL_MAX_ORDER);
  memset (dfb_m1, 0, sizeof(gdouble) * LOCAL_MAX_ORDER);

  kk = 0;
  for (k = 0; k < quadf->n; k++)
  {
    gdouble x = u + v * gsl_vector_get(quadf->inter_point, k);
    gdouble f_n = F->function (x, F->params);
    gint l;

//    printf ("Inter at (%g + %g * %g)[%g]=[%g]\n", u, v, gsl_vector_get(quadf->inter_point, k), x, f_n);
    GSL_REAL(*res) += f_n * gsl_vector_get (quadf->Re_b, k);
    GSL_IMAG(*res) += f_n * gsl_vector_get (quadf->Im_b, k);
    
    for (l = 0; l < quadf->order-1; l++)
    {
      gsl_vector_view dPn_1v = gsl_matrix_row (quadf->dPn_1,  l);
      gsl_vector_view dPn_m1v = gsl_matrix_row (quadf->dPn_m1,  l);
      gsl_vector_view err_dPn_1v = gsl_matrix_row (quadf->err.dPn_1,  l);
      gsl_vector_view err_dPn_m1v = gsl_matrix_row (quadf->err.dPn_m1,  l);

      dfa_1[l]  += f_n * gsl_vector_get (&dPn_1v.vector, k);
      dfa_m1[l] += f_n * gsl_vector_get (&dPn_m1v.vector, k);
      
      printf ("[%d, %d] -> %.15g | %.15g\n", k, l, gsl_vector_get (&dPn_1v.vector, k), gsl_vector_get (&dPn_m1v.vector, k));
      
      if ((k != quadf->order - 1) && (k != quadf->n - quadf->order))
      {
        dfb_1[l]  += f_n * gsl_vector_get (&err_dPn_1v.vector, kk);
        dfb_m1[l] += f_n * gsl_vector_get (&err_dPn_m1v.vector, kk);
      }
    }
    if ((k != quadf->order - 1) && (k != quadf->n - quadf->order))
      kk++;
  }

  psi_n = 1.0;
  kk = 0;
  for (k = quadf->order - 2; k >= 0; k--)
  {
    if (GSL_IS_EVEN(kk))
    {
      GSL_REAL(err_1)  += (dfa_1[k]  - dfb_1[k])  * ip[k % 4] * psi_n;
      GSL_REAL(err_m1) += (dfa_m1[k] - dfb_m1[k]) * ip[k % 4] * psi_n;
    }
    else
    {
      GSL_IMAG(err_1)  += (dfa_1[k]  - dfb_1[k])  * ip[k % 4] * psi_n;
      GSL_IMAG(err_m1) += (dfa_m1[k] - dfb_m1[k]) * ip[k % 4] * psi_n;
    }
    kk++;
    psi_n *= quadf->psi;
  }

  psi_n *= quadf->psi;

  *res = gsl_complex_mul(gsl_complex_polar (v, quadf->psi / quadf->range * u), *res);
  *err = (gsl_complex_abs (err_1) + gsl_complex_abs (err_m1)) / psi_n;
  
  return TRUE;
}
