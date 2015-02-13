/***************************************************************************
 *            magnus_iserles_ode.c
 *
 *  Thu Jun 18 14:01:28 2009
 *  Copyright  2009  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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
 * SECTION:magnus_iserles_ode
 * @title: NcmMIOde
 * @short_description: Magnus Iserles ode solver for fast oscillatory systems.
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/magnus_iserles_ode.h"

#include <math.h>

/**
 * ncm_magnus_iserles_ode_new: (skip)
 * @x0: FIXME
 * @g: a #NcmMIOdeFunction
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcmMIOde *
ncm_magnus_iserles_ode_new (long double x0, NcmMIOdeFunction *g)
{
  NcmMIOde *mi_ode = g_slice_new (NcmMIOde);

  mi_ode->h = 0.0;
  mi_ode->x0 = x0;
  mi_ode->xi = x0;
  mi_ode->g = g;
  mi_ode->g_x0 = g->func (x0, g->params);

  mi_ode->U = gsl_matrix_long_double_alloc (2, 2);
  mi_ode->A = gsl_matrix_long_double_alloc (2, 2);
  mi_ode->exp_A = gsl_matrix_long_double_alloc (2, 2);
  mi_ode->Omega = gsl_matrix_long_double_alloc (2, 2);
  mi_ode->exp_Omega = gsl_matrix_long_double_alloc (2, 2);

  gsl_matrix_long_double_set_identity (mi_ode->U);
  gsl_matrix_long_double_set_zero (mi_ode->A);
  
  return mi_ode;
}

/**
 * ncm_magnus_iserles_ode_step:
 * @mi_ode: a #NcmMIOde
 * @h: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gboolean
ncm_magnus_iserles_ode_step (NcmMIOde *mi_ode, long double h)
{
  long double Isin, Icos, I1;
  long double x0, xh_2, xh;
  long double g0, gh_2, gh;
  long double phi, phih, omega, omega2, psi, psi2, psi3;
  long double cos_psi, cos_phih;
  long double sin_psi, sin_phih;
  long double rb1, rb2, rb3, rbm;
  long double ib1, ib2, ib3, ibm;
  const long double a11 = 1.0 / 6.0;
  const long double a12 = 2.0 / 3.0;
  const long double a13 = 1.0 / 6.0;
  long double C_data[4];
/*
  gsl_matrix_view C_view = gsl_matrix_view_array (C_data, 2, 2);
  gsl_matrix *C = &C_view.matrix;
  gdouble D_data[4];
  gsl_matrix_view D_view = gsl_matrix_view_array (D_data, 2, 2);
  gsl_matrix *D = &D_view.matrix;
  gdouble expOmega_m1_data[4];
  gsl_matrix_view expOmega_m1_view = gsl_matrix_view_array (expOmega_m1_data, 2, 2);
  gsl_matrix *expOmega_m1 = &expOmega_m1_view.matrix;
  gint ret;
*/
  mi_ode->h = h;
  x0 = mi_ode->xi;
  xh_2 = x0 + h / 2.0;
  xh = x0 + h;
  
  g0 = mi_ode->g_x0;
  gh_2 = mi_ode->g->func (xh_2, mi_ode->g->params);
  gh = mi_ode->g->func (xh, mi_ode->g->params);

#define GS gh_2
  
  phi = sqrtl(fabsl(GS));
  phih = phi * h;

  cos_phih = cosl (phih);
  sin_phih = sinl (phih);

  gsl_matrix_long_double_set (mi_ode->exp_A, 0, 0, cos_phih);
  gsl_matrix_long_double_set (mi_ode->exp_A, 0, 1, sin_phih / phi);
  gsl_matrix_long_double_set (mi_ode->exp_A, 1, 0, -phi * sin_phih);
  gsl_matrix_long_double_set (mi_ode->exp_A, 1, 1, cos_phih);
  
  omega = 2.0 * phi;
  omega2 = omega * omega;
  psi = h * omega;
  psi2 = psi * psi;
  psi3 = psi2 * psi;
  cos_psi = cos_phih * cos_phih - sin_phih * sin_phih; //cos (psi);
  sin_psi = 2.0 * cos_phih * sin_phih; //sin (psi);

  rbm = -sin_psi / psi;
  ibm = (cos_psi - 1.0) / psi;
  
  rb1 = ((3.0 + cos_psi) * psi - 4.0 * sin_psi) / psi3;
  rb2 = (- 4.0 * psi * (1.0 + cos_psi) + 8.0 * sin_psi) / psi3;
  rb3 = (sin_psi * psi2 + (1.0 + 3.0 * cos_psi) * psi - 4.0 * sin_psi) / psi3;

  ib1 = (psi2 + sin_psi * psi - 4.0 * (1.0 - cos_psi)) / psi3;
  ib2 = (- 4.0 * psi * sin_psi + 8.0 * (1.0 - cos_psi)) / psi3;
  ib3 = (-cos_psi * psi2 + 3.0 * sin_psi * psi - 4.0 * (1.0 - cos_psi)) / psi3;

  Icos = h * (rb1 * g0 + rb2 * gh_2 + rb3 * gh + rbm * GS);
  Isin = h * (ib1 * g0 + ib2 * gh_2 + ib3 * gh + ibm * GS);
  I1   = h * (a11 * g0 + a12 * gh_2 + a13 * gh - GS);
/*
printf ("h: %.15Lg\n", h);
printf ("A: (% .15Lg, % .15Lg, % .15Lg)\n", x0, xh_2, xh);
printf ("G: (% .15Lg, % .15Lg, % .15Lg)\n", g0, gh_2, gh);
printf ("I: (% .15Lg, % .15Lg, % .15Lg)\n", Icos, Isin, I1);
*/
  {
    long double A = Isin / omega;
    long double A2 = A * A;
    long double B = (-Icos + I1) * 2.0 / omega2;
    long double C = -(Icos + I1) / 2.0;
    long double Delta2 = (A2 + B * C);
    long double Delta = sqrt(fabs(Delta2));
    long double sinn_Delta_Delta;
    long double cosn_Delta;
    //gdouble sinn_Delta_2_2;
//printf ("DAMN A B C: % .15Lg % .15Lg % .15Lg\n", A, B, C);

    if (Delta2 < 0)
    {
      gdouble sinn_Delta;
      sinn_Delta = sinl (Delta);
      sinn_Delta_Delta = sinn_Delta / Delta;
      cosn_Delta = cosl (Delta);
    }
    else
    {
      gdouble sinn_Delta;
      
      sinn_Delta = sinhl (Delta);
      sinn_Delta_Delta = sinn_Delta / Delta;
      cosn_Delta = coshl (Delta);
    }

    gsl_matrix_long_double_set (mi_ode->exp_Omega, 0, 0, cosn_Delta + A * sinn_Delta_Delta);
    gsl_matrix_long_double_set (mi_ode->exp_Omega, 0, 1, B * sinn_Delta_Delta);
    gsl_matrix_long_double_set (mi_ode->exp_Omega, 1, 0, C * sinn_Delta_Delta);
    gsl_matrix_long_double_set (mi_ode->exp_Omega, 1, 1, cosn_Delta - A * sinn_Delta_Delta);
  }
/*  
  gsl_matrix_memcpy (C, mi_ode->exp_A);  
  ret = gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, mi_ode->exp_A, expOmega_m1, 1.0, C);
  NCM_TEST_GSL_RESULT("ncm_magnus_iserles_ode_step", ret);
  {
    gdouble temp_u;
    temp_u     = gsl_matrix_get (C, 0, 0) * mi_ode->u + gsl_matrix_get (C, 0, 1) * mi_ode->up;
    mi_ode->up = gsl_matrix_get (C, 1, 0) * mi_ode->u + gsl_matrix_get (C, 1, 1) * mi_ode->up;
    mi_ode->u = temp_u;
  }
*/
  
#define EA(i,j) (gsl_matrix_long_double_get(mi_ode->exp_A,i-1,j-1))
#define EO(i,j) (gsl_matrix_long_double_get(mi_ode->exp_Omega,i-1,j-1))
#define U(i,j) (gsl_matrix_long_double_get(mi_ode->U,i-1,j-1))
  
  C_data[0] = (EA(1,1) * EO (1,1) + EA(1,2) * EO (2,1)) * U (1,1) + 
              (EA(1,1) * EO (1,2) + EA(1,2) * EO (2,2)) * U (2,1);
  C_data[1] = (EA(1,1) * EO (1,1) + EA(1,2) * EO (2,1)) * U (1,2) + 
              (EA(1,1) * EO (1,2) + EA(1,2) * EO (2,2)) * U (2,2);
  C_data[2] = (EA(2,1) * EO (1,1) + EA(2,2) * EO (2,1)) * U (1,1) + 
              (EA(2,1) * EO (1,2) + EA(2,2) * EO (2,2)) * U (2,1);
  C_data[3] = (EA(2,1) * EO (1,1) + EA(2,2) * EO (2,1)) * U (1,2) + 
              (EA(2,1) * EO (1,2) + EA(2,2) * EO (2,2)) * U (2,2);
  
  mi_ode->U->data[0] = C_data[0];
  mi_ode->U->data[1] = C_data[1];
  mi_ode->U->data[2] = C_data[2];
  mi_ode->U->data[3] = C_data[3];
  
  mi_ode->psi = psi;
  mi_ode->xi = xh;
  mi_ode->g_x0 = gh;
    
  return TRUE;
}

/**
 * ncm_magnus_iserles_ode_step_frac:
 * @mi_ode: a #NcmMIOde
 * @frac: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gboolean
ncm_magnus_iserles_ode_step_frac (NcmMIOde *mi_ode, long double frac)
{
  return ncm_magnus_iserles_ode_step (mi_ode, mi_ode->xi * frac);
}

/**
 * ncm_magnus_iserles_ode_eval_vec:
 * @mi_ode: a #NcmMIOde
 * @u_i: FIXME
 * @up_i: FIXME
 * @u: FIXME
 * @up: FIXME 
 *
 * FIXME
*/
void
ncm_magnus_iserles_ode_eval_vec (NcmMIOde *mi_ode, long double u_i, long double up_i, long double *u, long double *up)
{
  u[0]  = gsl_matrix_long_double_get (mi_ode->U, 0, 0) * u_i + 
    gsl_matrix_long_double_get (mi_ode->U, 0, 1) * up_i;
  up[0] = gsl_matrix_long_double_get (mi_ode->U, 1, 0) * u_i + 
    gsl_matrix_long_double_get (mi_ode->U, 1, 1) * up_i;
}
