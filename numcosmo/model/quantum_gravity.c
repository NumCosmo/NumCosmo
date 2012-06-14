/***************************************************************************
 *            quantum_gravity.c
 *
 *  Wed Jan 28 11:21:00 2009
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
 * SECTION:quantum_gravity
 * @title: Quantum Gravity Bouncing Model
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_exp.h>

typedef struct _NcHICosmoQG 
{
  gpointer cvode;
  gboolean initialized;
  N_Vector y;
  N_Vector yQ;
  gdouble reltol;
  gdouble abstol;
  gdouble lambda_i;
  gdouble lambda_d;
  gdouble lambda_b;
  gdouble eta_b;
  gdouble last_lambda;
  gdouble last_eta;
  GArray *lambda;
  GArray *eta;
  GArray *int_1_zeta2;
  GArray *cs2zeta2_int_1_zeta2;
  GArray *x;
  GArray *g;
  GArray *gbar;
  GArray *gbarbar;
  GArray *cs2;
  GArray *V;
  GArray *k_cross;
  gsl_interp *eta_lambda;
  gsl_interp *x_lambda;
  gsl_interp *gbar_lambda;
  gsl_interp *gbarbar_lambda;
  gsl_interp *int_1_zeta2_lambda;
  gsl_interp *cs2zeta2_int_1_zeta2_lambda;
  gsl_interp *lambda_k_cross;
  gsl_interp *cs2_lambda;
  gsl_interp *V_lambda;
  gsl_interp_accel *lambda_accel;
  gsl_interp_accel *eta_accel;
  gsl_interp_accel *z_accel;
  gsl_interp_accel *k_cross_accel;
  gboolean initialized_spline;
} NcHICosmoQG;

typedef struct _NcHICosmoQGScaleFactor
{
  NcmModel *model;
  gboolean first_order;
  gdouble sign;
  gdouble x_root_scale;
} NcHICosmoQGScaleFactor;

gboolean nc_hicosmo_qg_init_spline (NcHICosmoQG *qgint, NcmModel *model);

#define VECTOR   (model->params)
#define MACRO_H0 (ncm_vector_get (VECTOR, 0))
#define OMEGA_M  (ncm_vector_get (VECTOR, 1))
#define G_B      (ncm_vector_get (VECTOR, 2))
#define OMEGA_R  (ncm_vector_get (VECTOR, 3))
#define W        (ncm_vector_get (VECTOR, 4))
#define OMEGA_X  ((1.0 - OMEGA_R - OMEGA_M)*0.0)

NcHICosmo *
nc_hicosmo_qg_new (void)
{
	g_assert_not_reached ();
}


gdouble
nc_hicosmo_qg_cs2 (NcmModel *model, gdouble x, gpointer userdata)
{
  gdouble R = 3.0 * (1.0 + W) * OMEGA_M / (4.0 * OMEGA_R * x * pow(x, -3.0*W));
  return (1.0 + 3.0 * W * R) / (3.0 * (1.0 + R));
}

gdouble
nc_hicosmo_qg_dcs2 (NcmModel *model, gdouble x, gpointer userdata)
{
  gdouble R = 3.0 * (1.0 + W) * OMEGA_M / (4.0 * OMEGA_R * x * pow(x, -3.0*W));
  return (1.0 - 3.0 * W) / (3.0 * x) * ((1.0 - 3.0 * W)*R+3.0*W*R*R) / gsl_pow_2 (1.0 + R);
}

gdouble
nc_hicosmo_qg_beta (NcmModel *model, gdouble x, gpointer userdata)
{
  gdouble OM = 3.0 * (1.0 + W) * OMEGA_M * pow (x, 1.0 + 3.0 * W);
  gdouble OR = 4.0 * OMEGA_R * x * x;

  return (OM + OR) / 2.0;
}

gdouble
nc_hicosmo_qg_Omega_x2 (NcmModel *model, gdouble x, gpointer userdata)
{
  gdouble x3w = pow (x, 3.0 * W);
  return x * (x3w * OMEGA_M + x * OMEGA_R);
}

gdouble
nc_hicosmo_qg_dOmega_x2 (NcmModel *model, gdouble x, gpointer userdata)
{
  gdouble x3w = pow (x, 3.0 * W);
  return ((1.0 + 3.0 * W) * x3w * OMEGA_M + 2.0 * x * OMEGA_R);
}

gdouble
nc_hicosmo_qg_xbar (NcmModel *model, gdouble x, gpointer userdata)
{
  gdouble x3w = pow (x, 3.0 * W);
  return sqrt (x * (x3w * OMEGA_M + x * OMEGA_R));
}

gdouble
nc_hicosmo_qg_d2sqrtxxbarzeta_sqrtxxbarzeta (NcmModel *model, gdouble x, gpointer userdata)
{
  gdouble x3w = pow (x, 3.0 * W);
  gdouble x6w = x3w * x3w;
  gdouble x9w = x6w * x3w;
  gdouble x2 = x * x;
  gdouble x3 = x2 * x;
  gdouble w2 = W * W;
  gdouble one_pw = 1.0 + W;
  gdouble one_pw2 = gsl_pow_2 (one_pw);
  gdouble Om2 = OMEGA_M * OMEGA_M;
  gdouble Om3 = Om2 * OMEGA_M;
  gdouble Or = OMEGA_R;
  gdouble Or2 = Or * Or;
  gdouble Or3 = Or2 * Or;

  return ( x3w / x2 * OMEGA_M * 
          (81.0 * x9w * w2 * one_pw2 * (5.0 + 9.0 * (-2.0 + W) * W ) * Om3 + 
           36.0 * x * x6w * W * one_pw * (-2.0 + 3.0 * W) * (-1.0 + 3.0 * W) * (7.0 + 15.0 * W) * Om2 * Or + 
           12.0 * x2 * x3w * (1.0 + W) * (-1.0 + 3.0 * W) *
           (-28.0 + 3.0 * W * (28.0 + 9.0 * W * (-7.0 + W * (8.0 + 3.0 * W)))) * OMEGA_M * Or2
           - 32.0 * x3 * (-2.0 + 3.0 * W) * (-1.0 + 3.0 * W) * (-4.0 + 3.0 * W + 9.0 * w2) * Or3
           )
          ) / 
    (16.0 * gsl_pow_2 (x3w * OMEGA_M + x * Or) * 
     gsl_pow_2 (9.0 * x3w * W * (1.0 + W) * OMEGA_M + 4.0 * x * Or));
}

gdouble
nc_hicosmo_qg_cs2_xxbar2 (NcmModel *model, gdouble x, gpointer userdata)
{
  gdouble x3 = gsl_pow_3 (x);
  gdouble x4 = x3 * x;
  gdouble x3w = pow (x, 3.0 * W);
  gdouble R = 3.0 * (1.0 + W) * OMEGA_M / (4.0 * OMEGA_R * x / x3w);
  gdouble cs2 = (1.0 + 3.0 * W * R) / (3.0 * (1.0 + R));
  return cs2 / (x3 * x3w * OMEGA_M + x4 * OMEGA_R);
}

void
nc_hicosmo_qg_evolfunc (NcmModel *model, long double x, long double *x2d2sqrtxxbarzeta_sqrtxxbarzeta, long double *x2cs2_xxbar2)
{
  const long double w = W;
  const long double x3w = powl (x, 3.0L * w);
  const long double x6w = x3w * x3w;
  const long double x9w = x6w * x3w;
  const long double x2 = x * x;
  const long double x3 = x2 * x;
  const long double w2 = w * w;
  const long double Om = OMEGA_M;
  const long double one_pw = 1.0L + w;
  const long double one_pw2 = one_pw * one_pw;
  const long double Om2 = Om * Om;
  const long double Om3 = Om2 * Om;
  const long double Or = OMEGA_R * 1.0L;
  const long double Or2 = Or * Or;
  const long double Or3 = Or2 * Or;
  const long double one_R = (4.0L * Or * x / x3w) / 3.0L * one_pw * Om;
  const long double cs2 = (one_R + 3.0L * W) / (3.0L * (1.0L + one_R));
  const long double AAA = (x3w * Om + x * Or);
  const long double AAA2 = AAA * AAA;
  const long double BBB = (9.0L * x3w * W * one_pw * Om + 4.0L * x * Or);
  const long double BBB2 = BBB * BBB;
  
  *x2d2sqrtxxbarzeta_sqrtxxbarzeta = ( x3w * Om *
          (81.0L * x9w * w2 * one_pw2 * (5.0L + 9.0L * (-2.0L + W) * W ) * Om3 +
           36.0L * x * x6w * W * one_pw * (-2.0L + 3.0L * W) * (-1.0L + 3.0L * W) * (7.0L + 15.0L * W) * Om2 * Or +
           12.0L * x2 * x3w * one_pw * (-1.0L + 3.0L * W) *
           (-28.0L + 3.0L * W * (28.0L + 9.0L * W * (-7.0L + W * (8.0L + 3.0L * W)))) * Om * Or2
           - 32.0L * x3 * (-2.0L + 3.0L * W) * (-1.0L + 3.0L * W) * (-4.0L + 3.0L * W + 9.0L * w2) * Or3
           )
          ) /
    (16.0L * AAA2 * BBB2 );

  *x2cs2_xxbar2 = cs2 / (x * x3w * Om + x2 * Or);

}

gdouble
nc_hicosmo_qg_xxbarzeta2 (NcmModel *model, gdouble x, gpointer userdata)
{
  gdouble x3w = pow (x, 3.0 * W);
  
  return (3.0 * gsl_pow_2 (3.0 * x3w * (1.0 + W) * OMEGA_M + 4.0 * x * OMEGA_R)) /
   (2.0 * sqrt (x * (x3w * OMEGA_M + x * OMEGA_R)) *
     (9.0 * x3w * W * (1.0 + W) * OMEGA_M + 4.0 * x * OMEGA_R));
}

gdouble
nc_hicosmo_qg_dxxbarzeta2_xxbarzeta2 (NcmModel *model, gdouble x, gpointer userdata)
{
  gdouble x3w = pow (x, 3.0 * W);

  return ((-1.0 + 3.0 * W) * 
          ( 1.0 / x + OMEGA_R * ( 1.0 / (x3w * OMEGA_M + x * OMEGA_R) - 
                                 16.0 / ( 3.0 * x3w * (1.0 + W) * OMEGA_M + 4.0 * x * OMEGA_R) + 
                                 8.0 / ( 9.0 * x3w * W * (1.0 + W) * OMEGA_M + 4.0 * x * OMEGA_R)
                                 )
           )
          ) / 2.0;
}

gdouble
nc_hicosmo_qg_zeta (NcmModel *model, gdouble x, gpointer userdata)
{
  gdouble R = 3.0 * (1.0 + W) * OMEGA_M / (4.0 * OMEGA_R * x * pow(x, -3.0 * W));
  gdouble F = OMEGA_X / (OMEGA_R * gsl_pow_4 (x));

  return sqrt(6.0) / x * (1.0 + R) / sqrt ((1.0 + 4.0 * R / (3.0 * (1.0 + W)) + F) * (1.0 + 3.0 * W * R));
}

gdouble
nc_hicosmo_qg_dzeta_zeta (NcmModel *model, gdouble x, gpointer userdata)
{
  const gdouble R = 3.0 * (1.0 + W) * OMEGA_M / (4.0 * OMEGA_R * x * pow(x, -3.0 * W));
  const gdouble F = OMEGA_X / (OMEGA_R * gsl_pow_4 (x));
  const gdouble G1 = (6.0 * (1.0 + W) * F + 2.0 * (1.0 - 3.0 * W) * R) / (3.0 * (1.0 + W) * (1.0 + F) + 4.0 * R); 
  const gdouble G2 = 3.0 * W * (1.0 - 3.0 * W) * R / (2.0 * (1.0 + 3.0 * W * R));
  const gdouble cs2 = (1.0 + 3.0 * W * R) / (3.0 * (1.0 + R));

  return (3.0 * cs2 + G1 + G2 - 2.0) / x;
}

gdouble
nc_hicosmo_qg_dcs2zeta2_cs2zeta2 (NcmModel *model, gdouble x, gpointer userdata)
{
  const gdouble R = 3.0 * (1.0 + W) * OMEGA_M / (4.0 * OMEGA_R * x * pow(x, -3.0 * W));
  return -((6.0 * (1.0 + W) + R * (13.0 + 8.0 * R + 3.0 * (4.0 - 3.0 * W) * W) ) / ((1.0 + R) * x  * (3.0 + 4.0 * R + 3.0 * W)));
}

gdouble
nc_hicosmo_qg_dx2dcs2zeta2_cs2zeta2 (NcmModel *model, gdouble x, gpointer userdata)
{
  const gdouble R = 3.0 * (1.0 + W) * OMEGA_M / (4.0 * OMEGA_R * x * pow(x, -3.0 * W));
  const gdouble R2 = R  * R;
  const gdouble R3 = R2 * R;
  const gdouble R4 = R2 * R2;
  const gdouble onepw2 = (1.0 + W) * (1.0 + W);
  return (-32.0 * R4 - 18.0 * onepw2 + 3.0 * R * (1.0 + W) * (-28.0 + 9.0 * (-1.0 + W) * W * (1.0 + 3.0 * W)) + R2 * (-139.0 + 9.0 * W * (-19.0 + 3.0 * W * (1.0 + W))) - 4.0 * R3 * (26.0 + 9.0 * W * (3.0 + W * (-4.0 + 3.0 * W)))) / (gsl_pow_2 (1.0 + R) * gsl_pow_2 (3.0 + 4.0 * R + 3.0 * W));
}

gdouble
nc_hicosmo_qg_ddzeta_zeta (NcmModel *model, gdouble x, gpointer userdata)
{
  gdouble R = 3.0 * (1.0 + W) * OMEGA_M / (4.0 * OMEGA_R * x * pow(x, -3.0*W));
  gdouble F =  OMEGA_X / (OMEGA_R * gsl_pow_4 (x));
  gdouble G1 = (6.0 * (1.0 + W) * F + 2.0 * (1.0 - 3.0 * W) * R) / (3.0 * (1.0 + W) * (1.0 + F) + 4.0 * R); 
  gdouble G2 = 3.0 * W * (1.0 - 3.0 * W) * R / (2.0 * (1.0 + 3.0 * W * R));
  gdouble cs2 = (1.0 + 3.0 * W * R) / (3.0 * (1.0 + R));

  gdouble dcs2 = (1.0 - 3.0 * W) / (3.0 * x) * ((1.0 - 3.0 * W) * R + 3.0 * W * R * R) / gsl_pow_2 (1.0 + R);
  gdouble dG1 = (2.0 * G1 * G1 - (24.0 * (1.0 + W) * F + 2.0 * gsl_pow_2 (1.0 - 3.0 * W) * R) / (3.0 * (1.0 + W) * (1.0 + F) + 4.0 * R)) / x;
  gdouble dG2 = (2.0 * G2 * G2 - (1.0 - 3.0 * W) * G2) / x;

  gdouble dzeta_zeta = (3.0 * cs2 + G1 + G2 - 2.0) / x;
  
  return ((3.0 * dcs2 + dG1 + dG2)  - dzeta_zeta) / x + dzeta_zeta * dzeta_zeta;
}

gdouble
nc_hicosmo_qg_xddzeta_zeta_mxdzeta_zeta2_dzeta_zeta (NcmModel *model, gdouble x, gpointer userdata)
{
  gdouble R = 3.0 * (1.0 + W) * OMEGA_M / (4.0 * OMEGA_R * x * pow(x, -3.0*W));
  gdouble F =  OMEGA_X / (OMEGA_R * gsl_pow_4 (x));
  gdouble G1 = (6.0 * (1.0 + W)*F+2.0*(1.0-3.0*W)*R) / (3.0 * (1.0 + W)*(1.0 + F) + 4.0 * R); 
  gdouble G2 = 3.0 * W * (1.0 - 3.0 * W) * R / (2.0 * (1.0 + 3.0 * W * R));

  gdouble dcs2 = (1.0 - 3.0 * W) / (3.0 * x) * ((1.0 - 3.0 * W)*R+3.0*W*R*R) / gsl_pow_2 (1.0 + R);
  gdouble dG1 = (2.0 * G1*G1 - (24.0 * (1.0 + W) * F + 2.0 * gsl_pow_2 (1.0 - 3.0 * W) * R) / (3.0 * (1.0 + W)*(1.0 + F) + 4.0 * R)) / x;
  gdouble dG2 = (2.0 * G2*G2 - (1.0 - 3.0 * W) * G2) / x;

  return 3.0 * dcs2 + dG1 + dG2;
}

gdouble
nc_hicosmo_qg_gbar2 (NcmModel *model, gdouble x)
{
  gdouble x1_m3w = x * pow (x, -3.0 * W);
  gdouble x2 = x * x;
  gdouble x4 = x2 * x2;

  gdouble xb = exp(-G_B);
  gdouble xb2 = xb * xb;
  gdouble xb3 = xb2 * xb;
  gdouble xb3_m3w =  xb3 * pow (xb, -3.0 * W);
  
  gdouble xb6 = xb3 * xb3;

//  g_assert (x < xb*1e-16);

  return OMEGA_R * (1.0 - x2 / xb2) + OMEGA_M * (1.0 / x1_m3w - x2 / xb3_m3w) + OMEGA_X * (1.0 / x4 - x2 / xb6);
}

gdouble
nc_hicosmo_qg_gbarbar (NcmModel *model, gdouble x)
{
  gdouble x1_m3w = x * pow (x, -3.0 * W);
  gdouble x2 = x*x;
  gdouble x4 = x2*x2;

  gdouble xb = exp(-G_B);
  gdouble xb2 = xb*xb;
  gdouble xb3 = xb2*xb;
  gdouble xb3_m3w =  xb3 * pow (xb, -3.0 * W);
  gdouble xb6 = xb3*xb3;

  //g_assert (x < xb*1e-16);

  return OMEGA_R * x2 / xb2 + OMEGA_M * ((1.0 - 3.0 * W) / (x1_m3w * 2.0) + x2 / xb3_m3w) + OMEGA_X * (2.0 / x4 + x2 / xb6);
}

/**
 * nc_hicosmo_qg_h_to_R_matrix: (skip)
 * @model: FIXME
 * @x: FIXME
 * @T: FIXME
 * 
 * FIXME
 * 
 */
void
nc_hicosmo_qg_h_to_R_matrix (NcmModel *model, gdouble x, gsl_matrix *T)
{
  gdouble gbar2 = nc_hicosmo_qg_gbar2 (model, x);
  gdouble gbar = sqrt (gbar2);        
  gdouble F_x = nc_hicosmo_qg_xxbarzeta2 (model, x, NULL) / x;
  gdouble sqrt_x_F = sqrt (1.0 / F_x);
  gdouble xdF_F = x * nc_hicosmo_qg_dxxbarzeta2_xxbarzeta2 (model, x, NULL) / 1.0;
  
  g_assert ((T->size1 == T->size2) && (T->size1 == 2));

  gsl_matrix_set (T, 0, 0, sqrt_x_F);
  gsl_matrix_set (T, 0, 1, 0.0);
  gsl_matrix_set (T, 1, 0, -gbar * sqrt_x_F * (xdF_F - 1.0) / 2.0);
  gsl_matrix_set (T, 1, 1, gbar * sqrt_x_F);
}

/**
 * nc_hicosmo_qg_R_to_h_matrix: (skip)
 * @model: FIXME
 * @x: FIXME
 * @T: FIXME
 * 
 * FIXME
 * 
 */
void
nc_hicosmo_qg_R_to_h_matrix (NcmModel *model, gdouble x, gsl_matrix *T)
{
  nc_hicosmo_qg_h_to_R_matrix (model, x, T);
  {
    const gdouble det = gsl_matrix_get (T, 0, 0) * gsl_matrix_get (T, 1, 1) - gsl_matrix_get (T, 0, 1) * gsl_matrix_get (T, 1, 0);
    const gdouble T00 = gsl_matrix_get (T, 1, 1) / det;
    const gdouble T01 = -gsl_matrix_get (T, 0, 1) / det;
    const gdouble T10 = -gsl_matrix_get (T, 1, 0) / det;
    const gdouble T11 = gsl_matrix_get (T, 0, 0) / det;
    gsl_matrix_set (T, 0, 0, T00);
    gsl_matrix_set (T, 0, 1, T01);
    gsl_matrix_set (T, 1, 0, T10);
    gsl_matrix_set (T, 1, 1, T11);
  }
}

gdouble
nc_hicosmo_qg_d2a_n_deta2 (NcmModel *model, gdouble x)
{
  const gdouble x1_m3w = x * pow (x, -3.0 * W);
  const gdouble x2 = x*x;
  const gdouble x4 = x2*x2;
  const gdouble x5 = x4 * x;

  const gdouble xb = exp(-G_B);
  const gdouble xb2 = xb*xb;
  const gdouble xb3 = xb2*xb;
  const gdouble xb3_m3w = xb3 * pow (xb, -3.0 * W);
  const gdouble xb6 = xb3 * xb3;

  return x2 * 
    ( 
     OMEGA_R * x / xb2 + 
     OMEGA_M * ((1.0 - 3.0 * W) / (2.0 * x * x1_m3w) + x / xb3_m3w) + 
     OMEGA_X * (2.0 / x5 + x / xb6)
    );
}

gdouble
nc_hicosmo_qg_dd2a_n_deta2_da_n (NcmModel *model, gdouble x)
{
  const gdouble x1_3w = x * pow (x, 3.0 * W);
  const gdouble x2 = x*x;
  const gdouble x4 = x2*x2;

  const gdouble xb = exp(-G_B);
  const gdouble xb2 = xb * xb;
  const gdouble xb3 = xb2 * xb;
  const gdouble xb3_m3w = xb3 * pow (xb, -3.0 * W);
  const gdouble xb6 = xb3 * xb3;

  return
  - 3.0 * OMEGA_R * x4 / xb2 
  - OMEGA_M * (3.0 * W * (1.0 - 3.0 * W) * x1_3w / 2.0 + x4 / xb3_m3w)
  - OMEGA_X * 3.0 * (x4 / xb6 - 2.0 / x2);
}


gdouble
nc_hicosmo_qg_V (NcmModel *model, gdouble x, gpointer userdata)
{
  gdouble x2 = x*x;
  gdouble x3 = x2*x;
  gdouble x4 = x2*x2;
  gdouble x3w = pow (x, 3.0 * W);
  
  gdouble xb = exp(-G_B);
  gdouble xb2 = xb*xb;
  gdouble xb3 = xb2*xb;
  gdouble xb6 = xb3*xb3;
  gdouble xb3w = pow (xb, 3.0 * W);

  gdouble E = x2 * sqrt(OMEGA_R * (1.0 - x2 / xb2) + OMEGA_M*(x3w / x - x2  * xb3w / xb3) + OMEGA_X * (1.0 / x4 - x2 / xb6));
  gdouble E2 = E*E;
  gdouble V1 = E2 * nc_hicosmo_qg_ddzeta_zeta (model, x, NULL);
  gdouble V2 = 2.0 * E2 / x * nc_hicosmo_qg_dzeta_zeta (model, x, NULL);
  gdouble V3 = -(OMEGA_R * x2/xb2 + OMEGA_M * ((1.0 - 3.0 * W) * x3w / (x * 2.0) + x2  * xb3w / xb3) +OMEGA_X * (2.0 / x4 + x2 / xb6)) * 
    x3 * nc_hicosmo_qg_dzeta_zeta (model, x, NULL);

  if (gsl_finite (V1) && gsl_finite(V2))
    return V1 + V2 + V3;
  else
    return V3;
}

gdouble
nc_hicosmo_qg_lambda_x (NcmModel *model, gdouble x, gpointer userdata)
{
  gdouble sqrt_Omega_r = sqrt(OMEGA_R);
  gdouble sqrt_Omega_m = sqrt(OMEGA_M);
  gdouble one_m3omega = (1.0 - 3.0 * W);
  
  return 2.0 * asinh (pow(x, one_m3omega / 2.0) * sqrt_Omega_r / sqrt_Omega_m) / (one_m3omega * sqrt_Omega_r);
}

/**
 * nc_hicosmo_qg_past_sol: (skip)
 * @model: FIXME
 * @k: FIXME
 * @lambda: FIXME
 * @sol: FIXME
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
gboolean
nc_hicosmo_qg_past_sol (NcmModel *model, gdouble k, gdouble lambda, gsl_matrix *sol)
{
  gdouble k2 = k * k;
  gdouble one_p3w = 1.0 + 3.0 * W;
  gdouble A1 = (4.0 + one_p3w) / (2.0 * one_p3w);
  gdouble one_m3w = 1.0 - 3.0 * W;
  gdouble B = pow (lambda * one_m3w / 2.0, -2.0 * one_p3w / one_m3w);
  gdouble C = pow (OMEGA_M, -2.0 / one_m3w);
  gdouble argfunc = -k2 * W * B * C / gsl_pow_2 (one_p3w);
  gdouble s1 = gsl_sf_hyperg_0F1 (A1, argfunc);

  gdouble A2 = (-4.0 + 3.0 * one_p3w) / (2.0 * one_p3w);
  gdouble s2p = pow (lambda, 1.0 + 2.0 / one_m3w); 
  gdouble s2 = s2p * gsl_sf_hyperg_0F1 (A2, argfunc);

  gdouble ds1 = gsl_sf_hyperg_0F1 (A1 + 1.0, argfunc) * (-2.0 * one_p3w / one_m3w) * argfunc / (A1 * lambda);
  
  gdouble ds2 = 
    s2p * gsl_sf_hyperg_0F1 (A2 + 1.0, argfunc) * (-2.0 * one_p3w / one_m3w) * argfunc / (A2 * lambda) + 
    s2 * (1.0 + 2.0 / one_m3w) / lambda;

  gsl_matrix_set (sol, 0, 0, s1);
  gsl_matrix_set (sol, 0, 1, s2);
  
  gsl_matrix_set (sol, 1, 0, ds1);
  gsl_matrix_set (sol, 1, 1, ds2);
  
  return TRUE;
}

#define NC_C_Hb_H0 1.0
#define NC_C_SPEED 1.0e-5

gdouble
nc_hicosmo_qg_alphaprime2 (NcmModel *model, gdouble alpha, gpointer data)
{
  const gdouble alpha2 = alpha * alpha;
  const gdouble x = exp (-alpha2 / 2.0 - G_B);
  const gdouble x2 = x  * x;
  const gdouble x3w = pow (x, 3.0 * W);
  const gdouble three_1mW_2 = 3.0 * (1.0 - W) / 2.0;
  const gdouble speed = NC_C_SPEED;
  const gdouble Hb_H0 = NC_C_Hb_H0;
  const gdouble cc_fact = (alpha <= 0.0) ? (exp(alpha / speed) + Hb_H0) / (1.0 + exp (alpha / speed)) : (1.0 + Hb_H0 * exp(-alpha / speed)) / (1.0 + exp (-alpha / speed));
  const gdouble alphaprime2 = cc_fact * cc_fact *
    (
     OMEGA_R * x2 * gsl_sf_exprel (-alpha2) +
     OMEGA_M * x * x3w * three_1mW_2 * gsl_sf_exprel (-three_1mW_2 * alpha2) +
     OMEGA_X / x2 * 3.0 * gsl_sf_exprel (-3.0 * alpha2)
     );

  return alphaprime2;
}

gdouble
nc_hicosmo_qg_dalphaprime2_dalpha (NcmModel *model, gdouble alpha, gpointer data)
{
  const gdouble alpha2 = alpha * alpha;
  const gdouble alpha4 = alpha2 * alpha2;
  const gdouble x = exp (-alpha2 / 2.0 - G_B);
  const gdouble x2 = x  * x;
  const gdouble x3w = pow (x, 3.0 * W);
  const gdouble three_1mW_2 = 3.0 * (1.0 - W) / 2.0;
  const gdouble onep3W = 1.0 + 3.0 * W;
  const gdouble three_1mW_2alpha2 = three_1mW_2 * alpha2;
  gdouble er, em, ex, dalphaprime2_dalpha;
  if (alpha2 > -GSL_LOG_DBL_EPSILON + 3.0 * M_LN10)
  {
    er = (2.0 / alpha2 + 2.0 / alpha4);
    em = three_1mW_2 * (onep3W / three_1mW_2alpha2 + 2.0 / (alpha2 * three_1mW_2alpha2));
    ex = 3.0 * (-2.0 / (3.0 * alpha2) + 2.0 / (3.0 * alpha4));
  }
  else
  {
    er = (2.0 * gsl_sf_exprel (-alpha2) + exp (-alpha2) * gsl_sf_exprel_2 (alpha2));
    em = three_1mW_2 * (onep3W * gsl_sf_exprel (-three_1mW_2alpha2) + three_1mW_2 * exp (-three_1mW_2alpha2) * gsl_sf_exprel_2 (three_1mW_2alpha2));
    ex = 3.0 * (-2.0 * gsl_sf_exprel (-3.0 * alpha2) + 3.0 * exp (-3.0 * alpha2) * gsl_sf_exprel_2 (3.0 * alpha2));
  }

  {
    const gdouble speed = NC_C_SPEED;
    const gdouble Hb_H0 = NC_C_Hb_H0;
    const gdouble cc_fact = (alpha <= 0.0) ? (exp(alpha / speed) + Hb_H0) / (1.0 + exp (alpha / speed)) : (1.0 + Hb_H0 * exp(-alpha / speed)) / (1.0 + exp (-alpha / speed));    const gdouble dcc_fact2_dalpha = 2.0 / speed * cc_fact * (1.0 - Hb_H0) / gsl_pow_2 (2.0 * cosh (alpha / (2.0 * speed)));
    const gdouble DF = dcc_fact2_dalpha *
      (
        OMEGA_R * x2 * gsl_sf_exprel (-alpha2) +
        OMEGA_M * x * x3w * three_1mW_2 * gsl_sf_exprel (-three_1mW_2 * alpha2) +
        OMEGA_X / x2 * 3.0 * gsl_sf_exprel (-3.0 * alpha2)
        );
    
    dalphaprime2_dalpha = -alpha * cc_fact * cc_fact *
      (
        OMEGA_R * x2 * er +
        OMEGA_M * x * x3w * em +
        OMEGA_X / x2 * ex
        ) + DF;
  }

  return dalphaprime2_dalpha;
}

void
nc_hicosmo_qg_h_dust_sol (NcHICosmoQGMode *qgmode, long double x, long double ax_i, long double *h)
{
  NcmModel *model = qgmode->model;
  const long double k = qgmode->k;
  const long double x0 = ax_i;
  const long double x2 = x * x;
  const long double x3 = x2 * x;
  const long double w = W;
  const long double x3w = powl (x, 3.0L * w);
  const long double x03w = powl (x0, 3.0L * w);
  const long double Om = OMEGA_M;
  const long double F = (k / (1.0L + 3.0L * w)) * sqrtl (w / (Om));
  const long double f = F / sqrtl (x * x3w);
  const long double f0 = F / sqrtl (x0 * x03w) * 0.0;
  const long double f2 = f * f;
  const long double alpha_w = 3.0L * (1.0L - w) / (2.0L * (1.0L + 3.0L * w));
  const long double x3_m3w_4 = powl (x3 / x3w, 1.0L / 4.0L);
  const long double Om_w1_4 = powl (Om / w, 1.0L / 4.0L);
  const long double F2_1_3w = powl (F, 2.0L / (1.0L + 3.0L * w));
  const long double sqrt_pi = sqrtl (NC_C_PIL);
  const long double sin_alpha_pi = sinl (NC_C_PIL * alpha_w);
  const long double gamma_1_malpha = tgammal (1.0L - alpha_w);
  const long double gamma_1_palpha = tgammal (1.0L + alpha_w);
  const long double sqrt_k = sqrtl (k);
  long double resl, argp, argm;
  long double up, dup;
  long double um, dum;
  gsl_sf_result res;
  
  gsl_sf_hyperg_0F1_e (1.0 + alpha_w, -f2, &res);
  resl = res.val;
  up = 1.0 / (sqrt_k) * Om_w1_4 * F2_1_3w * sqrt_pi * resl / (x3_m3w_4 * gamma_1_palpha * sin_alpha_pi);
  gsl_sf_hyperg_0F1_e (1.0 - alpha_w, -f2, &res);
  resl = res.val;
  um = 1.0 / (sqrt_k) * Om_w1_4 * x3_m3w_4 * F * sqrt_pi * resl / (F2_1_3w * gamma_1_malpha * sin_alpha_pi);

  argp = -NC_C_PIL * alpha_w / 2.0L - NC_C_PIL / 4.0L - 2.0L * f0;
  argm =  NC_C_PIL * alpha_w / 2.0L - NC_C_PIL / 4.0L - 2.0L * f0;
  
  dup = ((-3.0L * (1.0L - w) / 4.0L) * up / x + (1.0L + 3.0L * w) * f2 / (x * (1.0L + alpha_w)) / (sqrt_k) * Om_w1_4 * F2_1_3w * sqrt_pi * gsl_sf_hyperg_0F1 (2.0L + alpha_w, -f2) / (x3_m3w_4 * gamma_1_palpha * sin_alpha_pi)) * x;
  dum = ((+3.0L * (1.0L - w) / 4.0L) * um / x + (1.0L + 3.0L * w) * f2 / (x * (1.0L - alpha_w)) / (sqrt_k) * Om_w1_4 * x3_m3w_4 * F * sqrt_pi * gsl_sf_hyperg_0F1 (2.0L - alpha_w, -f2) / (F2_1_3w * gamma_1_malpha * sin_alpha_pi)) * x;

  h[0] = -up * cosl (argp) + um * cosl (argm);      /* Re(h) */
  h[1] = -dup * cosl (argp) + dum * cosl (argm);    /* Re(hbar) */
  h[2] = -(-up * sinl (argp) + um * sinl (argm));   /* Im(h) */
  h[3] = -(-dup * sinl (argp) + dum * sinl (argm)); /* Im(hbar) */
  
  return;
}

void
nc_hicosmo_qg_R_dust_sol (NcHICosmoQGMode *qgmode, gdouble x, gdouble *R)
{
  NcmModel *model = qgmode->model; 
  const gdouble E0 = NC_C_Hb_H0;
  const gdouble sqrt_Omega_m = sqrt (OMEGA_M);
  const gdouble A = sqrt (M_PI / (E0 * (1.0 + 3.0 * W) * sqrt_Omega_m) * (2.0 * W) / (3.0 * (1.0 + W)));
  const gdouble alpha = 3.0 * (1.0 - W) / (2.0 * (1.0 + 3.0 * W));
  const gdouble nu_malpha = pow (x, 3.0 * (1.0 - W) / 4.0);
  const gdouble nu = pow (x, -(1.0 + 3.0 * W) / 2.0);
  const gdouble u = 2.0 * qgmode->k * sqrt(W) * nu / (E0 * (1.0 + 3.0 * W) * sqrt_Omega_m);
  const gdouble dnu_dx = - (1.0 + 3.0 * W) / 2.0 * pow (x, -3.0 * (1.0 + W) / 2.0);

  R[0] = A * nu_malpha * gsl_sf_bessel_Jnu (alpha, u);
  R[1] = -A * nu_malpha * gsl_sf_bessel_Jnu (alpha + 1.0, u) * u / nu * dnu_dx;
  
  R[2] = A * nu_malpha * gsl_sf_bessel_Ynu (alpha, u);
  R[3] = -A * nu_malpha * gsl_sf_bessel_Ynu (alpha + 1.0, u) * u / nu * dnu_dx;
}

#define LOCAL_MODEL_QG_SPLINE_FUNCTION(_func_,_param_) \
gdouble \
nc_hicosmo_qg_##_func_##_##_param_ (NcmModel *model, gdouble _param_, gint deriv) \
{ \
  NcHICosmoQG *qgint = (NcHICosmoQG *)model; \
  nc_hicosmo_qg_init_spline (qgint, model); \
  \
  if (qgint->_func_##_##_param_ == NULL) \
  { \
    qgint->_func_##_##_param_ = gsl_interp_alloc (gsl_interp_cspline, qgint->_param_->len); \
    gsl_interp_init (qgint->_func_##_##_param_,  \
                     &(g_array_index(qgint->_param_, gdouble, 0)),  \
                     &(g_array_index(qgint->_func_, gdouble, 0)),  \
                     qgint->_param_->len); \
  } \
\
  if (deriv == TRUE) \
    return gsl_interp_eval_deriv (qgint->_func_##_##_param_, \
                            &(g_array_index(qgint->_param_, gdouble, 0)), \
                            &(g_array_index(qgint->_func_, gdouble, 0)), \
                            _param_, \
                            qgint->_param_##_accel);\
  else if (deriv == FALSE) \
    return gsl_interp_eval (qgint->_func_##_##_param_, \
                            &(g_array_index(qgint->_param_, gdouble, 0)), \
                            &(g_array_index(qgint->_func_, gdouble, 0)), \
                            _param_, \
                            qgint->_param_##_accel);\
  else \
    return gsl_interp_eval_integ (qgint->_func_##_##_param_, \
                            &(g_array_index(qgint->_param_, gdouble, 0)), \
                            &(g_array_index(qgint->_func_, gdouble, 0)), \
                            qgint->lambda_i, _param_, \
                            qgint->_param_##_accel); \
}\


LOCAL_MODEL_QG_SPLINE_FUNCTION(eta, lambda)
LOCAL_MODEL_QG_SPLINE_FUNCTION(x, lambda)
LOCAL_MODEL_QG_SPLINE_FUNCTION(gbar, lambda)
LOCAL_MODEL_QG_SPLINE_FUNCTION(gbarbar, lambda)
LOCAL_MODEL_QG_SPLINE_FUNCTION(int_1_zeta2, lambda)
LOCAL_MODEL_QG_SPLINE_FUNCTION(cs2zeta2_int_1_zeta2, lambda)

LOCAL_MODEL_QG_SPLINE_FUNCTION(lambda, k_cross)

LOCAL_MODEL_QG_SPLINE_FUNCTION(cs2, lambda)
LOCAL_MODEL_QG_SPLINE_FUNCTION(V, lambda)

gdouble
nc_hicosmo_qg_get_lambda_f (NcmModel *model, gpointer userdata)
{
  NcHICosmoQG *qgint = (NcHICosmoQG *)model;
  nc_hicosmo_qg_init_spline (qgint, model);
  return qgint->last_lambda;
}

gdouble
nc_hicosmo_qg_get_lambda_i (NcmModel *model, gpointer userdata)
{
  NcHICosmoQG *qgint = (NcHICosmoQG *) model;
  nc_hicosmo_qg_init_spline (qgint, model);
  return qgint->lambda_i;
}

gdouble
nc_hicosmo_qg_get_lambda_d (NcmModel *model, gpointer userdata)
{
  NcHICosmoQG *qgint = (NcHICosmoQG *) model;
  nc_hicosmo_qg_init_spline (qgint, model);
  return qgint->lambda_d;
}


NcHICosmoQG *
nc_hicosmo_qgint_new ()
{
  NcHICosmoQG *qgint = g_slice_new (NcHICosmoQG);
  qgint->cvode = CVodeCreate(CV_BDF, CV_NEWTON);
  //qgint->cvode = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
  qgint->initialized = FALSE;
  qgint->y = N_VNew_Serial(2);
  qgint->yQ = N_VNew_Serial(2);

  qgint->lambda_accel = gsl_interp_accel_alloc ();
  qgint->eta_accel = gsl_interp_accel_alloc ();
  qgint->z_accel = gsl_interp_accel_alloc ();
  qgint->k_cross_accel = gsl_interp_accel_alloc ();

  qgint->eta_lambda = NULL;
  qgint->x_lambda = NULL;
  qgint->gbar_lambda = NULL;
  qgint->gbarbar_lambda = NULL;
  qgint->int_1_zeta2_lambda = NULL;
  qgint->cs2zeta2_int_1_zeta2_lambda = NULL;

  qgint->lambda_k_cross = NULL;

  qgint->cs2_lambda = NULL;
  qgint->V_lambda = NULL;
  
  qgint->lambda = NULL;
  qgint->eta = NULL;
  qgint->x = NULL;
  
  qgint->int_1_zeta2 = NULL;
  qgint->cs2zeta2_int_1_zeta2 = NULL;

  qgint->k_cross = NULL;
  
  qgint->g = NULL;
  qgint->gbar = NULL;
  qgint->gbarbar = NULL;

  qgint->cs2 = NULL;
  qgint->V = NULL;

  qgint->initialized_spline = FALSE;

  return qgint;
}

static gint 
scalefactor_step (realtype t, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHICosmoQGScaleFactor *qgscale = (NcHICosmoQGScaleFactor *)f_data;
  NcmModel *model = qgscale->model;
  const gdouble g = NV_Ith_S(y, 0);
  const gdouble expg = exp(g);
  const gdouble expg1_3w = exp (g * (1.0 - 3.0 * W));
  const gdouble exp2g = gsl_pow_2 (expg);
  const gdouble exp4g = gsl_pow_2 (exp2g);
  const gdouble exp2gb_2g = exp (2.0 * (G_B - g));
  const gdouble exp31mwgb_2g = exp (3.0 * (1.0 - W) * G_B - 2.0 * g);
  const gdouble exp6gb_2g = exp (6.0 * G_B - 2.0 * g);

//  printf ("%.15g %.15g %.15g\n", t, g, 1.0 / expg);
  
  if (qgscale->first_order)
  {
    NV_Ith_S(ydot, 0) = qgscale->sign * sqrt(nc_hicosmo_qg_gbar2 (model, 1.0 / expg));
    NV_Ith_S(ydot, 1) = 0.0;
  }
  else
  {
    NV_Ith_S(ydot, 0) = NV_Ith_S(y, 1);
    NV_Ith_S(ydot, 1) = OMEGA_R * exp2gb_2g + 
      OMEGA_M * (expg1_3w * (1.0 - 3.0 * W) / 2.0 + exp31mwgb_2g) + OMEGA_X * (2.0 * exp4g + exp6gb_2g);
  }
  return 0;
}

static gint 
scalefactor_step_J (gint N, realtype t, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHICosmoQGScaleFactor *qgscale = (NcHICosmoQGScaleFactor *)jac_data;
  NcmModel *model = qgscale->model;
  gdouble g = NV_Ith_S(y, 0);
  gdouble expg = exp(g);
  gdouble expg1_3w = exp (g * (1.0 - 3.0 * W));
  gdouble exp2g = gsl_pow_2 (expg);
  gdouble exp4g = gsl_pow_2 (exp2g);
  gdouble exp2gb_2g = exp (2.0 * (G_B-g));
  gdouble exp31mwgb_2g = exp (3.0 * (1.0 - W) * G_B - 2.0 * g);
  gdouble exp6gb_2g = exp (6.0 * G_B-2.0 * g);

  if (qgscale->first_order)
  {
    gdouble x = 1.0 / expg;
    DENSE_ELEM (J, 0, 0) = qgscale->sign * nc_hicosmo_qg_gbarbar (model, x) / sqrt(nc_hicosmo_qg_gbar2 (model, x));
    DENSE_ELEM (J, 0, 1) = 0.0;

    DENSE_ELEM (J, 1, 0) = 0.0;
    DENSE_ELEM (J, 1, 1) = 0.0;
  }
  else
  {
    DENSE_ELEM (J, 0, 0) = 0.0;
    DENSE_ELEM (J, 0, 1) = 1.0;

    DENSE_ELEM (J, 1, 0) = -2.0 * OMEGA_R * exp2gb_2g + OMEGA_M * (gsl_pow_2 (1.0 - 3.0 * W) * expg1_3w / 2.0 - 2.0 * exp31mwgb_2g) + OMEGA_X * (8.0 * exp4g - 2.0 * exp6gb_2g);
    DENSE_ELEM (J, 1, 1) = 0.0;
  }
  return 0;
}

gint 
scale_factor_root (gdouble t, N_Vector y, gdouble *gout, gpointer g_data)
{
  NcHICosmoQGScaleFactor *qgscale = (NcHICosmoQGScaleFactor *)g_data;
  gout[0] = NV_Ith_S(y, 0) - qgscale->x_root_scale;
  gout[1] = NV_Ith_S(y, 1);
  gout[2] = NV_Ith_S(y, 0) - 10.0 * M_LN10;
  return 0;
}

int 
scale_factor_time (gdouble t, N_Vector y, N_Vector yQdot, gpointer fQ_data)
{
  NcHICosmoQGScaleFactor *qgscale = (NcHICosmoQGScaleFactor *)fQ_data;
  NcmModel *model = qgscale->model;

  gdouble x = exp(-NV_Ith_S(y, 0));
  gdouble zeta = nc_hicosmo_qg_zeta (model, x, NULL);

  NV_Ith_S(yQdot, 0) = 1.0 / x;
  NV_Ith_S(yQdot, 1) = 1.0 / (x * zeta * zeta);
  return 0;
}

gboolean
nc_hicosmo_qgint_init (NcHICosmoQG *qgint, gdouble lambda, gdouble g, gdouble gbar, gdouble eta, gdouble int_1_zeta2)
{
  gint flag;

  NV_Ith_S(qgint->y, 0) = g;
  NV_Ith_S(qgint->y, 1) = gbar;
  
  NV_Ith_S(qgint->yQ, 0) = eta;
  NV_Ith_S(qgint->yQ, 1) = int_1_zeta2;
  qgint->reltol = 1e-11;
  qgint->abstol = 1e-250;

  if (!qgint->initialized)
  {
    flag = CVodeInit (qgint->cvode, &scalefactor_step, lambda, qgint->y); 
    CVODE_CHECK (&flag, "CVodeMalloc", 1, FALSE);
    CVodeQuadInit (qgint->cvode, scale_factor_time, qgint->yQ);
    qgint->initialized = TRUE;
  }
  else
  {
    flag = CVodeReInit (qgint->cvode, lambda, qgint->y);
    CVodeQuadReInit (qgint->cvode, qgint->yQ);
    CVODE_CHECK(&flag, "CVodeReInit", 1, FALSE);
  }
  
  return TRUE;
}

gboolean
nc_hicosmo_qgint_set_opts (NcHICosmoQG *qgint, NcHICosmoQGScaleFactor *qgscale)
{
  gint flag;

  flag = CVodeSStolerances (qgint->cvode, qgint->reltol, qgint->abstol);
  CVODE_CHECK(&flag, "CVodeSStolerances", 1, FALSE);
  
  flag = CVodeSetMaxNumSteps(qgint->cvode, 1000000);
  CVODE_CHECK(&flag, "CVodeSetMaxNumSteps", 1, FALSE);
  
  flag = CVodeSetUserData (qgint->cvode, qgscale);
  CVODE_CHECK(&flag, "CVodeSetFdata", 1, FALSE);
  
  if (TRUE)
  {
    flag = CVDense(qgint->cvode, 2);
    CVODE_CHECK(&flag, "CVDense", 1, FALSE);
    
    flag = CVDlsSetDenseJacFn (qgint->cvode, &scalefactor_step_J);
    CVODE_CHECK(&flag, "CVDlsSetDenseJacFn", 1, FALSE);
  }

//  flag = CVodeSetStabLimDet(qgint->cvode, TRUE);
//  CVODE_CHECK(&flag, "CVodeSetStabLimDet", 1, FALSE);
  
//  flag = CVodeSetMaxStep(qgint->cvode, 1.0); 
//  CVODE_CHECK(&flag, "CVodeSetMaxStep", 1, FALSE);        

//  flag = CVodeSetInitStep (qgint->cvode, 1e-15); 
//  CVODE_CHECK(&flag, "CVodeSetMaxStep", 1, FALSE);
  
  flag = CVodeRootInit (qgint->cvode, 3, scale_factor_root);
  CVODE_CHECK(&flag, "CVodeRootInit", 1, FALSE); 
  
//  flag = CVodeSetStopTime(qgint->cvode, qgint->tf);
//  CVODE_CHECK(&flag, "CVodeSetStopTime", 1, FALSE);

  return TRUE;
}

gboolean
nc_hicosmo_qg_init_spline (NcHICosmoQG *qgint, NcmModel *model)
{
  if (!qgint->initialized_spline)
  {
    gdouble lambda_i, temp;
    gdouble x_i = 1.0e-30;
    gdouble g_i = -log(x_i);
    gdouble gbar2 = nc_hicosmo_qg_gbar2 (model, x_i);
    gdouble gbar = -sqrt(gbar2);
    gdouble gbarbar = nc_hicosmo_qg_gbarbar (model, x_i);
    gdouble cs2, V, k_cross;
    NcHICosmoQGScaleFactor qgscale = {model, TRUE, -1.0};
    gboolean end_dust = FALSE;
    gdouble lambda_iQ;
    
    N_Vector dstep = N_VNew_Serial (2);

#define LOCAL_INIT_ARRAY(ArrAy) \
    if ((ArrAy) == NULL) ArrAy = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), 2000);  \
    else g_array_set_size (ArrAy, 0)
    
    LOCAL_INIT_ARRAY (qgint->lambda);
    LOCAL_INIT_ARRAY (qgint->eta);
    LOCAL_INIT_ARRAY (qgint->x);

    LOCAL_INIT_ARRAY (qgint->int_1_zeta2);
    LOCAL_INIT_ARRAY (qgint->cs2zeta2_int_1_zeta2);

    LOCAL_INIT_ARRAY (qgint->k_cross);

    LOCAL_INIT_ARRAY (qgint->g);
    LOCAL_INIT_ARRAY (qgint->gbar);
    LOCAL_INIT_ARRAY (qgint->gbarbar);

    LOCAL_INIT_ARRAY (qgint->cs2);
    LOCAL_INIT_ARRAY (qgint->V);
    /*
    qgscale.x_root_scale = -15.0 * M_LN10;
    qgscale.first_order = FALSE;
    nc_hicosmo_qgint_init (qgint, 0.0, G_B, 0.0, 0.0, 0.0);
    nc_hicosmo_qgint_set_opts (qgint, &qgscale);
    lambda_i = -1e-30;
    while (1)
    {
      gint flag;
      gdouble x, zeta, xdzeta_zeta, cs2zeta2_int_1_zeta2;

      flag = CVode (qgint->cvode, lambda_i * 1.001, qgint->y, &lambda_i, CV_ONE_STEP);
      CVODE_CHECK (&flag, "CVode[nc_hicosmo_qg_init_spline]", 1, FALSE);
      CVodeGetQuad (qgint->cvode, &lambda_iQ, qgint->yQ); 

      printf ("% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", lambda_i, NV_Ith_S (qgint->y, 0), NV_Ith_S (qgint->y, 1), NV_Ith_S (qgint->yQ, 0), exp(-NV_Ith_S (qgint->y, 0)));
      if (flag == CV_ROOT_RETURN)
      {
        gint roots[3];
        CVodeGetRootInfo (qgint->cvode, roots);
        if (roots[0])
          break;
      }
    }
    qgscale.x_root_scale = -2.0 * M_LN10;
    nc_hicosmo_qgint_init (qgint, 1.0, NV_Ith_S (qgint->y, 0), NV_Ith_S (qgint->y, 1), NV_Ith_S (qgint->yQ, 0), 0.0);
    nc_hicosmo_qgint_set_opts (qgint, &qgscale);
    lambda_i = 1.0;
    while (1)
    {
      gint flag;
      gdouble x, zeta, xdzeta_zeta, cs2zeta2_int_1_zeta2;

      flag = CVode (qgint->cvode, lambda_i * 1.001, qgint->y, &lambda_i, CV_ONE_STEP);
      CVODE_CHECK (&flag, "CVode[nc_hicosmo_qg_init_spline]", 1, FALSE);
      CVodeGetQuad (qgint->cvode, &lambda_iQ, qgint->yQ); 

      printf ("% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", lambda_i, NV_Ith_S (qgint->y, 0), NV_Ith_S (qgint->y, 1), NV_Ith_S (qgint->yQ, 0), exp(-NV_Ith_S (qgint->y, 0)));
      if (flag == CV_ROOT_RETURN)
      {
        gint roots[3];
        CVodeGetRootInfo (qgint->cvode, roots);
//        if (roots[2])
//          break;
      }
    }
*/

    
    
    qgint->lambda_i = nc_hicosmo_qg_lambda_x (model, x_i, NULL);
    g_array_append_val (qgint->lambda,  qgint->lambda_i);
    
    temp = 0.0; g_array_append_val (qgint->eta,  temp);
    temp = 0.0; g_array_append_val (qgint->int_1_zeta2, temp);
    temp = 0.0; g_array_append_val (qgint->cs2zeta2_int_1_zeta2, temp);

    g_array_append_val (qgint->x, x_i);
    g_array_append_val (qgint->g, g_i);
    g_array_append_val (qgint->gbar, gbar);
    g_array_append_val (qgint->gbarbar, gbarbar);

    cs2 = nc_hicosmo_qg_cs2 (model, x_i, NULL);
    V = nc_hicosmo_qg_V (model, x_i, NULL);
    g_array_append_val (qgint->cs2, cs2);
    g_array_append_val (qgint->V, V);

    k_cross = sqrt (V / cs2);
    g_array_append_val (qgint->k_cross, k_cross);






    
    nc_hicosmo_qgint_init (qgint, nc_hicosmo_qg_lambda_x (model, x_i, NULL), g_i, gbar, 0.0, 0.0);
    nc_hicosmo_qgint_set_opts (qgint, &qgscale);

    printf ("# INIT COND: %g %g %g %g %g\n", qgint->lambda_i, x_i, g_i, gbar, gbarbar);
    lambda_i = qgint->lambda_i;
    while (1)
    {
      gint flag;
      gdouble x, zeta, xdzeta_zeta, cs2zeta2_int_1_zeta2;

      flag = CVode (qgint->cvode, lambda_i * 1.001, qgint->y, &lambda_i, CV_NORMAL);
      CVODE_CHECK (&flag, "CVode[nc_hicosmo_qg_init_spline]", 1, FALSE);
      CVodeGetQuad (qgint->cvode, &lambda_iQ, qgint->yQ); 
      CVodeGetDky (qgint->cvode, lambda_i, 1, dstep); 

      g_array_append_val (qgint->lambda, lambda_i);
      g_array_append_val (qgint->int_1_zeta2, NV_Ith_S (qgint->yQ, 1));
      g_array_append_val (qgint->eta, NV_Ith_S (qgint->yQ, 0));
    
      x = exp(-NV_Ith_S (qgint->y, 0));
      g_array_append_val (qgint->x, x);
      g_array_append_val (qgint->g, NV_Ith_S (qgint->y, 0));
      g_array_append_val (qgint->gbar, NV_Ith_S (dstep, 0));
      g_array_append_val (qgint->gbarbar, NV_Ith_S (dstep, 1));

      cs2 = nc_hicosmo_qg_cs2 (model, x, NULL);
      V = nc_hicosmo_qg_V (model, x, NULL);
      g_array_append_val (qgint->cs2, cs2);
      g_array_append_val (qgint->V, V);

      zeta = nc_hicosmo_qg_zeta (model, x, NULL);
      cs2zeta2_int_1_zeta2 = cs2 * zeta * zeta / x * NV_Ith_S (qgint->yQ, 1);
      g_array_append_val (qgint->cs2zeta2_int_1_zeta2, cs2zeta2_int_1_zeta2);

      xdzeta_zeta = x * nc_hicosmo_qg_dzeta_zeta (model, x, NULL);
      
      if (!end_dust && (fabs(xdzeta_zeta+1) > 1e-14 || fabs((cs2-W)/W) > 1e-14))
      {
        qgint->lambda_d = lambda_i;
        end_dust = TRUE;
      }
      
      k_cross = sqrt(V / cs2);
      if (k_cross < 1.0)
        g_array_append_val (qgint->k_cross, k_cross);

//      printf ("%.15g %.15g %.15g %.15g %.15g\n", lambda_i, NV_Ith_S (qgint->yQ, 0), cs2, V, x);
//      printf ("%.15g %.15g %.15g %.15g\n", lambda_i, x, NV_Ith_S (qgint->y, 0), NV_Ith_S (qgint->yQ, 0));
//      printf ("%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n", lambda_i, x, NV_Ith_S (qgint->y, 0), NV_Ith_S (qgint->y, 1), sqrt(nc_hicosmo_qg_gbar2 (data, x)),fabs((fabs(NV_Ith_S (qgint->y, 1)) - sqrt(nc_hicosmo_qg_gbar2 (data, x)))/NV_Ith_S (qgint->y, 1)), NV_Ith_S (dstep, 1), nc_hicosmo_qg_cs2 (data, x), nc_hicosmo_qg_V (data, x));

      if (flag == CV_ROOT_RETURN)
      {
        gint roots[3];
        CVodeGetRootInfo (qgint->cvode, roots);

        if (roots[0] && qgscale.first_order && (qgscale.sign < 0.0))
        {
          qgscale.first_order = FALSE;
          nc_hicosmo_qgint_init (qgint, lambda_i, NV_Ith_S (qgint->y, 0), NV_Ith_S (dstep, 0), NV_Ith_S (qgint->yQ, 0), NV_Ith_S (qgint->yQ, 1));
          nc_hicosmo_qgint_set_opts (qgint, &qgscale);
          
          printf ("# [%.15g] g = %.15g, gbar = %.15g | %.15g (%.15f %.15f)\n", lambda_i, NV_Ith_S (qgint->y, 0), NV_Ith_S (qgint->y, 1), sqrt(nc_hicosmo_qg_gbar2 (model, 1.0)), fabs(nc_hicosmo_qg_gbar2 (model, 1.0) - gsl_pow_2 (NV_Ith_S (qgint->y, 1))), (gbar2 - nc_hicosmo_qg_gbar2 (model, 1.0))*1e-12);
        }
        
        if (roots[0] && (NV_Ith_S (qgint->y, 1) > 0.0))
        {
          qgscale.first_order = TRUE;
          qgscale.sign = 1.0;
          nc_hicosmo_qgint_init (qgint, lambda_i, NV_Ith_S (qgint->y, 0), 0.0, NV_Ith_S (qgint->yQ, 0), NV_Ith_S (qgint->yQ, 1));
          nc_hicosmo_qgint_set_opts (qgint, &qgscale);
          qgint->last_lambda = lambda_i;
          qgint->last_eta = NV_Ith_S (qgint->yQ, 0);
          printf ("# [%.15g] g = %.15g, gbar = %.15g | %.15g (%.15f %.15f)\n", lambda_i, NV_Ith_S (qgint->y, 0), NV_Ith_S (qgint->y, 1), sqrt(nc_hicosmo_qg_gbar2 (model, 1.0)), fabs(nc_hicosmo_qg_gbar2 (model, 1.0) - gsl_pow_2 (NV_Ith_S (qgint->y, 1))), (gbar2 - nc_hicosmo_qg_gbar2 (model, 1.0))*1e-12);
          break;
        }
        if (roots[1])
        {
          qgint->lambda_b = lambda_i;
          qgint->eta_b = NV_Ith_S (qgint->yQ, 0);
          qgint->last_lambda = 2.0 * lambda_i;
          printf ("# [%.15g] g = %.15g, gbar = %.15g | %.15g (%.15f %.15f)\n", lambda_i, NV_Ith_S (qgint->y, 0), NV_Ith_S (qgint->y, 1), sqrt(nc_hicosmo_qg_gbar2 (model, 1.0)), fabs(nc_hicosmo_qg_gbar2 (model, 1.0) - gsl_pow_2 (NV_Ith_S (qgint->y, 1))), (gbar2 - nc_hicosmo_qg_gbar2 (model, 1.0))*1e-12);
        }
        
        if (roots[2])
        {
          if (NV_Ith_S (qgint->y, 1) >= 0)
          {
            printf ("# [%.15g] g = %.15g, gbar = %.15g | %.15g (%.15f %.15f)\n", lambda_i, NV_Ith_S (qgint->y, 0), NV_Ith_S (qgint->y, 1), sqrt(nc_hicosmo_qg_gbar2 (model, 1.0)), fabs(nc_hicosmo_qg_gbar2 (model, 1.0) - gsl_pow_2 (NV_Ith_S (qgint->y, 1))), (gbar2 - nc_hicosmo_qg_gbar2 (model, 1.0))*1e-12);
            qgint->last_lambda = lambda_i;
            qgint->last_eta = NV_Ith_S (qgint->yQ, 0);
            break;
          }
        }
      }
    }

//    printf ("%.15g => %.15g\n", g_array_index(qgint->eta, gdouble, 0), g_array_index(qgint->cs2, gdouble, 0));
//    printf ("%.15g => %.15g\n", g_array_index(qgint->eta, gdouble, 1), g_array_index(qgint->cs2, gdouble, 1));
//    printf ("%.15g => %.15g\n", g_array_index(qgint->eta, gdouble, 2), g_array_index(qgint->cs2, gdouble, 2));
//    printf ("%.15g => %.15g\n", g_array_index(qgint->eta, gdouble, 3), g_array_index(qgint->cs2, gdouble, 3));
    N_VDestroy (dstep);
    qgint->initialized_spline = TRUE;
  }
  
  return TRUE;
}

void
nc_hicosmo_qg_max_z (NcmModel *model, gdouble *max, gdouble *trans)
{
  gdouble x1 = -123.0, x2 = -132.0, x3 = -312.0;
  gdouble G_B2;

  G_B2 = gsl_pow_2 (G_B);
  
  gsl_poly_solve_cubic (0.0, -G_B2, -OMEGA_M / OMEGA_R * G_B2, &x1, &x2, &x3);

  if (x3 == -312.0)
    *max = x1-1.0;
  else
    *max = x3-1.0;
  
  x1 = -123.0;
  x2 = -132.0; 
  x3 = -312.0;
  
  gsl_poly_solve_cubic (0.0, -G_B2 * 4.0 / 6.0, -OMEGA_M * 3.0 / (6.0 * OMEGA_R) * G_B2, &x1, &x2, &x3);
  
  if (x3 == -312.0)
    *trans = x1-1.0;
  else
    *trans = x3-1.0;

  printf ("# => %g %g\n", *max, *trans);

  return;
} 

gdouble
nc_hicosmo_qg_get_eta_b (NcmModel *model, gpointer userdata)
{
  return NC_C_c  * NC_C_MEGA_PARSEC / (sqrt(OMEGA_R) * exp(-G_B) * MACRO_H0 * 1.0e3);
}

gdouble
nc_hicosmo_qg_mode_EA2 (NcHICosmoQGMode *qgmode, gdouble x)
{
  gdouble cs2 = nc_hicosmo_qg_cs2 (qgmode->model, x, NULL);
  gdouble V = nc_hicosmo_qg_V (qgmode->model, x, NULL);
  gdouble EAk2 = cs2 * qgmode->k * qgmode->k - V;
  return EAk2;
}

#define _PERT_SYS_DIM 4
static gint Rk_step (realtype t, N_Vector y, N_Vector ydot, gpointer f_data);
static gint Rk_step_J (gint N, realtype t, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static gint Rak_step (realtype t, N_Vector y, N_Vector ydot, gpointer f_data);
static gint Rak_step_J (gint N, realtype t, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static gint hk_step (realtype lambda, N_Vector y, N_Vector ydot, gpointer f_data);
static gint hk_step_J (gint N, realtype t, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/**
 * nc_hicosmo_qg_pert_new: (skip)
 * @model: FIXME
 * @ax_i: FIXME
 * @x_i: FIXME
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
NcHICosmoQGMode *
nc_hicosmo_qg_pert_new (NcmModel *model, gdouble ax_i, gdouble x_i)
{
  NcHICosmoQGPertType type = NC_HICOSMO_QG_PERT_H;
  NcHICosmoQGMode *qgmode = g_slice_new (NcHICosmoQGMode);
  qgmode->cvode_R = CVodeCreate(CV_BDF, CV_NEWTON);
  qgmode->cvode_h = CVodeCreate(CV_BDF, CV_NEWTON);
  //qgmode->cvode = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
  qgmode->type = type;
  qgmode->init_R = FALSE;
  qgmode->init_h = FALSE;
  qgmode->y = N_VNew_Serial(_PERT_SYS_DIM);
  qgmode->model = model;
  qgmode->ax_i = ax_i;
  qgmode->x_i = x_i;

  switch (type)
  {
    case NC_HICOSMO_QG_PERT_CURVATURE:
      qgmode->cvode = qgmode->cvode_R;
      qgmode->f = Rk_step;
      qgmode->jac = Rk_step_J;
      qgmode->t_i = nc_hicosmo_qg_lambda_x (model, x_i, NULL);
      qgmode->initialized = &qgmode->init_R;
      break;
    case NC_HICOSMO_QG_PERT_H:
      qgmode->cvode = qgmode->cvode_h;
      qgmode->f = hk_step;
      qgmode->jac = hk_step_J;
      qgmode->t_i = log (x_i);
      qgmode->initialized = &qgmode->init_h;
      break;
    default:
      g_assert_not_reached ();
  }

  return qgmode;
}

gboolean
nc_hicosmo_qg_pert_set_opts (NcHICosmoQGMode *qgmode)
{
  gint flag;

  flag = CVodeSStolerances (qgmode->cvode, qgmode->reltol, qgmode->abstol);
  CVODE_CHECK(&flag, "CVodeSStolerances", 1, FALSE);
  
  flag = CVodeSetMaxNumSteps (qgmode->cvode, 1000000);
  CVODE_CHECK(&flag, "CVodeSetMaxNumSteps", 1, FALSE);

  flag = CVodeSetUserData (qgmode->cvode, qgmode);
  CVODE_CHECK(&flag, "CVodeSetUserData", 1, FALSE);

  if (TRUE)
  {
    flag = CVDense (qgmode->cvode, _PERT_SYS_DIM);
    CVODE_CHECK(&flag, "CVDense", 1, FALSE);

    flag = CVDlsSetDenseJacFn (qgmode->cvode, qgmode->jac);
    CVODE_CHECK (&flag, "CVDlsSetDenseJacFn", 1, FALSE);
  }

  return TRUE;
}

gboolean
nc_hicosmo_qg_pert_init (NcHICosmoQGMode *qgmode, gdouble k)
{
  long double h[4];
  gint flag;

  qgmode->k = k;
  nc_hicosmo_qg_h_dust_sol (qgmode, qgmode->x_i, qgmode->ax_i, h);
  switch (qgmode->type)
  {
    case NC_HICOSMO_QG_PERT_CURVATURE:
    {
      gdouble Tdata[4];
      gsl_matrix_view Tview = gsl_matrix_view_array (Tdata, 2, 2);
      gsl_matrix *T = &Tview.matrix;
      nc_hicosmo_qg_h_to_R_matrix (qgmode->model, qgmode->x_i, T);
      NV_Ith_S(qgmode->y, 0) = gsl_matrix_get (T, 0, 0) * h[0] + gsl_matrix_get (T, 0, 1) * h[1];
      NV_Ith_S(qgmode->y, 1) = gsl_matrix_get (T, 1, 0) * h[0] + gsl_matrix_get (T, 1, 1) * h[1];
      NV_Ith_S(qgmode->y, 2) = gsl_matrix_get (T, 0, 0) * h[2] + gsl_matrix_get (T, 0, 1) * h[3];
      NV_Ith_S(qgmode->y, 3) = gsl_matrix_get (T, 1, 0) * h[2] + gsl_matrix_get (T, 1, 1) * h[3];
      break;
    }
    case NC_HICOSMO_QG_PERT_H:
    {
      NV_Ith_S(qgmode->y, 0) = h[0];
      NV_Ith_S(qgmode->y, 1) = h[1];
      NV_Ith_S(qgmode->y, 2) = h[2];
      NV_Ith_S(qgmode->y, 3) = h[3];
    }
    default:
      break;
  }

  qgmode->reltol = 1e-11;
  qgmode->abstol = 0.0;

  qgmode->t = qgmode->t_i;
  if (!(*qgmode->initialized))
  {
    flag = CVodeInit (qgmode->cvode, qgmode->f, qgmode->t_i, qgmode->y);
    CVODE_CHECK(&flag, "CVodeMalloc", 1, FALSE);
    *qgmode->initialized = TRUE;
  }
  else
  {
    flag = CVodeReInit (qgmode->cvode, qgmode->t_i, qgmode->y);
    CVODE_CHECK(&flag, "CVodeReInit", 1, FALSE);
  }
  
  nc_hicosmo_qg_pert_set_opts (qgmode);
  return TRUE;
}

void
nc_hicosmo_qg_pert_R_to_h (NcHICosmoQGMode *qgmode, gdouble x, gdouble *R)
{
  gdouble tmp[4];
  gdouble Tdata[4];
  gsl_matrix_view Tview = gsl_matrix_view_array (Tdata, 2, 2);
  gsl_matrix *T = &Tview.matrix;
  nc_hicosmo_qg_R_to_h_matrix (qgmode->model, x, T);

  tmp[0] = gsl_matrix_get (T, 0, 0) * R[0] + gsl_matrix_get (T, 0, 1) * R[1];
  tmp[1] = gsl_matrix_get (T, 1, 0) * R[0] + gsl_matrix_get (T, 1, 1) * R[1];
  tmp[2] = gsl_matrix_get (T, 0, 0) * R[2] + gsl_matrix_get (T, 0, 1) * R[3];
  tmp[3] = gsl_matrix_get (T, 1, 0) * R[2] + gsl_matrix_get (T, 1, 1) * R[3];
  R[0] = tmp[0]; R[1] = tmp[1]; R[2] = tmp[2]; R[3] = tmp[3]; 
}

void
nc_hicosmo_qg_pert_h_to_R (NcHICosmoQGMode *qgmode, gdouble x, gdouble *h)
{
  gdouble tmp[4];
  gdouble Tdata[4];
  gsl_matrix_view Tview = gsl_matrix_view_array (Tdata, 2, 2);
  gsl_matrix *T = &Tview.matrix;
  nc_hicosmo_qg_h_to_R_matrix (qgmode->model, x, T);

  tmp[0] = gsl_matrix_get (T, 0, 0) * h[0] + gsl_matrix_get (T, 0, 1) * h[1];
  tmp[1] = gsl_matrix_get (T, 1, 0) * h[0] + gsl_matrix_get (T, 1, 1) * h[1];
  tmp[2] = gsl_matrix_get (T, 0, 0) * h[2] + gsl_matrix_get (T, 0, 1) * h[3];
  tmp[3] = gsl_matrix_get (T, 1, 0) * h[2] + gsl_matrix_get (T, 1, 1) * h[3];
  h[0] = tmp[0]; h[1] = tmp[1]; h[2] = tmp[2]; h[3] = tmp[3]; 
}

gboolean
nc_hicosmo_qg_pert_switch (NcHICosmoQGMode *qgmode, NcHICosmoQGPertType type, gdouble t)
{
  gint flag;

  switch (type)
  {
    case NC_HICOSMO_QG_PERT_CURVATURE:
    {
      qgmode->cvode = qgmode->cvode_R;
      qgmode->f = Rk_step;
      qgmode->jac = Rk_step_J;
      qgmode->initialized = &qgmode->init_R;        
      break;
    }
    case NC_HICOSMO_QG_PERT_H:
    {
      qgmode->cvode = qgmode->cvode_h;
      qgmode->f = hk_step;
      qgmode->jac = hk_step_J;
      qgmode->initialized = &qgmode->init_h;
      break;
    }
    default:
      g_assert_not_reached ();
      break;
  }
  qgmode->type = type;
  
  qgmode->reltol = 1e-11;
  qgmode->abstol = 1e-240;
  qgmode->t = t;

  if (!(*qgmode->initialized))
  {
    flag = CVodeInit (qgmode->cvode, qgmode->f, qgmode->t, qgmode->y);
    CVODE_CHECK(&flag, "CVodeMalloc", 1, FALSE);
    *qgmode->initialized = TRUE;
  }
  else
  {
    flag = CVodeReInit (qgmode->cvode, qgmode->t, qgmode->y);
    CVODE_CHECK(&flag, "CVodeReInit", 1, FALSE);
  }
  
  nc_hicosmo_qg_pert_set_opts (qgmode);
  return TRUE;
}

/*
gboolean
nc_hicosmo_qg_pert_evolve (NcHICosmoQGMode *qgmode)
{
  NcHICosmoQG *qgint = (NcHICosmoQG *)qgmode->cp->model->private_data;
  gint flag;
  gdouble xeq = NC_HICOSMO_Omega_m (qgmode->cp) / NC_HICOSMO_Omega_r (qgmode->cp);
  gdouble alpha_eq = log (xeq);
  gdouble lambda_eq;
  gdouble alpha, last_alpha, x, tmp[4], lambda, delta_phi, zeta;
  gint i;
  gdouble dalpha = 2.0 * fabs(qgmode->t) / 1.e5;

  nc_hicosmo_qg_init_spline (qgint, qgmode->cp);

  printf ("# Going to %.15e %.15e\n", exp(qgmode->t), xeq);

  flag = CVodeSetStopTime(qgmode->cvode, alpha_eq);
  CVODE_CHECK(&flag, "CVodeSetStopTime", 1, FALSE);

  last_alpha = qgmode->t;
  while (1)
  {
    flag = CVode (qgmode->cvode, alpha_eq, qgmode->y, &alpha, CV_ONE_STEP);
    CVODE_CHECK(&flag, "CVode[nc_hicosmo_qg_pert_evolve]", 1, FALSE);
    
    if (TRUE)
    {
      if (alpha >= last_alpha + dalpha)
      {
        last_alpha = alpha;
        x = exp (alpha);
        for (i = 0; i < 4; i++)
          tmp[i] = NV_Ith_S(qgmode->y, i);
        nc_hicosmo_qg_pert_h_to_R (qgmode, x, tmp);
        delta_phi = nc_hicosmo_qg_pert_powerspectrum (qgmode, x, tmp);
        zeta = nc_hicosmo_qg_zeta (qgmode->cp, x, NULL);
        printf ("% 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15Le\n", x, tmp[0], tmp[1], tmp[2], tmp[3], delta_phi, zeta, qgmode->k);
      }
    }
    if (alpha_eq == alpha)
      break;
  }

  x = exp (alpha);
  lambda_eq = nc_hicosmo_qg_lambda_x (qgmode->cp, exp (alpha), NULL);
  
  for (i = 0; i < 4; i++)
    tmp[i] = NV_Ith_S(qgmode->y, i);
  nc_hicosmo_qg_pert_h_to_R (qgmode, x, tmp);
  for (i = 0; i < 4; i++)
    NV_Ith_S(qgmode->y, i) = tmp[i];
  
  nc_hicosmo_qg_pert_switch (qgmode, NC_HICOSMO_QG_PERT_CURVATURE, lambda_eq);

  flag = CVodeSetStopTime(qgmode->cvode, qgint->last_lambda);
  CVODE_CHECK(&flag, "CVodeSetStopTime", 1, FALSE);

  while (1)
  {
    flag = CVode (qgmode->cvode, qgint->last_lambda, qgmode->y, &lambda, CV_ONE_STEP);
    CVODE_CHECK(&flag, "CVode[nc_hicosmo_qg_pert_evolve]", 1, FALSE);
    if (TRUE)
    {
      x = nc_hicosmo_qg_x_lambda (qgmode->cp, lambda, FALSE);
      alpha = log (x);
      if ( ((lambda <= qgint->lambda_b) && (alpha >= last_alpha + dalpha)) ||
           ((lambda >= qgint->lambda_b) && (alpha <= last_alpha - dalpha))
          )
      {
        last_alpha = alpha;
        delta_phi = nc_hicosmo_qg_pert_powerspectrum (qgmode, x, N_VGetArrayPointer (qgmode->y));
        zeta = nc_hicosmo_qg_zeta (qgmode->cp, x, NULL);
        printf ("% 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15Le\n", x, NV_Ith_S(qgmode->y, 0), NV_Ith_S(qgmode->y, 1), NV_Ith_S(qgmode->y, 2), NV_Ith_S(qgmode->y, 3), delta_phi, zeta, qgmode->k);
      }
    }
    if (qgint->last_lambda == lambda)
      break;
  }

  nc_hicosmo_qg_pert_switch (qgmode, NC_HICOSMO_QG_PERT_H, lambda_eq);
  return TRUE;
}
*/

#define OLD_CODE FALSE

/**
 * nc_hicosmo_qg_modefunc: (skip)
 * @model: FIXME
 * @k: FIXME
 * @x0: FIXME
 * @xf: FIXME
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
NcHICosmoQGMode *
nc_hicosmo_qg_modefunc (NcmModel *model, long double k, long double x0, long double xf)
{
  NcHICosmoQGMode *qgmode = g_slice_new (NcHICosmoQGMode);
  //qgmode->cvode = CVodeCreate (CV_BDF, CV_NEWTON);
  qgmode->cvode = CVodeCreate (CV_ADAMS, CV_FUNCTIONAL);
  qgmode->initialized = FALSE;
  qgmode->y = N_VNew_Serial(4);
  qgmode->yQ = N_VNew_Serial(4);
  qgmode->k = k;
  qgmode->alpha0 = logl(x0);
  qgmode->alphai = qgmode->alpha0;
  qgmode->alphaf = logl(xf);
  qgmode->model = model;
  qgmode->pw_spline = NULL;
  qgmode->initialized = &qgmode->init_h;
  qgmode->init_h = FALSE;

  return qgmode;
}

gboolean
nc_hicosmo_qg_pert_prepare_pw_spline (NcHICosmoQGMode *qgmode, gboolean verbose)
{
  gint i;
  GTimer *bench = NULL;
  gdouble part_time = 0.0;

  if (qgmode->pw_spline == NULL)
  {
    GArray *x = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), 100);
    GArray *y = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), 100);    
    g_array_set_size (x, 100);
    g_array_set_size (y, 100);
    qgmode->pw_spline = ncm_spline_gsl_new (gsl_interp_cspline);
		ncm_spline_set_array (qgmode->pw_spline, x ,y, FALSE);
    part_time = 0.0;
  }

  if (verbose)
  {
    bench = g_timer_new ();
    printf ("# Starting mode calculation\n");
  }
  else
    printf ("# Starting mode calculation\n#");
  
  for (i = 0; i < 100; i++)
  {
    qgmode->k = powl (10.0L, -4.0L + 8.0L / 99.0L * i);
    nc_hicosmo_qg_modefunc_init (qgmode);
    nc_hicosmo_qg_modefunc_evolve (qgmode);
		ncm_vector_set (qgmode->pw_spline->xv, i, log10(qgmode->k));
		ncm_vector_set (qgmode->pw_spline->yv, i, log10(nc_hicosmo_qg_pert_powerspectrum (qgmode, exp(qgmode->alphai), N_VGetArrayPointer (qgmode->y))));

//printf ("%.15g %.15g\n", qgmode->pw_spline->x[i], qgmode->pw_spline->y[i]);
    if (verbose)
    {
      printf ("# Mode[%d] = %.15Lg calculated, took %fs, total time %fs\n", i, qgmode->k, g_timer_elapsed (bench, NULL) - part_time, g_timer_elapsed (bench, NULL));
      part_time = g_timer_elapsed (bench, NULL);
    }
    else
      printf (".");
    fflush (stdout);
  }
  if (!verbose)
    printf ("\n");
//printf ("\n\n");
  ncm_spline_prepare (qgmode->pw_spline);
  return TRUE;
}

gdouble
nc_hicosmo_qg_pert_powerspectrum (NcHICosmoQGMode *qgmode, gdouble x, gdouble *R)
{
  NcmModel *model = qgmode->model;
  const gdouble Rbar = hypot (R[1], R[3]);
  const gdouble Rp = Rbar * x;
  const gdouble zeta = nc_hicosmo_qg_zeta (model, x, NULL);
  
  const gdouble gbar = sqrt(nc_hicosmo_qg_gbar2 (model, x));
  const gdouble delta_phi = gsl_pow_3 (x) * gbar * zeta * zeta * Rp / (sqrt (qgmode->k * M_PI) * NC_C_HUBBLE_RADIUS_PLANCK);

  if (FALSE)
  {
    const gdouble cs2 = nc_hicosmo_qg_cs2 (model, x, NULL);
    const gdouble cs = sqrt (cs2);
    const gdouble P1 = zeta * zeta * x * R[1];
    const gdouble P2 = zeta * zeta * x * R[3];
    const gdouble rho1 = hypot (P1 / zeta, R[0] * cs * qgmode->k * zeta);
    const gdouble rho2 = hypot (P2 / zeta, R[2] * cs * qgmode->k * zeta);
    const gdouble theta1 = atan (cs * qgmode->k * zeta * zeta * R[0] / P1);
    const gdouble theta2 = atan (cs * qgmode->k * zeta * zeta * R[2] / P2);
    const gdouble J1 = zeta * zeta * (x * x * R[1] * R[1] / (qgmode->k * cs) + qgmode->k * cs * R[0] * R[0]) / 2.0;
    const gdouble J2 = zeta * zeta * (x * x * R[3] * R[3] / (qgmode->k * cs) + qgmode->k * cs * R[2] * R[2]) / 2.0;
    
    printf ("% 12.8e % 20.15g %20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g ", x, J1, J2, theta1, rho1, theta2, rho2, R[0], zeta);
  }
  
  return delta_phi;
}

gboolean
nc_hicosmo_qg_modefunc_cvode_init (NcHICosmoQGMode *qgmode)
{
  gint flag;
  CVRhsFn f;
  if (OLD_CODE)
    f = &hk_step;
  else
    f = &Rak_step;
  
  if (!*qgmode->initialized)
  {
    flag = CVodeInit(qgmode->cvode, f, qgmode->alphai, qgmode->y);
    CVODE_CHECK(&flag, "CVodeMalloc", 1, FALSE);
    *qgmode->initialized = TRUE;
  }
  else
  {
    flag = CVodeReInit (qgmode->cvode, qgmode->alphai, qgmode->y);
    CVODE_CHECK(&flag, "CVodeReInit", 1, FALSE);
  }

  nc_hicosmo_qg_modefunc_set_opts (qgmode);
  return TRUE;
}

gboolean
nc_hicosmo_qg_modefunc_init (NcHICosmoQGMode *qgmode)
{
  NcmModel *model = qgmode->model;
  long double h[4];

  nc_hicosmo_qg_h_dust_sol (qgmode, exp (qgmode->alpha0), 1e-60L, h);
  NV_Ith_S(qgmode->y, 0) = h[0];
  NV_Ith_S(qgmode->y, 1) = h[1];
  NV_Ith_S(qgmode->y, 2) = h[2];
  NV_Ith_S(qgmode->y, 3) = h[3];

  if (!OLD_CODE)
  {
    gdouble alpha = -sqrt(-2.0 * (qgmode->alpha0 + G_B));
    gdouble alphap = sqrt(nc_hicosmo_qg_alphaprime2 (model, alpha, NULL));
    nc_hicosmo_qg_pert_h_to_R (qgmode, exp (qgmode->alpha0), N_VGetArrayPointer (qgmode->y));
    NV_Ith_S(qgmode->y, 1) *= exp (qgmode->alpha0) / alphap;
    NV_Ith_S(qgmode->y, 3) *= exp (qgmode->alpha0) / alphap;
  }
  
  printf ("# % .15g % .15Lg [% .15e % .15e % .15e % .15e] % .15g\n", exp (qgmode->alpha0), qgmode->k, NV_Ith_S(qgmode->y, 0), NV_Ith_S(qgmode->y, 1), NV_Ith_S(qgmode->y, 2), NV_Ith_S(qgmode->y, 3), hypot (NV_Ith_S(qgmode->y, 0), NV_Ith_S(qgmode->y, 2)) / exp (qgmode->alpha0));
  {
    gdouble alpha = -sqrt(-2.0 * (qgmode->alpha0 + G_B));
    gdouble R[4];
    gdouble x = exp (qgmode->alpha0);
    gint i;
    nc_hicosmo_qg_R_dust_sol (qgmode, x, R);
    printf ("# % .15g % .15Lg [% .15e % .15e % .15e % .15e] % .15g\n", x, qgmode->k,  R[0], R[1], R[2], R[3], hypot (R[0], R[2]) / x);
    R[1] *= - alpha * x;
    R[3] *= - alpha * x;
    for (i = 0; i < 4; i++)
      NV_Ith_S(qgmode->y, i) = R[i];
  }
  
  qgmode->reltol = 1e-13;
  qgmode->abstol = 0.0;

  return TRUE;
}

gboolean
nc_hicosmo_qg_modefunc_set_opts (NcHICosmoQGMode *qgmode)
{
  gint flag;

  flag = CVodeSStolerances (qgmode->cvode, qgmode->reltol, qgmode->abstol);
  CVODE_CHECK(&flag, "CVodeSStolerances", 1, FALSE);
  
  flag = CVodeSetMaxNumSteps(qgmode->cvode, 10000000);
  CVODE_CHECK(&flag, "CVodeSetMaxNumSteps", 1, FALSE);
  
  flag = CVodeSetUserData (qgmode->cvode, qgmode);
  CVODE_CHECK(&flag, "CVodeSetUserData", 1, FALSE);

  if (FALSE)
  {
    flag = CVDense(qgmode->cvode, 4);
    CVODE_CHECK(&flag, "CVDense", 1, FALSE);

    if (OLD_CODE)
    {
      flag = CVDlsSetDenseJacFn (qgmode->cvode, &hk_step_J);
      CVODE_CHECK(&flag, "CVDenseSetJacFn", 1, FALSE);
    }
    else
    {
      flag = CVDlsSetDenseJacFn (qgmode->cvode, &Rak_step_J);
      CVODE_CHECK(&flag, "CVDenseSetJacFn", 1, FALSE);
    }
  }

//  flag = CVodeSetStabLimDet(qgmode->cvode, TRUE);
//  CVODE_CHECK(&flag, "CVodeSetStabLimDet", 1, FALSE);
  
//  flag = CVodeSetMaxStep(qgmode->cvode, 10.0); 
//  CVODE_CHECK(&flag, "CVodeSetMaxStep", 1, FALSE);        
  
//  flag = CVodeRootInit (qgmode->cvode, 2, scale_factor_root, params);
//  CVODE_CHECK(&flag, "CVodeRootInit", 1, FALSE);        
  
//  flag = CVodeSetStopTime(qgmode->cvode, qgmode->tf);
//  CVODE_CHECK(&flag, "CVodeSetStopTime", 1, FALSE);
  
  return TRUE;
}

long double
mi_step_f (long double alpha, gpointer params)
{
  gdouble x = expl (alpha);
  NcHICosmoQGMode *qgmode = (NcHICosmoQGMode *)params;
  NcmModel *model = qgmode->model;
  long double x2d2sqrtxxbarzeta_sqrtxxbarzeta;// = nc_hicosmo_qg_d2sqrtxxbarzeta_sqrtxxbarzeta (data, x);
  long double x2cs2_xxbar2;// = nc_hicosmo_qg_cs2_xxbar2 (data, x);
  long double k2 = qgmode->k * qgmode->k;

  nc_hicosmo_qg_evolfunc (model, x, &x2d2sqrtxxbarzeta_sqrtxxbarzeta, &x2cs2_xxbar2);

  return (k2 * x2cs2_xxbar2 - x2d2sqrtxxbarzeta_sqrtxxbarzeta - 1.0 / 4.0);
}

gboolean
nc_hicosmo_qg_modefunc_evolve (NcHICosmoQGMode *qgmode)
{
  long double h[4], h1[4], last_x = -10000000.0;
  NcmMIOdeFunction F;
  NcmMIOde *mi_ode;
  gdouble tmp[4];
  gint flag, i;
  
  if (FALSE)
  {
    F.func = &mi_step_f;
    F.params = qgmode;
    mi_ode = ncm_magnus_iserles_ode_new (qgmode->alpha0, &F); /* FIXME LEAK!!! */

    do 
    {
      ncm_magnus_iserles_ode_step (mi_ode, 1e-5L);
      if (TRUE || (mi_ode->xi < 0 ? (mi_ode->xi > (last_x * 0.999)) : ((mi_ode->xi > (last_x * 1.001)))))
      {
        last_x = mi_ode->xi;
        if (TRUE)
        {
          nc_hicosmo_qg_h_dust_sol (qgmode, expl (mi_ode->xi), 1e-20L, h);
          ncm_magnus_iserles_ode_eval_vec (mi_ode, NV_Ith_S (qgmode->y, 0), NV_Ith_S (qgmode->y, 1), &h1[0], &h1[1]);
          ncm_magnus_iserles_ode_eval_vec (mi_ode, NV_Ith_S (qgmode->y, 2), NV_Ith_S (qgmode->y, 3), &h1[2], &h1[3]);
          for (i = 0; i < 4; i++) tmp[i] = h1[i];
          printf ("% .6Le % .6Le % .6Le % .6Le % .6Le % .6Le % .6Le % .6e %.15Le %.15Le %.15e %.15Le\n",
                  expl (mi_ode->xi), mi_ode->xi, mi_ode->psi,
                  h[0], h[1], h1[0], h1[1],
                  fabs(1.0 - h1[0] / h[0]), hypotl (h[0], h[2]), hypotl (h[1], h[3]),
                  nc_hicosmo_qg_pert_powerspectrum (qgmode, mi_ode->xi, tmp),
                  qgmode->k);
        }
      }
    } while (mi_ode->psi > 1.0e-1L && FALSE);

    qgmode->alphai = mi_ode->xi;
    ncm_magnus_iserles_ode_eval_vec (mi_ode, NV_Ith_S (qgmode->y, 0), NV_Ith_S (qgmode->y, 1), &h1[0], &h1[1]);
    ncm_magnus_iserles_ode_eval_vec (mi_ode, NV_Ith_S (qgmode->y, 2), NV_Ith_S (qgmode->y, 3), &h1[2], &h1[3]);

    NV_Ith_S (qgmode->y, 0) = h1[0];
    NV_Ith_S (qgmode->y, 1) = h1[1];
    NV_Ith_S (qgmode->y, 2) = h1[2];
    NV_Ith_S (qgmode->y, 3) = h1[3];
  }
  else
    qgmode->alphai = qgmode->alpha0;

  if (!OLD_CODE)
  {
    NcmModel *model = qgmode->model;
    qgmode->alphai = -sqrt(-2.0 * (qgmode->alphai + G_B));
    qgmode->alphaf = 6.0;//sqrt(-2.0 * (-15.0 * M_LN10 + G_B));
  }

  //qgmode->alphaf = 0.0;
  nc_hicosmo_qg_modefunc_cvode_init (qgmode);
  flag = CVodeSetStopTime(qgmode->cvode, qgmode->alphaf);
  CVODE_CHECK(&flag, "CVodeSetStopTime", 1, FALSE);
  printf ("# Ode solver changed at %.15g.\n", exp(qgmode->alphai));
  while (1)
  {
    gint flag;
    gdouble x_i;
    
    flag = CVode (qgmode->cvode, qgmode->alphaf, qgmode->y, &x_i, CV_ONE_STEP);
    CVODE_CHECK(&flag, "CVode[nc_hicosmo_qg_modefunc_evolve]", 1, FALSE);

    if (FALSE && (x_i < 0 ? (x_i > (last_x * 0.9999)) : ((x_i > (last_x * 1.0001)))))
    {
      last_x = x_i;
      if (TRUE)
      {
        NcmModel *model = qgmode->model;
        gdouble x = expl (x_i);
        gdouble F_x = nc_hicosmo_qg_xxbarzeta2 (model, x, NULL) / x;
        gdouble sqrt_x_F = sqrt (1.0 / F_x);
        gdouble R[4];
        gdouble RFo = hypot (NV_Ith_S (qgmode->y, 0), NV_Ith_S (qgmode->y, 2)) / NC_C_HUBBLE_RADIUS_PLANCK;
        //printf ("% .15g % .15g % .15g\n", exp (-x_i * x_i / 2.0 - G_B), x_i, gsl_pow_3 (qgmode->k) * gsl_pow_2(hypot (NV_Ith_S (qgmode->y, 0), NV_Ith_S (qgmode->y, 2))));

        nc_hicosmo_qg_h_dust_sol (qgmode, exp (x_i), 1e-20, h);
        if (OLD_CODE)
        {
          nc_hicosmo_qg_pert_h_to_R (qgmode, x, N_VGetArrayPointer (qgmode->y));
        }
        else
        {
          gdouble alphap = sqrt(nc_hicosmo_qg_alphaprime2 (model, x_i, NULL));
          gdouble hx = exp (-x_i * x_i / 2.0 - G_B);
          NV_Ith_S (qgmode->y, 1) *= alphap / hx;
          NV_Ith_S (qgmode->y, 3) *= alphap / hx;
        }
        nc_hicosmo_qg_R_dust_sol (qgmode, exp (-x_i * x_i / 2.0 - G_B), R);
        
        printf ("% .6e % .6e % .15e % .6e % .6e % .6e % .6e % .6e %.15e %.15e %.15e %.15Le\n", 
                exp (-x_i * x_i / 2.0 - G_B), x_i, RFo,
                R[0], -x_i * exp (-x_i * x_i / 2.0 - G_B) * R[1], NV_Ith_S (qgmode->y, 0), NV_Ith_S (qgmode->y, 1), 
                NV_Ith_S (qgmode->y, 2), NV_Ith_S (qgmode->y, 3), hypot (NV_Ith_S (qgmode->y, 0), NV_Ith_S (qgmode->y, 2)) * sqrt_x_F,
                nc_hicosmo_qg_pert_powerspectrum (qgmode, exp (-x_i * x_i / 2.0 - G_B), N_VGetArrayPointer (qgmode->y)),
                qgmode->k);
                
      }
    }

    if (x_i >= qgmode->alphaf)
    {
      qgmode->alphai = x_i;
      break;
    }
    
    if (FALSE)
    {
      glong nsteps, nfevals, 	nlinsetups, netfails;
      gint qlast, qcur;
 		 	gdouble hinused, hlast, hcur, tcur;

      CVodeGetIntegratorStats (qgmode->cvode, &nsteps, &nfevals, 	 
 		 	&nlinsetups, &netfails, &qlast, &qcur, 	 
 		 	&hinused, &hlast, &hcur, &tcur);

      printf ("# %ld %ld %ld %ld - %d %d - %.15e %.15e %.15e %.15e\n", nsteps, nfevals, 	nlinsetups, netfails, qlast, qcur, hinused, hlast, hcur, tcur);
    }
  }
  //printf ("\n\n");
  if (TRUE)
  {
    //printf ("% .15Lg %.15g\n", qgmode->k, hypot (NV_Ith_S (qgmode->y, 0), NV_Ith_S (qgmode->y, 2)));
    nc_hicosmo_qg_h_dust_sol (qgmode, exp (qgmode->alphai), 1e-20, h);
    if (OLD_CODE)
    {
      nc_hicosmo_qg_pert_h_to_R (qgmode, exp (qgmode->alphaf), N_VGetArrayPointer (qgmode->y));
    }
    else
    {
      NcmModel *model = qgmode->model;
      gdouble alphap = sqrt(nc_hicosmo_qg_alphaprime2 (model, qgmode->alphaf, NULL));
      gdouble hx = exp (-qgmode->alphaf * qgmode->alphaf / 2.0 - G_B);
      NV_Ith_S (qgmode->y, 1) *= alphap / hx;
      NV_Ith_S (qgmode->y, 3) *= alphap / hx;
    }
    
    printf ("%Lg %g %g %g %g %g %g %Lg %.15g\n", qgmode->k,
      NV_Ith_S (qgmode->y, 0), NV_Ith_S (qgmode->y, 1), 
      NV_Ith_S (qgmode->y, 2), NV_Ith_S (qgmode->y, 3), 
      hypot (NV_Ith_S (qgmode->y, 0), NV_Ith_S (qgmode->y, 2)), 
      hypot (NV_Ith_S (qgmode->y, 1), NV_Ith_S (qgmode->y, 3)),
      hypotl (h[0], h[2]),
      nc_hicosmo_qg_pert_powerspectrum (qgmode, qgmode->alphaf, N_VGetArrayPointer (qgmode->y))
      );
  }
  return TRUE;
}

static gint 
hk_step (realtype alpha, N_Vector y, N_Vector ydot, gpointer f_data)
{
  gdouble x = exp (alpha);
  NcHICosmoQGMode *qgmode = (NcHICosmoQGMode *)f_data;
  NcmModel *model = qgmode->model;
  long double x2d2sqrtxxbarzeta_sqrtxxbarzeta;// = nc_hicosmo_qg_d2sqrtxxbarzeta_sqrtxxbarzeta (data, x);
  long double x2cs2_xxbar2;// = nc_hicosmo_qg_cs2_xxbar2 (data, x);
  gdouble k2 = gsl_pow_2 (qgmode->k);

  nc_hicosmo_qg_evolfunc (model, x, &x2d2sqrtxxbarzeta_sqrtxxbarzeta, &x2cs2_xxbar2);
  
  NV_Ith_S(ydot, 0) = NV_Ith_S(y, 1);
  NV_Ith_S(ydot, 1) = (x2d2sqrtxxbarzeta_sqrtxxbarzeta + 1.0 / 4.0 - k2 * x2cs2_xxbar2) * NV_Ith_S(y, 0);
  
  NV_Ith_S(ydot, 2) = NV_Ith_S(y, 3);
  NV_Ith_S(ydot, 3) = (x2d2sqrtxxbarzeta_sqrtxxbarzeta + 1.0 / 4.0 - k2 * x2cs2_xxbar2) * NV_Ith_S(y, 2);
  
//  fprintf (stdout, "y [%.15e] %.15e %.15e %.15e %.15e\n", x, NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2), NV_Ith_S(y, 3));
//  fprintf (stdout, "dy[%.15e] %.15e %.15e %.15e %.15e\n", x, NV_Ith_S(ydot, 0), NV_Ith_S(ydot, 1), NV_Ith_S(ydot, 2), NV_Ith_S(ydot, 3));
  return 0;
}

static gint 
hk_step_J (gint N, realtype alpha, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  gdouble x = exp (alpha);
  NcHICosmoQGMode *qgmode = (NcHICosmoQGMode *)jac_data;
  NcmModel *model = qgmode->model;
  long double x2d2sqrtxxbarzeta_sqrtxxbarzeta;// = nc_hicosmo_qg_d2sqrtxxbarzeta_sqrtxxbarzeta (data, x);
  long double x2cs2_xxbar2;// = nc_hicosmo_qg_cs2_xxbar2 (data, x);
  gdouble k2 = gsl_pow_2 (qgmode->k);

  nc_hicosmo_qg_evolfunc (model, x, &x2d2sqrtxxbarzeta_sqrtxxbarzeta, &x2cs2_xxbar2);

  DENSE_ELEM (J, 0, 0) = 0.0;
  DENSE_ELEM (J, 0, 1) = 1.0;
  DENSE_ELEM (J, 0, 2) = 0.0;
  DENSE_ELEM (J, 0, 3) = 0.0;

  DENSE_ELEM (J, 1, 0) = (x2d2sqrtxxbarzeta_sqrtxxbarzeta + 1.0 / 4.0 - k2 * x2cs2_xxbar2);
  DENSE_ELEM (J, 1, 1) = 0.0;
  DENSE_ELEM (J, 1, 2) = 0.0;
  DENSE_ELEM (J, 1, 3) = 0.0;

  DENSE_ELEM (J, 2, 0) = 0.0;
  DENSE_ELEM (J, 2, 1) = 0.0;
  DENSE_ELEM (J, 2, 2) = 0.0;
  DENSE_ELEM (J, 2, 3) = 1.0;

  DENSE_ELEM (J, 3, 0) = 0.0;
  DENSE_ELEM (J, 3, 1) = 0.0;
  DENSE_ELEM (J, 3, 2) = (x2d2sqrtxxbarzeta_sqrtxxbarzeta + 1.0 / 4.0 - k2 * x2cs2_xxbar2);
  DENSE_ELEM (J, 3, 3) = 0.0;

  return 0;
}

static gint 
Rk_step (realtype lambda, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHICosmoQGMode *qgmode = (NcHICosmoQGMode *)f_data;
  gdouble x = nc_hicosmo_qg_x_lambda (qgmode->model, lambda, FALSE);
  gdouble gbar = nc_hicosmo_qg_gbar_lambda (qgmode->model, lambda, FALSE);
  gdouble x2 = x * x;
  gdouble cs2 = nc_hicosmo_qg_cs2 (qgmode->model, x, NULL);
//  gdouble dcs2 = nc_hicosmo_qg_dcs2 (qgmode->model, x, NULL);
  gdouble k2 = qgmode->k * qgmode->k;
  gdouble dzeta_zeta = nc_hicosmo_qg_dzeta_zeta (qgmode->model, x, NULL);
  gdouble cs2k2 = cs2 * k2;
  gdouble cs2k2_x2 =  cs2k2 / x2;
  gdouble x2dzeta_zeta_p1gbar = (2.0 * x * dzeta_zeta + 1.0) * gbar;

  NV_Ith_S(ydot, 0) = NV_Ith_S(y, 1);
  NV_Ith_S(ydot, 1) = x2dzeta_zeta_p1gbar * NV_Ith_S(y, 1)  - cs2k2_x2 * NV_Ith_S(y, 0);
  NV_Ith_S(ydot, 2) = NV_Ith_S (y, 3);
  NV_Ith_S(ydot, 3) = x2dzeta_zeta_p1gbar * NV_Ith_S(y, 3)  - cs2k2_x2 * NV_Ith_S(y, 2);

//  fprintf (stdout, "y [%.15e] %.15e %.15e %.15e %.15e\n", x, NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2), NV_Ith_S(y, 3));
//  fprintf (stdout, "dy[%.15e] %.15e %.15e %.15e %.15e\n", x, NV_Ith_S(ydot, 0), NV_Ith_S(ydot, 1), NV_Ith_S(ydot, 2), NV_Ith_S(ydot, 3));
 
  return 0;
}

static gint 
Rk_step_J (gint N, realtype lambda, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHICosmoQGMode *qgmode = (NcHICosmoQGMode *)jac_data;
  gdouble x = nc_hicosmo_qg_x_lambda (qgmode->model, lambda, FALSE);
  gdouble gbar = nc_hicosmo_qg_gbar_lambda (qgmode->model, lambda, FALSE);
  gdouble x2 = x*x;
  gdouble cs2 = nc_hicosmo_qg_cs2 (qgmode->model, x, NULL);
  gdouble k2 = qgmode->k * qgmode->k;
  gdouble dzeta_zeta = nc_hicosmo_qg_dzeta_zeta (qgmode->model, x, NULL);
  gdouble cs2k2 = cs2 * k2;
  gdouble cs2k2_x2 =  cs2k2 / x2;
  gdouble x2dzeta_zeta_p1gbar = (2.0 * x * dzeta_zeta + 1.0) * gbar;

  DENSE_ELEM (J, 0, 0) = 0.0;
  DENSE_ELEM (J, 0, 1) = 1.0;
  DENSE_ELEM (J, 0, 2) = 0.0;
  DENSE_ELEM (J, 0, 3) = 0.0;

  DENSE_ELEM (J, 1, 0) = -cs2k2_x2;
  DENSE_ELEM (J, 1, 1) = x2dzeta_zeta_p1gbar;
  DENSE_ELEM (J, 1, 2) = 0.0;
  DENSE_ELEM (J, 1, 3) = 0.0;

  DENSE_ELEM (J, 2, 0) = 0.0;
  DENSE_ELEM (J, 2, 1) = 0.0;
  DENSE_ELEM (J, 2, 2) = 0.0;
  DENSE_ELEM (J, 2, 3) = 1.0;

  DENSE_ELEM (J, 3, 0) = 0.0;
  DENSE_ELEM (J, 3, 1) = 0.0;
  DENSE_ELEM (J, 3, 2) = -cs2k2_x2;
  DENSE_ELEM (J, 3, 3) = x2dzeta_zeta_p1gbar;

  return 0;
}

static gint 
Rak_step (realtype alpha, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHICosmoQGMode *qgmode = (NcHICosmoQGMode *)f_data;
  NcmModel *model = qgmode->model;
  gdouble x = exp (-alpha * alpha / 2.0 - G_B);
  gdouble dalphap2_dalpha = nc_hicosmo_qg_dalphaprime2_dalpha (model, alpha, NULL);
  gdouble alphap2 = nc_hicosmo_qg_alphaprime2 (model, alpha, NULL);
  gdouble dzeta2_zeta2 = 2.0 * nc_hicosmo_qg_dzeta_zeta (model, x, NULL);
  gdouble cs2 = nc_hicosmo_qg_cs2 (model, x, NULL);
  gdouble k2 = qgmode->k * qgmode->k;  
  
  NV_Ith_S(ydot, 0) = NV_Ith_S(y, 1);
  NV_Ith_S(ydot, 1) = 
    (alpha * x * dzeta2_zeta2 - dalphap2_dalpha / (2.0 * alphap2)) * NV_Ith_S(y, 1) +
    -cs2 * k2 / alphap2 * NV_Ith_S(y, 0);  
  NV_Ith_S(ydot, 2) = NV_Ith_S(y, 3);
  NV_Ith_S(ydot, 3) = 
    (alpha * x * dzeta2_zeta2 - dalphap2_dalpha / (2.0 * alphap2)) * NV_Ith_S(y, 3) +
    -cs2 * k2 / alphap2 * NV_Ith_S(y, 2);
  
//  fprintf (stdout, "y [%.15e] %.15e %.15e %.15e %.15e\n", x, NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2), NV_Ith_S(y, 3));
//  fprintf (stdout, "dy[%.15e] %.15e %.15e %.15e %.15e\n", x, NV_Ith_S(ydot, 0), NV_Ith_S(ydot, 1), NV_Ith_S(ydot, 2), NV_Ith_S(ydot, 3));
  return 0;
}

static gint 
Rak_step_J (gint N, realtype alpha, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHICosmoQGMode *qgmode = (NcHICosmoQGMode *)jac_data;
  NcmModel *model = qgmode->model;
  gdouble x = exp (-alpha * alpha / 2.0 - G_B);
  gdouble dalphap2_dalpha = nc_hicosmo_qg_dalphaprime2_dalpha (model, alpha, NULL);
  gdouble alphap2 = nc_hicosmo_qg_alphaprime2 (model, alpha, NULL);
  gdouble dzeta2_zeta2 = 2.0 * nc_hicosmo_qg_dzeta_zeta (model, x, NULL);
  gdouble cs2 = nc_hicosmo_qg_cs2 (model, x, NULL);
  gdouble k2 = qgmode->k * qgmode->k;  

  DENSE_ELEM (J, 0, 0) = 0.0;
  DENSE_ELEM (J, 0, 1) = 1.0;
  DENSE_ELEM (J, 0, 2) = 0.0;
  DENSE_ELEM (J, 0, 3) = 0.0;

  DENSE_ELEM (J, 1, 0) = -cs2 * k2 / alphap2;
  DENSE_ELEM (J, 1, 1) = (alpha * x * dzeta2_zeta2 - dalphap2_dalpha / (2.0 * alphap2));
  DENSE_ELEM (J, 1, 2) = 0.0;
  DENSE_ELEM (J, 1, 3) = 0.0;

  DENSE_ELEM (J, 2, 0) = 0.0;
  DENSE_ELEM (J, 2, 1) = 0.0;
  DENSE_ELEM (J, 2, 2) = 0.0;
  DENSE_ELEM (J, 2, 3) = 1.0;

  DENSE_ELEM (J, 3, 0) = 0.0;
  DENSE_ELEM (J, 3, 1) = 0.0;
  DENSE_ELEM (J, 3, 2) = -cs2 * k2 / alphap2;
  DENSE_ELEM (J, 3, 3) = (alpha * x * dzeta2_zeta2 - dalphap2_dalpha / (2.0 * alphap2));

  return 0;
}
