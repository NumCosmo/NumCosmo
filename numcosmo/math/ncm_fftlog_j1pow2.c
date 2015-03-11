/***************************************************************************
 *            ncm_fftlog_j1pow2.c
 *
 *  Mon July 21 19:59:38 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_fftlog_j1pow2.c
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_fftlog_j1pow2
 * @title: NcmFftlogJ1pow2
 * @short_description: Logarithm fast fourier transform for a kernel with a spherical bessel of order one squared.
 *
 * Kernel $(j_1(x)/x)^2$.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fftlog_j1pow2.h"
#include "math/ncm_cfg.h"

#include <math.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_trig.h>

G_DEFINE_TYPE (NcmFftlogJ1pow2, ncm_fftlog_j1pow2, NCM_TYPE_FFTLOG);

static void
ncm_fftlog_j1pow2_init (NcmFftlogJ1pow2 *j1pow2)
{
}

static void
ncm_fftlog_j1pow2_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fftlog_j1pow2_parent_class)->finalize (object);
}

void _ncm_fftlog_j1pow2_get_Ym (NcmFftlog *fftlog);
void _ncm_fftlog_j1pow2_generate_Gr (NcmFftlog *fftlog);

static void
ncm_fftlog_j1pow2_class_init (NcmFftlogJ1pow2Class *klass)
{
  GObjectClass* object_class   = G_OBJECT_CLASS (klass);
  NcmFftlogClass *fftlog_class = NCM_FFTLOG_CLASS (klass);

  object_class->finalize = ncm_fftlog_j1pow2_finalize;

  fftlog_class->name        = "sphericalbessel1_x_pow_2";
  fftlog_class->ncomp       = 2;
  fftlog_class->get_Ym      = &_ncm_fftlog_j1pow2_get_Ym;
  fftlog_class->generate_Gr = &_ncm_fftlog_j1pow2_generate_Gr;
}

void 
_ncm_fftlog_j1pow2_get_Ym (NcmFftlog *fftlog)
{
  /*NcmFftlogJ1pow2 *j1pow2 = NCM_FFTLOG_J1POW2 (fftlog);*/
#ifdef NUMCOSMO_HAVE_FFTW3  
  gint i;
  for (i = -fftlog->N_2; i <= fftlog->N_2; i++)
  {
    const gint ii = (i < 0) ? i + fftlog->N : i;
    const gdouble a = 2.0 * M_PI / fftlog->Lk * i;
    const gdouble abs_a = fabs (a);
    const gdouble sign_a = a < 0 ? -1.0 : 1.0;
    complex double U = 36.0 * (a + I) / (cexp (M_LN2 * I * a) * (a * I - 5.0));

    gsl_sf_result lnr;
    gsl_sf_result arg;

    gsl_sf_lngamma_complex_e (-3.0, a, &lnr, &arg);
    U = (a != 0) ? U * sign_a * cexp (lnr.val + arg.val * I + gsl_sf_lnsinh (M_PI * abs_a * 0.5)) : 3.0 * M_PI / 5.0;

    fftlog->Ym[0][ii] = U * cexp (-2.0 * M_PI / fftlog->Lk * i * I * (fftlog->lnk0 + fftlog->lnr0));
    fftlog->Ym[1][ii] = -(1.0 + I * a) * fftlog->Ym[0][ii];
  }
#endif /* NUMCOSMO_HAVE_FFTW3 */
}

void 
_ncm_fftlog_j1pow2_generate_Gr (NcmFftlog *fftlog)
{
#ifdef NUMCOSMO_HAVE_FFTW3
  gint i, j = 0;
  
  for (i = -fftlog->N_2; i <= fftlog->N_2; i++)
  {
    const gint ii            = (i < 0) ? i + fftlog->N : i;
    const gdouble lnr        = fftlog->lnr0 + i * fftlog->Lk_N;
    
    const gdouble rGr        = creal (fftlog->Gr[0][ii]) / fftlog->N;
    const gdouble lnGr       = log (rGr) - lnr;
    const gdouble dlnGr_dlnr = creal (fftlog->Gr[1][ii]) / (rGr * fftlog->N);

/*printf ("%d % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", ii, lnr / M_LN10, lnGr, dlnGr_dlnr, creal (fftlog->Gr[0][ii]), creal (fftlog->Gr[1][ii]));*/
    ncm_vector_set (fftlog->lnr_vec, j, lnr);
    ncm_vector_set (fftlog->Gr_vec[0], j, lnGr);
    ncm_vector_set (fftlog->Gr_vec[1], j, dlnGr_dlnr);
    j++;
  }
#endif /* NUMCOSMO_HAVE_FFTW3 */
}

NcmFftlog *
ncm_fftlog_j1pow2_new (gdouble r0, gdouble k0, gdouble Lk, guint N)
{
  NcmFftlog *fftlog = g_object_new (NCM_TYPE_FFTLOG_J1POW2, 
                                    "r0", r0,
                                    "k0", k0,
                                    "Lk", Lk,
                                    "N", N,
                                    NULL);
  return fftlog;
}


