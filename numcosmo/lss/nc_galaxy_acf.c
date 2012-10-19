/***************************************************************************
 *            nc_galaxy_acf.c
 *
 *  Fri May 11 21:18:21 2012
 *  Copyright  2012 Fernando de Simoni & Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Fernando de Simoni & Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * SECTION:nc_galaxy_acf
 * @title: Galaxy Angular Corelation Function
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <gsl/gsl_sf_bessel.h>


G_DEFINE_TYPE (NcGalaxyAcf, nc_galaxy_acf, G_TYPE_OBJECT);

NcGalaxyAcf *
nc_galaxy_acf_new (NcGrowthFunc *gf, NcDistance *dist, NcTransferFunc *tf)
{
  NcGalaxyAcf *acf;
  acf = g_object_new (NC_TYPE_GALAXY_ACF, NULL);
  acf->gf = gf;      /* FIXME REF */
  acf->dist = dist;
  acf->tf = tf;
  acf->b = 2.0;
  acf->s = ncm_spline_cubic_notaknot_new ();
  {
	NcWindow *wp = nc_window_tophat_new ();
	acf->vp = nc_matter_var_new (NC_MATTER_VAR_FFT, wp, tf);
	nc_window_free (wp);
  }
  return acf;
}

typedef struct _NcmGalaxyAcfPsiKernel
{
  gdouble h;
  gdouble k;
  guint l;
  NcGalaxyAcf *acf;
  NcHICosmo *model;
  NcmSpline *jl;
} NcmGalaxyAcfPsiKernel;

static gdouble
ncm_galaxy_acf_psi_kernel_l (gdouble z, gpointer p)
{
  NcmGalaxyAcfPsiKernel *apk = (NcmGalaxyAcfPsiKernel *) p;
  const gdouble cd = NC_C_HUBBLE_RADIUS * nc_distance_comoving (apk->acf->dist, apk->model, z);
  const gdouble Dz = nc_growth_func_eval (apk->acf->gf, apk->model, z);
  const gdouble fz = nc_growth_func_eval_deriv (apk->acf->gf, apk->model, z);
  const gdouble x = apk->k * cd;
  const gdouble x2 = x * x;
  const gdouble sel_func = 1.0;
  const gdouble psi_kernel_l = Dz * sel_func * (apk->acf->b + fz * ((1.0 - (apk->l * (apk->l - 1)) / x2) ));
  return psi_kernel_l;
}

static gdouble
ncm_galaxy_acf_psi_kernel_lp1 (gdouble z, gpointer p)
{
  NcmGalaxyAcfPsiKernel *apk = (NcmGalaxyAcfPsiKernel *) p;
  const gdouble cd = NC_C_HUBBLE_RADIUS * nc_distance_comoving (apk->acf->dist, apk->model, z);
  const gdouble Dz = nc_growth_func_eval (apk->acf->gf, apk->model, z);
  const gdouble fz = nc_growth_func_eval_deriv (apk->acf->gf, apk->model, z);
  const gdouble x = apk->k * cd;
  const gdouble sel_func = 1.0;
  const gdouble psi_kernel_lp1 = -2.0 * Dz * sel_func * fz / x;
  return psi_kernel_lp1;
}

static gdouble
ncm_galaxy_acf_psi_kernel (gdouble z, gpointer p)
{
  NcmGalaxyAcfPsiKernel *apk = (NcmGalaxyAcfPsiKernel *) p;
  const gdouble cd = NC_C_HUBBLE_RADIUS * nc_distance_comoving (apk->acf->dist, apk->model, z);
  const gdouble x = apk->k * cd;
  const gdouble x2 = x * x;
  const gdouble sel_func = 1.0;
  const gdouble jl = gsl_sf_bessel_jl (apk->l, x);//ncm_sf_sbessel (apk->l, x);//
  const gdouble jlp1 = gsl_sf_bessel_jl (apk->l + 1, x); //ncm_sf_sbessel (apk->l + 1, x);//
  gdouble Dz, fz;
  nc_growth_func_eval_both (apk->acf->gf, apk->model, z, &Dz, &fz);

  fz *= -(1 + z) / Dz;

  return Dz * sel_func * (apk->acf->b * jl +
                          fz * (
                                ((apk->l * (1.0 - apk->l) + x2) * jl - 2 * x * jlp1) / x2
                                )
                          );
}

gdouble
ncm_galaxy_acf_psi (NcGalaxyAcf *acf, NcHICosmo *model, gdouble k, guint l)
{
  NcmGalaxyAcfPsiKernel psi_kernel;
  gsl_function F;
  F.function = &ncm_galaxy_acf_psi_kernel_l;
  F.function = &ncm_galaxy_acf_psi_kernel_lp1;
  F.function = &ncm_galaxy_acf_psi_kernel;
  F.params = &psi_kernel;

  psi_kernel.k = k;
  psi_kernel.l = l;
  psi_kernel.acf = acf;
  psi_kernel.model = model;
  //psi_kernel.jl = ncm_sf_sbessel_spline (l, 1e-5, xf, 1.0e-6);

  {
	gdouble res, err;
	nc_integral_locked_a_b (&F, 0.9, 1.0, 0.0, 1e-5, &res, &err);
	return res;
  }

}

static gdouble
ncm_galaxy_acf_psi_int (gdouble mu, gpointer p)
{
  static glong count = 0;
  NcmGalaxyAcfPsiKernel *apk = (NcmGalaxyAcfPsiKernel *) p;
  const gdouble k = exp (mu);
  const gdouble Pk = nc_transfer_func_matter_powerspectrum (apk->acf->tf, apk->model, k);
  const gdouble psi = ncm_galaxy_acf_psi (apk->acf, apk->model, k, apk->l);
printf ("% 20.15g % 20.15g\n", mu, psi * psi);
//printf ("% 20.15g % 20.15g\n", mu, k * k * k * Pk * psi * psi);
  //if ((count++) % 200 == 0)
    //printf(".");fflush (stdout);
  return k * k * k * Pk * psi * psi;
}

void
ncm_galaxy_acf_prepare_psi (NcGalaxyAcf *acf, NcHICosmo *model, guint l)
{
  NcmGalaxyAcfPsiKernel psi_kernel;
  gsl_function F;
  const gdouble sqrt_norma = nc_matter_var_sigma8_sqrtvar0 (acf->vp, model);
  const gdouble norma = sqrt_norma * sqrt_norma;
  F.function = &ncm_galaxy_acf_psi_int;
  F.params = &psi_kernel;

  psi_kernel.h = nc_hicosmo_h (model);
  psi_kernel.l = l;
  psi_kernel.acf = acf;
  psi_kernel.model = model;

  {
	gdouble res, err;
	printf("#");fflush (stdout);
	nc_integral_locked_a_b (&F,
	                        log (1.0e-10 / psi_kernel.h),
	                        log (1.0e3 / psi_kernel.h),
	                        0.0, 1e-6,
	                        &res, &err);
	//printf ("\n%u % 20.15g\n", l, norma * res);
	fflush (stdout);
	//exit (0);
	return res;
  }

/*
  ncm_spline_set_func (acf->s, NCM_SPLINE_FUNCTION_SPLINE, &F,
                       log (1.0e-12 * NC_C_HUBBLE_RADIUS),
                       log (1.0e1 * NC_C_HUBBLE_RADIUS),
                       0, 1e-6);
					   */
}

static void
nc_galaxy_acf_init (NcGalaxyAcf *nc_galaxy_acf)
{

}

static void
nc_galaxy_acf_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_acf_parent_class)->finalize (object);
}

static void
nc_galaxy_acf_class_init (NcGalaxyAcfClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->finalize = nc_galaxy_acf_finalize;
}

