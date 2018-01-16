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
 * @title: NcGalaxyAcf
 * @short_description: Galaxy angular correlation function.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_galaxy_acf.h"
#include "lss/nc_window_tophat.h"
#include "math/integral.h"
#include "math/ncm_spline_cubic_notaknot.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_bessel.h>
#endif /* NUMCOSMO_GIR_SCAN */

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
  NcHICosmo *cosmo;
  NcmSpline *jl;
} NcmGalaxyAcfPsiKernel;

static gdouble
nc_galaxy_acf_psi_kernel_l (gdouble z, gpointer p)
{
  NcmGalaxyAcfPsiKernel *apk = (NcmGalaxyAcfPsiKernel *) p;
  const gdouble cd = ncm_c_hubble_radius_hm1_Mpc () * nc_distance_comoving (apk->acf->dist, apk->cosmo, z);
  const gdouble Dz = nc_growth_func_eval (apk->acf->gf, apk->cosmo, z);
  const gdouble fz = nc_growth_func_eval_deriv (apk->acf->gf, apk->cosmo, z);
  const gdouble x = apk->k * cd;
  const gdouble x2 = x * x;
  const gdouble sel_func = 1.0;
  const gdouble psi_kernel_l = Dz * sel_func * (apk->acf->b + fz * ((1.0 - (apk->l * (apk->l - 1)) / x2) ));
  return psi_kernel_l;
}

static gdouble
nc_galaxy_acf_psi_kernel_lp1 (gdouble z, gpointer p)
{
  NcmGalaxyAcfPsiKernel *apk = (NcmGalaxyAcfPsiKernel *) p;
  const gdouble cd = ncm_c_hubble_radius_hm1_Mpc () * nc_distance_comoving (apk->acf->dist, apk->cosmo, z);
  const gdouble Dz = nc_growth_func_eval (apk->acf->gf, apk->cosmo, z);
  const gdouble fz = nc_growth_func_eval_deriv (apk->acf->gf, apk->cosmo, z);
  const gdouble x = apk->k * cd;
  const gdouble sel_func = 1.0;
  const gdouble psi_kernel_lp1 = -2.0 * Dz * sel_func * fz / x;
  return psi_kernel_lp1;
}

static gdouble
nc_galaxy_acf_psi_kernel (gdouble z, gpointer p)
{
  NcmGalaxyAcfPsiKernel *apk = (NcmGalaxyAcfPsiKernel *) p;
  const gdouble cd = ncm_c_hubble_radius_hm1_Mpc () * nc_distance_comoving (apk->acf->dist, apk->cosmo, z);
  const gdouble x = apk->k * cd;
  const gdouble x2 = x * x;
  const gdouble sel_func = 1.0;
  const gdouble jl = gsl_sf_bessel_jl (apk->l, x);//ncm_sf_sbessel (apk->l, x);//
  const gdouble jlp1 = gsl_sf_bessel_jl (apk->l + 1, x); //ncm_sf_sbessel (apk->l + 1, x);//
  gdouble Dz, fz;
  nc_growth_func_eval_both (apk->acf->gf, apk->cosmo, z, &Dz, &fz);

  fz *= -(1 + z) / Dz;

  return Dz * sel_func * (apk->acf->b * jl +
                          fz * (
                                ((apk->l * (1.0 - apk->l) + x2) * jl - 2 * x * jlp1) / x2
                                )
                          );
}

gdouble
nc_galaxy_acf_psi (NcGalaxyAcf *acf, NcHICosmo *cosmo, gdouble k, guint l)
{
  NcmGalaxyAcfPsiKernel psi_kernel;
  gsl_function F;
  F.function = &nc_galaxy_acf_psi_kernel_l;
  F.function = &nc_galaxy_acf_psi_kernel_lp1;
  F.function = &nc_galaxy_acf_psi_kernel;
  F.params = &psi_kernel;

  psi_kernel.k = k;
  psi_kernel.l = l;
  psi_kernel.acf = acf;
  psi_kernel.cosmo = cosmo;
  //psi_kernel.jl = ncm_sf_sbessel_spline (l, 1e-5, xf, 1.0e-6);

  {
	gdouble res, err;
	ncm_integral_locked_a_b (&F, 0.9, 1.0, 0.0, 1e-5, &res, &err);
	return res;
  }

}

static gdouble
nc_galaxy_acf_psi_int (gdouble mu, gpointer p)
{
//  static glong count = 0;
  NcmGalaxyAcfPsiKernel *apk = (NcmGalaxyAcfPsiKernel *) p;
  const gdouble k = exp (mu);
  const gdouble Pk = 1.0;/*nc_transfer_func_matter_powerspectrum (apk->acf->tf, apk->cosmo, k);*/
  const gdouble psi = nc_galaxy_acf_psi (apk->acf, apk->cosmo, k, apk->l);
printf ("% 20.15g % 20.15g\n", mu, psi * psi);
//printf ("% 20.15g % 20.15g\n", mu, k * k * k * Pk * psi * psi);
  //if ((count++) % 200 == 0)
    //printf(".");fflush (stdout);

  g_assert_not_reached ();
  
  return k * k * k * Pk * psi * psi;
}

void
nc_galaxy_acf_prepare_psi (NcGalaxyAcf *acf, NcHICosmo *cosmo, guint l)
{
  NcmGalaxyAcfPsiKernel psi_kernel;
  gsl_function F;
//  const gdouble sqrt_norma = nc_matter_var_sigma8_sqrtvar0 (acf->vp, cosmo);
//  const gdouble norma = sqrt_norma * sqrt_norma;
  F.function = &nc_galaxy_acf_psi_int;
  F.params = &psi_kernel;

  psi_kernel.h = nc_hicosmo_h (cosmo);
  psi_kernel.l = l;
  psi_kernel.acf = acf;
  psi_kernel.cosmo = cosmo;

  nc_transfer_func_prepare (acf->tf, cosmo);
  
  {
    gdouble res, err;
    printf("#");fflush (stdout);
    ncm_integral_locked_a_b (&F,
                             log (1.0e-10 / psi_kernel.h),
                             log (1.0e3 / psi_kernel.h),
                             0.0, 1e-6,
                             &res, &err);
    //printf ("\n%u % 20.15g\n", l, norma * res);
    fflush (stdout);
    //exit (0);
    return;
  }

/*
  ncm_spline_set_func (acf->s, NCM_SPLINE_FUNCTION_SPLINE, &F,
                       log (1.0e-12 * ncm_c_hubble_radius ()),
                       log (1.0e1 * ncm_c_hubble_radius ()),
                       0, 1e-6);
					   */
}

static void
nc_galaxy_acf_init (NcGalaxyAcf *nc_galaxy_acf)
{
  NCM_UNUSED (nc_galaxy_acf);
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

