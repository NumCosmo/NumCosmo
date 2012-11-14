/***************************************************************************
 *            test_nc_recomb.c
 *
 *  Wed November 14 11:38:56 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

void test_nc_recomb_seager_new (void);
void test_nc_recomb_seager_wmap_zstar (void);
void test_nc_recomb_seager_Xe_ini (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init ();
  /* ncm_cfg_enable_gsl_err_handler (); */

  g_test_add_func ("/numcosmo/nc_recomb/seager/new", &test_nc_recomb_seager_new);
  g_test_add_func ("/numcosmo/nc_recomb/seager/wmap/zstar", &test_nc_recomb_seager_wmap_zstar);
  g_test_add_func ("/numcosmo/nc_recomb/seager/wmap/Xe_ini", &test_nc_recomb_seager_Xe_ini);

  g_test_run ();
}

void
test_nc_recomb_seager_new ()
{
  NcRecomb *recomb = nc_recomb_seager_new ();
  g_assert (NC_IS_RECOMB (recomb));
  g_assert (NC_IS_RECOMB_SEAGER (recomb));
  nc_recomb_free (recomb);

  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    nc_recomb_free (recomb);
    exit (0);
  }
  g_test_trap_assert_failed ();

  recomb = nc_recomb_seager_new_full (1e-10, 2.2e9, 1e-5);
  g_assert (NC_IS_RECOMB (recomb));
  g_assert (NC_IS_RECOMB_SEAGER (recomb));

  g_assert (recomb->init_frac == 1e-10);
  g_assert (recomb->zi == 2.2e9);
  g_assert (recomb->prec == 1e-5);
  
  nc_recomb_free (recomb);

  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    nc_recomb_free (recomb);
    exit (0);
  }
  g_test_trap_assert_failed ();
  
  recomb = nc_recomb_new_from_name ("NcRecombSeager{'prec':<2e-7>}");
  g_assert (NC_IS_RECOMB (recomb));
  g_assert (NC_IS_RECOMB_SEAGER (recomb));
  nc_recomb_free (recomb);

  g_assert (recomb->prec == 2e-7);
  
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    nc_recomb_free (recomb);
    exit (0);
  }
  g_test_trap_assert_failed ();  
}

void
test_nc_recomb_seager_wmap_zstar (void)
{
  NcRecomb *recomb = nc_recomb_seager_new ();
  NcHICosmo *cosmo = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO, "NcHICosmoDEXcdm");
  nc_hicosmo_de_set_wmap5_params (cosmo);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_XCDM_W,  -1.0);

  nc_recomb_prepare_if_needed (recomb, cosmo);
  {
    const gdouble zstar = exp (-nc_recomb_tau_zstar (recomb, cosmo)) - 1.0;
    g_assert (ncm_cmp (zstar, 1088.9, 1e-4) == 0);
  }

  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_T_GAMMA0,  2.2250);
  nc_recomb_prepare_if_needed (recomb, cosmo);
  {
    const gdouble zstar = exp (-nc_recomb_tau_zstar (recomb, cosmo)) - 1.0;
    g_assert (ncm_cmp (zstar, 1323.2, 1e-4) == 0);
  }

  nc_recomb_free (recomb);
}

void 
test_nc_recomb_seager_Xe_ini (void)
{
  NcRecomb *recomb = nc_recomb_seager_new ();
  NcHICosmo *cosmo = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO, "NcHICosmoDEXcdm");
  nc_recomb_prepare_if_needed (recomb, cosmo);

  {
    const gdouble Xe_ini = nc_recomb_Xe (recomb, cosmo, recomb->lambdai);
    g_assert (ncm_cmp (Xe_ini, 1.0 + 2.0 * ncm_c_prim_XHe (), recomb->prec * 1.0e2) == 0);
  }
  
  nc_recomb_free (recomb);
}
