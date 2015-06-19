/***************************************************************************
 *            test_nc_cluster_pseudo_counts.c
 *
 *  Fri June 5 11:45:16 2015
 *  Copyright  2015  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2015 <pennalima@gmail.com>
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

#include <math.h>
#include <glib.h>
#include <glib-object.h>

#define _TEST_NC_CLUSTER_PSEUDO_COUNTS_DATA_POINTS 21

typedef struct _TestNcClusterPseudoCounts
{
  NcClusterPseudoCounts *cpc;
  NcClusterMass *mszl;
  NcHICosmo *cosmo;
  NcMatterVar *vp;
  gdouble z;
  gdouble Mobs[2];
  gdouble Mobs_params[2];
} TestNcClusterPseudoCounts;

void test_nc_cluster_pseudo_counts_new (TestNcClusterPseudoCounts *test, gconstpointer pdata);
void test_nc_cluster_pseudo_counts_1p2_integral (TestNcClusterPseudoCounts *test, gconstpointer pdata);
void test_nc_cluster_pseudo_counts_3d_integral (TestNcClusterPseudoCounts *test, gconstpointer pdata);
void test_nc_cluster_pseudo_counts_free (TestNcClusterPseudoCounts *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init ();
  ncm_cfg_enable_gsl_err_handler ();
  
  g_test_add ("/numcosmo/nc_cluster_pseudo_counts/1p2_integral/", TestNcClusterPseudoCounts, NULL, 
              &test_nc_cluster_pseudo_counts_new, 
              &test_nc_cluster_pseudo_counts_1p2_integral, 
              &test_nc_cluster_pseudo_counts_free);
 
  g_test_run ();
}

void _set_destroyed (gpointer b) { gboolean *destroyed = b; *destroyed = TRUE; }

void
test_nc_cluster_pseudo_counts_free (TestNcClusterPseudoCounts *test, gconstpointer pdata)
{  
  NC_CLUSTER_PSEUDO_COUNTS *cpc = NC_CLUSTER_PSEUDO_COUNTS (test->cpc);
  gboolean destroyed = FALSE;
  
  nc_matter_var_free (test->vp);
  nc_hicosmo_free (test->cosmo);
  nc_cluster_mass_free (test->mszl);
  
  g_object_set_data_full (G_OBJECT (cpc), "test-destroy", &destroyed, _set_destroyed);
  nc_cluster_pseudo_counts_free (cpc);
  g_assert (destroyed);
}

void
test_nc_cluster_pseudo_counts_new (TestNcClusterPseudoCounts *test, gconstpointer pdata)
{
  NcClusterPseudoCounts *cpc;
  NcClusterMass *mszl;
  NcHICosmo *cosmo;
  NcDistance *dist;
  NcWindow *wf; 
  NcTransferFunc *tf;
  NcMatterVar *vp;
  NcGrowthFunc *gf;
  NcMultiplicityFunc *mulf;
  NcMassFunction *mfp;
  gdouble m1, m2;
  gdouble z = g_test_rand_double_range (0.188, 0.890);

  cosmo = nc_hicosmo_new_from_name (NcHICosmo, "NcHICosmoDEXcdm");
  dist  = nc_distance_new (3.0);
  wf    = nc_window_new_from_name ("NcWindowTophat");
  tf    = nc_transfer_func_new_from_name ("NcTransferFuncEH");
  vp    = nc_matter_var_new (NcMatterVarStrategyFFT, wp, tf);
  gf    = nc_growth_func_new ();
  mulf  = nc_multiplicity_function_new ("NcMultiplicityFuncTinkerCrit{'Delta':<500.0>}");
  mfp   = nc_mass_function_new (dist, vp, gf, mulf);
  
  m1 = g_test_rand_double_range (1.235, 2.496); /* ln(M/M0), M0 = 10^14 h^-1 M_sun */
  m2 = m1 + 0.4;
  test->Mobs[0] = exp (m1) * 1.0e14;
  test->Mobs[1] = exp (m2) * 1.0e14; 
    
  cpc  = NC_CLUSTER_PSEUDO_COUNTS (nc_cluster_pseudo_counts_new (mfp));
  mszl = NC_CLUSTER_MASS_PLCL (nc_cluster_mass_new_from_name ("NcClusterMassPlCL{}"));

  g_assert (cpc != NULL);
  g_assert (mszl != NULL);
  test->cpc = NC_CLUSTER_PSEUDO_COUNTS (cpc);
  test->mszl = mszl;
  test->z = z;
  test->Mobs = Mobs;
  test->Mobs_params = Mobs_params;
  g_assert (NC_IS_CLUSTER_PSEUDO_COUNTS (cpc));

  nc_distance_free (dist);
  nc_window_free (wp);
  nc_transfer_func_free (tf);
  nc_growth_func_free (gf);
  nc_multiplicity_func_free (mulf);
  nc_mass_function_free (mfp);
}

void 
test_nc_cluster_mass_plcl_bounds_1p2_integral (TestNcClusterPseudoCounts *test, gconstpointer pdata)
{
  nc_matter_var_prepare (test->vp, test->cosmo);

  printf ("z = %.5g Msz = %.5g Ml = %.5g\n", test->z, test->Mobs[0], test->Mobs[1]);
}

