/***************************************************************************
 *            test_nc_halo_catalog_generator.c
 *
 *  Sat Jun 13 18:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_halo_catalog_generator.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

typedef struct _TestNcHaloCatalogGenerator
{
  NcClusterAbundance *cad;
  NcHaloCatalogGenerator *gen;
} TestNcHaloCatalogGenerator;

void test_nc_halo_catalog_generator_new (TestNcHaloCatalogGenerator *test, gconstpointer pdata);
void test_nc_halo_catalog_generator_free (TestNcHaloCatalogGenerator *test, gconstpointer pdata);

void test_nc_halo_catalog_generator_ref (TestNcHaloCatalogGenerator *test, gconstpointer pdata);
void test_nc_halo_catalog_generator_abundance (TestNcHaloCatalogGenerator *test, gconstpointer pdata);
void test_nc_halo_catalog_generator_footprint (TestNcHaloCatalogGenerator *test, gconstpointer pdata);
void test_nc_halo_catalog_generator_with_radius (TestNcHaloCatalogGenerator *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/halo_catalog_generator/ref", TestNcHaloCatalogGenerator, NULL,
              &test_nc_halo_catalog_generator_new,
              &test_nc_halo_catalog_generator_ref,
              &test_nc_halo_catalog_generator_free);

  g_test_add ("/nc/halo_catalog_generator/abundance", TestNcHaloCatalogGenerator, NULL,
              &test_nc_halo_catalog_generator_new,
              &test_nc_halo_catalog_generator_abundance,
              &test_nc_halo_catalog_generator_free);

  g_test_add ("/nc/halo_catalog_generator/footprint", TestNcHaloCatalogGenerator, NULL,
              &test_nc_halo_catalog_generator_new,
              &test_nc_halo_catalog_generator_footprint,
              &test_nc_halo_catalog_generator_free);

  g_test_add ("/nc/halo_catalog_generator/with_radius", TestNcHaloCatalogGenerator, NULL,
              &test_nc_halo_catalog_generator_new,
              &test_nc_halo_catalog_generator_with_radius,
              &test_nc_halo_catalog_generator_free);

  g_test_run ();

  return 0;
}

void
test_nc_halo_catalog_generator_new (TestNcHaloCatalogGenerator *test, gconstpointer pdata)
{
  NcDistance *dist         = nc_distance_new (3.0);
  NcTransferFunc *tf       = nc_transfer_func_eh_new ();
  NcPowspecML *ps_ml       = NC_POWSPEC_ML (nc_powspec_ml_transfer_new (tf));
  NcmPowspecFilter *psf    = ncm_powspec_filter_new (NCM_POWSPEC (ps_ml), NCM_POWSPEC_FILTER_TYPE_TOPHAT);
  NcMultiplicityFunc *mulf = NC_MULTIPLICITY_FUNC (nc_multiplicity_func_tinker_new_full (NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL, 200));
  NcHaloMassFunction *mfp  = nc_halo_mass_function_new (dist, psf, mulf);

  test->cad = nc_cluster_abundance_new (mfp, NULL);
  test->gen = nc_halo_catalog_generator_new (test->cad);

  g_assert_true (NC_IS_HALO_CATALOG_GENERATOR (test->gen));

  nc_halo_mass_function_free (mfp);
  nc_multiplicity_func_free (mulf);
  ncm_powspec_filter_free (psf);
  nc_powspec_ml_free (ps_ml);
  nc_transfer_func_free (tf);
  nc_distance_free (dist);
}

void
test_nc_halo_catalog_generator_free (TestNcHaloCatalogGenerator *test, gconstpointer pdata)
{
  nc_cluster_abundance_free (test->cad);
  NCM_TEST_FREE (nc_halo_catalog_generator_free, test->gen);
}

void
test_nc_halo_catalog_generator_ref (TestNcHaloCatalogGenerator *test, gconstpointer pdata)
{
  NcHaloCatalogGenerator *gen_ref = nc_halo_catalog_generator_ref (test->gen);

  g_assert_true (gen_ref == test->gen);

  nc_halo_catalog_generator_clear (&gen_ref);
  g_assert_null (gen_ref);

  g_assert_true (NC_IS_HALO_CATALOG_GENERATOR (test->gen));
}

void
test_nc_halo_catalog_generator_abundance (TestNcHaloCatalogGenerator *test, gconstpointer pdata)
{
  g_assert_true (nc_halo_catalog_generator_peek_abundance (test->gen) == test->cad);
}

void
test_nc_halo_catalog_generator_footprint (TestNcHaloCatalogGenerator *test, gconstpointer pdata)
{
  NcmSkyFootprintRectangular *rect = ncm_sky_footprint_rectangular_new (10.0, 40.0, -5.0, 25.0);

  g_assert_null (nc_halo_catalog_generator_peek_footprint (test->gen));

  nc_halo_catalog_generator_set_footprint (test->gen, NCM_SKY_FOOTPRINT (rect));
  g_assert_true (nc_halo_catalog_generator_peek_footprint (test->gen) == NCM_SKY_FOOTPRINT (rect));

  nc_halo_catalog_generator_set_footprint (test->gen, NULL);
  g_assert_null (nc_halo_catalog_generator_peek_footprint (test->gen));

  ncm_sky_footprint_rectangular_free (rect);
}

void
test_nc_halo_catalog_generator_with_radius (TestNcHaloCatalogGenerator *test, gconstpointer pdata)
{
  g_assert_false (nc_halo_catalog_generator_get_with_radius (test->gen));

  nc_halo_catalog_generator_set_with_radius (test->gen, TRUE);
  g_assert_true (nc_halo_catalog_generator_get_with_radius (test->gen));

  nc_halo_catalog_generator_set_with_radius (test->gen, FALSE);
  g_assert_false (nc_halo_catalog_generator_get_with_radius (test->gen));
}
