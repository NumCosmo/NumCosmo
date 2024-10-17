/***************************************************************************
 *            test_nc_galaxy_sd_shape.c
 *
 *  mon May 07 00:02:42 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Caio Lima de Oliveira 2024 <caiolimadeoliveira@pm.me>
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

typedef struct _TestNcGalaxySDShapeGauss
{
  NcGalaxySDObsRedshift *gsdor;
  NcGalaxySDPosition *gsdp;
  NcGalaxySDShape *gsds;
} TestNcGalaxySDShapeGauss;


static void test_nc_galaxy_sd_shape_gauss_new (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_free (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_serialize (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_model_id (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  /* g_test_set_nonfatal_assertions (); */

  g_test_add ("/nc/galaxy_sd_shape/gauss/serialize", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_serialize,
              &test_nc_galaxy_sd_shape_free);

  g_test_add ("/nc/galaxy_sd_shape/gauss/model_id", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_model_id,
              &test_nc_galaxy_sd_shape_free);

  g_test_run ();

  return 0;
}

static void
test_nc_galaxy_sd_shape_gauss_new (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDTrueRedshift *gsdtr = NC_GALAXY_SD_TRUE_REDSHIFT (nc_galaxy_sd_true_redshift_lsst_srd_new (0.0, 1100.0));
  NcGalaxySDObsRedshift *gsdor  = NC_GALAXY_SD_OBS_REDSHIFT (nc_galaxy_sd_obs_redshift_spec_new (gsdtr));
  NcGalaxySDPosition *gsdp      = NC_GALAXY_SD_POSITION (nc_galaxy_sd_position_flat_new (0.0, 360.0, -90.0, 90.0));
  NcGalaxySDShapeGauss *gsds    = nc_galaxy_sd_shape_gauss_new ();

  test->gsdor = gsdor;
  test->gsdp  = gsdp;
  test->gsds  = NC_GALAXY_SD_SHAPE (gsds);

  nc_galaxy_sd_true_redshift_free (gsdtr);

  g_assert_true (NC_IS_GALAXY_SD_SHAPE (gsds));
  g_assert_true (NC_IS_GALAXY_SD_SHAPE_GAUSS (gsds));
}

static void
test_nc_galaxy_sd_shape_free (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  nc_galaxy_sd_obs_redshift_free (test->gsdor);
  nc_galaxy_sd_position_free (test->gsdp);
  NCM_TEST_FREE (nc_galaxy_sd_shape_free, test->gsds);
}

static void
test_nc_galaxy_sd_shape_serialize (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDShape *gsds     = test->gsds;
  NcmSerialize *ser         = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  gchar *sdsg_ser           = ncm_serialize_to_string (ser, G_OBJECT (gsds), TRUE);
  NcGalaxySDShape *sdsg_dup = NC_GALAXY_SD_SHAPE (ncm_serialize_from_string (ser, sdsg_ser));

  g_assert_true (NC_IS_GALAXY_SD_SHAPE (sdsg_dup));

  ncm_serialize_free (ser);
  g_free (sdsg_ser);
  nc_galaxy_sd_shape_free (sdsg_dup);
}

static void
test_nc_galaxy_sd_shape_model_id (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcmMSet *model_set  = ncm_mset_empty_new ();
  NcmSerialize *ser   = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmModel *model_dup = ncm_model_dup (NCM_MODEL (test->gsds), ser);

  ncm_mset_set (model_set, model_dup, NULL);

  g_assert_true (NC_IS_GALAXY_SD_SHAPE (ncm_mset_peek (model_set, nc_galaxy_sd_shape_id ())));

  ncm_model_free (model_dup);
  ncm_mset_free (model_set);
  ncm_serialize_free (ser);
}

