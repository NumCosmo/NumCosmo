/***************************************************************************
 *            test_ncm_diff.c
 *
 *  Wed July 26 12:04:44 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_diff.c
 *
 * Copyright (C) 2017 - Sandro Dias Pinto Vitenti
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

void test_ncm_data_basic (void);
void test_ncm_data_funnel_basic (void);
void test_ncm_data_gaussmix2d_basic (void);
void test_ncm_data_rosenbrock_basic (void);
void test_ncm_dataset_basic (void);
void test_ncm_fftlog_basic (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/data/basic", test_ncm_data_basic);
  g_test_add_func ("/ncm/data_funnel/basic", test_ncm_data_funnel_basic);
  g_test_add_func ("/ncm/data_gaussmix2d/basic", test_ncm_data_gaussmix2d_basic);
  g_test_add_func ("/ncm/data_rosenbrock/basic", test_ncm_data_rosenbrock_basic);
  g_test_add_func ("/ncm/dataset/basic", test_ncm_dataset_basic);
  g_test_add_func ("/ncm/fftlog/basic", test_ncm_fftlog_basic);

  g_test_run ();
}

void
test_ncm_data_basic (void)
{
  NcmData *data = NCM_DATA (ncm_data_funnel_new ());
  NcmData *data2;

  g_assert (data != NULL);
  g_assert (NCM_IS_DATA (data));

  data2 = ncm_data_ref (data);
  ncm_data_clear (&data2);
  g_assert (data2 == NULL);

  g_assert (NCM_IS_DATA (data));

  NCM_TEST_FREE (ncm_data_free, data);
}

void
test_ncm_data_funnel_basic (void)
{
  NcmDataFunnel *funnel = ncm_data_funnel_new ();
  NcmDataFunnel *funnel2;

  g_assert (funnel != NULL);
  g_assert (NCM_IS_DATA_FUNNEL (funnel));

  funnel2 = ncm_data_funnel_ref (funnel);
  ncm_data_funnel_clear (&funnel2);
  g_assert (funnel2 == NULL);

  g_assert (NCM_IS_DATA_FUNNEL (funnel));

  NCM_TEST_FREE (ncm_data_funnel_free, funnel);
}

void
test_ncm_data_gaussmix2d_basic (void)
{
  NcmDataGaussMix2D *gm2d = ncm_data_gaussmix2d_new ();
  NcmDataGaussMix2D *gm2d2;

  g_assert (gm2d != NULL);
  g_assert (NCM_IS_DATA_GAUSSMIX2D (gm2d));

  gm2d2 = ncm_data_gaussmix2d_ref (gm2d);
  ncm_data_gaussmix2d_clear (&gm2d2);
  g_assert (gm2d2 == NULL);

  g_assert (NCM_IS_DATA_GAUSSMIX2D (gm2d));

  NCM_TEST_FREE (ncm_data_gaussmix2d_free, gm2d);
}

void
test_ncm_data_rosenbrock_basic (void)
{
  NcmDataRosenbrock *rosenbrock = ncm_data_rosenbrock_new ();
  NcmDataRosenbrock *rosenbrock2;

  g_assert (rosenbrock != NULL);
  g_assert (NCM_IS_DATA_ROSENBROCK (rosenbrock));

  rosenbrock2 = ncm_data_rosenbrock_ref (rosenbrock);
  ncm_data_rosenbrock_clear (&rosenbrock2);
  g_assert (rosenbrock2 == NULL);

  g_assert (NCM_IS_DATA_ROSENBROCK (rosenbrock));

  NCM_TEST_FREE (ncm_data_rosenbrock_free, rosenbrock);
}

void
test_ncm_dataset_basic (void)
{
  NcmDataset *dset = ncm_dataset_new ();
  NcmDataset *dset2;

  g_assert (dset != NULL);
  g_assert (NCM_IS_DATASET (dset));

  dset2 = ncm_dataset_ref (dset);
  ncm_dataset_clear (&dset2);
  g_assert (dset2 == NULL);

  g_assert (NCM_IS_DATASET (dset));

  NCM_TEST_FREE (ncm_dataset_free, dset);
}

void
test_ncm_fftlog_basic (void)
{
  NcmFftlog *fftlog = NCM_FFTLOG (ncm_fftlog_gausswin2_new (0.0, 0.0, 20.0, 1000));
  NcmFftlog *fftlog2;

  g_assert (fftlog != NULL);
  g_assert (NCM_IS_FFTLOG (fftlog));

  fftlog2 = ncm_fftlog_ref (fftlog);
  ncm_fftlog_clear (&fftlog2);
  g_assert (fftlog2 == NULL);

  g_assert (NCM_IS_FFTLOG (fftlog));

  NCM_TEST_FREE (ncm_fftlog_free, fftlog);
}

