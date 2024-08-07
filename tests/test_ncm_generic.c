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
void test_ncm_mpi_job_basic (void);
void test_ncm_mpi_job_test_basic (void);
void test_ncm_mpi_job_fit_basic (void);
void test_ncm_mpi_job_mcmc_basic (void);
void test_ncm_mpi_job_feval_basic (void);
void test_ncm_powspec_spline2d_basic (void);

void test_nc_de_cont_basic (void);
void test_nc_galaxy_sd_position_basic (void);
void test_nc_galaxy_sd_position_flat_basic (void);
void test_nc_galaxy_sd_position_lsst_srd_basic (void);
void test_nc_galaxy_sd_position_lsst_srd_y10_basic (void);
void test_nc_galaxy_sd_shape_basic (void);
void test_nc_galaxy_sd_shape_gauss_basic (void);
void test_nc_galaxy_sd_z_proxy_basic (void);
void test_nc_galaxy_sd_z_proxy_dirac_basic (void);
void test_nc_galaxy_sd_z_proxy_gauss_basic (void);
void test_nc_hicosmo_qgw_basic (void);
void test_nc_hipert_adiab_basic (void);
void test_nc_hipert_em_basic (void);
void test_nc_hipert_gw_basic (void);
void test_nc_hipert_two_fluids_basic (void);
void test_nc_hiprim_two_fluids_basic (void);

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
  g_test_add_func ("/ncm/mpi_job/basic", test_ncm_mpi_job_basic);
  g_test_add_func ("/ncm/mpi_job_test/basic", test_ncm_mpi_job_test_basic);
  g_test_add_func ("/ncm/mpi_job_fit/basic", test_ncm_mpi_job_fit_basic);
  g_test_add_func ("/ncm/mpi_job_mcmc/basic", test_ncm_mpi_job_mcmc_basic);
  g_test_add_func ("/ncm/mpi_job_feval/basic", test_ncm_mpi_job_feval_basic);
  g_test_add_func ("/ncm/powspec_spline2d/basic", test_ncm_powspec_spline2d_basic);

  g_test_add_func ("/nc/de_cont/basic", test_nc_de_cont_basic);
  g_test_add_func ("/nc/galaxy/sd_position/basic", test_nc_galaxy_sd_position_basic);
  g_test_add_func ("/nc/galaxy/sd_position_flat/basic", test_nc_galaxy_sd_position_flat_basic);
  g_test_add_func ("/nc/galaxy/sd_position_lsst_srd/basic", test_nc_galaxy_sd_position_lsst_srd_basic);
  g_test_add_func ("/nc/galaxy/sd_position_lsst_srd_y10/basic", test_nc_galaxy_sd_position_lsst_srd_y10_basic);
  g_test_add_func ("/nc/galaxy/sd_shape/basic", test_nc_galaxy_sd_shape_basic);
  g_test_add_func ("/nc/galaxy/sd_shape_gauss/basic", test_nc_galaxy_sd_shape_gauss_basic);
  g_test_add_func ("/nc/galaxy/sd_z_proxy/basic", test_nc_galaxy_sd_z_proxy_basic);
  g_test_add_func ("/nc/galaxy/sd_z_proxy_dirac/basic", test_nc_galaxy_sd_z_proxy_dirac_basic);
  g_test_add_func ("/nc/galaxy/sd_z_proxy_gauss/basic", test_nc_galaxy_sd_z_proxy_gauss_basic);
  g_test_add_func ("/nc/hicosmo/qgw/basic", test_nc_hicosmo_qgw_basic);

  g_test_add_func ("/nc/hipert/adiab/basic", test_nc_hipert_adiab_basic);
  g_test_add_func ("/nc/hipert/em/basic", test_nc_hipert_em_basic);
  g_test_add_func ("/nc/hipert/gw/basic", test_nc_hipert_gw_basic);
  g_test_add_func ("/nc/hipert/two_fluids/basic", test_nc_hipert_two_fluids_basic);
  g_test_add_func ("/nc/hiprim/two_fluids/basic", test_nc_hiprim_two_fluids_basic);

  g_test_run ();
}

void
test_ncm_data_basic (void)
{
  NcmData *data = NCM_DATA (ncm_data_funnel_new ());
  NcmData *data2;

  g_assert_true (data != NULL);
  g_assert_true (NCM_IS_DATA (data));

  data2 = ncm_data_ref (data);
  ncm_data_clear (&data2);
  g_assert_true (data2 == NULL);

  g_assert_true (NCM_IS_DATA (data));

  NCM_TEST_FREE (ncm_data_free, data);
}

void
test_ncm_data_funnel_basic (void)
{
  NcmDataFunnel *funnel = ncm_data_funnel_new ();
  NcmDataFunnel *funnel2;

  g_assert_true (funnel != NULL);
  g_assert_true (NCM_IS_DATA_FUNNEL (funnel));

  funnel2 = ncm_data_funnel_ref (funnel);
  ncm_data_funnel_clear (&funnel2);
  g_assert_true (funnel2 == NULL);

  g_assert_true (NCM_IS_DATA_FUNNEL (funnel));

  NCM_TEST_FREE (ncm_data_funnel_free, funnel);
}

void
test_ncm_data_gaussmix2d_basic (void)
{
  NcmDataGaussMix2D *gm2d = ncm_data_gaussmix2d_new ();
  NcmDataGaussMix2D *gm2d2;

  g_assert_true (gm2d != NULL);
  g_assert_true (NCM_IS_DATA_GAUSSMIX2D (gm2d));

  gm2d2 = ncm_data_gaussmix2d_ref (gm2d);
  ncm_data_gaussmix2d_clear (&gm2d2);
  g_assert_true (gm2d2 == NULL);

  g_assert_true (NCM_IS_DATA_GAUSSMIX2D (gm2d));

  NCM_TEST_FREE (ncm_data_gaussmix2d_free, gm2d);
}

void
test_ncm_data_rosenbrock_basic (void)
{
  NcmDataRosenbrock *rosenbrock = ncm_data_rosenbrock_new ();
  NcmDataRosenbrock *rosenbrock2;

  g_assert_true (rosenbrock != NULL);
  g_assert_true (NCM_IS_DATA_ROSENBROCK (rosenbrock));

  rosenbrock2 = ncm_data_rosenbrock_ref (rosenbrock);
  ncm_data_rosenbrock_clear (&rosenbrock2);
  g_assert_true (rosenbrock2 == NULL);

  g_assert_true (NCM_IS_DATA_ROSENBROCK (rosenbrock));

  NCM_TEST_FREE (ncm_data_rosenbrock_free, rosenbrock);
}

void
test_ncm_dataset_basic (void)
{
  NcmDataset *dset = ncm_dataset_new ();
  NcmDataset *dset2;

  g_assert_true (dset != NULL);
  g_assert_true (NCM_IS_DATASET (dset));

  dset2 = ncm_dataset_ref (dset);
  ncm_dataset_clear (&dset2);
  g_assert_true (dset2 == NULL);

  g_assert_true (NCM_IS_DATASET (dset));

  NCM_TEST_FREE (ncm_dataset_free, dset);
}

void
test_ncm_fftlog_basic (void)
{
  NcmFftlog *fftlog = NCM_FFTLOG (ncm_fftlog_gausswin2_new (0.0, 0.0, 20.0, 1000));
  NcmFftlog *fftlog2;

  g_assert_true (fftlog != NULL);
  g_assert_true (NCM_IS_FFTLOG (fftlog));

  fftlog2 = ncm_fftlog_ref (fftlog);
  ncm_fftlog_clear (&fftlog2);
  g_assert_true (fftlog2 == NULL);

  g_assert_true (NCM_IS_FFTLOG (fftlog));

  NCM_TEST_FREE (ncm_fftlog_free, fftlog);
}

void
test_ncm_mpi_job_basic (void)
{
  NcmMPIJob *mpi_job = NCM_MPI_JOB (ncm_mpi_job_test_new ());
  NcmMPIJob *mpi_job2;

  g_assert_true (mpi_job != NULL);
  g_assert_true (NCM_IS_MPI_JOB (mpi_job));

  mpi_job2 = ncm_mpi_job_ref (mpi_job);
  ncm_mpi_job_clear (&mpi_job2);
  g_assert_true (mpi_job2 == NULL);

  g_assert_true (NCM_IS_MPI_JOB (mpi_job));

  NCM_TEST_FREE (ncm_mpi_job_free, mpi_job);
}

void
test_ncm_mpi_job_test_basic (void)
{
  NcmMPIJobTest *mpi_job = ncm_mpi_job_test_new ();
  NcmMPIJobTest *mpi_job2;

  g_assert_true (mpi_job != NULL);
  g_assert_true (NCM_IS_MPI_JOB_TEST (mpi_job));

  mpi_job2 = ncm_mpi_job_test_ref (mpi_job);
  ncm_mpi_job_test_clear (&mpi_job2);
  g_assert_true (mpi_job2 == NULL);

  g_assert_true (NCM_IS_MPI_JOB_TEST (mpi_job));

  NCM_TEST_FREE (ncm_mpi_job_test_free, mpi_job);
}

void
test_ncm_mpi_job_fit_basic (void)
{
  NcmMSet *mset             = ncm_mset_empty_new ();
  NcmData *data             = NCM_DATA (ncm_data_funnel_new ());
  NcmDataset *dset          = ncm_dataset_new_list (data, NULL);
  NcmLikelihood *likelihood = ncm_likelihood_new (dset);
  NcmFit *fit               = ncm_fit_factory (NCM_FIT_TYPE_GSL_MMS, NULL, likelihood, mset, NCM_FIT_GRAD_NUMDIFF_FORWARD);
  NcmMPIJobFit *mpi_job     = ncm_mpi_job_fit_new (fit, NULL);
  NcmMPIJobFit *mpi_job2;

  g_assert_true (mpi_job != NULL);
  g_assert_true (NCM_IS_MPI_JOB_FIT (mpi_job));

  mpi_job2 = ncm_mpi_job_fit_ref (mpi_job);
  ncm_mpi_job_fit_clear (&mpi_job2);
  g_assert_true (mpi_job2 == NULL);

  g_assert_true (NCM_IS_MPI_JOB_FIT (mpi_job));

  NCM_TEST_FREE (ncm_mpi_job_fit_free, mpi_job);

  NCM_TEST_FREE (ncm_fit_free, fit);
  NCM_TEST_FREE (ncm_likelihood_free, likelihood);
  NCM_TEST_FREE (ncm_dataset_free, dset);
  NCM_TEST_FREE (ncm_data_free, data);
  NCM_TEST_FREE (ncm_mset_free, mset);
}

void
test_ncm_mpi_job_mcmc_basic (void)
{
  NcmMSet *mset             = ncm_mset_empty_new ();
  NcmData *data             = NCM_DATA (ncm_data_funnel_new ());
  NcmDataset *dset          = ncm_dataset_new_list (data, NULL);
  NcmLikelihood *likelihood = ncm_likelihood_new (dset);
  NcmFit *fit               = ncm_fit_factory (NCM_FIT_TYPE_GSL_MMS, NULL, likelihood, mset, NCM_FIT_GRAD_NUMDIFF_FORWARD);
  NcmMPIJobMCMC *mpi_job    = ncm_mpi_job_mcmc_new (fit, NULL);
  NcmMPIJobMCMC *mpi_job2;

  g_assert_true (mpi_job != NULL);
  g_assert_true (NCM_IS_MPI_JOB_MCMC (mpi_job));

  mpi_job2 = ncm_mpi_job_mcmc_ref (mpi_job);
  ncm_mpi_job_mcmc_clear (&mpi_job2);
  g_assert_true (mpi_job2 == NULL);

  g_assert_true (NCM_IS_MPI_JOB_MCMC (mpi_job));

  NCM_TEST_FREE (ncm_mpi_job_mcmc_free, mpi_job);

  NCM_TEST_FREE (ncm_fit_free, fit);
  NCM_TEST_FREE (ncm_likelihood_free, likelihood);
  NCM_TEST_FREE (ncm_dataset_free, dset);
  NCM_TEST_FREE (ncm_data_free, data);
  NCM_TEST_FREE (ncm_mset_free, mset);
}

void
test_ncm_mpi_job_feval_basic (void)
{
  NcmMSet *mset             = ncm_mset_empty_new ();
  NcmData *data             = NCM_DATA (ncm_data_funnel_new ());
  NcmDataset *dset          = ncm_dataset_new_list (data, NULL);
  NcmLikelihood *likelihood = ncm_likelihood_new (dset);
  NcmFit *fit               = ncm_fit_factory (NCM_FIT_TYPE_GSL_MMS, NULL, likelihood, mset, NCM_FIT_GRAD_NUMDIFF_FORWARD);
  NcmMPIJobFEval *mpi_job   = ncm_mpi_job_feval_new (fit, NULL);
  NcmMPIJobFEval *mpi_job2;

  g_assert_true (mpi_job != NULL);
  g_assert_true (NCM_IS_MPI_JOB_FEVAL (mpi_job));

  mpi_job2 = ncm_mpi_job_feval_ref (mpi_job);
  ncm_mpi_job_feval_clear (&mpi_job2);
  g_assert_true (mpi_job2 == NULL);

  g_assert_true (NCM_IS_MPI_JOB_FEVAL (mpi_job));

  NCM_TEST_FREE (ncm_mpi_job_feval_free, mpi_job);

  NCM_TEST_FREE (ncm_fit_free, fit);
  NCM_TEST_FREE (ncm_likelihood_free, likelihood);
  NCM_TEST_FREE (ncm_dataset_free, dset);
  NCM_TEST_FREE (ncm_data_free, data);
  NCM_TEST_FREE (ncm_mset_free, mset);
}

void
test_ncm_powspec_spline2d_basic (void)
{
  gdouble x[6]      = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  gdouble y[7]      = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  gdouble z[42]     = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, };
  NcmVector *xv     = ncm_vector_new_data_static (x, 6, 1);
  NcmVector *yv     = ncm_vector_new_data_static (y, 7, 1);
  NcmMatrix *zm     = ncm_matrix_new_data_static (z, 7, 6);
  NcmSpline *sc     = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  NcmSpline2d *sb2d = ncm_spline2d_bicubic_new (sc);

  ncm_spline2d_set (NCM_SPLINE2D (sb2d), xv, yv, zm, FALSE);

  {
    NcmPowspecSpline2d *ps_s2d = ncm_powspec_spline2d_new (NCM_SPLINE2D (sb2d));
    NcmPowspecSpline2d *ps_s2d2;

    g_assert_true (ps_s2d != NULL);
    g_assert_true (NCM_IS_POWSPEC_SPLINE2D (ps_s2d));

    ps_s2d2 = ncm_powspec_spline2d_ref (ps_s2d);
    ncm_powspec_spline2d_clear (&ps_s2d2);
    g_assert_true (ps_s2d2 == NULL);

    g_assert_true (NCM_IS_POWSPEC_SPLINE2D (ps_s2d));

    ncm_vector_free (xv);
    ncm_vector_free (yv);
    ncm_matrix_free (zm);
    ncm_spline_free (NCM_SPLINE (sc));

    NCM_TEST_FREE (ncm_powspec_spline2d_free, ps_s2d);
    NCM_TEST_FREE (ncm_spline2d_free, NCM_SPLINE2D (sb2d));
  }
}

void
test_nc_hicosmo_qgw_basic (void)
{
  NcHICosmoQGW *qgw = nc_hicosmo_qgw_new ();

  g_assert_true (qgw != NULL);
  g_assert_true (NC_IS_HICOSMO_QGW (qgw));

  NCM_TEST_FREE (nc_hicosmo_free, NC_HICOSMO (qgw));
}

void
test_nc_de_cont_basic (void)
{
  NcDECont *dec = nc_de_cont_new (0.3, 0.7, 0.01, 0.01);
  NcDECont *dec2;

  g_assert_true (dec != NULL);
  g_assert_true (NC_IS_DE_CONT (dec));

  dec2 = nc_de_cont_ref (dec);
  nc_de_cont_clear (&dec2);
  g_assert_true (dec2 == NULL);

  g_assert_true (NC_IS_DE_CONT (dec));

  nc_de_cont_set_k (dec, 0.1);
  g_assert_cmpfloat (nc_de_cont_get_k (dec), ==, 0.1);

  NCM_TEST_FREE (nc_de_cont_free, dec);
}

void
test_nc_galaxy_sd_position_basic (void)
{
  NcGalaxySDPosition *sdpos = NC_GALAXY_SD_POSITION (nc_galaxy_sd_position_flat_new (0, 2, 0, 2));
  NcGalaxySDPosition *sdpos2;

  g_assert_true (sdpos != NULL);
  g_assert_true (NC_IS_GALAXY_SD_POSITION (sdpos));

  sdpos2 = nc_galaxy_sd_position_ref (sdpos);
  nc_galaxy_sd_position_clear (&sdpos2);
  g_assert_true (sdpos2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_POSITION (sdpos));

  NCM_TEST_FREE (nc_galaxy_sd_position_free, sdpos);
}

void
test_nc_galaxy_sd_position_flat_basic (void)
{
  NcGalaxySDPositionFlat *gsdpf = nc_galaxy_sd_position_flat_new (0, 2, 0, 2);
  NcGalaxySDPositionFlat *gsdpf2;

  g_assert_true (gsdpf != NULL);
  g_assert_true (NC_IS_GALAXY_SD_POSITION_FLAT (gsdpf));

  gsdpf2 = nc_galaxy_sd_position_flat_ref (gsdpf);
  nc_galaxy_sd_position_flat_clear (&gsdpf2);
  g_assert_true (gsdpf2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_POSITION_FLAT (gsdpf));

  NCM_TEST_FREE (nc_galaxy_sd_position_flat_free, gsdpf);
}

void
test_nc_galaxy_sd_position_lsst_srd_basic (void)
{
  NcGalaxySDPositionLSSTSRD *gsdpls = nc_galaxy_sd_position_lsst_srd_new (0, 2, 0, 2);
  NcGalaxySDPositionLSSTSRD *gsdpls2;

  g_assert_true (gsdpls != NULL);
  g_assert_true (NC_IS_GALAXY_SD_POSITION_LSST_SRD (gsdpls));

  gsdpls2 = nc_galaxy_sd_position_lsst_srd_ref (gsdpls);
  nc_galaxy_sd_position_lsst_srd_clear (&gsdpls2);
  g_assert_true (gsdpls2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_POSITION_LSST_SRD (gsdpls));

  NCM_TEST_FREE (nc_galaxy_sd_position_lsst_srd_free, gsdpls);
}

void
test_nc_galaxy_sd_position_lsst_srd_y10_basic (void)
{
  NcGalaxySDPositionLSSTSRD *gsdplsy10 = nc_galaxy_sd_position_lsst_srd_new_y10 (0, 2, 0, 2);
  NcGalaxySDPositionLSSTSRD *gsdplsy102;

  g_assert_true (gsdplsy10 != NULL);
  g_assert_true (NC_IS_GALAXY_SD_POSITION_LSST_SRD (gsdplsy10));

  gsdplsy102 = nc_galaxy_sd_position_lsst_srd_ref (gsdplsy10);
  nc_galaxy_sd_position_lsst_srd_clear (&gsdplsy102);
  g_assert_true (gsdplsy102 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_POSITION_LSST_SRD (gsdplsy10));

  NCM_TEST_FREE (nc_galaxy_sd_position_lsst_srd_free, gsdplsy10);
}

void
test_nc_galaxy_sd_shape_basic (void)
{
  NcGalaxySDShape *gsds = NC_GALAXY_SD_SHAPE (nc_galaxy_sd_shape_gauss_new (0.0, 2.0, 0.05));
  NcGalaxySDShape *gsds2;

  g_assert_true (gsds != NULL);
  g_assert_true (NC_IS_GALAXY_SD_SHAPE (gsds));

  gsds2 = nc_galaxy_sd_shape_ref (gsds);
  nc_galaxy_sd_shape_clear (&gsds2);
  g_assert_true (gsds2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_SHAPE (gsds));

  NCM_TEST_FREE (nc_galaxy_sd_shape_free, gsds);
}

void
test_nc_galaxy_sd_shape_gauss_basic (void)
{
  NcGalaxySDShapeGauss *gsdsg = nc_galaxy_sd_shape_gauss_new (0.0, 2.0, 0.05);
  NcGalaxySDShapeGauss *gsdsg2;

  g_assert_true (gsdsg != NULL);
  g_assert_true (NC_IS_GALAXY_SD_SHAPE_GAUSS (gsdsg));

  gsdsg2 = nc_galaxy_sd_shape_gauss_ref (gsdsg);
  nc_galaxy_sd_shape_gauss_clear (&gsdsg2);
  g_assert_true (gsdsg2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_SHAPE_GAUSS (gsdsg));

  NCM_TEST_FREE (nc_galaxy_sd_shape_gauss_free, gsdsg);
}

void
test_nc_galaxy_sd_z_proxy_basic (void)
{
  NcGalaxySDZProxy *gsdzp = NC_GALAXY_SD_Z_PROXY (nc_galaxy_sd_z_proxy_dirac_new ());
  NcGalaxySDZProxy *gsdzp2;

  g_assert_true (gsdzp != NULL);
  g_assert_true (NC_IS_GALAXY_SD_Z_PROXY (gsdzp));

  gsdzp2 = nc_galaxy_sd_z_proxy_ref (gsdzp);
  nc_galaxy_sd_z_proxy_clear (&gsdzp2);
  g_assert_true (gsdzp2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_Z_PROXY (gsdzp));

  NCM_TEST_FREE (nc_galaxy_sd_z_proxy_free, gsdzp);
}

void
test_nc_galaxy_sd_z_proxy_dirac_basic (void)
{
  NcGalaxySDZProxyDirac *gsdzpd = nc_galaxy_sd_z_proxy_dirac_new ();
  NcGalaxySDZProxyDirac *gsdzpd2;

  g_assert_true (gsdzpd != NULL);
  g_assert_true (NC_IS_GALAXY_SD_Z_PROXY_DIRAC (gsdzpd));

  gsdzpd2 = nc_galaxy_sd_z_proxy_dirac_ref (gsdzpd);
  nc_galaxy_sd_z_proxy_dirac_clear (&gsdzpd2);
  g_assert_true (gsdzpd2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_Z_PROXY_DIRAC (gsdzpd));

  NCM_TEST_FREE (nc_galaxy_sd_z_proxy_dirac_free, gsdzpd);
}

void
test_nc_galaxy_sd_z_proxy_gauss_basic (void)
{
  NcGalaxySDZProxyGauss *gsdzpg = nc_galaxy_sd_z_proxy_gauss_new (0.0, 2.0, 0.05);
  NcGalaxySDZProxyGauss *gsdzpg2;

  g_assert_true (gsdzpg != NULL);
  g_assert_true (NC_IS_GALAXY_SD_Z_PROXY_GAUSS (gsdzpg));

  gsdzpg2 = nc_galaxy_sd_z_proxy_gauss_ref (gsdzpg);
  nc_galaxy_sd_z_proxy_gauss_clear (&gsdzpg2);
  g_assert_true (gsdzpg2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_Z_PROXY_GAUSS (gsdzpg));

  NCM_TEST_FREE (nc_galaxy_sd_z_proxy_gauss_free, gsdzpg);
}

void
test_nc_hipert_adiab_basic (void)
{
  NcHIPertAdiab *adiab = nc_hipert_adiab_new ();
  NcHIPertAdiab *adiab2;

  g_assert_true (adiab != NULL);
  g_assert_true (NC_IS_HIPERT_ADIAB (adiab));

  adiab2 = nc_hipert_adiab_ref (adiab);
  nc_hipert_adiab_clear (&adiab2);
  g_assert_true (adiab2 == NULL);

  g_assert_true (NC_IS_HIPERT_ADIAB (adiab));

  nc_hipert_adiab_set_k (adiab, 0.1);
  g_assert_cmpfloat (nc_hipert_adiab_get_k (adiab), ==, 0.1);

  NCM_TEST_FREE (nc_hipert_adiab_free, adiab);
}

void
test_nc_hipert_em_basic (void)
{
  NcHIPertEM *em = nc_hipert_em_new ();
  NcHIPertEM *em2;

  g_assert_true (em != NULL);
  g_assert_true (NC_IS_HIPERT_EM (em));

  em2 = nc_hipert_em_ref (em);
  nc_hipert_em_clear (&em2);
  g_assert_true (em2 == NULL);

  g_assert_true (NC_IS_HIPERT_EM (em));

  nc_hipert_em_set_k (em, 0.1);
  g_assert_cmpfloat (nc_hipert_em_get_k (em), ==, 0.1);

  NCM_TEST_FREE (nc_hipert_em_free, em);
}

void
test_nc_hipert_gw_basic (void)
{
  NcHIPertGW *gw = nc_hipert_gw_new ();
  NcHIPertGW *gw2;

  g_assert_true (gw != NULL);
  g_assert_true (NC_IS_HIPERT_GW (gw));

  gw2 = nc_hipert_gw_ref (gw);
  nc_hipert_gw_clear (&gw2);
  g_assert_true (gw2 == NULL);

  g_assert_true (NC_IS_HIPERT_GW (gw));

  nc_hipert_gw_set_k (gw, 0.1);
  g_assert_cmpfloat (nc_hipert_gw_get_k (gw), ==, 0.1);

  NCM_TEST_FREE (nc_hipert_gw_free, gw);
}

void
test_nc_hipert_two_fluids_basic (void)
{
  NcHIPertTwoFluids *tf = nc_hipert_two_fluids_new ();
  NcHIPertTwoFluids *tf2;

  g_assert_true (tf != NULL);
  g_assert_true (NC_IS_HIPERT_TWO_FLUIDS (tf));

  tf2 = nc_hipert_two_fluids_ref (tf);
  nc_hipert_two_fluids_clear (&tf2);
  g_assert_true (tf2 == NULL);

  g_assert_true (NC_IS_HIPERT_TWO_FLUIDS (tf));

  NCM_TEST_FREE (nc_hipert_two_fluids_free, tf);
}

void
test_nc_hiprim_two_fluids_basic (void)
{
  NcHIPrimTwoFluids *tf = nc_hiprim_two_fluids_new ();
  NcHIPrimTwoFluids *tf2;

  g_assert_true (tf != NULL);
  g_assert_true (NC_IS_HIPRIM_TWO_FLUIDS (tf));

  tf2 = nc_hiprim_two_fluids_ref (tf);
  nc_hiprim_two_fluids_clear (&tf2);
  g_assert_true (tf2 == NULL);

  g_assert_true (NC_IS_HIPRIM_TWO_FLUIDS (tf));

  NCM_TEST_FREE (nc_hiprim_two_fluids_free, tf);
}

