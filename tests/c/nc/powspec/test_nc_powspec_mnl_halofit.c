/***************************************************************************
 *            test_nc_powspec_mnl_halofit.c
 *
 *  Mon Oct 19 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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

#define TEST_NC_POWSPEC_MNL_HALOFIT_NZ 32
#define TEST_NC_POWSPEC_MNL_HALOFIT_NK 16
#define TEST_NC_POWSPEC_MNL_HALOFIT_NTHREADS 4
#define TEST_NC_POWSPEC_MNL_HALOFIT_NROUNDS 20

typedef struct _TestNcPowspecMNLHaloFit
{
  NcHICosmo *cosmo;
  NcPowspecML *ps_ml;
  NcPowspecMNLHaloFit *ps_nl;
  NcmVector *z;
  NcmVector *k;
  NcmMatrix *Pk_serial;
  NcmVector *Pk_z;
} TestNcPowspecMNLHaloFit;

typedef struct _TestNcPowspecMNLHaloFitThread
{
  TestNcPowspecMNLHaloFit *test;
  guint offset;
  gboolean vec;
  NcmMatrix *Pk;
  NcmVector *Pk_z;
} TestNcPowspecMNLHaloFitThread;

void test_nc_powspec_mnl_halofit_new (TestNcPowspecMNLHaloFit *test, gconstpointer pdata);
void test_nc_powspec_mnl_halofit_free (TestNcPowspecMNLHaloFit *test, gconstpointer pdata);

void test_nc_powspec_mnl_halofit_concurrent_eval (TestNcPowspecMNLHaloFit *test, gconstpointer pdata);
void test_nc_powspec_mnl_halofit_concurrent_eval_vec (TestNcPowspecMNLHaloFit *test, gconstpointer pdata);
void test_nc_powspec_mnl_halofit_reprepare (TestNcPowspecMNLHaloFit *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/nc/powspec/mnl/halofit/concurrent/eval", TestNcPowspecMNLHaloFit, NULL,
              &test_nc_powspec_mnl_halofit_new,
              &test_nc_powspec_mnl_halofit_concurrent_eval,
              &test_nc_powspec_mnl_halofit_free);

  g_test_add ("/nc/powspec/mnl/halofit/concurrent/eval_vec", TestNcPowspecMNLHaloFit, NULL,
              &test_nc_powspec_mnl_halofit_new,
              &test_nc_powspec_mnl_halofit_concurrent_eval_vec,
              &test_nc_powspec_mnl_halofit_free);

  g_test_add ("/nc/powspec/mnl/halofit/reprepare", TestNcPowspecMNLHaloFit, NULL,
              &test_nc_powspec_mnl_halofit_new,
              &test_nc_powspec_mnl_halofit_reprepare,
              &test_nc_powspec_mnl_halofit_free);

  g_test_run ();
}

/* Evaluates P(z_i, k_j) of the test's grid into @Pk with ncm_powspec_eval() or, with @vec,
 * ncm_powspec_eval_vec() into @Pk_z row by row, visiting the redshifts from @offset on so
 * that threads with different offsets evaluate different redshifts at the same time. */
static void
_test_nc_powspec_mnl_halofit_fill (TestNcPowspecMNLHaloFit *test, guint offset, gboolean vec, NcmMatrix *Pk, NcmVector *Pk_z)
{
  NcmPowspec *ps = NCM_POWSPEC (test->ps_nl);
  guint n, j;

  for (n = 0; n < TEST_NC_POWSPEC_MNL_HALOFIT_NZ; n++)
  {
    const guint i_z = (n + offset) % TEST_NC_POWSPEC_MNL_HALOFIT_NZ;
    const gdouble z = ncm_vector_get (test->z, i_z);

    if (vec)
    {
      ncm_powspec_eval_vec (ps, NCM_MODEL (test->cosmo), z, test->k, Pk_z);

      for (j = 0; j < TEST_NC_POWSPEC_MNL_HALOFIT_NK; j++)
        ncm_matrix_set (Pk, i_z, j, ncm_vector_get (Pk_z, j));
    }
    else
    {
      for (j = 0; j < TEST_NC_POWSPEC_MNL_HALOFIT_NK; j++)
        ncm_matrix_set (Pk, i_z, j, ncm_powspec_eval (ps, NCM_MODEL (test->cosmo), z, ncm_vector_get (test->k, j)));
    }
  }
}

static gpointer
_test_nc_powspec_mnl_halofit_thread (gpointer data)
{
  TestNcPowspecMNLHaloFitThread *th = (TestNcPowspecMNLHaloFitThread *) data;
  guint r;

  for (r = 0; r < TEST_NC_POWSPEC_MNL_HALOFIT_NROUNDS; r++)
    _test_nc_powspec_mnl_halofit_fill (th->test, th->offset + r, th->vec, th->Pk, th->Pk_z);

  return NULL;
}

/* Runs the threads over the grid and checks every value against the serial one, bit for bit. */
static void
_test_nc_powspec_mnl_halofit_concurrent (TestNcPowspecMNLHaloFit *test, gboolean vec)
{
  TestNcPowspecMNLHaloFitThread th[TEST_NC_POWSPEC_MNL_HALOFIT_NTHREADS];
  GThread *threads[TEST_NC_POWSPEC_MNL_HALOFIT_NTHREADS];
  guint t, i, j;

  for (t = 0; t < TEST_NC_POWSPEC_MNL_HALOFIT_NTHREADS; t++)
  {
    th[t].test   = test;
    th[t].offset = t * TEST_NC_POWSPEC_MNL_HALOFIT_NZ / TEST_NC_POWSPEC_MNL_HALOFIT_NTHREADS;
    th[t].vec    = vec;
    th[t].Pk     = ncm_matrix_new (TEST_NC_POWSPEC_MNL_HALOFIT_NZ, TEST_NC_POWSPEC_MNL_HALOFIT_NK);
    th[t].Pk_z   = ncm_vector_new (TEST_NC_POWSPEC_MNL_HALOFIT_NK);
    threads[t]   = g_thread_new ("halofit", &_test_nc_powspec_mnl_halofit_thread, &th[t]);
  }

  for (t = 0; t < TEST_NC_POWSPEC_MNL_HALOFIT_NTHREADS; t++)
  {
    g_thread_join (threads[t]);

    for (i = 0; i < TEST_NC_POWSPEC_MNL_HALOFIT_NZ; i++)
      for (j = 0; j < TEST_NC_POWSPEC_MNL_HALOFIT_NK; j++)
        g_assert_cmpfloat (ncm_matrix_get (th[t].Pk, i, j), ==, ncm_matrix_get (test->Pk_serial, i, j));

    ncm_matrix_free (th[t].Pk);
    ncm_vector_free (th[t].Pk_z);
  }
}

void
test_nc_powspec_mnl_halofit_new (TestNcPowspecMNLHaloFit *test, gconstpointer pdata)
{
  NcHIReion *reion   = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim *prim     = NC_HIPRIM (nc_hiprim_power_law_new ());
  NcTransferFunc *tf = NC_TRANSFER_FUNC (nc_transfer_func_eh_new ());
  guint i;

  test->cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new_full (reion, prim, NULL));
  test->ps_ml = NC_POWSPEC_ML (nc_powspec_ml_transfer_new (tf));
  test->ps_nl = nc_powspec_mnl_halofit_new (test->ps_ml, 3.0, 1.0e-4);
  test->z     = ncm_vector_new (TEST_NC_POWSPEC_MNL_HALOFIT_NZ);
  test->k     = ncm_vector_new (TEST_NC_POWSPEC_MNL_HALOFIT_NK);

  for (i = 0; i < TEST_NC_POWSPEC_MNL_HALOFIT_NZ; i++)
    ncm_vector_set (test->z, i, 2.5 * i / (TEST_NC_POWSPEC_MNL_HALOFIT_NZ - 1.0));

  for (i = 0; i < TEST_NC_POWSPEC_MNL_HALOFIT_NK; i++)
    ncm_vector_set (test->k, i, exp (log (1.0e-2) + log (1.0e3) * i / (TEST_NC_POWSPEC_MNL_HALOFIT_NK - 1.0)));

  ncm_powspec_prepare (NCM_POWSPEC (test->ps_nl), NCM_MODEL (test->cosmo));

  test->Pk_serial = ncm_matrix_new (TEST_NC_POWSPEC_MNL_HALOFIT_NZ, TEST_NC_POWSPEC_MNL_HALOFIT_NK);
  test->Pk_z      = ncm_vector_new (TEST_NC_POWSPEC_MNL_HALOFIT_NK);
  _test_nc_powspec_mnl_halofit_fill (test, 0, FALSE, test->Pk_serial, test->Pk_z);

  nc_transfer_func_free (tf);
  ncm_model_free (NCM_MODEL (reion));
  ncm_model_free (NCM_MODEL (prim));
}

void
test_nc_powspec_mnl_halofit_free (TestNcPowspecMNLHaloFit *test, gconstpointer pdata)
{
  ncm_matrix_free (test->Pk_serial);
  ncm_vector_free (test->Pk_z);
  ncm_vector_free (test->z);
  ncm_vector_free (test->k);
  NCM_TEST_FREE (ncm_powspec_free, NCM_POWSPEC (test->ps_nl));
  NCM_TEST_FREE (nc_powspec_ml_free, test->ps_ml);
  NCM_TEST_FREE (nc_hicosmo_free, test->cosmo);
}

void
test_nc_powspec_mnl_halofit_concurrent_eval (TestNcPowspecMNLHaloFit *test, gconstpointer pdata)
{
  _test_nc_powspec_mnl_halofit_concurrent (test, FALSE);
}

void
test_nc_powspec_mnl_halofit_concurrent_eval_vec (TestNcPowspecMNLHaloFit *test, gconstpointer pdata)
{
  _test_nc_powspec_mnl_halofit_concurrent (test, TRUE);
}

/* After a parameter change and a new prepare, the cached coefficients of the previous
 * preparation are not used: the values agree with a halofit object prepared only for the
 * new parameters. */
void
test_nc_powspec_mnl_halofit_reprepare (TestNcPowspecMNLHaloFit *test, gconstpointer pdata)
{
  NcPowspecMNLHaloFit *fresh = nc_powspec_mnl_halofit_new (test->ps_ml, 3.0, 1.0e-4);
  NcmMatrix *Pk              = ncm_matrix_new (TEST_NC_POWSPEC_MNL_HALOFIT_NZ, TEST_NC_POWSPEC_MNL_HALOFIT_NK);
  NcmMatrix *Pk_fresh        = ncm_matrix_new (TEST_NC_POWSPEC_MNL_HALOFIT_NZ, TEST_NC_POWSPEC_MNL_HALOFIT_NK);
  NcPowspecMNLHaloFit *ps_nl = test->ps_nl;
  guint i, j;

  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_C,
                            1.1 * ncm_model_orig_param_get (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_C));
  ncm_powspec_prepare (NCM_POWSPEC (test->ps_nl), NCM_MODEL (test->cosmo));
  _test_nc_powspec_mnl_halofit_fill (test, 0, FALSE, Pk, test->Pk_z);

  ncm_powspec_prepare (NCM_POWSPEC (fresh), NCM_MODEL (test->cosmo));
  test->ps_nl = fresh;
  _test_nc_powspec_mnl_halofit_fill (test, 0, FALSE, Pk_fresh, test->Pk_z);
  test->ps_nl = ps_nl;

  for (i = 0; i < TEST_NC_POWSPEC_MNL_HALOFIT_NZ; i++)
  {
    for (j = 0; j < TEST_NC_POWSPEC_MNL_HALOFIT_NK; j++)
    {
      g_assert_cmpfloat (ncm_matrix_get (Pk, i, j), ==, ncm_matrix_get (Pk_fresh, i, j));
      g_assert_cmpfloat (ncm_matrix_get (Pk, i, j), !=, ncm_matrix_get (test->Pk_serial, i, j));
    }
  }

  ncm_matrix_free (Pk);
  ncm_matrix_free (Pk_fresh);
  NCM_TEST_FREE (ncm_powspec_free, NCM_POWSPEC (fresh));
}

