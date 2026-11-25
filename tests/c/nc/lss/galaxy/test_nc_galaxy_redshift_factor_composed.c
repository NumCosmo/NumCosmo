/***************************************************************************
 *            test_nc_galaxy_redshift_factor_composed.c
 *
 *  Wed November 25 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_galaxy_redshift_factor_composed.c
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
#include <gsl/gsl_integration.h>

/* The highest order of the moments, and the panels and points of the reference rule */
#define TEST_MOM_N 32
#define TEST_REF_PANELS 4000
#define TEST_REF_GL 20

typedef struct _TestNcGalaxyRedshiftFactorComposed
{
  NcGalaxyRedshiftFactor *gsdr;
  NcGalaxyRedshiftFactorData *data;
  NcmMSet *mset;
} TestNcGalaxyRedshiftFactorComposed;

/* A case: the LSST SRD population, the photo-z and its width, and the selection window */
typedef struct _TestComposedCase
{
  gboolean lens;
  gdouble zp;
  gdouble sigma0;
  gdouble zp_min;
  gdouble zp_max;
} TestComposedCase;

static const TestComposedCase test_cases[2] =
{
  {FALSE, 0.6, 0.05, 0.0, 20.0},
  {TRUE, 1.2, 0.08, 0.8, 1.7},
};

/* A scheme whose density is -1 on [0, 1], for the abort on a negative density */
#define TEST_TYPE_NC_GALAXY_REDSHIFT_FACTOR_NEG (test_nc_galaxy_redshift_factor_neg_get_type ())

G_DECLARE_FINAL_TYPE (TestNcGalaxyRedshiftFactorNeg, test_nc_galaxy_redshift_factor_neg, TEST, NC_GALAXY_REDSHIFT_FACTOR_NEG, NcGalaxyRedshiftFactor)

struct _TestNcGalaxyRedshiftFactorNeg
{
  NcGalaxyRedshiftFactor parent_instance;
};

G_DEFINE_TYPE (TestNcGalaxyRedshiftFactorNeg, test_nc_galaxy_redshift_factor_neg, NC_TYPE_GALAXY_REDSHIFT_FACTOR)

static gdouble
_test_neg_integ_f (gpointer callback_data, const gdouble z, NcGalaxyRedshiftFactorData *data)
{
  return -1.0;
}

/* The integrand has no closure data: nothing to free, nothing to copy */
static void
_test_neg_integ_data_free (gpointer callback_data)
{
}

static gpointer
_test_neg_integ_data_copy (gpointer callback_data)
{
  return NULL;
}

static NcGalaxyRedshiftFactorIntegrand *
_test_neg_integ (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, gboolean use_lnp)
{
  return nc_galaxy_redshift_factor_integrand_new (&_test_neg_integ_f, &_test_neg_integ_data_free, &_test_neg_integ_data_copy, NULL, NULL);
}

static void
_test_neg_ldata_row (NcGalaxyRedshiftFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
}

static void
_test_neg_ldata_required_columns (NcGalaxyRedshiftFactorData *data, GList **columns)
{
}

/* A galaxy with no data of its own: the fragment hooks do nothing */
static void
_test_neg_data_init (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data)
{
  data->ldata                  = g_new0 (gint, 1);
  data->ldata_destroy          = &g_free;
  data->ldata_read_row         = &_test_neg_ldata_row;
  data->ldata_write_row        = &_test_neg_ldata_row;
  data->ldata_required_columns = &_test_neg_ldata_required_columns;
}

static void
test_nc_galaxy_redshift_factor_neg_init (TestNcGalaxyRedshiftFactorNeg *gsdr)
{
}

static void
test_nc_galaxy_redshift_factor_neg_class_init (TestNcGalaxyRedshiftFactorNegClass *klass)
{
  NcGalaxyRedshiftFactorClass *gsdr_class = NC_GALAXY_REDSHIFT_FACTOR_CLASS (klass);

  gsdr_class->data_init = &_test_neg_data_init;
  gsdr_class->integ     = &_test_neg_integ;
}

void test_nc_galaxy_redshift_factor_composed_new (TestNcGalaxyRedshiftFactorComposed *test, gconstpointer pdata);
void test_nc_galaxy_redshift_factor_composed_free (TestNcGalaxyRedshiftFactorComposed *test, gconstpointer pdata);

void test_nc_galaxy_redshift_factor_composed_moments (TestNcGalaxyRedshiftFactorComposed *test, gconstpointer pdata);
void test_nc_galaxy_redshift_factor_composed_cheb_reltol (TestNcGalaxyRedshiftFactorComposed *test, gconstpointer pdata);

void test_nc_galaxy_redshift_factor_composed_empty (void);
void test_nc_galaxy_redshift_factor_composed_empty_subprocess (void);
void test_nc_galaxy_redshift_factor_composed_negative (void);
void test_nc_galaxy_redshift_factor_composed_negative_subprocess (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/galaxy_redshift_factor_composed/moments/source", TestNcGalaxyRedshiftFactorComposed, &test_cases[0],
              &test_nc_galaxy_redshift_factor_composed_new,
              &test_nc_galaxy_redshift_factor_composed_moments,
              &test_nc_galaxy_redshift_factor_composed_free);
  g_test_add ("/nc/galaxy_redshift_factor_composed/moments/lens", TestNcGalaxyRedshiftFactorComposed, &test_cases[1],
              &test_nc_galaxy_redshift_factor_composed_new,
              &test_nc_galaxy_redshift_factor_composed_moments,
              &test_nc_galaxy_redshift_factor_composed_free);
  g_test_add ("/nc/galaxy_redshift_factor_composed/cheb_reltol", TestNcGalaxyRedshiftFactorComposed, &test_cases[0],
              &test_nc_galaxy_redshift_factor_composed_new,
              &test_nc_galaxy_redshift_factor_composed_cheb_reltol,
              &test_nc_galaxy_redshift_factor_composed_free);

  g_test_add_func ("/nc/galaxy_redshift_factor_composed/empty", &test_nc_galaxy_redshift_factor_composed_empty);
  g_test_add_func ("/nc/galaxy_redshift_factor_composed/empty/subprocess", &test_nc_galaxy_redshift_factor_composed_empty_subprocess);
  g_test_add_func ("/nc/galaxy_redshift_factor_composed/negative", &test_nc_galaxy_redshift_factor_composed_negative);
  g_test_add_func ("/nc/galaxy_redshift_factor_composed/negative/subprocess", &test_nc_galaxy_redshift_factor_composed_negative_subprocess);

  g_test_run ();

  return 0;
}

/* The scheme of a case and a galaxy with its photo-z and width */
static void
_test_composed_setup (const TestComposedCase *tc, NcGalaxyRedshiftFactor **gsdr, NcmMSet **mset, NcGalaxyRedshiftFactorData **data)
{
  NcGalaxyRedshiftPopLSSTSRD *pop = tc->lens ? nc_galaxy_redshift_pop_lsst_srd_new_y10_lens () : nc_galaxy_redshift_pop_lsst_srd_new_y1_source ();
  NcGalaxyRedshiftObsGauss *obs   = nc_galaxy_redshift_obs_gauss_new ();
  GStrvBuilder *builder           = g_strv_builder_new ();
  NcGalaxyWLObs *wlobs;
  GList *columns, *l;
  GStrv cols;

  *mset = ncm_mset_empty_new ();
  ncm_mset_set (*mset, NCM_MODEL (pop), NULL);
  ncm_mset_set (*mset, NCM_MODEL (obs), NULL);

  *gsdr = NC_GALAXY_REDSHIFT_FACTOR (nc_galaxy_redshift_factor_composed_new (tc->zp_min, tc->zp_max));
  *data = nc_galaxy_redshift_factor_data_new (*gsdr, *mset);

  columns = nc_galaxy_redshift_factor_data_required_columns (*data);

  for (l = columns; l != NULL; l = g_list_next (l))
    g_strv_builder_add (builder, l->data);

  cols  = g_strv_builder_end (builder);
  wlobs = nc_galaxy_wl_obs_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, NC_WL_ELLIPTICITY_FRAME_CELESTIAL, 1, cols);
  nc_galaxy_wl_obs_set (wlobs, "z", 0, 0.0, NULL);
  nc_galaxy_wl_obs_set (wlobs, "zp", 0, tc->zp, NULL);
  nc_galaxy_wl_obs_set (wlobs, "sigma0", 0, tc->sigma0, NULL);
  nc_galaxy_redshift_factor_data_read_row (*data, wlobs, 0);
  nc_galaxy_redshift_factor_prepare (*gsdr, *mset);

  g_strfreev (cols);
  g_strv_builder_unref (builder);
  g_list_free_full (columns, g_free);
  nc_galaxy_wl_obs_free (wlobs);
  nc_galaxy_redshift_pop_lsst_srd_free (pop);
  nc_galaxy_redshift_obs_gauss_free (obs);
}

/*
 * The moments of the integrand on [a, b], k <= TEST_MOM_N, by the TEST_REF_GL-point
 * Gauss-Legendre rule (a table of GSL) on TEST_REF_PANELS equal panels.
 */
static void
_test_moments_ref (TestNcGalaxyRedshiftFactorComposed *test, const gdouble a, const gdouble b, gdouble *m)
{
  NcGalaxyRedshiftFactorIntegrand *integ = nc_galaxy_redshift_factor_integ (test->gsdr, test->mset, FALSE);
  gsl_integration_glfixed_table *glt     = gsl_integration_glfixed_table_alloc (TEST_REF_GL);
  const gdouble h                        = (b - a) / TEST_REF_PANELS;
  guint i, j, k;

  nc_galaxy_redshift_factor_integrand_prepare (integ, test->mset);

  for (k = 0; k <= TEST_MOM_N; k++)
    m[k] = 0.0;

  for (i = 0; i < TEST_REF_PANELS; i++)
  {
    for (j = 0; j < TEST_REF_GL; j++)
    {
      gdouble zq, wq;

      gsl_integration_glfixed_point (a + i * h, a + (i + 1) * h, j, &zq, &wq, glt);
      {
        const gdouble t  = GSL_MAX (-1.0, GSL_MIN (1.0, (2.0 * zq - (a + b)) / (b - a)));
        const gdouble pw = wq * nc_galaxy_redshift_factor_integrand_eval (integ, zq, test->data);

        for (k = 0; k <= TEST_MOM_N; k++)
          m[k] += pw * cos (k * acos (t));
      }
    }
  }

  gsl_integration_glfixed_table_free (glt);
  nc_galaxy_redshift_factor_integrand_free (integ);
}

void
test_nc_galaxy_redshift_factor_composed_new (TestNcGalaxyRedshiftFactorComposed *test, gconstpointer pdata)
{
  _test_composed_setup (pdata, &test->gsdr, &test->mset, &test->data);
}

void
test_nc_galaxy_redshift_factor_composed_free (TestNcGalaxyRedshiftFactorComposed *test, gconstpointer pdata)
{
  nc_galaxy_redshift_factor_data_unref (test->data);
  ncm_mset_free (test->mset);
  NCM_TEST_FREE (nc_galaxy_redshift_factor_free, test->gsdr);
}

/*
 * The default moments, from the panels of sqrt (p) at the default cheb-reltol, on the
 * integration range of the galaxy and on a part of it, against _test_moments_ref(). The
 * difference measured 1.6e-15 of m_0 at most; 1e-13 of m_0 is asserted.
 */
void
test_nc_galaxy_redshift_factor_composed_moments (TestNcGalaxyRedshiftFactorComposed *test, gconstpointer pdata)
{
  GArray *m = NULL;
  gdouble m_ref[TEST_MOM_N + 1];
  gdouble z_lo, z_hi;
  guint r, k;

  nc_galaxy_redshift_factor_get_integ_lim (test->gsdr, test->mset, test->data, &z_lo, &z_hi);

  for (r = 0; r < 2; r++)
  {
    const gdouble a = (r == 0) ? z_lo : z_lo + 0.3 * (z_hi - z_lo);
    const gdouble b = (r == 0) ? z_hi : z_lo + 0.7 * (z_hi - z_lo);

    nc_galaxy_redshift_factor_cheb_moments (test->gsdr, test->mset, test->data, a, b, TEST_MOM_N, &m);
    _test_moments_ref (test, a, b, m_ref);

    g_assert_cmpuint (m->len, ==, TEST_MOM_N + 1);

    for (k = 0; k <= TEST_MOM_N; k++)
      ncm_assert_cmpdouble_e (g_array_index (m, gdouble, k), ==, m_ref[k], 0.0, 1.0e-13 * m_ref[0]);
  }

  g_array_unref (m);
}

/* The default of cheb-reltol, its setter and its property */
void
test_nc_galaxy_redshift_factor_composed_cheb_reltol (TestNcGalaxyRedshiftFactorComposed *test, gconstpointer pdata)
{
  gdouble v;

  g_assert_cmpfloat (nc_galaxy_redshift_factor_get_cheb_reltol (test->gsdr), ==, 1.0e-10);

  nc_galaxy_redshift_factor_set_cheb_reltol (test->gsdr, 1.0e-8);
  g_assert_cmpfloat (nc_galaxy_redshift_factor_get_cheb_reltol (test->gsdr), ==, 1.0e-8);

  g_object_set (test->gsdr, "cheb-reltol", 1.0e-12, NULL);
  g_object_get (test->gsdr, "cheb-reltol", &v, NULL);
  g_assert_cmpfloat (v, ==, 1.0e-12);
}

void
test_nc_galaxy_redshift_factor_composed_empty (void)
{
  g_test_trap_subprocess ("/nc/galaxy_redshift_factor_composed/empty/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the interval [0.5, 0.5] is empty*");
}

void
test_nc_galaxy_redshift_factor_composed_empty_subprocess (void)
{
  NcGalaxyRedshiftFactor *gsdr;
  NcGalaxyRedshiftFactorData *data;
  NcmMSet *mset;
  GArray *m = NULL;

  _test_composed_setup (&test_cases[0], &gsdr, &mset, &data);
  nc_galaxy_redshift_factor_cheb_moments (gsdr, mset, data, 0.5, 0.5, 2, &m);
}

void
test_nc_galaxy_redshift_factor_composed_negative (void)
{
  g_test_trap_subprocess ("/nc/galaxy_redshift_factor_composed/negative/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*TestNcGalaxyRedshiftFactorNeg gives p(z) = -1 < 0 at z = *");
}

void
test_nc_galaxy_redshift_factor_composed_negative_subprocess (void)
{
  NcGalaxyRedshiftFactor *gsdr     = g_object_new (TEST_TYPE_NC_GALAXY_REDSHIFT_FACTOR_NEG, NULL);
  NcmMSet *mset                    = ncm_mset_empty_new ();
  NcGalaxyRedshiftFactorData *data = nc_galaxy_redshift_factor_data_new (gsdr, mset);
  GArray *m                        = NULL;

  nc_galaxy_redshift_factor_cheb_moments (gsdr, mset, data, 0.0, 1.0, 2, &m);
}

