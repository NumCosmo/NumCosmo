/***************************************************************************
 *            test_nc_galaxy_shape_factor_wtg.c
 *
 *  Fri Nov 27 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * test_nc_galaxy_shape_factor_wtg.c
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * The shear-scatter populations and the WtG shape factor: reference counting through
 * every wrapper, the properties, the acceptance of a population by a shape factor, and
 * the two kernels against values computed here.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

static gint
_compare_double (gconstpointer a, gconstpointer b)
{
  const gdouble x = *(const gdouble *) a;
  const gdouble y = *(const gdouble *) b;

  return (x > y) - (x < y);
}

void test_nc_galaxy_shape_factor_wtg_ref (void);
void test_nc_galaxy_shape_pop_shear_scatter_gauss_ref (void);
void test_nc_galaxy_shape_pop_shear_scatter_gauss_props (void);
void test_nc_galaxy_shape_pop_shear_scatter_voigt_ref (void);
void test_nc_galaxy_shape_pop_shear_scatter_base_ref (void);
void test_nc_galaxy_shape_pop_shear_scatter_include_noise (void);
void test_nc_galaxy_shape_factor_wtg_accepts_pop (void);
void test_nc_galaxy_shape_pop_shear_scatter_gauss_kernel (void);
void test_nc_galaxy_shape_pop_shear_scatter_voigt_kernel (void);
void test_nc_galaxy_shape_pop_shear_scatter_gen (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/nc/galaxy_shape_factor_wtg/ref", &test_nc_galaxy_shape_factor_wtg_ref);
  g_test_add_func ("/nc/galaxy_shape_factor_wtg/accepts_pop", &test_nc_galaxy_shape_factor_wtg_accepts_pop);
  g_test_add_func ("/nc/galaxy_shape_pop_shear_scatter/base/ref", &test_nc_galaxy_shape_pop_shear_scatter_base_ref);
  g_test_add_func ("/nc/galaxy_shape_pop_shear_scatter/include_noise", &test_nc_galaxy_shape_pop_shear_scatter_include_noise);
  g_test_add_func ("/nc/galaxy_shape_pop_shear_scatter/gauss/ref", &test_nc_galaxy_shape_pop_shear_scatter_gauss_ref);
  g_test_add_func ("/nc/galaxy_shape_pop_shear_scatter/gauss/props", &test_nc_galaxy_shape_pop_shear_scatter_gauss_props);
  g_test_add_func ("/nc/galaxy_shape_pop_shear_scatter/gauss/kernel", &test_nc_galaxy_shape_pop_shear_scatter_gauss_kernel);
  g_test_add_func ("/nc/galaxy_shape_pop_shear_scatter/voigt/ref", &test_nc_galaxy_shape_pop_shear_scatter_voigt_ref);
  g_test_add_func ("/nc/galaxy_shape_pop_shear_scatter/voigt/kernel", &test_nc_galaxy_shape_pop_shear_scatter_voigt_kernel);
  g_test_add_func ("/nc/galaxy_shape_pop_shear_scatter/gen", &test_nc_galaxy_shape_pop_shear_scatter_gen);

  g_test_run ();

  return 0;
}

void
test_nc_galaxy_shape_factor_wtg_ref (void)
{
  NcGalaxyShapeFactorWtG *wtg = nc_galaxy_shape_factor_wtg_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  NcGalaxyShapeFactorWtG *ref = nc_galaxy_shape_factor_wtg_ref (wtg);

  g_assert_true (NC_IS_GALAXY_SHAPE_FACTOR_WTG (wtg));
  g_assert_true (NC_IS_GALAXY_SHAPE_FACTOR (wtg));
  g_assert_true (ref == wtg);
  g_assert_cmpint (nc_galaxy_shape_factor_get_ellip_conv (NC_GALAXY_SHAPE_FACTOR (wtg)), ==, NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);

  nc_galaxy_shape_factor_wtg_free (ref);
  nc_galaxy_shape_factor_wtg_clear (&wtg);
  g_assert_null (wtg);
}

void
test_nc_galaxy_shape_pop_shear_scatter_gauss_ref (void)
{
  NcGalaxyShapePopShearScatterGauss *gauss = nc_galaxy_shape_pop_shear_scatter_gauss_new (FALSE);
  NcGalaxyShapePopShearScatterGauss *ref   = nc_galaxy_shape_pop_shear_scatter_gauss_ref (gauss);

  g_assert_true (NC_IS_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS (gauss));
  g_assert_true (NC_IS_GALAXY_SHAPE_POP_SHEAR_SCATTER (gauss));
  g_assert_true (NC_IS_GALAXY_SHAPE_POP (gauss));
  g_assert_true (ref == gauss);

  nc_galaxy_shape_pop_shear_scatter_gauss_free (ref);
  nc_galaxy_shape_pop_shear_scatter_gauss_clear (&gauss);
  g_assert_null (gauss);
}

void
test_nc_galaxy_shape_pop_shear_scatter_voigt_ref (void)
{
  NcGalaxyShapePopShearScatterVoigt *voigt = nc_galaxy_shape_pop_shear_scatter_voigt_new ();
  NcGalaxyShapePopShearScatterVoigt *ref   = nc_galaxy_shape_pop_shear_scatter_voigt_ref (voigt);

  g_assert_true (NC_IS_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT (voigt));
  g_assert_true (NC_IS_GALAXY_SHAPE_POP_SHEAR_SCATTER (voigt));
  g_assert_true (ref == voigt);

  nc_galaxy_shape_pop_shear_scatter_voigt_free (ref);
  nc_galaxy_shape_pop_shear_scatter_voigt_clear (&voigt);
  g_assert_null (voigt);
}

/* The abstract class's own wrappers, reached through a concrete instance. */
void
test_nc_galaxy_shape_pop_shear_scatter_base_ref (void)
{
  NcGalaxyShapePopShearScatterVoigt *voigt = nc_galaxy_shape_pop_shear_scatter_voigt_new ();
  NcGalaxyShapePopShearScatter *base       = NC_GALAXY_SHAPE_POP_SHEAR_SCATTER (voigt);
  NcGalaxyShapePopShearScatter *ref        = nc_galaxy_shape_pop_shear_scatter_ref (base);

  g_assert_true (ref == base);

  nc_galaxy_shape_pop_shear_scatter_free (ref);

  ref = nc_galaxy_shape_pop_shear_scatter_ref (base);
  nc_galaxy_shape_pop_shear_scatter_clear (&ref);
  g_assert_null (ref);

  nc_galaxy_shape_pop_shear_scatter_voigt_clear (&voigt);
}

void
test_nc_galaxy_shape_pop_shear_scatter_gauss_props (void)
{
  NcGalaxyShapePopShearScatterGauss *fixed   = nc_galaxy_shape_pop_shear_scatter_gauss_new (FALSE);
  NcGalaxyShapePopShearScatterGauss *catalog = nc_galaxy_shape_pop_shear_scatter_gauss_new (TRUE);
  gboolean from_catalog                      = FALSE;
  gboolean include_noise                     = FALSE;

  g_assert_false (nc_galaxy_shape_pop_shear_scatter_gauss_get_sigma_from_catalog (fixed));
  g_assert_true (nc_galaxy_shape_pop_shear_scatter_gauss_get_sigma_from_catalog (catalog));

  g_object_get (catalog, "sigma-from-catalog", &from_catalog, "include-noise", &include_noise, NULL);
  g_assert_true (from_catalog);
  g_assert_true (include_noise);

  g_object_get (fixed, "sigma-from-catalog", &from_catalog, NULL);
  g_assert_false (from_catalog);

  /* A catalog-width population needs the e_rms column, a fixed-width one nothing. */
  {
    NcGalaxyShapePopData *data_c = nc_galaxy_shape_pop_data_new (NC_GALAXY_SHAPE_POP (catalog));
    NcGalaxyShapePopData *data_f = nc_galaxy_shape_pop_data_new (NC_GALAXY_SHAPE_POP (fixed));
    GList *cols_c                = nc_galaxy_shape_pop_data_required_columns (data_c);
    GList *cols_f                = nc_galaxy_shape_pop_data_required_columns (data_f);

    g_assert_cmpuint (g_list_length (cols_c), ==, 1);
    g_assert_cmpstr (cols_c->data, ==, NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS_COL_E_RMS);
    g_assert_null (cols_f);

    g_list_free_full (cols_c, g_free);
    nc_galaxy_shape_pop_data_unref (data_c);
    nc_galaxy_shape_pop_data_unref (data_f);
  }

  nc_galaxy_shape_pop_shear_scatter_gauss_clear (&fixed);
  nc_galaxy_shape_pop_shear_scatter_gauss_clear (&catalog);
}

void
test_nc_galaxy_shape_pop_shear_scatter_include_noise (void)
{
  NcGalaxyShapePopShearScatterGauss *gauss = nc_galaxy_shape_pop_shear_scatter_gauss_new (FALSE);
  NcGalaxyShapePopShearScatter *base       = NC_GALAXY_SHAPE_POP_SHEAR_SCATTER (gauss);
  NcGalaxyShapePopData *data               = nc_galaxy_shape_pop_data_new (NC_GALAXY_SHAPE_POP (gauss));
  const gdouble sigma                      = 0.3;
  const gdouble var_noise                  = 0.04;
  guint64 pkey0;

  ncm_model_param_set_by_name (NCM_MODEL (gauss), "sigma", sigma, NULL);
  nc_galaxy_shape_pop_prepare (NC_GALAXY_SHAPE_POP (gauss), data);

  g_assert_true (nc_galaxy_shape_pop_shear_scatter_get_include_noise (base));

  /* With the noise: variance sigma^2 + var_noise */
  {
    const gdouble var = sigma * sigma + var_noise;

    ncm_assert_cmpdouble_e (nc_galaxy_shape_pop_shear_scatter_eval_ln_p (base, data, 0.0, var_noise, 1.0),
                            ==, -0.5 * log (2.0 * M_PI * var), 1.0e-14, 0.0);
  }

  pkey0 = ncm_model_state_get_pkey (NCM_MODEL (gauss));
  nc_galaxy_shape_pop_shear_scatter_set_include_noise (base, FALSE);
  g_assert_false (nc_galaxy_shape_pop_shear_scatter_get_include_noise (base));
  /* a change of the setting must reach the caches that hash the model */
  g_assert_cmpuint (ncm_model_state_get_pkey (NCM_MODEL (gauss)), !=, pkey0);

  /* Without it: variance sigma^2 */
  ncm_assert_cmpdouble_e (nc_galaxy_shape_pop_shear_scatter_eval_ln_p (base, data, 0.0, var_noise, 1.0),
                          ==, -0.5 * log (2.0 * M_PI * sigma * sigma), 1.0e-14, 0.0);

  /* setting the same value again changes nothing */
  pkey0 = ncm_model_state_get_pkey (NCM_MODEL (gauss));
  nc_galaxy_shape_pop_shear_scatter_set_include_noise (base, FALSE);
  g_assert_cmpuint (ncm_model_state_get_pkey (NCM_MODEL (gauss)), ==, pkey0);

  nc_galaxy_shape_pop_data_unref (data);
  nc_galaxy_shape_pop_shear_scatter_gauss_clear (&gauss);
}

void
test_nc_galaxy_shape_factor_wtg_accepts_pop (void)
{
  NcGalaxyShapeFactor *wtg                 = NC_GALAXY_SHAPE_FACTOR (nc_galaxy_shape_factor_wtg_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET));
  NcGalaxyShapeFactor *var_add             = NC_GALAXY_SHAPE_FACTOR (nc_galaxy_shape_factor_var_add_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET));
  NcGalaxyShapePopShearScatterGauss *gauss = nc_galaxy_shape_pop_shear_scatter_gauss_new (FALSE);
  NcGalaxyShapePopShearScatterVoigt *voigt = nc_galaxy_shape_pop_shear_scatter_voigt_new ();
  NcGalaxyShapePopGauss *disc              = nc_galaxy_shape_pop_gauss_new ();

  g_assert_true (nc_galaxy_shape_factor_accepts_pop (wtg, NC_GALAXY_SHAPE_POP (gauss)));
  g_assert_true (nc_galaxy_shape_factor_accepts_pop (wtg, NC_GALAXY_SHAPE_POP (voigt)));
  g_assert_false (nc_galaxy_shape_factor_accepts_pop (wtg, NC_GALAXY_SHAPE_POP (disc)));

  g_assert_false (nc_galaxy_shape_factor_accepts_pop (var_add, NC_GALAXY_SHAPE_POP (gauss)));
  g_assert_false (nc_galaxy_shape_factor_accepts_pop (var_add, NC_GALAXY_SHAPE_POP (voigt)));
  g_assert_true (nc_galaxy_shape_factor_accepts_pop (var_add, NC_GALAXY_SHAPE_POP (disc)));

  nc_galaxy_shape_factor_clear (&wtg);
  nc_galaxy_shape_factor_clear (&var_add);
  nc_galaxy_shape_pop_shear_scatter_gauss_clear (&gauss);
  nc_galaxy_shape_pop_shear_scatter_voigt_clear (&voigt);
  nc_galaxy_shape_pop_gauss_clear (&disc);
}

/* The width of a catalog-width Gaussian is e_rms * eps_to_g, whatever sigma holds. */
void
test_nc_galaxy_shape_pop_shear_scatter_gauss_kernel (void)
{
  NcGalaxyShapePopShearScatterGauss *gauss = nc_galaxy_shape_pop_shear_scatter_gauss_new (TRUE);
  NcGalaxyShapePopShearScatter *base       = NC_GALAXY_SHAPE_POP_SHEAR_SCATTER (gauss);
  NcGalaxyShapePopData *data               = nc_galaxy_shape_pop_data_new (NC_GALAXY_SHAPE_POP (gauss));
  const gdouble e_rms                      = 0.36;
  const gdouble eps_to_g                   = 1.0 / (2.0 * 0.85);
  const gdouble var_noise                  = 0.03 * 0.03;
  const gdouble delta                      = 0.07;
  const gdouble sigma_g                    = e_rms * eps_to_g;
  const gdouble var                        = sigma_g * sigma_g + var_noise;

  ncm_model_param_set_by_name (NCM_MODEL (gauss), "sigma", 0.9, NULL); /* unused in this mode */
  data->e_rms = e_rms;
  nc_galaxy_shape_pop_prepare (NC_GALAXY_SHAPE_POP (gauss), data);

  ncm_assert_cmpdouble_e (nc_galaxy_shape_pop_shear_scatter_eval_ln_p (base, data, delta, var_noise, eps_to_g),
                          ==, -0.5 * delta * delta / var - 0.5 * log (2.0 * M_PI * var), 1.0e-14, 0.0);
  ncm_assert_cmpdouble_e (nc_galaxy_shape_pop_shear_scatter_eval_p (base, data, delta, var_noise, eps_to_g),
                          ==, exp (-0.5 * delta * delta / var) / sqrt (2.0 * M_PI * var), 1.0e-14, 0.0);

  nc_galaxy_shape_pop_data_unref (data);
  nc_galaxy_shape_pop_shear_scatter_gauss_clear (&gauss);
}

/* V (0) = exp (x^2) erfc (x) / (sigma_t sqrt (2 pi)), x = Gamma / (sigma_t sqrt (2)). */
void
test_nc_galaxy_shape_pop_shear_scatter_voigt_kernel (void)
{
  NcGalaxyShapePopShearScatterVoigt *voigt = nc_galaxy_shape_pop_shear_scatter_voigt_new ();
  NcGalaxyShapePopShearScatter *base       = NC_GALAXY_SHAPE_POP_SHEAR_SCATTER (voigt);
  NcGalaxyShapePopData *data               = nc_galaxy_shape_pop_data_new (NC_GALAXY_SHAPE_POP (voigt));
  const gdouble sigma                      = 0.22;
  const gdouble gamma                      = 0.03;
  const gdouble var_noise                  = 0.05 * 0.05;
  const gdouble sigma_t                    = sqrt (sigma * sigma + var_noise);
  const gdouble x                          = gamma / (M_SQRT2 * sigma_t);

  ncm_model_param_set_by_name (NCM_MODEL (voigt), "sigma", sigma, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (voigt), "gamma", gamma, NULL);
  nc_galaxy_shape_pop_prepare (NC_GALAXY_SHAPE_POP (voigt), data);

  ncm_assert_cmpdouble_e (nc_galaxy_shape_pop_shear_scatter_eval_p (base, data, 0.0, var_noise, 1.0),
                          ==, exp (x * x) * erfc (x) / (sigma_t * sqrt (2.0 * M_PI)), 1.0e-12, 0.0);

  /* even in delta, and with heavier wings than the Gaussian of the same core */
  ncm_assert_cmpdouble_e (nc_galaxy_shape_pop_shear_scatter_eval_p (base, data, 0.4, var_noise, 1.0),
                          ==, nc_galaxy_shape_pop_shear_scatter_eval_p (base, data, -0.4, var_noise, 1.0), 1.0e-14, 0.0);
  g_assert_cmpfloat (nc_galaxy_shape_pop_shear_scatter_eval_p (base, data, 1.5, var_noise, 1.0), >,
                     exp (-0.5 * 1.5 * 1.5 / (sigma_t * sigma_t)) / (sigma_t * sqrt (2.0 * M_PI)));

  nc_galaxy_shape_pop_data_unref (data);
  nc_galaxy_shape_pop_shear_scatter_voigt_clear (&voigt);
}

/* The Gaussian draws have the kernel's variance; the Voigt draws its median and quartiles. */
void
test_nc_galaxy_shape_pop_shear_scatter_gen (void)
{
  NcmRNG *rng                              = ncm_rng_seeded_new (NULL, 20261127);
  NcGalaxyShapePopShearScatterGauss *gauss = nc_galaxy_shape_pop_shear_scatter_gauss_new (FALSE);
  NcGalaxyShapePopShearScatterVoigt *voigt = nc_galaxy_shape_pop_shear_scatter_voigt_new ();
  NcGalaxyShapePopData *dg                 = nc_galaxy_shape_pop_data_new (NC_GALAXY_SHAPE_POP (gauss));
  NcGalaxyShapePopData *dv                 = nc_galaxy_shape_pop_data_new (NC_GALAXY_SHAPE_POP (voigt));
  const gdouble sigma                      = 0.25;
  const gdouble gamma                      = 0.02;
  const gdouble var_noise                  = 0.04 * 0.04;
  const guint n                            = 20000;
  GArray *draws                            = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), n);
  gdouble sum2                             = 0.0;
  guint i;

  ncm_model_param_set_by_name (NCM_MODEL (gauss), "sigma", sigma, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (voigt), "sigma", sigma, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (voigt), "gamma", gamma, NULL);
  nc_galaxy_shape_pop_prepare (NC_GALAXY_SHAPE_POP (gauss), dg);
  nc_galaxy_shape_pop_prepare (NC_GALAXY_SHAPE_POP (voigt), dv);

  for (i = 0; i < n; i++)
  {
    const gdouble d = nc_galaxy_shape_pop_shear_scatter_gen_delta (NC_GALAXY_SHAPE_POP_SHEAR_SCATTER (gauss), dg, rng, var_noise, 1.0);

    sum2 += d * d;
  }

  {
    const gdouble var = sigma * sigma + var_noise;

    /* variance of n Gaussian draws: relative error about sqrt (2 / n) = 1% */
    g_assert_cmpfloat (fabs (sum2 / n / var - 1.0), <, 0.05);
  }

  for (i = 0; i < n; i++)
  {
    const gdouble d = nc_galaxy_shape_pop_shear_scatter_gen_delta (NC_GALAXY_SHAPE_POP_SHEAR_SCATTER (voigt), dv, rng, var_noise, 1.0);

    g_array_append_val (draws, d);
  }

  g_array_sort (draws, &_compare_double);

  {
    const gdouble *s = (const gdouble *) draws->data;

    /* the Voigt distribution is symmetric about zero: median 0, and the first and third
     * quartiles are mirror images of each other */
    const gdouble q1 = s[n / 4];
    const gdouble q3 = s[3 * n / 4];

    g_assert_cmpfloat (fabs (s[n / 2]), <, 0.01);
    g_assert_cmpfloat (fabs (q1 + q3), <, 0.01);
    g_assert_cmpfloat (q3, >, 0.0);
  }

  g_array_unref (draws);
  nc_galaxy_shape_pop_data_unref (dg);
  nc_galaxy_shape_pop_data_unref (dv);
  nc_galaxy_shape_pop_shear_scatter_gauss_clear (&gauss);
  nc_galaxy_shape_pop_shear_scatter_voigt_clear (&voigt);
  ncm_rng_free (rng);
}

