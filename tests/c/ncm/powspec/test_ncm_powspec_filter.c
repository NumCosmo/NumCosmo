/***************************************************************************
 *            test_ncm_powspec_filter.c
 *
 *  Mon September 28 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
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

#define TEST_PL_AMPLITUDE (1.0e-9)

/*
 * A tabulated power law P = A k^n on k in [1e-6, 1e3] and z in [z0, z0 + 1], which
 * NcmPowspecSpline2d continues exactly past the table. For it
 * sigma^2 is proportional to R^-(n+3) with either window.
 */
static NcmPowspec *
_test_power_law_new (const gdouble n, const gdouble z0)
{
  const guint Nk   = 3000;
  const guint Nz   = 10;
  NcmVector *zv    = ncm_vector_new (Nz);
  NcmVector *lnkv  = ncm_vector_new (Nk);
  NcmMatrix *lnP   = ncm_matrix_new (Nk, Nz);
  NcmSpline2d *s2d = NCM_SPLINE2D (ncm_spline2d_bicubic_notaknot_new ());
  NcmPowspec *ps;
  guint i, j;

  for (j = 0; j < Nz; j++)
    ncm_vector_set (zv, j, z0 + j / (Nz - 1.0));

  for (i = 0; i < Nk; i++)
  {
    const gdouble lnk = log (1.0e-6) + log (1.0e9) * i / (Nk - 1.0);

    ncm_vector_set (lnkv, i, lnk);

    for (j = 0; j < Nz; j++)
      ncm_matrix_set (lnP, i, j, log (TEST_PL_AMPLITUDE) + n * lnk);
  }

  ncm_spline2d_set (s2d, zv, lnkv, lnP, TRUE);
  ps = NCM_POWSPEC (ncm_powspec_spline2d_new (s2d));

  ncm_vector_free (zv);
  ncm_vector_free (lnkv);
  ncm_matrix_free (lnP);
  ncm_spline2d_free (s2d);

  return ps;
}

static NcmPowspecFilter *
_test_filter_new (NcmPowspec *ps, NcmPowspecFilterType type, const guint nderivs, const gdouble reltol)
{
  NcmPowspecFilter *psf = ncm_powspec_filter_new (ps, type);

  ncm_powspec_filter_require_nderivs (psf, nderivs);
  ncm_powspec_filter_set_reltol (psf, reltol);
  ncm_powspec_filter_set_best_lnr0 (psf);
  ncm_powspec_filter_prepare (psf, NULL);

  return psf;
}

/*
 * For the power laws n = -2 and -1.5, d ln sigma^2 / d ln r = -(n + 3) and the second
 * log-derivative vanishes, for R in [1, 20]; at reltol 1e-6 the largest deviation
 * measured is 5.6e-8 on the first and 2.0e-7 on the second.
 */
static void
test_ncm_powspec_filter_power_law (void)
{
  const NcmPowspecFilterType types[] = {NCM_POWSPEC_FILTER_TYPE_TOPHAT, NCM_POWSPEC_FILTER_TYPE_GAUSS};
  const gdouble ns[]                 = {-2.0, -1.5};
  guint t, a, i;

  for (t = 0; t < G_N_ELEMENTS (types); t++)
  {
    for (a = 0; a < G_N_ELEMENTS (ns); a++)
    {
      NcmPowspec *ps        = _test_power_law_new (ns[a], 0.0);
      NcmPowspecFilter *psf = _test_filter_new (ps, types[t], 2, 1.0e-6);

      for (i = 0; i < 9; i++)
      {
        const gdouble lnr = log (20.0) * i / 8.0;

        ncm_assert_cmpdouble_e (ncm_powspec_filter_eval_dnlnvar_dlnrn (psf, 0.0, lnr, 1), ==, -(ns[a] + 3.0), 0.0, 3.0e-7);
        ncm_assert_cmpdouble_e (ncm_powspec_filter_eval_dnlnvar_dlnrn (psf, 0.0, lnr, 2), ==, 0.0, 0.0, 1.0e-6);
      }

      ncm_powspec_filter_free (psf);
      ncm_powspec_free (ps);
    }
  }
}

/*
 * The Gaussian window W = exp (-k^2 R^2 / 2) on P = A k^n gives
 * sigma^2 = A Gamma ((3 + n) / 2) / (4 pi^2 R^(3+n)); within reltol 1e-6 for R in
 * [1, 20] (measured 1.1e-7).
 */
static void
test_ncm_powspec_filter_gauss_closed_form (void)
{
  const gdouble ns[] = {-2.0, -1.5};
  guint a, i;

  for (a = 0; a < G_N_ELEMENTS (ns); a++)
  {
    NcmPowspec *ps        = _test_power_law_new (ns[a], 0.0);
    NcmPowspecFilter *psf = _test_filter_new (ps, NCM_POWSPEC_FILTER_TYPE_GAUSS, 1, 1.0e-6);

    for (i = 0; i < 9; i++)
    {
      const gdouble lnr   = log (20.0) * i / 8.0;
      const gdouble truth = TEST_PL_AMPLITUDE * tgamma (0.5 * (3.0 + ns[a])) / (4.0 * M_PI * M_PI * exp ((3.0 + ns[a]) * lnr));

      ncm_assert_cmpdouble_e (ncm_powspec_filter_eval_var_lnr (psf, 0.0, lnr), ==, truth, 1.0e-6, 0.0);
    }

    ncm_powspec_filter_free (psf);
    ncm_powspec_free (ps);
  }
}

/*
 * The top-hat variance of NcmPowspecAnalytic (BBKS, LCDM growth) on k in [1e-6, 1e2]
 * against the Arb table of the quadrature over the table, at R = 1, 8 and 50 Mpc where
 * the continuation past kmax contributes little; at reltol 1e-6, measured 3.7e-8 (R = 1)
 * and 4.3e-8 (R = 50).
 */
static void
test_ncm_powspec_filter_arb (void)
{
  NcmSerialize *ser   = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  gchar *path         = ncm_cfg_get_data_filename ("truth_tables/powspec/ncm_powspec_analytic_integrals.bin", TRUE);
  NcmObjDictStr *dict = ncm_serialize_dict_str_from_binfile (ser, path);
  NcmVector *k_range  = NCM_VECTOR (ncm_obj_dict_str_peek (dict, "k_range"));
  NcmMatrix *var      = NCM_MATRIX (ncm_obj_dict_str_peek (dict, "var"));
  NcmPowspec *ps      = NCM_POWSPEC (ncm_powspec_analytic_new (NCM_POWSPEC_ANALYTIC_SHAPE_BBKS, NCM_POWSPEC_ANALYTIC_GROWTH_LCDM));
  NcmPowspecFilter *psf;
  guint i;

  ncm_powspec_set_kmin (ps, ncm_vector_get (k_range, 0));
  ncm_powspec_set_kmax (ps, ncm_vector_get (k_range, 1));
  ncm_powspec_set_zi (ps, 0.0);
  ncm_powspec_set_zf (ps, 1.0);

  psf = _test_filter_new (ps, NCM_POWSPEC_FILTER_TYPE_TOPHAT, 1, 1.0e-6);

  for (i = 0; i < ncm_matrix_nrows (var); i++)
  {
    const gdouble z = ncm_matrix_get (var, i, 0);
    const gdouble R = ncm_matrix_get (var, i, 1);

    if (R >= 1.0)
      ncm_assert_cmpdouble_e (ncm_powspec_filter_eval_var (psf, z, R), ==, ncm_matrix_get (var, i, 2), 1.0e-7, 0.0);
  }

  ncm_powspec_filter_free (psf);
  ncm_powspec_free (ps);
  ncm_obj_dict_str_unref (dict);
  ncm_serialize_free (ser);
  g_free (path);
}

/* The evaluation forms agree with one another at a point. */
static void
test_ncm_powspec_filter_eval_forms (void)
{
  NcmPowspec *ps        = _test_power_law_new (-1.5, 0.0);
  NcmPowspecFilter *psf = _test_filter_new (ps, NCM_POWSPEC_FILTER_TYPE_TOPHAT, 2, 1.0e-6);
  const gdouble rs[]    = {0.5, 3.0, 11.0};
  guint i;

  for (i = 0; i < G_N_ELEMENTS (rs); i++)
  {
    const gdouble r     = rs[i];
    const gdouble lnr   = log (r);
    const gdouble var   = ncm_powspec_filter_eval_var_lnr (psf, 0.0, lnr);
    const gdouble dvar  = ncm_powspec_filter_eval_dnvar_dlnrn (psf, 0.0, lnr, 1);
    const gdouble d2var = ncm_powspec_filter_eval_dnvar_dlnrn (psf, 0.0, lnr, 2);

    ncm_assert_cmpdouble_e (ncm_powspec_filter_eval_var (psf, 0.0, r), ==, var, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_powspec_filter_eval_sigma (psf, 0.0, r), ==, sqrt (var), 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_powspec_filter_eval_sigma_lnr (psf, 0.0, lnr), ==, sqrt (var), 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_powspec_filter_eval_lnvar_lnr (psf, 0.0, lnr), ==, log (var), 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_powspec_filter_eval_dnvar_dlnrn (psf, 0.0, lnr, 0), ==, var, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_powspec_filter_eval_dvar_dlnr (psf, 0.0, lnr), ==, dvar, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_powspec_filter_eval_dlnvar_dlnr (psf, 0.0, lnr), ==, dvar / var, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_powspec_filter_eval_dlnvar_dr (psf, 0.0, lnr), ==, dvar / var / r, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_powspec_filter_eval_dnlnvar_dlnrn (psf, 0.0, lnr, 0), ==, log (var), 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_powspec_filter_eval_dnlnvar_dlnrn (psf, 0.0, lnr, 1), ==, dvar / var, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_powspec_filter_eval_dnlnvar_dlnrn (psf, 0.0, lnr, 2), ==, d2var / var - gsl_pow_2 (dvar / var), 1.0e-12, 0.0);
  }

  ncm_assert_cmpdouble_e (ncm_powspec_filter_volume_rm3 (psf), ==, 4.0 * M_PI / 3.0, 1.0e-15, 0.0);
  g_assert_true (ncm_powspec_filter_peek_powspec (psf) == ps);
  g_assert_cmpint (ncm_powspec_filter_get_filter_type (psf), ==, NCM_POWSPEC_FILTER_TYPE_TOPHAT);

  ncm_powspec_filter_free (psf);
  ncm_powspec_free (ps);

  ps  = _test_power_law_new (-1.5, 0.0);
  psf = _test_filter_new (ps, NCM_POWSPEC_FILTER_TYPE_GAUSS, 1, 1.0e-3);
  ncm_assert_cmpdouble_e (ncm_powspec_filter_volume_rm3 (psf), ==, pow (2.0 * M_PI, 1.5), 1.0e-15, 0.0);

  ncm_powspec_filter_free (psf);
  ncm_powspec_free (ps);
}

/*
 * require_nderivs never lowers the order; set_nderivs does, and the kept orders
 * evaluate within reltol after preparing again: the calibration also follows the
 * derivatives, so dropping orders can end at another size (measured 5.2e-9).
 */
static void
test_ncm_powspec_filter_nderivs (void)
{
  NcmPowspec *ps        = _test_power_law_new (-1.5, 0.0);
  NcmPowspecFilter *psf = ncm_powspec_filter_new (ps, NCM_POWSPEC_FILTER_TYPE_TOPHAT);
  const gdouble lnr     = log (5.0);
  gdouble var_before;
  guint nderivs;

  g_assert_cmpuint (ncm_powspec_filter_get_nderivs (psf), ==, 1);
  ncm_powspec_filter_require_nderivs (psf, 3);
  g_assert_cmpuint (ncm_powspec_filter_get_nderivs (psf), ==, 3);
  ncm_powspec_filter_require_nderivs (psf, 2);
  g_assert_cmpuint (ncm_powspec_filter_get_nderivs (psf), ==, 3);

  ncm_powspec_filter_set_reltol (psf, 1.0e-6);
  ncm_powspec_filter_prepare (psf, NULL);
  var_before = ncm_powspec_filter_eval_dnvar_dlnrn (psf, 0.0, lnr, 0);

  ncm_powspec_filter_set_nderivs (psf, 1);
  g_assert_cmpuint (ncm_powspec_filter_get_nderivs (psf), ==, 1);
  ncm_powspec_filter_prepare (psf, NULL);

  ncm_assert_cmpdouble_e (ncm_powspec_filter_eval_dnvar_dlnrn (psf, 0.0, lnr, 0), ==, var_before, 1.0e-6, 0.0);
  ncm_assert_cmpdouble_e (ncm_powspec_filter_eval_dnlnvar_dlnrn (psf, 0.0, lnr, 1), ==, -1.5, 0.0, 3.0e-7);

  g_object_set (psf, "nderivs", 2, NULL);
  g_object_get (psf, "nderivs", &nderivs, NULL);
  g_assert_cmpuint (nderivs, ==, 2);

  ncm_powspec_filter_free (psf);
  ncm_powspec_free (ps);
}

/*
 * The setters read back; require_zi/zf only widen the range; a new reltol recalibrates
 * to the grid a fresh filter at that reltol has; no grid before the first prepare.
 */
static void
test_ncm_powspec_filter_settings (void)
{
  NcmPowspec *ps        = _test_power_law_new (-1.5, 0.0);
  NcmPowspecFilter *psf = ncm_powspec_filter_new (ps, NCM_POWSPEC_FILTER_TYPE_TOPHAT);
  NcmPowspecFilter *fresh;
  guint N_k, N_z, N_k_fresh, N_z_fresh, max_k_knots;
  gdouble zi, zf;

  ncm_powspec_filter_get_nknots (psf, &N_k, &N_z);
  g_assert_cmpuint (N_k, ==, 0);
  g_assert_cmpuint (N_z, ==, 0);

  ncm_powspec_filter_set_max_k_knots (psf, 12345);
  ncm_powspec_filter_set_max_z_knots (psf, 321);
  g_assert_cmpuint (ncm_powspec_filter_get_max_k_knots (psf), ==, 12345);
  g_assert_cmpuint (ncm_powspec_filter_get_max_z_knots (psf), ==, 321);
  g_object_get (psf, "max-k-knots", &max_k_knots, NULL);
  g_assert_cmpuint (max_k_knots, ==, 12345);

  ncm_powspec_filter_set_reltol (psf, 1.0e-8);
  ncm_powspec_filter_set_reltol_z (psf, 1.0e-7);
  g_assert_cmpfloat (ncm_powspec_filter_get_reltol (psf), ==, 1.0e-8);
  g_assert_cmpfloat (ncm_powspec_filter_get_reltol_z (psf), ==, 1.0e-7);

  ncm_powspec_filter_set_zi (psf, 0.5);
  ncm_powspec_filter_set_zf (psf, 0.8);
  ncm_powspec_filter_require_zi (psf, 0.1);
  ncm_powspec_filter_require_zf (psf, 0.9);
  ncm_powspec_filter_require_zi (psf, 0.7);
  ncm_powspec_filter_require_zf (psf, 0.6);
  g_object_get (psf, "zi", &zi, "zf", &zf, NULL);
  g_assert_cmpfloat (zi, ==, 0.1);
  g_assert_cmpfloat (zf, ==, 0.9);

  ncm_powspec_filter_set_reltol (psf, 1.0e-2);
  ncm_powspec_filter_prepare (psf, NULL);
  ncm_powspec_filter_get_nknots (psf, &N_k, &N_z);

  ncm_powspec_filter_set_reltol (psf, 1.0e-8);
  ncm_powspec_filter_prepare_if_needed (psf, NULL);
  ncm_powspec_filter_get_nknots (psf, &N_k_fresh, &N_z_fresh);
  g_assert_cmpuint (N_k_fresh, >, N_k);

  fresh = ncm_powspec_filter_new (ps, NCM_POWSPEC_FILTER_TYPE_TOPHAT);
  ncm_powspec_filter_set_max_k_knots (fresh, 12345);
  ncm_powspec_filter_set_max_z_knots (fresh, 321);
  ncm_powspec_filter_set_reltol (fresh, 1.0e-8);
  ncm_powspec_filter_set_reltol_z (fresh, 1.0e-7);
  ncm_powspec_filter_set_zi (fresh, 0.1);
  ncm_powspec_filter_set_zf (fresh, 0.9);
  ncm_powspec_filter_prepare (fresh, NULL);
  ncm_powspec_filter_get_nknots (fresh, &N_k, &N_z);
  g_assert_cmpuint (N_k, ==, N_k_fresh);
  g_assert_cmpuint (N_z, ==, N_z_fresh);
  /* The same grid; the transforms may round differently (measured 1.7e-15). */
  ncm_assert_cmpdouble_e (ncm_powspec_filter_eval_var (fresh, 0.3, 8.0), ==, ncm_powspec_filter_eval_var (psf, 0.3, 8.0), 1.0e-13, 0.0);

  ncm_powspec_filter_free (fresh);
  ncm_powspec_filter_free (psf);
  ncm_powspec_free (ps);
}

/*
 * Before a prepare the r range is ln r0 -+ L/2; after it, the grid's own ends, within
 * one knot of that estimate and spanning (N_k - 1) L / N_k.
 */
static void
test_ncm_powspec_filter_r_range (void)
{
  NcmPowspec *ps        = _test_power_law_new (-1.5, 0.0);
  NcmPowspecFilter *psf = _test_filter_new (ps, NCM_POWSPEC_FILTER_TYPE_TOPHAT, 1, 1.0e-6);
  const gdouble L       = log (ncm_powspec_get_kmax (ps) / ncm_powspec_get_kmin (ps));
  const gdouble r_min   = ncm_powspec_filter_get_r_min (psf);
  const gdouble r_max   = ncm_powspec_filter_get_r_max (psf);
  gdouble lnr0;
  guint N_k, N_z;

  ncm_powspec_filter_get_nknots (psf, &N_k, &N_z);
  ncm_assert_cmpdouble_e (log (r_max / r_min), ==, (N_k - 1.0) * L / N_k, 1.0e-12, 0.0);

  g_object_get (psf, "lnr0", &lnr0, NULL);
  ncm_powspec_filter_set_lnr0 (psf, lnr0 + log (10.0));
  ncm_assert_cmpdouble_e (log (ncm_powspec_filter_get_r_min (psf)), ==, lnr0 + log (10.0) - 0.5 * L, 1.0e-14, 0.0);
  ncm_assert_cmpdouble_e (log (ncm_powspec_filter_get_r_max (psf)), ==, lnr0 + log (10.0) + 0.5 * L, 1.0e-14, 0.0);

  ncm_powspec_filter_prepare (psf, NULL);
  ncm_powspec_filter_get_nknots (psf, &N_k, &N_z);
  g_assert_cmpfloat (fabs (log (ncm_powspec_filter_get_r_min (psf) / (10.0 * r_min))), <, L / N_k);
  g_assert_cmpfloat (fabs (log (ncm_powspec_filter_get_r_max (psf) / (10.0 * r_max))), <, L / N_k);

  ncm_powspec_filter_set_best_lnr0 (psf);
  ncm_powspec_filter_prepare (psf, NULL);
  ncm_assert_cmpdouble_e (ncm_powspec_filter_get_r_min (psf), ==, r_min, 1.0e-12, 0.0);
  ncm_assert_cmpdouble_e (ncm_powspec_filter_get_r_max (psf), ==, r_max, 1.0e-12, 0.0);

  ncm_powspec_filter_free (psf);
  ncm_powspec_free (ps);
}

/* A table starting above z = 0: the calibration runs at zi. */
static void
test_ncm_powspec_filter_zi (void)
{
  NcmPowspec *ps        = _test_power_law_new (-1.5, 0.5);
  NcmPowspecFilter *psf = _test_filter_new (ps, NCM_POWSPEC_FILTER_TYPE_TOPHAT, 1, 1.0e-6);

  ncm_assert_cmpdouble_e (ncm_powspec_filter_eval_dnlnvar_dlnrn (psf, 0.7, log (5.0), 1), ==, -1.5, 0.0, 3.0e-7);

  ncm_powspec_filter_free (psf);
  ncm_powspec_free (ps);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/powspec_filter/power_law", &test_ncm_powspec_filter_power_law);
  g_test_add_func ("/ncm/powspec_filter/gauss_closed_form", &test_ncm_powspec_filter_gauss_closed_form);
  g_test_add_func ("/ncm/powspec_filter/arb", &test_ncm_powspec_filter_arb);
  g_test_add_func ("/ncm/powspec_filter/eval_forms", &test_ncm_powspec_filter_eval_forms);
  g_test_add_func ("/ncm/powspec_filter/nderivs", &test_ncm_powspec_filter_nderivs);
  g_test_add_func ("/ncm/powspec_filter/settings", &test_ncm_powspec_filter_settings);
  g_test_add_func ("/ncm/powspec_filter/r_range", &test_ncm_powspec_filter_r_range);
  g_test_add_func ("/ncm/powspec_filter/zi", &test_ncm_powspec_filter_zi);

  g_test_run ();
}

