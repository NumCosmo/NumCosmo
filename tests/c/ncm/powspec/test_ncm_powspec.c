/***************************************************************************
 *            test_ncm_powspec.c
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

/* NcmPowspecAnalytic at its defaults: BBKS shape, LCDM growth. */
static NcmPowspec *
_test_powspec_new (void)
{
  return NCM_POWSPEC (ncm_powspec_analytic_new (NCM_POWSPEC_ANALYTIC_SHAPE_BBKS,
                                                NCM_POWSPEC_ANALYTIC_GROWTH_LCDM));
}

/* The base-class implementations, called on a subclass that overrides them. */
static NcmPowspecClass *
_test_powspec_base_class (void)
{
  return NCM_POWSPEC_CLASS (g_type_class_peek (NCM_TYPE_POWSPEC));
}

/*
 * var_tophat_R, sigma_tophat_R, corr3d and sproj against the Arb table from
 * tests/tools/make_powspec_analytic_truth_table.py, at reltol 1e-9 over the
 * table's k range. The largest error measured is 4e-11 (var, R = 50 Mpc); sproj
 * at ell = 20 aborts on GSL roundoff from reltol 1e-10.
 */
static void
test_ncm_powspec_integrals (void)
{
  const gdouble reltol = 1.0e-9;
  NcmSerialize *ser    = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  gchar *path          = ncm_cfg_get_data_filename ("truth_tables/powspec/ncm_powspec_analytic_integrals.bin", TRUE);
  NcmObjDictStr *dict  = ncm_serialize_dict_str_from_binfile (ser, path);
  NcmVector *k_range   = NCM_VECTOR (ncm_obj_dict_str_peek (dict, "k_range"));
  NcmMatrix *var       = NCM_MATRIX (ncm_obj_dict_str_peek (dict, "var"));
  NcmMatrix *xi        = NCM_MATRIX (ncm_obj_dict_str_peek (dict, "xi"));
  NcmMatrix *sproj     = NCM_MATRIX (ncm_obj_dict_str_peek (dict, "sproj"));
  NcmPowspec *ps       = _test_powspec_new ();
  guint i;

  ncm_powspec_set_kmin (ps, ncm_vector_get (k_range, 0));
  ncm_powspec_set_kmax (ps, ncm_vector_get (k_range, 1));

  g_assert_cmpuint (ncm_matrix_nrows (var), ==, 5);
  g_assert_cmpuint (ncm_matrix_nrows (xi), ==, 3);
  g_assert_cmpuint (ncm_matrix_nrows (sproj), ==, 2);

  for (i = 0; i < ncm_matrix_nrows (var); i++)
  {
    const gdouble z     = ncm_matrix_get (var, i, 0);
    const gdouble R     = ncm_matrix_get (var, i, 1);
    const gdouble truth = ncm_matrix_get (var, i, 2);

    ncm_assert_cmpdouble_e (ncm_powspec_var_tophat_R (ps, NULL, reltol, z, R), ==, truth, reltol, 0.0);
    ncm_assert_cmpdouble_e (ncm_powspec_sigma_tophat_R (ps, NULL, reltol, z, R), ==, sqrt (truth), reltol, 0.0);
  }

  for (i = 0; i < ncm_matrix_nrows (xi); i++)
  {
    const gdouble z     = ncm_matrix_get (xi, i, 0);
    const gdouble r     = ncm_matrix_get (xi, i, 1);
    const gdouble truth = ncm_matrix_get (xi, i, 2);

    ncm_assert_cmpdouble_e (ncm_powspec_corr3d (ps, NULL, reltol, z, r), ==, truth, reltol, 0.0);
  }

  for (i = 0; i < ncm_matrix_nrows (sproj); i++)
  {
    const gint ell      = (gint) ncm_matrix_get (sproj, i, 0);
    const gdouble z1    = ncm_matrix_get (sproj, i, 1);
    const gdouble z2    = ncm_matrix_get (sproj, i, 2);
    const gdouble xi1   = ncm_matrix_get (sproj, i, 3);
    const gdouble xi2   = ncm_matrix_get (sproj, i, 4);
    const gdouble truth = ncm_matrix_get (sproj, i, 5);

    ncm_assert_cmpdouble_e (ncm_powspec_sproj (ps, NULL, reltol, ell, z1, z2, xi1, xi2), ==, truth, reltol, 0.0);
  }

  ncm_powspec_free (ps);
  ncm_obj_dict_str_unref (dict);
  ncm_serialize_free (ser);
  g_free (path);
}

/*
 * The default finite-difference derivatives against the closed form, on both z
 * stencils (one-sided below z = 2e-3). Largest errors measured: 1.1e-11 in z,
 * 8.7e-13 in k.
 */
static void
test_ncm_powspec_default_deriv (void)
{
  NcmPowspecClass *base = _test_powspec_base_class ();
  NcmPowspec *ps        = _test_powspec_new ();
  const gdouble zs[]    = {0.0, 1.0e-3, 0.5, 2.0};
  const gdouble ks[]    = {1.0e-4, 5.0e-2, 1.0};
  guint i, j;

  for (i = 0; i < G_N_ELEMENTS (zs); i++)
  {
    for (j = 0; j < G_N_ELEMENTS (ks); j++)
    {
      ncm_assert_cmpdouble_e (base->deriv_z (ps, NULL, zs[i], ks[j]), ==, ncm_powspec_deriv_z (ps, NULL, zs[i], ks[j]), 3.0e-11, 0.0);
      ncm_assert_cmpdouble_e (base->deriv_k (ps, NULL, zs[i], ks[j]), ==, ncm_powspec_deriv_k (ps, NULL, zs[i], ks[j]), 3.0e-12, 0.0);
    }
  }

  ncm_powspec_free (ps);
}

/*
 * The default spline meets NcmPowspec:reltol over z in [0, 2] and k in [1e-5, 10].
 * Measured maximum errors: 7.6e-4, 3.0e-6 and 1.4e-8 for reltol 1e-3, 1e-5, 1e-7.
 */
static void
test_ncm_powspec_default_spline_2d (void)
{
  NcmPowspecClass *base = _test_powspec_base_class ();
  NcmPowspec *ps        = _test_powspec_new ();
  const gdouble rels[]  = {1.0e-3, 1.0e-5, 1.0e-7};
  guint r;

  ncm_powspec_set_zi (ps, 0.0);
  ncm_powspec_set_zf (ps, 2.0);
  ncm_powspec_set_kmin (ps, 1.0e-5);
  ncm_powspec_set_kmax (ps, 10.0);

  for (r = 0; r < G_N_ELEMENTS (rels); r++)
  {
    NcmSpline2d *s2d;
    gdouble err = 0.0;
    guint i, j;

    ncm_powspec_set_reltol_spline (ps, rels[r]);
    s2d = base->get_spline_2d (ps, NULL);

    for (i = 0; i <= 12; i++)
    {
      const gdouble z = 2.0 * i / 12.0;

      for (j = 0; j < 400; j++)
      {
        const gdouble k = 1.0e-5 * pow (1.0e6, j / 399.0);

        err = GSL_MAX (err, fabs (ncm_spline2d_eval (s2d, z, k) / ncm_powspec_eval (ps, NULL, z, k) - 1.0));
      }
    }

    g_assert_cmpfloat (err, <, rels[r]);

    ncm_spline2d_free (s2d);
  }

  ncm_powspec_free (ps);
}

/* The default eval_vec is eval at each k. */
static void
test_ncm_powspec_default_eval_vec (void)
{
  NcmPowspecClass *base = _test_powspec_base_class ();
  NcmPowspec *ps        = _test_powspec_new ();
  NcmVector *k          = ncm_vector_new (20);
  NcmVector *Pk         = ncm_vector_new (20);
  guint i;

  for (i = 0; i < 20; i++)
    ncm_vector_set (k, i, 1.0e-4 * pow (1.0e5, i / 19.0));

  base->eval_vec (ps, NULL, 0.7, k, Pk);

  for (i = 0; i < 20; i++)
    g_assert_cmpfloat (ncm_vector_get (Pk, i), ==, ncm_powspec_eval (ps, NULL, 0.7, ncm_vector_get (k, i)));

  ncm_vector_free (k);
  ncm_vector_free (Pk);
  ncm_powspec_free (ps);
}

/* require_* move the range only outwards; the properties round-trip. */
static void
test_ncm_powspec_range (void)
{
  NcmPowspec *ps = _test_powspec_new ();

  ncm_powspec_set_zi (ps, 0.5);
  ncm_powspec_set_zf (ps, 1.0);
  ncm_powspec_set_kmin (ps, 1.0e-3);
  ncm_powspec_set_kmax (ps, 1.0);
  ncm_powspec_set_reltol_spline (ps, 1.0e-6);

  ncm_powspec_require_zi (ps, 0.7);
  ncm_powspec_require_zf (ps, 0.8);
  ncm_powspec_require_kmin (ps, 1.0e-2);
  ncm_powspec_require_kmax (ps, 0.5);

  g_assert_cmpfloat (ncm_powspec_get_zi (ps), ==, 0.5);
  g_assert_cmpfloat (ncm_powspec_get_zf (ps), ==, 1.0);
  g_assert_cmpfloat (ncm_powspec_get_kmin (ps), ==, 1.0e-3);
  g_assert_cmpfloat (ncm_powspec_get_kmax (ps), ==, 1.0);

  ncm_powspec_require_zi (ps, 0.1);
  ncm_powspec_require_zf (ps, 3.0);
  ncm_powspec_require_kmin (ps, 1.0e-5);
  ncm_powspec_require_kmax (ps, 20.0);

  g_assert_cmpfloat (ncm_powspec_get_zi (ps), ==, 0.1);
  g_assert_cmpfloat (ncm_powspec_get_zf (ps), ==, 3.0);
  g_assert_cmpfloat (ncm_powspec_get_kmin (ps), ==, 1.0e-5);
  g_assert_cmpfloat (ncm_powspec_get_kmax (ps), ==, 20.0);
  g_assert_cmpfloat (ncm_powspec_get_reltol_spline (ps), ==, 1.0e-6);

  {
    gdouble zi, zf, kmin, kmax, reltol;

    g_object_get (ps, "zi", &zi, "zf", &zf, "kmin", &kmin, "kmax", &kmax, "reltol", &reltol, NULL);
    g_assert_cmpfloat (zi, ==, 0.1);
    g_assert_cmpfloat (zf, ==, 3.0);
    g_assert_cmpfloat (kmin, ==, 1.0e-5);
    g_assert_cmpfloat (kmax, ==, 20.0);
    g_assert_cmpfloat (reltol, ==, 1.0e-6);
  }

  g_assert_nonnull (ncm_powspec_peek_model_ctrl (ps));

  ncm_powspec_free (ps);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/powspec/integrals", &test_ncm_powspec_integrals);
  g_test_add_func ("/ncm/powspec/default/deriv", &test_ncm_powspec_default_deriv);
  g_test_add_func ("/ncm/powspec/default/spline_2d", &test_ncm_powspec_default_spline_2d);
  g_test_add_func ("/ncm/powspec/default/eval_vec", &test_ncm_powspec_default_eval_vec);
  g_test_add_func ("/ncm/powspec/range", &test_ncm_powspec_range);

  g_test_run ();
}

