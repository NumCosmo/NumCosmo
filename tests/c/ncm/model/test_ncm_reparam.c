/***************************************************************************
 *            test_ncm_reparam.c
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

#define TEST_DIM 5

/* T = 1 + a random matrix with entries in [-0.3, 0.3], invertible, and v random. */
static NcmReparamLinear *
_test_reparam_new (NcmRNG *rng)
{
  NcmMatrix *T = ncm_matrix_new (TEST_DIM, TEST_DIM);
  NcmVector *v = ncm_vector_new (TEST_DIM);
  NcmReparamLinear *relin;
  guint i, j;

  for (i = 0; i < TEST_DIM; i++)
  {
    ncm_vector_set (v, i, ncm_rng_uniform_gen (rng, -1.0, 1.0));

    for (j = 0; j < TEST_DIM; j++)
      ncm_matrix_set (T, i, j, (i == j ? 1.0 : 0.0) + ncm_rng_uniform_gen (rng, -0.3, 0.3));
  }

  relin = ncm_reparam_linear_new (TEST_DIM, T, v);

  ncm_matrix_free (T);
  ncm_vector_free (v);

  return relin;
}

/* old2new gives T p + v, and new2old recovers p to rounding. */
static void
test_ncm_reparam_linear_round_trip (void)
{
  NcmRNG *rng             = ncm_rng_seeded_new (NULL, 20260928);
  NcmModel *model         = NCM_MODEL (ncm_model_mvnd_new (TEST_DIM));
  NcmReparamLinear *relin = _test_reparam_new (rng);
  NcmReparam *reparam     = NCM_REPARAM (relin);
  NcmVector *p            = ncm_model_orig_params_peek_vector (model);
  NcmVector *pn           = ncm_reparam_peek_params (reparam);
  NcmMatrix *T;
  NcmVector *v, *p0;
  guint i, j;

  g_object_get (relin, "matrix", &T, "vector", &v, NULL);
  g_assert_cmpuint (ncm_reparam_get_length (reparam), ==, TEST_DIM);

  for (i = 0; i < TEST_DIM; i++)
    ncm_vector_set (p, i, ncm_rng_uniform_gen (rng, -2.0, 2.0));

  p0 = ncm_vector_dup (p);

  ncm_reparam_old2new (reparam, model);

  for (i = 0; i < TEST_DIM; i++)
  {
    gdouble Tp = ncm_vector_get (v, i);

    for (j = 0; j < TEST_DIM; j++)
      Tp += ncm_matrix_get (T, i, j) * ncm_vector_get (p0, j);

    ncm_assert_cmpdouble_e (ncm_vector_get (pn, i), ==, Tp, 1.0e-14, 1.0e-15);
  }

  ncm_vector_set_zero (p);
  ncm_reparam_new2old (reparam, model);

  for (i = 0; i < TEST_DIM; i++)
    ncm_assert_cmpdouble_e (ncm_vector_get (p, i), ==, ncm_vector_get (p0, i), 1.0e-13, 1.0e-14);

  ncm_vector_free (p0);
  ncm_vector_free (v);
  ncm_matrix_free (T);
  ncm_reparam_free (reparam);
  ncm_model_free (model);
  ncm_rng_free (rng);
}

/* Descriptions are found by name, also after replacing one and after a serialization
 * round trip; the compatible type round-trips. */
static void
test_ncm_reparam_param_desc (void)
{
  NcmRNG *rng             = ncm_rng_seeded_new (NULL, 1);
  NcmReparamLinear *relin = _test_reparam_new (rng);
  NcmReparam *reparam     = NCM_REPARAM (relin);
  NcmSerialize *ser       = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcmReparam *dup;
  NcmSParam *sp;
  guint i;

  g_assert_null (ncm_reparam_peek_param_desc (reparam, 1));
  g_assert_null (ncm_reparam_get_param_desc (reparam, 1));

  ncm_reparam_set_param_desc_full (reparam, 1, "a", "a", -1.0, 1.0, 0.1, 0.0, 0.0, NCM_PARAM_TYPE_FREE);
  ncm_reparam_set_param_desc_full (reparam, 3, "b", "b", -1.0, 1.0, 0.1, 0.0, 0.0, NCM_PARAM_TYPE_FREE);

  g_assert_true (ncm_reparam_index_from_name (reparam, "a", &i));
  g_assert_cmpuint (i, ==, 1);
  g_assert_false (ncm_reparam_index_from_name (reparam, "c", &i));
  g_assert_cmpuint (i, ==, G_MAXUINT);

  /* Replacing the description of 1 drops its old name; the same name again is allowed. */
  ncm_reparam_set_param_desc_full (reparam, 1, "c", "c", -1.0, 1.0, 0.1, 0.0, 0.0, NCM_PARAM_TYPE_FREE);
  ncm_reparam_set_param_desc_full (reparam, 1, "c", "c", -2.0, 2.0, 0.1, 0.0, 0.0, NCM_PARAM_TYPE_FREE);
  g_assert_false (ncm_reparam_index_from_name (reparam, "a", &i));
  g_assert_true (ncm_reparam_index_from_name (reparam, "c", &i));
  g_assert_cmpuint (i, ==, 1);

  sp = ncm_reparam_get_param_desc (reparam, 1);
  g_assert_true (sp == ncm_reparam_peek_param_desc (reparam, 1));
  g_assert_cmpfloat (ncm_sparam_get_upper_bound (sp), ==, 2.0);
  ncm_sparam_free (sp);

  ncm_reparam_set_compat_type (reparam, NCM_TYPE_MODEL_MVND);
  g_assert_true (ncm_reparam_get_compat_type (reparam) == NCM_TYPE_MODEL_MVND);

  dup = NCM_REPARAM (ncm_serialize_dup_obj (ser, G_OBJECT (reparam)));
  g_assert_true (ncm_reparam_index_from_name (dup, "c", &i));
  g_assert_cmpuint (i, ==, 1);
  g_assert_true (ncm_reparam_index_from_name (dup, "b", &i));
  g_assert_cmpuint (i, ==, 3);
  g_assert_true (ncm_reparam_get_compat_type (dup) == NCM_TYPE_MODEL_MVND);

  ncm_reparam_free (dup);
  ncm_serialize_free (ser);
  ncm_reparam_free (reparam);
  ncm_rng_free (rng);
}

static void
test_ncm_reparam_duplicate_name (void)
{
  g_test_trap_subprocess ("/ncm/reparam/duplicate_name/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*already describes parameter 0*");
}

static void
test_ncm_reparam_duplicate_name_subprocess (void)
{
  NcmRNG *rng             = ncm_rng_seeded_new (NULL, 1);
  NcmReparamLinear *relin = _test_reparam_new (rng);

  ncm_reparam_set_param_desc_full (NCM_REPARAM (relin), 0, "a", "a", -1.0, 1.0, 0.1, 0.0, 0.0, NCM_PARAM_TYPE_FREE);
  ncm_reparam_set_param_desc_full (NCM_REPARAM (relin), 2, "a", "a", -1.0, 1.0, 0.1, 0.0, 0.0, NCM_PARAM_TYPE_FREE);
}

static void
test_ncm_reparam_singular (void)
{
  g_test_trap_subprocess ("/ncm/reparam/singular/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*singular*");
}

static void
test_ncm_reparam_singular_subprocess (void)
{
  NcmMatrix *T = ncm_matrix_new (2, 2);
  NcmVector *v = ncm_vector_new (2);

  ncm_matrix_set_all (T, 1.0);
  ncm_vector_set_zero (v);

  ncm_reparam_linear_new (2, T, v);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/reparam/linear/round_trip", &test_ncm_reparam_linear_round_trip);
  g_test_add_func ("/ncm/reparam/param_desc", &test_ncm_reparam_param_desc);
  g_test_add_func ("/ncm/reparam/duplicate_name", &test_ncm_reparam_duplicate_name);
  g_test_add_func ("/ncm/reparam/duplicate_name/subprocess", &test_ncm_reparam_duplicate_name_subprocess);
  g_test_add_func ("/ncm/reparam/singular", &test_ncm_reparam_singular);
  g_test_add_func ("/ncm/reparam/singular/subprocess", &test_ncm_reparam_singular_subprocess);

  g_test_run ();
}

