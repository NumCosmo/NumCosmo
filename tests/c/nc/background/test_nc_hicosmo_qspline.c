/***************************************************************************
 *            test_nc_hicosmo_qspline.c
 *
 *  Wed Oct 21 15:10:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_hicosmo_qspline.c
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

#define TEST_NKNOTS 10
#define TEST_ZF 2.0

void test_nc_hicosmo_qspline_lp_ls (void);

static void _test_nc_hicosmo_qspline_ls_check (NcHICosmoQSpline *model, NcmMSet *mset, NcmSpline *weight, const guint p, const gdouble sigma);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/nc/hicosmo_qspline/lp_ls", &test_nc_hicosmo_qspline_lp_ls);

  g_test_run ();
}

/*
 * The residual functions lp_q2_ls and wlp_q2_ls reproduce the squares of lp_q2 and wlp_q2 for
 * p = 2 and 10, and a Gaussian prior on them is the Gaussian prior on the norm.
 */
void
test_nc_hicosmo_qspline_lp_ls (void)
{
  const guint ps[]          = {2, 10};
  NcmSpline *s              = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  NcHICosmoQSpline *qspline = nc_hicosmo_qspline_new (s, TEST_NKNOTS, TEST_ZF);
  NcmMSet *mset             = ncm_mset_new (qspline, NULL, NULL);
  NcmVector *xv             = ncm_vector_new (9);
  NcmVector *yv             = ncm_vector_new (9);
  NcmSpline *weight;
  guint i;

  for (i = 0; i < TEST_NKNOTS; i++)
    ncm_model_orig_vparam_set (NCM_MODEL (qspline), NC_HICOSMO_QSPLINE_Q, i, -0.5 + 0.4 * sin (1.3 * i) + 0.03 * i);

  for (i = 0; i < 9; i++)
  {
    const gdouble z = TEST_ZF * i / 8.0;

    ncm_vector_set (xv, i, z);
    ncm_vector_set (yv, i, 1.0 + 0.5 * z + 0.2 * z * z);
  }

  weight = NCM_SPLINE (ncm_spline_cubic_notaknot_new_full (xv, yv, TRUE));

  for (i = 0; i < G_N_ELEMENTS (ps); i++)
  {
    _test_nc_hicosmo_qspline_ls_check (qspline, mset, NULL, ps[i], 0.7);
    _test_nc_hicosmo_qspline_ls_check (qspline, mset, weight, ps[i], 0.7);
  }

  ncm_spline_free (weight);
  ncm_vector_free (xv);
  ncm_vector_free (yv);
  ncm_mset_free (mset);
  NCM_TEST_FREE (nc_hicosmo_free, NC_HICOSMO (qspline));
  ncm_spline_free (s);
}

static void
_test_nc_hicosmo_qspline_ls_check (NcHICosmoQSpline *model, NcmMSet *mset, NcmSpline *weight, const guint p, const gdouble sigma)
{
  const gdouble x      = p;
  const guint len      = nc_hicosmo_qspline_d2_lp_residuals_len (model, weight, p);
  NcmMSetFuncList *fls = ncm_mset_func_list_new_dim ((weight == NULL) ? "NcHICosmoQSpline:lp_q2_ls" : "NcHICosmoQSpline:wlp_q2_ls", G_OBJECT (weight), len);
  NcmMSetFuncList *fsc = ncm_mset_func_list_new ((weight == NULL) ? "NcHICosmoQSpline:lp_q2" : "NcHICosmoQSpline:wlp_q2", G_OBJECT (weight));
  NcmPrior *pls        = NCM_PRIOR (ncm_prior_gauss_func_new (NCM_MSET_FUNC (fls), 0.0, sigma, x));
  NcmPrior *psc        = NCM_PRIOR (ncm_prior_gauss_func_new (NCM_MSET_FUNC (fsc), 0.0, sigma, x));
  GArray *f            = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), len);
  const gdouble norm   = (weight == NULL) ? nc_hicosmo_qspline_lp_norm (model, NCM_SPLINE_CURVATURE_D2, p) : nc_hicosmo_qspline_weighted_lp_norm (model, NCM_SPLINE_CURVATURE_D2, p, weight);
  gdouble x_f          = x;
  gdouble f2           = 0.0;
  gdouble prior_f2     = 0.0;
  guint k;

  g_assert_cmpuint (ncm_mset_func_get_dim (NCM_MSET_FUNC (pls)), ==, len);
  g_array_set_size (f, len);

  ncm_mset_func_eval (NCM_MSET_FUNC (fls), mset, &x_f, &g_array_index (f, gdouble, 0));

  for (k = 0; k < len; k++)
    f2 += gsl_pow_2 (g_array_index (f, gdouble, k));

  /* The residuals reproduce N_p^2: to the 1e-10 of the unweighted norm, and to the 2e-10 of
   * the weighted norm, whose numerator and denominator are integrals to 1e-10 each. */
  ncm_assert_cmpdouble_e (f2, ==, norm * norm, (weight == NULL) ? 1.0e-10 : 2.0e-10, 0.0);

  /* The prior on the residuals is the Gaussian prior on N_p, (N_p / sigma)^2. */
  ncm_mset_func_eval (NCM_MSET_FUNC (pls), mset, NULL, &g_array_index (f, gdouble, 0));

  for (k = 0; k < len; k++)
    prior_f2 += gsl_pow_2 (g_array_index (f, gdouble, k));

  ncm_assert_cmpdouble_e (prior_f2, ==, gsl_pow_2 (ncm_mset_func_eval0 (NCM_MSET_FUNC (psc), mset)),
                          (weight == NULL) ? 1.0e-10 : 2.0e-10, 0.0);

  g_array_unref (f);
  ncm_prior_free (pls);
  ncm_prior_free (psc);
  ncm_mset_func_list_free (fls);
  ncm_mset_func_list_free (fsc);
}

