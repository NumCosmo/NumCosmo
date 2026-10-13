/***************************************************************************
 *            autoknots_stress.c
 *
 *  Sat September 26 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * autoknots_stress.c
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

/*
 * Monte Carlo stress test of the knot placement of NcmSplineFunc on random functions, the
 * tests of Vitenti et al. (2025), https://doi.org/10.1016/j.ascom.2025.100970. Each
 * realization draws the parameters of the base function (--type) from --pdf with both
 * columns of the parameter table set to --p1 and --p2, places the knots with --ftype and
 * compares the spline with the function on a grid of --ngrid points; see
 * ncm_spline_func_test.c for the base functions and the statistics.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include "ncm_spline_func_test.h"

gint
main (gint argc, gchar *argv[])
{
  gchar *type_str        = g_strdup ("polynomial");
  gchar *ftype_str       = g_strdup ("function-spline");
  gchar *pdf_str         = g_strdup ("normal");
  gchar *output          = NULL;
  gint npar              = 7;
  gdouble p1             = 0.0;
  gdouble p2             = 1.0;
  gint nsim              = 100;
  gint ngrid             = 10000;
  gint64 seed            = 1;
  gdouble rel_error      = 1.0e-8;
  gdouble scale          = 0.0;
  gdouble xi             = 0.0;
  gdouble xf             = 1.0;
  GError *error          = NULL;
  GOptionEntry entries[] = {
    {"type",      't', 0, G_OPTION_ARG_STRING, &type_str,  "Base function: polynomial, polynomial-pos, cosine or exp-sinc", NULL},
    {"ftype",     'f', 0, G_OPTION_ARG_STRING, &ftype_str, "Knot placement, a NcmSplineFuncType nick", NULL},
    {"pdf",       'p', 0, G_OPTION_ARG_STRING, &pdf_str,   "Parameter distribution: flat or normal", NULL},
    {"npar",      'n', 0, G_OPTION_ARG_INT,    &npar,      "Number of parameters of the base function", NULL},
    {"p1",          0, 0, G_OPTION_ARG_DOUBLE, &p1,        "Lower limit (flat) or mean (normal) of every parameter", NULL},
    {"p2",          0, 0, G_OPTION_ARG_DOUBLE, &p2,        "Upper limit (flat) or standard deviation (normal) of every parameter", NULL},
    {"nsim",      'N', 0, G_OPTION_ARG_INT,    &nsim,      "Number of realizations", NULL},
    {"ngrid",     'g', 0, G_OPTION_ARG_INT,    &ngrid,     "Number of comparison points", NULL},
    {"seed",      's', 0, G_OPTION_ARG_INT64,  &seed,      "Random seed", NULL},
    {"rel-error", 'r', 0, G_OPTION_ARG_DOUBLE, &rel_error, "Relative tolerance", NULL},
    {"scale",       0, 0, G_OPTION_ARG_DOUBLE, &scale,     "Scale of the function values", NULL},
    {"xi",          0, 0, G_OPTION_ARG_DOUBLE, &xi,        "Lower limit", NULL},
    {"xf",          0, 0, G_OPTION_ARG_DOUBLE, &xf,        "Upper limit", NULL},
    {"output",    'o', 0, G_OPTION_ARG_STRING, &output,    "File for the statistics of every realization", NULL},
    { NULL, 0, 0, 0, NULL, NULL, NULL }
  };
  GOptionContext *context = g_option_context_new ("- stress test of NcmSplineFunc knot placement");
  NcmSplineFuncTest *sft;
  const GEnumValue *type_v, *ftype_v, *pdf_v;

  ncm_cfg_init ();

  g_option_context_add_main_entries (context, entries, NULL);

  if (!g_option_context_parse (context, &argc, &argv, &error))
    g_error ("autoknots_stress: %s", error->message);

  type_v  = ncm_cfg_get_enum_by_id_name_nick (NCM_TYPE_SPLINE_FUNC_TEST_TYPE, type_str);
  ftype_v = ncm_cfg_get_enum_by_id_name_nick (NCM_TYPE_SPLINE_FUNC_TYPE, ftype_str);
  pdf_v   = ncm_cfg_get_enum_by_id_name_nick (NCM_TYPE_SPLINE_FUNC_TEST_TYPE_PDF, pdf_str);

  if ((type_v == NULL) || (type_v->value == NCM_SPLINE_FUNC_TEST_TYPE_USER))
    g_error ("autoknots_stress: unknown base function `%s'.", type_str);

  if (ftype_v == NULL)
    g_error ("autoknots_stress: unknown knot placement `%s'.", ftype_str);

  if (pdf_v == NULL)
    g_error ("autoknots_stress: unknown parameter distribution `%s'.", pdf_str);

  if ((npar <= 0) || (nsim <= 0) || (ngrid <= 0))
    g_error ("autoknots_stress: npar, nsim and ngrid must be positive.");

  sft = ncm_spline_func_test_new ();

  ncm_spline_func_test_set_type (sft, type_v->value);
  ncm_spline_func_test_set_seed (sft, seed);
  ncm_spline_func_test_set_ngrid (sft, ngrid);
  ncm_spline_func_test_set_rel_error (sft, rel_error);
  ncm_spline_func_test_set_scale (sft, scale);
  ncm_spline_func_test_set_xi (sft, xi);
  ncm_spline_func_test_set_xf (sft, xf);
  ncm_spline_func_test_set_params_info_all (sft, npar, p1, p2);

  ncm_spline_func_test_prepare (sft, ftype_v->value, pdf_v->value);

  if (output != NULL)
    ncm_spline_func_test_monte_carlo_and_save_to_txt (sft, nsim, output);
  else
    ncm_spline_func_test_monte_carlo (sft, nsim);

  ncm_spline_func_test_log_vals_mc_stats (sft, NULL);

  ncm_spline_func_test_unref (sft);
  g_option_context_free (context);
  g_free (type_str);
  g_free (ftype_str);
  g_free (pdf_str);
  g_free (output);

  return 0;
}

