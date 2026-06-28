/***************************************************************************
 *            test_ncm_template.c
 *
 *  Thu May 10 15:13:31 2012
 *  Copyright  2023  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2023 <vitenti@uel.br>
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


typedef struct _TestNcmTemplate
{
  NcmTemplate *template;
  guint ntests;
} TestNcmTemplate;

void test_ncm_template_new (TestNcmTemplate *test, gconstpointer pdata);
void test_ncm_template_free (TestNcmTemplate *test, gconstpointer pdata);

void test_ncm_template_sanity (TestNcmTemplate *test, gconstpointer pdata);
void test_ncm_template_method1 (TestNcmTemplate *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/ncm/template/sanity", TestNcmTemplate, NULL,
              &test_ncm_template_new,
              &test_ncm_template_sanity,
              &test_ncm_template_free);

  g_test_add ("/ncm/template/sanity", TestNcmTemplate, NULL,
              &test_ncm_template_new,
              &test_ncm_template_method1,
              &test_ncm_template_free);

  g_test_run ();
}

void
test_ncm_template_new (TestNcmTemplate *test, gconstpointer pdata)
{
  NcmTemplate *template = ncm_template_new ();

  g_assert_true (NCM_IS_TEMPLATE (wf));

  test->template = template;
}

void
test_ncm_template_free (TestNcmTemplate *test, gconstpointer pdata)
{
  NcmTemplate *template = test->template;

  NCM_TEST_FREE (ncm_template_free, template);
}

void
test_ncm_template_sanity (TestNcmTemplate *test, gconstpointer pdata)
{
  g_assert_true (NCM_IS_TEMPLATE (test->template));

  ncm_template_ref (test->template);
  ncm_template_free (test->template);

  {
    NcmTemplate *template = ncm_template_ref (test->template);

    ncm_template_clear (&template);

    g_assert_true (template == NULL);
  }
}

void
test_ncm_template_method (TestNcmTemplate *test, gconstpointer pdata)
{
  gint arg1          = 1;
  gint ret1_expected = 1;
  gint ret1;

  ret1 = ncm_template_method1 (test->template, arg1);

  g_assert_true (ret1 == ret1_expected);
}

