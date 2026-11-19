/***************************************************************************
 *            test_nc_data_planck_lkl.c
 *
 *  Thu Nov 19 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_data_planck_lkl.c
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

#include <locale.h>
#include <glib.h>
#include <glib-object.h>

void test_nc_data_planck_lkl_decimal_comma_locale (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/nc/data_planck_lkl/decimal_comma_locale",
                   &test_nc_data_planck_lkl_decimal_comma_locale);

  g_test_run ();
}

/*
 * Builds the clik low-ell EE likelihood with LC_NUMERIC set to a locale whose decimal
 * separator is a comma (issue #484: clik crashed while parsing its data files), then
 * checks that the caller's LC_NUMERIC is unchanged. Skipped when no such locale is
 * installed.
 */
void
test_nc_data_planck_lkl_decimal_comma_locale (void)
{
  const gchar *candidates[] = {
    "fr_FR.UTF-8", "de_DE.UTF-8", "pt_BR.UTF-8", "es_ES.UTF-8", "it_IT.UTF-8", NULL
  };
  gchar *saved   = g_strdup (setlocale (LC_NUMERIC, NULL));
  gboolean found = FALSE;
  guint i;

  for (i = 0; candidates[i] != NULL; i++)
  {
    if ((setlocale (LC_NUMERIC, candidates[i]) != NULL) && (localeconv ()->decimal_point[0] == ','))
    {
      found = TRUE;
      break;
    }
  }

  if (!found)
  {
    setlocale (LC_NUMERIC, saved);
    g_free (saved);
    g_test_skip ("no locale with a decimal comma is installed");

    return;
  }

  {
    NcHIPertBoltzmannCBE *cbe = nc_hipert_boltzmann_cbe_new ();
    NcDataPlanckLKL *plik     = nc_data_planck_lkl_full_new_id (NC_DATA_PLANCK_LKL_TYPE_BASELINE_18_LOWL_EE,
                                                                NC_HIPERT_BOLTZMANN (cbe));

    g_assert_cmpuint (ncm_data_get_length (NCM_DATA (plik)), >, 0);
    g_assert_cmpstr (localeconv ()->decimal_point, ==, ",");

    g_clear_object (&plik);
    nc_hipert_boltzmann_cbe_clear (&cbe);
  }

  setlocale (LC_NUMERIC, saved);
  g_free (saved);
}

