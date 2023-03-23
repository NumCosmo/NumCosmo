/***************************************************************************
 *            test_nc_galaxy_acf.c
 *
 *  Fri May 11 21:18:21 2012
 *  Copyright  2012 Fernando de Simoni & Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Fernando de Simoni & Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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

typedef struct _TestNcGalaxyAcf
{
  NcGalaxyAcf *acf;
} TestNcGalaxyAcf;

void test_nc_galaxy_acf_new (TestNcGalaxyAcf *test, gconstpointer pdata);
void test_nc_galaxy_acf_free (TestNcGalaxyAcf *test, gconstpointer pdata);

void test_nc_galaxy_acf_traps (TestNcGalaxyAcf *test, gconstpointer pdata);
void test_nc_galaxy_acf_invalid_st (TestNcGalaxyAcf *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();
  
  g_test_set_nonfatal_assertions ();
  
  g_test_add ("/nc/galaxy/acf/traps", TestNcGalaxyAcf, NULL,
              &test_nc_galaxy_acf_new,
              &test_nc_galaxy_acf_traps,
              &test_nc_galaxy_acf_free);
#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_add ("/nc/galaxy/acf/invalid/st/subprocess", TestNcGalaxyAcf, NULL,
              &test_nc_galaxy_acf_new,
              &test_nc_galaxy_acf_invalid_st,
              &test_nc_galaxy_acf_free);
#endif
  g_test_run ();
}

void
test_nc_galaxy_acf_new (TestNcGalaxyAcf *test, gconstpointer pdata)
{
}

void
test_nc_galaxy_acf_free (TestNcGalaxyAcf *test, gconstpointer pdata)
{
}

void
test_nc_galaxy_acf_traps (TestNcGalaxyAcf *test, gconstpointer pdata)
{
#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_trap_subprocess ("/nc/galaxy/acf/invalid/st/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}

void
test_nc_galaxy_acf_invalid_st (TestNcGalaxyAcf *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}

