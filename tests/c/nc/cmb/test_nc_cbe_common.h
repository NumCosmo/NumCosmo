/***************************************************************************
 *            test_nc_cbe.c
 *
 *  Thu January 05 19:23:54 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2017 <vitenti@uel.br>
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
 * Models and checks shared by the NcCBE test executables.
 *
 * The checks divide sharply by cost: Cls and calc_ps run CLASS deep enough to
 * take 12.3 s and 4.4 s per model, while the rest are 0.4 s or less. Measured
 * per line, those two are 91% of this file's runtime for about 2% of the lines
 * only they reach, so they run as their own acceptance-tier executable and the
 * cheap checks stay in the instrumented lane. Splitting by check rather than
 * sharing a prepared state keeps each subtest's object its own -- calc_ps
 * mutates it with nc_cbe_set_calc_transfer().
 */

#ifndef _TEST_NC_CBE_COMMON_H
#define _TEST_NC_CBE_COMMON_H

#include <glib.h>
#include <numcosmo/numcosmo.h>

G_BEGIN_DECLS

typedef struct _TestNcCBE
{
  NcCBE *cbe;
  NcHICosmo *cosmo;
  guint ntests;
} TestNcCBE;

typedef struct _TestNcCBEFunc
{
  void (*func) (TestNcCBE *test, gconstpointer pdata);

  const gchar *name;
  gconstpointer pdata;
} TestNcCBEFunc;

void test_nc_cbe_free (TestNcCBE *test, gconstpointer pdata);

void test_nc_cbe_sanity (TestNcCBE *test, gconstpointer pdata);
void test_nc_cbe_compare_bg (TestNcCBE *test, gconstpointer pdata);
void test_nc_cbe_serialize (TestNcCBE *test, gconstpointer pdata);
void test_nc_cbe_calc_ps (TestNcCBE *test, gconstpointer pdata);
void test_nc_cbe_prec (TestNcCBE *test, gconstpointer pdata);
void test_nc_cbe_thermodyn (TestNcCBE *test, gconstpointer pdata);
void test_nc_cbe_Cls (TestNcCBE *test, gconstpointer pdata);
void test_nc_cbe_traps (TestNcCBE *test, gconstpointer pdata);

/**
 * test_nc_cbe_main:
 * @tests: the checks this executable runs, each against every model
 * @n_tests: how many
 * @with_traps: register the trap and its subprocess as well
 *
 * Registers @tests across every cosmology and runs them.
 */
gint test_nc_cbe_main (gint argc, gchar *argv[], const TestNcCBEFunc *tests,
                       guint n_tests, gboolean with_traps);

G_END_DECLS

#endif /* _TEST_NC_CBE_COMMON_H */

