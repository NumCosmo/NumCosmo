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
 * NcCBE -- the cheap checks; the CLASS-heavy ones are in test_nc_cbe_class.c.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include "test_nc_cbe_common.h"

static const TestNcCBEFunc tests[] = {
  { test_nc_cbe_sanity, "sanity", NULL },
  { test_nc_cbe_compare_bg, "compare_bg", NULL },
  { test_nc_cbe_serialize, "serialize", NULL },
  { test_nc_cbe_prec, "precision", NULL },
  { test_nc_cbe_thermodyn, "Thermodyn", NULL },
};

gint
main (gint argc, gchar *argv[])
{
  return test_nc_cbe_main (argc, argv, tests, G_N_ELEMENTS (tests), TRUE);
}

