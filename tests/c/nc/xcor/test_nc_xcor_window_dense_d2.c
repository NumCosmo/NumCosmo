/***************************************************************************
 *            test_nc_xcor_window_dense_d2.c
 *
 *  Tue Oct 13 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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
 * The Chebyshev closures against the dense certified windows at Bessel-derivative
 * order 2; the checks are in test_nc_xcor_window_common.c.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include "test_nc_xcor_window_common.h"

gint
main (gint argc, gchar *argv[])
{
  NcmObjDictStr *dense;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();
  g_test_set_nonfatal_assertions ();

  dense = test_nc_xcor_window_dense_table ();

  test_nc_xcor_window_register (dense, NULL, FALSE, 2, "/nc/xcor/window_dense");

  g_test_run ();

  ncm_obj_dict_str_unref (dense);
}

