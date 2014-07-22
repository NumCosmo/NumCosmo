/***************************************************************************
 *            test_nc_transfer_func.c
 *
 *  Wed May 16 21:32:07 2012
 *  Copyright  2012  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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

#include <glib.h>
#include <glib/gprintf.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>

void test_nc_transfer_func_new_bbks (void);
void test_nc_transfer_func_new_eh (void);
void test_nc_transfer_func_eval (void);
void test_nc_transfer_func_matter_powerspectrum (void);
void test_nc_transfer_func_free (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init ();
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/numcosmo/nc_transfer_func/bbks/new", &test_nc_transfer_func_new_bbks);
  //g_test_add_func ("/numcosmo/nc_transfer_func/bbks/eval", &test_nc_transfer_func_eval);
  //g_test_add_func ("/numcosmo/nc_transfer_func/bbks/matter_power", &test_nc_transfer_func_matter_powerspectrum);
  g_test_add_func ("/numcosmo/nc_transfer_func/bbks/free", &test_nc_transfer_func_free);

  g_test_add_func ("/numcosmo/nc_transfer_func/eh/new", &test_nc_transfer_func_new_eh);
  //g_test_add_func ("/numcosmo/nc_transfer_func/eh/eval", &test_nc_transfer_func_eval);
  //g_test_add_func ("/numcosmo/nc_transfer_func/eh/matter_power", &test_nc_transfer_func_matter_powerspectrum);
  g_test_add_func ("/numcosmo/nc_transfer_func/eh/free", &test_nc_transfer_func_free);

  g_test_run ();
}

NcTransferFunc *tf = NULL;
NcHICosmoLCDM *model = NULL;

void
test_nc_transfer_func_new_bbks (void)
{
  tf = nc_transfer_func_bbks_new ();
  g_assert (NC_IS_TRANSFER_FUNC (tf));
  g_assert (NC_IS_TRANSFER_FUNC_BBKS (tf));

  test_nc_transfer_func_free ();

#if !((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 30))
  tf = nc_transfer_func_new_from_name ("NcTransferFuncBBKS");
#else
  tf = nc_transfer_func_bbks_new ();
#endif

  model = nc_hicosmo_lcdm_new ();
  g_assert (NC_IS_TRANSFER_FUNC (tf));
  g_assert (NC_IS_TRANSFER_FUNC_BBKS (tf));  
}

void
test_nc_transfer_func_new_eh (void)
{
  tf = nc_transfer_func_eh_new ();
  model = nc_hicosmo_lcdm_new ();
  g_assert (NC_IS_TRANSFER_FUNC (tf));
  g_assert (NC_IS_TRANSFER_FUNC_EH (tf));

  test_nc_transfer_func_free ();

#if !((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 30))
  tf = nc_transfer_func_new_from_name ("NcTransferFuncEH");
#else
  tf = nc_transfer_func_eh_new ();
#endif

  g_assert (NC_IS_TRANSFER_FUNC (tf));
  g_assert (NC_IS_TRANSFER_FUNC_EH (tf));  
}

void
test_nc_transfer_func_free (void)
{
  NCM_TEST_FREE (nc_transfer_func_free, tf);
}

void 
test_nc_transfer_func_eval (void)
{
  gint i;
  gdouble tot = 0.0;

  for (i = 0; i < 100; i++)
  {
    gdouble kh = 1000.0 / 99.0 * i;
    gdouble T = nc_transfer_func_eval (tf, NC_HICOSMO (model), kh);
    tot += T;
  }
}
