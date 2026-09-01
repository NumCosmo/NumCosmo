/***************************************************************************
 *            test_ncm_fit.c
 *
 *  Sun February 04 16:02:57 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2018 <vitenti@uel.br>
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
 * gsl:ls. Every check lives in test_ncm_fit_common.c; this file names the
 * algorithm and nothing else.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>
#include "test_ncm_fit_common.h"

static const TestNcmFitAlgo algo = {
  .lib       = "gsl",
  .algo      = "ls",
  .fit_type  = NCM_FIT_TYPE_GSL_LS,
  .algo_str  = NULL,
  .fit_gtype = ncm_fit_gsl_ls_get_type,
  .max_dim   = TEST_NCM_FIT_DIM,
  .max_iter  = 10000000,
};

gint
main (gint argc, gchar *argv[])
{
  return test_ncm_fit_main (argc, argv, &algo);
}

