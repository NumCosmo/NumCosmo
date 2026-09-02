/***************************************************************************
 *            test_ncm_fit_common.h
 *
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
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
 * Checks shared by every NcmFit algorithm.
 *
 * The checks below are the same for all algorithms -- they take the fixture and
 * never mention which optimizer built it -- so they live here once, and each
 * algorithm gets a small executable that names itself and calls
 * test_ncm_fit_main(). One executable per algorithm keeps a failure attributable
 * without reading the whole log, gives each a timeout its own size, and lets
 * `meson test --slice` balance them by count.
 *
 * The algorithm reaches the checks as GLib fixture data: g_test_add() passes its
 * `tdata` to setup, test and teardown, so one setup function serves every
 * algorithm from the descriptor below.
 */

#ifndef _TEST_NCM_FIT_COMMON_H
#define _TEST_NCM_FIT_COMMON_H

#include <glib.h>
#include <numcosmo/numcosmo.h>

G_BEGIN_DECLS

#define TEST_NCM_FIT_FISHER_COV_RELTOL 1.0e-2
#define TEST_NCM_FIT_FISHER_BIAS_RELTOL 3.0e-1
#define TEST_NCM_FIT_DIM 10
#define TEST_NCM_FIT_DIM_SMALL 5

typedef struct _TestNcmFit
{
  NcmFit *fit;
  NcmRNG *rng;
  NcmDataGaussCovMVND *data_mvnd;
  guint ntests;
} TestNcmFit;

/**
 * TestNcmFitAlgo:
 * @lib: first path component, e.g. "gsl"
 * @algo: second path component, e.g. "mm_conjugate_fr"
 * @fit_type: which #NcmFit implementation ncm_fit_factory() should build
 * @algo_str: the algorithm name that implementation takes, or NULL
 * @fit_gtype: expected GType of the resulting fit, asserted after construction
 * @max_dim: upper bound of the randomly drawn problem dimension
 * @max_iter: iteration cap handed to the fit
 *
 * One algorithm's worth of what used to be macro arguments.
 */
typedef struct _TestNcmFitAlgo
{
  const gchar *lib;
  const gchar *algo;
  NcmFitType fit_type;
  const gchar *algo_str;

  GType (*fit_gtype) (void);

  guint max_dim;
  guint max_iter;
} TestNcmFitAlgo;

gint test_ncm_fit_main (gint argc, gchar *argv[], const TestNcmFitAlgo *algo);

G_END_DECLS

#endif /* _TEST_NCM_FIT_COMMON_H */

