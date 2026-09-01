/***************************************************************************
 *            test_ncm_stats_dist_common.h
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
 * The NcmStatsDist checks fuse two claims. Building an estimator -- adding observations,
 * preparing, fitting the interpolation with its cross-validation, evaluating, sampling,
 * serializing -- is mechanics that a small sample exercises completely. Requiring the
 * estimator to reproduce the distribution it was built from is statistics: a
 * Jensen-Shannon divergence over hundreds of draws, a covariance recovered to 50%, a
 * sampling frequency matched over ten million draws.
 *
 * The mechanics are what the instrumented lane needs and they are cheap at a small
 * sample; the statistical claims need the full sample and cost, in one binary, three
 * quarters of its runtime. So the same checks run in two modes, chosen by main().
 */

#ifndef _TEST_NCM_STATS_DIST_COMMON_H
#define _TEST_NCM_STATS_DIST_COMMON_H

#include <glib.h>
#include <numcosmo/numcosmo.h>

G_BEGIN_DECLS

/**
 * TestNcmStatsDistMode:
 * @TEST_NCM_STATS_DIST_MECHANICS: small samples, no statistical assertion
 * @TEST_NCM_STATS_DIST_DIVERGENCE: full samples, distribution recovery asserted
 */
typedef enum _TestNcmStatsDistMode
{
  TEST_NCM_STATS_DIST_MECHANICS,
  TEST_NCM_STATS_DIST_DIVERGENCE,
} TestNcmStatsDistMode;

gint test_ncm_stats_dist_main (gint argc, gchar *argv[], TestNcmStatsDistMode mode);

G_END_DECLS

#endif /* _TEST_NCM_STATS_DIST_COMMON_H */

