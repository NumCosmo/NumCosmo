/***************************************************************************
 *            test_ncm_fit_esmcmc_common.h
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
 * The ESMCMC checks come in two halves that were fused in one check each.
 *
 * Driving a chain -- start_run(), run(), trim_by_type(), validate(), end_run() -- is
 * mechanics, and a short chain exercises all of it. Requiring the sample covariance to
 * reproduce the data covariance is a statistical claim, converging as 1/sqrt(n), and the
 * checks pay for it with a loop that doubles the chain length up to fifteen times. That
 * makes their cost unbounded and their outcome noisy, neither of which belongs in an
 * instrumented lane that has to stay fast and deterministic.
 *
 * So the same checks run in two modes. %TEST_NCM_FIT_ESMCMC_MECHANICS runs a short fixed
 * chain and asserts the mechanics; %TEST_NCM_FIT_ESMCMC_RECOVERY runs the retry loop and
 * asserts the recovery. The mode is chosen once by each executable's main().
 */

#ifndef _TEST_NCM_FIT_ESMCMC_COMMON_H
#define _TEST_NCM_FIT_ESMCMC_COMMON_H

#include <glib.h>
#include <numcosmo/numcosmo.h>

G_BEGIN_DECLS

/**
 * TestNcmFitESMCMCMode:
 * @TEST_NCM_FIT_ESMCMC_MECHANICS: short fixed chains, no statistical assertion
 * @TEST_NCM_FIT_ESMCMC_RECOVERY: full chains with retries, covariance recovery asserted
 */
typedef enum _TestNcmFitESMCMCMode
{
  TEST_NCM_FIT_ESMCMC_MECHANICS,
  TEST_NCM_FIT_ESMCMC_RECOVERY,
} TestNcmFitESMCMCMode;

gint test_ncm_fit_esmcmc_main (gint argc, gchar *argv[], TestNcmFitESMCMCMode mode);

G_END_DECLS

#endif /* _TEST_NCM_FIT_ESMCMC_COMMON_H */

