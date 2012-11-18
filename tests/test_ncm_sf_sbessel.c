/***************************************************************************
 *            test_ncm_sf_sbessel.c
 *
 *  Tue July 03 13:35:29 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#include <gsl/gsl_sf_bessel.h>

#define NTOT 10000
#define XMAX 8.0
#define L 1200

gint
main (gint argc, gchar *argv[])
{
  if (FALSE)
  {
	GTimer *bench = g_timer_new ();
	gdouble total = 0.0;
	gdouble time_elap[NTOT];
	guint i, j;

	memset (time_elap, 0, sizeof (gdouble) * NTOT);

	g_test_init (&argc, &argv, NULL);
	ncm_cfg_init ();
	//ncm_cfg_enable_gsl_err_handler ();

	for (j = 430; j <= L; j++)
	{
	  printf ("# L = %u\n", j);
	  for (i = 0; i < NTOT; i++)
	  {
		const gdouble x = pow (10.0, -XMAX * 0.5 + XMAX / (NTOT - 1.0) * i);
		g_timer_start (bench);
		{
		  const gdouble ncm_jl = ncm_sf_sbessel (j, x);
		  const gdouble gsl_jl = gsl_sf_bessel_jl (j, x);
		  const gdouble err = fabs ((ncm_jl - gsl_jl) / ncm_jl);
		  if (ncm_jl > 1.0e-250 && gsl_jl > 1.0e-250)
			total = GSL_MAX (err, total);
		  time_elap[i] = (j * time_elap[i] + g_timer_elapsed (bench, NULL)) / (j + 1.0);

		  if (total > 1e-6)
		  {
			printf ("%u % 20.15g % 20.15g % 20.15g %e %e\n", j, x, ncm_jl, gsl_jl, err, total);
			exit (0);
		  }
		}
	  }
	}
	printf ("# TOTAL % 20.15g\n", total);

	for (i = 0; i < NTOT && FALSE; i++)
	{
	  const gdouble x = XMAX / (NTOT - 1.0) * i;
	  printf ("% 20.15g %e\n", x, time_elap[i]);
	}
  }
}

