/***************************************************************************
 *            test_nc_galaxy_acf.c
 *
 *  Fri May 11 21:18:21 2012
 *  Copyright  2012 Fernando de Simoni & Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Fernando de Simoni & Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init ();
  //ncm_cfg_enable_gsl_err_handler ();

  if (FALSE)
  {
	NcHICosmo *xcdm = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
	NcDistance *dist = nc_distance_new (1.2);
	NcGrowthFunc *gf = nc_growth_func_new ();
	NcTransferFunc *tf = nc_transfer_func_eh_new ();
	NcGalaxyAcf *acf = nc_galaxy_acf_new (gf, dist, tf);
	//guint l = 1;
	gint i;

	nc_distance_prepare (dist, xcdm);
	nc_growth_func_prepare (gf, NC_HICOSMO (xcdm));
	nc_transfer_func_prepare (tf, NC_HICOSMO (xcdm));
//ncm_model_params_log_all (NCM_MODEL (xcdm));

//	printf ("%u\n", ncm_vector_len (acf->s->xv));

	for (i = 0; i < 100; i++)
	{
	  ncm_galaxy_acf_prepare_psi (acf, xcdm, i);
	}
  }

  return 0;
}
