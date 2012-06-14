/***************************************************************************
 *            covariance.c
 *
 *  Sun Oct 18 14:35:28 2009
 *  Copyright  2009  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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
 
/**
 * SECTION:covariance
 * @title: Perturbation Covariance
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <glib.h>
#include <gsl/gsl_math.h>

NcLinearPert *gpert = NULL;
/*
static void
cov_direct (const int *ndim, const double x[], const int *ncomp, double f[])
{
  gint i;
  gpert->k = 1000.0 * x[0];
  printf ("# Go %f\n", gpert->k);
//  nc_pert_set_init_adiabatic (gpert);
//  nc_pert_linear_evolve_to_end (gpert);

  for (i = 0; i <= gpert->lmax; i++)
  {
//    f[i] = gsl_pow_2 (nc_pert_linear_get_theta (gpert, i)) / (gpert->k);
  }
}
*/
void
nc_pert_cov_direct (NcLinearPert *pert)
{
//  gint nregions, neval, fail;
//  gsl_vector *vintegral = gsl_vector_alloc (pert->lmax + 1);
//  gsl_vector *verror = gsl_vector_alloc (pert->lmax + 1);
//  gsl_vector *vprob = gsl_vector_alloc (pert->lmax + 1);
/*  
  gint i;
  
  gpert = pert;
  
  printf ("# Fail = %d\n", fail);
  for (i = 0; i <= pert->lmax; i++)
    printf ("%d %.15g\n", i, gsl_vector_get (vintegral, i));

  */
}
