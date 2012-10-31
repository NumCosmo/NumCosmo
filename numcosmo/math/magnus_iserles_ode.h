/***************************************************************************
 *            magnus_iserles_ode.h
 *
 *  Thu Jun 18 14:01:28 2009
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

#ifndef _NC_MAGNUS_ISERLES_ODE_H
#define _NC_MAGNUS_ISERLES_ODE_H

#include <glib.h>
#include <glib-object.h>
#include <gsl/gsl_matrix_long_double.h>

G_BEGIN_DECLS

typedef struct _NcmMIOdeFunction NcmMIOdeFunction;

/**
 * NcmMIOdeFunction:
 *
 * FIXME
 */
struct _NcmMIOdeFunction
{
  /*< private >*/
  long double (*func) (long double x, gpointer params);
  gpointer params;
};

typedef struct _NcmMIOde NcmMIOde;

/**
 * NcmMIOde:
 *
 * FIXME
 */
struct _NcmMIOde
{
  /*< private >*/
  long double x0;
  long double xi;
  long double h;
  long double g_x0;
  long double psi;
  gsl_matrix_long_double *A;
  gsl_matrix_long_double *exp_A;
  gsl_matrix_long_double *Omega;
  gsl_matrix_long_double *exp_Omega;
  gsl_matrix_long_double *U;
  NcmMIOdeFunction *g;
};

NcmMIOde *ncm_magnus_iserles_ode_new (long double x0, NcmMIOdeFunction *g);
gboolean ncm_magnus_iserles_ode_step (NcmMIOde *mi_ode, long double h);
gboolean ncm_magnus_iserles_ode_step_frac (NcmMIOde *mi_ode, long double frac);
void ncm_magnus_iserles_ode_eval_vec (NcmMIOde *mi_ode, long double u_i, long double up_i, long double *u, long double *up);

G_END_DECLS

#endif /* _NC_MAGNUS_ISERLES_ODE_H */


