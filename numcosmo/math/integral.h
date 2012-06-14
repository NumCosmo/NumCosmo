/***************************************************************************
 *            integral.h
 *
 *  Wed Aug 13 20:38:13 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
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

#ifndef _NC_INTEGRAL_H
#define _NC_INTEGRAL_H

#include <glib.h>

G_BEGIN_DECLS

gsl_integration_workspace **nc_integral_get_workspace ();

typedef struct _NcIntegrand2dim NcIntegrand2dim;

/**
 * NcIntegrand2dim:
 *
 * FIXME
 */
struct _NcIntegrand2dim
{
  /*< private >*/
	gpointer userdata;
	gdouble (*f) (gdouble x, gdouble y, gpointer userdata);
};

typedef struct _NcIntegralFixed NcIntegralFixed;

/**
 * NcIntegralFixed:
 *
 * FIXME
 */
struct _NcIntegralFixed
{
  /*< private >*/
  gdouble xl;
  gdouble xu;
  gdouble *int_nodes;
  gulong n_nodes;
  gulong rule_n;
#ifdef HAVE_GSL_GLF
  gsl_integration_glfixed_table *glt;
#endif /* HAVE_GSL_GLF */
};

gint nc_integral_locked_a_b (gsl_function *F, gdouble a, gdouble b, gdouble abstol, gdouble reltol, gdouble *result, gdouble *error);
gint nc_integral_locked_a_inf (gsl_function *F, gdouble a, gdouble abstol, gdouble reltol, gdouble *result, gdouble *error);

gint nc_integral_cached_0_x (NcFunctionCache *cache, gsl_function *F, gdouble x, gdouble *result, gdouble *error);
gint nc_integral_cached_x_inf (NcFunctionCache *cache, gsl_function *F, gdouble x, gdouble *result, gdouble *error);

gboolean ncm_integrate_2dim (NcIntegrand2dim *integ, gdouble xi, gdouble yi, gdouble xf, gdouble yf, gdouble epsrel, gdouble epsabs, gdouble *result, gdouble *error);

NcIntegralFixed *nc_integral_fixed_new (gulong n_nodes, gulong rule_n, gdouble xl, gdouble xu);
void nc_integral_fixed_free (NcIntegralFixed *intf);
void nc_integral_fixed_calc_nodes (NcIntegralFixed *intf, gsl_function *F);
gdouble nc_integral_fixed_nodes_eval (NcIntegralFixed *intf);
gdouble nc_integral_fixed_integ_mult (NcIntegralFixed *intf, gsl_function *F);
gdouble nc_integral_fixed_integ_posdef_mult (NcIntegralFixed *intf, gsl_function *F, gdouble max, gdouble reltol);

G_END_DECLS

#endif /* _NC_INTEGRAL_H */
