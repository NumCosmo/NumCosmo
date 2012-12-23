/***************************************************************************
 *            linear_ncm_gsl_odeiv2.c
 *
 *  Thu Nov 12 22:20:48 2009
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "perturbations/linear.h"

#ifdef HAVE_GSL_ODEIV2
#include <gsl/gsl_odeiv2.h>

#include "linear_internal.h"

typedef struct _GSLOde2Data
{
  const gsl_odeiv2_step_type *s_type;
  gsl_odeiv2_step *s;
  gsl_odeiv2_step *sa;
  gsl_odeiv2_control *c;
  gsl_odeiv2_evolve *e;
  gsl_odeiv2_driver *d;
  gsl_odeiv2_system sys;
  gsl_vector *yi;
  gsl_vector *y;
  gsl_vector *abstol;
} GSLOde2Data;

#define GSLODE2_DATA(a) ((GSLOde2Data *)(a))

static gpointer ncm_gsl_odeiv2_create (NcLinearPert *pert);
static void ncm_gsl_odeiv2_init (NcLinearPert *pert);
static void ncm_gsl_odeiv2_set_opts (NcLinearPert *pert);
static void ncm_gsl_odeiv2_reset (NcLinearPert *pert);
static void ncm_gsl_odeiv2_end_tight_coupling (NcLinearPert *pert);
static gboolean ncm_gsl_odeiv2_evol_step (NcLinearPert *pert, gdouble g);
static gboolean ncm_gsl_odeiv2_evol (NcLinearPert *pert, gdouble g);
static gboolean ncm_gsl_odeiv2_update_los (NcLinearPert *pert);
static void ncm_gsl_odeiv2_get_sources (NcLinearPert *pert, gdouble *S0, gdouble *S1, gdouble *S2);
static void ncm_gsl_odeiv2_free (NcLinearPert *pert);
static void ncm_gsl_odeiv2_print_stats (NcLinearPert *pert);
static gdouble ncm_gsl_odeiv2_get_z (NcLinearPert *pert);
static gdouble ncm_gsl_odeiv2_get_phi (NcLinearPert *pert);
static gdouble ncm_gsl_odeiv2_get_c0 (NcLinearPert *pert);
static gdouble ncm_gsl_odeiv2_get_b0 (NcLinearPert *pert);
static gdouble ncm_gsl_odeiv2_get_c1 (NcLinearPert *pert);
static gdouble ncm_gsl_odeiv2_get_b1 (NcLinearPert *pert);
static gdouble ncm_gsl_odeiv2_get_theta2 (NcLinearPert *pert);
static gdouble ncm_gsl_odeiv2_get (NcLinearPert *pert, guint n);
static gdouble ncm_gsl_odeiv2_get_theta (NcLinearPert *pert, guint n);
static gdouble ncm_gsl_odeiv2_get_theta_p (NcLinearPert *pert, guint n);

static NcLinearPertOdeSolver _ncm_gsl_odeiv2_solver = {
  &ncm_gsl_odeiv2_create,
  &ncm_gsl_odeiv2_init,
  &ncm_gsl_odeiv2_set_opts,
  &ncm_gsl_odeiv2_reset,
  &ncm_gsl_odeiv2_evol_step,
  &ncm_gsl_odeiv2_evol,
  &ncm_gsl_odeiv2_update_los,
  &ncm_gsl_odeiv2_get_sources,
  &ncm_gsl_odeiv2_free,
  &ncm_gsl_odeiv2_print_stats,
  &ncm_gsl_odeiv2_get_z,
  &ncm_gsl_odeiv2_get_phi,
  &ncm_gsl_odeiv2_get_c0,
  &ncm_gsl_odeiv2_get_b0,
  &ncm_gsl_odeiv2_get_c1,
  &ncm_gsl_odeiv2_get_b1,
  &ncm_gsl_odeiv2_get,
  &ncm_gsl_odeiv2_get_theta,
  &ncm_gsl_odeiv2_get_theta_p,
  NULL
};
NcLinearPertOdeSolver *ncm_gsl_odeiv2_solver = &_ncm_gsl_odeiv2_solver;

static gint ncm_gsl_odeiv2_step (gdouble lambda, const gdouble y[], gdouble ydot[], gpointer params);
static gint ncm_gsl_odeiv2_band_J (gdouble lambda, const gdouble y[], gdouble *J, gdouble dfdt[], gpointer user_data);

static gpointer
ncm_gsl_odeiv2_create (NcLinearPert *pert)
{
  GSLOde2Data *data = g_slice_new (GSLOde2Data);

  data->yi = gsl_vector_alloc (pert->sys_size);
  data->y = gsl_vector_alloc (pert->sys_size);
  data->abstol = gsl_vector_alloc (pert->sys_size);

//  data->s_type = gsl_odeiv2_step_rk2;
//  data->s_type = gsl_odeiv2_step_rk4;
  data->s_type = gsl_odeiv2_step_rkf45;
//  data->s_type = gsl_odeiv2_step_rkck;
//  data->s_type = gsl_odeiv2_step_rk8pd;
//  data->s_type = gsl_odeiv2_step_rk2imp;
//  data->s_type = gsl_odeiv2_step_rk4imp;
  //data->s_type = gsl_odeiv2_step_bsimp;
  //data->s_type = gsl_odeiv2_step_msadams;
  //data->s_type = gsl_odeiv2_step_msbdf;

  data->s = gsl_odeiv2_step_alloc (data->s_type, pert->sys_size);
  data->sa = gsl_odeiv2_step_alloc (gsl_odeiv2_step_rkf45, pert->sys_size);
  data->e = gsl_odeiv2_evolve_alloc (pert->sys_size);
	data->d = NULL;
  data->c = NULL;

  data->sys.function = &ncm_gsl_odeiv2_step;
  data->sys.jacobian = &ncm_gsl_odeiv2_band_J;
  data->sys.dimension = pert->sys_size;
  data->sys.params = NULL;

  return data;
}

static void
ncm_gsl_odeiv2_set_opts (NcLinearPert *pert)
{
  GSLOde2Data *data = GSLODE2_DATA (pert->solver->data);
  data->sys.params = pert;

  return;
}

static void
ncm_gsl_odeiv2_reset (NcLinearPert *pert)
{
  GSLOde2Data *data = GSLODE2_DATA (pert->solver->data);

	pert->pws->dlambda = 1e-1;

  gsl_odeiv2_evolve_reset (data->e);

  if (data->c != NULL)
    gsl_odeiv2_control_free (data->c);
  data->c = gsl_odeiv2_control_scaled_new (1.0, pert->reltol, 1.0, 0.0, data->abstol->data, pert->sys_size);

	if (data->d != NULL)
		gsl_odeiv2_driver_free (data->d);

	data->d = gsl_odeiv2_driver_alloc_scaled_new (&data->sys, data->s_type, pert->pws->dlambda, 1.0, pert->reltol, 1.0, 0.0, data->abstol->data);

	gsl_odeiv2_control_set_driver (data->c, data->d);

  ncm_gsl_odeiv2_set_opts (pert);

  return;
}

static void
change_alg (NcLinearPert *pert)
{
  GSLOde2Data *data = GSLODE2_DATA (pert->solver->data);
  gsl_odeiv2_step *tmp = data->s;
  data->s = data->sa;
  data->sa = tmp;
  ncm_gsl_odeiv2_reset (pert);
}

static gboolean
ncm_gsl_odeiv2_evol_step (NcLinearPert *pert, gdouble lambda)
{
  GSLOde2Data *data = GSLODE2_DATA (pert->solver->data);
  gint status;

  status = gsl_odeiv2_evolve_apply (data->e, data->c, data->s, &data->sys, &pert->pws->lambda, lambda, &pert->pws->dlambda, data->y->data);
  if (status != GSL_SUCCESS)
    g_error ("Argg!! %d\n", status);

  if (pert->pws->tight_coupling && pert->pws->tight_coupling_end)
    ncm_gsl_odeiv2_end_tight_coupling (pert);

  if (pert->pws->lambda > -4.5)
    change_alg (pert);

  if (pert->pws->lambda == lambda)
    return TRUE;
  else
    return FALSE;
}

static gboolean
ncm_gsl_odeiv2_evol (NcLinearPert *pert, gdouble lambda)
{
  GSLOde2Data *data = GSLODE2_DATA (pert->solver->data);
	gdouble flambda = lambda;

	//g = 1e-3;
  while (flambda > pert->pws->lambda)
  {
		printf ("# Aqui\n"); fflush (stdout);
    gint status = gsl_odeiv2_evolve_apply (data->e, data->c, data->s, &data->sys, &pert->pws->lambda, lambda, &pert->pws->dlambda, data->y->data);
		//gint status = gsl_odeiv2_driver_apply (data->d, &pert->pws->g, g, data->y->data);
		printf ("# Aqui\n"); fflush (stdout);
		//g *= 1.1;
		printf ("%d % 20.15g % 20.15g % 20.15g\n", status, pert->pws->lambda, lambda, pert->pws->dlambda);
    if (pert->pws->tight_coupling && pert->pws->tight_coupling_end)
      ncm_gsl_odeiv2_end_tight_coupling (pert);
  }
  if (pert->pws->lambda == lambda)
    return TRUE;
  else
    return FALSE;
}

static gboolean
ncm_gsl_odeiv2_update_los (NcLinearPert *pert)
{
  g_assert_not_reached ();
  return TRUE;
}


static void
ncm_gsl_odeiv2_free (NcLinearPert *pert)
{
//  GSLOde2Data *data = GSLODE2_DATA (pert->solver->data);
}

static void
ncm_gsl_odeiv2_print_stats (NcLinearPert *pert)
{
//  GSLOde2Data *data = GSLODE2_DATA (pert->solver->data);
}

#define LINEAR_VECTOR_PREPARE gdouble *y = GSLODE2_DATA(pert->solver->data)->y->data
#define LINEAR_VEC_COMP(v,i) ((v)[(i)])
#define LINEAR_MATRIX_E(M,i,j) ((M)[i * pert->sys_size + j])
#define LINEAR_VECTOR_SET_ALL(v,c,n) do {gint _i_i; for (_i_i = 0; _i_i < (n); _i_i++) (v)[_i_i] = (c); } while (FALSE)
#define LINEAR_VEC_ABSTOL(pert) (GSLODE2_DATA(pert->solver->data)->abstol->data)
#define LINEAR_VEC_LOS_THETA(pert) (GSLODE2_DATA(pert->solver->data)->y->data)
#define LINEAR_STEP_RET gint
#define LINEAR_NAME_SUFFIX(base) ncm_gsl_odeiv2_##base
#define LINEAR_STEP_PARAMS gdouble lambda, const gdouble y[], gdouble ydot[], gpointer user_data
#define LINEAR_JAC_PARAMS gdouble lambda, const gdouble y[], gdouble *J, gdouble dfdt[], gpointer user_data
#define LINEAR_STEP_RET_VAL return GSL_SUCCESS

#include "linear_generic.c"

#undef LINEAR_VECTOR_PREPARE
#undef LINEAR_VECTOR_SET_ALL
#undef LINEAR_VEC_COMP
#undef LINEAR_STEP_RET
#undef LINEAR_NAME_SUFFIX
#undef LINEAR_STEP_PARAMS
#undef LINEAR_STEP_RET_VAL

#endif /* HAVE_GSL_ODEIV2 */
