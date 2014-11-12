/***************************************************************************
 *            linear_gsl_ode.c
 *
 *  Thu Nov 12 15:56:10 2009
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

#include <gsl/gsl_odeiv.h>
#include "linear_internal.h"

typedef struct _GSLOdeData
{
  const gsl_odeiv_step_type *s_type;
  gsl_odeiv_step *s;
  gsl_odeiv_control *c;
  gsl_odeiv_evolve *e;
  gsl_odeiv_system sys;
  gsl_vector *yi;
  gsl_vector *y;  
  gsl_vector *abstol;  
} GSLOdeData;

#define GSLODE_DATA(a) ((GSLOdeData *)(a))

static gpointer gsl_ode_create (NcLinearPert *pert);
static void gsl_ode_init (NcLinearPert *pert);
static void gsl_ode_set_opts (NcLinearPert *pert);
static void gsl_ode_reset (NcLinearPert *pert);
static void gsl_ode_end_tight_coupling (NcLinearPert *pert);
static gboolean gsl_ode_evol_step (NcLinearPert *pert, gdouble g);
static gboolean gsl_ode_evol (NcLinearPert *pert, gdouble g);
static gboolean gsl_ode_update_los (NcLinearPert *pert);
static void gsl_ode_get_sources (NcLinearPert *pert, gdouble *S0, gdouble *S1, gdouble *S2);
static void gsl_ode_free (NcLinearPert *pert);
static void gsl_ode_print_stats (NcLinearPert *pert);
static gdouble gsl_ode_get_z (NcLinearPert *pert);
static gdouble gsl_ode_get_phi (NcLinearPert *pert);
static gdouble gsl_ode_get_c0 (NcLinearPert *pert);
static gdouble gsl_ode_get_b0 (NcLinearPert *pert);
static gdouble gsl_ode_get_c1 (NcLinearPert *pert);
static gdouble gsl_ode_get_b1 (NcLinearPert *pert);
static gdouble gsl_ode_get_theta2 (NcLinearPert *pert);
static gdouble gsl_ode_get (NcLinearPert *pert, guint n);
static gdouble gsl_ode_get_theta (NcLinearPert *pert, guint n);
static gdouble gsl_ode_get_theta_p (NcLinearPert *pert, guint n);

static NcLinearPertOdeSolver _gsl_ode_solver = { 
  &gsl_ode_create,
  &gsl_ode_init,
  &gsl_ode_set_opts,
  &gsl_ode_reset,
  &gsl_ode_evol_step,
  &gsl_ode_evol,
  &gsl_ode_update_los,
  &gsl_ode_get_sources,
  &gsl_ode_free,
  &gsl_ode_print_stats,
  &gsl_ode_get_z,
  &gsl_ode_get_phi,
  &gsl_ode_get_c0,
  &gsl_ode_get_b0,
  &gsl_ode_get_c1,
  &gsl_ode_get_b1,
  &gsl_ode_get,
  &gsl_ode_get_theta,
  &gsl_ode_get_theta_p,
  NULL,
  NULL,
  NULL,
};
NcLinearPertOdeSolver *gsl_ode_solver = &_gsl_ode_solver;

static gint gsl_ode_step (gdouble lambda, const gdouble y[], gdouble ydot[], gpointer params);
static gint gsl_ode_band_J (gdouble lambda, const gdouble y[], gdouble *J, gdouble dfdt[], gpointer user_data);

static gpointer 
gsl_ode_create (NcLinearPert *pert)
{
  GSLOdeData *data = g_slice_new (GSLOdeData);

  data->yi = gsl_vector_alloc (pert->sys_size);
  data->y = gsl_vector_alloc (pert->sys_size);
  data->abstol = gsl_vector_alloc (pert->sys_size);

//  data->s_type = gsl_odeiv_step_rk2; 
//  data->s_type = gsl_odeiv_step_rk4; 
  data->s_type = gsl_odeiv_step_rkf45;
//  data->s_type = gsl_odeiv_step_rkck; 
//  data->s_type = gsl_odeiv_step_rk8pd;
//  data->s_type = gsl_odeiv_step_rk2imp;
//  data->s_type = gsl_odeiv_step_rk4imp;
//  data->s_type = gsl_odeiv_step_bsimp;
//  data->s_type = gsl_odeiv_step_gear1;
//  data->s_type = gsl_odeiv_step_gear2;
  
  data->s = gsl_odeiv_step_alloc (data->s_type, pert->sys_size);
  data->e = gsl_odeiv_evolve_alloc (pert->sys_size);
  data->c = NULL;
  
  data->sys.function = &gsl_ode_step;
  data->sys.jacobian = &gsl_ode_band_J;
  data->sys.dimension = pert->sys_size;
  data->sys.params = NULL;
  
  return data;
}

static void
gsl_ode_set_opts (NcLinearPert *pert)
{
  GSLOdeData *data = GSLODE_DATA (pert->solver->data);
  data->sys.params = pert;
  
  return;
}

static void
gsl_ode_reset (NcLinearPert *pert)
{
  GSLOdeData *data = GSLODE_DATA (pert->solver->data);  

  gsl_odeiv_evolve_reset (data->e);

  if (data->c != NULL)
    gsl_odeiv_control_free (data->c);
  data->c = gsl_odeiv_control_scaled_new (1.0, pert->reltol, 1.0, 0.0, data->abstol->data, pert->sys_size);

  pert->pws->dlambda = 1e-10;
  
  gsl_ode_set_opts (pert);

  return;
}

static gboolean
gsl_ode_evol_step (NcLinearPert *pert, gdouble lambda)
{
  GSLOdeData *data = GSLODE_DATA (pert->solver->data);

  gint status = gsl_odeiv_evolve_apply (data->e, data->c, data->s, &data->sys, &pert->pws->lambda, lambda, &pert->pws->dlambda, data->y->data);
  if (status != GSL_SUCCESS)
    g_error ("Argg!!");
  
  if (pert->pws->tight_coupling && pert->pws->tight_coupling_end)
    gsl_ode_end_tight_coupling (pert); 

  if (pert->pws->lambda == lambda)
    return TRUE;
  else
    return FALSE;
}

static gboolean
gsl_ode_evol (NcLinearPert *pert, gdouble lambda)
{
  GSLOdeData *data = GSLODE_DATA (pert->solver->data);

  while (lambda > pert->pws->lambda)
  {
    gsl_odeiv_evolve_apply (data->e, data->c, data->s, &data->sys, &pert->pws->lambda, lambda, &pert->pws->dlambda, data->y->data);
    
    if (pert->pws->tight_coupling && pert->pws->tight_coupling_end)
      gsl_ode_end_tight_coupling (pert);
  }
  if (pert->pws->lambda == lambda)
    return TRUE;
  else
    return FALSE;
}

static gboolean
gsl_ode_update_los (NcLinearPert *pert)
{
  NCM_UNUSED (pert);
  g_assert_not_reached ();
  return TRUE;
}

static void 
gsl_ode_free (NcLinearPert *pert)
{
  NCM_UNUSED (pert);
  //  GSLOdeData *data = GSLODE_DATA (pert->solver->data);
}

static void 
gsl_ode_print_stats (NcLinearPert *pert)
{
  NCM_UNUSED (pert);
//  GSLOdeData *data = GSLODE_DATA (pert->solver->data);
}

#define LINEAR_VECTOR_PREPARE gdouble *y = GSLODE_DATA(pert->solver->data)->y->data
#define LINEAR_VEC_COMP(v,i) ((v)[(i)])
#define LINEAR_MATRIX_E(M,i,j) ((M)[i * pert->sys_size + j])
#define LINEAR_VECTOR_SET_ALL(v,c,n) G_STMT_START {guint _i_i; for (_i_i = 0; _i_i < (n); _i_i++) (v)[_i_i] = (c); } G_STMT_END
#define LINEAR_VEC_ABSTOL(pert) (GSLODE_DATA(pert->solver->data)->abstol->data)
#define LINEAR_VEC_LOS_THETA(pert) (GSLODE_DATA(pert->solver->data)->y->data)
#define LINEAR_STEP_RET gint
#define LINEAR_NAME_SUFFIX(base) gsl_ode_##base
#define LINEAR_STEP_PARAMS gdouble lambda, const gdouble y[], gdouble ydot[], gpointer user_data
#define LINEAR_JAC_PARAMS gdouble lambda, const gdouble y[], gdouble *J, gdouble dfdt[], gpointer user_data
#define LINEAR_STEP_RET_VAL return 0
#define LINEAR_JAC_UNUSED G_STMT_START { NCM_UNUSED (y); NCM_UNUSED (dfdt); } G_STMT_END

#include "linear_generic.c"

#undef LINEAR_VECTOR_PREPARE
#undef LINEAR_VECTOR_SET_ALL
#undef LINEAR_VEC_COMP
#undef LINEAR_STEP_RET
#undef LINEAR_NAME_SUFFIX
#undef LINEAR_STEP_PARAMS
#undef LINEAR_STEP_RET_VAL
