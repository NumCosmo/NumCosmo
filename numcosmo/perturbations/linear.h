/***************************************************************************
 *            linear.h
 *
 *  Sat Oct 25 21:02:53 2008
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

#ifndef _NC_LINEAR_H
#define _NC_LINEAR_H

#include <glib.h>

G_BEGIN_DECLS

/**
 * NcLinearPertVars:
 * @NC_PERT_B0: FIXME
 * @NC_PERT_THETA0: FIXME  
 * @NC_PERT_C0: FIXME
 * @NC_PERT_PHI: FIXME
 * @NC_PERT_B1: FIXME
 * @NC_PERT_THETA1: FIXME
 * @NC_PERT_C1: FIXME
 * @NC_PERT_THETA2: FIXME
 * @NC_PERT_THETA_P0: FIXME
 * @NC_PERT_THETA_P1: FIXME
 * @NC_PERT_THETA_P2: FIXME
 *
 * FIXME
 */  
typedef enum _NcLinearPertVars
{
  NC_PERT_B0 = 0,
  NC_PERT_THETA0,  
  NC_PERT_C0,
  NC_PERT_PHI,
  NC_PERT_B1,
  NC_PERT_THETA1,
  NC_PERT_C1,
  NC_PERT_THETA2,
  NC_PERT_THETA_P0,
  NC_PERT_THETA_P1,
  NC_PERT_THETA_P2
} NcLinearPertVars;

#define NC_PERTURBATION_BASE_SIZE (NC_PERT_THETA2 + 1)

#define NC_PERT_dB0     NC_PERT_B0
#define NC_PERT_V       NC_PERT_C1
#define NC_PERT_T       NC_PERT_B1
#define NC_PERT_dTHETA0 NC_PERT_THETA0
#define NC_PERT_U       NC_PERT_THETA1

extern gint _itheta_table[3];
extern gint _itheta_p_table[3];

#define NC_PERT_THETA(n)   ((n <= 2) ? (_itheta_table[n])   : (NC_PERT_THETA_P2 + 1) + (2*(n-3)))
#define NC_PERT_THETA_P(n) ((n <= 2) ? (_itheta_p_table[n]) : (NC_PERT_THETA_P2 + 1) + (2*(n-3)+1))

typedef struct _NcLinearPertWorkSpace
{
  struct _NcLinearPert *pert;
  gboolean tight_coupling;
  gboolean tight_coupling_end;
  gdouble g_int;
  gdouble g;
  gdouble dg;
  gdouble k;
} NcLinearPertWorkSpace;

typedef struct _NcLinearPert NcLinearPert;

/**
 * NcLinearPert:
 * 
 * FIXME
 */
struct _NcLinearPert
{
  /*< private >*/
  NcmModel *model;
  NcThermodynRecomb *recomb;
  NcScaleFactor *a;
  struct _NcLinearPertOdeSolver *solver;
  NcLinearPertWorkSpace *pws;
  gdouble eta0;
  gdouble g0;
  gdouble gf;
  gdouble g_opt_cutoff;
  gdouble g_rec;
  gdouble g_rec_10m2_max[2];
  gdouble abstol;
  gdouble reltol;
  gdouble tc_abstol;
  gdouble tc_reltol;
  guint lmax;
  guint sys_size;
};

#define NC_PERTURBATIONS_G2X(g0,g) (exp ((g0) - (g)))
#define NC_PERTURBATIONS_X2G(g0,x) (-log (x) + (g0))

typedef struct _NcLinearPertOdeSolver NcLinearPertOdeSolver;

/**
 * NcLinearPertOdeSolver:
 * 
 * FIXME
 */
struct _NcLinearPertOdeSolver
{
  /*< private >*/
  gpointer (*create) (NcLinearPert *pert);
  void (*init) (NcLinearPert *pert);
  void (*set_opts) (NcLinearPert *pert);
  void (*reset) (NcLinearPert *pert);
  gboolean (*evol_step) (NcLinearPert *pert, gdouble g);
  gboolean (*evol) (NcLinearPert *pert, gdouble g);
  gboolean (*update_los) (NcLinearPert *pert);
  void (*get_sources) (NcLinearPert *pert, gdouble *S0, gdouble *S1, gdouble *S2);
  void (*free) (NcLinearPert *pert);
  void (*print_stats) (NcLinearPert *pert);
  gdouble (*get_z) (NcLinearPert *pert);
  gdouble (*get_phi) (NcLinearPert *pert);
  gdouble (*get_c0) (NcLinearPert *pert);
  gdouble (*get_b0) (NcLinearPert *pert);
  gdouble (*get_c1) (NcLinearPert *pert);
  gdouble (*get_b1) (NcLinearPert *pert);
  gdouble (*get) (NcLinearPert *pert, guint n);
  gdouble (*get_theta) (NcLinearPert *pert, gint n);
  gdouble (*get_theta_p) (NcLinearPert *pert, gint n);
  gdouble (*get_los_theta) (NcLinearPert *pert, gint n);
	void (*print_all) (NcLinearPert *pert);
  gpointer data;
};

/**
 * NcLinearPertSplineTypes:
 * @NC_LINEAR_PERTURBATIONS_SPLINE_SOURCES: FIXME
 * @NC_LINEAR_PERTURBATIONS_SPLINE_PHI: FIXME
 * @NC_LINEAR_PERTURBATIONS_SPLINE_THETA0: FIXME
 * @NC_LINEAR_PERTURBATIONS_SPLINE_C0: FIXME
 * @NC_LINEAR_PERTURBATIONS_SPLINE_B0: FIXME
 * @NC_LINEAR_PERTURBATIONS_SPLINE_THETA1: FIXME
 * @NC_LINEAR_PERTURBATIONS_SPLINE_C1: FIXME
 * @NC_LINEAR_PERTURBATIONS_SPLINE_B1: FIXME
 * @NC_LINEAR_PERTURBATIONS_SPLINE_THETA2: FIXME
 * 
 * FIXME
 */
typedef enum _NcLinearPertSplineTypes
{
  NC_LINEAR_PERTURBATIONS_SPLINE_SOURCES = 1 << NC_PERTURBATION_BASE_SIZE,
  NC_LINEAR_PERTURBATIONS_SPLINE_PHI     = 1 << NC_PERT_PHI,
  NC_LINEAR_PERTURBATIONS_SPLINE_THETA0  = 1 << NC_PERT_THETA0,
  NC_LINEAR_PERTURBATIONS_SPLINE_C0      = 1 << NC_PERT_C0,
  NC_LINEAR_PERTURBATIONS_SPLINE_B0      = 1 << NC_PERT_B0,
  NC_LINEAR_PERTURBATIONS_SPLINE_THETA1  = 1 << NC_PERT_THETA1,
  NC_LINEAR_PERTURBATIONS_SPLINE_C1      = 1 << NC_PERT_C1,
  NC_LINEAR_PERTURBATIONS_SPLINE_B1      = 1 << NC_PERT_B1,
  NC_LINEAR_PERTURBATIONS_SPLINE_THETA2  = 1 << NC_PERT_THETA2,
} NcLinearPertSplineTypes;

#define NC_LINEAR_PERTURBATIONS_SPLINE_ALL (~0)

typedef struct _NcLinearPertSplines NcLinearPertSplines;

/**
 * NcLinearPertSplines:
 */
struct _NcLinearPertSplines
{
  /*< private >*/
  NcLinearPert *pert;
  gulong n_deta, n_evol;
  gdouble k0, k1;
  NcLinearPertSplineTypes types;
  NcmVector *ga;
  NcmVector *Sg_data[3];
  NcmVector *ka;
  NcmMatrix *Sk_data[3];
  NcmVector *sdata[NC_PERTURBATION_BASE_SIZE];
  NcmSpline **Sk[3];
  NcmSpline *Sg[3];
  NcmSpline *s[NC_PERTURBATION_BASE_SIZE];
  NcmSpline *Nc;
};

#define NC_LINEAR_PERTURBATIONS_GET_SPLINE(pspline,n) ((pspline)->s[n])
#define NC_LINEAR_PERTURBATIONS(p) ((NcLinearPert *)(p))
#define NC_PERTURBATION_START_X (1.0e12)

typedef struct _NcLinearPertTF NcLinearPertTF;

/**
 * NcLinearPertTF:
 * 
 * FIXME
 */
struct _NcLinearPertTF
{
  /*< private >*/
  NcLinearPert *pert;
  gdouble logk0;
  gdouble logk1;
  gulong np;
  NcmSpline *logPhi_logk;
};

NcLinearPert *nc_pert_linear_new (NcmModel *model, guint lmax, gdouble tc_reltol, gdouble reltol, gdouble tc_abstol, gdouble abstol);
NcLinearPertSplines *nc_pert_linear_splines_new (NcLinearPert *pert, NcLinearPertSplineTypes types, gulong n_deta, gulong n_evol, gdouble k0, gdouble k1);
void nc_pert_linear_prepare_splines (NcLinearPertSplines *pspline);

void nc_pert_linear_free (NcLinearPert * pert);
void nc_pert_linear_splines_free (NcLinearPertSplines *pspline);

gboolean nc_pert_linear_spline_set_source_at (NcLinearPertSplines *pspline, gdouble k);
gboolean nc_pert_linear_calc_Nc_spline (NcLinearPertSplines *pspline, NcmSpline *pw_spline, GArray *los_table, gulong n_interp);
gdouble nc_pert_linear_los_integrate (NcLinearPertSplines *pspline, glong l, gdouble k);
GArray *nc_pert_linear_create_los_table (gint lmax_los, gint *los_ini, gint *los_step);

extern gint _nc_default_los_init[];
extern gint _nc_default_los_step[];

#define nc_pert_get_default_los_table(lmax) (nc_pert_linear_create_los_table (lmax, _nc_default_los_init, _nc_default_los_step))

extern NcLinearPertOdeSolver *cvodes_solver;
#ifdef HAVE_GSL_ODEIV2
extern NcLinearPertOdeSolver *ncm_gsl_odeiv2_solver;
#endif /* HAVE_GSL_ODEIV2 */

NcLinearPertTF *nc_pert_transfer_function_new (NcLinearPert *pert, gdouble k0, gdouble k1, gulong np);
void nc_pert_transfer_function_prepare (NcLinearPertTF *perttf);
gdouble nc_pert_transfer_function_get (NcLinearPertTF *perttf, gdouble kh);

G_END_DECLS

#endif /* LINEAR_H */
