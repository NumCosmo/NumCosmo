/***************************************************************************
 *            recomb.h
 *
 *  Sun Oct  5 20:40:46 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
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

#ifndef _NC_RECOMB_H
#define _NC_RECOMB_H

#include <glib.h>
#include <cvodes/cvodes.h>             /* prototypes for CVODE fcts. and consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., and macros */
#include <cvodes/cvodes_dense.h>       /* prototype for CVDense */
#include <cvodes/cvodes_band.h>        /* prototype for CVBand */
#include <sundials/sundials_band.h>
#include <sundials/sundials_dense.h> /* definitions DenseMat DENSE_ELEM */   
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_spline.h>

G_BEGIN_DECLS

/**
 * NcThermodynRecombType:
 * @NC_THERMODYN_RECOMB_SINGLE: FIXME
 * @NC_THERMODYN_RECOMB_SINGLE_Tm: FIXME
 * @NC_THERMODYN_RECOMB_FULL: FIXME
 * 
 * FIXME
 */
typedef enum _NcThermodynRecombType 
{
  NC_THERMODYN_RECOMB_SINGLE = 0,
  NC_THERMODYN_RECOMB_SINGLE_Tm,
  NC_THERMODYN_RECOMB_FULL
} NcThermodynRecombType;

typedef struct _NcThermodynRecomb NcThermodynRecomb;

/**
 * NcThermodynRecomb:
 * 
 * FIXME
 */
struct _NcThermodynRecomb
{
  /*< private >*/
  gpointer cvode;
  NcHICosmo *model;
  gboolean init;
  guint n;
  gdouble x0;
  gdouble x_rec;
  gdouble x_rec_10m2_max[2];
  gdouble x_opt_cutoff;
  N_Vector y0;
  N_Vector y;
  gdouble abstol;
  gdouble reltol;
  CVRhsFn ion;
  CVDlsDenseJacFn ion_J;
  NcFunctionCache *H_ion_fraction;
  NcFunctionCache *optical_depth;
  gdouble He_single_ionized_x;
  gsl_interp_accel *dtau_dx_accel;
  gsl_interp_accel *tau_accel;
  gsl_spline *tau_spline;
  gsl_spline *log_g_spline;
  gsl_spline *dtau_dx_spline;
  gsl_spline *dtau_dxR_spline;
  gboolean init_spline;
};

gdouble nc_thermodyn_H_ionization_saha (NcHICosmo *model, gdouble x);
gdouble nc_thermodyn_HeI_ionization_saha (NcHICosmo *model, gdouble x);
gdouble nc_thermodyn_HeII_ionization_saha (NcHICosmo *model, gdouble x);

gdouble nc_thermodyn_HeII_ionization_saha_x (NcHICosmo *model, gdouble frac);

NcThermodynRecomb *nc_thermodyn_recomb_new (NcHICosmo *model, NcThermodynRecombType type);
void nc_thermodyn_recomb_free (NcThermodynRecomb *recomb);
void nc_thermodyn_recomb_dtau_dx_init_spline (NcThermodynRecomb *recomb);
gboolean nc_thermodyn_recomb_reinit (NcThermodynRecomb *recomb);
gboolean nc_thermodyn_recomb_reset (NcThermodynRecomb *recomb);
gboolean nc_thermodyn_recomb_evolve (NcThermodynRecomb *recomb, gdouble x);

gdouble nc_thermodyn_recomb_get_Xe (NcThermodynRecomb *recomb);
gdouble nc_thermodyn_recomb_get_Xe_at (NcThermodynRecomb *recomb, gdouble x);

gdouble nc_thermodyn_recomb_dtau_dx (NcThermodynRecomb *recomb, gdouble x);
gdouble nc_thermodyn_recomb_taubar (NcThermodynRecomb *recomb, gdouble x);
gdouble nc_thermodyn_recomb_taubbar (NcThermodynRecomb *recomb, gdouble x);
gdouble nc_thermodyn_recomb_taubbbar (NcThermodynRecomb *recomb, gdouble x);
void nc_thermodyn_recomb_taubbar_val (NcThermodynRecomb *recomb, gdouble x, gdouble *taubar, gdouble *taubbar);
gdouble nc_thermodyn_recomb_optical_depth (NcThermodynRecomb *recomb, gdouble x);
gdouble nc_thermodyn_recomb_optical_depth_x0_x1 (NcThermodynRecomb *recomb, gdouble x0, gdouble x1);
gdouble nc_thermodyn_recomb_optical_depth_R_x0_x1 (NcThermodynRecomb *recomb, gdouble x0, gdouble x1);
gdouble nc_thermodyn_recomb_g (NcThermodynRecomb *recomb, gdouble x);
gdouble nc_thermodyn_recomb_gbar (NcThermodynRecomb *recomb, gdouble x);
gdouble nc_thermodyn_recomb_gbbar (NcThermodynRecomb *recomb, gdouble x);

gdouble nc_thermodyn_recomb_peak_x (NcThermodynRecomb *recomb);
gdouble nc_thermodyn_recomb_log_g (NcThermodynRecomb *recomb, gdouble x);

gdouble nc_thermodyn_recomb_H_case_B (NcHICosmo *model, gdouble Tm);
gdouble nc_thermodyn_recomb_H_case_B_dTm (NcHICosmo *model, gdouble Tm);

gdouble nc_thermodyn_recomb_HeI_case_B (NcHICosmo *model, gdouble Tm);
gdouble nc_thermodyn_recomb_HeI_case_B_dTm (NcHICosmo *model, gdouble Tm);

gdouble nc_thermodyn_H_ionization_rate (NcHICosmo *model, gdouble XH, gdouble Tm, gdouble XHeII, gdouble x);
gdouble nc_thermodyn_HeII_ionization_rate (NcHICosmo *model, gdouble Xp, gdouble Tm, gdouble XHeII, gdouble x);

void nc_thermodyn_H_ionization_rate_grad (NcHICosmo *model, gdouble Xp, gdouble Tm, gdouble XHeII, gdouble x, gsl_vector *grad);
void nc_thermodyn_HeII_ionization_rate_grad (NcHICosmo *model, gdouble Xp, gdouble Tm, gdouble XHeII, gdouble x, gsl_vector *grad);

G_END_DECLS

#endif /* RECOMB_H */
