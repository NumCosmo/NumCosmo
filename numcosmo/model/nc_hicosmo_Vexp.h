/***************************************************************************
 *            nc_hicosmo_Vexp.h
 *
 *  Fri October 28 13:27:25 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2016 <sandro@isoftware.com.br>
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

#ifndef _NC_HICOSMO_VEXP_H_
#define _NC_HICOSMO_VEXP_H_

#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_hicosmo.h>
#include <nvector/nvector_serial.h>

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_VEXP             (nc_hicosmo_Vexp_get_type ())
#define NC_HICOSMO_VEXP(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_VEXP, NcHICosmoVexp))
#define NC_HICOSMO_VEXP_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_VEXP, NcHICosmoVexpClass))
#define NC_IS_HICOSMO_VEXP(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_VEXP))
#define NC_IS_HICOSMO_VEXP_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_VEXP))
#define NC_HICOSMO_VEXP_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_VEXP, NcHICosmoVexpClass))

typedef struct _NcHICosmoVexpClass NcHICosmoVexpClass;
typedef struct _NcHICosmoVexp NcHICosmoVexp;

struct _NcHICosmoVexpClass
{
  /*< private >*/
  NcHICosmoClass parent_class;
};

/**
 * NcHICosmoVexpParams:
 * @NC_HICOSMO_VEXP_H0: FIXME
 * @NC_HICOSMO_VEXP_OMEGA_C: FIXME
 * @NC_HICOSMO_VEXP_OMEGA_L: FIXME
 * @NC_HICOSMO_VEXP_SIGMA_PHI: FIXME
 * @NC_HICOSMO_VEXP_D_PHI: FIXME
 * @NC_HICOSMO_VEXP_ALPHA_B: FIXME
 * @NC_HICOSMO_VEXP_X_B: FIXME
 *
 * FIXME
 * 
 */
typedef enum _NcHICosmoVexpParams
{
  NC_HICOSMO_VEXP_H0 = 0,
  NC_HICOSMO_VEXP_OMEGA_C,
  NC_HICOSMO_VEXP_OMEGA_L,
  NC_HICOSMO_VEXP_SIGMA_PHI,
  NC_HICOSMO_VEXP_D_PHI,
  NC_HICOSMO_VEXP_ALPHA_B,
  NC_HICOSMO_VEXP_X_B,        /*< private >*/
  NC_HICOSMO_VEXP_SPARAM_LEN, /*< skip >*/
} NcHICosmoVexpParams;

struct _NcHICosmoVexp
{
  /*< private >*/
  NcHICosmo parent_instance;
  gpointer cvode_qt;
  gpointer cvode_clp;
  gpointer cvode_clm;
  gboolean qt_init;
  gboolean clm_init;
  gboolean clp_init;
  gboolean glue_de;
  N_Vector y_qt;
  N_Vector ydot_qt;
  N_Vector y_cl;
  gint cl_bc, cl_be;
  gdouble RH_lp;
  gdouble alpha_b;
  gdouble a_0de;
  gdouble a_0c, a_0e;
  gdouble qc, qe;
  gdouble Ec, Ee;
  gdouble alpha_qc, alpha_qe;
  gdouble alpha_0c, alpha_0e;
  gdouble tau_x0;
  gdouble c1c, c1e;
  gdouble c2c, c2e;
  GArray *evol_c;
  GArray *evol_e;
  NcmSpline *xtau_s;
  NcmSpline *lnN_s;
  NcmSpline *tau_dlnx_dtau_s;
  NcmSpline *E2_s;
};

GType nc_hicosmo_Vexp_get_type (void) G_GNUC_CONST;

NcHICosmoVexp *nc_hicosmo_Vexp_new (void);

gdouble nc_hicosmo_Vexp_tau_min (NcHICosmoVexp *Vexp);
gdouble nc_hicosmo_Vexp_tau_max (NcHICosmoVexp *Vexp);

#define NC_HICOSMO_VEXP_DEFAULT_H0 (70.0)
#define NC_HICOSMO_VEXP_DEFAULT_OMEGA_C (0.25)
#define NC_HICOSMO_VEXP_DEFAULT_OMEGA_L (0.75)
#define NC_HICOSMO_VEXP_DEFAULT_SIGMA_PHI (0.4)
#define NC_HICOSMO_VEXP_DEFAULT_D_PHI (-0.3)
#define NC_HICOSMO_VEXP_DEFAULT_ALPHA_B (0.1)
#define NC_HICOSMO_VEXP_DEFAULT_X_B (1.0e30)

G_END_DECLS

#endif /* _NC_HICOSMO_VEXP_H_ */
