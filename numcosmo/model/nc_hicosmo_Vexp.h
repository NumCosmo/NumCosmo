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

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_VEXP             (nc_hicosmo_Vexp_get_type ())
#define NC_HICOSMO_VEXP(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_VEXP, NcHICosmoVexp))
#define NC_HICOSMO_VEXP_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_VEXP, NcHICosmoVexpClass))
#define NC_IS_HICOSMO_VEXP(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_VEXP))
#define NC_IS_HICOSMO_VEXP_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_VEXP))
#define NC_HICOSMO_VEXP_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_VEXP, NcHICosmoVexpClass))

typedef struct _NcHICosmoVexpClass NcHICosmoVexpClass;
typedef struct _NcHICosmoVexp NcHICosmoVexp;
typedef struct _NcHICosmoVexpPrivate NcHICosmoVexpPrivate;

struct _NcHICosmoVexpClass
{
  /*< private >*/
  NcHICosmoClass parent_class;
};

/**
 * NcHICosmoVexpSParams:
 * @NC_HICOSMO_VEXP_H0: Hubble constant
 * @NC_HICOSMO_VEXP_OMEGA_C: scalar field energy density when in a dust-like phase
 * @NC_HICOSMO_VEXP_OMEGA_L: scalar field energy density when in a dark energy-like phase
 * @NC_HICOSMO_VEXP_SIGMA_PHI: standard deviation of the scalar field wave function
 * @NC_HICOSMO_VEXP_D_PHI: mean of the scalar field wave function distribution
 * @NC_HICOSMO_VEXP_ALPHA_B: lorgarithm base e of the scale factor at the bounce
 * @NC_HICOSMO_VEXP_X_B: ratio of the scale factor today and at the bounce 
 *
 * FIXME
 * 
 */
typedef enum /*< enum,underscore_name=NC_HICOSMO_VEXP_SPARAMS >*/
{
  NC_HICOSMO_VEXP_H0 = 0,
  NC_HICOSMO_VEXP_OMEGA_C,
  NC_HICOSMO_VEXP_OMEGA_L,
  NC_HICOSMO_VEXP_SIGMA_PHI,
  NC_HICOSMO_VEXP_D_PHI,
  NC_HICOSMO_VEXP_ALPHA_B,
  NC_HICOSMO_VEXP_X_B,        
  /* < private > */
  NC_HICOSMO_VEXP_SPARAM_LEN, /*< skip >*/
} NcHICosmoVexpSParams;

struct _NcHICosmoVexp
{
  /*< private >*/
  NcHICosmo parent_instance;
  NcHICosmoVexpPrivate *priv;
};

GType nc_hicosmo_Vexp_get_type (void) G_GNUC_CONST;

NcHICosmoVexp *nc_hicosmo_Vexp_new (void);

gdouble nc_hicosmo_Vexp_tau_min (NcHICosmoVexp *Vexp);
gdouble nc_hicosmo_Vexp_tau_max (NcHICosmoVexp *Vexp);

gdouble nc_hicosmo_Vexp_tau_qt_c (NcHICosmoVexp *Vexp);
gdouble nc_hicosmo_Vexp_tau_qt_e (NcHICosmoVexp *Vexp);

gdouble nc_hicosmo_Vexp_xbe (NcHICosmoVexp *Vexp);
gdouble nc_hicosmo_Vexp_xbc (NcHICosmoVexp *Vexp);

gdouble nc_hicosmo_Vexp_x_tau (NcHICosmoVexp *Vexp, const gdouble tau);

gdouble nc_hicosmo_Vexp_tau_xe (NcHICosmoVexp *Vexp, const gdouble xe);
gdouble nc_hicosmo_Vexp_tau_xc (NcHICosmoVexp *Vexp, const gdouble xc);

gdouble nc_hicosmo_Vexp_alpha (NcHICosmoVexp *Vexp, const gdouble tau);
gdouble nc_hicosmo_Vexp_phi (NcHICosmoVexp *Vexp, const gdouble tau);
gdouble nc_hicosmo_Vexp_Ricci_scale (NcHICosmoVexp *Vexp, const gdouble tau);
void nc_hicosmo_Vexp_x_y (NcHICosmoVexp *Vexp, const gdouble tau, gdouble *x, gdouble *y);

#define NC_HICOSMO_VEXP_DEFAULT_H0 (70.0)
#define NC_HICOSMO_VEXP_DEFAULT_OMEGA_C (0.25)
#define NC_HICOSMO_VEXP_DEFAULT_OMEGA_L (0.75)
#define NC_HICOSMO_VEXP_DEFAULT_SIGMA_PHI (0.4)
#define NC_HICOSMO_VEXP_DEFAULT_D_PHI (-0.3)
#define NC_HICOSMO_VEXP_DEFAULT_ALPHA_B (0.1)
#define NC_HICOSMO_VEXP_DEFAULT_X_B (1.0e30)

#define NC_HICOSMO_VEXP_DEBUG_EVOL_QT (FALSE)
#define NC_HICOSMO_VEXP_DEBUG_EVOL_CL (FALSE)

G_END_DECLS

#endif /* _NC_HICOSMO_VEXP_H_ */
