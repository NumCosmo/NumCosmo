/***************************************************************************
 *            nc_hicosmo_qgrw.h
 *
 *  Wed June 04 10:04:37 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>>
 ****************************************************************************/
/*
 * nc_hicosmo_qgrw.h
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_HICOSMO_QGRW_H_
#define _NC_HICOSMO_QGRW_H_

#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/perturbations/nc_hipert_itwo_fluids.h>

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_QGRW             (nc_hicosmo_qgrw_get_type ())
#define NC_HICOSMO_QGRW(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_QGRW, NcHICosmoQGRW))
#define NC_HICOSMO_QGRW_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_QGRW, NcHICosmoQGRWClass))
#define NC_IS_HICOSMO_QGRW(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_QGRW))
#define NC_IS_HICOSMO_QGRW_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_QGRW))
#define NC_HICOSMO_QGRW_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_QGRW, NcHICosmoQGRWClass))

#define NC_TYPE_HIPERT_WKB_QGRW_ZETA     (nc_hipert_wkb_qgrw_zeta_get_type ())

typedef struct _NcHICosmoQGRWClass NcHICosmoQGRWClass;
typedef struct _NcHICosmoQGRW NcHICosmoQGRW;

/**
 * NcHICosmoQGRWSParams:
 * @NC_HICOSMO_QGRW_H0: Hubble constant.
 * @NC_HICOSMO_QGRW_OMEGA_R: Radiation density at $a_0$.
 * @NC_HICOSMO_QGRW_OMEGA_W: $w$-fluid density at $a_0$.
 * @NC_HICOSMO_QGRW_W: $w$-fluid equation of state.
 * @NC_HICOSMO_QGRW_X_B: Redshift at the bounce.
 *
 * Parameter of the Quantum Gravity Radiation W model.
 * 
 */
typedef enum /*< enum,underscore_name=NC_HICOSMO_QGRW_SPARAMS >*/
{
  NC_HICOSMO_QGRW_H0 = 0,
  NC_HICOSMO_QGRW_OMEGA_R,
  NC_HICOSMO_QGRW_OMEGA_W,
  NC_HICOSMO_QGRW_W,  
  NC_HICOSMO_QGRW_X_B,    
  /* < private > */
  NC_HICOSMO_QGRW_SPARAM_LEN, /*< skip >*/
} NcHICosmoQGRWSParams;

/**
 * NC_HICOSMO_QGRW_DEFAULT_H0: (value 73.0)
 * 
 * Default value for $H_0$.
 */ 
#define NC_HICOSMO_QGRW_DEFAULT_H0      ncm_c_hubble_cte_wmap ()

/**
 * NC_HICOSMO_QGRW_DEFAULT_OMEGA_R: (value 1.0e-5)
 * 
 * Default $\Omega_{r0}$.
 */
#define NC_HICOSMO_QGRW_DEFAULT_OMEGA_R (1.0e-5)

/**
 * NC_HICOSMO_QGRW_DEFAULT_OMEGA_W: (value 0.99999)
 * 
 * Default $\Omega_{w0}$.
 */
#define NC_HICOSMO_QGRW_DEFAULT_OMEGA_W (1.0 - NC_HICOSMO_QGRW_DEFAULT_OMEGA_R)

/**
 * NC_HICOSMO_QGRW_DEFAULT_W: (value 1.0e-12)
 * 
 * Default $w$.
 */
#define NC_HICOSMO_QGRW_DEFAULT_W       (1.0e-12)

/**
 * NC_HICOSMO_QGRW_DEFAULT_OMEGA_X_B: (value 1.0e30)
 * 
 * Default $x_b$.
 */
#define NC_HICOSMO_QGRW_DEFAULT_X_B     (1.0e30)

struct _NcHICosmoQGRWClass
{
  NcHICosmoClass parent_class;
  /*< private >*/
};

struct _NcHICosmoQGRW
{
  NcHICosmo parent_instance;
  /*< private >*/
  NcHIPertITwoFluidsEOM eom_two_fluids;
  NcHIPertITwoFluidsTV tv_two_fluids;
};

GType nc_hicosmo_qgrw_get_type (void) G_GNUC_CONST;
GType nc_hipert_wkb_qgrw_zeta_get_type (void) G_GNUC_CONST;

NcHICosmoQGRW *nc_hicosmo_qgrw_new (void);

G_END_DECLS

#endif /* _NC_HICOSMO_QGRW_H_ */
