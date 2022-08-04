/***************************************************************************
 *            nc_hicosmo_sfb.h
 *
 *  Wed June 04 10:04:37 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>>
 ****************************************************************************/
/*
 * nc_hicosmo_sfb.h
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

#ifndef _NC_HICOSMO_SFB_H_
#define _NC_HICOSMO_SFB_H_

#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/perturbations/nc_hipert_adiab.h>

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_SFB             (nc_hicosmo_sfb_get_type ())
#define NC_HICOSMO_SFB(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_SFB, NcHICosmoSFB))
#define NC_HICOSMO_SFB_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_SFB, NcHICosmoSFBClass))
#define NC_IS_HICOSMO_SFB(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_SFB))
#define NC_IS_HICOSMO_SFB_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_SFB))
#define NC_HICOSMO_SFB_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_SFB, NcHICosmoSFBClass))


typedef struct _NcHICosmoSFBClass NcHICosmoSFBClass;
typedef struct _NcHICosmoSFB NcHICosmoSFB;

/**
 * NcHICosmoSFBSParams:
 * @NC_HICOSMO_SFB_H0: Hubble constant.
 * @NC_HICOSMO_SFB_OMEGA_R: Radiation density at $a_0$.
 * @NC_HICOSMO_SFB_OMEGA_W: $w$-fluid density at $a_0$.
 * @NC_HICOSMO_SFB_W: $w$-fluid equation of state.
 * @NC_HICOSMO_SFB_X_B: Redshift at the bounce.
 *
 * 
 */
typedef enum /*< enum,underscore_name=NC_HICOSMO_SFB_SPARAMS >*/
{
  NC_HICOSMO_SFB_H0 = 0,
  NC_HICOSMO_SFB_OMEGA_R,
  NC_HICOSMO_SFB_OMEGA_W,
  NC_HICOSMO_SFB_W,  
  NC_HICOSMO_SFB_X_B,    
  NC_HICOSMO_SFB_TAU_B,
  /* < private > */
  NC_HICOSMO_SFB_SPARAM_LEN, /*< skip >*/
} NcHICosmoSFBSParams;

/**
 * NC_HICOSMO_SFB_DEFAULT_H0: (value 73.0)
 * 
 * Default value for $H_0$.
 */ 
#define NC_HICOSMO_SFB_DEFAULT_H0      ncm_c_hubble_cte_planck6_base ()

/**
 * NC_HICOSMO_SFB_DEFAULT_OMEGA_R: (value 1.0e-5)
 * 
 * Default $\Omega_{r0}$.
 */
#define NC_HICOSMO_SFB_DEFAULT_OMEGA_R (1.0e-5)

/**
 * NC_HICOSMO_SFB_DEFAULT_OMEGA_W: (value 0.9999)
 * 
 * Default $\Omega_{w0}$.
 */
#define NC_HICOSMO_SFB_DEFAULT_OMEGA_W (1 - NC_HICOSMO_SFB_DEFAULT_OMEGA_R)

/**
 * NC_HICOSMO_SFB_DEFAULT_W: (value 0.33)
 * 
 * Default $w$.
 */
#define NC_HICOSMO_SFB_DEFAULT_W       (0.33)

/**
 * NC_HICOSMO_SFB_DEFAULT_OMEGA_X_B: (value 1.0e30)
 * 
 * Default $x_b$.
 */
#define NC_HICOSMO_SFB_DEFAULT_X_B     (1.0e25)

/**
 * NC_HICOSMO_SFB_DEFAULT_TAU_B: (value 0.0 )
 *
 * Default $tau_b$.
 */
#define NC_HICOSMO_SFB_DEFAULT_TAU_B     (1 / (NC_HICOSMO_SFB_DEFAULT_X_B * pow(NC_HICOSMO_SFB_DEFAULT_OMEGA_R, 0.5)))


struct _NcHICosmoSFBClass
{
  NcHICosmoClass parent_class;
  /*< private >*/
};

struct _NcHICosmoSFB
{
  NcHICosmo parent_instance;
  /*< private >*/
  /*NcHIPertITwoFluidsEOM eom_two_fluids;
  NcHIPertITwoFluidsTV tv_two_fluids;*/
};

GType nc_hicosmo_sfb_get_type (void) G_GNUC_CONST;

NcHICosmoSFB *nc_hicosmo_sfb_new (void);

G_END_DECLS

#endif /* _NC_HICOSMO_SFB_H_ */
