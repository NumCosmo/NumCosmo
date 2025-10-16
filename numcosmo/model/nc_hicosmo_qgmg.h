/***************************************************************************
 *            nc_hicosmo_qgmg.h
 *
 *  Thu October 16 10:26:33 2025
 *  Copyright  2025  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>>
 ****************************************************************************/
/*
 * nc_hicosmo_qgmg.h
 * Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_HICOSMO_QGMG_H_
#define _NC_HICOSMO_QGMG_H_

#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/perturbations/nc_hipert_itwo_fluids.h>

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_QGMG             (nc_hicosmo_qgmg_get_type ())
#define NC_HICOSMO_QGMG(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_QGMG, NcHICosmoQGMG))
#define NC_HICOSMO_QGMG_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_QGMG, NcHICosmoQGMGClass))
#define NC_IS_HICOSMO_QGMG(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_QGMG))
#define NC_IS_HICOSMO_QGMG_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_QGMG))
#define NC_HICOSMO_QGMG_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_QGMG, NcHICosmoQGMGClass))

typedef struct _NcHICosmoQGMGClass NcHICosmoQGMGClass;
typedef struct _NcHICosmoQGMG NcHICosmoQGMG;

/**
 * NcHICosmoQGMGSParams:
 * @NC_HICOSMO_QGMG_H0: Hubble constant.
 * @NC_HICOSMO_QGMG_OMEGA_R: Radiation density at $a_0$.
 * @NC_HICOSMO_QGMG_OMEGA_W: $w$-fluid density at $a_0$.
 * @NC_HICOSMO_QGMG_W: $w$-fluid equation of state.
 * @NC_HICOSMO_QGMG_X_B: Redshift at the bounce.
 *
 * Parameter of the Quantum Gravity Radiation W model.
 *
 */
typedef enum /*< enum,underscore_name=NC_HICOSMO_QGMG_SPARAMS >*/
{
  NC_HICOSMO_QGMG_H0 = 0,
  NC_HICOSMO_QGMG_OMEGA_R,
  NC_HICOSMO_QGMG_OMEGA_W,
  NC_HICOSMO_QGMG_W,
  NC_HICOSMO_QGMG_X_B,
  /* < private > */
  NC_HICOSMO_QGMG_SPARAM_LEN, /*< skip >*/
} NcHICosmoQGMGSParams;

/**
 * NC_HICOSMO_QGMG_DEFAULT_H0: (value 73.0)
 *
 * Default value for $H_0$.
 */
#define NC_HICOSMO_QGMG_DEFAULT_H0      ncm_c_hubble_cte_planck6_base ()

/**
 * NC_HICOSMO_QGMG_DEFAULT_OMEGA_R: (value 1.0e-5)
 *
 * Default $\Omega_{r0}$.
 */
#define NC_HICOSMO_QGMG_DEFAULT_OMEGA_R (1.0e-5)

/**
 * NC_HICOSMO_QGMG_DEFAULT_OMEGA_W: (value 0.99999)
 *
 * Default $\Omega_{w0}$.
 */
#define NC_HICOSMO_QGMG_DEFAULT_OMEGA_W (1.0 - NC_HICOSMO_QGMG_DEFAULT_OMEGA_R)

/**
 * NC_HICOSMO_QGMG_DEFAULT_W: (value 1.0e-12)
 *
 * Default $w$.
 */
#define NC_HICOSMO_QGMG_DEFAULT_W       (1.0e-12)

/**
 * NC_HICOSMO_QGMG_DEFAULT_OMEGA_X_B: (value 1.0e30)
 *
 * Default $x_b$.
 */
#define NC_HICOSMO_QGMG_DEFAULT_X_B     (1.0e30)

struct _NcHICosmoQGMGClass
{
  NcHICosmoClass parent_class;
  /*< private >*/
};

struct _NcHICosmoQGMG
{
  NcHICosmo parent_instance;
  /*< private >*/
  NcHIPertITwoFluidsEOM eom_two_fluids;
  NcHIPertITwoFluidsWKB wkb_two_fluids;
  NcHIPertITwoFluidsTV tv_two_fluids;
};

GType nc_hicosmo_qgmg_get_type (void) G_GNUC_CONST;

NcHICosmoQGMG *nc_hicosmo_qgmg_new (void);

G_END_DECLS

#endif /* _NC_HICOSMO_QGMG_H_ */

