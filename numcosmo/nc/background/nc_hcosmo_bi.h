/***************************************************************************
 *            nc_hcosmo_bi.h
 *
 *  Sat Oct 17 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hcosmo_bi.h
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_HCOSMO_BI_H_
#define _NC_HCOSMO_BI_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/model/ncm_model.h>
#include <numcosmo/ncm/model/ncm_mset.h>
#include <numcosmo/ncm/algebra/ncm_quaternion.h>

G_BEGIN_DECLS

#define NC_TYPE_HCOSMO_BI (nc_hcosmo_bi_get_type ())

G_DECLARE_FINAL_TYPE (NcHCosmoBI, nc_hcosmo_bi, NC, HCOSMO_BI, NcmModel)

/**
 * NcHCosmoBISParams:
 * @NC_HCOSMO_BI_H1: Hubble rate of axis 1 today $H_1$ [km/s/Mpc]
 * @NC_HCOSMO_BI_H2: Hubble rate of axis 2 today $H_2$ [km/s/Mpc]
 * @NC_HCOSMO_BI_H3: Hubble rate of axis 3 today $H_3$ [km/s/Mpc]
 * @NC_HCOSMO_BI_OMEGA_M: matter density parameter today $\Omega_{m0}$
 * @NC_HCOSMO_BI_U1: orientation parameter $u_1 \in [0, 1]$
 * @NC_HCOSMO_BI_U2: orientation parameter $u_2 \in [0, 1]$
 * @NC_HCOSMO_BI_U3: orientation parameter $u_3 \in [0, 1]$
 *
 * #NcHCosmoBI model parameters.
 *
 */
typedef enum /*< enum,underscore_name=NC_HCOSMO_BI_SPARAMS,prefix=NC_HCOSMO_BI >*/
{
  NC_HCOSMO_BI_H1,
  NC_HCOSMO_BI_H2,
  NC_HCOSMO_BI_H3,
  NC_HCOSMO_BI_OMEGA_M,
  NC_HCOSMO_BI_U1,
  NC_HCOSMO_BI_U2,
  NC_HCOSMO_BI_U3,
  /* < private > */
  NC_HCOSMO_BI_SPARAM_LEN, /*< skip >*/
} NcHCosmoBISParams;

NCM_MSET_MODEL_DECLARE_ID (nc_hcosmo_bi);

NcHCosmoBI *nc_hcosmo_bi_new (void);
NcHCosmoBI *nc_hcosmo_bi_ref (NcHCosmoBI *cosmo);

void nc_hcosmo_bi_free (NcHCosmoBI *cosmo);
void nc_hcosmo_bi_clear (NcHCosmoBI **cosmo);

gdouble nc_hcosmo_bi_H0 (NcHCosmoBI *cosmo);
gdouble nc_hcosmo_bi_RH_Mpc (NcHCosmoBI *cosmo);
gdouble nc_hcosmo_bi_Omega_s0 (NcHCosmoBI *cosmo);
gdouble nc_hcosmo_bi_Omega_L0 (NcHCosmoBI *cosmo);
gdouble nc_hcosmo_bi_age (NcHCosmoBI *cosmo);

void nc_hcosmo_bi_ln_a (NcHCosmoBI *cosmo, const gdouble t_lb, gdouble *ln_a1, gdouble *ln_a2, gdouble *ln_a3);
void nc_hcosmo_bi_E (NcHCosmoBI *cosmo, const gdouble t_lb, gdouble *E1, gdouble *E2, gdouble *E3);

void nc_hcosmo_bi_get_quaternion (NcHCosmoBI *cosmo, NcmQuaternion *q);

#define NC_HCOSMO_BI_DEFAULT_H1 (70.0)
#define NC_HCOSMO_BI_DEFAULT_H2 (70.0)
#define NC_HCOSMO_BI_DEFAULT_H3 (70.0)
#define NC_HCOSMO_BI_DEFAULT_OMEGA_M (0.3)
#define NC_HCOSMO_BI_DEFAULT_U1 (1.0)
#define NC_HCOSMO_BI_DEFAULT_U2 (0.0)
#define NC_HCOSMO_BI_DEFAULT_U3 (0.0)

#define NC_HCOSMO_BI_DEFAULT_PARAMS_ABSTOL (0.0)

G_END_DECLS

#endif /* _NC_HCOSMO_BI_H_ */

