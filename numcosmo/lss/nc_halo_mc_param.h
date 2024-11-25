/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_halo_mc_param.h
 *
 *  Sat Oct 12 11:51:53 2024
 *  Copyright  2024  Mariana Penna-Lima
 *  <pennalima@unb.br>
 ****************************************************************************/
/*
 * nc_halo_mc_param.h
 * Copyright (C) 2024 Mariana Penna-Lima <pennalima@unb.br>
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NC_HALO_MC_PARAM_H_
#define _NC_HALO_MC_PARAM_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_halo_mass_summary.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_MC_PARAM (nc_halo_mc_param_get_type ())

G_DECLARE_FINAL_TYPE (NcHaloMCParam, nc_halo_mc_param, NC, HALO_MC_PARAM, NcHaloMassSummary)

/**
 * NcHaloMCParamSParams:
 * @NC_HALO_MC_PARAM_C_DELTA: concentration parameter $r_\Delta$
 * @NC_HALO_MC_PARAM_LOG10M_DELTA: halo mass $\log_{10}(M_\Delta)$
 *
 * Fundamental parametrization of the profile $\rho(r)$.
 * Both mass and concentration are parameters.
 *
 */
typedef enum /*< enum,underscore_name=NC_HALO_MC_PARAM_SPARAMS >*/
{
  NC_HALO_MC_PARAM_C_DELTA,
  NC_HALO_MC_PARAM_LOG10M_DELTA,
  /* < private > */
  NC_HALO_MC_PARAM_SPARAM_LEN, /*< skip >*/
} NcHaloMCParamSParams;

#define NC_HALO_MC_PARAM_LOCAL_SPARAM_LEN (NC_HALO_MC_PARAM_SPARAM_LEN - 0)

NcHaloMCParam *nc_halo_mc_param_new (const NcHaloMassSummaryMassDef mdef, const gdouble Delta);
NcHaloMCParam *nc_halo_mc_param_ref (NcHaloMCParam *hmcp);

void nc_halo_mc_param_free (NcHaloMCParam *hmcp);
void nc_halo_mc_param_clear (NcHaloMCParam **hmcp);

#define NC_HALO_MC_PARAM_DEFAULT_C_DELTA (4.0)
#define NC_HALO_MC_PARAM_DEFAULT_LOG10M_DELTA (log10 (2.0e14))
#define NC_HALO_MC_PARAM_DEFAULT_PARAMS_ABSTOL (0.0)

G_END_DECLS

#endif /* _NC_HALO_MC_PARAM_H_ */

