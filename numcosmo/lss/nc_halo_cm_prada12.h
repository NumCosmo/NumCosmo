/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_halo_cm_prada12.h
 *
 *  Mon Dec 16 15:01:55 2024
 *  Copyright  2024  Mariana Penna-Lima <pennalima@unb.br>, Thais Mikami Ornellas <thais.ornellas@uel.br> 
 ****************************************************************************/
/*
 * nc_halo_cm_prada12.h
 * Copyright (C) 2024 Mariana Penna-Lima <pennalima@unb.br>, Thais Mikami Ornellas <thais.ornellas@uel.br>
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

#ifndef _NC_HALO_CM_PRADA12_H_
#define _NC_HALO_CM_PRADA12_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_halo_mass_summary.h>
#include <numcosmo/math/ncm_powspec_filter.h>
#include <numcosmo/lss/nc_halo_mass_function.h>
#include <numcosmo/model/nc_hicosmo_de.h>


G_BEGIN_DECLS

#define NC_TYPE_HALO_CM_PRADA12 (nc_halo_cm_prada12_get_type ())

G_DECLARE_FINAL_TYPE (NcHaloCMPrada12, nc_halo_cm_prada12, NC, HALO_CM_PRADA12, NcHaloMassSummary)

/**
 * NcHaloCMPrada12SParams:
 * @NC_HALO_CM_PRADA12_LOG10M_DELTA: halo mass $\log_{10}(M_\Delta)$
 *
 * Fundamental parametrization of the profile $\rho(r)$.
 * The halo mass is a paremeter while the concentration is given by the 
 * Prada et al. (2012) concentration-mass relation.
 *
 */
typedef enum /*< enum,underscore_name=NC_HALO_CM_PRADA12_SPARAMS >*/
{
  NC_HALO_CM_PRADA12_LOG10M_DELTA,
  /* < private > */
  NC_HALO_CM_PRADA12_SPARAM_LEN, /*< skip >*/
} NcHaloCMPrada12SParams;

#define NC_HALO_CM_PRADA12_LOCAL_SPARAM_LEN (NC_HALO_CM_PRADA12_SPARAM_LEN - 0)

NcHaloCMPrada12 *nc_halo_cm_prada12_new (const NcHaloMassSummaryMassDef mdef, const gdouble Delta);
NcHaloCMPrada12 *nc_halo_cm_prada12_ref (NcHaloCMPrada12 *hcmp);

void nc_halo_cm_prada12_free (NcHaloCMPrada12 *hcmp);
void nc_halo_cm_prada12_clear (NcHaloCMPrada12 **hcmp);

#define NC_HALO_CM_PRADA12_DEFAULT_LOG10M_DELTA (log10 (2.0e14))
#define NC_HALO_CM_PRADA12_DEFAULT_PARAMS_ABSTOL (0.0)

G_END_DECLS

#endif /* _NC_HALO_CM_PRADA12_H_ */

