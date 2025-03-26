/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_halo_cm_duffy08.h
 *
 *  Thu Dec 05 09:42:15 2024
 *  Copyright  2024  Mariana Penna-Lima <pennalima@unb.br>, Thais Mikami Ornellas <thais.ornellas@uel.br> 
 ****************************************************************************/
/*
 * nc_halo_cm_duffy08.h
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

#ifndef _NC_HALO_CM_DUFFY08_H_
#define _NC_HALO_CM_DUFFY08_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_halo_mass_summary.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_CM_DUFFY08 (nc_halo_cm_duffy08_get_type ())

G_DECLARE_FINAL_TYPE (NcHaloCMDuffy08, nc_halo_cm_duffy08, NC, HALO_CM_DUFFY08, NcHaloMassSummary)

/**
 * NcHaloCMDuffy08SParams:
 * @NC_HALO_CM_DUFFY08_LOG10M_DELTA: halo mass $\log_{10}(M_\Delta)$
 *
 * Fundamental parametrization of the profile $\rho(r)$.
 * The halo mass is a paremeter while the concentration is given by the 
 * Duffy et al. (2008) concentration-mass relation.
 *
 */
typedef enum /*< enum,underscore_name=NC_HALO_CM_DUFFY08_SPARAMS >*/
{
  NC_HALO_CM_DUFFY08_LOG10M_DELTA,
  /* < private > */
  NC_HALO_CM_DUFFY08_SPARAM_LEN, /*< skip >*/
} NcHaloCMDuffy08SParams;

#define NC_HALO_CM_DUFFY08_LOCAL_SPARAM_LEN (NC_HALO_CM_DUFFY08_SPARAM_LEN - 0)

NcHaloCMDuffy08 *nc_halo_cm_duffy08_new (const NcHaloMassSummaryMassDef mdef, const gdouble Delta);
NcHaloCMDuffy08 *nc_halo_cm_duffy08_ref (NcHaloCMDuffy08 *hcmd);

void nc_halo_cm_duffy08_free (NcHaloCMDuffy08 *hcmd);
void nc_halo_cm_duffy08_clear (NcHaloCMDuffy08 **hcmd);

#define NC_HALO_CM_DUFFY08_DEFAULT_LOG10M_DELTA (log10 (2.0e14))
#define NC_HALO_CM_DUFFY08_DEFAULT_PARAMS_ABSTOL (0.0)

G_END_DECLS

#endif /* _NC_HALO_CM_DUFFY08_H_ */

