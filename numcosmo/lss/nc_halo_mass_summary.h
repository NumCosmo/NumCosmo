/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_halo_mass_summary.h
 *
 *  Thu Oct 10 14:12:27 2024
 *  Copyright  2024  Mariana Penna-Lima
 *  <pennalima@unb.br>
 ****************************************************************************/
/*
 * nc_halo_mass_summary.h
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
#ifndef _NC_HALO_MASS_SUMMARY_H_
#define _NC_HALO_MASS_SUMMARY_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/nc_hicosmo.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_MASS_SUMMARY (nc_halo_mass_summary_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcHaloMassSummary, nc_halo_mass_summary, NC, HALO_MASS_SUMMARY, NcmModel)

struct _NcHaloMassSummaryClass
{
  /*< private >*/
  NcmModelClass parent_class;

  gdouble (*mass) (NcHaloMassSummary *hms);
  gdouble (*concentration) (NcHaloMassSummary *hms, NcHICosmo *cosmo);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[15];
};

/**
 * NcHaloMassSummaryMassDef:
 * @NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN: halo mass defined in terms of the mean density $\rho_\mathrm{bg} = \rho_m(z)$
 * @NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL: halo mass defined in terms of the critical density $\rho_\mathrm{bg} = \rho_\mathrm{crit}(z)$
 * @NC_HALO_MASS_SUMMARY_MASS_DEF_VIRIAL: halo mass defined in terms of virial overdensity times the critical density $\rho_\mathrm{bg} = \rho_\mathrm{crit
 *
 * Spherical overdensity halo mass: $$M_\Delta = \frac{4\pi}{3} \Delta \rho_\mathrm{bg} r_\Delta^3,$$
 * where $\rho_\mathrm{bg}$ is the background density of the universe at redshift z, $\rho_\mathrm{bg} (z)$.
 * For @NC_HALO_MASS_SUMMARY_MASS_DEF_VIRIAL the virial overdensity is defined as:
 * \begin{equation}\label{def:DVir}
 * \Delta_\mathrm{Vir} = 18 \pi^2 + 82 x - 39 x^2, \quad x \equiv \Omega_m(z) - 1.
 * \end{equation}
 *
 */
typedef enum _NcHaloMassSummaryMassDef
{
  NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN = 0,
  NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL,
  NC_HALO_MASS_SUMMARY_MASS_DEF_VIRIAL,
  /* < private > */
  NC_HALO_MASS_SUMMARY_MASS_DEF_LEN, /*< skip >*/
} NcHaloMassSummaryMassDef;


NCM_MSET_MODEL_DECLARE_ID (nc_halo_mass_summary);

NcHaloMassSummary *nc_halo_mass_summary_ref (NcHaloMassSummary *hms);

void nc_halo_mass_summary_free (NcHaloMassSummary *hms);
void nc_halo_mass_summary_clear (NcHaloMassSummary **hms);

void nc_halo_mass_summary_set_Delta (NcHaloMassSummary *hms, const gdouble Delta);
gdouble nc_halo_mass_summary_get_Delta (NcHaloMassSummary *hms);

gdouble nc_halo_mass_summary_mass (NcHaloMassSummary *hms);
gdouble nc_halo_mass_summary_concentration (NcHaloMassSummary *hms, NcHICosmo *cosmo);

gdouble nc_halo_mass_summary_Delta (NcHaloMassSummary *hms, NcHICosmo *cosmo, const gdouble z);
gdouble nc_halo_mass_summary_rho_bg (NcHaloMassSummary *hms, NcHICosmo *cosmo, const gdouble z);
gdouble nc_halo_mass_summary_Delta_rho_bg (NcHaloMassSummary *hms, NcHICosmo *cosmo, const gdouble z);

G_END_DECLS

#endif /* _NC_HALO_MASS_SUMMARY_H_ */

