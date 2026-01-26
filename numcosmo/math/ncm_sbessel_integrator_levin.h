/***************************************************************************
 *            ncm_sbessel_integrator_levin.h
 *
 *  Sat January 25 00:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_sbessel_integrator_levin.h
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

#ifndef _NCM_SBESSEL_INTEGRATOR_LEVIN_H_
#define _NCM_SBESSEL_INTEGRATOR_LEVIN_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_sbessel_integrator.h>

G_BEGIN_DECLS

#define NCM_TYPE_SBESSEL_INTEGRATOR_LEVIN (ncm_sbessel_integrator_levin_get_type ())

G_DECLARE_FINAL_TYPE (NcmSBesselIntegratorLevin, ncm_sbessel_integrator_levin, NCM, SBESSEL_INTEGRATOR_LEVIN, NcmSBesselIntegrator)

NcmSBesselIntegratorLevin *ncm_sbessel_integrator_levin_new (guint lmin, guint lmax);
NcmSBesselIntegratorLevin *ncm_sbessel_integrator_levin_ref (NcmSBesselIntegratorLevin *sbilv);
void ncm_sbessel_integrator_levin_free (NcmSBesselIntegratorLevin *sbilv);
void ncm_sbessel_integrator_levin_clear (NcmSBesselIntegratorLevin **sbilv);

void ncm_sbessel_integrator_levin_set_n_panels (NcmSBesselIntegratorLevin *sbilv, guint n_panels);
guint ncm_sbessel_integrator_levin_get_n_panels (NcmSBesselIntegratorLevin *sbilv);

void ncm_sbessel_integrator_levin_set_min_order (NcmSBesselIntegratorLevin *sbilv, guint min_order);
guint ncm_sbessel_integrator_levin_get_min_order (NcmSBesselIntegratorLevin *sbilv);

void ncm_sbessel_integrator_levin_set_max_order (NcmSBesselIntegratorLevin *sbilv, guint max_order);
guint ncm_sbessel_integrator_levin_get_max_order (NcmSBesselIntegratorLevin *sbilv);

void ncm_sbessel_integrator_levin_set_reltol (NcmSBesselIntegratorLevin *sbilv, gdouble reltol);
gdouble ncm_sbessel_integrator_levin_get_reltol (NcmSBesselIntegratorLevin *sbilv);

G_END_DECLS

#endif /* _NCM_SBESSEL_INTEGRATOR_LEVIN_H_ */

