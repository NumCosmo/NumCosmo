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
#include <numcosmo/ncm/specfunc/ncm_sbessel_integrator.h>

G_BEGIN_DECLS

#define NCM_TYPE_SBESSEL_INTEGRATOR_LEVIN (ncm_sbessel_integrator_levin_get_type ())

G_DECLARE_FINAL_TYPE (NcmSBesselIntegratorLevin, ncm_sbessel_integrator_levin, NCM, SBESSEL_INTEGRATOR_LEVIN, NcmSBesselIntegrator)

NcmSBesselIntegratorLevin *ncm_sbessel_integrator_levin_new (guint ell_min, guint ell_max);
NcmSBesselIntegratorLevin *ncm_sbessel_integrator_levin_new_full (guint ell_min, guint ell_max, gdouble y_knots_min, gdouble y_knots_max, guint n_knots, guint ell_cache_max, gdouble reltol, guint cheb_min_order, gdouble cheb_reltol);
NcmSBesselIntegratorLevin *ncm_sbessel_integrator_levin_ref (NcmSBesselIntegratorLevin *sbilv);
void ncm_sbessel_integrator_levin_free (NcmSBesselIntegratorLevin *sbilv);
void ncm_sbessel_integrator_levin_clear (NcmSBesselIntegratorLevin **sbilv);

void ncm_sbessel_integrator_levin_set_max_order (NcmSBesselIntegratorLevin *sbilv, guint max_order);
guint ncm_sbessel_integrator_levin_get_max_order (NcmSBesselIntegratorLevin *sbilv);

void ncm_sbessel_integrator_levin_set_reltol (NcmSBesselIntegratorLevin *sbilv, gdouble reltol);
gdouble ncm_sbessel_integrator_levin_get_reltol (NcmSBesselIntegratorLevin *sbilv);

void ncm_sbessel_integrator_levin_set_cheb_min_order (NcmSBesselIntegratorLevin *sbilv, guint cheb_min_order);
guint ncm_sbessel_integrator_levin_get_cheb_min_order (NcmSBesselIntegratorLevin *sbilv);

void ncm_sbessel_integrator_levin_set_cheb_reltol (NcmSBesselIntegratorLevin *sbilv, gdouble cheb_reltol);
gdouble ncm_sbessel_integrator_levin_get_cheb_reltol (NcmSBesselIntegratorLevin *sbilv);

gdouble ncm_sbessel_integrator_levin_get_y_knots_min (NcmSBesselIntegratorLevin *sbilv);
gdouble ncm_sbessel_integrator_levin_get_y_knots_max (NcmSBesselIntegratorLevin *sbilv);
guint ncm_sbessel_integrator_levin_get_n_knots (NcmSBesselIntegratorLevin *sbilv);
guint ncm_sbessel_integrator_levin_get_ell_cache_max (NcmSBesselIntegratorLevin *sbilv);

#define NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_Y_KNOTS_MIN (1.0e-4)
#define NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_Y_KNOTS_MAX (1.0e6)
#define NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_N_KNOTS (21)
#define NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_ELL_CACHE_MAX (1200)
#define NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_RELTOL (1.0e-13)
#define NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_CHEB_MIN_ORDER (2)
#define NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_CHEB_RELTOL (1.0e-8)


G_END_DECLS

#endif /* _NCM_SBESSEL_INTEGRATOR_LEVIN_H_ */

