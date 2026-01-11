/***************************************************************************
 *            ncm_sbessel_integrator_fftl.h
 *
 *  Thu January 09 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_sbessel_integrator_fftl.h
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

#ifndef _NCM_SBESSEL_INTEGRATOR_FFTL_H_
#define _NCM_SBESSEL_INTEGRATOR_FFTL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_sbessel_integrator.h>

G_BEGIN_DECLS

#define NCM_TYPE_SBESSEL_INTEGRATOR_FFTL (ncm_sbessel_integrator_fftl_get_type ())

G_DECLARE_FINAL_TYPE (NcmSBesselIntegratorFFTL, ncm_sbessel_integrator_fftl, NCM, SBESSEL_INTEGRATOR_FFTL, NcmSBesselIntegrator)

NcmSBesselIntegratorFFTL *ncm_sbessel_integrator_fftl_new (guint lmin, guint lmax);
NcmSBesselIntegratorFFTL *ncm_sbessel_integrator_fftl_ref (NcmSBesselIntegratorFFTL *sbilf);
void ncm_sbessel_integrator_fftl_free (NcmSBesselIntegratorFFTL *sbilf);
void ncm_sbessel_integrator_fftl_clear (NcmSBesselIntegratorFFTL **sbilf);

G_END_DECLS

#endif /* _NCM_SBESSEL_INTEGRATOR_FFTL_H_ */
