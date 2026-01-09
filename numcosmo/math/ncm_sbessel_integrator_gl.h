/***************************************************************************
 *            ncm_sbessel_integrator_gl.h
 *
 *  Thu January 09 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_sbessel_integrator_gl.h
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

#ifndef _NCM_SBESSEL_INTEGRATOR_GL_H_
#define _NCM_SBESSEL_INTEGRATOR_GL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_sbessel_integrator.h>

G_BEGIN_DECLS

#define NCM_TYPE_SBESSEL_INTEGRATOR_GL (ncm_sbessel_integrator_gl_get_type ())

G_DECLARE_FINAL_TYPE (NcmSBesselIntegratorGL, ncm_sbessel_integrator_gl, NCM, SBESSEL_INTEGRATOR_GL, NcmSBesselIntegrator)

NcmSBesselIntegratorGL *ncm_sbessel_integrator_gl_new (guint lmin, guint lmax);
NcmSBesselIntegratorGL *ncm_sbessel_integrator_gl_ref (NcmSBesselIntegratorGL *sbigl);
void ncm_sbessel_integrator_gl_free (NcmSBesselIntegratorGL *sbigl);
void ncm_sbessel_integrator_gl_clear (NcmSBesselIntegratorGL **sbigl);

void ncm_sbessel_integrator_gl_set_npts (NcmSBesselIntegratorGL *sbigl, guint npts);
guint ncm_sbessel_integrator_gl_get_npts (NcmSBesselIntegratorGL *sbigl);
void ncm_sbessel_integrator_gl_set_margin (NcmSBesselIntegratorGL *sbigl, gdouble margin);
gdouble ncm_sbessel_integrator_gl_get_margin (NcmSBesselIntegratorGL *sbigl);
void ncm_sbessel_integrator_gl_set_nosc (NcmSBesselIntegratorGL *sbigl, gdouble nosc);
gdouble ncm_sbessel_integrator_gl_get_nosc (NcmSBesselIntegratorGL *sbigl);

G_END_DECLS

#endif /* _NCM_SBESSEL_INTEGRATOR_GL_H_ */

