/***************************************************************************
 *            ncm_spectral_panels.h
 *
 *  Sat Oct 10 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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

#ifndef _NCM_SPECTRAL_PANELS_H_
#define _NCM_SPECTRAL_PANELS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/ncm/algebra/ncm_matrix.h>
#include <numcosmo/ncm/algebra/ncm_spectral.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPECTRAL_PANELS (ncm_spectral_panels_get_type ())

G_DECLARE_FINAL_TYPE (NcmSpectralPanels, ncm_spectral_panels, NCM, SPECTRAL_PANELS, GObject)

NcmSpectralPanels *ncm_spectral_panels_new (guint fdim);
NcmSpectralPanels *ncm_spectral_panels_ref (NcmSpectralPanels *sp);

void ncm_spectral_panels_free (NcmSpectralPanels *sp);
void ncm_spectral_panels_clear (NcmSpectralPanels **sp);

void ncm_spectral_panels_append (NcmSpectralPanels *sp, gdouble a, gdouble b, NcmMatrix *coeffs);
void ncm_spectral_panels_reset (NcmSpectralPanels *sp);

guint ncm_spectral_panels_get_fdim (NcmSpectralPanels *sp);
guint ncm_spectral_panels_get_n_panels (NcmSpectralPanels *sp);
void ncm_spectral_panels_get_range (NcmSpectralPanels *sp, gdouble *a, gdouble *b);
const gdouble *ncm_spectral_panels_peek_panel (NcmSpectralPanels *sp, guint i, gdouble *a, gdouble *b, guint *N);
NcmMatrix *ncm_spectral_panels_get_panel (NcmSpectralPanels *sp, guint i, gdouble *a, gdouble *b);
guint ncm_spectral_panels_find_panel (NcmSpectralPanels *sp, gdouble x);

void ncm_spectral_panels_eval (NcmSpectralPanels *sp, gdouble x, NcmVector *fval);
gdouble ncm_spectral_panels_eval_row (NcmSpectralPanels *sp, guint r, gdouble x);
gboolean ncm_spectral_panels_restrict (NcmSpectralPanels *sp, NcmSpectral *spectral, gdouble a, gdouble b, NcmMatrix **coeffs);

G_END_DECLS

#endif /* _NCM_SPECTRAL_PANELS_H_ */

