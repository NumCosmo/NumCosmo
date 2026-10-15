/***************************************************************************
 *            ncm_spectral_panels_product.h
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

#ifndef _NCM_SPECTRAL_PANELS_PRODUCT_H_
#define _NCM_SPECTRAL_PANELS_PRODUCT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/ncm/algebra/ncm_spectral_panels.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPECTRAL_PANELS_PRODUCT (ncm_spectral_panels_product_get_type ())

G_DECLARE_FINAL_TYPE (NcmSpectralPanelsProduct, ncm_spectral_panels_product, NCM, SPECTRAL_PANELS_PRODUCT, GObject)

NcmSpectralPanelsProduct *ncm_spectral_panels_product_new (guint p);
NcmSpectralPanelsProduct *ncm_spectral_panels_product_ref (NcmSpectralPanelsProduct *prod);

void ncm_spectral_panels_product_free (NcmSpectralPanelsProduct *prod);
void ncm_spectral_panels_product_clear (NcmSpectralPanelsProduct **prod);

guint ncm_spectral_panels_product_get_p (NcmSpectralPanelsProduct *prod);

void ncm_spectral_panels_product_integrate (NcmSpectralPanelsProduct *prod, NcmSpectralPanels *P, NcmVector *s1, NcmSpectralPanels *Q, NcmVector *s2, NcmVector *result, NcmVector *err);
void ncm_spectral_panels_product_integrate_rows (NcmSpectralPanelsProduct *prod, NcmSpectralPanels *P, guint o1, NcmVector *s1, NcmSpectralPanels *Q, guint o2, NcmVector *s2, NcmVector *result, NcmVector *err);

G_END_DECLS

#endif /* _NCM_SPECTRAL_PANELS_PRODUCT_H_ */

