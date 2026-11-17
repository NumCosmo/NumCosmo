/***************************************************************************
 *            ncm_spectral_polar_panels.h
 *
 *  Wed Nov 11 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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

#ifndef _NCM_SPECTRAL_POLAR_PANELS_H_
#define _NCM_SPECTRAL_POLAR_PANELS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_matrix.h>
#include <numcosmo/ncm/core/ncm_obj_array.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPECTRAL_POLAR_PANELS (ncm_spectral_polar_panels_get_type ())
#define NCM_TYPE_SPECTRAL_POLAR_PANELS_ARRAY (ncm_spectral_polar_panels_array_get_type ())

G_DECLARE_FINAL_TYPE (NcmSpectralPolarPanels, ncm_spectral_polar_panels, NCM, SPECTRAL_POLAR_PANELS, GObject)
G_DECLARE_FINAL_TYPE (NcmSpectralPolarPanelsArray, ncm_spectral_polar_panels_array, NCM, SPECTRAL_POLAR_PANELS_ARRAY, GObject)

NcmSpectralPolarPanels *ncm_spectral_polar_panels_new (void);
NcmSpectralPolarPanels *ncm_spectral_polar_panels_ref (NcmSpectralPolarPanels *spp);

void ncm_spectral_polar_panels_free (NcmSpectralPolarPanels *spp);
void ncm_spectral_polar_panels_clear (NcmSpectralPolarPanels **spp);

void ncm_spectral_polar_panels_reset (NcmSpectralPolarPanels *spp);
void ncm_spectral_polar_panels_set_edges (NcmSpectralPolarPanels *spp, GArray *r_edges, GArray *theta_edges);
void ncm_spectral_polar_panels_append (NcmSpectralPolarPanels *spp, gdouble ra, gdouble rb, gdouble ta, gdouble tb, NcmMatrix *coeffs);
void ncm_spectral_polar_panels_append_data (NcmSpectralPolarPanels *spp, gdouble ra, gdouble rb, gdouble ta, gdouble tb, guint l_r, guint l_theta, const gdouble *c);

guint ncm_spectral_polar_panels_get_n_panels (NcmSpectralPolarPanels *spp);
void ncm_spectral_polar_panels_get_range (NcmSpectralPolarPanels *spp, gdouble *ra, gdouble *rb, gdouble *ta, gdouble *tb);
const gdouble *ncm_spectral_polar_panels_peek_panel (NcmSpectralPolarPanels *spp, guint i, gdouble *ra, gdouble *rb, gdouble *ta, gdouble *tb, guint *n_r, guint *n_theta);
void ncm_spectral_polar_panels_trim (NcmSpectralPolarPanels *spp, gdouble tol);
void ncm_spectral_polar_panels_trim_modes (guint n, NcmSpectralPolarPanels **spp, gdouble tol);
void ncm_spectral_polar_panels_combine (NcmSpectralPolarPanels *out, guint n, NcmSpectralPolarPanels **spp, const gdouble *w);
void ncm_spectral_polar_panels_peek_edges (NcmSpectralPolarPanels *spp, GArray **r_edges, GArray **theta_edges);
NcmMatrix *ncm_spectral_polar_panels_get_panel (NcmSpectralPolarPanels *spp, guint i, gdouble *ra, gdouble *rb, gdouble *ta, gdouble *tb);
guint ncm_spectral_polar_panels_find_panel (NcmSpectralPolarPanels *spp, gdouble r, gdouble theta);
guint ncm_spectral_polar_panels_find_panel_near (NcmSpectralPolarPanels *spp, guint hint, gdouble r, gdouble theta);

gdouble ncm_spectral_polar_panels_eval (NcmSpectralPolarPanels *spp, gdouble r, gdouble theta);
gdouble ncm_spectral_polar_panels_eval_on (NcmSpectralPolarPanels *spp, guint i, gdouble r, gdouble theta);
void ncm_spectral_polar_panels_eval_on_array (NcmSpectralPolarPanels *spp, guint i, guint n, const gdouble *r, const gdouble *theta, gdouble *out);

typedef struct _NcmSpectralPolarPanelsComb NcmSpectralPolarPanelsComb;

NcmSpectralPolarPanelsComb *ncm_spectral_polar_panels_comb_new (void);
void ncm_spectral_polar_panels_comb_free (NcmSpectralPolarPanelsComb *comb);
void ncm_spectral_polar_panels_comb_reset (NcmSpectralPolarPanelsComb *comb);
void ncm_spectral_polar_panels_eval_on_array_combine (guint n_sets, NcmSpectralPolarPanels **spp, const gdouble *w, NcmSpectralPolarPanelsComb *comb, guint i, guint n, const gdouble *r, const gdouble *theta, gdouble *out);

typedef struct _NcmSpectralPolarPanelsPack NcmSpectralPolarPanelsPack;

gsize ncm_spectral_polar_panels_pack_get_size (const NcmSpectralPolarPanelsPack *pack);
void ncm_spectral_polar_panels_pack_get_range (const NcmSpectralPolarPanelsPack *pack, gdouble *ra, gdouble *rb, gdouble *ta, gdouble *tb);
guint ncm_spectral_polar_panels_pack_find_panel_near (const NcmSpectralPolarPanelsPack *pack, guint hint, gdouble r, gdouble theta);
void ncm_spectral_polar_panels_pack_eval_on_array_combine (const NcmSpectralPolarPanelsPack *pack, const gdouble *w, NcmSpectralPolarPanelsComb *comb, guint i, guint n, const gdouble *r, const gdouble *theta, gdouble *out);

NcmSpectralPolarPanelsArray *ncm_spectral_polar_panels_array_new (void);
NcmSpectralPolarPanelsArray *ncm_spectral_polar_panels_array_ref (NcmSpectralPolarPanelsArray *spa);

void ncm_spectral_polar_panels_array_free (NcmSpectralPolarPanelsArray *spa);
void ncm_spectral_polar_panels_array_clear (NcmSpectralPolarPanelsArray **spa);

void ncm_spectral_polar_panels_array_add (NcmSpectralPolarPanelsArray *spa, NcmSpectralPolarPanels *spp);
guint ncm_spectral_polar_panels_array_len (NcmSpectralPolarPanelsArray *spa);
NcmSpectralPolarPanels *ncm_spectral_polar_panels_array_peek (NcmSpectralPolarPanelsArray *spa, guint i);
NcmSpectralPolarPanels **ncm_spectral_polar_panels_array_peek_data (NcmSpectralPolarPanelsArray *spa);
void ncm_spectral_polar_panels_array_pack (NcmSpectralPolarPanelsArray *spa);
const NcmSpectralPolarPanelsPack *ncm_spectral_polar_panels_array_peek_pack (NcmSpectralPolarPanelsArray *spa);

G_END_DECLS

#endif /* _NCM_SPECTRAL_POLAR_PANELS_H_ */

