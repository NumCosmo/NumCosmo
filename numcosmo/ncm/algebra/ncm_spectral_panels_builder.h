/***************************************************************************
 *            ncm_spectral_panels_builder.h
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

#ifndef _NCM_SPECTRAL_PANELS_BUILDER_H_
#define _NCM_SPECTRAL_PANELS_BUILDER_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/ncm/algebra/ncm_matrix.h>
#include <numcosmo/ncm/algebra/ncm_spectral_panels.h>

G_BEGIN_DECLS

/**
 * NcmSpectralPanelsF:
 * @user_data: (closure): user data
 * @i: index of the continuous interval holding @x
 * @x: the point
 * @fval: $F(x)$, filled in place
 *
 * Vector-valued function expanded in panels. The breakpoints split the range into
 * continuous intervals, numbered from 0; a panel lies in one of them, and $F$ is called
 * with its index, so at a breakpoint it returns the limit from the side of that panel.
 */
typedef void (*NcmSpectralPanelsF) (gpointer user_data, guint i, gdouble x, NcmVector *fval);

/**
 * NcmSpectralPanelsAccept:
 * @user_data: (closure): user data
 * @a: left edge of the panel
 * @b: right edge of the panel
 * @c: the coefficients of the panel, rows of @N
 * @N: number of coefficients in each row
 * @tda: distance between the starts of consecutive rows of @c
 *
 * Further acceptance test of a converged panel.
 *
 * Returns: whether the panel is accepted.
 */
typedef gboolean (*NcmSpectralPanelsAccept) (gpointer user_data, gdouble a, gdouble b, const gdouble *c, guint N, guint tda);

#define NCM_TYPE_SPECTRAL_PANELS_BUILDER (ncm_spectral_panels_builder_get_type ())

G_DECLARE_FINAL_TYPE (NcmSpectralPanelsBuilder, ncm_spectral_panels_builder, NCM, SPECTRAL_PANELS_BUILDER, GObject)

NcmSpectralPanelsBuilder *ncm_spectral_panels_builder_new (void);
NcmSpectralPanelsBuilder *ncm_spectral_panels_builder_ref (NcmSpectralPanelsBuilder *builder);

void ncm_spectral_panels_builder_free (NcmSpectralPanelsBuilder *builder);
void ncm_spectral_panels_builder_clear (NcmSpectralPanelsBuilder **builder);

void ncm_spectral_panels_builder_set_level_min (NcmSpectralPanelsBuilder *builder, guint level_min);
void ncm_spectral_panels_builder_set_level_cap (NcmSpectralPanelsBuilder *builder, guint level_cap);
void ncm_spectral_panels_builder_set_min_panels_per_efold (NcmSpectralPanelsBuilder *builder, gdouble min_panels_per_efold);
void ncm_spectral_panels_builder_set_min_panel_fraction (NcmSpectralPanelsBuilder *builder, gdouble min_panel_fraction);

guint ncm_spectral_panels_builder_get_level_min (NcmSpectralPanelsBuilder *builder);
guint ncm_spectral_panels_builder_get_level_cap (NcmSpectralPanelsBuilder *builder);
gdouble ncm_spectral_panels_builder_get_min_panels_per_efold (NcmSpectralPanelsBuilder *builder);
gdouble ncm_spectral_panels_builder_get_min_panel_fraction (NcmSpectralPanelsBuilder *builder);

void ncm_spectral_panels_builder_build (NcmSpectralPanelsBuilder *builder, NcmSpectralPanelsF F, gpointer user_data, gdouble a, gdouble b, NcmVector *breakpoints, NcmMatrix *f_minus, NcmMatrix *f_plus, gdouble reltol, gdouble abstol, NcmSpectralPanels *sp);
void ncm_spectral_panels_builder_build_full (NcmSpectralPanelsBuilder *builder, NcmSpectralPanelsF F, NcmSpectralPanelsAccept accept, gpointer user_data, gdouble a, gdouble b, NcmVector *breakpoints, NcmMatrix *f_minus, NcmMatrix *f_plus, gdouble reltol, gdouble abstol, NcmSpectralPanels *sp);

G_END_DECLS

#endif /* _NCM_SPECTRAL_PANELS_BUILDER_H_ */

