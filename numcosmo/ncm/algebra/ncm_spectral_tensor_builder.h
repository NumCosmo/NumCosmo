/***************************************************************************
 *            ncm_spectral_tensor_builder.h
 *
 *  Thu Nov 19 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
 * Copyright (C) Caio Lima de Oliveira 2026 <caiolimadeoliveira@pm.me>
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

#ifndef _NCM_SPECTRAL_TENSOR_BUILDER_H_
#define _NCM_SPECTRAL_TENSOR_BUILDER_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_spectral_tensor.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPECTRAL_TENSOR_BUILDER (ncm_spectral_tensor_builder_get_type ())
#define NCM_SPECTRAL_TENSOR_BUILDER_ERROR (ncm_spectral_tensor_builder_error_quark ())

G_DECLARE_FINAL_TYPE (NcmSpectralTensorBuilder, ncm_spectral_tensor_builder, NCM, SPECTRAL_TENSOR_BUILDER, GObject)

/**
 * NcmSpectralTensorBuilderError:
 * @NCM_SPECTRAL_TENSOR_BUILDER_ERROR_MAX_DEPTH: a split beyond #NcmSpectralTensorBuilder:max-depth
 *
 * Errors of ncm_spectral_tensor_builder_build().
 */
typedef enum _NcmSpectralTensorBuilderError /*< prefix=NCM_SPECTRAL_TENSOR_BUILDER_ERROR >*/
{
  NCM_SPECTRAL_TENSOR_BUILDER_ERROR_MAX_DEPTH,
} NcmSpectralTensorBuilderError;

/**
 * NcmSpectralTensorF:
 * @n_vars: number of variables
 * @n: the number of points of each variable, @n_vars entries
 * @x: the points of each variable, @x[i] holding @n[i] values
 * @fval: the values of $F$ at the $\prod_i n_i$ points of the tensor product of the points
 *   of all variables, the last variable running fastest, filled in place
 * @user_data: user data
 *
 * The function a #NcmSpectralTensorBuilder tabulates, evaluated on a tensor product of
 * points. The arrays are valid during the call only.
 */
typedef void (*NcmSpectralTensorF) (guint n_vars, const guint *n, const gdouble * const *x, gdouble *fval, gpointer user_data);

GQuark ncm_spectral_tensor_builder_error_quark (void);

NcmSpectralTensorBuilder *ncm_spectral_tensor_builder_new (void);
NcmSpectralTensorBuilder *ncm_spectral_tensor_builder_ref (NcmSpectralTensorBuilder *builder);
void ncm_spectral_tensor_builder_free (NcmSpectralTensorBuilder *builder);
void ncm_spectral_tensor_builder_clear (NcmSpectralTensorBuilder **builder);

void ncm_spectral_tensor_builder_set_level_min (NcmSpectralTensorBuilder *builder, guint level_min);
void ncm_spectral_tensor_builder_set_level_cap (NcmSpectralTensorBuilder *builder, guint level_cap);
void ncm_spectral_tensor_builder_set_max_depth (NcmSpectralTensorBuilder *builder, guint max_depth);
void ncm_spectral_tensor_builder_set_trim_frac (NcmSpectralTensorBuilder *builder, gdouble trim_frac);
guint ncm_spectral_tensor_builder_get_level_min (NcmSpectralTensorBuilder *builder);
guint ncm_spectral_tensor_builder_get_level_cap (NcmSpectralTensorBuilder *builder);
guint ncm_spectral_tensor_builder_get_max_depth (NcmSpectralTensorBuilder *builder);
gdouble ncm_spectral_tensor_builder_get_trim_frac (NcmSpectralTensorBuilder *builder);

void ncm_spectral_tensor_builder_build (NcmSpectralTensorBuilder *builder, NcmSpectralTensorF F, NcmSpectralTensor *st, const guint *level_min, const guint *level_cap, gdouble reltol, gdouble abstol, gpointer user_data, guint64 *n_eval, GError **error);

G_END_DECLS

#endif /* _NCM_SPECTRAL_TENSOR_BUILDER_H_ */

