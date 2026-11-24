/***************************************************************************
 *            ncm_spectral_tensor.h
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

#ifndef _NCM_SPECTRAL_TENSOR_H_
#define _NCM_SPECTRAL_TENSOR_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPECTRAL_TENSOR (ncm_spectral_tensor_get_type ())
#define NCM_SPECTRAL_TENSOR_MAX_VARS (16)

G_DECLARE_FINAL_TYPE (NcmSpectralTensor, ncm_spectral_tensor, NCM, SPECTRAL_TENSOR, GObject)

/**
 * NcmSpectralTensorNodeKind:
 * @NCM_SPECTRAL_TENSOR_NODE_KIND_EMPTY: an empty node, not divided and with no series
 * @NCM_SPECTRAL_TENSOR_NODE_KIND_SPLIT: a split node, divided along one variable
 * @NCM_SPECTRAL_TENSOR_NODE_KIND_PANEL: a panel node, not divided and with a series
 *
 * The kinds of node of the tree of an #NcmSpectralTensor. Only a split node has children.
 */
typedef enum _NcmSpectralTensorNodeKind /*< prefix=NCM_SPECTRAL_TENSOR_NODE_KIND >*/
{
  NCM_SPECTRAL_TENSOR_NODE_KIND_EMPTY,
  NCM_SPECTRAL_TENSOR_NODE_KIND_SPLIT,
  NCM_SPECTRAL_TENSOR_NODE_KIND_PANEL,
} NcmSpectralTensorNodeKind;

/**
 * NcmSpectralTensorPanelFormat:
 * @NCM_SPECTRAL_TENSOR_PANEL_FORMAT_DENSE: the panel stores its coefficients
 * @NCM_SPECTRAL_TENSOR_PANEL_FORMAT_TUCKER: the panel stores the factors and the core of a
 *   Tucker decomposition of its coefficients
 *
 * The storage formats of a panel of an #NcmSpectralTensor.
 */
typedef enum _NcmSpectralTensorPanelFormat /*< prefix=NCM_SPECTRAL_TENSOR_PANEL_FORMAT >*/
{
  NCM_SPECTRAL_TENSOR_PANEL_FORMAT_DENSE,
  NCM_SPECTRAL_TENSOR_PANEL_FORMAT_TUCKER,
} NcmSpectralTensorPanelFormat;

NcmSpectralTensor *ncm_spectral_tensor_new (guint n_vars);
NcmSpectralTensor *ncm_spectral_tensor_ref (NcmSpectralTensor *st);
void ncm_spectral_tensor_free (NcmSpectralTensor *st);
void ncm_spectral_tensor_clear (NcmSpectralTensor **st);

void ncm_spectral_tensor_set_domain (NcmSpectralTensor *st, guint n_lower, const gdouble *lower, guint n_upper, const gdouble *upper);
void ncm_spectral_tensor_reset (NcmSpectralTensor *st);
guint ncm_spectral_tensor_split (NcmSpectralTensor *st, guint node, guint i, guint n_edges, const gdouble *edges);
guint ncm_spectral_tensor_split_empty (NcmSpectralTensor *st, guint i, guint n_edges, const gdouble *edges);
guint ncm_spectral_tensor_set_panel (NcmSpectralTensor *st, guint node, guint n_K, const guint *n_coeffs, gsize len, const gdouble *coeffs);
guint ncm_spectral_tensor_set_panel_tucker (NcmSpectralTensor *st, guint node, guint n_K, const guint *n_coeffs, guint n_r, const guint *ranks, gsize len, const gdouble *data);

guint ncm_spectral_tensor_get_n_vars (NcmSpectralTensor *st);
guint ncm_spectral_tensor_get_n_nodes (NcmSpectralTensor *st);
guint ncm_spectral_tensor_get_n_empty (NcmSpectralTensor *st);
guint ncm_spectral_tensor_get_n_panels (NcmSpectralTensor *st);
void ncm_spectral_tensor_get_range (NcmSpectralTensor *st, guint i, gdouble *a, gdouble *b);

NcmSpectralTensorNodeKind ncm_spectral_tensor_get_node_kind (NcmSpectralTensor *st, guint node);
void ncm_spectral_tensor_peek_node_bounds (NcmSpectralTensor *st, guint node, const gdouble **lower, const gdouble **upper);
void ncm_spectral_tensor_get_node_bounds (NcmSpectralTensor *st, guint node, GArray **lower, GArray **upper);
GArray *ncm_spectral_tensor_get_node_split (NcmSpectralTensor *st, guint node, guint *i, guint *first_child);
guint ncm_spectral_tensor_get_node_panel (NcmSpectralTensor *st, guint node);

const gdouble *ncm_spectral_tensor_peek_panel (NcmSpectralTensor *st, guint p, const gdouble **lower, const gdouble **upper, const guint **n_coeffs);
GArray *ncm_spectral_tensor_get_panel (NcmSpectralTensor *st, guint p, GArray **lower, GArray **upper, GArray **n_coeffs);
NcmSpectralTensorPanelFormat ncm_spectral_tensor_get_panel_format (NcmSpectralTensor *st, guint p);
GArray *ncm_spectral_tensor_get_panel_tucker (NcmSpectralTensor *st, guint p, GArray **ranks);

guint ncm_spectral_tensor_find_panel (NcmSpectralTensor *st, guint n_vars, const gdouble *x);
gdouble ncm_spectral_tensor_eval (NcmSpectralTensor *st, guint n_vars, const gdouble *x);
gdouble ncm_spectral_tensor_eval_on (NcmSpectralTensor *st, guint p, guint n_vars, const gdouble *x);
void ncm_spectral_tensor_eval_array (NcmSpectralTensor *st, GArray *x, GArray **fval);

G_END_DECLS

#endif /* _NCM_SPECTRAL_TENSOR_H_ */

