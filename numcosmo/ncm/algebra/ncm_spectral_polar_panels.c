/***************************************************************************
 *            ncm_spectral_polar_panels.c
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

/**
 * NcmSpectralPolarPanels:
 *
 * A function of a radius and an angle as Chebyshev panels.
 *
 * The radial edges $r_0 < r_1 < \dots < r_{n_r}$ and the angular edges
 * $\theta_0 < \theta_1 < \dots < \theta_{n_\theta}$ split
 * $[r_0, r_{n_r}] \times [\theta_0, \theta_{n_\theta}]$ into the $n_r n_\theta$ initial
 * panels $[r_i, r_{i+1}] \times [\theta_j, \theta_{j+1}]$. An initial panel holds one panel
 * or the panels that came from bisecting it, which cover it without overlap. On a panel
 * $[r_a, r_b] \times [\theta_a, \theta_b]$ with $K_r$ radial and $K_\theta$ angular
 * coefficients the function is
 * $$F (r, \theta) = \sum_{k = 0}^{K_r - 1} \sum_{l = 0}^{K_\theta - 1} c^{kl} \, T_k (s_r) \, T_l (s_\theta),$$
 * with $s_r = (2 r - r_a - r_b) / (r_b - r_a)$ and
 * $s_\theta = (2 \theta - \theta_a - \theta_b) / (\theta_b - \theta_a)$. The coefficients
 * of a panel are stored by angular index: row $l$ holds
 * $c^{0l}, c^{1l}, \dots, c^{(K_r - 1) l}$, and the numbers may differ between panels. A
 * panel is appended with $K_r = 2^{l_r} + 1$ and $K_\theta = 2^{l_\theta} + 1$ for its
 * levels $l_r$ and $l_\theta$; ncm_spectral_polar_panels_trim() then drops the trailing
 * coefficients that do not matter, which leaves any numbers.
 *
 * ncm_spectral_polar_build() fills the object from a function. The same can be done by
 * hand: ncm_spectral_polar_panels_set_edges(), then ncm_spectral_polar_panels_append() for
 * each panel, in the order of the initial panels (ascending $i n_\theta + j$) with the
 * panels of one initial panel one after the other.
 *
 * ncm_spectral_polar_panels_eval() finds the panel of a point and evaluates the sum
 * above. A point on an edge belongs to the panel that starts there. The range is not
 * checked: outside of it the series of the nearest panel is extrapolated. The angle is
 * not wrapped, so it must be given inside $[\theta_0, \theta_{n_\theta}]$.
 *
 * ncm_spectral_polar_panels_reset() keeps the arrays, so an object refilled at every
 * build stops allocating once it has grown. Reading and evaluating do not change the
 * object, so one object may be evaluated in several threads.
 *
 * See <a href="../../theory/ncm/algebra/spectral_polar.html">Spectral Integration on a
 * Disk</a> for the formulas.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/algebra/ncm_spectral_polar_panels.h"
#include "ncm/algebra/ncm_spectral_polar_panels_private.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <string.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_PANELS,
};

G_DEFINE_TYPE (NcmSpectralPolarPanels, ncm_spectral_polar_panels, G_TYPE_OBJECT)

static void
ncm_spectral_polar_panels_init (NcmSpectralPolarPanels *spp)
{
  spp->r_edges     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  spp->theta_edges = g_array_new (FALSE, FALSE, sizeof (gdouble));
  spp->first       = g_array_new (FALSE, FALSE, sizeof (guint));
  spp->panels      = g_array_new (FALSE, FALSE, sizeof (NcmSpectralPolarPanelsPanel));
  spp->coeffs      = g_array_new (FALSE, FALSE, sizeof (gdouble));
  spp->init_panels = NULL;
}

static void
_ncm_spectral_polar_panels_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSpectralPolarPanels *spp = NCM_SPECTRAL_POLAR_PANELS (object);

  g_return_if_fail (NCM_IS_SPECTRAL_POLAR_PANELS (object));

  switch (prop_id)
  {
    case PROP_PANELS:
      spp->init_panels = g_value_dup_variant (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static GVariant *_ncm_spectral_polar_panels_get_variant (NcmSpectralPolarPanels *spp);

static void
_ncm_spectral_polar_panels_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSpectralPolarPanels *spp = NCM_SPECTRAL_POLAR_PANELS (object);

  g_return_if_fail (NCM_IS_SPECTRAL_POLAR_PANELS (object));

  switch (prop_id)
  {
    case PROP_PANELS:
      g_value_take_variant (value, _ncm_spectral_polar_panels_get_variant (spp));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

/*
 * The edges and the panels as a (adada(dddduuad)), each panel with its numbers of radial
 * and angular coefficients. The arrays are copied: reset and append rewrite them in place,
 * and a variant must not change after it is made.
 */
static GVariant *
_ncm_spectral_polar_panels_get_variant (NcmSpectralPolarPanels *spp)
{
  GVariantBuilder builder;
  GVariant *r_edges, *theta_edges;
  guint i;

  r_edges     = g_variant_new_fixed_array (G_VARIANT_TYPE ("d"), spp->r_edges->data, spp->r_edges->len, sizeof (gdouble));
  theta_edges = g_variant_new_fixed_array (G_VARIANT_TYPE ("d"), spp->theta_edges->data, spp->theta_edges->len, sizeof (gdouble));

  g_variant_builder_init (&builder, G_VARIANT_TYPE ("a(dddduuad)"));

  for (i = 0; i < spp->panels->len; i++)
  {
    const NcmSpectralPolarPanelsPanel *panel = &g_array_index (spp->panels, NcmSpectralPolarPanelsPanel, i);
    const gsize len                          = (gsize) panel->n_r * panel->n_theta;
    GVariant *coeffs                         = g_variant_new_fixed_array (G_VARIANT_TYPE ("d"),
                                                                          &g_array_index (spp->coeffs, gdouble, panel->offset),
                                                                          len, sizeof (gdouble));

    g_variant_builder_add (&builder, "(dddduu@ad)", panel->ra, panel->rb, panel->ta, panel->tb,
                           panel->n_r, panel->n_theta, coeffs);
  }

  return g_variant_ref_sink (g_variant_new ("(@ad@ad@a(dddduuad))", r_edges, theta_edges, g_variant_builder_end (&builder)));
}

static void _ncm_spectral_polar_panels_set_edges_data (NcmSpectralPolarPanels *spp, const gdouble *r_edges, guint n_r_edges, const gdouble *theta_edges, guint n_theta_edges);
static void _ncm_spectral_polar_panels_append_sized (NcmSpectralPolarPanels *spp, gdouble ra, gdouble rb, gdouble ta, gdouble tb, guint n_r, guint n_theta, const gdouble *c);

/*
 * The edges and the panels given as a property, added in order so that each is checked as
 * an appended one. A panel holds its numbers of coefficients K_r and K_theta; a panel
 * written by an earlier version holds its levels instead, which is recognized by the
 * length of its coefficients: (2^a + 1) (2^b + 1) is never a b.
 */
static void
_ncm_spectral_polar_panels_constructed (GObject *object)
{
  NcmSpectralPolarPanels *spp = NCM_SPECTRAL_POLAR_PANELS (object);

  /* Chain up : start */
  G_OBJECT_CLASS (ncm_spectral_polar_panels_parent_class)->constructed (object);

  if (spp->init_panels != NULL)
  {
    GVariant *r_edges_var      = g_variant_get_child_value (spp->init_panels, 0);
    GVariant *theta_edges_var  = g_variant_get_child_value (spp->init_panels, 1);
    GVariant *panels_var       = g_variant_get_child_value (spp->init_panels, 2);
    const guint n_panels       = g_variant_n_children (panels_var);
    gsize n_r_edges            = 0;
    gsize n_theta_edges        = 0;
    const gdouble *r_edges     = g_variant_get_fixed_array (r_edges_var, &n_r_edges, sizeof (gdouble));
    const gdouble *theta_edges = g_variant_get_fixed_array (theta_edges_var, &n_theta_edges, sizeof (gdouble));
    guint i;

    if ((n_r_edges > 0) || (n_theta_edges > 0) || (n_panels > 0))
      _ncm_spectral_polar_panels_set_edges_data (spp, r_edges, n_r_edges, theta_edges, n_theta_edges);

    for (i = 0; i < n_panels; i++)
    {
      GVariant *panel_var = g_variant_get_child_value (panels_var, i);
      GVariant *coeffs_var;
      gdouble ra, rb, ta, tb;
      guint u_r, u_theta;
      const gdouble *c;
      gsize len = 0;

      g_variant_get (panel_var, "(dddduu@ad)", &ra, &rb, &ta, &tb, &u_r, &u_theta, &coeffs_var);
      c = g_variant_get_fixed_array (coeffs_var, &len, sizeof (gdouble));

      if ((u_r > 0) && (u_theta > 0) && (len == (gsize) u_r * u_theta))
        _ncm_spectral_polar_panels_append_sized (spp, ra, rb, ta, tb, u_r, u_theta, c);
      else if ((u_r <= NCM_SPECTRAL_POLAR_PANELS_LEVEL_MAX) && (u_theta <= NCM_SPECTRAL_POLAR_PANELS_LEVEL_MAX) &&
               (len == (gsize) ((1u << u_r) + 1) * ((1u << u_theta) + 1)))
        ncm_spectral_polar_panels_append_data (spp, ra, rb, ta, tb, u_r, u_theta, c);
      else
        g_error ("ncm_spectral_polar_panels: panel %u has %zu coefficients, neither %u by %u nor at the levels %u and %u.",
                 i, len, u_theta, u_r, u_r, u_theta);

      g_variant_unref (coeffs_var);
      g_variant_unref (panel_var);
    }

    g_variant_unref (r_edges_var);
    g_variant_unref (theta_edges_var);
    g_variant_unref (panels_var);
    g_clear_pointer (&spp->init_panels, g_variant_unref);
  }
}

static void
_ncm_spectral_polar_panels_finalize (GObject *object)
{
  NcmSpectralPolarPanels *spp = NCM_SPECTRAL_POLAR_PANELS (object);

  g_array_unref (spp->r_edges);
  g_array_unref (spp->theta_edges);
  g_array_unref (spp->first);
  g_array_unref (spp->panels);
  g_array_unref (spp->coeffs);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spectral_polar_panels_parent_class)->finalize (object);
}

static void
ncm_spectral_polar_panels_class_init (NcmSpectralPolarPanelsClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_spectral_polar_panels_set_property;
  object_class->get_property = &_ncm_spectral_polar_panels_get_property;
  object_class->constructed  = &_ncm_spectral_polar_panels_constructed;
  object_class->finalize     = &_ncm_spectral_polar_panels_finalize;

  /**
   * NcmSpectralPolarPanels:panels:
   *
   * The edges and the panels as a #GVariant of type `(adada(dddduuad))`: the radial edges,
   * the angular edges, and for each panel $r_a$, $r_b$, $\theta_a$, $\theta_b$, $K_r$,
   * $K_\theta$ and its $K_\theta K_r$ coefficients by angular index. The three arrays are
   * empty when no edges are set. A panel given with its levels $l_r$ and $l_\theta$ and
   * $(2^{l_\theta} + 1) (2^{l_r} + 1)$ coefficients, as earlier versions wrote it, is read
   * as well.
   */
  g_object_class_install_property (object_class,
                                   PROP_PANELS,
                                   g_param_spec_variant ("panels",
                                                         NULL,
                                                         "Edges, panels and coefficients",
                                                         G_VARIANT_TYPE ("(adada(dddduuad))"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_spectral_polar_panels_new:
 *
 * Creates a panel set with no edges and no panels.
 *
 * Returns: (transfer full): a new #NcmSpectralPolarPanels.
 */
NcmSpectralPolarPanels *
ncm_spectral_polar_panels_new (void)
{
  return g_object_new (NCM_TYPE_SPECTRAL_POLAR_PANELS, NULL);
}

/**
 * ncm_spectral_polar_panels_ref:
 * @spp: a #NcmSpectralPolarPanels
 *
 * Increases the reference count of @spp by one.
 *
 * Returns: (transfer full): @spp.
 */
NcmSpectralPolarPanels *
ncm_spectral_polar_panels_ref (NcmSpectralPolarPanels *spp)
{
  return g_object_ref (spp);
}

/**
 * ncm_spectral_polar_panels_free:
 * @spp: a #NcmSpectralPolarPanels
 *
 * Decreases the reference count of @spp by one.
 */
void
ncm_spectral_polar_panels_free (NcmSpectralPolarPanels *spp)
{
  g_object_unref (spp);
}

/**
 * ncm_spectral_polar_panels_clear:
 * @spp: a #NcmSpectralPolarPanels
 *
 * Decreases the reference count of *@spp by one and sets *@spp to %NULL.
 */
void
ncm_spectral_polar_panels_clear (NcmSpectralPolarPanels **spp)
{
  g_clear_object (spp);
}

/**
 * ncm_spectral_polar_panels_reset:
 * @spp: a #NcmSpectralPolarPanels
 *
 * Removes the edges and every panel, keeping the arrays for reuse, so that refilling @spp
 * allocates nothing until it holds more edges or more coefficients than before.
 */
void
ncm_spectral_polar_panels_reset (NcmSpectralPolarPanels *spp)
{
  g_array_set_size (spp->r_edges, 0);
  g_array_set_size (spp->theta_edges, 0);
  g_array_set_size (spp->first, 0);
  g_array_set_size (spp->panels, 0);
  g_array_set_size (spp->coeffs, 0);
}

/**
 * ncm_spectral_polar_panels_set_edges:
 * @spp: a #NcmSpectralPolarPanels
 * @r_edges: (element-type gdouble): the radial edges $r_0, \dots, r_{n_r}$
 * @theta_edges: (element-type gdouble): the angular edges $\theta_0, \dots, \theta_{n_\theta}$
 *
 * Resets @spp and sets the edges of its initial panels, copying both arrays. Each array
 * must be ascending with two entries or more; otherwise this is an error.
 */
void
ncm_spectral_polar_panels_set_edges (NcmSpectralPolarPanels *spp, GArray *r_edges, GArray *theta_edges)
{
  _ncm_spectral_polar_panels_set_edges_data (spp, (const gdouble *) r_edges->data, r_edges->len,
                                             (const gdouble *) theta_edges->data, theta_edges->len);
}

static void
_ncm_spectral_polar_panels_set_edges_data (NcmSpectralPolarPanels *spp, const gdouble *r_edges, guint n_r_edges, const gdouble *theta_edges, guint n_theta_edges)
{
  guint i;

  if ((n_r_edges < 2) || (n_theta_edges < 2))
    g_error ("ncm_spectral_polar_panels_set_edges: %u radial and %u angular edges, two or more of each are needed.",
             n_r_edges, n_theta_edges);

  for (i = 0; i + 1 < n_r_edges; i++)
  {
    if (!(r_edges[i] < r_edges[i + 1]))
      g_error ("ncm_spectral_polar_panels_set_edges: the radial edges are not ascending at index %u.", i);
  }

  for (i = 0; i + 1 < n_theta_edges; i++)
  {
    if (!(theta_edges[i] < theta_edges[i + 1]))
      g_error ("ncm_spectral_polar_panels_set_edges: the angular edges are not ascending at index %u.", i);
  }

  ncm_spectral_polar_panels_reset (spp);

  g_array_append_vals (spp->r_edges, r_edges, n_r_edges);
  g_array_append_vals (spp->theta_edges, theta_edges, n_theta_edges);
  g_array_set_size (spp->first, (n_r_edges - 1) * (n_theta_edges - 1) + 1);
  memset (spp->first->data, 0, sizeof (guint) * spp->first->len);
}

static gboolean _ncm_spectral_polar_panels_level (guint n_nodes, guint *level);

/**
 * ncm_spectral_polar_panels_append:
 * @spp: a #NcmSpectralPolarPanels
 * @ra: first radial edge of the panel
 * @rb: last radial edge of the panel
 * @ta: first angular edge of the panel
 * @tb: last angular edge of the panel
 * @coeffs: the $(2^{l_\theta} + 1)$ by $(2^{l_r} + 1)$ coefficients of the panel, row $l$ holding $c^{0l}, \dots, c^{2^{l_r} l}$
 *
 * Adds the panel $[r_a, r_b] \times [\theta_a, \theta_b]$ after the last one, copying
 * @coeffs. The edges must be set, the panel must be inside one initial panel, no panel of
 * a later initial panel may be present, and the numbers of rows and of columns of @coeffs
 * must each be a power of two plus one; otherwise this is an error. The panels of one
 * initial panel must cover it without overlap, which is not checked.
 */
void
ncm_spectral_polar_panels_append (NcmSpectralPolarPanels *spp, gdouble ra, gdouble rb, gdouble ta, gdouble tb, NcmMatrix *coeffs)
{
  const guint n_rows = ncm_matrix_nrows (coeffs);
  const guint n_cols = ncm_matrix_ncols (coeffs);
  guint l_r, l_theta;

  if (!_ncm_spectral_polar_panels_level (n_cols, &l_r) || !_ncm_spectral_polar_panels_level (n_rows, &l_theta))
    g_error ("ncm_spectral_polar_panels_append: the coefficients are %u by %u, each must be a power of two plus one.",
             n_rows, n_cols);

  if (ncm_matrix_tda (coeffs) == n_cols)
  {
    ncm_spectral_polar_panels_append_data (spp, ra, rb, ta, tb, l_r, l_theta, ncm_matrix_data (coeffs));
  }
  else
  {
    NcmMatrix *packed = ncm_matrix_dup (coeffs);

    ncm_spectral_polar_panels_append_data (spp, ra, rb, ta, tb, l_r, l_theta, ncm_matrix_data (packed));
    ncm_matrix_free (packed);
  }
}

/* Whether n_nodes is 2^level + 1 for a level up to the largest one, which is then set */
static gboolean
_ncm_spectral_polar_panels_level (guint n_nodes, guint *level)
{
  guint l;

  for (l = 0; l <= NCM_SPECTRAL_POLAR_PANELS_LEVEL_MAX; l++)
  {
    if (n_nodes == (1u << l) + 1)
    {
      *level = l;

      return TRUE;
    }
  }

  return FALSE;
}

/**
 * ncm_spectral_polar_panels_append_data: (skip)
 * @spp: a #NcmSpectralPolarPanels
 * @ra: first radial edge of the panel
 * @rb: last radial edge of the panel
 * @ta: first angular edge of the panel
 * @tb: last angular edge of the panel
 * @l_r: radial level of the panel
 * @l_theta: angular level of the panel
 * @c: the coefficients of the panel, $2^{l_\theta} + 1$ contiguous rows of $2^{l_r} + 1$
 *
 * As ncm_spectral_polar_panels_append(), for coefficients stored in a buffer.
 */
void
ncm_spectral_polar_panels_append_data (NcmSpectralPolarPanels *spp, gdouble ra, gdouble rb, gdouble ta, gdouble tb, guint l_r, guint l_theta, const gdouble *c)
{
  if ((l_r > NCM_SPECTRAL_POLAR_PANELS_LEVEL_MAX) || (l_theta > NCM_SPECTRAL_POLAR_PANELS_LEVEL_MAX))
    g_error ("ncm_spectral_polar_panels_append: the levels %u and %u are above %u.", l_r, l_theta,
             NCM_SPECTRAL_POLAR_PANELS_LEVEL_MAX);

  _ncm_spectral_polar_panels_append_sized (spp, ra, rb, ta, tb, (1u << l_r) + 1, (1u << l_theta) + 1, c);
}

/* Appends a panel with n_r radial and n_theta angular coefficients, after the checks of
 * ncm_spectral_polar_panels_append() */
static void
_ncm_spectral_polar_panels_append_sized (NcmSpectralPolarPanels *spp, gdouble ra, gdouble rb, gdouble ta, gdouble tb, guint n_r, guint n_theta, const gdouble *c)
{
  NcmSpectralPolarPanelsPanel panel = { ra, rb, ta, tb, n_r, n_theta, spp->coeffs->len };
  guint *first                      = (guint *) spp->first->data;
  const guint n_cells               = spp->first->len - 1;
  guint i, j, cell, q;
  gsize len;

  if (spp->r_edges->len < 2)
    g_error ("ncm_spectral_polar_panels_append: the edges are not set.");

  if (!(ra < rb) || !(ta < tb))
    g_error ("ncm_spectral_polar_panels_append: the panel [%.17g, %.17g] x [%.17g, %.17g] is empty.", ra, rb, ta, tb);

  i    = _ncm_spectral_polar_panels_edge_index ((const gdouble *) spp->r_edges->data, spp->r_edges->len, 0.5 * (ra + rb));
  j    = _ncm_spectral_polar_panels_edge_index ((const gdouble *) spp->theta_edges->data, spp->theta_edges->len, 0.5 * (ta + tb));
  cell = i * (spp->theta_edges->len - 1) + j;

  if ((ra < g_array_index (spp->r_edges, gdouble, i)) || (rb > g_array_index (spp->r_edges, gdouble, i + 1)) ||
      (ta < g_array_index (spp->theta_edges, gdouble, j)) || (tb > g_array_index (spp->theta_edges, gdouble, j + 1)))
    g_error ("ncm_spectral_polar_panels_append: the panel [%.17g, %.17g] x [%.17g, %.17g] is not inside one initial panel.",
             ra, rb, ta, tb);

  if (first[cell + 1] != spp->panels->len)
    g_error ("ncm_spectral_polar_panels_append: the panel [%.17g, %.17g] x [%.17g, %.17g] is in the initial panel %u, "
             "and a later initial panel has panels already.", ra, rb, ta, tb, cell);

  len = (gsize) n_r * n_theta;
  g_array_append_vals (spp->coeffs, c, len);
  g_array_append_val (spp->panels, panel);

  for (q = cell + 1; q <= n_cells; q++)
    first[q] = spp->panels->len;
}

/**
 * ncm_spectral_polar_panels_get_n_panels:
 * @spp: a #NcmSpectralPolarPanels
 *
 * Returns: the number of panels.
 */
guint
ncm_spectral_polar_panels_get_n_panels (NcmSpectralPolarPanels *spp)
{
  return spp->panels->len;
}

/**
 * ncm_spectral_polar_panels_get_range:
 * @spp: a #NcmSpectralPolarPanels
 * @ra: (out): the first radial edge
 * @rb: (out): the last radial edge
 * @ta: (out): the first angular edge
 * @tb: (out): the last angular edge
 *
 * Gets the range covered by the initial panels; the edges of @spp must be set.
 */
void
ncm_spectral_polar_panels_get_range (NcmSpectralPolarPanels *spp, gdouble *ra, gdouble *rb, gdouble *ta, gdouble *tb)
{
  g_assert_cmpuint (spp->r_edges->len, >, 1);

  *ra = g_array_index (spp->r_edges, gdouble, 0);
  *rb = g_array_index (spp->r_edges, gdouble, spp->r_edges->len - 1);
  *ta = g_array_index (spp->theta_edges, gdouble, 0);
  *tb = g_array_index (spp->theta_edges, gdouble, spp->theta_edges->len - 1);
}

/**
 * ncm_spectral_polar_panels_peek_panel: (skip)
 * @spp: a #NcmSpectralPolarPanels
 * @i: panel index
 * @ra: (out): first radial edge of the panel
 * @rb: (out): last radial edge of the panel
 * @ta: (out): first angular edge of the panel
 * @tb: (out): last angular edge of the panel
 * @n_r: (out): number of radial coefficients of the panel
 * @n_theta: (out): number of angular coefficients of the panel
 *
 * Gets panel @i, which must be below ncm_spectral_polar_panels_get_n_panels(). Its
 * @n_theta rows of @n_r coefficients are contiguous, and the returned pointer is valid
 * until @spp is changed.
 *
 * Returns: (transfer none): the coefficients of panel @i.
 */
const gdouble *
ncm_spectral_polar_panels_peek_panel (NcmSpectralPolarPanels *spp, guint i, gdouble *ra, gdouble *rb, gdouble *ta, gdouble *tb, guint *n_r, guint *n_theta)
{
  const NcmSpectralPolarPanelsPanel *panel;

  g_assert_cmpuint (i, <, spp->panels->len);

  panel = &g_array_index (spp->panels, NcmSpectralPolarPanelsPanel, i);

  *ra      = panel->ra;
  *rb      = panel->rb;
  *ta      = panel->ta;
  *tb      = panel->tb;
  *n_r     = panel->n_r;
  *n_theta = panel->n_theta;

  return &g_array_index (spp->coeffs, gdouble, panel->offset);
}

/**
 * ncm_spectral_polar_panels_get_panel:
 * @spp: a #NcmSpectralPolarPanels
 * @i: panel index
 * @ra: (out): first radial edge of the panel
 * @rb: (out): last radial edge of the panel
 * @ta: (out): first angular edge of the panel
 * @tb: (out): last angular edge of the panel
 *
 * Gets a copy of panel @i, which must be below ncm_spectral_polar_panels_get_n_panels().
 *
 * Returns: (transfer full): the $K_\theta$ by $K_r$ coefficients of panel @i.
 */
NcmMatrix *
ncm_spectral_polar_panels_get_panel (NcmSpectralPolarPanels *spp, guint i, gdouble *ra, gdouble *rb, gdouble *ta, gdouble *tb)
{
  guint n_r, n_t;
  const gdouble *c  = ncm_spectral_polar_panels_peek_panel (spp, i, ra, rb, ta, tb, &n_r, &n_t);
  NcmMatrix *coeffs = ncm_matrix_new (n_t, n_r);

  memcpy (ncm_matrix_data (coeffs), c, sizeof (gdouble) * n_t * n_r);

  return coeffs;
}

/*
 * The leading block of the n_theta by n_r coefficients c with the fewest coefficients
 * whose complement has l1 norm at most tol: *keep_theta rows of *keep_r columns. sum is
 * scratch of n_theta n_r entries, filled with the l1 norms of the leading blocks.
 */
static void
_ncm_spectral_polar_panels_trim_block (const gdouble *c, guint n_r, guint n_theta, gdouble tol, gdouble *sum, guint *keep_r, guint *keep_theta)
{
  gdouble total;
  guint best_r = n_r;
  guint best_t = n_theta;
  guint k, l;

  /* sum[l n_r + k]: l1 norm of rows 0..l, columns 0..k */
  for (l = 0; l < n_theta; l++)
  {
    gdouble row = 0.0;

    for (k = 0; k < n_r; k++)
    {
      row             += fabs (c[(gsize) l * n_r + k]);
      sum[l * n_r + k] = row + ((l > 0) ? sum[(l - 1) * n_r + k] : 0.0);
    }
  }

  total = sum[(n_theta - 1) * n_r + (n_r - 1)];

  for (l = 0; l < n_theta; l++)
  {
    for (k = 0; k < n_r; k++)
    {
      if (total - sum[l * n_r + k] <= tol)
      {
        if ((gsize) (l + 1) * (k + 1) < (gsize) best_t * best_r)
        {
          best_t = l + 1;
          best_r = k + 1;
        }

        break;
      }
    }
  }

  *keep_r     = best_r;
  *keep_theta = best_t;
}

/**
 * ncm_spectral_polar_panels_trim:
 * @spp: a #NcmSpectralPolarPanels
 * @tol: largest l1 norm of the coefficients dropped from a panel
 *
 * Cuts each panel to the leading block of its coefficients, the lowest radial and angular
 * indices, with the fewest coefficients among those whose dropped coefficients have l1
 * norm at most @tol. Since $|T_k| \le 1$ on the panel, the function changes by at most
 * @tol at any point of it. The coefficients of the panels that remain are packed into an
 * array of their size, so the object takes less memory and its evaluation fewer
 * operations; the serialized object keeps the trimmed panels. Objects serialized before
 * panels could be trimmed (each panel given by its levels) are read as they were, untrimmed.
 */
void
ncm_spectral_polar_panels_trim (NcmSpectralPolarPanels *spp, gdouble tol)
{
  GArray *scratch = g_array_new (FALSE, FALSE, sizeof (gdouble));
  gdouble *c      = (gdouble *) spp->coeffs->data;
  gsize cursor    = 0;
  guint i;

  g_assert_cmpfloat (tol, >=, 0.0);

  for (i = 0; i < spp->panels->len; i++)
  {
    NcmSpectralPolarPanelsPanel *panel = &g_array_index (spp->panels, NcmSpectralPolarPanelsPanel, i);
    const gdouble *src                 = c + panel->offset;
    guint keep_r, keep_t, l;

    g_array_set_size (scratch, (gsize) panel->n_r * panel->n_theta);
    _ncm_spectral_polar_panels_trim_block (src, panel->n_r, panel->n_theta, tol, (gdouble *) scratch->data, &keep_r, &keep_t);

    /* The block moves to cursor <= offset with a row stride keep_r <= n_r, so every row
     * goes to an address at or below the one it is read from */
    for (l = 0; l < keep_t; l++)
      memmove (c + cursor + (gsize) l * keep_r, src + (gsize) l * panel->n_r, sizeof (gdouble) * keep_r);

    panel->offset  = cursor;
    panel->n_r     = keep_r;
    panel->n_theta = keep_t;
    cursor        += (gsize) keep_r * keep_t;
  }

  /* A new array of the trimmed size: g_array_set_size() would keep the allocation */
  {
    GArray *packed = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), cursor);

    g_array_append_vals (packed, spp->coeffs->data, cursor);
    g_array_unref (spp->coeffs);
    spp->coeffs = packed;
  }

  g_array_unref (scratch);
}

/* Index sort of the l1 norms, ascending (insertion: a few tens of sets) */
static void
_ncm_spectral_polar_panels_sort_l1 (const gdouble *l1, guint n, guint *idx)
{
  guint i, j;

  for (i = 0; i < n; i++)
    idx[i] = i;

  for (i = 1; i < n; i++)
  {
    const guint k = idx[i];

    for (j = i; (j > 0) && (l1[idx[j - 1]] > l1[k]); j--)
      idx[j] = idx[j - 1];

    idx[j] = k;
  }
}

/**
 * ncm_spectral_polar_panels_trim_modes: (skip)
 * @n: the number of panel sets
 * @spp: (array length=n): the panel sets, with the same panels and, in each panel, the same numbers of coefficients
 * @tol: the largest change of $\sum_q w_q F_q$ at any point, for any weights with $|w_q| \le 1$
 *
 * Cuts the @n sets, the terms $F_q$ of an expansion combined with weights of modulus at most
 * one (ncm_spectral_polar_panels_combine()), so that the combination changes by at most @tol
 * at any point. In each panel the terms with the smallest l1 norms are dropped while the sum
 * of their norms stays within @tol / 2 (each keeps one coefficient, zero), and each term kept
 * is cut to its leading block as in ncm_spectral_polar_panels_trim(), with an equal share of
 * the other half. The sets keep their panels; their numbers of coefficients then differ.
 */
void
ncm_spectral_polar_panels_trim_modes (guint n, NcmSpectralPolarPanels **spp, gdouble tol)
{
  GArray *scratch = g_array_new (FALSE, FALSE, sizeof (gdouble));
  gdouble *l1     = g_new (gdouble, n);
  guint *idx      = g_new (guint, n);
  guint *keep_r   = g_new (guint, n);
  guint *keep_t   = g_new (guint, n);
  gsize *cursor   = g_new0 (gsize, n);
  guint i, q;

  g_assert_cmpuint (n, >, 0);
  g_assert_cmpfloat (tol, >=, 0.0);

  for (q = 1; q < n; q++)
    g_assert_cmpuint (spp[q]->panels->len, ==, spp[0]->panels->len);

  for (i = 0; i < spp[0]->panels->len; i++)
  {
    const NcmSpectralPolarPanelsPanel *p0 = &g_array_index (spp[0]->panels, NcmSpectralPolarPanelsPanel, i);
    const guint n_r                       = p0->n_r;
    const guint n_theta                   = p0->n_theta;
    const gsize len                       = (gsize) n_r * n_theta;
    gdouble dropped                       = 0.0;
    guint n_drop                          = 0;
    guint k, l;

    g_array_set_size (scratch, len);

    for (q = 0; q < n; q++)
    {
      const NcmSpectralPolarPanelsPanel *pq = &g_array_index (spp[q]->panels, NcmSpectralPolarPanelsPanel, i);
      const gdouble *c                      = (const gdouble *) spp[q]->coeffs->data + pq->offset;

      g_assert_cmpuint (pq->n_r, ==, n_r);
      g_assert_cmpuint (pq->n_theta, ==, n_theta);

      l1[q] = 0.0;

      for (k = 0; k < len; k++)
        l1[q] += fabs (c[k]);
    }

    _ncm_spectral_polar_panels_sort_l1 (l1, n, idx);

    while ((n_drop < n) && (dropped + l1[idx[n_drop]] <= 0.5 * tol))
    {
      dropped += l1[idx[n_drop]];
      n_drop++;
    }

    for (k = 0; k < n; k++)
    {
      const guint qq = idx[k];

      if (k < n_drop)
      {
        keep_r[qq] = 0;
        keep_t[qq] = 0;
      }
      else
      {
        const NcmSpectralPolarPanelsPanel *pq = &g_array_index (spp[qq]->panels, NcmSpectralPolarPanelsPanel, i);

        _ncm_spectral_polar_panels_trim_block ((const gdouble *) spp[qq]->coeffs->data + pq->offset, n_r, n_theta,
                                               0.5 * tol / (n - n_drop), (gdouble *) scratch->data, &keep_r[qq], &keep_t[qq]);
      }
    }

    for (q = 0; q < n; q++)
    {
      NcmSpectralPolarPanelsPanel *pq = &g_array_index (spp[q]->panels, NcmSpectralPolarPanelsPanel, i);
      gdouble *c                      = (gdouble *) spp[q]->coeffs->data;
      const gsize offset              = pq->offset;

      if (keep_r[q] == 0)
      {
        /* A dropped term keeps one coefficient, zero */
        c[cursor[q]] = 0.0;
        pq->n_r      = 1;
        pq->n_theta  = 1;
        pq->offset   = cursor[q];
        cursor[q]   += 1;
        continue;
      }

      /* As in ncm_spectral_polar_panels_trim(): every row moves to an address at or below its own */
      for (l = 0; l < keep_t[q]; l++)
        memmove (c + cursor[q] + (gsize) l * keep_r[q], c + offset + (gsize) l * n_r, sizeof (gdouble) * keep_r[q]);

      pq->offset  = cursor[q];
      pq->n_r     = keep_r[q];
      pq->n_theta = keep_t[q];
      cursor[q]  += (gsize) keep_r[q] * keep_t[q];
    }
  }

  for (q = 0; q < n; q++)
    g_array_set_size (spp[q]->coeffs, cursor[q]);

  g_array_unref (scratch);
  g_free (l1);
  g_free (idx);
  g_free (keep_r);
  g_free (keep_t);
  g_free (cursor);
}

/**
 * ncm_spectral_polar_panels_combine: (skip)
 * @out: the panel set to fill
 * @n: the number of panel sets
 * @spp: (array length=n): the panel sets, with the same panels
 * @w: (array length=n): the weights
 *
 * Sets @out to the expansion of $\sum_q w_q F_q$, with @spp[q] that of $F_q$: the panels of
 * the sets, each with the largest numbers of radial and angular coefficients of the sets
 * there, and $\sum_q w_q c_q$ with the coefficients of each set added to the leading block
 * they fill. The layout of @out is made only when it differs, so that a set combined again
 * with other weights only writes its coefficients.
 */
void
ncm_spectral_polar_panels_combine (NcmSpectralPolarPanels *out, guint n, NcmSpectralPolarPanels **spp, const gdouble *w)
{
  const guint n_panels = spp[0]->panels->len;
  gboolean same_layout = (out->panels->len == n_panels) && (out->r_edges->len == spp[0]->r_edges->len) &&
                         (out->theta_edges->len == spp[0]->theta_edges->len);
  gdouble *c;
  guint i, q;

  g_assert_cmpuint (n, >, 0);

  for (q = 1; q < n; q++)
    g_assert_cmpuint (spp[q]->panels->len, ==, n_panels);

  /* The largest block of each panel */
  for (i = 0; same_layout && (i < n_panels); i++)
  {
    const NcmSpectralPolarPanelsPanel *po = &g_array_index (out->panels, NcmSpectralPolarPanelsPanel, i);
    guint n_r                             = 0;
    guint n_theta                         = 0;

    for (q = 0; q < n; q++)
    {
      const NcmSpectralPolarPanelsPanel *pq = &g_array_index (spp[q]->panels, NcmSpectralPolarPanelsPanel, i);

      n_r     = MAX (n_r, pq->n_r);
      n_theta = MAX (n_theta, pq->n_theta);
    }

    same_layout = (po->n_r == n_r) && (po->n_theta == n_theta);
  }

  if (!same_layout)
  {
    gsize cursor = 0;

    g_array_set_size (out->r_edges, 0);
    g_array_append_vals (out->r_edges, spp[0]->r_edges->data, spp[0]->r_edges->len);
    g_array_set_size (out->theta_edges, 0);
    g_array_append_vals (out->theta_edges, spp[0]->theta_edges->data, spp[0]->theta_edges->len);
    g_array_set_size (out->first, 0);
    g_array_append_vals (out->first, spp[0]->first->data, spp[0]->first->len);
    g_array_set_size (out->panels, 0);
    g_array_append_vals (out->panels, spp[0]->panels->data, n_panels);

    for (i = 0; i < n_panels; i++)
    {
      NcmSpectralPolarPanelsPanel *po = &g_array_index (out->panels, NcmSpectralPolarPanelsPanel, i);

      po->n_r     = 0;
      po->n_theta = 0;

      for (q = 0; q < n; q++)
      {
        const NcmSpectralPolarPanelsPanel *pq = &g_array_index (spp[q]->panels, NcmSpectralPolarPanelsPanel, i);

        po->n_r     = MAX (po->n_r, pq->n_r);
        po->n_theta = MAX (po->n_theta, pq->n_theta);
      }

      po->offset = cursor;
      cursor    += (gsize) po->n_r * po->n_theta;
    }

    g_array_set_size (out->coeffs, cursor);
  }

  c = (gdouble *) out->coeffs->data;
  memset (c, 0, sizeof (gdouble) * out->coeffs->len);

  for (i = 0; i < n_panels; i++)
  {
    const NcmSpectralPolarPanelsPanel *po = &g_array_index (out->panels, NcmSpectralPolarPanelsPanel, i);
    gdouble *co                           = c + po->offset;

    for (q = 0; q < n; q++)
    {
      const NcmSpectralPolarPanelsPanel *pq = &g_array_index (spp[q]->panels, NcmSpectralPolarPanelsPanel, i);
      const gdouble *cq                     = (const gdouble *) spp[q]->coeffs->data + pq->offset;
      const gdouble wq                      = w[q];
      guint l, k;

      if (wq == 0.0)
        continue;

      for (l = 0; l < pq->n_theta; l++)
      {
        gdouble *row_o       = co + (gsize) l * po->n_r;
        const gdouble *row_q = cq + (gsize) l * pq->n_r;

        #pragma omp simd

        for (k = 0; k < pq->n_r; k++)
          row_o[k] += wq * row_q[k];
      }
    }
  }
}

/**
 * ncm_spectral_polar_panels_peek_edges: (skip)
 * @spp: a #NcmSpectralPolarPanels
 * @r_edges: (out) (transfer none): the radial edges of the initial panels
 * @theta_edges: (out) (transfer none): their angular edges
 *
 * Peeks the edges set by ncm_spectral_polar_panels_set_edges().
 */
void
ncm_spectral_polar_panels_peek_edges (NcmSpectralPolarPanels *spp, GArray **r_edges, GArray **theta_edges)
{
  *r_edges     = spp->r_edges;
  *theta_edges = spp->theta_edges;
}

/**
 * ncm_spectral_polar_panels_find_panel:
 * @spp: a #NcmSpectralPolarPanels
 * @r: a radius
 * @theta: an angle
 *
 * Finds the panel holding $(r, \theta)$: the initial panel by bisection over the radial
 * and over the angular edges, then a search over the panels of that initial panel, which
 * must have at least one. A point on an edge belongs to the panel that starts there, and a
 * point outside the range to the nearest panel.
 *
 * Returns: the index of the panel.
 */
guint
ncm_spectral_polar_panels_find_panel (NcmSpectralPolarPanels *spp, gdouble r, gdouble theta)
{
  return _ncm_spectral_polar_panels_find_panel_in ((const gdouble *) spp->r_edges->data, spp->r_edges->len,
                                                   (const gdouble *) spp->theta_edges->data, spp->theta_edges->len,
                                                   (const guint *) spp->first->data,
                                                   (const NcmSpectralPolarPanelsPanel *) spp->panels->data, r, theta);
}

/**
 * ncm_spectral_polar_panels_find_panel_near:
 * @spp: a #NcmSpectralPolarPanels
 * @hint: a panel index, usually the result of the previous search
 * @r: a radius
 * @theta: an angle
 *
 * Finds the panel holding $(r, \theta)$, the same panel as
 * ncm_spectral_polar_panels_find_panel(), first checking whether the point is inside
 * panel @hint, $r_a \le r < r_b$ and $\theta_a \le \theta < \theta_b$. Points evaluated
 * in sequence along a curve usually fall in the same panel, and then the bisections are
 * skipped. Any @hint is accepted: one that is not a panel index or that does not hold the
 * point leads to the full search.
 *
 * Returns: the index of the panel.
 */
guint
ncm_spectral_polar_panels_find_panel_near (NcmSpectralPolarPanels *spp, guint hint, gdouble r, gdouble theta)
{
  if (hint < spp->panels->len)
  {
    const NcmSpectralPolarPanelsPanel *panel = &g_array_index (spp->panels, NcmSpectralPolarPanelsPanel, hint);

    /* Inside the box of the panel the full search returns it: the panels of an initial
     * panel tile it, and the rules for the sides only apply outside it. */
    if ((r >= panel->ra) && (r < panel->rb) && (theta >= panel->ta) && (theta < panel->tb))
      return hint;
  }

  return ncm_spectral_polar_panels_find_panel (spp, r, theta);
}

/**
 * ncm_spectral_polar_panels_eval:
 * @spp: a #NcmSpectralPolarPanels
 * @r: a radius
 * @theta: an angle
 *
 * Evaluates the function at $(r, \theta)$ on the panel given by
 * ncm_spectral_polar_panels_find_panel(). Outside the range of the edges the series of the
 * nearest panel is extrapolated.
 *
 * Returns: $F (r, \theta)$.
 */
gdouble
ncm_spectral_polar_panels_eval (NcmSpectralPolarPanels *spp, gdouble r, gdouble theta)
{
  return ncm_spectral_polar_panels_eval_on (spp, ncm_spectral_polar_panels_find_panel (spp, r, theta), r, theta);
}

/**
 * ncm_spectral_polar_panels_eval_on_array: (skip)
 * @spp: a #NcmSpectralPolarPanels
 * @i: panel index
 * @n: number of points
 * @r: (array length=n): the radii
 * @theta: (array length=n): the angles
 * @out: (array length=n): the values
 *
 * Evaluates the series of panel @i at the @n points ($r_j$, $\theta_j$), the values of
 * ncm_spectral_polar_panels_eval_on() at each point within rounding: the two do the same
 * operations in the same order but may contract them into fma differently, so they agree
 * within rounding of the sum of the moduli of the coefficients (1e-13 of it in the tests),
 * not bit for bit. The points are taken several at a time, each coefficient being loaded
 * once for all of them, so that evaluating points that share a panel together costs less
 * than evaluating them one by one. @out may be @r or @theta. @i must be below
 * ncm_spectral_polar_panels_get_n_panels(), which is not checked.
 */
void
ncm_spectral_polar_panels_eval_on_array (NcmSpectralPolarPanels *spp, guint i, guint n, const gdouble *r, const gdouble *theta, gdouble *out)
{
  const NcmSpectralPolarPanelsPanel *panel = &g_array_index (spp->panels, NcmSpectralPolarPanelsPanel, i);
  const gdouble *c                         = &g_array_index (spp->coeffs, gdouble, panel->offset);

  _ncm_spectral_polar_panels_eval_points (panel, c, n, r, theta, out);
}

/**
 * ncm_spectral_polar_panels_eval_on:
 * @spp: a #NcmSpectralPolarPanels
 * @i: panel index
 * @r: a radius
 * @theta: an angle
 *
 * Evaluates the series of panel @i at $(r, \theta)$ by Clenshaw recurrences, the radial
 * one for each angular index and the angular one over their results. @i must be below
 * ncm_spectral_polar_panels_get_n_panels(), which is not checked.
 *
 * Returns: the series of panel @i at $(r, \theta)$.
 */
gdouble
ncm_spectral_polar_panels_eval_on (NcmSpectralPolarPanels *spp, guint i, gdouble r, gdouble theta)
{
  const NcmSpectralPolarPanelsPanel *panel = &g_array_index (spp->panels, NcmSpectralPolarPanelsPanel, i);

  return _ncm_spectral_polar_panels_eval_on_c (panel, &g_array_index (spp->coeffs, gdouble, panel->offset), r, theta);
}

