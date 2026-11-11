/***************************************************************************
 *            ncm_spectral_polar_panels.c
 *
 *  Wed Nov 11 2026
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
 * $[r_a, r_b] \times [\theta_a, \theta_b]$ of levels $l_r$ and $l_\theta$ the function is
 * $$F (r, \theta) = \sum_{k = 0}^{2^{l_r}} \sum_{l = 0}^{2^{l_\theta}} c^{kl} \, T_k (s_r) \, T_l (s_\theta),$$
 * with $s_r = (2 r - r_a - r_b) / (r_b - r_a)$ and
 * $s_\theta = (2 \theta - \theta_a - \theta_b) / (\theta_b - \theta_a)$. The coefficients
 * of a panel are stored by angular index: row $l$ holds
 * $c^{0l}, c^{1l}, \dots, c^{2^{l_r} l}$, and the levels may differ between panels.
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
 * Disc</a> for the formulas.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/algebra/ncm_spectral_polar_panels.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <string.h>
#endif /* NUMCOSMO_GIR_SCAN */

/* Largest level of a panel in either direction */
#define NCM_SPECTRAL_POLAR_PANELS_LEVEL_MAX (20)

enum
{
  PROP_0,
  PROP_PANELS,
};

/* A panel, its levels and the position of its coefficients in the coefficient array */
typedef struct _NcmSpectralPolarPanelsPanel
{
  gdouble ra;
  gdouble rb;
  gdouble ta;
  gdouble tb;
  guint l_r;
  guint l_theta;
  gsize offset;
} NcmSpectralPolarPanelsPanel;

struct _NcmSpectralPolarPanels
{
  /*< private >*/
  GObject parent_instance;
  GArray *r_edges;       /* gdouble, the radial edges of the initial panels */
  GArray *theta_edges;   /* gdouble, their angular edges */
  GArray *first;         /* guint, n_r n_theta + 1 entries: initial panel c holds the panels
                          * first[c], ..., first[c + 1] - 1 */
  GArray *panels;        /* NcmSpectralPolarPanelsPanel */
  GArray *coeffs;        /* gdouble, the coefficients of every panel one after the other */
  GVariant *init_panels; /* the construct-only property, consumed by constructed() */
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
 * The edges and the panels as a (adada(dddduuad)). The arrays are copied: reset and append
 * rewrite them in place, and a variant must not change after it is made.
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
    const gsize len                          = (gsize) ((1u << panel->l_r) + 1) * ((1u << panel->l_theta) + 1);
    GVariant *coeffs                         = g_variant_new_fixed_array (G_VARIANT_TYPE ("d"),
                                                                          &g_array_index (spp->coeffs, gdouble, panel->offset),
                                                                          len, sizeof (gdouble));

    g_variant_builder_add (&builder, "(dddduu@ad)", panel->ra, panel->rb, panel->ta, panel->tb,
                           panel->l_r, panel->l_theta, coeffs);
  }

  return g_variant_ref_sink (g_variant_new ("(@ad@ad@a(dddduuad))", r_edges, theta_edges, g_variant_builder_end (&builder)));
}

static void _ncm_spectral_polar_panels_set_edges_data (NcmSpectralPolarPanels *spp, const gdouble *r_edges, guint n_r_edges, const gdouble *theta_edges, guint n_theta_edges);

/* The edges and the panels given as a property, added in order so that each is checked as an appended one */
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
      guint l_r, l_theta;
      const gdouble *c;
      gsize len = 0;

      g_variant_get (panel_var, "(dddduu@ad)", &ra, &rb, &ta, &tb, &l_r, &l_theta, &coeffs_var);
      c = g_variant_get_fixed_array (coeffs_var, &len, sizeof (gdouble));

      if ((l_r > NCM_SPECTRAL_POLAR_PANELS_LEVEL_MAX) || (l_theta > NCM_SPECTRAL_POLAR_PANELS_LEVEL_MAX) ||
          (len != (gsize) ((1u << l_r) + 1) * ((1u << l_theta) + 1)))
        g_error ("ncm_spectral_polar_panels: panel %u has %zu coefficients at the levels %u and %u.", i, len, l_r, l_theta);

      ncm_spectral_polar_panels_append_data (spp, ra, rb, ta, tb, l_r, l_theta, c);
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
   * the angular edges, and for each panel $r_a$, $r_b$, $\theta_a$, $\theta_b$, $l_r$,
   * $l_\theta$ and its coefficients by angular index. The three arrays are empty when no
   * edges are set.
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

static guint _ncm_spectral_polar_panels_edge_index (const GArray *edges, gdouble x);

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
  NcmSpectralPolarPanelsPanel panel = { ra, rb, ta, tb, l_r, l_theta, spp->coeffs->len };
  guint *first                      = (guint *) spp->first->data;
  const guint n_cells               = spp->first->len - 1;
  guint i, j, cell, q;
  gsize len;

  if (spp->r_edges->len < 2)
    g_error ("ncm_spectral_polar_panels_append: the edges are not set.");

  if (!(ra < rb) || !(ta < tb))
    g_error ("ncm_spectral_polar_panels_append: the panel [%.17g, %.17g] x [%.17g, %.17g] is empty.", ra, rb, ta, tb);

  if ((l_r > NCM_SPECTRAL_POLAR_PANELS_LEVEL_MAX) || (l_theta > NCM_SPECTRAL_POLAR_PANELS_LEVEL_MAX))
    g_error ("ncm_spectral_polar_panels_append: the levels %u and %u are above %u.", l_r, l_theta,
             NCM_SPECTRAL_POLAR_PANELS_LEVEL_MAX);

  i    = _ncm_spectral_polar_panels_edge_index (spp->r_edges, 0.5 * (ra + rb));
  j    = _ncm_spectral_polar_panels_edge_index (spp->theta_edges, 0.5 * (ta + tb));
  cell = i * (spp->theta_edges->len - 1) + j;

  if ((ra < g_array_index (spp->r_edges, gdouble, i)) || (rb > g_array_index (spp->r_edges, gdouble, i + 1)) ||
      (ta < g_array_index (spp->theta_edges, gdouble, j)) || (tb > g_array_index (spp->theta_edges, gdouble, j + 1)))
    g_error ("ncm_spectral_polar_panels_append: the panel [%.17g, %.17g] x [%.17g, %.17g] is not inside one initial panel.",
             ra, rb, ta, tb);

  if (first[cell + 1] != spp->panels->len)
    g_error ("ncm_spectral_polar_panels_append: the panel [%.17g, %.17g] x [%.17g, %.17g] is in the initial panel %u, "
             "and a later initial panel has panels already.", ra, rb, ta, tb, cell);

  len = (gsize) ((1u << l_r) + 1) * ((1u << l_theta) + 1);
  g_array_append_vals (spp->coeffs, c, len);
  g_array_append_val (spp->panels, panel);

  for (q = cell + 1; q <= n_cells; q++)
    first[q] = spp->panels->len;
}

/*
 * The index i of the interval [edges_i, edges_{i+1}) that holds x, by bisection. A point
 * below the first edge gives the first interval, and a point on or above the last edge the
 * last one.
 */
static guint
_ncm_spectral_polar_panels_edge_index (const GArray *edges, gdouble x)
{
  const gdouble *e = (const gdouble *) edges->data;
  guint lo         = 0;
  guint hi         = edges->len - 2;

  while (lo < hi)
  {
    const guint mid = (lo + hi + 1) / 2;

    if (x >= e[mid])
      lo = mid;
    else
      hi = mid - 1;
  }

  return lo;
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
 * @l_r: (out): radial level of the panel
 * @l_theta: (out): angular level of the panel
 *
 * Gets panel @i, which must be below ncm_spectral_polar_panels_get_n_panels(). Its
 * $2^{l_\theta} + 1$ rows of $2^{l_r} + 1$ coefficients are contiguous, and the returned
 * pointer is valid until @spp is changed.
 *
 * Returns: (transfer none): the coefficients of panel @i.
 */
const gdouble *
ncm_spectral_polar_panels_peek_panel (NcmSpectralPolarPanels *spp, guint i, gdouble *ra, gdouble *rb, gdouble *ta, gdouble *tb, guint *l_r, guint *l_theta)
{
  const NcmSpectralPolarPanelsPanel *panel;

  g_assert_cmpuint (i, <, spp->panels->len);

  panel = &g_array_index (spp->panels, NcmSpectralPolarPanelsPanel, i);

  *ra      = panel->ra;
  *rb      = panel->rb;
  *ta      = panel->ta;
  *tb      = panel->tb;
  *l_r     = panel->l_r;
  *l_theta = panel->l_theta;

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
 * Returns: (transfer full): the $(2^{l_\theta} + 1)$ by $(2^{l_r} + 1)$ coefficients of panel @i.
 */
NcmMatrix *
ncm_spectral_polar_panels_get_panel (NcmSpectralPolarPanels *spp, guint i, gdouble *ra, gdouble *rb, gdouble *ta, gdouble *tb)
{
  guint l_r, l_theta;
  const gdouble *c  = ncm_spectral_polar_panels_peek_panel (spp, i, ra, rb, ta, tb, &l_r, &l_theta);
  const guint n_r   = (1u << l_r) + 1;
  const guint n_t   = (1u << l_theta) + 1;
  NcmMatrix *coeffs = ncm_matrix_new (n_t, n_r);

  memcpy (ncm_matrix_data (coeffs), c, sizeof (gdouble) * n_t * n_r);

  return coeffs;
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
  const guint i      = _ncm_spectral_polar_panels_edge_index (spp->r_edges, r);
  const guint j      = _ncm_spectral_polar_panels_edge_index (spp->theta_edges, theta);
  const guint cell   = i * (spp->theta_edges->len - 1) + j;
  const guint *first = (const guint *) spp->first->data;
  const guint p0     = first[cell];
  const guint p1     = first[cell + 1];

  g_assert_cmpuint (p0, <, p1);

  if (p0 + 1 == p1)
    return p0;

  {
    const NcmSpectralPolarPanelsPanel *panels = (const NcmSpectralPolarPanelsPanel *) spp->panels->data;
    const gdouble cell_ra                     = g_array_index (spp->r_edges, gdouble, i);
    const gdouble cell_rb                     = g_array_index (spp->r_edges, gdouble, i + 1);
    const gdouble cell_ta                     = g_array_index (spp->theta_edges, gdouble, j);
    const gdouble cell_tb                     = g_array_index (spp->theta_edges, gdouble, j + 1);
    guint p;

    /* A panel on a side of the initial panel also takes the points beyond that side */
    for (p = p0; p + 1 < p1; p++)
    {
      const NcmSpectralPolarPanelsPanel *panel = &panels[p];

      if (((r >= panel->ra) || (panel->ra == cell_ra)) && ((r < panel->rb) || (panel->rb == cell_rb)) &&
          ((theta >= panel->ta) || (panel->ta == cell_ta)) && ((theta < panel->tb) || (panel->tb == cell_tb)))
        return p;
    }

    return p;
  }
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
  const gdouble *c                         = &g_array_index (spp->coeffs, gdouble, panel->offset);
  const guint n                            = 1u << panel->l_r;
  const guint m                            = 1u << panel->l_theta;
  const gdouble sr                         = (2.0 * r - panel->ra - panel->rb) / (panel->rb - panel->ra);
  const gdouble st                         = (2.0 * theta - panel->ta - panel->tb) / (panel->tb - panel->ta);
  const gdouble sr2                        = 2.0 * sr;
  const gdouble st2                        = 2.0 * st;
  gdouble B1                               = 0.0;
  gdouble B2                               = 0.0;
  guint k, l;

  for (l = m; l > 0; l--)
  {
    const gdouble *row = c + (gsize) l * (n + 1);
    gdouble b1         = 0.0;
    gdouble b2         = 0.0;

    for (k = n; k > 0; k--)
    {
      const gdouble b0 = sr2 * b1 - b2 + row[k];

      b2 = b1;
      b1 = b0;
    }

    {
      const gdouble B0 = st2 * B1 - B2 + (sr * b1 - b2 + row[0]);

      B2 = B1;
      B1 = B0;
    }
  }

  {
    gdouble b1 = 0.0;
    gdouble b2 = 0.0;

    for (k = n; k > 0; k--)
    {
      const gdouble b0 = sr2 * b1 - b2 + c[k];

      b2 = b1;
      b1 = b0;
    }

    return st * B1 - B2 + (sr * b1 - b2 + c[0]);
  }
}

