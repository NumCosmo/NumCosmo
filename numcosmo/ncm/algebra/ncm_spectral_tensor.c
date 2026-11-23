/***************************************************************************
 *            ncm_spectral_tensor.c
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

/**
 * NcmSpectralTensor:
 *
 * Tensor-product Chebyshev series of a function of several variables on a divided domain.
 *
 * A smooth scalar function of several variables that is expensive to evaluate can be
 * tabulated once and recomputed from the table. An #NcmSpectralTensor is such a table: it
 * divides the domain of the function into panels and stores, for each panel, the
 * coefficients of a tensor product of Chebyshev series. The value of the function at a
 * point is approximated by the sum of the series of the panel that contains the point. The
 * class #NcmSpectralTensorBuilder computes the table from the function, with an estimated
 * error below a requested tolerance; the table can also be filled by hand.
 *
 * The table approximates a real function $F (x_0, \dots, x_{n - 1})$ of $n$ variables, the
 * property #NcmSpectralTensor:n-vars, on the domain
 * $$\Omega \equiv [a_0, b_0] \times [a_1, b_1] \times \dots \times [a_{n - 1}, b_{n - 1}].$$
 * The table divides $\Omega$ recursively, and the subdomains form a tree. Each node of the
 * tree has a subdomain, the Cartesian product of one interval per variable. The interval
 * of a node in variable $i$ is $[p_i, q_i]$. The subdomain of the root is $\Omega$, so for
 * the root $p_i = a_i$ and $q_i = b_i$. A node has one of the three kinds of
 * #NcmSpectralTensorNodeKind, and only a split node has children. A split node divides its
 * subdomain along one variable $i$: with the interior edges
 * $p_i < e_1 < \dots < e_{m - 1} < q_i$, its $m$ children have the intervals
 * $[p_i, e_1], [e_1, e_2], \dots, [e_{m - 1}, q_i]$ in variable $i$ and the intervals of
 * the node in the other variables. A panel node is not divided and holds a series. An
 * empty node is not divided and holds no series. The panels of a table with no empty node
 * cover $\Omega$ and intersect only on their boundaries.
 *
 * A panel is the subdomain of a panel node,
 * $$\mathcal{P} \equiv [p_0, q_0] \times [p_1, q_1] \times \dots \times [p_{n - 1}, q_{n - 1}],$$
 * together with its series. For $(x_0, \dots, x_{n - 1}) \in \mathcal{P}$ the table
 * approximates
 * $$F (x_0, \dots, x_{n - 1}) \approx \sum_{k_0 = 0}^{K_0 - 1} \dots \sum_{k_{n - 1} = 0}^{K_{n - 1} - 1} c^{k_0 \dots k_{n - 1}} \prod_{i = 0}^{n - 1} T_{k_i} (s_i),$$
 * where $T_k$ is the Chebyshev polynomial of the first kind and degree $k$,
 * $s_i = (2 x_i - p_i - q_i) / (q_i - p_i)$, $K_i \ge 1$ is the number of coefficients of
 * the panel in variable $i$, and $c^{k_0 \dots k_{n - 1}}$ are its coefficients. The
 * numbers $K_i$ may differ between variables and between panels. The coefficients of a
 * panel are stored with the first variable running fastest: the coefficient
 * $c^{k_0 \dots k_{n - 1}}$ is at position $k_0 + K_0 (k_1 + K_1 (k_2 + \dots))$.
 *
 * Nodes are numbered from 0, the root, in the order they are created, and the children of
 * a split node have consecutive numbers. A panel node also has a panel number, from 0, in
 * the order the series are set. To fill the table by hand,
 * ncm_spectral_tensor_set_domain() sets $\Omega$ and makes the root an empty node,
 * ncm_spectral_tensor_split() makes an empty node a split node, and
 * ncm_spectral_tensor_set_panel() makes an empty node a panel node. One call of
 * ncm_spectral_tensor_split_empty() per variable divides $\Omega$ into a grid.
 *
 * The function ncm_spectral_tensor_eval() finds the panel that contains a point and sums
 * its series. The search starts at the root and, at each split node, continues at the
 * child whose interval contains the coordinate of the point in the variable of the split
 * node. A point on an interior edge belongs to the child whose interval starts at that
 * edge. The point is not checked to be in $\Omega$: a point outside $\Omega$ reaches the
 * panel that contains the point obtained by clamping each coordinate $x_i$ to
 * $[a_i, b_i]$, and the series of that panel is extrapolated. The tolerance of the table
 * does not apply to an extrapolated value, and the evaluation does not report it. A caller
 * whose points may be outside $\Omega$ compares each coordinate with the range given by
 * ncm_spectral_tensor_get_range() before the call, and computes the function directly for
 * a point outside. A search that reaches an empty node aborts with an error. The function
 * ncm_spectral_tensor_eval_array() evaluates several points in one call and sums
 * consecutive points on the same panel together, which takes less time per point.
 *
 * The functions ncm_spectral_tensor_set_domain() and ncm_spectral_tensor_reset() keep the
 * arrays, so a table refilled many times stops allocating once it has grown. Reading and
 * evaluating do not change the table, so one table may be evaluated in several threads.
 *
 * See <a href="../../theory/ncm/algebra/spectral_tensor.html">Tensor-Product Chebyshev
 * Tables</a> for the formulas.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/algebra/ncm_spectral_tensor.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <string.h>
#endif /* NUMCOSMO_GIR_SCAN */

/* Largest number of coefficients of a panel in one variable, 2^20 + 1 */
#define NCM_SPECTRAL_TENSOR_N_COEFFS_MAX ((1u << 20) + 1)

enum
{
  PROP_0,
  PROP_N_VARS,
  PROP_TREE,
};

/*
 * A node. A split node has its variable, its n_edges interior edges at edges_offset in
 * edges and the index of its first child in index; a panel node has its panel index in
 * index.
 */
typedef struct _NcmSpectralTensorNode
{
  guint kind;
  guint index;
  guint var;
  guint n_edges;
  gsize edges_offset;
} NcmSpectralTensorNode;

/* A panel: its panel node and the position of its coefficients in coeffs */
typedef struct _NcmSpectralTensorPanel
{
  guint node;
  gsize offset;
} NcmSpectralTensorPanel;

struct _NcmSpectralTensor
{
  /*< private >*/
  GObject parent_instance;
  guint n_vars;
  guint n_empty;       /* number of empty nodes */
  GArray *nodes;       /* NcmSpectralTensorNode; no node until the domain is set */
  GArray *node_bounds; /* gdouble, 2 n_vars per node: its lower ends, then its upper ends */
  GArray *splits;      /* guint, the split nodes, in creation order */
  GArray *edges;       /* gdouble, the interior edges of all split nodes in sequence */
  GArray *panels;      /* NcmSpectralTensorPanel, in creation order */
  GArray *n_coeffs;    /* guint, n_vars per panel */
  GArray *coeffs;      /* gdouble, the coefficients of every panel one after the other */
  GVariant *init_tree; /* the construct-only property, consumed by constructed() */
};

G_DEFINE_TYPE (NcmSpectralTensor, ncm_spectral_tensor, G_TYPE_OBJECT)

static void
ncm_spectral_tensor_init (NcmSpectralTensor *st)
{
  st->n_vars      = 0;
  st->n_empty     = 0;
  st->nodes       = g_array_new (FALSE, FALSE, sizeof (NcmSpectralTensorNode));
  st->node_bounds = g_array_new (FALSE, FALSE, sizeof (gdouble));
  st->splits      = g_array_new (FALSE, FALSE, sizeof (guint));
  st->edges       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  st->panels      = g_array_new (FALSE, FALSE, sizeof (NcmSpectralTensorPanel));
  st->n_coeffs    = g_array_new (FALSE, FALSE, sizeof (guint));
  st->coeffs      = g_array_new (FALSE, FALSE, sizeof (gdouble));
  st->init_tree   = NULL;
}

static void
_ncm_spectral_tensor_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSpectralTensor *st = NCM_SPECTRAL_TENSOR (object);

  g_return_if_fail (NCM_IS_SPECTRAL_TENSOR (object));

  switch (prop_id)
  {
    case PROP_N_VARS:
      st->n_vars = g_value_get_uint (value);
      break;
    case PROP_TREE:
      st->init_tree = g_value_dup_variant (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static GVariant *_ncm_spectral_tensor_get_variant (NcmSpectralTensor *st);

static void
_ncm_spectral_tensor_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSpectralTensor *st = NCM_SPECTRAL_TENSOR (object);

  g_return_if_fail (NCM_IS_SPECTRAL_TENSOR (object));

  switch (prop_id)
  {
    case PROP_N_VARS:
      g_value_set_uint (value, st->n_vars);
      break;
    case PROP_TREE:
      g_value_take_variant (value, _ncm_spectral_tensor_get_variant (st));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static gsize _ncm_spectral_tensor_block_len (const guint *n_coeffs, guint n_vars);

/*
 * The tree as a (adada(uuad)a(uauad)): the lower and the upper ends of the domain, the
 * split nodes in creation order and the panel nodes in creation order. The arrays are
 * copied: the functions that fill the table rewrite them in place, and a variant must not
 * change after it is made. The four arrays are empty when the domain is not set.
 */
static GVariant *
_ncm_spectral_tensor_get_variant (NcmSpectralTensor *st)
{
  const guint n_vars = st->n_vars;
  const guint n_dom  = (st->nodes->len > 0) ? n_vars : 0;
  GVariantBuilder splits_builder;
  GVariantBuilder panels_builder;
  guint q;

  g_variant_builder_init (&splits_builder, G_VARIANT_TYPE ("a(uuad)"));

  for (q = 0; q < st->splits->len; q++)
  {
    const guint node                = g_array_index (st->splits, guint, q);
    const NcmSpectralTensorNode *nd = &g_array_index (st->nodes, NcmSpectralTensorNode, node);

    g_variant_builder_add (&splits_builder, "(uu@ad)", node, nd->var,
                           g_variant_new_fixed_array (G_VARIANT_TYPE ("d"), &g_array_index (st->edges, gdouble, nd->edges_offset),
                                                      nd->n_edges, sizeof (gdouble)));
  }

  g_variant_builder_init (&panels_builder, G_VARIANT_TYPE ("a(uauad)"));

  for (q = 0; q < st->panels->len; q++)
  {
    const NcmSpectralTensorPanel *pn = &g_array_index (st->panels, NcmSpectralTensorPanel, q);
    const guint *K                   = &g_array_index (st->n_coeffs, guint, (gsize) n_vars * q);

    g_variant_builder_add (&panels_builder, "(u@au@ad)", pn->node,
                           g_variant_new_fixed_array (G_VARIANT_TYPE ("u"), K, n_vars, sizeof (guint)),
                           g_variant_new_fixed_array (G_VARIANT_TYPE ("d"), &g_array_index (st->coeffs, gdouble, pn->offset),
                                                      _ncm_spectral_tensor_block_len (K, n_vars), sizeof (gdouble)));
  }

  return g_variant_ref_sink (g_variant_new ("(@ad@ad@a(uuad)@a(uauad))",
                                            g_variant_new_fixed_array (G_VARIANT_TYPE ("d"), st->node_bounds->data, n_dom, sizeof (gdouble)),
                                            g_variant_new_fixed_array (G_VARIANT_TYPE ("d"), st->node_bounds->data + sizeof (gdouble) * n_vars, n_dom, sizeof (gdouble)),
                                            g_variant_builder_end (&splits_builder), g_variant_builder_end (&panels_builder)));
}

/*
 * The tree given as a property: the domain, then the split nodes and the panel nodes in
 * their creation order, each through the public function, so that node and panel indices
 * are those of the original and every entry is checked as one made by hand.
 */
static void
_ncm_spectral_tensor_constructed (GObject *object)
{
  NcmSpectralTensor *st = NCM_SPECTRAL_TENSOR (object);

  /* Chain up : start */
  G_OBJECT_CLASS (ncm_spectral_tensor_parent_class)->constructed (object);

  if (st->init_tree != NULL)
  {
    GVariant *lower_var  = g_variant_get_child_value (st->init_tree, 0);
    GVariant *upper_var  = g_variant_get_child_value (st->init_tree, 1);
    GVariant *splits_var = g_variant_get_child_value (st->init_tree, 2);
    GVariant *panels_var = g_variant_get_child_value (st->init_tree, 3);
    const guint n_splits = g_variant_n_children (splits_var);
    const guint n_panels = g_variant_n_children (panels_var);
    gsize n_lower        = 0;
    gsize n_upper        = 0;
    const gdouble *lower = g_variant_get_fixed_array (lower_var, &n_lower, sizeof (gdouble));
    const gdouble *upper = g_variant_get_fixed_array (upper_var, &n_upper, sizeof (gdouble));
    guint q;

    if ((n_lower > 0) || (n_upper > 0) || (n_splits > 0) || (n_panels > 0))
      ncm_spectral_tensor_set_domain (st, n_lower, lower, n_upper, upper);

    for (q = 0; q < n_splits; q++)
    {
      GVariant *split_var = g_variant_get_child_value (splits_var, q);
      GVariant *edges_var;
      guint node, var;
      gsize n_edges = 0;
      const gdouble *edges;

      g_variant_get (split_var, "(uu@ad)", &node, &var, &edges_var);
      edges = g_variant_get_fixed_array (edges_var, &n_edges, sizeof (gdouble));

      ncm_spectral_tensor_split (st, node, var, n_edges, edges);

      g_variant_unref (edges_var);
      g_variant_unref (split_var);
    }

    for (q = 0; q < n_panels; q++)
    {
      GVariant *panel_var = g_variant_get_child_value (panels_var, q);
      GVariant *K_var, *c_var;
      guint node;
      gsize n_K = 0;
      gsize len = 0;
      const guint *K;
      const gdouble *c;

      g_variant_get (panel_var, "(u@au@ad)", &node, &K_var, &c_var);
      K = g_variant_get_fixed_array (K_var, &n_K, sizeof (guint));
      c = g_variant_get_fixed_array (c_var, &len, sizeof (gdouble));

      ncm_spectral_tensor_set_panel (st, node, n_K, K, len, c);

      g_variant_unref (K_var);
      g_variant_unref (c_var);
      g_variant_unref (panel_var);
    }

    g_variant_unref (lower_var);
    g_variant_unref (upper_var);
    g_variant_unref (splits_var);
    g_variant_unref (panels_var);
    g_clear_pointer (&st->init_tree, g_variant_unref);
  }
}

static void
_ncm_spectral_tensor_finalize (GObject *object)
{
  NcmSpectralTensor *st = NCM_SPECTRAL_TENSOR (object);

  g_array_unref (st->nodes);
  g_array_unref (st->node_bounds);
  g_array_unref (st->splits);
  g_array_unref (st->edges);
  g_array_unref (st->panels);
  g_array_unref (st->n_coeffs);
  g_array_unref (st->coeffs);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spectral_tensor_parent_class)->finalize (object);
}

static void
ncm_spectral_tensor_class_init (NcmSpectralTensorClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_spectral_tensor_set_property;
  object_class->get_property = &_ncm_spectral_tensor_get_property;
  object_class->constructed  = &_ncm_spectral_tensor_constructed;
  object_class->finalize     = &_ncm_spectral_tensor_finalize;

  /**
   * NcmSpectralTensor:n-vars:
   *
   * The number of variables $n$ of $F$, from 1 to %NCM_SPECTRAL_TENSOR_MAX_VARS.
   */
  g_object_class_install_property (object_class,
                                   PROP_N_VARS,
                                   g_param_spec_uint ("n-vars",
                                                      NULL,
                                                      "Number of variables",
                                                      1, NCM_SPECTRAL_TENSOR_MAX_VARS, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSpectralTensor:tree:
   *
   * The table as a #GVariant of type `(adada(uuad)a(uauad))`: the lower ends $a_i$ and the
   * upper ends $b_i$ of the domain, then for each split node, in creation order, its node,
   * its variable and its interior edges, then for each panel node, in creation order, its
   * node, its numbers of coefficients $K_i$ and its coefficients, the first variable
   * running fastest.
   * The four arrays are empty when the domain is not set.
   */
  g_object_class_install_property (object_class,
                                   PROP_TREE,
                                   g_param_spec_variant ("tree",
                                                         NULL,
                                                         "Domain, split nodes and panel nodes",
                                                         G_VARIANT_TYPE ("(adada(uuad)a(uauad))"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_spectral_tensor_new:
 * @n_vars: number of variables of $F$
 *
 * Creates a table of a function of @n_vars variables, with no domain.
 *
 * Returns: (transfer full): a new #NcmSpectralTensor.
 */
NcmSpectralTensor *
ncm_spectral_tensor_new (guint n_vars)
{
  return g_object_new (NCM_TYPE_SPECTRAL_TENSOR,
                       "n-vars", n_vars,
                       NULL);
}

/**
 * ncm_spectral_tensor_ref:
 * @st: a #NcmSpectralTensor
 *
 * Increases the reference count of @st by one.
 *
 * Returns: (transfer full): @st.
 */
NcmSpectralTensor *
ncm_spectral_tensor_ref (NcmSpectralTensor *st)
{
  return g_object_ref (st);
}

/**
 * ncm_spectral_tensor_free:
 * @st: a #NcmSpectralTensor
 *
 * Decreases the reference count of @st by one.
 */
void
ncm_spectral_tensor_free (NcmSpectralTensor *st)
{
  g_object_unref (st);
}

/**
 * ncm_spectral_tensor_clear:
 * @st: a #NcmSpectralTensor
 *
 * Decreases the reference count of *@st by one and sets *@st to %NULL.
 */
void
ncm_spectral_tensor_clear (NcmSpectralTensor **st)
{
  g_clear_object (st);
}

/* Leaves the tree with the root only, empty; its bounds, the first 2 n_vars values of node_bounds, are kept */
static void
_ncm_spectral_tensor_root_only (NcmSpectralTensor *st)
{
  const NcmSpectralTensorNode root = {NCM_SPECTRAL_TENSOR_NODE_KIND_EMPTY, 0, 0, 0, 0};

  g_array_set_size (st->nodes, 0);
  g_array_append_val (st->nodes, root);
  g_array_set_size (st->node_bounds, 2 * st->n_vars);
  g_array_set_size (st->splits, 0);
  g_array_set_size (st->edges, 0);
  g_array_set_size (st->panels, 0);
  g_array_set_size (st->n_coeffs, 0);
  g_array_set_size (st->coeffs, 0);

  st->n_empty = 1;
}

/**
 * ncm_spectral_tensor_set_domain:
 * @st: a #NcmSpectralTensor
 * @n_lower: number of entries of @lower
 * @lower: (array length=n_lower): the lower ends $a_i$ of the domain
 * @n_upper: number of entries of @upper
 * @upper: (array length=n_upper): the upper ends $b_i$ of the domain
 *
 * To start filling a table, sets its domain $\Omega$ and removes every split node and
 * every panel node: the tree is then the root, node 0, an empty node. The call is an error
 * unless @n_lower and @n_upper are #NcmSpectralTensor:n-vars and every $a_i < b_i$.
 */
void
ncm_spectral_tensor_set_domain (NcmSpectralTensor *st, guint n_lower, const gdouble *lower, guint n_upper, const gdouble *upper)
{
  guint i;

  if ((n_lower != st->n_vars) || (n_upper != st->n_vars))
    g_error ("ncm_spectral_tensor_set_domain: %u lower ends and %u upper ends for %u variables.", n_lower, n_upper, st->n_vars);

  for (i = 0; i < st->n_vars; i++)
  {
    if (!(lower[i] < upper[i]))
      g_error ("ncm_spectral_tensor_set_domain: the range [%.17g, %.17g] of variable %u is empty.", lower[i], upper[i], i);
  }

  g_array_set_size (st->node_bounds, 0);
  g_array_append_vals (st->node_bounds, lower, st->n_vars);
  g_array_append_vals (st->node_bounds, upper, st->n_vars);

  _ncm_spectral_tensor_root_only (st);
}

/**
 * ncm_spectral_tensor_reset:
 * @st: a #NcmSpectralTensor
 *
 * Before a refill on the same domain, removes every split node and every panel node: the
 * tree is then the root, node 0, an empty node. A table with no domain is left as it is.
 * The arrays are kept for reuse, so refilling @st allocates nothing until it holds more
 * nodes or more coefficients than before.
 */
void
ncm_spectral_tensor_reset (NcmSpectralTensor *st)
{
  if (st->nodes->len > 0)
    _ncm_spectral_tensor_root_only (st);
}

/* Aborts unless node is an empty node of st; func is the function called, for the message */
static void
_ncm_spectral_tensor_check_empty (NcmSpectralTensor *st, const gchar *func, guint node)
{
  if (node >= st->nodes->len)
    g_error ("%s: node %u, the table has %u nodes.", func, node, st->nodes->len);

  if (g_array_index (st->nodes, NcmSpectralTensorNode, node).kind != NCM_SPECTRAL_TENSOR_NODE_KIND_EMPTY)
    g_error ("%s: node %u is not empty.", func, node);
}

/**
 * ncm_spectral_tensor_split:
 * @st: a #NcmSpectralTensor
 * @node: node index
 * @i: variable index
 * @n_edges: number of interior edges
 * @edges: (array length=n_edges): the interior edges $e_1, \dots, e_{m - 1}$
 *
 * To divide a subdomain, makes the empty node @node a split node of variable @i with
 * @n_edges + 1 children, copying @edges. The children are new empty nodes with consecutive
 * indices, in the order of their intervals. The call is an error unless @node is an empty
 * node, @i is below #NcmSpectralTensor:n-vars, @n_edges is 1 or more, and the edges are
 * ascending and strictly inside the interval of @node in variable @i.
 *
 * Returns: the index of the first child.
 */
guint
ncm_spectral_tensor_split (NcmSpectralTensor *st, guint node, guint i, guint n_edges, const gdouble *edges)
{
  const guint n_vars                = st->n_vars;
  const NcmSpectralTensorNode empty = {NCM_SPECTRAL_TENSOR_NODE_KIND_EMPTY, 0, 0, 0, 0};
  gdouble lower[NCM_SPECTRAL_TENSOR_MAX_VARS];
  gdouble upper[NCM_SPECTRAL_TENSOR_MAX_VARS];
  gsize edges_offset;
  guint first_child;
  gdouble a, b;
  guint k;

  _ncm_spectral_tensor_check_empty (st, "ncm_spectral_tensor_split", node);

  if (i >= n_vars)
    g_error ("ncm_spectral_tensor_split: variable %u, the table has %u variables.", i, n_vars);

  if (n_edges < 1)
    g_error ("ncm_spectral_tensor_split: a split needs one interior edge or more.");

  if ((guint64) st->nodes->len + n_edges + 1 > G_MAXUINT)
    g_error ("ncm_spectral_tensor_split: more than %u nodes.", G_MAXUINT);

  /* The bounds of the node are copied: appending the children may move node_bounds */
  memcpy (lower, &g_array_index (st->node_bounds, gdouble, (gsize) 2 * n_vars * node), sizeof (gdouble) * n_vars);
  memcpy (upper, &g_array_index (st->node_bounds, gdouble, (gsize) 2 * n_vars * node + n_vars), sizeof (gdouble) * n_vars);

  a = lower[i];
  b = upper[i];

  for (k = 0; k < n_edges; k++)
  {
    if (!(((k == 0) ? a : edges[k - 1]) < edges[k]) || !(edges[k] < b))
      g_error ("ncm_spectral_tensor_split: the interior edges of variable %u must be ascending and strictly inside "
               "(%.17g, %.17g), edge %u is %.17g.", i, a, b, k, edges[k]);
  }

  edges_offset = st->edges->len;
  first_child  = st->nodes->len;

  g_array_append_vals (st->edges, edges, n_edges);
  g_array_append_val (st->splits, node);

  {
    NcmSpectralTensorNode *nd = &g_array_index (st->nodes, NcmSpectralTensorNode, node);

    nd->kind         = NCM_SPECTRAL_TENSOR_NODE_KIND_SPLIT;
    nd->index        = first_child;
    nd->var          = i;
    nd->n_edges      = n_edges;
    nd->edges_offset = edges_offset;
  }

  for (k = 0; k <= n_edges; k++)
  {
    lower[i] = (k == 0) ? a : g_array_index (st->edges, gdouble, edges_offset + k - 1);
    upper[i] = (k == n_edges) ? b : g_array_index (st->edges, gdouble, edges_offset + k);

    g_array_append_val (st->nodes, empty);
    g_array_append_vals (st->node_bounds, lower, n_vars);
    g_array_append_vals (st->node_bounds, upper, n_vars);
  }

  st->n_empty += n_edges;

  return first_child;
}

/**
 * ncm_spectral_tensor_split_empty:
 * @st: a #NcmSpectralTensor
 * @i: variable index
 * @n_edges: number of edges
 * @edges: (array length=n_edges): edges of variable @i
 *
 * To divide the domain as a grid, one call per variable, splits every node that is empty
 * when the function is called along variable @i, as ncm_spectral_tensor_split() with the
 * entries of @edges that are strictly inside the interval of the node in that variable. A
 * node with no such entry is not split, so @edges may include the ends of the domain. The
 * call is an error unless the domain is set, @i is below #NcmSpectralTensor:n-vars and the
 * edges are ascending.
 *
 * Returns: the number of nodes split.
 */
guint
ncm_spectral_tensor_split_empty (NcmSpectralTensor *st, guint i, guint n_edges, const gdouble *edges)
{
  const guint n_vars  = st->n_vars;
  const guint n_nodes = st->nodes->len;
  guint n_split       = 0;
  guint node, k;

  if (n_nodes == 0)
    g_error ("ncm_spectral_tensor_split_empty: the domain is not set.");

  if (i >= n_vars)
    g_error ("ncm_spectral_tensor_split_empty: variable %u, the table has %u variables.", i, n_vars);

  for (k = 0; k + 1 < n_edges; k++)
  {
    if (!(edges[k] < edges[k + 1]))
      g_error ("ncm_spectral_tensor_split_empty: the edges of variable %u are not ascending at index %u.", i, k);
  }

  for (node = 0; node < n_nodes; node++)
  {
    if (g_array_index (st->nodes, NcmSpectralTensorNode, node).kind == NCM_SPECTRAL_TENSOR_NODE_KIND_EMPTY)
    {
      const gdouble a = g_array_index (st->node_bounds, gdouble, (gsize) 2 * n_vars * node + i);
      const gdouble b = g_array_index (st->node_bounds, gdouble, (gsize) 2 * n_vars * node + n_vars + i);
      guint k0        = 0;
      guint k1;

      while ((k0 < n_edges) && !(edges[k0] > a))
        k0++;

      k1 = k0;

      while ((k1 < n_edges) && (edges[k1] < b))
        k1++;

      if (k1 > k0)
      {
        ncm_spectral_tensor_split (st, node, i, k1 - k0, edges + k0);
        n_split++;
      }
    }
  }

  return n_split;
}

/**
 * ncm_spectral_tensor_set_panel:
 * @st: a #NcmSpectralTensor
 * @node: node index
 * @n_K: number of entries of @n_coeffs
 * @n_coeffs: (array length=n_K): the numbers of coefficients $K_i$ of the panel
 * @len: number of coefficients of the panel, $\prod_i K_i$
 * @coeffs: (array length=len): the coefficients of the panel, the first variable running fastest
 *
 * To give a subdomain its series, makes the empty node @node a panel node, copying
 * @coeffs. The call is an error unless @node is an empty node, @n_K is
 * #NcmSpectralTensor:n-vars, every $K_i \ge 1$ and @len is $\prod_i K_i$.
 *
 * Returns: the index of the panel.
 */
guint
ncm_spectral_tensor_set_panel (NcmSpectralTensor *st, guint node, guint n_K, const guint *n_coeffs, gsize len, const gdouble *coeffs)
{
  const NcmSpectralTensorPanel pn = {node, st->coeffs->len};
  guint i;

  _ncm_spectral_tensor_check_empty (st, "ncm_spectral_tensor_set_panel", node);

  if (n_K != st->n_vars)
    g_error ("ncm_spectral_tensor_set_panel: %u numbers of coefficients for %u variables.", n_K, st->n_vars);

  for (i = 0; i < n_K; i++)
  {
    if ((n_coeffs[i] < 1) || (n_coeffs[i] > NCM_SPECTRAL_TENSOR_N_COEFFS_MAX))
      g_error ("ncm_spectral_tensor_set_panel: %u coefficients in variable %u, outside [1, %u].",
               n_coeffs[i], i, NCM_SPECTRAL_TENSOR_N_COEFFS_MAX);
  }

  if (len != _ncm_spectral_tensor_block_len (n_coeffs, n_K))
    g_error ("ncm_spectral_tensor_set_panel: %zu coefficients, the numbers of coefficients give %zu.",
             len, _ncm_spectral_tensor_block_len (n_coeffs, n_K));

  g_array_append_vals (st->coeffs, coeffs, len);
  g_array_append_vals (st->n_coeffs, n_coeffs, n_K);
  g_array_append_val (st->panels, pn);

  g_array_index (st->nodes, NcmSpectralTensorNode, node).kind  = NCM_SPECTRAL_TENSOR_NODE_KIND_PANEL;
  g_array_index (st->nodes, NcmSpectralTensorNode, node).index = st->panels->len - 1;

  st->n_empty--;

  return st->panels->len - 1;
}

/* The number of coefficients of a panel, prod_i n_coeffs[i] */
static gsize
_ncm_spectral_tensor_block_len (const guint *n_coeffs, guint n_vars)
{
  gsize len = 1;
  guint i;

  for (i = 0; i < n_vars; i++)
    len *= n_coeffs[i];

  return len;
}

/**
 * ncm_spectral_tensor_get_n_vars:
 * @st: a #NcmSpectralTensor
 *
 * Returns: the number of variables, #NcmSpectralTensor:n-vars.
 */
guint
ncm_spectral_tensor_get_n_vars (NcmSpectralTensor *st)
{
  return st->n_vars;
}

/**
 * ncm_spectral_tensor_get_n_nodes:
 * @st: a #NcmSpectralTensor
 *
 * Returns: the number of nodes, 0 when the domain is not set.
 */
guint
ncm_spectral_tensor_get_n_nodes (NcmSpectralTensor *st)
{
  return st->nodes->len;
}

/**
 * ncm_spectral_tensor_get_n_empty:
 * @st: a #NcmSpectralTensor
 *
 * A table with a domain and no empty node can be evaluated at every point.
 *
 * Returns: the number of empty nodes.
 */
guint
ncm_spectral_tensor_get_n_empty (NcmSpectralTensor *st)
{
  return st->n_empty;
}

/**
 * ncm_spectral_tensor_get_n_panels:
 * @st: a #NcmSpectralTensor
 *
 * Returns: the number of panels.
 */
guint
ncm_spectral_tensor_get_n_panels (NcmSpectralTensor *st)
{
  return st->panels->len;
}

/**
 * ncm_spectral_tensor_get_range:
 * @st: a #NcmSpectralTensor
 * @i: variable index, below #NcmSpectralTensor:n-vars
 * @a: (out): the lower end $a_i$ of the domain
 * @b: (out): the upper end $b_i$ of the domain
 *
 * Gets the range of variable @i in the domain, which must be set.
 */
void
ncm_spectral_tensor_get_range (NcmSpectralTensor *st, guint i, gdouble *a, gdouble *b)
{
  g_assert_cmpuint (i, <, st->n_vars);
  g_assert_cmpuint (st->nodes->len, >, 0);

  *a = g_array_index (st->node_bounds, gdouble, i);
  *b = g_array_index (st->node_bounds, gdouble, st->n_vars + i);
}

/**
 * ncm_spectral_tensor_get_node_kind:
 * @st: a #NcmSpectralTensor
 * @node: node index, below ncm_spectral_tensor_get_n_nodes()
 *
 * Returns: the kind of node @node.
 */
NcmSpectralTensorNodeKind
ncm_spectral_tensor_get_node_kind (NcmSpectralTensor *st, guint node)
{
  g_assert_cmpuint (node, <, st->nodes->len);

  return g_array_index (st->nodes, NcmSpectralTensorNode, node).kind;
}

/**
 * ncm_spectral_tensor_peek_node_bounds: (skip)
 * @st: a #NcmSpectralTensor
 * @node: node index, below ncm_spectral_tensor_get_n_nodes()
 * @lower: (out): the lower ends of the subdomain of the node, one per variable
 * @upper: (out): the upper ends of the subdomain of the node, one per variable
 *
 * To read a subdomain without copies, gets the bounds of node @node. The pointers are valid
 * until @st is changed.
 */
void
ncm_spectral_tensor_peek_node_bounds (NcmSpectralTensor *st, guint node, const gdouble **lower, const gdouble **upper)
{
  g_assert_cmpuint (node, <, st->nodes->len);

  *lower = &g_array_index (st->node_bounds, gdouble, (gsize) 2 * st->n_vars * node);
  *upper = *lower + st->n_vars;
}

/**
 * ncm_spectral_tensor_get_node_bounds:
 * @st: a #NcmSpectralTensor
 * @node: node index, below ncm_spectral_tensor_get_n_nodes()
 * @lower: (out) (transfer full) (element-type gdouble): the lower ends of the subdomain of the node
 * @upper: (out) (transfer full) (element-type gdouble): the upper ends of the subdomain of the node
 *
 * As ncm_spectral_tensor_peek_node_bounds(), with copies in new arrays.
 */
void
ncm_spectral_tensor_get_node_bounds (NcmSpectralTensor *st, guint node, GArray **lower, GArray **upper)
{
  const gdouble *lo, *up;

  ncm_spectral_tensor_peek_node_bounds (st, node, &lo, &up);

  *lower = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), st->n_vars);
  *upper = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), st->n_vars);

  g_array_append_vals (*lower, lo, st->n_vars);
  g_array_append_vals (*upper, up, st->n_vars);
}

/**
 * ncm_spectral_tensor_get_node_split:
 * @st: a #NcmSpectralTensor
 * @node: node index of a split node
 * @i: (out): the variable of the split node
 * @first_child: (out): the index of its first child
 *
 * To walk the tree, gets the split node @node: its variable, its first child and a copy of
 * its interior edges; it has one child more than interior edges. A @node that is not a
 * split node is an error.
 *
 * Returns: (transfer full) (element-type gdouble): the interior edges of the split node.
 */
GArray *
ncm_spectral_tensor_get_node_split (NcmSpectralTensor *st, guint node, guint *i, guint *first_child)
{
  const NcmSpectralTensorNode *nd;
  GArray *edges;

  if ((node >= st->nodes->len) || (g_array_index (st->nodes, NcmSpectralTensorNode, node).kind != NCM_SPECTRAL_TENSOR_NODE_KIND_SPLIT))
    g_error ("ncm_spectral_tensor_get_node_split: node %u is not a split node.", node);

  nd    = &g_array_index (st->nodes, NcmSpectralTensorNode, node);
  edges = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), nd->n_edges);

  g_array_append_vals (edges, &g_array_index (st->edges, gdouble, nd->edges_offset), nd->n_edges);

  *i           = nd->var;
  *first_child = nd->index;

  return edges;
}

/**
 * ncm_spectral_tensor_get_node_panel:
 * @st: a #NcmSpectralTensor
 * @node: node index of a panel node
 *
 * To walk the tree, gets the panel number of the panel node @node. A @node that is not a
 * panel node is an error.
 *
 * Returns: the index of the panel.
 */
guint
ncm_spectral_tensor_get_node_panel (NcmSpectralTensor *st, guint node)
{
  if ((node >= st->nodes->len) || (g_array_index (st->nodes, NcmSpectralTensorNode, node).kind != NCM_SPECTRAL_TENSOR_NODE_KIND_PANEL))
    g_error ("ncm_spectral_tensor_get_node_panel: node %u is not a panel node.", node);

  return g_array_index (st->nodes, NcmSpectralTensorNode, node).index;
}

/**
 * ncm_spectral_tensor_peek_panel: (skip)
 * @st: a #NcmSpectralTensor
 * @p: panel index
 * @lower: (out): the lower ends $p_i$ of the panel, one per variable
 * @upper: (out): the upper ends $q_i$ of the panel, one per variable
 * @n_coeffs: (out): the numbers of coefficients $K_i$ of the panel, one per variable
 *
 * To read a panel without copies, gets panel @p, which must be below
 * ncm_spectral_tensor_get_n_panels(). Its $\prod_i K_i$ coefficients are contiguous, the
 * first variable running fastest. The returned pointers are valid until @st is changed.
 *
 * Returns: (transfer none): the coefficients of panel @p.
 */
const gdouble *
ncm_spectral_tensor_peek_panel (NcmSpectralTensor *st, guint p, const gdouble **lower, const gdouble **upper, const guint **n_coeffs)
{
  const guint n_vars = st->n_vars;
  const NcmSpectralTensorPanel *pn;

  g_assert_cmpuint (p, <, st->panels->len);

  pn = &g_array_index (st->panels, NcmSpectralTensorPanel, p);

  *lower    = &g_array_index (st->node_bounds, gdouble, (gsize) 2 * n_vars * pn->node);
  *upper    = *lower + n_vars;
  *n_coeffs = &g_array_index (st->n_coeffs, guint, (gsize) n_vars * p);

  return &g_array_index (st->coeffs, gdouble, pn->offset);
}

/**
 * ncm_spectral_tensor_get_panel:
 * @st: a #NcmSpectralTensor
 * @p: panel index
 * @lower: (out) (transfer full) (element-type gdouble): the lower ends $p_i$ of the panel
 * @upper: (out) (transfer full) (element-type gdouble): the upper ends $q_i$ of the panel
 * @n_coeffs: (out) (transfer full) (element-type guint): the numbers of coefficients $K_i$ of the panel
 *
 * As ncm_spectral_tensor_peek_panel(), with copies in new arrays. The coefficients have
 * the first variable running fastest.
 *
 * Returns: (transfer full) (element-type gdouble): the coefficients of panel @p.
 */
GArray *
ncm_spectral_tensor_get_panel (NcmSpectralTensor *st, guint p, GArray **lower, GArray **upper, GArray **n_coeffs)
{
  const guint n_vars = st->n_vars;
  const gdouble *lo, *up;
  const guint *K;
  const gdouble *c = ncm_spectral_tensor_peek_panel (st, p, &lo, &up, &K);
  const gsize len  = _ncm_spectral_tensor_block_len (K, n_vars);
  GArray *coeffs   = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), len);

  *lower    = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), n_vars);
  *upper    = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), n_vars);
  *n_coeffs = g_array_sized_new (FALSE, FALSE, sizeof (guint), n_vars);

  g_array_append_vals (*lower, lo, n_vars);
  g_array_append_vals (*upper, up, n_vars);
  g_array_append_vals (*n_coeffs, K, n_vars);
  g_array_append_vals (coeffs, c, len);

  return coeffs;
}

/**
 * ncm_spectral_tensor_find_panel:
 * @st: a #NcmSpectralTensor
 * @n_vars: number of variables, #NcmSpectralTensor:n-vars
 * @x: (array length=n_vars): the point
 *
 * In a filled table, finds the panel that contains @x. The search starts at the root and,
 * at each split node, counts by bisection the interior edges at or below the coordinate of
 * @x in the variable of the split node, and continues at the child of that position. A
 * point on an interior edge belongs to the child whose interval starts at that edge, and a
 * point outside the domain reaches the panel that contains the point obtained by clamping
 * each coordinate to the domain. The call is an error
 * when @n_vars is not #NcmSpectralTensor:n-vars, the domain is not set, or the search
 * reaches an empty node.
 *
 * Returns: the index of the panel.
 */
guint
ncm_spectral_tensor_find_panel (NcmSpectralTensor *st, guint n_vars, const gdouble *x)
{
  const NcmSpectralTensorNode *nodes = (const NcmSpectralTensorNode *) st->nodes->data;
  const gdouble *edges               = (const gdouble *) st->edges->data;
  guint node                         = 0;

  if (n_vars != st->n_vars)
    g_error ("ncm_spectral_tensor_find_panel: a point of %u values for %u variables.", n_vars, st->n_vars);

  if (st->nodes->len == 0)
    g_error ("ncm_spectral_tensor_find_panel: the domain is not set.");

  while (nodes[node].kind == NCM_SPECTRAL_TENSOR_NODE_KIND_SPLIT)
  {
    const NcmSpectralTensorNode *nd = &nodes[node];
    const gdouble *e                = edges + nd->edges_offset;
    const gdouble xi                = x[nd->var];
    guint lo                        = 0;
    guint hi                        = nd->n_edges;

    /* The number of interior edges at or below xi, the index of the child */
    while (lo < hi)
    {
      const guint mid = (lo + hi) / 2;

      if (xi >= e[mid])
        lo = mid + 1;
      else
        hi = mid;
    }

    node = nd->index + lo;
  }

  if (nodes[node].kind != NCM_SPECTRAL_TENSOR_NODE_KIND_PANEL)
    g_error ("ncm_spectral_tensor_find_panel: the point is in the empty node %u.", node);

  return nodes[node].index;
}

/**
 * ncm_spectral_tensor_eval:
 * @st: a #NcmSpectralTensor
 * @n_vars: number of variables, #NcmSpectralTensor:n-vars
 * @x: (array length=n_vars): the point
 *
 * In a filled table, evaluates the function at @x on the panel given by
 * ncm_spectral_tensor_find_panel(), with its errors. Outside the domain the series of the
 * panel that contains the point obtained by clamping each coordinate to the domain is
 * extrapolated, without an error and outside the tolerance of the table: the caller
 * compares the point with the range of ncm_spectral_tensor_get_range() before the call.
 *
 * Returns: $F (x)$.
 */
gdouble
ncm_spectral_tensor_eval (NcmSpectralTensor *st, guint n_vars, const gdouble *x)
{
  return ncm_spectral_tensor_eval_on (st, ncm_spectral_tensor_find_panel (st, n_vars, x), n_vars, x);
}

/* The series sum_{k < K} c[k] T_k (s) by the Clenshaw recurrence */
static inline __attribute__ ((always_inline)) gdouble

_ncm_spectral_tensor_clenshaw_1 (const gdouble *c, guint K, gdouble s)
{
  const gdouble s2 = 2.0 * s;
  gdouble b1       = 0.0;
  gdouble b2       = 0.0;
  guint k;

  for (k = K - 1; k > 0; k--)
  {
    const gdouble b0 = s2 * b1 - b2 + c[k];

    b2 = b1;
    b1 = b0;
  }

  return s * b1 - b2 + c[0];
}

/*
 * The series of a panel in the variables 0, ..., d, d >= 1, at s[0..d]: the Clenshaw
 * recurrence in variable d, whose coefficient k is the series in the variables below d of
 * the block at c + k stride[d]. Used for five variables or more, and for a panel with more
 * than NCM_SPECTRAL_TENSOR_BASIS_K_MAX coefficients in a variable.
 */
static gdouble
_ncm_spectral_tensor_clenshaw (const gdouble *c, const guint *K, const gsize *stride, const gdouble *s, guint d)
{
  const gdouble s2 = 2.0 * s[d];
  gdouble b1       = 0.0;
  gdouble b2       = 0.0;
  guint k;

  for (k = K[d] - 1; k > 0; k--)
  {
    const gdouble *ck = c + (gsize) k * stride[d];
    const gdouble a   = (d == 1) ? _ncm_spectral_tensor_clenshaw_1 (ck, K[0], s[0]) : _ncm_spectral_tensor_clenshaw (ck, K, stride, s, d - 1);
    const gdouble b0  = s2 * b1 - b2 + a;

    b2 = b1;
    b1 = b0;
  }

  {
    const gdouble a = (d == 1) ? _ncm_spectral_tensor_clenshaw_1 (c, K[0], s[0]) : _ncm_spectral_tensor_clenshaw (c, K, stride, s, d - 1);

    return s[d] * b1 - b2 + a;
  }
}

/*
 * Largest number of coefficients in a variable of the sums with the values T_k (s_i): that
 * of a panel of NcmSpectralTensorBuilder, 2^10 + 1. The values are on the stack.
 */
#define NCM_SPECTRAL_TENSOR_BASIS_K_MAX (1025)

/*
 * Largest number of values T_k (s) computed by the product formulas or by four chains; above
 * it the three-term recurrence is used. Measured error of the values near s = +-1, against
 * the three-term recurrence: equal up to 17 values, twice at 33, five times at 129. It is
 * also the largest row length with its own copy of the sums.
 */
#define NCM_SPECTRAL_TENSOR_BASIS_K_FAST (33)

/* T_0 (s), ..., T_{K - 1} (s) by the three-term recurrence */
static inline __attribute__ ((always_inline)) void

_ncm_spectral_tensor_cheb_T (gdouble *T, const guint K, const gdouble s)
{
  const gdouble s2 = 2.0 * s;
  guint k;

  T[0] = 1.0;

  if (K > 1)
    T[1] = s;

  for (k = 2; k < K; k++)
    T[k] = s2 * T[k - 1] - T[k - 2];
}

/*
 * T_0 (s), ..., T_{K - 1} (s) by the product formulas T_2k = 2 T_k^2 - 1 and
 * T_2k+1 = 2 T_k T_k+1 - s: the value k depends on about log2 (k) earlier ones. For K a
 * constant at the call site, where the loop unrolls.
 */
static inline __attribute__ ((always_inline)) void

_ncm_spectral_tensor_cheb_T_prod (gdouble *T, const guint K, const gdouble s)
{
  guint k;

  T[0] = 1.0;

  if (K > 1)
    T[1] = s;

  for (k = 2; k < K; k++)
  {
    const guint h = k / 2;

    T[k] = (k & 1) ? 2.0 * T[h] * T[h + 1] - s : 2.0 * T[h] * T[h] - 1.0;
  }
}

/*
 * T_0 (s), ..., T_{K - 1} (s) in four chains held in registers, T_k+4 = 2 T_4 T_k - T_k-4
 * with T_-k = T_k for the first step: four values per dependent step. Writes up to
 * 4 ceil (K / 4) values, at least four, so T has room for three values more than K.
 */
static inline __attribute__ ((always_inline)) void

_ncm_spectral_tensor_cheb_T_chains (gdouble *T, const guint K, const gdouble s)
{
  const gdouble s2 = 2.0 * s;
  gdouble a0       = 1.0;
  gdouble a1       = s;
  gdouble a2       = s2 * s - 1.0;
  gdouble a3       = s2 * a2 - s;
  const gdouble T4 = 2.0 * a2 * a2 - 1.0;
  const gdouble m4 = 2.0 * T4;
  gdouble b0       = T4;
  gdouble b1       = m4 * a1 - a3;
  gdouble b2       = m4 * a2 - a2;
  gdouble b3       = m4 * a3 - a1;
  guint k;

  T[0] = a0;
  T[1] = a1;
  T[2] = a2;
  T[3] = a3;

  for (k = 4; k < K; k += 4)
  {
    const gdouble c0 = m4 * b0 - a0;
    const gdouble c1 = m4 * b1 - a1;
    const gdouble c2 = m4 * b2 - a2;
    const gdouble c3 = m4 * b3 - a3;

    T[k]     = b0;
    T[k + 1] = b1;
    T[k + 2] = b2;
    T[k + 3] = b3;

    a0 = b0;
    a1 = b1;
    a2 = b2;
    a3 = b3;
    b0 = c0;
    b1 = c1;
    b2 = c2;
    b3 = c3;
  }
}

/* T_0 (s), ..., T_{K - 1} (s) for K known at run time: four chains up to K_FAST values */
static inline __attribute__ ((always_inline)) void

_ncm_spectral_tensor_cheb_T_any (gdouble *T, const guint K, const gdouble s)
{
  if (K <= NCM_SPECTRAL_TENSOR_BASIS_K_FAST)
    _ncm_spectral_tensor_cheb_T_chains (T, K, s);
  else
    _ncm_spectral_tensor_cheb_T (T, K, s);
}

/*
 * The series of a block of two variables, sum_k1 T1[k1] (row k1 . T0), from the values
 * T0[k0] = T_k0 (s_0) and T1[k1] = T_k1 (s_1). A row is a dot product with T0, four rows in
 * one pass over the columns: the dot products are independent, and a pass is one loop.
 * With K0 a constant at the call site the loops have no remainder.
 */
static inline __attribute__ ((always_inline)) gdouble

_ncm_spectral_tensor_basis_2 (const gdouble *c, const guint K0, const guint K1, const gdouble *T0, const gdouble *T1)
{
  gdouble acc = 0.0;
  guint j, k;

  for (k = 0; k + 4 <= K1; k += 4)
  {
    const gdouble *r0 = c + (gsize) k * K0;
    const gdouble *r1 = r0 + K0;
    const gdouble *r2 = r1 + K0;
    const gdouble *r3 = r2 + K0;
    gdouble d0        = 0.0;
    gdouble d1        = 0.0;
    gdouble d2        = 0.0;
    gdouble d3        = 0.0;

    #pragma omp simd reduction(+:d0,d1,d2,d3)

    for (j = 0; j < K0; j++)
    {
      const gdouble t = T0[j];

      d0 += r0[j] * t;
      d1 += r1[j] * t;
      d2 += r2[j] * t;
      d3 += r3[j] * t;
    }

    acc += T1[k] * d0 + T1[k + 1] * d1 + T1[k + 2] * d2 + T1[k + 3] * d3;
  }

  for ( ; k < K1; k++)
  {
    const gdouble *row = c + (gsize) k * K0;
    gdouble d          = 0.0;

    #pragma omp simd reduction(+:d)

    for (j = 0; j < K0; j++)
      d += row[j] * T0[j];

    acc += T1[k] * d;
  }

  return acc;
}

typedef void (*NcmSpectralTensorChebFunc) (gdouble *T, const gdouble s);
typedef gdouble (*NcmSpectralTensorBlockFunc) (const gdouble *c, const guint K0, const guint K1, const gdouble *T0, const gdouble *T1);

/*
 * The copies for the row lengths k0 = 1, ..., NCM_SPECTRAL_TENSOR_BASIS_K_FAST: the values
 * T_0 (s), ..., T_{k0 - 1} (s) by the product formulas, and
 * _ncm_spectral_tensor_basis_2() with K0 = k0, which ignores its argument K0. The tables
 * have the copies of row length k0 at index k0 - 1.
 */
/* *INDENT-OFF* */
#define NCM_SPECTRAL_TENSOR_BASIS_K0_LIST(X) \
  X (1) X (2) X (3) X (4) X (5) X (6) X (7) X (8) X (9) X (10) X (11) \
  X (12) X (13) X (14) X (15) X (16) X (17) X (18) X (19) X (20) X (21) X (22) \
  X (23) X (24) X (25) X (26) X (27) X (28) X (29) X (30) X (31) X (32) X (33)

#define NCM_SPECTRAL_TENSOR_BASIS_FIXED(k0) \
  static void \
  _ncm_spectral_tensor_cheb_T_ ## k0 (gdouble *T, const gdouble s) \
  { \
    _ncm_spectral_tensor_cheb_T_prod (T, k0, s); \
  } \
  static gdouble \
  _ncm_spectral_tensor_block_ ## k0 (const gdouble *c, const guint K0, const guint K1, const gdouble *T0, const gdouble *T1) \
  { \
    return _ncm_spectral_tensor_basis_2 (c, k0, K1, T0, T1); \
  }

#define NCM_SPECTRAL_TENSOR_BASIS_CHEB_ENTRY(k0) &_ncm_spectral_tensor_cheb_T_ ## k0,
#define NCM_SPECTRAL_TENSOR_BASIS_BLOCK_ENTRY(k0) &_ncm_spectral_tensor_block_ ## k0,

NCM_SPECTRAL_TENSOR_BASIS_K0_LIST (NCM_SPECTRAL_TENSOR_BASIS_FIXED)

static const NcmSpectralTensorChebFunc _ncm_spectral_tensor_cheb_T_fixed[NCM_SPECTRAL_TENSOR_BASIS_K_FAST] = {
  NCM_SPECTRAL_TENSOR_BASIS_K0_LIST (NCM_SPECTRAL_TENSOR_BASIS_CHEB_ENTRY)
};

static const NcmSpectralTensorBlockFunc _ncm_spectral_tensor_block_fixed[NCM_SPECTRAL_TENSOR_BASIS_K_FAST] = {
  NCM_SPECTRAL_TENSOR_BASIS_K0_LIST (NCM_SPECTRAL_TENSOR_BASIS_BLOCK_ENTRY)
};
/* *INDENT-ON* */

/* _ncm_spectral_tensor_basis_2() for a row length above NCM_SPECTRAL_TENSOR_BASIS_K_FAST */
static gdouble
_ncm_spectral_tensor_block_any (const gdouble *c, const guint K0, const guint K1, const gdouble *T0, const gdouble *T1)
{
  return _ncm_spectral_tensor_basis_2 (c, K0, K1, T0, T1);
}

/*
 * The coefficients of panel p, with its numbers of coefficients in *K and the Chebyshev
 * variables of x on it in s. Always inlined, so that n_vars is a constant in each caller.
 */
static inline __attribute__ ((always_inline)) const gdouble *

_ncm_spectral_tensor_panel_at (NcmSpectralTensor *st, guint p, const guint n_vars, const gdouble *x, const guint **K, gdouble *s)
{
  const NcmSpectralTensorPanel *pn = &g_array_index (st->panels, NcmSpectralTensorPanel, p);
  const gdouble *lower             = &g_array_index (st->node_bounds, gdouble, (gsize) 2 * n_vars * pn->node);
  const gdouble *upper             = lower + n_vars;
  guint i;

  for (i = 0; i < n_vars; i++)
    s[i] = (2.0 * x[i] - lower[i] - upper[i]) / (upper[i] - lower[i]);

  *K = &g_array_index (st->n_coeffs, guint, (gsize) n_vars * p);

  return &g_array_index (st->coeffs, gdouble, pn->offset);
}

/* The series of panel p at x by nested Clenshaw recurrences, for any number of variables */
static gdouble
_ncm_spectral_tensor_eval_on_clenshaw (NcmSpectralTensor *st, guint p, guint n_vars, const gdouble *x)
{
  const guint *K;
  gdouble s[NCM_SPECTRAL_TENSOR_MAX_VARS];
  gsize stride[NCM_SPECTRAL_TENSOR_MAX_VARS];
  const gdouble *c = _ncm_spectral_tensor_panel_at (st, p, n_vars, x, &K, s);
  guint i;

  if (n_vars == 1)
    return _ncm_spectral_tensor_clenshaw_1 (c, K[0], s[0]);

  stride[0] = 1;

  for (i = 1; i < n_vars; i++)
    stride[i] = stride[i - 1] * K[i - 1];

  return _ncm_spectral_tensor_clenshaw (c, K, stride, s, n_vars - 1);
}

/*
 * The series of panel p of a table of n_vars variables, two to four, at x. The values
 * T_k (s_i) of every variable are computed first, those of the first variable by the copy of
 * its row length. The series is then the sum over the blocks of the last variables of
 * T_k (s_i) times the series of the block, the blocks of two variables summed by the copy of
 * the row length. A variable with more than NCM_SPECTRAL_TENSOR_BASIS_K_MAX coefficients
 * sends the panel to the Clenshaw recurrences. Always inlined, n_vars a constant.
 */
static inline __attribute__ ((always_inline)) gdouble

_ncm_spectral_tensor_eval_on_basis (NcmSpectralTensor *st, guint p, const guint n_vars, const gdouble *x)
{
  const guint *K;
  gdouble s[4];
  const gdouble *c = _ncm_spectral_tensor_panel_at (st, p, n_vars, x, &K, s);
  gdouble T0[NCM_SPECTRAL_TENSOR_BASIS_K_MAX + 3];
  gdouble T1[NCM_SPECTRAL_TENSOR_BASIS_K_MAX + 3];
  NcmSpectralTensorBlockFunc block;
  guint K0, i;

  for (i = 0; i < n_vars; i++)
  {
    if (K[i] > NCM_SPECTRAL_TENSOR_BASIS_K_MAX)
      return _ncm_spectral_tensor_eval_on_clenshaw (st, p, n_vars, x);
  }

  K0 = K[0];

  if (K0 <= NCM_SPECTRAL_TENSOR_BASIS_K_FAST)
  {
    _ncm_spectral_tensor_cheb_T_fixed[K0 - 1](T0, s[0]);
    block = _ncm_spectral_tensor_block_fixed[K0 - 1];
  }
  else
  {
    _ncm_spectral_tensor_cheb_T (T0, K0, s[0]);
    block = &_ncm_spectral_tensor_block_any;
  }

  _ncm_spectral_tensor_cheb_T_any (T1, K[1], s[1]);

  if (n_vars == 2)
  {
    return block (c, K0, K[1], T0, T1);
  }
  else
  {
    const gsize stride2 = (gsize) K0 * K[1];
    gdouble T2[NCM_SPECTRAL_TENSOR_BASIS_K_MAX + 3];
    gdouble acc2 = 0.0;
    guint k2;

    _ncm_spectral_tensor_cheb_T_any (T2, K[2], s[2]);

    if (n_vars == 3)
    {
      for (k2 = 0; k2 < K[2]; k2++)
        acc2 += T2[k2] * block (c + (gsize) k2 * stride2, K0, K[1], T0, T1);

      return acc2;
    }
    else
    {
      const gsize stride3 = stride2 * K[2];
      gdouble T3[NCM_SPECTRAL_TENSOR_BASIS_K_MAX + 3];
      gdouble acc3 = 0.0;
      guint k3;

      _ncm_spectral_tensor_cheb_T_any (T3, K[3], s[3]);

      for (k3 = 0; k3 < K[3]; k3++)
      {
        const gdouble *c3 = c + (gsize) k3 * stride3;

        acc2 = 0.0;

        for (k2 = 0; k2 < K[2]; k2++)
          acc2 += T2[k2] * block (c3 + (gsize) k2 * stride2, K0, K[1], T0, T1);

        acc3 += T3[k3] * acc2;
      }

      return acc3;
    }
  }
}

/**
 * ncm_spectral_tensor_eval_on:
 * @st: a #NcmSpectralTensor
 * @p: panel index
 * @n_vars: number of variables, #NcmSpectralTensor:n-vars
 * @x: (array length=n_vars): the point
 *
 * For a point whose panel is known, evaluates the series of panel @p at @x. With one
 * variable the series is summed by the Clenshaw recurrence. With two to four variables and
 * at most 1025 coefficients in each, it is summed from the values $T_k (s_i)$, computed by
 * recurrences of the Chebyshev polynomials. Otherwise it is summed by Clenshaw recurrences,
 * nested from the last variable, outermost, to the first. The index @p must be below
 * ncm_spectral_tensor_get_n_panels() and @n_vars must be #NcmSpectralTensor:n-vars;
 * neither is checked.
 *
 * Returns: the series of panel @p at @x.
 */
gdouble
ncm_spectral_tensor_eval_on (NcmSpectralTensor *st, guint p, guint n_vars, const gdouble *x)
{
  switch (n_vars)
  {
    case 1:
    {
      const guint *K;
      gdouble s[1];
      const gdouble *c = _ncm_spectral_tensor_panel_at (st, p, 1, x, &K, s);

      return _ncm_spectral_tensor_clenshaw_1 (c, K[0], s[0]);
    }

    case 2:
      return _ncm_spectral_tensor_eval_on_basis (st, p, 2, x);

    case 3:
      return _ncm_spectral_tensor_eval_on_basis (st, p, 3, x);

    case 4:
      return _ncm_spectral_tensor_eval_on_basis (st, p, 4, x);

    default:
      return _ncm_spectral_tensor_eval_on_clenshaw (st, p, n_vars, x);
  }
}

/*
 * Largest number of coefficients in a variable of the sums of several points at once: that
 * of a panel at the default level-cap of NcmSpectralTensorBuilder, 2^7 + 1. The values
 * T_k (s_i) of the points of a pass are on the stack.
 */
#define NCM_SPECTRAL_TENSOR_LANES_K_MAX (129)

/*
 * The series of a panel of n_vars variables, one to four, at M points at once: s[i][j] is
 * the Chebyshev variable i of point j. With one variable, M Clenshaw recurrences. With two
 * to four, the values T[i][k][j] = T_k (s_i) of point j by the three-term recurrence, then
 * the sums with every coefficient loaded once and multiplied into the M points: the loops
 * over the points are over M entries, a constant at the call site, as n_vars is. Every sum
 * runs from the last index to the first, so that the terms of a decaying series are added
 * from the smallest: the rounding error is then that of the Clenshaw recurrences.
 */
static inline __attribute__ ((always_inline)) void

_ncm_spectral_tensor_lanes (const gdouble *c, const guint n_vars, const guint M, const guint *K, gdouble s[4][8], gdouble *out)
{
  guint i, j, k;

  if (n_vars == 1)
  {
    gdouble b1[8];
    gdouble b2[8];

    for (j = 0; j < M; j++)
      b1[j] = b2[j] = 0.0;

    for (k = K[0] - 1; k > 0; k--)
    {
      const gdouble ck = c[k];

      for (j = 0; j < M; j++)
      {
        const gdouble b0 = 2.0 * s[0][j] * b1[j] - b2[j] + ck;

        b2[j] = b1[j];
        b1[j] = b0;
      }
    }

    for (j = 0; j < M; j++)
      out[j] = s[0][j] * b1[j] - b2[j] + c[0];
  }
  else
  {
    const guint K0 = K[0];
    const guint K1 = K[1];
    const guint K2 = (n_vars > 2) ? K[2] : 1;
    const guint K3 = (n_vars > 3) ? K[3] : 1;
    gdouble T[4][NCM_SPECTRAL_TENSOR_LANES_K_MAX][8];
    gdouble acc3[8];
    guint k1, k2, k3;

    for (i = 0; i < n_vars; i++)
    {
      for (j = 0; j < M; j++)
      {
        T[i][0][j] = 1.0;
        T[i][1][j] = s[i][j];
      }

      for (k = 2; k < K[i]; k++)
      {
        for (j = 0; j < M; j++)
          T[i][k][j] = 2.0 * s[i][j] * T[i][k - 1][j] - T[i][k - 2][j];
      }
    }

    for (j = 0; j < M; j++)
      acc3[j] = 0.0;

    for (k3 = K3; k3-- > 0; )
    {
      gdouble acc2[8];

      for (j = 0; j < M; j++)
        acc2[j] = 0.0;

      for (k2 = K2; k2-- > 0; )
      {
        const gdouble *block = c + ((gsize) k3 * K2 + k2) * K1 * K0;
        gdouble acc1[8];

        for (j = 0; j < M; j++)
          acc1[j] = 0.0;

        for (k1 = K1; k1-- > 0; )
        {
          const gdouble *row = block + (gsize) k1 * K0;
          gdouble u[8];

          for (j = 0; j < M; j++)
            u[j] = 0.0;

          for (k = K0; k-- > 0; )
          {
            const gdouble ck = row[k];

            for (j = 0; j < M; j++)
              u[j] += ck * T[0][k][j];
          }

          for (j = 0; j < M; j++)
            acc1[j] += T[1][k1][j] * u[j];
        }

        if (n_vars > 2)
          for (j = 0; j < M; j++)
            acc2[j] += T[2][k2][j] * acc1[j];

        else
          for (j = 0; j < M; j++)
            acc2[j] = acc1[j];
      }

      if (n_vars > 3)
        for (j = 0; j < M; j++)
          acc3[j] += T[3][k3][j] * acc2[j];

      else
        for (j = 0; j < M; j++)
          acc3[j] = acc2[j];
    }

    for (j = 0; j < M; j++)
      out[j] = acc3[j];
  }
}

/*
 * The series of panel p at the n consecutive points of x, for n_vars from one to four and at
 * most NCM_SPECTRAL_TENSOR_LANES_K_MAX coefficients in each variable: passes of 8 points,
 * then one of 4, then the points left one at a time. Always inlined, n_vars a constant.
 */
static inline __attribute__ ((always_inline)) void

_ncm_spectral_tensor_eval_on_run (NcmSpectralTensor *st, guint p, const guint n_vars, guint n, const gdouble *x, gdouble *fval)
{
  const NcmSpectralTensorPanel *pn = &g_array_index (st->panels, NcmSpectralTensorPanel, p);
  const gdouble *lower             = &g_array_index (st->node_bounds, gdouble, (gsize) 2 * n_vars * pn->node);
  const gdouble *upper             = lower + n_vars;
  const guint *K                   = &g_array_index (st->n_coeffs, guint, (gsize) n_vars * p);
  const gdouble *c                 = &g_array_index (st->coeffs, gdouble, pn->offset);
  gdouble s[4][8];
  guint i, j;

  while (n >= 8)
  {
    for (i = 0; i < n_vars; i++)
    {
      for (j = 0; j < 8; j++)
        s[i][j] = (2.0 * x[n_vars * j + i] - lower[i] - upper[i]) / (upper[i] - lower[i]);
    }

    _ncm_spectral_tensor_lanes (c, n_vars, 8, K, s, fval);

    x    += 8 * n_vars;
    fval += 8;
    n    -= 8;
  }

  if (n >= 4)
  {
    for (i = 0; i < n_vars; i++)
    {
      for (j = 0; j < 4; j++)
        s[i][j] = (2.0 * x[n_vars * j + i] - lower[i] - upper[i]) / (upper[i] - lower[i]);
    }

    _ncm_spectral_tensor_lanes (c, n_vars, 4, K, s, fval);

    x    += 4 * n_vars;
    fval += 4;
    n    -= 4;
  }

  for (j = 0; j < n; j++)
    fval[j] = ncm_spectral_tensor_eval_on (st, p, n_vars, &x[n_vars * j]);
}

/* _ncm_spectral_tensor_eval_on_run() for the number of variables of st, one copy for each */
static void
_ncm_spectral_tensor_eval_on_array (NcmSpectralTensor *st, guint p, guint n, const gdouble *x, gdouble *fval)
{
  const guint n_vars = st->n_vars;
  const guint *K     = &g_array_index (st->n_coeffs, guint, (gsize) n_vars * p);
  gboolean lanes     = (n_vars <= 4) && (n >= 4);
  guint i, j;

  for (i = 0; lanes && (i < n_vars); i++)
    lanes = (K[i] <= NCM_SPECTRAL_TENSOR_LANES_K_MAX);

  if (!lanes)
  {
    for (j = 0; j < n; j++)
      fval[j] = ncm_spectral_tensor_eval_on (st, p, n_vars, &x[n_vars * j]);

    return;
  }

  switch (n_vars)
  {
    case 1:
      _ncm_spectral_tensor_eval_on_run (st, p, 1, n, x, fval);
      break;

    case 2:
      _ncm_spectral_tensor_eval_on_run (st, p, 2, n, x, fval);
      break;

    case 3:
      _ncm_spectral_tensor_eval_on_run (st, p, 3, n, x, fval);
      break;

    default:
      _ncm_spectral_tensor_eval_on_run (st, p, 4, n, x, fval);
      break;
  }
}

/**
 * ncm_spectral_tensor_eval_array:
 * @st: a #NcmSpectralTensor
 * @x: (element-type gdouble): the points, the variables of a point together
 * @fval: (out callee-allocates) (transfer full) (element-type gdouble): the values at the points
 *
 * In a filled table, evaluates the function at several points: point $j$ is the
 * #NcmSpectralTensor:n-vars entries of @x that start at $j$ times that number, and its
 * value is entry $j$ of @fval. A C caller may pass in @fval an array to reuse, or a pointer
 * to %NULL for a new one. The panel of each point is the one of
 * ncm_spectral_tensor_find_panel(), with its errors, and a point outside the domain is
 * extrapolated as in ncm_spectral_tensor_eval(). Consecutive points on the same panel
 * are summed together, 8 points at a time, which takes less time per point than
 * ncm_spectral_tensor_eval(): points ordered along a line, or grouped by panel, gain the
 * most. The values agree with those of ncm_spectral_tensor_eval() to rounding, not bit by
 * bit. The call is an error when the length of @x is not a multiple of
 * #NcmSpectralTensor:n-vars.
 */
void
ncm_spectral_tensor_eval_array (NcmSpectralTensor *st, GArray *x, GArray **fval)
{
  const guint n_vars   = st->n_vars;
  const guint n_points = x->len / n_vars;
  const gdouble *xd    = &g_array_index (x, gdouble, 0);
  gdouble *fd;
  guint j = 0;
  guint p;

  if (x->len != n_points * n_vars)
    g_error ("ncm_spectral_tensor_eval_array: %u values are not points of %u variables.", x->len, n_vars);

  if (*fval == NULL)
    *fval = g_array_new (FALSE, FALSE, sizeof (gdouble));

  g_array_set_size (*fval, n_points);

  if (n_points == 0)
    return;

  fd = &g_array_index (*fval, gdouble, 0);
  p  = ncm_spectral_tensor_find_panel (st, n_vars, xd);

  while (j < n_points)
  {
    guint p_next = p;
    guint run    = 1;

    /* The run of points after j on the panel of j, and the panel of the point that ends it */
    while (j + run < n_points)
    {
      p_next = ncm_spectral_tensor_find_panel (st, n_vars, &xd[(gsize) n_vars * (j + run)]);

      if (p_next != p)
        break;

      run++;
    }

    _ncm_spectral_tensor_eval_on_array (st, p, run, &xd[(gsize) n_vars * j], &fd[j]);

    j += run;
    p  = p_next;
  }
}

