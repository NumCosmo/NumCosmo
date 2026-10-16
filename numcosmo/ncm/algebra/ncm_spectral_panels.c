/***************************************************************************
 *            ncm_spectral_panels.c
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

/**
 * NcmSpectralPanels:
 *
 * A vector function on an interval as contiguous Chebyshev panels.
 *
 * The edges $x_0 < x_1 < \dots < x_n$ split $[x_0, x_n]$ into $n$ panels. On panel $i$,
 * $[x_i, x_{i+1}]$, each entry $F_r$ of $F$ is a Chebyshev series in the variable $s \in [-1,
 * 1]$ of ncm_spectral_x_to_s(), with its coefficients in row $r$ of an
 * #NcmSpectralPanels:fdim by $N_i$ block, stored row after row; $N_i$ may differ between panels. A point on an
 * interior edge belongs to the panel that starts there, so where $F$ jumps at an edge its
 * value there is the limit from the right.
 *
 * Panels are added in ascending order by ncm_spectral_panels_append(). Each panel's
 * coefficients live in an array that ncm_spectral_panels_reset() keeps for reuse, so a
 * panel set rebuilt at every evaluation of a likelihood stops allocating once it has grown. Reading
 * and evaluating do not change the object.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/algebra/ncm_spectral_panels.h"
#include <math.h>
#include <float.h>
#include <string.h>

enum
{
  PROP_0,
  PROP_FDIM,
  PROP_PANELS,
};

struct _NcmSpectralPanels
{
  /*< private >*/
  GObject parent_instance;
  guint fdim;
  guint n_panels;
  GArray *edges;         /* gdouble, n_panels + 1 entries once a panel is added */
  GPtrArray *pool;       /* GArray of gdouble, the first n_panels hold the panels' fdim rows of N_i */
  GVariant *init_panels; /* the construct-only property, consumed by constructed() */
};

G_DEFINE_TYPE (NcmSpectralPanels, ncm_spectral_panels, G_TYPE_OBJECT)

static void
ncm_spectral_panels_init (NcmSpectralPanels *sp)
{
  sp->fdim        = 0;
  sp->n_panels    = 0;
  sp->edges       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  sp->pool        = g_ptr_array_new_with_free_func ((GDestroyNotify) g_array_unref);
  sp->init_panels = NULL;
}

static void
_ncm_spectral_panels_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSpectralPanels *sp = NCM_SPECTRAL_PANELS (object);

  g_return_if_fail (NCM_IS_SPECTRAL_PANELS (object));

  switch (prop_id)
  {
    case PROP_FDIM:
      sp->fdim = g_value_get_uint (value);
      break;
    case PROP_PANELS:
      sp->init_panels = g_value_dup_variant (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static GVariant *_ncm_spectral_panels_get_variant (NcmSpectralPanels *sp);

static void
_ncm_spectral_panels_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSpectralPanels *sp = NCM_SPECTRAL_PANELS (object);

  g_return_if_fail (NCM_IS_SPECTRAL_PANELS (object));

  switch (prop_id)
  {
    case PROP_FDIM:
      g_value_set_uint (value, sp->fdim);
      break;
    case PROP_PANELS:
      g_value_take_variant (value, _ncm_spectral_panels_get_variant (sp));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

/* The edges and the coefficients of every panel as a (adaad). The arrays are copied, not
 * shared as by ncm_cfg_array_to_variant(): reset and append rewrite them in place, and a
 * variant must not change after it is made. */
static GVariant *
_ncm_spectral_panels_get_variant (NcmSpectralPanels *sp)
{
  GVariantBuilder builder;
  GVariant *edges;
  guint i;

  edges = g_variant_new_fixed_array (G_VARIANT_TYPE ("d"), sp->edges->data, sp->edges->len, sizeof (gdouble));

  g_variant_builder_init (&builder, G_VARIANT_TYPE ("aad"));

  for (i = 0; i < sp->n_panels; i++)
  {
    GArray *coeffs = g_ptr_array_index (sp->pool, i);

    g_variant_builder_add_value (&builder, g_variant_new_fixed_array (G_VARIANT_TYPE ("d"), coeffs->data, coeffs->len, sizeof (gdouble)));
  }

  return g_variant_ref_sink (g_variant_new ("(@ad@aad)", edges, g_variant_builder_end (&builder)));
}

/* The panels given as a property, added in order so that each is checked as an appended one */
static void
_ncm_spectral_panels_constructed (GObject *object)
{
  NcmSpectralPanels *sp = NCM_SPECTRAL_PANELS (object);

  /* Chain up : start */
  G_OBJECT_CLASS (ncm_spectral_panels_parent_class)->constructed (object);

  if (sp->init_panels != NULL)
  {
    GVariant *edges_var  = g_variant_get_child_value (sp->init_panels, 0);
    GVariant *panels_var = g_variant_get_child_value (sp->init_panels, 1);
    const guint n_panels = g_variant_n_children (panels_var);
    gsize n_edges        = 0;
    const gdouble *edges = g_variant_get_fixed_array (edges_var, &n_edges, sizeof (gdouble));
    guint i;

    if (((n_edges > 0) || (n_panels > 0)) && (n_edges != n_panels + 1))
      g_error ("ncm_spectral_panels: %zu edges for %u panels, expected %u.", n_edges, n_panels, n_panels + 1);

    for (i = 0; i < n_panels; i++)
    {
      GVariant *coeffs_var = g_variant_get_child_value (panels_var, i);
      gsize len            = 0;
      const gdouble *c     = g_variant_get_fixed_array (coeffs_var, &len, sizeof (gdouble));

      if (len % sp->fdim != 0)
        g_error ("ncm_spectral_panels: panel %u has %zu coefficients, not a multiple of fdim = %u.", i, len, sp->fdim);

      ncm_spectral_panels_append_data (sp, edges[i], edges[i + 1], c, len / sp->fdim, len / sp->fdim);
      g_variant_unref (coeffs_var);
    }

    g_variant_unref (edges_var);
    g_variant_unref (panels_var);
    g_clear_pointer (&sp->init_panels, g_variant_unref);
  }
}

static void
_ncm_spectral_panels_finalize (GObject *object)
{
  NcmSpectralPanels *sp = NCM_SPECTRAL_PANELS (object);

  g_array_unref (sp->edges);
  g_ptr_array_unref (sp->pool);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spectral_panels_parent_class)->finalize (object);
}

static void
ncm_spectral_panels_class_init (NcmSpectralPanelsClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_spectral_panels_set_property;
  object_class->get_property = &_ncm_spectral_panels_get_property;
  object_class->constructed  = &_ncm_spectral_panels_constructed;
  object_class->finalize     = &_ncm_spectral_panels_finalize;

  /**
   * NcmSpectralPanels:fdim:
   *
   * The dimension of the values of $F$, the number of rows of every panel.
   */
  g_object_class_install_property (object_class,
                                   PROP_FDIM,
                                   g_param_spec_uint ("fdim",
                                                      NULL,
                                                      "Dimension of the values",
                                                      1, G_MAXUINT, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSpectralPanels:panels:
   *
   * The panels as a #GVariant of type `(adaad)`: the edges in ascending order, then for
   * each panel its #NcmSpectralPanels:fdim rows of coefficients one after the other. Both
   * arrays are empty when there are no panels.
   */
  g_object_class_install_property (object_class,
                                   PROP_PANELS,
                                   g_param_spec_variant ("panels",
                                                         NULL,
                                                         "Panel edges and coefficients",
                                                         G_VARIANT_TYPE ("(adaad)"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

static gdouble _ncm_spectral_panels_clenshaw (const gdouble *c, guint N, const gdouble s);

/**
 * ncm_spectral_panels_new:
 * @fdim: dimension of the values of $F$
 *
 * Creates a panel set with no panels.
 *
 * Returns: (transfer full): a new #NcmSpectralPanels.
 */
NcmSpectralPanels *
ncm_spectral_panels_new (guint fdim)
{
  return g_object_new (NCM_TYPE_SPECTRAL_PANELS,
                       "fdim", fdim,
                       NULL);
}

/**
 * ncm_spectral_panels_ref:
 * @sp: a #NcmSpectralPanels
 *
 * Increases the reference count of @sp by one.
 *
 * Returns: (transfer full): @sp.
 */
NcmSpectralPanels *
ncm_spectral_panels_ref (NcmSpectralPanels *sp)
{
  return g_object_ref (sp);
}

/**
 * ncm_spectral_panels_free:
 * @sp: a #NcmSpectralPanels
 *
 * Decreases the reference count of @sp by one.
 */
void
ncm_spectral_panels_free (NcmSpectralPanels *sp)
{
  g_object_unref (sp);
}

/**
 * ncm_spectral_panels_clear:
 * @sp: a #NcmSpectralPanels
 *
 * Decreases the reference count of *@sp by one and sets *@sp to %NULL.
 */
void
ncm_spectral_panels_clear (NcmSpectralPanels **sp)
{
  g_clear_object (sp);
}

/**
 * ncm_spectral_panels_append:
 * @sp: a #NcmSpectralPanels
 * @a: left edge of the panel
 * @b: right edge of the panel
 * @coeffs: the #NcmSpectralPanels:fdim by $N$ coefficients of the panel
 *
 * Adds the panel $[a, b]$ after the last one, copying @coeffs. It must start where the last
 * panel ends, @a must be below @b, and @coeffs must have #NcmSpectralPanels:fdim rows;
 * otherwise this is an error.
 */
void
ncm_spectral_panels_append (NcmSpectralPanels *sp, gdouble a, gdouble b, NcmMatrix *coeffs)
{
  if (ncm_matrix_nrows (coeffs) != sp->fdim)
    g_error ("ncm_spectral_panels_append: the coefficients have %u rows, expected fdim = %u.",
             ncm_matrix_nrows (coeffs), sp->fdim);

  ncm_spectral_panels_append_data (sp, a, b, ncm_matrix_data (coeffs), ncm_matrix_ncols (coeffs), ncm_matrix_tda (coeffs));
}

/**
 * ncm_spectral_panels_append_data: (skip)
 * @sp: a #NcmSpectralPanels
 * @a: left edge of the panel
 * @b: right edge of the panel
 * @c: the coefficients of the panel, #NcmSpectralPanels:fdim rows of @N
 * @N: number of coefficients in each row
 * @tda: distance between the starts of consecutive rows of @c
 *
 * As ncm_spectral_panels_append(), for coefficients stored in a buffer.
 */
void
ncm_spectral_panels_append_data (NcmSpectralPanels *sp, gdouble a, gdouble b, const gdouble *c, guint N, guint tda)
{
  GArray *coeffs;
  guint r;

  if (N == 0)
    g_error ("ncm_spectral_panels_append: the panel [%.17g, %.17g] has no coefficients.", a, b);

  if (!(a < b))
    g_error ("ncm_spectral_panels_append: the panel [%.17g, %.17g] is empty.", a, b);

  if (sp->edges->len == 0)
  {
    g_array_append_val (sp->edges, a);
  }
  else
  {
    const gdouble last = g_array_index (sp->edges, gdouble, sp->edges->len - 1);

    if (a != last)
      g_error ("ncm_spectral_panels_append: the panel starts at %.17g, not at the last edge %.17g.",
               a, last);
  }

  g_array_append_val (sp->edges, b);

  if (sp->n_panels == sp->pool->len)
    g_ptr_array_add (sp->pool, g_array_new (FALSE, FALSE, sizeof (gdouble)));

  coeffs = g_ptr_array_index (sp->pool, sp->n_panels);
  g_array_set_size (coeffs, (gsize) sp->fdim * N);

  for (r = 0; r < sp->fdim; r++)
    memcpy (&g_array_index (coeffs, gdouble, (gsize) r * N), c + (gsize) r * tda, sizeof (gdouble) * N);

  sp->n_panels++;
}

/**
 * ncm_spectral_panels_reset:
 * @sp: a #NcmSpectralPanels
 *
 * Removes every panel, keeping #NcmSpectralPanels:fdim and the arrays of the panels for
 * reuse, so that refilling @sp allocates nothing until it holds more panels, or larger
 * ones, than before.
 */
void
ncm_spectral_panels_reset (NcmSpectralPanels *sp)
{
  g_array_set_size (sp->edges, 0);
  sp->n_panels = 0;
}

/**
 * ncm_spectral_panels_get_fdim:
 * @sp: a #NcmSpectralPanels
 *
 * Returns: the dimension of the values of $F$.
 */
guint
ncm_spectral_panels_get_fdim (NcmSpectralPanels *sp)
{
  return sp->fdim;
}

/**
 * ncm_spectral_panels_get_n_panels:
 * @sp: a #NcmSpectralPanels
 *
 * Returns: the number of panels.
 */
guint
ncm_spectral_panels_get_n_panels (NcmSpectralPanels *sp)
{
  return sp->n_panels;
}

/**
 * ncm_spectral_panels_get_range:
 * @sp: a #NcmSpectralPanels
 * @a: (out): the first edge
 * @b: (out): the last edge
 *
 * Gets the interval covered by the panels; @sp must have at least one panel.
 */
void
ncm_spectral_panels_get_range (NcmSpectralPanels *sp, gdouble *a, gdouble *b)
{
  g_assert_cmpuint (sp->edges->len, >, 1);

  *a = g_array_index (sp->edges, gdouble, 0);
  *b = g_array_index (sp->edges, gdouble, sp->edges->len - 1);
}

/**
 * ncm_spectral_panels_peek_panel: (skip)
 * @sp: a #NcmSpectralPanels
 * @i: panel index
 * @a: (out): left edge of the panel
 * @b: (out): right edge of the panel
 * @N: (out): number of coefficients in each row
 *
 * Gets panel @i, which must be below ncm_spectral_panels_get_n_panels(). Its
 * #NcmSpectralPanels:fdim rows are contiguous, row $r$ starting $r N$ entries after the
 * returned pointer, which is valid until @sp is changed.
 *
 * Returns: (transfer none): the coefficients of panel @i.
 */
const gdouble *
ncm_spectral_panels_peek_panel (NcmSpectralPanels *sp, guint i, gdouble *a, gdouble *b, guint *N)
{
  GArray *coeffs;

  g_assert_cmpuint (i, <, sp->n_panels);

  coeffs = g_ptr_array_index (sp->pool, i);

  *a = g_array_index (sp->edges, gdouble, i);
  *b = g_array_index (sp->edges, gdouble, i + 1);
  *N = coeffs->len / sp->fdim;

  return (const gdouble *) coeffs->data;
}

/**
 * ncm_spectral_panels_get_panel:
 * @sp: a #NcmSpectralPanels
 * @i: panel index
 * @a: (out): left edge of the panel
 * @b: (out): right edge of the panel
 *
 * Gets a copy of panel @i, which must be below ncm_spectral_panels_get_n_panels().
 *
 * Returns: (transfer full): the #NcmSpectralPanels:fdim by $N$ coefficients of panel @i.
 */
NcmMatrix *
ncm_spectral_panels_get_panel (NcmSpectralPanels *sp, guint i, gdouble *a, gdouble *b)
{
  guint N;
  const gdouble *c  = ncm_spectral_panels_peek_panel (sp, i, a, b, &N);
  NcmMatrix *coeffs = ncm_matrix_new (sp->fdim, N);
  guint r;

  for (r = 0; r < sp->fdim; r++)
    memcpy (ncm_matrix_ptr (coeffs, r, 0), c + (gsize) r * N, sizeof (gdouble) * N);

  return coeffs;
}

/**
 * ncm_spectral_panels_find_panel:
 * @sp: a #NcmSpectralPanels
 * @x: a point
 *
 * Finds the panel holding @x by bisection over the edges; @sp must have at least one
 * panel. A point on an interior edge belongs to the panel that starts there, a point
 * below the first edge to the first panel, and a point on or above the last edge to the
 * last panel.
 *
 * Returns: the index of the panel.
 */
guint
ncm_spectral_panels_find_panel (NcmSpectralPanels *sp, gdouble x)
{
  const gdouble *edges = (const gdouble *) sp->edges->data;
  guint lo             = 0;
  guint hi;

  g_assert_cmpuint (sp->edges->len, >, 1);

  hi = sp->edges->len - 2;

  while (lo < hi)
  {
    const guint mid = (lo + hi + 1) / 2;

    if (x >= edges[mid])
      lo = mid;
    else
      hi = mid - 1;
  }

  return lo;
}

/**
 * ncm_spectral_panels_eval:
 * @sp: a #NcmSpectralPanels
 * @x: a point
 * @fval: the values of $F$ at @x
 *
 * Evaluates every entry of $F$ at @x into @fval, of length #NcmSpectralPanels:fdim, on
 * the panel given by ncm_spectral_panels_find_panel(). Outside the range of the panels the
 * series of the nearest end panel is extrapolated.
 */
void
ncm_spectral_panels_eval (NcmSpectralPanels *sp, gdouble x, NcmVector *fval)
{
  const guint i = ncm_spectral_panels_find_panel (sp, x);
  gdouble a, b;
  guint N, r;
  const gdouble *c = ncm_spectral_panels_peek_panel (sp, i, &a, &b, &N);
  const gdouble s  = ncm_spectral_x_to_s (a, b, x);

  for (r = 0; r < sp->fdim; r++)
    ncm_vector_set (fval, r, _ncm_spectral_panels_clenshaw (c + (gsize) r * N, N, s));
}

/**
 * ncm_spectral_panels_eval_row:
 * @sp: a #NcmSpectralPanels
 * @r: row index
 * @x: a point
 *
 * Evaluates the entry $F_r$ at @x, as ncm_spectral_panels_eval() does.
 *
 * Returns: $F_r(x)$.
 */
gdouble
ncm_spectral_panels_eval_row (NcmSpectralPanels *sp, guint r, gdouble x)
{
  return ncm_spectral_panels_eval_row_on (sp, ncm_spectral_panels_find_panel (sp, x), r, x);
}

/**
 * ncm_spectral_panels_eval_row_on:
 * @sp: a #NcmSpectralPanels
 * @i: panel index
 * @r: row index
 * @x: a point
 *
 * Evaluates the entry $F_r$ at @x on panel @i, whatever panel holds @x. On an edge of
 * panel @i this is the limit from inside it, where ncm_spectral_panels_eval_row() takes
 * the panel that starts at the edge.
 *
 * Returns: the series of row @r of panel @i at @x.
 */
gdouble
ncm_spectral_panels_eval_row_on (NcmSpectralPanels *sp, guint i, guint r, gdouble x)
{
  gdouble a, b;
  guint N;
  const gdouble *c = ncm_spectral_panels_peek_panel (sp, i, &a, &b, &N);

  return _ncm_spectral_panels_clenshaw (c + (gsize) r * N, N, ncm_spectral_x_to_s (a, b, x));
}

/* The series of N coefficients c at s in [-1, 1] by the Clenshaw recurrence */
static gdouble
_ncm_spectral_panels_clenshaw (const gdouble *c, guint N, const gdouble s)
{
  const gdouble two_s = 2.0 * s;
  gdouble b_1         = 0.0;
  gdouble b_2         = 0.0;
  gint n;

  for (n = (gint) N - 1; n >= 1; n--)
  {
    const gdouble b_0 = two_s * b_1 - b_2 + c[n];

    b_2 = b_1;
    b_1 = b_0;
  }

  return s * b_1 - b_2 + c[0];
}

/**
 * ncm_spectral_panels_restrict:
 * @sp: a #NcmSpectralPanels
 * @spectral: a #NcmSpectral
 * @a: left end of the subinterval
 * @b: right end of the subinterval
 * @coeffs: (out) (transfer full): the coefficients on [@a, @b]
 *
 * Computes the Chebyshev coefficients of every entry of $F$ on [@a, @b], which must lie
 * inside one panel, see ncm_spectral_panels_restrict_rows().
 *
 * Returns: whether [@a, @b] lies inside one panel; when it does not, @coeffs is unchanged.
 */
gboolean
ncm_spectral_panels_restrict (NcmSpectralPanels *sp, NcmSpectral *spectral, gdouble a, gdouble b, NcmMatrix **coeffs)
{
  return ncm_spectral_panels_restrict_rows (sp, spectral, a, b, 0, sp->fdim, coeffs);
}

/**
 * ncm_spectral_panels_restrict_rows:
 * @sp: a #NcmSpectralPanels
 * @spectral: a #NcmSpectral
 * @a: left end of the subinterval
 * @b: right end of the subinterval
 * @r0: the first row
 * @n: the number of rows
 * @coeffs: (out) (transfer full): the coefficients on [@a, @b], @n rows
 *
 * Computes the Chebyshev coefficients of the entries [@r0, @r0 + @n) of $F$ on [@a, @b],
 * which must lie inside one panel; row $i$ of @coeffs is entry @r0 + $i$. An end that
 * edges computed in another variable put a rounding outside the panel, within $10^{-12}$
 * of its width or four units in the last place of the larger edge, is clamped to the
 * edge. The restriction of a polynomial to a subinterval is a polynomial of the same
 * degree, so the coefficients follow exactly from those of the panel by
 * ncm_spectral_chebyshev_rebase_rows_data(), with the scratch space of @spectral; when
 * [@a, @b] is the panel they are its own. A matrix of the right shape in @coeffs is
 * reused, and a new one is allocated otherwise. @a must be below @b, and the rows inside
 * #NcmSpectralPanels:fdim.
 *
 * Returns: whether [@a, @b] lies inside one panel; when it does not, @coeffs is unchanged.
 */
gboolean
ncm_spectral_panels_restrict_rows (NcmSpectralPanels *sp, NcmSpectral *spectral, gdouble a, gdouble b, guint r0, guint n, NcmMatrix **coeffs)
{
  const guint i = ncm_spectral_panels_find_panel (sp, 0.5 * (a + b));
  gdouble pa, pb;
  guint N;
  const gdouble *c = ncm_spectral_panels_peek_panel (sp, i, &pa, &pb, &N);

  const gdouble slack = MAX (1.0e-12 * (pb - pa), 4.0 * DBL_EPSILON * MAX (fabs (pa), fabs (pb)));

  g_assert_cmpfloat (a, <, b);

  if (r0 + n > sp->fdim)
    g_error ("ncm_spectral_panels_restrict_rows: rows [%u, %u) are not inside the %u rows of the panel set.",
             r0, r0 + n, sp->fdim);

  if ((a < pa - slack) || (b > pb + slack))
    return FALSE;

  a = MAX (a, pa);
  b = MIN (b, pb);
  c = c + (gsize) r0 * N;

  if ((*coeffs != NULL) && ((ncm_matrix_nrows (*coeffs) != n) || (ncm_matrix_ncols (*coeffs) != N)))
    ncm_matrix_clear (coeffs);

  if (*coeffs == NULL)
    *coeffs = ncm_matrix_new (n, N);

  if ((a == pa) && (b == pb))
  {
    guint r;

    for (r = 0; r < n; r++)
      memcpy (ncm_matrix_ptr (*coeffs, r, 0), c + (gsize) r * N, sizeof (gdouble) * N);
  }
  else
  {
    ncm_spectral_chebyshev_rebase_rows_data (spectral, c, n, N, N, pa, pb, a, b, *coeffs);
  }

  return TRUE;
}

static gint
_ncm_spectral_panels_cmp_gdouble (gconstpointer a, gconstpointer b)
{
  const gdouble x = *(const gdouble *) a;
  const gdouble y = *(const gdouble *) b;

  return (x > y) - (x < y);
}

/**
 * ncm_spectral_panels_new_sum_rows: (skip)
 * @spectral: a #NcmSpectral, for the scratch of the rebase
 * @n: the number of rows
 * @sps: the panel set of each row
 * @rows: the row of each in its panel set
 * @scales: the scale $s_r$ of each row: row $r$ is a function of $x = k / s_r$
 * @x_min: the lower end of each row's support in its $x$, the row being zero outside
 * @x_max: the upper end of each row's support in its $x$
 * @abstol: trailing coefficients of a cell's sum below this are dropped, keeping at least
 *   one; 0 keeps all
 *
 * Builds the sum of @n scaled rows as one row in $k$. The partition is the union of the
 * rows' support ends and of their panel edges inside the supports, all in $k$; on each
 * cell every row whose support contains it is rebased onto the cell with
 * ncm_spectral_chebyshev_rebase_rows_data(), which is exact for a polynomial, and the
 * coefficients are added. A cell no row covers holds zero. A row with an empty support
 * ($x_\mathrm{min} \ge x_\mathrm{max}$) is skipped. Edges closer than 16 ulps are merged,
 * which moves a row's edge by that much. Dropping the coefficients below @abstol changes
 * a cell's sum by at most their number times @abstol.
 *
 * Returns: (transfer full): a one-row #NcmSpectralPanels over the union of the supports,
 *   without panels when every support is empty
 */
NcmSpectralPanels *
ncm_spectral_panels_new_sum_rows (NcmSpectral *spectral, guint n, NcmSpectralPanels * const *sps, const guint *rows, const gdouble *scales, const gdouble *x_min, const gdouble *x_max, gdouble abstol)
{
  NcmSpectralPanels *out = ncm_spectral_panels_new (1);
  GArray *edges          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *sum            = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GPtrArray *rebased     = g_ptr_array_new (); /* the rebase output of each N met, indexed by N */
  guint r, e, j;

  for (r = 0; r < n; r++)
  {
    NcmSpectralPanels *sp = sps[r];
    const gdouble s       = scales[r];
    const gdouble lo      = s * x_min[r];
    const gdouble hi      = s * x_max[r];
    guint p;

    g_assert_cmpuint (rows[r], <, sp->fdim);
    g_assert_cmpfloat (s, >, 0.0);

    if (!(x_min[r] < x_max[r]))
      continue;

    g_array_append_val (edges, lo);
    g_array_append_val (edges, hi);

    for (p = 1; p + 1 < sp->edges->len; p++)
    {
      const gdouble k = s * g_array_index (sp->edges, gdouble, p);

      if ((k > lo) && (k < hi))
        g_array_append_val (edges, k);
    }
  }

  if (edges->len > 0)
  {
    guint kept = 0;

    g_array_sort (edges, &_ncm_spectral_panels_cmp_gdouble);

    for (e = 0; e < edges->len; e++)
    {
      const gdouble k = g_array_index (edges, gdouble, e);

      if ((kept == 0) || (k - g_array_index (edges, gdouble, kept - 1) > 16.0 * DBL_EPSILON * fabs (k)))
        g_array_index (edges, gdouble, kept++) = k;
    }

    g_array_set_size (edges, kept);
  }

  for (e = 1; e < edges->len; e++)
  {
    const gdouble a   = g_array_index (edges, gdouble, e - 1);
    const gdouble b   = g_array_index (edges, gdouble, e);
    const gdouble mid = 0.5 * (a + b);

    g_array_set_size (sum, 0);

    for (r = 0; r < n; r++)
    {
      NcmSpectralPanels *sp = sps[r];
      const gdouble s       = scales[r];
      const gdouble x       = mid / s;
      const gdouble *c;
      NcmMatrix *m;
      gdouble pa, pb;
      guint p, N;

      if (!(x_min[r] < x_max[r]) || (x <= x_min[r]) || (x >= x_max[r]))
        continue;

      p = ncm_spectral_panels_find_panel (sp, x);
      c = ncm_spectral_panels_peek_panel (sp, p, &pa, &pb, &N) + (gsize) rows[r] * N;

      if (rebased->len <= N)
        g_ptr_array_set_size (rebased, N + 1);

      if (g_ptr_array_index (rebased, N) == NULL)
        g_ptr_array_index (rebased, N) = ncm_matrix_new (1, N);

      m = g_ptr_array_index (rebased, N);
      ncm_spectral_chebyshev_rebase_rows_data (spectral, c, 1, N, N, pa, pb, a / s, b / s, m);

      if (sum->len < N)
      {
        const guint old = sum->len;

        g_array_set_size (sum, N);

        for (j = old; j < N; j++)
          g_array_index (sum, gdouble, j) = 0.0;
      }

      for (j = 0; j < N; j++)
        g_array_index (sum, gdouble, j) += ncm_matrix_get (m, 0, j);
    }

    if (sum->len == 0)
    {
      const gdouble zero = 0.0;

      ncm_spectral_panels_append_data (out, a, b, &zero, 1, 1);
      continue;
    }

    if (abstol > 0.0)
    {
      guint keep = sum->len;

      while ((keep > 1) && (fabs (g_array_index (sum, gdouble, keep - 1)) < abstol))
        keep--;

      g_array_set_size (sum, keep);
    }

    ncm_spectral_panels_append_data (out, a, b, (const gdouble *) sum->data, sum->len, sum->len);
  }

  for (j = 0; j < rebased->len; j++)
    if (g_ptr_array_index (rebased, j) != NULL)
      ncm_matrix_free (g_ptr_array_index (rebased, j));

  g_array_unref (edges);
  g_array_unref (sum);
  g_ptr_array_unref (rebased);

  return out;
}

