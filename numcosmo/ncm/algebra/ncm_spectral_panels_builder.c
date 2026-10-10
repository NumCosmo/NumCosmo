/***************************************************************************
 *            ncm_spectral_panels_builder.c
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
 * NcmSpectralPanelsBuilder:
 *
 * Adaptive construction of the Chebyshev panels of a vector function.
 *
 * ncm_spectral_panels_builder_build() expands $F$ on $[a, b]$ into an #NcmSpectralPanels.
 * The breakpoints, the points where $F$ or a derivative may jump, split $[a, b]$ into
 * continuous intervals, and no panel straddles one. With
 * #NcmSpectralPanelsBuilder:min-panels-per-efold nonzero each interval is first split at
 * the geometric grid $x = e^{j / p}$, $j$ integer and $p$ that property, anchored at 1 so
 * that panel sets of the same variable share their grid edges. Each initial panel is then
 * expanded on nested Chebyshev-Lobatto grids from level #NcmSpectralPanelsBuilder:level-min
 * up to #NcmSpectralPanelsBuilder:level-cap, as
 * ncm_spectral_compute_chebyshev_coeffs_batch_adaptive_cap() does, and bisected when it
 * does not converge; the halves reuse the values of $F$ at their ends.
 *
 * A build takes its scratch space from a pool, so one builder may serve builds in
 * several threads; its properties must not change during a build.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/algebra/ncm_spectral_panels_builder.h"
#include "ncm/algebra/ncm_spectral.h"
#include "ncm/core/ncm_memory_pool.h"

#include <math.h>

enum
{
  PROP_0,
  PROP_LEVEL_MIN,
  PROP_LEVEL_CAP,
  PROP_MIN_PANELS_PER_EFOLD,
  PROP_MIN_PANEL_FRACTION,
};

struct _NcmSpectralPanelsBuilder
{
  /*< private >*/
  GObject parent_instance;
  guint level_min;
  guint level_cap;
  gdouble min_panels_per_efold;
  gdouble min_panel_fraction;
  NcmMemoryPool *work_pool; /* NcmSpectralPanelsBuilderWork, one per concurrent build */
};

/* The scratch space of one build */
typedef struct _NcmSpectralPanelsBuilderWork
{
  NcmSpectral *spectral;
  guint fdim;
  GPtrArray *ends;      /* NcmMatrix, per split depth the 3 by fdim F at a, the midpoint and b */
  GPtrArray *ends_rows; /* NcmVector, the three row views of each ends matrix */
  GArray *edges;        /* NcmSpectralPanelsBuilderEdge, the initial edges */
  NcmVector *f_a;       /* the known values at the ends of an initial panel */
  NcmVector *f_b;
} NcmSpectralPanelsBuilderWork;

/* An initial edge, and the index of the breakpoint it is, or -1 */
typedef struct _NcmSpectralPanelsBuilderEdge
{
  gdouble x;
  gint bp;
} NcmSpectralPanelsBuilderEdge;

/* What the recursive split needs, for one build */
typedef struct _NcmSpectralPanelsBuilderCall
{
  NcmSpectralPanelsF F;
  NcmSpectralPanelsAccept accept;
  gpointer user_data;
  guint interval;
  guint level_min;
  guint level_cap;
  gdouble min_panel_fraction;
  gdouble reltol;
  gdouble abstol;
  NcmSpectralPanelsBuilderWork *work;
  NcmSpectralPanels *sp;
} NcmSpectralPanelsBuilderCall;

G_DEFINE_TYPE (NcmSpectralPanelsBuilder, ncm_spectral_panels_builder, G_TYPE_OBJECT)

static gpointer _ncm_spectral_panels_builder_work_new (gpointer userdata);
static void _ncm_spectral_panels_builder_work_free (gpointer p);

static void
ncm_spectral_panels_builder_init (NcmSpectralPanelsBuilder *builder)
{
  builder->level_min            = 0;
  builder->level_cap            = 0;
  builder->min_panels_per_efold = 0.0;
  builder->min_panel_fraction   = 0.0;
  builder->work_pool            = ncm_memory_pool_new (_ncm_spectral_panels_builder_work_new, NULL,
                                                       _ncm_spectral_panels_builder_work_free);
}

static void
_ncm_spectral_panels_builder_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSpectralPanelsBuilder *builder = NCM_SPECTRAL_PANELS_BUILDER (object);

  g_return_if_fail (NCM_IS_SPECTRAL_PANELS_BUILDER (object));

  switch (prop_id)
  {
    case PROP_LEVEL_MIN:
      ncm_spectral_panels_builder_set_level_min (builder, g_value_get_uint (value));
      break;
    case PROP_LEVEL_CAP:
      ncm_spectral_panels_builder_set_level_cap (builder, g_value_get_uint (value));
      break;
    case PROP_MIN_PANELS_PER_EFOLD:
      ncm_spectral_panels_builder_set_min_panels_per_efold (builder, g_value_get_double (value));
      break;
    case PROP_MIN_PANEL_FRACTION:
      ncm_spectral_panels_builder_set_min_panel_fraction (builder, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_spectral_panels_builder_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSpectralPanelsBuilder *builder = NCM_SPECTRAL_PANELS_BUILDER (object);

  g_return_if_fail (NCM_IS_SPECTRAL_PANELS_BUILDER (object));

  switch (prop_id)
  {
    case PROP_LEVEL_MIN:
      g_value_set_uint (value, builder->level_min);
      break;
    case PROP_LEVEL_CAP:
      g_value_set_uint (value, builder->level_cap);
      break;
    case PROP_MIN_PANELS_PER_EFOLD:
      g_value_set_double (value, builder->min_panels_per_efold);
      break;
    case PROP_MIN_PANEL_FRACTION:
      g_value_set_double (value, builder->min_panel_fraction);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_spectral_panels_builder_finalize (GObject *object)
{
  NcmSpectralPanelsBuilder *builder = NCM_SPECTRAL_PANELS_BUILDER (object);

  ncm_memory_pool_free (builder->work_pool, TRUE);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spectral_panels_builder_parent_class)->finalize (object);
}

static void
ncm_spectral_panels_builder_class_init (NcmSpectralPanelsBuilderClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_spectral_panels_builder_set_property;
  object_class->get_property = &_ncm_spectral_panels_builder_get_property;
  object_class->finalize     = &_ncm_spectral_panels_builder_finalize;

  /**
   * NcmSpectralPanelsBuilder:level-min:
   *
   * The level at which the expansion of a panel starts, below
   * #NcmSpectralPanelsBuilder:level-cap.
   */
  g_object_class_install_property (object_class,
                                   PROP_LEVEL_MIN,
                                   g_param_spec_uint ("level-min",
                                                      NULL,
                                                      "Starting level of a panel",
                                                      0, 29, 3,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSpectralPanelsBuilder:level-cap:
   *
   * The highest level tried on a panel before it is bisected.
   */
  g_object_class_install_property (object_class,
                                   PROP_LEVEL_CAP,
                                   g_param_spec_uint ("level-cap",
                                                      NULL,
                                                      "Highest level of a panel",
                                                      1, 30, 5,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSpectralPanelsBuilder:min-panels-per-efold:
   *
   * The minimum number of panels per e-fold between breakpoints, 0 to disable. It avoids
   * premature acceptance, a first grid too coarse to see a feature, and chains of failed
   * splits on hard intervals. A nonzero value requires $a > 0$.
   */
  g_object_class_install_property (object_class,
                                   PROP_MIN_PANELS_PER_EFOLD,
                                   g_param_spec_double ("min-panels-per-efold",
                                                        NULL,
                                                        "Minimum panels per e-fold",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSpectralPanelsBuilder:min-panel-fraction:
   *
   * The width, as a fraction of $\max(|a|, |b|)$ of the panel $[a, b]$, below which a
   * panel is no longer bisected: its expansion up to #NcmSpectralPanelsBuilder:level-cap
   * must converge, or the build is an error.
   */
  g_object_class_install_property (object_class,
                                   PROP_MIN_PANEL_FRACTION,
                                   g_param_spec_double ("min-panel-fraction",
                                                        NULL,
                                                        "Minimum panel width fraction",
                                                        G_MINDOUBLE, 1.0, 1.0e-6,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

static gpointer
_ncm_spectral_panels_builder_work_new (gpointer userdata)
{
  NcmSpectralPanelsBuilderWork *work = g_new0 (NcmSpectralPanelsBuilderWork, 1);

  work->spectral  = ncm_spectral_new ();
  work->fdim      = 0;
  work->ends      = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_matrix_free);
  work->ends_rows = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_vector_free);
  work->edges     = g_array_new (FALSE, FALSE, sizeof (NcmSpectralPanelsBuilderEdge));
  work->f_a       = NULL;
  work->f_b       = NULL;

  return work;
}

static void
_ncm_spectral_panels_builder_work_free (gpointer p)
{
  NcmSpectralPanelsBuilderWork *work = p;

  ncm_spectral_free (work->spectral);
  g_ptr_array_unref (work->ends_rows);
  g_ptr_array_unref (work->ends);
  g_array_unref (work->edges);
  ncm_vector_clear (&work->f_a);
  ncm_vector_clear (&work->f_b);

  g_free (work);
}

/**
 * ncm_spectral_panels_builder_new:
 *
 * Creates a builder with the default properties.
 *
 * Returns: (transfer full): a new #NcmSpectralPanelsBuilder.
 */
NcmSpectralPanelsBuilder *
ncm_spectral_panels_builder_new (void)
{
  return g_object_new (NCM_TYPE_SPECTRAL_PANELS_BUILDER, NULL);
}

/**
 * ncm_spectral_panels_builder_ref:
 * @builder: a #NcmSpectralPanelsBuilder
 *
 * Increases the reference count of @builder by one.
 *
 * Returns: (transfer full): @builder.
 */
NcmSpectralPanelsBuilder *
ncm_spectral_panels_builder_ref (NcmSpectralPanelsBuilder *builder)
{
  return g_object_ref (builder);
}

/**
 * ncm_spectral_panels_builder_free:
 * @builder: a #NcmSpectralPanelsBuilder
 *
 * Decreases the reference count of @builder by one.
 */
void
ncm_spectral_panels_builder_free (NcmSpectralPanelsBuilder *builder)
{
  g_object_unref (builder);
}

/**
 * ncm_spectral_panels_builder_clear:
 * @builder: a #NcmSpectralPanelsBuilder
 *
 * Decreases the reference count of *@builder by one and sets *@builder to %NULL.
 */
void
ncm_spectral_panels_builder_clear (NcmSpectralPanelsBuilder **builder)
{
  g_clear_object (builder);
}

/**
 * ncm_spectral_panels_builder_set_level_min:
 * @builder: a #NcmSpectralPanelsBuilder
 * @level_min: starting level of a panel
 *
 * Sets #NcmSpectralPanelsBuilder:level-min.
 */
void
ncm_spectral_panels_builder_set_level_min (NcmSpectralPanelsBuilder *builder, guint level_min)
{
  builder->level_min = level_min;
}

/**
 * ncm_spectral_panels_builder_set_level_cap:
 * @builder: a #NcmSpectralPanelsBuilder
 * @level_cap: highest level of a panel
 *
 * Sets #NcmSpectralPanelsBuilder:level-cap.
 */
void
ncm_spectral_panels_builder_set_level_cap (NcmSpectralPanelsBuilder *builder, guint level_cap)
{
  builder->level_cap = level_cap;
}

/**
 * ncm_spectral_panels_builder_set_min_panels_per_efold:
 * @builder: a #NcmSpectralPanelsBuilder
 * @min_panels_per_efold: minimum panels per e-fold
 *
 * Sets #NcmSpectralPanelsBuilder:min-panels-per-efold.
 */
void
ncm_spectral_panels_builder_set_min_panels_per_efold (NcmSpectralPanelsBuilder *builder, gdouble min_panels_per_efold)
{
  builder->min_panels_per_efold = min_panels_per_efold;
}

/**
 * ncm_spectral_panels_builder_set_min_panel_fraction:
 * @builder: a #NcmSpectralPanelsBuilder
 * @min_panel_fraction: minimum panel width fraction
 *
 * Sets #NcmSpectralPanelsBuilder:min-panel-fraction.
 */
void
ncm_spectral_panels_builder_set_min_panel_fraction (NcmSpectralPanelsBuilder *builder, gdouble min_panel_fraction)
{
  builder->min_panel_fraction = min_panel_fraction;
}

/**
 * ncm_spectral_panels_builder_get_level_min:
 * @builder: a #NcmSpectralPanelsBuilder
 *
 * Returns: #NcmSpectralPanelsBuilder:level-min.
 */
guint
ncm_spectral_panels_builder_get_level_min (NcmSpectralPanelsBuilder *builder)
{
  return builder->level_min;
}

/**
 * ncm_spectral_panels_builder_get_level_cap:
 * @builder: a #NcmSpectralPanelsBuilder
 *
 * Returns: #NcmSpectralPanelsBuilder:level-cap.
 */
guint
ncm_spectral_panels_builder_get_level_cap (NcmSpectralPanelsBuilder *builder)
{
  return builder->level_cap;
}

/**
 * ncm_spectral_panels_builder_get_min_panels_per_efold:
 * @builder: a #NcmSpectralPanelsBuilder
 *
 * Returns: #NcmSpectralPanelsBuilder:min-panels-per-efold.
 */
gdouble
ncm_spectral_panels_builder_get_min_panels_per_efold (NcmSpectralPanelsBuilder *builder)
{
  return builder->min_panels_per_efold;
}

/**
 * ncm_spectral_panels_builder_get_min_panel_fraction:
 * @builder: a #NcmSpectralPanelsBuilder
 *
 * Returns: #NcmSpectralPanelsBuilder:min-panel-fraction.
 */
gdouble
ncm_spectral_panels_builder_get_min_panel_fraction (NcmSpectralPanelsBuilder *builder)
{
  return builder->min_panel_fraction;
}

/**
 * ncm_spectral_panels_builder_build:
 * @builder: a #NcmSpectralPanelsBuilder
 * @F: (scope call) (closure user_data): function to expand
 * @user_data: user data for @F
 * @a: left end of the range
 * @b: right end of the range
 * @breakpoints: (nullable): points where $F$ may jump, ascending and inside $(a, b)$
 * @f_minus: (nullable): $F$ at each breakpoint from the left, one row per breakpoint
 * @f_plus: (nullable): $F$ at each breakpoint from the right, one row per breakpoint
 * @reltol: relative tolerance of a panel
 * @abstol: absolute tolerance of a panel, in units of $F$
 * @sp: the panel set to fill
 *
 * Resets @sp and fills it with the panels of $F$ on [@a, @b], with
 * #NcmSpectralPanels:fdim of @sp the dimension of the values of $F$. A panel is accepted
 * as ncm_spectral_compute_chebyshev_coeffs_batch_adaptive_cap() accepts an expansion,
 * with @reltol and @abstol. @f_minus and @f_plus are given together or not at all, with
 * #NcmSpectralPanels:fdim columns; a known value is used instead of evaluating $F$ at the
 * end of a panel.
 */
void
ncm_spectral_panels_builder_build (NcmSpectralPanelsBuilder *builder, NcmSpectralPanelsF F, gpointer user_data,
                                   gdouble a, gdouble b, NcmVector *breakpoints, NcmMatrix *f_minus, NcmMatrix *f_plus,
                                   gdouble reltol, gdouble abstol, NcmSpectralPanels *sp)
{
  ncm_spectral_panels_builder_build_full (builder, F, NULL, user_data, a, b, breakpoints, f_minus, f_plus,
                                          reltol, abstol, sp);
}

static void _ncm_spectral_panels_builder_check (NcmSpectralPanelsBuilder *builder, guint fdim, gdouble a, gdouble b, NcmVector *breakpoints, NcmMatrix *f_minus, NcmMatrix *f_plus);
static void _ncm_spectral_panels_builder_work_prepare (NcmSpectralPanelsBuilderWork *work, guint fdim, guint level_cap);
static void _ncm_spectral_panels_builder_initial_edges (NcmSpectralPanelsBuilder *builder, NcmSpectralPanelsBuilderWork *work, gdouble a, gdouble b, NcmVector *breakpoints);
static void _ncm_spectral_panels_builder_copy_row (NcmMatrix *m, guint row, NcmVector *v);
static void _ncm_spectral_panels_builder_split (NcmSpectralPanelsBuilderCall *call, gdouble a, gdouble b, guint depth, NcmVector *f_a, NcmVector *f_b);

/**
 * ncm_spectral_panels_builder_build_full: (skip)
 * @builder: a #NcmSpectralPanelsBuilder
 * @F: function to expand
 * @accept: (nullable): further acceptance test of a converged panel
 * @user_data: user data for @F and @accept
 * @a: left end of the range
 * @b: right end of the range
 * @breakpoints: (nullable): points where $F$ may jump, ascending and inside $(a, b)$
 * @f_minus: (nullable): $F$ at each breakpoint from the left, one row per breakpoint
 * @f_plus: (nullable): $F$ at each breakpoint from the right, one row per breakpoint
 * @reltol: relative tolerance of a panel
 * @abstol: absolute tolerance of a panel, in units of $F$
 * @sp: the panel set to fill
 *
 * As ncm_spectral_panels_builder_build(), with a converged panel kept only when @accept
 * returns %TRUE and bisected otherwise. A panel narrower than
 * #NcmSpectralPanelsBuilder:min-panel-fraction is kept without @accept.
 */
void
ncm_spectral_panels_builder_build_full (NcmSpectralPanelsBuilder *builder, NcmSpectralPanelsF F, NcmSpectralPanelsAccept accept,
                                        gpointer user_data, gdouble a, gdouble b, NcmVector *breakpoints,
                                        NcmMatrix *f_minus, NcmMatrix *f_plus, gdouble reltol, gdouble abstol,
                                        NcmSpectralPanels *sp)
{
  const guint fdim                   = ncm_spectral_panels_get_fdim (sp);
  NcmSpectralPanelsBuilderWork **ptr = ncm_memory_pool_get (builder->work_pool);
  NcmSpectralPanelsBuilderWork *work = *ptr;
  NcmSpectralPanelsBuilderCall call  = {
    .F                  = F,
    .accept             = accept,
    .user_data          = user_data,
    .interval           = 0,
    .level_min          = builder->level_min,
    .level_cap          = builder->level_cap,
    .min_panel_fraction = builder->min_panel_fraction,
    .reltol             = reltol,
    .abstol             = abstol,
    .work               = work,
    .sp                 = sp
  };
  guint k;

  _ncm_spectral_panels_builder_check (builder, fdim, a, b, breakpoints, f_minus, f_plus);
  _ncm_spectral_panels_builder_work_prepare (work, fdim, builder->level_cap);
  _ncm_spectral_panels_builder_initial_edges (builder, work, a, b, breakpoints);

  ncm_spectral_panels_reset (sp);

  /* An initial edge that is a breakpoint has its values on each side only when they are
   * given; a grid edge has none, the first level evaluates it. */
  for (k = 0; k + 1 < work->edges->len; k++)
  {
    const NcmSpectralPanelsBuilderEdge *e0 = &g_array_index (work->edges, NcmSpectralPanelsBuilderEdge, k);
    const NcmSpectralPanelsBuilderEdge *e1 = &g_array_index (work->edges, NcmSpectralPanelsBuilderEdge, k + 1);
    NcmVector *f_a                         = NULL;
    NcmVector *f_b                         = NULL;

    if (e0->bp >= 0)
    {
      call.interval = e0->bp + 1;

      if (f_plus != NULL)
      {
        _ncm_spectral_panels_builder_copy_row (f_plus, e0->bp, work->f_a);
        f_a = work->f_a;
      }
    }

    if ((e1->bp >= 0) && (f_minus != NULL))
    {
      _ncm_spectral_panels_builder_copy_row (f_minus, e1->bp, work->f_b);
      f_b = work->f_b;
    }

    _ncm_spectral_panels_builder_split (&call, e0->x, e1->x, 0, f_a, f_b);
  }

  ncm_memory_pool_return (ptr);
}

/* The arguments of a build, checked before anything is evaluated */
static void
_ncm_spectral_panels_builder_check (NcmSpectralPanelsBuilder *builder, guint fdim, gdouble a, gdouble b,
                                    NcmVector *breakpoints, NcmMatrix *f_minus, NcmMatrix *f_plus)
{
  const guint n_bp = (breakpoints != NULL) ? ncm_vector_len (breakpoints) : 0;
  guint j;

  if (builder->level_min >= builder->level_cap)
    g_error ("ncm_spectral_panels_builder_build: level-min %u is not below level-cap %u.",
             builder->level_min, builder->level_cap);

  if (!(a < b))
    g_error ("ncm_spectral_panels_builder_build: the range [%.17g, %.17g] is empty.", a, b);

  if ((builder->min_panels_per_efold > 0.0) && !(a > 0.0))
    g_error ("ncm_spectral_panels_builder_build: min-panels-per-efold %g needs a > 0, got a = %.17g.",
             builder->min_panels_per_efold, a);

  for (j = 0; j < n_bp; j++)
  {
    const gdouble x    = ncm_vector_get (breakpoints, j);
    const gdouble left = (j == 0) ? a : ncm_vector_get (breakpoints, j - 1);

    if (!((x > left) && (x < b)))
      g_error ("ncm_spectral_panels_builder_build: breakpoint %u at %.17g is not ascending inside (%.17g, %.17g).",
               j, x, a, b);
  }

  if ((f_minus == NULL) != (f_plus == NULL))
    g_error ("ncm_spectral_panels_builder_build: f_minus and f_plus must be given together.");

  if (f_minus != NULL)
    if ((ncm_matrix_nrows (f_minus) != n_bp) || (ncm_matrix_ncols (f_minus) != fdim) ||
        (ncm_matrix_nrows (f_plus) != n_bp) || (ncm_matrix_ncols (f_plus) != fdim))
      g_error ("ncm_spectral_panels_builder_build: f_minus and f_plus must be %u by %u.", n_bp, fdim);
}

/* Scratch of the dimension of this build; the ends of earlier depths are kept when it is the same */
static void
_ncm_spectral_panels_builder_work_prepare (NcmSpectralPanelsBuilderWork *work, guint fdim, guint level_cap)
{
  if (work->fdim != fdim)
  {
    g_ptr_array_set_size (work->ends_rows, 0);
    g_ptr_array_set_size (work->ends, 0);
    ncm_vector_clear (&work->f_a);
    ncm_vector_clear (&work->f_b);

    work->f_a  = ncm_vector_new (fdim);
    work->f_b  = ncm_vector_new (fdim);
    work->fdim = fdim;
  }

  if (ncm_spectral_get_max_level (work->spectral) < level_cap)
    ncm_spectral_set_max_level (work->spectral, level_cap);
}

static gint
_ncm_spectral_panels_builder_cmp_edge (gconstpointer pa, gconstpointer pb)
{
  const gdouble a = ((const NcmSpectralPanelsBuilderEdge *) pa)->x;
  const gdouble b = ((const NcmSpectralPanelsBuilderEdge *) pb)->x;

  return (a > b) - (a < b);
}

/*
 * The initial edges, ascending from a to b: the breakpoints, and with min-panels-per-efold
 * p nonzero the grid exp (j / p) inside the range. A grid point closer than a quarter of
 * the grid spacing, in ln x, to a or b or a breakpoint is dropped.
 */
static void
_ncm_spectral_panels_builder_initial_edges (NcmSpectralPanelsBuilder *builder, NcmSpectralPanelsBuilderWork *work,
                                            gdouble a, gdouble b, NcmVector *breakpoints)
{
  const guint n_bp                  = (breakpoints != NULL) ? ncm_vector_len (breakpoints) : 0;
  NcmSpectralPanelsBuilderEdge edge = {a, -1};
  guint j;

  g_array_set_size (work->edges, 0);
  g_array_append_val (work->edges, edge);

  for (j = 0; j < n_bp; j++)
  {
    edge.x  = ncm_vector_get (breakpoints, j);
    edge.bp = j;
    g_array_append_val (work->edges, edge);
  }

  edge.x  = b;
  edge.bp = -1;
  g_array_append_val (work->edges, edge);

  if (builder->min_panels_per_efold > 0.0)
  {
    const gdouble h       = 1.0 / builder->min_panels_per_efold;
    const glong j_min     = (glong) floor (log (a) / h) + 1;
    const glong j_max     = (glong) ceil (log (b) / h) - 1;
    const guint n_initial = work->edges->len;
    glong jj;
    guint i;

    for (jj = j_min; jj <= j_max; jj++)
    {
      const gdouble x = exp (jj * h);
      gboolean near   = FALSE;

      for (i = 0; i < n_initial; i++)
      {
        if (fabs (log (x / g_array_index (work->edges, NcmSpectralPanelsBuilderEdge, i).x)) < 0.25 * h)
        {
          near = TRUE;
          break;
        }
      }

      if (!near)
      {
        edge.x  = x;
        edge.bp = -1;
        g_array_append_val (work->edges, edge);
      }
    }

    g_array_sort (work->edges, _ncm_spectral_panels_builder_cmp_edge);
  }
}

static void
_ncm_spectral_panels_builder_copy_row (NcmMatrix *m, guint row, NcmVector *v)
{
  guint r;

  for (r = 0; r < ncm_vector_len (v); r++)
    ncm_vector_set (v, r, ncm_matrix_get (m, row, r));
}

static void _ncm_spectral_panels_builder_F (gpointer user_data, gdouble x, NcmVector *fval);

/*
 * Expands on [a, b] and appends the panel to the panel set, bisecting where the expansion
 * up to level-cap does not converge or accept refuses it. A split hands rows 0 and 1 of
 * the ends of its depth to the lower half and rows 1 and 2 to the upper half; the lower
 * half's own splits write one depth further down, so the rows survive for the upper half.
 */
static void
_ncm_spectral_panels_builder_split (NcmSpectralPanelsBuilderCall *call, gdouble a, gdouble b, guint depth,
                                    NcmVector *f_a, NcmVector *f_b)
{
  NcmSpectralPanelsBuilderWork *work = call->work;
  NcmMatrix *ends;
  const gdouble *c;
  guint level, tda;

  if (work->ends->len <= depth)
  {
    NcmMatrix *m = ncm_matrix_new (3, work->fdim);
    guint r;

    g_ptr_array_add (work->ends, m);

    for (r = 0; r < 3; r++)
      g_ptr_array_add (work->ends_rows, ncm_matrix_get_row (m, r));
  }

  ends = g_ptr_array_index (work->ends, depth);

  level = ncm_spectral_compute_chebyshev_coeffs_batch_adaptive_cap (work->spectral, _ncm_spectral_panels_builder_F, work->fdim,
                                                                    a, b, call->level_min, call->level_cap,
                                                                    call->reltol, call->abstol, FALSE, f_a, f_b, ends,
                                                                    NULL, call);

  if (level > 0)
  {
    const guint N = (1 << level) + 1;

    c = ncm_spectral_peek_batch_coeffs (work->spectral, &tda);

    if ((call->accept == NULL) || call->accept (call->user_data, a, b, c, N, tda))
    {
      ncm_spectral_panels_append_data (call->sp, a, b, c, N, tda);

      return;
    }
  }

  /* A panel that does not converge however far it is split is not a resolution
   * problem; refusing to bisect past a fraction of its scale turns a hang into an
   * error that names the tolerances. */
  if ((b - a) < call->min_panel_fraction * MAX (fabs (a), fabs (b)))
  {
    level = ncm_spectral_compute_chebyshev_coeffs_batch_adaptive_cap (work->spectral, _ncm_spectral_panels_builder_F, work->fdim,
                                                                      a, b, call->level_min, call->level_cap,
                                                                      call->reltol, call->abstol, TRUE, f_a, f_b, NULL,
                                                                      NULL, call);
    c = ncm_spectral_peek_batch_coeffs (work->spectral, &tda);
    ncm_spectral_panels_append_data (call->sp, a, b, c, (1 << level) + 1, tda);

    return;
  }

  {
    const gdouble mid = 0.5 * (a + b);
    NcmVector *f_lo   = g_ptr_array_index (work->ends_rows, 3 * depth + 0);
    NcmVector *f_mid  = g_ptr_array_index (work->ends_rows, 3 * depth + 1);
    NcmVector *f_hi   = g_ptr_array_index (work->ends_rows, 3 * depth + 2);

    _ncm_spectral_panels_builder_split (call, a, mid, depth + 1, f_lo, f_mid);
    _ncm_spectral_panels_builder_split (call, mid, b, depth + 1, f_mid, f_hi);
  }
}

/* F of the build on the continuous interval of the panel being expanded */
static void
_ncm_spectral_panels_builder_F (gpointer user_data, gdouble x, NcmVector *fval)
{
  NcmSpectralPanelsBuilderCall *call = user_data;

  call->F (call->user_data, call->interval, x, fval);
}

