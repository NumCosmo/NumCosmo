/***************************************************************************
 *            ncm_spectral_tensor_builder.c
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
 * NcmSpectralTensorBuilder:
 *
 * Adaptive Chebyshev expansion of a function into an #NcmSpectralTensor.
 *
 * A table of a function agrees with it only where its panels and numbers of coefficients
 * follow the variation of the function, which is not known in advance. An
 * #NcmSpectralTensorBuilder computes them: it fills the empty nodes of the table of a
 * function $F (x_0, \dots, x_{n - 1})$, given by a callback, refining each node until an
 * estimate of the error of its series is within a tolerance. The domain of the table must
 * be set before the construction, and nodes may be split beforehand, for example at a kink
 * of $F$.
 *
 * The variables of the table are those the callback receives: the builder applies no
 * change of variable. A function that varies on a logarithmic scale in one of its
 * arguments is tabulated in the logarithm of that argument, and the callback converts the
 * points it receives. For example, a function of a width $\sigma \in [0.02, 1.5]$ that
 * varies most at small $\sigma$ is tabulated in $x = \ln \sigma$, and the callback
 * computes $\sigma = e^x$. The caller evaluates the table at the same variables.
 *
 * The function ncm_spectral_tensor_builder_build() expands $F$ on each empty node of the
 * table, in the order of the node numbers. On a node with the interval $[p_i, q_i]$ in
 * variable $i$, the points of variable $i$ are the $N_i + 1$ Chebyshev-Lobatto points of
 * $[p_i, q_i]$, where $N_i = 2^{l_i}$ and the level $l_i$ starts at the level-min of the
 * variable. The builder evaluates $F$ on the grid formed by the tensor product of the
 * points of all variables, and a Chebyshev transform along each variable gives the
 * coefficients $c^{k_0 \dots k_{n - 1}}$. The series passes the test when
 * $$D \equiv \sum_i D_i \le \varepsilon \equiv \max (\varepsilon_\mathrm{rel} V, \varepsilon_\mathrm{abs}),$$
 * where $V$ is the largest $|F|$ on the current grid of the node,
 * $\varepsilon_\mathrm{rel}$ and $\varepsilon_\mathrm{abs}$ are the relative and absolute
 * tolerances of the construction, the second in the units of $F$, and
 * $$D_i \equiv 2 \sum_{k_i > N_i / 2} \sum_{k_{i'}, i' \ne i} |c^{k_0 \dots k_{n - 1}}|$$
 * is the exact $\ell_1$ change of the coefficients when the grid of variable $i$ is
 * replaced by the grid of the level below, which keeps every other point. The sum $D$ is
 * the estimate of the error of the series; it is not a bound on the distance to $F$. A
 * peak of $F$ narrower than the spacing of the grid and located between two points is
 * absent from the values of $F$ on the grid, and the series can pass the test without
 * representing it. When the position and the width of such a peak are known, the domain is
 * split there before the construction.
 *
 * A node whose series fails the test is refined in the variable $i$ with the largest $D_i$
 * among the variables whose level-min is below their level-cap, the first of them on a
 * tie. The level $l_i$ is raised by one: the values of $F$ already computed are kept, and
 * $F$ is evaluated only at the new points. The node is split at the midpoint of
 * $[p_i, q_i]$ instead when $l_i$ equals the level-cap of the variable, or is one below it
 * and the value of $D_i$ predicted for level $l_i + 1$ from a geometric decay of its
 * coefficients is above $\varepsilon$. Both children are expanded from the level-min of
 * every variable. A variable whose level-min equals its level-cap has a fixed number of
 * points: it is never chosen for refinement, and its $D_i$ is part of $D$.
 *
 * The tolerance $\varepsilon$ must be above the error of the values of $F$. The sum $D$
 * adds the absolute values of many coefficients, and each coefficient carries the error of
 * the values given by the callback, which is at least their rounding, of order
 * $10^{-16} |F|$. A node whose $\varepsilon$ is below the value of $D$ produced by that
 * error does not pass the test at any level. It is bisected until
 * #NcmSpectralTensorBuilder:max-depth is reached, and the build ends with the error
 * %NCM_SPECTRAL_TENSOR_BUILDER_ERROR_MAX_DEPTH. A build that fails in this way may
 * evaluate $F$ many times before the error. For a function of large $|F|$ the limit is
 * reached first by the absolute tolerance: with $|F|$ of order $10^4$,
 * $\varepsilon_\mathrm{abs} = 10^{-12}$ is below the rounding of the values and the build
 * fails, while $\varepsilon_\mathrm{rel} = 10^{-12}$ gives $\varepsilon = 10^{-12} V$ on
 * each node and the build ends.
 *
 * A node whose series passes the test becomes a panel node, which stores all its
 * coefficients by default. With #NcmSpectralTensorBuilder:trim-frac $f > 0$, slices of
 * highest index are removed one at a time, each time the slice with the smallest $\ell_1$
 * norm, while the $\ell_1$ norm $S$ of all the coefficients removed satisfies
 * $S \le f \varepsilon$. Since $|T_k| \le 1$, the stored series differs from the full
 * series by at most $S$ at every point of the panel, so the estimated error of the stored
 * series is at most $D + S \le (1 + f) \varepsilon$.
 *
 * A build takes its scratch space from a pool, so one builder may be used by builds in
 * several threads. Its properties must not change during a build.
 *
 * See <a href="../../theory/ncm/algebra/spectral_tensor.html">Tensor-Product Chebyshev
 * Tables</a> for the formulas.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/algebra/ncm_spectral_tensor_builder.h"
#include "ncm/algebra/ncm_spectral.h"
#include "ncm/algebra/ncm_spectral_private.h"
#include "ncm/core/ncm_memory_pool.h"
#include "ncm/core/ncm_util.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <string.h>
#include <fftw3.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

/* Largest level of a variable, 2^10 + 1 points */
#define NCM_SPECTRAL_TENSOR_BUILDER_LEVEL_MAX (10)

enum
{
  PROP_0,
  PROP_LEVEL_MIN,
  PROP_LEVEL_CAP,
  PROP_MAX_DEPTH,
  PROP_TRIM_FRAC,
};

struct _NcmSpectralTensorBuilder
{
  /*< private >*/
  GObject parent_instance;
  guint level_min;
  guint level_cap;
  guint max_depth;
  gdouble trim_frac;
  NcmMemoryPool *work_pool; /* NcmSpectralTensorBuilderWork, one per concurrent build */
};

/* A buffer from fftw_malloc and its capacity; growing it discards its contents */
typedef struct _NcmSpectralTensorBuilderBuf
{
  gdouble *data;
  gsize cap;
} NcmSpectralTensorBuilderBuf;

/* The scratch space of one build */
typedef struct _NcmSpectralTensorBuilderWork
{
  NcmSpectral *spectral;                                 /* the transforms from node values to Chebyshev coefficients */
  GArray *xc[NCM_SPECTRAL_TENSOR_BUILDER_LEVEL_MAX + 1]; /* xc[c]: cos (pi q / 2^c), q = 0, ..., 2^c, NULL until used */
  GArray *nodes[NCM_SPECTRAL_TENSOR_MAX_VARS];           /* the nodes of each variable at its current level */
  GArray *b[NCM_SPECTRAL_TENSOR_MAX_VARS];               /* the l1 masses of the coefficients per index of each variable */
  GArray *new_nodes;                                     /* the nodes a doubling adds */
  GArray *fval;                                          /* the values the callback returns for a doubling */
  GArray *block;                                         /* the stored block of a trimmed panel */
  NcmSpectralTensorBuilderBuf grid;                      /* the values on the grid of the current levels, last variable fastest */
  NcmSpectralTensorBuilderBuf grid_new;                  /* the grid a doubling writes */
  NcmSpectralTensorBuilderBuf t1;                        /* the intermediate transforms */
  NcmSpectralTensorBuilderBuf t2;
  NcmSpectralTensorBuilderBuf coeffs; /* the coefficients, first variable fastest */
} NcmSpectralTensorBuilderWork;

/* The arguments and the state of one build */
typedef struct _NcmSpectralTensorBuilderCall
{
  NcmSpectralTensorF F;
  gpointer user_data;
  NcmSpectralTensor *st;
  guint n_vars;
  guint level_min[NCM_SPECTRAL_TENSOR_MAX_VARS];
  guint level_cap[NCM_SPECTRAL_TENSOR_MAX_VARS];
  gdouble reltol;
  gdouble abstol;
  guint max_depth;
  gdouble trim_frac;
  guint64 n_eval;
  gboolean failed;
  gdouble failed_lower[NCM_SPECTRAL_TENSOR_MAX_VARS];
  gdouble failed_upper[NCM_SPECTRAL_TENSOR_MAX_VARS];
  guint failed_depth;
  NcmSpectralTensorBuilderWork *work;
} NcmSpectralTensorBuilderCall;

/* *INDENT-OFF* */
G_DEFINE_QUARK (ncm-spectral-tensor-builder-error, ncm_spectral_tensor_builder_error)
/* *INDENT-ON* */
G_DEFINE_TYPE (NcmSpectralTensorBuilder, ncm_spectral_tensor_builder, G_TYPE_OBJECT)

static gpointer _ncm_spectral_tensor_builder_work_new (gpointer userdata);
static void _ncm_spectral_tensor_builder_work_free (gpointer p);

static void
ncm_spectral_tensor_builder_init (NcmSpectralTensorBuilder *builder)
{
  builder->level_min = 0;
  builder->level_cap = 0;
  builder->max_depth = 0;
  builder->trim_frac = 0.0;
  builder->work_pool = ncm_memory_pool_new (&_ncm_spectral_tensor_builder_work_new, NULL, &_ncm_spectral_tensor_builder_work_free);
}

static void
_ncm_spectral_tensor_builder_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSpectralTensorBuilder *builder = NCM_SPECTRAL_TENSOR_BUILDER (object);

  g_return_if_fail (NCM_IS_SPECTRAL_TENSOR_BUILDER (object));

  switch (prop_id)
  {
    case PROP_LEVEL_MIN:
      ncm_spectral_tensor_builder_set_level_min (builder, g_value_get_uint (value));
      break;
    case PROP_LEVEL_CAP:
      ncm_spectral_tensor_builder_set_level_cap (builder, g_value_get_uint (value));
      break;
    case PROP_MAX_DEPTH:
      ncm_spectral_tensor_builder_set_max_depth (builder, g_value_get_uint (value));
      break;
    case PROP_TRIM_FRAC:
      ncm_spectral_tensor_builder_set_trim_frac (builder, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_spectral_tensor_builder_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSpectralTensorBuilder *builder = NCM_SPECTRAL_TENSOR_BUILDER (object);

  g_return_if_fail (NCM_IS_SPECTRAL_TENSOR_BUILDER (object));

  switch (prop_id)
  {
    case PROP_LEVEL_MIN:
      g_value_set_uint (value, builder->level_min);
      break;
    case PROP_LEVEL_CAP:
      g_value_set_uint (value, builder->level_cap);
      break;
    case PROP_MAX_DEPTH:
      g_value_set_uint (value, builder->max_depth);
      break;
    case PROP_TRIM_FRAC:
      g_value_set_double (value, builder->trim_frac);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_spectral_tensor_builder_finalize (GObject *object)
{
  NcmSpectralTensorBuilder *builder = NCM_SPECTRAL_TENSOR_BUILDER (object);

  ncm_memory_pool_free (builder->work_pool, TRUE);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spectral_tensor_builder_parent_class)->finalize (object);
}

static void
ncm_spectral_tensor_builder_class_init (NcmSpectralTensorBuilderClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_spectral_tensor_builder_set_property;
  object_class->get_property = &_ncm_spectral_tensor_builder_get_property;
  object_class->finalize     = &_ncm_spectral_tensor_builder_finalize;

  /**
   * NcmSpectralTensorBuilder:level-min:
   *
   * The starting level of every variable when ncm_spectral_tensor_builder_build() is
   * given no levels: $2^l + 1$ points.
   */
  g_object_class_install_property (object_class,
                                   PROP_LEVEL_MIN,
                                   g_param_spec_uint ("level-min",
                                                      NULL,
                                                      "Starting level",
                                                      2, 9, 3,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSpectralTensorBuilder:level-cap:
   *
   * The largest level of every variable when ncm_spectral_tensor_builder_build() is given
   * no levels; a node that needs more is split at its midpoint.
   */
  g_object_class_install_property (object_class,
                                   PROP_LEVEL_CAP,
                                   g_param_spec_uint ("level-cap",
                                                      NULL,
                                                      "Largest level",
                                                      3, NCM_SPECTRAL_TENSOR_BUILDER_LEVEL_MAX, 7,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSpectralTensorBuilder:max-depth:
   *
   * The largest number of bisections below a node that was empty when the build started,
   * along any path; a node that needs more fails the build with
   * %NCM_SPECTRAL_TENSOR_BUILDER_ERROR_MAX_DEPTH.
   */
  g_object_class_install_property (object_class,
                                   PROP_MAX_DEPTH,
                                   g_param_spec_uint ("max-depth",
                                                      NULL,
                                                      "Largest number of bisections below a starting node",
                                                      0, G_MAXUINT, 30,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSpectralTensorBuilder:trim-frac:
   *
   * The fraction $f \in [0, 1)$ of the tolerance $\varepsilon$ of a node that its stored
   * coefficients may give away: slices of highest index of the coefficients of a panel node
   * are removed while the $\ell_1$ norm of everything removed is at most $f \varepsilon$.
   * With 0, every coefficient is stored.
   */
  g_object_class_install_property (object_class,
                                   PROP_TRIM_FRAC,
                                   g_param_spec_double ("trim-frac",
                                                        NULL,
                                                        "Fraction of the tolerance given to the removed coefficients",
                                                        0.0, 1.0, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

static gpointer
_ncm_spectral_tensor_builder_work_new (gpointer userdata)
{
  NcmSpectralTensorBuilderWork *work = g_new0 (NcmSpectralTensorBuilderWork, 1);
  guint i;

  work->spectral  = ncm_spectral_new ();
  work->new_nodes = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->fval      = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->block     = g_array_new (FALSE, FALSE, sizeof (gdouble));

  for (i = 0; i <= NCM_SPECTRAL_TENSOR_BUILDER_LEVEL_MAX; i++)
    work->xc[i] = NULL;

  for (i = 0; i < NCM_SPECTRAL_TENSOR_MAX_VARS; i++)
  {
    work->nodes[i] = g_array_new (FALSE, FALSE, sizeof (gdouble));
    work->b[i]     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  }

  return work;
}

static void
_ncm_spectral_tensor_builder_work_free (gpointer p)
{
  NcmSpectralTensorBuilderWork *work = p;
  guint i;

  ncm_spectral_free (work->spectral);
  g_array_unref (work->new_nodes);
  g_array_unref (work->fval);
  g_array_unref (work->block);

  for (i = 0; i <= NCM_SPECTRAL_TENSOR_BUILDER_LEVEL_MAX; i++)
    g_clear_pointer (&work->xc[i], g_array_unref);

  for (i = 0; i < NCM_SPECTRAL_TENSOR_MAX_VARS; i++)
  {
    g_array_unref (work->nodes[i]);
    g_array_unref (work->b[i]);
  }

  g_clear_pointer (&work->grid.data, fftw_free);
  g_clear_pointer (&work->grid_new.data, fftw_free);
  g_clear_pointer (&work->t1.data, fftw_free);
  g_clear_pointer (&work->t2.data, fftw_free);
  g_clear_pointer (&work->coeffs.data, fftw_free);

  g_free (work);
}

/**
 * ncm_spectral_tensor_builder_new:
 *
 * Creates a builder with the default properties.
 *
 * Returns: (transfer full): a new #NcmSpectralTensorBuilder.
 */
NcmSpectralTensorBuilder *
ncm_spectral_tensor_builder_new (void)
{
  return g_object_new (NCM_TYPE_SPECTRAL_TENSOR_BUILDER, NULL);
}

/**
 * ncm_spectral_tensor_builder_ref:
 * @builder: a #NcmSpectralTensorBuilder
 *
 * Increases the reference count of @builder by one.
 *
 * Returns: (transfer full): @builder.
 */
NcmSpectralTensorBuilder *
ncm_spectral_tensor_builder_ref (NcmSpectralTensorBuilder *builder)
{
  return g_object_ref (builder);
}

/**
 * ncm_spectral_tensor_builder_free:
 * @builder: a #NcmSpectralTensorBuilder
 *
 * Decreases the reference count of @builder by one.
 */
void
ncm_spectral_tensor_builder_free (NcmSpectralTensorBuilder *builder)
{
  g_object_unref (builder);
}

/**
 * ncm_spectral_tensor_builder_clear:
 * @builder: a #NcmSpectralTensorBuilder
 *
 * Decreases the reference count of *@builder by one and sets *@builder to %NULL.
 */
void
ncm_spectral_tensor_builder_clear (NcmSpectralTensorBuilder **builder)
{
  g_clear_object (builder);
}

/**
 * ncm_spectral_tensor_builder_set_level_min:
 * @builder: a #NcmSpectralTensorBuilder
 * @level_min: the starting level
 *
 * Sets #NcmSpectralTensorBuilder:level-min.
 */
void
ncm_spectral_tensor_builder_set_level_min (NcmSpectralTensorBuilder *builder, guint level_min)
{
  builder->level_min = level_min;
}

/**
 * ncm_spectral_tensor_builder_set_level_cap:
 * @builder: a #NcmSpectralTensorBuilder
 * @level_cap: the largest level
 *
 * Sets #NcmSpectralTensorBuilder:level-cap.
 */
void
ncm_spectral_tensor_builder_set_level_cap (NcmSpectralTensorBuilder *builder, guint level_cap)
{
  builder->level_cap = level_cap;
}

/**
 * ncm_spectral_tensor_builder_set_max_depth:
 * @builder: a #NcmSpectralTensorBuilder
 * @max_depth: the largest number of bisections
 *
 * Sets #NcmSpectralTensorBuilder:max-depth.
 */
void
ncm_spectral_tensor_builder_set_max_depth (NcmSpectralTensorBuilder *builder, guint max_depth)
{
  builder->max_depth = max_depth;
}

/**
 * ncm_spectral_tensor_builder_set_trim_frac:
 * @builder: a #NcmSpectralTensorBuilder
 * @trim_frac: fraction of the tolerance given to the removed coefficients, in $[0, 1)$
 *
 * Sets #NcmSpectralTensorBuilder:trim-frac. A value outside $[0, 1)$ aborts.
 */
void
ncm_spectral_tensor_builder_set_trim_frac (NcmSpectralTensorBuilder *builder, gdouble trim_frac)
{
  if (!((trim_frac >= 0.0) && (trim_frac < 1.0)))
    g_error ("ncm_spectral_tensor_builder_set_trim_frac: %g is outside [0, 1).", trim_frac);

  builder->trim_frac = trim_frac;
}

/**
 * ncm_spectral_tensor_builder_get_level_min:
 * @builder: a #NcmSpectralTensorBuilder
 *
 * Returns: #NcmSpectralTensorBuilder:level-min.
 */
guint
ncm_spectral_tensor_builder_get_level_min (NcmSpectralTensorBuilder *builder)
{
  return builder->level_min;
}

/**
 * ncm_spectral_tensor_builder_get_level_cap:
 * @builder: a #NcmSpectralTensorBuilder
 *
 * Returns: #NcmSpectralTensorBuilder:level-cap.
 */
guint
ncm_spectral_tensor_builder_get_level_cap (NcmSpectralTensorBuilder *builder)
{
  return builder->level_cap;
}

/**
 * ncm_spectral_tensor_builder_get_max_depth:
 * @builder: a #NcmSpectralTensorBuilder
 *
 * Returns: #NcmSpectralTensorBuilder:max-depth.
 */
guint
ncm_spectral_tensor_builder_get_max_depth (NcmSpectralTensorBuilder *builder)
{
  return builder->max_depth;
}

/**
 * ncm_spectral_tensor_builder_get_trim_frac:
 * @builder: a #NcmSpectralTensorBuilder
 *
 * Returns: #NcmSpectralTensorBuilder:trim-frac.
 */
gdouble
ncm_spectral_tensor_builder_get_trim_frac (NcmSpectralTensorBuilder *builder)
{
  return builder->trim_frac;
}

static void _ncm_spectral_tensor_builder_check (NcmSpectralTensorBuilderCall *call);
static void _ncm_spectral_tensor_builder_work_prepare (NcmSpectralTensorBuilderCall *call);
static void _ncm_spectral_tensor_builder_expand (NcmSpectralTensorBuilderCall *call, guint node, guint depth);

/**
 * ncm_spectral_tensor_builder_build: (skip)
 * @builder: a #NcmSpectralTensorBuilder
 * @F: (scope call) (closure user_data): the function $F$
 * @st: the table to fill, with its domain set
 * @level_min: (nullable): the starting level of each variable, or %NULL
 * @level_cap: (nullable): the largest level of each variable, or %NULL
 * @reltol: relative tolerance of a panel
 * @abstol: absolute tolerance of a panel, in units of $F$
 * @user_data: user data for @F
 * @n_eval: (out) (optional): number of values of $F$ computed
 * @error: a #GError
 *
 * For a table with empty nodes, expands $F$ on each of them as described in
 * #NcmSpectralTensorBuilder, so that every empty node becomes a panel node or a split node
 * whose descendants are split nodes and panel nodes. The split nodes and panel nodes @st
 * already has are kept. @level_min
 * and @level_cap have one entry per variable of @st; %NULL gives
 * #NcmSpectralTensorBuilder:level-min or #NcmSpectralTensorBuilder:level-cap to every
 * variable. A variable with equal entries has a fixed number of points. The call is an
 * error unless the following hold: the domain of @st is set, the levels satisfy
 * $1 \le l_\mathrm{min} \le l_\mathrm{cap} \le 10$ with $l_\mathrm{min} < l_\mathrm{cap}$
 * for at least one variable, and the tolerances @reltol and @abstol are non-negative and
 * not both zero.
 *
 * A node that needs a split after #NcmSpectralTensorBuilder:max-depth bisections, or whose
 * midpoint is not strictly inside its interval in floating point, sets @error to
 * %NCM_SPECTRAL_TENSOR_BUILDER_ERROR_MAX_DEPTH, and @st is left as by
 * ncm_spectral_tensor_reset(): its domain, with the root an empty node.
 */
void
ncm_spectral_tensor_builder_build (NcmSpectralTensorBuilder *builder, NcmSpectralTensorF F, NcmSpectralTensor *st,
                                   const guint *level_min, const guint *level_cap, gdouble reltol, gdouble abstol,
                                   gpointer user_data, guint64 *n_eval, GError **error)
{
  NcmSpectralTensorBuilderWork **ptr = ncm_memory_pool_get (builder->work_pool);
  NcmSpectralTensorBuilderCall call;
  guint n_nodes, node, i;

  call.F            = F;
  call.user_data    = user_data;
  call.st           = st;
  call.n_vars       = ncm_spectral_tensor_get_n_vars (st);
  call.reltol       = reltol;
  call.abstol       = abstol;
  call.max_depth    = builder->max_depth;
  call.trim_frac    = builder->trim_frac;
  call.n_eval       = 0;
  call.failed       = FALSE;
  call.failed_depth = 0;
  call.work         = *ptr;

  for (i = 0; i < call.n_vars; i++)
  {
    call.level_min[i] = (level_min != NULL) ? level_min[i] : builder->level_min;
    call.level_cap[i] = (level_cap != NULL) ? level_cap[i] : builder->level_cap;
  }

  _ncm_spectral_tensor_builder_check (&call);
  _ncm_spectral_tensor_builder_work_prepare (&call);

  /* The nodes that are empty now; the nodes an expansion adds are filled by it */
  n_nodes = ncm_spectral_tensor_get_n_nodes (st);

  for (node = 0; node < n_nodes; node++)
  {
    if (ncm_spectral_tensor_get_node_kind (st, node) == NCM_SPECTRAL_TENSOR_NODE_KIND_EMPTY)
      _ncm_spectral_tensor_builder_expand (&call, node, 0);
  }

  ncm_memory_pool_return (ptr);

  if (n_eval != NULL)
    *n_eval = call.n_eval;

  if (call.failed)
  {
    GString *panel = g_string_new (NULL);

    for (i = 0; i < call.n_vars; i++)
      g_string_append_printf (panel, "%s[% 22.15g, % 22.15g]", (i > 0) ? " x " : "", call.failed_lower[i], call.failed_upper[i]);

    ncm_spectral_tensor_reset (st);
    ncm_util_set_or_call_error (error, NCM_SPECTRAL_TENSOR_BUILDER_ERROR, NCM_SPECTRAL_TENSOR_BUILDER_ERROR_MAX_DEPTH,
                                "ncm_spectral_tensor_builder_build: the subdomain %s does not pass the test after %u bisections.",
                                panel->str, call.failed_depth);
    g_string_free (panel, TRUE);
  }
}

/* The arguments of a build, checked before anything is evaluated */
static void
_ncm_spectral_tensor_builder_check (NcmSpectralTensorBuilderCall *call)
{
  gboolean refinable = FALSE;
  guint i;

  if (!(call->reltol >= 0.0) || !(call->abstol >= 0.0) || ((call->reltol == 0.0) && (call->abstol == 0.0)))
    g_error ("ncm_spectral_tensor_builder_build: reltol %g and abstol %g must be non-negative and not both zero.",
             call->reltol, call->abstol);

  if (ncm_spectral_tensor_get_n_nodes (call->st) == 0)
    g_error ("ncm_spectral_tensor_builder_build: the domain of the table is not set.");

  for (i = 0; i < call->n_vars; i++)
  {
    if ((call->level_min[i] < 1) || (call->level_min[i] > call->level_cap[i]) ||
        (call->level_cap[i] > NCM_SPECTRAL_TENSOR_BUILDER_LEVEL_MAX))
      g_error ("ncm_spectral_tensor_builder_build: the levels %u and %u of variable %u must satisfy "
               "1 <= level-min <= level-cap <= %u.", call->level_min[i], call->level_cap[i], i,
               NCM_SPECTRAL_TENSOR_BUILDER_LEVEL_MAX);

    refinable = refinable || (call->level_min[i] < call->level_cap[i]);
  }

  if (!refinable)
    g_error ("ncm_spectral_tensor_builder_build: every variable has fixed nodes, level-min = level-cap.");
}

/* Ensures that buf holds len values; its contents are lost when it grows */
static void
_ncm_spectral_tensor_builder_reserve (NcmSpectralTensorBuilderBuf *buf, gsize len)
{
  if (len <= buf->cap)
    return;

  g_clear_pointer (&buf->data, fftw_free);
  buf->data = fftw_malloc (sizeof (gdouble) * len);
  buf->cap  = len;
}

/*
 * Sizes the point and coefficient-mass arrays for the level caps of the call and computes
 * the cosine tables of those caps: xc[c][q] = cos (pi q / 2^c), with the middle entry set
 * to 0.
 */
static void
_ncm_spectral_tensor_builder_work_prepare (NcmSpectralTensorBuilderCall *call)
{
  NcmSpectralTensorBuilderWork *work = call->work;
  guint i, k;

  for (i = 0; i < call->n_vars; i++)
  {
    const guint cap = call->level_cap[i];
    const guint np  = (1u << cap) + 1;

    g_array_set_size (work->nodes[i], np);
    g_array_set_size (work->b[i], np);

    if (work->xc[cap] == NULL)
    {
      gdouble *xc;

      work->xc[cap] = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), np);
      g_array_set_size (work->xc[cap], np);
      xc = &g_array_index (work->xc[cap], gdouble, 0);

      for (k = 0; k < np; k++)
        xc[k] = cos (M_PI * k / (np - 1.0));

      xc[(np - 1) / 2] = 0.0;
    }
  }
}

/* The points of variable i on [a, b] at level l */
static void
_ncm_spectral_tensor_builder_nodes (NcmSpectralTensorBuilderCall *call, guint i, gdouble a, gdouble b, guint l)
{
  const guint n      = 1u << l;
  const guint sh     = call->level_cap[i] - l;
  const gdouble mid  = 0.5 * (a + b);
  const gdouble half = 0.5 * (b - a);
  const gdouble *xc  = &g_array_index (call->work->xc[call->level_cap[i]], gdouble, 0);
  gdouble *nodes     = &g_array_index (call->work->nodes[i], gdouble, 0);
  guint q;

  for (q = 0; q <= n; q++)
    nodes[q] = mid + half * xc[q << sh];
}

/* Samples F on the full grid of the levels l into work->grid */
static void
_ncm_spectral_tensor_builder_sample (NcmSpectralTensorBuilderCall *call, const guint *l)
{
  NcmSpectralTensorBuilderWork *work = call->work;
  guint n[NCM_SPECTRAL_TENSOR_MAX_VARS];
  const gdouble *x[NCM_SPECTRAL_TENSOR_MAX_VARS];
  gsize total = 1;
  guint i;

  for (i = 0; i < call->n_vars; i++)
  {
    n[i]   = (1u << l[i]) + 1;
    x[i]   = &g_array_index (work->nodes[i], gdouble, 0);
    total *= n[i];
  }

  _ncm_spectral_tensor_builder_reserve (&work->grid, total);
  call->F (call->n_vars, n, x, work->grid.data, call->user_data);
  call->n_eval += total;
}

/*
 * Doubles the level of variable j on [a, b]: moves its points from q to 2 q, computes the
 * new ones at the odd positions, evaluates F at the new points of variable j combined with
 * every point of the other variables, and interleaves those values with the grid in
 * work->grid.
 */
static void
_ncm_spectral_tensor_builder_double (NcmSpectralTensorBuilderCall *call, guint *l, guint j, gdouble a, gdouble b)
{
  NcmSpectralTensorBuilderWork *work = call->work;
  const guint N                      = 1u << l[j];
  const guint sh                     = call->level_cap[j] - (l[j] + 1);
  const gdouble mid                  = 0.5 * (a + b);
  const gdouble half                 = 0.5 * (b - a);
  const gdouble *xc                  = &g_array_index (work->xc[call->level_cap[j]], gdouble, 0);
  gdouble *nodes                     = &g_array_index (work->nodes[j], gdouble, 0);
  guint n[NCM_SPECTRAL_TENSOR_MAX_VARS];
  const gdouble *x[NCM_SPECTRAL_TENSOR_MAX_VARS];
  gsize outer = 1;
  gsize inner = 1;
  gsize o, q, count;
  guint i;

  g_array_set_size (work->new_nodes, N);

  for (q = N; q > 0; q--)
    nodes[2 * q] = nodes[q];

  for (q = 0; q < N; q++)
  {
    const gdouble v = mid + half * xc[(2 * q + 1) << sh];

    nodes[2 * q + 1]                            = v;
    g_array_index (work->new_nodes, gdouble, q) = v;
  }

  for (i = 0; i < call->n_vars; i++)
  {
    if (i == j)
    {
      n[i] = N;
      x[i] = &g_array_index (work->new_nodes, gdouble, 0);
    }
    else
    {
      n[i] = (1u << l[i]) + 1;
      x[i] = &g_array_index (work->nodes[i], gdouble, 0);

      if (i < j)
        outer *= n[i];
      else
        inner *= n[i];
    }
  }

  count = outer * N * inner;
  g_array_set_size (work->fval, count);
  call->F (call->n_vars, n, x, &g_array_index (work->fval, gdouble, 0), call->user_data);
  call->n_eval += count;

  _ncm_spectral_tensor_builder_reserve (&work->grid_new, outer * (2 * N + 1) * inner);

  for (o = 0; o < outer; o++)
  {
    const gdouble *old = work->grid.data + o * (N + 1) * inner;
    const gdouble *fv  = &g_array_index (work->fval, gdouble, o * N * inner);
    gdouble *dst       = work->grid_new.data + o * (2 * N + 1) * inner;

    for (q = 0; q <= N; q++)
      memcpy (dst + 2 * q * inner, old + q * inner, sizeof (gdouble) * inner);

    for (q = 0; q < N; q++)
      memcpy (dst + (2 * q + 1) * inner, fv + q * inner, sizeof (gdouble) * inner);
  }

  {
    const NcmSpectralTensorBuilderBuf tmp = work->grid;

    work->grid     = work->grid_new;
    work->grid_new = tmp;
  }

  l[j]++;
}

/*
 * Transforms the grid of the levels l into its Chebyshev coefficients in work->coeffs, the
 * first variable running fastest, and fills work->b[i] with the l1 masses of the
 * coefficients per index of variable i. Returns the largest |F| on the grid.
 *
 * The grid has axis i for variable i, the last one contiguous. Each step transforms along
 * axis 0, which holds variable t, then moves that axis to position n - 1 - t, so after the
 * last step the axes hold the variables n - 1, ..., 0. With two variables this is one
 * transpose between the two transforms.
 */
static gdouble
_ncm_spectral_tensor_builder_coeffs (NcmSpectralTensorBuilderCall *call, const guint *l)
{
  NcmSpectralTensorBuilderWork *work = call->work;
  const guint n_vars                 = call->n_vars;
  guint dims[NCM_SPECTRAL_TENSOR_MAX_VARS];
  guint ax[NCM_SPECTRAL_TENSOR_MAX_VARS];
  guint k[NCM_SPECTRAL_TENSOR_MAX_VARS];
  const gdouble *src = work->grid.data;
  gdouble *c;
  gsize total = 1;
  gdouble vs  = 0.0;
  gsize q;
  guint i, t;

  for (i = 0; i < n_vars; i++)
  {
    dims[i] = (1u << l[i]) + 1;
    ax[i]   = i;
    total  *= dims[i];
  }

  /* The maximum does not depend on the order, so the reduction may be reordered */
  {
    const gdouble * restrict grid = work->grid.data;

    #pragma omp simd reduction(max:vs)

    for (q = 0; q < total; q++)
      vs = MAX (vs, fabs (grid[q]));
  }

  _ncm_spectral_tensor_builder_reserve (&work->t1, total);
  _ncm_spectral_tensor_builder_reserve (&work->t2, total);
  _ncm_spectral_tensor_builder_reserve (&work->coeffs, total);

  for (t = 0; t < n_vars; t++)
  {
    const gsize A = dims[t];

    if (t + 1 == n_vars)
    {
      _ncm_spectral_transform (work->spectral, A, total / A, src, work->coeffs.data);
      break;
    }

    _ncm_spectral_transform (work->spectral, A, total / A, src, work->t1.data);

    {
      const guint p = n_vars - 1 - t;
      gsize B       = 1;
      gsize C       = 1;
      gsize a, b;
      guint new_ax[NCM_SPECTRAL_TENSOR_MAX_VARS];

      for (i = 1; i <= p; i++)
        B *= dims[ax[i]];

      for (i = p + 1; i < n_vars; i++)
        C *= dims[ax[i]];

      /* out[b][a][c] = in[a][b][c] */
      if (C == 1)
      {
        for (a = 0; a < A; a++)
        {
          const gdouble *row = work->t1.data + a * B;

          for (b = 0; b < B; b++)
            work->t2.data[b * A + a] = row[b];
        }
      }
      else
      {
        for (a = 0; a < A; a++)
          for (b = 0; b < B; b++)
            memcpy (work->t2.data + (b * A + a) * C, work->t1.data + (a * B + b) * C, sizeof (gdouble) * C);
      }

      for (i = 0; i < p; i++)
        new_ax[i] = ax[i + 1];

      new_ax[p] = ax[0];

      for (i = p + 1; i < n_vars; i++)
        new_ax[i] = ax[i];

      memcpy (ax, new_ax, sizeof (guint) * n_vars);
    }

    src = work->t2.data;
  }

  for (i = 0; i < n_vars; i++)
  {
    memset (&g_array_index (work->b[i], gdouble, 0), 0, sizeof (gdouble) * dims[i]);
    k[i] = 0;
  }

  c = work->coeffs.data;

  for (q = 0; q < total; q++)
  {
    const gdouble ac = fabs (c[q]);

    for (i = 0; i < n_vars; i++)
      g_array_index (work->b[i], gdouble, k[i]) += ac;

    for (i = 0; i < n_vars; i++)
    {
      if (++k[i] < dims[i])
        break;

      k[i] = 0;
    }
  }

  return vs;
}

/* The changes D_i of the coefficients of the levels l, from work->b; returns their sum */
static gdouble
_ncm_spectral_tensor_builder_change (NcmSpectralTensorBuilderCall *call, const guint *l, gdouble *Di)
{
  gdouble D = 0.0;
  guint i, k;

  for (i = 0; i < call->n_vars; i++)
  {
    const guint N    = 1u << l[i];
    const gdouble *b = &g_array_index (call->work->b[i], gdouble, 0);
    gdouble d        = 0.0;

    for (k = N / 2 + 1; k <= N; k++)
      d += 2.0 * b[k];

    Di[i] = d;
    D    += d;
  }

  return D;
}

/*
 * The value the test of the next doubling is predicted to have, from the l1 masses
 * b[0..M] of the coefficients along one variable: with B_1 and B_2 the masses of the
 * bands (M/4, M/2] and (M/2, M], geometric decay gives B_2 / B_1 = q (1 + q) and a next
 * band of B_2 q^2 (1 + q^2), as _ncm_spectral_batch_cannot_converge() in ncm_spectral.c.
 */
static gdouble
_ncm_spectral_tensor_builder_predict (const gdouble *b, guint M)
{
  gdouble B1 = 0.0;
  gdouble B2 = 0.0;
  guint i;

  for (i = M / 4 + 1; i <= M / 2; i++)
    B1 += b[i];

  for (i = M / 2 + 1; i <= M; i++)
    B2 += b[i];

  if (B2 == 0.0)
    return 0.0;

  if (B1 == 0.0)
    return GSL_POSINF;

  {
    const gdouble ratio = B2 / B1;
    const gdouble q     = 0.5 * (sqrt (1.0 + 4.0 * ratio) - 1.0);

    return 2.0 * B2 * q * q * (1.0 + q * q);
  }
}

/*
 * Makes the node a panel node of call->st with the coefficients of the levels l: all of
 * them, or with call->trim_frac > 0 the leading block left after removing slices of highest
 * index, each time the slice with the smallest l1 norm, while the norm of everything removed
 * is at most budget = trim_frac tol.
 */
static void
_ncm_spectral_tensor_builder_set_panel (NcmSpectralTensorBuilderCall *call, guint node, const guint *l, gdouble budget)
{
  NcmSpectralTensorBuilderWork *work = call->work;
  const guint n_vars                 = call->n_vars;
  const gdouble *c                   = work->coeffs.data;
  guint dims[NCM_SPECTRAL_TENSOR_MAX_VARS];
  guint K[NCM_SPECTRAL_TENSOR_MAX_VARS];
  gsize stride[NCM_SPECTRAL_TENSOR_MAX_VARS];
  guint k[NCM_SPECTRAL_TENSOR_MAX_VARS];
  gboolean cut = FALSE;
  gsize total  = 1;
  gsize len, q;
  guint i;

  for (i = 0; i < n_vars; i++)
  {
    dims[i]   = (1u << l[i]) + 1;
    K[i]      = dims[i];
    stride[i] = total;
    total    *= dims[i];
  }

  if (call->trim_frac > 0.0)
  {
    gdouble S = 0.0;

    while (TRUE)
    {
      gdouble norm[NCM_SPECTRAL_TENSOR_MAX_VARS];
      gint best = -1;

      for (i = 0; i < n_vars; i++)
      {
        norm[i] = 0.0;
        k[i]    = 0;
      }

      /* The l1 norms of the last slice of the block in every variable */
      while (TRUE)
      {
        gsize pos = 0;

        for (i = 0; i < n_vars; i++)
          pos += k[i] * stride[i];

        for (i = 0; i < n_vars; i++)
        {
          if (k[i] + 1 == K[i])
            norm[i] += fabs (c[pos]);
        }

        for (i = 0; i < n_vars; i++)
        {
          if (++k[i] < K[i])
            break;

          k[i] = 0;
        }

        if (i == n_vars)
          break;
      }

      for (i = 0; i < n_vars; i++)
      {
        if ((K[i] > 1) && ((best < 0) || (norm[i] < norm[best])))
          best = i;
      }

      if ((best < 0) || !(S + norm[best] <= budget))
        break;

      S += norm[best];
      K[best]--;
      cut = TRUE;
    }
  }

  if (!cut)
  {
    ncm_spectral_tensor_set_panel (call->st, node, n_vars, K, total, c);

    return;
  }

  len = 1;

  for (i = 0; i < n_vars; i++)
  {
    len *= K[i];
    k[i] = 0;
  }

  g_array_set_size (work->block, len);

  for (q = 0; q < len; q++)
  {
    gsize pos = 0;

    for (i = 0; i < n_vars; i++)
      pos += k[i] * stride[i];

    g_array_index (work->block, gdouble, q) = c[pos];

    for (i = 0; i < n_vars; i++)
    {
      if (++k[i] < K[i])
        break;

      k[i] = 0;
    }
  }

  ncm_spectral_tensor_set_panel (call->st, node, n_vars, K, len, &g_array_index (work->block, gdouble, 0));
}

/*
 * Expands F on the empty node from the starting levels until the test passes, increasing
 * the level of the refinable variable with the largest change or, at its level cap,
 * splitting the node at its midpoint in that variable and expanding the two children. A
 * node whose series passes the test becomes a panel node of call->st. A node at max_depth
 * that needs a split, or whose midpoint is not inside its interval, sets call->failed,
 * after which nothing more is evaluated.
 */
static void
_ncm_spectral_tensor_builder_expand (NcmSpectralTensorBuilderCall *call, guint node, guint depth)
{
  const guint n_vars = call->n_vars;
  gdouble lower[NCM_SPECTRAL_TENSOR_MAX_VARS];
  gdouble upper[NCM_SPECTRAL_TENSOR_MAX_VARS];
  guint l[NCM_SPECTRAL_TENSOR_MAX_VARS];
  gdouble Di[NCM_SPECTRAL_TENSOR_MAX_VARS];
  guint i;

  if (call->failed)
    return;

  /* The bounds are copied: a split moves the arrays of the table */
  {
    const gdouble *lo, *up;

    ncm_spectral_tensor_peek_node_bounds (call->st, node, &lo, &up);
    memcpy (lower, lo, sizeof (gdouble) * n_vars);
    memcpy (upper, up, sizeof (gdouble) * n_vars);
  }

  for (i = 0; i < n_vars; i++)
  {
    l[i] = call->level_min[i];
    _ncm_spectral_tensor_builder_nodes (call, i, lower[i], upper[i], l[i]);
  }

  _ncm_spectral_tensor_builder_sample (call, l);

  while (TRUE)
  {
    const gdouble vs  = _ncm_spectral_tensor_builder_coeffs (call, l);
    const gdouble tol = MAX (call->reltol * vs, call->abstol);
    const gdouble D   = _ncm_spectral_tensor_builder_change (call, l, Di);
    gboolean split;
    guint j = n_vars;

    if (D <= tol)
    {
      _ncm_spectral_tensor_builder_set_panel (call, node, l, call->trim_frac * tol);

      return;
    }

    /* The refinable variable with the largest change, the first one on a tie */
    for (i = 0; i < n_vars; i++)
    {
      if ((call->level_min[i] < call->level_cap[i]) && ((j == n_vars) || (Di[i] > Di[j])))
        j = i;
    }

    split = (l[j] >= call->level_cap[j]) ||
            ((l[j] + 1 == call->level_cap[j]) &&
             (_ncm_spectral_tensor_builder_predict (&g_array_index (call->work->b[j], gdouble, 0), 1u << l[j]) > tol));

    if (split)
    {
      const gdouble mid = 0.5 * (lower[j] + upper[j]);
      guint child;

      if ((depth >= call->max_depth) || !(lower[j] < mid) || !(mid < upper[j]))
      {
        call->failed       = TRUE;
        call->failed_depth = depth;
        memcpy (call->failed_lower, lower, sizeof (gdouble) * n_vars);
        memcpy (call->failed_upper, upper, sizeof (gdouble) * n_vars);

        return;
      }

      child = ncm_spectral_tensor_split (call->st, node, j, 1, &mid);

      _ncm_spectral_tensor_builder_expand (call, child, depth + 1);
      _ncm_spectral_tensor_builder_expand (call, child + 1, depth + 1);

      return;
    }

    _ncm_spectral_tensor_builder_double (call, l, j, lower[j], upper[j]);
  }
}

