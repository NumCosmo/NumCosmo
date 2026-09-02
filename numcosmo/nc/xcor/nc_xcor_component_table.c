/***************************************************************************
 *            nc_xcor_component_table.c
 *
 *  Wed September 02 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_xcor_component_table.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 *
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
 * NcXcorComponentTable:
 *
 * One tabulated component of a radial window: $(\chi, W)$ samples plus what the
 * window is a window of.
 *
 * A tracer is in general a sum of terms with different Bessel weights and
 * $\ell$ prefactors: a galaxy-count tracer has a density term weighted by
 * $j_\ell$, a redshift-space distortion term weighted by $j_\ell''$ and a
 * magnification term weighted by $j_\ell/(k\chi)^2$ with an $\ell(\ell+1)$
 * prefactor. Each of those is one of these objects; #NcXcorKernelTable holds
 * the list and integrates every component over its own support.
 *
 * The samples are reconstructed with a #NcmSplineBSpline of the requested
 * order (degree 7 by default: a cubic caps a 2000-sample window at
 * $\sim 10^{-8}$ while degree 7 reaches $\sim 10^{-15}$ from the same data).
 * Leading and trailing runs of exact zeros are trimmed, keeping one zero on
 * each side so the reconstruction still goes to zero at the reported edge.
 *
 * The convention is #NcXcorKernelRadial's: $\chi$ in Mpc, $P(k)$ at $z = 0$
 * and all redshift dependence -- growth included -- carried by $W$.
 *
 * For the kinds that carry $1/(k\chi)^2$ the reported support never starts at
 * the origin. A lensing window is linear in $\chi$ there, so $W/\chi^2$ is a
 * $1/\chi$ the solver handles but the origin itself is not a valid lower limit;
 * when the table's support starts below %NC_XCOR_COMPONENT_TABLE_CHI_FLOOR
 * (0.01 Mpc) it is moved there. The window omitted below that point is
 * $\sim W'(0)\,\chi_\mathrm{min}^2/2$, negligible for any survey window. The
 * floor is a fixed distance, not a fraction of the sample spacing: the closure's
 * $k$ range for these kinds grows as the inverse of the lower limit, and a
 * 4000-sample lensing table with the floor at a thousandth of its first sample
 * (1.5 Mpc) cost seven times the same window tabulated on 1000 samples.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/xcor/nc_xcor_component_table.h"
#include "ncm/spline/ncm_spline_bspline.h"
#include "nc_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcXcorComponentTable
{
  /*< private >*/
  GObject parent_instance;
  NcmVector *chi;
  NcmVector *W;
  NcXcorKernelTableKind kind;
  guint order;
  gboolean normalize;
  NcmSpline *spline;
  gdouble chi_min;
  gdouble chi_max;
  gdouble norm;
};

enum
{
  PROP_0,
  PROP_CHI,
  PROP_W,
  PROP_KIND,
  PROP_ORDER,
  PROP_NORMALIZE,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcXcorComponentTable, nc_xcor_component_table, G_TYPE_OBJECT)

static void
nc_xcor_component_table_init (NcXcorComponentTable *xcct)
{
  xcct->chi       = NULL;
  xcct->W         = NULL;
  xcct->kind      = NC_XCOR_KERNEL_TABLE_KIND_DENSITY;
  xcct->order     = NCM_SPLINE_BSPLINE_DEFAULT_ORDER;
  xcct->normalize = TRUE;
  xcct->spline    = NULL;
  xcct->chi_min   = 0.0;
  xcct->chi_max   = 0.0;
  xcct->norm      = 1.0;
}

static void
nc_xcor_component_table_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorComponentTable *xcct = NC_XCOR_COMPONENT_TABLE (object);

  g_return_if_fail (NC_IS_XCOR_COMPONENT_TABLE (object));

  switch (prop_id)
  {
    case PROP_CHI:
      ncm_vector_clear (&xcct->chi);
      xcct->chi = g_value_dup_object (value);
      break;
    case PROP_W:
      ncm_vector_clear (&xcct->W);
      xcct->W = g_value_dup_object (value);
      break;
    case PROP_KIND:
      xcct->kind = g_value_get_enum (value);
      break;
    case PROP_ORDER:
      xcct->order = g_value_get_uint (value);
      break;
    case PROP_NORMALIZE:
      xcct->normalize = g_value_get_boolean (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_xcor_component_table_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorComponentTable *xcct = NC_XCOR_COMPONENT_TABLE (object);

  g_return_if_fail (NC_IS_XCOR_COMPONENT_TABLE (object));

  switch (prop_id)
  {
    case PROP_CHI:
      g_value_set_object (value, xcct->chi);
      break;
    case PROP_W:
      g_value_set_object (value, xcct->W);
      break;
    case PROP_KIND:
      g_value_set_enum (value, xcct->kind);
      break;
    case PROP_ORDER:
      g_value_set_uint (value, xcct->order);
      break;
    case PROP_NORMALIZE:
      g_value_set_boolean (value, xcct->normalize);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

/*
 * Trim leading and trailing runs of exact zeros. One zero sample is kept on each
 * side so the reconstruction still goes to zero at the reported edge rather than
 * being cut where the window is already non-zero.
 */
static void
_nc_xcor_component_table_trim (NcmVector *W, gsize *first, gsize *last)
{
  const gsize len = ncm_vector_len (W);
  gsize i         = 0;
  gsize j         = len - 1;

  while ((i < len - 1) && (ncm_vector_get (W, i) == 0.0))
    i++;

  while ((j > i) && (ncm_vector_get (W, j) == 0.0))
    j--;

  *first = (i > 0) ? i - 1 : 0;
  *last  = (j < len - 1) ? j + 1 : len - 1;
}

static gboolean
_nc_xcor_component_table_kind_has_inverse_square (NcXcorKernelTableKind kind)
{
  return (kind == NC_XCOR_KERNEL_TABLE_KIND_SHEAR) || (kind == NC_XCOR_KERNEL_TABLE_KIND_CONVERGENCE);
}

static void
_nc_xcor_component_table_build (NcXcorComponentTable *xcct)
{
  gsize first, last, n;
  NcmVector *chi_t, *W_t;

  if ((xcct->chi == NULL) || (xcct->W == NULL))
    g_error ("nc_xcor_component_table: both chi and W must be given.");

  if (ncm_vector_len (xcct->chi) != ncm_vector_len (xcct->W))
    g_error ("nc_xcor_component_table: chi and W have different lengths (%u != %u).",
             ncm_vector_len (xcct->chi), ncm_vector_len (xcct->W));

  if (ncm_vector_len (xcct->chi) < xcct->order)
    g_error ("nc_xcor_component_table: a degree %u reconstruction needs at least %u samples, got %u.",
             xcct->order - 1, xcct->order, ncm_vector_len (xcct->chi));

  _nc_xcor_component_table_trim (xcct->W, &first, &last);
  n = last - first + 1;

  if (n < xcct->order)
    g_error ("nc_xcor_component_table: the window's support holds %u samples, "
             "fewer than the %u a degree %u reconstruction needs.",
             (guint) n, xcct->order, xcct->order - 1);

  chi_t = ncm_vector_get_subvector (xcct->chi, first, n);
  W_t   = ncm_vector_get_subvector (xcct->W, first, n);

  ncm_spline_clear (&xcct->spline);
  xcct->spline  = NCM_SPLINE (ncm_spline_bspline_new_full (xcct->order, chi_t, W_t, TRUE));
  xcct->chi_min = ncm_vector_get (chi_t, 0);
  xcct->chi_max = ncm_vector_get (chi_t, n - 1);

  xcct->norm = xcct->normalize ?
               ncm_spline_eval_integ (xcct->spline, xcct->chi_min, xcct->chi_max) : 1.0;

  if (!gsl_finite (xcct->norm) || (xcct->norm == 0.0))
    g_error ("nc_xcor_component_table: the tabulated window integrates to %g over "
             "[%g, %g] Mpc; it cannot be normalized.", xcct->norm, xcct->chi_min, xcct->chi_max);

  /* See the class documentation: 1/(k chi)^2 kinds cannot start at the origin. */
  if (_nc_xcor_component_table_kind_has_inverse_square (xcct->kind) && (xcct->chi_min < NC_XCOR_COMPONENT_TABLE_CHI_FLOOR))
  {
    if (xcct->chi_max <= NC_XCOR_COMPONENT_TABLE_CHI_FLOOR)
      g_error ("nc_xcor_component_table: a window with a 1/(k chi)^2 weight must extend "
               "beyond %g Mpc.", NC_XCOR_COMPONENT_TABLE_CHI_FLOOR);

    xcct->chi_min = NC_XCOR_COMPONENT_TABLE_CHI_FLOOR;
  }

  ncm_vector_free (chi_t);
  ncm_vector_free (W_t);
}

static void
nc_xcor_component_table_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_xcor_component_table_parent_class)->constructed (object);

  _nc_xcor_component_table_build (NC_XCOR_COMPONENT_TABLE (object));
}

static void
nc_xcor_component_table_dispose (GObject *object)
{
  NcXcorComponentTable *xcct = NC_XCOR_COMPONENT_TABLE (object);

  ncm_vector_clear (&xcct->chi);
  ncm_vector_clear (&xcct->W);
  ncm_spline_clear (&xcct->spline);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_component_table_parent_class)->dispose (object);
}

static void
nc_xcor_component_table_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_component_table_parent_class)->finalize (object);
}

static void
nc_xcor_component_table_class_init (NcXcorComponentTableClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &nc_xcor_component_table_set_property;
  object_class->get_property = &nc_xcor_component_table_get_property;
  object_class->constructed  = &nc_xcor_component_table_constructed;
  object_class->dispose      = &nc_xcor_component_table_dispose;
  object_class->finalize     = &nc_xcor_component_table_finalize;

  /**
   * NcXcorComponentTable:chi:
   *
   * Comoving distances of the samples, in Mpc, strictly increasing.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHI,
                                   g_param_spec_object ("chi",
                                                        NULL,
                                                        "Sample comoving distances in Mpc",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorComponentTable:W:
   *
   * Window samples $W(\chi)$, one per entry of #NcXcorComponentTable:chi.
   */
  g_object_class_install_property (object_class,
                                   PROP_W,
                                   g_param_spec_object ("W",
                                                        NULL,
                                                        "Window samples",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorComponentTable:kind:
   *
   * What the window is a window of, which fixes the Bessel weight and the
   * factors applied on top of it.
   */
  g_object_class_install_property (object_class,
                                   PROP_KIND,
                                   g_param_spec_enum ("kind",
                                                      NULL,
                                                      "Window kind",
                                                      NC_TYPE_XCOR_KERNEL_TABLE_KIND,
                                                      NC_XCOR_KERNEL_TABLE_KIND_DENSITY,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorComponentTable:order:
   *
   * B-spline order of the reconstruction; the degree is one less. The default is
   * degree 7, not a cubic: a cubic caps a 2000-sample window at $\sim 10^{-8}$
   * while degree 7 reaches $\sim 10^{-15}$ from the same data.
   */
  g_object_class_install_property (object_class,
                                   PROP_ORDER,
                                   g_param_spec_uint ("order",
                                                      NULL,
                                                      "B-spline order of the reconstruction",
                                                      2, 10, NCM_SPLINE_BSPLINE_DEFAULT_ORDER,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorComponentTable:normalize:
   *
   * Whether to rescale the window to unit integral over its support. Turn it off
   * when the supplied table is already normalized in the same convention, or
   * when it is one term of a tracer whose relative amplitudes must be kept.
   */
  g_object_class_install_property (object_class,
                                   PROP_NORMALIZE,
                                   g_param_spec_boolean ("normalize",
                                                         NULL,
                                                         "Rescale to unit integral over the support",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_xcor_component_table_new:
 * @chi: a #NcmVector of sample comoving distances in Mpc
 * @W: a #NcmVector of window samples
 *
 * Creates a density component of default order, normalized to unit integral.
 *
 * Returns: (transfer full): a new #NcXcorComponentTable
 */
NcXcorComponentTable *
nc_xcor_component_table_new (NcmVector *chi, NcmVector *W)
{
  return g_object_new (NC_TYPE_XCOR_COMPONENT_TABLE,
                       "chi", chi,
                       "W", W,
                       NULL);
}

/**
 * nc_xcor_component_table_new_full:
 * @chi: a #NcmVector of sample comoving distances in Mpc
 * @W: a #NcmVector of window samples
 * @kind: a #NcXcorKernelTableKind
 * @order: B-spline order of the reconstruction
 * @normalize: whether to rescale the window to unit integral
 *
 * Creates a component with every option set.
 *
 * Returns: (transfer full): a new #NcXcorComponentTable
 */
NcXcorComponentTable *
nc_xcor_component_table_new_full (NcmVector *chi, NcmVector *W, NcXcorKernelTableKind kind, guint order, gboolean normalize)
{
  return g_object_new (NC_TYPE_XCOR_COMPONENT_TABLE,
                       "chi", chi,
                       "W", W,
                       "kind", kind,
                       "order", order,
                       "normalize", normalize,
                       NULL);
}

/**
 * nc_xcor_component_table_ref:
 * @xcct: a #NcXcorComponentTable
 *
 * Increases the reference count of @xcct by one.
 *
 * Returns: (transfer full): @xcct
 */
NcXcorComponentTable *
nc_xcor_component_table_ref (NcXcorComponentTable *xcct)
{
  return g_object_ref (xcct);
}

/**
 * nc_xcor_component_table_free:
 * @xcct: a #NcXcorComponentTable
 *
 * Decreases the reference count of @xcct by one.
 */
void
nc_xcor_component_table_free (NcXcorComponentTable *xcct)
{
  g_object_unref (xcct);
}

/**
 * nc_xcor_component_table_clear:
 * @xcct: a #NcXcorComponentTable
 *
 * Decreases the reference count of *@xcct by one and sets *@xcct to NULL.
 */
void
nc_xcor_component_table_clear (NcXcorComponentTable **xcct)
{
  g_clear_object (xcct);
}

/**
 * nc_xcor_component_table_get_kind:
 * @xcct: a #NcXcorComponentTable
 *
 * Returns: the window kind.
 */
NcXcorKernelTableKind
nc_xcor_component_table_get_kind (NcXcorComponentTable *xcct)
{
  return xcct->kind;
}

/**
 * nc_xcor_component_table_get_order:
 * @xcct: a #NcXcorComponentTable
 *
 * Returns: the B-spline order of the reconstruction.
 */
guint
nc_xcor_component_table_get_order (NcXcorComponentTable *xcct)
{
  return xcct->order;
}

/**
 * nc_xcor_component_table_get_normalize:
 * @xcct: a #NcXcorComponentTable
 *
 * Returns: whether the window is rescaled to unit integral.
 */
gboolean
nc_xcor_component_table_get_normalize (NcXcorComponentTable *xcct)
{
  return xcct->normalize;
}

/**
 * nc_xcor_component_table_get_norm:
 * @xcct: a #NcXcorComponentTable
 *
 * Gets the integral of the supplied table over its support, which the window
 * is divided by when #NcXcorComponentTable:normalize is set, or 1 otherwise.
 *
 * Returns: the normalization constant.
 */
gdouble
nc_xcor_component_table_get_norm (NcXcorComponentTable *xcct)
{
  return xcct->norm;
}

/**
 * nc_xcor_component_table_peek_spline:
 * @xcct: a #NcXcorComponentTable
 *
 * Gets the reconstruction of the supplied samples, before normalization.
 *
 * Returns: (transfer none): the #NcmSpline
 */
NcmSpline *
nc_xcor_component_table_peek_spline (NcXcorComponentTable *xcct)
{
  return xcct->spline;
}

/**
 * nc_xcor_component_table_peek_knots:
 * @xcct: a #NcXcorComponentTable
 *
 * Gets the abscissae the reconstruction is built on: the supplied samples with
 * the leading and trailing zero runs trimmed.
 *
 * Returns: (transfer none): the knots, in Mpc
 */
NcmVector *
nc_xcor_component_table_peek_knots (NcXcorComponentTable *xcct)
{
  return ncm_spline_peek_xv (xcct->spline);
}

/**
 * nc_xcor_component_table_set_samples:
 * @xcct: a #NcXcorComponentTable
 * @chi: a #NcmVector of sample comoving distances in Mpc
 * @W: a #NcmVector of window samples
 *
 * Replaces the tabulated samples and rebuilds the reconstruction, keeping the
 * kind, order and normalization. The table may change length. This is how a
 * sampler feeds a new cosmology's window into a kernel that stays registered
 * with its #NcXcorSolver, so the solver's per-block integrators survive the
 * step; the owning kernel must be told through nc_xcor_kernel_mark_outdated(),
 * which nc_xcor_kernel_table_replace_samples() does.
 */
void
nc_xcor_component_table_set_samples (NcXcorComponentTable *xcct, NcmVector *chi, NcmVector *W)
{
  ncm_vector_clear (&xcct->chi);
  ncm_vector_clear (&xcct->W);
  xcct->chi = ncm_vector_ref (chi);
  xcct->W   = ncm_vector_ref (W);

  _nc_xcor_component_table_build (xcct);
}

/**
 * nc_xcor_component_table_get_support:
 * @xcct: a #NcXcorComponentTable
 * @chi_min: (out): lower end of the support, in Mpc
 * @chi_max: (out): upper end of the support, in Mpc
 *
 * Gets the interval the component is integrated over. It is the trimmed table's
 * range, except that a $1/(k\chi)^2$ kind never starts at the origin (see the
 * class documentation).
 */
void
nc_xcor_component_table_get_support (NcXcorComponentTable *xcct, gdouble *chi_min, gdouble *chi_max)
{
  *chi_min = xcct->chi_min;
  *chi_max = xcct->chi_max;
}

/**
 * nc_xcor_component_table_eval_W:
 * @xcct: a #NcXcorComponentTable
 * @chi: comoving distance in Mpc, inside the support
 *
 * Evaluates the normalized window.
 *
 * Returns: $W(\chi)$ divided by the normalization constant.
 */
gdouble
nc_xcor_component_table_eval_W (NcXcorComponentTable *xcct, gdouble chi)
{
  return ncm_spline_eval (xcct->spline, chi) / xcct->norm;
}

/**
 * nc_xcor_component_table_eval_kernel_factor:
 * @xcct: a #NcXcorComponentTable
 * @chi: comoving distance in Mpc
 * @k: wave number in Mpc$^{-1}$
 *
 * Evaluates the $(\chi, k)$ factor the kind carries: $1/(k\chi)^2$ for
 * %NC_XCOR_KERNEL_TABLE_KIND_SHEAR and %NC_XCOR_KERNEL_TABLE_KIND_CONVERGENCE,
 * 1 otherwise.
 *
 * Returns: the factor.
 */
gdouble
nc_xcor_component_table_eval_kernel_factor (NcXcorComponentTable *xcct, gdouble chi, gdouble k)
{
  if (_nc_xcor_component_table_kind_has_inverse_square (xcct->kind))
  {
    const gdouble kchi = k * chi;

    return 1.0 / (kchi * kchi);
  }

  return 1.0;
}

/**
 * nc_xcor_component_table_eval_prefactor:
 * @xcct: a #NcXcorComponentTable
 * @l: multipole $\ell$
 *
 * Evaluates the $\ell$ prefactor the kind carries:
 * $\sqrt{(\ell+2)(\ell+1)\ell(\ell-1)}$ for %NC_XCOR_KERNEL_TABLE_KIND_SHEAR
 * (zero where the radicand is not positive), $\ell(\ell+1)$ for
 * %NC_XCOR_KERNEL_TABLE_KIND_CONVERGENCE, 1 otherwise.
 *
 * Returns: the prefactor.
 */
gdouble
nc_xcor_component_table_eval_prefactor (NcXcorComponentTable *xcct, gint l)
{
  switch (xcct->kind)
  {
    case NC_XCOR_KERNEL_TABLE_KIND_SHEAR:
    {
      const gdouble v = (l + 2.0) * (l + 1.0) * l * (l - 1.0);

      return (v > 0.0) ? sqrt (v) : 0.0;
    }
    case NC_XCOR_KERNEL_TABLE_KIND_CONVERGENCE:
      return l * (l + 1.0);

    default:
      return 1.0;
  }
}

/**
 * nc_xcor_component_table_get_bessel_deriv:
 * @xcct: a #NcXcorComponentTable
 *
 * Gets the derivative order of the spherical Bessel weight: 2 for
 * %NC_XCOR_KERNEL_TABLE_KIND_RSD, 0 otherwise.
 *
 * Returns: the Bessel-derivative order.
 */
guint
nc_xcor_component_table_get_bessel_deriv (NcXcorComponentTable *xcct)
{
  return (xcct->kind == NC_XCOR_KERNEL_TABLE_KIND_RSD) ? 2 : 0;
}

