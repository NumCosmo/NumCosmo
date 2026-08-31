/***************************************************************************
 *            nc_xcor_kernel_table.c
 *
 *  Sat August 30 2026
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_xcor_kernel_table.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 *
 * numcosmo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * numcosmo is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcXcorKernelTable:
 *
 * Radial window supplied as a table of $(\chi, W)$ samples.
 *
 * Every other #NcXcorKernelRadial is a formula. This one is the case a survey
 * pipeline actually hands you: the window already evaluated on a grid, with no
 * closed form behind it. It reconstructs that table with a #NcmSplineBSpline
 * and is otherwise an ordinary radial kernel, so #NcXcor and #NcXcorSolver
 * drive it exactly as they drive the analytic ones.
 *
 * The reconstruction order is the point. A cubic spline never reaches machine
 * precision at any sample density, so tabulated input silently caps everything
 * downstream; on a 2000-sample window a cubic reconstruction is accurate to
 * $\sim 10^{-8}$ while degree 7 reaches $\sim 10^{-15}$ from the same data.
 * The default is therefore %NCM_SPLINE_BSPLINE_DEFAULT_ORDER, degree 7, not a
 * cubic.
 *
 * The convention is #NcXcorKernelRadial's: $\chi$ in Mpc, $P(k)$ at $z=0$, and
 * all redshift dependence -- growth included -- carried by $W$. A table
 * exported from CCL's tracer kernels or from the N5K challenge already
 * satisfies it.
 *
 * Leading and trailing runs of exact zeros are trimmed from the support, so the
 * oscillatory integral runs over the interval the window actually occupies
 * rather than the whole tabulated range.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/xcor/nc_xcor_kernel_table.h"
#include "ncm/spline/ncm_spline_bspline.h"
#include "nc_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcXcorKernelTable
{
  /*< private >*/
  NcXcorKernelRadial parent_instance;

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

G_DEFINE_TYPE (NcXcorKernelTable, nc_xcor_kernel_table, NC_TYPE_XCOR_KERNEL_RADIAL)

static void
nc_xcor_kernel_table_init (NcXcorKernelTable *xckt)
{
  xckt->chi       = NULL;
  xckt->W         = NULL;
  xckt->kind      = NC_XCOR_KERNEL_TABLE_KIND_DENSITY;
  xckt->order     = NCM_SPLINE_BSPLINE_DEFAULT_ORDER;
  xckt->normalize = TRUE;
  xckt->spline    = NULL;
  xckt->chi_min   = 0.0;
  xckt->chi_max   = 0.0;
  xckt->norm      = 1.0;
}

static void
nc_xcor_kernel_table_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorKernelTable *xckt = NC_XCOR_KERNEL_TABLE (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_TABLE (object));

  switch (prop_id)
  {
    case PROP_CHI:
      ncm_vector_clear (&xckt->chi);
      xckt->chi = g_value_dup_object (value);
      break;
    case PROP_W:
      ncm_vector_clear (&xckt->W);
      xckt->W = g_value_dup_object (value);
      break;
    case PROP_KIND:
      xckt->kind = g_value_get_enum (value);
      break;
    case PROP_ORDER:
      xckt->order = g_value_get_uint (value);
      break;
    case PROP_NORMALIZE:
      xckt->normalize = g_value_get_boolean (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_xcor_kernel_table_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorKernelTable *xckt = NC_XCOR_KERNEL_TABLE (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_TABLE (object));

  switch (prop_id)
  {
    case PROP_CHI:
      g_value_set_object (value, xckt->chi);
      break;
    case PROP_W:
      g_value_set_object (value, xckt->W);
      break;
    case PROP_KIND:
      g_value_set_enum (value, xckt->kind);
      break;
    case PROP_ORDER:
      g_value_set_uint (value, xckt->order);
      break;
    case PROP_NORMALIZE:
      g_value_set_boolean (value, xckt->normalize);
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
_nc_xcor_kernel_table_trim (NcmVector *chi, NcmVector *W, gsize *first, gsize *last)
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

static void
nc_xcor_kernel_table_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_xcor_kernel_table_parent_class)->constructed (object);
  {
    NcXcorKernelTable *xckt = NC_XCOR_KERNEL_TABLE (object);
    gsize first, last, n;
    NcmVector *chi_t, *W_t;

    if ((xckt->chi == NULL) || (xckt->W == NULL))
      g_error ("nc_xcor_kernel_table_constructed: both chi and W must be given.");

    if (ncm_vector_len (xckt->chi) != ncm_vector_len (xckt->W))
      g_error ("nc_xcor_kernel_table_constructed: chi and W have different lengths (%u != %u).",
               ncm_vector_len (xckt->chi), ncm_vector_len (xckt->W));

    if (ncm_vector_len (xckt->chi) < xckt->order)
      g_error ("nc_xcor_kernel_table_constructed: a degree %u reconstruction needs at least %u samples, got %u.",
               xckt->order - 1, xckt->order, ncm_vector_len (xckt->chi));

    _nc_xcor_kernel_table_trim (xckt->chi, xckt->W, &first, &last);
    n = last - first + 1;

    if (n < xckt->order)
      g_error ("nc_xcor_kernel_table_constructed: the window's support holds %u samples, "
               "fewer than the %u a degree %u reconstruction needs.",
               (guint) n, xckt->order, xckt->order - 1);

    chi_t = ncm_vector_get_subvector (xckt->chi, first, n);
    W_t   = ncm_vector_get_subvector (xckt->W, first, n);

    xckt->spline  = NCM_SPLINE (ncm_spline_bspline_new_full (xckt->order, chi_t, W_t, TRUE));
    xckt->chi_min = ncm_vector_get (chi_t, 0);
    xckt->chi_max = ncm_vector_get (chi_t, n - 1);

    xckt->norm = xckt->normalize ?
                 ncm_spline_eval_integ (xckt->spline, xckt->chi_min, xckt->chi_max) : 1.0;

    if (!gsl_finite (xckt->norm) || (xckt->norm == 0.0))
      g_error ("nc_xcor_kernel_table_constructed: the tabulated window integrates to %g over "
               "[%g, %g] Mpc; it cannot be normalized.", xckt->norm, xckt->chi_min, xckt->chi_max);

    ncm_vector_free (chi_t);
    ncm_vector_free (W_t);
  }
}

static void
nc_xcor_kernel_table_dispose (GObject *object)
{
  NcXcorKernelTable *xckt = NC_XCOR_KERNEL_TABLE (object);

  ncm_vector_clear (&xckt->chi);
  ncm_vector_clear (&xckt->W);
  ncm_spline_clear (&xckt->spline);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_table_parent_class)->dispose (object);
}

static void
nc_xcor_kernel_table_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_table_parent_class)->finalize (object);
}

static guint _nc_xcor_kernel_table_get_n_comps (NcXcorKernelRadial *xcka);
static gdouble _nc_xcor_kernel_table_eval_W_comp (NcXcorKernelRadial *xcka, guint comp, gdouble chi);
static void _nc_xcor_kernel_table_get_comp_support (NcXcorKernelRadial *xcka, guint comp, gdouble *chi_min, gdouble *chi_max);
static gdouble _nc_xcor_kernel_table_eval_kernel_factor (NcXcorKernelRadial *xcka, NcHICosmo *cosmo, gdouble chi, gdouble k);
static gdouble _nc_xcor_kernel_table_eval_prefactor (NcXcorKernelRadial *xcka, NcHICosmo *cosmo, gint l);

static void
nc_xcor_kernel_table_class_init (NcXcorKernelTableClass *klass)
{
  GObjectClass *object_class            = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class            = NCM_MODEL_CLASS (klass);
  NcXcorKernelRadialClass *parent_class = NC_XCOR_KERNEL_RADIAL_CLASS (klass);

  object_class->constructed = &nc_xcor_kernel_table_constructed;
  object_class->dispose     = &nc_xcor_kernel_table_dispose;
  object_class->finalize    = &nc_xcor_kernel_table_finalize;
  model_class->set_property = &nc_xcor_kernel_table_set_property;
  model_class->get_property = &nc_xcor_kernel_table_get_property;

  ncm_model_class_set_name_nick (model_class, "Tabulated radial window", "XcorTable");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  /**
   * NcXcorKernelTable:chi:
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
   * NcXcorKernelTable:W:
   *
   * Window samples $W(\chi)$, one per entry of #NcXcorKernelTable:chi.
   */
  g_object_class_install_property (object_class,
                                   PROP_W,
                                   g_param_spec_object ("W",
                                                        NULL,
                                                        "Window samples",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelTable:kind:
   *
   * What the window is a window of, which fixes the factors applied on top of it.
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
   * NcXcorKernelTable:order:
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
   * NcXcorKernelTable:normalize:
   *
   * Whether to rescale the window to unit integral over its support. Turn it off
   * only when the supplied table is already normalized in the same convention.
   */
  g_object_class_install_property (object_class,
                                   PROP_NORMALIZE,
                                   g_param_spec_boolean ("normalize",
                                                         NULL,
                                                         "Rescale to unit integral over the support",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  parent_class->get_n_comps        = &_nc_xcor_kernel_table_get_n_comps;
  parent_class->eval_W_comp        = &_nc_xcor_kernel_table_eval_W_comp;
  parent_class->get_comp_support   = &_nc_xcor_kernel_table_get_comp_support;
  parent_class->eval_kernel_factor = &_nc_xcor_kernel_table_eval_kernel_factor;
  parent_class->eval_prefactor     = &_nc_xcor_kernel_table_eval_prefactor;

  ncm_model_class_add_impl_flag (model_class, NC_XCOR_KERNEL_IMPL_ALL);
}

static guint
_nc_xcor_kernel_table_get_n_comps (NcXcorKernelRadial *xcka)
{
  return 1;
}

static gdouble
_nc_xcor_kernel_table_eval_W_comp (NcXcorKernelRadial *xcka, guint comp, gdouble chi)
{
  NcXcorKernelTable *xckt = NC_XCOR_KERNEL_TABLE (xcka);

  return ncm_spline_eval (xckt->spline, chi) / xckt->norm;
}

static void
_nc_xcor_kernel_table_get_comp_support (NcXcorKernelRadial *xcka, guint comp, gdouble *chi_min, gdouble *chi_max)
{
  NcXcorKernelTable *xckt = NC_XCOR_KERNEL_TABLE (xcka);

  *chi_min = xckt->chi_min;
  *chi_max = xckt->chi_max;
}

static gdouble
_nc_xcor_kernel_table_eval_kernel_factor (NcXcorKernelRadial *xcka, NcHICosmo *cosmo, gdouble chi, gdouble k)
{
  NcXcorKernelTable *xckt = NC_XCOR_KERNEL_TABLE (xcka);

  if (xckt->kind == NC_XCOR_KERNEL_TABLE_KIND_SHEAR)
  {
    const gdouble kchi = k * chi;

    return 1.0 / (kchi * kchi);
  }

  return 1.0;
}

static gdouble
_nc_xcor_kernel_table_eval_prefactor (NcXcorKernelRadial *xcka, NcHICosmo *cosmo, gint l)
{
  NcXcorKernelTable *xckt = NC_XCOR_KERNEL_TABLE (xcka);

  if (xckt->kind == NC_XCOR_KERNEL_TABLE_KIND_SHEAR)
  {
    const gdouble v = (l + 2.0) * (l + 1.0) * l * (l - 1.0);

    return (v > 0.0) ? sqrt (v) : 0.0;
  }

  return 1.0;
}

/**
 * nc_xcor_kernel_table_new:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @chi: a #NcmVector of sample comoving distances in Mpc
 * @W: a #NcmVector of window samples
 *
 * Creates a new #NcXcorKernelTable of kind %NC_XCOR_KERNEL_TABLE_KIND_DENSITY,
 * reconstructed at the default order and normalized to unit integral.
 *
 * Returns: (transfer full): a new #NcXcorKernelTable
 */
NcXcorKernelTable *
nc_xcor_kernel_table_new (NcDistance *dist, NcmPowspec *ps, NcmVector *chi, NcmVector *W)
{
  NcXcorKernelTable *xckt = g_object_new (NC_TYPE_XCOR_KERNEL_TABLE,
                                          "dist", dist,
                                          "powspec", ps,
                                          "chi", chi,
                                          "W", W,
                                          NULL);

  return xckt;
}

/**
 * nc_xcor_kernel_table_new_full:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @chi: a #NcmVector of sample comoving distances in Mpc
 * @W: a #NcmVector of window samples
 * @kind: what the window is a window of
 * @order: B-spline order of the reconstruction, 2 to 10
 * @normalize: whether to rescale to unit integral over the support
 * @sbi: a #NcmSBesselIntegrator
 *
 * Creates a new #NcXcorKernelTable carrying @sbi, as nc_xcor_kernel_table_new()
 * does not. A #NcXcorKernel only accepts the non-Limber modes of
 * nc_xcor_kernel_set_l_limber() once it holds an integrator, so this is the
 * constructor to use for them.
 *
 * Returns: (transfer full): a new #NcXcorKernelTable
 */
NcXcorKernelTable *
nc_xcor_kernel_table_new_full (NcDistance *dist, NcmPowspec *ps, NcmVector *chi, NcmVector *W, NcXcorKernelTableKind kind, guint order, gboolean normalize, NcmSBesselIntegrator *sbi)
{
  NcXcorKernelTable *xckt = g_object_new (NC_TYPE_XCOR_KERNEL_TABLE,
                                          "dist", dist,
                                          "powspec", ps,
                                          "chi", chi,
                                          "W", W,
                                          "kind", kind,
                                          "order", order,
                                          "normalize", normalize,
                                          "integrator", sbi,
                                          NULL);

  return xckt;
}

/**
 * nc_xcor_kernel_table_get_kind:
 * @xckt: a #NcXcorKernelTable
 *
 * Returns: what the window is a window of.
 */
NcXcorKernelTableKind
nc_xcor_kernel_table_get_kind (NcXcorKernelTable *xckt)
{
  return xckt->kind;
}

/**
 * nc_xcor_kernel_table_get_order:
 * @xckt: a #NcXcorKernelTable
 *
 * Returns: the B-spline order of the reconstruction.
 */
guint
nc_xcor_kernel_table_get_order (NcXcorKernelTable *xckt)
{
  return xckt->order;
}

/**
 * nc_xcor_kernel_table_get_normalize:
 * @xckt: a #NcXcorKernelTable
 *
 * Returns: whether the window was rescaled to unit integral.
 */
gboolean
nc_xcor_kernel_table_get_normalize (NcXcorKernelTable *xckt)
{
  return xckt->normalize;
}

/**
 * nc_xcor_kernel_table_get_norm:
 * @xckt: a #NcXcorKernelTable
 *
 * Gets the divisor applied to the reconstruction, that is the integral of the
 * supplied table over its support, or 1 when #NcXcorKernelTable:normalize is
 * %FALSE.
 *
 * Returns: the normalization.
 */
gdouble
nc_xcor_kernel_table_get_norm (NcXcorKernelTable *xckt)
{
  return xckt->norm;
}

/**
 * nc_xcor_kernel_table_peek_spline:
 * @xckt: a #NcXcorKernelTable
 *
 * Gets the reconstruction itself, un-normalized.
 *
 * Returns: (transfer none): the #NcmSpline reconstructing the table.
 */
NcmSpline *
nc_xcor_kernel_table_peek_spline (NcXcorKernelTable *xckt)
{
  return xckt->spline;
}

/**
 * nc_xcor_kernel_table_peek_knots:
 * @xckt: a #NcXcorKernelTable
 *
 * Gets the reconstruction's breakpoints, that is the sample abscissae kept after
 * trimming. These are the natural panel edges for the radial quadrature: a panel
 * holding no interior breakpoint carries a forcing that is exactly a polynomial.
 *
 * Returns: (transfer none): the breakpoints, in Mpc.
 */
NcmVector *
nc_xcor_kernel_table_peek_knots (NcXcorKernelTable *xckt)
{
  return ncm_spline_peek_xv (xckt->spline);
}

