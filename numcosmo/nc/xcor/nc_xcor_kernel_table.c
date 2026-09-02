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
 * Radial window supplied as a table of $(\chi, W)$ samples, or as a list of
 * such tables.
 *
 * Every other #NcXcorKernelRadial is a formula. This one is the case a survey
 * pipeline actually hands you: the window already evaluated on a grid, with no
 * closed form behind it. Each table is one #NcXcorComponentTable, which
 * reconstructs its samples with a #NcmSplineBSpline and declares what the
 * window is a window of; this kernel integrates every component over its own
 * support and sums them, and is otherwise an ordinary radial kernel, so
 * #NcXcor and #NcXcorSolver drive it exactly as they drive the analytic ones.
 *
 * A tracer with several terms -- density plus redshift-space distortions plus
 * magnification for galaxy counts, shear plus intrinsic alignments for
 * lensing -- is one kernel with one component per term, each carrying its own
 * Bessel weight and $\ell$ prefactor through its #NcXcorComponentTable:kind.
 *
 * The single-table properties #NcXcorKernelTable:chi, #NcXcorKernelTable:W,
 * #NcXcorKernelTable:kind, #NcXcorKernelTable:order and
 * #NcXcorKernelTable:normalize build one component at construction; they are
 * the convenience for the one-term case and are consumed by it, so a kernel
 * always serializes through #NcXcorKernelTable:components.
 *
 * The convention is #NcXcorKernelRadial's: $\chi$ in Mpc, $P(k)$ at $z=0$, and
 * all redshift dependence -- growth included -- carried by $W$. A table
 * exported from CCL's tracer kernels or from the N5K challenge already
 * satisfies it.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/xcor/nc_xcor_kernel_table.h"
#include "ncm/spline/ncm_spline_bspline.h"
#include "ncm/core/ncm_obj_array.h"
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
  NcmObjArray *components;
};

enum
{
  PROP_0,
  PROP_CHI,
  PROP_W,
  PROP_KIND,
  PROP_ORDER,
  PROP_NORMALIZE,
  PROP_COMPONENTS,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcXcorKernelTable, nc_xcor_kernel_table, NC_TYPE_XCOR_KERNEL_RADIAL)

static void
nc_xcor_kernel_table_init (NcXcorKernelTable *xckt)
{
  xckt->chi        = NULL;
  xckt->W          = NULL;
  xckt->kind       = NC_XCOR_KERNEL_TABLE_KIND_DENSITY;
  xckt->order      = NCM_SPLINE_BSPLINE_DEFAULT_ORDER;
  xckt->normalize  = TRUE;
  xckt->components = NULL;
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
    case PROP_COMPONENTS:
      ncm_obj_array_clear (&xckt->components);
      xckt->components = g_value_dup_boxed (value);
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
      g_value_set_enum (value, nc_xcor_kernel_table_get_kind (xckt));
      break;
    case PROP_ORDER:
      g_value_set_uint (value, nc_xcor_kernel_table_get_order (xckt));
      break;
    case PROP_NORMALIZE:
      g_value_set_boolean (value, nc_xcor_kernel_table_get_normalize (xckt));
      break;
    case PROP_COMPONENTS:
      g_value_set_boxed (value, xckt->components);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

/*
 * Runs before the radial parent's constructed, which asks for the number of
 * components: the list must exist by then. The single-table properties are
 * consumed here, so the object always describes itself through the list.
 */
static void
nc_xcor_kernel_table_constructed (GObject *object)
{
  NcXcorKernelTable *xckt = NC_XCOR_KERNEL_TABLE (object);

  if (xckt->components == NULL)
  {
    NcXcorComponentTable *xcct;

    if ((xckt->chi == NULL) || (xckt->W == NULL))
      g_error ("nc_xcor_kernel_table_constructed: give either a list of components or both chi and W.");

    xcct             = nc_xcor_component_table_new_full (xckt->chi, xckt->W, xckt->kind, xckt->order, xckt->normalize);
    xckt->components = ncm_obj_array_new ();
    ncm_obj_array_add (xckt->components, G_OBJECT (xcct));
    nc_xcor_component_table_free (xcct);
  }
  else
  {
    guint i;

    if ((xckt->chi != NULL) || (xckt->W != NULL))
      g_error ("nc_xcor_kernel_table_constructed: give either a list of components or chi and W, not both.");

    if (ncm_obj_array_len (xckt->components) == 0)
      g_error ("nc_xcor_kernel_table_constructed: the list of components is empty.");

    for (i = 0; i < ncm_obj_array_len (xckt->components); i++)
      if (!NC_IS_XCOR_COMPONENT_TABLE (ncm_obj_array_peek (xckt->components, i)))
        g_error ("nc_xcor_kernel_table_constructed: component %u is a %s, not a NcXcorComponentTable.",
                 i, G_OBJECT_TYPE_NAME (ncm_obj_array_peek (xckt->components, i)));
  }

  ncm_vector_clear (&xckt->chi);
  ncm_vector_clear (&xckt->W);

  /* Chain up : end, the parent counts the components */
  G_OBJECT_CLASS (nc_xcor_kernel_table_parent_class)->constructed (object);
}

static void
nc_xcor_kernel_table_dispose (GObject *object)
{
  NcXcorKernelTable *xckt = NC_XCOR_KERNEL_TABLE (object);

  ncm_vector_clear (&xckt->chi);
  ncm_vector_clear (&xckt->W);
  ncm_obj_array_clear (&xckt->components);

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
static gdouble _nc_xcor_kernel_table_eval_kernel_factor (NcXcorKernelRadial *xcka, guint comp, NcHICosmo *cosmo, gdouble chi, gdouble k);
static gdouble _nc_xcor_kernel_table_eval_prefactor (NcXcorKernelRadial *xcka, guint comp, NcHICosmo *cosmo, gint l);
static guint _nc_xcor_kernel_table_get_comp_bessel_deriv (NcXcorKernelRadial *xcka, guint comp);

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
   * Comoving distances of the samples, in Mpc, strictly increasing. Builds a
   * single component with #NcXcorKernelTable:W at construction; %NULL once
   * the kernel exists.
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
   * Window samples $W(\chi)$, one per entry of #NcXcorKernelTable:chi. See
   * #NcXcorKernelTable:chi.
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
   * Kind of the single component built from #NcXcorKernelTable:chi and
   * #NcXcorKernelTable:W. Reads back the first component's kind.
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
   * B-spline order of the single component built from #NcXcorKernelTable:chi
   * and #NcXcorKernelTable:W; see #NcXcorComponentTable:order. Reads back the
   * first component's order.
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
   * Whether the single component built from #NcXcorKernelTable:chi and
   * #NcXcorKernelTable:W is rescaled to unit integral; see
   * #NcXcorComponentTable:normalize. Reads back the first component's flag.
   */
  g_object_class_install_property (object_class,
                                   PROP_NORMALIZE,
                                   g_param_spec_boolean ("normalize",
                                                         NULL,
                                                         "Rescale to unit integral over the support",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelTable:components:
   *
   * The tabulated components, a #NcmObjArray of #NcXcorComponentTable, one per
   * term of the tracer. Either this or the pair #NcXcorKernelTable:chi and
   * #NcXcorKernelTable:W is given at construction, not both.
   */
  g_object_class_install_property (object_class,
                                   PROP_COMPONENTS,
                                   g_param_spec_boxed ("components",
                                                       NULL,
                                                       "Tabulated components",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  parent_class->get_n_comps           = &_nc_xcor_kernel_table_get_n_comps;
  parent_class->eval_W_comp           = &_nc_xcor_kernel_table_eval_W_comp;
  parent_class->get_comp_support      = &_nc_xcor_kernel_table_get_comp_support;
  parent_class->eval_kernel_factor    = &_nc_xcor_kernel_table_eval_kernel_factor;
  parent_class->eval_prefactor        = &_nc_xcor_kernel_table_eval_prefactor;
  parent_class->get_comp_bessel_deriv = &_nc_xcor_kernel_table_get_comp_bessel_deriv;

  ncm_model_class_add_impl_flag (model_class, NC_XCOR_KERNEL_IMPL_ALL);
}

static NcXcorComponentTable *
_nc_xcor_kernel_table_comp (NcXcorKernelRadial *xcka, guint comp)
{
  NcXcorKernelTable *xckt = NC_XCOR_KERNEL_TABLE (xcka);

  return NC_XCOR_COMPONENT_TABLE (ncm_obj_array_peek (xckt->components, comp));
}

static guint
_nc_xcor_kernel_table_get_n_comps (NcXcorKernelRadial *xcka)
{
  NcXcorKernelTable *xckt = NC_XCOR_KERNEL_TABLE (xcka);

  return ncm_obj_array_len (xckt->components);
}

static gdouble
_nc_xcor_kernel_table_eval_W_comp (NcXcorKernelRadial *xcka, guint comp, gdouble chi)
{
  return nc_xcor_component_table_eval_W (_nc_xcor_kernel_table_comp (xcka, comp), chi);
}

static void
_nc_xcor_kernel_table_get_comp_support (NcXcorKernelRadial *xcka, guint comp, gdouble *chi_min, gdouble *chi_max)
{
  nc_xcor_component_table_get_support (_nc_xcor_kernel_table_comp (xcka, comp), chi_min, chi_max);
}

static gdouble
_nc_xcor_kernel_table_eval_kernel_factor (NcXcorKernelRadial *xcka, guint comp, NcHICosmo *cosmo, gdouble chi, gdouble k)
{
  return nc_xcor_component_table_eval_kernel_factor (_nc_xcor_kernel_table_comp (xcka, comp), chi, k);
}

static gdouble
_nc_xcor_kernel_table_eval_prefactor (NcXcorKernelRadial *xcka, guint comp, NcHICosmo *cosmo, gint l)
{
  return nc_xcor_component_table_eval_prefactor (_nc_xcor_kernel_table_comp (xcka, comp), l);
}

static guint
_nc_xcor_kernel_table_get_comp_bessel_deriv (NcXcorKernelRadial *xcka, guint comp)
{
  return nc_xcor_component_table_get_bessel_deriv (_nc_xcor_kernel_table_comp (xcka, comp));
}

/**
 * nc_xcor_kernel_table_new:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @chi: a #NcmVector of sample comoving distances in Mpc
 * @W: a #NcmVector of window samples
 *
 * Creates a one-component density kernel of default order, normalized to unit
 * integral, with no integrator set.
 *
 * Returns: (transfer full): a new #NcXcorKernelTable
 */
NcXcorKernelTable *
nc_xcor_kernel_table_new (NcDistance *dist, NcmPowspec *ps, NcmVector *chi, NcmVector *W)
{
  return g_object_new (NC_TYPE_XCOR_KERNEL_TABLE,
                       "dist", dist,
                       "powspec", ps,
                       "chi", chi,
                       "W", W,
                       NULL);
}

/**
 * nc_xcor_kernel_table_new_full:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @chi: a #NcmVector of sample comoving distances in Mpc
 * @W: a #NcmVector of window samples
 * @kind: a #NcXcorKernelTableKind
 * @order: B-spline order of the reconstruction
 * @normalize: whether to rescale the window to unit integral
 * @sbi: (nullable): a #NcmSBesselIntegrator
 *
 * Creates a one-component kernel with every option set.
 *
 * Returns: (transfer full): a new #NcXcorKernelTable
 */
NcXcorKernelTable *
nc_xcor_kernel_table_new_full (NcDistance *dist, NcmPowspec *ps, NcmVector *chi, NcmVector *W, NcXcorKernelTableKind kind, guint order, gboolean normalize, NcmSBesselIntegrator *sbi)
{
  return g_object_new (NC_TYPE_XCOR_KERNEL_TABLE,
                       "dist", dist,
                       "powspec", ps,
                       "chi", chi,
                       "W", W,
                       "kind", kind,
                       "order", order,
                       "normalize", normalize,
                       "integrator", sbi,
                       NULL);
}

/**
 * nc_xcor_kernel_table_new_from_components:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @components: a #NcmObjArray of #NcXcorComponentTable
 * @sbi: (nullable): a #NcmSBesselIntegrator
 *
 * Creates a kernel from a list of tabulated components, one per term of the
 * tracer.
 *
 * Returns: (transfer full): a new #NcXcorKernelTable
 */
NcXcorKernelTable *
nc_xcor_kernel_table_new_from_components (NcDistance *dist, NcmPowspec *ps, NcmObjArray *components, NcmSBesselIntegrator *sbi)
{
  return g_object_new (NC_TYPE_XCOR_KERNEL_TABLE,
                       "dist", dist,
                       "powspec", ps,
                       "components", components,
                       "integrator", sbi,
                       NULL);
}

/**
 * nc_xcor_kernel_table_get_n_components:
 * @xckt: a #NcXcorKernelTable
 *
 * Returns: the number of tabulated components.
 */
guint
nc_xcor_kernel_table_get_n_components (NcXcorKernelTable *xckt)
{
  return ncm_obj_array_len (xckt->components);
}

/**
 * nc_xcor_kernel_table_peek_component:
 * @xckt: a #NcXcorKernelTable
 * @i: component index
 *
 * Returns: (transfer none): the @i-th #NcXcorComponentTable
 */
NcXcorComponentTable *
nc_xcor_kernel_table_peek_component (NcXcorKernelTable *xckt, guint i)
{
  g_assert_cmpuint (i, <, ncm_obj_array_len (xckt->components));

  return NC_XCOR_COMPONENT_TABLE (ncm_obj_array_peek (xckt->components, i));
}

/**
 * nc_xcor_kernel_table_replace_samples:
 * @xckt: a #NcXcorKernelTable
 * @i: component index
 * @chi: a #NcmVector of sample comoving distances in Mpc
 * @W: a #NcmVector of window samples
 *
 * Replaces the samples of the @i-th component in place (see
 * nc_xcor_component_table_set_samples()) and marks the kernel outdated, so its
 * next nc_xcor_kernel_prepare_if_needed() prepares again. The kernel object,
 * and with it its registration in a #NcXcorSolver and the solver's per-block
 * integrators, is unchanged: this is the per-step update of a sampler whose
 * windows come from a tabulating code.
 */
void
nc_xcor_kernel_table_replace_samples (NcXcorKernelTable *xckt, guint i, NcmVector *chi, NcmVector *W)
{
  nc_xcor_component_table_set_samples (nc_xcor_kernel_table_peek_component (xckt, i), chi, W);
  nc_xcor_kernel_mark_outdated (NC_XCOR_KERNEL (xckt));
}

/**
 * nc_xcor_kernel_table_get_kind:
 * @xckt: a #NcXcorKernelTable
 *
 * Returns: the kind of the first component.
 */
NcXcorKernelTableKind
nc_xcor_kernel_table_get_kind (NcXcorKernelTable *xckt)
{
  return nc_xcor_component_table_get_kind (nc_xcor_kernel_table_peek_component (xckt, 0));
}

/**
 * nc_xcor_kernel_table_get_order:
 * @xckt: a #NcXcorKernelTable
 *
 * Returns: the B-spline order of the first component.
 */
guint
nc_xcor_kernel_table_get_order (NcXcorKernelTable *xckt)
{
  return nc_xcor_component_table_get_order (nc_xcor_kernel_table_peek_component (xckt, 0));
}

/**
 * nc_xcor_kernel_table_get_normalize:
 * @xckt: a #NcXcorKernelTable
 *
 * Returns: whether the first component is rescaled to unit integral.
 */
gboolean
nc_xcor_kernel_table_get_normalize (NcXcorKernelTable *xckt)
{
  return nc_xcor_component_table_get_normalize (nc_xcor_kernel_table_peek_component (xckt, 0));
}

/**
 * nc_xcor_kernel_table_get_norm:
 * @xckt: a #NcXcorKernelTable
 *
 * Returns: the normalization constant of the first component; see
 * nc_xcor_component_table_get_norm().
 */
gdouble
nc_xcor_kernel_table_get_norm (NcXcorKernelTable *xckt)
{
  return nc_xcor_component_table_get_norm (nc_xcor_kernel_table_peek_component (xckt, 0));
}

/**
 * nc_xcor_kernel_table_peek_spline:
 * @xckt: a #NcXcorKernelTable
 *
 * Returns: (transfer none): the reconstruction of the first component's samples
 */
NcmSpline *
nc_xcor_kernel_table_peek_spline (NcXcorKernelTable *xckt)
{
  return nc_xcor_component_table_peek_spline (nc_xcor_kernel_table_peek_component (xckt, 0));
}

/**
 * nc_xcor_kernel_table_peek_knots:
 * @xckt: a #NcXcorKernelTable
 *
 * Returns: (transfer none): the first component's knots, in Mpc
 */
NcmVector *
nc_xcor_kernel_table_peek_knots (NcXcorKernelTable *xckt)
{
  return nc_xcor_component_table_peek_knots (nc_xcor_kernel_table_peek_component (xckt, 0));
}

