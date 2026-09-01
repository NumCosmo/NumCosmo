/***************************************************************************
 *            nc_xcor_kernel_radial.c
 *
 *  Thu August 21 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
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
 * NcXcorKernelRadial:
 *
 * Base class for kernels whose radial window is defined as a function of
 * comoving distance.
 *
 * Every other #NcXcorKernel builds its window from splined physical input, so
 * its $C_\ell$ has no independently known value. A subclass of this one is
 * defined by a formula, so it is exactly evaluable outside NumCosmo and can
 * serve both as a test case with a known answer and as the shared input of a
 * cross-code benchmark.
 *
 * The window is a function of comoving distance in Mpc, normalized to unit
 * integral over its support,
 * \begin{equation}
 *   \int_{\chi_\mathrm{min}}^{\chi_\mathrm{max}} W(\chi)\,\mathrm{d}\chi = 1 ,
 * \end{equation}
 *
 * and the kernel entering the radial integral is
 * \begin{equation}
 *   K(\chi,k) = W(\chi)\sqrt{P(k, z=0)} .
 * \end{equation}
 *
 * The power spectrum is evaluated at $z = 0$, so all redshift dependence,
 * including growth, is carried by $W$. This is the convention of the N5K
 * challenge and of CCL's tracer kernels, and it keeps $z(\chi)$, itself a
 * spline, out of the exact path.
 *
 * Distances are in Mpc rather than the library's internal $R_H$, so the
 * resulting $C_\ell$ is directly comparable to a code working in Mpc with no
 * cosmology-dependent conversion factor.
 *
 * A subclass declares its components: the maximal intervals on which the
 * window is supported. Each becomes one #NcXcorKernelComponent, integrated
 * over exactly that interval and summed with the others, so no boundary is
 * ever interior to an integration domain.
 *
 * The split is by support, not by feature. A window with several bumps that
 * overlap into one continuous stretch is a single component, because it is
 * smooth throughout and cutting it would create edges rather than declare
 * them. Bumps whose supports are disjoint are separate components.
 *
 * A component must vanish outside the interval it reports. The truncation is
 * part of the definition rather than an implementation detail, since a
 * comparison against an external code is only meaningful if both truncate
 * identically.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/core/ncm_cfg.h"
#include "nc/xcor/nc_xcor_kernel_radial.h"
#include "nc/xcor/nc_xcor_kernel_component.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcXcorKernelRadialPrivate
{
  GPtrArray *comps;
  NcXcorKernelRadialKDep *kdep;
  gdouble z_min;
  gdouble z_max;
  gdouble z_mid;
} NcXcorKernelRadialPrivate;

enum
{
  PROP_0,
  PROP_SCALE_DEPENDENCE,
  PROP_SIZE,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcXcorKernelRadial, nc_xcor_kernel_radial, NC_TYPE_XCOR_KERNEL)

#define NC_XCOR_KERNEL_RADIAL_GET_PRIVATE(o) ((NcXcorKernelRadialPrivate *) nc_xcor_kernel_radial_get_instance_private ((NcXcorKernelRadial *) (o)))

/*
 * One per declared interval: W_comp(chi) sqrt(P(k, z=0)), with chi and k converted
 * from the library's internal units on entry and the result scaled by R_H so
 * that the internal radial integral reproduces the Mpc-normalized one.
 */

typedef struct _RadialComponentData
{
  NcXcorKernelRadial *xcka;
  NcmPowspec *ps;
  NcXcorKernelRadialKDep *kdep;
  guint comp;
} RadialComponentData;

#define _NC_XCOR_KERNEL_COMPONENT_RADIAL_GET_DATA(comp) \
        ((RadialComponentData *) ((guint8 *) (comp) + sizeof (NcXcorKernelComponent)))

static gdouble _radial_component_eval_kernel (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble xi, gdouble k);
static gdouble _radial_component_eval_prefactor (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble k, gint l);
static void _radial_component_get_limits (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble *xi_min, gdouble *xi_max, gdouble *k_min, gdouble *k_max);
static void _radial_component_data_clear (RadialComponentData *data);

NC_XCOR_KERNEL_COMPONENT_DEFINE_TYPE (NC, XCOR_KERNEL_COMPONENT_RADIAL,
                                      NcXcorKernelComponentRadial,
                                      nc_xcor_kernel_component_radial,
                                      _radial_component_eval_kernel,
                                      _radial_component_eval_prefactor,
                                      _radial_component_get_limits,
                                      RadialComponentData,
                                      _radial_component_data_clear)

static void
nc_xcor_kernel_radial_init (NcXcorKernelRadial *xcka)
{
  NcXcorKernelRadialPrivate * const self = NC_XCOR_KERNEL_RADIAL_GET_PRIVATE (xcka);

  self->comps = NULL;
  self->kdep  = NULL;
  self->z_min = GSL_NAN;
  self->z_max = GSL_NAN;
  self->z_mid = GSL_NAN;
}

static void
nc_xcor_kernel_radial_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_xcor_kernel_radial_parent_class)->constructed (object);
  {
    NcXcorKernelRadial *xcka               = NC_XCOR_KERNEL_RADIAL (object);
    NcXcorKernelRadialPrivate * const self = NC_XCOR_KERNEL_RADIAL_GET_PRIVATE (xcka);
    NcmPowspec *ps                         = nc_xcor_kernel_peek_powspec (NC_XCOR_KERNEL (xcka));
    const guint n_comps                    = nc_xcor_kernel_radial_get_n_comps (xcka);
    guint i;

    if ((n_comps == 0) || (n_comps > NC_XCOR_KERNEL_RADIAL_MAX_COMPS))
      g_error ("nc_xcor_kernel_radial_constructed: %s declares %u components, but a kernel "
               "carries between 1 and %d (NC_XCOR_KERNEL_RADIAL_MAX_COMPS).",
               G_OBJECT_TYPE_NAME (object), n_comps, NC_XCOR_KERNEL_RADIAL_MAX_COMPS);

    g_assert_null (self->comps);
    self->comps = g_ptr_array_new_full (n_comps, g_object_unref);

    for (i = 0; i < n_comps; i++)
    {
      NcXcorKernelComponent *comp = g_object_new (nc_xcor_kernel_component_radial_get_type (), NULL);
      RadialComponentData *data   = _NC_XCOR_KERNEL_COMPONENT_RADIAL_GET_DATA (comp);

      data->xcka = xcka;
      data->ps   = ps;
      data->kdep = self->kdep;
      data->comp = i;

      g_ptr_array_add (self->comps, comp);
    }
  }
}

static void
nc_xcor_kernel_radial_dispose (GObject *object)
{
  NcXcorKernelRadialPrivate * const self = NC_XCOR_KERNEL_RADIAL_GET_PRIVATE (object);

  g_clear_pointer (&self->comps, g_ptr_array_unref);
  nc_xcor_kernel_radial_kdep_clear (&self->kdep);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_radial_parent_class)->dispose (object);
}

static void
nc_xcor_kernel_radial_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_radial_parent_class)->finalize (object);
}

static void _nc_xcor_kernel_radial_get_z_range (NcXcorKernel *xclk, gdouble *zmin, gdouble *zmax, gdouble *zmid);
static gdouble _nc_xcor_kernel_radial_eval_limber_z (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l);
static gdouble _nc_xcor_kernel_radial_eval_limber_z_prefactor (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l);
static void _nc_xcor_kernel_radial_prepare (NcXcorKernel *xclk, NcHICosmo *cosmo);
static void _nc_xcor_kernel_radial_add_noise (NcXcorKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin);
static guint _nc_xcor_kernel_radial_obs_len (NcXcorKernel *xclk);
static guint _nc_xcor_kernel_radial_obs_params_len (NcXcorKernel *xclk);
static GPtrArray *_nc_xcor_kernel_radial_get_component_list (NcXcorKernel *xclk);
static gdouble _nc_xcor_kernel_radial_eval_kernel_factor_default (NcXcorKernelRadial *xcka, NcHICosmo *cosmo, gdouble chi, gdouble k);
static gdouble _nc_xcor_kernel_radial_eval_prefactor_default (NcXcorKernelRadial *xcka, NcHICosmo *cosmo, gint l);

static void
nc_xcor_kernel_radial_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorKernelRadialPrivate * const self = NC_XCOR_KERNEL_RADIAL_GET_PRIVATE (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_RADIAL (object));

  switch (prop_id)
  {
    case PROP_SCALE_DEPENDENCE:
      nc_xcor_kernel_radial_kdep_clear (&self->kdep);
      self->kdep = g_value_dup_object (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_xcor_kernel_radial_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorKernelRadialPrivate * const self = NC_XCOR_KERNEL_RADIAL_GET_PRIVATE (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_RADIAL (object));

  switch (prop_id)
  {
    case PROP_SCALE_DEPENDENCE:
      g_value_set_object (value, self->kdep);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_xcor_kernel_radial_class_init (NcXcorKernelRadialClass *klass)
{
  GObjectClass *object_class      = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class      = NCM_MODEL_CLASS (klass);
  NcXcorKernelClass *parent_class = NC_XCOR_KERNEL_CLASS (klass);

  object_class->constructed = &nc_xcor_kernel_radial_constructed;
  object_class->dispose     = &nc_xcor_kernel_radial_dispose;
  object_class->finalize    = &nc_xcor_kernel_radial_finalize;

  parent_class->get_z_range             = &_nc_xcor_kernel_radial_get_z_range;
  parent_class->eval_limber_z           = &_nc_xcor_kernel_radial_eval_limber_z;
  parent_class->eval_limber_z_prefactor = &_nc_xcor_kernel_radial_eval_limber_z_prefactor;
  parent_class->prepare                 = &_nc_xcor_kernel_radial_prepare;
  parent_class->add_noise               = &_nc_xcor_kernel_radial_add_noise;
  parent_class->obs_len                 = &_nc_xcor_kernel_radial_obs_len;
  parent_class->obs_params_len          = &_nc_xcor_kernel_radial_obs_params_len;
  parent_class->get_component_list      = &_nc_xcor_kernel_radial_get_component_list;

  klass->eval_kernel_factor = &_nc_xcor_kernel_radial_eval_kernel_factor_default;
  klass->eval_prefactor     = &_nc_xcor_kernel_radial_eval_prefactor_default;

  model_class->set_property = &nc_xcor_kernel_radial_set_property;
  model_class->get_property = &nc_xcor_kernel_radial_get_property;

  /* Installs the GObject-level property plumbing on this class. NcmModel routes
   * a property to the set_property of the class that *owns* its pspec, so the
   * one below is handled here and not by whichever subclass is instantiated. */
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  /**
   * NcXcorKernelRadial:scale-dependence:
   *
   * Optional #NcXcorKernelRadialKDep multiplying the radial integrand, or
   * %NULL for none. Any shape may carry any scale dependence, so the two are
   * varied independently rather than baked into separate kernel classes.
   */
  g_object_class_install_property (object_class,
                                   PROP_SCALE_DEPENDENCE,
                                   g_param_spec_object ("scale-dependence",
                                                        NULL,
                                                        "Scale-dependent factor multiplying the radial integrand",
                                                        NC_TYPE_XCOR_KERNEL_RADIAL_KDEP,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_xcor_kernel_radial_peek_kdep:
 * @xcka: a #NcXcorKernelRadial
 *
 * Gets the scale-dependent factor the kernel carries, if any.
 *
 * Returns: (transfer none) (nullable): the #NcXcorKernelRadialKDep, or %NULL.
 */
NcXcorKernelRadialKDep *
nc_xcor_kernel_radial_peek_kdep (NcXcorKernelRadial *xcka)
{
  return NC_XCOR_KERNEL_RADIAL_GET_PRIVATE (xcka)->kdep;
}

/**
 * nc_xcor_kernel_radial_get_n_comps: (virtual get_n_comps)
 * @xcka: a #NcXcorKernelRadial
 *
 * Gets the number of components the window is split into: its maximal intervals
 * of support. Each is integrated separately over exactly its own interval, so
 * that a boundary is never interior to an integration domain.
 *
 * Returns: the number of components, 1 to %NC_XCOR_KERNEL_RADIAL_MAX_COMPS.
 */
guint
nc_xcor_kernel_radial_get_n_comps (NcXcorKernelRadial *xcka)
{
  return NC_XCOR_KERNEL_RADIAL_GET_CLASS (xcka)->get_n_comps (xcka);
}

/**
 * nc_xcor_kernel_radial_eval_W_comp: (virtual eval_W_comp)
 * @xcka: a #NcXcorKernelRadial
 * @comp: component index
 * @chi: comoving distance $\chi$ in Mpc
 *
 * Evaluates component @comp of the radial window at @chi, zero outside the
 * interval that nc_xcor_kernel_radial_get_comp_support() reports for it. The
 * components sum to a window of unit integral.
 *
 * Returns: the component's contribution to $W(\chi)$, in $\mathrm{Mpc}^{-1}$.
 */
gdouble
nc_xcor_kernel_radial_eval_W_comp (NcXcorKernelRadial *xcka, guint comp, gdouble chi)
{
  return NC_XCOR_KERNEL_RADIAL_GET_CLASS (xcka)->eval_W_comp (xcka, comp, chi);
}

/**
 * nc_xcor_kernel_radial_get_comp_support: (virtual get_comp_support)
 * @xcka: a #NcXcorKernelRadial
 * @comp: component index
 * @chi_min: (out): lower end of the component's support, in Mpc
 * @chi_max: (out): upper end of the component's support, in Mpc
 *
 * Gets the interval outside which component @comp vanishes. It is exactly the
 * interval that component is integrated over, so the boundary is exact by
 * construction rather than something the quadrature has to locate.
 */
void
nc_xcor_kernel_radial_get_comp_support (NcXcorKernelRadial *xcka, guint comp, gdouble *chi_min, gdouble *chi_max)
{
  NC_XCOR_KERNEL_RADIAL_GET_CLASS (xcka)->get_comp_support (xcka, comp, chi_min, chi_max);
}

/**
 * nc_xcor_kernel_radial_eval_W:
 * @xcka: a #NcXcorKernelRadial
 * @chi: comoving distance $\chi$ in Mpc
 *
 * Evaluates the whole radial window, the sum of its components, normalized to unit
 * integral over its support and zero outside it.
 *
 * The cut-off is a step, not a taper: a truncated shape is finite at its own
 * edge and this returns exactly zero one ulp beyond it. Passing this function
 * directly to a quadrature over $[\chi_\mathrm{min}, \chi_\mathrm{max}]$ is
 * therefore unsafe, because a node placed just outside by rounding sees the
 * step and a spectral rule cannot fit it: #NcmSBesselIntegratorLevin aborts on
 * max-order instead of converging. Clamp $\chi$ into the support first, which
 * is what this class does internally on every integration path.
 *
 * Returns: the value of $W(\chi)$, in $\mathrm{Mpc}^{-1}$.
 */
gdouble
nc_xcor_kernel_radial_eval_W (NcXcorKernelRadial *xcka, gdouble chi)
{
  const guint n_comps = nc_xcor_kernel_radial_get_n_comps (xcka);
  gdouble W           = 0.0;
  guint i;

  for (i = 0; i < n_comps; i++)
    W += nc_xcor_kernel_radial_eval_W_comp (xcka, i, chi);

  return W;
}

/**
 * nc_xcor_kernel_radial_get_support:
 * @xcka: a #NcXcorKernelRadial
 * @chi_min: (out): lower end of the support, in Mpc
 * @chi_max: (out): upper end of the support, in Mpc
 *
 * Gets the hull of the components' intervals, outside which the whole window
 * vanishes. Components may be disjoint, so the hull can contain gaps where the
 * window is zero; it is what fixes the kernel's redshift range, not what
 * anything is integrated over.
 */
void
nc_xcor_kernel_radial_get_support (NcXcorKernelRadial *xcka, gdouble *chi_min, gdouble *chi_max)
{
  const guint n_comps = nc_xcor_kernel_radial_get_n_comps (xcka);
  guint i;

  nc_xcor_kernel_radial_get_comp_support (xcka, 0, chi_min, chi_max);

  for (i = 1; i < n_comps; i++)
  {
    gdouble lo, hi;

    nc_xcor_kernel_radial_get_comp_support (xcka, i, &lo, &hi);

    *chi_min = GSL_MIN (*chi_min, lo);
    *chi_max = GSL_MAX (*chi_max, hi);
  }
}

/*
 * Clamp roundoff outside a component's reported interval before evaluating it.
 *
 * A component is integrated over exactly the interval it reports, but the
 * caller reaches the ends of that interval indirectly -- the radial path
 * through y = k xi and back, the Limber path through z(chi) and back -- so the
 * argument can land just outside, where the component is exactly zero. That
 * step is a discontinuity the Levin panel's Chebyshev fit cannot resolve at any
 * order, and it aborts on the last panel of every k. Clamping onto the interval
 * changes no value the caller could have asked for.
 */
static gdouble
_nc_xcor_kernel_radial_eval_W_comp_clamped (NcXcorKernelRadial *xcka, guint comp, gdouble chi)
{
  gdouble chi_min, chi_max;

  nc_xcor_kernel_radial_get_comp_support (xcka, comp, &chi_min, &chi_max);

  return nc_xcor_kernel_radial_eval_W_comp (xcka, comp, GSL_MIN (GSL_MAX (chi, chi_min), chi_max));
}

/*
 * The z range follows from the support in comoving distance, which needs a
 * cosmology to convert; get_z_range() is not given one, so prepare() -- which
 * is -- caches it.
 */
static void
_nc_xcor_kernel_radial_get_z_range (NcXcorKernel *xclk, gdouble *zmin, gdouble *zmax, gdouble *zmid)
{
  NcXcorKernelRadialPrivate * const self = NC_XCOR_KERNEL_RADIAL_GET_PRIVATE (xclk);

  if (gsl_isnan (self->z_min))
    g_error ("_nc_xcor_kernel_radial_get_z_range: the redshift range of a radial kernel "
             "follows from its support in comoving distance and is only known once the kernel "
             "has been prepared; call nc_xcor_kernel_prepare() first.");

  *zmin = self->z_min;
  *zmax = self->z_max;
  *zmid = self->z_mid;
}

/*
 * The window per unit Hubble radius, times the two factors this class owes the
 * Limber integrand: its (chi, k) factor, and sqrt(P(k,0)/P(k,z)).
 *
 * The second exists because the conventions differ. This class carries the
 * growth in W and pairs it with P(k, 0) (see the class documentation), while
 * nc_xcor_limber_z.c multiplies the two kernels by P(k, z) -- the convention of
 * the physical kernels, whose W carries no growth. Each kernel returning
 * sqrt(P(k,0)/P(k,z)) makes the product come out as P(k, 0) without the shared
 * integrand having to know which kind of kernel it holds.
 */
static gdouble
_nc_xcor_kernel_radial_eval_limber_z (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l)
{
  NcXcorKernelRadial *xcka               = NC_XCOR_KERNEL_RADIAL (xclk);
  NcXcorKernelRadialPrivate * const self = NC_XCOR_KERNEL_RADIAL_GET_PRIVATE (xcka);
  const gdouble RH_Mpc                   = nc_hicosmo_RH_Mpc (cosmo);
  const gdouble chi                      = xck->xi_z * RH_Mpc;
  const guint n_comps                    = nc_xcor_kernel_radial_get_n_comps (xcka);
  gdouble W                              = 0.0;
  guint i;

  /* The Limber path has one integration domain for the whole kernel, so a
   * component is used only where chi falls in its own interval; outside, it
   * really is zero and must stay so. */
  for (i = 0; i < n_comps; i++)
  {
    gdouble lo, hi;

    nc_xcor_kernel_radial_get_comp_support (xcka, i, &lo, &hi);

    if ((chi >= lo) && (chi <= hi))
      W += nc_xcor_kernel_radial_eval_W_comp (xcka, i, chi);
  }

  {
    /* Limber fixes k = (l + 1/2) / chi, so the (chi, k) factor is evaluated
     * there rather than being left out of this branch. */
    const gdouble k = (l + 0.5) / chi;
    const gdouble f = nc_xcor_kernel_radial_eval_kernel_factor (xcka, cosmo, chi, k);

    /* This class carries the growth in W and pairs it with P(k, 0), but the
     * Limber integrand multiplies the two kernels by P(k, z). Each kernel hands
     * back sqrt(P(k,0)/P(k,z)) so the product restores P(k, 0). */
    NcmPowspec *ps       = nc_xcor_kernel_peek_powspec (xclk);
    const gdouble P_z    = ncm_powspec_eval (ps, NCM_MODEL (cosmo), z, k);
    const gdouble P_0    = ncm_powspec_eval (ps, NCM_MODEL (cosmo), 0.0, k);
    const gdouble growth = (P_z > 0.0) ? sqrt (P_0 / P_z) : 0.0;

    /* The scale dependence is part of the kernel, so it belongs on both
     * branches. Limber knows k, so there is nothing here it cannot evaluate. */
    const gdouble g = (self->kdep != NULL) ? nc_xcor_kernel_radial_kdep_eval (self->kdep, chi, k) : 1.0;

    return RH_Mpc * W * f * g * growth;
  }
}

static gdouble
_nc_xcor_kernel_radial_eval_limber_z_prefactor (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l)
{
  return nc_xcor_kernel_radial_eval_prefactor (NC_XCOR_KERNEL_RADIAL (xclk), cosmo, l);
}

static void
_nc_xcor_kernel_radial_prepare (NcXcorKernel *xclk, NcHICosmo *cosmo)
{
  NcXcorKernelRadialPrivate * const self = NC_XCOR_KERNEL_RADIAL_GET_PRIVATE (xclk);
  NcDistance *dist                       = nc_xcor_kernel_peek_dist (xclk);
  NcmPowspec *ps                         = nc_xcor_kernel_peek_powspec (xclk);
  const gdouble RH_Mpc                   = nc_hicosmo_RH_Mpc (cosmo);
  gdouble chi_min, chi_max;

  nc_distance_prepare_if_needed (dist, cosmo);
  ncm_powspec_prepare_if_needed (ps, NCM_MODEL (cosmo));

  nc_xcor_kernel_radial_get_support (NC_XCOR_KERNEL_RADIAL (xclk), &chi_min, &chi_max);

  self->z_min = nc_distance_inv_comoving (dist, cosmo, chi_min / RH_Mpc);
  self->z_max = nc_distance_inv_comoving (dist, cosmo, chi_max / RH_Mpc);
  self->z_mid = nc_distance_inv_comoving (dist, cosmo, 0.5 * (chi_min + chi_max) / RH_Mpc);

  {
    guint i;

    for (i = 0; i < self->comps->len; i++)
      nc_xcor_kernel_component_prepare (g_ptr_array_index (self->comps, i), cosmo);
  }
}

static void
_nc_xcor_kernel_radial_add_noise (NcXcorKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin)
{
  /* A radial kernel describes a window, not a survey: it carries no noise. */
}

static guint
_nc_xcor_kernel_radial_obs_len (NcXcorKernel *xclk)
{
  return 1;
}

static guint
_nc_xcor_kernel_radial_obs_params_len (NcXcorKernel *xclk)
{
  return 0;
}

static GPtrArray *
_nc_xcor_kernel_radial_get_component_list (NcXcorKernel *xclk)
{
  NcXcorKernelRadialPrivate * const self = NC_XCOR_KERNEL_RADIAL_GET_PRIVATE (xclk);
  GPtrArray *comp_list                   = g_ptr_array_new_full (self->comps->len, g_object_unref);
  guint i;

  for (i = 0; i < self->comps->len; i++)
    g_ptr_array_add (comp_list, nc_xcor_kernel_component_ref (g_ptr_array_index (self->comps, i)));

  return comp_list;
}

static gdouble
_nc_xcor_kernel_radial_eval_kernel_factor_default (NcXcorKernelRadial *xcka, NcHICosmo *cosmo, gdouble chi, gdouble k)
{
  return 1.0;
}

static gdouble
_nc_xcor_kernel_radial_eval_prefactor_default (NcXcorKernelRadial *xcka, NcHICosmo *cosmo, gint l)
{
  return 1.0;
}

/**
 * nc_xcor_kernel_radial_eval_kernel_factor: (virtual eval_kernel_factor)
 * @xcka: a #NcXcorKernelRadial
 * @cosmo: a #NcHICosmo
 * @chi: comoving distance $\chi$ in Mpc
 * @k: wave number in Mpc$^{-1}$
 *
 * Evaluates a multiplicative factor of $(\chi, k)$ applied to the kernel on top
 * of $W(\chi)\sqrt{P(k,0)}$. It must not depend on $\ell$: one closure serves a
 * whole block of multipoles, so anything $\ell$-dependent belongs in
 * nc_xcor_kernel_radial_eval_prefactor() instead.
 *
 * Shear is the motivating case, where the kernel carries $1/(k\chi)^2$.
 *
 * Defaults to 1.
 *
 * Returns: the factor at @chi and @k.
 */
gdouble
nc_xcor_kernel_radial_eval_kernel_factor (NcXcorKernelRadial *xcka, NcHICosmo *cosmo, gdouble chi, gdouble k)
{
  return NC_XCOR_KERNEL_RADIAL_GET_CLASS (xcka)->eval_kernel_factor (xcka, cosmo, chi, k);
}

/**
 * nc_xcor_kernel_radial_eval_prefactor: (virtual eval_prefactor)
 * @xcka: a #NcXcorKernelRadial
 * @cosmo: a #NcHICosmo
 * @l: multipole $\ell$
 *
 * Evaluates a multiplicative factor that depends only on $\ell$, applied to the
 * whole kernel. Shear contributes
 * $\sqrt{(\ell+2)(\ell+1)\ell(\ell-1)}$ here.
 *
 * Defaults to 1.
 *
 * Returns: the prefactor at @l.
 */
gdouble
nc_xcor_kernel_radial_eval_prefactor (NcXcorKernelRadial *xcka, NcHICosmo *cosmo, gint l)
{
  return NC_XCOR_KERNEL_RADIAL_GET_CLASS (xcka)->eval_prefactor (xcka, cosmo, l);
}

/*
 * Component implementation
 */

static void
_radial_component_data_clear (RadialComponentData *data)
{
  /* Weak references to the parent kernel and to what it already owns. */
}

static gdouble
_radial_component_eval_kernel (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble xi, gdouble k)
{
  RadialComponentData *data = _NC_XCOR_KERNEL_COMPONENT_RADIAL_GET_DATA (comp);
  const gdouble RH_Mpc      = nc_hicosmo_RH_Mpc (cosmo);
  const gdouble chi         = xi * RH_Mpc;
  const gdouble k_Mpc       = k / RH_Mpc;
  const gdouble W           = _nc_xcor_kernel_radial_eval_W_comp_clamped (data->xcka, data->comp, chi);
  const gdouble powspec     = ncm_powspec_eval (data->ps, NCM_MODEL (cosmo), 0.0, k_Mpc);
  const gdouble g           = (data->kdep != NULL) ? nc_xcor_kernel_radial_kdep_eval (data->kdep, chi, k_Mpc) : 1.0;
  const gdouble f           = nc_xcor_kernel_radial_eval_kernel_factor (data->xcka, cosmo, chi, k_Mpc);

  return RH_Mpc * W * g * f * sqrt (powspec);
}

static gdouble
_radial_component_eval_prefactor (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble k, gint l)
{
  RadialComponentData *data = _NC_XCOR_KERNEL_COMPONENT_RADIAL_GET_DATA (comp);

  return nc_xcor_kernel_radial_eval_prefactor (data->xcka, cosmo, l);
}

static void
_radial_component_get_limits (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble *xi_min, gdouble *xi_max, gdouble *k_min, gdouble *k_max)
{
  RadialComponentData *data = _NC_XCOR_KERNEL_COMPONENT_RADIAL_GET_DATA (comp);
  const gdouble RH_Mpc      = nc_hicosmo_RH_Mpc (cosmo);
  gdouble chi_min, chi_max;

  ncm_powspec_prepare_if_needed (data->ps, NCM_MODEL (cosmo));
  nc_xcor_kernel_radial_get_comp_support (data->xcka, data->comp, &chi_min, &chi_max);

  *xi_min = chi_min / RH_Mpc;
  *xi_max = chi_max / RH_Mpc;
  *k_min  = ncm_powspec_get_kmin (data->ps) * RH_Mpc;
  *k_max  = ncm_powspec_get_kmax (data->ps) * RH_Mpc;
}

