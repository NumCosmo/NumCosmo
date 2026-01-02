/***************************************************************************
 *            nc_xcor_kernel.c
 *
 *  Tue July 14 12:00:00 2015
 *  Copyright  2015  Cyrille Doux
 *  <cdoux@apc.in2p3.fr>
 *  Sat December 27 20:21:01 2025
 *  Copyright  2025  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2015 Cyrille Doux <cdoux@apc.in2p3.fr>
 * Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcXcorKernel:
 *
 * Base object for the kernels of projected observables used in cross-correlations.
 *
 * The projected field and its kernel are linked by
 * \begin{equation}
 * $A(\hat{\mathbf{n}}) = \int_0^\infty dz \ W^A(z) \ \delta(\chi(z)\hat{\mathbf{n}}, z)$
 * \end{equation}
 * where $\delta$ is the matter density field.
 *
 * Kernels also implement the noise power spectrum.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_integrate.h"
#include "math/ncm_memory_pool.h"
#include "math/ncm_cfg.h"
#include "math/ncm_serialize.h"
#include "math/ncm_powspec.h"
#include "nc_distance.h"
#include "xcor/nc_xcor_kernel.h"
#include "xcor/nc_xcor.h"
#include "nc_enum_types.h"

typedef struct _NcXcorKernelPrivate
{
  /*< private >*/
  NcmModel parent_instance;
  gdouble cons_factor;
  gdouble zmin, zmax, zmid;

  NcDistance *dist;
  NcmPowspec *ps;

  gdouble (*eval_kernel_func) (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble k, const NcXcorKinetic *xck, gint l);
  void (*get_k_range_func) (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l, gdouble *kmin, gdouble *kmax);

  NcXcorKernelIntegMethod integ_method;
} NcXcorKernelPrivate;

enum
{
  PROP_0,
  PROP_ZMIN,
  PROP_ZMAX,
  PROP_DIST,
  PROP_POWSPEC,
  PROP_INTEG_METHOD,
  PROP_SIZE,
};


G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcXcorKernel, nc_xcor_kernel, NCM_TYPE_MODEL)
G_DEFINE_BOXED_TYPE (NcXcorKinetic, nc_xcor_kinetic, nc_xcor_kinetic_copy, nc_xcor_kinetic_free)

static void
nc_xcor_kernel_init (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  self->cons_factor      = 0.0;
  self->zmin             = 0.0;
  self->zmax             = 0.0;
  self->zmid             = 0.0;
  self->dist             = NULL;
  self->ps               = NULL;
  self->eval_kernel_func = NULL;
  self->get_k_range_func = NULL;
  self->integ_method     = NC_XCOR_KERNEL_INTEG_METHOD_LEN;
}

static void
_nc_xcor_kernel_dispose (GObject *object)
{
  NcXcorKernel *xclk        = NC_XCOR_KERNEL (object);
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  nc_distance_clear (&self->dist);
  ncm_powspec_clear (&self->ps);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_parent_class)->dispose (object);
}

static void
_nc_xcor_kernel_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_parent_class)->finalize (object);
}

static void
_nc_xcor_kernel_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_xcor_kernel_parent_class)->constructed (object);
  {
    NcXcorKernel *xclk        = NC_XCOR_KERNEL (object);
    NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

    if (self->eval_kernel_func == NULL)
      g_error ("nc_xcor_kernel_constructed: eval_kernel_func was not set by subclass implementation. "
               "Subclasses must call nc_xcor_kernel_set_eval_kernel_func() in their constructed method.");

    if (self->get_k_range_func == NULL)
      g_error ("nc_xcor_kernel_constructed: get_k_range_func was not set by subclass implementation. "
               "Subclasses must call nc_xcor_kernel_set_get_k_range_func() in their constructed method.");

    if (self->dist == NULL)
      g_error ("nc_xcor_kernel_constructed: dist property was not set. "
               "The 'dist' property must be provided at construction time.");

    if (self->ps == NULL)
      g_error ("nc_xcor_kernel_constructed: powspec property was not set. "
               "The 'powspec' property must be provided at construction time.");

    if (self->integ_method == NC_XCOR_KERNEL_INTEG_METHOD_LEN)
      g_error ("nc_xcor_kernel_constructed: integ-method property was not set. "
               "The 'integ-method' property must be provided at construction time.");

    nc_distance_compute_inv_comoving (self->dist, TRUE);
  }
}

static void
_nc_xcor_kernel_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorKernel *xclk        = NC_XCOR_KERNEL (object);
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  g_return_if_fail (NC_IS_XCOR_KERNEL (object));

  switch (prop_id)
  {
    case PROP_ZMIN:
      self->zmin = g_value_get_double (value);
      break;
    case PROP_ZMAX:
      self->zmax = g_value_get_double (value);
      break;
    case PROP_DIST:
      nc_distance_clear (&self->dist);
      self->dist = g_value_dup_object (value);
      break;
    case PROP_POWSPEC:
      ncm_powspec_clear (&self->ps);
      self->ps = g_value_dup_object (value);
      break;
    case PROP_INTEG_METHOD:
      self->integ_method = g_value_get_enum (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_xcor_kernel_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorKernel *xclk        = NC_XCOR_KERNEL (object);
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  g_return_if_fail (NC_IS_XCOR_KERNEL (object));

  switch (prop_id)
  {
    case PROP_ZMIN:
      g_value_set_double (value, self->zmin);
      break;
    case PROP_ZMAX:
      g_value_set_double (value, self->zmax);
      break;
    case PROP_DIST:
      g_value_set_object (value, self->dist);
      break;
    case PROP_POWSPEC:
      g_value_set_object (value, self->ps);
      break;
    case PROP_INTEG_METHOD:
      g_value_set_enum (value, self->integ_method);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

NCM_MSET_MODEL_REGISTER_ID (nc_xcor_kernel, NC_TYPE_XCOR_KERNEL);

static void
nc_xcor_kernel_class_init (NcXcorKernelClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_xcor_kernel_set_property;
  model_class->get_property = &_nc_xcor_kernel_get_property;
  object_class->constructed = &_nc_xcor_kernel_constructed;
  object_class->dispose     = &_nc_xcor_kernel_dispose;
  object_class->finalize    = &_nc_xcor_kernel_finalize;

  ncm_model_class_set_name_nick (model_class, "Cross-correlation Kernels", "xcor-kernel");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  ncm_model_class_check_params_info (NCM_MODEL_CLASS (klass));

  g_object_class_install_property (object_class,
                                   PROP_ZMIN,
                                   g_param_spec_double ("zmin",
                                                        NULL,
                                                        "Minimum redshift",
                                                        0.0, 1e5, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ZMAX,
                                   g_param_spec_double ("zmax",
                                                        NULL,
                                                        "Maximum redshift",
                                                        0.0, 1e5, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("dist",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_POWSPEC,
                                   g_param_spec_object ("powspec",
                                                        NULL,
                                                        "Power spectrum object",
                                                        NCM_TYPE_POWSPEC,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_INTEG_METHOD,
                                   g_param_spec_enum ("integ-method",
                                                      NULL,
                                                      "Integration method",
                                                      NC_TYPE_XCOR_KERNEL_INTEG_METHOD,
                                                      NC_XCOR_KERNEL_INTEG_METHOD_LIMBER,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  ncm_mset_model_register_id (model_class, "NcXcorKernel", "Cross-correlation Kernels",
                              NULL, TRUE, NCM_MSET_MODEL_MAIN);
}

/**
 * nc_xcor_kernel_ref:
 * @xclk: a #NcXcorKernel
 *
 * Increases the reference count of @xclk by one.
 *
 * Returns: (transfer full): @xclk.
 */
NcXcorKernel *
nc_xcor_kernel_ref (NcXcorKernel *xclk)
{
  return g_object_ref (xclk);
}

/**
 * nc_xcor_kernel_free:
 * @xclk: a #NcXcorKernel
 *
 * Decreases the reference count of @xclk by one. If the reference count
 * reaches zero, the object is freed.
 *
 */
void
nc_xcor_kernel_free (NcXcorKernel *xclk)
{
  g_object_unref (xclk);
}

/**
 * nc_xcor_kernel_clear:
 * @xclk: a #NcXcorKernel
 *
 * Atomically decrements the reference count of @xclk by one.
 * If the reference count drops to zero, all memory allocated by @xclk is
 * released. @xclk is set to NULL after being freed.
 *
 */
void
nc_xcor_kernel_clear (NcXcorKernel **xclk)
{
  g_clear_object (xclk);
}

/**
 * nc_xcor_kinetic_copy:
 * @xck: a #NcXcorKinetic
 *
 * Creates a copy of @xck.
 *
 * Returns: (transfer full): a new #NcXcorKinetic copy of @xck.
 */
NcXcorKinetic *
nc_xcor_kinetic_copy (NcXcorKinetic *xck)
{
  NcXcorKinetic *xck_copy = g_new (NcXcorKinetic, 1);

  xck_copy[0] = xck[0];

  return xck_copy;
}

/**
 * nc_xcor_kinetic_free:
 * @xck: a #NcXcorKinetic
 *
 * Frees @xck.
 *
 */
void
nc_xcor_kinetic_free (NcXcorKinetic *xck)
{
  g_free (xck);
}

/**
 * nc_xcor_kernel_obs_len: (virtual obs_len)
 * @xclk: a #NcXcorKernel
 *
 * Gets the number of observables required by this kernel.
 *
 * Returns: the number of observables
 */
guint
nc_xcor_kernel_obs_len (NcXcorKernel *xclk)
{
  return NC_XCOR_KERNEL_GET_CLASS (xclk)->obs_len (xclk);
}

/**
 * nc_xcor_kernel_obs_params_len: (virtual obs_params_len)
 * @xclk: a #NcXcorKernel
 *
 * Gets the number of parameters needed to describe the observables
 * for this kernel (e.g., measurement uncertainties, systematic parameters).
 *
 * Returns: the number of observable parameters
 */
guint
nc_xcor_kernel_obs_params_len (NcXcorKernel *xclk)
{
  return NC_XCOR_KERNEL_GET_CLASS (xclk)->obs_params_len (xclk);
}

/**
 * nc_xcor_kernel_set_z_range:
 * @xclk: a #NcXcorKernel
 * @zmin: minimum redshift
 * @zmax: maximum redshift
 * @zmid: mid redshift
 *
 * Set the redshift range of the kernel.
 *
 */
void
nc_xcor_kernel_set_z_range (NcXcorKernel *xclk, gdouble zmin, gdouble zmax, gdouble zmid)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  self->zmin = zmin;
  self->zmax = zmax;
  self->zmid = zmid;
}

/**
 * nc_xcor_kernel_get_z_range:
 * @xclk: a #NcXcorKernel
 * @zmin: (out): minimum redshift
 * @zmax: (out): maximum redshift
 * @zmid: (out) (allow-none): mid redshift
 *
 * Get the redshift range of the kernel.
 *
 */
void
nc_xcor_kernel_get_z_range (NcXcorKernel *xclk, gdouble *zmin, gdouble *zmax, gdouble *zmid)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  *zmin = self->zmin;
  *zmax = self->zmax;
  *zmid = self->zmid;
}

/**
 * nc_xcor_kernel_set_const_factor:
 * @xclk: a #NcXcorKernel
 * @cf: a #gdouble
 *
 * Set the constant factor of the kernel.
 *
 */
void
nc_xcor_kernel_set_const_factor (NcXcorKernel *xclk, gdouble cf)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  self->cons_factor = cf;
}

/**
 * nc_xcor_kernel_get_const_factor:
 * @xclk: a #NcXcorKernel
 *
 * Get the constant factor of the kernel.
 *
 * Returns: the constant factor.
 */
gdouble
nc_xcor_kernel_get_const_factor (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return self->cons_factor;
}

/**
 * nc_xcor_kernel_peek_dist:
 * @xclk: a #NcXcorKernel
 *
 * Peeks the distance object from the kernel. This method is intended
 * for use by subclass implementations.
 *
 * Returns: (transfer none): the distance object.
 */
NcDistance *
nc_xcor_kernel_peek_dist (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return self->dist;
}

/**
 * nc_xcor_kernel_peek_powspec:
 * @xclk: a #NcXcorKernel
 *
 * Peeks the power spectrum object from the kernel. This method is intended
 * for use by subclass implementations.
 *
 * Returns: (transfer none): the power spectrum object.
 */
NcmPowspec *
nc_xcor_kernel_peek_powspec (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return self->ps;
}

/**
 * nc_xcor_kernel_set_eval_kernel_func: (skip)
 * @xclk: a #NcXcorKernel
 * @eval_kernel_func: (scope notified): function pointer to evaluate the kernel
 *
 * Sets the function pointer that will be used to evaluate the kernel.
 * This method should only be called by subclass implementations during
 * the constructed phase to set the appropriate kernel evaluation function
 * based on the integration method.
 *
 */
void
nc_xcor_kernel_set_eval_kernel_func (NcXcorKernel *xclk, NcXcorKernelEvalFunc eval_kernel_func)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  self->eval_kernel_func = eval_kernel_func;
}

/**
 * nc_xcor_kernel_set_get_k_range_func: (skip)
 * @xclk: a #NcXcorKernel
 * @get_k_range_func: (scope notified): function pointer to get the k range
 *
 * Sets the function pointer that will be used to get the wavenumber range.
 * This method should only be called by subclass implementations during
 * the constructed phase to set the appropriate k range function
 * based on the integration method.
 *
 */
void
nc_xcor_kernel_set_get_k_range_func (NcXcorKernel *xclk, NcXcorKernelGetKRangeFunc get_k_range_func)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  self->get_k_range_func = get_k_range_func;
}

/**
 * nc_xcor_kernel_eval_kernel:
 * @xclk: a #NcXcorKernel
 * @cosmo: a #NcHICosmo
 * @k: wavenumber
 * @xck: a #NcXcorKinetic
 * @l: multipole
 *
 * Evaluates the kernel at wavenumber @k and multipole @l by calling the
 * function pointer set by subclasses. This method provides optimized
 * dispatch without virtual method overhead.
 *
 * Returns: the kernel evaluation result
 */
gdouble
nc_xcor_kernel_eval_kernel (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble k, const NcXcorKinetic *xck, gint l)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  g_assert (self->eval_kernel_func != NULL);

  return self->eval_kernel_func (xclk, cosmo, k, xck, l);
}

/**
 * nc_xcor_kernel_get_k_range:
 * @xclk: a #NcXcorKernel
 * @cosmo: a #NcHICosmo
 * @l: multipole
 * @kmin: (out): minimum wavenumber
 * @kmax: (out): maximum wavenumber
 *
 * Gets the wavenumber range for the kernel by calling the
 * function pointer set by subclasses. This method provides optimized
 * dispatch without virtual method overhead.
 *
 */
void
nc_xcor_kernel_get_k_range (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l, gdouble *kmin, gdouble *kmax)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  g_assert (self->get_k_range_func != NULL);

  self->get_k_range_func (xclk, cosmo, l, kmin, kmax);
}

/**
 * nc_xcor_kernel_get_integ_method:
 * @xclk: a #NcXcorKernel
 *
 * Gets the integration method for the kernel.
 *
 * Returns: the #NcXcorKernelIntegMethod
 */
NcXcorKernelIntegMethod
nc_xcor_kernel_get_integ_method (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return self->integ_method;
}

/**
 * nc_xcor_kernel_eval_radial_weight: (virtual eval_radial_weight)
 * @xclk: a #NcXcorKernel
 * @cosmo: a #NcHICosmo
 * @z: a #gdouble
 * @xck: a #NcXcorKinetic
 * @l: a #gint
 *
 * Evaluates the Limber kernel at redshift @z for multipole @l.
 * The kinetic quantities (comoving distance and Hubble parameter) are
 * provided in @xck. Returns zero if @z is outside the kernel's redshift range.
 *
 * Returns: the kernel value $W(z,\ell)$
 */
gdouble
nc_xcor_kernel_eval_radial_weight (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  if ((self->zmin <= z) && (self->zmax >= z))
    return NC_XCOR_KERNEL_GET_CLASS (xclk)->eval_radial_weight (xclk, cosmo, z, xck, l);
  else
    return 0.0;
}

/**
 * nc_xcor_kernel_eval_radial_weight_full:
 * @xclk: a #NcXcorKernel
 * @cosmo: a #NcHICosmo
 * @z: a #gdouble
 * @dist: a #NcDistance
 * @l: a #gint
 *
 * Evaluates the Limber kernel at redshift @z for multipole @l, including
 * the normalization factor. This function computes the kinetic quantities
 * internally using @dist and applies the kernel's constant factor.
 *
 * Returns: the normalized kernel value $c \times W(z,\ell)$
 */
gdouble
nc_xcor_kernel_eval_radial_weight_full (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, NcDistance *dist, gint l)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);
  const gdouble xi_z        = nc_distance_comoving (dist, cosmo, z); /* in units of Hubble radius */
  const gdouble E_z         = nc_hicosmo_E (cosmo, z);
  const NcXcorKinetic xck   = { xi_z, E_z };

  if ((self->zmin <= z) && (self->zmax >= z))
    return NC_XCOR_KERNEL_GET_CLASS (xclk)->eval_radial_weight (xclk, cosmo, z, &xck, l) * self->cons_factor;
  else
    return 0.0;
}

/**
 * nc_xcor_kernel_add_noise: (virtual add_noise)
 * @xclk: a #NcXcorKernel
 * @vp1: a #NcmVector
 * @vp2: a #NcmVector
 * @lmin: a #guint
 *
 * vp2 = vp1 + noise spectrum
 *
 */
void
nc_xcor_kernel_add_noise (NcXcorKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin)
{
  NC_XCOR_KERNEL_GET_CLASS (xclk)->add_noise (xclk, vp1, vp2, lmin);
}

/**
 * nc_xcor_kernel_prepare: (virtual prepare)
 * @xclk: a #NcXcorKernel
 * @cosmo: a NcHICosmo
 *
 * Prepares the kernel for evaluation with the given cosmological model.
 * This may involve precomputing quantities that depend on @cosmo but not
 * on redshift or multipole.
 *
 */
void
nc_xcor_kernel_prepare (NcXcorKernel *xclk, NcHICosmo *cosmo)
{
  return NC_XCOR_KERNEL_GET_CLASS (xclk)->prepare (xclk, cosmo);
}

static void
_nc_xcor_kernel_log_all_models_go (GType model_type, guint n)
{
  guint nc, i, j;
  GType *models = g_type_children (model_type, &nc);

  for (i = 0; i < nc; i++)
  {
    guint ncc;
    GType *modelsc = g_type_children (models[i], &ncc);

    g_message ("#  ");

    for (j = 0; j < n; j++)
      g_message (" ");

    g_message ("%s\n", g_type_name (models[i]));

    if (ncc)
      _nc_xcor_kernel_log_all_models_go (models[i], n + 2);

    g_free (modelsc);
  }

  g_free (models);
}

/**
 * nc_xcor_kernel_log_all_models:
 *
 * Logs all registered #NcXcorLimberKernel subclasses to the message log.
 * This is useful for debugging and discovering available kernel implementations.
 *
 */
void
nc_xcor_kernel_log_all_models (void)
{
  g_message ("# Registred NcXcorKernel:%s are:\n",
             g_type_name (NC_TYPE_XCOR_KERNEL));
  _nc_xcor_kernel_log_all_models_go (NC_TYPE_XCOR_KERNEL, 0);
}

