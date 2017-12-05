/***************************************************************************
 *            nc_xcor_limber_kernel.c
 *
 *  Tue July 14 12:00:00 2015
 *  Copyright  2015  Cyrille Doux
 *  <cdoux@apc.in2p3.fr>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2015 Cyrille Doux <cdoux@apc.in2p3.fr>
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
  * SECTION:nc_xcor_limber_kernel
  * @title: NcXcorLimberKernel
  * @short_description: Abstract object for the kernels of projected observables used in cross-correlations.
  *
  * The projected field and its kernel are linked by
  * \begin{equation}
  * $A(\hat{\mathbf{n}}) = \int_0^\infty dz \ W^A(z) \ \delta(\chi(z)\hat{\mathbf{n}}, z)$
  * \end{equation}
  * where $\delta$ is the matter density field.
  *
  * Kernels also implement the nosie power spectrum.
  */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/integral.h"
#include "math/memory_pool.h"
#include "math/ncm_cfg.h"
#include "math/ncm_serialize.h"
#include "xcor/nc_xcor_limber_kernel.h"
#include "xcor/nc_xcor.h"

G_DEFINE_ABSTRACT_TYPE (NcXcorLimberKernel, nc_xcor_limber_kernel, NCM_TYPE_MODEL);

enum
{
  PROP_0,
  PROP_ZMIN,
  PROP_ZMAX,
  PROP_SIZE,
};

static void
nc_xcor_limber_kernel_init (NcXcorLimberKernel* xclk)
{
  xclk->cons_factor = 0.0;
  xclk->zmin        = 0.0;
  xclk->zmax        = 0.0;
}

static void
_nc_xcor_limber_kernel_dispose (GObject* object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_limber_kernel_parent_class)->dispose (object);
}

static void
_nc_xcor_limber_kernel_finalize (GObject* object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_limber_kernel_parent_class)->finalize (object);
}

static void
_nc_xcor_limber_kernel_set_property (GObject* object, guint prop_id, const GValue* value, GParamSpec* pspec)
{
  NcXcorLimberKernel* xclk = NC_XCOR_LIMBER_KERNEL (object);
  g_return_if_fail (NC_IS_XCOR_LIMBER_KERNEL (object));

  switch (prop_id)
  {
    case PROP_ZMIN:
      xclk->zmin = g_value_get_double (value);
      break;
    case PROP_ZMAX:
      xclk->zmax = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_xcor_limber_kernel_get_property (GObject* object, guint prop_id, GValue* value, GParamSpec* pspec)
{
  NcXcorLimberKernel* xclk = NC_XCOR_LIMBER_KERNEL (object);
  g_return_if_fail (NC_IS_XCOR_LIMBER_KERNEL (object));

  switch (prop_id)
  {
    case PROP_ZMIN:
      g_value_set_double (value, xclk->zmin);
      break;
    case PROP_ZMAX:
      g_value_set_double (value, xclk->zmax);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

NCM_MSET_MODEL_REGISTER_ID (nc_xcor_limber_kernel, NC_TYPE_XCOR_LIMBER_KERNEL);

static void nc_xcor_limber_kernel_class_init (NcXcorLimberKernelClass* klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmModelClass* model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_xcor_limber_kernel_set_property;
  model_class->get_property = &_nc_xcor_limber_kernel_get_property;
  object_class->dispose     = &_nc_xcor_limber_kernel_dispose;
  object_class->finalize    = &_nc_xcor_limber_kernel_finalize;

  ncm_model_class_set_name_nick (model_class, "Cross-correlation Limber Kernels", "xcor-limber-kernel");
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

  ncm_mset_model_register_id (model_class, "NcXcorLimberKernel", "Cross-correlation Limber Kernels",
                              NULL, TRUE, NCM_MSET_MODEL_MAIN);
}

/**
 * nc_xcor_limber_kernel_new_from_name:
 * @xcor_name: string which specifies the type of the observable
 *
 * This function returns a new #NcXcorLimberKernel whose type is defined by
 * @xcor_name.
 *
 * Returns: (transfer full): A new #NcXcorLimberKernel.
   */
NcXcorLimberKernel *
nc_xcor_limber_kernel_new_from_name (gchar* xcor_name)
{
  GObject* obj = ncm_serialize_global_from_string (xcor_name);
  GType xcor_type = G_OBJECT_TYPE (obj);

  if (!g_type_is_a (xcor_type, NC_TYPE_XCOR_LIMBER_KERNEL))
    g_error ("nc_xcor_limber_kernel_new_from_name: NcXcorLimberKernel %s do not "
             "descend from %s.",
             xcor_name, g_type_name (NC_TYPE_XCOR_LIMBER_KERNEL));

  return NC_XCOR_LIMBER_KERNEL (obj);
}

/**
 * nc_xcor_limber_kernel_ref:
 * @xclk: a #NcXcorLimberKernel
 *
 * FIXME
 *
 * Returns: (transfer full): @xclk.
   */
NcXcorLimberKernel *
nc_xcor_limber_kernel_ref (NcXcorLimberKernel* xclk)
{
  return g_object_ref (xclk);
}

/**
 * nc_xcor_limber_kernel_free:
 * @xclk: a #NcXcorLimberKernel
 *
 * FIXME
 *
 */
void
nc_xcor_limber_kernel_free (NcXcorLimberKernel* xclk)
{
  g_object_unref (xclk);
}

/**
 * nc_xcor_limber_kernel_clear:
 * @xclk: a #NcXcorLimberKernel
 *
 * FIXME
 *
 */
void
nc_xcor_limber_kernel_clear (NcXcorLimberKernel **xclk)
{
  g_clear_object (xclk);
}

/**
 * nc_xcor_limber_kernel_obs_len: (virtual obs_len)
 * @xclk: a #NcXcorLimberKernel
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint nc_xcor_limber_kernel_obs_len (NcXcorLimberKernel *xclk)
{
  return NC_XCOR_LIMBER_KERNEL_GET_CLASS (xclk)->obs_len (xclk);
}

/**
 * nc_xcor_limber_kernel_obs_params_len: (virtual obs_params_len)
 * @xclk: a #NcXcorLimberKernel
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint
nc_xcor_limber_kernel_obs_params_len (NcXcorLimberKernel *xclk)
{
  return NC_XCOR_LIMBER_KERNEL_GET_CLASS (xclk)->obs_params_len (xclk);
}

/**
 * nc_xcor_limber_kernel_eval: (virtual eval)
 * @xclk: a #NcXcorLimberKernel
 * @cosmo: a #NcHICosmo
 * @z: a #gdouble
 * @xck: a #NcXcorKinetic
 * @l: a #gint
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_xcor_limber_kernel_eval (NcXcorLimberKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l)
{
  if ((xclk->zmin <= z) && (xclk->zmax >= z))
    return NC_XCOR_LIMBER_KERNEL_GET_CLASS (xclk)->eval (xclk, cosmo, z, xck, l);
  else
    return 0.0;
}

/**
 * nc_xcor_limber_kernel_eval_full:
 * @xclk: a #NcXcorLimberKernel
 * @cosmo: a #NcHICosmo
 * @z: a #gdouble
 * @dist: a #NcDistance
 * @l: a #gint
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_xcor_limber_kernel_eval_full (NcXcorLimberKernel *xclk, NcHICosmo *cosmo, gdouble z, NcDistance *dist, gint l)
{
  const gdouble xi_z = nc_distance_comoving (dist, cosmo, z); // in units of Hubble radius
  const gdouble E_z = nc_hicosmo_E (cosmo, z);
  const NcXcorKinetic xck = { xi_z, E_z };
  if ((xclk->zmin <= z) && (xclk->zmax >= z))
    return NC_XCOR_LIMBER_KERNEL_GET_CLASS (xclk)->eval (xclk, cosmo, z, &xck, l) * xclk->cons_factor;
  else
    return 0.0;
}

/**
 * nc_xcor_limber_kernel_add_noise: (virtual add_noise)
 * @xclk: a #NcXcorLimberKernel
 * @vp1: a #NcmVector
 * @vp2: a #NcmVector
 * @lmin: a #guint
 *
 * vp2 = vp1 + noise spectrum
 *
 */
void
nc_xcor_limber_kernel_add_noise (NcXcorLimberKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin)
{
  NC_XCOR_LIMBER_KERNEL_GET_CLASS (xclk)-> add_noise (xclk, vp1, vp2, lmin);
}

/**
 * nc_xcor_limber_kernel_prepare: (virtual prepare)
 * @xclk: a #NcXcorLimberKernel
 * @cosmo: a NcHICosmo
 *
 * FIXME
 *
 */
void
nc_xcor_limber_kernel_prepare (NcXcorLimberKernel *xclk, NcHICosmo *cosmo)
{
  return NC_XCOR_LIMBER_KERNEL_GET_CLASS (xclk)->prepare (xclk, cosmo);
}

static void
_nc_xcor_limber_kernel_log_all_models_go (GType model_type, guint n)
{
  guint nc, i, j;
  GType *models = g_type_children (model_type, &nc);
  for (i = 0; i < nc; i++)
  {
    guint ncc;
    GType* modelsc = g_type_children (models[i], &ncc);

    g_message ("#  ");
    for (j = 0; j < n; j++)
      g_message (" ");
    g_message ("%s\n", g_type_name (models[i]));
    if (ncc)
      _nc_xcor_limber_kernel_log_all_models_go (models[i], n + 2);

    g_free (modelsc);
  }
  g_free (models);
}

/**
 * nc_xcor_limber_kernel_log_all_models:
 *
 * FIXME
 *
 */
void
nc_xcor_limber_kernel_log_all_models (void)
{
  g_message ("# Registred NcXcorLimberKernel:%s are:\n",
             g_type_name (NC_TYPE_XCOR_LIMBER_KERNEL));
  _nc_xcor_limber_kernel_log_all_models_go (NC_TYPE_XCOR_LIMBER_KERNEL, 0);
}
