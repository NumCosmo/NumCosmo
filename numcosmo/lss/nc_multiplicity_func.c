/***************************************************************************
 *            nc_multiplicity_func.c
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/

/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
 * SECTION:nc_multiplicity_func
 * @title: NcMultiplicityFunc
 * @short_description: Dark matter halo multiplicity function.
 *
 * The  multiplicity function comprises information about the non-linear regime
 * of halo (structure) formation. The mass function can be written as
 * \begin{equation}\label{def:multip}
 * \frac{dn (M,z)}{d\ln M} = - \frac{\rho_m (z)}{M} f(\sigma_R, z) \frac{1}{\sigma_R} \frac{d\sigma_R}{d\ln M},
 * \end{equation}
 * where $\rho_m(z)$ is the mean matter density at redshift $z$, $f(\sigma_R, z)$ is the multiplicity function,
 * and $\sigma_R$ is the variance of the linear density contrast filtered on the length scale $R$ associated to the
 * mass $M$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_multiplicity_func.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"
#include "numcosmo/nc_enum_types.h"

struct _NcMultiplicityFuncPrivate
{
  gint place_holder;
};

enum
{
  PROP_0,
  PROP_MDEF,
  PROP_DELTA,
  PROP_SIZE,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcMultiplicityFunc, nc_multiplicity_func, G_TYPE_OBJECT);

static void
nc_multiplicity_func_init (NcMultiplicityFunc *mulf)
{
  NcMultiplicityFuncPrivate * const self = mulf->priv = nc_multiplicity_func_get_instance_private (mulf);

  self->place_holder = 0;
}

static void
_nc_multiplicity_func_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcMultiplicityFunc *mulf = NC_MULTIPLICITY_FUNC (object);

  /* NcMultiplicityFuncPrivate * const self = mulf->priv; */

  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC (object));

  switch (prop_id)
  {
    case PROP_MDEF:
      NC_MULTIPLICITY_FUNC_GET_CLASS (mulf)->set_mdef (mulf, g_value_get_enum (value));
      break;
    case PROP_DELTA:
      NC_MULTIPLICITY_FUNC_GET_CLASS (mulf)->set_Delta (mulf, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcMultiplicityFunc *mulf = NC_MULTIPLICITY_FUNC (object);

  /* NcMultiplicityFuncPrivate * const self = mulf->priv; */

  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC (object));

  switch (prop_id)
  {
    case PROP_MDEF:
      g_value_set_enum (value, nc_multiplicity_func_get_mdef (mulf));
      break;
    case PROP_DELTA:
      g_value_set_double (value, nc_multiplicity_func_get_Delta (mulf));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_multiplicity_func_parent_class)->finalize (object);
}

static void
_nc_multiplicity_func_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef)
{
  g_error ("method set_mdef not implemented by %s.", G_OBJECT_TYPE_NAME (mulf));
}

static void
_nc_multiplicity_func_set_Delta (NcMultiplicityFunc *mulf, gdouble Delta)
{
  g_error ("method set_Delta not implemented by %s.", G_OBJECT_TYPE_NAME (mulf));
}

static NcMultiplicityFuncMassDef
_nc_multiplicity_func_get_mdef (NcMultiplicityFunc *mulf)
{
  g_error ("method get_mdef not implemented by %s.", G_OBJECT_TYPE_NAME (mulf));

  return -1;
}

static gdouble
_nc_multiplicity_func_get_Delta (NcMultiplicityFunc *mulf)
{
  g_error ("method get_Delta not implemented by %s.", G_OBJECT_TYPE_NAME (mulf));

  return -1;
}

static gdouble
_nc_multiplicity_func_get_matter_Delta (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble z)
{
  g_error ("method get_matter_Delta not implemented by %s.", G_OBJECT_TYPE_NAME (mulf));

  return -1;
}

static gdouble
_nc_multiplicity_func_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)
{
  g_error ("method eval not implemented by %s.", G_OBJECT_TYPE_NAME (mulf));

  return 0.0;
}

static gboolean _nc_multiplicity_func_has_correction_factor (NcMultiplicityFunc *mulf);

static gdouble
_nc_multiplicity_func_correction_factor (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z, gdouble lnM)
{
  g_error ("method correction factor not implemented by %s.", G_OBJECT_TYPE_NAME (mulf));

  return 0.0;
}

static void
nc_multiplicity_func_class_init (NcMultiplicityFuncClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_multiplicity_func_set_property;
  object_class->get_property = &_nc_multiplicity_func_get_property;
  object_class->finalize     = &_nc_multiplicity_func_finalize;

  /**
   * NcMultiplicityFunc:mass-def:
   *
   * It refers to the halo finder used to obtain the multiplicity function (e.g., SO and FoF), and
   * the background density $\rho_\mathrm{bg}$ used in the mass definition \eqref{eq:mrr}.
   * See the enumerator #NcMultiplicityFuncMassDef for more details about the
   * background density definition.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_MDEF,
                                   g_param_spec_enum ("mass-def",
                                                      NULL,
                                                      "Mass definition",
                                                      NC_TYPE_MULTIPLICITY_FUNC_MASS_DEF, NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


/**
 * NcMultiplicityFunc:mass-Delta:
 *
 *
 *
 */
  g_object_class_install_property (object_class,
                                   PROP_DELTA,
                                   g_param_spec_double ("Delta",
                                                        NULL,
                                                        "Delta",
                                                        1.0, G_MAXDOUBLE, 200.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));



  klass->set_mdef              = &_nc_multiplicity_func_set_mdef;
  klass->set_Delta             = &_nc_multiplicity_func_set_Delta;
  klass->get_mdef              = &_nc_multiplicity_func_get_mdef;
  klass->get_Delta             = &_nc_multiplicity_func_get_Delta;
  klass->get_matter_Delta      = &_nc_multiplicity_func_get_matter_Delta;
  klass->eval                  = &_nc_multiplicity_func_eval;
  klass->has_correction_factor = &_nc_multiplicity_func_has_correction_factor;
  klass->correction_factor     = &_nc_multiplicity_func_correction_factor;
}

static gboolean
_nc_multiplicity_func_has_correction_factor (NcMultiplicityFunc *mulf)
{
  return FALSE;
}

/**
 * nc_multiplicity_func_free:
 * @mulf: a #NcMultiplicityFunc
 *
 * Atomically decrements the reference count of @mulf by one. If the reference count drops to 0,
 * all memory allocated by @mulf is released.
 *
 */
void
nc_multiplicity_func_free (NcMultiplicityFunc *mulf)
{
  g_object_unref (mulf);
}

/**
 * nc_multiplicity_func_clear:
 * @mulf: a #NcMultiplicityFunc
 *
 * Atomically decrements the reference count of @mulf by one. If the reference count drops to 0,
 * all memory allocated by @mulf is released. Set pointer to NULL.
 *
 */
void
nc_multiplicity_func_clear (NcMultiplicityFunc **mulf)
{
  g_clear_object (mulf);
}

/**
 * nc_multiplicity_func_set_mdef: (virtual set_mdef)
 * @mulf: a #NcMultiplicityFunc
 * @mdef: a #NcMultiplicityFuncMassDef
 *
 * Sets the mass definition to mdef.
 *
 */
void
nc_multiplicity_func_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef)
{
  NC_MULTIPLICITY_FUNC_GET_CLASS (mulf)->set_mdef (mulf, mdef);
}

/**
 * nc_multiplicity_func_set_Delta: (virtual set_Delta)
 * @mulf: a #NcMultiplicityFunc
 * @Delta: the $Delta$ value
 *
 * Sets the Delta for mean or matter mass definition to mdef.
 *
 */
void
nc_multiplicity_func_set_Delta (NcMultiplicityFunc *mulf, gdouble Delta)
{
  NcMultiplicityFuncMassDef mdef = nc_multiplicity_func_get_mdef (mulf);

  switch (mdef)
  {
    case NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN:
      NC_MULTIPLICITY_FUNC_GET_CLASS (mulf)->set_Delta (mulf, Delta);
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL:
      NC_MULTIPLICITY_FUNC_GET_CLASS (mulf)->set_Delta (mulf, Delta);
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_VIRIAL:
      g_error ("NcMultiplicityFuncMassDefVirial does not support Delta def");
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_FOF:
      g_error ("NcMultiplicityFuncMassDefFof does not support Delta def");
      break;
    default:
      g_assert_not_reached ();
      break;
  }
}

/**
 * nc_multiplicity_func_get_mdef: (virtual get_mdef)
 * @mulf: a #NcMultiplicityFunc
 *
 * Gets the mass definition.
 *
 * Returns: mdef.
 */
NcMultiplicityFuncMassDef
nc_multiplicity_func_get_mdef (NcMultiplicityFunc *mulf)
{
  return NC_MULTIPLICITY_FUNC_GET_CLASS (mulf)->get_mdef (mulf);
}

/**
 * nc_multiplicity_func_get_Delta: (virtual get_Delta)
 * @mulf: a #NcMultiplicityFunc
 *
 * Gets the Delta Value.
 *
 * Returns: Delta.
 */
gdouble
nc_multiplicity_func_get_Delta (NcMultiplicityFunc *mulf)
{
  return NC_MULTIPLICITY_FUNC_GET_CLASS (mulf)->get_Delta (mulf);
}

/**
 * nc_multiplicity_func_get_matter_Delta: (virtual get_matter_Delta)
 * @mulf: a #NcMultiplicityFunc
 * @cosmo: a #NcHICosmo
 * @z: the redshift $z$
 *
 * Gets the Delta Value with the matter density definition.
 *
 * Returns: Delta.
 */

gdouble
nc_multiplicity_func_get_matter_Delta (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble z)
{
  NcMultiplicityFuncMassDef mdef = nc_multiplicity_func_get_mdef (mulf);

  switch (mdef)
  {
    case NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN:

      return nc_multiplicity_func_get_Delta (mulf);

      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL:
    {
      const gdouble E2      = nc_hicosmo_E2 (cosmo, z);
      const gdouble Omega_m = nc_hicosmo_E2Omega_m (cosmo, z) / E2;
      const gdouble Delta_z = nc_multiplicity_func_get_Delta (mulf) / Omega_m;

      return Delta_z;

      break;
    }
    case NC_MULTIPLICITY_FUNC_MASS_DEF_VIRIAL:
      g_error ("NcMultiplicityFuncMassDefVirial does not support Delta def");
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_FOF:
      g_error ("NcMultiplicityFuncMassDefFof does not support Delta def");
      break;
    default:
      g_assert_not_reached ();
      break;
  }
  return 0.0;
}

/**
 * nc_multiplicity_func_eval: (virtual eval)
 * @mulf: a #NcMultiplicityFunc
 * @cosmo: a #NcHICosmo
 * @sigma: standard fluctuation of the matter density contrast
 * @z: redshift
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_multiplicity_func_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)
{
  return NC_MULTIPLICITY_FUNC_GET_CLASS (mulf)->eval (mulf, cosmo, sigma, z);
}

/**
 * nc_multiplicity_func_has_correction_factor: (virtual has_correction_factor)
 * @mulf: a #NcMultiplicityFunc
 *
 * Checks if the multiplicity function has a correction factor, e.g., NcMultiplicityFuncBocquet (when using
 * matter density).
 *
 * Returns: a gboolean
 */
gboolean
nc_multiplicity_func_has_correction_factor (NcMultiplicityFunc *mulf)
{
  return NC_MULTIPLICITY_FUNC_GET_CLASS (mulf)->has_correction_factor (mulf);
}

/**
 * nc_multiplicity_func_correction_factor: (virtual correction_factor)
 * @mulf: a #NcMultiplicityFunc
 * @cosmo: a #NcHICosmo
 * @sigma: standard fluctuation of the matter density contrast
 * @z: redshift
 * @lnM: logarithm base e of the mass
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_multiplicity_func_correction_factor (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z, gdouble lnM)
{
  return NC_MULTIPLICITY_FUNC_GET_CLASS (mulf)->correction_factor (mulf, cosmo, sigma, z, lnM);
}

