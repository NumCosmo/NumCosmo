/***************************************************************************
 *            nc_multiplicity_func_crocce.c
 *
 *  Wed Feb 15 13:36:09 2017
 *  Copyright  2017  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2017 <pennalima@gmail.com>
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
 * SECTION:nc_multiplicity_func_crocce
 * @title: NcMultiplicityFuncCrocce
 * @short_description: Dark matter halo -- Crocce multiplicity function.
 *
 * Dark matter halo multipliciticy function fitted for the MICE simulations. 
 * They used friends of friends algorithm, FoF(0.2). See reference arXiv:0907.0019.
 * 
 * $$f_{\textrm{MICE}} (\sigma, z) = A(z) \left[ \sigma^{-a(z)} + b(z) \right] e^{\left[ - \frac{c(z)}{\sigma^2}  \right]}$$,
 * where $A(z) = 0.58 (1+z)^{-0.13}$, $a(z) = 1.37(1+z)^{-0.15}$, $b(z) = 0.3(1+z)^{-0.084}$, 
 * and $c(z) = 1.036(1+z)^{-0.024}$. 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_multiplicity_func_crocce.h"

struct _NcMultiplicityFuncCroccePrivate
{
  NcMultiplicityFuncMassDef mdef;
  gdouble A0;
  gdouble a0;
  gdouble b0;
  gdouble c0;
};

enum
{
  PROP_0,
  PROP_A0,
  PROP_A1,
  PROP_B0,
  PROP_C0,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcMultiplicityFuncCrocce, nc_multiplicity_func_crocce, NC_TYPE_MULTIPLICITY_FUNC);

static void
nc_multiplicity_func_crocce_init (NcMultiplicityFuncCrocce *mc)
{
  NcMultiplicityFuncCroccePrivate * const self = mc->priv = nc_multiplicity_func_crocce_get_instance_private (mc);

  self->mdef    = NC_MULTIPLICITY_FUNC_MASS_DEF_LEN;
  self->A0 = 0.0;
  self->a0 = 0.0;
  self->b0 = 0.0;
  self->c0 = 0.0;
}

static void
_nc_multiplicity_func_crocce_set_property (GObject * object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcMultiplicityFuncCrocce *mc = NC_MULTIPLICITY_FUNC_CROCCE (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_CROCCE (object));

  switch (prop_id)
  {
    case PROP_A0:
      nc_multiplicity_func_crocce_set_A0 (mc, g_value_get_double (value));
      break;
    case PROP_A1:
      nc_multiplicity_func_crocce_set_a0 (mc, g_value_get_double (value));
      break;
    case PROP_B0:
      nc_multiplicity_func_crocce_set_b0 (mc, g_value_get_double (value));
      break;
    case PROP_C0:
      nc_multiplicity_func_crocce_set_c0 (mc, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_crocce_get_property (GObject * object, guint prop_id, GValue * value, GParamSpec * pspec)
{
  NcMultiplicityFuncCrocce *mc = NC_MULTIPLICITY_FUNC_CROCCE (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_CROCCE (object));

  switch (prop_id)
  {
    case PROP_A0:
      g_value_set_double (value, nc_multiplicity_func_crocce_get_A0 (mc));
      break;
    case PROP_A1:
      g_value_set_double (value, nc_multiplicity_func_crocce_get_a0 (mc));
      break;
    case PROP_B0:
      g_value_set_double (value, nc_multiplicity_func_crocce_get_b0 (mc));
      break;
    case PROP_C0:
      g_value_set_double (value, nc_multiplicity_func_crocce_get_c0 (mc));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_crocce_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_multiplicity_func_crocce_parent_class)->finalize (object);
}

static void _nc_multiplicity_func_crocce_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef); 
static NcMultiplicityFuncMassDef _nc_multiplicity_func_crocce_get_mdef (NcMultiplicityFunc *mulf);
static gdouble _nc_multiplicity_func_crocce_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z);

static void
nc_multiplicity_func_crocce_class_init (NcMultiplicityFuncCrocceClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass* parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);

  object_class->set_property = &_nc_multiplicity_func_crocce_set_property;
  object_class->get_property = &_nc_multiplicity_func_crocce_get_property;
  object_class->finalize     = &_nc_multiplicity_func_crocce_finalize;

  /**
   * NcMultiplicityFuncCrocce:A0:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_A0,
                                   g_param_spec_double ("A0",
                                                        NULL,
                                                        "A0",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.58,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcMultiplicityFuncCrocce:a0:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_A1,
                                   g_param_spec_double ("a0",
                                                        NULL,
                                                        "a0",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.37,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcMultiplicityFuncCrocce:b0:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_B0,
                                   g_param_spec_double ("b0",
                                                        NULL,
                                                        "b0",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcMultiplicityFuncCrocce:c0:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_C0,
                                   g_param_spec_double ("c0",
                                                        NULL,
                                                        "c0",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.036,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  parent_class->set_mdef = &_nc_multiplicity_func_crocce_set_mdef;
  parent_class->get_mdef = &_nc_multiplicity_func_crocce_get_mdef;
  parent_class->eval     = &_nc_multiplicity_func_crocce_eval;
}

static void 
_nc_multiplicity_func_crocce_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef)
{
  NcMultiplicityFuncCrocce *mc = NC_MULTIPLICITY_FUNC_CROCCE (mulf);
  NcMultiplicityFuncCroccePrivate * const self = mc->priv;

  switch (mdef)
  {
    case NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN:
      g_error ("NcMultiplicityFuncCrocce does not support critical mass def");
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL:
      g_error ("NcMultiplicityFuncCrocce does not support critical mass def");
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_VIRIAL:
      g_error ("NcMultiplicityFuncCrocce does not support virial mass def");
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_FOF:
      /* nothing to do */
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  self->mdef = mdef;
}

static NcMultiplicityFuncMassDef 
_nc_multiplicity_func_crocce_get_mdef (NcMultiplicityFunc *mulf)
{
  NcMultiplicityFuncCrocce *mc = NC_MULTIPLICITY_FUNC_CROCCE (mulf);
  NcMultiplicityFuncCroccePrivate * const self = mc->priv;

  return self->mdef;
}

static gdouble
_nc_multiplicity_func_crocce_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)   
{
  NcMultiplicityFuncCrocce *mc = NC_MULTIPLICITY_FUNC_CROCCE (mulf);
  NcMultiplicityFuncCroccePrivate * const self = mc->priv;
  
  const gdouble A = self->A0 * pow (1.0 + z, -0.13);
  const gdouble a = self->a0 * pow (1.0 + z, -0.15);
  const gdouble b = self->b0 * pow (1.0 + z, -0.084);
  const gdouble c = self->c0 * pow (1.0 + z, -0.024);    
  const gdouble f_crocce = A * (pow (sigma, -a) + b) * exp (- c / (sigma * sigma));     

  NCM_UNUSED (cosmo);
  
  return f_crocce;
}

/**
 * nc_multiplicity_func_crocce_new:
 *   
 * FIXME
 *
 * Returns: A new #NcMultiplicityFuncCrocce.
 */
NcMultiplicityFuncCrocce *
nc_multiplicity_func_crocce_new (void)
{
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_CROCCE,
                       "mass-def", NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN,
                       NULL);
}

/**
 * nc_multiplicity_func_crocce_ref:
 * @mc: a #NcMultiplicityFuncCrocce
 *
 * Increases the reference count of @mc by one.
 *
 * Returns: (transfer full): @mc
 */
NcMultiplicityFuncCrocce *
nc_multiplicity_func_crocce_ref (NcMultiplicityFuncCrocce *mc)
{
  return g_object_ref (mc);
}

/**
 * nc_multiplicity_func_crocce_free:
 * @mc: a #NcMultiplicityFuncCrocce
 *
 * Atomically decrements the reference count of @mc by one. If the reference count drops to 0,
 * all memory allocated by @mc is released.
 *
 */
void
nc_multiplicity_func_crocce_free (NcMultiplicityFuncCrocce *mc)
{
  g_object_unref (mc);
}

/**
 * nc_multiplicity_func_crocce_clear:
 * @mc: a #NcMultiplicityFuncCrocce
 *
 * Atomically decrements the reference count of @mc by one. If the reference count drops to 0,
 * all memory allocated by @mc is released. Set the pointer to NULL;
 *
 */
void
nc_multiplicity_func_crocce_clear (NcMultiplicityFuncCrocce **mc)
{
  g_clear_object (mc);
}

/**
 * nc_multiplicity_func_crocce_set_A0:
 * @mc: a #NcMultiplicityFuncCrocce.
 * @A0: value of #NcMultiplicityFuncCrocce:A0.
 *
 * Sets the value @A0 to the #NcMultiplicityFuncCrocce:A0 property.
 *
 */
void
nc_multiplicity_func_crocce_set_A0 (NcMultiplicityFuncCrocce *mc, gdouble A0)
{
  NcMultiplicityFuncCroccePrivate * const self = mc->priv;

  g_assert (A0 >= 0);

  self->A0 = A0;
}

/**
 * nc_multiplicity_func_crocce_get_A0:
 * @mc: a #NcMultiplicityFuncCrocce.
 *
 * Returns: the value of #NcMultiplicityFuncCrocce:A0 property.
 */
gdouble
nc_multiplicity_func_crocce_get_A0 (const NcMultiplicityFuncCrocce *mc)
{
  NcMultiplicityFuncCroccePrivate * const self = mc->priv;

  return self->A0;
}

/**
 * nc_multiplicity_func_crocce_set_a0:
 * @mc: a #NcMultiplicityFuncCrocce.
 * @a0: value of #NcMultiplicityFuncCrocce:a0.
 *
 * Sets the value @a0 to the #NcMultiplicityFuncCrocce:a0 property.
 *
 */
void
nc_multiplicity_func_crocce_set_a0 (NcMultiplicityFuncCrocce *mc, gdouble a0)
{
  NcMultiplicityFuncCroccePrivate * const self = mc->priv;

  g_assert (a0 >= 0);

  self->a0 = a0;
}

/**
 * nc_multiplicity_func_crocce_get_a0:
 * @mc: a #NcMultiplicityFuncCrocce.
 *
 * Returns: the value of #NcMultiplicityFuncCrocce:a0 property.
 */
gdouble
nc_multiplicity_func_crocce_get_a0 (const NcMultiplicityFuncCrocce *mc)
{
  NcMultiplicityFuncCroccePrivate * const self = mc->priv;

  return self->a0;
}

/**
 * nc_multiplicity_func_crocce_set_b0:
 * @mc: a #NcMultiplicityFuncCrocce.
 * @b0: value of #NcMultiplicityFuncCrocce:b0.
 *
 * Sets the value @b0 to the #NcMultiplicityFuncCrocce:b0 property.
 *
 */
void
nc_multiplicity_func_crocce_set_b0 (NcMultiplicityFuncCrocce *mc, gdouble b0)
{
  NcMultiplicityFuncCroccePrivate * const self = mc->priv;

  g_assert (b0 >= 0);

  self->b0 = b0;
}

/**
 * nc_multiplicity_func_crocce_get_b0:
 * @mc: a #NcMultiplicityFuncCrocce.
 *
 * Returns: the value of #NcMultiplicityFuncCrocce:b0 property.
 */
gdouble
nc_multiplicity_func_crocce_get_b0 (const NcMultiplicityFuncCrocce *mc)
{
  NcMultiplicityFuncCroccePrivate * const self = mc->priv;

  return self->b0;
}

/**
 * nc_multiplicity_func_crocce_set_c0:
 * @mc: a #NcMultiplicityFuncCrocce.
 * @c0: value of #NcMultiplicityFuncCrocce:c0.
 *
 * Sets the value @c0 to the #NcMultiplicityFuncCrocce:c0 property.
 *
 */
void
nc_multiplicity_func_crocce_set_c0 (NcMultiplicityFuncCrocce *mc, gdouble c0)
{
  NcMultiplicityFuncCroccePrivate * const self = mc->priv;

  g_assert (c0 >= 0);

  self->c0 = c0;
}

/**
 * nc_multiplicity_func_crocce_get_c0:
 * @mc: a #NcMultiplicityFuncCrocce.
 *
 * Returns: the value of #NcMultiplicityFuncCrocce:c0 property.
 */
gdouble
nc_multiplicity_func_crocce_get_c0 (const NcMultiplicityFuncCrocce *mc)
{
  NcMultiplicityFuncCroccePrivate * const self = mc->priv;

  return self->c0;
}

