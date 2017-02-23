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

G_DEFINE_TYPE (NcMultiplicityFuncCrocce, nc_multiplicity_func_crocce, NC_TYPE_MULTIPLICITY_FUNC);

enum
{
  PROP_0,
  PROP_A0,
  PROP_A1,
  PROP_B0,
  PROP_C0,
};

/**
 * nc_multiplicity_func_crocce_new:
 * @A0: FIXME
 * @a0: FIXME 
 * @b0: FIXME
 * @c0: FIXME 
 *   
 * FIXME
 *
 * Returns: A new #NcMultiplicityFunc.
 */
NcMultiplicityFunc *
nc_multiplicity_func_crocce_new (gdouble A0, gdouble a0, gdouble b0, gdouble c0)
{
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_CROCCE,
                       "A0", A0,
                       "a0", a0,
                       "b0", b0, 
                       "c0", c0,
                       NULL);
}

static gdouble
_nc_multiplicity_func_crocce_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)   
{
  NcMultiplicityFuncCrocce *mulf_crocce = NC_MULTIPLICITY_FUNC_CROCCE (mulf);
  const gdouble A = mulf_crocce->A0 * pow(1.0 + z, -0.13);
  const gdouble a = mulf_crocce->a0 * pow(1.0 + z, -0.15);
  const gdouble b = mulf_crocce->b0 * pow(1.0 + z, -0.084);
  const gdouble c = mulf_crocce->c0 * pow(1.0 + z, -0.024);	
  const gdouble f_crocce = A * (pow(sigma, -a) + b) * exp(-c / (sigma * sigma));	

  NCM_UNUSED (cosmo);
  
  //	printf ("%.15g %.15g %.15g %.15g | %.15g %.15g %.15g %.15g\n", mulf_crocce->A0, mulf_crocce->a0, mulf_crocce->b0, mulf_crocce->c0,
  //	        A, a, b, c);
  return f_crocce;
}

/**
 * nc_multiplicity_func_crocce_set_A0:
 * @mulf_crocce: a #NcMultiplicityFuncCrocce.
 * @A0: value of #NcMultiplicityFuncCrocce:A0.
 *
 * Sets the value @A0 to the #NcMultiplicityFuncCrocce:A0 property.
 *
 */
void
nc_multiplicity_func_crocce_set_A0 (NcMultiplicityFuncCrocce *mulf_crocce, gdouble A0)
{
  g_assert (A0 >= 0);
  mulf_crocce->A0 = A0;
}

/**
 * nc_multiplicity_func_crocce_get_A0:
 * @mulf_crocce: a #NcMultiplicityFuncCrocce.
 *
 * Returns: the value of #NcMultiplicityFuncCrocce:A0 property.
 */
gdouble
nc_multiplicity_func_crocce_get_A0 (const NcMultiplicityFuncCrocce *mulf_crocce)
{
  return mulf_crocce->A0;
}

/**
 * nc_multiplicity_func_crocce_set_a0:
 * @mulf_crocce: a #NcMultiplicityFuncCrocce.
 * @a0: value of #NcMultiplicityFuncCrocce:a0.
 *
 * Sets the value @a0 to the #NcMultiplicityFuncCrocce:a0 property.
 *
 */
void
nc_multiplicity_func_crocce_set_a0 (NcMultiplicityFuncCrocce *mulf_crocce, gdouble a0)
{
  g_assert (a0 >= 0);
  mulf_crocce->a0 = a0;
}

/**
 * nc_multiplicity_func_crocce_get_a0:
 * @mulf_crocce: a #NcMultiplicityFuncCrocce.
 *
 * Returns: the value of #NcMultiplicityFuncCrocce:a0 property.
 */
gdouble
nc_multiplicity_func_crocce_get_a0 (const NcMultiplicityFuncCrocce *mulf_crocce)
{
  return mulf_crocce->a0;
}

/**
 * nc_multiplicity_func_crocce_set_b0:
 * @mulf_crocce: a #NcMultiplicityFuncCrocce.
 * @b0: value of #NcMultiplicityFuncCrocce:b0.
 *
 * Sets the value @b0 to the #NcMultiplicityFuncCrocce:b0 property.
 *
 */
void
nc_multiplicity_func_crocce_set_b0 (NcMultiplicityFuncCrocce *mulf_crocce, gdouble b0)
{
  g_assert (b0 >= 0);
  mulf_crocce->b0 = b0;
}

/**
 * nc_multiplicity_func_crocce_get_b0:
 * @mulf_crocce: a #NcMultiplicityFuncCrocce.
 *
 * Returns: the value of #NcMultiplicityFuncCrocce:b0 property.
 */
gdouble
nc_multiplicity_func_crocce_get_b0 (const NcMultiplicityFuncCrocce *mulf_crocce)
{
  return mulf_crocce->b0;
}

/**
 * nc_multiplicity_func_crocce_set_c0:
 * @mulf_crocce: a #NcMultiplicityFuncCrocce.
 * @c0: value of #NcMultiplicityFuncCrocce:c0.
 *
 * Sets the value @c0 to the #NcMultiplicityFuncCrocce:c0 property.
 *
 */
void
nc_multiplicity_func_crocce_set_c0 (NcMultiplicityFuncCrocce *mulf_crocce, gdouble c0)
{
  g_assert (c0 >= 0);
  mulf_crocce->c0 = c0;
}

/**
 * nc_multiplicity_func_crocce_get_c0:
 * @mulf_crocce: a #NcMultiplicityFuncCrocce.
 *
 * Returns: the value of #NcMultiplicityFuncCrocce:c0 property.
 */
gdouble
nc_multiplicity_func_crocce_get_c0 (const NcMultiplicityFuncCrocce *mulf_crocce)
{
  return mulf_crocce->c0;
}

static void
nc_multiplicity_func_crocce_init (NcMultiplicityFuncCrocce *mulf_crocce)
{
  /* TODO: Add initialization code here */
  mulf_crocce->A0 = 0.58;
  mulf_crocce->a0 = 1.37;
  mulf_crocce->b0 = 0.3;
  mulf_crocce->c0 = 1.036;
}

static void
_nc_multiplicity_func_crocce_finalize (GObject *object)
{
  /* TODO: Add deinitalization code here */

  G_OBJECT_CLASS (nc_multiplicity_func_crocce_parent_class)->finalize (object);
}

static void
_nc_multiplicity_func_crocce_set_property (GObject * object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcMultiplicityFuncCrocce *mulf_crocce = NC_MULTIPLICITY_FUNC_CROCCE (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_CROCCE (object));

  switch (prop_id)
  {
	case PROP_A0:
	  mulf_crocce->A0 = g_value_get_double (value);
	  break;
	case PROP_A1:
	  mulf_crocce->a0 = g_value_get_double (value);
	  break;
	case PROP_B0:
	  mulf_crocce->b0 = g_value_get_double (value);
	  break;
	case PROP_C0:
	  mulf_crocce->c0 = g_value_get_double (value);
	  break;
	default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_crocce_get_property (GObject * object, guint prop_id, GValue * value, GParamSpec * pspec)
{
  NcMultiplicityFuncCrocce *mulf_crocce = NC_MULTIPLICITY_FUNC_CROCCE (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_CROCCE (object));

  switch (prop_id)
  {
	case PROP_A0:
	  g_value_set_double (value, mulf_crocce->A0);
	  break;
	case PROP_A1:
	  g_value_set_double (value, mulf_crocce->a0);
	  break;
	case PROP_B0:
	  g_value_set_double (value, mulf_crocce->b0);
	  break;
	case PROP_C0:
	  g_value_set_double (value, mulf_crocce->c0);
	  break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_multiplicity_func_crocce_class_init (NcMultiplicityFuncCrocceClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass* parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);

  parent_class->eval = &_nc_multiplicity_func_crocce_eval;

  object_class->finalize = _nc_multiplicity_func_crocce_finalize;
  object_class->set_property = _nc_multiplicity_func_crocce_set_property;
  object_class->get_property = _nc_multiplicity_func_crocce_get_property;

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
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
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
 
}

