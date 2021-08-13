/***************************************************************************
 *            nc_multiplicity_func_warren.c
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
 * SECTION:nc_multiplicity_func_warren
 * @title: NcMultiplicityFuncWarren
 * @short_description: Dark matter halo -- Warren multiplicity function.
 *
 * FIXME
 * Reference: astro-ph/0506395 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_multiplicity_func_warren.h"

struct _NcMultiplicityFuncWarrenPrivate
{
  NcMultiplicityFuncMassDef mdef;
  gdouble A;
  gdouble a;
  gdouble b;
  gdouble c;
};

enum
{
  PROP_0,
  PROP_A,
  PROP_A1,
  PROP_B,
  PROP_C,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcMultiplicityFuncWarren, nc_multiplicity_func_warren, NC_TYPE_MULTIPLICITY_FUNC);

static void
nc_multiplicity_func_warren_init (NcMultiplicityFuncWarren *mw)
{
  NcMultiplicityFuncWarrenPrivate * const self = mw->priv = nc_multiplicity_func_warren_get_instance_private (mw);

  self->mdef    = NC_MULTIPLICITY_FUNC_MASS_DEF_LEN;
  self->A = 0.0;
  self->a = 0.0;
  self->b = 0.0;
  self->c = 0.0;
}

static void
_nc_multiplicity_func_warren_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcMultiplicityFuncWarren *mw = NC_MULTIPLICITY_FUNC_WARREN (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_WARREN (object));

  switch (prop_id)
  {
    case PROP_A:
      nc_multiplicity_func_warren_set_A (mw, g_value_get_double (value));
      break;
    case PROP_A1:
      nc_multiplicity_func_warren_set_a (mw, g_value_get_double (value));
      break;
    case PROP_B:
      nc_multiplicity_func_warren_set_b (mw, g_value_get_double (value));
      break;
    case PROP_C:
      nc_multiplicity_func_warren_set_c (mw, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_warren_get_property (GObject * object, guint prop_id, GValue * value, GParamSpec * pspec)
{
  NcMultiplicityFuncWarren *mw = NC_MULTIPLICITY_FUNC_WARREN (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_WARREN (object));

  switch (prop_id)
  {
    case PROP_A:
      g_value_set_double (value, nc_multiplicity_func_warren_get_A (mw));
      break;
    case PROP_A1:
      g_value_set_double (value, nc_multiplicity_func_warren_get_a (mw));
      break;
    case PROP_B:
      g_value_set_double (value, nc_multiplicity_func_warren_get_b (mw));
      break;
    case PROP_C:
      g_value_set_double (value, nc_multiplicity_func_warren_get_c (mw));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_warren_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_multiplicity_func_warren_parent_class)->finalize (object);
}

// _NC_MULTIPLICITY_FUNCTION_WARREN_DATASET_0506395 = {0.7234, 1.625, 0.2538, 1.1982};

static void _nc_multiplicity_func_warren_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef); 
static NcMultiplicityFuncMassDef _nc_multiplicity_func_warren_get_mdef (NcMultiplicityFunc *mulf);
static gdouble _nc_multiplicity_func_warren_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z);

static void
nc_multiplicity_func_warren_class_init (NcMultiplicityFuncWarrenClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass* parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);

  object_class->set_property = &_nc_multiplicity_func_warren_set_property;
  object_class->get_property = &_nc_multiplicity_func_warren_get_property;
  object_class->finalize     = &_nc_multiplicity_func_warren_finalize;
  
  /**
   * NcMultiplicityFuncWarren:A:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_A,
                                   g_param_spec_double ("A",
                                                        NULL,
                                                        "A",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.7234,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcMultiplicityFuncWarren:a:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_A1,
                                   g_param_spec_double ("a",
                                                        NULL,
                                                        "a",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.625,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcMultiplicityFuncWarren:b:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_B,
                                   g_param_spec_double ("b",
                                                        NULL,
                                                        "b",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.2538,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcMultiplicityFuncWarren:c:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_C,
                                   g_param_spec_double ("c",
                                                        NULL,
                                                        "c",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.1982,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  parent_class->set_mdef = &_nc_multiplicity_func_warren_set_mdef;
  parent_class->get_mdef = &_nc_multiplicity_func_warren_get_mdef;
  parent_class->eval     = &_nc_multiplicity_func_warren_eval;
}

static void 
_nc_multiplicity_func_warren_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef)
{
  NcMultiplicityFuncWarren *mw = NC_MULTIPLICITY_FUNC_WARREN (mulf);
  NcMultiplicityFuncWarrenPrivate * const self = mw->priv;

  switch (mdef)
  {
    case NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN:
      /* nothing to do */
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL:
      g_error ("NcMultiplicityFuncWarren does not support critical mass def");
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_VIRIAL:
      g_error ("NcMultiplicityFuncWarren does not support virial mass def");
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_FOF:
      g_error ("NcMultiplicityFuncWarren does not support fof mass def");
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  self->mdef = mdef;
}

static NcMultiplicityFuncMassDef 
_nc_multiplicity_func_warren_get_mdef (NcMultiplicityFunc *mulf)
{
  NcMultiplicityFuncWarren *mw = NC_MULTIPLICITY_FUNC_WARREN (mulf);
  NcMultiplicityFuncWarrenPrivate * const self = mw->priv;

  return self->mdef;
}

static gdouble
_nc_multiplicity_func_warren_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)   /* $f(\sigma)$ Warren: astro-ph/0506395 */
{
  NcMultiplicityFuncWarren *mw = NC_MULTIPLICITY_FUNC_WARREN (mulf);
  NcMultiplicityFuncWarrenPrivate * const self = mw->priv;
  
  gdouble f_Warren = self->A * (pow(sigma, - self->a) + self->b) * exp(-(self->c)/ (sigma * sigma) );

  NCM_UNUSED (cosmo);
  NCM_UNUSED (z);

  return f_Warren;
}

/**
 * nc_multiplicity_func_warren_new:
 *   
 * FIXME
 *
 * Returns: A new #NcMultiplicityFuncWarren.
 */
NcMultiplicityFuncWarren *
nc_multiplicity_func_warren_new (void)
{
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_WARREN,
                       "mass-def", NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN,
                       NULL);
}

/**
 * nc_multiplicity_func_warren_ref:
 * @mw: a #NcMultiplicityFuncWarren
 *
 * Increases the reference count of @mw by one.
 *
 * Returns: (transfer full): @mw
 */
NcMultiplicityFuncWarren *
nc_multiplicity_func_warren_ref (NcMultiplicityFuncWarren *mw)
{
  return g_object_ref (mw);
}

/**
 * nc_multiplicity_func_warren_free:
 * @mw: a #NcMultiplicityFuncWarren
 *
 * Atomically decrements the reference count of @mw by one. If the reference count drops to 0,
 * all memory allocated by @mw is released.
 *
 */
void
nc_multiplicity_func_warren_free (NcMultiplicityFuncWarren *mw)
{
  g_object_unref (mw);
}

/**
 * nc_multiplicity_func_warren_clear:
 * @mw: a #NcMultiplicityFuncWarren
 *
 * Atomically decrements the reference count of @mw by one. If the reference count drops to 0,
 * all memory allocated by @mw is released. Set the pointer to NULL;
 *
 */
void
nc_multiplicity_func_warren_clear (NcMultiplicityFuncWarren **mw)
{
  g_clear_object (mw);
}

/**
 * nc_multiplicity_func_warren_set_A:
 * @mw: a #NcMultiplicityFuncWarren.
 * @A: value of #NcMultiplicityFuncWarren:A.
 *
 * Sets the value @A to the #NcMultiplicityFuncWarren:A property.
 *
 */
void
nc_multiplicity_func_warren_set_A (NcMultiplicityFuncWarren *mw, gdouble A)
{
  NcMultiplicityFuncWarrenPrivate * const self = mw->priv;

  g_assert (A >= 0);

  self->A = A;
}

/**
 * nc_multiplicity_func_warren_get_A:
 * @mw: a #NcMultiplicityFuncWarren.
 *
 * Returns: the value of #NcMultiplicityFuncWarren:A property.
 */
gdouble
nc_multiplicity_func_warren_get_A (const NcMultiplicityFuncWarren *mw)
{
  NcMultiplicityFuncWarrenPrivate * const self = mw->priv;

  return self->A;
}

/**
 * nc_multiplicity_func_warren_set_a:
 * @mw: a #NcMultiplicityFuncWarren.
 * @a: value of #NcMultiplicityFuncWarren:a.
 *
 * Sets the value @a to the #NcMultiplicityFuncWarren:a property.
 *
 */
void
nc_multiplicity_func_warren_set_a (NcMultiplicityFuncWarren *mw, gdouble a)
{
  NcMultiplicityFuncWarrenPrivate * const self = mw->priv;

  g_assert (a >= 0);

  self->a = a;
}

/**
 * nc_multiplicity_func_warren_get_a:
 * @mw: a #NcMultiplicityFuncWarren.
 *
 * Returns: the value of #NcMultiplicityFuncWarren:a property.
 */
gdouble
nc_multiplicity_func_warren_get_a (const NcMultiplicityFuncWarren *mw)
{
  NcMultiplicityFuncWarrenPrivate * const self = mw->priv;

  return self->a;
}

/**
 * nc_multiplicity_func_warren_set_b:
 * @mw: a #NcMultiplicityFuncWarren.
 * @b: value of #NcMultiplicityFuncWarren:b.
 *
 * Sets the value @b to the #NcMultiplicityFuncWarren:b property.
 *
 */
void
nc_multiplicity_func_warren_set_b (NcMultiplicityFuncWarren *mw, gdouble b)
{
  NcMultiplicityFuncWarrenPrivate * const self = mw->priv;

  g_assert (b >= 0);

  self->b = b;
}

/**
 * nc_multiplicity_func_warren_get_b:
 * @mw: a #NcMultiplicityFuncWarren.
 *
 * Returns: the value of #NcMultiplicityFuncWarren:b property.
 */
gdouble
nc_multiplicity_func_warren_get_b (const NcMultiplicityFuncWarren *mw)
{
  NcMultiplicityFuncWarrenPrivate * const self = mw->priv;

  return self->b;
}

/**
 * nc_multiplicity_func_warren_set_c:
 * @mw: a #NcMultiplicityFuncWarren.
 * @c: value of #NcMultiplicityFuncWarren:c.
 *
 * Sets the value @c to the #NcMultiplicityFuncWarren:c property.
 *
 */
void
nc_multiplicity_func_warren_set_c (NcMultiplicityFuncWarren *mw, gdouble c)
{
  NcMultiplicityFuncWarrenPrivate * const self = mw->priv;

  g_assert (c >= 0);

  self->c = c;
}

/**
 * nc_multiplicity_func_warren_get_c:
 * @mw: a #NcMultiplicityFuncWarren.
 *
 * Returns: the value of #NcMultiplicityFuncWarren:c property.
 */
gdouble
nc_multiplicity_func_warren_get_c (const NcMultiplicityFuncWarren *mw)
{
  NcMultiplicityFuncWarrenPrivate * const self = mw->priv;

  return self->c;
}