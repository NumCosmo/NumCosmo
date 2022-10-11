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
  gdouble Delta;
};

enum
{
  PROP_0,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcMultiplicityFuncWarren, nc_multiplicity_func_warren, NC_TYPE_MULTIPLICITY_FUNC);

static void
nc_multiplicity_func_warren_init (NcMultiplicityFuncWarren *mw)
{
  NcMultiplicityFuncWarrenPrivate * const self = mw->priv = nc_multiplicity_func_warren_get_instance_private (mw);

  self->mdef    = NC_MULTIPLICITY_FUNC_MASS_DEF_LEN;
  self->Delta = 0.0;
}

static void
_nc_multiplicity_func_warren_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  /* NcMultiplicityFuncWarren *mw = NC_MULTIPLICITY_FUNC_WARREN (object); */
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_WARREN (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_warren_get_property (GObject * object, guint prop_id, GValue * value, GParamSpec * pspec)
{
  /*NcMultiplicityFuncWarren *mw = NC_MULTIPLICITY_FUNC_WARREN (object); */
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_WARREN (object));

  switch (prop_id)
  {
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

static void _nc_multiplicity_func_warren_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef);
static void _nc_multiplicity_func_warren_set_Delta (NcMultiplicityFunc *mulf, gdouble Delta); 
static NcMultiplicityFuncMassDef _nc_multiplicity_func_warren_get_mdef (NcMultiplicityFunc *mulf);
static gdouble _nc_multiplicity_func_warren_get_Delta (NcMultiplicityFunc *mulf); 
static gdouble _nc_multiplicity_func_warren_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z);

static void
nc_multiplicity_func_warren_class_init (NcMultiplicityFuncWarrenClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass* parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);

  object_class->set_property = &_nc_multiplicity_func_warren_set_property;
  object_class->get_property = &_nc_multiplicity_func_warren_get_property;
  object_class->finalize     = &_nc_multiplicity_func_warren_finalize;
  
  parent_class->set_mdef = &_nc_multiplicity_func_warren_set_mdef;
  parent_class->set_Delta = &_nc_multiplicity_func_warren_set_Delta;
  parent_class->get_mdef = &_nc_multiplicity_func_warren_get_mdef;
  parent_class->get_Delta = &_nc_multiplicity_func_warren_get_Delta;
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
  /* NcMultiplicityFuncWarren *mw = NC_MULTIPLICITY_FUNC_WARREN (mulf);
  NcMultiplicityFuncWarrenPrivate * const self = mw->priv; */
  
  gdouble f_Warren = 0.7234 * (pow(sigma, - 1.625) + 0.2538) * exp(-1.1982 / (sigma * sigma) );

  NCM_UNUSED (mulf);
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
 * nc_multiplicity_func_warren_set_Delta:
 * @mt: a #NcMultiplicityFuncWarren.
 * @Delta: value of #NcMultiplicityFuncWarren:Delta.
 *
 * Sets the value @Delta to the #NcMultiplicityFuncWarren:Delta property.
 *
 */
static void
_nc_multiplicity_func_warren_set_Delta (NcMultiplicityFunc *mulf, gdouble Delta)
{
  NcMultiplicityFuncWarren *mw = NC_MULTIPLICITY_FUNC_WARREN (mulf);
  NcMultiplicityFuncWarrenPrivate * const self = mw->priv;

  self->Delta = Delta;
}

/**
 * nc_multiplicity_func_warren_get_Delta:
 * @mt: a #NcMultiplicityFuncWarren.
 *
 * Returns: the value of #NcMultiplicityFuncWarren:Delta property.
 */
gdouble
_nc_multiplicity_func_warren_get_Delta (NcMultiplicityFunc *mulf)
{
  NcMultiplicityFuncWarren *mw = NC_MULTIPLICITY_FUNC_WARREN (mulf);
  NcMultiplicityFuncWarrenPrivate * const self = mw->priv;
  
  return self->Delta;
}