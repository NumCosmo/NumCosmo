/***************************************************************************
 *            nc_multiplicity_func_watson.c
 *
 *  Sat Sep 11 17:45:13 2021
 *  Copyright  2021  Cinthia Nunes de Lima / Mariana Penna Lima
 *  <cinthia.n.lima@uel.br> / <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Cinthia Nunes de Lima <cinthia.n.lima@uel.br>, Mariana Penna Lima 2012 <pennalima@gmail.com>
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
 * SECTION:nc_multiplicity_func_watson
 * @title: NcMultiplicityFuncWatson
 * @short_description: Dark matter halo -- Watson multiplicity function.
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_multiplicity_func_watson.h"

struct _NcMultiplicityFuncWatsonPrivate
{
  NcMultiplicityFuncMassDef mdef;
};

enum
{
  PROP_0,
  PROP_SIZE
};

G_DEFINE_TYPE_WITH_PRIVATE (NcMultiplicityFuncWatson, nc_multiplicity_func_watson, NC_TYPE_MULTIPLICITY_FUNC);

static void
nc_multiplicity_func_watson_init (NcMultiplicityFuncWatson *mwat)
{
  NcMultiplicityFuncWatsonPrivate * const self = mwat->priv = nc_multiplicity_func_watson_get_instance_private (mwat);

  self->mdef    = NC_MULTIPLICITY_FUNC_MASS_DEF_LEN;
}

static void
_nc_multiplicity_func_watson_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  /* NcMultiplicityFuncWatson *mwat = NC_MULTIPLICITY_FUNC_WATSON (object); */
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_WATSON (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_watson_get_property (GObject * object, guint prop_id, GValue * value, GParamSpec * pspec)
{
  /* NcMultiplicityFuncWatson *mwat = NC_MULTIPLICITY_FUNC_WATSON (object); */
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_WATSON (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_watson_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_multiplicity_func_watson_parent_class)->finalize (object);
}

static void _nc_multiplicity_func_watson_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef); 
static NcMultiplicityFuncMassDef _nc_multiplicity_func_watson_get_mdef (NcMultiplicityFunc *mulf);
static gdouble _nc_multiplicity_func_watson_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z);

// _NC_MULTIPLICITY_FUNCTION_WATSON_DATASET_FOF_0005260 = {0.315, 0.0, 0.61, 0.0, 3.8, 0.0};

static void
nc_multiplicity_func_watson_class_init (NcMultiplicityFuncWatsonClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass* parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);

  object_class->set_property = _nc_multiplicity_func_watson_set_property;
  object_class->get_property = _nc_multiplicity_func_watson_get_property;
  object_class->finalize     = _nc_multiplicity_func_watson_finalize;

  parent_class->set_mdef = &_nc_multiplicity_func_watson_set_mdef;
  parent_class->get_mdef = &_nc_multiplicity_func_watson_get_mdef;
  parent_class->eval     = &_nc_multiplicity_func_watson_eval;
}

static void 
_nc_multiplicity_func_watson_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef)
{
  NcMultiplicityFuncWatson *mwat = NC_MULTIPLICITY_FUNC_WATSON (mulf);
  NcMultiplicityFuncWatsonPrivate * const self = mwat->priv;

  switch (mdef)
  {
    case NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN:
      g_error ("NcMultiplicityFuncWatson does not support fof mass def");
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL:
      g_error ("NcMultiplicityFuncWatson does not support critical mass def");
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_VIRIAL:
      g_error ("NcMultiplicityFuncWatson does not support virial mass def");
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
_nc_multiplicity_func_watson_get_mdef (NcMultiplicityFunc *mulf)
{
  NcMultiplicityFuncWatson *mwat = NC_MULTIPLICITY_FUNC_WATSON (mulf);
  NcMultiplicityFuncWatsonPrivate * const self = mwat->priv;

  return self->mdef;
}

static gdouble
_nc_multiplicity_func_watson_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)
{
  /* NcMultiplicityFuncWatson *mwat = NC_MULTIPLICITY_FUNC_WATSON (mulf);
  NcMultiplicityFuncWatsonPrivate * const self = mwat->priv; */
  
  gdouble f_Watson = 0.315 * exp(-pow(fabs(-log(sigma) + 0.61), 3.8));

  NCM_UNUSED (mulf);
  NCM_UNUSED (cosmo);
  NCM_UNUSED (z);

  return f_Watson;
}

/**
 * nc_multiplicity_func_watson_new:
 *   
 * FIXME
 *
 * Returns: A new #NcMultiplicityFuncWatson.
 */
NcMultiplicityFuncWatson *
nc_multiplicity_func_watson_new (void)
{
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_WATSON,
                       "mass-def", NC_MULTIPLICITY_FUNC_MASS_DEF_FOF,
                       NULL);
}

/**
 * nc_multiplicity_func_watson_ref:
 * @mwat: a #NcMultiplicityFuncWatson
 *
 * Increases the reference count of @mwat by one.
 *
 * Returns: (transfer full): @mwat
 */
NcMultiplicityFuncWatson *
nc_multiplicity_func_watson_ref (NcMultiplicityFuncWatson *mwat)
{
  return g_object_ref (mwat);
}

/**
 * nc_multiplicity_func_watson_free:
 * @mwat: a #NcMultiplicityFuncWatson
 *
 * Atomically decrements the reference count of @mwat by one. If the reference count drops to 0,
 * all memory allocated by @mwat is released.
 *
 */
void
nc_multiplicity_func_watson_free (NcMultiplicityFuncWatson *mwat)
{
  g_object_unref (mwat);
}

/**
 * nc_multiplicity_func_watson_clear:
 * @mwat: a #NcMultiplicityFuncWatson
 *
 * Atomically decrements the reference count of @mwat by one. If the reference count drops to 0,
 * all memory allocated by @mwat is released. Set the pointer to NULL;
 *
 */
void
nc_multiplicity_func_watson_clear (NcMultiplicityFuncWatson **mwat)
{
  g_clear_object (mwat);
}

