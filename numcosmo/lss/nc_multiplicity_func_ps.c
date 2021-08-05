/***************************************************************************
 *            nc_multiplicity_func_ps.c
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
 * SECTION:nc_multiplicity_func_ps
 * @title: NcMultiplicityFuncPS
 * @short_description: Dark matter halo -- Press-Schechter multiplicity function.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_multiplicity_func_ps.h"

struct _NcMultiplicityFuncPSPrivate
{
  NcMultiplicityFuncMassDef mdef;
  gdouble delta_c;
};

enum
{
  PROP_0,
  PROP_DELTA_C,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcMultiplicityFuncPS, nc_multiplicity_func_ps, NC_TYPE_MULTIPLICITY_FUNC);

static void
nc_multiplicity_func_ps_init (NcMultiplicityFuncPS *mps)
{
  NcMultiplicityFuncPSPrivate * const self = mps->priv = nc_multiplicity_func_ps_get_instance_private (mps);

  self->mdef    = NC_MULTIPLICITY_FUNC_MASS_DEF_LEN;
  self->delta_c = 0.0;
}

static void
_nc_multiplicity_func_ps_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcMultiplicityFuncPS *mps = NC_MULTIPLICITY_FUNC_PS (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_PS (object));

  switch (prop_id)
  {
    case PROP_DELTA_C:
      nc_multiplicity_func_ps_set_delta_c (mps, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_ps_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcMultiplicityFuncPS *mps = NC_MULTIPLICITY_FUNC_PS (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_PS (object));

  switch (prop_id)
  {
    case PROP_DELTA_C:
      g_value_set_double (value, nc_multiplicity_func_ps_get_delta_c (mps));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_ps_finalize (GObject *object)
{

  G_OBJECT_CLASS (nc_multiplicity_func_ps_parent_class)->finalize (object);
}

static void _nc_multiplicity_func_ps_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef); 
static NcMultiplicityFuncMassDef _nc_multiplicity_func_ps_get_mdef (NcMultiplicityFunc *mulf);
static gdouble _nc_multiplicity_func_ps_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z);

static void
nc_multiplicity_func_ps_class_init (NcMultiplicityFuncPSClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass* parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);

  object_class->set_property = &_nc_multiplicity_func_ps_set_property;
  object_class->get_property = &_nc_multiplicity_func_ps_get_property;
  object_class->finalize     = &_nc_multiplicity_func_ps_finalize;

  /**
   * NcMultiplicityFuncPS:critical_delta:
   *
   * Critical matter density contrast.
   */
  g_object_class_install_property (object_class,
                                   PROP_DELTA_C,
                                   g_param_spec_double ("critical-delta",
                                                        NULL,
                                                        "Critical delta",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.6864701998411454502,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  parent_class->set_mdef = &_nc_multiplicity_func_ps_set_mdef;
  parent_class->get_mdef = &_nc_multiplicity_func_ps_get_mdef;
  parent_class->eval     = &_nc_multiplicity_func_ps_eval;       
}

static void 
_nc_multiplicity_func_ps_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef)
{
  NcMultiplicityFuncPS *mps = NC_MULTIPLICITY_FUNC_PS (mulf);
  NcMultiplicityFuncPSPrivate * const self = mps->priv;

  switch (mdef)
  {
    case NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN:
      /* nothing to do */
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL:
      g_error ("NcMultiplicityFuncPS does not support critical mass def");
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_VIRIAL:
      g_error ("NcMultiplicityFuncPS does not support virial mass def");
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_FOF:
      g_error ("NcMultiplicityFuncPS does not support fof mass def");
      break;
    default:
      g_assert_not_reached ();
      break;
  }
  
  self->mdef = mdef;
}

static NcMultiplicityFuncMassDef 
_nc_multiplicity_func_ps_get_mdef (NcMultiplicityFunc *mulf)
{
  NcMultiplicityFuncPS *mps = NC_MULTIPLICITY_FUNC_PS (mulf);
  NcMultiplicityFuncPSPrivate * const self = mps->priv;

  return self->mdef;
}


static gdouble
_nc_multiplicity_func_ps_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)       /* f(\sigma) - Press \& Schechter (PS) */
{
  NcMultiplicityFuncPS *mps = NC_MULTIPLICITY_FUNC_PS (mulf);
  NcMultiplicityFuncPSPrivate * const self = mps->priv;

  const gdouble c1   = sqrt (2.0 / M_PI);     /* c1 = \sqrt{\frac{2}{\pi}} */
  const gdouble x    = self->delta_c / sigma; /* \delta_c \sigma^{-1} */
  const gdouble x2   = x * x;
  const gdouble f_PS = c1 * x * exp(-0.5 * x2);

  return f_PS;
}

/**
 * nc_multiplicity_func_ps_new:
 *
 * FIXME
 *
 * Returns: A new #NcMultiplicityFunc.
 */
NcMultiplicityFuncPS *
nc_multiplicity_func_ps_new (void)
{
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_PS,
                       "mass-def", NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN,
                       NULL);
}

/**
 * nc_multiplicity_func_ps_ref:
 * @mps: a #NcMultiplicityFuncPS
 *
 * Increases the reference count of @mps by one.
 *
 * Returns: (transfer full): @mps
 */
NcMultiplicityFuncPS *
nc_multiplicity_func_ps_ref (NcMultiplicityFuncPS *mps)
{
  return g_object_ref (mps);
}

/**
 * nc_multiplicity_func_ps_free:
 * @mps: a #NcMultiplicityFuncPS
 *
 * Atomically decrements the reference count of @mps by one. If the reference count drops to 0,
 * all memory allocated by @mps is released.
 *
 */
void
nc_multiplicity_func_ps_free (NcMultiplicityFuncPS *mps)
{
  g_object_unref (mps);
}

/**
 * nc_multiplicity_func_ps_clear:
 * @mps: a #NcMultiplicityFuncPS
 *
 * Atomically decrements the reference count of @mps by one. If the reference count drops to 0,
 * all memory allocated by @mps is released. Set the pointer to NULL;
 *
 */
void
nc_multiplicity_func_ps_clear (NcMultiplicityFuncPS **mps)
{
  g_clear_object (mps);
}

/**
 * nc_multiplicity_func_ps_set_delta_c:
 * @mps: a #NcMultiplicityFuncPS.
 * @delta_c: value of #NcMultiplicityFuncPS:critical-delta.
 *
 * Sets the value @delta_c to the #NcMultiplicityFuncPS:critical-delta property.
 *
 */
void
nc_multiplicity_func_ps_set_delta_c (NcMultiplicityFuncPS *mps, gdouble delta_c)
{
  NcMultiplicityFuncPSPrivate * const self = mps->priv;
  
  g_assert (delta_c >= 0);
  
  self->delta_c = delta_c;
}

/**
 * nc_multiplicity_func_ps_get_delta_c:
 * @mps: a #NcMultiplicityFuncPS.
 *
 * Returns: the value of #NcMultiplicityFuncPS:critical-delta property.
 */
gdouble
nc_multiplicity_func_ps_get_delta_c (const NcMultiplicityFuncPS *mps)
{
  NcMultiplicityFuncPSPrivate * const self = mps->priv;

  return self->delta_c;
}

