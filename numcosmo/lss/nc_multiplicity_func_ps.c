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

G_DEFINE_TYPE (NcMultiplicityFuncPS, nc_multiplicity_func_ps, NC_TYPE_MULTIPLICITY_FUNC);

enum
{
  PROP_0,
  PROP_DELTA_C
};

/**
 * nc_multiplicity_func_ps_new:
 * @delta_c: FIXME
 *
 * FIXME
 *
 * Returns: A new #NcMultiplicityFunc.
 */
NcMultiplicityFunc *
nc_multiplicity_func_ps_new (gdouble delta_c)
{
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_PS,
                       "critical-delta", delta_c,
                       NULL);
}

static gdouble
_nc_multiplicity_func_ps_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)       /* f(\sigma) - Press \& Schechter (PS) */
{
  NcMultiplicityFuncPS *mulf_ps = NC_MULTIPLICITY_FUNC_PS (mulf);
  const gdouble c1 = sqrt (2.0 / M_PI);               /* c1 = \sqrt{\frac{2}{\pi}} */
  gdouble x = mulf_ps->delta_c / sigma;        /* \delta_c \sigma^{-1} */
  gdouble x2 = x * x;
  gdouble f_PS = c1 * x * exp(-x2 / 2.0);

  NCM_UNUSED (cosmo);
  NCM_UNUSED (z);
  
  return f_PS;
}

/**
 * nc_multiplicity_func_ps_set_delta_c:
 * @mulf_ps: a #NcMultiplicityFuncPS.
 * @delta_c: value of #NcMultiplicityFuncPS:critical-delta.
 *
 * Sets the value @delta_c to the #NcMultiplicityFuncPS:critical-delta property.
 *
 */
void
nc_multiplicity_func_ps_set_delta_c (NcMultiplicityFuncPS *mulf_ps, gdouble delta_c)
{
  g_assert (delta_c >= 0);
  mulf_ps->delta_c = delta_c;
}

/**
 * nc_multiplicity_func_ps_get_delta_c:
 * @mulf_ps: a #NcMultiplicityFuncPS.
 *
 * Returns: the value of #NcMultiplicityFuncPS:critical_delta property.
 */
gdouble
nc_multiplicity_func_ps_get_delta_c (const NcMultiplicityFuncPS *mulf_ps)
{
  return mulf_ps->delta_c;
}

static void
nc_multiplicity_func_ps_init (NcMultiplicityFuncPS *mulf_ps)
{
  /* TODO: Add initialization code here */
  mulf_ps->delta_c = 1.686;
}

static void
_nc_multiplicity_func_ps_finalize (GObject *object)
{
  /* TODO: Add deinitalization code here */

  G_OBJECT_CLASS (nc_multiplicity_func_ps_parent_class)->finalize (object);
}

static void
_nc_multiplicity_func_ps_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcMultiplicityFuncPS *mulf_ps = NC_MULTIPLICITY_FUNC_PS (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_PS (object));

  switch (prop_id)
  {
    case PROP_DELTA_C:
      mulf_ps->delta_c = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_ps_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcMultiplicityFuncPS *mulf_ps = NC_MULTIPLICITY_FUNC_PS (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_PS (object));

  switch (prop_id)
  {
    case PROP_DELTA_C:
      g_value_set_double (value, mulf_ps->delta_c);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_multiplicity_func_ps_class_init (NcMultiplicityFuncPSClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass* parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);

  parent_class->eval = &_nc_multiplicity_func_ps_eval;

  object_class->finalize = _nc_multiplicity_func_ps_finalize;
  object_class->set_property = _nc_multiplicity_func_ps_set_property;
  object_class->get_property = _nc_multiplicity_func_ps_get_property;

  /**
   * NcMultiplicityFuncPS:critical_delta:
   *
   * FIXME
   */
  g_object_class_install_property (object_class,
                                   PROP_DELTA_C,
                                   g_param_spec_double ("critical-delta",
                                                        NULL,
                                                        "Critical delta",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.686,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

