/***************************************************************************
 *            nc_multiplicity_func_jenkins.c
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
 * SECTION:nc_multiplicity_func_jenkins
 * @title: NcMultiplicityFuncJenkins
 * @short_description: Dark matter halo -- Jenkins multiplicity function.
 *
 * FIXME
 * Reference: astro-ph/0005260
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_multiplicity_func_jenkins.h"

struct _NcMultiplicityFuncJenkinsPrivate
{
  NcMultiplicityFuncMassDef mdef;
  gdouble A; 
  gdouble B;
  gdouble epsilon;
};

enum
{
  PROP_0,
  PROP_A,
  PROP_B,
  PROP_EPSILON,
  PROP_SIZE
};

G_DEFINE_TYPE_WITH_PRIVATE (NcMultiplicityFuncJenkins, nc_multiplicity_func_jenkins, NC_TYPE_MULTIPLICITY_FUNC);

static void
nc_multiplicity_func_jenkins_init (NcMultiplicityFuncJenkins *mj)
{
  NcMultiplicityFuncJenkinsPrivate * const self = mj->priv = nc_multiplicity_func_jenkins_get_instance_private (mj);

  self->mdef    = NC_MULTIPLICITY_FUNC_MASS_DEF_LEN;
  self->A       = 0.0;
  self->B       = 0.0;
  self->epsilon = 0.0;
}

static void
_nc_multiplicity_func_jenkins_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcMultiplicityFuncJenkins *mj = NC_MULTIPLICITY_FUNC_JENKINS (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_JENKINS (object));

  switch (prop_id)
  {
    case PROP_A:
      nc_multiplicity_func_jenkins_set_A (mj, g_value_get_double (value));
      break;
    case PROP_B:
      nc_multiplicity_func_jenkins_set_B (mj, g_value_get_double (value));
      break;
    case PROP_EPSILON:
      nc_multiplicity_func_jenkins_set_epsilon (mj, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_jenkins_get_property (GObject * object, guint prop_id, GValue * value, GParamSpec * pspec)
{
  NcMultiplicityFuncJenkins *mj = NC_MULTIPLICITY_FUNC_JENKINS (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_JENKINS (object));

  switch (prop_id)
  {
    case PROP_A:
      g_value_set_double (value, nc_multiplicity_func_jenkins_get_A (mj));
      break;
    case PROP_B:
      g_value_set_double (value, nc_multiplicity_func_jenkins_get_B (mj));
      break;
    case PROP_EPSILON:
      g_value_set_double (value, nc_multiplicity_func_jenkins_get_epsilon (mj));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_jenkins_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_multiplicity_func_jenkins_parent_class)->finalize (object);
}

static void _nc_multiplicity_func_jenkins_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef); 
static NcMultiplicityFuncMassDef _nc_multiplicity_func_jenkins_get_mdef (NcMultiplicityFunc *mulf);
static gdouble _nc_multiplicity_func_jenkins_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z);

// _NC_MULTIPLICITY_FUNCTION_JENKINS_DATASET_FOF_0005260 = {0.315, 0.0, 0.61, 0.0, 3.8, 0.0};

static void
nc_multiplicity_func_jenkins_class_init (NcMultiplicityFuncJenkinsClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass* parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);

  object_class->set_property = _nc_multiplicity_func_jenkins_set_property;
  object_class->get_property = _nc_multiplicity_func_jenkins_get_property;
  object_class->finalize     = _nc_multiplicity_func_jenkins_finalize;

  /**
   * NcMultiplicityFuncJenkins:A:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_A,
                                   g_param_spec_double ("A",
                                                        NULL,
                                                        "A",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.315,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcMultiplicityFuncJenkins:B:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_B,
                                   g_param_spec_double ("B",
                                                        NULL,
                                                        "B",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.61,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

 /**
   * NcMultiplicityFuncJenkins:epsilon:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_EPSILON,
                                   g_param_spec_double ("epsilon",
                                                        NULL,
                                                        "epsilon",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 3.8,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  parent_class->set_mdef = &_nc_multiplicity_func_jenkins_set_mdef;
  parent_class->get_mdef = &_nc_multiplicity_func_jenkins_get_mdef;
  parent_class->eval     = &_nc_multiplicity_func_jenkins_eval;
}

static void 
_nc_multiplicity_func_jenkins_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef)
{
  NcMultiplicityFuncJenkins *mj = NC_MULTIPLICITY_FUNC_JENKINS (mulf);
  NcMultiplicityFuncJenkinsPrivate * const self = mj->priv;

  switch (mdef)
  {
    case NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN:
      g_error ("NcMultiplicityFuncJenkins does not support fof mass def");
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL:
      g_error ("NcMultiplicityFuncJenkins does not support critical mass def");
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_VIRIAL:
      g_error ("NcMultiplicityFuncJenkins does not support virial mass def");
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
_nc_multiplicity_func_jenkins_get_mdef (NcMultiplicityFunc *mulf)
{
  NcMultiplicityFuncJenkins *mj = NC_MULTIPLICITY_FUNC_JENKINS (mulf);
  NcMultiplicityFuncJenkinsPrivate * const self = mj->priv;

  return self->mdef;
}

static gdouble
_nc_multiplicity_func_jenkins_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)
{
  NcMultiplicityFuncJenkins *mj = NC_MULTIPLICITY_FUNC_JENKINS (mulf);
  NcMultiplicityFuncJenkinsPrivate * const self = mj->priv;
  
  gdouble f_Jenkins = self->A * exp(-pow(fabs(-log(sigma) + self->B), self->epsilon));

  NCM_UNUSED (cosmo);
  NCM_UNUSED (z);

  return f_Jenkins;
}

/**
 * nc_multiplicity_func_jenkins_new:
 *   
 * FIXME
 *
 * Returns: A new #NcMultiplicityFuncJenkins.
 */
NcMultiplicityFuncJenkins *
nc_multiplicity_func_jenkins_new (void)
{
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_JENKINS,
                       "mass-def", NC_MULTIPLICITY_FUNC_MASS_DEF_FOF,
                       NULL);
}

/**
 * nc_multiplicity_func_jenkins_ref:
 * @mj: a #NcMultiplicityFuncJenkins
 *
 * Increases the reference count of @mj by one.
 *
 * Returns: (transfer full): @mj
 */
NcMultiplicityFuncJenkins *
nc_multiplicity_func_jenkins_ref (NcMultiplicityFuncJenkins *mj)
{
  return g_object_ref (mj);
}

/**
 * nc_multiplicity_func_jenkins_free:
 * @mj: a #NcMultiplicityFuncJenkins
 *
 * Atomically decrements the reference count of @mj by one. If the reference count drops to 0,
 * all memory allocated by @mj is released.
 *
 */
void
nc_multiplicity_func_jenkins_free (NcMultiplicityFuncJenkins *mj)
{
  g_object_unref (mj);
}

/**
 * nc_multiplicity_func_jenkins_clear:
 * @mj: a #NcMultiplicityFuncJenkins
 *
 * Atomically decrements the reference count of @mj by one. If the reference count drops to 0,
 * all memory allocated by @mj is released. Set the pointer to NULL;
 *
 */
void
nc_multiplicity_func_jenkins_clear (NcMultiplicityFuncJenkins **mj)
{
  g_clear_object (mj);
}

/**
 * nc_multiplicity_func_jenkins_set_A:
 * @mj: a #NcMultiplicityFuncJenkins.
 * @A: value of #NcMultiplicityFuncJenkins:A.
 *
 * Sets the value @A to the #NcMultiplicityFuncJenkins:A property.
 *
 */
void
nc_multiplicity_func_jenkins_set_A (NcMultiplicityFuncJenkins *mj, gdouble A)
{
  NcMultiplicityFuncJenkinsPrivate * const self = mj->priv;

  g_assert (A >= 0);

  self->A = A;
}

/**
 * nc_multiplicity_func_jenkins_get_A:
 * @mj: a #NcMultiplicityFuncJenkins.
 *
 * Returns: the value of #NcMultiplicityFuncJenkins:A property.
 */
gdouble
nc_multiplicity_func_jenkins_get_A (const NcMultiplicityFuncJenkins *mj)
{
  NcMultiplicityFuncJenkinsPrivate * const self = mj->priv;

  return self->A;
}

/**
 * nc_multiplicity_func_jenkins_set_B:
 * @mj: a #NcMultiplicityFuncJenkins.
 * @B: value of #NcMultiplicityFuncJenkins:B.
 *
 * Sets the value @B to the #NcMultiplicityFuncJenkins:B property.
 *
 */
void
nc_multiplicity_func_jenkins_set_B (NcMultiplicityFuncJenkins *mj, gdouble B)
{
  NcMultiplicityFuncJenkinsPrivate * const self = mj->priv;

  g_assert (B >= 0);

  self->B = B;
}

/**
 * nc_multiplicity_func_jenkins_get_B:
 * @mj: a #NcMultiplicityFuncJenkins.
 *
 * Returns: the value of #NcMultiplicityFuncJenkins:B property.
 */
gdouble
nc_multiplicity_func_jenkins_get_B (const NcMultiplicityFuncJenkins *mj)
{
  NcMultiplicityFuncJenkinsPrivate * const self = mj->priv;

  return self->B;
}

/**
 * nc_multiplicity_func_jenkins_set_epsilon:
 * @mj: a #NcMultiplicityFuncJenkins.
 * @epsilon: value of #NcMultiplicityFuncJenkins:epsilon.
 *
 * Sets the value @epsilon to the #NcMultiplicityFuncJenkins:epsilon property.
 *
 */
void
nc_multiplicity_func_jenkins_set_epsilon (NcMultiplicityFuncJenkins *mj, gdouble epsilon)
{
  NcMultiplicityFuncJenkinsPrivate * const self = mj->priv;

  g_assert (epsilon >= 0);

  self->epsilon = epsilon;
}

/**
 * nc_multiplicity_func_jenkins_get_epsilon:
 * @mj: a #NcMultiplicityFuncJenkins.
 *
 * Returns: the value of #NcMultiplicityFuncJenkins:epsilon property.
 */
gdouble
nc_multiplicity_func_jenkins_get_epsilon (const NcMultiplicityFuncJenkins *mj)
{
  NcMultiplicityFuncJenkinsPrivate * const self = mj->priv;

  return self->epsilon;
}

