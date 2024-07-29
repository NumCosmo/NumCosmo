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
 * Computes the multiplicity function of dark matter halos using the Jenkins et al. (2001) model.
 * Jenkins et al. (2001) [arXiv:astro-ph/0005260] is a parametrization of the Press-Schechter multiplicity function.
 *
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
  gdouble Delta;
};

enum
{
  PROP_0,
  PROP_SIZE
};

G_DEFINE_TYPE_WITH_PRIVATE (NcMultiplicityFuncJenkins, nc_multiplicity_func_jenkins, NC_TYPE_MULTIPLICITY_FUNC)

static void
nc_multiplicity_func_jenkins_init (NcMultiplicityFuncJenkins *mj)
{
  NcMultiplicityFuncJenkinsPrivate * const self = mj->priv = nc_multiplicity_func_jenkins_get_instance_private (mj);

  self->mdef = NC_MULTIPLICITY_FUNC_MASS_DEF_LEN;
}

static void
_nc_multiplicity_func_jenkins_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  /* NcMultiplicityFuncJenkins *mj = NC_MULTIPLICITY_FUNC_JENKINS (object); */
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_JENKINS (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_jenkins_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /* NcMultiplicityFuncJenkins *mj = NC_MULTIPLICITY_FUNC_JENKINS (object); */
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_JENKINS (object));

  switch (prop_id)
  {
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
static void _nc_multiplicity_func_jenkins_set_Delta (NcMultiplicityFunc *mulf, gdouble Delta);
static double _nc_multiplicity_func_jenkins_get_Delta (NcMultiplicityFunc *mulf);
static gdouble _nc_multiplicity_func_jenkins_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z);

/* _NC_MULTIPLICITY_FUNCTION_JENKINS_DATASET_FOF_0005260 = {0.315, 0.0, 0.61, 0.0, 3.8, 0.0}; */

static void
nc_multiplicity_func_jenkins_class_init (NcMultiplicityFuncJenkinsClass *klass)
{
  GObjectClass *object_class            = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass *parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);

  object_class->set_property = _nc_multiplicity_func_jenkins_set_property;
  object_class->get_property = _nc_multiplicity_func_jenkins_get_property;
  object_class->finalize     = _nc_multiplicity_func_jenkins_finalize;

  parent_class->set_mdef  = &_nc_multiplicity_func_jenkins_set_mdef;
  parent_class->get_mdef  = &_nc_multiplicity_func_jenkins_get_mdef;
  parent_class->set_Delta = &_nc_multiplicity_func_jenkins_set_Delta;
  parent_class->get_Delta = &_nc_multiplicity_func_jenkins_get_Delta;
  parent_class->eval      = &_nc_multiplicity_func_jenkins_eval;
}

static void
_nc_multiplicity_func_jenkins_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef)
{
  NcMultiplicityFuncJenkins *mj                 = NC_MULTIPLICITY_FUNC_JENKINS (mulf);
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
  NcMultiplicityFuncJenkins *mj                 = NC_MULTIPLICITY_FUNC_JENKINS (mulf);
  NcMultiplicityFuncJenkinsPrivate * const self = mj->priv;

  return self->mdef;
}

static void
_nc_multiplicity_func_jenkins_set_Delta (NcMultiplicityFunc *mulf, gdouble Delta)
{
  NcMultiplicityFuncJenkins *mj                 = NC_MULTIPLICITY_FUNC_JENKINS (mulf);
  NcMultiplicityFuncJenkinsPrivate * const self = mj->priv;

  self->Delta = Delta;
}

static gdouble
_nc_multiplicity_func_jenkins_get_Delta (NcMultiplicityFunc *mulf)
{
  NcMultiplicityFuncJenkins *mj                 = NC_MULTIPLICITY_FUNC_JENKINS (mulf);
  NcMultiplicityFuncJenkinsPrivate * const self = mj->priv;

  return self->Delta;
}

static gdouble
_nc_multiplicity_func_jenkins_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)
{
  /* NcMultiplicityFuncJenkins *mj = NC_MULTIPLICITY_FUNC_JENKINS (mulf);
   *  NcMultiplicityFuncJenkinsPrivate * const self = mj->priv; */

  gdouble f_Jenkins = 0.315 * exp (-pow (fabs (-log (sigma) + 0.61), 3.8));

  NCM_UNUSED (mulf);
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

