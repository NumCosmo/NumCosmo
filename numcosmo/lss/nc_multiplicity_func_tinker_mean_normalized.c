/***************************************************************************
 *            nc_multiplicity_func_tinker_mean_normalized.c
 *
 *  Thu July 03 12:29:55 2014
 *  Copyright  2014  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_multiplicity_func_tinker_mean_normalized.c
 * Copyright (C) 2014 Mariana Penna Lima <pennalima@gmail.com>
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
 * SECTION:nc_multiplicity_func_tinker_mean_normalized
 * @title: NcMultiplicityFuncTinkerMeanNormalized
 * @short_description: Dark matter halo -- Tinker normalized multiplicity function mean matter density.
 *
 * Normalized multiplicity function that has to be used in the halo model approach.
 * 
 * FIXME
 * Reference: arXiv:1001.3162
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_multiplicity_func_tinker_mean_normalized.h"

struct _NcMultiplicityFuncTinkerMeanNormalizedPrivate
{
  NcMultiplicityFuncMassDef mdef;
  gdouble Delta;
  guint int_Delta;
  gdouble alpha;
  gdouble beta;
  gdouble phi;
  gdouble eta;
  gdouble gamma;
};

enum
{
  PROP_0,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcMultiplicityFuncTinkerMeanNormalized, nc_multiplicity_func_tinker_mean_normalized, NC_TYPE_MULTIPLICITY_FUNC);

static void
nc_multiplicity_func_tinker_mean_normalized_init (NcMultiplicityFuncTinkerMeanNormalized *mt10)
{
  NcMultiplicityFuncTinkerMeanNormalizedPrivate * const self = mt10->priv = nc_multiplicity_func_tinker_mean_normalized_get_instance_private (mt10);

  self->mdef      = NC_MULTIPLICITY_FUNC_MASS_DEF_LEN;  
  self->Delta     = 0.0;
  self->int_Delta = 0;
  self->alpha     = 0.0;
  self->beta      = 0.0;
  self->phi       = 0.0;
  self->eta       = 0.0;
  self->gamma     = 0.0;
}

static void
_nc_multiplicity_func_tinker_mean_normalized_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  /*NcMultiplicityFuncTinkerMeanNormalized *mt10 = NC_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED (object);*/
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_tinker_mean_normalized_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /*NcMultiplicityFuncTinkerMeanNormalized *mt10 = NC_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED (object);*/
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_tinker_mean_normalized_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_multiplicity_func_tinker_mean_normalized_parent_class)->finalize (object);
}

static void _nc_multiplicity_func_tinker_mean_normalized_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef);
static void _nc_multiplicity_func_tinker_mean_normalized_set_Delta (NcMultiplicityFunc *mulf, gdouble Delta);
static NcMultiplicityFuncMassDef _nc_multiplicity_func_tinker_mean_normalized_get_mdef (NcMultiplicityFunc *mulf);
static gdouble _nc_multiplicity_func_tinker_mean_normalized_get_Delta (NcMultiplicityFunc *mulf);
static gdouble _nc_multiplicity_func_tinker_mean_normalized_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z);

static void
nc_multiplicity_func_tinker_mean_normalized_class_init (NcMultiplicityFuncTinkerMeanNormalizedClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass* parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);

  object_class->set_property = &_nc_multiplicity_func_tinker_mean_normalized_set_property;
  object_class->get_property = &_nc_multiplicity_func_tinker_mean_normalized_get_property;
  object_class->finalize     = &_nc_multiplicity_func_tinker_mean_normalized_finalize;



  parent_class->set_mdef = &_nc_multiplicity_func_tinker_mean_normalized_set_mdef;
  parent_class->set_Delta = &_nc_multiplicity_func_tinker_mean_normalized_set_Delta;
  parent_class->get_mdef = &_nc_multiplicity_func_tinker_mean_normalized_get_mdef;
  parent_class->get_Delta = &_nc_multiplicity_func_tinker_mean_normalized_get_Delta;
  parent_class->eval     = &_nc_multiplicity_func_tinker_mean_normalized_eval;
}

static void 
_nc_multiplicity_func_tinker_mean_normalized_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef)
{
  NcMultiplicityFuncTinkerMeanNormalized *mt10 = NC_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED (mulf);
  NcMultiplicityFuncTinkerMeanNormalizedPrivate * const self = mt10->priv;

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
_nc_multiplicity_func_tinker_mean_normalized_get_mdef (NcMultiplicityFunc *mulf)
{
  NcMultiplicityFuncTinkerMeanNormalized *mt10 = NC_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED (mulf);
  NcMultiplicityFuncTinkerMeanNormalizedPrivate * const self = mt10->priv;

  return self->mdef;
}

static gdouble
_nc_multiplicity_func_tinker_mean_normalized_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)   /* $g(\sigma) = \nu x f(\nu)$ Tinker: Eq. 8 1001.3162 */
{
  NcMultiplicityFuncTinkerMeanNormalized *mt10 = NC_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED (mulf);
  NcMultiplicityFuncTinkerMeanNormalizedPrivate * const self = mt10->priv;

  const gdouble nu = 1.686 / sigma;
  gdouble f_Tinker_mean_normalized;
  
  NCM_UNUSED (cosmo);

  if (self->Delta == 200.0)
  {
    const gdouble alpha_z = self->alpha;
    const gdouble beta_z  = self->beta * pow(1.0 + z, 0.2);
    const gdouble phi_z   = self->phi * pow(1.0 + z, -0.08);
    const gdouble eta_z   = self->eta * pow(1.0 + z, 0.27);
    const gdouble gamma_z = self->gamma * pow(1.0 + z, -0.01);
    f_Tinker_mean_normalized = alpha_z * (1.0 + pow(beta_z * nu, -2.0 * phi_z)) 
                               * pow(nu, 2.0 * eta_z) * exp(-gamma_z * nu * nu / 2.0) * nu;
  }    
  else 
  {
    f_Tinker_mean_normalized = nu * self->alpha * (1.0 + pow(self->beta * nu, -2.0 * self->phi)) 
                               * pow(nu, 2.0 * self->eta) * exp(-self->gamma * nu * nu / 2.0);
  }
  
    return f_Tinker_mean_normalized;
}

/**
 * nc_multiplicity_func_tinker_mean_normalized_new:
 *   
 * FIXME
 *
 * Returns: A new #NcMultiplicityFuncTinkerMeanNormalized.
 */
NcMultiplicityFuncTinkerMeanNormalized *
nc_multiplicity_func_tinker_mean_normalized_new (void)
{
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED,
                       "mass-def", NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN,
                       NULL);
}

/**
 * nc_multiplicity_func_tinker_mean_normalized_ref:
 * @mt10: a #NcMultiplicityFuncTinkerMeanNormalized
 *
 * Increases the reference count of @mt10 by one.
 *
 * Returns: (transfer full): @mt10
 */
NcMultiplicityFuncTinkerMeanNormalized *
nc_multiplicity_func_tinker_mean_normalized_ref (NcMultiplicityFuncTinkerMeanNormalized *mt10)
{
  return g_object_ref (mt10);
}

/**
 * nc_multiplicity_func_tinker_mean_normalized_free:
 * @mt10: a #NcMultiplicityFuncTinkerMeanNormalized
 *
 * Atomically decrements the reference count of @mt10 by one. If the reference count drops to 0,
 * all memory allocated by @mt10 is released.
 *
 */
void
nc_multiplicity_func_tinker_mean_normalized_free (NcMultiplicityFuncTinkerMeanNormalized *mt10)
{
  g_object_unref (mt10);
}

/**
 * nc_multiplicity_func_tinker_mean_normalized_clear:
 * @mt10: a #NcMultiplicityFuncTinkerMeanNormalized
 *
 * Atomically decrements the reference count of @mt10 by one. If the reference count drops to 0,
 * all memory allocated by @mt10 is released. Set the pointer to NULL;
 *
 */
void
nc_multiplicity_func_tinker_mean_normalized_clear (NcMultiplicityFuncTinkerMeanNormalized **mt10)
{
  g_clear_object (mt10);
}

/**
 * nc_multiplicity_func_tinker_mean_normalized_set_Delta:
 * @mt10: a #NcMultiplicityFuncTinkerMeanNormalized.
 * @Delta: value of #NcMultiplicityFuncTinkerMeanNormalized:Delta.
 *
 * Sets the value @Delta to the #NcMultiplicityFuncTinkerMeanNormalized:Delta property.
 *
 */
static void
_nc_multiplicity_func_tinker_mean_normalized_set_Delta (NcMultiplicityFunc *mulf, gdouble Delta)
{
  NcMultiplicityFuncTinkerMeanNormalized *mt10 = NC_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED (mulf);
  NcMultiplicityFuncTinkerMeanNormalizedPrivate * const self = mt10->priv;
  
  const guint int_Delta = Delta;
  g_assert (Delta >= 0);
  g_assert (Delta == int_Delta);
  
  self->Delta = Delta;
  self->int_Delta = int_Delta;

  switch (self->int_Delta)
  {
    case 200:
      self->alpha =  0.368;
      self->beta  =  0.589;
      self->phi   =  -0.729;
      self->eta   = -0.243;
      self->gamma = 0.864;
      break;
    case 300:
      self->alpha = 0.363;
      self->beta  = 0.585;
      self->phi   = -0.789;
      self->eta   = -0.261;
      self->gamma = 0.922;
      break;
    case 400:
      self->alpha = 0.385;
      self->beta  = 0.544;
      self->phi   = -0.910;
      self->eta   = -0.261;
      self->gamma = 0.987;
      break;  
    case 600:
     self->alpha = 0.389;
     self->beta  = 0.543;
     self->phi   = -1.05;
     self->eta   = -0.273;
     self->gamma = 1.09;
     break;
    case 800:
     self->alpha = 0.393;
     self->beta  = 0.564;
     self->phi   = -1.20;
     self->eta   = -0.278;
     self->gamma = 1.20;
     break;
    case 1200:
     self->alpha = 0.365;
     self->beta  = 0.623;
     self->phi   = -1.26;
     self->eta   = -0.301;
     self->gamma = 1.34;
     break;
    case 1600:
     self->alpha = 0.379;
     self->beta  = 0.637;
     self->phi   = -1.45;
     self->eta   = -0.301;
     self->gamma = 1.50;
     break;
    case 2400:
     self->alpha = 0.355;
     self->beta  = 0.673;
     self->phi   = -1.50;
     self->eta   = -0.319;
     self->gamma = 1.68;
     break;
    case 3200:
     self->alpha = 0.327;
     self->beta  = 0.702;
     self->phi   = -1.49;
     self->eta   = -0.336;
     self->gamma = 1.81;
     break; 
    default:
      g_error ("NcMultiplicityFuncTinkerMeanNormalized: Delta == %u not supported.", self->int_Delta);
    break;
  }
  
}

/**
 * nc_multiplicity_func_tinker_mean_normalized_get_Delta:
 * @mt10: a #NcMultiplicityFuncTinkerMeanNormalized.
 *
 * Returns: the value of #NcMultiplicityFuncTinkerMeanNormalized:Delta property.
 */
static gdouble
_nc_multiplicity_func_tinker_mean_normalized_get_Delta (NcMultiplicityFunc *mulf)
{
  NcMultiplicityFuncTinkerMeanNormalized *mt10 = NC_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED (mulf);
  NcMultiplicityFuncTinkerMeanNormalizedPrivate * const self = mt10->priv;

  return self->Delta;
}


