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
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_multiplicity_func_tinker_mean_normalized.h"

G_DEFINE_TYPE (NcMultiplicityFuncTinkerMeanNormalized, nc_multiplicity_func_tinker_mean_normalized, NC_TYPE_MULTIPLICITY_FUNC);

enum
{
  PROP_0,
  PROP_DELTA
};

/**
 * nc_multiplicity_func_tinker_mean_normalized_new:
 * @Delta: FIXME 
 *   
 * FIXME
 *
 * Returns: A new #NcMultiplicityFunc.
 */
NcMultiplicityFunc *
nc_multiplicity_func_tinker_mean_normalized_new (gdouble Delta)
{
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED,
                       "Delta", Delta,
                       NULL);
}

static gdouble
_nc_multiplicity_func_tinker_mean_normalized_eval (NcMultiplicityFunc *mulf, NcHICosmo *model, gdouble sigma, gdouble z)   /* $g(\sigma) = \nu x f(\nu)$ Tinker: Eq. 8 1001.3162 */
{
  NcMultiplicityFuncTinkerMeanNormalized *mtmn = NC_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED (mulf);
  const gdouble nu = 1.686 / sigma;
  gdouble f_Tinker_mean_normalized;
  
  NCM_UNUSED (model);

  if (mtmn->Delta == 200.0)
  {
    const gdouble alpha_z = mtmn->alpha;
    const gdouble beta_z  = mtmn->beta * pow(1.0 + z, 0.2);
    const gdouble phi_z   = mtmn->phi * pow(1.0 + z, -0.08);
    const gdouble eta_z   = mtmn->eta * pow(1.0 + z, 0.27);
    const gdouble gamma_z = mtmn->gamma * pow(1.0 + z, -0.01);
    f_Tinker_mean_normalized = alpha_z * (1.0 + pow(beta_z * nu, -2.0 * phi_z)) 
                               * pow(nu, 2.0 * eta_z) * exp(-gamma_z * nu * nu / 2.0) * nu;
  }    
  else 
  {
    f_Tinker_mean_normalized = nu * mtmn->alpha * (1.0 + pow(mtmn->beta * nu, -2.0 * mtmn->phi)) 
                               * pow(nu, 2.0 * mtmn->eta) * exp(-mtmn->gamma * nu * nu / 2.0);
  }
  
    return f_Tinker_mean_normalized;
}

/**
 * nc_multiplicity_func_tinker_mean_normalized_set_Delta:
 * @mtmn: a #NcMultiplicityFuncTinkerMeanNormalized.
 * @Delta: value of #NcMultiplicityFuncTinkerMeanNormalized:Delta.
 *
 * Sets the value @Delta to the #NcMultiplicityFuncTinkerMeanNormalized:Delta property.
 *
 */
void
nc_multiplicity_func_tinker_mean_normalized_set_Delta (NcMultiplicityFuncTinkerMeanNormalized *mtmn, gdouble Delta)
{
  const guint int_Delta = mtmn->Delta;
  g_assert (Delta >= 0);
  g_assert (Delta == int_Delta);
  
  mtmn->Delta = Delta;
  mtmn->int_Delta = int_Delta;

  switch (mtmn->int_Delta)
  {
    case 200:
      mtmn->alpha =  0.368;
      mtmn->beta  =  0.589;
      mtmn->phi   =  -0.729;
      mtmn->eta   = -0.243;
      mtmn->gamma = 0.864;
      break;
    case 300:
      mtmn->alpha = 0.363;
      mtmn->beta  = 0.585;
      mtmn->phi   = -0.789;
      mtmn->eta   = -0.261;
      mtmn->gamma = 0.922;
      break;
    case 400:
      mtmn->alpha = 0.385;
      mtmn->beta  = 0.544;
      mtmn->phi   = -0.910;
      mtmn->eta   = -0.261;
      mtmn->gamma = 0.987;
      break;  
    case 600:
     mtmn->alpha = 0.389;
     mtmn->beta  = 0.543;
     mtmn->phi   = -1.05;
     mtmn->eta   = -0.273;
     mtmn->gamma = 1.09;
     break;
    case 800:
     mtmn->alpha = 0.393;
     mtmn->beta  = 0.564;
     mtmn->phi   = -1.20;
     mtmn->eta   = -0.278;
     mtmn->gamma = 1.20;
     break;
    case 1200:
     mtmn->alpha = 0.365;
     mtmn->beta  = 0.623;
     mtmn->phi   = -1.26;
     mtmn->eta   = -0.301;
     mtmn->gamma = 1.34;
     break;
    case 1600:
     mtmn->alpha = 0.379;
     mtmn->beta  = 0.637;
     mtmn->phi   = -1.45;
     mtmn->eta   = -0.301;
     mtmn->gamma = 1.50;
     break;
    case 2400:
     mtmn->alpha = 0.355;
     mtmn->beta  = 0.673;
     mtmn->phi   = -1.50;
     mtmn->eta   = -0.319;
     mtmn->gamma = 1.68;
     break;
    case 3200:
     mtmn->alpha = 0.327;
     mtmn->beta  = 0.702;
     mtmn->phi   = -1.49;
     mtmn->eta   = -0.336;
     mtmn->gamma = 1.81;
     break; 
    default:
      g_error ("NcMultiplicityFuncTinkerMeanNormalized: Delta == %u not supported.", mtmn->int_Delta);
    break;
  }
  
}

/**
 * nc_multiplicity_func_tinker_mean_normalized_get_Delta:
 * @mtmn: a #NcMultiplicityFuncTinkerMeanNormalized.
 *
 * Returns: the value of #NcMultiplicityFuncTinkerMeanNormalized:Delta property.
 */
gdouble
nc_multiplicity_func_tinker_mean_normalized_get_Delta (const NcMultiplicityFuncTinkerMeanNormalized *mtmn)
{
  return mtmn->Delta;
}

static void
nc_multiplicity_func_tinker_mean_normalized_init (NcMultiplicityFuncTinkerMeanNormalized *mtmn)
{
  mtmn->Delta     = 200.0;
  mtmn->int_Delta = 200;
}

static void
_nc_multiplicity_func_tinker_mean_normalized_finalize (GObject *object)
{
  /* TODO: Add deinitalization code here */

  G_OBJECT_CLASS (nc_multiplicity_func_tinker_mean_normalized_parent_class)->finalize (object);
}

static void
_nc_multiplicity_func_tinker_mean_normalized_set_property (GObject * object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcMultiplicityFuncTinkerMeanNormalized *mtmn = NC_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED (object));

  switch (prop_id)
  {
    case PROP_DELTA:
      nc_multiplicity_func_tinker_mean_normalized_set_Delta (mtmn, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_tinker_mean_normalized_get_property (GObject * object, guint prop_id, GValue * value, GParamSpec * pspec)
{
  NcMultiplicityFuncTinkerMeanNormalized *mtmn = NC_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED (object));

  switch (prop_id)
  {
    case PROP_DELTA:
      g_value_set_double (value, mtmn->Delta);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_multiplicity_func_tinker_mean_normalized_class_init (NcMultiplicityFuncTinkerMeanNormalizedClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass* parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);

  parent_class->eval = &_nc_multiplicity_func_tinker_mean_normalized_eval;
                      
  object_class->finalize = _nc_multiplicity_func_tinker_mean_normalized_finalize;
  object_class->set_property = _nc_multiplicity_func_tinker_mean_normalized_set_property;
  object_class->get_property = _nc_multiplicity_func_tinker_mean_normalized_get_property;

  /**
   * NcMultiplicityFuncTinkerMeanNormalized:Delta:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_DELTA,
                                   g_param_spec_double ("Delta",
                                                        NULL,
                                                        "Delta",
                                                        200.0, 3200.0, 200.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}


