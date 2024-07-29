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
 *
 * FIXME
 * Reference: arXiv:1212.0095v4
 */ 

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_multiplicity_func_watson.h"

struct _NcMultiplicityFuncWatsonPrivate
{
  NcMultiplicityFuncMassDef mdef;
  gdouble (*eval) (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z); 
  gdouble Delta;
};

enum
{
  PROP_0,
  PROP_SIZE
};

G_DEFINE_TYPE_WITH_PRIVATE (NcMultiplicityFuncWatson, nc_multiplicity_func_watson, NC_TYPE_MULTIPLICITY_FUNC)

static void
nc_multiplicity_func_watson_init (NcMultiplicityFuncWatson *mwat)
{
  NcMultiplicityFuncWatsonPrivate * const self = mwat->priv = nc_multiplicity_func_watson_get_instance_private (mwat);

  self->mdef    = NC_MULTIPLICITY_FUNC_MASS_DEF_LEN;
}

static void
_nc_multiplicity_func_watson_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  /*NcMultiplicityFuncWatson *mwat = NC_MULTIPLICITY_FUNC_WATSON (object);*/
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
  /*NcMultiplicityFuncWatson *mwat = NC_MULTIPLICITY_FUNC_WATSON (object);*/
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
static void _nc_multiplicity_func_watson_set_Delta (NcMultiplicityFunc *mulf, gdouble Delta); 
static NcMultiplicityFuncMassDef _nc_multiplicity_func_watson_get_mdef (NcMultiplicityFunc *mulf);
static gdouble _nc_multiplicity_func_watson_get_Delta (NcMultiplicityFunc *mulf); 
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
  parent_class->set_Delta = &_nc_multiplicity_func_watson_set_Delta;
  parent_class->get_mdef = &_nc_multiplicity_func_watson_get_mdef;
  parent_class->get_Delta = &_nc_multiplicity_func_watson_get_Delta;
  parent_class->eval     = &_nc_multiplicity_func_watson_eval;
    
}

static gdouble
_nc_multiplicity_func_watson_fof_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)
{
  /*NcMultiplicityFuncWatson *mwat = NC_MULTIPLICITY_FUNC_WATSON (mulf); */
  /*NcMultiplicityFuncWatsonPrivate * const self = mwat->priv; */
  const gdouble A     = 0.282;
  const gdouble alpha = 2.163;
  const gdouble beta  = 1.406;
  const gdouble gamma = 1.210;
  
  gdouble   f_Watson = A * ( pow(beta / sigma, alpha) + 1.0) * exp (- gamma / (sigma * sigma) );
  
  NCM_UNUSED (mulf);
  NCM_UNUSED (cosmo);
  NCM_UNUSED (z);

  return f_Watson;
}

static gdouble
_nc_multiplicity_func_watson_mean_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)
{
  NcMultiplicityFuncWatson *mwat = NC_MULTIPLICITY_FUNC_WATSON (mulf);
  NcMultiplicityFuncWatsonPrivate * const self = mwat->priv; 
  const gdouble Omega_m = nc_hicosmo_E2Omega_m (cosmo, z);
  const gdouble Delta_178 = self->Delta / 178.0;
  
  gdouble A, alpha, beta, gamma;
  gdouble f_Watson;
   
  
  if (z == 0)
  {
	A = 0.194;
	alpha = 1.805;
	beta = 2.267;
	gamma = 1.287;
  }
  
  else if (z>=6)
  {
	A = 0.563;
	alpha = 3.810;
	beta = 0.874;
	gamma = 1.453;
  }
  
  else
  {
	A = Omega_m * (1.097 * pow(1 + z, -3.216) + 0.074);
	alpha = Omega_m * (5.907 * pow(1 + z, -3.058) + 2.349) ;
	beta = Omega_m * (3.136 * pow(1 + z, -3.599) + 2.344);
	gamma = 1.318;
	  
  }

  if (self->Delta == 178) 
  {
    f_Watson = A * ( pow(beta / sigma, alpha) + 1.0) * exp (- gamma / (sigma * sigma) );
 	
  }
  
  else
  {
	gdouble f_Watson_178 = A * ( pow(beta / sigma, alpha) + 1.0) * exp (- gamma / (sigma * sigma) );
 		
 	const gdouble C = exp ( 0.023 * ( Delta_178  - 1.0 ) );
	const gdouble d = - 0.456 * Omega_m - 0.139;
	const gdouble p = 0.072;
	const gdouble q = 2.130;
 		 	
	const gdouble Gamma = (C * pow(Delta_178, d)* exp (p * (1 - Delta_178)/ pow( sigma, q)));
 		  
	f_Watson = f_Watson_178 * Gamma;
	
  }
	
	return f_Watson;
  
}

static void 
_nc_multiplicity_func_watson_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef)
{
  NcMultiplicityFuncWatson *mwat = NC_MULTIPLICITY_FUNC_WATSON (mulf);
  NcMultiplicityFuncWatsonPrivate * const self = mwat->priv;

  switch (mdef)
  {
    case NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN:
      self->eval = &_nc_multiplicity_func_watson_mean_eval;
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL:
      g_error ("NcMultiplicityFuncWatson does not support critical mass def");
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_VIRIAL:
      g_error ("NcMultiplicityFuncWatson does not support virial mass def");
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_FOF:
      self->eval = &_nc_multiplicity_func_watson_fof_eval;
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
_nc_multiplicity_func_watson_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)   /* $f(\sigma)$ watson: astro-ph/0803.2706 */
{
  NcMultiplicityFuncWatson *mwat = NC_MULTIPLICITY_FUNC_WATSON (mulf);
  NcMultiplicityFuncWatsonPrivate * const self = mwat->priv;
  
  return self->eval(mulf, cosmo, sigma, z);
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

/**
 * nc_multiplicity_func_watson_set_Delta:
 * @mwat: a #NcMultiplicityFuncWatson.
 * @Delta: value of #NcMultiplicityFuncWatson:Delta.
 *
 * Sets the value @Delta to the #NcMultiplicityFuncWatson:Delta property.
 *
 */
void
_nc_multiplicity_func_watson_set_Delta (NcMultiplicityFunc *mulf, gdouble Delta)
{
  NcMultiplicityFuncWatson *mwat = NC_MULTIPLICITY_FUNC_WATSON (mulf);
  NcMultiplicityFuncWatsonPrivate * const self = mwat->priv;

  g_assert (Delta >= 0);

  self->Delta = Delta;
}

/**
 * nc_multiplicity_func_watson_get_Delta:
 * @mwat: a #NcMultiplicityFuncWatson.
 *
 * Returns: the value of #NcMultiplicityFuncWatson:Delta property.
 */
gdouble
_nc_multiplicity_func_watson_get_Delta (NcMultiplicityFunc *mulf)
{
  NcMultiplicityFuncWatson *mwat = NC_MULTIPLICITY_FUNC_WATSON (mulf);
  NcMultiplicityFuncWatsonPrivate * const self = mwat->priv;
  
  return self->Delta;
}


