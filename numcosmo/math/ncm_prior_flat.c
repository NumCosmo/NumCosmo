/***************************************************************************
 *            ncm_prior_flat.c
 *
 *  Wed August 03 16:58:19 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_prior_flat.c
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * SECTION:ncm_prior_flat
 * @title: NcmPriorFlat
 * @short_description: A flat prior for NcmLikelihood
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_prior_flat.h"

enum
{
  PROP_0,
  PROP_X_LOW,
  PROP_X_UPP,
  PROP_S,
  PROP_VARIABLE,
};

G_DEFINE_TYPE (NcmPriorFlat, ncm_prior_flat, NCM_TYPE_PRIOR)

static void
ncm_prior_flat_init (NcmPriorFlat *pf)
{
  pf->x_low = 0.0;
  pf->x_upp = 0.0;
  pf->s     = 0.0;
}

static void
_ncm_prior_flat_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmPriorFlat *pf = NCM_PRIOR_FLAT (object);
  g_return_if_fail (NCM_IS_PRIOR_FLAT (object));

  switch (prop_id)
  {
    case PROP_X_LOW:
      pf->x_low = g_value_get_double (value);
      break;
    case PROP_X_UPP:
      pf->x_upp = g_value_get_double (value);
      break;
    case PROP_S:
      pf->s = g_value_get_double (value);
      break;
    case PROP_VARIABLE:
      pf->var = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_prior_flat_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmPriorFlat *pf = NCM_PRIOR_FLAT (object);
  g_return_if_fail (NCM_IS_PRIOR_FLAT (object));

  switch (prop_id)
  {
    case PROP_X_LOW:
      g_value_set_double (value, pf->x_low);
      break;
    case PROP_X_UPP:
      g_value_set_double (value, pf->x_upp);
      break;
    case PROP_S:
      g_value_set_double (value, pf->s);
      break;
    case PROP_VARIABLE:
      g_value_set_double (value, pf->var);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_prior_flat_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_prior_flat_parent_class)->finalize (object);
}

static void _ncm_prior_flat_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res);

static void
ncm_prior_flat_class_init (NcmPriorFlatClass *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcmMSetFuncClass *mset_func_class = NCM_MSET_FUNC_CLASS (klass);

  object_class->set_property = &_ncm_prior_flat_set_property;
  object_class->get_property = &_ncm_prior_flat_get_property;
  object_class->finalize     = &_ncm_prior_flat_finalize;

  g_object_class_install_property (object_class,
                                   PROP_X_LOW,
                                   g_param_spec_double ("x-low",
                                                        NULL,
                                                        "lower limit",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_X_UPP,
                                   g_param_spec_double ("x-upp",
                                                        NULL,
                                                        "upper limit",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_S,
                                   g_param_spec_double ("scale",
                                                        NULL,
                                                        "border scale",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_VARIABLE,
                                   g_param_spec_double ("variable",
                                                        NULL,
                                                        "variable",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  NCM_PRIOR_CLASS (klass)->is_m2lnL = FALSE;
  mset_func_class->eval        = &_ncm_prior_flat_eval;
}

static void 
_ncm_prior_flat_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res)
{
  NcmPriorFlat *pf       = NCM_PRIOR_FLAT (func);
  const gdouble huge_exp = 20.0;
  const gdouble mean     = NCM_PRIOR_FLAT_GET_CLASS (pf)->mean (pf, mset);

  res[0] = 
    exp (2.0 * huge_exp / pf->s * ((pf->x_low - mean) + pf->s / 2.0)) +
    exp (2.0 * huge_exp / pf->s * ((mean - pf->x_upp) + pf->s / 2.0))
    ;
}

/**
 * ncm_prior_flat_ref:
 * @pf: a #NcmPriorFlat
 * 
 * Increases the reference count of @pf atomically.
 * 
 * Returns: (transfer full): @pf.
 */
NcmPriorFlat *
ncm_prior_flat_ref (NcmPriorFlat *pf)
{
  return g_object_ref (pf);
}

/**
 * ncm_prior_flat_free:
 * @pf: a #NcmPriorFlat
 * 
 * Decreases the reference count of @pf atomically.
 * 
 */
void 
ncm_prior_flat_free (NcmPriorFlat *pf)
{
  g_object_unref (pf);
}

/**
 * ncm_prior_flat_clear:
 * @pf: a #NcmPriorFlat
 * 
 * Decreases the reference count of *@pf and sets *@pf to NULL.
 * 
 */
void 
ncm_prior_flat_clear (NcmPriorFlat **pf)
{
  g_clear_object (pf);
}
