/***************************************************************************
 *            ncm_calc.c
 *
 *  Mon March 21 15:32:40 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_calc.c
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_calc
 * @title: NcmData
 * @short_description: Abstract class for implementing calculator objects.
 * 
 * Base class describing calculator objects, i.e., objetcs that use
 * one or more #NcmModel to perform calculations.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_calc.h"
#include "math/ncm_model_ctrl.h"

enum
{
  PROP_0,
  PROP_RELTOL,
  PROP_ABSTOL
};

G_DEFINE_TYPE (NcmCalc, ncm_calc, G_TYPE_OBJECT);

static void
ncm_calc_init (NcmCalc *calc)
{
  calc->reltol = 0.0;
  calc->abstol = 0.0;
  calc->ctrl   = g_ptr_array_new ();

  g_ptr_array_set_free_func (calc->ctrl, (GDestroyNotify) ncm_model_ctrl_free);
}

static void
ncm_calc_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmCalc *calc = NCM_CALC (object);
  g_return_if_fail (NCM_IS_CALC (object));

  switch (prop_id)
  {
    case PROP_RELTOL:
      ncm_calc_set_reltol (calc, g_value_get_double (value));
      break;
    case PROP_ABSTOL:
      ncm_calc_set_abstol (calc, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_calc_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmCalc *calc = NCM_CALC (object);
  g_return_if_fail (NCM_IS_CALC (object));

  switch (prop_id)
  {
    case PROP_RELTOL:
      g_value_set_double (value, ncm_calc_get_reltol (calc));
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, ncm_calc_get_abstol (calc));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_calc_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_calc_parent_class)->constructed (object);
  {
    NcmCalc *calc            = NCM_CALC (object);
    NcmCalcClass *calc_class = NCM_CALC_GET_CLASS (calc);
    const guint len = calc_class->dep_list->len;
    guint i;
    
    for (i = 0; i < len; i++)
    {
      NcmModelCtrl *ctrl_i = ncm_model_ctrl_new (NULL);
      g_ptr_array_add (calc->ctrl, ctrl_i);
    }
  }
}

static void
ncm_calc_dispose (GObject *object)
{
  NcmCalc *calc = NCM_CALC (object);

  g_clear_pointer (&calc->ctrl, g_ptr_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_calc_parent_class)->dispose (object);
}

static void
ncm_calc_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_calc_parent_class)->finalize (object);
}

static void
ncm_calc_class_init (NcmCalcClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->constructed  = ncm_calc_constructed;
  object_class->set_property = ncm_calc_set_property;
  object_class->get_property = ncm_calc_get_property;
  object_class->dispose      = ncm_calc_dispose;
  object_class->finalize     = ncm_calc_finalize;

  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance",
                                                        0.0, 1.0, NCM_CALC_DEFAULT_RELTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_ABSTOL,
                                   g_param_spec_double ("abstol",
                                                        NULL,
                                                        "Absolute tolerance",
                                                        0.0, G_MAXDOUBLE, NCM_CALC_DEFAULT_ABSTOL,
                                                        G_PARAM_READABLE | G_PARAM_WRITABLE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

void _ncm_calc_class_cdep (GType *t) { t[0] = G_TYPE_INVALID; }

/**
 * ncm_calc_class_set_num_dep:
 * @calc_class: a #NcmCalcClass
 * @ndep: number of dependencies #NcmModel
 * 
 * Sets the number of dependencies #NcmModel
 * 
 */
void 
ncm_calc_class_set_num_dep (NcmCalcClass *calc_class, guint ndep)
{
  g_assert (NCM_CALC_CLASS (ncm_calc_parent_class)->dep_list == calc_class->dep_list);
  
  if (ndep > NCM_CALC_MAX_DEPS)
    g_error ("ncm_calc_class_set_num_dep: NcmCalc with more than %u dependencies are not supported.", NCM_CALC_MAX_DEPS);
  
  calc_class->dep_list = g_array_new (TRUE, TRUE, sizeof (GType));
  g_array_set_size (calc_class->dep_list, ndep);
  calc_class->ndep = ndep;

  g_array_set_clear_func (calc_class->dep_list, (GDestroyNotify) &_ncm_calc_class_cdep);
}

/**
 * ncm_calc_class_set_dep:
 * @calc_class: a #NcmCalcClass
 * @p: position of the #NcmModel
 * @dep_model: a #GType
 * 
 * Sets the @p-th #NcmModel to be of GType @dep_model.
 * 
 */
void 
ncm_calc_class_set_dep (NcmCalcClass *calc_class, guint p, GType dep_model)
{
  guint i;
  g_assert (g_array_index (calc_class->dep_list, GType, p) == G_TYPE_INVALID);
  g_assert (g_type_is_a (dep_model, NCM_TYPE_MODEL));
  
  for (i = 0; i < calc_class->dep_list->len; i++)
  {
    GType dep_i = g_array_index (calc_class->dep_list, GType, p);
    if (dep_i != G_TYPE_INVALID)
    {
      if (g_type_is_a (dep_i, dep_model) || g_type_is_a (dep_model, dep_i))
        g_error ("ncm_calc_class_set_dep: Model of type `%s', position %u: `%s' was already set.",
                 g_type_name (dep_model), i, g_type_name (dep_i));
    }
  }
}

/**
 * ncm_calc_class_check:
 * @calc_class: a #NcmCalcClass
 * 
 * Checks if all dependencies were set consistently.
 * 
 */
void
ncm_calc_class_check (NcmCalcClass *calc_class)
{
  guint i;
  for (i = 0; i < calc_class->dep_list->len; i++)
  {
    GType dep_i = g_array_index (calc_class->dep_list, GType, i);
    g_assert (g_type_is_a (dep_i, NCM_TYPE_MODEL));
  }
}

/**
 * ncm_calc_prepare_array:
 * @calc: a #NcmCalc
 * @ma: (in) (array zero-terminated=1) (element-type NcmModel): array of #NcmModel
 * 
 * Prepares @calc using the models in @ma.
 * 
 */
void 
ncm_calc_prepare_array (NcmCalc *calc, NcmModel **ma)
{
  NcmCalcClass *calc_class = NCM_CALC_GET_CLASS (calc);
  guint order[NCM_CALC_MAX_DEPS];
  guint mask = 0;
  guint i;

  for (i = 0; i < calc_class->dep_list->len; i++)
  {
    NcmModel *model;
    GType dep_i = g_array_index (calc_class->dep_list, GType, i);
    gboolean found = FALSE;
    NcmModelID main_model = ncm_model_type_main_model (dep_i);
    guint j = 0;

    if (main_model != -1)
    {
      GType submodel_dep_i = dep_i;
      gboolean sub_found = FALSE;

      dep_i = ncm_mset_get_type_by_id (main_model);

      while ((model = ma[j]) != NULL)
      {
        if (g_type_is_a (G_OBJECT_TYPE (model), dep_i))
        {
          const guint submodel_len = ncm_model_get_submodel_len (model);
          guint k;
          found = TRUE;
          for (k = 0; k < submodel_len; k++)
          {
            NcmModel *submodel = ncm_model_peek_submodel (model, k);
            if (g_type_is_a (G_OBJECT_TYPE (submodel), submodel_dep_i))
            {
              sub_found = TRUE;
            }
          }
        }
        j++;
      }
      if (!found)
        g_error ("ncm_calc_prepare_array: cannot find the main NcmModel type `%s' (of required submodel `%s') in the model array.",
                 g_type_name (dep_i), g_type_name (submodel_dep_i));
      if (!sub_found)
        g_error ("ncm_calc_prepare_array: main NcmModel type `%s' found in the model array, but it does not contain the required submodel `%s'.",
                 g_type_name (dep_i), g_type_name (submodel_dep_i));
    }
    else
    {
      while ((model = ma[j]) != NULL)
      {
        guint mask_j = 1 << j;
        if (g_type_is_a (G_OBJECT_TYPE (model), dep_i) && !(mask & mask_j))
        {
          mask = mask | mask_j;
          found = TRUE;
          order[i] = j;
        }      
        j++;
      }
    }

    if (!found)
      g_error ("ncm_calc_prepare_array: cannot find NcmModel type `%s' in the model array.",
               g_type_name (dep_i));
  }

  switch (calc_class->dep_list->len)
  {
    case 0:
      ((NcmCalcPrepare0) calc_class->prepare) (calc);
      break;
    case 1:
      ((NcmCalcPrepare1) calc_class->prepare) (calc, ma[order[0]]);
      break;
    case 2:
      ((NcmCalcPrepare2) calc_class->prepare) (calc, ma[order[0]], ma[order[1]]);
      break;
    case 3:
      ((NcmCalcPrepare3) calc_class->prepare) (calc, ma[order[0]], ma[order[1]], ma[order[2]]);
      break;
    case 4:
      ((NcmCalcPrepare4) calc_class->prepare) (calc, ma[order[0]], ma[order[1]], ma[order[2]], ma[order[3]]);
      break;
    case 5:
      ((NcmCalcPrepare5) calc_class->prepare) (calc, ma[order[0]], ma[order[1]], ma[order[2]], ma[order[3]], ma[order[4]]);
      break;
    case 6:
      ((NcmCalcPrepare6) calc_class->prepare) (calc, ma[order[0]], ma[order[1]], ma[order[2]], ma[order[3]], ma[order[4]], ma[order[5]]);
      break;
    default:
      g_assert_not_reached ();
      break;
  }
}

/**
 * ncm_calc_prepare_if_needed_array:
 * @calc: a #NcmCalc
 * @ma: (in) (array zero-terminated=1) (element-type NcmModel): array of #NcmModel
 * 
 * Prepares @calc using the models in @ma.
 * 
 */
void 
ncm_calc_prepare_if_needed_array (NcmCalc *calc, NcmModel **ma)
{
  

}

/**
 * ncm_calc_prepare_if_needed_vargs: (skip)
 * @calc: a #NcmCalc
 * @...: list of #NcmModel
 * 
 * Prepares @calc using the models in @... .
 * 
 */
void 
ncm_calc_prepare_if_needed_vargs (NcmCalc *calc, ...)
{

}

/**
 * ncm_calc_set_reltol:
 * @calc: a #NcmCalc
 * @reltol: the relative tolerance
 * 
 * Sets the relative tolerance to @reltol.
 * 
 */
void 
ncm_calc_set_reltol (NcmCalc *calc, const gdouble reltol)
{
  calc->reltol = reltol;
}

/**
 * ncm_calc_set_abstol:
 * @calc: a #NcmCalc
 * @abstol: the absolute tolerance
 * 
 * Sets the absolute tolerance to @abstol.
 * 
 */
void 
ncm_calc_set_abstol (NcmCalc *calc, const gdouble abstol)
{
  calc->abstol = abstol;
}

/**
 * ncm_calc_get_reltol:
 * @calc: a #NcmCalc
 * 
 * Gets the relative tolerance.
 * 
 * Returns: the relative tolerance.
 */
gdouble 
ncm_calc_get_reltol (NcmCalc *calc)
{
  return calc->reltol;
}

/**
 * ncm_calc_get_abstol:
 * @calc: a #NcmCalc
 * 
 * Gets the absolute tolerance.
 * 
 * Returns: the relative tolerance.
 */
gdouble 
ncm_calc_get_abstol (NcmCalc *calc)
{
  return calc->abstol;
}
