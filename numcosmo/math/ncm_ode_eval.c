/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_ode_eval.c
 *
 *  Thu December 13 11:13:24 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_ode_eval.c
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_ode_eval
 * @title: NcmODEEval
 * @short_description: Abstract class for ODE system evaluation
 *
 * This class implement an abstract interface between the ODE system
 * and the evaluation of $\mathrm{d}f$ and $J$. 
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_ode_eval.h"

struct _NcmODEEvalPrivate
{
  guint sys_size;
};

enum
{
  PROP_0,
  PROP_SYS_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmODEEval, ncm_ode_eval, G_TYPE_OBJECT);

static void
ncm_ode_eval_init (NcmODEEval *ode_eval)
{
  NcmODEEvalPrivate * const self = ode_eval->priv = ncm_ode_eval_get_instance_private (ode_eval);

  self->sys_size = 0;
}

static void
_ncm_ode_eval_finalize (GObject *object)
{

  
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_ode_eval_parent_class)->finalize (object);
}

static void
_ncm_ode_eval_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NCM_IS_ODE_EVAL (object));

  switch (prop_id)
  {
    case PROP_SYS_SIZE:
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_ode_eval_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NCM_IS_ODE_EVAL (object));

  switch (prop_id)
  {
    case PROP_SYS_SIZE:
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static gint _ncm_ode_eval_df (NcmODEEval *ode_eval, const guint sys_size, const gdouble t, const gdouble * restrict f, gdouble * restrict df) { g_error ("_ncm_ode_eval_df: not implemented by `%s'.", G_OBJECT_TYPE_NAME (ode_eval)); return 1; }
static gint _ncm_ode_eval_J_dense (NcmODEEval *ode_eval, const guint sys_size, const gdouble t, const gdouble * restrict f, gdouble ** restrict J_col) { g_error ("_ncm_ode_eval_J_dense: not implemented by `%s'.", G_OBJECT_TYPE_NAME (ode_eval)); return 1; }

static void
ncm_ode_eval_class_init (NcmODEEvalClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_ode_eval_set_property;
  object_class->get_property = &_ncm_ode_eval_get_property;
  object_class->finalize     = &_ncm_ode_eval_finalize;

  g_object_class_install_property (object_class,
                                   PROP_SYS_SIZE,
                                   g_param_spec_uint ("sys-size",
                                                      NULL,
                                                      "ODE system size",
                                                      1, G_MAXUINT, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->df      = &_ncm_ode_eval_df;
  klass->J_dense = &_ncm_ode_eval_J_dense;
}

/**
 * ncm_ode_eval_ref:
 * @ode_eval: a #NcmODEEval
 *
 * Increase the reference of @ode_eval by one.
 *
 * Returns: (transfer full): @ode_eval.
 */
NcmODEEval *
ncm_ode_eval_ref (NcmODEEval *ode_eval)
{
  return g_object_ref (ode_eval);
}

/**
 * ncm_ode_eval_free:
 * @ode_eval: a #NcmODEEval
 *
 * Decrease the reference count of @ode_eval by one.
 *
 */
void
ncm_ode_eval_free (NcmODEEval *ode_eval)
{
  g_object_unref (ode_eval);
}

/**
 * ncm_ode_eval_clear:
 * @ode_eval: a #NcmODEEval
 *
 * Decrease the reference count of @ode_eval by one, and sets the pointer *@ode_eval to
 * NULL.
 *
 */
void
ncm_ode_eval_clear (NcmODEEval **ode_eval)
{
  g_clear_object (ode_eval);
}

/**
 * ncm_ode_eval_df: (virtual df)
 * @ode_eval: a #NcmODEEval
 * @sys_size: ODE system size 
 * @t: the current time $t$
 * @f: (array length=sys_size): ODE system current state $f$
 * @df: (inout) (array length=sys_size): Vector to hold the time derivatives $\mathrm{d}f$
 * 
 * Computes the time derivatives of the ODE system in @df using the 
 * current state in @f.
 * 
 * Return: status
 */
/**
 * ncm_ode_eval_J_dense: (virtual J_dense)
 * @ode_eval: a #NcmODEEval
 * @sys_size: ODE system size 
 * @t: the current time $t$
 * @f: (array length=sys_size): ODE system current state $f$
 * @J_col: (array length=sys_size): Array containing the Jacobian columns 
 *
 * Computes the jacobian matrix $J$ of the ODE system in @J_col using the 
 * current state in @f.
 * 
 * Return: status
 */
