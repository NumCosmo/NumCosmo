/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_ode.c
 *
 *  Wed December 12 11:03:10 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_ode.c
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
 * SECTION:ncm_ode
 * @title: NcmODE
 * @short_description: Abstract class for ODE solvers
 *
 * This class determine the methods needed to implement a ODE solver.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_ode.h"

struct _NcmODEPrivate
{
  guint sys_size;
};

enum
{
  PROP_0,
  PROP_SYS_SIZE
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmODE, ncm_ode, G_TYPE_OBJECT);

static void
ncm_ode_init (NcmODE *ode)
{
  NcmODEPrivate * const self = ode->priv = ncm_ode_get_instance_private (ode);

  self->sys_size = 0;
}

static void
_ncm_ode_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmODE *ode = NCM_ODE (object);  
  g_return_if_fail (NCM_IS_ODE (object));

  switch (prop_id)
  {
    case PROP_SYS_SIZE:
      ncm_ode_set_sys_size (ode, g_value_get_uint (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_ode_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmODE *ode = NCM_ODE (object);
  g_return_if_fail (NCM_IS_ODE (object));

  switch (prop_id)
  {
    case PROP_SYS_SIZE:
      g_value_set_uint (value, ncm_ode_get_sys_size (ode));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_ode_dispose (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_ode_parent_class)->dispose (object);
}

static void
_ncm_ode_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_ode_parent_class)->finalize (object);
}

void _ncm_ode_set_sys_size (NcmODE *ode, guint sys_size) { g_error ("_ncm_ode_set_sys_size: not implemented by `%s'.", G_OBJECT_TYPE_NAME (ode)); }

static void
ncm_ode_class_init (NcmODEClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_ode_set_property;
  object_class->get_property = &_ncm_ode_get_property;
  object_class->dispose      = &_ncm_ode_dispose;
  object_class->finalize     = &_ncm_ode_finalize;

  g_object_class_install_property (object_class,
                                   PROP_SYS_SIZE,
                                   g_param_spec_uint ("sys-size",
                                                      NULL,
                                                      "ODE system size",
                                                      1, G_MAXUINT, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->set_sys_size = &_ncm_ode_set_sys_size;
}

/**
 * ncm_ode_ref:
 * @ode: a #NcmODE
 *
 * Increase the reference of @ode by one.
 *
 * Returns: (transfer full): @ode.
 */
NcmODE *
ncm_ode_ref (NcmODE *ode)
{
  return g_object_ref (ode);
}

/**
 * ncm_ode_free:
 * @ode: a #NcmODE
 *
 * Decrease the reference count of @ode by one.
 *
 */
void
ncm_ode_free (NcmODE *ode)
{
  g_object_unref (ode);
}

/**
 * ncm_ode_clear:
 * @ode: a #NcmODE
 *
 * Decrease the reference count of @ode by one, and sets the pointer *@ode to
 * NULL.
 *
 */
void
ncm_ode_clear (NcmODE **ode)
{
  g_clear_object (ode);
}

/**
 * ncm_ode_set_sys_size: (virtual set_sys_size)
 * @ode: a #NcmODE
 * @sys_size: system size
 *
 * Sets the ODE system size to @sys_size.
 * 
 */
void 
ncm_ode_set_sys_size (NcmODE *ode, guint sys_size)
{
  NcmODEPrivate * const self = ode->priv;

  self->sys_size = sys_size;

  NCM_ODE_GET_CLASS (ode)->set_sys_size (ode, sys_size);
}

/**
 * ncm_ode_get_sys_size:
 * @ode: a #NcmODE
 *
 * Gets the current ODE system size.
 * 
 * Returns: current ODE system size.
 */
guint 
ncm_ode_get_sys_size (NcmODE *ode)
{
  NcmODEPrivate * const self = ode->priv;
  return self->sys_size;
}
