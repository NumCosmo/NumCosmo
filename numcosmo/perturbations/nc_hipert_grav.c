/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_hipert_grav.c
 *
 *  Thu October 12 14:32:22 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_grav.c
 * Copyright (C) 2017 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_hipert_grav
 * @title: NcHIPertComp
 * @short_description: Abstract class describing a general first order gravitation theory.
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "perturbations/nc_hipert_grav.h"
#include "nc_enum_types.h"

struct _NcHIPertGravPrivate
{
  gint a;
};

enum
{
  PROP_0,
  PROP_GAUGE
};

G_DEFINE_TYPE (NcHIPertGrav, nc_hipert_grav, G_TYPE_OBJECT);

static void
nc_hipert_grav_init (NcHIPertGrav *grav)
{
  grav->priv = G_TYPE_INSTANCE_GET_PRIVATE (grav, NC_TYPE_HIPERT_GRAV, NcHIPertGravPrivate);
}

static void
_nc_hipert_grav_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_HIPERT_GRAV (object));

  switch (prop_id)
  {
    case PROP_GAUGE:
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_grav_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_HIPERT_GRAV (object));

  switch (prop_id)
  {
    case PROP_GAUGE:
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_grav_dispose (GObject *object)
{

  G_OBJECT_CLASS (nc_hipert_grav_parent_class)->dispose (object);
}

static void
_nc_hipert_grav_finalize (GObject *object)
{

  G_OBJECT_CLASS (nc_hipert_grav_parent_class)->finalize (object);
}

static void
nc_hipert_grav_class_init (NcHIPertGravClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  g_type_class_add_private (klass, sizeof (NcHIPertGravPrivate));

  object_class->set_property = &_nc_hipert_grav_set_property;
  object_class->get_property = &_nc_hipert_grav_get_property;
  object_class->dispose      = &_nc_hipert_grav_dispose;
  object_class->finalize     = &_nc_hipert_grav_finalize;

  g_object_class_install_property (object_class,
                                   PROP_GAUGE,
                                   g_param_spec_enum ("gauge",
                                                      NULL,
                                                      "gauge",
                                                      NC_TYPE_HIPERT_GRAV_GAUGE,
                                                      NC_HIPERT_GRAV_GAUGE_SYNCHRONOUS,
                                                      G_PARAM_READABLE | G_PARAM_WRITABLE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}


