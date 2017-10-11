/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_hipert_comp.c
 *
 *  Wed October 11 15:54:13 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_comp.c
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
 * SECTION:nc_hipert_comp
 * @title: NcHIPertComp
 * @short_description: Abstract class describing a general perturbation compoment.
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "perturbations/nc_hipert_comp.h"

struct _NcHIPertCompPrivate
{
  gint a;
};

enum
{
  PROP_0,
  PROP_DIM
};

G_DEFINE_TYPE (NcHIPertComp, nc_hipert_comp, G_TYPE_OBJECT);

static void
nc_hipert_comp_init (NcHIPertComp *comp)
{
  comp->priv = G_TYPE_INSTANCE_GET_PRIVATE (nc_hipert_comp, NC_TYPE_HIPERT_COMP, NcHIPertCompPrivate);

}

static void
nc_hipert_comp_finalize (GObject *object)
{

  G_OBJECT_CLASS (nc_hipert_comp_parent_class)->finalize (object);
}

static void
nc_hipert_comp_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_HIPERT_COMP (object));

  switch (prop_id)
  {
    case PROP_DIM:
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hipert_comp_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_HIPERT_COMP (object));

  switch (prop_id)
  {
    case PROP_DIM:
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hipert_comp_class_init (NcHIPertCompClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  g_type_class_add_private (klass, sizeof (NcHIPertCompPrivate));

  object_class->finalize = nc_hipert_comp_finalize;
  object_class->set_property = nc_hipert_comp_set_property;
  object_class->get_property = nc_hipert_comp_get_property;

  g_object_class_install_property (object_class,
                                   PROP_DIM,
                                   g_param_spec_uint ("dim",
                                                      NULL,
                                                      "dimension",
                                                      0, 
                                                      G_MAXUINT,
                                                      2,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}


