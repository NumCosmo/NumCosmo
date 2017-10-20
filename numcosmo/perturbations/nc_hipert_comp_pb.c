/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_hipert_comp_pb.c
 *
 *  Fri October 13 11:10:24 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_comp_pb.c
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
 * SECTION:nc_hipert_comp_pb
 * @title: NcHIPertCompPB
 * @short_description: Photon-Baryon plasma compoment
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "perturbations/nc_hipert_comp_pb.h"

struct _NcHIPertCompPBPrivate
{
  gint a;
};

enum
{
  PROP_0,
  PROP_BLA
};

G_DEFINE_TYPE (NcHIPertCompPB, nc_hipert_comp_pb, NC_TYPE_HIPERT_COMP);

static void
nc_hipert_comp_pb_init (NcHIPertCompPB *pb)
{
  pb->priv = G_TYPE_INSTANCE_GET_PRIVATE (pb, NC_TYPE_HIPERT_COMP_PB, NcHIPertCompPBPrivate);
}

static void
_nc_hipert_comp_pb_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_HIPERT_COMP_PB (object));

  switch (prop_id)
  {
    case PROP_BLA:
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_comp_pb_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_HIPERT_COMP_PB (object));

  switch (prop_id)
  {
    case PROP_BLA:
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_comp_pb_dispose (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_comp_pb_parent_class)->dispose (object);
}

static void
_nc_hipert_comp_pb_finalize (GObject *object)
{
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_comp_pb_parent_class)->finalize (object);
}

NC_HIPERT_BG_VAR_ID_FUNC_IMPL (nc_hipert_comp_pb, NcHIPertCompPB);

static void
nc_hipert_comp_pb_class_init (NcHIPertCompPBClass *klass)
{
  GObjectClass* object_class    = G_OBJECT_CLASS (klass);

  g_type_class_add_private (klass, sizeof (NcHIPertCompPBPrivate));

  object_class->set_property = &_nc_hipert_comp_pb_set_property;
  object_class->get_property = &_nc_hipert_comp_pb_get_property;
  object_class->dispose      = &_nc_hipert_comp_pb_dispose;
  object_class->finalize     = &_nc_hipert_comp_pb_finalize;

  g_object_class_install_property (object_class,
                                   PROP_BLA,
                                   g_param_spec_int ("bla",
                                                     NULL,
                                                     "ha",
                                                     G_MININT,
                                                     G_MAXINT,
                                                     0,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  nc_hipert_bg_var_class_register_id ("NcHIPertCompPB", 
                                      "First order Photon-Baryons background variables", 
                                      NULL,
                                      0);
}

/**
 * nc_hipert_comp_pb_new:
 *
 * Creates a new #NcHIPertCompPB.
 *
 * Returns: (transfer full): the newly created #NcHIPertCompPB.
 */
NcHIPertCompPB *
nc_hipert_comp_pb_new (void)
{
  NcHIPertCompPB *pb = g_object_new (NC_TYPE_HIPERT_COMP_PB,
                                     NULL);

  return pb;
}

/**
 * nc_hipert_comp_pb_ref:
 * @pb: a #NcHIPertCompPB
 *
 * Increases the reference count of @pb.
 *
 * Returns: (transfer full): @pb.
 */
NcHIPertCompPB *
nc_hipert_comp_pb_ref (NcHIPertCompPB *pb)
{
  return g_object_ref (pb);
}

/**
 * nc_hipert_comp_pb_free:
 * @pb: a #NcHIPertCompPB
 *
 * Decreases the reference count of @pb.
 *
 */
void 
nc_hipert_comp_pb_free (NcHIPertCompPB *pb)
{
  g_object_unref (pb);
}

/**
 * nc_hipert_comp_pb_clear:
 * @pb: a #NcHIPertCompPB
 *
 * Decreases the reference count of *@pb and sets the pointer *@pb to NULL.
 *
 */
void 
nc_hipert_comp_pb_clear (NcHIPertCompPB **pb)
{
  g_clear_object (pb);
}

