/***************************************************************************
 *            nc_hipert_first_order.c
 *
 *  Mon October 09 16:58:16 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_first_order.c
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
 * SECTION:nc_hipert_first_order
 * @title: NcHIPertFirstOrder
 * @short_description: Base object for implementing first order perturbation in a Friedmann background.
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "perturbations/nc_hipert_first_order.h"

struct _NcHIPertFirstOrderPrivate
{
  gint a;
};

enum
{
  PROP_0,
  PROP_NCOMP
};

G_DEFINE_TYPE (NcHIPertFirstOrder, nc_hipert_first_order, NC_TYPE_HIPERT_BOLTZMANN);

static void
nc_hipert_first_order_init (NcHIPertFirstOrder *nc_hipert_first_order)
{
  nc_hipert_first_order->priv = G_TYPE_INSTANCE_GET_PRIVATE (nc_hipert_first_order, NC_TYPE_HIPERT_FIRST_ORDER, NcHIPertFirstOrderPrivate);
}

static void
_nc_hipert_first_order_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_HIPERT_FIRST_ORDER (object));

  switch (prop_id)
  {
    case PROP_NCOMP:
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_first_order_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_HIPERT_FIRST_ORDER (object));

  switch (prop_id)
  {
    case PROP_NCOMP:
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_first_order_dispose (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_first_order_parent_class)->dispose (object);
}

static void
_nc_hipert_first_order_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_first_order_parent_class)->finalize (object);
}

static void
nc_hipert_first_order_class_init (NcHIPertFirstOrderClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  g_type_class_add_private (klass, sizeof (NcHIPertFirstOrderPrivate));

  object_class->set_property = &_nc_hipert_first_order_set_property;
  object_class->get_property = &_nc_hipert_first_order_get_property;
  object_class->dispose      = &_nc_hipert_first_order_dispose;
  object_class->finalize     = &_nc_hipert_first_order_finalize;

  g_object_class_install_property (object_class,
                                   PROP_NCOMP,
                                   g_param_spec_uint ("ncomp",
                                                      NULL,
                                                      "Number of components",
                                                      0, 
                                                      G_MAXUINT,
                                                      1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}
