/***************************************************************************
 *            nc_powspec_ml.c
 *
 *  Thu February 18 12:32:13 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_powspec_ml.c
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
 * SECTION:nc_powspec_ml
 * @title: NcPowspecML
 * @short_description: Abstrac class for linear matter power spectrum implementation.
 *
 * This module comprises the set of functions to compute a power spectrum and
 * derived quantities.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_powspec_ml.h"

G_DEFINE_ABSTRACT_TYPE (NcPowspecML, nc_powspec_ml, NCM_TYPE_POWSPEC);

static void
nc_powspec_ml_init (NcPowspecML *nc_powspec_ml)
{
}

static void
nc_powspec_ml_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_powspec_ml_parent_class)->finalize (object);
}

static void
nc_powspec_ml_class_init (NcPowspecMLClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->finalize = nc_powspec_ml_finalize;
}
