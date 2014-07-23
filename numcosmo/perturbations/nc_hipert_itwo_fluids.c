/***************************************************************************
 *            nc_hipert_itwo_fluids.c
 *
 *  Tue July 22 17:36:57 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_itwo_fluids.c
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_hipert_itwo_fluids
 * @title: Two Fluids Perturbation Interface
 * @short_description: Perturbation interface for two fluids system 
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_hipert_itwo_fluids.h"

G_DEFINE_INTERFACE (NcHIPertITwoFluids, nc_hipert_itwo_fluids, 0);
G_DEFINE_BOXED_TYPE (NcHIPertITwoFluidsEOM, nc_hipert_itwo_fluids_eom, nc_hipert_itwo_fluids_eom_dup, nc_hipert_itwo_fluids_eom_free);

static void
nc_hipert_itwo_fluids_default_init (NcHIPertITwoFluidsInterface *iface)
{
  iface->nuA2          = NULL;
  iface->dmzetanuA_nuA = NULL;
  iface->nuB2          = NULL;
  iface->dmSnuB_nuB    = NULL;
  iface->eom           = NULL;
}

/**
 * nc_hipert_itwo_fluids_eom_dup:
 * @tf_eom: a #NcHIPertITwoFluidsEOM.
 *
 * Duplicates @tf_eom.
 * 
 * Returns: (transfer full): a copy of @tf_eom.
 */
NcHIPertITwoFluidsEOM *
nc_hipert_itwo_fluids_eom_dup (NcHIPertITwoFluidsEOM *tf_eom)
{
  NcHIPertITwoFluidsEOM *tf_eom_dup = g_new (NcHIPertITwoFluidsEOM, 1);
  *tf_eom_dup = *tf_eom;
  return tf_eom_dup;
}

/**
 * nc_hipert_itwo_fluids_eom_free:
 * @tf_eom: a #NcHIPertITwoFluidsEOM.
 *
 * Frees @tf_eom.
 * 
 */
void
nc_hipert_itwo_fluids_eom_free (NcHIPertITwoFluidsEOM *tf_eom)
{
  g_free (tf_eom);
}
