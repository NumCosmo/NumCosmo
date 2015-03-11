/***************************************************************************
 *            nc_hipert_iadiab.c
 *
 *  Fri July 18 15:17:11 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_iadiab.c
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
 * SECTION:nc_hipert_iadiab
 * @title: NcHIPertIAdiab
 * @short_description: Perturbation interface for adiabatic mode only.
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "perturbations/nc_hipert_iadiab.h"

G_DEFINE_INTERFACE (NcHIPertIAdiab, nc_hipert_iadiab, 0);
G_DEFINE_BOXED_TYPE (NcHIPertIAdiabEOM, nc_hipert_iadiab_eom, nc_hipert_iadiab_eom_dup, nc_hipert_iadiab_eom_free);

static void
nc_hipert_iadiab_default_init (NcHIPertIAdiabInterface *iface)
{
  iface->nuA2          = NULL;
  iface->dmzetanuA_nuA = NULL;
  iface->eom           = NULL;
}

/**
 * nc_hipert_iadiab_eom_dup:
 * @adiab_eom: a #NcHIPertIAdiabEOM.
 *
 * Duplicates @adiab_eom.
 * 
 * Returns: (transfer full): a copy of @adiab.
 */
NcHIPertIAdiabEOM *
nc_hipert_iadiab_eom_dup (NcHIPertIAdiabEOM *adiab_eom)
{
  NcHIPertIAdiabEOM *adiab_eom_dup = g_new (NcHIPertIAdiabEOM, 1);
  *adiab_eom_dup = *adiab_eom;
  return adiab_eom_dup;
}

/**
 * nc_hipert_iadiab_eom_free:
 * @adiab_eom: a #NcHIPertIAdiabEOM.
 *
 * Frees @adiab_eom.
 * 
 */
void
nc_hipert_iadiab_eom_free (NcHIPertIAdiabEOM *adiab_eom)
{
  g_free (adiab_eom);
}
