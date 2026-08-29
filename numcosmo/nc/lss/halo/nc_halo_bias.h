/***************************************************************************
 *            nc_halo_bias.h
 *
 *  Tue June 28 15:41:57 2011
 *  Copyright  2011  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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

#ifndef _NC_HALO_BIAS_H_
#define _NC_HALO_BIAS_H_

#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/halo/nc_halo_mass_function.h>
#include <numcosmo/nc/background/nc_hicosmo.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_BIAS             (nc_halo_bias_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcHaloBias, nc_halo_bias, NC, HALO_BIAS, GObject)

struct _NcHaloBiasClass
{
  /*< private >*/
  GObjectClass parent_class;

  gdouble (*eval) (NcHaloBias *bias, NcHICosmo *cosmo, gdouble sigma, gdouble lnM, gdouble z);

  /* Padding to allow adding up to 17 more virtual functions without breaking ABI. */
  gpointer padding[17];
};

gdouble nc_halo_bias_eval (NcHaloBias *bias, NcHICosmo *cosmo, gdouble sigma, gdouble lnM, gdouble z);
void nc_halo_bias_free (NcHaloBias *bias);
void nc_halo_bias_clear (NcHaloBias **bias);

NcHaloMassFunction *nc_halo_bias_peek_mass_function (NcHaloBias *bias);

gdouble nc_halo_bias_integrand (NcHaloBias *mbiasf, NcHICosmo *cosmo, gdouble lnM, gdouble z);

G_END_DECLS

#endif /* _NC_HALO_BIAS_H_ */

