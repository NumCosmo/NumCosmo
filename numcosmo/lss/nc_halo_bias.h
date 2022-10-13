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
#include <numcosmo/lss/nc_halo_mass_function.h>
#include <numcosmo/nc_hicosmo.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_BIAS             (nc_halo_bias_get_type ())
#define NC_HALO_BIAS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HALO_BIAS, NcHaloBias))
#define NC_HALO_BIAS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HALO_BIAS, NcHaloBiasClass))
#define NC_IS_HALO_BIAS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HALO_BIAS))
#define NC_IS_HALO_BIAS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HALO_BIAS))
#define NC_HALO_BIAS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HALO_BIAS, NcHaloBiasClass))

typedef struct _NcHaloBiasClass NcHaloBiasClass;
typedef struct _NcHaloBias NcHaloBias;

struct _NcHaloBiasClass
{
  /*< private >*/
  GObjectClass parent_class;
  gdouble (*eval) (NcHaloBias *bias, NcHICosmo *cosmo, gdouble sigma, gdouble z); 
};

struct _NcHaloBias
{
  /*< private >*/
  GObject parent_instance;
  NcHaloMassFunction *mfp;
};

GType nc_halo_bias_get_type (void) G_GNUC_CONST;

gdouble nc_halo_bias_eval (NcHaloBias *bias, NcHICosmo *cosmo, gdouble sigma, gdouble z);
void nc_halo_bias_free (NcHaloBias *bias);
void nc_halo_bias_clear (NcHaloBias **bias);

gdouble nc_halo_bias_integrand (NcHaloBias *mbiasf, NcHICosmo *cosmo, gdouble lnM, gdouble z);

G_END_DECLS

#endif /* _NC_HALO_BIAS_H_ */
