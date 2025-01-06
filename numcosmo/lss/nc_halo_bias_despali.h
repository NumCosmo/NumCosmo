/***************************************************************************
 *            nc_halo_bias_despali.h
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
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

#ifndef _NC_HALO_BIAS_DESPALI_H_
#define _NC_HALO_BIAS_DESPALI_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_halo_bias.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_BIAS_DESPALI             (nc_halo_bias_despali_get_type ())
#define NC_HALO_BIAS_DESPALI(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HALO_BIAS_DESPALI, NcHaloBiasDespali))
#define NC_HALO_BIAS_DESPALI_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HALO_BIAS_DESPALI, NcHaloBiasDespaliClass))
#define NC_IS_HALO_BIAS_DESPALI(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HALO_BIAS_DESPALI))
#define NC_IS_HALO_BIAS_DESPALI_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HALO_BIAS_DESPALI))
#define NC_HALO_BIAS_DESPALI_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HALO_BIAS_DESPALI, NcHaloBiasDespaliClass))

typedef struct _NcHaloBiasDespaliClass NcHaloBiasDespaliClass;
typedef struct _NcHaloBiasDespali NcHaloBiasDespali;
struct _NcHaloBiasDespaliClass
{
  /*< private >*/
  NcHaloBiasClass parent_class;
};

struct _NcHaloBiasDespali
{
  /*< private >*/
  NcHaloBias parent_instance;
  gboolean eo;
  gboolean cmf;
};

GType nc_halo_bias_despali_get_type (void) G_GNUC_CONST;

NcHaloBiasDespali *nc_halo_bias_despali_new (NcHaloMassFunction *mfp);
NcHaloBiasDespali *nc_halo_bias_despali_new_full (NcHaloMassFunction *mfp, gboolean EO, gboolean CMF);
NcHaloBiasDespali *nc_halo_bias_despali_ref (NcHaloBiasDespali *biasf_despali);

void nc_halo_bias_despali_free (NcHaloBiasDespali *biasf_despali);
void nc_halo_bias_despali_clear (NcHaloBiasDespali **biasf_despali);

gdouble nc_halo_bias_despali_delta_c (NcHaloBiasDespali *biasf_despali , NcHICosmo *cosmo ,gdouble z);

gdouble nc_halo_bias_despali_delta_vir (NcHaloBiasDespali *biasf_despali , NcHICosmo *cosmo ,gdouble z);
void nc_halo_bias_despali_set_eo (NcHaloBiasDespali *biasf_despali, gboolean on);
gboolean nc_halo_bias_despali_get_eo (NcHaloBiasDespali *biasf_despali);
void nc_halo_bias_despali_set_cmf (NcHaloBiasDespali *biasf_despali, gboolean on);
gboolean nc_halo_bias_despali_get_cmf (NcHaloBiasDespali *biasf_despali);

G_END_DECLS

#endif /* _NC_HALO_BIAS_DESPALI_H_ */
