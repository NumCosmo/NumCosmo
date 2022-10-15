/***************************************************************************
 *            nc_halo_bias_tinker.h
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

#ifndef _NC_HALO_BIAS_TINKER_H_
#define _NC_HALO_BIAS_TINKER_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_halo_bias.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_BIAS_TINKER             (nc_halo_bias_tinker_get_type ())
#define NC_HALO_BIAS_TINKER(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HALO_BIAS_TINKER, NcHaloBiasTinker))
#define NC_HALO_BIAS_TINKER_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HALO_BIAS_TINKER, NcHaloBiasTinkerClass))
#define NC_IS_HALO_BIAS_TINKER(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HALO_BIAS_TINKER))
#define NC_IS_HALO_BIAS_TINKER_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HALO_BIAS_TINKER))
#define NC_HALO_BIAS_TINKER_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HALO_BIAS_TINKER, NcHaloBiasTinkerClass))

typedef struct _NcHaloBiasTinkerClass NcHaloBiasTinkerClass;
typedef struct _NcHaloBiasTinker NcHaloBiasTinker;

struct _NcHaloBiasTinkerClass
{
  /*< private >*/
  NcHaloBiasClass parent_class;
};

struct _NcHaloBiasTinker
{
  /*< private >*/
  NcHaloBias parent_instance;
  gdouble delta_c;
  gdouble B;
  gdouble b;
  gdouble c;
};

GType nc_halo_bias_tinker_get_type (void) G_GNUC_CONST;

NcHaloBiasTinker *nc_halo_bias_tinker_new (NcHaloMassFunction *mfp);
NcHaloBiasTinker *nc_halo_bias_tinker_new_full (NcHaloMassFunction *mfp, gdouble delta_c, gdouble B, gdouble b, gdouble c);
void nc_halo_bias_tinker_set_delta_c (NcHaloBiasTinker *biasf_tinker, gdouble delta_c);
gdouble nc_halo_bias_tinker_get_delta_c (const NcHaloBiasTinker *biasf_tinker);
void nc_halo_bias_tinker_set_B (NcHaloBiasTinker *biasf_tinker, gdouble B);
gdouble nc_halo_bias_tinker_get_B (const NcHaloBiasTinker *biasf_tinker);
void nc_halo_bias_tinker_set_b (NcHaloBiasTinker *biasf_tinker, gdouble b);
gdouble nc_halo_bias_tinker_get_b (const NcHaloBiasTinker *biasf_tinker);
void nc_halo_bias_tinker_set_c (NcHaloBiasTinker *biasf_tinker, gdouble c);
gdouble nc_halo_bias_tinker_get_c (const NcHaloBiasTinker *biasf_tinker);

G_END_DECLS

#endif /* _NC_HALO_BIAS_TINKER_H_ */

