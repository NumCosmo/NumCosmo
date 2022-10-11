/***************************************************************************
 *            nc_halo_bias_st_spher.h
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

#ifndef _NC_HALO_BIAS_ST_SPHER_H_
#define _NC_HALO_BIAS_ST_SPHER_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_halo_bias.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_BIAS_ST_SPHER             (nc_halo_bias_st_spher_get_type ())
#define NC_HALO_BIAS_ST_SPHER(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HALO_BIAS_ST_SPHER, NcHaloBiasSTSpher))
#define NC_HALO_BIAS_ST_SPHER_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HALO_BIAS_ST_SPHER, NcHaloBiasSTSpherClass))
#define NC_IS_HALO_BIAS_ST_SPHER(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HALO_BIAS_ST_SPHER))
#define NC_IS_HALO_BIAS_ST_SPHER_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HALO_BIAS_ST_SPHER))
#define NC_HALO_BIAS_ST_SPHER_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HALO_BIAS_ST_SPHER, NcHaloBiasSTSpherClass))

typedef struct _NcHaloBiasSTSpherClass NcHaloBiasSTSpherClass;
typedef struct _NcHaloBiasSTSpher NcHaloBiasSTSpher;

struct _NcHaloBiasSTSpherClass
{
  /*< private >*/
  NcHaloBiasClass parent_class;
};

struct _NcHaloBiasSTSpher
{
  /*< private >*/
  NcHaloBias parent_instance;
  gdouble delta_c;
  gdouble a;
  gdouble p;
};

GType nc_halo_bias_st_spher_get_type (void) G_GNUC_CONST;

NcHaloBias *nc_halo_bias_st_spher_new (gdouble delta_c, gdouble a, gdouble p);
void nc_halo_bias_st_spher_set_delta_c (NcHaloBiasSTSpher *biasf_st_spher, gdouble delta_c);
gdouble nc_halo_bias_st_spher_get_delta_c (const NcHaloBiasSTSpher *biasf_st_spher);
void nc_halo_bias_st_spher_set_a (NcHaloBiasSTSpher *biasf_st_spher, gdouble a);
gdouble nc_halo_bias_st_spher_get_a (const NcHaloBiasSTSpher *biasf_st_spher);
void nc_halo_bias_st_spher_set_p (NcHaloBiasSTSpher *biasf_st_spher, gdouble p);
gdouble nc_halo_bias_st_spher_get_p (const NcHaloBiasSTSpher *biasf_st_spher);

G_END_DECLS

#endif /* _NC_HALO_BIAS_ST_SPHER_H_ */
