/***************************************************************************
 *            nc_hicosmo_qconst.h
 *
 *  Mon Aug 11 19:58:22 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef _NC_HICOSMO_QCONST_H_
#define _NC_HICOSMO_QCONST_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_hicosmo.h>

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_QCONST             (nc_hicosmo_qconst_get_type ())
#define NC_HICOSMO_QCONST(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_QCONST, NcHICosmoQConst))
#define NC_HICOSMO_QCONST_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_QCONST, NcHICosmoQConstClass))
#define NC_IS_HICOSMO_QCONST(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_QCONST))
#define NC_IS_HICOSMO_QCONST_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_QCONST))
#define NC_HICOSMO_QCONST_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_QCONST, NcHICosmoQConstClass))

typedef struct _NcHICosmoQConstClass NcHICosmoQConstClass;
typedef struct _NcHICosmoQConst NcHICosmoQConst;

/**
 * NcHICosmoQConstSParams:
 * @NC_HICOSMO_QCONST_H0: Hubble constant
 * @NC_HICOSMO_QCONST_OMEGA_T: total energy density of the universe
 * @NC_HICOSMO_QCONST_CD: comoving distance in units of the Hubble radius today
 * @NC_HICOSMO_QCONST_E: normalized Hubble function at Z1
 * @NC_HICOSMO_QCONST_Q: constant vlaue of the deceleration function
 * @NC_HICOSMO_QCONST_Z1: initial redshift
 *
 */
typedef enum /*< enum,underscore_name=NC_HICOSMO_QCONST_SPARAMS >*/
{
  NC_HICOSMO_QCONST_H0 = 0,
  NC_HICOSMO_QCONST_OMEGA_T,
  NC_HICOSMO_QCONST_CD,
  NC_HICOSMO_QCONST_E,
  NC_HICOSMO_QCONST_Q,
  NC_HICOSMO_QCONST_Z1,         
  /* < private > */
  NC_HICOSMO_QCONST_SPARAM_LEN, /*< skip >*/
} NcHICosmoQConstSParams;

#define NC_HICOSMO_QCONST_DEFAULT_H0      ncm_c_hubble_cte_wmap ()
#define NC_HICOSMO_QCONST_DEFAULT_OMEGA_T ( 1.0)
#define NC_HICOSMO_QCONST_DEFAULT_CD      ( 0.0)
#define NC_HICOSMO_QCONST_DEFAULT_E       ( 1.0)
#define NC_HICOSMO_QCONST_DEFAULT_Q       (-0.5)
#define NC_HICOSMO_QCONST_DEFAULT_Z1      ( 0.0)

struct _NcHICosmoQConstClass
{
  /*< private >*/
  NcHICosmoClass parent_class;
};

struct _NcHICosmoQConst
{
  /*< private >*/
  NcHICosmo parent_instance;
};

GType nc_hicosmo_qconst_get_type (void) G_GNUC_CONST;

NcHICosmoQConst *nc_hicosmo_qconst_new (void);

G_END_DECLS

#endif /* _NC_HICOSMO_QCONST_H_ */
