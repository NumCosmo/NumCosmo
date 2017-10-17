/***************************************************************************
 *            nc_hicosmo_qlinear.h
 *
 *  Mon Aug 11 20:00:24 2008
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

#ifndef _NC_HICOSMO_QLINEAR_H_
#define _NC_HICOSMO_QLINEAR_H_

#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_hicosmo.h>

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_QLINEAR             (nc_hicosmo_qlinear_get_type ())
#define NC_HICOSMO_QLINEAR(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_QLINEAR, NcHICosmoQLinear))
#define NC_HICOSMO_QLINEAR_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_QLINEAR, NcHICosmoQLinearClass))
#define NC_IS_HICOSMO_QLINEAR(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_QLINEAR))
#define NC_IS_HICOSMO_QLINEAR_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_QLINEAR))
#define NC_HICOSMO_QLINEAR_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_QLINEAR, NcHICosmoQLinearClass))

typedef struct _NcHICosmoQLinearClass NcHICosmoQLinearClass;
typedef struct _NcHICosmoQLinear NcHICosmoQLinear;

/**
 * NcHICosmoQLinearSParams:
 * @NC_HICOSMO_QLINEAR_H0: Hubble constant
 * @NC_HICOSMO_QLINEAR_OMEGA_T: total energy density of the universe
 * @NC_HICOSMO_QLINEAR_CD: comoving distance in units of the Hubble radius today
 * @NC_HICOSMO_QLINEAR_E: normalized Hubble function at Z1
 * @NC_HICOSMO_QLINEAR_Q: q-intercept term of the linear function
 * @NC_HICOSMO_QLINEAR_QP: slope of the linear function
 * @NC_HICOSMO_QLINEAR_Z1: initial redshift
 *
 */
typedef enum /*< enum,underscore_name=NC_HICOSMO_QLINEAR_SPARAMS >*/
{
  NC_HICOSMO_QLINEAR_H0 = 0,
  NC_HICOSMO_QLINEAR_OMEGA_T,
  NC_HICOSMO_QLINEAR_CD,
  NC_HICOSMO_QLINEAR_E,
  NC_HICOSMO_QLINEAR_Q,
  NC_HICOSMO_QLINEAR_QP,
  NC_HICOSMO_QLINEAR_Z1,         
  /* < private > */
  NC_HICOSMO_QLINEAR_SPARAM_LEN, /*< skip >*/
} NcHICosmoQLinearSParams;

#define NC_HICOSMO_QLINEAR_DEFAULT_H0      ncm_c_hubble_cte_wmap ()
#define NC_HICOSMO_QLINEAR_DEFAULT_OMEGA_T ( 1.0)
#define NC_HICOSMO_QLINEAR_DEFAULT_CD      ( 0.0)
#define NC_HICOSMO_QLINEAR_DEFAULT_E       ( 1.0)
#define NC_HICOSMO_QLINEAR_DEFAULT_Q       (-0.5)
#define NC_HICOSMO_QLINEAR_DEFAULT_QP      ( 1.0)
#define NC_HICOSMO_QLINEAR_DEFAULT_Z1      ( 0.0)

struct _NcHICosmoQLinearClass
{
  /*< private >*/
  NcHICosmoClass parent_class;
};

struct _NcHICosmoQLinear
{
  /*< private >*/
  NcHICosmo parent_instance;
};

GType nc_hicosmo_qlinear_get_type (void) G_GNUC_CONST;

NcHICosmoQLinear *nc_hicosmo_qlinear_new (void);
gdouble nc_hicosmo_qlinear_dE (gdouble z2, gdouble z1, gdouble q, gdouble qp);

G_END_DECLS

#endif /* _NC_HICOSMO_QLINEAR_H_ */
