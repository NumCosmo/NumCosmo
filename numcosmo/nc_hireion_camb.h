/***************************************************************************
 *            nc_hireion_camb.h
 *
 *  Thu December 10 11:56:27 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hireion_camb.h
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_HIREION_CAMB_H_
#define _NC_HIREION_CAMB_H_

#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model_ctrl.h>
#include <numcosmo/nc_hireion.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NC_TYPE_HIREION_CAMB             (nc_hireion_camb_get_type ())
#define NC_HIREION_CAMB(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIREION_CAMB, NcHIReionCamb))
#define NC_HIREION_CAMB_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIREION_CAMB, NcHIReionCambClass))
#define NC_IS_HIREION_CAMB(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIREION_CAMB))
#define NC_IS_HIREION_CAMB_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIREION_CAMB))
#define NC_HIREION_CAMB_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIREION_CAMB, NcHIReionCambClass))

typedef struct _NcHIReionCambClass NcHIReionCambClass;
typedef struct _NcHIReionCamb NcHIReionCamb;

struct _NcHIReionCambClass
{
  /*< private >*/
  NcHIReionClass parent_class;
};

/**
 * NcHIReionCambSParams:
 * @NC_HIREION_CAMB_HII_HEII_Z: FIXME
 * @NC_HIREION_CAMB_HEIII_Z: FIXME
 *
 * FIXME
 * 
 */
typedef enum /*< enum,underscore_name=NC_HIREION_CAMB_SPARAMS >*/
{
  NC_HIREION_CAMB_HII_HEII_Z = 0,
  NC_HIREION_CAMB_HEIII_Z,    
  /* < private > */
  NC_HIREION_CAMB_SPARAM_LEN, /*< skip >*/
} NcHIReionCambSParams;


struct _NcHIReionCamb
{
  /*< private >*/
  NcHIReion parent_instance;
  gdouble HII_HeII_reion_delta;
  gdouble HeIII_reion_delta;
  gdouble HII_HeII_reion_expo;
  gdouble HII_HeII_reion_delta_eff;
  gdouble HII_HeII_reion_x_pow_expo;
  gboolean HEII_reionized;
  gsl_root_fsolver *fsol;
  NcmModelCtrl *tau_ctrl;
};

GType nc_hireion_camb_get_type (void) G_GNUC_CONST;

NcHIReionCamb *nc_hireion_camb_new (void);
gdouble nc_hireion_camb_calc_z_from_tau (NcHIReionCamb *reion_camb, NcHICosmo *cosmo, const gdouble tau);
void nc_hireion_camb_set_z_from_tau (NcHIReionCamb *reion_camb, NcHICosmo *cosmo, const gdouble tau);
void nc_hireion_camb_z_to_tau (NcHIReionCamb *reion_camb, NcHICosmo *cosmo);

#define NC_HIREION_CAMB_DEFAULT_HII_HEII_REION_DELTA (0.5)
#define NC_HIREION_CAMB_DEFAULT_HEIII_REION_DELTA    (0.5)
#define NC_HIREION_CAMB_DEFAULT_HII_HEII_REION_EXPO  (1.5)

#define NC_HIREION_CAMB_DEFAULT_HII_HEII_Z (13.0)
#define NC_HIREION_CAMB_DEFAULT_HEIII_Z    (3.5)

G_END_DECLS

#endif /* _NC_HIREION_CAMB_H_ */
