/***************************************************************************
 *            nc_hipert_adiabatic.h
 *
 *  Tue June 03 17:20:57 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_adiabatic.h
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

#ifndef _NC_HIPERT_ADIABATIC_H_
#define _NC_HIPERT_ADIABATIC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_ode_spline.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/nc_hipert.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_ADIABATIC             (nc_hipert_adiabatic_get_type ())
#define NC_HIPERT_ADIABATIC(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT_ADIABATIC, NcHIPertAdiabatic))
#define NC_HIPERT_ADIABATIC_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPERT_ADIABATIC, NcHIPertAdiabaticClass))
#define NC_IS_HIPERT_ADIABATIC(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT_ADIABATIC))
#define NC_IS_HIPERT_ADIABATIC_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPERT_ADIABATIC))
#define NC_HIPERT_ADIABATIC_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPERT_ADIABATIC, NcHIPertAdiabaticClass))

typedef struct _NcHIPertAdiabaticClass NcHIPertAdiabaticClass;
typedef struct _NcHIPertAdiabatic NcHIPertAdiabatic;

struct _NcHIPertAdiabaticClass
{
  NcHIPertClass parent_class;
};

typedef enum _NcHIPertAdiabaticVars
{
  NC_HIPERT_ADIABATIC_RE_ZETA = 0,
  NC_HIPERT_ADIABATIC_IM_ZETA,
  NC_HIPERT_ADIABATIC_RE_PZETA,
  NC_HIPERT_ADIABATIC_IM_PZETA,
} NcHIPertAdiabaticVars;

struct _NcHIPertAdiabatic
{
  NcHIPert parent_instance;
  /*< private >*/
  NcmOdeSpline *wkb_phase;
};

GType nc_hipert_adiabatic_get_type (void) G_GNUC_CONST;

NcHIPertAdiabatic *nc_hipert_adiabatic_new (void);
NcHIPertAdiabatic *nc_hipert_adiabatic_ref (NcHIPertAdiabatic *pa);
void nc_hipert_adiabatic_free (NcHIPertAdiabatic *pa);
void nc_hipert_adiabatic_clear (NcHIPertAdiabatic **pa);

void nc_hipert_adiabatic_prepare_wkb (NcHIPertAdiabatic *pa, NcHICosmo *cosmo, gdouble alpha_i, gdouble alpha_f);

gdouble nc_hipert_adiabatic_mass_zeta (NcHIPertAdiabatic *pa, NcHICosmo *cosmo, gdouble alpha);
gdouble nc_hipert_adiabatic_freq2_zeta (NcHIPertAdiabatic *pa, NcHICosmo *cosmo, gdouble k, gdouble alpha);

void nc_hipert_adiabatic_zeta_wkb (NcHIPertAdiabatic *pa, NcHICosmo *cosmo, gdouble k, gdouble alpha, gdouble *Re_zeta, gdouble *Im_zeta);
void nc_hipert_adiabatic_zeta_Pzeta_wkb (NcHIPertAdiabatic *pa, NcHICosmo *cosmo, gdouble k, gdouble alpha, gdouble *Re_zeta, gdouble *Im_zeta, gdouble *Re_Pzeta, gdouble *Im_Pzeta);
void nc_hipert_adiabatic_v_wkb (NcHIPertAdiabatic *pa, NcHICosmo *cosmo, gdouble k, gdouble alpha, gdouble *Re_v, gdouble *Im_v);

void nc_hipert_adiabatic_set_init_cond (NcHIPertAdiabatic *pa, NcHICosmo *cosmo, gdouble k, gdouble alphai, gdouble Re_zeta, gdouble Im_zeta, gdouble Re_Pzeta, gdouble Im_Pzeta);
void nc_hipert_adiabatic_set_init_cond_wkb (NcHIPertAdiabatic *pa, NcHICosmo *cosmo, gdouble k, gdouble alphai);
void nc_hipert_adiabatic_evolve (NcHIPertAdiabatic *pa, NcHICosmo *cosmo, gdouble alphaf);
void nc_hipert_adiabatic_get_values (NcHIPertAdiabatic *pa, gdouble *alphai, gdouble *Re_zeta, gdouble *Im_zeta, gdouble *Re_Pzeta, gdouble *Im_Pzeta);

G_END_DECLS

#endif /* _NC_HIPERT_ADIABATIC_H_ */
