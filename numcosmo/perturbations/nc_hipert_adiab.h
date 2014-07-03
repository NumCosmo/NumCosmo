/***************************************************************************
 *            nc_hipert_adiab.h
 *
 *  Tue June 03 17:20:57 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_adiab.h
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

#ifndef _NC_HIPERT_ADIAB_H_
#define _NC_HIPERT_ADIAB_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_ode_spline.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/nc_hipert.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_ADIAB             (nc_hipert_adiab_get_type ())
#define NC_HIPERT_ADIAB(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT_ADIAB, NcHIPertAdiab))
#define NC_HIPERT_ADIAB_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPERT_ADIAB, NcHIPertAdiabClass))
#define NC_IS_HIPERT_ADIAB(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT_ADIAB))
#define NC_IS_HIPERT_ADIAB_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPERT_ADIAB))
#define NC_HIPERT_ADIAB_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPERT_ADIAB, NcHIPertAdiabClass))

typedef struct _NcHIPertAdiabClass NcHIPertAdiabClass;
typedef struct _NcHIPertAdiab NcHIPertAdiab;

struct _NcHIPertAdiabClass
{
  /*< private >*/
  NcHIPertClass parent_class;
};

/**
 * NcHIPertAdiabVars:
 * @NC_HIPERT_ADIAB_RE_ZETA: $\text{Re}(\zeta)$
 * @NC_HIPERT_ADIAB_IM_ZETA: $\text{Im}(\zeta)$
 * @NC_HIPERT_ADIAB_RE_PZETA: $\text{Re}(P_\zeta)$
 * @NC_HIPERT_ADIAB_IM_PZETA: $\text{Im}(P_\zeta)$
 * 
 * Perturbation variables enumerator.
 * 
 */
typedef enum _NcHIPertAdiabVars
{
  NC_HIPERT_ADIAB_RE_ZETA = 0,
  NC_HIPERT_ADIAB_IM_ZETA,
  NC_HIPERT_ADIAB_RE_PZETA,
  NC_HIPERT_ADIAB_IM_PZETA,
} NcHIPertAdiabVars;

struct _NcHIPertAdiab
{
  /*< private >*/
  NcHIPert parent_instance;
  NcmOdeSpline *wkb_phase;
};

GType nc_hipert_adiab_get_type (void) G_GNUC_CONST;

NcHIPertAdiab *nc_hipert_adiab_new (void);
NcHIPertAdiab *nc_hipert_adiab_ref (NcHIPertAdiab *pa);
void nc_hipert_adiab_free (NcHIPertAdiab *pa);
void nc_hipert_adiab_clear (NcHIPertAdiab **pa);

void nc_hipert_adiab_prepare_wkb (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alpha_i, gdouble alpha_f);

void nc_hipert_adiab_wkb_zeta (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alpha, gdouble *Re_zeta, gdouble *Im_zeta);
void nc_hipert_adiab_wkb_zeta_Pzeta (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alpha, gdouble *Re_zeta, gdouble *Im_zeta, gdouble *Re_Pzeta, gdouble *Im_Pzeta);
void nc_hipert_adiab_wkb_v (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alpha, gdouble *Re_v, gdouble *Im_v);
gdouble nc_hipert_adiab_wkb_maxtime (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alpha0, gdouble alpha1);

void nc_hipert_adiab_set_init_cond (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alphai, gdouble Re_zeta, gdouble Im_zeta, gdouble Re_Pzeta, gdouble Im_Pzeta);
void nc_hipert_adiab_set_init_cond_wkb (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alphai);
void nc_hipert_adiab_evolve (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alphaf);
void nc_hipert_adiab_get_values (NcHIPertAdiab *pa, gdouble *alphai, gdouble *Re_zeta, gdouble *Im_zeta, gdouble *Re_Pzeta, gdouble *Im_Pzeta);

G_END_DECLS

#endif /* _NC_HIPERT_ADIAB_H_ */
