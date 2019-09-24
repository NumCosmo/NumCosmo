/***************************************************************************
 *            nc_hipert_itwo_fluids.h
 *
 *  Tue July 22 17:37:11 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_itwo_fluids.h
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

#ifndef _NC_HIPERT_ITWO_FLUIDS_H_
#define _NC_HIPERT_ITWO_FLUIDS_H_

#include <glib-object.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/perturbations/nc_hipert_wkb.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_ITWO_FLUIDS               (nc_hipert_itwo_fluids_get_type ())
#define NC_HIPERT_ITWO_FLUIDS(obj)               (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT_ITWO_FLUIDS, NcHIPertITwoFluids))
#define NC_IS_HIPERT_ITWO_FLUIDS(obj)            (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT_ITWO_FLUIDS))
#define NC_HIPERT_ITWO_FLUIDS_GET_INTERFACE(obj) (G_TYPE_INSTANCE_GET_INTERFACE ((obj), NC_TYPE_HIPERT_ITWO_FLUIDS, NcHIPertITwoFluidsInterface))

typedef struct _NcHIPertITwoFluids NcHIPertITwoFluids;
typedef struct _NcHIPertITwoFluidsInterface NcHIPertITwoFluidsInterface;
typedef struct _NcHIPertITwoFluidsEOM NcHIPertITwoFluidsEOM;
typedef struct _NcHIPertITwoFluidsTV NcHIPertITwoFluidsTV;

typedef NcHIPertITwoFluidsEOM *(*NcHIPertITwoFluidsFuncEOM) (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);
typedef NcHIPertITwoFluidsTV *(*NcHIPertITwoFluidsFuncTV) (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);

struct _NcHIPertITwoFluidsInterface
{
  /*< private >*/
  GTypeInterface parent;
  NcHIPertITwoFluidsFuncEOM eom;
  NcHIPertITwoFluidsFuncTV tv;
};

/**
 * NcHIPertITwoFluidsEOM:
 * 
 * FIXME
 * 
 */
struct _NcHIPertITwoFluidsEOM
{
  /*< private >*/
  guint64 skey;
  gdouble alpha;
  gdouble k;
  gdouble nu1;
  gdouble nu2;
  gdouble gammabar11;
  gdouble gammabar22;
  gdouble gammabar12;
  gdouble taubar;
  gdouble m_zeta;
  gdouble m_s;
  gdouble mnu2_zeta;
  gdouble mnu2_s;
  gdouble y;
};

/**
 * NcHIPertITwoFluidsVars:
 * @NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R: FIXME
 * @NC_HIPERT_ITWO_FLUIDS_VARS_S_R:  FIXME
 * @NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R: FIXME
 * @NC_HIPERT_ITWO_FLUIDS_VARS_PS_R:  FIXME
 * @NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I: FIXME
 * @NC_HIPERT_ITWO_FLUIDS_VARS_S_I:  FIXME
 * @NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I: FIXME
 * @NC_HIPERT_ITWO_FLUIDS_VARS_PS_I:  FIXME
 * 
 * FIXME
 * 
 */
typedef enum /*< enum,underscore_name=NC_HIPERT_ITWO_FLUIDS_VARS >*/
{
  NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R = 0,
  NC_HIPERT_ITWO_FLUIDS_VARS_S_R,
  NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R,
  NC_HIPERT_ITWO_FLUIDS_VARS_PS_R,
  NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I,
  NC_HIPERT_ITWO_FLUIDS_VARS_S_I,
  NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I,
  NC_HIPERT_ITWO_FLUIDS_VARS_PS_I,  
  /* < private > */
  NC_HIPERT_ITWO_FLUIDS_VARS_LEN,   /*< skip >*/
} NcHIPertITwoFluidsVars;

#define NC_HIPERT_ITWO_FLUIDS_VARS_Q_R1 NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_R
#define NC_HIPERT_ITWO_FLUIDS_VARS_Q_R2 NC_HIPERT_ITWO_FLUIDS_VARS_S_R
#define NC_HIPERT_ITWO_FLUIDS_VARS_P_R1 NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_R
#define NC_HIPERT_ITWO_FLUIDS_VARS_P_R2 NC_HIPERT_ITWO_FLUIDS_VARS_PS_R
#define NC_HIPERT_ITWO_FLUIDS_VARS_Q_I1 NC_HIPERT_ITWO_FLUIDS_VARS_ZETA_I
#define NC_HIPERT_ITWO_FLUIDS_VARS_Q_I2 NC_HIPERT_ITWO_FLUIDS_VARS_S_I
#define NC_HIPERT_ITWO_FLUIDS_VARS_P_I1 NC_HIPERT_ITWO_FLUIDS_VARS_PZETA_I
#define NC_HIPERT_ITWO_FLUIDS_VARS_P_I2 NC_HIPERT_ITWO_FLUIDS_VARS_PS_I

/**
 * NcHIPertITwoFluidsTV:
 * 
 * FIXME
 * 
 */
struct _NcHIPertITwoFluidsTV
{
  /*< private >*/
  guint64 skey;
  gdouble alpha;
  gdouble k;
  gdouble zeta[NC_HIPERT_ITWO_FLUIDS_VARS_LEN / 2];
  gdouble s[NC_HIPERT_ITWO_FLUIDS_VARS_LEN / 2];
  gdouble Pzeta[NC_HIPERT_ITWO_FLUIDS_VARS_LEN / 2];
  gdouble Ps[NC_HIPERT_ITWO_FLUIDS_VARS_LEN / 2];
};

GType nc_hipert_itwo_fluids_eom_get_type (void) G_GNUC_CONST;
GType nc_hipert_itwo_fluids_tv_get_type (void) G_GNUC_CONST;
GType nc_hipert_itwo_fluids_get_type (void) G_GNUC_CONST;

NcHIPertITwoFluidsEOM *nc_hipert_itwo_fluids_eom_dup (NcHIPertITwoFluidsEOM *tf_eom);
void nc_hipert_itwo_fluids_eom_free (NcHIPertITwoFluidsEOM *tf_eom);

NcHIPertITwoFluidsTV *nc_hipert_itwo_fluids_tv_dup (NcHIPertITwoFluidsTV *tf_tv);
void nc_hipert_itwo_fluids_tv_free (NcHIPertITwoFluidsTV *tf_tv);

NCM_INLINE NcHIPertITwoFluidsEOM *nc_hipert_itwo_fluids_eom_eval (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);
NCM_INLINE NcHIPertITwoFluidsTV *nc_hipert_itwo_fluids_tv_eval (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);

G_END_DECLS

#endif /* _NC_HIPERT_ITWO_FLUIDS_H_ */

#ifndef _NC_HIPERT_ITWO_FLUIDS_INLINE_H_
#define _NC_HIPERT_ITWO_FLUIDS_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE NcHIPertITwoFluidsEOM *
nc_hipert_itwo_fluids_eom_eval (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k)
{
  return NC_HIPERT_ITWO_FLUIDS_GET_INTERFACE (itf)->eom (itf, alpha, k);
}

NCM_INLINE NcHIPertITwoFluidsTV *
nc_hipert_itwo_fluids_tv_eval (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k)
{
  return NC_HIPERT_ITWO_FLUIDS_GET_INTERFACE (itf)->tv (itf, alpha, k);
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_HIPERT_ITWO_FLUIDS_INLINE_H_ */
