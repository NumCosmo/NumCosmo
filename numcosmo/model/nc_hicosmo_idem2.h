/***************************************************************************
 *            nc_hicosmo_idem2.h
 *
 *  Wed March 08 13:49:14 2017
 *  Copyright  2017  Rodrigo vom Marttens
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Rodrigo vom Marttens 2017 <rodrigovonmarttens@gmail.com>
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

#ifndef _NC_HICOSMO_IDEM2_H_
#define _NC_HICOSMO_IDEM2_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/math/ncm_likelihood.h>
#include <numcosmo/math/ncm_spline2d.h>
#include <numcosmo/math/ncm_integral1d_ptr.h>

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_IDEM2             (nc_hicosmo_idem2_get_type ())
#define NC_HICOSMO_IDEM2(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_IDEM2, NcHICosmoIDEM2))
#define NC_HICOSMO_IDEM2_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_IDEM2, NcHICosmoIDEM2Class))
#define NC_IS_HICOSMO_IDEM2(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_IDEM2))
#define NC_IS_HICOSMO_IDEM2_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_IDEM2))
#define NC_HICOSMO_IDEM2_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_IDEM2, NcHICosmoIDEM2Class))

typedef struct _NcHICosmoIDEM2Class NcHICosmoIDEM2Class;
typedef struct _NcHICosmoIDEM2 NcHICosmoIDEM2;
typedef struct _NcHICosmoIDEM2Private NcHICosmoIDEM2Private;

typedef gdouble (*NcHICosmoIDEM2Func1) (NcHICosmoIDEM2 *cosmo_idem2, gdouble z);

/**
 * NcHICosmoIDEM2SParams:
 * @NC_HICOSMO_IDEM2_H0: FIXME
 * @NC_HICOSMO_IDEM2_OMEGA_C: FIXME
 * @NC_HICOSMO_IDEM2_OMEGA_X: FIXME
 * @NC_HICOSMO_IDEM2_T_GAMMA0: FIXME
 * @NC_HICOSMO_IDEM2_HE_YP: FIXME
 * @NC_HICOSMO_IDEM2_ENNU: FIXME
 * @NC_HICOSMO_IDEM2_OMEGA_B: FIXME
 * @NC_HICOSMO_IDEM2_GAMMA: FIXME
 *
 * FIXME
 *
 */
typedef enum /*< enum,underscore_name=NC_HICOSMO_DE_IDEM2_SPARAMS >*/
{
  NC_HICOSMO_IDEM2_H0 = 0,
  NC_HICOSMO_IDEM2_OMEGA_C,
  NC_HICOSMO_IDEM2_OMEGA_X,
  NC_HICOSMO_IDEM2_T_GAMMA0,
  NC_HICOSMO_IDEM2_HE_YP,
  NC_HICOSMO_IDEM2_ENNU,
  NC_HICOSMO_IDEM2_OMEGA_B,    
  NC_HICOSMO_IDEM2_GAMMA,    
  /* < private > */
  NC_HICOSMO_IDEM2_SPARAM_LEN, /*< skip >*/
} NcHICosmoIDEM2SParams;

/**
 * NcHICosmoIDEM2VParams:
 * @NC_HICOSMO_IDEM2_MASSNU_M: FIXME
 * @NC_HICOSMO_IDEM2_MASSNU_T: FIXME
 * @NC_HICOSMO_IDEM2_MASSNU_MU: FIXME
 * @NC_HICOSMO_IDEM2_MASSNU_G: FIXME
 * 
 * FIXME
 *
 */
typedef enum /*< enum,underscore_name=NC_HICOSMO_DE_IDEM2_VPARAMS >*/
{
  NC_HICOSMO_IDEM2_MASSNU_M = 0,
  NC_HICOSMO_IDEM2_MASSNU_T,
  NC_HICOSMO_IDEM2_MASSNU_MU,
  NC_HICOSMO_IDEM2_MASSNU_G,
  /* < private > */
  NC_HICOSMO_IDEM2_VPARAM_LEN, /*< skip >*/
} NcHICosmoIDEM2VParams;

#define NC_HICOSMO_IDEM2_DEFAULT_H0        ncm_c_hubble_cte_wmap ()
#define NC_HICOSMO_IDEM2_DEFAULT_OMEGA_C   (0.2568)
#define NC_HICOSMO_IDEM2_DEFAULT_OMEGA_X   (0.70)
#define NC_HICOSMO_IDEM2_DEFAULT_OMEGA_B   (0.0432)
#define NC_HICOSMO_IDEM2_DEFAULT_T_GAMMA0  (2.7245)
#define NC_HICOSMO_IDEM2_DEFAULT_HE_YP     (0.24)
#define NC_HICOSMO_IDEM2_DEFAULT_ENNU      (3.046)
#define NC_HICOSMO_IDEM2_DEFAULT_GAMMA     (0.0)
#define NC_HICOSMO_IDEM2_DEFAULT_NU_MASS   (1.0e-5)
#define NC_HICOSMO_IDEM2_DEFAULT_NU_T      (0.71611)
#define NC_HICOSMO_IDEM2_DEFAULT_NU_MU     (0.0)
#define NC_HICOSMO_IDEM2_DEFAULT_NU_G      (1.0)

struct _NcHICosmoIDEM2Class
{
  /*< private >*/
  NcHICosmoClass parent_class;
};

struct _NcHICosmoIDEM2
{
  /*< private >*/
  NcHICosmo parent_instance;
  NcHICosmoIDEM2Private *priv;
};

GType nc_hicosmo_idem2_get_type (void) G_GNUC_CONST;

void nc_hicosmo_idem2_omega_x2omega_k (NcHICosmoIDEM2 *cosmo_idem2);
void nc_hicosmo_idem2_cmb_params (NcHICosmoIDEM2 *cosmo_idem2);

/***********************************************************************/
/* Reparam CMB                                                         */
/***********************************************************************/

#define NC_TYPE_HICOSMO_IDEM2_REPARAM_CMB             (nc_hicosmo_idem2_reparam_cmb_get_type ())
#define NC_HICOSMO_IDEM2_REPARAM_CMB(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_IDEM2_REPARAM_CMB, NcHICosmoIDEM2ReparamCMB))
#define NC_HICOSMO_IDEM2_REPARAM_CMB_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_IDEM2_REPARAM_CMB, NcHICosmoIDEM2ReparamCMBClass))
#define NC_IS_HICOSMO_IDEM2_REPARAM_CMB(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_IDEM2_REPARAM_CMB))
#define NC_IS_HICOSMO_IDEM2_REPARAM_CMB_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_IDEM2_REPARAM_CMB))
#define NC_HICOSMO_IDEM2_REPARAM_CMB_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_IDEM2_REPARAM_CMB, NcHICosmoIDEM2ReparamCMBClass))

typedef struct _NcHICosmoIDEM2ReparamCMBClass NcHICosmoIDEM2ReparamCMBClass;
typedef struct _NcHICosmoIDEM2ReparamCMB NcHICosmoIDEM2ReparamCMB;

struct _NcHICosmoIDEM2ReparamCMBClass
{
  /*< private >*/
  NcmReparamClass parent_class;
};

struct _NcHICosmoIDEM2ReparamCMB
{
  /*< private >*/
  NcmReparam parent_instance;
};

GType nc_hicosmo_idem2_reparam_cmb_get_type (void) G_GNUC_CONST;

NcHICosmoIDEM2ReparamCMB *nc_hicosmo_idem2_reparam_cmb_new (guint length);

/***********************************************************************/
/* Reparam Omega_x -> Omega_k                                          */
/***********************************************************************/

#define NC_TYPE_HICOSMO_IDEM2_REPARAM_OK             (nc_hicosmo_idem2_reparam_ok_get_type ())
#define NC_HICOSMO_IDEM2_REPARAM_OK(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_IDEM2_REPARAM_OK, NcHICosmoIDEM2ReparamOk))
#define NC_HICOSMO_IDEM2_REPARAM_OK_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_IDEM2_REPARAM_OK, NcHICosmoIDEM2ReparamOkClass))
#define NC_IS_HICOSMO_IDEM2_REPARAM_OK(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_IDEM2_REPARAM_OK))
#define NC_IS_HICOSMO_IDEM2_REPARAM_OK_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_IDEM2_REPARAM_OK))
#define NC_HICOSMO_IDEM2_REPARAM_OK_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_IDEM2_REPARAM_OK, NcHICosmoIDEM2ReparamOkClass))

typedef struct _NcHICosmoIDEM2ReparamOkClass NcHICosmoIDEM2ReparamOkClass;
typedef struct _NcHICosmoIDEM2ReparamOk NcHICosmoIDEM2ReparamOk;

struct _NcHICosmoIDEM2ReparamOkClass
{
  /*< private >*/
  NcmReparamClass parent_class;
};

struct _NcHICosmoIDEM2ReparamOk
{
  /*< private >*/
  NcmReparam parent_instance;
};

GType nc_hicosmo_idem2_reparam_ok_get_type (void) G_GNUC_CONST;

NcHICosmoIDEM2ReparamOk *nc_hicosmo_idem2_reparam_ok_new (guint length);

G_END_DECLS

#endif /* _NC_HICOSMO_IDEM2_H_ */
