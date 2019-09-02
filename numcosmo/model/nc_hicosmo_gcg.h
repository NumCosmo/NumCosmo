/***************************************************************************
 *            nc_hicosmo_gcg.h
 *
 *  Wed March 08 13:49:14 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2017 <sandro@isoftware.com.br>
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

#ifndef _NC_HICOSMO_GCG_H_
#define _NC_HICOSMO_GCG_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/math/ncm_likelihood.h>
#include <numcosmo/math/ncm_spline2d.h>
#include <numcosmo/math/ncm_integral1d_ptr.h>

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_GCG             (nc_hicosmo_gcg_get_type ())
#define NC_HICOSMO_GCG(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_GCG, NcHICosmoGCG))
#define NC_HICOSMO_GCG_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_GCG, NcHICosmoGCGClass))
#define NC_IS_HICOSMO_GCG(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_GCG))
#define NC_IS_HICOSMO_GCG_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_GCG))
#define NC_HICOSMO_GCG_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_GCG, NcHICosmoGCGClass))

typedef struct _NcHICosmoGCGClass NcHICosmoGCGClass;
typedef struct _NcHICosmoGCG NcHICosmoGCG;
typedef struct _NcHICosmoGCGPrivate NcHICosmoGCGPrivate;

typedef gdouble (*NcHICosmoGCGFunc1) (NcHICosmoGCG *cosmo_gcg, gdouble z);

/**
 * NcHICosmoGCGSParams:
 * @NC_HICOSMO_GCG_H0: FIXME
 * @NC_HICOSMO_GCG_OMEGA_C: FIXME
 * @NC_HICOSMO_GCG_OMEGA_X: FIXME
 * @NC_HICOSMO_GCG_T_GAMMA0: FIXME
 * @NC_HICOSMO_GCG_HE_YP: FIXME
 * @NC_HICOSMO_GCG_ENNU: FIXME
 * @NC_HICOSMO_GCG_OMEGA_B: FIXME
 * @NC_HICOSMO_GCG_GAMMA: FIXME
 *
 * FIXME
 *
 */
typedef enum /*< enum,underscore_name=NC_HICOSMO_DE_GCG_SPARAMS >*/
{
  NC_HICOSMO_GCG_H0 = 0,
  NC_HICOSMO_GCG_OMEGA_C,
  NC_HICOSMO_GCG_OMEGA_X,
  NC_HICOSMO_GCG_T_GAMMA0,
  NC_HICOSMO_GCG_HE_YP,
  NC_HICOSMO_GCG_ENNU,
  NC_HICOSMO_GCG_OMEGA_B,    
  NC_HICOSMO_GCG_GAMMA,   
  /* < private > */
  NC_HICOSMO_GCG_SPARAM_LEN, /*< skip >*/
} NcHICosmoGCGSParams;

/**
 * NcHICosmoGCGVParams:
 * @NC_HICOSMO_GCG_MASSNU_M: FIXME
 * @NC_HICOSMO_GCG_MASSNU_T: FIXME
 * @NC_HICOSMO_GCG_MASSNU_MU: FIXME
 * @NC_HICOSMO_GCG_MASSNU_G: FIXME
 * 
 * FIXME
 *
 */
typedef enum /*< enum,underscore_name=NC_HICOSMO_DE_GCG_VPARAMS >*/
{
  NC_HICOSMO_GCG_MASSNU_M = 0,
  NC_HICOSMO_GCG_MASSNU_T,
  NC_HICOSMO_GCG_MASSNU_MU,
  NC_HICOSMO_GCG_MASSNU_G,   
  /* < private > */
  NC_HICOSMO_GCG_VPARAM_LEN, /*< skip >*/
} NcHICosmoGCGVParams;

#define NC_HICOSMO_GCG_DEFAULT_H0        ncm_c_hubble_cte_wmap ()
#define NC_HICOSMO_GCG_DEFAULT_OMEGA_C   (0.2568)
#define NC_HICOSMO_GCG_DEFAULT_OMEGA_X   (0.70)
#define NC_HICOSMO_GCG_DEFAULT_OMEGA_B   (0.0432)
#define NC_HICOSMO_GCG_DEFAULT_T_GAMMA0  (2.7245)
#define NC_HICOSMO_GCG_DEFAULT_HE_YP     (0.24)
#define NC_HICOSMO_GCG_DEFAULT_ENNU      (3.046)
#define NC_HICOSMO_GCG_DEFAULT_GAMMA     (0.0)
#define NC_HICOSMO_GCG_DEFAULT_NU_MASS   (1.0e-5)
#define NC_HICOSMO_GCG_DEFAULT_NU_T      (0.71611)
#define NC_HICOSMO_GCG_DEFAULT_NU_MU     (0.0)
#define NC_HICOSMO_GCG_DEFAULT_NU_G      (1.0)

struct _NcHICosmoGCGClass
{
  /*< private >*/
  NcHICosmoClass parent_class;
};

struct _NcHICosmoGCG
{
  /*< private >*/
  NcHICosmo parent_instance;
  NcHICosmoGCGPrivate *priv;
};

GType nc_hicosmo_gcg_get_type (void) G_GNUC_CONST;

void nc_hicosmo_gcg_omega_x2omega_k (NcHICosmoGCG *cosmo_gcg);
void nc_hicosmo_gcg_cmb_params (NcHICosmoGCG *cosmo_gcg);

/***********************************************************************/
/* Reparam CMB                                                         */
/***********************************************************************/

#define NC_TYPE_HICOSMO_GCG_REPARAM_CMB             (nc_hicosmo_gcg_reparam_cmb_get_type ())
#define NC_HICOSMO_GCG_REPARAM_CMB(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_GCG_REPARAM_CMB, NcHICosmoGCGReparamCMB))
#define NC_HICOSMO_GCG_REPARAM_CMB_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_GCG_REPARAM_CMB, NcHICosmoGCGReparamCMBClass))
#define NC_IS_HICOSMO_GCG_REPARAM_CMB(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_GCG_REPARAM_CMB))
#define NC_IS_HICOSMO_GCG_REPARAM_CMB_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_GCG_REPARAM_CMB))
#define NC_HICOSMO_GCG_REPARAM_CMB_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_GCG_REPARAM_CMB, NcHICosmoGCGReparamCMBClass))

typedef struct _NcHICosmoGCGReparamCMBClass NcHICosmoGCGReparamCMBClass;
typedef struct _NcHICosmoGCGReparamCMB NcHICosmoGCGReparamCMB;

struct _NcHICosmoGCGReparamCMBClass
{
  /*< private >*/
  NcmReparamClass parent_class;
};

struct _NcHICosmoGCGReparamCMB
{
  /*< private >*/
  NcmReparam parent_instance;
};

GType nc_hicosmo_gcg_reparam_cmb_get_type (void) G_GNUC_CONST;

NcHICosmoGCGReparamCMB *nc_hicosmo_gcg_reparam_cmb_new (guint length);

/***********************************************************************/
/* Reparam Omega_x -> Omega_k                                          */
/***********************************************************************/

#define NC_TYPE_HICOSMO_GCG_REPARAM_OK             (nc_hicosmo_gcg_reparam_ok_get_type ())
#define NC_HICOSMO_GCG_REPARAM_OK(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_GCG_REPARAM_OK, NcHICosmoGCGReparamOk))
#define NC_HICOSMO_GCG_REPARAM_OK_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_GCG_REPARAM_OK, NcHICosmoGCGReparamOkClass))
#define NC_IS_HICOSMO_GCG_REPARAM_OK(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_GCG_REPARAM_OK))
#define NC_IS_HICOSMO_GCG_REPARAM_OK_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_GCG_REPARAM_OK))
#define NC_HICOSMO_GCG_REPARAM_OK_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_GCG_REPARAM_OK, NcHICosmoGCGReparamOkClass))

typedef struct _NcHICosmoGCGReparamOkClass NcHICosmoGCGReparamOkClass;
typedef struct _NcHICosmoGCGReparamOk NcHICosmoGCGReparamOk;

struct _NcHICosmoGCGReparamOkClass
{
  /*< private >*/
  NcmReparamClass parent_class;
};

struct _NcHICosmoGCGReparamOk
{
  /*< private >*/
  NcmReparam parent_instance;
};

GType nc_hicosmo_gcg_reparam_ok_get_type (void) G_GNUC_CONST;

NcHICosmoGCGReparamOk *nc_hicosmo_gcg_reparam_ok_new (guint length);

G_END_DECLS

#endif /* _NC_HICOSMO_GCG_H_ */
