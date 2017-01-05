/***************************************************************************
 *            nc_hicosmo_de.h
 *
 *  Mon Aug 11 19:54:00 2008
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

#ifndef _NC_HICOSMO_DE_H_
#define _NC_HICOSMO_DE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/math/ncm_likelihood.h>
#include <numcosmo/math/ncm_spline2d.h>
#include <numcosmo/math/ncm_integral1d_ptr.h>

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_DE             (nc_hicosmo_de_get_type ())
#define NC_HICOSMO_DE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_DE, NcHICosmoDE))
#define NC_HICOSMO_DE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_DE, NcHICosmoDEClass))
#define NC_IS_HICOSMO_DE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_DE))
#define NC_IS_HICOSMO_DE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_DE))
#define NC_HICOSMO_DE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_DE, NcHICosmoDEClass))

typedef struct _NcHICosmoDEClass NcHICosmoDEClass;
typedef struct _NcHICosmoDE NcHICosmoDE;
typedef struct _NcHICosmoDEPrivate NcHICosmoDEPrivate;

/**
 * NcHICosmoDEImpl:
 * @NC_HICOSMO_DE_IMPL_E2Omega_de: FIXME
 * @NC_HICOSMO_DE_IMPL_dE2Omega_de_dz: FIXME
 * @NC_HICOSMO_DE_IMPL_d2E2Omega_de_dz2: FIXME
 * @NC_HICOSMO_DE_IMPL_w_de: FIXME
 *
 * FIXME
 *
 */
typedef enum _NcHICosmoDEImpl
{
  NC_HICOSMO_DE_IMPL_E2Omega_de = NC_HICOSMO_IMPL_LAST,
  NC_HICOSMO_DE_IMPL_dE2Omega_de_dz,
  NC_HICOSMO_DE_IMPL_d2E2Omega_de_dz2,
  NC_HICOSMO_DE_IMPL_w_de, /*< private >*/
  NC_HICOSMO_DE_IMPL_LAST, /*< skip >*/
} NcHICosmoDEImpl;

typedef gdouble (*NcHICosmoDEFunc1) (NcHICosmoDE *cosmo_de, gdouble z);

/**
 * NcHICosmoDEParams:
 * @NC_HICOSMO_DE_H0: FIXME
 * @NC_HICOSMO_DE_OMEGA_C: FIXME
 * @NC_HICOSMO_DE_OMEGA_X: FIXME
 * @NC_HICOSMO_DE_T_GAMMA0: FIXME
 * @NC_HICOSMO_DE_HE_YP: FIXME
 * @NC_HICOSMO_DE_ENNU: FIXME
 * @NC_HICOSMO_DE_OMEGA_B: FIXME
 *
 * FIXME
 *
 */
typedef enum _NcHICosmoDEParams
{
  NC_HICOSMO_DE_H0 = 0,
  NC_HICOSMO_DE_OMEGA_C,
  NC_HICOSMO_DE_OMEGA_X,
  NC_HICOSMO_DE_T_GAMMA0,
  NC_HICOSMO_DE_HE_YP,
  NC_HICOSMO_DE_ENNU,
  NC_HICOSMO_DE_OMEGA_B,    /*< private >*/
  NC_HICOSMO_DE_SPARAM_LEN, /*< skip >*/
} NcHICosmoDEParams;

/**
 * NcHICosmoDEVParams:
 * @NC_HICOSMO_DE_MASSNU_M: FIXME
 * @NC_HICOSMO_DE_MASSNU_T: FIXME
 *
 * FIXME
 *
 */
typedef enum _NcHICosmoDEVParams
{
  NC_HICOSMO_DE_MASSNU_M = 0,
  NC_HICOSMO_DE_MASSNU_T,
  NC_HICOSMO_DE_MASSNU_MU,
  NC_HICOSMO_DE_MASSNU_G,   /*< private >*/
  NC_HICOSMO_DE_VPARAM_LEN, /*< skip >*/
} NcHICosmoDEVParams;

#define NC_HICOSMO_DE_DEFAULT_H0        ncm_c_hubble_cte_wmap ()
#define NC_HICOSMO_DE_DEFAULT_OMEGA_C   (0.2568)
#define NC_HICOSMO_DE_DEFAULT_OMEGA_X   (0.70)
#define NC_HICOSMO_DE_DEFAULT_OMEGA_B   (0.0432)
#define NC_HICOSMO_DE_DEFAULT_T_GAMMA0  (2.7245)
#define NC_HICOSMO_DE_DEFAULT_HE_YP     (0.24)
#define NC_HICOSMO_DE_DEFAULT_ENNU      (3.046)
#define NC_HICOSMO_DE_DEFAULT_NU_MASS   (1.0e-5)
#define NC_HICOSMO_DE_DEFAULT_NU_T      (0.71611)
#define NC_HICOSMO_DE_DEFAULT_NU_MU     (0.0)
#define NC_HICOSMO_DE_DEFAULT_NU_G      (1.0)

struct _NcHICosmoDEClass
{
  /*< private >*/
  NcHICosmoClass parent_class;
  NcHICosmoDEFunc1 E2Omega_de;
  NcHICosmoDEFunc1 dE2Omega_de_dz;
  NcHICosmoDEFunc1 d2E2Omega_de_dz2;
  NcHICosmoDEFunc1 w_de;
};

struct _NcHICosmoDE
{
  /*< private >*/
  NcHICosmo parent_instance;
  NcHICosmoDEPrivate *priv;
};

GType nc_hicosmo_de_get_type (void) G_GNUC_CONST;

void nc_hicosmo_de_set_wmap5_params (NcHICosmoDE *cosmo_de);
void nc_hicosmo_de_omega_x2omega_k (NcHICosmoDE *cosmo_de);
void nc_hicosmo_de_cmb_params (NcHICosmoDE *cosmo_de);
void nc_hicosmo_de_new_add_bbn (NcmLikelihood *lh);

void nc_hicosmo_de_set_E2Omega_de_impl (NcHICosmoDEClass *cosmo_de_class, NcHICosmoDEFunc1 f);
void nc_hicosmo_de_set_dE2Omega_de_dz_impl (NcHICosmoDEClass *cosmo_de_class, NcHICosmoDEFunc1 f);
void nc_hicosmo_de_set_d2E2Omega_de_dz2_impl (NcHICosmoDEClass *cosmo_de_class, NcHICosmoDEFunc1 f);
void nc_hicosmo_de_set_w_de_impl (NcHICosmoDEClass *cosmo_de_class, NcHICosmoDEFunc1 f);

G_INLINE_FUNC gdouble nc_hicosmo_de_E2Omega_de (NcHICosmoDE *cosmo_de, gdouble z);
G_INLINE_FUNC gdouble nc_hicosmo_de_dE2Omega_de_dz (NcHICosmoDE *cosmo_de, gdouble z);
G_INLINE_FUNC gdouble nc_hicosmo_de_d2E2Omega_de_dz2 (NcHICosmoDE *cosmo_de, gdouble z);
G_INLINE_FUNC gdouble nc_hicosmo_de_w_de (NcHICosmoDE *cosmo_de, gdouble z);
G_INLINE_FUNC gdouble nc_hicosmo_de_E2Omega_de_onepw (NcHICosmoDE *cosmo_de, gdouble z);

G_END_DECLS

#endif /* _NC_HICOSMO_DE_H_ */

#ifndef _NC_HICOSMO_DE_INLINE_H_
#define _NC_HICOSMO_DE_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

NCM_MODEL_FUNC1_IMPL (NC_HICOSMO_DE,NcHICosmoDE,nc_hicosmo_de,E2Omega_de,z)
NCM_MODEL_FUNC1_IMPL (NC_HICOSMO_DE,NcHICosmoDE,nc_hicosmo_de,dE2Omega_de_dz,z)
NCM_MODEL_FUNC1_IMPL (NC_HICOSMO_DE,NcHICosmoDE,nc_hicosmo_de,d2E2Omega_de_dz2,z)
NCM_MODEL_FUNC1_IMPL (NC_HICOSMO_DE,NcHICosmoDE,nc_hicosmo_de,w_de,z)

G_INLINE_FUNC gdouble
nc_hicosmo_de_E2Omega_de_onepw (NcHICosmoDE *cosmo_de, gdouble z)
{
  return nc_hicosmo_de_E2Omega_de (cosmo_de, z) * (1.0 + nc_hicosmo_de_w_de (cosmo_de, z));
}

G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_HICOSMO_DE_INLINE_H_ */
