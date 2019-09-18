/***************************************************************************
 *            nc_hiprim.h
 *
 *  Tue October 27 12:12:11 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hiprim.h
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

#ifndef _NC_HIPRIM_H_
#define _NC_HIPRIM_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/math/ncm_c.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_mset_func.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPRIM             (nc_hiprim_get_type ())
#define NC_HIPRIM(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPRIM, NcHIPrim))
#define NC_HIPRIM_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPRIM, NcHIPrimClass))
#define NC_IS_HIPRIM(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPRIM))
#define NC_IS_HIPRIM_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPRIM))
#define NC_HIPRIM_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPRIM, NcHIPrimClass))

typedef struct _NcHIPrimClass NcHIPrimClass;
/*typedef struct _NcHIPrim NcHIPrim;*/ /* already defined in nc_hicosmo.h! */

/**
 * NcHIPrimImpl:
 * @NC_HIPRIM_IMPL_lnSA_powspec_lnk: Logarithm of the Scalar Adiabatic power spectrum as a function of $\ln(k\mathrm{Mpc})$
 * @NC_HIPRIM_IMPL_lnT_powspec_lnk: Logarithm of the Tensor power spectrum as a function of $\ln(k\mathrm{Mpc})$
 *
 * Methods to be implementd by every primordial model.
 */
typedef enum /*< flags,underscore_name=NC_HIPRIM_IMPL >*/
{
  NC_HIPRIM_IMPL_lnSA_powspec_lnk = 1 << 0,
  NC_HIPRIM_IMPL_lnT_powspec_lnk  = 1 << 1, 
  /* < private > */
  NC_HIPRIM_IMPL_LAST             = 1 << 2, /*< skip >*/
} NcHIPrimImpl;

typedef gdouble (*NcHIPrimFunc1) (NcHIPrim *prim, gdouble lnk);

struct _NcHIPrimClass
{
  /*< private >*/
  NcmModelClass parent_class;
  NcHIPrimFunc1 lnSA_powspec_lnk;
  NcHIPrimFunc1 lnT_powspec_lnk;
  gdouble (*testee) (NcHIPrim *prim, gdouble x);
};

struct _NcHIPrim
{
  /*< private >*/
  NcmModel parent_instance;
  gdouble k_pivot;
  gdouble lnk_pivot;
};

GType nc_hiprim_get_type (void) G_GNUC_CONST;

NCM_MSET_MODEL_DECLARE_ID (nc_hiprim);

NcHIPrim *nc_hiprim_new_from_name (GType parent_type, gchar *prim_name);
NcHIPrim *nc_hiprim_ref (NcHIPrim *prim);
void nc_hiprim_free (NcHIPrim *prim);
void nc_hiprim_clear (NcHIPrim **prim);

void nc_hiprim_log_all_models (GType parent);

void nc_hiprim_set_k_pivot (NcHIPrim *prim, gdouble k_pivot);
gdouble nc_hiprim_get_k_pivot (NcHIPrim *prim);
gdouble nc_hiprim_get_lnk_pivot (NcHIPrim *prim);

NCM_INLINE gdouble nc_hiprim_lnSA_powspec_lnk (NcHIPrim *prim, const gdouble lnk);
NCM_INLINE gdouble nc_hiprim_lnT_powspec_lnk (NcHIPrim *prim, const gdouble lnk);

NCM_INLINE gdouble nc_hiprim_SA_powspec_k (NcHIPrim *prim, const gdouble k);
NCM_INLINE gdouble nc_hiprim_T_powspec_k (NcHIPrim *prim, const gdouble k);

NCM_INLINE gdouble nc_hiprim_SA_Ampl (NcHIPrim *prim);
NCM_INLINE gdouble nc_hiprim_T_Ampl (NcHIPrim *prim);
NCM_INLINE gdouble nc_hiprim_T_SA_ratio (NcHIPrim *prim);

void nc_hiprim_set_lnSA_powspec_lnk_impl (NcHIPrimClass *model_class, NcHIPrimFunc1 f);
void nc_hiprim_set_lnT_powspec_lnk_impl (NcHIPrimClass *model_class, NcHIPrimFunc1 f);

#define NC_HIPRIM_DEFAULT_K_PIVOT (0.05) /* In units of 1/Mpc */
#define NC_HIPRIM_DEFAULT_PARAMS_RELTOL (1e-7)
#define NC_HIPRIM_DEFAULT_PARAMS_ABSTOL (0.0)

G_END_DECLS

#endif /* _NC_HIPRIM_H_ */

#ifndef _NC_HIPRIM_INLINE_H_
#define _NC_HIPRIM_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_MODEL_FUNC1_IMPL (NC_HIPRIM,NcHIPrim,nc_hiprim,lnSA_powspec_lnk,lnk)
NCM_MODEL_FUNC1_IMPL (NC_HIPRIM,NcHIPrim,nc_hiprim,lnT_powspec_lnk,lnk)

NCM_INLINE gdouble
nc_hiprim_SA_powspec_k (NcHIPrim *prim, const gdouble k)
{
  return exp (nc_hiprim_lnSA_powspec_lnk (prim, log (k)));
}

NCM_INLINE gdouble
nc_hiprim_T_powspec_k (NcHIPrim *prim, const gdouble k)
{
  return exp (nc_hiprim_lnT_powspec_lnk (prim, log (k)));
}

NCM_INLINE gdouble
nc_hiprim_SA_Ampl (NcHIPrim *prim)
{
  return nc_hiprim_SA_powspec_k (prim, prim->k_pivot);
}

NCM_INLINE gdouble nc_hiprim_T_Ampl (NcHIPrim *prim)
{
  return nc_hiprim_T_powspec_k (prim, prim->k_pivot);
}

NCM_INLINE gdouble
nc_hiprim_T_SA_ratio (NcHIPrim *prim)
{
  return exp (nc_hiprim_lnT_powspec_lnk (prim, prim->lnk_pivot) - nc_hiprim_lnSA_powspec_lnk (prim, prim->lnk_pivot));
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_HIPRIM_INLINE_H_ */
