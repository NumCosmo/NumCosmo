/***************************************************************************
 *            nc_hiprim_expc.h
 *
 *  Sun May 08 18:13:08 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hiprim_expc.h
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_HIPRIM_EXPC_H_
#define _NC_HIPRIM_EXPC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_c.h>
#include <numcosmo/nc_hiprim.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPRIM_EXPC             (nc_hiprim_expc_get_type ())
#define NC_HIPRIM_EXPC(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPRIM_EXPC, NcHIPrimExpc))
#define NC_HIPRIM_EXPC_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPRIM_EXPC, NcHIPrimExpcClass))
#define NC_IS_HIPRIM_EXPC(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPRIM_EXPC))
#define NC_IS_HIPRIM_EXPC_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPRIM_EXPC))
#define NC_HIPRIM_EXPC_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPRIM_EXPC, NcHIPrimExpcClass))

typedef struct _NcHIPrimExpcClass NcHIPrimExpcClass;
typedef struct _NcHIPrimExpc NcHIPrimExpc;

/**
 * NcHIPrimExpcSParams:
 * @NC_HIPRIM_EXPC_LN10E10ASA: Amplitude of the adiabatic scalar mode $\ln(10^10A_{SA})$
 * @NC_HIPRIM_EXPC_N_SA: Adiabatic scalar spectral index
 * @NC_HIPRIM_EXPC_LAMBDAC: Exp parameter $\lambda_c$
 * @NC_HIPRIM_EXPC_LNKC: Exp parameter $\ln(k_c)$
 * @NC_HIPRIM_EXPC_C: Exp cut parameter $c$
 * @NC_HIPRIM_EXPC_T_SA_RATIO: Tensor-to-scalar ratio
 * @NC_HIPRIM_EXPC_N_T: Tensor spectral index
 * 
 * FIXME
 * 
 */
typedef enum /*< enum,underscore_name=NC_HIPRIM_EXPC_SPARAMS >*/
{
  NC_HIPRIM_EXPC_LN10E10ASA,
  NC_HIPRIM_EXPC_N_SA,
  NC_HIPRIM_EXPC_LAMBDAC,
  NC_HIPRIM_EXPC_LNKC,
  NC_HIPRIM_EXPC_C,          
  NC_HIPRIM_EXPC_T_SA_RATIO,
  NC_HIPRIM_EXPC_N_T,        
  /* < private > */
  NC_HIPRIM_EXPC_SPARAM_LEN, /*< skip >*/
} NcHIPrimExpcSParams;

struct _NcHIPrimExpcClass
{
  /*< private >*/
  NcHIPrimClass parent_class;
};

struct _NcHIPrimExpc
{
  /*< private >*/
  NcHIPrim parent_instance;
};

GType nc_hiprim_expc_get_type (void) G_GNUC_CONST;

NcHIPrimExpc *nc_hiprim_expc_new (void);

#define NC_HIPRIM_EXPC_DEFAULT_LN10E10ASA (3.179)
#define NC_HIPRIM_EXPC_DEFAULT_N_SA (0.9742)
#define NC_HIPRIM_EXPC_DEFAULT_LAMBDAC (0.5)
#define NC_HIPRIM_EXPC_DEFAULT_LNKC (-7.98)
#define NC_HIPRIM_EXPC_DEFAULT_C (0.5)
#define NC_HIPRIM_EXPC_DEFAULT_T_SA_RATIO (0.2)
#define NC_HIPRIM_EXPC_DEFAULT_N_T (0.0)

G_END_DECLS

#endif /* _NC_HIPRIM_EXPC_H_ */
