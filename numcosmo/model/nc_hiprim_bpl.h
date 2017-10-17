/***************************************************************************
 *            nc_hiprim_bpl.h
 *
 *  Sun May 08 18:56:20 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hiprim_bpl.h
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

#ifndef _NC_HIPRIM_BPL_H_
#define _NC_HIPRIM_BPL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_c.h>
#include <numcosmo/nc_hiprim.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPRIM_BPL             (nc_hiprim_bpl_get_type ())
#define NC_HIPRIM_BPL(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPRIM_BPL, NcHIPrimBPL))
#define NC_HIPRIM_BPL_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPRIM_BPL, NcHIPrimBPLClass))
#define NC_IS_HIPRIM_BPL(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPRIM_BPL))
#define NC_IS_HIPRIM_BPL_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPRIM_BPL))
#define NC_HIPRIM_BPL_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPRIM_BPL, NcHIPrimBPLClass))

typedef struct _NcHIPrimBPLClass NcHIPrimBPLClass;
typedef struct _NcHIPrimBPL NcHIPrimBPL;

/**
 * NcHIPrimBPLSParams:
 * @NC_HIPRIM_BPL_LN10E10ASA: Amplitude of the adiabatic scalar mode $\ln(10^10A_{SA})$
 * @NC_HIPRIM_BPL_N_SA: Adiabatic scalar spectral index
 * @NC_HIPRIM_BPL_DELTA: Exp parameter $\delta$
 * @NC_HIPRIM_BPL_LNKB: Exp parameter $\ln(k_b)$
 * @NC_HIPRIM_BPL_T_SA_RATIO: Tensor-to-scalar ratio
 * @NC_HIPRIM_BPL_N_T: Tensor spectral index
 * 
 * FIXME
 * 
 */
typedef enum /*< enum,underscore_name=NC_HIPRIM_BPL_SPARAMS >*/
{
  NC_HIPRIM_BPL_LN10E10ASA,
  NC_HIPRIM_BPL_N_SA,
  NC_HIPRIM_BPL_DELTA,
  NC_HIPRIM_BPL_LNKB,       
  NC_HIPRIM_BPL_T_SA_RATIO,
  NC_HIPRIM_BPL_N_T,        
  /* < private > */
  NC_HIPRIM_BPL_SPARAM_LEN, /*< skip >*/
} NcHIPrimBPLSParams;

struct _NcHIPrimBPLClass
{
  /*< private >*/
  NcHIPrimClass parent_class;
};

struct _NcHIPrimBPL
{
  /*< private >*/
  NcHIPrim parent_instance;
};

GType nc_hiprim_bpl_get_type (void) G_GNUC_CONST;

NcHIPrimBPL *nc_hiprim_bpl_new (void);

#define NC_HIPRIM_BPL_DEFAULT_LN10E10ASA (3.179)
#define NC_HIPRIM_BPL_DEFAULT_N_SA (0.9742)
#define NC_HIPRIM_BPL_DEFAULT_DELTA (1.14)
#define NC_HIPRIM_BPL_DEFAULT_LNKB (-7.55)
#define NC_HIPRIM_BPL_DEFAULT_T_SA_RATIO (0.2)
#define NC_HIPRIM_BPL_DEFAULT_N_T (0.0)

G_END_DECLS

#endif /* _NC_HIPRIM_BPL_H_ */
