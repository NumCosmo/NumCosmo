/***************************************************************************
 *            nc_hiprim_sbpl.h
 *
 *  Thu July 20 13:06:34 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hiprim_sbpl.h
 * Copyright (C) 2017 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_HIPRIM_SBPL_H_
#define _NC_HIPRIM_SBPL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_c.h>
#include <numcosmo/nc_hiprim.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPRIM_SBPL             (nc_hiprim_sbpl_get_type ())
#define NC_HIPRIM_SBPL(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPRIM_SBPL, NcHIPrimSBPL))
#define NC_HIPRIM_SBPL_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPRIM_SBPL, NcHIPrimSBPLClass))
#define NC_IS_HIPRIM_SBPL(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPRIM_SBPL))
#define NC_IS_HIPRIM_SBPL_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPRIM_SBPL))
#define NC_HIPRIM_SBPL_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPRIM_SBPL, NcHIPrimSBPLClass))

typedef struct _NcHIPrimSBPLClass NcHIPrimSBPLClass;
typedef struct _NcHIPrimSBPL NcHIPrimSBPL;

/**
 * NcHIPrimSBPLSParams:
 * @NC_HIPRIM_SBPL_LN10E10ASA: Amplitude of the adiabatic scalar mode $\ln(10^10A_{SA})$
 * @NC_HIPRIM_SBPL_N_SA: Adiabatic scalar spectral index
 * @NC_HIPRIM_SBPL_DELTA: Exp parameter $\delta$
 * @NC_HIPRIM_SBPL_RA: Amplitude ratio $R_{A_{SA}}$
 * @NC_HIPRIM_SBPL_LNKB: Exp parameter $\ln(k_b)$
 * @NC_HIPRIM_SBPL_LAMBDA: Exp parameter $\lambda$
 * @NC_HIPRIM_SBPL_T_SA_RATIO: Tensor-to-scalar ratio
 * @NC_HIPRIM_SBPL_N_T: Tensor spectral index
 * 
 * FIXME
 * 
 */
typedef enum /*< enum,underscore_name=NC_HIPRIM_SBPL_SPARAMS >*/
{
  NC_HIPRIM_SBPL_LN10E10ASA,
  NC_HIPRIM_SBPL_N_SA,
  NC_HIPRIM_SBPL_DELTA,
  NC_HIPRIM_SBPL_RA, 
  NC_HIPRIM_SBPL_LNKB,
  NC_HIPRIM_SBPL_LAMBDA,
  NC_HIPRIM_SBPL_T_SA_RATIO,
  NC_HIPRIM_SBPL_N_T,        
  /* < private > */
  NC_HIPRIM_SBPL_SPARAM_LEN, /*< skip >*/
} NcHIPrimSBPLSParams;

struct _NcHIPrimSBPLClass
{
  /*< private >*/
  NcHIPrimClass parent_class;
};

struct _NcHIPrimSBPL
{
  /*< private >*/
  NcHIPrim parent_instance;
};

GType nc_hiprim_sbpl_get_type (void) G_GNUC_CONST;

NcHIPrimSBPL *nc_hiprim_sbpl_new (void);

#define NC_HIPRIM_SBPL_DEFAULT_LN10E10ASA (3.179)
#define NC_HIPRIM_SBPL_DEFAULT_N_SA (0.9742)
#define NC_HIPRIM_SBPL_DEFAULT_DELTA (0.0)
#define NC_HIPRIM_SBPL_DEFAULT_RA (0.8)
#define NC_HIPRIM_SBPL_DEFAULT_LNKB (-7.55)
#define NC_HIPRIM_SBPL_DEFAULT_LAMBDA (10.0)
#define NC_HIPRIM_SBPL_DEFAULT_T_SA_RATIO (0.2)
#define NC_HIPRIM_SBPL_DEFAULT_N_T (0.0)

G_END_DECLS

#endif /* _NC_HIPRIM_SBPL_H_ */
