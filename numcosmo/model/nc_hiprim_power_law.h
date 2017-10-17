/***************************************************************************
 *            nc_hiprim_power_law.h
 *
 *  Tue October 27 14:14:03 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hiprim_power_law.h
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

#ifndef _NC_HIPRIM_POWER_LAW_H_
#define _NC_HIPRIM_POWER_LAW_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_c.h>
#include <numcosmo/nc_hiprim.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPRIM_POWER_LAW             (nc_hiprim_power_law_get_type ())
#define NC_HIPRIM_POWER_LAW(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPRIM_POWER_LAW, NcHIPrimPowerLaw))
#define NC_HIPRIM_POWER_LAW_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPRIM_POWER_LAW, NcHIPrimPowerLawClass))
#define NC_IS_HIPRIM_POWER_LAW(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPRIM_POWER_LAW))
#define NC_IS_HIPRIM_POWER_LAW_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPRIM_POWER_LAW))
#define NC_HIPRIM_POWER_LAW_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPRIM_POWER_LAW, NcHIPrimPowerLawClass))

typedef struct _NcHIPrimPowerLawClass NcHIPrimPowerLawClass;
typedef struct _NcHIPrimPowerLaw NcHIPrimPowerLaw;

/**
 * NcHIPrimPowerLawSParams:
 * @NC_HIPRIM_POWER_LAW_LN10E10ASA: Amplitude of the adiabatic scalar mode $\ln(10^{10}\mathcal{A}_\mathrm{s})$
 * @NC_HIPRIM_POWER_LAW_T_SA_RATIO: Tensor-to-scalar ratio $r$
 * @NC_HIPRIM_POWER_LAW_N_SA: Adiabatic scalar spectral index $n_s$
 * @NC_HIPRIM_POWER_LAW_N_T: Tensor spectral index $n_T$
 *
 * Primordial adiabatic scalar power spectrum:
 * $$ \mathcal{P}_{SA}(k) = \mathcal{A}_\mathrm{s}\left(\frac{k}{k_\star}\right)^{n_s -1 }.$$
 * 
 * Primordial tensor power spectrum:
 * $$ \mathcal{P}_T(k) = r \mathcal{A}_\mathrm{s} \left(\frac{k}{k_\star}\right)^{n_T -1 }.$$
 *
 */
typedef enum /*< enum,underscore_name=NC_HIPRIM_POWER_LAW_SPARAMS >*/
{
  NC_HIPRIM_POWER_LAW_LN10E10ASA,
  NC_HIPRIM_POWER_LAW_T_SA_RATIO,
  NC_HIPRIM_POWER_LAW_N_SA,
  NC_HIPRIM_POWER_LAW_N_T, 
  /* < private > */
  NC_HIPRIM_POWER_LAW_SPARAM_LEN, /*< skip >*/
} NcHIPrimPowerLawSParams;

struct _NcHIPrimPowerLawClass
{
  /*< private >*/
  NcHIPrimClass parent_class;
};

struct _NcHIPrimPowerLaw
{
  /*< private >*/
  NcHIPrim parent_instance;

};

GType nc_hiprim_power_law_get_type (void) G_GNUC_CONST;

#define NC_HIPRIM_POWER_LAW_DEFAULT_LN10E10ASA (3.179)
#define NC_HIPRIM_POWER_LAW_DEFAULT_T_SA_RATIO (0.2)
#define NC_HIPRIM_POWER_LAW_DEFAULT_N_SA (0.9742)
#define NC_HIPRIM_POWER_LAW_DEFAULT_N_T (0.0)

NcHIPrimPowerLaw *nc_hiprim_power_law_new (void);

G_END_DECLS

#endif /* _NC_HIPRIM_POWER_LAW_H_ */
