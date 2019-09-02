/***************************************************************************
 *            nc_planck_fi_cor_ttteee.h
 *
 *  Fri April 22 14:47:48 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_planck_fi_cor_ttteee.h
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

#ifndef _NC_PLANCK_FI_COR_TTTEEE_H_
#define _NC_PLANCK_FI_COR_TTTEEE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_likelihood.h>
#include <numcosmo/nc_planck_fi_cor_tt.h>

G_BEGIN_DECLS

#define NC_TYPE_PLANCK_FI_COR_TTTEEE             (nc_planck_fi_cor_ttteee_get_type ())
#define NC_PLANCK_FI_COR_TTTEEE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_PLANCK_FI_COR_TTTEEE, NcPlanckFICorTTTEEE))
#define NC_PLANCK_FI_COR_TTTEEE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_PLANCK_FI_COR_TTTEEE, NcPlanckFICorTTTEEEClass))
#define NC_IS_PLANCK_FI_COR_TTTEEE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_PLANCK_FI_COR_TTTEEE))
#define NC_IS_PLANCK_FI_COR_TTTEEE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_PLANCK_FI_COR_TTTEEE))
#define NC_PLANCK_FI_COR_TTTEEE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_PLANCK_FI_COR_TTTEEE, NcPlanckFICorTTTEEEClass))

typedef struct _NcPlanckFICorTTTEEEClass NcPlanckFICorTTTEEEClass;
typedef struct _NcPlanckFICorTTTEEE NcPlanckFICorTTTEEE;

/**
 * NcPlanckFICorTTTEEESParams:
 * @NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_100: EE amplitude of Galactic dust power at $\ell = 500$ at $100\,$GHz (in $\mu\mathrm{K}^2$)
 * @NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_100_143: EE amplitude of Galactic dust power at $\ell = 500$ at $100 \times 143\,$GHz (in $\mu\mathrm{K}^2$)
 * @NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_100_217: EE amplitude of Galactic dust power at $\ell = 500$ at $100 \times 217\,$GHz (in $\mu\mathrm{K}^2$)
 * @NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_143: EE amplitude of Galactic dust power at $\ell = 500$ at $143\,$GHz (in $\mu\mathrm{K}^2$)
 * @NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_143_217: EE amplitude of Galactic dust power at $\ell = 500$ at $143 \times 217\,$GHz (in $\mu\mathrm{K}^2$)
 * @NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_217: EE amplitude of Galactic dust power at $\ell = 500$ at $217\,$GHz (in $\mu\mathrm{K}^2$)
 * @NC_PLANCK_FI_COR_TTTEEE_galf_EE_index: the dust EE template slope
 * @NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_100: TE amplitude of Galactic dust power at $\ell = 500$ at $100\,$GHz (in $\mu\mathrm{K}^2$)
 * @NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_100_143: TE amplitude of Galactic dust power at $\ell = 500$ at $100 \times 143\,$GHz (in $\mu\mathrm{K}^2$)
 * @NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_100_217: TE amplitude of Galactic dust power at $\ell = 500$ at $100 \times 217\,$GHz (in $\mu\mathrm{K}^2$)
 * @NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_143: TE amplitude of Galactic dust power at $\ell = 500$ at $143\,$GHz (in $\mu\mathrm{K}^2$)
 * @NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_143_217: TE amplitude of Galactic dust power at $\ell = 500$ at $143 \times 217\,$GHz (in $\mu\mathrm{K}^2$)
 * @NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_217: TE amplitude of Galactic dust power at $\ell = 500$ at $217\,$GHz (in $\mu\mathrm{K}^2$)
 * @NC_PLANCK_FI_COR_TTTEEE_galf_TE_index: the dust TE template slope
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_0T_0E: beam-leakage parameter, $\epsilon_0$, $100\times100$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_0T_0E: beam-leakage parameter, $\epsilon_1$, $100\times100$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_0T_0E: beam-leakage parameter, $\epsilon_2$, $100\times100$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_0T_0E: beam-leakage parameter, $\epsilon_3$, $100\times100$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_0T_0E: beam-leakage parameter, $\epsilon_4$, $100\times100$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_0T_1E: beam-leakage parameter, $\epsilon_0$, $100\times143$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_0T_1E: beam-leakage parameter, $\epsilon_1$, $100\times143$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_0T_1E: beam-leakage parameter, $\epsilon_2$, $100\times143$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_0T_1E: beam-leakage parameter, $\epsilon_3$, $100\times143$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_0T_1E: beam-leakage parameter, $\epsilon_4$, $100\times143$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_0T_2E: beam-leakage parameter, $\epsilon_0$, $100\times217$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_0T_2E: beam-leakage parameter, $\epsilon_1$, $100\times217$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_0T_2E: beam-leakage parameter, $\epsilon_2$, $100\times217$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_0T_2E: beam-leakage parameter, $\epsilon_3$, $100\times217$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_0T_2E: beam-leakage parameter, $\epsilon_4$, $100\times217$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_1T_1E: beam-leakage parameter, $\epsilon_0$, $143\times143$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_1T_1E: beam-leakage parameter, $\epsilon_1$, $143\times143$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_1T_1E: beam-leakage parameter, $\epsilon_2$, $143\times143$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_1T_1E: beam-leakage parameter, $\epsilon_3$, $143\times143$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_1T_1E: beam-leakage parameter, $\epsilon_4$, $143\times143$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_1T_2E: beam-leakage parameter, $\epsilon_0$, $143\times217$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_1T_2E: beam-leakage parameter, $\epsilon_1$, $143\times217$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_1T_2E: beam-leakage parameter, $\epsilon_2$, $143\times217$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_1T_2E: beam-leakage parameter, $\epsilon_3$, $143\times217$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_1T_2E: beam-leakage parameter, $\epsilon_4$, $143\times217$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_2T_2E: beam-leakage parameter, $\epsilon_0$, $217\times217$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_2T_2E: beam-leakage parameter, $\epsilon_1$, $217\times217$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_2T_2E: beam-leakage parameter, $\epsilon_2$, $217\times217$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_2T_2E: beam-leakage parameter, $\epsilon_3$, $217\times217$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_2T_2E: beam-leakage parameter, $\epsilon_4$, $217\times217$ TE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_0E_0E: beam-leakage parameter, $\epsilon_0$, $100\times100$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_0E_0E: beam-leakage parameter, $\epsilon_1$, $100\times100$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_0E_0E: beam-leakage parameter, $\epsilon_2$, $100\times100$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_0E_0E: beam-leakage parameter, $\epsilon_3$, $100\times100$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_0E_0E: beam-leakage parameter, $\epsilon_4$, $100\times100$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_0E_1E: beam-leakage parameter, $\epsilon_0$, $100\times143$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_0E_1E: beam-leakage parameter, $\epsilon_1$, $100\times143$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_0E_1E: beam-leakage parameter, $\epsilon_2$, $100\times143$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_0E_1E: beam-leakage parameter, $\epsilon_3$, $100\times143$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_0E_1E: beam-leakage parameter, $\epsilon_4$, $100\times143$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_0E_2E: beam-leakage parameter, $\epsilon_0$, $100\times217$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_0E_2E: beam-leakage parameter, $\epsilon_1$, $100\times217$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_0E_2E: beam-leakage parameter, $\epsilon_2$, $100\times217$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_0E_2E: beam-leakage parameter, $\epsilon_3$, $100\times217$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_0E_2E: beam-leakage parameter, $\epsilon_4$, $100\times217$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_1E_1E: beam-leakage parameter, $\epsilon_0$, $143\times143$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_1E_1E: beam-leakage parameter, $\epsilon_1$, $143\times143$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_1E_1E: beam-leakage parameter, $\epsilon_2$, $143\times143$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_1E_1E: beam-leakage parameter, $\epsilon_3$, $143\times143$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_1E_1E: beam-leakage parameter, $\epsilon_4$, $143\times143$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_1E_2E: beam-leakage parameter, $\epsilon_0$, $143\times217$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_1E_2E: beam-leakage parameter, $\epsilon_1$, $143\times217$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_1E_2E: beam-leakage parameter, $\epsilon_2$, $143\times217$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_1E_2E: beam-leakage parameter, $\epsilon_3$, $143\times217$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_1E_2E: beam-leakage parameter, $\epsilon_4$, $143\times217$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_2E_2E: beam-leakage parameter, $\epsilon_0$, $217\times217$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_2E_2E: beam-leakage parameter, $\epsilon_1$, $217\times217$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_2E_2E: beam-leakage parameter, $\epsilon_2$, $217\times217$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_2E_2E: beam-leakage parameter, $\epsilon_3$, $217\times217$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_2E_2E: beam-leakage parameter, $\epsilon_4$, $217\times217$ EE
 * @NC_PLANCK_FI_COR_TTTEEE_calib_100P: Calibration of the 100 EE spectra
 * @NC_PLANCK_FI_COR_TTTEEE_calib_143P: Calibration of the 143 EE spectra
 * @NC_PLANCK_FI_COR_TTTEEE_calib_217P: Calibration of the 217 EE spectra
 * @NC_PLANCK_FI_COR_TTTEEE_A_pol: Calibration of the polarization relative to the temperature
 * 
 * Planck Foregound and Instrument parameters, compatible with 2013 and 2015
 * releases (see [Planck 2015 results XI (2015)][XPlanckCollaboration2015a]).
 * 
 */
typedef enum /*< enum,underscore_name=NC_PLANCK_FI_COR_TTTEEE_SPARAMS >*/
{
  NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_100 = NC_PLANCK_FI_COR_TT_SPARAM_LEN,
  NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_100_143,
  NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_100_217,
  NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_143,
  NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_143_217,
  NC_PLANCK_FI_COR_TTTEEE_galf_EE_A_217,
  NC_PLANCK_FI_COR_TTTEEE_galf_EE_index,
  NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_100,
  NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_100_143,
  NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_100_217,
  NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_143,
  NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_143_217,
  NC_PLANCK_FI_COR_TTTEEE_galf_TE_A_217,
  NC_PLANCK_FI_COR_TTTEEE_galf_TE_index,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_0T_0E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_0T_0E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_0T_0E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_0T_0E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_0T_0E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_0T_1E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_0T_1E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_0T_1E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_0T_1E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_0T_1E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_0T_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_0T_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_0T_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_0T_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_0T_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_1T_1E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_1T_1E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_1T_1E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_1T_1E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_1T_1E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_1T_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_1T_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_1T_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_1T_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_1T_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_2T_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_2T_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_2T_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_2T_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_2T_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_0E_0E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_0E_0E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_0E_0E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_0E_0E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_0E_0E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_0E_1E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_0E_1E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_0E_1E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_0E_1E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_0E_1E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_0E_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_0E_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_0E_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_0E_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_0E_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_1E_1E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_1E_1E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_1E_1E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_1E_1E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_1E_1E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_1E_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_1E_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_1E_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_1E_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_1E_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_0_2E_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_1_2E_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_2_2E_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_3_2E_2E,
  NC_PLANCK_FI_COR_TTTEEE_bleak_epsilon_4_2E_2E,
  NC_PLANCK_FI_COR_TTTEEE_calib_100P,
  NC_PLANCK_FI_COR_TTTEEE_calib_143P,
  NC_PLANCK_FI_COR_TTTEEE_calib_217P,
  NC_PLANCK_FI_COR_TTTEEE_A_pol,      
  /* < private > */
  NC_PLANCK_FI_COR_TTTEEE_SPARAM_LEN, /*< skip >*/
} NcPlanckFICorTTTEEESParams;

#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_EE_A_100     ( 0.060)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_EE_A_100_143 ( 0.050)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_EE_A_100_217 ( 0.110)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_EE_A_143     ( 0.10)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_EE_A_143_217 ( 0.240)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_EE_A_217     ( 0.72)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_EE_index     (-2.4)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_TE_A_100     ( 0.140)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_TE_A_100_143 ( 0.120)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_TE_A_100_217 ( 0.30)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_TE_A_143     ( 0.240)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_TE_A_143_217 ( 0.60)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_TE_A_217     ( 1.80)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_galf_TE_index     (-2.4)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_0T_0E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_0T_0E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_0T_0E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_0T_0E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_0T_0E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_0T_1E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_0T_1E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_0T_1E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_0T_1E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_0T_1E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_0T_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_0T_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_0T_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_0T_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_0T_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_1T_1E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_1T_1E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_1T_1E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_1T_1E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_1T_1E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_1T_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_1T_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_1T_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_1T_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_1T_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_2T_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_2T_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_2T_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_2T_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_2T_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_0E_0E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_0E_0E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_0E_0E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_0E_0E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_0E_0E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_0E_1E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_0E_1E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_0E_1E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_0E_1E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_0E_1E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_0E_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_0E_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_0E_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_0E_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_0E_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_1E_1E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_1E_1E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_1E_1E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_1E_1E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_1E_1E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_1E_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_1E_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_1E_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_1E_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_1E_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_2E_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_2E_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_2E_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_2E_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_2E_2E (0.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_calib_100P (1.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_calib_143P (1.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_calib_217P (1.0)
#define NC_PLANCK_FI_COR_TTTEEE_DEFAULT_A_pol      (1.0)

struct _NcPlanckFICorTTTEEEClass
{
  /*< private >*/
  NcPlanckFICorTTClass parent_class;
};

struct _NcPlanckFICorTTTEEE
{
  /*< private >*/
  NcPlanckFICorTT parent_instance;
};

GType nc_planck_fi_cor_ttteee_get_type (void) G_GNUC_CONST;

void nc_planck_fi_cor_ttteee_add_galf_priors (NcmLikelihood *lh, NcmVector *mean, NcmVector *sigma);
void nc_planck_fi_cor_ttteee_add_calib_priors (NcmLikelihood *lh, NcmVector *mean, NcmVector *sigma);
void nc_planck_fi_cor_ttteee_add_sz_prior (NcmLikelihood *lh, gdouble f_tSZ, gdouble mean, gdouble sigma);

void nc_planck_fi_cor_ttteee_add_default_galf_priors (NcmLikelihood *lh);
void nc_planck_fi_cor_ttteee_add_default_calib_priors (NcmLikelihood *lh);
void nc_planck_fi_cor_ttteee_add_default_sz_prior (NcmLikelihood *lh);

void nc_planck_fi_cor_ttteee_add_all_default_priors (NcmLikelihood *lh);

G_END_DECLS

#endif /* _NC_PLANCK_FI_COR_TTTEEE_H_ */
