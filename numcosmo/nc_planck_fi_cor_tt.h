/***************************************************************************
 *            nc_planck_fi_cor_tt.h
 *
 *  Thu October 22 16:22:20 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_planck_fi_cor_tt.h
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

#ifndef _NC_PLANCK_FI_COR_TT_H_
#define _NC_PLANCK_FI_COR_TT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_likelihood.h>
#include <numcosmo/nc_planck_fi.h>

G_BEGIN_DECLS

#define NC_TYPE_PLANCK_FI_COR_TT             (nc_planck_fi_cor_tt_get_type ())
#define NC_PLANCK_FI_COR_TT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_PLANCK_FI_COR_TT, NcPlanckFICorTT))
#define NC_PLANCK_FI_COR_TT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_PLANCK_FI_COR_TT, NcPlanckFICorTTClass))
#define NC_IS_PLANCK_FI_COR_TT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_PLANCK_FI_COR_TT))
#define NC_IS_PLANCK_FI_COR_TT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_PLANCK_FI_COR_TT))
#define NC_PLANCK_FI_COR_TT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_PLANCK_FI_COR_TT, NcPlanckFICorTTClass))

typedef struct _NcPlanckFICorTTClass NcPlanckFICorTTClass;
typedef struct _NcPlanckFICorTT NcPlanckFICorTT;

/**
 * NcPlanckFICorTTSParams:
 * @NC_PLANCK_FI_COR_TT_A_cib_217: Contribution of CIB power to $\mathcal{D}^{217}_{3000}$ at the Planck CMB frequency for $217\,$GHz (in $\mu\mathrm{K}^2$)
 * @NC_PLANCK_FI_COR_TT_cib_index: The effective slope of the CIB spectrum
 * @NC_PLANCK_FI_COR_TT_xi_sz_cib: Correlation coefficient between the CIB and tSZ
 * @NC_PLANCK_FI_COR_TT_A_sz: Contribution of tSZ to $\mathcal{D}_{3000}^{143\times 143}$ at $143\,$GHz (in $\mu\mathrm{K}^2$)
 * @NC_PLANCK_FI_COR_TT_ps_A_100_100: Contribution of Poisson point-source power to $\mathcal{D}^{100\times 100}_{3000}$ for Planck (in $\mu\mathrm{K}^2$)
 * @NC_PLANCK_FI_COR_TT_ps_A_143_143: Contribution of Poisson point-source power to $\mathcal{D}^{143\times 143}_{3000}$ for Planck (in $\mu\mathrm{K}^2$)
 * @NC_PLANCK_FI_COR_TT_ps_A_143_217: Contribution of Poisson point-source power to $\mathcal{D}^{143\times 217}_{3000}$ for Planck (in $\mu\mathrm{K}^2$)
 * @NC_PLANCK_FI_COR_TT_ps_A_217_217: Contribution of Poisson point-source power to $\mathcal{D}^{217\times 217}_{3000}$ for Planck (in $\mu\mathrm{K}^2$)
 * @NC_PLANCK_FI_COR_TT_ksz_norm: Contribution of kSZ to $\mathcal{D}_{3000}$ (in $\mu\mathrm{K}^2$)
 * @NC_PLANCK_FI_COR_TT_gal545_A_100: Amplitude of Galactic dust power at $\ell=200$ at $100\,$GHz (in $\mu\mathrm{K}^2$)
 * @NC_PLANCK_FI_COR_TT_gal545_A_143: Amplitude of Galactic dust power at $\ell=200$ at $143\,$GHz (in $\mu\mathrm{K}^2$)
 * @NC_PLANCK_FI_COR_TT_gal545_A_143_217: Amplitude of Galactic dust power at $\ell=200$ at $143\times 217\,$GHz (in $\mu\mathrm{K}^2$)
 * @NC_PLANCK_FI_COR_TT_gal545_A_217: Amplitude of Galactic dust power at $\ell=200$ at $217\,$GHz (in $\mu\mathrm{K}^2$)
 * @NC_PLANCK_FI_COR_TT_calib_100T: Power spectrum calibration for the $100\,$GHz
 * @NC_PLANCK_FI_COR_TT_calib_217T: Power spectrum calibration for the $217\,$GHz
 * @NC_PLANCK_FI_COR_TT_A_planck: Absolute map calibration for Planck
 *
 * Planck Foregound and Instrument parameters, compatible with 2013 and 2015
 * releases (see [Planck 2015 results XI (2015)][XPlanckCollaboration2015a]).
 *
 */
typedef enum /*< enum,underscore_name=NC_PLANCK_FI_COR_TT_SPARAMS >*/
{
  NC_PLANCK_FI_COR_TT_A_cib_217 = 0,
  NC_PLANCK_FI_COR_TT_cib_index,
  NC_PLANCK_FI_COR_TT_xi_sz_cib,
  NC_PLANCK_FI_COR_TT_A_sz,
  NC_PLANCK_FI_COR_TT_ps_A_100_100,
  NC_PLANCK_FI_COR_TT_ps_A_143_143,
  NC_PLANCK_FI_COR_TT_ps_A_143_217,
  NC_PLANCK_FI_COR_TT_ps_A_217_217,
  NC_PLANCK_FI_COR_TT_ksz_norm,
  NC_PLANCK_FI_COR_TT_gal545_A_100,
  NC_PLANCK_FI_COR_TT_gal545_A_143,
  NC_PLANCK_FI_COR_TT_gal545_A_143_217,
  NC_PLANCK_FI_COR_TT_gal545_A_217,
  NC_PLANCK_FI_COR_TT_calib_100T,
  NC_PLANCK_FI_COR_TT_calib_217T,
  NC_PLANCK_FI_COR_TT_A_planck,   
  /* < private > */
  NC_PLANCK_FI_COR_TT_SPARAM_LEN, /*< skip >*/
} NcPlanckFICorTTSParams;

#define NC_PLANCK_FI_COR_TT_DEFAULT_A_cib_217        (100.0)
#define NC_PLANCK_FI_COR_TT_DEFAULT_cib_index        (-1.3)
#define NC_PLANCK_FI_COR_TT_DEFAULT_xi_sz_cib        (0.5)
#define NC_PLANCK_FI_COR_TT_DEFAULT_A_sz             (5.0)
#define NC_PLANCK_FI_COR_TT_DEFAULT_ps_A_100_100     (200.0)
#define NC_PLANCK_FI_COR_TT_DEFAULT_ps_A_143_143     (200.0)
#define NC_PLANCK_FI_COR_TT_DEFAULT_ps_A_143_217     (200.0)
#define NC_PLANCK_FI_COR_TT_DEFAULT_ps_A_217_217     (200.0)
#define NC_PLANCK_FI_COR_TT_DEFAULT_ksz_norm         (5.0)
#define NC_PLANCK_FI_COR_TT_DEFAULT_gal545_A_100     (7.0)
#define NC_PLANCK_FI_COR_TT_DEFAULT_gal545_A_143     (9.0)
#define NC_PLANCK_FI_COR_TT_DEFAULT_gal545_A_143_217 (21.0)
#define NC_PLANCK_FI_COR_TT_DEFAULT_gal545_A_217     (80.0)
#define NC_PLANCK_FI_COR_TT_DEFAULT_calib_100T       (0.9990004)
#define NC_PLANCK_FI_COR_TT_DEFAULT_calib_217T       (0.99501)
#define NC_PLANCK_FI_COR_TT_DEFAULT_A_planck         (1.0)

struct _NcPlanckFICorTTClass
{
  /*< private >*/
  NcPlanckFIClass parent_class;
};

struct _NcPlanckFICorTT
{
  /*< private >*/
  NcPlanckFI parent_instance;
};

GType nc_planck_fi_cor_tt_get_type (void) G_GNUC_CONST;

void nc_planck_fi_cor_tt_add_gal_priors (NcmLikelihood *lh, NcmVector *mean, NcmVector *sigma);
void nc_planck_fi_cor_tt_add_calib_priors (NcmLikelihood *lh, NcmVector *mean, NcmVector *sigma);
void nc_planck_fi_cor_tt_add_sz_prior (NcmLikelihood *lh, gdouble f_tSZ, gdouble mean, gdouble sigma);

void nc_planck_fi_cor_tt_add_default_gal_priors (NcmLikelihood *lh);
void nc_planck_fi_cor_tt_add_default_calib_priors (NcmLikelihood *lh);
void nc_planck_fi_cor_tt_add_default_sz_prior (NcmLikelihood *lh);

void nc_planck_fi_cor_tt_add_all_default_priors (NcmLikelihood *lh);

G_END_DECLS

#endif /* _NC_PLANCK_FI_COR_TT_H_ */
