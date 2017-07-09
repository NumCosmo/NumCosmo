/***************************************************************************
 *            nc_recomb_seager.h
 *
 *  Mon November 05 18:28:36 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_RECOMB_SEAGER_H_
#define _NC_RECOMB_SEAGER_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_recomb.h>

G_BEGIN_DECLS

#define NC_TYPE_RECOMB_SEAGER             (nc_recomb_seager_get_type ())
#define NC_RECOMB_SEAGER(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_RECOMB_SEAGER, NcRecombSeager))
#define NC_RECOMB_SEAGER_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_RECOMB_SEAGER, NcRecombSeagerClass))
#define NC_IS_RECOMB_SEAGER(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_RECOMB_SEAGER))
#define NC_IS_RECOMB_SEAGER_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_RECOMB_SEAGER))
#define NC_RECOMB_SEAGER_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_RECOMB_SEAGER, NcRecombSeagerClass))

typedef struct _NcRecombSeagerClass NcRecombSeagerClass;
typedef struct _NcRecombSeager NcRecombSeager;

struct _NcRecombSeagerClass
{
  /*< private >*/
  NcRecombClass parent_class;
};

/**
 * NcRecombSeagerOpt:
 * @NC_RECOM_SEAGER_OPT_HII_FUDGE: Includes fudge factor in the case_B recombination fitting formulas.
 * @NC_RECOM_SEAGER_OPT_HII_FUDGE_GAUSS_COR: Includes gaussian correction in the case_B recombination fitting formulas.
 * @NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1: Includes Sobolev scape probability for the $2p\,{}^1\\!P_{1} \to 1s\,{}^1\\!S_{0}$.
 * @NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1_CO: Also includes the continum opacity effect due to H.
 * @NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_3P012: Includes Sobolev scape probability for the $2p\,{}^3\\!P_{0,1,2} \to 1s\,{}^1\\!S_{0}$.
 * @NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_3P012_CO: Also includes the continum opacity effect due to H.
 * @NC_RECOM_SEAGER_OPT_ALL: All options.
 *
 * FIXME
 *
 */
typedef enum _NcRecombSeagerOpt
{
  NC_RECOM_SEAGER_OPT_HII_FUDGE             = 1 << 0,
  NC_RECOM_SEAGER_OPT_HII_FUDGE_GAUSS_COR   = 1 << 1,
  NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1      = 1 << 2,
  NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_1P1_CO   = 1 << 3,
  NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_3P012    = 1 << 4,
  NC_RECOM_SEAGER_OPT_HEII_SOBOLEV_3P012_CO = 1 << 5,
  NC_RECOM_SEAGER_OPT_ALL                   = (1 << 6) - 1, /*< private >*/
  NC_RECOM_SEAGER_OPT_LEN,                                  /*< skip >*/
} NcRecombSeagerOpt;

typedef gdouble (*NcRecombSeagerKHI2p2Pmean) (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble H);
typedef gdouble (*NcRecombSeagerKHeI2p) (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble XHI, const gdouble T, const gdouble XHeI, const gdouble H, const gdouble n_H);
typedef void (*NcRecombSeagerKHeI2pGrad) (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble x, const gdouble XHI, const gdouble T, const gdouble XHeI, const gdouble H, const gdouble n_H, gdouble grad[3]);

struct _NcRecombSeager
{
  /*< private >*/
  NcRecomb parent_instance;
  NcRecombSeagerOpt opts;
  NcRecombSeagerKHI2p2Pmean K_HI_2p_2Pmean;
  NcRecombSeagerKHeI2p KX_HeI_2p_1P1;
  NcRecombSeagerKHeI2p KX_HeI_2p_3Pmean;
  NcRecombSeagerKHeI2pGrad KX_HeI_2p_1P1_grad;
  NcRecombSeagerKHeI2pGrad KX_HeI_2p_3Pmean_grad;
  gdouble H_fudge;
  gdouble AGauss1, AGauss2;
  gdouble zGauss1, zGauss2;
  gdouble wGauss1, wGauss2;
  gdouble A2P_s;
  gdouble A2P_t;
  gdouble sigma_He_2P_s;
  gdouble sigma_He_2P_t;
  gdouble Pb, Pb_t;
  gdouble Qb, Qb_t;
  gpointer cvode;
  gboolean init;
  N_Vector y0;
  N_Vector y;
  N_Vector abstol;
  guint n;
  NcmSpline *Xe_s;
  NcmSpline *Xe_reion_s;
  NcmSpline *Xe_recomb_s;
  NcmSpline *XHII_s;
  NcmSpline *XHeII_s;
};

GType nc_recomb_seager_get_type (void) G_GNUC_CONST;

NcRecombSeager *nc_recomb_seager_new (void);
NcRecombSeager *nc_recomb_seager_new_full (gdouble init_frac, gdouble zi, gdouble prec);
NcRecombSeager *nc_recomb_seager_ref (NcRecombSeager *recomb_seager);
void nc_recomb_seager_free (NcRecombSeager *recomb_seager);
void nc_recomb_seager_clear (NcRecombSeager **recomb_seager);

void nc_recomb_seager_set_options (NcRecombSeager *recomb_seager, NcRecombSeagerOpt opts);
void nc_recomb_seager_set_switch (NcRecombSeager *recomb_seager, guint H_switch, guint He_switch);
NcRecombSeagerOpt nc_recomb_seager_get_options (NcRecombSeager *recomb_seager);

gdouble nc_recomb_seager_pequignot_HI_case_B (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble Tm);
gdouble nc_recomb_seager_pequignot_HI_case_B_dTm (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble Tm);
gdouble nc_recomb_seager_hummer_HeI_case_B (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble Tm);
gdouble nc_recomb_seager_hummer_HeI_case_B_dTm (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble Tm);
gdouble nc_recomb_seager_hummer_HeI_case_B_trip (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble Tm);
gdouble nc_recomb_seager_hummer_HeI_case_B_trip_dTm (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble Tm);

#define NC_RECOMB_SEAGER_HUMMER_HEI_CASE_B_T1 (pow (10.0, 5.114))
#define NC_RECOMB_SEAGER_HUMMER_HEI_CASE_B_T2 (pow (10.0, 0.477121))
#define NC_RECOMB_SEAGER_HUMMER_HEI_CASE_B_P (0.711)
#define NC_RECOMB_SEAGER_HUMMER_HEI_CASE_B_Q (pow (10.0, -16.744))
#define NC_RECOMB_SEAGER_HUMMER_HEI_CASE_B_P_TRIP (0.761)
#define NC_RECOMB_SEAGER_HUMMER_HEI_CASE_B_Q_TRIP (pow (10.0, -16.306))

G_END_DECLS

#endif /* _NC_RECOMB_SEAGER_H_ */
