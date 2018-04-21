/***************************************************************************
 *            nc_hipert_boltzmann.h
 *
 *  Sat Oct 25 21:02:53 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_boltzmann.h
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

#ifndef _NC_HIPERT_BOLTZMANN_H_
#define _NC_HIPERT_BOLTZMANN_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model_ctrl.h>
#include <numcosmo/nc_hiprim.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/perturbations/nc_hipert.h>
#include <numcosmo/nc_recomb.h>
#include <numcosmo/nc_scalefactor.h>
#include <numcosmo/data/nc_data_cmb.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_BOLTZMANN             (nc_hipert_boltzmann_get_type ())
#define NC_HIPERT_BOLTZMANN(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT_BOLTZMANN, NcHIPertBoltzmann))
#define NC_HIPERT_BOLTZMANN_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPERT_BOLTZMANN, NcHIPertBoltzmannClass))
#define NC_IS_HIPERT_BOLTZMANN(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT_BOLTZMANN))
#define NC_IS_HIPERT_BOLTZMANN_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPERT_BOLTZMANN))
#define NC_HIPERT_BOLTZMANN_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPERT_BOLTZMANN, NcHIPertBoltzmannClass))

typedef struct _NcHIPertBoltzmannClass NcHIPertBoltzmannClass;
typedef struct _NcHIPertBoltzmann NcHIPertBoltzmann;

typedef void (*NcHIPertBoltzmannCreate) (NcHIPertBoltzmann *pb, NcHICosmo *cosmo);
typedef void (*NcHIPertBoltzmannPrepare) (NcHIPertBoltzmann *pb, NcHICosmo *cosmo);
typedef void (*NcHIPertBoltzmannConf) (NcHIPertBoltzmann *pb);
typedef void (*NcHIPertBoltzmannEvol) (NcHIPertBoltzmann *pb, gdouble g);
typedef gboolean (*NcHIPertBoltzmannTest) (NcHIPertBoltzmann *pb);
typedef void (*NcHIPertBoltzmannSources) (NcHIPertBoltzmann *pb, gdouble *S0, gdouble *S1, gdouble *S2);
typedef gdouble (*NcHIPertBoltzmannGet) (NcHIPertBoltzmann *pb);
typedef gdouble (*NcHIPertBoltzmannGetN) (NcHIPertBoltzmann *pb, guint n);
typedef void (*NcHIPertBoltzmannGetCl) (NcHIPertBoltzmann *pb, NcmVector *Cls);

struct _NcHIPertBoltzmannClass
{
  /*< private >*/
  NcHIPertClass parent_class;
  NcHIPertBoltzmannCreate init;
  NcHIPertBoltzmannConf set_opts;
  NcHIPertBoltzmannConf reset;
  NcHIPertBoltzmannEvol evol_step;
  NcHIPertBoltzmannEvol evol;
  NcHIPertBoltzmannPrepare prepare;
  NcHIPertBoltzmannPrepare prepare_if_needed;
  NcHIPertBoltzmannSources get_sources;
  NcHIPertBoltzmannConf print_stats;
  NcHIPertBoltzmannGet get_z;
  NcHIPertBoltzmannGet get_phi;
  NcHIPertBoltzmannGet get_c0;
  NcHIPertBoltzmannGet get_b0;
  NcHIPertBoltzmannGet get_c1;
  NcHIPertBoltzmannGet get_b1;
  NcHIPertBoltzmannGetN get;
  NcHIPertBoltzmannGetN get_theta;
  NcHIPertBoltzmannGetN get_theta_p;
  NcHIPertBoltzmannGetN get_los_theta;
  NcHIPertBoltzmannGetCl get_PHIPHI_Cls;
  NcHIPertBoltzmannGetCl get_TT_Cls;
  NcHIPertBoltzmannGetCl get_EE_Cls;
  NcHIPertBoltzmannGetCl get_BB_Cls;
  NcHIPertBoltzmannGetCl get_TE_Cls;
  NcHIPertBoltzmannGetCl get_TB_Cls;
  NcHIPertBoltzmannGetCl get_EB_Cls;
  NcHIPertBoltzmannConf print_all;
  gpointer data;
};

/**
 * NcHIPertBoltzmannVars:
 * @NC_HIPERT_BOLTZMANN_B0: FIXME
 * @NC_HIPERT_BOLTZMANN_THETA0: FIXME
 * @NC_HIPERT_BOLTZMANN_C0: FIXME
 * @NC_HIPERT_BOLTZMANN_PHI: FIXME
 * @NC_HIPERT_BOLTZMANN_B1: FIXME
 * @NC_HIPERT_BOLTZMANN_THETA1: FIXME
 * @NC_HIPERT_BOLTZMANN_C1: FIXME
 * @NC_HIPERT_BOLTZMANN_THETA2: FIXME
 * @NC_HIPERT_BOLTZMANN_THETA_P0: FIXME
 * @NC_HIPERT_BOLTZMANN_THETA_P1: FIXME
 * @NC_HIPERT_BOLTZMANN_THETA_P2: FIXME
 *
 * FIXME
 */
typedef enum /*< enum,underscore_name=NC_HIPERT_BOLTZMANN_VARS  >*/
{
  NC_HIPERT_BOLTZMANN_B0 = 0,
  NC_HIPERT_BOLTZMANN_THETA0,
  NC_HIPERT_BOLTZMANN_C0,
  NC_HIPERT_BOLTZMANN_PHI,
  NC_HIPERT_BOLTZMANN_B1,
  NC_HIPERT_BOLTZMANN_THETA1,
  NC_HIPERT_BOLTZMANN_C1,
  NC_HIPERT_BOLTZMANN_THETA2,
  NC_HIPERT_BOLTZMANN_THETA_P0,
  NC_HIPERT_BOLTZMANN_THETA_P1,
  NC_HIPERT_BOLTZMANN_THETA_P2, 
  /* < private > */
  NC_HIPERT_BOLTZMANN_LEN,      /*< skip >*/
} NcHIPertBoltzmannVars;

#define NC_HIPERT_BOLTZMANN_BASE_SIZE (NC_HIPERT_BOLTZMANN_THETA2 + 1)

#define NC_HIPERT_BOLTZMANN_dB0     NC_HIPERT_BOLTZMANN_B0
#define NC_HIPERT_BOLTZMANN_V       NC_HIPERT_BOLTZMANN_C1
#define NC_HIPERT_BOLTZMANN_T       NC_HIPERT_BOLTZMANN_B1
#define NC_HIPERT_BOLTZMANN_dTHETA0 NC_HIPERT_BOLTZMANN_THETA0
#define NC_HIPERT_BOLTZMANN_U       NC_HIPERT_BOLTZMANN_THETA1

extern guint _itheta_table[3];
extern guint _itheta_p_table[3];

#define NC_HIPERT_BOLTZMANN_THETA(n)   ((n <= 2) ? (_itheta_table[n])   : (NC_HIPERT_BOLTZMANN_THETA_P2 + 1) + (2 * (n - 3)))
#define NC_HIPERT_BOLTZMANN_THETA_P(n) ((n <= 2) ? (_itheta_p_table[n]) : (NC_HIPERT_BOLTZMANN_THETA_P2 + 1) + (2 * (n - 3) + 1))

struct _NcHIPertBoltzmann
{
  /*< private >*/
  NcHIPert parent_instance;
  NcRecomb *recomb;
  NcHICosmo *cosmo;
  NcScalefactor *a;
  gdouble eta0;
  gdouble lambdai;
  gdouble lambdaf;
  gdouble lambda_opt_cutoff;
  gdouble lambda_rec;
  gdouble lambda_rec_10m2_max[2];
  gdouble lambda;
  NcDataCMBDataType target_Cls;
  gboolean calc_transfer;
  gboolean use_lensed_Cls;
  gboolean use_tensor;
  guint PHIPHI_lmax, TT_lmax, EE_lmax, BB_lmax, TE_lmax, TB_lmax, EB_lmax;
  gboolean tight_coupling;
  NcmModelCtrl *ctrl_cosmo;
  NcmModelCtrl *ctrl_prim;
};

GType nc_hipert_boltzmann_get_type (void) G_GNUC_CONST;

NcHIPertBoltzmann *nc_hipert_boltzmann_ref (NcHIPertBoltzmann *pb);
void nc_hipert_boltzmann_free (NcHIPertBoltzmann *pb);
void nc_hipert_boltzmann_clear (NcHIPertBoltzmann **pb);

void nc_hipert_boltzmann_set_recomb (NcHIPertBoltzmann *pb, NcRecomb *recomb);

void nc_hipert_boltzmann_set_target_Cls (NcHIPertBoltzmann *pb, NcDataCMBDataType tCls);
void nc_hipert_boltzmann_append_target_Cls (NcHIPertBoltzmann *pb, NcDataCMBDataType tCls);
NcDataCMBDataType nc_hipert_boltzmann_get_target_Cls (NcHIPertBoltzmann *pb);

void nc_hipert_boltzmann_set_calc_transfer (NcHIPertBoltzmann *pb, gboolean calc_transfer);
gboolean nc_hipert_boltzmann_get_calc_transfer (NcHIPertBoltzmann *pb);

void nc_hipert_boltzmann_set_lensed_Cls (NcHIPertBoltzmann *pb, gboolean use_lensed_Cls);
gboolean nc_hipert_boltzmann_lensed_Cls (NcHIPertBoltzmann *pb);

void nc_hipert_boltzmann_set_tensor (NcHIPertBoltzmann *pb, gboolean use_tensor);
gboolean nc_hipert_boltzmann_tensor (NcHIPertBoltzmann *pb);

void nc_hipert_boltzmann_set_PHIPHI_lmax (NcHIPertBoltzmann *pb, guint lmax);
void nc_hipert_boltzmann_set_TT_lmax (NcHIPertBoltzmann *pb, guint lmax);
void nc_hipert_boltzmann_set_EE_lmax (NcHIPertBoltzmann *pb, guint lmax);
void nc_hipert_boltzmann_set_BB_lmax (NcHIPertBoltzmann *pb, guint lmax);
void nc_hipert_boltzmann_set_TE_lmax (NcHIPertBoltzmann *pb, guint lmax);
void nc_hipert_boltzmann_set_TB_lmax (NcHIPertBoltzmann *pb, guint lmax);
void nc_hipert_boltzmann_set_EB_lmax (NcHIPertBoltzmann *pb, guint lmax);

guint nc_hipert_boltzmann_get_PHIPHI_lmax (NcHIPertBoltzmann *pb);
guint nc_hipert_boltzmann_get_TT_lmax (NcHIPertBoltzmann *pb);
guint nc_hipert_boltzmann_get_EE_lmax (NcHIPertBoltzmann *pb);
guint nc_hipert_boltzmann_get_BB_lmax (NcHIPertBoltzmann *pb);
guint nc_hipert_boltzmann_get_TE_lmax (NcHIPertBoltzmann *pb);
guint nc_hipert_boltzmann_get_TB_lmax (NcHIPertBoltzmann *pb);
guint nc_hipert_boltzmann_get_EB_lmax (NcHIPertBoltzmann *pb);

void nc_hipert_boltzmann_prepare (NcHIPertBoltzmann *pb, NcHICosmo *cosmo);
void nc_hipert_boltzmann_prepare_if_needed (NcHIPertBoltzmann *pb, NcHICosmo *cosmo);

void nc_hipert_boltzmann_get_PHIPHI_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls);
void nc_hipert_boltzmann_get_TT_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls);
void nc_hipert_boltzmann_get_EE_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls);
void nc_hipert_boltzmann_get_BB_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls);
void nc_hipert_boltzmann_get_TE_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls);
void nc_hipert_boltzmann_get_TB_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls);
void nc_hipert_boltzmann_get_EB_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls);

#define NC_HIPERT_BOLTZMANN_LAMBDA2X(lambda) (exp (-(lambda)))
#define NC_HIPERT_BOLTZMANN_X2LAMBDA(x) (-log (x))

G_END_DECLS

#endif /* _NC_HIPERT_BOLTZMANN_H_ */
