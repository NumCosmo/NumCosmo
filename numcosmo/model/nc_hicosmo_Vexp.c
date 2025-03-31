/***************************************************************************
 *            nc_hicosmo_Vexp.c
 *
 *  Fri October 28 13:27:53 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2016 <vitenti@uel.br>
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

/**
 * NcHICosmoVexp:
 *
 * Single scalar field with an exponential potential.
 *
 * Bounce cosmological model assuming a single scalar field with an exponential
 * potential. For details see [Bacalhau et al. (2017)][XBacalhau2017].
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hicosmo_Vexp.h"
#include "perturbations/nc_hipert_adiab.h"
#include "perturbations/nc_hipert_gw.h"
#include "perturbations/nc_hipert_em.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "nc_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <nvector/nvector_serial.h>
#include <cvode/cvode.h>
#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#define SUN_DENSE_ACCESS SM_ELEMENT_D

#include <sundials/sundials_types.h>

#endif /* NUMCOSMO_GIR_SCAN */

static void nc_hipert_iadiab_interface_init (NcHIPertIAdiabInterface *iface);
static void nc_hipert_igw_interface_init (NcHIPertIGWInterface *iface);
static void nc_hipert_iem_interface_init (NcHIPertIEMInterface *iface);

struct _NcHICosmoVexpPrivate
{
  gpointer cvode_qt;
  gpointer cvode_clp;
  gpointer cvode_clm;
  gboolean qt_init;
  gboolean clm_init;
  gboolean clp_init;
  gboolean glue_de;
  gboolean set_xb_max;
  SUNContext sunctx;
  N_Vector y_qt;
  N_Vector ydot_qt;
  N_Vector y_cl;
  SUNMatrix A;
  SUNLinearSolver LS;
  gint cl_bc, cl_be;
  gdouble RH_lp;
  gdouble alpha_b;
  gdouble a_0de;
  gdouble a_0c, a_0e;
  gdouble qc, qe;
  gdouble Ec, Ee;
  gdouble alpha_qc, alpha_qe;
  gdouble alpha_0c, alpha_0e;
  gdouble tau_x0;
  gdouble tau_x0_i;
  gdouble tau_x0_f;
  gdouble tau_qt_c;
  gdouble tau_qt_e;
  gdouble c1c, c1e;
  gdouble c2c, c2e;
  GArray *evol_c;
  GArray *evol_e;
  NcmSpline *E2_s;
  NcmSpline *lnqc_mtau;
  NcmSpline *lnqe_tau;
  NcmSpline *phi_tau;
  NcHICosmoVexpEMCoupling em_coupling;
  gdouble (*em_eval_xi) (NcHICosmoVexp *Vexp, const gdouble tau, const gdouble k);
  gdouble (*em_eval_F1) (NcHICosmoVexp *Vexp, const gdouble tau, const gdouble k);
  gdouble (*em_eval_nu) (NcHICosmoVexp *Vexp, const gdouble tau, const gdouble k);
  gdouble (*em_eval_m) (NcHICosmoVexp *Vexp, const gdouble tau, const gdouble k);
};

G_DEFINE_TYPE_WITH_CODE (NcHICosmoVexp, nc_hicosmo_Vexp, NC_TYPE_HICOSMO,
                         G_IMPLEMENT_INTERFACE (NC_TYPE_HIPERT_IADIAB,
                                                nc_hipert_iadiab_interface_init)
                         G_IMPLEMENT_INTERFACE (NC_TYPE_HIPERT_IGW,
                                                nc_hipert_igw_interface_init)
                         G_IMPLEMENT_INTERFACE (NC_TYPE_HIPERT_IEM,
                                                nc_hipert_iem_interface_init)
                         G_ADD_PRIVATE (NcHICosmoVexp)
                        );

enum
{
  PROP_0,
  PROP_GLUE_DE,
  PROP_SET_XB_MAX,
  PROP_EM_COUPLING,
  PROP_SIZE,
};

typedef struct _NcHICosmoVexpState
{
  gdouble tau;
  gdouble phi;
} NcHICosmoVexpState;

static void
nc_hicosmo_Vexp_init (NcHICosmoVexp *Vexp)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  if (SUNContext_Create (SUN_COMM_NULL, &self->sunctx))
    g_error ("ERROR: SUNContext_Create failed\n");

  self->cvode_qt   = CVodeCreate (CV_BDF, self->sunctx);
  self->cvode_clp  = CVodeCreate (CV_BDF, self->sunctx);
  self->cvode_clm  = CVodeCreate (CV_BDF, self->sunctx);
  self->qt_init    = FALSE;
  self->clp_init   = FALSE;
  self->clm_init   = FALSE;
  self->glue_de    = FALSE;
  self->set_xb_max = FALSE;
  self->y_qt       = N_VNew_Serial (2, self->sunctx);
  self->ydot_qt    = N_VNew_Serial (2, self->sunctx);
  self->y_cl       = N_VNew_Serial (2, self->sunctx);
  self->A          = SUNDenseMatrix (2, 2, self->sunctx);
  self->LS         = SUNLinSol_Dense (self->y_qt, self->A, self->sunctx);

  NCM_CVODE_CHECK ((gpointer) self->A, "SUNDenseMatrix", 0, );
  NCM_CVODE_CHECK ((gpointer) self->LS, "SUNDenseLinearSolver", 0, );

  self->RH_lp    = 0.0;
  self->alpha_b  = 0.0;
  self->a_0de    = 0.0;
  self->a_0c     = 0.0;
  self->a_0e     = 0.0;
  self->qc       = 0.0;
  self->qe       = 0.0;
  self->Ec       = 0.0;
  self->Ee       = 0.0;
  self->alpha_qc = 0.0;
  self->alpha_qe = 0.0;
  self->alpha_0c = 0.0;
  self->alpha_0e = 0.0;
  self->tau_x0   = 0.0;
  self->tau_x0_i = 0.0;
  self->tau_x0_f = 0.0;
  self->tau_qt_c = 0.0;
  self->tau_qt_e = 0.0;
  self->c1c      = 0.0;
  self->c1e      = 0.0;
  self->c2c      = 0.0;
  self->c2e      = 0.0;
  self->evol_c   = g_array_sized_new (TRUE, TRUE, sizeof (NcHICosmoVexpState), 1000);
  self->evol_e   = g_array_sized_new (TRUE, TRUE, sizeof (NcHICosmoVexpState), 1000);

  self->E2_s      = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  self->lnqc_mtau = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  self->lnqe_tau  = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  self->phi_tau   = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());

  {
    GArray *mtau  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    GArray *tau   = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    GArray *tau_q = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    GArray *mz    = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);

    GArray *lnqc = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    GArray *lnqe = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    GArray *phi  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    GArray *E2   = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);

    g_array_set_size (mtau,  1000);
    g_array_set_size (tau,   1000);
    g_array_set_size (tau_q, 1000);
    g_array_set_size (mz,    1000);
    g_array_set_size (lnqc,  1000);
    g_array_set_size (lnqe,  1000);
    g_array_set_size (phi,   1000);
    g_array_set_size (E2,    1000);

    {
      NcmVector *mtau_v  = ncm_vector_new_array (mtau);
      NcmVector *tau_v   = ncm_vector_new_array (tau);
      NcmVector *tau_q_v = ncm_vector_new_array (tau_q);
      NcmVector *mz_v    = ncm_vector_new_array (mz);

      NcmVector *lnqc_v = ncm_vector_new_array (lnqc);
      NcmVector *lnqe_v = ncm_vector_new_array (lnqe);
      NcmVector *phi_v  = ncm_vector_new_array (phi);
      NcmVector *E2_v   = ncm_vector_new_array (E2);

      ncm_spline_set (self->lnqc_mtau, mtau_v, lnqc_v, FALSE);
      ncm_spline_set (self->lnqe_tau,   tau_v, lnqe_v, FALSE);
      ncm_spline_set (self->phi_tau,  tau_q_v,  phi_v, FALSE);
      ncm_spline_set (self->E2_s,        mz_v,   E2_v, FALSE);

      g_array_unref (mtau);
      g_array_unref (tau);
      g_array_unref (tau_q);
      g_array_unref (mz);

      g_array_unref (lnqc);
      g_array_unref (lnqe);
      g_array_unref (phi);
      g_array_unref (E2);

      ncm_vector_clear (&mtau_v);
      ncm_vector_clear (&tau_v);
      ncm_vector_clear (&tau_q_v);
      ncm_vector_clear (&mz_v);

      ncm_vector_clear (&lnqc_v);
      ncm_vector_clear (&lnqe_v);
      ncm_vector_clear (&phi_v);
      ncm_vector_clear (&E2_v);
    }
  }

  self->em_coupling = NC_HICOSMO_VEXP_EM_COUPLING_INVALID;
  self->em_eval_xi  = NULL;
  self->em_eval_F1  = NULL;
  self->em_eval_nu  = NULL;
  self->em_eval_m   = NULL;
}

static void
_nc_hicosmo_Vexp_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHICosmoVexp *Vexp               = NC_HICOSMO_VEXP (object);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  g_return_if_fail (NC_IS_HICOSMO_VEXP (object));

  switch (prop_id)
  {
    case PROP_GLUE_DE:
      g_value_set_boolean (value, self->glue_de);
      break;
    case PROP_SET_XB_MAX:
      g_value_set_boolean (value, self->set_xb_max);
      break;
    case PROP_EM_COUPLING:
      g_value_set_enum (value, nc_hicosmo_Vexp_get_em_coupling (Vexp));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_hicosmo_Vexp_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHICosmoVexp *Vexp               = NC_HICOSMO_VEXP (object);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  g_return_if_fail (NC_IS_HICOSMO_VEXP (object));

  switch (prop_id)
  {
    case PROP_GLUE_DE:
      self->glue_de = g_value_get_boolean (value);
      ncm_model_state_mark_outdated (NCM_MODEL (Vexp));
      break;
    case PROP_SET_XB_MAX:
      self->set_xb_max = g_value_get_boolean (value);
      ncm_model_state_mark_outdated (NCM_MODEL (Vexp));
      break;
    case PROP_EM_COUPLING:
      nc_hicosmo_Vexp_set_em_coupling (Vexp, g_value_get_enum (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_hicosmo_Vexp_dispose (GObject *object)
{
  NcHICosmoVexp *Vexp               = NC_HICOSMO_VEXP (object);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  ncm_spline_clear (&self->E2_s);

  ncm_spline_clear (&self->lnqc_mtau);
  ncm_spline_clear (&self->lnqe_tau);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_Vexp_parent_class)->dispose (object);
}

static void
_nc_hicosmo_Vexp_finalize (GObject *object)
{
  NcHICosmoVexp *Vexp               = NC_HICOSMO_VEXP (object);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  if (self->cvode_qt != NULL)
  {
    CVodeFree (&self->cvode_qt);
    self->cvode_qt = NULL;
  }

  if (self->cvode_clp != NULL)
  {
    CVodeFree (&self->cvode_clp);
    self->cvode_clp = NULL;
  }

  if (self->cvode_clm != NULL)
  {
    CVodeFree (&self->cvode_clm);
    self->cvode_clm = NULL;
  }

  if (self->y_qt != NULL)
  {
    N_VDestroy (self->y_qt);
    self->y_qt = NULL;
  }

  if (self->ydot_qt != NULL)
  {
    N_VDestroy (self->ydot_qt);
    self->ydot_qt = NULL;
  }

  if (self->y_cl != NULL)
  {
    N_VDestroy (self->y_cl);
    self->y_cl = NULL;
  }

  if (self->A != NULL)
  {
    SUNMatDestroy (self->A);
    self->A = NULL;
  }

  if (self->LS != NULL)
  {
    SUNLinSolFree (self->LS);
    self->LS = NULL;
  }

  g_clear_pointer (&self->evol_c, g_array_unref);
  g_clear_pointer (&self->evol_e, g_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_Vexp_parent_class)->finalize (object);
}

static gdouble _nc_hicosmo_Vexp_H0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_Vexp_xb (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_Vexp_E2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_Vexp_dE2_dz (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_Vexp_d2E2_dz2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_Vexp_Omega_t0 (NcHICosmo *cosmo);

static void
nc_hicosmo_Vexp_class_init (NcHICosmoVexpClass *klass)
{
  GObjectClass *object_class   = G_OBJECT_CLASS (klass);
  NcHICosmoClass *parent_class = NC_HICOSMO_CLASS (klass);
  NcmModelClass *model_class   = NCM_MODEL_CLASS (klass);

  object_class->dispose     = &_nc_hicosmo_Vexp_dispose;
  object_class->finalize    = &_nc_hicosmo_Vexp_finalize;
  model_class->set_property = &_nc_hicosmo_Vexp_set_property;
  model_class->get_property = &_nc_hicosmo_Vexp_get_property;

  ncm_model_class_set_name_nick (model_class, "V_\\exp", "Vexp");
  ncm_model_class_add_params (model_class, NC_HICOSMO_VEXP_SPARAM_LEN, 0, PROP_SIZE);

  /* Set H_0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_VEXP_H0, "H_0", "H0",
                              10.0, 500.0, 1.0,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_VEXP_DEFAULT_H0,
                              NCM_PARAM_TYPE_FIXED);
  /* Set Omega_c0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_VEXP_OMEGA_C, "\\Omega_{c0}", "Omegac",
                              1e-8,  10.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_VEXP_DEFAULT_OMEGA_C,
                              NCM_PARAM_TYPE_FREE);
  /* Set Omega_L0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_VEXP_OMEGA_L, "\\Omega_{\\Lambda0}", "OmegaL",
                              1e-8,  1.0e100, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_VEXP_DEFAULT_OMEGA_L,
                              NCM_PARAM_TYPE_FREE);
  /* Set sigmaphi param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_VEXP_SIGMA_PHI, "\\sigma_{\\phi}", "sigmaphi",
                              1e-8,  100.0e10, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_VEXP_DEFAULT_SIGMA_PHI,
                              NCM_PARAM_TYPE_FREE);
  /* Set cphi param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_VEXP_D_PHI, "d_\\phi", "dphi",
                              -10.0, 10.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_VEXP_DEFAULT_D_PHI,
                              NCM_PARAM_TYPE_FIXED);
  /* Set alpha_0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_VEXP_ALPHA_B, "\\alpha_b", "alphab",
                              0.0,  G_MAXDOUBLE, 1.0e-5,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_VEXP_DEFAULT_ALPHA_B,
                              NCM_PARAM_TYPE_FIXED);
  /* Set xb param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_VEXP_X_B, "x_b", "xb",
                              1.0e10,  1.0e60, 1.0e20,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_VEXP_DEFAULT_X_B,
                              NCM_PARAM_TYPE_FIXED);

  /* Set em_b param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_VEXP_EM_ALPHA, "\\alpha_{\\mathrm{em}}", "alphaem",
                              1.0e-1,  50.0, 1.0e-1,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_VEXP_DEFAULT_EM_ALPHA,
                              NCM_PARAM_TYPE_FIXED);
  /* Set em_beta param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_VEXP_EM_BETA, "\\beta_\\mathrm{em}", "betaem",
                              1.0e-10,  1.0e4, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_VEXP_DEFAULT_EM_BETA,
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  g_object_class_install_property (object_class,
                                   PROP_GLUE_DE,
                                   g_param_spec_boolean ("glue-de",
                                                         NULL,
                                                         "Whether to glue to a DE phase",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SET_XB_MAX,
                                   g_param_spec_boolean ("set-xb-max",
                                                         NULL,
                                                         "Whether to use max xb allowed by the matching",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_EM_COUPLING,
                                   g_param_spec_enum ("em-coupling",
                                                      NULL,
                                                      "Electromagnetic coupling",
                                                      NC_TYPE_HICOSMO_VEXP_EM_COUPLING,
                                                      NC_HICOSMO_VEXP_EM_COUPLING_NONE,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  nc_hicosmo_set_H0_impl        (parent_class, &_nc_hicosmo_Vexp_H0);
  nc_hicosmo_set_xb_impl        (parent_class, &_nc_hicosmo_Vexp_xb);
  nc_hicosmo_set_E2_impl        (parent_class, &_nc_hicosmo_Vexp_E2);
  nc_hicosmo_set_Omega_t0_impl  (parent_class, &_nc_hicosmo_Vexp_Omega_t0);

  nc_hicosmo_set_dE2_dz_impl    (parent_class, &_nc_hicosmo_Vexp_dE2_dz);
  nc_hicosmo_set_d2E2_dz2_impl  (parent_class, &_nc_hicosmo_Vexp_d2E2_dz2);
}

void _nc_hicosmo_Vexp_init_x_alpha_recur (NcmVector *alpha_x_series, NcmVector *x_alpha_series, const gdouble pf, guint nmax);

static gdouble
_nc_hicosmo_Vexp_init_x_alpha_series (const gdouble dalpha)
{
  const gdouble an[25] = {
    +2.1213203435596425732,
    -3.1819805153394638598,
    +0.0,
    +14.318912319027587369,
    -38.661063261374485897,
    +30.069715869957933475,
    +133.47271840236429655,
    -562.31135761409854934,
    +852.95970820407459509,
    +879.60822680923629430,
    -7833.1017094058379035,
    +17969.469996896215007,
    -5132.5309707529267568,
    -97212.641401182025488,
    +325940.66821809623956,
    -370774.67123243092981,
    -949452.26167748074539,
    +5.2556176004654939455e6,
    -9.9203647608405914131e6,
    -3.5322017625266955972e6,
    +7.4749084074478288480e7,
    -2.0578628366558384105e8,
    +1.4591930291931796152e8,
    +8.8839574285216097773e8,
    -3.6743882793601701195e9
  };
  register gint i;
  register gdouble res = an[24];

  for (i = 23; i >= 0; i--)
  {
#ifdef HAVE_FMA
    res = fma (dalpha, res, an[i]);
#else
    res = dalpha * res + an[i];
#endif
  }

  return dalpha * res;
}

#define VECTOR    (NCM_MODEL (cosmo))
#define MACRO_H0  (ncm_model_orig_param_get (VECTOR, NC_HICOSMO_VEXP_H0))
#define OMEGA_C   (ncm_model_orig_param_get (VECTOR, NC_HICOSMO_VEXP_OMEGA_C))
#define OMEGA_L   (ncm_model_orig_param_get (VECTOR, NC_HICOSMO_VEXP_OMEGA_L))
#define SIGMA_PHI (ncm_model_orig_param_get (VECTOR, NC_HICOSMO_VEXP_SIGMA_PHI))
#define D_PHI     (ncm_model_orig_param_get (VECTOR, NC_HICOSMO_VEXP_D_PHI))
#define ALPHA_B   (ncm_model_orig_param_get (VECTOR, NC_HICOSMO_VEXP_ALPHA_B))
#define X_B       (ncm_model_orig_param_get (VECTOR, NC_HICOSMO_VEXP_X_B))
#define EM_ALPHA  (ncm_model_orig_param_get (VECTOR, NC_HICOSMO_VEXP_EM_ALPHA))
#define EM_BETA   (ncm_model_orig_param_get (VECTOR, NC_HICOSMO_VEXP_EM_BETA))

#define LAMBDAp (1.0 + 1.0 / M_SQRT2)
#define LAMBDAm (1.0 - 1.0 / M_SQRT2)

static gdouble _nc_hicosmo_Vexp_adiab_eval_xi (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_Vexp_adiab_eval_F1 (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_Vexp_adiab_eval_nu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_Vexp_adiab_eval_m (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_Vexp_adiab_eval_unit (NcHIPertIAdiab *iad);
static gdouble _nc_hicosmo_Vexp_adiab_eval_x (NcHIPertIAdiab *iad, const gdouble tau);
static gdouble _nc_hicosmo_Vexp_adiab_eval_p2Psi (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_Vexp_adiab_eval_p2drho (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);

static void
nc_hipert_iadiab_interface_init (NcHIPertIAdiabInterface *iface)
{
  iface->eval_xi     = &_nc_hicosmo_Vexp_adiab_eval_xi;
  iface->eval_F1     = &_nc_hicosmo_Vexp_adiab_eval_F1;
  iface->eval_nu     = &_nc_hicosmo_Vexp_adiab_eval_nu;
  iface->eval_m      = &_nc_hicosmo_Vexp_adiab_eval_m;
  iface->eval_unit   = &_nc_hicosmo_Vexp_adiab_eval_unit;
  iface->eval_x      = &_nc_hicosmo_Vexp_adiab_eval_x;
  iface->eval_p2Psi  = &_nc_hicosmo_Vexp_adiab_eval_p2Psi;
  iface->eval_p2drho = &_nc_hicosmo_Vexp_adiab_eval_p2drho;
}

static gdouble _nc_hicosmo_Vexp_gw_eval_xi (NcHIPertIGW *igw, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_Vexp_gw_eval_F1 (NcHIPertIGW *igw, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_Vexp_gw_eval_nu (NcHIPertIGW *igw, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_Vexp_gw_eval_m (NcHIPertIGW *igw, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_Vexp_gw_eval_unit (NcHIPertIGW *igw);
static gdouble _nc_hicosmo_Vexp_gw_eval_x (NcHIPertIGW *igw, const gdouble tau);

static void
nc_hipert_igw_interface_init (NcHIPertIGWInterface *iface)
{
  iface->eval_xi   = &_nc_hicosmo_Vexp_gw_eval_xi;
  iface->eval_F1   = &_nc_hicosmo_Vexp_gw_eval_F1;
  iface->eval_nu   = &_nc_hicosmo_Vexp_gw_eval_nu;
  iface->eval_m    = &_nc_hicosmo_Vexp_gw_eval_m;
  iface->eval_unit = &_nc_hicosmo_Vexp_gw_eval_unit;
  iface->eval_x    = &_nc_hicosmo_Vexp_gw_eval_x;
}

static gdouble _nc_hicosmo_Vexp_em_eval_xi (NcHIPertIEM *iem, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_Vexp_em_eval_F1 (NcHIPertIEM *iem, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_Vexp_em_eval_nu (NcHIPertIEM *iem, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_Vexp_em_eval_m (NcHIPertIEM *iem, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_Vexp_em_eval_unit (NcHIPertIEM *iem);
static gdouble _nc_hicosmo_Vexp_em_eval_x (NcHIPertIEM *iem, const gdouble tau);
static gdouble _nc_hicosmo_Vexp_em_eval_lapse (NcHIPertIEM *iem, const gdouble tau);

static void
nc_hipert_iem_interface_init (NcHIPertIEMInterface *iface)
{
  iface->eval_xi    = &_nc_hicosmo_Vexp_em_eval_xi;
  iface->eval_F1    = &_nc_hicosmo_Vexp_em_eval_F1;
  iface->eval_nu    = &_nc_hicosmo_Vexp_em_eval_nu;
  iface->eval_m     = &_nc_hicosmo_Vexp_em_eval_m;
  iface->eval_unit  = &_nc_hicosmo_Vexp_em_eval_unit;
  iface->eval_x     = &_nc_hicosmo_Vexp_em_eval_x;
  iface->eval_lapse = &_nc_hicosmo_Vexp_em_eval_lapse;
}

static gdouble
_nc_hicosmo_Vexp_q (const gdouble epsilon, const gint cl_b)
{
  if (cl_b > 0)
    return fabs (epsilon) / (LAMBDAm - fabs (epsilon));
  else
    return fabs (epsilon) / (LAMBDAp - fabs (epsilon));
}

static gdouble
_nc_hicosmo_Vexp_x_q (const gdouble q, const gint cl_b)
{
  if (cl_b > 0)
    return 1.0 - LAMBDAm / (1.0 / q + 1.0);
  else
    return -1.0 + LAMBDAp / (1.0 / q + 1.0);
}

static gdouble
_nc_hicosmo_Vexp_1mx_q (const gdouble q, const gint cl_b)
{
  if (cl_b > 0)
    return LAMBDAm / (1.0 / q + 1.0);
  else
    return 2.0 - LAMBDAp / (1.0 / q + 1.0);
}

static gdouble
_nc_hicosmo_Vexp_1px_q (const gdouble q, const gint cl_b)
{
  if (cl_b > 0)
    return 2.0 - LAMBDAm / (1.0 / q + 1.0);
  else
    return LAMBDAp / (1.0 / q + 1.0);
}

static gdouble
_nc_hicosmo_Vexp_1sqrt2_x_q (const gdouble q, const gint cl_b)
{
  if (cl_b > 0)
    return -LAMBDAm * (1.0 / (1.0 + q));
  else
    return +LAMBDAp * (1.0 / (1.0 + q));
}

static gdouble
_nc_hicosmo_Vexp_ac_a0c_q (NcHICosmoVexp *Vexp, const gdouble q, const gint cl_b)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);
  const gdouble onemx               = _nc_hicosmo_Vexp_1mx_q (q, cl_b);
  const gdouble onepx               = _nc_hicosmo_Vexp_1px_q (q, cl_b);
  const gdouble onesqrt2mx          = _nc_hicosmo_Vexp_1sqrt2_x_q (q, cl_b);

  return cbrt (self->c2c * pow (onemx, LAMBDAp) * pow (onepx, LAMBDAm) / gsl_pow_2 (onesqrt2mx));
}

static gdouble
_nc_hicosmo_Vexp_ae_a0e_q (NcHICosmoVexp *Vexp, const gdouble q, const gint cl_b)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);
  const gdouble onemx               = _nc_hicosmo_Vexp_1mx_q (q, cl_b);
  const gdouble onepx               = _nc_hicosmo_Vexp_1px_q (q, cl_b);
  const gdouble onesqrt2mx          = _nc_hicosmo_Vexp_1sqrt2_x_q (q, cl_b);

  return cbrt (self->c2e * pow (onemx, LAMBDAp) * pow (onepx, LAMBDAm) / gsl_pow_2 (onesqrt2mx));
}

static gdouble
_nc_hicosmo_Vexp_Hc_H0_q (NcHICosmoVexp *Vexp, const gdouble q, const gint cl_b)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);
  const gdouble onemx               = _nc_hicosmo_Vexp_1mx_q (q, cl_b);
  const gdouble onepx               = _nc_hicosmo_Vexp_1px_q (q, cl_b);
  const gdouble onesqrt2mx          = _nc_hicosmo_Vexp_1sqrt2_x_q (q, cl_b);

  return -sqrt (self->c1c) * pow (onemx, -LAMBDAp) * pow (onepx, -LAMBDAm) * fabs (onesqrt2mx) / self->c2c;
}

static gdouble
_nc_hicosmo_Vexp_He_H0_q (NcHICosmoVexp *Vexp, const gdouble q, const gint cl_b)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);
  const gdouble onemx               = _nc_hicosmo_Vexp_1mx_q (q, cl_b);
  const gdouble onepx               = _nc_hicosmo_Vexp_1px_q (q, cl_b);
  const gdouble onesqrt2mx          = _nc_hicosmo_Vexp_1sqrt2_x_q (q, cl_b);

  return +sqrt (self->c1e) * pow (onemx, -LAMBDAp) * pow (onepx, -LAMBDAm) * fabs (onesqrt2mx) / self->c2e;
}

static gdouble
_nc_hicosmo_Vexp_cl_dlnx_dalpha (NcHICosmoVexp *Vexp, const gdouble q, const gint cl_b)
{
  const gdouble x          = _nc_hicosmo_Vexp_x_q (q, cl_b);
  const gdouble onemx      = _nc_hicosmo_Vexp_1mx_q (q, cl_b);
  const gdouble onepx      = _nc_hicosmo_Vexp_1px_q (q, cl_b);
  const gdouble onesqrt2mx = _nc_hicosmo_Vexp_1sqrt2_x_q (q, cl_b);

  return 3.0 * onepx * onemx * onesqrt2mx / x;
}

static void
_nc_hicosmo_Vexp_dalpha_dtau_dphi_dtQ (NcHICosmoVexp *Vexp, const gdouble alpha, const gdouble phi, gdouble *dalpha, gdouble *dphi)
{
  NcHICosmo *cosmo = NC_HICOSMO (Vexp);

  const gdouble sigma        = SIGMA_PHI;
  const gdouble sigma2       = sigma * sigma;
  const gdouble d            = D_PHI;
  const gdouble arg_tri      = 2.0 * alpha * d;
  const gdouble arg_hyp      = sigma2 * alpha * phi;
  const gdouble cosh_arg_hyp = cosh (arg_hyp);
  const gdouble tanh_arg_hyp = tanh (arg_hyp);
  const gdouble sin_arg_tri  = sin (arg_tri);
  const gdouble cos_arg_tri  = cos (arg_tri);
  const gdouble sin_cosh     = sin_arg_tri / cosh_arg_hyp;
  const gdouble cos_cosh     = cos_arg_tri / cosh_arg_hyp;

  const gdouble dphi_num   = -alpha * sigma2 * sin_cosh + 2.0 * d * cos_cosh + 2.0 * d;
  const gdouble dalpha_num = phi * sigma2 * sin_cosh + 2.0 * d * tanh_arg_hyp;
  const gdouble den        = 2.0 * (cos_cosh + 1.0);

  dphi[0]   = dphi_num / den;
  dalpha[0] = dalpha_num / den;

  return;
}

static void
_nc_hicosmo_Vexp_H_x (NcHICosmoVexp *Vexp, const gdouble tQ, const gdouble alpha, const gdouble phi, gdouble *H, gdouble *x)
{
  NcHICosmo *cosmo = NC_HICOSMO (Vexp);

  const gdouble sigma        = SIGMA_PHI;
  const gdouble sigma2       = sigma * sigma;
  const gdouble d            = D_PHI;
  const gdouble arg_tri      = 2.0 * alpha * d;
  const gdouble arg_hyp      = sigma2 * alpha * phi;
  const gdouble cosh_arg_hyp = cosh (arg_hyp);
  const gdouble tanh_arg_hyp = tanh (arg_hyp);
  const gdouble sin_arg_tri  = sin (arg_tri);
  const gdouble cos_arg_tri  = cos (arg_tri);
  const gdouble sin_cosh     = sin_arg_tri / cosh_arg_hyp;
  const gdouble cos_cosh     = cos_arg_tri / cosh_arg_hyp;

  const gdouble dphi_num   = -alpha * sigma2 * sin_cosh + 2.0 * d * cos_cosh + 2.0 * d;
  const gdouble dalpha_num = phi * sigma2 * sin_cosh + 2.0 * d * tanh_arg_hyp;
  const gdouble den        = 2.0 * (cos_cosh + 1.0);

  H[0] = dalpha_num / (den * exp (3.0 * alpha));
  x[0] = dphi_num / dalpha_num;

  return;
}

static gdouble
_nc_hicosmo_Vexp_epsilon (NcHICosmoVexp *Vexp, const gdouble tQ, const gdouble alpha, const gdouble phi)
{
  NcHICosmo *cosmo = NC_HICOSMO (Vexp);

  const gdouble sigma        = SIGMA_PHI;
  const gdouble sigma2       = sigma * sigma;
  const gdouble d            = D_PHI;
  const gdouble arg_tri      = 2.0 * alpha * d;
  const gdouble arg_hyp      = sigma2 * alpha * phi;
  const gdouble cosh_arg_hyp = cosh (arg_hyp);
  const gdouble tanh_arg_hyp = tanh (arg_hyp);
  const gdouble sin_arg_tri  = sin (arg_tri);
  const gdouble cos_arg_tri  = cos (arg_tri);
  const gdouble sin_cosh     = sin_arg_tri / cosh_arg_hyp;
  const gdouble cos_cosh     = cos_arg_tri / cosh_arg_hyp;

  const gdouble dphi_num   = -alpha * sigma2 * sin_cosh + 2.0 * d * cos_cosh + 2.0 * d;
  const gdouble dalpha_num = phi * sigma2 * sin_cosh + 2.0 * d * tanh_arg_hyp;

  const gdouble x  = dphi_num / dalpha_num;
  const gdouble x2 = x * x;
  gdouble x2m1;

  if (x2 > 1.1)
  {
    x2m1 = x2 - 1.0;
  }
  else
  {
    const gdouble epsilon_num = (-alpha * sigma2 * sin_cosh + 2.0 * d * cos_cosh) * tanh_arg_hyp - phi * sigma2 * sin_cosh;
    const gdouble epsilon_dem = (phi * sigma2 * sin_cosh + 2.0 * d * tanh_arg_hyp) * tanh_arg_hyp;
    const gdouble epsilon     = epsilon_num / epsilon_dem;

    x2m1 = (epsilon * epsilon + 2.0 * epsilon / tanh_arg_hyp + 1.0 / gsl_pow_2 (tanh_arg_hyp * cosh_arg_hyp));
  }

  return -ncm_util_sqrt1px_m1 (x2m1);
}

static gdouble
_nc_hicosmo_Vexp_qt_dlnx_dalpha (NcHICosmoVexp *Vexp, const gdouble tQ, const gdouble alpha, const gdouble phi)
{
  NcHICosmo *cosmo           = NC_HICOSMO (Vexp);
  const gdouble sigma        = SIGMA_PHI;
  const gdouble sigma2       = sigma * sigma;
  const gdouble d            = D_PHI;
  const gdouble arg_tri      = 2.0 * alpha * d;
  const gdouble arg_hyp      = sigma2 * alpha * phi;
  const gdouble cosh_arg_hyp = cosh (arg_hyp);
  const gdouble tanh_arg_hyp = tanh (arg_hyp);
  const gdouble sin_arg_tri  = sin (arg_tri);
  const gdouble cos_arg_tri  = cos (arg_tri);
  const gdouble sin_cosh     = sin_arg_tri / cosh_arg_hyp;
  const gdouble cos_cosh     = cos_arg_tri / cosh_arg_hyp;

  const gdouble dphi_num   = -alpha * sigma2 * sin_cosh + 2.0 * d * cos_cosh + 2.0 * d;
  const gdouble dalpha_num = phi * sigma2 * sin_cosh + 2.0 * d * tanh_arg_hyp;

  const gdouble x  = dphi_num / dalpha_num;
  const gdouble x2 = x * x;
  gdouble onemx2;

  if (x2 > 1.1)
  {
    onemx2 = 1.0 - x2;
  }
  else
  {
    const gdouble e_num = (-alpha * sigma2 * sin_cosh + 2.0 * d * cos_cosh) * tanh_arg_hyp - phi * sigma2 * sin_cosh;
    const gdouble e_dem = (phi * sigma2 * sin_cosh + 2.0 * d * tanh_arg_hyp) * tanh_arg_hyp;
    const gdouble e     = e_num / e_dem;

    onemx2 = -(e * e + 2.0 * e / tanh_arg_hyp + 1.0 / gsl_pow_2 (tanh_arg_hyp * cosh_arg_hyp));
  }

  {
    const gdouble alpha2 = alpha * alpha;
    const gdouble f1     = +onemx2 * sigma2 * (alpha + phi / x);
    const gdouble f2     = -2.0 * d * sigma2 * alpha * cos_cosh * onemx2 / dphi_num;
    const gdouble f3     = -sin_cosh * (4.0 * d * d + sigma2 * (1.0 + x2 * (1.0 + alpha2 * sigma2)) + sigma2 * sigma2 * phi * (2.0 * x * alpha + phi)) / dphi_num;

    return (f1 + f2 + f3);
  }
}

static gdouble
_nc_hicosmo_Vexp_qt_phi_tau (NcHICosmoVexp *Vexp, const gdouble tau)
{
  NcHICosmo *cosmo                  = NC_HICOSMO (Vexp);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  const gdouble sigma        = SIGMA_PHI;
  const gdouble sigma2       = sigma * sigma;
  const gdouble d            = D_PHI;
  const gdouble arg_tri      = self->alpha_b * d;
  const gdouble sin_2arg_tri = sin (2.0 * arg_tri);
  const gdouble cos_arg_tri  = cos (arg_tri);
  const gdouble tan_arg_tri  = tan (arg_tri);
  const gdouble t1           = sigma2 * (2.0 * d * self->alpha_b + sin_2arg_tri);
  const gdouble t2           = 2.0 * d - self->alpha_b * sigma2 * tan_arg_tri;

  return tau * fabs (cos_arg_tri) * t2 * sqrt (2.0 / (t1 * t2));
}

static gdouble
_nc_hicosmo_Vexp_a_0_by_ea (NcHICosmoVexp *Vexp, const gdouble tQ, const gdouble epsilon, const gdouble a, const gint cl_b)
{
  NcHICosmo *cosmo                  = NC_HICOSMO (Vexp);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  const gdouble d    = D_PHI;
  const gdouble Oc   = OMEGA_C;
  const gdouble LpLm = pow (LAMBDAp, LAMBDAm);
  const gdouble LmLp = pow (LAMBDAm, LAMBDAp);

  if (cl_b > 0)
  {
    const gdouble twoLm = pow (2.0, LAMBDAm);
    const gdouble a_0   = (pow (fabs (epsilon), LAMBDAp / 3.0) / a) * cbrt (gsl_pow_2 (d * self->RH_lp) * twoLm / (LpLm * LmLp * Oc));

    return a_0;
  }
  else
  {
    if (self->glue_de)
    {
      const gdouble twoLp = pow (2.0, LAMBDAp);
      const gdouble a_0   = (a / pow (fabs (epsilon), LAMBDAm / 3.0)) * cbrt (2.0 * gsl_pow_2 (LAMBDAp) / twoLp);

      return a_0;
    }
    else
    {
      const gdouble twoLp = pow (2.0, LAMBDAp);
      const gdouble a_0   = (pow (fabs (epsilon), LAMBDAm / 3.0) / a) * cbrt (gsl_pow_2 (d * self->RH_lp) * twoLp / (LpLm * LmLp * Oc));

      return a_0;
    }
  }
}

static gint
_nc_hicosmo_Vexp_qt_f (sunrealtype tQ, N_Vector y_qt, N_Vector ydot_qt, gpointer f_data)
{
  NcHICosmoVexp *Vexp               = NC_HICOSMO_VEXP (f_data);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  const gdouble alpha = NV_Ith_S (y_qt, 0) + self->alpha_b;
  const gdouble phi   = NV_Ith_S (y_qt, 1);
  gdouble dalpha, dphi;

  _nc_hicosmo_Vexp_dalpha_dtau_dphi_dtQ (Vexp, alpha, phi, &dalpha, &dphi);

  NV_Ith_S (ydot_qt, 0) = dalpha;
  NV_Ith_S (ydot_qt, 1) = dphi;

  return 0;
}

static void
_nc_hicosmo_Vexp_qt_ddfdv (NcHICosmoVexp *Vexp, NcHICosmo *cosmo, const gdouble alpha, const gdouble phi, gdouble J[2][2])
{
  const gdouble sigma        = SIGMA_PHI;
  const gdouble sigma2       = sigma * sigma;
  const gdouble d            = D_PHI;
  const gdouble arg_tri      = 2.0 * alpha * d;
  const gdouble arg_hyp      = sigma2 * alpha * phi;
  const gdouble cosh_arg_hyp = cosh (arg_hyp);
  const gdouble tanh_arg_hyp = tanh (arg_hyp);
  const gdouble sin_arg_tri  = sin (arg_tri);
  const gdouble cos_arg_tri  = cos (arg_tri);
  const gdouble sin_cosh     = sin_arg_tri / cosh_arg_hyp;
  const gdouble cos_cosh     = cos_arg_tri / cosh_arg_hyp;

  const gdouble dphi_num   = -alpha * sigma2 * sin_cosh + 2.0 * d * cos_cosh + 2.0 * d;
  const gdouble dalpha_num = phi * sigma2 * sin_cosh + 2.0 * d * tanh_arg_hyp;
  const gdouble den        = 2.0 * (cos_cosh + 1.0);

  const gdouble dalpha = dalpha_num / den;
  const gdouble dphi   = dphi_num / den;

  J[0][0] = d * sigma2 * phi - dalpha * (-2.0 * d * sin_cosh + sigma2 * phi * tanh_arg_hyp) / (cos_cosh + 1.0);
  J[0][1] = ((sigma2 * sin_cosh + 2.0 * d * alpha * sigma2) * cos_cosh + sigma2 * sin_cosh - phi * sigma2 * sin_cosh * sigma2 * alpha * tanh_arg_hyp + 2.0 * d * alpha * sigma2 / gsl_pow_2 (cosh_arg_hyp)) / (den * (cos_cosh + 1.0));

  J[1][0] = -dphi * (-2.0 * d * sin_cosh + sigma2 * phi * tanh_arg_hyp) / (cos_cosh + 1.0) + (-sigma2 * sin_cosh - 4.0 * d * d * sin_cosh - 2.0 * d * alpha * sigma2 * cos_cosh + 2.0 * d * sigma2 * phi * tanh_arg_hyp) / den;
  J[1][1] = sigma2 * sigma2 * alpha * alpha * sin_cosh * tanh_arg_hyp / (2.0 * gsl_pow_2 (cos_cosh + 1.0));
}

static gdouble
_nc_hicosmo_Vexp_qt_Ricci_scale (NcHICosmoVexp *Vexp, NcHICosmo *cosmo, const gdouble alpha, const gdouble phi)
{
  gdouble ddfdv[2][2];
  gdouble ddalpha, dalpha, dphi;

  _nc_hicosmo_Vexp_qt_ddfdv (Vexp, cosmo, alpha, phi, ddfdv);
  _nc_hicosmo_Vexp_dalpha_dtau_dphi_dtQ (Vexp, alpha, phi, &dalpha, &dphi);

  ddalpha = ddfdv[0][0] * dalpha + ddfdv[0][1] * dphi;

  return exp (3.0 * alpha) / sqrt (6.0 * fabs (ddalpha - dalpha * dalpha));
}

static gint
_nc_hicosmo_Vexp_qt_J (sunrealtype tQ, N_Vector y_qt, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHICosmoVexp *Vexp               = NC_HICOSMO_VEXP (jac_data);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);
  NcHICosmo *cosmo                  = NC_HICOSMO (Vexp);

  const gdouble alpha = NV_Ith_S (y_qt, 0) + self->alpha_b;
  const gdouble phi   = NV_Ith_S (y_qt, 1);
  gdouble ddfdv[2][2];

  _nc_hicosmo_Vexp_qt_ddfdv (Vexp, cosmo, alpha, phi, ddfdv);

  SUN_DENSE_ACCESS (J, 0, 0) = ddfdv[0][0];
  SUN_DENSE_ACCESS (J, 0, 1) = ddfdv[0][1];

  SUN_DENSE_ACCESS (J, 1, 0) = ddfdv[1][0];
  SUN_DENSE_ACCESS (J, 1, 1) = ddfdv[1][1];

  return 0;
}

static gint
_nc_hicosmo_Vexp_qt_root (sunrealtype tQ, N_Vector y_qt, sunrealtype *gout, gpointer user_data)
{
  NcHICosmoVexp *Vexp               = NC_HICOSMO_VEXP (user_data);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);
  NcHICosmo *cosmo                  = NC_HICOSMO (Vexp);

  const gdouble alpha   = NV_Ith_S (y_qt, 0) + self->alpha_b;
  const gdouble phi     = NV_Ith_S (y_qt, 1);
  const gdouble epsilon = _nc_hicosmo_Vexp_epsilon (Vexp, tQ, alpha, phi);
  const gdouble x_b     = X_B;
  const gdouble a_b     = exp (self->alpha_b);
  const gdouble a_0p    = x_b * a_b;
  const gdouble a_0m    = self->glue_de ? self->a_0de : a_0p;

  gdouble H_lp, x;

  _nc_hicosmo_Vexp_H_x (Vexp, tQ, alpha, phi, &H_lp, &x);

  gout[0] = log (fabs (epsilon * x) * 1.0e4);
  gout[1] = log (_nc_hicosmo_Vexp_a_0_by_ea (Vexp, tQ, epsilon, exp (alpha), +1) / a_0p);
  gout[2] = log (_nc_hicosmo_Vexp_a_0_by_ea (Vexp, tQ, epsilon, exp (alpha), -1) / a_0m);
  gout[3] = log (fabs (H_lp * self->RH_lp * x));

  return 0;
}

static gint
_nc_hicosmo_Vexp_clp_f (sunrealtype tQ, N_Vector y_clp, N_Vector ydot_clp, gpointer f_data)
{
  const gdouble q = exp (NV_Ith_S (y_clp, 0));
  const gdouble E = NV_Ith_S (y_clp, 1);
  const gdouble x = _nc_hicosmo_Vexp_x_q (q, +1);

  NV_Ith_S (ydot_clp, 0) = +6.0 * (LAMBDAm + 0.25 * q) / (1.0 + q);
  NV_Ith_S (ydot_clp, 1) = -3.0 * x * x * E;

  return 0;
}

static gint
_nc_hicosmo_Vexp_clm_f (sunrealtype tQ, N_Vector y_clm, N_Vector ydot_clm, gpointer f_data)
{
  const gdouble q = exp (NV_Ith_S (y_clm, 0));
  const gdouble E = NV_Ith_S (y_clm, 1);
  const gdouble x = _nc_hicosmo_Vexp_x_q (q, -1);

  NV_Ith_S (ydot_clm, 0) = +6.0 * (LAMBDAp + 0.25 * q) / (1.0 + q);
  NV_Ith_S (ydot_clm, 1) = -3.0 * x * x * E;

  return 0;
}

static gint
_nc_hicosmo_Vexp_clm_root (sunrealtype tQ, N_Vector y_cl, sunrealtype *gout, gpointer user_data)
{
  const gdouble q = exp (NV_Ith_S (y_cl, 0));

  gout[0] = _nc_hicosmo_Vexp_x_q (q, -1);
  gout[1] = _nc_hicosmo_Vexp_x_q (q, -1) - 0.1;
  gout[2] = _nc_hicosmo_Vexp_x_q (q, -1) + 0.1;

  return 0;
}

#define RELTOL (1.0e-14)

static gdouble
_nc_hicosmo_Vexp_q_smalltq_ralpha (NcHICosmoVexp *Vexp, const gdouble tQ)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);
  NcHICosmo *cosmo                  = NC_HICOSMO (Vexp);

  const gdouble sigma    = SIGMA_PHI;
  const gdouble sigma2   = sigma * sigma;
  const gdouble d        = D_PHI;
  const gdouble tQ2      = tQ * tQ;
  const gdouble arg_trig = d * self->alpha_b;

  return 0.125 * tQ2 * sigma2 * (d * self->alpha_b / gsl_pow_2 (cos (arg_trig)) + tan (arg_trig)) * (2.0 * d - self->alpha_b * sigma2 *  tan (arg_trig));
}

static gdouble
_nc_hicosmo_Vexp_q_smalltq_phi (NcHICosmoVexp *Vexp, const gdouble tQ)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);
  NcHICosmo *cosmo                  = NC_HICOSMO (Vexp);

  const gdouble sigma    = SIGMA_PHI;
  const gdouble sigma2   = sigma * sigma;
  const gdouble d        = D_PHI;
  const gdouble arg_trig = d * self->alpha_b;

  return 0.5 * tQ * (2.0 * d - self->alpha_b * sigma2 *  tan (arg_trig));
}

#define NC_HICOSMO_VEXP_START_T (1.0e-80)
#define NC_HICOSMO_VEXP_START_T_SAVE (NC_HICOSMO_VEXP_START_T * 1.0e10)

static void
_nc_hicosmo_Vexp_init_qt (NcHICosmoVexp *Vexp, const gdouble direction)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);
  const gdouble tQ_i                = NC_HICOSMO_VEXP_START_T * direction;
  gint flag;

  NV_Ith_S (self->y_qt, 0) = _nc_hicosmo_Vexp_q_smalltq_ralpha (Vexp, tQ_i);
  NV_Ith_S (self->y_qt, 1) = _nc_hicosmo_Vexp_q_smalltq_phi (Vexp, tQ_i);

  if (!self->qt_init)
  {
    flag = CVodeInit (self->cvode_qt, &_nc_hicosmo_Vexp_qt_f, tQ_i, self->y_qt);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    flag = CVodeSetMaxStep (self->cvode_qt, G_MAXUINT32);
    NCM_CVODE_CHECK (&flag, "CVodeSetMaxStep", 1, );

    flag = CVodeSStolerances (self->cvode_qt, RELTOL, 0.0);
    NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );

    flag = CVodeSetUserData (self->cvode_qt, Vexp);
    NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

    flag = CVodeSetLinearSolver (self->cvode_qt, self->LS, self->A);
    NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, );

    flag = CVodeSetJacFn (self->cvode_qt, &_nc_hicosmo_Vexp_qt_J);
    NCM_CVODE_CHECK (&flag, "CVodeSetJacFn", 1, );

    flag = CVodeRootInit (self->cvode_qt, 4, &_nc_hicosmo_Vexp_qt_root);
    NCM_CVODE_CHECK (&flag, "CVodeRootInit", 1, );

    flag = CVodeSetInitStep (self->cvode_qt, tQ_i * 1.0e-4);
    NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );

    flag = CVodeSetMaxStep (self->cvode_qt, 1.0);
    NCM_CVODE_CHECK (&flag, "CVodeSetMaxStep", 1, );
  }
  else
  {
    flag = CVodeReInit (self->cvode_qt, tQ_i, self->y_qt);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    flag = CVodeRootInit (self->cvode_qt, 4, &_nc_hicosmo_Vexp_qt_root);
    NCM_CVODE_CHECK (&flag, "CVodeRootInit", 1, );

    flag = CVodeSetInitStep (self->cvode_qt, tQ_i * 1.0e-4);
    NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );
  }
}

static void
_nc_hicosmo_Vexp_init_clp (NcHICosmoVexp *Vexp, gdouble alpha_q, gdouble q0, gdouble E0)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);
  gint flag;

  NV_Ith_S (self->y_cl, 0) = log (q0);
  NV_Ith_S (self->y_cl, 1) = E0;

  if (!self->clp_init)
  {
    flag = CVodeInit (self->cvode_clp, &_nc_hicosmo_Vexp_clp_f, alpha_q, self->y_cl);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    flag = CVodeSetMaxStep (self->cvode_clp, G_MAXUINT32);
    NCM_CVODE_CHECK (&flag, "CVodeSetMaxStep", 1, );

    flag = CVodeSStolerances (self->cvode_clp, RELTOL, 0.0);
    NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );

    flag = CVodeSetUserData (self->cvode_clp, Vexp);
    NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

    flag = CVodeSetLinearSolver (self->cvode_clp, self->LS, self->A);
    NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, );
  }
  else
  {
    flag = CVodeReInit (self->cvode_clp, alpha_q, self->y_cl);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
  }
}

static void
_nc_hicosmo_Vexp_init_clm (NcHICosmoVexp *Vexp, gdouble alpha_q, gdouble q0, gdouble E0)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);
  gint flag;

  NV_Ith_S (self->y_cl, 0) = log (q0);
  NV_Ith_S (self->y_cl, 1) = E0;

  if (!self->clm_init)
  {
    flag = CVodeInit (self->cvode_clm, &_nc_hicosmo_Vexp_clm_f, alpha_q, self->y_cl);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    flag = CVodeSetMaxStep (self->cvode_clm, G_MAXUINT32);
    NCM_CVODE_CHECK (&flag, "CVodeSetMaxStep", 1, );

    flag = CVodeSStolerances (self->cvode_clm, RELTOL, 0.0);
    NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );

    flag = CVodeSetUserData (self->cvode_clm, Vexp);
    NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

    flag = CVodeSetLinearSolver (self->cvode_clm, self->LS, self->A);
    NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, );

    flag = CVodeRootInit (self->cvode_clm, 3, &_nc_hicosmo_Vexp_clm_root);
    NCM_CVODE_CHECK (&flag, "CVodeRootInit", 1, );
  }
  else
  {
    flag = CVodeReInit (self->cvode_clm, alpha_q, self->y_cl);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    flag = CVodeRootInit (self->cvode_clm, 3, &_nc_hicosmo_Vexp_clm_root);
    NCM_CVODE_CHECK (&flag, "CVodeRootInit", 1, );
  }
}

static void
_nc_hicosmo_Vexp_evolve_qt (NcHICosmoVexp *Vexp, gdouble tQ_f)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);
  NcHICosmo *cosmo                  = NC_HICOSMO (Vexp);

  const gdouble tau_sign = GSL_SIGN (tQ_f);
  GArray *earray         = tQ_f > 0.0 ? self->evol_e : self->evol_c;
  gboolean root1_found   = FALSE;
  gint flag;

  _nc_hicosmo_Vexp_init_qt (Vexp, GSL_SIGN (tQ_f));

  flag = CVodeSetStopTime (self->cvode_qt, tQ_f);
  NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

  while (TRUE)
  {
    gdouble tQ;
    gboolean root_found;

    flag = CVode (self->cvode_qt, tQ_f, self->y_qt, &tQ, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "CVode", 1, );

    root_found = (flag == CV_ROOT_RETURN);

    {
      const gdouble ralpha       = NV_Ith_S (self->y_qt, 0);
      const gdouble alpha        = ralpha + self->alpha_b;
      const gdouble phi          = NV_Ith_S (self->y_qt, 1);
      const gdouble tau          = tau_sign * sqrt (2.0 * ralpha);
      const NcHICosmoVexpState s = {tau, phi};

      if (fabs (tQ) > NC_HICOSMO_VEXP_START_T_SAVE)
        g_array_append_val (earray, s);

      if (NC_HICOSMO_VEXP_DEBUG_EVOL_QT)
        printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g %e %e\n",
                tQ,
                tau,
                ralpha,
                phi,
                _nc_hicosmo_Vexp_q_smalltq_ralpha (Vexp, tQ),
                _nc_hicosmo_Vexp_q_smalltq_phi (Vexp, tQ),
                ralpha / _nc_hicosmo_Vexp_q_smalltq_ralpha (Vexp, tQ) - 1.0,
                phi / _nc_hicosmo_Vexp_q_smalltq_phi (Vexp, tQ) - 1.0
               );

      if (root_found)
      {
        gint rootsfound[4] = {0, 0, 0, 0};
        gdouble H_lp, x, set_xb_max;

        _nc_hicosmo_Vexp_H_x (Vexp, tQ, alpha, phi, &H_lp, &x);

        flag = CVodeGetRootInfo (self->cvode_qt, rootsfound);
        NCM_CVODE_CHECK (&flag, "CVodeGetRootInfo", 1, );

        if (NC_HICOSMO_VEXP_DEBUG_EVOL_QT)
          printf ("ROOT: %d %d %d %d\n", rootsfound[0], rootsfound[1], rootsfound[2], rootsfound[3]);

        if (rootsfound[3] != 0)
          g_error ("_nc_hicosmo_Vexp_evolve_qt: quantum phase is too long, unable to match with the classical phase.");

        if (rootsfound[0] != 0)
          root1_found = TRUE;

        set_xb_max = self->set_xb_max && ((x > 0.0) || ((x < 0.0) && !self->glue_de));

        if (
          (!set_xb_max && (((x > 0.0) && (rootsfound[1] != 0)) || ((x < 0.0) && (rootsfound[2] != 0))))
          ||
          (set_xb_max && root1_found)
           )
        {
          const gint cl_b       = GSL_SIGN (x);
          const gdouble epsilon = _nc_hicosmo_Vexp_epsilon (Vexp, tQ, alpha, phi);
          const gdouble a_q     = exp (alpha);
          const gdouble LpLm    = pow (LAMBDAp, LAMBDAm);
          const gdouble LmLp    = pow (LAMBDAm, LAMBDAp);
          const gdouble d       = D_PHI;

          const gdouble a_0     = _nc_hicosmo_Vexp_a_0_by_ea (Vexp, tQ, epsilon, a_q, cl_b);
          const gdouble alpha_0 = log (a_0);
          const gdouble alpha_q = alpha;
          const gdouble q       = _nc_hicosmo_Vexp_q (epsilon, cl_b);
          const gdouble E       = H_lp * self->RH_lp;
          const gdouble c1      = (cl_b > 0) ? gsl_pow_2 (LAMBDAm * d * self->RH_lp) / gsl_pow_6 (a_0) : gsl_pow_2 (LAMBDAp * d * self->RH_lp) / gsl_pow_6 (a_0);
          const gdouble c2      = c1 / (OMEGA_C * LpLm * LmLp);

          if (set_xb_max && (tQ_f < 0.0))
            ncm_model_orig_param_set (VECTOR, NC_HICOSMO_VEXP_X_B, a_0 / exp (self->alpha_b));

          if (!root1_found)
            g_warning ("_nc_hicosmo_Vexp_evolve_qt: imperfect match of the classical solution.");

          if (tQ_f < 0.0)
          {
            self->cl_bc    = cl_b;
            self->a_0c     = a_0;
            self->alpha_0c = alpha_0;
            self->alpha_qc = alpha_q;
            self->qc       = q;
            self->Ec       = E;
            self->c1c      = c1;
            self->c2c      = c2;
            self->tau_qt_c = tau;

            if ((cl_b < 0) && self->glue_de)
            {
              self->c1c = 0.5 * OMEGA_L;
              self->c2c = 0.5;

              /*
               *  self->alpha_0c = log (2.0 * gsl_pow_2 (LAMBDAp * d * self->RH_lp) / OMEGA_L) / 6.0;
               *  self->a_0c     = exp (self->alpha_0c);
               */
            }
          }
          else
          {
            self->cl_be    = cl_b;
            self->a_0e     = a_0;
            self->alpha_0e = alpha_0;
            self->alpha_qe = alpha_q;
            self->qe       = q;
            self->Ee       = E;
            self->c1e      = c1;
            self->c2e      = c2;
            self->tau_qt_e = tau;

            if ((cl_b < 0) && self->glue_de)
            {
              self->c1e = 0.5 * OMEGA_L;
              self->c2e = 0.5;

              /*
               *  self->alpha_0e = log (2.0 * gsl_pow_2 (LAMBDAp * d * self->RH_lp) / OMEGA_L) / 6.0;
               *  self->a_0e     = exp (self->alpha_0c);
               */
            }
          }

          break;
        }
      }

      if (tQ == tQ_f)
        break;
    }
  }
}

static void
_nc_hicosmo_Vexp_evolve_cl_c (NcHICosmoVexp *Vexp)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);
  NcmVector *lnqc_mtau_xv           = ncm_spline_peek_xv (self->lnqc_mtau);
  NcmVector *lnqc_mtau_yv           = ncm_spline_peek_yv (self->lnqc_mtau);
  const gdouble alpha_f             = 200.0 + self->alpha_0c;
  GArray *mtau_a                    = ncm_vector_get_array (lnqc_mtau_xv);
  GArray *lnqc_a                    = ncm_vector_get_array (lnqc_mtau_yv);
  gboolean found_lower              = FALSE;
  gpointer cvode;
  gint flag;

  if (self->cl_bc > 0)
  {
    _nc_hicosmo_Vexp_init_clp (Vexp, self->alpha_qc, self->qc, self->Ec);
    cvode = self->cvode_clp;
  }
  else
  {
    _nc_hicosmo_Vexp_init_clm (Vexp, self->alpha_qc, self->qc, self->Ec);
    cvode = self->cvode_clm;
  }

  flag = CVodeSetStopTime (cvode, alpha_f);
  NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

  g_array_set_size (mtau_a, 0);
  g_array_set_size (lnqc_a, 0);

  while (TRUE)
  {
    gdouble alpha;

    flag = CVode (cvode, alpha_f, self->y_cl, &alpha, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "CVode", 1, );

    {
      const gdouble lnq  = NV_Ith_S (self->y_cl, 0);
      const gdouble En   = NV_Ith_S (self->y_cl, 1);
      const gdouble q    = exp (lnq);
      const gdouble tau  = GSL_SIGN (En) * sqrt (2.0 * (alpha - self->alpha_b));
      const gdouble mtau = -tau;

      if (flag == CV_ROOT_RETURN)
      {
        gint rootsfound[3] = {0, 0, 0};

        flag = CVodeGetRootInfo (cvode, rootsfound);
        NCM_CVODE_CHECK (&flag, "CVodeGetRootInfo", 1, );

        if (rootsfound[0])
        {
          self->tau_x0 = tau;

          if (rootsfound[1] || rootsfound[2])
            g_error ("_nc_hicosmo_Vexp_evolve_cl_c: unknown error!");
        }

        if (rootsfound[1])
        {
          if (found_lower)
          {
            self->tau_x0_i = tau;
          }
          else
          {
            self->tau_x0_f = tau;
            found_lower    = TRUE;
          }

          if (rootsfound[0] || rootsfound[2])
            g_error ("_nc_hicosmo_Vexp_evolve_cl_c: unknown error!");
        }

        if (rootsfound[2])
        {
          if (found_lower)
          {
            self->tau_x0_i = tau;
          }
          else
          {
            self->tau_x0_f = tau;
            found_lower    = TRUE;
          }

          if (rootsfound[0] || rootsfound[1])
            g_error ("_nc_hicosmo_Vexp_evolve_cl_c: unknown error!");
        }
      }

      g_array_append_val (mtau_a, mtau);
      g_array_append_val (lnqc_a,  lnq);

      if (NC_HICOSMO_VEXP_DEBUG_EVOL_CL)
      {
        const gdouble x    = _nc_hicosmo_Vexp_x_q (q, self->cl_bc);
        const gdouble a_a0 = _nc_hicosmo_Vexp_ac_a0c_q (Vexp, q, self->cl_bc);
        const gdouble E    = _nc_hicosmo_Vexp_Hc_H0_q (Vexp, q, self->cl_bc);

        printf ("# C % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15e % 21.15e % 21.15e % 21.15e\n",
                tau, alpha, lnq, x, exp (alpha - self->alpha_0c), a_a0, exp (alpha - self->alpha_0c) / a_a0 - 1.0, E, En, E / En - 1.0);
      }
    }

    if (alpha == alpha_f)
      break;
  }

  ncm_spline_set_array (self->lnqc_mtau, mtau_a, lnqc_a, TRUE);

  g_array_unref (mtau_a);
  g_array_unref (lnqc_a);
}

static void
_nc_hicosmo_Vexp_evolve_cl_e (NcHICosmoVexp *Vexp)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);
  NcmVector *lnqe_tau_xv            = ncm_spline_peek_xv (self->lnqe_tau);
  NcmVector *lnqe_tau_yv            = ncm_spline_peek_yv (self->lnqe_tau);
  NcmVector *E2_s_xv                = ncm_spline_peek_xv (self->E2_s);
  NcmVector *E2_s_yv                = ncm_spline_peek_yv (self->E2_s);
  const gdouble alpha_f             = 200.0 + self->alpha_0e;
  GArray *tau_a                     = ncm_vector_get_array (lnqe_tau_xv);
  GArray *lnqe_a                    = ncm_vector_get_array (lnqe_tau_yv);
  GArray *mz_a                      = ncm_vector_get_array (E2_s_xv);
  GArray *E2_a                      = ncm_vector_get_array (E2_s_yv);
  gdouble last_mz                   = GSL_NEGINF;
  gboolean found_lower              = FALSE;
  gpointer cvode;
  gint flag;

  if (self->cl_be > 0)
  {
    _nc_hicosmo_Vexp_init_clp (Vexp, self->alpha_qe, self->qe, self->Ee);
    cvode = self->cvode_clp;
  }
  else
  {
    _nc_hicosmo_Vexp_init_clm (Vexp, self->alpha_qe, self->qe, self->Ee);
    cvode = self->cvode_clm;
  }

  flag = CVodeSetStopTime (cvode, alpha_f);
  NCM_CVODE_CHECK (&flag, "CVodeRootInit", 1, );

  g_array_set_size (tau_a,  0);
  g_array_set_size (lnqe_a, 0);

  g_array_set_size (mz_a, 0);
  g_array_set_size (E2_a, 0);

  while (TRUE)
  {
    gdouble alpha;

    flag = CVode (cvode, alpha_f, self->y_cl, &alpha, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "CVode", 1, );

    {
      const gdouble lnq    = NV_Ith_S (self->y_cl, 0);
      const gdouble q      = exp (lnq);
      const gdouble En     = NV_Ith_S (self->y_cl, 1);
      const gdouble tau    = GSL_SIGN (En) * sqrt (2.0 * (alpha - self->alpha_b));
      const gdouble lna_a0 = alpha - self->alpha_0e;
      const gdouble mz     = -expm1 (-lna_a0);

      if (flag == CV_ROOT_RETURN)
      {
        gint rootsfound[3] = {0, 0, 0};

        flag = CVodeGetRootInfo (cvode, rootsfound);
        NCM_CVODE_CHECK (&flag, "CVodeGetRootInfo", 1, );

        if (rootsfound[0])
        {
          self->tau_x0 = tau;

          if (rootsfound[1] || rootsfound[2])
            g_error ("_nc_hicosmo_Vexp_evolve_cl_c: unknown error!");
        }

        if (rootsfound[1])
        {
          if (found_lower)
          {
            self->tau_x0_f = tau;
          }
          else
          {
            self->tau_x0_i = tau;
            found_lower    = TRUE;
          }

          if (rootsfound[0] || rootsfound[2])
            g_error ("_nc_hicosmo_Vexp_evolve_cl_c: unknown error!");
        }

        if (rootsfound[2])
        {
          if (found_lower)
          {
            self->tau_x0_f = tau;
          }
          else
          {
            self->tau_x0_i = tau;
            found_lower    = TRUE;
          }

          if (rootsfound[0] || rootsfound[1])
            g_error ("_nc_hicosmo_Vexp_evolve_cl_c: unknown error!");
        }
      }

      g_array_append_val (tau_a,  tau);
      g_array_append_val (lnqe_a, lnq);

      if ((mz - fabs (mz) * 1.0e-7 > last_mz) && (mz < 0.1))
      {
        const gdouble E2 = En * En;

        g_array_append_val (mz_a, mz);
        g_array_append_val (E2_a, E2);

        last_mz = mz;
      }

      if (NC_HICOSMO_VEXP_DEBUG_EVOL_CL)
      {
        const gdouble x    = _nc_hicosmo_Vexp_x_q (q, self->cl_be);
        const gdouble a_a0 = _nc_hicosmo_Vexp_ae_a0e_q (Vexp, q, self->cl_be);
        const gdouble E    = _nc_hicosmo_Vexp_He_H0_q (Vexp, q, self->cl_be);

        printf ("# E % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15e % 21.15e % 21.15e % 21.15e\n",
                tau, alpha, lnq, x, exp (alpha - self->alpha_0e), a_a0, exp (alpha - self->alpha_0e) / a_a0 - 1.0, E, En, E / En - 1.0);
      }
    }

    if (alpha == alpha_f)
      break;
  }

  ncm_spline_set_array (self->lnqe_tau, tau_a, lnqe_a, TRUE);
  ncm_spline_set_array (self->E2_s,      mz_a,   E2_a, TRUE);

  g_array_unref (tau_a);
  g_array_unref (lnqe_a);

  g_array_unref (mz_a);
  g_array_unref (E2_a);
}

static void
_nc_hicosmo_Vexp_prepare (NcHICosmoVexp *Vexp)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);
  NcHICosmo *cosmo                  = NC_HICOSMO (Vexp);
  NcmVector *phi_tau_xv             = ncm_spline_peek_xv (self->phi_tau);
  NcmVector *phi_tau_yv             = ncm_spline_peek_yv (self->phi_tau);
  const gdouble tQ_f                = 1.0e20;

  if (!ncm_model_state_is_update (NCM_MODEL (Vexp)))
  {
    const gdouble d = D_PHI;

    self->RH_lp   = nc_hicosmo_RH_planck (cosmo);
    self->alpha_b = ALPHA_B;
    self->a_0de   = pow (2.0 * gsl_pow_2 (d * LAMBDAp * self->RH_lp) / OMEGA_L, 1.0 / 6.0);

    g_array_set_size (self->evol_c, 0);
    g_array_set_size (self->evol_e, 0);

    _nc_hicosmo_Vexp_evolve_qt (Vexp, -tQ_f);
    _nc_hicosmo_Vexp_evolve_qt (Vexp, +tQ_f);
    _nc_hicosmo_Vexp_evolve_cl_c (Vexp);
    _nc_hicosmo_Vexp_evolve_cl_e (Vexp);

    {
      glong i;

      GArray *tau_qt_a = ncm_vector_get_array (phi_tau_xv);
      GArray *phi_a    = ncm_vector_get_array (phi_tau_yv);

      g_array_set_size (tau_qt_a, 0);
      g_array_set_size (phi_a,    0);

      for (i = self->evol_c->len - 1; i >= 0; i--)
      {
        const NcHICosmoVexpState s = g_array_index (self->evol_c, NcHICosmoVexpState, i);
        const gdouble phi_o_tau    = s.phi / s.tau;

        g_array_append_val (tau_qt_a, s.tau);
        g_array_append_val (phi_a,    phi_o_tau);
      }

      for (i = 0; i < self->evol_e->len; i++)
      {
        const NcHICosmoVexpState s = g_array_index (self->evol_e, NcHICosmoVexpState, i);
        const gdouble phi_o_tau    = s.phi / s.tau;

        g_array_append_val (tau_qt_a, s.tau);
        g_array_append_val (phi_a,    phi_o_tau);
      }

      ncm_spline_set_array (self->phi_tau, tau_qt_a, phi_a, TRUE);

      g_array_unref (tau_qt_a);
      g_array_unref (phi_a);
    }

    ncm_model_state_set_update (NCM_MODEL (Vexp));
  }
  else
  {
    return;
  }
}

/****************************************************************************
 * Normalized Hubble function
 ****************************************************************************/

static gdouble
_nc_hicosmo_Vexp_E2 (NcHICosmo *cosmo, gdouble z)
{
  NcHICosmoVexp *Vexp               = NC_HICOSMO_VEXP (cosmo);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);

  return ncm_spline_eval (self->E2_s, -z);
}

/****************************************************************************
 * Normalized Hubble function redshift derivative
 ****************************************************************************/

static gdouble
_nc_hicosmo_Vexp_dE2_dz (NcHICosmo *cosmo, gdouble z)
{
  NcHICosmoVexp *Vexp               = NC_HICOSMO_VEXP (cosmo);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);

  return -ncm_spline_eval_deriv (self->E2_s, -z);
}

static gdouble
_nc_hicosmo_Vexp_d2E2_dz2 (NcHICosmo *cosmo, gdouble z)
{
  NcHICosmoVexp *Vexp               = NC_HICOSMO_VEXP (cosmo);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);

  return ncm_spline_eval_deriv2 (self->E2_s, -z);
}

/****************************************************************************
 * Simple functions
 ****************************************************************************/
static gdouble
_nc_hicosmo_Vexp_H0 (NcHICosmo *cosmo)
{
  return MACRO_H0;
}

static gdouble
_nc_hicosmo_Vexp_Omega_t0 (NcHICosmo *cosmo)
{
  return OMEGA_C + OMEGA_L;
}

static gdouble
_nc_hicosmo_Vexp_xb (NcHICosmo *cosmo)
{
  return X_B;
}

#define _NC_HICOSMO_VEXP_PHIA NC_HICOSMO_VEXP_START_T

/*
 * Interface implementation -- NcHIPertIAdiab
 */

static gdouble
_nc_hicosmo_Vexp_adiab_eval_xi (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmoVexp *Vexp               = NC_HICOSMO_VEXP (iad);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);

  {
    const gdouble tau2    = tau * tau;
    const gdouble lna_a0e = 0.5 * tau2 + self->alpha_b - self->alpha_0e;

    gdouble mnu;
    guint branch = 0;

    if (tau > self->tau_qt_c)
    {
      branch++;

      if (tau > self->tau_qt_e)
        branch++;
    }

    switch (branch)
    {
      case 0: /* Classical contraction */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqc_mtau, -tau);
        const gdouble q   = exp (lnq);
        const gdouble x   = _nc_hicosmo_Vexp_x_q (q, self->cl_bc);

        if (fabs (x) < 1.0e-1)
        {
          const gdouble dalpha = 0.5 * (tau - self->tau_x0) * (tau + self->tau_x0);
          const gdouble x_s    = _nc_hicosmo_Vexp_init_x_alpha_series (dalpha);

          mnu = k * x_s * x_s * exp (2.0 * lna_a0e);
        }
        else
        {
          mnu = k * x * x * exp (2.0 * lna_a0e);
        }

        break;
      }
      case 1: /* Quantum phase */
      {
        const gdouble alpha = 0.5 * tau * tau + self->alpha_b;
        const gdouble phi   = (fabs (tau) < _NC_HICOSMO_VEXP_PHIA) ? _nc_hicosmo_Vexp_qt_phi_tau (Vexp, tau) : ncm_spline_eval (self->phi_tau, tau) * tau;
        gdouble H_lp, xq;

        _nc_hicosmo_Vexp_H_x (Vexp, 0.0, alpha, phi, &H_lp, &xq);

        mnu = k * xq * xq * exp (2.0 * lna_a0e);

        break;
      }
      case 2: /* Classical expansion */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqe_tau, tau);
        const gdouble q   = exp (lnq);
        const gdouble x   = _nc_hicosmo_Vexp_x_q (q, self->cl_be);

        if (fabs (x) < 1.0e-1)
        {
          const gdouble dalpha = 0.5 * (tau - self->tau_x0) * (tau + self->tau_x0);
          const gdouble x_s    = _nc_hicosmo_Vexp_init_x_alpha_series (dalpha);

          mnu = k * x_s * x_s * exp (2.0 * lna_a0e);
        }
        else
        {
          mnu = k * x * x * exp (2.0 * lna_a0e);
        }

        break;
      }
      default:
        g_assert_not_reached ();
        break;
    }

    return log (mnu);
  }
}

static gdouble
_nc_hicosmo_Vexp_adiab_eval_F1 (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmoVexp *Vexp               = NC_HICOSMO_VEXP (iad);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);

  {
    const gdouble tau2    = tau * tau;
    const gdouble lna_a0e = 0.5 * tau2 + self->alpha_b - self->alpha_0e;
    guint branch          = 0;
    gdouble dlnmnu, nu;

    if (tau > self->tau_qt_c)
    {
      branch++;

      if (tau > self->tau_qt_e)
        branch++;
    }

    switch (branch)
    {
      case 0: /* Classical contraction */
      {
        const gdouble lnq         = ncm_spline_eval (self->lnqc_mtau, -tau);
        const gdouble q           = exp (lnq);
        const gdouble dlnx_dalpha = _nc_hicosmo_Vexp_cl_dlnx_dalpha (Vexp, q, self->cl_bc);
        const gdouble E           = _nc_hicosmo_Vexp_Hc_H0_q (Vexp, q, self->cl_bc);

        nu     = exp (-lna_a0e) * k * tau / E;
        dlnmnu = 2.0 * (dlnx_dalpha + 1.0) * tau;

        break;
      }
      case 1: /* Quantum phase */
      {
        const gdouble alpha = 0.5 * tau * tau + self->alpha_b;
        const gdouble phi   = (fabs (tau) < _NC_HICOSMO_VEXP_PHIA) ? _nc_hicosmo_Vexp_qt_phi_tau (Vexp, tau) : ncm_spline_eval (self->phi_tau, tau) * tau;
        gdouble H_lp, xq;

        _nc_hicosmo_Vexp_H_x (Vexp, 0.0, alpha, phi, &H_lp, &xq);

        {
          const gdouble dlnx_dalpha = _nc_hicosmo_Vexp_qt_dlnx_dalpha (Vexp, 0.0, alpha, phi);

          nu     = exp (-lna_a0e) * k * tau / (H_lp * self->RH_lp);
          dlnmnu = 2.0 * (dlnx_dalpha + 1.0) * tau;
        }

        break;
      }
      case 2: /* Classical expansion */
      {
        const gdouble lnq         = ncm_spline_eval (self->lnqe_tau, tau);
        const gdouble q           = exp (lnq);
        const gdouble dlnx_dalpha = _nc_hicosmo_Vexp_cl_dlnx_dalpha (Vexp, q, self->cl_be);
        const gdouble E           = _nc_hicosmo_Vexp_He_H0_q (Vexp, q, self->cl_be);

        nu     = exp (-lna_a0e) * k * tau / E;
        dlnmnu = 2.0 * (dlnx_dalpha + 1.0) * tau;

        break;
      }
      default:
        g_assert_not_reached ();
        break;
    }

    return dlnmnu / (2.0 * nu);
  }
}

static gdouble
_nc_hicosmo_Vexp_adiab_eval_nu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmoVexp *Vexp               = NC_HICOSMO_VEXP (iad);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);

  {
    const gdouble tau2    = tau * tau;
    const gdouble lna_a0e = 0.5 * tau2 + self->alpha_b - self->alpha_0e;
    gdouble nu;
    guint branch = 0;

    if (tau > self->tau_qt_c)
    {
      branch++;

      if (tau > self->tau_qt_e)
        branch++;
    }

    switch (branch)
    {
      case 0: /* Classical contraction */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqc_mtau, -tau);
        const gdouble q   = exp (lnq);
        const gdouble E   = _nc_hicosmo_Vexp_Hc_H0_q (Vexp, q, self->cl_bc);

        nu = exp (-lna_a0e) * k * tau / E;

        break;
      }
      case 1: /* Quantum phase */
      {
        const gdouble alpha = 0.5 * tau * tau + self->alpha_b;
        const gdouble phi   = (fabs (tau) < _NC_HICOSMO_VEXP_PHIA) ? _nc_hicosmo_Vexp_qt_phi_tau (Vexp, tau) : ncm_spline_eval (self->phi_tau, tau) * tau;
        gdouble H_lp, xq;

        _nc_hicosmo_Vexp_H_x (Vexp, 0.0, alpha, phi, &H_lp, &xq);

        nu = exp (-lna_a0e) * k * tau / (H_lp * self->RH_lp);

        break;
      }
      case 2: /* Classical expansion */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqe_tau, tau);
        const gdouble q   = exp (lnq);
        const gdouble E   = _nc_hicosmo_Vexp_He_H0_q (Vexp, q, self->cl_be);

        nu = exp (-lna_a0e) * k * tau / E;

        break;
      }
      default:
        g_assert_not_reached ();
        break;
    }

    return nu;
  }
}

static gdouble
_nc_hicosmo_Vexp_adiab_eval_m (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmoVexp *Vexp               = NC_HICOSMO_VEXP (iad);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);

  {
    const gdouble tau2    = tau * tau;
    const gdouble lna_a0e = 0.5 * tau2 + self->alpha_b - self->alpha_0e;

    gdouble m;
    guint branch = 0;

    if (tau > self->tau_qt_c)
    {
      branch++;

      if (tau > self->tau_qt_e)
        branch++;
    }

    switch (branch)
    {
      case 0: /* Classical contraction */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqc_mtau, -tau);
        const gdouble q   = exp (lnq);
        const gdouble x   = _nc_hicosmo_Vexp_x_q (q, self->cl_bc);
        const gdouble E   = _nc_hicosmo_Vexp_Hc_H0_q (Vexp, q, self->cl_bc);

        if (fabs (x) < 1.0e-1)
        {
          const gdouble dalpha = 0.5 * (tau - self->tau_x0) * (tau + self->tau_x0);
          const gdouble x_s    = _nc_hicosmo_Vexp_init_x_alpha_series (dalpha);

          m = x_s * x_s * exp (3.0 * lna_a0e) * E / tau;
        }
        else
        {
          m = x * x * exp (3.0 * lna_a0e) * E / tau;
        }

        break;
      }
      case 1: /* Quantum phase */
      {
        const gdouble alpha = 0.5 * tau * tau + self->alpha_b;
        const gdouble phi   = (fabs (tau) < _NC_HICOSMO_VEXP_PHIA) ? _nc_hicosmo_Vexp_qt_phi_tau (Vexp, tau) : ncm_spline_eval (self->phi_tau, tau) * tau;
        gdouble H_lp, xq;

        _nc_hicosmo_Vexp_H_x (Vexp, 0.0, alpha, phi, &H_lp, &xq);

        m = xq * xq * exp (3.0 * lna_a0e) * (H_lp * self->RH_lp) / tau;

        break;
      }
      case 2: /* Classical expansion */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqe_tau, tau);
        const gdouble q   = exp (lnq);
        const gdouble x   = _nc_hicosmo_Vexp_x_q (q, self->cl_be);
        const gdouble E   = _nc_hicosmo_Vexp_He_H0_q (Vexp, q, self->cl_be);

        if (fabs (x) < 1.0e-1)
        {
          const gdouble dalpha = 0.5 * (tau - self->tau_x0) * (tau + self->tau_x0);
          const gdouble x_s    = _nc_hicosmo_Vexp_init_x_alpha_series (dalpha);

          m = x_s * x_s * exp (3.0 * lna_a0e) * E / tau;
        }
        else
        {
          m = x * x * exp (3.0 * lna_a0e) * E / tau;
        }

        break;
      }
      default:
        g_assert_not_reached ();
        break;
    }

    return m;
  }
}

static gdouble
_nc_hicosmo_Vexp_adiab_eval_unit (NcHIPertIAdiab *iad)
{
  NcHICosmo *cosmo     = NC_HICOSMO (iad);
  const gdouble RH_lp  = nc_hicosmo_RH_planck (cosmo);
  const gdouble factor = sqrt (8.0 * ncm_c_pi () / 6.0);


  return factor / RH_lp;
}

static gdouble
_nc_hicosmo_Vexp_adiab_eval_x (NcHIPertIAdiab *iad, const gdouble tau)
{
  return nc_hicosmo_Vexp_xe_tau (NC_HICOSMO_VEXP (iad), tau);
}

static gdouble
_nc_hicosmo_Vexp_adiab_eval_p2Psi (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmoVexp *Vexp = NC_HICOSMO_VEXP (iad);
  const gdouble E     = nc_hicosmo_Vexp_E_tau (Vexp, tau);
  const gdouble x     = _nc_hicosmo_Vexp_adiab_eval_x (iad, tau);

  g_assert_not_reached ();
  /* Needs a double check */

  return -1.0 / (6.0 * E * x * k * k);
}

static gdouble
_nc_hicosmo_Vexp_adiab_eval_p2drho (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  g_assert_not_reached ();
}

/*
 * Interface implementation -- NcHIPertIGW
 */

static gdouble
_nc_hicosmo_Vexp_gw_eval_xi (NcHIPertIGW *igw, const gdouble tau, const gdouble k)
{
  NcHICosmoVexp *Vexp               = NC_HICOSMO_VEXP (igw);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);

  {
    const gdouble tau2    = tau * tau;
    const gdouble lna_a0e = 0.5 * tau2 + self->alpha_b - self->alpha_0e;

    return 2.0 * lna_a0e + log (k);
  }
}

static gdouble
_nc_hicosmo_Vexp_gw_eval_F1 (NcHIPertIGW *igw, const gdouble tau, const gdouble k)
{
  NcHICosmoVexp *Vexp               = NC_HICOSMO_VEXP (igw);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);

  {
    const gdouble tau2    = tau * tau;
    const gdouble lna_a0e = 0.5 * tau2 + self->alpha_b - self->alpha_0e;
    gdouble nu;
    guint branch = 0;

    if (tau > self->tau_qt_c)
    {
      branch++;

      if (tau > self->tau_qt_e)
        branch++;
    }

    switch (branch)
    {
      case 0: /* Classical contraction */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqc_mtau, -tau);
        const gdouble q   = exp (lnq);
        const gdouble E   = _nc_hicosmo_Vexp_Hc_H0_q (Vexp, q, self->cl_bc);

        nu = exp (-lna_a0e) * k * tau / E;

        break;
      }
      case 1: /* Quantum phase */
      {
        const gdouble alpha = 0.5 * tau * tau + self->alpha_b;
        const gdouble phi   = (fabs (tau) < _NC_HICOSMO_VEXP_PHIA) ? _nc_hicosmo_Vexp_qt_phi_tau (Vexp, tau) : (ncm_spline_eval (self->phi_tau, tau) * tau);
        gdouble H_lp, xq;

        _nc_hicosmo_Vexp_H_x (Vexp, 0.0, alpha, phi, &H_lp, &xq);

        nu = exp (-lna_a0e) * k * tau / (H_lp * self->RH_lp);

        break;
      }
      case 2: /* Classical expansion */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqe_tau, tau);
        const gdouble q   = exp (lnq);
        const gdouble E   = _nc_hicosmo_Vexp_He_H0_q (Vexp, q, self->cl_be);

        nu = exp (-lna_a0e) * k * tau / E;

        break;
      }
      default:
        g_assert_not_reached ();
        break;
    }

    return tau / nu; /* dxi = 2 * tau, F1 = dxi / (2nu). */
  }
}

static gdouble
_nc_hicosmo_Vexp_gw_eval_nu (NcHIPertIGW *igw, const gdouble tau, const gdouble k)
{
  NcHICosmoVexp *Vexp               = NC_HICOSMO_VEXP (igw);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);

  {
    const gdouble tau2    = tau * tau;
    const gdouble lna_a0e = 0.5 * tau2 + self->alpha_b - self->alpha_0e;
    gdouble nu;
    guint branch = 0;

    if (tau > self->tau_qt_c)
    {
      branch++;

      if (tau > self->tau_qt_e)
        branch++;
    }

    switch (branch)
    {
      case 0: /* Classical contraction */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqc_mtau, -tau);
        const gdouble q   = exp (lnq);
        const gdouble E   = _nc_hicosmo_Vexp_Hc_H0_q (Vexp, q, self->cl_bc);

        nu = exp (-lna_a0e) * k * tau / E;

        break;
      }
      case 1: /* Quantum phase */
      {
        const gdouble alpha = 0.5 * tau * tau + self->alpha_b;
        const gdouble phi   = (fabs (tau) < _NC_HICOSMO_VEXP_PHIA) ? _nc_hicosmo_Vexp_qt_phi_tau (Vexp, tau) : (ncm_spline_eval (self->phi_tau, tau) * tau);
        gdouble H_lp, xq;

        _nc_hicosmo_Vexp_H_x (Vexp, 0.0, alpha, phi, &H_lp, &xq);

        nu = exp (-lna_a0e) * k * tau / (H_lp * self->RH_lp);

        break;
      }
      case 2: /* Classical expansion */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqe_tau, tau);
        const gdouble q   = exp (lnq);
        const gdouble E   = _nc_hicosmo_Vexp_He_H0_q (Vexp, q, self->cl_be);

        nu = exp (-lna_a0e) * k * tau / E;

        break;
      }
      default:
        g_assert_not_reached ();
        break;
    }

    return nu;
  }
}

static gdouble
_nc_hicosmo_Vexp_gw_eval_m (NcHIPertIGW *igw, const gdouble tau, const gdouble k)
{
  NcHICosmoVexp *Vexp               = NC_HICOSMO_VEXP (igw);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);

  {
    const gdouble tau2    = tau * tau;
    const gdouble lna_a0e = 0.5 * tau2 + self->alpha_b - self->alpha_0e;
    guint branch          = 0;
    gdouble m;

    if (tau > self->tau_qt_c)
    {
      branch++;

      if (tau > self->tau_qt_e)
        branch++;
    }

    switch (branch)
    {
      case 0: /* Classical contraction */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqc_mtau, -tau);
        const gdouble q   = exp (lnq);
        const gdouble E   = _nc_hicosmo_Vexp_Hc_H0_q (Vexp, q, self->cl_bc);

        m = exp (3.0 * lna_a0e) * E / tau;

        break;
      }
      case 1: /* Quantum phase */
      {
        const gdouble alpha = 0.5 * tau * tau + self->alpha_b;
        const gdouble phi   = (fabs (tau) < _NC_HICOSMO_VEXP_PHIA) ? _nc_hicosmo_Vexp_qt_phi_tau (Vexp, tau) : (ncm_spline_eval (self->phi_tau, tau) * tau);
        gdouble H_lp, xq;

        _nc_hicosmo_Vexp_H_x (Vexp, 0.0, alpha, phi, &H_lp, &xq);

        m = exp (3.0 * lna_a0e) * (H_lp * self->RH_lp) / tau;

        break;
      }
      case 2: /* Classical expansion */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqe_tau, tau);
        const gdouble q   = exp (lnq);
        const gdouble E   = _nc_hicosmo_Vexp_He_H0_q (Vexp, q, self->cl_be);

        m = exp (3.0 * lna_a0e) * E / tau;

        break;
      }
      default:
        g_assert_not_reached ();
        break;
    }

    return m;
  }
}

static gdouble
_nc_hicosmo_Vexp_gw_eval_unit (NcHIPertIGW *igw)
{
  NcHICosmo *cosmo     = NC_HICOSMO (igw);
  const gdouble RH_lp  = nc_hicosmo_RH_planck (cosmo);
  const gdouble factor = sqrt (32.0 * ncm_c_pi ());

  return factor / RH_lp;
}

static gdouble
_nc_hicosmo_Vexp_gw_eval_x (NcHIPertIGW *igw, const gdouble tau)
{
  return nc_hicosmo_Vexp_xe_tau (NC_HICOSMO_VEXP (igw), tau);
}

/*
 * Interface implementation -- NcHIPertIGW
 */

static gdouble
_nc_hicosmo_Vexp_em_eval_xi (NcHIPertIEM *iem, const gdouble tau, const gdouble k)
{
  NcHICosmoVexp *Vexp               = NC_HICOSMO_VEXP (iem);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);

  return self->em_eval_xi (Vexp, tau, k);
}

static gdouble
_nc_hicosmo_Vexp_em_eval_F1 (NcHIPertIEM *iem, const gdouble tau, const gdouble k)
{
  NcHICosmoVexp *Vexp               = NC_HICOSMO_VEXP (iem);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);

  return self->em_eval_F1 (Vexp, tau, k);
}

static gdouble
_nc_hicosmo_Vexp_em_eval_nu (NcHIPertIEM *iem, const gdouble tau, const gdouble k)
{
  NcHICosmoVexp *Vexp               = NC_HICOSMO_VEXP (iem);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);

  return self->em_eval_nu (Vexp, tau, k);
}

static gdouble
_nc_hicosmo_Vexp_em_eval_m (NcHIPertIEM *iem, const gdouble tau, const gdouble k)
{
  NcHICosmoVexp *Vexp               = NC_HICOSMO_VEXP (iem);
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);

  return self->em_eval_m (Vexp, tau, k);
}

static gdouble
_nc_hicosmo_Vexp_em_eval_unit (NcHIPertIEM *iem)
{
  NcHICosmo *cosmo         = NC_HICOSMO (iem);
  const gdouble hbar_c_mu0 = ncm_c_hbar () * ncm_c_c () * ncm_c_magnetic_constant ();
  const gdouble RH         = nc_hicosmo_RH_Mpc (cosmo) * ncm_c_Mpc ();

  /* Factor in Gauss (1.0e-4 Tesla) */

  return sqrt (hbar_c_mu0) * 1.0e4 / gsl_pow_2 (RH);
}

static gdouble
_nc_hicosmo_Vexp_em_eval_x (NcHIPertIEM *iem, const gdouble tau)
{
  return nc_hicosmo_Vexp_xe_tau (NC_HICOSMO_VEXP (iem), tau);
}

static gdouble _nc_hicosmo_Vexp_em_eval_m_none (NcHICosmoVexp *Vexp, const gdouble tau, const gdouble k);

static gdouble
_nc_hicosmo_Vexp_em_eval_lapse (NcHIPertIEM *iem, const gdouble tau)
{
  return 1.0 / (4.0 * _nc_hicosmo_Vexp_em_eval_m_none (NC_HICOSMO_VEXP (iem), tau, 0.0));
}

/*
 * EM implementation that is independent of the specific form of the
 * potential. This is the default implementation.
 */

static gdouble
_nc_hicosmo_Vexp_em_eval_nu_all (NcHICosmoVexp *Vexp, const gdouble tau, const gdouble k)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);

  {
    const gdouble tau2    = tau * tau;
    const gdouble lna_a0e = 0.5 * tau2 + self->alpha_b - self->alpha_0e;
    gdouble nu;
    guint branch = 0;

    if (tau > self->tau_qt_c)
    {
      branch++;

      if (tau > self->tau_qt_e)
        branch++;
    }

    switch (branch)
    {
      case 0: /* Classical contraction */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqc_mtau, -tau);
        const gdouble q   = exp (lnq);
        const gdouble E   = _nc_hicosmo_Vexp_Hc_H0_q (Vexp, q, self->cl_bc);

        nu = exp (-lna_a0e) * k * tau / E;

        break;
      }
      case 1: /* Quantum phase */
      {
        const gdouble alpha = 0.5 * tau * tau + self->alpha_b;
        const gdouble phi   = (fabs (tau) < _NC_HICOSMO_VEXP_PHIA) ? _nc_hicosmo_Vexp_qt_phi_tau (Vexp, tau) : (ncm_spline_eval (self->phi_tau, tau) * tau);
        gdouble H_lp, xq;

        _nc_hicosmo_Vexp_H_x (Vexp, 0.0, alpha, phi, &H_lp, &xq);

        nu = exp (-lna_a0e) * k * tau / (H_lp * self->RH_lp);

        break;
      }
      case 2: /* Classical expansion */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqe_tau, tau);
        const gdouble q   = exp (lnq);
        const gdouble E   = _nc_hicosmo_Vexp_He_H0_q (Vexp, q, self->cl_be);

        nu = exp (-lna_a0e) * k * tau / E;

        break;
      }
      default:
        g_assert_not_reached ();
        break;
    }

    return nu;
  }
}

/*
 * EM implementation that is independent of the specific form of the
 * potential.
 */

static gdouble
_nc_hicosmo_Vexp_em_eval_xi_none (NcHICosmoVexp *Vexp, const gdouble tau, const gdouble k)
{
  return log (0.25 * k);
}

static gdouble
_nc_hicosmo_Vexp_em_eval_F1_none (NcHICosmoVexp *Vexp, const gdouble tau, const gdouble k)
{
  return 0.0;
}

static gdouble
_nc_hicosmo_Vexp_em_eval_m_none (NcHICosmoVexp *Vexp, const gdouble tau, const gdouble k)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);
  const gdouble F                   = 0.25;

  {
    const gdouble tau2    = tau * tau;
    const gdouble lna_a0e = 0.5 * tau2 + self->alpha_b - self->alpha_0e;
    gdouble N;
    guint branch = 0;

    if (tau > self->tau_qt_c)
    {
      branch++;

      if (tau > self->tau_qt_e)
        branch++;
    }

    switch (branch)
    {
      case 0: /* Classical contraction */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqc_mtau, -tau);
        const gdouble q   = exp (lnq);
        const gdouble E   = _nc_hicosmo_Vexp_Hc_H0_q (Vexp, q, self->cl_bc);

        N = exp (-lna_a0e) * tau / E;

        break;
      }
      case 1: /* Quantum phase */
      {
        const gdouble alpha = 0.5 * tau * tau + self->alpha_b;
        const gdouble phi   = (fabs (tau) < _NC_HICOSMO_VEXP_PHIA) ? _nc_hicosmo_Vexp_qt_phi_tau (Vexp, tau) : (ncm_spline_eval (self->phi_tau, tau) * tau);
        gdouble H_lp, xq;

        _nc_hicosmo_Vexp_H_x (Vexp, 0.0, alpha, phi, &H_lp, &xq);

        N = exp (-lna_a0e) * tau / (H_lp * self->RH_lp);

        break;
      }
      case 2: /* Classical expansion */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqe_tau, tau);
        const gdouble q   = exp (lnq);
        const gdouble E   = _nc_hicosmo_Vexp_He_H0_q (Vexp, q, self->cl_be);

        N = exp (-lna_a0e) * tau / E;

        break;
      }
      default:
        g_assert_not_reached ();
        break;
    }

    return F / N;
  }
}

static gdouble
_nc_hicosmo_Vexp_em_eval_xi_cauchy (NcHICosmoVexp *Vexp, const gdouble tau, const gdouble k)
{
  NcHICosmo *cosmo   = NC_HICOSMO (Vexp);
  const gdouble beta = EM_BETA;
  const gdouble B    = exp (gsl_pow_2 (EM_ALPHA));
  const gdouble phi  = nc_hicosmo_Vexp_phi (Vexp, tau);

  return log (0.25 * k)  + log1p (4.0 * B / (1.0 + gsl_pow_2 (phi / beta)));
}

static gdouble
_nc_hicosmo_Vexp_em_eval_F1_cauchy (NcHICosmoVexp *Vexp, const gdouble tau, const gdouble k)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);
  NcHICosmo *cosmo                  = NC_HICOSMO (Vexp);
  const gdouble beta                = EM_BETA;
  const gdouble B                   = exp (gsl_pow_2 (EM_ALPHA));
  const gdouble phi                 = nc_hicosmo_Vexp_phi (Vexp, tau);
  const gdouble t1                  = -8.0 * phi * B / (beta * beta) * 1.0 / gsl_pow_2 (1.0 + gsl_pow_2 (phi / beta)) * 1.0 / (1.0 + 4.0 * B / (1.0 + gsl_pow_2 (phi / beta)));
  gdouble x;

  {
    guint branch = 0;

    if (tau > self->tau_qt_c)
    {
      branch++;

      if (tau > self->tau_qt_e)
        branch++;
    }

    switch (branch)
    {
      case 0: /* Classical contraction */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqc_mtau, -tau);
        const gdouble q   = exp (lnq);
        const gdouble xl  = _nc_hicosmo_Vexp_x_q (q, self->cl_bc);

        if (fabs (xl) < 1.0e-1)
        {
          const gdouble dalpha = 0.5 * (tau - self->tau_x0) * (tau + self->tau_x0);
          const gdouble x_s    = _nc_hicosmo_Vexp_init_x_alpha_series (dalpha);

          x = x_s;
        }
        else
        {
          x = xl;
        }

        break;
      }
      case 1: /* Quantum phase */
      {
        const gdouble alpha = 0.5 * tau * tau + self->alpha_b;
        const gdouble phi   = (fabs (tau) < _NC_HICOSMO_VEXP_PHIA) ? _nc_hicosmo_Vexp_qt_phi_tau (Vexp, tau) : ncm_spline_eval (self->phi_tau, tau) * tau;
        gdouble H_lp, xq;

        _nc_hicosmo_Vexp_H_x (Vexp, 0.0, alpha, phi, &H_lp, &xq);

        x = xq;

        break;
      }
      case 2: /* Classical expansion */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqe_tau, tau);
        const gdouble q   = exp (lnq);
        const gdouble xl  = _nc_hicosmo_Vexp_x_q (q, self->cl_be);

        if (fabs (xl) < 1.0e-1)
        {
          const gdouble dalpha = 0.5 * (tau - self->tau_x0) * (tau + self->tau_x0);
          const gdouble x_s    = _nc_hicosmo_Vexp_init_x_alpha_series (dalpha);

          x = x_s;
        }
        else
        {
          x = xl;
        }

        break;
      }
      default:
        g_assert_not_reached ();
        break;
    }
  }

  {
    const gdouble nu = _nc_hicosmo_Vexp_em_eval_nu_all (Vexp, tau, k);

    return t1 * tau * x / (2.0 * nu);
  }
}

static gdouble
_nc_hicosmo_Vexp_em_eval_m_cauchy (NcHICosmoVexp *Vexp, const gdouble tau, const gdouble k)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);
  NcHICosmo *cosmo                  = NC_HICOSMO (Vexp);
  const gdouble beta                = EM_BETA;
  const gdouble B                   = exp (gsl_pow_2 (EM_ALPHA));
  const gdouble phi                 = nc_hicosmo_Vexp_phi (Vexp, tau);
  const gdouble F                   = 0.25 + B / (1.0 + gsl_pow_2 (phi / beta));

  {
    const gdouble tau2    = tau * tau;
    const gdouble lna_a0e = 0.5 * tau2 + self->alpha_b - self->alpha_0e;
    gdouble N;
    guint branch = 0;

    if (tau > self->tau_qt_c)
    {
      branch++;

      if (tau > self->tau_qt_e)
        branch++;
    }

    switch (branch)
    {
      case 0: /* Classical contraction */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqc_mtau, -tau);
        const gdouble q   = exp (lnq);
        const gdouble E   = _nc_hicosmo_Vexp_Hc_H0_q (Vexp, q, self->cl_bc);

        N = exp (-lna_a0e) * tau / E;

        break;
      }
      case 1: /* Quantum phase */
      {
        const gdouble alpha = 0.5 * tau * tau + self->alpha_b;
        const gdouble phi   = (fabs (tau) < _NC_HICOSMO_VEXP_PHIA) ? _nc_hicosmo_Vexp_qt_phi_tau (Vexp, tau) : (ncm_spline_eval (self->phi_tau, tau) * tau);
        gdouble H_lp, xq;

        _nc_hicosmo_Vexp_H_x (Vexp, 0.0, alpha, phi, &H_lp, &xq);

        N = exp (-lna_a0e) * tau / (H_lp * self->RH_lp);

        break;
      }
      case 2: /* Classical expansion */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqe_tau, tau);
        const gdouble q   = exp (lnq);
        const gdouble E   = _nc_hicosmo_Vexp_He_H0_q (Vexp, q, self->cl_be);

        N = exp (-lna_a0e) * tau / E;

        break;
      }
      default:
        g_assert_not_reached ();
        break;
    }

    return F / N;
  }
}

static gdouble
_nc_hicosmo_Vexp_em_eval_xi_gauss (NcHICosmoVexp *Vexp, const gdouble tau, const gdouble k)
{
  NcHICosmo *cosmo       = NC_HICOSMO (Vexp);
  const gdouble beta     = EM_BETA;
  const gdouble em_alpha = EM_ALPHA;
  const gdouble phi      = nc_hicosmo_Vexp_phi (Vexp, tau);
  const gdouble exponent = gsl_pow_2 (em_alpha) - gsl_pow_2 (phi / beta);

  return log (0.25 * k)  + log1p (4.0 * exp (exponent));
}

static gdouble
_nc_hicosmo_Vexp_em_eval_F1_gauss (NcHICosmoVexp *Vexp, const gdouble tau, const gdouble k)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);
  NcHICosmo *cosmo                  = NC_HICOSMO (Vexp);
  const gdouble beta                = EM_BETA;
  const gdouble em_alpha            = EM_ALPHA;
  const gdouble phi                 = nc_hicosmo_Vexp_phi (Vexp, tau);
  const gdouble exponent            = gsl_pow_2 (em_alpha) - gsl_pow_2 (phi / beta);
  const gdouble t1                  = -8.0 * phi / (beta * beta) * exp (exponent) / (1.0 + 4.0 * exp (exponent));
  gdouble x;

  {
    guint branch = 0;

    if (tau > self->tau_qt_c)
    {
      branch++;

      if (tau > self->tau_qt_e)
        branch++;
    }

    switch (branch)
    {
      case 0: /* Classical contraction */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqc_mtau, -tau);
        const gdouble q   = exp (lnq);
        const gdouble xl  = _nc_hicosmo_Vexp_x_q (q, self->cl_bc);

        if (fabs (xl) < 1.0e-1)
        {
          const gdouble dalpha = 0.5 * (tau - self->tau_x0) * (tau + self->tau_x0);
          const gdouble x_s    = _nc_hicosmo_Vexp_init_x_alpha_series (dalpha);

          x = x_s;
        }
        else
        {
          x = xl;
        }

        break;
      }
      case 1: /* Quantum phase */
      {
        const gdouble alpha = 0.5 * tau * tau + self->alpha_b;
        const gdouble phi   = (fabs (tau) < _NC_HICOSMO_VEXP_PHIA) ? _nc_hicosmo_Vexp_qt_phi_tau (Vexp, tau) : ncm_spline_eval (self->phi_tau, tau) * tau;
        gdouble H_lp, xq;

        _nc_hicosmo_Vexp_H_x (Vexp, 0.0, alpha, phi, &H_lp, &xq);

        x = xq;

        break;
      }
      case 2: /* Classical expansion */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqe_tau, tau);
        const gdouble q   = exp (lnq);
        const gdouble xl  = _nc_hicosmo_Vexp_x_q (q, self->cl_be);

        if (fabs (xl) < 1.0e-1)
        {
          const gdouble dalpha = 0.5 * (tau - self->tau_x0) * (tau + self->tau_x0);
          const gdouble x_s    = _nc_hicosmo_Vexp_init_x_alpha_series (dalpha);

          x = x_s;
        }
        else
        {
          x = xl;
        }

        break;
      }
      default:
        g_assert_not_reached ();
        break;
    }
  }

  {
    const gdouble nu = _nc_hicosmo_Vexp_em_eval_nu_all (Vexp, tau, k);

    return t1 * tau * x / (2.0 * nu);
  }
}

static gdouble
_nc_hicosmo_Vexp_em_eval_m_gauss (NcHICosmoVexp *Vexp, const gdouble tau, const gdouble k)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);
  NcHICosmo *cosmo                  = NC_HICOSMO (Vexp);
  const gdouble beta                = EM_BETA;
  const gdouble em_alpha            = EM_ALPHA;
  const gdouble phi                 = nc_hicosmo_Vexp_phi (Vexp, tau);
  const gdouble exponent            = gsl_pow_2 (em_alpha) - gsl_pow_2 (phi / beta);
  const gdouble F                   = 0.25 + exp (exponent);

  {
    const gdouble tau2    = tau * tau;
    const gdouble lna_a0e = 0.5 * tau2 + self->alpha_b - self->alpha_0e;
    gdouble N;
    guint branch = 0;

    if (tau > self->tau_qt_c)
    {
      branch++;

      if (tau > self->tau_qt_e)
        branch++;
    }

    switch (branch)
    {
      case 0: /* Classical contraction */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqc_mtau, -tau);
        const gdouble q   = exp (lnq);
        const gdouble E   = _nc_hicosmo_Vexp_Hc_H0_q (Vexp, q, self->cl_bc);

        N = exp (-lna_a0e) * tau / E;

        break;
      }
      case 1: /* Quantum phase */
      {
        const gdouble alpha = 0.5 * tau * tau + self->alpha_b;
        const gdouble phi   = (fabs (tau) < _NC_HICOSMO_VEXP_PHIA) ? _nc_hicosmo_Vexp_qt_phi_tau (Vexp, tau) : (ncm_spline_eval (self->phi_tau, tau) * tau);
        gdouble H_lp, xq;

        _nc_hicosmo_Vexp_H_x (Vexp, 0.0, alpha, phi, &H_lp, &xq);

        N = exp (-lna_a0e) * tau / (H_lp * self->RH_lp);

        break;
      }
      case 2: /* Classical expansion */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqe_tau, tau);
        const gdouble q   = exp (lnq);
        const gdouble E   = _nc_hicosmo_Vexp_He_H0_q (Vexp, q, self->cl_be);

        N = exp (-lna_a0e) * tau / E;

        break;
      }
      default:
        g_assert_not_reached ();
        break;
    }

    return F / N;
  }
}

/**
 * nc_hicosmo_Vexp_new:
 *
 * This function instantiates a new object of type #NcHICosmoVexp.
 *
 * Returns: A new #NcHICosmoVexp
 */
NcHICosmoVexp *
nc_hicosmo_Vexp_new (void)
{
  NcHICosmoVexp *Vexp = g_object_new (NC_TYPE_HICOSMO_VEXP, NULL);

  return Vexp;
}

/**
 * nc_hicosmo_Vexp_set_em_coupling:
 * @Vexp: a #NcHICosmoVexp
 * @coupling: a #NcHICosmoVexpEMCoupling
 *
 * Set the electromagnetic coupling.
 */
void
nc_hicosmo_Vexp_set_em_coupling (NcHICosmoVexp *Vexp, const NcHICosmoVexpEMCoupling coupling)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  if (self->em_coupling == coupling)
    return;

  switch (coupling)
  {
    case NC_HICOSMO_VEXP_EM_COUPLING_NONE:
      self->em_eval_F1 = _nc_hicosmo_Vexp_em_eval_F1_none;
      self->em_eval_nu = _nc_hicosmo_Vexp_em_eval_nu_all;
      self->em_eval_m  = _nc_hicosmo_Vexp_em_eval_m_none;
      self->em_eval_xi = _nc_hicosmo_Vexp_em_eval_xi_none;
      break;
    case NC_HICOSMO_VEXP_EM_COUPLING_CAUCHY:
      self->em_eval_F1 = _nc_hicosmo_Vexp_em_eval_F1_cauchy;
      self->em_eval_nu = _nc_hicosmo_Vexp_em_eval_nu_all;
      self->em_eval_m  = _nc_hicosmo_Vexp_em_eval_m_cauchy;
      self->em_eval_xi = _nc_hicosmo_Vexp_em_eval_xi_cauchy;
      break;
    case NC_HICOSMO_VEXP_EM_COUPLING_GAUSS:
      self->em_eval_F1 = _nc_hicosmo_Vexp_em_eval_F1_gauss;
      self->em_eval_nu = _nc_hicosmo_Vexp_em_eval_nu_all;
      self->em_eval_m  = _nc_hicosmo_Vexp_em_eval_m_gauss;
      self->em_eval_xi = _nc_hicosmo_Vexp_em_eval_xi_gauss;
      break;
    default:
      g_error ("Unknown electromagnetic coupling");

      return;
  }

  self->em_coupling = coupling;
  ncm_model_state_mark_outdated (NCM_MODEL (Vexp));
}

/**
 * nc_hicosmo_Vexp_get_em_coupling:
 * @Vexp: a #NcHICosmoVexp
 *
 * Get the electromagnetic coupling.
 *
 * Returns: a #NcHICosmoVexpEMCoupling
 */
NcHICosmoVexpEMCoupling
nc_hicosmo_Vexp_get_em_coupling (NcHICosmoVexp *Vexp)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  return self->em_coupling;
}

/**
 * nc_hicosmo_Vexp_tau_min:
 * @Vexp: a #NcHICosmoVexp
 *
 * The minimum value of the time variable suitable to describe the bounce, $\tau_{min}$.
 *
 * Returns: $\tau_{min}$
 */
gdouble
nc_hicosmo_Vexp_tau_min (NcHICosmoVexp *Vexp)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);
  {
    NcmVector *lnqc_mtau_xv = ncm_spline_peek_xv (self->lnqc_mtau);

    return -ncm_vector_get (lnqc_mtau_xv, ncm_vector_len (lnqc_mtau_xv) - 1);
  }
}

/**
 * nc_hicosmo_Vexp_tau_max:
 * @Vexp: a #NcHICosmoVexp
 *
 * The maximum value of the time variable suitable to describe the bounce, $\tau_{max}$.
 *
 * Returns: $\tau_{max}$
 */
gdouble
nc_hicosmo_Vexp_tau_max (NcHICosmoVexp *Vexp)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);
  {
    NcmVector *lnqe_tau_xv = ncm_spline_peek_xv (self->lnqe_tau);

    return ncm_vector_get (lnqe_tau_xv, ncm_vector_len (lnqe_tau_xv) - 1);
  }
}

/**
 * nc_hicosmo_Vexp_tau_qt_c:
 * @Vexp: a #NcHICosmoVexp
 *
 * Value of the time $\tau$ when the quantum regime begins during the contraction
 * phase.
 *
 * Returns: $\tau_\mathrm{qc}$.
 */
gdouble
nc_hicosmo_Vexp_tau_qt_c (NcHICosmoVexp *Vexp)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);

  return self->tau_qt_c;
}

/**
 * nc_hicosmo_Vexp_tau_qt_e:
 * @Vexp: a #NcHICosmoVexp
 *
 * Value of the time $\tau$ when the quantum regime ends during the expanding phase.
 *
 * Returns: $\tau_\mathrm{qe}$.
 */
gdouble
nc_hicosmo_Vexp_tau_qt_e (NcHICosmoVexp *Vexp)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);

  return self->tau_qt_e;
}

/**
 * nc_hicosmo_Vexp_xbe:
 * @Vexp: a #NcHICosmoVexp
 *
 * Computes the value of $a_{0\mathrm{e}} / a_b$ where $a_{0\mathrm{e}}$ is the value
 * of the scale factor at the scale when the Hubble parameter is equal to $\Omega_i
 * H_0$. Where $\Omega_i$ is the density parameter of the cold dark matter if $d_\phi >
 * 0$ or the density parameter of the cosmological constant if $d_\phi < 0$.
 *
 * Returns: $a_{0\mathrm{e}} / a_b$.
 */
gdouble
nc_hicosmo_Vexp_xbe (NcHICosmoVexp *Vexp)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);

  return exp (-self->alpha_b + self->alpha_0e);
}

/**
 * nc_hicosmo_Vexp_xbc:
 * @Vexp: a #NcHICosmoVexp
 *
 * Computes the value of $a_{0\mathrm{c}} / a_b$ where $a_{0\mathrm{c}}$ is the value
 * of the scale factor at the scale when the Hubble parameter is equal to $\Omega_i
 * H_0$. Where $\Omega_i$ is the density parameter of the cold dark matter if $d_\phi <
 * 0$ or the density parameter of the cosmological constant if $d_\phi > 0$.
 *
 * Returns: $a_{0\mathrm{c}} / a_b$.
 */
gdouble
nc_hicosmo_Vexp_xbc (NcHICosmoVexp *Vexp)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);

  return exp (-self->alpha_b + self->alpha_0c);
}

/**
 * nc_hicosmo_Vexp_alpha_0e:
 * @Vexp: a #NcHICosmoVexp
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hicosmo_Vexp_alpha_0e (NcHICosmoVexp *Vexp)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);

  return self->alpha_0e;
}

/**
 * nc_hicosmo_Vexp_alpha_0c:
 * @Vexp: a #NcHICosmoVexp
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hicosmo_Vexp_alpha_0c (NcHICosmoVexp *Vexp)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);

  return self->alpha_0c;
}

/**
 * nc_hicosmo_Vexp_xe_tau:
 * @Vexp: a #NcHICosmoVexp
 * @tau: $\tau$
 *
 * Computes the value of $x_e$ at time $\tau$. Note that all interface implementations of
 * perturbations use the observables computed at $x_\mathrm{e} = 1$, for example, the
 * physical wave number $k$ at $\tau$ is $k_\mathrm{phys} = k x_\mathrm{e}$.
 *
 * Returns: $x_\mathrm{e}$.
 */
gdouble
nc_hicosmo_Vexp_xe_tau (NcHICosmoVexp *Vexp, const gdouble tau)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);
  {
    const gdouble alpha = self->alpha_b + 0.5 * tau * tau;

    return exp (self->alpha_0e - alpha);
  }
}

/**
 * nc_hicosmo_Vexp_xc_tau:
 * @Vexp: a #NcHICosmoVexp
 * @tau: $\tau$
 *
 * Computes the value of $x_c$ at time $\tau$. Note that all interface implementations of
 * perturbations use the observables computed at $x_\mathrm{e} = 1$, see nc_hicosmo_Vexp_xe_tau().
 *
 *
 * Returns: $x_mathrm{c}$.
 */
gdouble
nc_hicosmo_Vexp_xc_tau (NcHICosmoVexp *Vexp, const gdouble tau)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);
  {
    const gdouble alpha = self->alpha_b + 0.5 * tau * tau;

    return exp (self->alpha_0c - alpha);
  }
}

/**
 * nc_hicosmo_Vexp_tau_xe:
 * @Vexp: a #NcHICosmoVexp
 * @xe: $x_e$
 *
 * FIXME
 *
 * Returns: $\tau$.
 */
gdouble
nc_hicosmo_Vexp_tau_xe (NcHICosmoVexp *Vexp, const gdouble xe)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);
  {
    const gdouble log_xe = log (xe);
    const gdouble alpha  = self->alpha_0e - log_xe;

    /*
     * We need to test for alpha < self->alpha_b, the fabs is to avoid domain
     * problems when alpha is approximately equal to self->alpha_b.
     */
    return +sqrt (2.0 * fabs (alpha - self->alpha_b));
  }
}

/**
 * nc_hicosmo_Vexp_tau_xc:
 * @Vexp: a #NcHICosmoVexp
 * @xc: $x_c$
 *
 * FIXME
 *
 * Returns: $\tau$.
 */
gdouble
nc_hicosmo_Vexp_tau_xc (NcHICosmoVexp *Vexp, const gdouble xc)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);
  {
    const gdouble log_xc = log (xc);
    const gdouble alpha  = self->alpha_0c - log_xc;

    /*
     * We need to test for alpha < self->alpha_b, the fabs is to avoid domain
     * problems when alpha is approximately equal to self->alpha_b.
     */
    return -sqrt (2.0 * fabs (alpha - self->alpha_b));
  }
}

/**
 * nc_hicosmo_Vexp_alpha:
 * @Vexp: a #NcHICosmoVexp
 * @tau: $\tau$
 *
 * Computes $\alpha = \ln a$, where $a$ is the scale factor, at time @tau.
 *
 * Returns: $\alpha(\tau)$.
 */
gdouble
nc_hicosmo_Vexp_alpha (NcHICosmoVexp *Vexp, const gdouble tau)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);

  return 0.5 * tau * tau + self->alpha_b;
}

/**
 * nc_hicosmo_Vexp_phi:
 * @Vexp: a #NcHICosmoVexp
 * @tau: $\tau$
 *
 * Computes the scalar field $\phi$ at time @tau.
 *
 * Returns: $\phi(\tau)$.
 */
gdouble
nc_hicosmo_Vexp_phi (NcHICosmoVexp *Vexp, const gdouble tau)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);
  {
    guint branch = 0;

    if (tau > self->tau_qt_c)
    {
      branch++;

      if (tau > self->tau_qt_e)
        branch++;
    }

    switch (branch)
    {
      case 0: /* Classical contraction */
      {
        const gdouble tau0      = self->tau_qt_c;
        const gdouble phi0      = ncm_spline_eval (self->phi_tau, tau0) * tau0;
        const gdouble lnq0      = ncm_spline_eval (self->lnqc_mtau, -tau0);
        const gdouble q0        = exp (lnq0);
        const gdouble E0        = _nc_hicosmo_Vexp_Hc_H0_q (Vexp, q0, self->cl_bc);
        const gdouble sqrt1mx02 = sqrt (_nc_hicosmo_Vexp_1mx_q (q0, self->cl_bc) * _nc_hicosmo_Vexp_1px_q (q0, self->cl_bc));
        const gdouble L0        = fabs (E0) * sqrt1mx02 * exp (3.0 * phi0 / M_SQRT2);

        const gdouble lnq      = ncm_spline_eval (self->lnqc_mtau, -tau);
        const gdouble q        = exp (lnq);
        const gdouble E        = _nc_hicosmo_Vexp_Hc_H0_q (Vexp, q, self->cl_bc);
        const gdouble sqrt1mx2 = sqrt (_nc_hicosmo_Vexp_1mx_q (q, self->cl_bc) * _nc_hicosmo_Vexp_1px_q (q, self->cl_bc));

        return M_SQRT2 * log (L0 / fabs (E * sqrt1mx2)) / 3.0;

        break;
      }
      case 1: /* Quantum phase */
      {
        const gdouble phi = (fabs (tau) < _NC_HICOSMO_VEXP_PHIA) ? _nc_hicosmo_Vexp_qt_phi_tau (Vexp, tau) : ncm_spline_eval (self->phi_tau, tau) * tau;

        return phi;

        break;
      }
      case 2: /* Classical expansion */
      {
        const gdouble tau0      = self->tau_qt_e;
        const gdouble phi0      = ncm_spline_eval (self->phi_tau, tau0) * tau0;
        const gdouble lnq0      = ncm_spline_eval (self->lnqe_tau, tau0);
        const gdouble q0        = exp (lnq0);
        const gdouble E0        = _nc_hicosmo_Vexp_He_H0_q (Vexp, q0, self->cl_be);
        const gdouble sqrt1mx02 = sqrt (_nc_hicosmo_Vexp_1mx_q (q0, self->cl_be) * _nc_hicosmo_Vexp_1px_q (q0, self->cl_be));
        const gdouble L0        = fabs (E0) * sqrt1mx02 * exp (3.0 * phi0 / M_SQRT2);

        const gdouble lnq      = ncm_spline_eval (self->lnqe_tau, tau);
        const gdouble q        = exp (lnq);
        const gdouble E        = _nc_hicosmo_Vexp_He_H0_q (Vexp, q, self->cl_be);
        const gdouble sqrt1mx2 = sqrt (_nc_hicosmo_Vexp_1mx_q (q, self->cl_be) * _nc_hicosmo_Vexp_1px_q (q, self->cl_be));

        return M_SQRT2 * log (L0 / fabs (E * sqrt1mx2)) / 3.0;

        break;
      }
      default:
        g_assert_not_reached ();
        break;
    }
  }
}

/**
 * nc_hicosmo_Vexp_Ricci_scale:
 * @Vexp: a #NcHICosmoVexp
 * @tau: $\tau$
 *
 * FIXME
 *
 * Returns: $L_R(\tau) / \ell_\mathrm{P}$.
 */
gdouble
nc_hicosmo_Vexp_Ricci_scale (NcHICosmoVexp *Vexp, const gdouble tau)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);
  {
    guint branch = 0;

    if (tau > self->tau_qt_c)
    {
      branch++;

      if (tau > self->tau_qt_e)
        branch++;
    }

    switch (branch)
    {
      case 0: /* Classical contraction */
      {
        const gdouble lnq      = ncm_spline_eval (self->lnqc_mtau, -tau);
        const gdouble q        = exp (lnq);
        const gdouble E        = _nc_hicosmo_Vexp_Hc_H0_q (Vexp, q, self->cl_bc);
        const gdouble two_3mx2 = _nc_hicosmo_Vexp_1mx_q (q, self->cl_bc) * _nc_hicosmo_Vexp_1px_q (q, self->cl_bc) - 1.0 / 3.0;

        return self->RH_lp / ((sqrt (18.0 * fabs (two_3mx2)) * fabs (E)));

        break;
      }
      case 1: /* Quantum phase */
      {
        const gdouble alpha = 0.5 * tau * tau + self->alpha_b;
        const gdouble phi   = (fabs (tau) < _NC_HICOSMO_VEXP_PHIA) ? _nc_hicosmo_Vexp_qt_phi_tau (Vexp, tau) : ncm_spline_eval (self->phi_tau, tau) * tau;

        return _nc_hicosmo_Vexp_qt_Ricci_scale (Vexp, NC_HICOSMO (Vexp), alpha, phi);

        break;
      }
      case 2: /* Classical expansion */
      {
        const gdouble lnq      = ncm_spline_eval (self->lnqe_tau, tau);
        const gdouble q        = exp (lnq);
        const gdouble E        = _nc_hicosmo_Vexp_He_H0_q (Vexp, q, self->cl_be);
        const gdouble two_3mx2 = _nc_hicosmo_Vexp_1mx_q (q, self->cl_be) * _nc_hicosmo_Vexp_1px_q (q, self->cl_be) - 1.0 / 3.0;

        return self->RH_lp / ((sqrt (18.0 * fabs (two_3mx2)) * fabs (E)));

        break;
      }
      default:
        g_assert_not_reached ();
        break;
    }
  }
}

/**
 * nc_hicosmo_Vexp_E_tau:
 * @Vexp: a #NcHICosmoVexp
 * @tau: $\tau$
 *
 * Computes the value of $E(\tau)$, where $E(\tau) = H(\tau) / H_0$.
 *
 * Returns: $E(\tau)$.
 */
gdouble
nc_hicosmo_Vexp_E_tau (NcHICosmoVexp *Vexp, const gdouble tau)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);
  {
    guint branch = 0;

    if (tau > self->tau_qt_c)
    {
      branch++;

      if (tau > self->tau_qt_e)
        branch++;
    }

    switch (branch)
    {
      case 0: /* Classical contraction */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqc_mtau, -tau);
        const gdouble q   = exp (lnq);
        const gdouble E   = _nc_hicosmo_Vexp_Hc_H0_q (Vexp, q, self->cl_bc);

        return E;

        break;
      }
      case 1: /* Quantum phase */
      {
        const gdouble alpha = 0.5 * tau * tau + self->alpha_b;
        const gdouble phi   = (fabs (tau) < _NC_HICOSMO_VEXP_PHIA) ? _nc_hicosmo_Vexp_qt_phi_tau (Vexp, tau) : (ncm_spline_eval (self->phi_tau, tau) * tau);
        gdouble H_lp, xq;

        _nc_hicosmo_Vexp_H_x (Vexp, 0.0, alpha, phi, &H_lp, &xq);

        return (H_lp * self->RH_lp);

        break;
      }
      case 2: /* Classical expansion */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqe_tau, tau);
        const gdouble q   = exp (lnq);
        const gdouble E   = _nc_hicosmo_Vexp_He_H0_q (Vexp, q, self->cl_be);

        return E;

        break;
      }
      default:
        g_assert_not_reached ();

        return 0.0;

        break;
    }
  }
}

/**
 * nc_hicosmo_Vexp_x_y:
 * @Vexp: a #NcHICosmoVexp
 * @tau: $\tau$
 * @x: (out): the value of $x(\tau)$
 * @y: (out): the value of $y(\tau)$
 *
 * FIXME
 *
 */
void
nc_hicosmo_Vexp_x_y (NcHICosmoVexp *Vexp, const gdouble tau, gdouble *x, gdouble *y)
{
  NcHICosmoVexpPrivate * const self = nc_hicosmo_Vexp_get_instance_private (Vexp);

  _nc_hicosmo_Vexp_prepare (Vexp);
  {
    guint branch = 0;

    if (tau > self->tau_qt_c)
    {
      branch++;

      if (tau > self->tau_qt_e)
        branch++;
    }

    switch (branch)
    {
      case 0: /* Classical contraction */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqc_mtau, -tau);
        const gdouble q   = exp (lnq);
        const gdouble xl  = _nc_hicosmo_Vexp_x_q (q, self->cl_bc);

        if (fabs (xl) < 1.0e-1)
        {
          const gdouble dalpha = 0.5 * (tau - self->tau_x0) * (tau + self->tau_x0);
          const gdouble x_s    = _nc_hicosmo_Vexp_init_x_alpha_series (dalpha);

          x[0] = x_s;
        }
        else
        {
          x[0] = xl;
        }

        break;
      }
      case 1: /* Quantum phase */
      {
        const gdouble alpha = 0.5 * tau * tau + self->alpha_b;
        const gdouble phi   = (fabs (tau) < _NC_HICOSMO_VEXP_PHIA) ? _nc_hicosmo_Vexp_qt_phi_tau (Vexp, tau) : ncm_spline_eval (self->phi_tau, tau) * tau;
        gdouble H_lp, xq;

        _nc_hicosmo_Vexp_H_x (Vexp, 0.0, alpha, phi, &H_lp, &xq);

        x[0] = xq;

        break;
      }
      case 2: /* Classical expansion */
      {
        const gdouble lnq = ncm_spline_eval (self->lnqe_tau, tau);
        const gdouble q   = exp (lnq);
        const gdouble xl  = _nc_hicosmo_Vexp_x_q (q, self->cl_be);

        if (fabs (xl) < 1.0e-1)
        {
          const gdouble dalpha = 0.5 * (tau - self->tau_x0) * (tau + self->tau_x0);
          const gdouble x_s    = _nc_hicosmo_Vexp_init_x_alpha_series (dalpha);

          x[0] = x_s;
        }
        else
        {
          x[0] = xl;
        }

        break;
      }
      default:
        g_assert_not_reached ();
        break;
    }

    y[0] = GSL_SIGN (tau) * sqrt (fabs (1.0 - x[0] * x[0]));

    return;
  }
}

