/***************************************************************************
 *            nc_hipert_boltzmann_std.c (linear_generic.c)
 *
 *  Thu Nov 12 12:37:02 2009
 *  Copyright  2009  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hipert_boltzmann_std.c
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

/**
 * NcHIPertBoltzmannStd:
 *
 * Perturbations object for standard Boltzmann hierarchy model.
 *
 * This object is a subclass of #NcHIPertBoltzmann and is designed to implement the
 * perturbations for the standard Boltzmann hierarchy model.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_hipert_boltzmann_std.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>

#define SUN_DENSE_ACCESS SM_ELEMENT_D
#define SUN_BAND_ACCESS SM_ELEMENT_D
#endif /* NUMCOSMO_GIR_SCAN */

#include "perturbations/nc_hipert_private.h"

enum
{
  PROP_0,
  PROP_LMAX,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcHIPertBoltzmannStd, nc_hipert_boltzmann_std, NC_TYPE_HIPERT_BOLTZMANN)

static void
nc_hipert_boltzmann_std_init (NcHIPertBoltzmannStd *pbs)
{
}

static void
nc_hipert_boltzmann_std_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_boltzmann_std_parent_class)->finalize (object);
}

static void
nc_hipert_boltzmann_std_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  /*NcHIPertBoltzmannStd *pbs = NC_HIPERT_BOLTZMANN_STD (object);*/
  g_return_if_fail (NC_IS_HIPERT_BOLTZMANN_STD (object));

  switch (prop_id)
  {
    case PROP_LMAX:
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hipert_boltzmann_std_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /*NcHIPertBoltzmannStd *pbs = NC_HIPERT_BOLTZMANN_STD (object);*/
  g_return_if_fail (NC_IS_HIPERT_BOLTZMANN_STD (object));

  switch (prop_id)
  {
    case PROP_LMAX:
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void _nc_hipert_boltzmann_std_init (NcHIPertBoltzmann *pb, NcHICosmo *cosmo);
static void _nc_hipert_boltzmann_std_set_opts (NcHIPertBoltzmann *pb);
static void _nc_hipert_boltzmann_std_reset (NcHIPertBoltzmann *pb);
static void _nc_hipert_boltzmann_std_evol_step (NcHIPertBoltzmann *pb, gdouble lambda);
static void _nc_hipert_boltzmann_std_evol (NcHIPertBoltzmann *pb, gdouble g);
static void _nc_hipert_boltzmann_std_get_sources (NcHIPertBoltzmann *pb, gdouble *S0, gdouble *S1, gdouble *S2);
static void _nc_hipert_boltzmann_std_print_stats (NcHIPertBoltzmann *pb);
static gdouble _nc_hipert_boltzmann_std_get_z (NcHIPertBoltzmann *pb);
static gdouble _nc_hipert_boltzmann_std_get_phi (NcHIPertBoltzmann *pb);
static gdouble _nc_hipert_boltzmann_std_get_c0 (NcHIPertBoltzmann *pb);
static gdouble _nc_hipert_boltzmann_std_get_b0 (NcHIPertBoltzmann *pb);
static gdouble _nc_hipert_boltzmann_std_get_c1 (NcHIPertBoltzmann *pb);
static gdouble _nc_hipert_boltzmann_std_get_b1 (NcHIPertBoltzmann *pb);
static gdouble _nc_hipert_boltzmann_std_get (NcHIPertBoltzmann *pb, guint n);
static gdouble _nc_hipert_boltzmann_std_get_theta (NcHIPertBoltzmann *pb, guint n);
static gdouble _nc_hipert_boltzmann_std_get_theta_p (NcHIPertBoltzmann *pb, guint n);
static gdouble _nc_hipert_boltzmann_std_get_los_theta (NcHIPertBoltzmann *pb, guint n);
static void _nc_hipert_boltzmann_std_print_all (NcHIPertBoltzmann *pb);

static void
nc_hipert_boltzmann_std_class_init (NcHIPertBoltzmannStdClass *klass)
{
  GObjectClass *object_class       = G_OBJECT_CLASS (klass);
  NcHIPertBoltzmannClass *pb_class = NC_HIPERT_BOLTZMANN_CLASS (klass);

  object_class->finalize     = nc_hipert_boltzmann_std_finalize;
  object_class->set_property = nc_hipert_boltzmann_std_set_property;
  object_class->get_property = nc_hipert_boltzmann_std_get_property;

  g_object_class_install_property (object_class,
                                   PROP_LMAX,
                                   g_param_spec_uint ("l-maxa",
                                                      NULL,
                                                      "Last multipole",
                                                      2, G_MAXUINT32, 2,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  pb_class->init          = &_nc_hipert_boltzmann_std_init;
  pb_class->set_opts      = &_nc_hipert_boltzmann_std_set_opts;
  pb_class->reset         = &_nc_hipert_boltzmann_std_reset;
  pb_class->evol_step     = &_nc_hipert_boltzmann_std_evol_step;
  pb_class->evol          = &_nc_hipert_boltzmann_std_evol;
  pb_class->get_sources   = &_nc_hipert_boltzmann_std_get_sources;
  pb_class->print_stats   = &_nc_hipert_boltzmann_std_print_stats;
  pb_class->get_z         = &_nc_hipert_boltzmann_std_get_z;
  pb_class->get_phi       = &_nc_hipert_boltzmann_std_get_phi;
  pb_class->get_c0        = &_nc_hipert_boltzmann_std_get_c0;
  pb_class->get_b0        = &_nc_hipert_boltzmann_std_get_b0;
  pb_class->get_c1        = &_nc_hipert_boltzmann_std_get_c1;
  pb_class->get_b1        = &_nc_hipert_boltzmann_std_get_b1;
  pb_class->get           = &_nc_hipert_boltzmann_std_get;
  pb_class->get_theta     = &_nc_hipert_boltzmann_std_get_theta;
  pb_class->get_theta_p   = &_nc_hipert_boltzmann_std_get_theta_p;
  pb_class->get_los_theta = &_nc_hipert_boltzmann_std_get_los_theta;
  pb_class->print_all     = &_nc_hipert_boltzmann_std_print_all;
}

#define _NC_PHI (NV_Ith_S (pert->priv->y, NC_HIPERT_BOLTZMANN_PHI))
#define _NC_C0 (NV_Ith_S (pert->priv->y, NC_HIPERT_BOLTZMANN_C0))

#define _NC_B0 (NV_Ith_S (pert->priv->y, NC_HIPERT_BOLTZMANN_B0))
#define _NC_dB0 (NV_Ith_S (pert->priv->y, NC_HIPERT_BOLTZMANN_dB0))

#define _NC_THETA0 (NV_Ith_S (pert->priv->y, NC_HIPERT_BOLTZMANN_THETA0))
#define _NC_dTHETA0 (NV_Ith_S (pert->priv->y, NC_HIPERT_BOLTZMANN_dTHETA0))

#define _NC_C1 (NV_Ith_S (pert->priv->y, NC_HIPERT_BOLTZMANN_C1))
#define _NC_V (NV_Ith_S (pert->priv->y, NC_HIPERT_BOLTZMANN_V))

#define _NC_U (NV_Ith_S (pert->priv->y, NC_HIPERT_BOLTZMANN_U))
#define _NC_THETA1 (NV_Ith_S (pert->priv->y, NC_HIPERT_BOLTZMANN_THETA1))

#define _NC_T (NV_Ith_S (pert->priv->y, NC_HIPERT_BOLTZMANN_T))
#define _NC_B1 (NV_Ith_S (pert->priv->y, NC_HIPERT_BOLTZMANN_B1))

#define _NC_THETA2 (NV_Ith_S (pert->priv->y, NC_HIPERT_BOLTZMANN_THETA2))
#define _NC_THETA(n) (NV_Ith_S (pert->priv->y, NC_HIPERT_BOLTZMANN_THETA (n)))

#define _NC_THETA_P0 (NV_Ith_S (pert->priv->y, NC_HIPERT_BOLTZMANN_THETA_P0))
#define _NC_THETA_P1 (NV_Ith_S (pert->priv->y, NC_HIPERT_BOLTZMANN_THETA_P1))
#define _NC_THETA_P2 (NV_Ith_S (pert->priv->y, NC_HIPERT_BOLTZMANN_THETA_P2))
#define _NC_THETA_P(n) (NV_Ith_S (pert->priv->y, NC_HIPERT_BOLTZMANN_THETA_P (n)))

#define _NC_DPHI (NV_Ith_S (ydot, NC_HIPERT_BOLTZMANN_PHI))
#define _NC_DC0 (NV_Ith_S (ydot, NC_HIPERT_BOLTZMANN_C0))

#define _NC_DB0 (NV_Ith_S (ydot, NC_HIPERT_BOLTZMANN_B0))
#define _NC_DdB0 (NV_Ith_S (ydot, NC_HIPERT_BOLTZMANN_dB0))

#define _NC_DTHETA0 (NV_Ith_S (ydot, NC_HIPERT_BOLTZMANN_THETA0))
#define _NC_DdTHETA0 (NV_Ith_S (ydot, NC_HIPERT_BOLTZMANN_dTHETA0))

#define _NC_DV (NV_Ith_S (ydot, NC_HIPERT_BOLTZMANN_V))
#define _NC_DC1 (NV_Ith_S (ydot, NC_HIPERT_BOLTZMANN_C1))

#define _NC_DU (NV_Ith_S (ydot, NC_HIPERT_BOLTZMANN_U))
#define _NC_DTHETA1 (NV_Ith_S (ydot, NC_HIPERT_BOLTZMANN_THETA1))

#define _NC_DT (NV_Ith_S (ydot, NC_HIPERT_BOLTZMANN_T))
#define _NC_DB1 (NV_Ith_S (ydot, NC_HIPERT_BOLTZMANN_B1))

#define _NC_DTHETA2 (NV_Ith_S (ydot, NC_HIPERT_BOLTZMANN_THETA2))
#define _NC_DTHETA(n) (NV_Ith_S (ydot, NC_HIPERT_BOLTZMANN_THETA (n)))

#define _NC_DTHETA_P0 (NV_Ith_S (ydot, NC_HIPERT_BOLTZMANN_THETA_P0))
#define _NC_DTHETA_P1 (NV_Ith_S (ydot, NC_HIPERT_BOLTZMANN_THETA_P1))
#define _NC_DTHETA_P2 (NV_Ith_S (ydot, NC_HIPERT_BOLTZMANN_THETA_P2))
#define _NC_DTHETA_P(n) (NV_Ith_S (ydot, NC_HIPERT_BOLTZMANN_THETA_P (n)))

static void
_nc_hipert_boltzmann_std_init (NcHIPertBoltzmann *pb, NcHICosmo *cosmo)
{
  NcHIPert *pert             = NC_HIPERT (pb);
  const gdouble x            = NC_HIPERT_BOLTZMANN_LAMBDA2X (pb->lambdai);
  const gdouble z            = x - 1.0;
  const gdouble E2           = nc_hicosmo_E2 (cosmo, z);
  const gdouble E            = sqrt (E2);
  const gdouble k_E          = nc_hipert_get_mode_k (pert) / E;
  const gdouble kx_E         = x * k_E;
  const gdouble tauprime     = nc_recomb_dtau_dlambda (pb->recomb, cosmo, pb->lambdai);
  const gdouble kx_Etauprime = kx_E / tauprime;
  gdouble theta1;
  guint i;

  nc_hicosmo_clear (&pb->cosmo);
  pb->cosmo = nc_hicosmo_ref (cosmo);

  N_VConst (0.0, pert->priv->y);

  _NC_PHI    = 1.0;
  _NC_C0     = 3.0 * _NC_PHI / 2.0;
  _NC_THETA0 = _NC_B0 = _NC_C0;

  _NC_THETA1 = _NC_B1 = _NC_C1 = theta1 = -kx_E / 6.0 * _NC_PHI;

  _NC_THETA2   = -8.0 / 15.0 * kx_Etauprime * theta1;
  _NC_THETA_P0 =  5.0 /  4.0 * _NC_THETA2;
  _NC_THETA_P2 =  1.0 /  4.0 * _NC_THETA2;

  _NC_THETA_P1 = -1.0 / 4.0 * kx_Etauprime * _NC_THETA2;

  for (i = 3; i <= pb->TT_lmax; i++)
  {
    const gdouble l  = i;
    const gdouble f1 = -l * kx_Etauprime / (2.0 * l + 1.0);

    _NC_THETA (i)   = f1 * _NC_THETA (i - 1);
    _NC_THETA_P (i) = f1 * _NC_THETA_P (i - 1);
  }

  return;
}

static gint _nc_hipert_boltzmann_std_step (sunrealtype lambda, N_Vector y, N_Vector ydot, gpointer user_data);
static gint _nc_hipert_boltzmann_std_band_J (sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static void
_nc_hipert_boltzmann_std_set_opts (NcHIPertBoltzmann *pb)
{
  NcHIPert *pert = NC_HIPERT (pb);
  gint flag, i;

  N_VConst (pert->priv->abstol, pert->priv->vec_abstol);

  for (i = 0; i <= NC_HIPERT_BOLTZMANN_THETA_P2; i++)
    NV_Ith_S (pert->priv->vec_abstol, i) = 0.0;

  flag = CVodeSVtolerances (pert->priv->cvode, pert->priv->reltol, pert->priv->vec_abstol);
  NCM_CVODE_CHECK (&flag, "CVodeSVtolerances", 1, );

  flag = CVodeSetMaxNumSteps (pert->priv->cvode, 1000000);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

  flag = CVodeSetUserData (pert->priv->cvode, pert);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

  flag = CVodeSetLinearSolver (pert->priv->cvode, pert->priv->LS, pert->priv->A);
  NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, );

  flag = CVodeSetJacFn (pert->priv->cvode, &_nc_hipert_boltzmann_std_band_J);
  NCM_CVODE_CHECK (&flag, "CVodeSetJacFn", 1, );

  flag = CVodeSetStopTime (pert->priv->cvode, pb->lambdaf);
  NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

  return;
}

static void
_nc_hipert_boltzmann_std_reset (NcHIPertBoltzmann *pb)
{
  gint flag;
  NcHIPert *pert = NC_HIPERT (pb);

  if (!pert->priv->cvode_init)
  {
    flag = CVodeInit (pert->priv->cvode, &_nc_hipert_boltzmann_std_step, pb->lambdai, pert->priv->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
    pert->priv->cvode_init = TRUE;
  }
  else
  {
    flag = CVodeReInit (pert->priv->cvode, pb->lambdai, pert->priv->y);
    NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );
  }

  _nc_hipert_boltzmann_std_set_opts (pb);

  return;
}

static void
_nc_hipert_boltzmann_std_evol_step (NcHIPertBoltzmann *pb, gdouble lambda)
{
  NcHIPert *pert = NC_HIPERT (pb);
  gint flag;
  gdouble lambdai;

  flag = CVode (pert->priv->cvode, lambda, pert->priv->y, &lambdai, CV_ONE_STEP);
  NCM_CVODE_CHECK (&flag, "CVode", 1, );
}

static void
_nc_hipert_boltzmann_std_evol (NcHIPertBoltzmann *pb, gdouble lambda)
{
  NcHIPert *pert = NC_HIPERT (pb);
  gint flag;
  gdouble lambdai = pb->lambdai;

  while (lambdai < lambda)
  {
    flag = CVode (pert->priv->cvode, lambda, pert->priv->y, &lambdai, CV_NORMAL);
    NCM_CVODE_CHECK (&flag, "CVode", 1, );
  }
}

static void
_nc_hipert_boltzmann_std_get_sources (NcHIPertBoltzmann *pb, gdouble *S0, gdouble *S1, gdouble *S2)
{
}

static void
_nc_hipert_boltzmann_std_print_stats (NcHIPertBoltzmann *pb)
{
}

static gdouble
_nc_hipert_boltzmann_std_get_z (NcHIPertBoltzmann *pb)
{
  return 0.0;
}

static gdouble
_nc_hipert_boltzmann_std_get_phi (NcHIPertBoltzmann *pb)
{
  return 0.0;
}

static gdouble
_nc_hipert_boltzmann_std_get_c0 (NcHIPertBoltzmann *pb)
{
  return 0.0;
}

static gdouble
_nc_hipert_boltzmann_std_get_b0 (NcHIPertBoltzmann *pb)
{
  return 0.0;
}

static gdouble
_nc_hipert_boltzmann_std_get_c1 (NcHIPertBoltzmann *pb)
{
  return 0.0;
}

static gdouble
_nc_hipert_boltzmann_std_get_b1 (NcHIPertBoltzmann *pb)
{
  return 0.0;
}

static gdouble
_nc_hipert_boltzmann_std_get (NcHIPertBoltzmann *pb, guint n)
{
  return 0.0;
}

static gdouble
_nc_hipert_boltzmann_std_get_theta (NcHIPertBoltzmann *pb, guint n)
{
  return 0.0;
}

static gdouble
_nc_hipert_boltzmann_std_get_theta_p (NcHIPertBoltzmann *pb, guint n)
{
  return 0.0;
}

static gdouble
_nc_hipert_boltzmann_std_get_los_theta (NcHIPertBoltzmann *pb, guint n)
{
  return 0.0;
}

static void
_nc_hipert_boltzmann_std_print_all (NcHIPertBoltzmann *pb)
{
}

/**
 * nc_hipert_boltzmann_std_new:
 * @recomb: a #NcRecomb.
 * @lmax: FIXME
 *
 * FIXME
 *
 * Returns: (transfer full): a new #NcHIPertBoltzmannStd object.
 */
NcHIPertBoltzmannStd *
nc_hipert_boltzmann_std_new (NcRecomb *recomb, guint lmax)
{
  NcHIPertBoltzmannStd *pbs = g_object_new (NC_TYPE_HIPERT_BOLTZMANN_STD,
                                            "recomb", recomb,
                                            "l-max", lmax,
                                            NULL);

  NC_HIPERT_BOLTZMANN (pbs)->lambdai = NC_HIPERT_BOLTZMANN_X2LAMBDA (1.0e9);
  NC_HIPERT_BOLTZMANN (pbs)->lambdaf = NC_HIPERT_BOLTZMANN_X2LAMBDA (1.0);
  nc_hipert_set_stiff_solver (NC_HIPERT (pbs), TRUE);

  return pbs;
}

static gint
_nc_hipert_boltzmann_std_step (sunrealtype lambda, N_Vector y, N_Vector ydot, gpointer user_data)
{
  NcHIPertBoltzmann *pb  = NC_HIPERT_BOLTZMANN (user_data);
  NcHIPert *pert         = NC_HIPERT (pb);
  NcHICosmo *cosmo       = pb->cosmo;
  const guint lmax       = pb->TT_lmax;
  const gdouble Omega_r0 = nc_hicosmo_Omega_r0 (cosmo);
  const gdouble Omega_b0 = nc_hicosmo_Omega_b0 (cosmo);
  const gdouble Omega_c0 = nc_hicosmo_Omega_c0 (cosmo);
  const gdouble Omega_m0 = nc_hicosmo_Omega_m0 (cosmo);
  const gdouble R0       = 4.0 * Omega_r0 / (3.0 * Omega_b0);
  const gdouble x        = NC_HIPERT_BOLTZMANN_LAMBDA2X (lambda);
  const gdouble R        = R0 * x;
  const gdouble x2       = x * x;
  const gdouble x3       = x2 * x;
  const gdouble k        = nc_hipert_get_mode_k (pert);
  const gdouble k2       = k * k;
  const gdouble E2       = nc_hicosmo_E2 (cosmo, x - 1.0);
  const gdouble E        = sqrt (E2);
  const gdouble k_E      = k / E;
  const gdouble kx_E     = x * k_E;
  const gdouble kx_3E    = kx_E / 3.0;
  const gdouble k2x_3E   = k * kx_3E;
  const gdouble k2x_3E2  = k2x_3E / E;
  const gdouble k2x2_3E2 = x * k2x_3E2;
  const gdouble psi      = -_NC_PHI - 12.0 * x2 / k2 * Omega_r0 * _NC_THETA2;
  const gdouble PI       = _NC_THETA2 + _NC_THETA_P0 + _NC_THETA_P2;
  const gdouble dErm2_dx = (3.0 * Omega_m0 * x2 + 4.0 * Omega_r0 * x3);
  const gdouble taubar   = nc_recomb_dtau_dlambda (pb->recomb, cosmo, lambda);
  guint i;

  if (pb->tight_coupling)
  {
    _NC_DPHI = psi - k2x2_3E2 * _NC_PHI - x / (2.0 * E2) *
               (
      dErm2_dx * (_NC_PHI - _NC_C0)
      - (3.0 * Omega_b0 * _NC_dB0 * x2 + 4.0 * Omega_r0 * x3 * _NC_dTHETA0)
               );

    _NC_DC0      = -kx_E * _NC_V + k2x2_3E2 * (_NC_C0 - _NC_PHI);
    _NC_DdTHETA0 = -kx_E * (R * _NC_U + _NC_T) / (R + 1.0);
    _NC_DdB0     = -kx_E * R * (_NC_U - _NC_T) / (R + 1.0);

    _NC_DV = -_NC_V - kx_3E * (kx_E * _NC_V - k2x2_3E2 * _NC_C0 + x3 * (3.0 * Omega_b0 * _NC_dB0 + 4.0 * Omega_r0 * x * _NC_dTHETA0) / (2.0 * E2));
    _NC_DU = kx_3E * (_NC_dTHETA0 - 2.0 * _NC_THETA2) + _NC_V;
    _NC_DT = kx_3E * (_NC_dTHETA0 - 2.0 * _NC_THETA2) + _NC_V + R * (_NC_U - _NC_T) / (R + 1.0) + taubar * (1.0 + R) * _NC_T;

    _NC_DTHETA2 = kx_E * (2.0 * ((R * _NC_U + _NC_T) / (R + 1.0) + _NC_V - kx_3E * (_NC_C0 - _NC_PHI)) - 3.0 * _NC_THETA (3)) / 5.0 + taubar * (_NC_THETA2 - PI / 10.0);
  }
  else
  {
    gdouble Delta = 0.0;

    _NC_DPHI = psi - k2x2_3E2 * _NC_PHI + x / (2.0 * E2) *
               (
      (3.0 * (Omega_c0 * _NC_C0 * x2 + Omega_b0 * _NC_B0 * x2) + 4.0 * Omega_r0 * x3 * _NC_THETA0)
      - dErm2_dx * _NC_PHI
               );

    _NC_DC0 = -kx_E * _NC_C1;
    _NC_DB0 = -kx_E * _NC_B1;

    _NC_DTHETA0 = -kx_E * _NC_THETA1;

    if (TRUE)
    {
      const gdouble R0x_onepR0x = R0 * x / (R0 * x + 1.0);
      const gdouble taunp       = (taubar * R0 + 1.0) * x - R0x_onepR0x;
      const gdouble Sigma       = _NC_THETA1 + _NC_B1 / (R0 * x);

      Delta = -1.0 / taunp * (R0x_onepR0x * Sigma + kx_3E * (_NC_THETA0 - 2.0 * _NC_THETA2 - _NC_PHI));


      printf ("% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n",
              -log (NC_HIPERT_BOLTZMANN_LAMBDA2X (lambda)),
              (_NC_THETA0 - _NC_PHI + _NC_THETA1 / kx_3E),
              _NC_THETA0, -_NC_PHI, _NC_THETA1 / kx_3E
             );

/*      printf ("% 20.15g % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e\n",
 *             -log (NC_HIPERT_BOLTZMANN_LAMBDA2X (lambda)),
 *             _NC_THETA0, _NC_THETA1,
 *             fabs ((_NC_B1 - _NC_THETA1)/_NC_THETA1),
 *             fabs ((_NC_C1 - _NC_THETA1)/_NC_THETA1),
 *             fabs ((_NC_C0 - _NC_THETA0)/_NC_THETA0),
 *             fabs ((_NC_B0 - _NC_THETA0)/_NC_THETA0), kx_E);
 *//*printf ("# % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e\n", R0x_onepR0x * Sigma, kx_3E * _NC_THETA0, - kx_3E * _NC_PHI, - kx_3E * 2.0 * _NC_THETA2, Delta); */

/*      printf ("% 20.15g % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e\n",
 *             -log (NC_HIPERT_BOLTZMANN_LAMBDA2X (lambda)),
 *             _NC_THETA0, _NC_B0, _NC_C0, _NC_THETA1, _NC_B1, _NC_C1, kx_E, nc_hipert_get_mode_k (pert));
 */
    }
    else
    {
      Delta = _NC_THETA1 - _NC_B1;
    }

    /*printf ("% 20.15g % 20.15e % 20.15g % 20.15g\n", lambda, NC_HIPERT_BOLTZMANN_LAMBDA2X (lambda), Delta, _NC_THETA1 - _NC_B1); */
    Delta = _NC_THETA1 - _NC_B1;

    _NC_DC1 = -_NC_C1 + kx_3E * psi;
    _NC_DB1 = -_NC_B1 + kx_3E * psi - taubar * R0 * x * Delta;

    _NC_DTHETA1 = kx_3E * (_NC_THETA0 - 2.0 * _NC_THETA2 - _NC_PHI + psi) + taubar * Delta;
    _NC_DTHETA2 = kx_E * (2.0 * _NC_THETA1 - 3.0 * _NC_THETA (3)) / 5.0 + taubar * (_NC_THETA2 - PI / 10.0);
  }

  /*printf ("% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", */
  /*        lambda, _NC_PHI, _NC_C0, _NC_B0, _NC_THETA0, _NC_C1, _NC_B1, _NC_THETA1); */
  /*printf ("% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", */
  /*        lambda, _NC_DPHI, _NC_DC0, _NC_DB0, _NC_DTHETA0, _NC_DC1, _NC_DB1, _NC_DTHETA1); */
/*  printf ("% 20.15g % 20.15g\n", _NC_DTHETA1, taubar * (_NC_THETA1 - _NC_B1)); */

  _NC_DTHETA_P0 = -kx_E * _NC_THETA_P1 + taubar * (_NC_THETA_P0 - PI / 2.0);
  _NC_DTHETA_P1 = kx_3E * (_NC_THETA_P0 - 2.0 * _NC_THETA_P2) + taubar * (_NC_THETA_P1);
  _NC_DTHETA_P2 = kx_E * (2.0 * _NC_THETA_P1 - 3.0 * _NC_THETA_P (3)) / 5.0 + taubar * (_NC_THETA_P2 - PI / 10.0);

  for (i = 3; i < lmax; i++)
  {
    const gdouble l  = i;
    const gdouble f1 = kx_E / (2.0 * l + 1.0);

    _NC_DTHETA (i)   = f1 * (l * _NC_THETA (i - 1)   - (1.0 + l) * _NC_THETA (i + 1))  + taubar * _NC_THETA (i);
    _NC_DTHETA_P (i) = f1 * (l * _NC_THETA_P (i - 1) - (1.0 + l) * _NC_THETA_P (i + 1)) + taubar * _NC_THETA_P (i);
    /*printf ("%u % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", i, _NC_THETA (i), _NC_DTHETA (i), f1 * l * _NC_THETA (i - 1), - f1 * (1.0 + l) * _NC_THETA (i + 1), taubar * _NC_THETA(i)); */
  }

  i = lmax; /* Not really needed but ... */
  {
    const gdouble l  = lmax;
    const gdouble f1 = kx_E / (2.0 * l + 1.0);

    _NC_DTHETA (i)   = f1 * l * _NC_THETA (i - 1)   + taubar * _NC_THETA (i);
    _NC_DTHETA_P (i) = f1 * l * _NC_THETA_P (i - 1) + taubar * _NC_THETA_P (i);
    /*printf ("%u % 20.15g % 20.15g\n", i, _NC_THETA (i), _NC_DTHETA (i)); */
  }

  if (FALSE)
  {
    const gdouble eta  = nc_scalefactor_eval_eta_x (pb->a, x);
    const gdouble keta = nc_hipert_get_mode_k (pert) * eta;

    guint n = 16;

    g_message ("%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n", lambda, k, eta, _NC_THETA (n + 1), (2.0 * n + 1.0) * _NC_THETA (n) / keta - _NC_THETA (n - 1), (2.0 * n + 1.0) * _NC_THETA (n) / keta, _NC_THETA (n - 1),
               -(n + 1.0) * kx_E * _NC_THETA (n) / (taubar * (2.0 * n + 3.0)), _NC_THETA (n), taubar);
    g_debug ("%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n", lambda, k, eta, _NC_THETA_P (n + 1), (2.0 * n + 1.0) * _NC_THETA_P (n) / keta - _NC_THETA_P (n - 1), (2.0 * n + 1.0) * _NC_THETA_P (n) / keta, _NC_THETA_P (n - 1),
             -(n + 1.0) * kx_E * _NC_THETA_P (n) / (taubar * (2.0 * n + 3.0)), _NC_THETA_P (n), taubar);
  }

  if (FALSE)
  {
    g_message ("%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g ",
               lambda,
               _NC_DPHI / _NC_PHI,
               _NC_DC0 / _NC_C0,
               _NC_DB0 / _NC_B0,
               _NC_DTHETA0 / _NC_THETA0,
               _NC_DC1 / _NC_C1,
               _NC_DB1 / _NC_B1,
               _NC_DTHETA1 / _NC_THETA1,
               _NC_DTHETA2 / _NC_THETA2,
               _NC_DTHETA_P0 / _NC_THETA_P0,
               _NC_DTHETA_P1 / _NC_THETA_P1,
               _NC_DTHETA_P2 / _NC_THETA_P2);

    for (i = 3; i <= 16; i++)
      g_message ("%.15g ", _NC_DTHETA (i) / _NC_THETA (i));

    g_message ("%.15g\n", taubar);
  }

  return 0;
}

#define _NC_BAND_ELEM(J, i, j) SUN_BAND_ACCESS ((J), (gint) (i), (gint) (j))

static gint
_nc_hipert_boltzmann_std_band_J (sunrealtype lambda, N_Vector y, N_Vector fy, SUNMatrix J, void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHIPertBoltzmann *pb  = NC_HIPERT_BOLTZMANN (jac_data);
  NcHIPert *pert         = NC_HIPERT (pb);
  NcHICosmo *cosmo       = pb->cosmo;
  const guint lmax       = pb->TT_lmax;
  const gdouble Omega_r0 = nc_hicosmo_Omega_r0 (cosmo);
  const gdouble Omega_b0 = nc_hicosmo_Omega_b0 (cosmo);
  const gdouble Omega_c0 = nc_hicosmo_Omega_c0 (cosmo);
  const gdouble Omega_m0 = nc_hicosmo_Omega_m0 (cosmo);
  const gdouble R0       = 4.0 * Omega_r0 / (3.0 * Omega_b0);
  const gdouble x        = NC_HIPERT_BOLTZMANN_LAMBDA2X (lambda);
  const gdouble x2       = x * x;
  const gdouble x3       = x2 * x;
  const gdouble x4       = x3 * x;
  const gdouble k        = nc_hipert_get_mode_k (pert);
  const gdouble k2       = k * k;
  const gdouble E2       = nc_hicosmo_E2 (cosmo, x - 1.0);
  const gdouble E        = sqrt (E2);
  const gdouble k_E      = k / E;
  const gdouble kx_E     = x * k_E;
  const gdouble kx_3E    = kx_E / 3.0;
  const gdouble k2x_3E   = k * kx_3E;
  const gdouble k2x_3E2  = k2x_3E / E;
  const gdouble k2x2_3E2 = x * k2x_3E2;
  const gdouble dErm2_dx = (3.0 * Omega_m0 * x2 + 4.0 * Omega_r0 * x3);
  const gdouble taubar   = nc_recomb_dtau_dlambda (pb->recomb, cosmo, lambda);
  guint i;

  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_PHI, NC_HIPERT_BOLTZMANN_PHI)    = -1.0 - k2x2_3E2 - x * dErm2_dx / (2.0 * E2);
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_PHI, NC_HIPERT_BOLTZMANN_C0)     = 3.0 * Omega_c0 * x3 / (2.0 * E2);
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_PHI, NC_HIPERT_BOLTZMANN_B0)     = 3.0 * Omega_b0 * x3 / (2.0 * E2);
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_PHI, NC_HIPERT_BOLTZMANN_THETA0) = 2.0 * x4 * Omega_r0 / E2;
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_PHI, NC_HIPERT_BOLTZMANN_THETA2) = -12.0 * x2 / k2 * Omega_r0;

  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_C0, NC_HIPERT_BOLTZMANN_C1) = -kx_E;

  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_B0, NC_HIPERT_BOLTZMANN_B1) = -kx_E;

  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA0, NC_HIPERT_BOLTZMANN_THETA1) = -kx_E;

  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_C1, NC_HIPERT_BOLTZMANN_PHI)    = -kx_3E;
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_C1, NC_HIPERT_BOLTZMANN_THETA2) = -12.0 * x3 / (E * k) * Omega_r0;
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_C1, NC_HIPERT_BOLTZMANN_C1)     = -1.0;

  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_B1, NC_HIPERT_BOLTZMANN_PHI)    = -kx_3E;
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_B1, NC_HIPERT_BOLTZMANN_THETA1) = -taubar * R0 * x;
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_B1, NC_HIPERT_BOLTZMANN_THETA2) = -12.0 * x3 / (E * k) * Omega_r0;
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_B1, NC_HIPERT_BOLTZMANN_B1)     = -1.0 + taubar * R0 * x;

  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA1, NC_HIPERT_BOLTZMANN_PHI)    = -2.0 * kx_E / 3.0;
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA1, NC_HIPERT_BOLTZMANN_THETA0) = kx_E / 3.0;
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA1, NC_HIPERT_BOLTZMANN_B1)     = -taubar;
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA1, NC_HIPERT_BOLTZMANN_THETA1) = taubar;
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA1, NC_HIPERT_BOLTZMANN_THETA2) = -2.0 * kx_E / 3.0 - 4.0 * x3 / (E * k) * Omega_r0;

  /* _NC_DTHETA_P2 = kx_E * (2.0 * _NC_THETA_P1 - 3.0 * _NC_THETA_P(3)) / 5.0 + taubar * (_NC_THETA_P2 - PI / 10.0); */
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA2, NC_HIPERT_BOLTZMANN_THETA1)    = 2.0 * kx_E / 5.0;
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA2, NC_HIPERT_BOLTZMANN_THETA2)    = taubar * 9.0 / 10.0;
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA2, NC_HIPERT_BOLTZMANN_THETA (3)) = -3.0 * kx_E / 5.0;
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA2, NC_HIPERT_BOLTZMANN_THETA_P0)  = -taubar / 10.0;
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA2, NC_HIPERT_BOLTZMANN_THETA_P2)  = -taubar / 10.0;

  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA_P0, NC_HIPERT_BOLTZMANN_THETA2)   = -taubar / 2.0;
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA_P0, NC_HIPERT_BOLTZMANN_THETA_P0) =  taubar / 2.0;
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA_P0, NC_HIPERT_BOLTZMANN_THETA_P1) =  -kx_E;
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA_P0, NC_HIPERT_BOLTZMANN_THETA_P2) = -taubar / 2.0;

  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA_P1, NC_HIPERT_BOLTZMANN_THETA_P0) = kx_E / 3.0;
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA_P1, NC_HIPERT_BOLTZMANN_THETA_P1) = taubar;
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA_P1, NC_HIPERT_BOLTZMANN_THETA_P2) = -2.0 * kx_E / 3.0;

  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA_P2, NC_HIPERT_BOLTZMANN_THETA2)      = -taubar / 10.0;
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA_P2, NC_HIPERT_BOLTZMANN_THETA_P0)    = -taubar / 10.0;
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA_P2, NC_HIPERT_BOLTZMANN_THETA_P1)    = 2.0 * kx_E / 5.0;
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA_P2, NC_HIPERT_BOLTZMANN_THETA_P2)    = 9.0 * taubar / 10.0;
  _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA_P2, NC_HIPERT_BOLTZMANN_THETA_P (3)) = -3.0 * kx_E / 5.0;

  for (i = 3; i < lmax; i++)
  {
    const gdouble l  = i;
    const gdouble f1 = kx_E / (2.0 * l + 1.0);
    const gdouble f2 = f1 * (1.0 + l);

    /* _NC_DTHETA(i) = f1 * i * _NC_THETA(i-1) + taubar * _NC_THETA(i); */
    _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA (i), NC_HIPERT_BOLTZMANN_THETA (i - 1)) = f1 * l;
    _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA (i), NC_HIPERT_BOLTZMANN_THETA (i))     = taubar;
    _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA (i), NC_HIPERT_BOLTZMANN_THETA (i + 1)) = -f2;

    _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA_P (i), NC_HIPERT_BOLTZMANN_THETA_P (i - 1)) = f1 * l;
    _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA_P (i), NC_HIPERT_BOLTZMANN_THETA_P (i))     = taubar;
    _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA_P (i), NC_HIPERT_BOLTZMANN_THETA_P (i + 1)) = -f2;
  }

  i = lmax; /* Not really needed but ... */
  {
    const gdouble l  = i;
    const gdouble f1 = kx_E / (2.0 * i + 1.0);

    _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA (i), NC_HIPERT_BOLTZMANN_THETA (i - 1)) = f1 * l;
    _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA (i), NC_HIPERT_BOLTZMANN_THETA (i))     = taubar;

    _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA_P (i), NC_HIPERT_BOLTZMANN_THETA_P (i - 1)) = f1 * l;
    _NC_BAND_ELEM (J, NC_HIPERT_BOLTZMANN_THETA_P (i), NC_HIPERT_BOLTZMANN_THETA_P (i))     = taubar;
  }

  return 0;
}

