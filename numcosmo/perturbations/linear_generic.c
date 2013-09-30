/***************************************************************************
 *            linear_generic.c
 *
 *  Thu Nov 12 12:37:02 2009
 *  Copyright  2009  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com>
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

static void
LINEAR_NAME_SUFFIX (init) (NcLinearPert *pert)
{
  LINEAR_VECTOR_PREPARE;
  const gdouble x = NC_PERTURBATIONS_LAMBDA2X (pert->lambdai);
  const gdouble z = x - 1.0;
  const gdouble E2 = nc_hicosmo_E2 (pert->cosmo, z);
  const gdouble E = sqrt (E2);
  const gdouble Omega_r = nc_hicosmo_Omega_r (pert->cosmo);
  const gdouble Omega_b = nc_hicosmo_Omega_b (pert->cosmo);
  const gdouble R0 = 4.0 * Omega_r / (3.0 * Omega_b);
  const gdouble taubar = nc_recomb_dtau_dlambda (pert->recomb, pert->cosmo, pert->lambdai);
  const gdouble taudot = E / x * taubar;
  const gdouble keta = pert->pws->k / (sqrt (Omega_r) * x);
  const gdouble keta3 = gsl_pow_3 (keta);
  const gdouble k_taudot = pert->pws->k / taudot;
  gdouble theta1;
  guint i;

  pert->pws->lambda_int = pert->lambdai;
  pert->pws->lambda     = pert->lambdai;

  LINEAR_VECTOR_SET_ALL (y, 0.0, pert->sys_size);

  pert->pws->tight_coupling = FALSE;
  pert->pws->tight_coupling_end = FALSE;

  _NC_PHI = 1.0;
  _NC_C0 = 3.0 * _NC_PHI / 2.0;

  if (pert->pws->tight_coupling)
  {
    _NC_V = keta3 * _NC_C0 / 36.0;
    _NC_U = _NC_V / 3.0;
    _NC_T = -(4.0 * _NC_V) / (3.0 * taubar * R0 * x);
    _NC_dB0 = _NC_dTHETA0 = -keta * _NC_V / 12.0;

    theta1 = -keta / 6.0 * _NC_PHI;
  }
  else if (FALSE)
  {
    _NC_THETA0 = _NC_B0 = _NC_C0;
    _NC_THETA1 = _NC_B1 = _NC_C1 = theta1 = -keta / 6.0 * _NC_PHI;
  }
  else
  {
    _NC_THETA0 = _NC_B0 = _NC_C0;
    _NC_THETA1 = _NC_C1 = theta1 = -keta / 6.0 * _NC_PHI;
    _NC_B1 = 1e-80;
    _NC_THETA0 = 1e-80;
  }

  _NC_THETA2 = -8.0 / 15.0 * k_taudot * theta1;
  _NC_THETA_P0 = 5.0 / 4.0 * _NC_THETA2;
  _NC_THETA_P1 = -1.0 / 4.0 * k_taudot * _NC_THETA2;
  _NC_THETA_P2 = 1.0 / 4.0 * _NC_THETA2;

  for (i = 3; i <= pert->lmax; i++)
  {
    const gdouble f1 = i * 1.0 / (2.0 * i + 1.0);
    _NC_THETA(i) = f1 * _NC_THETA(i-1);
    _NC_THETA_P(i) = f1 * _NC_THETA_P(i-1);
  }

#ifdef SIMUL_LOS_INT
  {
    guint last_l = 2;
    gdouble thetal = _NC_THETA2;
    for (i = 0; i < pert->los_table->len; i++)
    {
      guint l = g_array_index(pert->los_table, guint, i);
      if (l == 0)
        LINEAR_VEC_COMP(LINEAR_VEC_LOS_THETA(pert), i) = _NC_THETA0;
      else if (l == 1)
        LINEAR_VEC_COMP(LINEAR_VEC_LOS_THETA(pert), i) = _NC_THETA1;
      else if (l == 2)
        LINEAR_VEC_COMP(LINEAR_VEC_LOS_THETA(pert), i) = _NC_THETA2;
      else
      {
        guint j;
        for (j = last_l + 1; j <= l; j++)
        {
          const gdouble f1 = j * 1.0 / (2.0 * j + 1.0);
          thetal *= f1;
          last_l = j;
        }
        LINEAR_VEC_COMP(LINEAR_VEC_LOS_THETA(pert), i) = thetal;
      }
    }
  }
#endif

  LINEAR_NAME_SUFFIX (reset) (pert);

  return;
}

static void
LINEAR_NAME_SUFFIX (end_tight_coupling) (NcLinearPert *pert)
{
  LINEAR_VECTOR_PREPARE;
  const gdouble x = NC_PERTURBATIONS_LAMBDA2X (pert->pws->lambda);
  const gdouble Omega_r = nc_hicosmo_Omega_r (pert->cosmo);
  const gdouble Omega_b = nc_hicosmo_Omega_b (pert->cosmo);
  const gdouble R0 = 4.0 * Omega_r / (3.0 * Omega_b);
  const gdouble R = R0 * x;
  const gdouble E2 = nc_hicosmo_E2 (pert->cosmo, x-1.0);
  const gdouble kx_3E = pert->pws->k * x / (3.0 * sqrt (E2));
//  const gdouble taubar = nc_recomb_taubar (pert->recomb, lambda);

  if (FALSE)
    return;

  if (!pert->pws->tight_coupling)
    return;

  printf ("# tight coupling ended at %f => %f\n", pert->pws->lambda, exp (-pert->pws->lambda));

  pert->pws->tight_coupling = FALSE;
  pert->pws->tight_coupling_end = FALSE;

  {
    gdouble b0, theta0, c1, b1, theta1;
    b0 = _NC_dB0 + _NC_C0;
    theta0 = _NC_dTHETA0 + _NC_C0;

    c1 = _NC_V - kx_3E * (_NC_C0 - _NC_PHI);
    b1 = (R * (_NC_U - _NC_T)) / (R + 1.0) + _NC_V - kx_3E * (_NC_C0 - _NC_PHI);
    theta1 = (R * _NC_U + _NC_T) / (R + 1.0) + _NC_V - kx_3E * (_NC_C0 - _NC_PHI);

    _NC_B0 = b0;
    _NC_THETA0 = theta0;
    _NC_C1 = c1;
    _NC_B1 = b1;
    _NC_THETA1 = theta1;
  }

  LINEAR_VEC_COMP (LINEAR_VEC_ABSTOL(pert), NC_PERT_T) = 0.0;
#ifdef SIMUL_LOS_INT
  LINEAR_NAME_SUFFIX(update_los) (pert);
#endif

  {
    guint i;
    pert->reltol = 1e-4;
    LINEAR_VECTOR_SET_ALL (LINEAR_VEC_ABSTOL(pert), 1e-13, pert->sys_size);
    for (i = 0; i <= NC_PERT_THETA_P2; i++)
      LINEAR_VEC_COMP (LINEAR_VEC_ABSTOL(pert), i) = 0.0;
  }

  LINEAR_NAME_SUFFIX (reset) (pert);

  return;
}

gdouble
LINEAR_NAME_SUFFIX (get_k_H) (NcLinearPert *pert)
{
  const gdouble k = pert->pws->k;
  const gdouble c_H0 = ncm_c_c () / (nc_hicosmo_H0 (pert->cosmo) * 1.0e3);
  return k / c_H0;
}

gdouble
LINEAR_NAME_SUFFIX (get_z) (NcLinearPert *pert)
{
  return NC_PERTURBATIONS_LAMBDA2X (pert->pws->lambda) - 1.0;
}

gdouble
LINEAR_NAME_SUFFIX (get) (NcLinearPert *pert, guint n)
{
  switch (n)
  {
    case NC_PERT_PHI:
      return LINEAR_NAME_SUFFIX(get_phi) (pert);
    case NC_PERT_THETA0:
      return LINEAR_NAME_SUFFIX(get_theta) (pert, 0);
    case NC_PERT_B0:
      return LINEAR_NAME_SUFFIX(get_b0) (pert);
    case NC_PERT_C0:
      return LINEAR_NAME_SUFFIX(get_c0) (pert);
    case NC_PERT_THETA1:
      return LINEAR_NAME_SUFFIX(get_theta) (pert, 1);
    case NC_PERT_B1:
      return LINEAR_NAME_SUFFIX(get_b1) (pert);
    case NC_PERT_C1:
      return LINEAR_NAME_SUFFIX(get_c1) (pert);
    case NC_PERT_THETA2:
      return LINEAR_NAME_SUFFIX(get_theta2) (pert);
    default:
      g_error ("get n = %u not implemented", n);
  }
}

gdouble
LINEAR_NAME_SUFFIX (get_phi) (NcLinearPert *pert)
{
  LINEAR_VECTOR_PREPARE;
  return _NC_PHI;
}

gdouble
LINEAR_NAME_SUFFIX (get_c0) (NcLinearPert *pert)
{
  LINEAR_VECTOR_PREPARE;
  return _NC_C0 - _NC_PHI;
}

gdouble
LINEAR_NAME_SUFFIX (get_b0) (NcLinearPert *pert)
{
  LINEAR_VECTOR_PREPARE;
  if (pert->pws->tight_coupling)
    return _NC_dB0 + (_NC_C0 - _NC_PHI);
  else
    return _NC_B0 - _NC_PHI;
}

gdouble
LINEAR_NAME_SUFFIX (get_c1) (NcLinearPert *pert)
{
  LINEAR_VECTOR_PREPARE;
  if (pert->pws->tight_coupling)
  {
    const gdouble x = NC_PERTURBATIONS_LAMBDA2X (pert->pws->lambda);
    const gdouble E2 = nc_hicosmo_E2 (pert->cosmo, x-1.0);
    const gdouble kx_3E = pert->pws->k * x / (3.0 * sqrt (E2));

    return _NC_V - kx_3E * (_NC_C0 - _NC_PHI);
  }
  else
    return _NC_C1;
}

gdouble
LINEAR_NAME_SUFFIX (get_b1) (NcLinearPert *pert)
{
  LINEAR_VECTOR_PREPARE;
  if (pert->pws->tight_coupling)
  {
    const gdouble x = NC_PERTURBATIONS_LAMBDA2X (pert->pws->lambda);
    const gdouble Omega_r = nc_hicosmo_Omega_r (pert->cosmo);
    const gdouble Omega_b = nc_hicosmo_Omega_b (pert->cosmo);
    const gdouble R0 = 4.0 * Omega_r / (3.0 * Omega_b);
    const gdouble R = R0 * x;
    const gdouble E2 = nc_hicosmo_E2 (pert->cosmo, x-1.0);
    const gdouble kx_3E = pert->pws->k * x / (3.0 * sqrt (E2));

    return (R * (_NC_U - _NC_T)) / (R + 1.0) + _NC_V - kx_3E * (_NC_C0 - _NC_PHI);
  }
  else
    return _NC_B1;
}

gdouble
LINEAR_NAME_SUFFIX (get_theta2) (NcLinearPert *pert)
{
  LINEAR_VECTOR_PREPARE;
  return _NC_THETA(2);
}

gdouble
LINEAR_NAME_SUFFIX (get_theta) (NcLinearPert *pert, guint l)
{
  LINEAR_VECTOR_PREPARE;
  if (pert->pws->tight_coupling && l < 2)
  {
    const gdouble x = NC_PERTURBATIONS_LAMBDA2X (pert->pws->lambda);
    const gdouble Omega_r = nc_hicosmo_Omega_r (pert->cosmo);
    const gdouble Omega_b = nc_hicosmo_Omega_b (pert->cosmo);
    const gdouble R0 = 4.0 * Omega_r / (3.0 * Omega_b);
    const gdouble R = R0 * x;
    const gdouble E2 = nc_hicosmo_E2 (pert->cosmo, x-1.0);
    const gdouble kx_3E = pert->pws->k * x / (3.0 * sqrt (E2));

    if (l == 0)
      return _NC_dTHETA0 + (_NC_C0 - _NC_PHI);
    else if (l == 1)
      return (R * _NC_U + _NC_T) / (R + 1.0) + _NC_V - kx_3E * (_NC_C0 - _NC_PHI);
    else
      return _NC_THETA(l);
  }
  else
    if (l == 0)
      return _NC_THETA0 - _NC_PHI;
    else
      return _NC_THETA(l);
}

gdouble
LINEAR_NAME_SUFFIX (get_theta_p) (NcLinearPert *pert, guint l)
{
  LINEAR_VECTOR_PREPARE;
  return _NC_THETA_P (l);
}

void
LINEAR_NAME_SUFFIX (print_all) (NcLinearPert *pert)
{
	gdouble lambda = pert->pws->lambda;
	guint i;
	LINEAR_VECTOR_PREPARE;
  {
    g_message ("%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g ", lambda, _NC_PHI, _NC_C0, _NC_B0, _NC_THETA0, _NC_THETA_P0, _NC_C1, _NC_B1, _NC_THETA1, _NC_THETA_P1, _NC_THETA2, _NC_THETA_P2);
    for (i = 3; i <= 16; i++)
      g_message ("%.15g %.15g ", _NC_THETA (i), _NC_THETA_P (i));
		g_message ("\n");
  }
}

gdouble
LINEAR_NAME_SUFFIX (get_los_theta) (NcLinearPert *pert, guint l)
{
  return LINEAR_VEC_COMP(LINEAR_VEC_LOS_THETA(pert), l);
}

static void
LINEAR_NAME_SUFFIX (get_sources) (NcLinearPert *pert, gdouble *S0, gdouble *S1, gdouble *S2)
{
  LINEAR_VECTOR_PREPARE;
  const gdouble x = NC_PERTURBATIONS_LAMBDA2X (pert->pws->lambda);
  const gdouble taubar = nc_recomb_dtau_dlambda (pert->recomb, pert->cosmo, pert->pws->lambda);
  gdouble opt = nc_recomb_tau (pert->recomb, pert->cosmo, pert->pws->lambda);
  gdouble tau_log_abs_taubar = -opt + log(fabs(taubar));

  if (tau_log_abs_taubar > GSL_LOG_DBL_MIN + 0.01)
  {
    const gdouble Omega_r = nc_hicosmo_Omega_r (pert->cosmo);
    const gdouble Omega_b = nc_hicosmo_Omega_b (pert->cosmo);
    const gdouble Omega_c = nc_hicosmo_Omega_c (pert->cosmo);
    const gdouble Omega_m = nc_hicosmo_Omega_m (pert->cosmo);
    const gdouble x2 = x*x;
    const gdouble x3 = x2*x;
    const gdouble k = pert->pws->k;
    const gdouble k2 = k*k;
    const gdouble E2 = nc_hicosmo_E2 (pert->cosmo, x-1.0);
    const gdouble E = sqrt(E2);
    const gdouble k_E = k / E;
    const gdouble kx_E = x * k_E;
    const gdouble kx_3E = kx_E / 3.0;
    const gdouble k2x_3E = k * kx_3E;
    const gdouble k2x_3E2 = k2x_3E / E;
    const gdouble k2x2_3E2 = x * k2x_3E2;
    const gdouble dErm2_dx = (3.0 * Omega_m * x2 + 4.0 * Omega_r * x3);
    const gdouble psi = -_NC_PHI - 12.0 * x2 / k2 * Omega_r * _NC_THETA2;
    const gdouble PI = _NC_THETA2 + _NC_THETA_P0 + _NC_THETA_P2;
    const gdouble exp_tau = exp (-opt);
    const gdouble taubar_exp_tau = GSL_SIGN(taubar) * exp(tau_log_abs_taubar);
    const gdouble R0 = 4.0 * Omega_r / (3.0 * Omega_b);
    const gdouble R = R0 * x;
    gdouble dphi, b1, theta0;

    if (pert->pws->tight_coupling)
    {
      dphi = psi - k2x2_3E2 * _NC_PHI - x / (2.0 * E2) *
      (
        dErm2_dx * (_NC_PHI - _NC_C0)
        -(3.0 * Omega_b * _NC_dB0 * x2 + 4.0 * Omega_r * x3 * _NC_dTHETA0)
        );
      theta0 = _NC_dTHETA0 + (_NC_C0 - _NC_PHI);
      b1 = (R * (_NC_U - _NC_T)) / (R + 1.0) + _NC_V - kx_3E * (_NC_C0 - _NC_PHI);
    }
    else
    {
      dphi = psi - k2x2_3E2 * _NC_PHI - x / (2.0 * E2) *
      (
        dErm2_dx * _NC_PHI
        -(3.0 * (Omega_c * _NC_C0 * x2 + Omega_b * _NC_B0 * x2) + 4.0 * Omega_r * x3 * _NC_THETA0)
        );
      theta0 = _NC_THETA0 - _NC_PHI;
      b1 = _NC_B1;
    }

    *S0 = E / x * (-exp_tau * dphi - taubar_exp_tau * (theta0 + PI / 4.0));
    *S1 = E / x * (exp_tau * kx_E * psi - taubar_exp_tau * 3.0 * b1);
    *S2 = E / x * (-taubar_exp_tau * PI * 3.0 / 4.0);
  }
  else
    *S0 = *S1 = *S2 = 0.0;

  return;
}

static LINEAR_STEP_RET
LINEAR_NAME_SUFFIX (step) (LINEAR_STEP_PARAMS)
{
  NcLinearPert *pert = NC_LINEAR_PERTURBATIONS (user_data);
  const guint lmax = pert->lmax;
  const gdouble Omega_r = nc_hicosmo_Omega_r (pert->cosmo);
  const gdouble Omega_b = nc_hicosmo_Omega_b (pert->cosmo);
  const gdouble Omega_c = nc_hicosmo_Omega_c (pert->cosmo);
  const gdouble Omega_m = nc_hicosmo_Omega_m (pert->cosmo);
  const gdouble R0 = 4.0 * Omega_r / (3.0 * Omega_b);
  const gdouble x = NC_PERTURBATIONS_LAMBDA2X (lambda);
  const gdouble R = R0 * x;
  const gdouble x2 = x*x;
  const gdouble x3 = x2*x;
  const gdouble k = pert->pws->k;
  const gdouble k2 = k*k;
  const gdouble E2 = nc_hicosmo_E2 (pert->cosmo, x-1.0);
  const gdouble E = sqrt(E2);
  const gdouble k_E = k / E;
  const gdouble kx_E = x * k_E;
  const gdouble kx_3E = kx_E / 3.0;
  const gdouble k2x_3E = k * kx_3E;
  const gdouble k2x_3E2 = k2x_3E / E;
  const gdouble k2x2_3E2 = x * k2x_3E2;
  const gdouble psi = -_NC_PHI - 12.0 * x2 / k2 * Omega_r * _NC_THETA2;
  const gdouble PI = _NC_THETA2 + _NC_THETA_P0 + _NC_THETA_P2;
  const gdouble dErm2_dx = (3.0 * Omega_m * x2 + 4.0 * Omega_r * x3);
  const gdouble taubar = nc_recomb_dtau_dlambda (pert->recomb, pert->cosmo, lambda);
  guint i;

  if (fabs(_NC_V * _NC_TIGHT_COUPLING_END) > kx_3E * _NC_C0)
    pert->pws->tight_coupling_end = TRUE;

  if (pert->pws->tight_coupling)
  {
    _NC_DPHI = psi - k2x2_3E2 * _NC_PHI - x / (2.0 * E2) *
      (
        dErm2_dx * (_NC_PHI - _NC_C0)
        -(3.0 * Omega_b * _NC_dB0 * x2 + 4.0 * Omega_r * x3 * _NC_dTHETA0)
        );

    _NC_DC0 = -kx_E * _NC_V + k2x2_3E2 * (_NC_C0 - _NC_PHI);
    _NC_DdTHETA0 = -kx_E * (R * _NC_U + _NC_T) / (R + 1.0);
    _NC_DdB0 = -kx_E * R * (_NC_U - _NC_T) / (R + 1.0);

    _NC_DV = -_NC_V - kx_3E * (kx_E * _NC_V - k2x2_3E2 * _NC_C0 + x3 * (3.0 * Omega_b * _NC_dB0 + 4.0 * Omega_r * x * _NC_dTHETA0) / (2.0 * E2));
    _NC_DU = kx_3E * (_NC_dTHETA0 - 2.0 * _NC_THETA2) + _NC_V;
    _NC_DT = kx_3E * (_NC_dTHETA0 - 2.0 * _NC_THETA2) + _NC_V + R * (_NC_U - _NC_T) / (R + 1.0) + taubar * (1.0 + R) * _NC_T;

    _NC_DTHETA2 = kx_E * (2.0 * ((R * _NC_U + _NC_T) / (R + 1.0) + _NC_V - kx_3E * (_NC_C0 - _NC_PHI)) - 3.0 * _NC_THETA(3)) / 5.0 + taubar * (_NC_THETA2 - PI / 10.0);
  }
  else if (FALSE)
  {
    _NC_DPHI = psi - k2x2_3E2 * _NC_PHI + x / (2.0 * E2) *
      (
        (3.0 * (Omega_c * _NC_C0 * x2 + Omega_b * _NC_B0 * x2) + 4.0 * Omega_r * x3 * _NC_THETA0)
        -dErm2_dx * _NC_PHI
        );

    _NC_DC0 = -kx_E * _NC_C1;
    _NC_DB0 = -kx_E * _NC_B1;
    
    _NC_DTHETA0 = -kx_E * _NC_THETA1;

    _NC_DC1 = -_NC_C1 + kx_3E * psi;
    _NC_DB1 = -_NC_B1 + kx_3E * psi + taubar * R0 * x * (_NC_B1 - _NC_THETA1);
    
    _NC_DTHETA1 = kx_3E * (_NC_THETA0 - 2.0 * _NC_THETA2 - _NC_PHI + psi) + taubar * (_NC_THETA1 - _NC_B1);
    _NC_DTHETA2 = kx_E * (2.0 * _NC_THETA1 - 3.0 * _NC_THETA(3)) / 5.0 + taubar * (_NC_THETA2 - PI / 10.0);
  }
  else
  {
    const gdouble theta0 = _NC_THETA0 + _NC_PHI - _NC_THETA1 / kx_3E;
    const gdouble d3E_kx = 1.0 / kx_3E - 3.0 * nc_hicosmo_dE2_dz (pert->cosmo, x) / (2.0 * E * k);
      
    _NC_DPHI = psi - k2x2_3E2 * _NC_PHI + x / (2.0 * E2) *
      (
        (3.0 * (Omega_c * _NC_C0 * x2 + Omega_b * _NC_B0 * x2) + 4.0 * Omega_r * x3 * theta0)
        -dErm2_dx * _NC_PHI
        );

    _NC_DC0 = -kx_E * _NC_C1;
    _NC_DB0 = -kx_E * (_NC_B1 + _NC_THETA1);

    _NC_DC1 = -_NC_C1 + kx_3E * psi;
    _NC_DB1 = -_NC_B1 - kx_3E * (_NC_THETA0 - 2.0 * _NC_THETA2) + taubar * (R0 * x + 1.0) * _NC_B1;
    
    _NC_DTHETA1 = _NC_THETA1 + kx_3E * (_NC_THETA0 - 2.0 * _NC_THETA2 + psi) - taubar * _NC_B1;

    _NC_DTHETA2 = kx_E * (2.0 * _NC_THETA1 - 3.0 * _NC_THETA(3)) / 5.0 + taubar * (_NC_THETA2 - PI / 10.0);

    _NC_DTHETA0 = -kx_E * _NC_THETA1 - _NC_DPHI + _NC_DTHETA1 / kx_3E + d3E_kx * _NC_THETA1;
    
fprintf (stderr, "# % 1.7e % 1.7e % 1.7e % 1.7e % 1.7e % 1.7e % 1.7e % 1.7e\n", 
         lambda, x, _NC_THETA0, _NC_B1, _NC_THETA1, _NC_DTHETA0, _NC_DB1, _NC_DTHETA1);
fprintf (stderr, "# % 1.7e % 1.7e % 1.7e\n", 
         -_NC_B1, - kx_3E * (_NC_THETA0 - 2.0 * _NC_THETA2), + taubar * (R0 * x + 1.0) * _NC_B1);

  }


  /*
printf ("%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n",
        g,
        kx_3E * _NC_THETA0,
        - 2.0 * kx_3E * _NC_THETA2,
        - kx_3E * _NC_PHI,
        + kx_3E * psi,
        taubar * _NC_THETA1,
        - taubar * _NC_B1,
        _NC_THETA1, _NC_DTHETA1);
	*/
//printf ("% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", lambda, _NC_PHI, _NC_DPHI, _NC_DC0, _NC_C0);

  _NC_DTHETA_P0 = -kx_E * _NC_THETA_P1 + taubar * (_NC_THETA_P0 - PI / 2.0);
  _NC_DTHETA_P1 = kx_3E * (_NC_THETA_P0 - 2.0 * _NC_THETA_P2) + taubar * (_NC_THETA_P1);
  _NC_DTHETA_P2 = kx_E * (2.0 * _NC_THETA_P1 - 3.0 * _NC_THETA_P(3)) / 5.0 + taubar * (_NC_THETA_P2 - PI / 10.0);

  for (i = 3; i < lmax; i++)
  {
    const gdouble f1 = kx_E / (2.0 * i + 1.0);
    _NC_DTHETA(i) =
      f1 * (
            i * _NC_THETA(i-1) -
            (1.0 + i) * _NC_THETA (i+1)
            ) +
      taubar * _NC_THETA(i);
    _NC_DTHETA_P(i) =
      f1 *
      (
       i * _NC_THETA_P(i-1) -
       (1.0 + i) * _NC_THETA_P (i+1)) +
      taubar * _NC_THETA_P(i);
  }

  i = lmax; /* Not really needed but ... */
  {
    const gdouble f1 = kx_E / (2.0 * i + 1.0);
#ifdef _NC_USE_CUTOFF
    if (_NC_USE_CUTOFF_TEST)
    {
      const gdouble keta = pert->pws->k * nc_scale_factor_t_x (pert->a, x);
      _NC_DTHETA(i) = kx_E * (_NC_THETA(i-1) - (i + 1.0) * _NC_THETA(i) / keta) + taubar * _NC_THETA(i);
      _NC_DTHETA_P(i) = kx_E * (_NC_THETA_P(i-1) - (i + 1.0) * _NC_THETA_P(i) / keta) + taubar * _NC_THETA_P(i);
    }
    else
#endif
    {
      _NC_DTHETA(i) =
        f1 * i * _NC_THETA(i-1) +
        taubar * _NC_THETA(i);
      _NC_DTHETA_P(i) =
        f1 * i * _NC_THETA_P(i-1) +
        taubar * _NC_THETA_P(i);
    }
  }

  if (FALSE)
  {
    const gdouble eta = nc_scale_factor_t_x (pert->a, x);
    const gdouble keta = pert->pws->k * eta;

    guint n = 16;
    g_message ("%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n", lambda, k, eta, _NC_THETA(n + 1), (2.0 * n+ 1.0) * _NC_THETA(n) / keta - _NC_THETA(n - 1), (2.0 * n+ 1.0) * _NC_THETA(n) / keta, _NC_THETA(n - 1),
            - (n + 1.0) * kx_E * _NC_THETA(n) / (taubar * (2.0 * n + 3.0)), _NC_THETA(n), taubar);
    g_debug ("%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n", lambda, k, eta, _NC_THETA_P(n + 1), (2.0 * n+ 1.0) * _NC_THETA_P(n) / keta - _NC_THETA_P(n - 1), (2.0 * n+ 1.0) * _NC_THETA_P(n) / keta, _NC_THETA_P(n - 1),
            - (n + 1.0) * kx_E * _NC_THETA_P(n) / (taubar * (2.0 * n + 3.0)), _NC_THETA_P(n), taubar);
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
      g_message ("%.15g ", _NC_DTHETA(i) / _NC_THETA(i));
    g_message ("%.15g\n", taubar);
  }

  LINEAR_STEP_RET_VAL;
}

gint
LINEAR_NAME_SUFFIX (band_J) (LINEAR_JAC_PARAMS)
{
  NcLinearPert *pert = (NcLinearPert *)user_data;
  const guint lmax = pert->lmax;
  const gdouble Omega_r = nc_hicosmo_Omega_r (pert->cosmo);
  const gdouble Omega_b = nc_hicosmo_Omega_b (pert->cosmo);
  const gdouble Omega_c = nc_hicosmo_Omega_c (pert->cosmo);
  const gdouble Omega_m = nc_hicosmo_Omega_m (pert->cosmo);
  const gdouble R0 = 4.0 * Omega_r / (3.0 * Omega_b);
  const gdouble x = NC_PERTURBATIONS_LAMBDA2X (lambda);
  const gdouble R = R0 * x;
  const gdouble x2 = x*x;
  const gdouble x3 = x2*x;
  const gdouble x4 = x3*x;
  const gdouble k = pert->pws->k;
  const gdouble k2 = k*k;
  const gdouble E2 = nc_hicosmo_E2 (pert->cosmo, x - 1.0);
  const gdouble E = sqrt(E2);
  const gdouble k_E = k / E;
  const gdouble kx_E = x * k_E;
  const gdouble kx_3E = kx_E / 3.0;
  const gdouble k2x_3E = k * kx_3E;
  const gdouble k2x_3E2 = k2x_3E / E;
  const gdouble k2x2_3E2 = x * k2x_3E2;
  const gdouble dErm2_dx = (3.0 * Omega_m * x2 + 4.0 * Omega_r * x3);
  const gdouble taubar = nc_recomb_dtau_dlambda (pert->recomb, pert->cosmo, lambda);
  guint i;

  LINEAR_JAC_UNUSED;

  if (pert->pws->tight_coupling)
  {
    /*
     * const gdouble psi = -_NC_PHI - 12.0 * x2 / k2 * Omega_r * _NC_THETA2;
     * _NC_DPHI = psi - k2x2_3E2 * _NC_PHI - x / (2.0 * E2) *
     *  (
       *    dErm2_dx * (_NC_PHI - _NC_C0)
     *   -(3.0 * Omega_b * _NC_dB0 * x2 + 4.0 * Omega_r * x3 * _NC_dTHETA0)
     *   );
     */
    LINEAR_MATRIX_E (J, NC_PERT_PHI, NC_PERT_PHI)      = -1.0 - k2x2_3E2 - x * dErm2_dx / (2.0 * E2);
    LINEAR_MATRIX_E (J, NC_PERT_PHI, NC_PERT_C0)       = +x * dErm2_dx / (2.0 * E2);
    LINEAR_MATRIX_E (J, NC_PERT_PHI, NC_PERT_dB0)      = +3.0 * Omega_b * x3 / (2.0 * E2);
    LINEAR_MATRIX_E (J, NC_PERT_PHI, NC_PERT_dTHETA0)  = +2.0 * x4 * Omega_r / E2;
    LINEAR_MATRIX_E (J, NC_PERT_PHI, NC_PERT_THETA2)   = -12.0 * x2 / k2 * Omega_r;

    /* _NC_DC0 = -kx_E * _NC_V + k2x2_3E2 * (_NC_C0 - _NC_PHI); */
    LINEAR_MATRIX_E (J, NC_PERT_C0, NC_PERT_PHI)       = -k2x2_3E2;
    LINEAR_MATRIX_E (J, NC_PERT_C0, NC_PERT_C0)        = +k2x2_3E2;
    LINEAR_MATRIX_E (J, NC_PERT_C0, NC_PERT_V)         = -kx_E;

    /* _NC_DdB0 = -kx_E * R * (_NC_U - _NC_T) / (R + 1.0); */
    LINEAR_MATRIX_E (J, NC_PERT_dB0, NC_PERT_U)        = -kx_E * R / (R + 1.0);
    LINEAR_MATRIX_E (J, NC_PERT_dB0, NC_PERT_T)        = +kx_E * R / (R + 1.0);

    /* _NC_DdTHETA0 = -kx_E * (R * _NC_U + _NC_T) / (R + 1.0); */
    LINEAR_MATRIX_E (J, NC_PERT_dTHETA0, NC_PERT_U)    = -kx_E * R / (R + 1.0);
    LINEAR_MATRIX_E (J, NC_PERT_dTHETA0, NC_PERT_T)    = -kx_E / (R + 1.0);

    /* _NC_DV = -_NC_V - kx_3E * (kx_E * _NC_V - k2x2_3E2 * _NC_C0 + x3 * (3.0 * Omega_b * _NC_dB0 + 4.0 * Omega_r * x * _NC_dTHETA0) / (2.0 * E2)); */
    LINEAR_MATRIX_E (J, NC_PERT_V, NC_PERT_dB0)        = -kx_E * x3 *  Omega_b / (2.0 * E2);
    LINEAR_MATRIX_E (J, NC_PERT_V, NC_PERT_dTHETA0)    = -kx_3E * x4 * 2.0 * Omega_r / E2;
    LINEAR_MATRIX_E (J, NC_PERT_V, NC_PERT_C0)         = +kx_3E * k2x2_3E2;
    LINEAR_MATRIX_E (J, NC_PERT_V, NC_PERT_V)          = -1.0 - k2x2_3E2;

    /* _NC_DU = kx_3E * (_NC_dTHETA0 - 2.0 * _NC_THETA2) + _NC_V; */
    LINEAR_MATRIX_E (J, NC_PERT_U, NC_PERT_dTHETA0)    = +kx_3E;
    LINEAR_MATRIX_E (J, NC_PERT_U, NC_PERT_V)          = +1.0;
    LINEAR_MATRIX_E (J, NC_PERT_U, NC_PERT_THETA2)     = -2.0 * kx_3E;

    /* _NC_DT = kx_3E * (_NC_dTHETA0 - 2.0 * _NC_THETA2) + _NC_V + R * (_NC_U - _NC_T) / (R + 1.0) + taubar * (1.0 + R) * _NC_T;  */
    LINEAR_MATRIX_E (J, NC_PERT_T, NC_PERT_dTHETA0)    = +kx_3E;
    LINEAR_MATRIX_E (J, NC_PERT_T, NC_PERT_V)          = +1.0;
    LINEAR_MATRIX_E (J, NC_PERT_T, NC_PERT_U)          = +R / (R + 1.0);
    LINEAR_MATRIX_E (J, NC_PERT_T, NC_PERT_T)          = -R / (R + 1.0) + (1.0 + R0 * x) * taubar;
    LINEAR_MATRIX_E (J, NC_PERT_T, NC_PERT_THETA2)     = -2.0 * kx_3E;

    /*
     * _NC_DTHETA2 = kx_E * (2.0 * ((R * _NC_U + _NC_T) / (R + 1.0) + _NC_V - kx_3E * (_NC_C0 - _NC_PHI)) -
     * 3.0 * _NC_THETA(3)) / 5.0 + taubar * (_NC_THETA2 - PI / 10.0);
     */
    LINEAR_MATRIX_E (J, NC_PERT_THETA2, NC_PERT_PHI)   = +2.0 * k2x2_3E2 / 5.0;
    LINEAR_MATRIX_E (J, NC_PERT_THETA2, NC_PERT_C0)    = -2.0 * k2x2_3E2 / 5.0;
    LINEAR_MATRIX_E (J, NC_PERT_THETA2, NC_PERT_V)     = +2.0 * kx_E / 5.0;
    LINEAR_MATRIX_E (J, NC_PERT_THETA2, NC_PERT_U)     = +2.0 * kx_E / 5.0 * R / (1.0 + R);
    LINEAR_MATRIX_E (J, NC_PERT_THETA2, NC_PERT_T)     = +2.0 * kx_E / 5.0 / (1.0 + R);

    LINEAR_MATRIX_E (J, NC_PERT_THETA2, NC_PERT_THETA2)   = taubar * 9.0 / 10.0;
    LINEAR_MATRIX_E (J, NC_PERT_THETA2, NC_PERT_THETA(3)) = -3.0 * kx_E / 5.0;
    LINEAR_MATRIX_E (J, NC_PERT_THETA2, NC_PERT_THETA_P0) = -taubar / 10.0;
    LINEAR_MATRIX_E (J, NC_PERT_THETA2, NC_PERT_THETA_P2) = -taubar / 10.0;
  }
  else
  {
    LINEAR_MATRIX_E (J, NC_PERT_PHI, NC_PERT_PHI)         = -1.0 - k2x2_3E2 - x * dErm2_dx / (2.0 * E2);
    LINEAR_MATRIX_E (J, NC_PERT_PHI, NC_PERT_C0)          = 3.0 * Omega_c * x3 / (2.0 * E2);
    LINEAR_MATRIX_E (J, NC_PERT_PHI, NC_PERT_B0)          = 3.0 * Omega_b * x3 / (2.0 * E2);
    LINEAR_MATRIX_E (J, NC_PERT_PHI, NC_PERT_THETA0)      = 2.0 * x4 * Omega_r / E2;
    LINEAR_MATRIX_E (J, NC_PERT_PHI, NC_PERT_THETA2)      = -12.0 * x2 / k2 * Omega_r;

    LINEAR_MATRIX_E (J, NC_PERT_C0, NC_PERT_C1)           = -kx_E;

    LINEAR_MATRIX_E (J, NC_PERT_B0, NC_PERT_B1)           = -kx_E;

    LINEAR_MATRIX_E (J, NC_PERT_THETA0, NC_PERT_THETA1)   = -kx_E;

    LINEAR_MATRIX_E (J, NC_PERT_C1, NC_PERT_PHI)          = -kx_3E;
    LINEAR_MATRIX_E (J, NC_PERT_C1, NC_PERT_THETA2)       = -12.0 * x3 / (E*k) * Omega_r;
    LINEAR_MATRIX_E (J, NC_PERT_C1, NC_PERT_C1)           = -1.0;

    LINEAR_MATRIX_E (J, NC_PERT_B1, NC_PERT_PHI)          = -kx_3E;
    LINEAR_MATRIX_E (J, NC_PERT_B1, NC_PERT_THETA1)       = -taubar * R0 * x;
    LINEAR_MATRIX_E (J, NC_PERT_B1, NC_PERT_THETA2)       = -12.0 * x3 / (E*k) * Omega_r;
    LINEAR_MATRIX_E (J, NC_PERT_B1, NC_PERT_B1)           = -1.0 + taubar * R0 * x;

    LINEAR_MATRIX_E (J, NC_PERT_THETA1, NC_PERT_PHI)      = -2.0 * kx_E / 3.0;
    LINEAR_MATRIX_E (J, NC_PERT_THETA1, NC_PERT_THETA0)   = kx_E / 3.0;
    LINEAR_MATRIX_E (J, NC_PERT_THETA1, NC_PERT_B1)       = -taubar;
    LINEAR_MATRIX_E (J, NC_PERT_THETA1, NC_PERT_THETA1)   = taubar;
    LINEAR_MATRIX_E (J, NC_PERT_THETA1, NC_PERT_THETA2)   = -2.0 * kx_E / 3.0 - 4.0 * x3 / (E*k) * Omega_r;

    /* _NC_DTHETA_P2 = kx_E * (2.0 * _NC_THETA_P1 - 3.0 * _NC_THETA_P(3)) / 5.0 + taubar * (_NC_THETA_P2 - PI / 10.0); */
    LINEAR_MATRIX_E (J, NC_PERT_THETA2, NC_PERT_THETA1)   = 2.0 * kx_E / 5.0;
    LINEAR_MATRIX_E (J, NC_PERT_THETA2, NC_PERT_THETA2)   = taubar * 9.0 / 10.0;
    LINEAR_MATRIX_E (J, NC_PERT_THETA2, NC_PERT_THETA(3)) = -3.0 * kx_E / 5.0;
    LINEAR_MATRIX_E (J, NC_PERT_THETA2, NC_PERT_THETA_P0) = -taubar / 10.0;
    LINEAR_MATRIX_E (J, NC_PERT_THETA2, NC_PERT_THETA_P2) = -taubar / 10.0;
  }

  LINEAR_MATRIX_E (J, NC_PERT_THETA_P0, NC_PERT_THETA2)     = -taubar / 2.0;
  LINEAR_MATRIX_E (J, NC_PERT_THETA_P0, NC_PERT_THETA_P0)   =  taubar / 2.0;
  LINEAR_MATRIX_E (J, NC_PERT_THETA_P0, NC_PERT_THETA_P1)   =  -kx_E;
  LINEAR_MATRIX_E (J, NC_PERT_THETA_P0, NC_PERT_THETA_P2)   = -taubar / 2.0;

  LINEAR_MATRIX_E (J, NC_PERT_THETA_P1, NC_PERT_THETA_P0)   = kx_E / 3.0;
  LINEAR_MATRIX_E (J, NC_PERT_THETA_P1, NC_PERT_THETA_P1)   = taubar;
  LINEAR_MATRIX_E (J, NC_PERT_THETA_P1, NC_PERT_THETA_P2)   = -2.0 * kx_E / 3.0;

  LINEAR_MATRIX_E (J, NC_PERT_THETA_P2, NC_PERT_THETA2)     = -taubar / 10.0;
  LINEAR_MATRIX_E (J, NC_PERT_THETA_P2, NC_PERT_THETA_P0)   = -taubar / 10.0;
  LINEAR_MATRIX_E (J, NC_PERT_THETA_P2, NC_PERT_THETA_P1)   = 2.0 * kx_E / 5.0;
  LINEAR_MATRIX_E (J, NC_PERT_THETA_P2, NC_PERT_THETA_P2)   = 9.0 * taubar / 10.0;
  LINEAR_MATRIX_E (J, NC_PERT_THETA_P2, NC_PERT_THETA_P(3)) = -3.0 * kx_E / 5.0;

  for (i = 3; i < lmax; i++)
  {
    const gdouble f1 = kx_E / (2.0 * i + 1.0);
    const gdouble f2 = f1 * (1.0 + i);

    /* _NC_DTHETA(i) = f1 * i * _NC_THETA(i-1) + taubar * _NC_THETA(i); */
    LINEAR_MATRIX_E (J, NC_PERT_THETA(i), NC_PERT_THETA(i-1)) = f1 * i;
    LINEAR_MATRIX_E (J, NC_PERT_THETA(i), NC_PERT_THETA(i))   = taubar;
    LINEAR_MATRIX_E (J, NC_PERT_THETA(i), NC_PERT_THETA(i+1)) = -f2;

    LINEAR_MATRIX_E (J, NC_PERT_THETA_P(i), NC_PERT_THETA_P(i-1)) = f1 * i;
    LINEAR_MATRIX_E (J, NC_PERT_THETA_P(i), NC_PERT_THETA_P(i))   = taubar;
    LINEAR_MATRIX_E (J, NC_PERT_THETA_P(i), NC_PERT_THETA_P(i+1)) = -f2;
  }

  i = lmax; /* Not really needed but ... */
  {
    const gdouble f1 = kx_E / (2.0 * i + 1.0);

#ifdef _NC_USE_CUTOFF
    if (_NC_USE_CUTOFF_TEST)
    {
      const gdouble keta = pert->pws->k * nc_scale_factor_t_x (pert->a, x);

      LINEAR_MATRIX_E (J, NC_PERT_THETA(i), NC_PERT_THETA(i-1)) = kx_E;
      LINEAR_MATRIX_E (J, NC_PERT_THETA(i), NC_PERT_THETA(i))   = -(i + 1.0) * kx_E / keta + taubar;

      LINEAR_MATRIX_E (J, NC_PERT_THETA_P(i), NC_PERT_THETA_P(i-1)) = kx_E;
      LINEAR_MATRIX_E (J, NC_PERT_THETA_P(i), NC_PERT_THETA_P(i))   = -(i + 1.0) * kx_E / keta + taubar;
    }
    else
#endif
    {
      LINEAR_MATRIX_E (J, NC_PERT_THETA(i), NC_PERT_THETA(i-1)) = f1 * i;
      LINEAR_MATRIX_E (J, NC_PERT_THETA(i), NC_PERT_THETA(i))   = taubar;

      LINEAR_MATRIX_E (J, NC_PERT_THETA_P(i), NC_PERT_THETA_P(i-1)) = f1 * i;
      LINEAR_MATRIX_E (J, NC_PERT_THETA_P(i), NC_PERT_THETA_P(i))   = taubar;
    }
  }

  return 0;
}
