/***************************************************************************
 *            nc_cbe_precision.c
 *
 *  Sun October 25 20:44:43 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_cbe_precision.c
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

/**
 * SECTION:nc_cbe_precision
 * @title: NcCBEPrecision
 * @short_description: CLASS (Cosmic Linear Anisotropy Solving System) backend for perturbations
 *
 * This object provides a front-end for CLASS precision structure.
 *
 * If you use this object please cite: [Blas (2011) CLASS II][XBlas2011],
 * see also:
 * - [Lesgourgues (2011) CLASS I][XLesgourgues2011],
 * - [Lesgourgues (2011) CLASS III][XLesgourgues2011a],
 * - [Lesgourgues (2011) CLASS IV][XLesgourgues2011b] and
 * - [CLASS website](http://class-code.net/).
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */

/*
 * It must be include before anything else, several symbols clash
 * with the default includes.
 */
#include "class/include/class.h"

#include "build_cfg.h"

#include "nc_cbe_precision.h"
#include "math/ncm_cfg.h"

struct _NcCBEPrecisionPrivate
{
  struct precision ppr;
};

enum
{
  PROP_0,
  PROP_A_INI_A_0,
  PROP_BACK_INT_STEP,
  PROP_BACK_TOL,
  PROP_INITIAL_OMEGA_R_TOL,
  PROP_M_NCDM_TOL,
  PROP_NCDM_TOL,
  PROP_NCDM_SYNCHRONOUS_TOL,
  PROP_NCDM_NEWTONIAN_TOL,
  PROP_NCDM_BG_TOL,
  PROP_NCDM_INITIAL_W_TOL,
  PROP_SAFE_PHI_SCF,
  PROP_TAU_EQ_TOL,
  PROP_SBBN_FILE,
  PROP_RECFAST_Z_INI,
  PROP_RECFAST_NZ0,
  PROP_THERMO_INTEGRATION_TOL,
  PROP_RECFAST_HE_SWITCH,
  PROP_RECFAST_FUDGE_HE,
  PROP_RECFAST_H_SWITCH,
  PROP_RECFAST_FUDGE_H,
  PROP_RECFAST_DELTA_FUDGE_H,
  PROP_RECFAST_A_GAUSS_1,
  PROP_RECFAST_A_GAUSS_2,
  PROP_RECFAST_Z_GAUSS_1,
  PROP_RECFAST_Z_GAUSS_2,
  PROP_RECFAST_W_GAUSS_1,
  PROP_RECFAST_W_GAUSS_2,
  PROP_RECFAST_Z_HE_1,
  PROP_RECFAST_DELTA_Z_HE_1,
  PROP_RECFAST_Z_HE_2,
  PROP_RECFAST_DELTA_Z_HE_2,
  PROP_RECFAST_Z_HE_3,
  PROP_RECFAST_DELTA_Z_HE_3,
  PROP_RECFAST_X_HE0_TRIGGER,
  PROP_RECFAST_X_HE0_TRIGGER2,
  PROP_RECFAST_X_HE0_TRIGGER_DELTA,
  PROP_RECFAST_X_H0_TRIGGER,
  PROP_RECFAST_X_H0_TRIGGER2,
  PROP_RECFAST_X_H0_TRIGGER_DELTA,
  PROP_RECFAST_H_FRAC,
  PROP_HYREC_ALPHA_INF_FILE,
  PROP_HYREC_R_INF_FILE,
  PROP_HYREC_TWO_PHOTON_TABLES_FILE,
  PROP_REION_Z_START_MAX,
  PROP_REION_SAMPLING,
  PROP_REION_OPT_DEPTH_TOL,
  PROP_REION_START_FACTOR,
  PROP_THERMO_RATE_SMOOTHING_RADIUS,
  PROP_EVOLVER,
  PROP_K_MIN_TAU0,
  PROP_K_MAX_TAU0_OVER_L_MAX,
  PROP_K_STEP_SUB,
  PROP_K_STEP_SUPER,
  PROP_K_STEP_TRANSITION,
  PROP_K_STEP_SUPER_REDUCTION,
  PROP_K_PER_DECADE_FOR_PK,
  PROP_K_PER_DECADE_FOR_BAO,
  PROP_K_BAO_CENTER,
  PROP_K_BAO_WIDTH,
  PROP_START_SMALL_K_AT_TAU_C_OVER_TAU_H,
  PROP_START_LARGE_K_AT_TAU_H_OVER_TAU_K,
  PROP_TIGHT_COUPLING_TRIGGER_TAU_C_OVER_TAU_H,
  PROP_TIGHT_COUPLING_TRIGGER_TAU_C_OVER_TAU_K,
  PROP_START_SOURCES_AT_TAU_C_OVER_TAU_H,
  PROP_TIGHT_COUPLING_APPROXIMATION,
  PROP_L_MAX_G,
  PROP_L_MAX_POL_G,
  PROP_L_MAX_DR,
  PROP_L_MAX_UR,
  PROP_L_MAX_NCDM,
  PROP_L_MAX_G_TEN,
  PROP_L_MAX_POL_G_TEN,
  PROP_CURVATURE_INI,
  PROP_ENTROPY_INI,
  PROP_GW_INI,
  PROP_PERTURB_INTEGRATION_STEPSIZE,
  PROP_TOL_TAU_APPROX,
  PROP_TOL_PERTURB_INTEGRATION,
  PROP_PERTURB_SAMPLING_STEPSIZE,
  PROP_RADIATION_STREAMING_APPROXIMATION,
  PROP_RADIATION_STREAMING_TRIGGER_TAU_OVER_TAU_K,
  PROP_RADIATION_STREAMING_TRIGGER_TAU_C_OVER_TAU,
  PROP_UR_FLUID_APPROXIMATION,
  PROP_UR_FLUID_TRIGGER_TAU_OVER_TAU_K,
  PROP_NCDM_FLUID_APPROXIMATION,
  PROP_NCDM_FLUID_TRIGGER_TAU_OVER_TAU_K,
  PROP_NEGLECT_CMB_SOURCES_BELOW_VISIBILITY,
  PROP_K_PER_DECADE_PRIMORDIAL,
  PROP_PRIMORDIAL_INFLATION_RATIO_MIN,
  PROP_PRIMORDIAL_INFLATION_RATIO_MAX,
  PROP_PRIMORDIAL_INFLATION_PHI_INI_MAXIT,
  PROP_PRIMORDIAL_INFLATION_PT_STEPSIZE,
  PROP_PRIMORDIAL_INFLATION_BG_STEPSIZE,
  PROP_PRIMORDIAL_INFLATION_TOL_INTEGRATION,
  PROP_PRIMORDIAL_INFLATION_ATTRACTOR_PRECISION_PIVOT,
  PROP_PRIMORDIAL_INFLATION_ATTRACTOR_PRECISION_INITIAL,
  PROP_PRIMORDIAL_INFLATION_ATTRACTOR_MAXIT,
  PROP_PRIMORDIAL_INFLATION_TOL_CURVATURE,
  PROP_PRIMORDIAL_INFLATION_AH_INI_TARGET,
  PROP_PRIMORDIAL_INFLATION_END_DPHI,
  PROP_PRIMORDIAL_INFLATION_END_LOGSTEP,
  PROP_PRIMORDIAL_INFLATION_SMALL_EPSILON,
  PROP_PRIMORDIAL_INFLATION_SMALL_EPSILON_TOL,
  PROP_PRIMORDIAL_INFLATION_EXTRA_EFOLDS,
  PROP_L_LOGSTEP,
  PROP_L_LINSTEP,
  PROP_HYPER_X_MIN,
  PROP_HYPER_SAMPLING_FLAT,
  PROP_HYPER_SAMPLING_CURVED_LOW_NU,
  PROP_HYPER_SAMPLING_CURVED_HIGH_NU,
  PROP_HYPER_NU_SAMPLING_STEP,
  PROP_HYPER_PHI_MIN_ABS,
  PROP_HYPER_X_TOL,
  PROP_HYPER_FLAT_APPROXIMATION_NU,
  PROP_Q_LINSTEP,
  PROP_Q_LOGSTEP_SPLINE,
  PROP_Q_LOGSTEP_OPEN,
  PROP_Q_LOGSTEP_TRAPZD,
  PROP_Q_NUMSTEP_TRANSITION,
  PROP_TRANSFER_NEGLECT_DELTA_K_S_T0,
  PROP_TRANSFER_NEGLECT_DELTA_K_S_T1,
  PROP_TRANSFER_NEGLECT_DELTA_K_S_T2,
  PROP_TRANSFER_NEGLECT_DELTA_K_S_E,
  PROP_TRANSFER_NEGLECT_DELTA_K_V_T1,
  PROP_TRANSFER_NEGLECT_DELTA_K_V_T2,
  PROP_TRANSFER_NEGLECT_DELTA_K_V_E,
  PROP_TRANSFER_NEGLECT_DELTA_K_V_B,
  PROP_TRANSFER_NEGLECT_DELTA_K_T_T2,
  PROP_TRANSFER_NEGLECT_DELTA_K_T_E,
  PROP_TRANSFER_NEGLECT_DELTA_K_T_B,
  PROP_TRANSFER_NEGLECT_LATE_SOURCE,
  PROP_L_SWITCH_LIMBER,
  PROP_L_SWITCH_LIMBER_FOR_NC_LOCAL_OVER_Z,
  PROP_L_SWITCH_LIMBER_FOR_NC_LOS_OVER_Z,
  PROP_SELECTION_CUT_AT_SIGMA,
  PROP_SELECTION_SAMPLING,
  PROP_SELECTION_SAMPLING_BESSEL,
  PROP_SELECTION_SAMPLING_BESSEL_LOS,
  PROP_SELECTION_TOPHAT_EDGE,
  PROP_HALOFIT_MIN_K_NONLINEAR,
  PROP_HALOFIT_MIN_K_MAX,
  PROP_HALOFIT_K_PER_DECADE,
  PROP_HALOFIT_SIGMA_PRECISION,
  PROP_HALOFIT_SIGMA_TOL,
  PROP_PK_EQ_Z_MAX,
  PROP_PK_EQ_TOL,
  PROP_ACCURATE_LENSING,
  PROP_NUM_MU_MINUS_LMAX,
  PROP_DELTA_L_MAX,
  PROP_TOL_GAUSS_LEGENDRE,
  PROP_SMALLEST_ALLOWED_VARIATION,
  PROP_SIZE
};

G_DEFINE_TYPE_WITH_PRIVATE (NcCBEPrecision, nc_cbe_precision, G_TYPE_OBJECT);

static void
nc_cbe_precision_init (NcCBEPrecision* cbe_prec)
{
  cbe_prec->priv = nc_cbe_precision_get_instance_private (cbe_prec);
}

static void
nc_cbe_precision_set_property (GObject* object, guint prop_id, const GValue* value, GParamSpec* pspec)
{
  NcCBEPrecision* cbe_prec = NC_CBE_PRECISION (object);
  g_return_if_fail (NC_IS_CBE_PRECISION (object));

  switch (prop_id)
  {
    case PROP_A_INI_A_0:
      cbe_prec->priv->ppr.a_ini_over_a_today_default = g_value_get_double (value);
      break;
    case PROP_BACK_INT_STEP:
      cbe_prec->priv->ppr.back_integration_stepsize  = g_value_get_double (value);
      break;
    case PROP_BACK_TOL:
      cbe_prec->priv->ppr.tol_background_integration = g_value_get_double (value);
      break;
    case PROP_INITIAL_OMEGA_R_TOL:
      cbe_prec->priv->ppr.tol_initial_Omega_r        = g_value_get_double (value);
      break;
    case PROP_M_NCDM_TOL:
      cbe_prec->priv->ppr.tol_M_ncdm                 = g_value_get_double (value);
      break;
    case PROP_NCDM_TOL:
      cbe_prec->priv->ppr.tol_ncdm                   = g_value_get_double (value);
      break;
    case PROP_NCDM_SYNCHRONOUS_TOL:
      cbe_prec->priv->ppr.tol_ncdm_synchronous       = g_value_get_double (value);
      break;
    case PROP_NCDM_NEWTONIAN_TOL:
      cbe_prec->priv->ppr.tol_ncdm_newtonian         = g_value_get_double (value);
      break;
    case PROP_NCDM_BG_TOL:
      cbe_prec->priv->ppr.tol_ncdm_bg                = g_value_get_double (value);
      break;
    case PROP_NCDM_INITIAL_W_TOL:
      cbe_prec->priv->ppr.tol_ncdm_initial_w         = g_value_get_double (value);
      break;
    case PROP_SAFE_PHI_SCF:
      cbe_prec->priv->ppr.safe_phi_scf               = g_value_get_double (value);
      break;
    case PROP_TAU_EQ_TOL:
      cbe_prec->priv->ppr.tol_tau_eq                 = g_value_get_double (value);
      break;
    case PROP_SBBN_FILE:
    {
      gchar *sBBN_file    = g_value_dup_string (value);
      guint sBBN_file_len = strlen (sBBN_file);

      g_assert_cmpuint (sBBN_file_len, <= ,_FILENAMESIZE_);

      if (!g_file_test (sBBN_file, G_FILE_TEST_EXISTS))
      {
        g_warning ("nc_cbe_precision_set_property: file `%s' not found, using default.", sBBN_file);
        g_free (sBBN_file);

        sBBN_file     = ncm_cfg_get_data_filename ("class_data"G_DIR_SEPARATOR_S"bbn"G_DIR_SEPARATOR_S"sBBN.dat", TRUE);
        sBBN_file_len = strlen (sBBN_file);
        g_assert_cmpuint (sBBN_file_len, <= ,_FILENAMESIZE_);
      }

      memcpy (cbe_prec->priv->ppr.sBBN_file, sBBN_file, sBBN_file_len + 1);
      g_free (sBBN_file);

      break;
    }
    case PROP_RECFAST_Z_INI:
      cbe_prec->priv->ppr.recfast_z_initial          = g_value_get_double (value);
      break;
    case PROP_RECFAST_NZ0:
      cbe_prec->priv->ppr.recfast_Nz0                = g_value_get_int (value);
      break;
    case PROP_THERMO_INTEGRATION_TOL:
      cbe_prec->priv->ppr.tol_thermo_integration     = g_value_get_double (value);
      break;
    case PROP_RECFAST_HE_SWITCH:
      cbe_prec->priv->ppr.recfast_Heswitch           = g_value_get_int (value);
      break;
    case PROP_RECFAST_FUDGE_HE:
      cbe_prec->priv->ppr.recfast_fudge_He           = g_value_get_double (value);
      break;
    case PROP_RECFAST_H_SWITCH:
      cbe_prec->priv->ppr.recfast_Hswitch            = g_value_get_int (value);
      break;
    case PROP_RECFAST_FUDGE_H:
      cbe_prec->priv->ppr.recfast_fudge_H            = g_value_get_double (value);
      break;
    case PROP_RECFAST_DELTA_FUDGE_H:
      cbe_prec->priv->ppr.recfast_delta_fudge_H      = g_value_get_double (value);
      break;
    case PROP_RECFAST_A_GAUSS_1:
      cbe_prec->priv->ppr.recfast_AGauss1            = g_value_get_double (value);
      break;
    case PROP_RECFAST_A_GAUSS_2:
      cbe_prec->priv->ppr.recfast_AGauss2            = g_value_get_double (value);
      break;
    case PROP_RECFAST_Z_GAUSS_1:
      cbe_prec->priv->ppr.recfast_zGauss1            = g_value_get_double (value);
      break;
    case PROP_RECFAST_Z_GAUSS_2:
      cbe_prec->priv->ppr.recfast_zGauss2            = g_value_get_double (value);
      break;
    case PROP_RECFAST_W_GAUSS_1:
      cbe_prec->priv->ppr.recfast_wGauss1            = g_value_get_double (value);
      break;
    case PROP_RECFAST_W_GAUSS_2:
      cbe_prec->priv->ppr.recfast_wGauss2            = g_value_get_double (value);
      break;
    case PROP_RECFAST_Z_HE_1:
      cbe_prec->priv->ppr.recfast_z_He_1             = g_value_get_double (value);
      break;
    case PROP_RECFAST_DELTA_Z_HE_1:
      cbe_prec->priv->ppr.recfast_delta_z_He_1       = g_value_get_double (value);
      break;
    case PROP_RECFAST_Z_HE_2:
      cbe_prec->priv->ppr.recfast_z_He_2             = g_value_get_double (value);
      break;
    case PROP_RECFAST_DELTA_Z_HE_2:
      cbe_prec->priv->ppr.recfast_delta_z_He_2       = g_value_get_double (value);
      break;
    case PROP_RECFAST_Z_HE_3:
      cbe_prec->priv->ppr.recfast_z_He_3             = g_value_get_double (value);
      break;
    case PROP_RECFAST_DELTA_Z_HE_3:
      cbe_prec->priv->ppr.recfast_delta_z_He_3       = g_value_get_double (value);
      break;
    case PROP_RECFAST_X_HE0_TRIGGER:
      cbe_prec->priv->ppr.recfast_x_He0_trigger      = g_value_get_double (value);
      break;
    case PROP_RECFAST_X_HE0_TRIGGER2:
      cbe_prec->priv->ppr.recfast_x_He0_trigger2     = g_value_get_double (value);
      break;
    case PROP_RECFAST_X_HE0_TRIGGER_DELTA:
      cbe_prec->priv->ppr.recfast_x_He0_trigger_delta = g_value_get_double (value);
      break;
    case PROP_RECFAST_X_H0_TRIGGER:
      cbe_prec->priv->ppr.recfast_x_H0_trigger       = g_value_get_double (value);
      break;
    case PROP_RECFAST_X_H0_TRIGGER2:
      cbe_prec->priv->ppr.recfast_x_H0_trigger2      = g_value_get_double (value);
      break;
    case PROP_RECFAST_X_H0_TRIGGER_DELTA:
      cbe_prec->priv->ppr.recfast_x_H0_trigger_delta = g_value_get_double (value);
      break;
    case PROP_RECFAST_H_FRAC:
      cbe_prec->priv->ppr.recfast_H_frac             = g_value_get_double (value);
      break;
    case PROP_HYREC_ALPHA_INF_FILE:
    {
      gchar *Alpha_inf_file    = g_value_dup_string (value);
      guint Alpha_inf_file_len = strlen (Alpha_inf_file);

      g_assert_cmpuint (Alpha_inf_file_len, <= ,_FILENAMESIZE_);

      if (!g_file_test (Alpha_inf_file, G_FILE_TEST_EXISTS))
      {
        g_warning ("nc_cbe_precision_set_property: file `%s' not found, using default.", Alpha_inf_file);
        g_free (Alpha_inf_file);

        Alpha_inf_file     = ncm_cfg_get_data_filename ("class_data"G_DIR_SEPARATOR_S"hyrec"G_DIR_SEPARATOR_S"Alpha_inf.dat", TRUE);
        Alpha_inf_file_len = strlen (Alpha_inf_file);
        g_assert_cmpuint (Alpha_inf_file_len, <= ,_FILENAMESIZE_);
      }

      memcpy (cbe_prec->priv->ppr.hyrec_Alpha_inf_file, Alpha_inf_file, Alpha_inf_file_len + 1);
      g_free (Alpha_inf_file);

      break;
    }
    case PROP_HYREC_R_INF_FILE:
    {
      gchar *R_inf_file    = g_value_dup_string (value);
      guint R_inf_file_len = strlen (R_inf_file);

      g_assert_cmpuint (R_inf_file_len, <= ,_FILENAMESIZE_);

      if (!g_file_test (R_inf_file, G_FILE_TEST_EXISTS))
      {
        g_warning ("nc_cbe_precision_set_property: file `%s' not found, using default.", R_inf_file);
        g_free (R_inf_file);

        R_inf_file     = ncm_cfg_get_data_filename ("class_data"G_DIR_SEPARATOR_S"hyrec"G_DIR_SEPARATOR_S"R_inf.dat", TRUE);
        R_inf_file_len = strlen (R_inf_file);

        g_assert_cmpuint (R_inf_file_len, <= ,_FILENAMESIZE_);
      }

      memcpy (cbe_prec->priv->ppr.hyrec_R_inf_file, R_inf_file, R_inf_file_len + 1);
      g_free (R_inf_file);

      break;
    }
    case PROP_HYREC_TWO_PHOTON_TABLES_FILE:
    {
      gchar *two_photon_tables_file    = g_value_dup_string (value);
      guint two_photon_tables_file_len = strlen (two_photon_tables_file);

      g_assert_cmpuint (two_photon_tables_file_len, <= ,_FILENAMESIZE_);

      if (!g_file_test (two_photon_tables_file, G_FILE_TEST_EXISTS))
      {
        g_warning ("nc_cbe_precision_set_property: file `%s' not found, using default.", two_photon_tables_file);
        g_free (two_photon_tables_file);

        two_photon_tables_file     = ncm_cfg_get_data_filename ("class_data"G_DIR_SEPARATOR_S"hyrec"G_DIR_SEPARATOR_S"two_photon_tables.dat", TRUE);
        two_photon_tables_file_len = strlen (two_photon_tables_file);

        g_assert_cmpuint (two_photon_tables_file_len, <= ,_FILENAMESIZE_);
      }

      memcpy (cbe_prec->priv->ppr.hyrec_two_photon_tables_file, two_photon_tables_file, two_photon_tables_file_len + 1);

      g_free (two_photon_tables_file);
      break;
    }
    case PROP_REION_Z_START_MAX:
      cbe_prec->priv->ppr.reionization_z_start_max   = g_value_get_double (value);
      break;
    case PROP_REION_SAMPLING:
      cbe_prec->priv->ppr.reionization_sampling      = g_value_get_double (value);
      break;
    case PROP_REION_OPT_DEPTH_TOL:
      cbe_prec->priv->ppr.reionization_optical_depth_tol = g_value_get_double (value);
      break;
    case PROP_REION_START_FACTOR:
      cbe_prec->priv->ppr.reionization_start_factor  = g_value_get_double (value);
      break;
    case PROP_THERMO_RATE_SMOOTHING_RADIUS:
      cbe_prec->priv->ppr.thermo_rate_smoothing_radius = g_value_get_int (value);
      break;
    case PROP_EVOLVER:
      cbe_prec->priv->ppr.evolver                    = g_value_get_int (value);
      break;
    case PROP_K_MIN_TAU0:
      cbe_prec->priv->ppr.k_min_tau0                 = g_value_get_double (value);
      break;
    case PROP_K_MAX_TAU0_OVER_L_MAX:
      cbe_prec->priv->ppr.k_max_tau0_over_l_max      = g_value_get_double (value);
      break;
    case PROP_K_STEP_SUB:
      cbe_prec->priv->ppr.k_step_sub                 = g_value_get_double (value);
      break;
    case PROP_K_STEP_SUPER:
      cbe_prec->priv->ppr.k_step_super               = g_value_get_double (value);
      break;
    case PROP_K_STEP_TRANSITION:
      cbe_prec->priv->ppr.k_step_transition          = g_value_get_double (value);
      break;
    case PROP_K_STEP_SUPER_REDUCTION:
      cbe_prec->priv->ppr.k_step_super_reduction     = g_value_get_double (value);
      break;
    case PROP_K_PER_DECADE_FOR_PK:
      cbe_prec->priv->ppr.k_per_decade_for_pk        = g_value_get_double (value);
      break;
    case PROP_K_PER_DECADE_FOR_BAO:
      cbe_prec->priv->ppr.k_per_decade_for_bao       = g_value_get_double (value);
      break;
    case PROP_K_BAO_CENTER:
      cbe_prec->priv->ppr.k_bao_center               = g_value_get_double (value);
      break;
    case PROP_K_BAO_WIDTH:
      cbe_prec->priv->ppr.k_bao_width                = g_value_get_double (value);
      break;
    case PROP_START_SMALL_K_AT_TAU_C_OVER_TAU_H:
      cbe_prec->priv->ppr.start_small_k_at_tau_c_over_tau_h       = g_value_get_double (value);
      break;
    case PROP_START_LARGE_K_AT_TAU_H_OVER_TAU_K:
      cbe_prec->priv->ppr.start_large_k_at_tau_h_over_tau_k       = g_value_get_double (value);
      break;
    case PROP_TIGHT_COUPLING_TRIGGER_TAU_C_OVER_TAU_H:
      cbe_prec->priv->ppr.tight_coupling_trigger_tau_c_over_tau_h = g_value_get_double (value);
      break;
    case PROP_TIGHT_COUPLING_TRIGGER_TAU_C_OVER_TAU_K:
      cbe_prec->priv->ppr.tight_coupling_trigger_tau_c_over_tau_k = g_value_get_double (value);
      break;
    case PROP_START_SOURCES_AT_TAU_C_OVER_TAU_H:
      cbe_prec->priv->ppr.start_sources_at_tau_c_over_tau_h       = g_value_get_double (value);
      break;
    case PROP_TIGHT_COUPLING_APPROXIMATION:
      cbe_prec->priv->ppr.tight_coupling_approximation            = g_value_get_int (value);
      break;
    case PROP_L_MAX_G:
      cbe_prec->priv->ppr.l_max_g                    = g_value_get_int (value);
      break;
    case PROP_L_MAX_POL_G:
      cbe_prec->priv->ppr.l_max_pol_g                = g_value_get_int (value);
      break;
    case PROP_L_MAX_DR:
      cbe_prec->priv->ppr.l_max_dr                   = g_value_get_int (value);
      break;
    case PROP_L_MAX_UR:
      cbe_prec->priv->ppr.l_max_ur                   = g_value_get_int (value);
      break;
    case PROP_L_MAX_NCDM:
      cbe_prec->priv->ppr.l_max_ncdm                 = g_value_get_int (value);
      break;
    case PROP_L_MAX_G_TEN:
      cbe_prec->priv->ppr.l_max_g_ten                = g_value_get_int (value);
      break;
    case PROP_L_MAX_POL_G_TEN:
      cbe_prec->priv->ppr.l_max_pol_g_ten            = g_value_get_int (value);
      break;
    case PROP_CURVATURE_INI:
      cbe_prec->priv->ppr.curvature_ini              = g_value_get_double (value);
      break;
    case PROP_ENTROPY_INI:
      cbe_prec->priv->ppr.entropy_ini                = g_value_get_double (value);
      break;
    case PROP_GW_INI:
      cbe_prec->priv->ppr.gw_ini                     = g_value_get_double (value);
      break;
    case PROP_PERTURB_INTEGRATION_STEPSIZE:
      cbe_prec->priv->ppr.perturb_integration_stepsize = g_value_get_double (value);
      break;
    case PROP_TOL_TAU_APPROX:
      cbe_prec->priv->ppr.tol_tau_approx             = g_value_get_double (value);
      break;
    case PROP_TOL_PERTURB_INTEGRATION:
      cbe_prec->priv->ppr.tol_perturb_integration    = g_value_get_double (value);
      break;
    case PROP_PERTURB_SAMPLING_STEPSIZE:
      cbe_prec->priv->ppr.perturb_sampling_stepsize  = g_value_get_double (value);
      break;
    case PROP_RADIATION_STREAMING_APPROXIMATION:
      cbe_prec->priv->ppr.radiation_streaming_approximation = g_value_get_int (value);
      break;
    case PROP_RADIATION_STREAMING_TRIGGER_TAU_OVER_TAU_K:
      cbe_prec->priv->ppr.radiation_streaming_trigger_tau_over_tau_k = g_value_get_double (value);
      break;
    case PROP_RADIATION_STREAMING_TRIGGER_TAU_C_OVER_TAU:
      cbe_prec->priv->ppr.radiation_streaming_trigger_tau_c_over_tau = g_value_get_double (value);
      break;
    case PROP_UR_FLUID_APPROXIMATION:
      cbe_prec->priv->ppr.ur_fluid_approximation     = g_value_get_int (value);
      break;
    case PROP_UR_FLUID_TRIGGER_TAU_OVER_TAU_K:
      cbe_prec->priv->ppr.ur_fluid_trigger_tau_over_tau_k = g_value_get_double (value);
      break;
    case PROP_NCDM_FLUID_APPROXIMATION:
      cbe_prec->priv->ppr.ncdm_fluid_approximation   = g_value_get_int (value);
      break;
    case PROP_NCDM_FLUID_TRIGGER_TAU_OVER_TAU_K:
      cbe_prec->priv->ppr.ncdm_fluid_trigger_tau_over_tau_k = g_value_get_double (value);
      break;
    case PROP_NEGLECT_CMB_SOURCES_BELOW_VISIBILITY:
      cbe_prec->priv->ppr.neglect_CMB_sources_below_visibility = g_value_get_double (value);
      break;
    case PROP_K_PER_DECADE_PRIMORDIAL:
      cbe_prec->priv->ppr.k_per_decade_primordial    = g_value_get_double (value);
      break;
    case PROP_PRIMORDIAL_INFLATION_RATIO_MIN:
      cbe_prec->priv->ppr.primordial_inflation_ratio_min = g_value_get_double (value);
      break;
    case PROP_PRIMORDIAL_INFLATION_RATIO_MAX:
      cbe_prec->priv->ppr.primordial_inflation_ratio_max = g_value_get_double (value);
      break;
    case PROP_PRIMORDIAL_INFLATION_PHI_INI_MAXIT:
      cbe_prec->priv->ppr.primordial_inflation_phi_ini_maxit = g_value_get_int (value);
      break;
    case PROP_PRIMORDIAL_INFLATION_PT_STEPSIZE:
      cbe_prec->priv->ppr.primordial_inflation_pt_stepsize = g_value_get_double (value);
      break;
    case PROP_PRIMORDIAL_INFLATION_BG_STEPSIZE:
      cbe_prec->priv->ppr.primordial_inflation_bg_stepsize = g_value_get_double (value);
      break;
    case PROP_PRIMORDIAL_INFLATION_TOL_INTEGRATION:
      cbe_prec->priv->ppr.primordial_inflation_tol_integration = g_value_get_double (value);
      break;
    case PROP_PRIMORDIAL_INFLATION_ATTRACTOR_PRECISION_PIVOT:
      cbe_prec->priv->ppr.primordial_inflation_attractor_precision_pivot = g_value_get_double (value);
      break;
    case PROP_PRIMORDIAL_INFLATION_ATTRACTOR_PRECISION_INITIAL:
      cbe_prec->priv->ppr.primordial_inflation_attractor_precision_initial = g_value_get_double (value);
      break;
    case PROP_PRIMORDIAL_INFLATION_ATTRACTOR_MAXIT:
      cbe_prec->priv->ppr.primordial_inflation_attractor_maxit = g_value_get_int (value);
      break;
    case PROP_PRIMORDIAL_INFLATION_TOL_CURVATURE:
      cbe_prec->priv->ppr.primordial_inflation_tol_curvature = g_value_get_double (value);
      break;
    case PROP_PRIMORDIAL_INFLATION_AH_INI_TARGET:
      cbe_prec->priv->ppr.primordial_inflation_aH_ini_target = g_value_get_double (value);
      break;
    case PROP_PRIMORDIAL_INFLATION_END_DPHI:
      cbe_prec->priv->ppr.primordial_inflation_end_dphi = g_value_get_double (value);
      break;
    case PROP_PRIMORDIAL_INFLATION_END_LOGSTEP:
      cbe_prec->priv->ppr.primordial_inflation_end_logstep = g_value_get_double (value);
      break;
    case PROP_PRIMORDIAL_INFLATION_SMALL_EPSILON:
      cbe_prec->priv->ppr.primordial_inflation_small_epsilon = g_value_get_double (value);
      break;
    case PROP_PRIMORDIAL_INFLATION_SMALL_EPSILON_TOL:
      cbe_prec->priv->ppr.primordial_inflation_small_epsilon_tol = g_value_get_double (value);
      break;
    case PROP_PRIMORDIAL_INFLATION_EXTRA_EFOLDS:
      cbe_prec->priv->ppr.primordial_inflation_extra_efolds = g_value_get_double (value);
      break;
    case PROP_L_LOGSTEP:
      cbe_prec->priv->ppr.l_logstep                  = g_value_get_double (value);
      break;
    case PROP_L_LINSTEP:
      cbe_prec->priv->ppr.l_linstep                  = g_value_get_int (value);
      break;
    case PROP_HYPER_X_MIN:
      cbe_prec->priv->ppr.hyper_x_min                = g_value_get_double (value);
      break;
    case PROP_HYPER_SAMPLING_FLAT:
      cbe_prec->priv->ppr.hyper_sampling_flat        = g_value_get_double (value);
      break;
    case PROP_HYPER_SAMPLING_CURVED_LOW_NU:
      cbe_prec->priv->ppr.hyper_sampling_curved_low_nu = g_value_get_double (value);
      break;
    case PROP_HYPER_SAMPLING_CURVED_HIGH_NU:
      cbe_prec->priv->ppr.hyper_sampling_curved_high_nu = g_value_get_double (value);
      break;
    case PROP_HYPER_NU_SAMPLING_STEP:
      cbe_prec->priv->ppr.hyper_nu_sampling_step     = g_value_get_double (value);
      break;
    case PROP_HYPER_PHI_MIN_ABS:
      cbe_prec->priv->ppr.hyper_phi_min_abs          = g_value_get_double (value);
      break;
    case PROP_HYPER_X_TOL:
      cbe_prec->priv->ppr.hyper_x_tol                = g_value_get_double (value);
      break;
    case PROP_HYPER_FLAT_APPROXIMATION_NU:
      cbe_prec->priv->ppr.hyper_flat_approximation_nu = g_value_get_double (value);
      break;
    case PROP_Q_LINSTEP:
      cbe_prec->priv->ppr.q_linstep                  = g_value_get_double (value);
      break;
    case PROP_Q_LOGSTEP_SPLINE:
      cbe_prec->priv->ppr.q_logstep_spline           = g_value_get_double (value);
      break;
    case PROP_Q_LOGSTEP_OPEN:
      cbe_prec->priv->ppr.q_logstep_open             = g_value_get_double (value);
      break;
    case PROP_Q_LOGSTEP_TRAPZD:
      cbe_prec->priv->ppr.q_logstep_trapzd           = g_value_get_double (value);
      break;
    case PROP_Q_NUMSTEP_TRANSITION:
      cbe_prec->priv->ppr.q_numstep_transition       = g_value_get_double (value);
      break;
    case PROP_TRANSFER_NEGLECT_DELTA_K_S_T0:
      cbe_prec->priv->ppr.transfer_neglect_delta_k_S_t0 = g_value_get_double (value);
      break;
    case PROP_TRANSFER_NEGLECT_DELTA_K_S_T1:
      cbe_prec->priv->ppr.transfer_neglect_delta_k_S_t1 = g_value_get_double (value);
      break;
    case PROP_TRANSFER_NEGLECT_DELTA_K_S_T2:
      cbe_prec->priv->ppr.transfer_neglect_delta_k_S_t2 = g_value_get_double (value);
      break;
    case PROP_TRANSFER_NEGLECT_DELTA_K_S_E:
      cbe_prec->priv->ppr.transfer_neglect_delta_k_S_e  = g_value_get_double (value);
      break;
    case PROP_TRANSFER_NEGLECT_DELTA_K_V_T1:
      cbe_prec->priv->ppr.transfer_neglect_delta_k_V_t1 = g_value_get_double (value);
      break;
    case PROP_TRANSFER_NEGLECT_DELTA_K_V_T2:
      cbe_prec->priv->ppr.transfer_neglect_delta_k_V_t2 = g_value_get_double (value);
      break;
    case PROP_TRANSFER_NEGLECT_DELTA_K_V_E:
      cbe_prec->priv->ppr.transfer_neglect_delta_k_V_e  = g_value_get_double (value);
      break;
    case PROP_TRANSFER_NEGLECT_DELTA_K_V_B:
      cbe_prec->priv->ppr.transfer_neglect_delta_k_V_b  = g_value_get_double (value);
      break;
    case PROP_TRANSFER_NEGLECT_DELTA_K_T_T2:
      cbe_prec->priv->ppr.transfer_neglect_delta_k_T_t2 = g_value_get_double (value);
      break;
    case PROP_TRANSFER_NEGLECT_DELTA_K_T_E:
      cbe_prec->priv->ppr.transfer_neglect_delta_k_T_e = g_value_get_double (value);
      break;
    case PROP_TRANSFER_NEGLECT_DELTA_K_T_B:
      cbe_prec->priv->ppr.transfer_neglect_delta_k_T_b = g_value_get_double (value);
      break;
    case PROP_TRANSFER_NEGLECT_LATE_SOURCE:
      cbe_prec->priv->ppr.transfer_neglect_late_source = g_value_get_double (value);
      break;
    case PROP_L_SWITCH_LIMBER:
      cbe_prec->priv->ppr.l_switch_limber            = g_value_get_double (value);
      break;
    case PROP_L_SWITCH_LIMBER_FOR_NC_LOCAL_OVER_Z:
      cbe_prec->priv->ppr.l_switch_limber_for_nc_local_over_z = g_value_get_double (value);
      break;
    case PROP_L_SWITCH_LIMBER_FOR_NC_LOS_OVER_Z:
      cbe_prec->priv->ppr.l_switch_limber_for_nc_los_over_z   = g_value_get_double (value);
      break;
    case PROP_SELECTION_CUT_AT_SIGMA:
      cbe_prec->priv->ppr.selection_cut_at_sigma     = g_value_get_double (value);
      break;
    case PROP_SELECTION_SAMPLING:
      cbe_prec->priv->ppr.selection_sampling         = g_value_get_double (value);
      break;
    case PROP_SELECTION_SAMPLING_BESSEL:
      cbe_prec->priv->ppr.selection_sampling_bessel  = g_value_get_double (value);
      break;
    case PROP_SELECTION_SAMPLING_BESSEL_LOS:
      cbe_prec->priv->ppr.selection_sampling_bessel_los = g_value_get_double (value);
      break;
    case PROP_SELECTION_TOPHAT_EDGE:
      cbe_prec->priv->ppr.selection_tophat_edge      = g_value_get_double (value);
      break;
    case PROP_HALOFIT_MIN_K_NONLINEAR:
      cbe_prec->priv->ppr.halofit_min_k_nonlinear    = g_value_get_double (value);
      break;
    case PROP_HALOFIT_MIN_K_MAX:
      cbe_prec->priv->ppr.halofit_min_k_max          = g_value_get_double (value);
      break;
    case PROP_HALOFIT_K_PER_DECADE:
      cbe_prec->priv->ppr.halofit_k_per_decade       = g_value_get_double (value);
      break;
    case PROP_HALOFIT_SIGMA_PRECISION:
      cbe_prec->priv->ppr.halofit_sigma_precision    = g_value_get_double (value);
      break;
    case PROP_HALOFIT_SIGMA_TOL:
      cbe_prec->priv->ppr.halofit_tol_sigma          = g_value_get_double (value);
      break;
    case PROP_PK_EQ_Z_MAX:
      cbe_prec->priv->ppr.pk_eq_z_max                = g_value_get_double (value);
      break;
    case PROP_PK_EQ_TOL:
      cbe_prec->priv->ppr.pk_eq_tol                  = g_value_get_double (value);
      break;
    case PROP_ACCURATE_LENSING:
      cbe_prec->priv->ppr.accurate_lensing           = g_value_get_int (value);
      break;
    case PROP_NUM_MU_MINUS_LMAX:
      cbe_prec->priv->ppr.num_mu_minus_lmax          = g_value_get_int (value);
      break;
    case PROP_DELTA_L_MAX:
      cbe_prec->priv->ppr.delta_l_max                = g_value_get_int (value);
      break;
    case PROP_SMALLEST_ALLOWED_VARIATION:
      cbe_prec->priv->ppr.smallest_allowed_variation = g_value_get_double (value);
      break;
    case PROP_TOL_GAUSS_LEGENDRE:
      cbe_prec->priv->ppr.tol_gauss_legendre         = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_cbe_precision_get_property (GObject* object, guint prop_id, GValue* value, GParamSpec* pspec)
{
  NcCBEPrecision* cbe_prec = NC_CBE_PRECISION (object);
  g_return_if_fail (NC_IS_CBE_PRECISION (object));

  switch (prop_id)
  {
    case PROP_A_INI_A_0:
      g_value_set_double (value, cbe_prec->priv->ppr.a_ini_over_a_today_default);
      break;
    case PROP_BACK_INT_STEP:
      g_value_set_double (value, cbe_prec->priv->ppr.back_integration_stepsize);
      break;
    case PROP_BACK_TOL:
      g_value_set_double (value, cbe_prec->priv->ppr.tol_background_integration);
      break;
    case PROP_INITIAL_OMEGA_R_TOL:
      g_value_set_double (value, cbe_prec->priv->ppr.tol_initial_Omega_r);
      break;
    case PROP_M_NCDM_TOL:
      g_value_set_double (value, cbe_prec->priv->ppr.tol_M_ncdm);
      break;
    case PROP_NCDM_TOL:
      g_value_set_double (value, cbe_prec->priv->ppr.tol_ncdm);
      break;
    case PROP_NCDM_SYNCHRONOUS_TOL:
      g_value_set_double (value, cbe_prec->priv->ppr.tol_ncdm_synchronous);
      break;
    case PROP_NCDM_NEWTONIAN_TOL:
      g_value_set_double (value, cbe_prec->priv->ppr.tol_ncdm_newtonian);
      break;
    case PROP_NCDM_BG_TOL:
      g_value_set_double (value, cbe_prec->priv->ppr.tol_ncdm_bg);
      break;
    case PROP_NCDM_INITIAL_W_TOL:
      g_value_set_double (value, cbe_prec->priv->ppr.tol_ncdm_initial_w);
      break;
    case PROP_SAFE_PHI_SCF:
      g_value_set_double (value, cbe_prec->priv->ppr.safe_phi_scf);
      break;
    case PROP_TAU_EQ_TOL:
      g_value_set_double (value, cbe_prec->priv->ppr.tol_tau_eq);
      break;
    case PROP_SBBN_FILE:
      g_value_set_string (value, cbe_prec->priv->ppr.sBBN_file);
      break;
    case PROP_RECFAST_Z_INI:
      g_value_set_double (value, cbe_prec->priv->ppr.recfast_z_initial);
      break;
    case PROP_RECFAST_NZ0:
      g_value_set_int (value, cbe_prec->priv->ppr.recfast_Nz0);
      break;
    case PROP_THERMO_INTEGRATION_TOL:
      g_value_set_double (value, cbe_prec->priv->ppr.tol_thermo_integration);
      break;
    case PROP_RECFAST_HE_SWITCH:
      g_value_set_int (value, cbe_prec->priv->ppr.recfast_Heswitch);
      break;
    case PROP_RECFAST_FUDGE_HE:
      g_value_set_double (value, cbe_prec->priv->ppr.recfast_fudge_He);
      break;
    case PROP_RECFAST_H_SWITCH:
      g_value_set_int (value, cbe_prec->priv->ppr.recfast_Hswitch);
      break;
    case PROP_RECFAST_FUDGE_H:
      g_value_set_double (value, cbe_prec->priv->ppr.recfast_fudge_H);
      break;
    case PROP_RECFAST_DELTA_FUDGE_H:
      g_value_set_double (value, cbe_prec->priv->ppr.recfast_delta_fudge_H);
      break;
    case PROP_RECFAST_A_GAUSS_1:
      g_value_set_double (value, cbe_prec->priv->ppr.recfast_AGauss1);
      break;
    case PROP_RECFAST_A_GAUSS_2:
      g_value_set_double (value, cbe_prec->priv->ppr.recfast_AGauss2);
      break;
    case PROP_RECFAST_Z_GAUSS_1:
      g_value_set_double (value, cbe_prec->priv->ppr.recfast_zGauss1);
      break;
    case PROP_RECFAST_Z_GAUSS_2:
      g_value_set_double (value, cbe_prec->priv->ppr.recfast_zGauss2);
      break;
    case PROP_RECFAST_W_GAUSS_1:
      g_value_set_double (value, cbe_prec->priv->ppr.recfast_wGauss1);
      break;
    case PROP_RECFAST_W_GAUSS_2:
      g_value_set_double (value, cbe_prec->priv->ppr.recfast_wGauss2);
      break;
    case PROP_RECFAST_Z_HE_1:
      g_value_set_double (value, cbe_prec->priv->ppr.recfast_z_He_1);
      break;
    case PROP_RECFAST_DELTA_Z_HE_1:
      g_value_set_double (value, cbe_prec->priv->ppr.recfast_delta_z_He_1);
      break;
    case PROP_RECFAST_Z_HE_2:
      g_value_set_double (value, cbe_prec->priv->ppr.recfast_z_He_2);
      break;
    case PROP_RECFAST_DELTA_Z_HE_2:
      g_value_set_double (value, cbe_prec->priv->ppr.recfast_delta_z_He_2);
      break;
    case PROP_RECFAST_Z_HE_3:
      g_value_set_double (value, cbe_prec->priv->ppr.recfast_z_He_3);
      break;
    case PROP_RECFAST_DELTA_Z_HE_3:
      g_value_set_double (value, cbe_prec->priv->ppr.recfast_delta_z_He_3);
      break;
    case PROP_RECFAST_X_HE0_TRIGGER:
      g_value_set_double (value, cbe_prec->priv->ppr.recfast_x_He0_trigger);
      break;
    case PROP_RECFAST_X_HE0_TRIGGER2:
      g_value_set_double (value, cbe_prec->priv->ppr.recfast_x_He0_trigger2);
      break;
    case PROP_RECFAST_X_HE0_TRIGGER_DELTA:
      g_value_set_double (value, cbe_prec->priv->ppr.recfast_x_He0_trigger_delta);
      break;
    case PROP_RECFAST_X_H0_TRIGGER:
      g_value_set_double (value, cbe_prec->priv->ppr.recfast_x_H0_trigger);
      break;
    case PROP_RECFAST_X_H0_TRIGGER2:
      g_value_set_double (value, cbe_prec->priv->ppr.recfast_x_H0_trigger2);
      break;
    case PROP_RECFAST_X_H0_TRIGGER_DELTA:
      g_value_set_double (value, cbe_prec->priv->ppr.recfast_x_H0_trigger_delta);
      break;
    case PROP_RECFAST_H_FRAC:
      g_value_set_double (value, cbe_prec->priv->ppr.recfast_H_frac);
      break;
    case PROP_HYREC_ALPHA_INF_FILE:
      g_value_set_string (value, cbe_prec->priv->ppr.hyrec_Alpha_inf_file);
      break;
    case PROP_HYREC_R_INF_FILE:
      g_value_set_string (value, cbe_prec->priv->ppr.hyrec_R_inf_file);
      break;
    case PROP_HYREC_TWO_PHOTON_TABLES_FILE:
      g_value_set_string (value, cbe_prec->priv->ppr.hyrec_two_photon_tables_file);
      break;
    case PROP_REION_Z_START_MAX:
      g_value_set_double (value, cbe_prec->priv->ppr.reionization_z_start_max);
      break;
    case PROP_REION_SAMPLING:
      g_value_set_double (value, cbe_prec->priv->ppr.reionization_sampling);
      break;
    case PROP_REION_OPT_DEPTH_TOL:
      g_value_set_double (value, cbe_prec->priv->ppr.reionization_optical_depth_tol);
      break;
    case PROP_REION_START_FACTOR:
      g_value_set_double (value, cbe_prec->priv->ppr.reionization_start_factor);
      break;
    case PROP_THERMO_RATE_SMOOTHING_RADIUS:
      g_value_set_int (value, cbe_prec->priv->ppr.thermo_rate_smoothing_radius);
      break;
    case PROP_EVOLVER:
      g_value_set_int (value, cbe_prec->priv->ppr.evolver);
      break;
    case PROP_K_MIN_TAU0:
      g_value_set_double (value, cbe_prec->priv->ppr.k_min_tau0);
      break;
    case PROP_K_MAX_TAU0_OVER_L_MAX:
      g_value_set_double (value, cbe_prec->priv->ppr.k_max_tau0_over_l_max);
      break;
    case PROP_K_STEP_SUB:
      g_value_set_double (value, cbe_prec->priv->ppr.k_step_sub);
      break;
    case PROP_K_STEP_SUPER:
      g_value_set_double (value, cbe_prec->priv->ppr.k_step_super);
      break;
    case PROP_K_STEP_TRANSITION:
      g_value_set_double (value, cbe_prec->priv->ppr.k_step_transition);
      break;
    case PROP_K_STEP_SUPER_REDUCTION:
      g_value_set_double (value, cbe_prec->priv->ppr.k_step_super_reduction);
      break;
    case PROP_K_PER_DECADE_FOR_PK:
      g_value_set_double (value, cbe_prec->priv->ppr.k_per_decade_for_pk);
      break;
    case PROP_K_PER_DECADE_FOR_BAO:
      g_value_set_double (value, cbe_prec->priv->ppr.k_per_decade_for_bao);
      break;
    case PROP_K_BAO_CENTER:
      g_value_set_double (value, cbe_prec->priv->ppr.k_bao_center);
      break;
    case PROP_K_BAO_WIDTH:
      g_value_set_double (value, cbe_prec->priv->ppr.k_bao_width);
      break;
    case PROP_START_SMALL_K_AT_TAU_C_OVER_TAU_H:
      g_value_set_double (value, cbe_prec->priv->ppr.start_small_k_at_tau_c_over_tau_h);
      break;
    case PROP_START_LARGE_K_AT_TAU_H_OVER_TAU_K:
      g_value_set_double (value, cbe_prec->priv->ppr.start_large_k_at_tau_h_over_tau_k);
      break;
    case PROP_TIGHT_COUPLING_TRIGGER_TAU_C_OVER_TAU_H:
      g_value_set_double (value, cbe_prec->priv->ppr.tight_coupling_trigger_tau_c_over_tau_h);
      break;
    case PROP_TIGHT_COUPLING_TRIGGER_TAU_C_OVER_TAU_K:
      g_value_set_double (value, cbe_prec->priv->ppr.tight_coupling_trigger_tau_c_over_tau_k);
      break;
    case PROP_START_SOURCES_AT_TAU_C_OVER_TAU_H:
      g_value_set_double (value, cbe_prec->priv->ppr.start_sources_at_tau_c_over_tau_h);
      break;
    case PROP_TIGHT_COUPLING_APPROXIMATION:
      g_value_set_int (value, cbe_prec->priv->ppr.tight_coupling_approximation);
      break;
    case PROP_L_MAX_G:
      g_value_set_int (value, cbe_prec->priv->ppr.l_max_g);
      break;
    case PROP_L_MAX_POL_G:
      g_value_set_int (value, cbe_prec->priv->ppr.l_max_pol_g);
      break;
    case PROP_L_MAX_DR:
      g_value_set_int (value, cbe_prec->priv->ppr.l_max_dr);
      break;
    case PROP_L_MAX_UR:
      g_value_set_int (value, cbe_prec->priv->ppr.l_max_ur);
      break;
    case PROP_L_MAX_NCDM:
      g_value_set_int (value, cbe_prec->priv->ppr.l_max_ncdm);
      break;
    case PROP_L_MAX_G_TEN:
      g_value_set_int (value, cbe_prec->priv->ppr.l_max_g_ten);
      break;
    case PROP_L_MAX_POL_G_TEN:
      g_value_set_int (value, cbe_prec->priv->ppr.l_max_pol_g_ten);
      break;
    case PROP_CURVATURE_INI:
      g_value_set_double (value, cbe_prec->priv->ppr.curvature_ini);
      break;
    case PROP_ENTROPY_INI:
      g_value_set_double (value, cbe_prec->priv->ppr.entropy_ini);
      break;
    case PROP_GW_INI:
      g_value_set_double (value, cbe_prec->priv->ppr.gw_ini);
      break;
    case PROP_PERTURB_INTEGRATION_STEPSIZE:
      g_value_set_double (value, cbe_prec->priv->ppr.perturb_integration_stepsize);
      break;
    case PROP_TOL_TAU_APPROX:
      g_value_set_double (value, cbe_prec->priv->ppr.tol_tau_approx);
      break;
    case PROP_TOL_PERTURB_INTEGRATION:
      g_value_set_double (value, cbe_prec->priv->ppr.tol_perturb_integration);
      break;
    case PROP_PERTURB_SAMPLING_STEPSIZE:
      g_value_set_double (value, cbe_prec->priv->ppr.perturb_sampling_stepsize);
      break;
    case PROP_RADIATION_STREAMING_APPROXIMATION:
      g_value_set_int (value, cbe_prec->priv->ppr.radiation_streaming_approximation);
      break;
    case PROP_RADIATION_STREAMING_TRIGGER_TAU_OVER_TAU_K:
      g_value_set_double (value, cbe_prec->priv->ppr.radiation_streaming_trigger_tau_over_tau_k);
      break;
    case PROP_RADIATION_STREAMING_TRIGGER_TAU_C_OVER_TAU:
      g_value_set_double (value, cbe_prec->priv->ppr.radiation_streaming_trigger_tau_c_over_tau);
      break;
    case PROP_UR_FLUID_APPROXIMATION:
      g_value_set_int (value, cbe_prec->priv->ppr.ur_fluid_approximation);
      break;
    case PROP_UR_FLUID_TRIGGER_TAU_OVER_TAU_K:
      g_value_set_double (value, cbe_prec->priv->ppr.ur_fluid_trigger_tau_over_tau_k);
      break;
    case PROP_NCDM_FLUID_APPROXIMATION:
      g_value_set_int (value, cbe_prec->priv->ppr.ncdm_fluid_approximation);
      break;
    case PROP_NCDM_FLUID_TRIGGER_TAU_OVER_TAU_K:
      g_value_set_double (value, cbe_prec->priv->ppr.ncdm_fluid_trigger_tau_over_tau_k);
      break;
    case PROP_NEGLECT_CMB_SOURCES_BELOW_VISIBILITY:
      g_value_set_double (value, cbe_prec->priv->ppr.neglect_CMB_sources_below_visibility);
      break;
    case PROP_K_PER_DECADE_PRIMORDIAL:
      g_value_set_double (value, cbe_prec->priv->ppr.k_per_decade_primordial);
      break;
    case PROP_PRIMORDIAL_INFLATION_RATIO_MIN:
      g_value_set_double (value, cbe_prec->priv->ppr.primordial_inflation_ratio_min);
      break;
    case PROP_PRIMORDIAL_INFLATION_RATIO_MAX:
      g_value_set_double (value, cbe_prec->priv->ppr.primordial_inflation_ratio_max);
      break;
    case PROP_PRIMORDIAL_INFLATION_PHI_INI_MAXIT:
      g_value_set_int (value, cbe_prec->priv->ppr.primordial_inflation_phi_ini_maxit);
      break;
    case PROP_PRIMORDIAL_INFLATION_PT_STEPSIZE:
      g_value_set_double (value, cbe_prec->priv->ppr.primordial_inflation_pt_stepsize);
      break;
    case PROP_PRIMORDIAL_INFLATION_BG_STEPSIZE:
      g_value_set_double (value, cbe_prec->priv->ppr.primordial_inflation_bg_stepsize);
      break;
    case PROP_PRIMORDIAL_INFLATION_TOL_INTEGRATION:
      g_value_set_double (value, cbe_prec->priv->ppr.primordial_inflation_tol_integration);
      break;
    case PROP_PRIMORDIAL_INFLATION_ATTRACTOR_PRECISION_PIVOT:
      g_value_set_double (value, cbe_prec->priv->ppr.primordial_inflation_attractor_precision_pivot);
      break;
    case PROP_PRIMORDIAL_INFLATION_ATTRACTOR_PRECISION_INITIAL:
      g_value_set_double (value, cbe_prec->priv->ppr.primordial_inflation_attractor_precision_initial);
      break;
    case PROP_PRIMORDIAL_INFLATION_ATTRACTOR_MAXIT:
      g_value_set_int (value, cbe_prec->priv->ppr.primordial_inflation_attractor_maxit);
      break;
    case PROP_PRIMORDIAL_INFLATION_TOL_CURVATURE:
      g_value_set_double (value, cbe_prec->priv->ppr.primordial_inflation_tol_curvature);
      break;
    case PROP_PRIMORDIAL_INFLATION_AH_INI_TARGET:
      g_value_set_double (value, cbe_prec->priv->ppr.primordial_inflation_aH_ini_target);
      break;
    case PROP_PRIMORDIAL_INFLATION_END_DPHI:
      g_value_set_double (value, cbe_prec->priv->ppr.primordial_inflation_end_dphi);
      break;
    case PROP_PRIMORDIAL_INFLATION_END_LOGSTEP:
      g_value_set_double (value, cbe_prec->priv->ppr.primordial_inflation_end_logstep);
      break;
    case PROP_PRIMORDIAL_INFLATION_SMALL_EPSILON:
      g_value_set_double (value, cbe_prec->priv->ppr.primordial_inflation_small_epsilon);
      break;
    case PROP_PRIMORDIAL_INFLATION_SMALL_EPSILON_TOL:
      g_value_set_double (value, cbe_prec->priv->ppr.primordial_inflation_small_epsilon_tol);
      break;
    case PROP_PRIMORDIAL_INFLATION_EXTRA_EFOLDS:
      g_value_set_double (value, cbe_prec->priv->ppr.primordial_inflation_extra_efolds);
      break;
    case PROP_L_LOGSTEP:
      g_value_set_double (value, cbe_prec->priv->ppr.l_logstep);
      break;
    case PROP_L_LINSTEP:
      g_value_set_int (value, cbe_prec->priv->ppr.l_linstep);
      break;
    case PROP_HYPER_X_MIN:
      g_value_set_double (value, cbe_prec->priv->ppr.hyper_x_min);
      break;
    case PROP_HYPER_SAMPLING_FLAT:
      g_value_set_double (value, cbe_prec->priv->ppr.hyper_sampling_flat);
      break;
    case PROP_HYPER_SAMPLING_CURVED_LOW_NU:
      g_value_set_double (value, cbe_prec->priv->ppr.hyper_sampling_curved_low_nu);
      break;
    case PROP_HYPER_SAMPLING_CURVED_HIGH_NU:
      g_value_set_double (value, cbe_prec->priv->ppr.hyper_sampling_curved_high_nu);
      break;
    case PROP_HYPER_NU_SAMPLING_STEP:
      g_value_set_double (value, cbe_prec->priv->ppr.hyper_nu_sampling_step);
      break;
    case PROP_HYPER_PHI_MIN_ABS:
      g_value_set_double (value, cbe_prec->priv->ppr.hyper_phi_min_abs);
      break;
    case PROP_HYPER_X_TOL:
      g_value_set_double (value, cbe_prec->priv->ppr.hyper_x_tol);
      break;
    case PROP_HYPER_FLAT_APPROXIMATION_NU:
      g_value_set_double (value, cbe_prec->priv->ppr.hyper_flat_approximation_nu);
      break;
    case PROP_Q_LINSTEP:
      g_value_set_double (value, cbe_prec->priv->ppr.q_linstep);
      break;
    case PROP_Q_LOGSTEP_SPLINE:
      g_value_set_double (value, cbe_prec->priv->ppr.q_logstep_spline);
      break;
    case PROP_Q_LOGSTEP_OPEN:
      g_value_set_double (value, cbe_prec->priv->ppr.q_logstep_open);
      break;
    case PROP_Q_LOGSTEP_TRAPZD:
      g_value_set_double (value, cbe_prec->priv->ppr.q_logstep_trapzd);
      break;
    case PROP_Q_NUMSTEP_TRANSITION:
      g_value_set_double (value, cbe_prec->priv->ppr.q_numstep_transition);
      break;
    case PROP_TRANSFER_NEGLECT_DELTA_K_S_T0:
      g_value_set_double (value, cbe_prec->priv->ppr.transfer_neglect_delta_k_S_t0);
      break;
    case PROP_TRANSFER_NEGLECT_DELTA_K_S_T1:
      g_value_set_double (value, cbe_prec->priv->ppr.transfer_neglect_delta_k_S_t1);
      break;
    case PROP_TRANSFER_NEGLECT_DELTA_K_S_T2:
      g_value_set_double (value, cbe_prec->priv->ppr.transfer_neglect_delta_k_S_t2);
      break;
    case PROP_TRANSFER_NEGLECT_DELTA_K_S_E:
      g_value_set_double (value, cbe_prec->priv->ppr.transfer_neglect_delta_k_S_e);
      break;
    case PROP_TRANSFER_NEGLECT_DELTA_K_V_T1:
      g_value_set_double (value, cbe_prec->priv->ppr.transfer_neglect_delta_k_V_t1);
      break;
    case PROP_TRANSFER_NEGLECT_DELTA_K_V_T2:
      g_value_set_double (value, cbe_prec->priv->ppr.transfer_neglect_delta_k_V_t2);
      break;
    case PROP_TRANSFER_NEGLECT_DELTA_K_V_E:
      g_value_set_double (value, cbe_prec->priv->ppr.transfer_neglect_delta_k_V_e);
      break;
    case PROP_TRANSFER_NEGLECT_DELTA_K_V_B:
      g_value_set_double (value, cbe_prec->priv->ppr.transfer_neglect_delta_k_V_b);
      break;
    case PROP_TRANSFER_NEGLECT_DELTA_K_T_T2:
      g_value_set_double (value, cbe_prec->priv->ppr.transfer_neglect_delta_k_T_t2);
      break;
    case PROP_TRANSFER_NEGLECT_DELTA_K_T_E:
      g_value_set_double (value, cbe_prec->priv->ppr.transfer_neglect_delta_k_T_e);
      break;
    case PROP_TRANSFER_NEGLECT_DELTA_K_T_B:
      g_value_set_double (value, cbe_prec->priv->ppr.transfer_neglect_delta_k_T_b);
      break;
    case PROP_TRANSFER_NEGLECT_LATE_SOURCE:
      g_value_set_double (value, cbe_prec->priv->ppr.transfer_neglect_late_source);
      break;
    case PROP_L_SWITCH_LIMBER:
      g_value_set_double (value, cbe_prec->priv->ppr.l_switch_limber);
      break;
    case PROP_L_SWITCH_LIMBER_FOR_NC_LOCAL_OVER_Z:
      g_value_set_double (value, cbe_prec->priv->ppr.l_switch_limber_for_nc_local_over_z);
      break;
    case PROP_L_SWITCH_LIMBER_FOR_NC_LOS_OVER_Z:
      g_value_set_double (value, cbe_prec->priv->ppr.l_switch_limber_for_nc_los_over_z);
      break;
    case PROP_SELECTION_CUT_AT_SIGMA:
      g_value_set_double (value, cbe_prec->priv->ppr.selection_cut_at_sigma);
      break;
    case PROP_SELECTION_SAMPLING:
      g_value_set_double (value, cbe_prec->priv->ppr.selection_sampling);
      break;
    case PROP_SELECTION_SAMPLING_BESSEL:
      g_value_set_double (value, cbe_prec->priv->ppr.selection_sampling_bessel);
      break;
    case PROP_SELECTION_SAMPLING_BESSEL_LOS:
      g_value_set_double (value, cbe_prec->priv->ppr.selection_sampling_bessel_los);
      break;
    case PROP_SELECTION_TOPHAT_EDGE:
      g_value_set_double (value, cbe_prec->priv->ppr.selection_tophat_edge);
      break;
    case PROP_HALOFIT_MIN_K_NONLINEAR:
      g_value_set_double (value, cbe_prec->priv->ppr.halofit_min_k_nonlinear);
      break;
    case PROP_HALOFIT_MIN_K_MAX:
      g_value_set_double (value, cbe_prec->priv->ppr.halofit_min_k_max);
      break;
    case PROP_HALOFIT_K_PER_DECADE:
      g_value_set_double (value, cbe_prec->priv->ppr.halofit_k_per_decade);
      break;
    case PROP_HALOFIT_SIGMA_PRECISION:
      g_value_set_double (value, cbe_prec->priv->ppr.halofit_sigma_precision);
      break;
    case PROP_HALOFIT_SIGMA_TOL:
      g_value_set_double (value, cbe_prec->priv->ppr.halofit_tol_sigma);
      break;
    case PROP_PK_EQ_Z_MAX:
      g_value_set_double (value, cbe_prec->priv->ppr.pk_eq_z_max);
      break;
    case PROP_PK_EQ_TOL:
      g_value_set_double (value, cbe_prec->priv->ppr.pk_eq_tol);
      break;
    case PROP_ACCURATE_LENSING:
      g_value_set_int (value, cbe_prec->priv->ppr.accurate_lensing);
      break;
    case PROP_NUM_MU_MINUS_LMAX:
      g_value_set_int (value, cbe_prec->priv->ppr.num_mu_minus_lmax);
      break;
    case PROP_DELTA_L_MAX:
      g_value_set_int (value, cbe_prec->priv->ppr.delta_l_max);
      break;
    case PROP_SMALLEST_ALLOWED_VARIATION:
      g_value_set_double (value, cbe_prec->priv->ppr.smallest_allowed_variation);
      break;
    case PROP_TOL_GAUSS_LEGENDRE:
      g_value_set_double (value, cbe_prec->priv->ppr.tol_gauss_legendre);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_cbe_precision_finalize (GObject* object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_cbe_precision_parent_class)->finalize (object);
}

static void
nc_cbe_precision_class_init (NcCBEPrecisionClass* klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &nc_cbe_precision_set_property;
  object_class->get_property = &nc_cbe_precision_get_property;
  object_class->finalize     = &nc_cbe_precision_finalize;

  /*
   * Background related parameters.
   */
  g_object_class_install_property (object_class,
                                   PROP_A_INI_A_0,
                                   g_param_spec_double ("a-ini-over-a-today-default",
                                                        NULL,
                                                        "default initial value of scale factor in background integration, in units of scale factor today",
                                                        0.0, G_MAXDOUBLE, 1.0e-14,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_BACK_INT_STEP,
                                   g_param_spec_double ("back-integration-stepsize",
                                                        NULL,
                                                        "default step d tau in background integration, in units of conformal Hubble time ($d \\tau$ = back_integration_stepsize / aH )",
                                                        0.0, G_MAXDOUBLE, 7.0e-3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_BACK_TOL,
                                   g_param_spec_double ("tol-background-integration",
                                                        NULL,
                                                        "parameter controlling precision of background integration",
                                                        0.0, G_MAXDOUBLE, 1.0e-2,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_INITIAL_OMEGA_R_TOL,
                                   g_param_spec_double ("tol-initial-Omega-r",
                                                        NULL,
                                                        "parameter controlling how deep inside radiation domination must the initial time be chosen",
                                                        0.0, G_MAXDOUBLE, 1.0e-4,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_M_NCDM_TOL,
                                   g_param_spec_double ("tol-M-ncdm",
                                                        NULL,
                                                        "parameter controlling relative precision of ncdm mass for given ncdm current density",
                                                        0.0, G_MAXDOUBLE, 1.0e-7,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NCDM_TOL,
                                   g_param_spec_double ("tol-ncdm",
                                                        NULL,
                                                        "parameter controlling relative precision of integrals over ncdm phase-space distribution during perturbation calculation",
                                                        0.0, G_MAXDOUBLE, 1.0e-3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NCDM_SYNCHRONOUS_TOL,
                                   g_param_spec_double ("tol-ncdm-synchronous",
                                                        NULL,
                                                        "parameter controlling relative precision of integrals over ncdm phase-space distribution during perturbation calculation - synchronous",
                                                        0.0, G_MAXDOUBLE, 1.0e-3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NCDM_NEWTONIAN_TOL,
                                   g_param_spec_double ("tol-ncdm-newtonian",
                                                        NULL,
                                                        "parameter controlling relative precision of integrals over ncdm phase-space distribution during perturbation calculation - newtonian",
                                                        0.0, G_MAXDOUBLE, 1.0e-5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NCDM_BG_TOL,
                                   g_param_spec_double ("tol-ncdm-bg",
                                                        NULL,
                                                        "parameter controlling relative precision of integrals over ncdm phase-space distribution during background evolution",
                                                        0.0, G_MAXDOUBLE, 1.0e-5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NCDM_INITIAL_W_TOL,
                                   g_param_spec_double ("tol-ncdm-initial-w",
                                                        NULL,
                                                        "parameter controlling how relativistic must non-cold relics be at initial time",
                                                        0.0, G_MAXDOUBLE, 1.0e-3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SAFE_PHI_SCF,
                                   g_param_spec_double ("safe-phi-scf",
                                                        NULL,
                                                        "parameter controlling the initial scalar field in background functions",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0, /* Undefined in CLASS: FIXME */
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TAU_EQ_TOL,
                                   g_param_spec_double ("tol-tau-eq",
                                                        NULL,
                                                        "parameter controlling precision with which tau_eq (conformal time at radiation/matter equality) is found (units: Mpc)",
                                                        0.0, G_MAXDOUBLE, 1.0e-6,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /*
   * Thermodynamics related parameters
   */
  {
    gchar *sBBN_file = ncm_cfg_get_data_filename ("class_data" G_DIR_SEPARATOR_S "bbn" G_DIR_SEPARATOR_S "sBBN_2017.dat", TRUE);
    g_object_class_install_property (object_class,
                                     PROP_SBBN_FILE,
                                     g_param_spec_string ("sBBN-file",
                                                          NULL,
                                                          "SBBN filename",
                                                          sBBN_file,
                                                          G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
    g_free (sBBN_file);
  }
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_Z_INI,
                                   g_param_spec_double ("recfast-z-initial",
                                                        NULL,
                                                        "initial redshift in recfast",
                                                        0.0, G_MAXDOUBLE, 1.0e4,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_NZ0,
                                   g_param_spec_int ("recfast-Nz0",
                                                     NULL,
                                                     "number of integration steps",
                                                     0, G_MAXINT, 20000,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_THERMO_INTEGRATION_TOL,
                                   g_param_spec_double ("tol-thermo-integration",
                                                        NULL,
                                                        "precision of each integration step",
                                                        0.0, G_MAXDOUBLE, 1.0e-2,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_HE_SWITCH,
                                   g_param_spec_int ("recfast-Heswitch",
                                                     NULL,
                                                     "Recfast He switch",
                                                     0, G_MAXINT, 6,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_FUDGE_HE,
                                   g_param_spec_double ("recfast-fudge-He",
                                                        NULL,
                                                        "Recfast fudge He",
                                                        0.0, G_MAXDOUBLE, 0.86,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_H_SWITCH,
                                   g_param_spec_int ("recfast-Hswitch",
                                                     NULL,
                                                     "recfast 1.5 switching parameter",
                                                     0, G_MAXINT, _TRUE_,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_FUDGE_H,
                                   g_param_spec_double ("recfast-fudge-H",
                                                        NULL,
                                                        "H fudge factor when recfast_Hswitch set to false (v1.4 fudging)",
                                                        0.0, G_MAXDOUBLE, 1.14,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_DELTA_FUDGE_H,
                                   g_param_spec_double ("recfast-delta-fudge-H",
                                                        NULL,
                                                        "correction to H fudge factor in v1.5",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, -0.015,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_A_GAUSS_1,
                                   g_param_spec_double ("recfast-AGauss1",
                                                        NULL,
                                                        "Amplitude of 1st Gaussian",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, -0.14,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_A_GAUSS_2,
                                   g_param_spec_double ("recfast-AGauss2",
                                                        NULL,
                                                        "Amplitude of 2st Gaussian",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.079,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_Z_GAUSS_1,
                                   g_param_spec_double ("recfast-zGauss1",
                                                        NULL,
                                                        "ln(1+z) of 1st Gaussian",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 7.28,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_Z_GAUSS_2,
                                   g_param_spec_double ("recfast-zGauss2",
                                                        NULL,
                                                        "ln(1+z) of 2st Gaussian",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 6.73,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_W_GAUSS_1,
                                   g_param_spec_double ("recfast-wGauss1",
                                                        NULL,
                                                        "Width of 2st Gaussian",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.18,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_W_GAUSS_2,
                                   g_param_spec_double ("recfast-wGauss2",
                                                        NULL,
                                                        "Width of 2st Gaussian",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.33,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_Z_HE_1,
                                   g_param_spec_double ("recfast-z-He-1",
                                                        NULL,
                                                        "down to which redshift Helium fully ionized",
                                                        0.0, G_MAXDOUBLE, 8000.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_DELTA_Z_HE_1,
                                   g_param_spec_double ("recfast-delta-z-He-1",
                                                        NULL,
                                                        "z range over which transition is smoothed",
                                                        0.0, G_MAXDOUBLE, 50.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_Z_HE_2,
                                   g_param_spec_double ("recfast-z-He-2",
                                                        NULL,
                                                        "down to which redshift first Helium recombination not complete",
                                                        0.0, G_MAXDOUBLE, 5000.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_DELTA_Z_HE_2,
                                   g_param_spec_double ("recfast-delta-z-He-2",
                                                        NULL,
                                                        "z range over which transition is smoothed",
                                                        0.0, G_MAXDOUBLE, 100.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_Z_HE_3,
                                   g_param_spec_double ("recfast-z-He-3",
                                                        NULL,
                                                        "down to which redshift Helium singly ionized",
                                                        0.0, G_MAXDOUBLE, 3500.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_DELTA_Z_HE_3,
                                   g_param_spec_double ("recfast-delta-z-He-3",
                                                        NULL,
                                                        "z range over which transition is smoothed",
                                                        0.0, G_MAXDOUBLE, 50.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_X_HE0_TRIGGER,
                                   g_param_spec_double ("recfast-x-He0-trigger",
                                                        NULL,
                                                        "value below which recfast uses the full equation for Helium",
                                                        0.0, G_MAXDOUBLE, 0.995,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_X_HE0_TRIGGER2,
                                   g_param_spec_double ("recfast-x-He0-trigger2",
                                                        NULL,
                                                        "a second threshold used in derivative routine",
                                                        0.0, G_MAXDOUBLE, 0.995,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_X_HE0_TRIGGER_DELTA,
                                   g_param_spec_double ("recfast-x-He0-trigger-delta",
                                                        NULL,
                                                        "x_He range over which transition is smoothed",
                                                        0.0, G_MAXDOUBLE, 0.05,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_X_H0_TRIGGER,
                                   g_param_spec_double ("recfast-x-H0-trigger",
                                                        NULL,
                                                        "value below which recfast uses the full equation for Hydrogen",
                                                        0.0, G_MAXDOUBLE, 0.995,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_X_H0_TRIGGER2,
                                   g_param_spec_double ("recfast-x-H0-trigger2",
                                                        NULL,
                                                        "a second threshold used in derivative routine",
                                                        0.0, G_MAXDOUBLE, 0.995,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_X_H0_TRIGGER_DELTA,
                                   g_param_spec_double ("recfast-x-H0-trigger-delta",
                                                        NULL,
                                                        "x_H range over which transition is smoothed",
                                                        0.0, G_MAXDOUBLE, 0.05,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_H_FRAC,
                                   g_param_spec_double ("recfast-H-frac",
                                                        NULL,
                                                        "governs time at which full equation of evolution for Tmat is used",
                                                        0.0, G_MAXDOUBLE, 1.0e-3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  {
    gchar *hyrec_Alpha_inf_file = ncm_cfg_get_data_filename ("class_data" G_DIR_SEPARATOR_S "hyrec" G_DIR_SEPARATOR_S "Alpha_inf.dat", TRUE);
    g_object_class_install_property (object_class,
                                     PROP_HYREC_ALPHA_INF_FILE,
                                     g_param_spec_string ("hyrec-Alpha-inf-file",
                                                          NULL,
                                                          "Hyrec Alpha inf file",
                                                          hyrec_Alpha_inf_file,
                                                          G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
    g_free (hyrec_Alpha_inf_file);
  }
  {
    gchar *hyrec_R_inf_file = ncm_cfg_get_data_filename ("class_data" G_DIR_SEPARATOR_S "hyrec" G_DIR_SEPARATOR_S "R_inf.dat", TRUE);
    g_object_class_install_property (object_class,
                                     PROP_HYREC_R_INF_FILE,
                                     g_param_spec_string ("hyrec-R-inf-file",
                                                          NULL,
                                                          "Hyrec R inf file",
                                                          hyrec_R_inf_file,
                                                          G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
    g_free (hyrec_R_inf_file);
  }
  {
    gchar *hyrec_two_photon_tables_file = ncm_cfg_get_data_filename ("class_data" G_DIR_SEPARATOR_S "hyrec" G_DIR_SEPARATOR_S "two_photon_tables.dat", TRUE);
    g_object_class_install_property (object_class,
                                     PROP_HYREC_TWO_PHOTON_TABLES_FILE,
                                     g_param_spec_string ("hyrec-two-photon-tables-file",
                                                          NULL,
                                                          "Hyrec two photon tables file",
                                                          hyrec_two_photon_tables_file,
                                                          G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
    g_free (hyrec_two_photon_tables_file);
  }
  /* for reionization */
  g_object_class_install_property (object_class,
                                   PROP_REION_Z_START_MAX,
                                   g_param_spec_double ("reionization-z-start-max",
                                                        NULL,
                                                        "maximum redshift at which reionization should start. If not, return an error",
                                                        0.0, G_MAXDOUBLE, 50.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_REION_SAMPLING,
                                   g_param_spec_double ("reionization-sampling",
                                                        NULL,
                                                        "control stepsize in z during reionization",
                                                        0.0, G_MAXDOUBLE, 5.0e-2,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_REION_OPT_DEPTH_TOL,
                                   g_param_spec_double ("reionization-optical-depth-tol",
                                                        NULL,
                                                        "fractional error on optical_depth",
                                                        0.0, G_MAXDOUBLE, 1.0e-4,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_REION_START_FACTOR,
                                   g_param_spec_double ("reionization-start-factor",
                                                        NULL,
                                                        "parameter for CAMB-like parametrization",
                                                        0.0, G_MAXDOUBLE, 8.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /* general */
  g_object_class_install_property (object_class,
                                   PROP_THERMO_RATE_SMOOTHING_RADIUS,
                                   g_param_spec_int ("thermo-rate-smoothing-radius",
                                                     NULL,
                                                     "plays a minor (almost aesthetic) role in the definition of the variation rate of thermodynamical quantities",
                                                     0, G_MAXINT, 50,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /*
   * Perturbations parameters
   */
  g_object_class_install_property (object_class,
                                   PROP_EVOLVER,
                                   g_param_spec_int ("evolver",
                                                     NULL,
                                                     "which type of evolver for integrating perturbations (Runge-Kutta? Stiff?...)",
                                                     0, G_MAXINT, ndf15,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_K_MIN_TAU0,
                                   g_param_spec_double ("k-min-tau0",
                                                        NULL,
                                                        "number defining k_min for the computation of Cl's and P(k)'s (dimensionless): (k_min tau_0), usually chosen much smaller than one",
                                                        0.0, G_MAXDOUBLE, 0.1,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_K_MAX_TAU0_OVER_L_MAX,
                                   g_param_spec_double ("k-max-tau0-over-l-max",
                                                        NULL,
                                                        "number defining k_max for the computation of Cl's (dimensionless): (k_max tau_0)/l_max, usually chosen around two (very relevant for accuracy of lensed ClTT at highest l's)",
                                                        0.0, G_MAXDOUBLE, 2.4,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_K_STEP_SUB,
                                   g_param_spec_double ("k-step-sub",
                                                        NULL,
                                                        "step in k space, in units of one period of acoustic oscillation at decoupling, for scales inside sound horizon at decoupling",
                                                        0.0, G_MAXDOUBLE, 0.05,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_K_STEP_SUPER,
                                   g_param_spec_double ("k-step-super",
                                                        NULL,
                                                        "step in k space, in units of one period of acoustic oscillation at decoupling, for scales above sound horizon at decoupling",
                                                        0.0, G_MAXDOUBLE, 0.002,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_K_STEP_TRANSITION,
                                   g_param_spec_double ("k-step-transition",
                                                        NULL,
                                                        "dimensionless number regulating the transition from 'sub' steps to 'super' steps. Decrease for more precision",
                                                        0.0, G_MAXDOUBLE, 0.2,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_K_STEP_SUPER_REDUCTION,
                                   g_param_spec_double ("k-step-super-reduction",
                                                        NULL,
                                                        "the step k_step_super is reduced by this amount in the k-->0 limit (below scale of Hubble and/or curvature radius)",
                                                        0.0, G_MAXDOUBLE, 0.1,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_K_PER_DECADE_FOR_PK,
                                   g_param_spec_double ("k-per-decade-for-pk",
                                                        NULL,
                                                        "if values needed between kmax inferred from k_oscillations and k_kmax_for_pk, this gives the number of k per decade outside the BAO region",
                                                        0.0, G_MAXDOUBLE, 10.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_K_PER_DECADE_FOR_BAO,
                                   g_param_spec_double ("k-per-decade-for-bao",
                                                        NULL,
                                                        "if values needed between kmax inferred from k_oscillations and k_kmax_for_pk, this gives the number of k per decade inside the BAO region (for finer sampling)",
                                                        0.0, G_MAXDOUBLE, 70.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_K_BAO_CENTER,
                                   g_param_spec_double ("k-bao-center",
                                                        NULL,
                                                        "in ln(k) space, the central value of the BAO region where sampling is finer is defined as k_rec times this number (recommended: 3, i.e. finest sampling near 3rd BAO peak)",
                                                        0.0, G_MAXDOUBLE, 3.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_K_BAO_WIDTH,
                                   g_param_spec_double ("k-bao-width",
                                                        NULL,
                                                        "in ln(k) space, width of the BAO region where sampling is finer: this number gives roughly the number of BAO oscillations well resolved on both sides of the central value (recommended: 4, i.e. finest sampling from before first up to 3+4=7th peak)",
                                                        0.0, G_MAXDOUBLE, 4.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_START_SMALL_K_AT_TAU_C_OVER_TAU_H,
                                   g_param_spec_double ("start-small-k-at-tau-c-over-tau-h",
                                                        NULL,
                                                        "largest wavelengths start being sampled when universe is sufficiently opaque. This is quantified in terms of the ratio of thermo to hubble time scales, $\tau_c/\tau_H$. Start when start_largek_at_tau_c_over_tau_h equals this ratio. Decrease this value to start integrating the wavenumbers earlier in time.",
                                                        0.0, G_MAXDOUBLE, 0.0015,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_START_LARGE_K_AT_TAU_H_OVER_TAU_K,
                                   g_param_spec_double ("start-large-k-at-tau-h-over-tau-k",
                                                        NULL,
                                                        "largest wavelengths start being sampled when mode is sufficiently outside Hibble scale. This is quantified in terms of the ratio of hubble time scale to wavenumber time scale, $\tau_h/\tau_k$ wich is roughly equal to (k*tau). Start when this ratio equals start_large_k_at_tau_k_over_tau_h. Decrease this value to start integrating the wavenumbers earlier in time.",
                                                        0.0, G_MAXDOUBLE, 0.07,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TIGHT_COUPLING_TRIGGER_TAU_C_OVER_TAU_H,
                                   g_param_spec_double ("tight-coupling-trigger-tau-c-over-tau-h",
                                                        NULL,
                                                        "when to switch off tight-coupling approximation: first condition: $\\tau_c/\\tau_H$ > tight_coupling_trigger_tau_c_over_tau_h. Decrease this value to switch off earlier in time.  If this number is larger than start_sources_at_tau_c_over_tau_h, the code returns an error, because the source computation requires tight-coupling to be switched off.",
                                                        0.0, G_MAXDOUBLE, 0.015,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TIGHT_COUPLING_TRIGGER_TAU_C_OVER_TAU_K,
                                   g_param_spec_double ("tight-coupling-trigger-tau-c-over-tau-k",
                                                        NULL,
                                                        "when to switch off tight-coupling approximation: second condition: $\\tau_c/\\tau_k \\equiv k \\tau_c$ < tight_coupling_trigger_tau_c_over_tau_k. Decrease this value to switch off earlier in time.",
                                                        0.0, G_MAXDOUBLE, 0.01,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_START_SOURCES_AT_TAU_C_OVER_TAU_H,
                                   g_param_spec_double ("start-sources-at-tau-c-over-tau-h",
                                                        NULL,
                                                        "sources start being sampled when universe is sufficiently opaque. This is quantified in terms of the ratio of thermo to hubble time scales, $\\tau_c/\\tau_H$. Start when start_sources_at_tau_c_over_tau_h equals this ratio. Decrease this value to start sampling the sources earlier in time.",
                                                        0.0, G_MAXDOUBLE, 0.008,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TIGHT_COUPLING_APPROXIMATION,
                                   g_param_spec_int ("tight-coupling-approximation",
                                                     NULL,
                                                     "Tight coupling approximation scheme",
                                                     0, G_MAXINT, (gint)compromise_CLASS,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_L_MAX_G,
                                   g_param_spec_int ("l-max-g",
                                                     NULL,
                                                     "number of momenta in Boltzmann hierarchy for photon temperature (scalar)",
                                                     4, G_MAXINT, 12,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_L_MAX_POL_G,
                                   g_param_spec_int ("l-max-pol-g",
                                                     NULL,
                                                     "number of momenta in Boltzmann hierarchy for photon polarisation (scalar)",
                                                     4, G_MAXINT, 10,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_L_MAX_DR,
                                   g_param_spec_int ("l-max-dr",
                                                     NULL,
                                                     "number of momenta in Boltzmann hierarchy for decay radiation",
                                                     4, G_MAXINT, 17,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_L_MAX_UR,
                                   g_param_spec_int ("l-max-ur",
                                                     NULL,
                                                     "number of momenta in Boltzmann hierarchy for relativistic neutrino/relics (scalar)",
                                                     4, G_MAXINT, 17,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_L_MAX_NCDM,
                                   g_param_spec_int ("l-max-ncdm",
                                                     NULL,
                                                     "number of momenta in Boltzmann hierarchy for relativistic neutrino/relics (scalar)",
                                                     4, G_MAXINT, 17,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_L_MAX_G_TEN,
                                   g_param_spec_int ("l-max-g-ten",
                                                     NULL,
                                                     "number of momenta in Boltzmann hierarchy for photon temperature (tensor)",
                                                     4, G_MAXINT, 5,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_L_MAX_POL_G_TEN,
                                   g_param_spec_int ("l-max-pol-g-ten",
                                                     NULL,
                                                     "number of momenta in Boltzmann hierarchy for photon polarisation (tensor)",
                                                     4, G_MAXINT, 5,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_CURVATURE_INI,
                                   g_param_spec_double ("curvature-ini",
                                                        NULL,
                                                        "initial curvature; used to fix adiabatic initial conditions; must remain fixed to one as long as the primordial adiabatic spectrum stands for the curvature power spectrum",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ENTROPY_INI,
                                   g_param_spec_double ("entropy-ini",
                                                        NULL,
                                                        "initial entropy; used to fix isocurvature initial conditions; must remain fixed to one as long as the primordial isocurvature spectrum stands for an entropy power spectrum",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_GW_INI,
                                   g_param_spec_double ("gw-ini",
                                                        NULL,
                                                        "initial condition for tensor metric perturbation h",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PERTURB_INTEGRATION_STEPSIZE,
                                   g_param_spec_double ("perturb-integration-stepsize",
                                                        NULL,
                                                        "default step $d \\tau$ in perturbation integration, in units of the timescale involved in the equations (usally, the min of $1/k$, $1/aH$, $1/\\dot{\\kappa}$)",
                                                        0.0, G_MAXDOUBLE, 0.5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TOL_TAU_APPROX,
                                   g_param_spec_double ("tol-tau-approx",
                                                        NULL,
                                                        "precision with which the code should determine (by bisection) the times at which sources start being sampled, and at which approximations must be switched on/off (units of Mpc)",
                                                        0.0, G_MAXDOUBLE, 1.0e-10,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TOL_PERTURB_INTEGRATION,
                                   g_param_spec_double ("tol-perturb-integration",
                                                        NULL,
                                                        "control parameter for the precision of the perturbation integration",
                                                        0.0, G_MAXDOUBLE, 1.0e-5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PERTURB_SAMPLING_STEPSIZE,
                                   g_param_spec_double ("perturb-sampling-stepsize",
                                                        NULL,
                                                        "default step $d \\tau$ for sampling the source function, in units of the timescale involved in the sources: $(\\dot{\\kappa}- \\ddot{\\kappa}/\\dot{\\kappa})^{-1}$",
                                                        0.0, G_MAXDOUBLE, 0.10,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RADIATION_STREAMING_APPROXIMATION,
                                   g_param_spec_int ("radiation-streaming-approximation",
                                                     NULL,
                                                     "method for switching off photon perturbations",
                                                     0, G_MAXINT, rsa_MD_with_reio,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RADIATION_STREAMING_TRIGGER_TAU_OVER_TAU_K,
                                   g_param_spec_double ("radiation-streaming-trigger-tau-over-tau-k",
                                                        NULL,
                                                        "when to switch off photon perturbations, ie when to switch on photon free-streaming approximation (keep density and thtau, set shear and higher momenta to zero): first condition: $k \tau$ > radiation_streaming_trigger_tau_h_over_tau_k",
                                                        0.0, G_MAXDOUBLE, 45.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RADIATION_STREAMING_TRIGGER_TAU_C_OVER_TAU,
                                   g_param_spec_double ("radiation-streaming-trigger-tau-c-over-tau",
                                                        NULL,
                                                        "when to switch off photon perturbations, ie when to switch on photon free-streaming approximation (keep density and theta, set shear and higher momenta to zero): second condition:",
                                                        0.0, G_MAXDOUBLE, 5.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_UR_FLUID_APPROXIMATION,
                                   g_param_spec_int ("ur-fluid-approximation",
                                                     NULL,
                                                     "UR fluid approximation scheme",
                                                     0, G_MAXINT, ufa_CLASS,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_UR_FLUID_TRIGGER_TAU_OVER_TAU_K,
                                   g_param_spec_double ("ur-fluid-trigger-tau-over-tau-k",
                                                        NULL,
                                                        "when to switch off ur (massless neutrinos / ultra-relativistic relics) fluid approximation",
                                                        0.0, G_MAXDOUBLE, 30.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NCDM_FLUID_APPROXIMATION,
                                   g_param_spec_int ("ncdm-fluid-approximation",
                                                     NULL,
                                                     "NCDM fluid approximation scheme",
                                                     0, G_MAXINT, ncdmfa_CLASS,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NCDM_FLUID_TRIGGER_TAU_OVER_TAU_K,
                                   g_param_spec_double ("ncdm-fluid-trigger-tau-over-tau-k",
                                                        NULL,
                                                        "when to switch off ncdm (massive neutrinos / non-cold relics) fluid approximation",
                                                        0.0, G_MAXDOUBLE, 31.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NEGLECT_CMB_SOURCES_BELOW_VISIBILITY,
                                   g_param_spec_double ("neglect-CMB-sources-below-visibility",
                                                        NULL,
                                                        "neglect CMB sources below visibility",
                                                        0.0, G_MAXDOUBLE, 1.0e-3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /*
   * Primordial spectra parameters
   */
  g_object_class_install_property (object_class,
                                   PROP_K_PER_DECADE_PRIMORDIAL,
                                   g_param_spec_double ("k-per-decade-primordial",
                                                        NULL,
                                                        "logarithmic sampling for primordial spectra (number of points per decade in k space)",
                                                        0.0, G_MAXDOUBLE, 10.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PRIMORDIAL_INFLATION_RATIO_MIN,
                                   g_param_spec_double ("primordial-inflation-ratio-min",
                                                        NULL,
                                                        "primordial inflation ratio min",
                                                        0.0, G_MAXDOUBLE, 100.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PRIMORDIAL_INFLATION_RATIO_MAX,
                                   g_param_spec_double ("primordial-inflation-ratio-max",
                                                        NULL,
                                                        "primordial inflation ratio max",
                                                        0.0, G_MAXDOUBLE, 1.0 / 50.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PRIMORDIAL_INFLATION_PHI_INI_MAXIT,
                                   g_param_spec_int ("primordial-inflation-phi-ini-maxit",
                                                     NULL,
                                                     "primordial inflation phi ini maxit",
                                                     0, G_MAXINT, 10000,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PRIMORDIAL_INFLATION_PT_STEPSIZE,
                                   g_param_spec_double ("primordial-inflation-pt-stepsize",
                                                        NULL,
                                                        "primordial inflation pt stepsize",
                                                        0.0, G_MAXDOUBLE, 0.01,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PRIMORDIAL_INFLATION_BG_STEPSIZE,
                                   g_param_spec_double ("primordial-inflation-bg-stepsize",
                                                        NULL,
                                                        "primordial inflation bg stepsize",
                                                        0.0, G_MAXDOUBLE, 0.005,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PRIMORDIAL_INFLATION_TOL_INTEGRATION,
                                   g_param_spec_double ("primordial-inflation-tol-integration",
                                                        NULL,
                                                        "primordial inflation tol integration",
                                                        0.0, G_MAXDOUBLE, 1.0e-3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PRIMORDIAL_INFLATION_ATTRACTOR_PRECISION_PIVOT,
                                   g_param_spec_double ("primordial-inflation-attractor-precision-pivot",
                                                        NULL,
                                                        "primordial inflation attractor_precision_pivot",
                                                        0.0, G_MAXDOUBLE, 0.001,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PRIMORDIAL_INFLATION_ATTRACTOR_PRECISION_INITIAL,
                                   g_param_spec_double ("primordial-inflation-attractor-precision-initial",
                                                        NULL,
                                                        "primordial inflation attractor_precision_initial",
                                                        0.0, G_MAXDOUBLE, 0.1,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PRIMORDIAL_INFLATION_ATTRACTOR_MAXIT,
                                   g_param_spec_int ("primordial-inflation-attractor-maxit",
                                                     NULL,
                                                     "primordial inflation attractor_maxit",
                                                     0, G_MAXINT, 10,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PRIMORDIAL_INFLATION_TOL_CURVATURE,
                                   g_param_spec_double ("primordial-inflation-tol-curvature",
                                                        NULL,
                                                        "primordial inflation tol curvature",
                                                        0.0, G_MAXDOUBLE, 1.0e-3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PRIMORDIAL_INFLATION_AH_INI_TARGET,
                                   g_param_spec_double ("primordial-inflation-aH-ini-target",
                                                        NULL,
                                                        "primordial inflation aH ini target",
                                                        0.0, G_MAXDOUBLE, 0.9,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PRIMORDIAL_INFLATION_END_DPHI,
                                   g_param_spec_double ("primordial-inflation-end-dphi",
                                                        NULL,
                                                        "primordial inflation end dphi",
                                                        0.0, G_MAXDOUBLE, 1.0e-10,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PRIMORDIAL_INFLATION_END_LOGSTEP,
                                   g_param_spec_double ("primordial-inflation-end-logstep",
                                                        NULL,
                                                        "primordial inflation end logstep",
                                                        0.0, G_MAXDOUBLE, 10.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PRIMORDIAL_INFLATION_SMALL_EPSILON,
                                   g_param_spec_double ("primordial-inflation-small-epsilon",
                                                        NULL,
                                                        "primordial inflation small epsilon",
                                                        0.0, G_MAXDOUBLE, 0.1,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PRIMORDIAL_INFLATION_SMALL_EPSILON_TOL,
                                   g_param_spec_double ("primordial-inflation-small-epsilon-tol",
                                                        NULL,
                                                        "primordial inflation small epsilon tol",
                                                        0.0, G_MAXDOUBLE, 0.01,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PRIMORDIAL_INFLATION_EXTRA_EFOLDS,
                                   g_param_spec_double ("primordial-inflation-extra-efolds",
                                                        NULL,
                                                        "primordial inflation extra efolds",
                                                        0.0, G_MAXDOUBLE, 2.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /*
   * Transfer functions parameters
   */
  g_object_class_install_property (object_class,
                                   PROP_L_LOGSTEP,
                                   g_param_spec_double ("l-logstep",
                                                        NULL,
                                                        "maximum spacing of values of l over which Bessel and transfer functions are sampled (so, spacing becomes linear instead of logarithmic at some point)",
                                                        0.0, G_MAXDOUBLE, 1.12,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_L_LINSTEP,
                                   g_param_spec_int ("l-linstep",
                                                     NULL,
                                                     "factor for logarithmic spacing of values of l over which bessel and transfer functions are sampled",
                                                     0, G_MAXINT, 40,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_HYPER_X_MIN,
                                   g_param_spec_double ("hyper-x-min",
                                                        NULL,
                                                        "hyper x min",
                                                        0.0, G_MAXDOUBLE, 1.e-5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_HYPER_SAMPLING_FLAT,
                                   g_param_spec_double ("hyper-sampling-flat",
                                                        NULL,
                                                        "hyper sampling flat",
                                                        0.0, G_MAXDOUBLE, 8.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_HYPER_SAMPLING_CURVED_LOW_NU,
                                   g_param_spec_double ("hyper-sampling-curved-low-nu",
                                                        NULL,
                                                        "hyper sampling_curved_low_nu",
                                                        0.0, G_MAXDOUBLE, 7.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_HYPER_SAMPLING_CURVED_HIGH_NU,
                                   g_param_spec_double ("hyper-sampling-curved-high-nu",
                                                        NULL,
                                                        "hyper sampling_curved_high_nu",
                                                        0.0, G_MAXDOUBLE, 3.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_HYPER_NU_SAMPLING_STEP,
                                   g_param_spec_double ("hyper-nu-sampling-step",
                                                        NULL,
                                                        "hyper nu sampling step",
                                                        0.0, G_MAXDOUBLE, 1000.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_HYPER_PHI_MIN_ABS,
                                   g_param_spec_double ("hyper-phi-min-abs",
                                                        NULL,
                                                        "hyper phi min abs",
                                                        0.0, G_MAXDOUBLE, 1.0e-10,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_HYPER_X_TOL,
                                   g_param_spec_double ("hyper-x-tol",
                                                        NULL,
                                                        "hyper x tol",
                                                        0.0, G_MAXDOUBLE, 1.0e-4,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_HYPER_FLAT_APPROXIMATION_NU,
                                   g_param_spec_double ("hyper-flat-approximation-nu",
                                                        NULL,
                                                        "hyper flat approximation nu",
                                                        0.0, G_MAXDOUBLE, 4000.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_Q_LINSTEP,
                                   g_param_spec_double ("q-linstep",
                                                        NULL,
                                                        "asymptotic linear sampling step in q space, in units of 2pi/r_a(tau_rec) (comoving angular diameter distance to recombination)",
                                                        0.0, G_MAXDOUBLE, 0.45,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_Q_LOGSTEP_SPLINE,
                                   g_param_spec_double ("q-logstep-spline",
                                                        NULL,
                                                        "initial logarithmic sampling step in q space, in units of 2pi/r_a(tau_rec) (comoving angular diameter distance to recombination)",
                                                        0.0, G_MAXDOUBLE, 170.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_Q_LOGSTEP_OPEN,
                                   g_param_spec_double ("q-logstep-open",
                                                        NULL,
                                                        "in open models, the value of q_logstep_spline must be decreased according to curvature. Increasing this number will make the calculation more accurate for large positive Omega_k0",
                                                        0.0, G_MAXDOUBLE, 6.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_Q_LOGSTEP_TRAPZD,
                                   g_param_spec_double ("q-logstep-trapzd",
                                                        NULL,
                                                        "initial logarithmic sampling step in q space, in units of 2pi/r_a(tau_rec) (comoving angular diameter distance to recombination), in the case of small q's in the closed case, for which one must used trapezoidal integration instead of spline (the number of q's for which this is the case decreases with curvature and vanishes in the flat limit)",
                                                        0.0, G_MAXDOUBLE, 20.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_Q_NUMSTEP_TRANSITION,
                                   g_param_spec_double ("q-numstep-transition",
                                                        NULL,
                                                        "number of steps for the transition from q_logstep_trapzd steps to q_logstep_spline steps (transition must be smooth for spline)",
                                                        0.0, G_MAXDOUBLE, 250.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TRANSFER_NEGLECT_DELTA_K_S_T0,
                                   g_param_spec_double ("transfer-neglect-delta-k-S-t0",
                                                        NULL,
                                                        "range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero",
                                                        0.0, G_MAXDOUBLE, 0.15,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TRANSFER_NEGLECT_DELTA_K_S_T1,
                                   g_param_spec_double ("transfer-neglect-delta-k-S-t1",
                                                        NULL,
                                                        "range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero",
                                                        0.0, G_MAXDOUBLE, 0.04,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TRANSFER_NEGLECT_DELTA_K_S_T2,
                                   g_param_spec_double ("transfer-neglect-delta-k-S-t2",
                                                        NULL,
                                                        "range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero",
                                                        0.0, G_MAXDOUBLE, 0.15,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TRANSFER_NEGLECT_DELTA_K_S_E,
                                   g_param_spec_double ("transfer-neglect-delta-k-S-e",
                                                        NULL,
                                                        "range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero",
                                                        0.0, G_MAXDOUBLE, 0.11,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TRANSFER_NEGLECT_DELTA_K_V_T1,
                                   g_param_spec_double ("transfer-neglect-delta-k-V-t1",
                                                        NULL,
                                                        "range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TRANSFER_NEGLECT_DELTA_K_V_T2,
                                   g_param_spec_double ("transfer-neglect-delta-k-V-t2",
                                                        NULL,
                                                        "range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TRANSFER_NEGLECT_DELTA_K_V_E,
                                   g_param_spec_double ("transfer-neglect-delta-k-V-e",
                                                        NULL,
                                                        "range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TRANSFER_NEGLECT_DELTA_K_V_B,
                                   g_param_spec_double ("transfer-neglect-delta-k-V-b",
                                                        NULL,
                                                        "range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TRANSFER_NEGLECT_DELTA_K_T_T2,
                                   g_param_spec_double ("transfer-neglect-delta-k-T-t2",
                                                        NULL,
                                                        "range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero",
                                                        0.0, G_MAXDOUBLE, 0.2,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TRANSFER_NEGLECT_DELTA_K_T_E,
                                   g_param_spec_double ("transfer-neglect-delta-k-T-e",
                                                        NULL,
                                                        "range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero",
                                                        0.0, G_MAXDOUBLE, 0.25,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TRANSFER_NEGLECT_DELTA_K_T_B,
                                   g_param_spec_double ("transfer-neglect-delta-k-T-b",
                                                        NULL,
                                                        "range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero",
                                                        0.0, G_MAXDOUBLE, 0.1,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TRANSFER_NEGLECT_LATE_SOURCE,
                                   g_param_spec_double ("transfer-neglect-late-source",
                                                        NULL,
                                                        "range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero",
                                                        0.0, G_MAXDOUBLE, 400.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_L_SWITCH_LIMBER,
                                   g_param_spec_double ("l-switch-limber",
                                                        NULL,
                                                        "when to use the Limber approximation for project gravitational potential cl's",
                                                        0.0, G_MAXDOUBLE, 10.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_L_SWITCH_LIMBER_FOR_NC_LOCAL_OVER_Z,
                                   g_param_spec_double ("l-switch-limber-for-nc-local-over-z",
                                                        NULL,
                                                        "when to use the Limber approximation for local number count contributions to cl's (relative to central redshift of each bin)",
                                                        0.0, G_MAXDOUBLE, 100.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_L_SWITCH_LIMBER_FOR_NC_LOS_OVER_Z,
                                   g_param_spec_double ("l-switch-limber-for-nc-los-over-z",
                                                        NULL,
                                                        "when to use the Limber approximation for number count contributions to cl's integrated along the line-of-sight (relative to central redshift of each bin)",
                                                        0.0, G_MAXDOUBLE, 30.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SELECTION_CUT_AT_SIGMA,
                                   g_param_spec_double ("selection-cut-at-sigma",
                                                        NULL,
                                                        "in sigma units, where to cut gaussian selection functions",
                                                        0.0, G_MAXDOUBLE, 5.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SELECTION_SAMPLING,
                                   g_param_spec_double ("selection-sampling",
                                                        NULL,
                                                        "controls sampling of integral over time when selection functions vary quicker than Bessel functions. Increase for better sampling.",
                                                        0.0, G_MAXDOUBLE, 50.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SELECTION_SAMPLING_BESSEL,
                                   g_param_spec_double ("selection-sampling-bessel",
                                                        NULL,
                                                        "controls sampling of integral over time when selection functions vary slower than Bessel functions. Increase for better sampling",
                                                        0.0, G_MAXDOUBLE, 20.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SELECTION_SAMPLING_BESSEL_LOS,
                                   g_param_spec_double ("selection-sampling-bessel-los",
                                                        NULL,
                                                        "controls sampling of integral over time when selection functions vary slower than Bessel functions. This parameter is specific to number counts contributions to Cl integrated along the line of sight. Increase for better sampling",
                                                        0.0, G_MAXDOUBLE, 20.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SELECTION_TOPHAT_EDGE,
                                   g_param_spec_double ("selection-tophat-edge",
                                                        NULL,
                                                        "controls how smooth are the edge of top-hat window function (<<1 for very sharp, 0.1 for sharp)",
                                                        0.0, G_MAXDOUBLE, 0.1,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /*
   * Nonlinear module parameters
   */
  g_object_class_install_property (object_class,
                                   PROP_HALOFIT_MIN_K_NONLINEAR,
                                   g_param_spec_double ("halofit-min-k-nonlinear",
                                                        NULL,
                                                        "value of k in 1/Mpc above which non-linear corrections will be computed",
                                                        0.0, G_MAXDOUBLE, 1.0e-4,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_HALOFIT_MIN_K_MAX,
                                   g_param_spec_double ("halofit-min-k-max",
                                                        NULL,
                                                        "when halofit is used, k_max must be at least equal to this value (otherwise halofit could not find the scale of non-linearity)",
                                                        0.0, G_MAXDOUBLE, 5.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_HALOFIT_K_PER_DECADE,
                                   g_param_spec_double ("halofit-k-per-decade",
                                                        NULL,
                                                        "halofit needs to evalute integrals (linear power spectrum times some kernels). They are sampled using this logarithmic step size.",
                                                        0.0, G_MAXDOUBLE, 80.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_HALOFIT_SIGMA_PRECISION,
                                   g_param_spec_double ("halofit-sigma-precision",
                                                        NULL,
                                                        "a smaller value will lead to a more precise halofit result at the highest requested redshift, at the expense of requiring a larger k_max",
                                                        0.0, G_MAXDOUBLE, 0.05,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_HALOFIT_SIGMA_TOL,
                                   g_param_spec_double ("halofit-tol-sigma",
                                                        NULL,
                                                        "tolerance required on sigma(R) when matching the condition sigma(R_nl)=1, whcih defines the wavenumber of non-linearity, k_nl=1./R_nl",
                                                        0.0, G_MAXDOUBLE, 1.0e-6,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PK_EQ_Z_MAX,
                                   g_param_spec_double ("pk-eq-z-max",
                                                        NULL,
                                                        "Maximum z until which the Pk_equal method of 0810.0190 and 1601.07230 is used",
                                                        0.0, G_MAXDOUBLE, 5.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PK_EQ_TOL,
                                   g_param_spec_double ("pk-eq-tol",
                                                        NULL,
                                                        "tolerance for finding the equivalent models of the pk_equal method",
                                                        0.0, G_MAXDOUBLE, 1.0e-7,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /*
   * parameter related to lensing
   */
  g_object_class_install_property (object_class,
                                   PROP_ACCURATE_LENSING,
                                   g_param_spec_int ("accurate-lensing",
                                                     NULL,
                                                     "switch between Gauss-Legendre quadrature integration and simple quadrature on a subdomain of angles",
                                                     0.0, G_MAXINT, _FALSE_,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NUM_MU_MINUS_LMAX,
                                   g_param_spec_int ("num-mu-minus-lmax",
                                                     NULL,
                                                     "difference between num_mu and l_max, increase for more precision",
                                                     0.0, G_MAXINT, 70,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_DELTA_L_MAX,
                                   g_param_spec_int ("delta-l-max",
                                                     NULL,
                                                     "difference between l_max in unlensed and lensed spectra",
                                                     0.0, G_MAXINT, 500,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /*
   * general precision parameters
   */

  g_object_class_install_property (object_class,
                                   PROP_SMALLEST_ALLOWED_VARIATION,
                                   g_param_spec_double ("smallest-allowed-variation",
                                                        NULL,
                                                        "machine-dependent, defined by the implementation",
                                                        0.0, G_MAXDOUBLE, DBL_EPSILON,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TOL_GAUSS_LEGENDRE,
                                   g_param_spec_double ("tol-gauss-legendre",
                                                        NULL,
                                                        "tolerance with which quadrature points are found: must be very small for an accurate integration (if not entered manually, set automatically to match implementation precision)",
                                                        0.0, G_MAXDOUBLE, DBL_EPSILON,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}


/**
 * nc_cbe_precision_new: (constructor)
 *
 * FIXME
 *
 * Returns: (transfer full): a new #NcCBEPrecision
 */
NcCBEPrecision*
nc_cbe_precision_new (void)
{
  NcCBEPrecision* cbe_prec = g_object_new (NC_TYPE_CBE_PRECISION,
                                           NULL);
  return cbe_prec;
}

/**
 * nc_cbe_precision_ref:
 * @cbe_prec: a #NcCBEPrecision.
 *
 * Increases the reference count of @cbe_prec.
 *
 * Returns: (transfer full): @cbe_prec.
 */
NcCBEPrecision*
nc_cbe_precision_ref (NcCBEPrecision* cbe_prec)
{
  return g_object_ref (cbe_prec);
}

/**
 * nc_cbe_precision_free:
 * @cbe_prec: a #NcCBEPrecision.
 *
 * Decreases the reference count of @cbe_prec.
 *
 */
void nc_cbe_precision_free (NcCBEPrecision* cbe_prec)
{
  g_object_unref (cbe_prec);
}

/**
 * nc_cbe_precision_clear:
 * @cbe_prec: a #NcCBEPrecision.
 *
 * Decreases the reference count of *@cbe_prec and sets *@cbe_prec to NULL.
 *
 */
void nc_cbe_precision_clear (NcCBEPrecision** cbe_prec)
{
  g_clear_object (cbe_prec);
}

/**
 * nc_cbe_precision_assert_default:
 * @cbe_prec: a #NcCBEPrecision.
 *
 * Check agaist CLASS default values.
 *
 */
void 
nc_cbe_precision_assert_default (NcCBEPrecision* cbe_prec)
{
  struct precision ppr;
  input_default_precision (&ppr);

#define _CMP_DBL(name) g_assert_cmpfloat (ppr.name, ==, cbe_prec->priv->ppr.name)
#define _CMP_STR(name)                                        \
G_STMT_START                                                  \
{                                                             \
  gchar *s1 = g_path_get_basename (ppr.name);                 \
  gchar *s2 = g_path_get_basename (cbe_prec->priv->ppr.name); \
  g_assert_cmpstr (s1, ==, s2);                               \
  /*printf ("`%s'\n", ppr.name);*/                            \
  g_free (s1);                                                \
  g_free (s2);                                                \
}                                                             \
G_STMT_END

  _CMP_STR (sBBN_file);
  _CMP_STR (hyrec_Alpha_inf_file);
  _CMP_STR (hyrec_R_inf_file);
  _CMP_STR (hyrec_two_photon_tables_file);

  _CMP_DBL (a_ini_over_a_today_default);
  _CMP_DBL (back_integration_stepsize);
  _CMP_DBL (tol_background_integration);
  _CMP_DBL (tol_initial_Omega_r);
  _CMP_DBL (tol_M_ncdm);
  _CMP_DBL (tol_ncdm);
  _CMP_DBL (tol_ncdm_synchronous);
  _CMP_DBL (tol_ncdm_newtonian);
  _CMP_DBL (tol_ncdm_bg);
  _CMP_DBL (tol_ncdm_initial_w);
  _CMP_DBL (safe_phi_scf); /* CHECK: FIXME */
  _CMP_DBL (tol_tau_eq);
  _CMP_DBL (recfast_z_initial);
  _CMP_DBL (recfast_Nz0);
  _CMP_DBL (tol_thermo_integration);
  _CMP_DBL (recfast_Heswitch);
  _CMP_DBL (recfast_fudge_He);
  _CMP_DBL (recfast_Hswitch);
  _CMP_DBL (recfast_fudge_H);
  _CMP_DBL (recfast_delta_fudge_H);
  _CMP_DBL (recfast_AGauss1);
  _CMP_DBL (recfast_AGauss2);
  _CMP_DBL (recfast_zGauss1);
  _CMP_DBL (recfast_zGauss2);
  _CMP_DBL (recfast_wGauss1);
  _CMP_DBL (recfast_wGauss2);
  _CMP_DBL (recfast_z_He_1);
  _CMP_DBL (recfast_delta_z_He_1);
  _CMP_DBL (recfast_z_He_2);
  _CMP_DBL (recfast_delta_z_He_2);
  _CMP_DBL (recfast_z_He_3);
  _CMP_DBL (recfast_delta_z_He_3);
  _CMP_DBL (recfast_x_He0_trigger);
  _CMP_DBL (recfast_x_He0_trigger2);
  _CMP_DBL (recfast_x_He0_trigger_delta);
  _CMP_DBL (recfast_x_H0_trigger);
  _CMP_DBL (recfast_x_H0_trigger2);
  _CMP_DBL (recfast_x_H0_trigger_delta);
  _CMP_DBL (recfast_H_frac);
  _CMP_DBL (reionization_z_start_max);
  _CMP_DBL (reionization_sampling);
  _CMP_DBL (reionization_optical_depth_tol);
  _CMP_DBL (reionization_start_factor);
  _CMP_DBL (thermo_rate_smoothing_radius);
  _CMP_DBL (evolver);
  _CMP_DBL (k_min_tau0);
  _CMP_DBL (k_max_tau0_over_l_max);
  _CMP_DBL (k_step_sub);
  _CMP_DBL (k_step_super);
  _CMP_DBL (k_step_transition);
  _CMP_DBL (k_step_super_reduction);
  _CMP_DBL (k_per_decade_for_pk);
  _CMP_DBL (k_per_decade_for_bao);
  _CMP_DBL (k_bao_center);
  _CMP_DBL (k_bao_width);
  _CMP_DBL (start_small_k_at_tau_c_over_tau_h);
  _CMP_DBL (start_large_k_at_tau_h_over_tau_k);
  _CMP_DBL (tight_coupling_trigger_tau_c_over_tau_h);
  _CMP_DBL (tight_coupling_trigger_tau_c_over_tau_k);
  _CMP_DBL (start_sources_at_tau_c_over_tau_h);
  _CMP_DBL (tight_coupling_approximation);
  _CMP_DBL (l_max_g);
  _CMP_DBL (l_max_pol_g);
  _CMP_DBL (l_max_dr);
  _CMP_DBL (l_max_ur);
  _CMP_DBL (l_max_ncdm);
  _CMP_DBL (l_max_g_ten);
  _CMP_DBL (l_max_pol_g_ten);
  _CMP_DBL (curvature_ini);
  _CMP_DBL (entropy_ini);
  _CMP_DBL (gw_ini);
  _CMP_DBL (perturb_integration_stepsize);
  _CMP_DBL (tol_tau_approx);
  _CMP_DBL (tol_perturb_integration);
  _CMP_DBL (perturb_sampling_stepsize);
  _CMP_DBL (radiation_streaming_approximation);
  _CMP_DBL (radiation_streaming_trigger_tau_over_tau_k);
  _CMP_DBL (radiation_streaming_trigger_tau_c_over_tau);
  _CMP_DBL (ur_fluid_approximation);
  _CMP_DBL (ur_fluid_trigger_tau_over_tau_k);
  _CMP_DBL (ncdm_fluid_approximation);
  _CMP_DBL (ncdm_fluid_trigger_tau_over_tau_k);
  _CMP_DBL (neglect_CMB_sources_below_visibility);
  _CMP_DBL (k_per_decade_primordial);
  _CMP_DBL (primordial_inflation_ratio_min);
  _CMP_DBL (primordial_inflation_ratio_max);
  _CMP_DBL (primordial_inflation_phi_ini_maxit);
  _CMP_DBL (primordial_inflation_pt_stepsize);
  _CMP_DBL (primordial_inflation_bg_stepsize);
  _CMP_DBL (primordial_inflation_tol_integration);
  _CMP_DBL (primordial_inflation_attractor_precision_pivot);
  _CMP_DBL (primordial_inflation_attractor_precision_initial);
  _CMP_DBL (primordial_inflation_attractor_maxit);
  _CMP_DBL (primordial_inflation_tol_curvature);
  _CMP_DBL (primordial_inflation_aH_ini_target);
  _CMP_DBL (primordial_inflation_end_dphi);
  _CMP_DBL (primordial_inflation_end_logstep);
  _CMP_DBL (primordial_inflation_small_epsilon);
  _CMP_DBL (primordial_inflation_small_epsilon_tol);
  _CMP_DBL (primordial_inflation_extra_efolds);
  _CMP_DBL (l_logstep);
  _CMP_DBL (l_linstep);
  _CMP_DBL (hyper_x_min);
  _CMP_DBL (hyper_sampling_flat);
  _CMP_DBL (hyper_sampling_curved_low_nu);
  _CMP_DBL (hyper_sampling_curved_high_nu);
  _CMP_DBL (hyper_nu_sampling_step);
  _CMP_DBL (hyper_phi_min_abs);
  _CMP_DBL (hyper_x_tol);
  _CMP_DBL (hyper_flat_approximation_nu);
  _CMP_DBL (q_linstep);
  _CMP_DBL (q_logstep_spline);
  _CMP_DBL (q_logstep_open);
  _CMP_DBL (q_logstep_trapzd);
  _CMP_DBL (q_numstep_transition);
  _CMP_DBL (transfer_neglect_delta_k_S_t0);
  _CMP_DBL (transfer_neglect_delta_k_S_t1);
  _CMP_DBL (transfer_neglect_delta_k_S_t2);
  _CMP_DBL (transfer_neglect_delta_k_S_e);
  _CMP_DBL (transfer_neglect_delta_k_V_t1);
  _CMP_DBL (transfer_neglect_delta_k_V_t2);
  _CMP_DBL (transfer_neglect_delta_k_V_e);
  _CMP_DBL (transfer_neglect_delta_k_V_b);
  _CMP_DBL (transfer_neglect_delta_k_T_t2);
  _CMP_DBL (transfer_neglect_delta_k_T_e);
  _CMP_DBL (transfer_neglect_delta_k_T_b);
  _CMP_DBL (transfer_neglect_late_source);
  _CMP_DBL (l_switch_limber);
  _CMP_DBL (l_switch_limber_for_nc_local_over_z);
  _CMP_DBL (l_switch_limber_for_nc_los_over_z);
  _CMP_DBL (selection_cut_at_sigma);
  _CMP_DBL (selection_sampling);
  _CMP_DBL (selection_sampling_bessel);
  _CMP_DBL (selection_sampling_bessel_los);
  _CMP_DBL (selection_tophat_edge);
  _CMP_DBL (halofit_min_k_nonlinear);
  _CMP_DBL (halofit_min_k_max);
  _CMP_DBL (halofit_k_per_decade);
  _CMP_DBL (halofit_sigma_precision);
  _CMP_DBL (halofit_tol_sigma);
  _CMP_DBL (pk_eq_z_max);
  _CMP_DBL (pk_eq_tol);
  _CMP_DBL (accurate_lensing);
  _CMP_DBL (num_mu_minus_lmax);
  _CMP_DBL (delta_l_max);
  _CMP_DBL (smallest_allowed_variation);
  _CMP_DBL (tol_gauss_legendre);
}
