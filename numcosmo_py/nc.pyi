import typing

from gi.repository import GLib
from gi.repository import GObject
from gi.repository import NumCosmoMath

T = typing.TypeVar("T")

CLUSTER_MASS_ASCASO_DEFAULT_MU_P0: float = 3.19
CLUSTER_MASS_ASCASO_DEFAULT_MU_P1: int = 0
CLUSTER_MASS_ASCASO_DEFAULT_MU_P2: int = 0
CLUSTER_MASS_ASCASO_DEFAULT_PARAMS_ABSTOL: float = 0.0
CLUSTER_MASS_ASCASO_DEFAULT_SIGMA_P0: float = 0.33
CLUSTER_MASS_ASCASO_DEFAULT_SIGMA_P1: int = 0
CLUSTER_MASS_ASCASO_DEFAULT_SIGMA_P2: float = 0.0
CLUSTER_MASS_BENSON_DEFAULT_A_SZ: float = 5.58
CLUSTER_MASS_BENSON_DEFAULT_B_SZ: float = 1.32
CLUSTER_MASS_BENSON_DEFAULT_C_SZ: float = 0.87
CLUSTER_MASS_BENSON_DEFAULT_D_SZ: float = 0.24
CLUSTER_MASS_BENSON_DEFAULT_PARAMS_ABSTOL: float = 0.0
CLUSTER_MASS_BENSON_M_LOWER_BOUND: float = 10000000000000.0
CLUSTER_MASS_BENSON_XI_ZETA_DIST_CUT: float = 2.0
CLUSTER_MASS_BENSON_XRAY_DEFAULT_A_X: float = 5.77
CLUSTER_MASS_BENSON_XRAY_DEFAULT_B_X: float = 0.57
CLUSTER_MASS_BENSON_XRAY_DEFAULT_C_X: float = 0.4
CLUSTER_MASS_BENSON_XRAY_DEFAULT_D_X: float = 0.12
CLUSTER_MASS_BENSON_XRAY_DEFAULT_PARAMS_ABSTOL: float = 0.0
CLUSTER_MASS_LNNORMAL_DEFAULT_BIAS: float = 0.0
CLUSTER_MASS_LNNORMAL_DEFAULT_PARAMS_ABSTOL: float = 0.0
CLUSTER_MASS_LNNORMAL_DEFAULT_SIGMA: float = 0.04
CLUSTER_MASS_PLCL_DEFAULT_A_L: float = 0.9
CLUSTER_MASS_PLCL_DEFAULT_A_SZ: float = 1.0
CLUSTER_MASS_PLCL_DEFAULT_B_L: float = 0.0
CLUSTER_MASS_PLCL_DEFAULT_B_SZ: float = 0.2
CLUSTER_MASS_PLCL_DEFAULT_COR: float = 0.5
CLUSTER_MASS_PLCL_DEFAULT_PARAMS_ABSTOL: float = 0.0
CLUSTER_MASS_PLCL_DEFAULT_SD_L: float = 0.2
CLUSTER_MASS_PLCL_DEFAULT_SD_SZ: float = 0.3
CLUSTER_MASS_PLCL_MCL: int = 1
CLUSTER_MASS_PLCL_MPL: int = 0
CLUSTER_MASS_PLCL_SD_CL: int = 1
CLUSTER_MASS_PLCL_SD_PL: int = 0
CLUSTER_MASS_VANDERLINDE_DEFAULT_A_SZ: float = 6.01
CLUSTER_MASS_VANDERLINDE_DEFAULT_B_SZ: float = 1.31
CLUSTER_MASS_VANDERLINDE_DEFAULT_C_SZ: float = 1.6
CLUSTER_MASS_VANDERLINDE_DEFAULT_D_SZ: float = 0.21
CLUSTER_MASS_VANDERLINDE_DEFAULT_PARAMS_ABSTOL: float = 0.0
CLUSTER_PHOTOZ_GAUSS_BIAS: int = 0
CLUSTER_PHOTOZ_GAUSS_SIGMA: int = 1
CLUSTER_PSEUDO_COUNTS_DEFAULT_DELTAZ: float = 0.99
CLUSTER_PSEUDO_COUNTS_DEFAULT_LNMCUT: float = 33.0
CLUSTER_PSEUDO_COUNTS_DEFAULT_PARAMS_ABSTOL: float = 0.0
CLUSTER_PSEUDO_COUNTS_DEFAULT_SD_MCUT: float = 0.206
CLUSTER_PSEUDO_COUNTS_DEFAULT_ZMIN: float = 0.188
CLUSTER_REDSHIFT_PHOTOZ_GAUSS_GLOBAL_DEFAULT_BIAS: float = 0.0
CLUSTER_REDSHIFT_PHOTOZ_GAUSS_GLOBAL_DEFAULT_PARAMS_ABSTOL: float = 0.0
CLUSTER_REDSHIFT_PHOTOZ_GAUSS_GLOBAL_DEFAULT_SIGMA0: float = 0.03
DATA_BAO_RDV_LEN: int = 1
DATA_CLUSTER_PSEUDO_COUNTS_RESAMPLE_MAX_TRIES: int = 100000
DATA_SNIA_COV_LEN: int = 1
DATA_SNIA_SIMPLE_LEN: int = 1
DATA_XCOR_DL: int = 10
DATA_XCOR_MAX: int = 5
GALAXY_SD_OBS_REDSHIFT_COL_Z: str = r"z"
GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA: str = r"sigma_z"
GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA0: str = r"sigma_0"
GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_ZP: str = r"zp"
GALAXY_SD_POSITION_COL_DEC: str = r"dec"
GALAXY_SD_POSITION_COL_RA: str = r"ra"
GALAXY_SD_SHAPE_COL_COORD: str = r"coord"
GALAXY_SD_SHAPE_COL_EPSILON_INT_1: str = r"epsilon_int_1"
GALAXY_SD_SHAPE_COL_EPSILON_INT_2: str = r"epsilon_int_2"
GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_1: str = r"epsilon_obs_1"
GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_2: str = r"epsilon_obs_2"
GALAXY_SD_SHAPE_GAUSS_COL_STD_NOISE: str = r"std_noise"
GALAXY_SD_SHAPE_GAUSS_DEFAULT_PARAMS_ABSTOL: float = 0.0
GALAXY_SD_SHAPE_GAUSS_DEFAULT_SIGMA: float = 0.3
GALAXY_SD_SHAPE_GAUSS_HSC_COL_C1: str = r"c1"
GALAXY_SD_SHAPE_GAUSS_HSC_COL_C2: str = r"c2"
GALAXY_SD_SHAPE_GAUSS_HSC_COL_EPSILON_OBS_1: str = r"epsilon_obs_1"
GALAXY_SD_SHAPE_GAUSS_HSC_COL_EPSILON_OBS_2: str = r"epsilon_obs_2"
GALAXY_SD_SHAPE_GAUSS_HSC_COL_M: str = r"m"
GALAXY_SD_SHAPE_GAUSS_HSC_COL_STD_NOISE: str = r"std_noise"
GALAXY_SD_SHAPE_GAUSS_HSC_COL_STD_SHAPE: str = r"std_shape"
GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_ALPHA: float = 0.78
GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_BETA: float = 2.0
GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_PARAMS_ABSTOL: float = 0.0
GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_Z0: float = 0.13
GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_Z_HIGH: float = 20.0
GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_Z_LOW: float = 0.0
GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_ALPHA: float = 0.68
GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_BETA: float = 2.0
GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_Z0: float = 0.11
HALO_CM_DUFFY08_DEFAULT_PARAMS_ABSTOL: float = 0.0
HALO_CM_DUFFY08_LOCAL_SPARAM_LEN: int = 1
HALO_CM_KLYPIN11_DEFAULT_PARAMS_ABSTOL: float = 0.0
HALO_CM_KLYPIN11_LOCAL_SPARAM_LEN: int = 1
HALO_CM_PARAM_DEFAULT_C_DELTA: float = 4.0
HALO_CM_PARAM_DEFAULT_PARAMS_ABSTOL: float = 0.0
HALO_CM_PARAM_LOCAL_SPARAM_LEN: int = 2
HALO_DENSITY_PROFILE_DK14_DEFAULT_BETA: float = 4.0
HALO_DENSITY_PROFILE_DK14_DEFAULT_PARAMS_ABSTOL: float = 0.0
HALO_DENSITY_PROFILE_DK14_DEFAULT_RT: float = 1.0
HALO_DENSITY_PROFILE_EINASTO_DEFAULT_ALPHA: float = 0.25
HALO_DENSITY_PROFILE_EINASTO_DEFAULT_PARAMS_ABSTOL: float = 0.0
HALO_DENSITY_PROFILE_EINASTO_LOCAL_SPARAM_LEN: int = 1
HALO_POSITION_DEFAULT_DEC: float = 0.0
HALO_POSITION_DEFAULT_PARAMS_ABSTOL: float = 0.0
HALO_POSITION_DEFAULT_RA: float = 0.0
HALO_POSITION_DEFAULT_Z: float = 0.2
HICOSMO_DEFAULT_PARAMS_ABSTOL: float = 0.0
HICOSMO_DEFAULT_PARAMS_RELTOL: float = 0.0
HICOSMO_DE_CPL_DEFAULT_W0: float = 1.0
HICOSMO_DE_CPL_DEFAULT_W1: float = 0.0
HICOSMO_DE_CPL_N: int = 9
HICOSMO_DE_DEFAULT_ENNU: float = 3.046
HICOSMO_DE_DEFAULT_HE_YP: float = 0.24
HICOSMO_DE_DEFAULT_NU_G: float = 1.0
HICOSMO_DE_DEFAULT_NU_MASS: float = 1e-05
HICOSMO_DE_DEFAULT_NU_MU: float = 0.0
HICOSMO_DE_DEFAULT_NU_T: float = 0.71611
HICOSMO_DE_DEFAULT_OMEGA_B: float = 0.0432
HICOSMO_DE_DEFAULT_OMEGA_C: float = 0.2568
HICOSMO_DE_DEFAULT_OMEGA_X: float = 0.7
HICOSMO_DE_DEFAULT_T_GAMMA0: float = 2.7245
HICOSMO_DE_JBP_DEFAULT_W0: float = 1.0
HICOSMO_DE_JBP_DEFAULT_W1: float = 0.0
HICOSMO_DE_WSPLINE_DEFAULT_W0: float = 1.0
HICOSMO_DE_WSPLINE_N: int = 5
HICOSMO_DE_XCDM_DEFAULT_W0: float = 1.0
HICOSMO_DE_XCDM_N: int = 8
HICOSMO_GCG_DEFAULT_ENNU: float = 3.046
HICOSMO_GCG_DEFAULT_GAMMA: float = 0.0
HICOSMO_GCG_DEFAULT_HE_YP: float = 0.24
HICOSMO_GCG_DEFAULT_NU_G: float = 1.0
HICOSMO_GCG_DEFAULT_NU_MASS: float = 1e-05
HICOSMO_GCG_DEFAULT_NU_MU: float = 0.0
HICOSMO_GCG_DEFAULT_NU_T: float = 0.71611
HICOSMO_GCG_DEFAULT_OMEGA_B: float = 0.0432
HICOSMO_GCG_DEFAULT_OMEGA_C: float = 0.2568
HICOSMO_GCG_DEFAULT_OMEGA_X: float = 0.7
HICOSMO_GCG_DEFAULT_T_GAMMA0: float = 2.7245
HICOSMO_IDEM2_DEFAULT_ENNU: float = 3.046
HICOSMO_IDEM2_DEFAULT_GAMMA: float = 0.0
HICOSMO_IDEM2_DEFAULT_HE_YP: float = 0.24
HICOSMO_IDEM2_DEFAULT_NU_G: float = 1.0
HICOSMO_IDEM2_DEFAULT_NU_MASS: float = 1e-05
HICOSMO_IDEM2_DEFAULT_NU_MU: float = 0.0
HICOSMO_IDEM2_DEFAULT_NU_T: float = 0.71611
HICOSMO_IDEM2_DEFAULT_OMEGA_B: float = 0.0432
HICOSMO_IDEM2_DEFAULT_OMEGA_C: float = 0.2568
HICOSMO_IDEM2_DEFAULT_OMEGA_X: float = 0.7
HICOSMO_IDEM2_DEFAULT_T_GAMMA0: float = 2.7245
HICOSMO_OMEGA_K0_LIMIT: float = 0.0
HICOSMO_QCONST_DEFAULT_CD: float = 0.0
HICOSMO_QCONST_DEFAULT_E: float = 1.0
HICOSMO_QCONST_DEFAULT_OMEGA_T: float = 1.0
HICOSMO_QCONST_DEFAULT_Q: float = 0.5
HICOSMO_QCONST_DEFAULT_Z1: float = 0.0
HICOSMO_QGRW_DEFAULT_OMEGA_R: float = 1e-05
HICOSMO_QGRW_DEFAULT_OMEGA_W: int = 0
HICOSMO_QGRW_DEFAULT_W: float = 1e-12
HICOSMO_QGRW_DEFAULT_X_B: float = 1e30
HICOSMO_QGW_DEFAULT_OMEGA_W: float = 1.0
HICOSMO_QGW_DEFAULT_W: float = 1e-05
HICOSMO_QGW_DEFAULT_X_B: float = 1e30
HICOSMO_QLINEAR_DEFAULT_CD: float = 0.0
HICOSMO_QLINEAR_DEFAULT_E: float = 1.0
HICOSMO_QLINEAR_DEFAULT_OMEGA_T: float = 1.0
HICOSMO_QLINEAR_DEFAULT_Q: float = 0.5
HICOSMO_QLINEAR_DEFAULT_QP: float = 1.0
HICOSMO_QLINEAR_DEFAULT_Z1: float = 0.0
HICOSMO_QRBF_DEFAULT_AS_DRAG: float = 0.035
HICOSMO_QRBF_DEFAULT_H0: float = 73.0
HICOSMO_QRBF_DEFAULT_OMEGA_T: float = 1.0
HICOSMO_QRBF_DEFAULT_RBF_CENTERS: float = 1.0
HICOSMO_QRBF_DEFAULT_RBF_CENTERS_LEN: int = 3
HICOSMO_QRBF_DEFAULT_RBF_COEFFS: float = 0.5
HICOSMO_QRBF_DEFAULT_RBF_COEFFS_LEN: int = 3
HICOSMO_QRBF_DEFAULT_RBF_H: float = 0.5
HICOSMO_QSPLINE_CONT_PRIOR_ABSTOL: int = 0
HICOSMO_QSPLINE_CONT_PRIOR_LNSIGMA: int = 0
HICOSMO_QSPLINE_DEFAULT_AS_DRAG: float = 0.035
HICOSMO_QSPLINE_DEFAULT_OMEGA_T: float = 1.0
HICOSMO_QSPLINE_DEFAULT_Q: float = 0.5
HICOSMO_QSPLINE_DEFAULT_Q_LEN: int = 3
HICOSMO_VEXP_DEBUG_EVOL_CL: bool = False
HICOSMO_VEXP_DEBUG_EVOL_QT: bool = False
HICOSMO_VEXP_DEFAULT_ALPHA_B: float = 0.1
HICOSMO_VEXP_DEFAULT_D_PHI: float = 0.3
HICOSMO_VEXP_DEFAULT_EM_ALPHA: float = 13.6
HICOSMO_VEXP_DEFAULT_EM_BETA: float = 0.1
HICOSMO_VEXP_DEFAULT_H0: float = 70.0
HICOSMO_VEXP_DEFAULT_OMEGA_C: float = 0.25
HICOSMO_VEXP_DEFAULT_OMEGA_L: float = 0.75
HICOSMO_VEXP_DEFAULT_SIGMA_PHI: float = 0.4
HICOSMO_VEXP_DEFAULT_X_B: float = 1e30
HIPERT_BG_VAR_DEFAULT_ZF: float = 1000000000.0
HIPERT_BOLTZMANN_BASE_SIZE: int = 8
HIPRIM_ATAN_DEFAULT_C2: float = 0.5
HIPRIM_ATAN_DEFAULT_C3: float = 1.0
HIPRIM_ATAN_DEFAULT_LAMBDA: float = 1.0
HIPRIM_ATAN_DEFAULT_LN10E10ASA: float = 3.179
HIPRIM_ATAN_DEFAULT_LNKC: float = 5.3
HIPRIM_ATAN_DEFAULT_N_SA: float = 0.9742
HIPRIM_ATAN_DEFAULT_N_T: float = 0.0
HIPRIM_ATAN_DEFAULT_T_SA_RATIO: float = 0.2
HIPRIM_BPL_DEFAULT_DELTA: float = 1.14
HIPRIM_BPL_DEFAULT_LN10E10ASA: float = 3.179
HIPRIM_BPL_DEFAULT_LNKB: float = 7.55
HIPRIM_BPL_DEFAULT_N_SA: float = 0.9742
HIPRIM_BPL_DEFAULT_N_T: float = 0.0
HIPRIM_BPL_DEFAULT_T_SA_RATIO: float = 0.2
HIPRIM_DEFAULT_K_PIVOT: float = 0.05
HIPRIM_DEFAULT_PARAMS_ABSTOL: float = 0.0
HIPRIM_DEFAULT_PARAMS_RELTOL: float = 0.0
HIPRIM_EXPC_DEFAULT_C: float = 0.5
HIPRIM_EXPC_DEFAULT_LAMBDAC: float = 0.5
HIPRIM_EXPC_DEFAULT_LN10E10ASA: float = 3.179
HIPRIM_EXPC_DEFAULT_LNKC: float = 7.98
HIPRIM_EXPC_DEFAULT_N_SA: float = 0.9742
HIPRIM_EXPC_DEFAULT_N_T: float = 0.0
HIPRIM_EXPC_DEFAULT_T_SA_RATIO: float = 0.2
HIPRIM_POWER_LAW_DEFAULT_LN10E10ASA: float = 3.179
HIPRIM_POWER_LAW_DEFAULT_N_SA: float = 0.9742
HIPRIM_POWER_LAW_DEFAULT_N_T: float = 0.0
HIPRIM_POWER_LAW_DEFAULT_T_SA_RATIO: float = 0.2
HIPRIM_SBPL_DEFAULT_DELTA: float = 0.0
HIPRIM_SBPL_DEFAULT_LAMBDA: float = 10.0
HIPRIM_SBPL_DEFAULT_LN10E10ASA: float = 3.179
HIPRIM_SBPL_DEFAULT_LNKB: float = 7.55
HIPRIM_SBPL_DEFAULT_N_SA: float = 0.9742
HIPRIM_SBPL_DEFAULT_N_T: float = 0.0
HIPRIM_SBPL_DEFAULT_RA: float = 0.8
HIPRIM_SBPL_DEFAULT_T_SA_RATIO: float = 0.2
HIPRIM_TWO_FLUIDS_DEFAULT_LN10E10ASA: float = 3.179
HIPRIM_TWO_FLUIDS_DEFAULT_LNK0: int = 0
HIPRIM_TWO_FLUIDS_DEFAULT_LNW: int = 0
HIPRIM_TWO_FLUIDS_DEFAULT_N_T: float = 0.0
HIPRIM_TWO_FLUIDS_DEFAULT_T_SA_RATIO: float = 0.2
HIREION_CAMB_DEFAULT_HEIII_REION_DELTA: float = 0.5
HIREION_CAMB_DEFAULT_HEIII_Z: float = 3.5
HIREION_CAMB_DEFAULT_HII_HEII_REION_DELTA: float = 0.5
HIREION_CAMB_DEFAULT_HII_HEII_REION_EXPO: float = 1.5
HIREION_CAMB_DEFAULT_HII_HEII_Z: float = 13.0
HIREION_DEFAULT_PARAMS_ABSTOL: float = 0.0
MULTIPLICITY_FUNC_DELTA_C0: float = 1.68647
PLANCK_FI_COR_TTTEEE_DEFAULT_A_cnoise_e2e_100_100_EE: float = 1.0
PLANCK_FI_COR_TTTEEE_DEFAULT_A_cnoise_e2e_143_143_EE: float = 1.0
PLANCK_FI_COR_TTTEEE_DEFAULT_A_cnoise_e2e_217_217_EE: float = 1.0
PLANCK_FI_COR_TTTEEE_DEFAULT_A_pol: float = 1.0
PLANCK_FI_COR_TTTEEE_DEFAULT_A_sbpx_100_100_EE: float = 1.0
PLANCK_FI_COR_TTTEEE_DEFAULT_A_sbpx_100_143_EE: float = 1.0
PLANCK_FI_COR_TTTEEE_DEFAULT_A_sbpx_100_217_EE: float = 1.0
PLANCK_FI_COR_TTTEEE_DEFAULT_A_sbpx_143_143_EE: float = 1.0
PLANCK_FI_COR_TTTEEE_DEFAULT_A_sbpx_143_217_EE: float = 1.0
PLANCK_FI_COR_TTTEEE_DEFAULT_A_sbpx_217_217_EE: float = 1.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_0E_0E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_0E_1E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_0E_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_0T_0E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_0T_1E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_0T_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_1E_1E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_1E_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_1T_1E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_1T_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_2E_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_0_2T_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_0E_0E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_0E_1E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_0E_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_0T_0E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_0T_1E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_0T_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_1E_1E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_1E_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_1T_1E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_1T_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_2E_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_1_2T_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_0E_0E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_0E_1E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_0E_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_0T_0E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_0T_1E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_0T_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_1E_1E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_1E_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_1T_1E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_1T_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_2E_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_2_2T_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_0E_0E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_0E_1E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_0E_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_0T_0E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_0T_1E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_0T_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_1E_1E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_1E_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_1T_1E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_1T_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_2E_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_3_2T_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_0E_0E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_0E_1E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_0E_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_0T_0E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_0T_1E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_0T_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_1E_1E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_1E_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_1T_1E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_1T_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_2E_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_bleak_epsilon_4_2T_2E: float = 0.0
PLANCK_FI_COR_TTTEEE_DEFAULT_calib_100P: float = 1.0
PLANCK_FI_COR_TTTEEE_DEFAULT_calib_143P: float = 1.0
PLANCK_FI_COR_TTTEEE_DEFAULT_calib_217P: float = 1.0
PLANCK_FI_COR_TTTEEE_DEFAULT_galf_EE_A_100: float = 0.06
PLANCK_FI_COR_TTTEEE_DEFAULT_galf_EE_A_100_143: float = 0.05
PLANCK_FI_COR_TTTEEE_DEFAULT_galf_EE_A_100_217: float = 0.11
PLANCK_FI_COR_TTTEEE_DEFAULT_galf_EE_A_143: float = 0.1
PLANCK_FI_COR_TTTEEE_DEFAULT_galf_EE_A_143_217: float = 0.24
PLANCK_FI_COR_TTTEEE_DEFAULT_galf_EE_A_217: float = 0.72
PLANCK_FI_COR_TTTEEE_DEFAULT_galf_EE_index: float = 2.4
PLANCK_FI_COR_TTTEEE_DEFAULT_galf_TE_A_100: float = 0.14
PLANCK_FI_COR_TTTEEE_DEFAULT_galf_TE_A_100_143: float = 0.12
PLANCK_FI_COR_TTTEEE_DEFAULT_galf_TE_A_100_217: float = 0.3
PLANCK_FI_COR_TTTEEE_DEFAULT_galf_TE_A_143: float = 0.24
PLANCK_FI_COR_TTTEEE_DEFAULT_galf_TE_A_143_217: float = 0.6
PLANCK_FI_COR_TTTEEE_DEFAULT_galf_TE_A_217: float = 1.8
PLANCK_FI_COR_TTTEEE_DEFAULT_galf_TE_index: float = 2.4
PLANCK_FI_COR_TT_DEFAULT_A_cib_217: float = 100.0
PLANCK_FI_COR_TT_DEFAULT_A_planck: float = 1.0
PLANCK_FI_COR_TT_DEFAULT_A_sbpx_100_100_TT: float = 1.0
PLANCK_FI_COR_TT_DEFAULT_A_sbpx_143_143_TT: float = 1.0
PLANCK_FI_COR_TT_DEFAULT_A_sbpx_143_217_TT: float = 1.0
PLANCK_FI_COR_TT_DEFAULT_A_sbpx_217_217_TT: float = 1.0
PLANCK_FI_COR_TT_DEFAULT_A_sz: float = 5.0
PLANCK_FI_COR_TT_DEFAULT_calib_100T: float = 0.999
PLANCK_FI_COR_TT_DEFAULT_calib_217T: float = 0.99501
PLANCK_FI_COR_TT_DEFAULT_cib_index: float = 1.3
PLANCK_FI_COR_TT_DEFAULT_gal545_A_100: float = 7.0
PLANCK_FI_COR_TT_DEFAULT_gal545_A_143: float = 9.0
PLANCK_FI_COR_TT_DEFAULT_gal545_A_143_217: float = 21.0
PLANCK_FI_COR_TT_DEFAULT_gal545_A_217: float = 80.0
PLANCK_FI_COR_TT_DEFAULT_ksz_norm: float = 5.0
PLANCK_FI_COR_TT_DEFAULT_ps_A_100_100: float = 200.0
PLANCK_FI_COR_TT_DEFAULT_ps_A_143_143: float = 200.0
PLANCK_FI_COR_TT_DEFAULT_ps_A_143_217: float = 200.0
PLANCK_FI_COR_TT_DEFAULT_ps_A_217_217: float = 200.0
PLANCK_FI_COR_TT_DEFAULT_xi_sz_cib: float = 0.5
PLANCK_FI_DEFAULT_PARAMS_ABSTOL: float = 0.0
POWSPEC_ML_CBE_INTERN_KMAX: float = 10.0
POWSPEC_ML_CBE_INTERN_KMIN: float = 1e-05
POWSPEC_MNL_HALOFIT_F1aPOW: float = 0.0732
POWSPEC_MNL_HALOFIT_F1bPOW: float = 0.0307
POWSPEC_MNL_HALOFIT_F2aPOW: float = 0.1423
POWSPEC_MNL_HALOFIT_F2bPOW: float = 0.0585
POWSPEC_MNL_HALOFIT_F3aPOW: float = 0.0725
POWSPEC_MNL_HALOFIT_F3bPOW: float = 0.0743
POWSPEC_MNL_HALOFIT_LOGRMIN: float = 35.0
RECOMB_SEAGER_HUMMER_HEI_CASE_B_P: float = 0.711
RECOMB_SEAGER_HUMMER_HEI_CASE_B_P_TRIP: float = 0.761
RECOMB_STARTING_X: float = 1000000000000.0
REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_A: float = 0.0
REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_B: float = 0.0
REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_C: float = 0.0
REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_PARAMS_ABSTOL: float = 0.0
REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_VGAMMA: float = 0.05
REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_VSIGMA: float = 0.3
REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_XP: float = 0.2
SCALEFACTOR_DEFAULT_A0: float = 1.0
SCALEFACTOR_DEFAULT_ABSTOL: float = 0.0
SCALEFACTOR_DEFAULT_RELTOL: float = 0.0
SCALEFACTOR_DEFAULT_ZF: float = 100000000000000.0
SCALEFACTOR_MIN_ETA_STEP: float = 0.0
SCALEFACTOR_OMEGA_K_ZERO: float = 0.0
SNIA_DIST_COV_DEFAULT_ALPHA: float = 0.145
SNIA_DIST_COV_DEFAULT_BETA: float = 3.16
SNIA_DIST_COV_DEFAULT_M1: float = 19.168613
SNIA_DIST_COV_DEFAULT_M2: float = 19.185613
SNIA_DIST_COV_DEFAULT_MU: float = 18.0
SNIA_DIST_COV_DEFAULT_PARAMS_ABSTOL: float = 0.0
SNIA_DIST_COV_LNSIGMA_INT_DEFAULT_LEN: int = 4
SNIA_DIST_COV_MU_DEFAULT_LEN: int = 0
WINDOW_VOLUME_GAUSSIAN: int = 0
WINDOW_VOLUME_TOPHAT: int = 0
WL_SURFACE_MASS_DENSITY_DEFAULT_PARAMS_ABSTOL: float = 0.0
WL_SURFACE_MASS_DENSITY_DEFAULT_PCC: float = 0.8
WL_SURFACE_MASS_DENSITY_DEFAULT_ROFF: float = 1.0
XCOR_LIMBER_KERNEL_CMB_ISW_DEFAULT_PARAMS_ABSTOL: float = 0.0
XCOR_LIMBER_KERNEL_CMB_LENSING_DEFAULT_PARAMS_ABSTOL: float = 0.0
XCOR_LIMBER_KERNEL_GAL_BIAS_DEFAULT_LEN: int = 1
XCOR_LIMBER_KERNEL_GAL_DEFAULT_BIAS: float = 1.0
XCOR_LIMBER_KERNEL_GAL_DEFAULT_MAG_BIAS: float = 0.4
XCOR_LIMBER_KERNEL_GAL_DEFAULT_NOISE_BIAS: float = 0.0
XCOR_LIMBER_KERNEL_GAL_DEFAULT_PARAMS_ABSTOL: float = 0.0
XCOR_LIMBER_KERNEL_GAL_G_FUNC_LEN: int = 200
XCOR_LIMBER_KERNEL_WEAK_LENSING_DEFAULT_PARAMS_ABSTOL: float = 0.0
XCOR_PRECISION: float = 1e-06
_lock = ...  # FIXME Constant
_namespace: str = r"NumCosmo"
_version: str = r"1.0"

def bias_mean_prepare(cad: ClusterAbundance, cosmo: HICosmo) -> None: ...
def bias_mean_val(
    cad: ClusterAbundance, cosmo: HICosmo, lnMl: float, lnMu: float, z: float
) -> float: ...
def ca_mean_bias(
    cad: ClusterAbundance, cosmo: HICosmo, lnM: float, z: float
) -> float: ...
def ca_mean_bias_Mobs_denominator(
    cad: ClusterAbundance, cosmo: HICosmo, lnMobs: float, z: float
) -> float: ...
def ca_mean_bias_Mobs_numerator(
    cad: ClusterAbundance, cosmo: HICosmo, lnMobs: float, z: float
) -> float: ...
def ca_mean_bias_denominator(
    cad: ClusterAbundance, cosmo: HICosmo, lnM: float, z: float
) -> float: ...
def ca_mean_bias_numerator(
    cad: ClusterAbundance, cosmo: HICosmo, lnM: float, z: float
) -> float: ...
def data_bao_create(dist: Distance, id: DataBaoId) -> NumCosmoMath.Data: ...
def data_cmb_create(dist: Distance, id: DataCMBId) -> NumCosmoMath.Data: ...
def data_snia_cov_error_quark() -> int: ...
def halo_density_profile_nfw_class_set_ni(num: bool) -> None: ...

class CBE(GObject.Object):
    r"""
    :Constructors:

    ::

        CBE(**properties)
        new() -> NumCosmo.CBE
        prec_file_new(prec_filename:str) -> NumCosmo.CBE
        prec_new(cbe_prec:NumCosmo.CBEPrecision) -> NumCosmo.CBE

    Object NcCBE

    Properties from NcCBE:
      precision -> NcCBEPrecision: precision
        CLASS precision object
      target-Cls -> NcDataCMBDataType: target-Cls
        Target Cls to calculate
      calc-transfer -> gboolean: calc-transfer
        Whether to calculate the transfer function
      use-lensed-Cls -> gboolean: use-lensed-Cls
        Whether to use lensed Cls
      use-tensor -> gboolean: use-tensor
        Whether to use tensor contributions
      use-thermodyn -> gboolean: use-thermodyn
        Whether to use the thermodynamics module
      scalar-lmax -> guint: scalar-lmax
        Scalar modes l_max
      vector-lmax -> guint: vector-lmax
        Vector modes l_max
      tensor-lmax -> guint: tensor-lmax
        Tensor modes l_max
      matter-pk-maxz -> gdouble: matter-pk-maxz
        Maximum redshift for matter Pk
      matter-pk-maxk -> gdouble: matter-pk-maxk
        Maximum mode k for matter Pk
      use-ppf -> gboolean: use-ppf
        Whether to use PPF
      verbosity -> guint: verbosity
        Verbosity

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        calc_transfer: bool
        matter_pk_maxk: float
        matter_pk_maxz: float
        precision: CBEPrecision
        scalar_lmax: int
        target_Cls: DataCMBDataType
        tensor_lmax: int
        use_lensed_Cls: bool
        use_ppf: bool
        use_tensor: bool
        use_thermodyn: bool
        vector_lmax: int
        verbosity: int

    props: Props = ...
    parent_instance: GObject.Object = ...
    prec: CBEPrecision = ...
    priv: CBEPrivate = ...
    a: Scalefactor = ...
    bg_verbose: int = ...
    thermo_verbose: int = ...
    pert_verbose: int = ...
    transfer_verbose: int = ...
    prim_verbose: int = ...
    spectra_verbose: int = ...
    nonlin_verbose: int = ...
    lensing_verbose: int = ...
    target_Cls: DataCMBDataType = ...
    calc_transfer: bool = ...
    use_lensed_Cls: bool = ...
    use_tensor: bool = ...
    use_thermodyn: bool = ...
    scalar_lmax: int = ...
    vector_lmax: int = ...
    tensor_lmax: int = ...
    ctrl_cosmo: NumCosmoMath.ModelCtrl = ...
    ctrl_prim: NumCosmoMath.ModelCtrl = ...
    call: typing.Callable[[CBE, HICosmo], None] = ...
    free: typing.Callable[[CBE], None] = ...
    allocated: bool = ...
    thermodyn_prepared: bool = ...
    def __init__(
        self,
        calc_transfer: bool = ...,
        matter_pk_maxk: float = ...,
        matter_pk_maxz: float = ...,
        precision: CBEPrecision = ...,
        scalar_lmax: int = ...,
        target_Cls: DataCMBDataType = ...,
        tensor_lmax: int = ...,
        use_lensed_Cls: bool = ...,
        use_ppf: bool = ...,
        use_tensor: bool = ...,
        use_thermodyn: bool = ...,
        vector_lmax: int = ...,
        verbosity: int = ...,
    ) -> None: ...
    @staticmethod
    def clear(cbe: CBE) -> None: ...
    def compare_bg(self, cosmo: HICosmo, log_cmp: bool) -> float: ...
    def get_all_Cls(
        self,
        PHIPHI_Cls: NumCosmoMath.Vector,
        TT_Cls: NumCosmoMath.Vector,
        EE_Cls: NumCosmoMath.Vector,
        BB_Cls: NumCosmoMath.Vector,
        TE_Cls: NumCosmoMath.Vector,
    ) -> None: ...
    def get_matter_ps(self) -> NumCosmoMath.Spline2d: ...
    def get_max_matter_pk_k(self) -> float: ...
    def get_max_matter_pk_z(self) -> float: ...
    def get_scalar_lmax(self) -> int: ...
    def get_sigma8(self) -> float: ...
    def get_target_Cls(self) -> DataCMBDataType: ...
    def get_tensor_lmax(self) -> int: ...
    def get_vector_lmax(self) -> int: ...
    def lensed_Cls(self) -> bool: ...
    @classmethod
    def new(cls) -> CBE: ...
    def peek_precision(self) -> CBEPrecision: ...
    @classmethod
    def prec_file_new(cls, prec_filename: str) -> CBE: ...
    @classmethod
    def prec_new(cls, cbe_prec: CBEPrecision) -> CBE: ...
    def prepare(self, cosmo: HICosmo) -> None: ...
    def prepare_if_needed(self, cosmo: HICosmo) -> None: ...
    def ref(self) -> CBE: ...
    def set_calc_transfer(self, calc_transfer: bool) -> None: ...
    def set_lensed_Cls(self, use_lensed_Cls: bool) -> None: ...
    def set_max_matter_pk_k(self, kmax: float) -> None: ...
    def set_max_matter_pk_z(self, zmax: float) -> None: ...
    def set_precision(self, cbe_prec: CBEPrecision) -> None: ...
    def set_scalar_lmax(self, scalar_lmax: int) -> None: ...
    def set_target_Cls(self, target_Cls: DataCMBDataType) -> None: ...
    def set_tensor(self, use_tensor: bool) -> None: ...
    def set_tensor_lmax(self, tensor_lmax: int) -> None: ...
    def set_thermodyn(self, use_thermodyn: bool) -> None: ...
    def set_vector_lmax(self, vector_lmax: int) -> None: ...
    def tensor(self) -> bool: ...
    def thermodyn(self) -> bool: ...
    def thermodyn_get_Xe(self) -> NumCosmoMath.Spline: ...
    def thermodyn_prepare(self, cosmo: HICosmo) -> None: ...
    def thermodyn_prepare_if_needed(self, cosmo: HICosmo) -> None: ...
    def thermodyn_v_tau_max_z(self) -> float: ...
    def thermodyn_z_d(self) -> float: ...
    def use_ppf(self, use_ppf: bool) -> None: ...

class CBEClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        CBEClass()
    """

    parent_class: GObject.ObjectClass = ...

class CBEPrecision(GObject.Object):
    r"""
    :Constructors:

    ::

        CBEPrecision(**properties)
        new() -> NumCosmo.CBEPrecision

    Object NcCBEPrecision

    Properties from NcCBEPrecision:
      a-ini-over-a-today-default -> gdouble: a-ini-over-a-today-default
        default initial value of scale factor in background integration, in units of scale factor today
      back-integration-stepsize -> gdouble: back-integration-stepsize
        default step d tau in background integration, in units of conformal Hubble time ($d \tau$ = back_integration_stepsize / aH )
      tol-background-integration -> gdouble: tol-background-integration
        parameter controlling precision of background integration
      tol-initial-Omega-r -> gdouble: tol-initial-Omega-r
        parameter controlling how deep inside radiation domination must the initial time be chosen
      tol-M-ncdm -> gdouble: tol-M-ncdm
        parameter controlling relative precision of ncdm mass for given ncdm current density
      tol-ncdm -> gdouble: tol-ncdm
        parameter controlling relative precision of integrals over ncdm phase-space distribution during perturbation calculation
      tol-ncdm-synchronous -> gdouble: tol-ncdm-synchronous
        parameter controlling relative precision of integrals over ncdm phase-space distribution during perturbation calculation - synchronous
      tol-ncdm-newtonian -> gdouble: tol-ncdm-newtonian
        parameter controlling relative precision of integrals over ncdm phase-space distribution during perturbation calculation - newtonian
      tol-ncdm-bg -> gdouble: tol-ncdm-bg
        parameter controlling relative precision of integrals over ncdm phase-space distribution during background evolution
      tol-ncdm-initial-w -> gdouble: tol-ncdm-initial-w
        parameter controlling how relativistic must non-cold relics be at initial time
      safe-phi-scf -> gdouble: safe-phi-scf
        parameter controlling the initial scalar field in background functions
      tol-tau-eq -> gdouble: tol-tau-eq
        parameter controlling precision with which tau_eq (conformal time at radiation/matter equality) is found (units: Mpc)
      sBBN-file -> gchararray: sBBN-file
        SBBN filename
      recfast-z-initial -> gdouble: recfast-z-initial
        initial redshift in recfast
      recfast-Nz0 -> gint: recfast-Nz0
        number of integration steps
      tol-thermo-integration -> gdouble: tol-thermo-integration
        precision of each integration step
      recfast-Heswitch -> gint: recfast-Heswitch
        Recfast He switch
      recfast-fudge-He -> gdouble: recfast-fudge-He
        Recfast fudge He
      recfast-Hswitch -> gint: recfast-Hswitch
        recfast 1.5 switching parameter
      recfast-fudge-H -> gdouble: recfast-fudge-H
        H fudge factor when recfast_Hswitch set to false (v1.4 fudging)
      recfast-delta-fudge-H -> gdouble: recfast-delta-fudge-H
        correction to H fudge factor in v1.5
      recfast-AGauss1 -> gdouble: recfast-AGauss1
        Amplitude of 1st Gaussian
      recfast-AGauss2 -> gdouble: recfast-AGauss2
        Amplitude of 2st Gaussian
      recfast-zGauss1 -> gdouble: recfast-zGauss1
        ln(1+z) of 1st Gaussian
      recfast-zGauss2 -> gdouble: recfast-zGauss2
        ln(1+z) of 2st Gaussian
      recfast-wGauss1 -> gdouble: recfast-wGauss1
        Width of 2st Gaussian
      recfast-wGauss2 -> gdouble: recfast-wGauss2
        Width of 2st Gaussian
      recfast-z-He-1 -> gdouble: recfast-z-He-1
        down to which redshift Helium fully ionized
      recfast-delta-z-He-1 -> gdouble: recfast-delta-z-He-1
        z range over which transition is smoothed
      recfast-z-He-2 -> gdouble: recfast-z-He-2
        down to which redshift first Helium recombination not complete
      recfast-delta-z-He-2 -> gdouble: recfast-delta-z-He-2
        z range over which transition is smoothed
      recfast-z-He-3 -> gdouble: recfast-z-He-3
        down to which redshift Helium singly ionized
      recfast-delta-z-He-3 -> gdouble: recfast-delta-z-He-3
        z range over which transition is smoothed
      recfast-x-He0-trigger -> gdouble: recfast-x-He0-trigger
        value below which recfast uses the full equation for Helium
      recfast-x-He0-trigger2 -> gdouble: recfast-x-He0-trigger2
        a second threshold used in derivative routine
      recfast-x-He0-trigger-delta -> gdouble: recfast-x-He0-trigger-delta
        x_He range over which transition is smoothed
      recfast-x-H0-trigger -> gdouble: recfast-x-H0-trigger
        value below which recfast uses the full equation for Hydrogen
      recfast-x-H0-trigger2 -> gdouble: recfast-x-H0-trigger2
        a second threshold used in derivative routine
      recfast-x-H0-trigger-delta -> gdouble: recfast-x-H0-trigger-delta
        x_H range over which transition is smoothed
      recfast-H-frac -> gdouble: recfast-H-frac
        governs time at which full equation of evolution for Tmat is used
      hyrec-Alpha-inf-file -> gchararray: hyrec-Alpha-inf-file
        Hyrec Alpha inf file
      hyrec-R-inf-file -> gchararray: hyrec-R-inf-file
        Hyrec R inf file
      hyrec-two-photon-tables-file -> gchararray: hyrec-two-photon-tables-file
        Hyrec two photon tables file
      reionization-z-start-max -> gdouble: reionization-z-start-max
        maximum redshift at which reionization should start. If not, return an error
      reionization-sampling -> gdouble: reionization-sampling
        control stepsize in z during reionization
      reionization-optical-depth-tol -> gdouble: reionization-optical-depth-tol
        fractional error on optical_depth
      reionization-start-factor -> gdouble: reionization-start-factor
        parameter for CAMB-like parametrization
      thermo-rate-smoothing-radius -> gint: thermo-rate-smoothing-radius
        plays a minor (almost aesthetic) role in the definition of the variation rate of thermodynamical quantities
      evolver -> gint: evolver
        which type of evolver for integrating perturbations (Runge-Kutta? Stiff?...)
      k-min-tau0 -> gdouble: k-min-tau0
        number defining k_min for the computation of Cl's and P(k)'s (dimensionless): (k_min tau_0), usually chosen much smaller than one
      k-max-tau0-over-l-max -> gdouble: k-max-tau0-over-l-max
        number defining k_max for the computation of Cl's (dimensionless): (k_max tau_0)/l_max, usually chosen around two (very relevant for accuracy of lensed ClTT at highest l's)
      k-step-sub -> gdouble: k-step-sub
        step in k space, in units of one period of acoustic oscillation at decoupling, for scales inside sound horizon at decoupling
      k-step-super -> gdouble: k-step-super
        step in k space, in units of one period of acoustic oscillation at decoupling, for scales above sound horizon at decoupling
      k-step-transition -> gdouble: k-step-transition
        dimensionless number regulating the transition from 'sub' steps to 'super' steps. Decrease for more precision
      k-step-super-reduction -> gdouble: k-step-super-reduction
        the step k_step_super is reduced by this amount in the k-->0 limit (below scale of Hubble and/or curvature radius)
      k-per-decade-for-pk -> gdouble: k-per-decade-for-pk
        if values needed between kmax inferred from k_oscillations and k_kmax_for_pk, this gives the number of k per decade outside the BAO region
      k-per-decade-for-bao -> gdouble: k-per-decade-for-bao
        if values needed between kmax inferred from k_oscillations and k_kmax_for_pk, this gives the number of k per decade inside the BAO region (for finer sampling)
      k-bao-center -> gdouble: k-bao-center
        in ln(k) space, the central value of the BAO region where sampling is finer is defined as k_rec times this number (recommended: 3, i.e. finest sampling near 3rd BAO peak)
      k-bao-width -> gdouble: k-bao-width
        in ln(k) space, width of the BAO region where sampling is finer: this number gives roughly the number of BAO oscillations well resolved on both sides of the central value (recommended: 4, i.e. finest sampling from before first up to 3+4=7th peak)
      start-small-k-at-tau-c-over-tau-h -> gdouble: start-small-k-at-tau-c-over-tau-h
        largest wavelengths start being sampled when universe is sufficiently opaque. This is quantified in terms of the ratio of thermo to hubble time scales, $	au_c/	au_H$. Start when start_largek_at_tau_c_over_tau_h equals this ratio. Decrease this value to start integrating the wavenumbers earlier in time.
      start-large-k-at-tau-h-over-tau-k -> gdouble: start-large-k-at-tau-h-over-tau-k
        largest wavelengths start being sampled when mode is sufficiently outside Hibble scale. This is quantified in terms of the ratio of hubble time scale to wavenumber time scale, $	au_h/	au_k$ wich is roughly equal to (k*tau). Start when this ratio equals start_large_k_at_tau_k_over_tau_h. Decrease this value to start integrating the wavenumbers earlier in time.
      tight-coupling-trigger-tau-c-over-tau-h -> gdouble: tight-coupling-trigger-tau-c-over-tau-h
        when to switch off tight-coupling approximation: first condition: $\tau_c/\tau_H$ > tight_coupling_trigger_tau_c_over_tau_h. Decrease this value to switch off earlier in time.  If this number is larger than start_sources_at_tau_c_over_tau_h, the code returns an error, because the source computation requires tight-coupling to be switched off.
      tight-coupling-trigger-tau-c-over-tau-k -> gdouble: tight-coupling-trigger-tau-c-over-tau-k
        when to switch off tight-coupling approximation: second condition: $\tau_c/\tau_k \equiv k \tau_c$ < tight_coupling_trigger_tau_c_over_tau_k. Decrease this value to switch off earlier in time.
      start-sources-at-tau-c-over-tau-h -> gdouble: start-sources-at-tau-c-over-tau-h
        sources start being sampled when universe is sufficiently opaque. This is quantified in terms of the ratio of thermo to hubble time scales, $\tau_c/\tau_H$. Start when start_sources_at_tau_c_over_tau_h equals this ratio. Decrease this value to start sampling the sources earlier in time.
      tight-coupling-approximation -> gint: tight-coupling-approximation
        Tight coupling approximation scheme
      l-max-g -> gint: l-max-g
        number of momenta in Boltzmann hierarchy for photon temperature (scalar)
      l-max-pol-g -> gint: l-max-pol-g
        number of momenta in Boltzmann hierarchy for photon polarisation (scalar)
      l-max-dr -> gint: l-max-dr
        number of momenta in Boltzmann hierarchy for decay radiation
      l-max-ur -> gint: l-max-ur
        number of momenta in Boltzmann hierarchy for relativistic neutrino/relics (scalar)
      l-max-ncdm -> gint: l-max-ncdm
        number of momenta in Boltzmann hierarchy for relativistic neutrino/relics (scalar)
      l-max-g-ten -> gint: l-max-g-ten
        number of momenta in Boltzmann hierarchy for photon temperature (tensor)
      l-max-pol-g-ten -> gint: l-max-pol-g-ten
        number of momenta in Boltzmann hierarchy for photon polarisation (tensor)
      curvature-ini -> gdouble: curvature-ini
        initial curvature; used to fix adiabatic initial conditions; must remain fixed to one as long as the primordial adiabatic spectrum stands for the curvature power spectrum
      entropy-ini -> gdouble: entropy-ini
        initial entropy; used to fix isocurvature initial conditions; must remain fixed to one as long as the primordial isocurvature spectrum stands for an entropy power spectrum
      gw-ini -> gdouble: gw-ini
        initial condition for tensor metric perturbation h
      perturb-integration-stepsize -> gdouble: perturb-integration-stepsize
        default step $d \tau$ in perturbation integration, in units of the timescale involved in the equations (usally, the min of $1/k$, $1/aH$, $1/\dot{\kappa}$)
      tol-tau-approx -> gdouble: tol-tau-approx
        precision with which the code should determine (by bisection) the times at which sources start being sampled, and at which approximations must be switched on/off (units of Mpc)
      tol-perturb-integration -> gdouble: tol-perturb-integration
        control parameter for the precision of the perturbation integration
      perturb-sampling-stepsize -> gdouble: perturb-sampling-stepsize
        default step $d \tau$ for sampling the source function, in units of the timescale involved in the sources: $(\dot{\kappa}- \ddot{\kappa}/\dot{\kappa})^{-1}$
      radiation-streaming-approximation -> gint: radiation-streaming-approximation
        method for switching off photon perturbations
      radiation-streaming-trigger-tau-over-tau-k -> gdouble: radiation-streaming-trigger-tau-over-tau-k
        when to switch off photon perturbations, ie when to switch on photon free-streaming approximation (keep density and thtau, set shear and higher momenta to zero): first condition: $k 	au$ > radiation_streaming_trigger_tau_h_over_tau_k
      radiation-streaming-trigger-tau-c-over-tau -> gdouble: radiation-streaming-trigger-tau-c-over-tau
        when to switch off photon perturbations, ie when to switch on photon free-streaming approximation (keep density and theta, set shear and higher momenta to zero): second condition:
      ur-fluid-approximation -> gint: ur-fluid-approximation
        UR fluid approximation scheme
      ur-fluid-trigger-tau-over-tau-k -> gdouble: ur-fluid-trigger-tau-over-tau-k
        when to switch off ur (massless neutrinos / ultra-relativistic relics) fluid approximation
      ncdm-fluid-approximation -> gint: ncdm-fluid-approximation
        NCDM fluid approximation scheme
      ncdm-fluid-trigger-tau-over-tau-k -> gdouble: ncdm-fluid-trigger-tau-over-tau-k
        when to switch off ncdm (massive neutrinos / non-cold relics) fluid approximation
      neglect-CMB-sources-below-visibility -> gdouble: neglect-CMB-sources-below-visibility
        neglect CMB sources below visibility
      k-per-decade-primordial -> gdouble: k-per-decade-primordial
        logarithmic sampling for primordial spectra (number of points per decade in k space)
      primordial-inflation-ratio-min -> gdouble: primordial-inflation-ratio-min
        primordial inflation ratio min
      primordial-inflation-ratio-max -> gdouble: primordial-inflation-ratio-max
        primordial inflation ratio max
      primordial-inflation-phi-ini-maxit -> gint: primordial-inflation-phi-ini-maxit
        primordial inflation phi ini maxit
      primordial-inflation-pt-stepsize -> gdouble: primordial-inflation-pt-stepsize
        primordial inflation pt stepsize
      primordial-inflation-bg-stepsize -> gdouble: primordial-inflation-bg-stepsize
        primordial inflation bg stepsize
      primordial-inflation-tol-integration -> gdouble: primordial-inflation-tol-integration
        primordial inflation tol integration
      primordial-inflation-attractor-precision-pivot -> gdouble: primordial-inflation-attractor-precision-pivot
        primordial inflation attractor_precision_pivot
      primordial-inflation-attractor-precision-initial -> gdouble: primordial-inflation-attractor-precision-initial
        primordial inflation attractor_precision_initial
      primordial-inflation-attractor-maxit -> gint: primordial-inflation-attractor-maxit
        primordial inflation attractor_maxit
      primordial-inflation-tol-curvature -> gdouble: primordial-inflation-tol-curvature
        primordial inflation tol curvature
      primordial-inflation-aH-ini-target -> gdouble: primordial-inflation-aH-ini-target
        primordial inflation aH ini target
      primordial-inflation-end-dphi -> gdouble: primordial-inflation-end-dphi
        primordial inflation end dphi
      primordial-inflation-end-logstep -> gdouble: primordial-inflation-end-logstep
        primordial inflation end logstep
      primordial-inflation-small-epsilon -> gdouble: primordial-inflation-small-epsilon
        primordial inflation small epsilon
      primordial-inflation-small-epsilon-tol -> gdouble: primordial-inflation-small-epsilon-tol
        primordial inflation small epsilon tol
      primordial-inflation-extra-efolds -> gdouble: primordial-inflation-extra-efolds
        primordial inflation extra efolds
      l-logstep -> gdouble: l-logstep
        maximum spacing of values of l over which Bessel and transfer functions are sampled (so, spacing becomes linear instead of logarithmic at some point)
      l-linstep -> gint: l-linstep
        factor for logarithmic spacing of values of l over which bessel and transfer functions are sampled
      hyper-x-min -> gdouble: hyper-x-min
        hyper x min
      hyper-sampling-flat -> gdouble: hyper-sampling-flat
        hyper sampling flat
      hyper-sampling-curved-low-nu -> gdouble: hyper-sampling-curved-low-nu
        hyper sampling_curved_low_nu
      hyper-sampling-curved-high-nu -> gdouble: hyper-sampling-curved-high-nu
        hyper sampling_curved_high_nu
      hyper-nu-sampling-step -> gdouble: hyper-nu-sampling-step
        hyper nu sampling step
      hyper-phi-min-abs -> gdouble: hyper-phi-min-abs
        hyper phi min abs
      hyper-x-tol -> gdouble: hyper-x-tol
        hyper x tol
      hyper-flat-approximation-nu -> gdouble: hyper-flat-approximation-nu
        hyper flat approximation nu
      q-linstep -> gdouble: q-linstep
        asymptotic linear sampling step in q space, in units of 2pi/r_a(tau_rec) (comoving angular diameter distance to recombination)
      q-logstep-spline -> gdouble: q-logstep-spline
        initial logarithmic sampling step in q space, in units of 2pi/r_a(tau_rec) (comoving angular diameter distance to recombination)
      q-logstep-open -> gdouble: q-logstep-open
        in open models, the value of q_logstep_spline must be decreased according to curvature. Increasing this number will make the calculation more accurate for large positive Omega_k0
      q-logstep-trapzd -> gdouble: q-logstep-trapzd
        initial logarithmic sampling step in q space, in units of 2pi/r_a(tau_rec) (comoving angular diameter distance to recombination), in the case of small q's in the closed case, for which one must used trapezoidal integration instead of spline (the number of q's for which this is the case decreases with curvature and vanishes in the flat limit)
      q-numstep-transition -> gdouble: q-numstep-transition
        number of steps for the transition from q_logstep_trapzd steps to q_logstep_spline steps (transition must be smooth for spline)
      transfer-neglect-delta-k-S-t0 -> gdouble: transfer-neglect-delta-k-S-t0
        range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero
      transfer-neglect-delta-k-S-t1 -> gdouble: transfer-neglect-delta-k-S-t1
        range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero
      transfer-neglect-delta-k-S-t2 -> gdouble: transfer-neglect-delta-k-S-t2
        range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero
      transfer-neglect-delta-k-S-e -> gdouble: transfer-neglect-delta-k-S-e
        range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero
      transfer-neglect-delta-k-V-t1 -> gdouble: transfer-neglect-delta-k-V-t1
        range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero
      transfer-neglect-delta-k-V-t2 -> gdouble: transfer-neglect-delta-k-V-t2
        range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero
      transfer-neglect-delta-k-V-e -> gdouble: transfer-neglect-delta-k-V-e
        range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero
      transfer-neglect-delta-k-V-b -> gdouble: transfer-neglect-delta-k-V-b
        range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero
      transfer-neglect-delta-k-T-t2 -> gdouble: transfer-neglect-delta-k-T-t2
        range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero
      transfer-neglect-delta-k-T-e -> gdouble: transfer-neglect-delta-k-T-e
        range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero
      transfer-neglect-delta-k-T-b -> gdouble: transfer-neglect-delta-k-T-b
        range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero
      transfer-neglect-late-source -> gdouble: transfer-neglect-late-source
        range of k values (in 1/Mpc) taken into account in transfer function: for l < (k-delta_k)*tau0, ie for k > (l/tau0 + delta_k), the transfer function is set to zero
      l-switch-limber -> gdouble: l-switch-limber
        when to use the Limber approximation for project gravitational potential cl's
      l-switch-limber-for-nc-local-over-z -> gdouble: l-switch-limber-for-nc-local-over-z
        when to use the Limber approximation for local number count contributions to cl's (relative to central redshift of each bin)
      l-switch-limber-for-nc-los-over-z -> gdouble: l-switch-limber-for-nc-los-over-z
        when to use the Limber approximation for number count contributions to cl's integrated along the line-of-sight (relative to central redshift of each bin)
      selection-cut-at-sigma -> gdouble: selection-cut-at-sigma
        in sigma units, where to cut gaussian selection functions
      selection-sampling -> gdouble: selection-sampling
        controls sampling of integral over time when selection functions vary quicker than Bessel functions. Increase for better sampling.
      selection-sampling-bessel -> gdouble: selection-sampling-bessel
        controls sampling of integral over time when selection functions vary slower than Bessel functions. Increase for better sampling
      selection-sampling-bessel-los -> gdouble: selection-sampling-bessel-los
        controls sampling of integral over time when selection functions vary slower than Bessel functions. This parameter is specific to number counts contributions to Cl integrated along the line of sight. Increase for better sampling
      selection-tophat-edge -> gdouble: selection-tophat-edge
        controls how smooth are the edge of top-hat window function (<<1 for very sharp, 0.1 for sharp)
      halofit-min-k-nonlinear -> gdouble: halofit-min-k-nonlinear
        value of k in 1/Mpc above which non-linear corrections will be computed
      halofit-min-k-max -> gdouble: halofit-min-k-max
        when halofit is used, k_max must be at least equal to this value (otherwise halofit could not find the scale of non-linearity)
      halofit-k-per-decade -> gdouble: halofit-k-per-decade
        halofit needs to evalute integrals (linear power spectrum times some kernels). They are sampled using this logarithmic step size.
      halofit-sigma-precision -> gdouble: halofit-sigma-precision
        a smaller value will lead to a more precise halofit result at the highest requested redshift, at the expense of requiring a larger k_max
      halofit-tol-sigma -> gdouble: halofit-tol-sigma
        tolerance required on sigma(R) when matching the condition sigma(R_nl)=1, whcih defines the wavenumber of non-linearity, k_nl=1./R_nl
      pk-eq-z-max -> gdouble: pk-eq-z-max
        Maximum z until which the Pk_equal method of 0810.0190 and 1601.07230 is used
      pk-eq-tol -> gdouble: pk-eq-tol
        tolerance for finding the equivalent models of the pk_equal method
      accurate-lensing -> gint: accurate-lensing
        switch between Gauss-Legendre quadrature integration and simple quadrature on a subdomain of angles
      num-mu-minus-lmax -> gint: num-mu-minus-lmax
        difference between num_mu and l_max, increase for more precision
      delta-l-max -> gint: delta-l-max
        difference between l_max in unlensed and lensed spectra
      tol-gauss-legendre -> gdouble: tol-gauss-legendre
        tolerance with which quadrature points are found: must be very small for an accurate integration (if not entered manually, set automatically to match implementation precision)
      smallest-allowed-variation -> gdouble: smallest-allowed-variation
        machine-dependent, defined by the implementation

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        a_ini_over_a_today_default: float
        accurate_lensing: int
        back_integration_stepsize: float
        curvature_ini: float
        delta_l_max: int
        entropy_ini: float
        evolver: int
        gw_ini: float
        halofit_k_per_decade: float
        halofit_min_k_max: float
        halofit_min_k_nonlinear: float
        halofit_sigma_precision: float
        halofit_tol_sigma: float
        hyper_flat_approximation_nu: float
        hyper_nu_sampling_step: float
        hyper_phi_min_abs: float
        hyper_sampling_curved_high_nu: float
        hyper_sampling_curved_low_nu: float
        hyper_sampling_flat: float
        hyper_x_min: float
        hyper_x_tol: float
        hyrec_Alpha_inf_file: str
        hyrec_R_inf_file: str
        hyrec_two_photon_tables_file: str
        k_bao_center: float
        k_bao_width: float
        k_max_tau0_over_l_max: float
        k_min_tau0: float
        k_per_decade_for_bao: float
        k_per_decade_for_pk: float
        k_per_decade_primordial: float
        k_step_sub: float
        k_step_super: float
        k_step_super_reduction: float
        k_step_transition: float
        l_linstep: int
        l_logstep: float
        l_max_dr: int
        l_max_g: int
        l_max_g_ten: int
        l_max_ncdm: int
        l_max_pol_g: int
        l_max_pol_g_ten: int
        l_max_ur: int
        l_switch_limber: float
        l_switch_limber_for_nc_local_over_z: float
        l_switch_limber_for_nc_los_over_z: float
        ncdm_fluid_approximation: int
        ncdm_fluid_trigger_tau_over_tau_k: float
        neglect_CMB_sources_below_visibility: float
        num_mu_minus_lmax: int
        perturb_integration_stepsize: float
        perturb_sampling_stepsize: float
        pk_eq_tol: float
        pk_eq_z_max: float
        primordial_inflation_aH_ini_target: float
        primordial_inflation_attractor_maxit: int
        primordial_inflation_attractor_precision_initial: float
        primordial_inflation_attractor_precision_pivot: float
        primordial_inflation_bg_stepsize: float
        primordial_inflation_end_dphi: float
        primordial_inflation_end_logstep: float
        primordial_inflation_extra_efolds: float
        primordial_inflation_phi_ini_maxit: int
        primordial_inflation_pt_stepsize: float
        primordial_inflation_ratio_max: float
        primordial_inflation_ratio_min: float
        primordial_inflation_small_epsilon: float
        primordial_inflation_small_epsilon_tol: float
        primordial_inflation_tol_curvature: float
        primordial_inflation_tol_integration: float
        q_linstep: float
        q_logstep_open: float
        q_logstep_spline: float
        q_logstep_trapzd: float
        q_numstep_transition: float
        radiation_streaming_approximation: int
        radiation_streaming_trigger_tau_c_over_tau: float
        radiation_streaming_trigger_tau_over_tau_k: float
        recfast_AGauss1: float
        recfast_AGauss2: float
        recfast_H_frac: float
        recfast_Heswitch: int
        recfast_Hswitch: int
        recfast_Nz0: int
        recfast_delta_fudge_H: float
        recfast_delta_z_He_1: float
        recfast_delta_z_He_2: float
        recfast_delta_z_He_3: float
        recfast_fudge_H: float
        recfast_fudge_He: float
        recfast_wGauss1: float
        recfast_wGauss2: float
        recfast_x_H0_trigger: float
        recfast_x_H0_trigger_delta: float
        recfast_x_H0_trigger2: float
        recfast_x_He0_trigger: float
        recfast_x_He0_trigger_delta: float
        recfast_x_He0_trigger2: float
        recfast_z_He_1: float
        recfast_z_He_2: float
        recfast_z_He_3: float
        recfast_z_initial: float
        recfast_zGauss1: float
        recfast_zGauss2: float
        reionization_optical_depth_tol: float
        reionization_sampling: float
        reionization_start_factor: float
        reionization_z_start_max: float
        sBBN_file: str
        safe_phi_scf: float
        selection_cut_at_sigma: float
        selection_sampling: float
        selection_sampling_bessel: float
        selection_sampling_bessel_los: float
        selection_tophat_edge: float
        smallest_allowed_variation: float
        start_large_k_at_tau_h_over_tau_k: float
        start_small_k_at_tau_c_over_tau_h: float
        start_sources_at_tau_c_over_tau_h: float
        thermo_rate_smoothing_radius: int
        tight_coupling_approximation: int
        tight_coupling_trigger_tau_c_over_tau_h: float
        tight_coupling_trigger_tau_c_over_tau_k: float
        tol_M_ncdm: float
        tol_background_integration: float
        tol_gauss_legendre: float
        tol_initial_Omega_r: float
        tol_ncdm: float
        tol_ncdm_bg: float
        tol_ncdm_initial_w: float
        tol_ncdm_newtonian: float
        tol_ncdm_synchronous: float
        tol_perturb_integration: float
        tol_tau_approx: float
        tol_tau_eq: float
        tol_thermo_integration: float
        transfer_neglect_delta_k_S_e: float
        transfer_neglect_delta_k_S_t0: float
        transfer_neglect_delta_k_S_t1: float
        transfer_neglect_delta_k_S_t2: float
        transfer_neglect_delta_k_T_b: float
        transfer_neglect_delta_k_T_e: float
        transfer_neglect_delta_k_T_t2: float
        transfer_neglect_delta_k_V_b: float
        transfer_neglect_delta_k_V_e: float
        transfer_neglect_delta_k_V_t1: float
        transfer_neglect_delta_k_V_t2: float
        transfer_neglect_late_source: float
        ur_fluid_approximation: int
        ur_fluid_trigger_tau_over_tau_k: float

    props: Props = ...
    parent_instance: GObject.Object = ...
    priv: CBEPrecisionPrivate = ...
    def __init__(
        self,
        a_ini_over_a_today_default: float = ...,
        accurate_lensing: int = ...,
        back_integration_stepsize: float = ...,
        curvature_ini: float = ...,
        delta_l_max: int = ...,
        entropy_ini: float = ...,
        evolver: int = ...,
        gw_ini: float = ...,
        halofit_k_per_decade: float = ...,
        halofit_min_k_max: float = ...,
        halofit_min_k_nonlinear: float = ...,
        halofit_sigma_precision: float = ...,
        halofit_tol_sigma: float = ...,
        hyper_flat_approximation_nu: float = ...,
        hyper_nu_sampling_step: float = ...,
        hyper_phi_min_abs: float = ...,
        hyper_sampling_curved_high_nu: float = ...,
        hyper_sampling_curved_low_nu: float = ...,
        hyper_sampling_flat: float = ...,
        hyper_x_min: float = ...,
        hyper_x_tol: float = ...,
        hyrec_Alpha_inf_file: str = ...,
        hyrec_R_inf_file: str = ...,
        hyrec_two_photon_tables_file: str = ...,
        k_bao_center: float = ...,
        k_bao_width: float = ...,
        k_max_tau0_over_l_max: float = ...,
        k_min_tau0: float = ...,
        k_per_decade_for_bao: float = ...,
        k_per_decade_for_pk: float = ...,
        k_per_decade_primordial: float = ...,
        k_step_sub: float = ...,
        k_step_super: float = ...,
        k_step_super_reduction: float = ...,
        k_step_transition: float = ...,
        l_linstep: int = ...,
        l_logstep: float = ...,
        l_max_dr: int = ...,
        l_max_g: int = ...,
        l_max_g_ten: int = ...,
        l_max_ncdm: int = ...,
        l_max_pol_g: int = ...,
        l_max_pol_g_ten: int = ...,
        l_max_ur: int = ...,
        l_switch_limber: float = ...,
        l_switch_limber_for_nc_local_over_z: float = ...,
        l_switch_limber_for_nc_los_over_z: float = ...,
        ncdm_fluid_approximation: int = ...,
        ncdm_fluid_trigger_tau_over_tau_k: float = ...,
        neglect_CMB_sources_below_visibility: float = ...,
        num_mu_minus_lmax: int = ...,
        perturb_integration_stepsize: float = ...,
        perturb_sampling_stepsize: float = ...,
        pk_eq_tol: float = ...,
        pk_eq_z_max: float = ...,
        primordial_inflation_aH_ini_target: float = ...,
        primordial_inflation_attractor_maxit: int = ...,
        primordial_inflation_attractor_precision_initial: float = ...,
        primordial_inflation_attractor_precision_pivot: float = ...,
        primordial_inflation_bg_stepsize: float = ...,
        primordial_inflation_end_dphi: float = ...,
        primordial_inflation_end_logstep: float = ...,
        primordial_inflation_extra_efolds: float = ...,
        primordial_inflation_phi_ini_maxit: int = ...,
        primordial_inflation_pt_stepsize: float = ...,
        primordial_inflation_ratio_max: float = ...,
        primordial_inflation_ratio_min: float = ...,
        primordial_inflation_small_epsilon: float = ...,
        primordial_inflation_small_epsilon_tol: float = ...,
        primordial_inflation_tol_curvature: float = ...,
        primordial_inflation_tol_integration: float = ...,
        q_linstep: float = ...,
        q_logstep_open: float = ...,
        q_logstep_spline: float = ...,
        q_logstep_trapzd: float = ...,
        q_numstep_transition: float = ...,
        radiation_streaming_approximation: int = ...,
        radiation_streaming_trigger_tau_c_over_tau: float = ...,
        radiation_streaming_trigger_tau_over_tau_k: float = ...,
        recfast_AGauss1: float = ...,
        recfast_AGauss2: float = ...,
        recfast_H_frac: float = ...,
        recfast_Heswitch: int = ...,
        recfast_Hswitch: int = ...,
        recfast_Nz0: int = ...,
        recfast_delta_fudge_H: float = ...,
        recfast_delta_z_He_1: float = ...,
        recfast_delta_z_He_2: float = ...,
        recfast_delta_z_He_3: float = ...,
        recfast_fudge_H: float = ...,
        recfast_fudge_He: float = ...,
        recfast_wGauss1: float = ...,
        recfast_wGauss2: float = ...,
        recfast_x_H0_trigger: float = ...,
        recfast_x_H0_trigger_delta: float = ...,
        recfast_x_H0_trigger2: float = ...,
        recfast_x_He0_trigger: float = ...,
        recfast_x_He0_trigger_delta: float = ...,
        recfast_x_He0_trigger2: float = ...,
        recfast_z_He_1: float = ...,
        recfast_z_He_2: float = ...,
        recfast_z_He_3: float = ...,
        recfast_z_initial: float = ...,
        recfast_zGauss1: float = ...,
        recfast_zGauss2: float = ...,
        reionization_optical_depth_tol: float = ...,
        reionization_sampling: float = ...,
        reionization_start_factor: float = ...,
        reionization_z_start_max: float = ...,
        sBBN_file: str = ...,
        safe_phi_scf: float = ...,
        selection_cut_at_sigma: float = ...,
        selection_sampling: float = ...,
        selection_sampling_bessel: float = ...,
        selection_sampling_bessel_los: float = ...,
        selection_tophat_edge: float = ...,
        smallest_allowed_variation: float = ...,
        start_large_k_at_tau_h_over_tau_k: float = ...,
        start_small_k_at_tau_c_over_tau_h: float = ...,
        start_sources_at_tau_c_over_tau_h: float = ...,
        thermo_rate_smoothing_radius: int = ...,
        tight_coupling_approximation: int = ...,
        tight_coupling_trigger_tau_c_over_tau_h: float = ...,
        tight_coupling_trigger_tau_c_over_tau_k: float = ...,
        tol_M_ncdm: float = ...,
        tol_background_integration: float = ...,
        tol_gauss_legendre: float = ...,
        tol_initial_Omega_r: float = ...,
        tol_ncdm: float = ...,
        tol_ncdm_bg: float = ...,
        tol_ncdm_initial_w: float = ...,
        tol_ncdm_newtonian: float = ...,
        tol_ncdm_synchronous: float = ...,
        tol_perturb_integration: float = ...,
        tol_tau_approx: float = ...,
        tol_tau_eq: float = ...,
        tol_thermo_integration: float = ...,
        transfer_neglect_delta_k_S_e: float = ...,
        transfer_neglect_delta_k_S_t0: float = ...,
        transfer_neglect_delta_k_S_t1: float = ...,
        transfer_neglect_delta_k_S_t2: float = ...,
        transfer_neglect_delta_k_T_b: float = ...,
        transfer_neglect_delta_k_T_e: float = ...,
        transfer_neglect_delta_k_T_t2: float = ...,
        transfer_neglect_delta_k_V_b: float = ...,
        transfer_neglect_delta_k_V_e: float = ...,
        transfer_neglect_delta_k_V_t1: float = ...,
        transfer_neglect_delta_k_V_t2: float = ...,
        transfer_neglect_late_source: float = ...,
        ur_fluid_approximation: int = ...,
        ur_fluid_trigger_tau_over_tau_k: float = ...,
    ) -> None: ...
    def assert_default(self) -> None: ...
    @staticmethod
    def clear(cbe_prec: CBEPrecision) -> None: ...
    def free(self) -> None: ...
    @classmethod
    def new(cls) -> CBEPrecision: ...
    def ref(self) -> CBEPrecision: ...

class CBEPrecisionClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        CBEPrecisionClass()
    """

    parent_class: GObject.ObjectClass = ...

class CBEPrecisionPrivate(GObject.GPointer): ...
class CBEPrivate(GObject.GPointer): ...

class ClusterAbundance(GObject.Object):
    r"""
    :Constructors:

    ::

        ClusterAbundance(**properties)
        new(mfp:NumCosmo.HaloMassFunction, mbiasf:NumCosmo.HaloBias=None) -> NumCosmo.ClusterAbundance
        nodist_new(mfp:NumCosmo.HaloMassFunction, mbiasf:NumCosmo.HaloBias=None) -> NumCosmo.ClusterAbundance

    Object NcClusterAbundance

    Properties from NcClusterAbundance:
      halo-mass-function -> NcHaloMassFunction: halo-mass-function
        Mass Function
      mean-bias -> NcHaloBias: mean-bias
        Mean Halo Bias Function

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        halo_mass_function: HaloMassFunction
        mean_bias: HaloBias

    props: Props = ...
    parent_instance: GObject.Object = ...
    mfp: HaloMassFunction = ...
    mbiasf: HaloBias = ...
    N: typing.Callable[
        [ClusterAbundance, HICosmo, ClusterRedshift, ClusterMass], float
    ] = ...
    intp_d2N: typing.Callable[
        [ClusterAbundance, HICosmo, ClusterRedshift, ClusterMass, float, float], float
    ] = ...
    intp_d2N_bias: typing.Callable[
        [
            ClusterAbundance,
            HICosmo,
            ClusterRedshift,
            ClusterMass,
            float,
            float,
            float,
            float,
        ],
        float,
    ] = ...
    norma: float = ...
    log_norma: float = ...
    lnMi: float = ...
    lnMf: float = ...
    zi: float = ...
    zf: float = ...
    lnM_epsilon: float = ...
    z_epsilon: float = ...
    optimize: bool = ...
    purity: int = ...
    sd_lnM: int = ...
    dbdlnM: NumCosmoMath.Spline2d = ...
    inv_z: NumCosmoMath.Spline = ...
    inv_lnM: NumCosmoMath.Spline = ...
    inv_lnM_z: NumCosmoMath.Spline2d = ...
    rng: int = ...
    ctrl_cosmo: NumCosmoMath.ModelCtrl = ...
    ctrl_reion: NumCosmoMath.ModelCtrl = ...
    ctrl_z: NumCosmoMath.ModelCtrl = ...
    ctrl_m: NumCosmoMath.ModelCtrl = ...
    def __init__(
        self, halo_mass_function: HaloMassFunction = ..., mean_bias: HaloBias = ...
    ) -> None: ...
    @staticmethod
    def clear(cad: ClusterAbundance) -> None: ...
    def d2n(
        self,
        cosmo: HICosmo,
        clusterz: ClusterRedshift,
        clusterm: ClusterMass,
        lnM: float,
        z: float,
    ) -> float: ...
    def free(self) -> None: ...
    def intp_bin_d2n(
        self,
        cosmo: HICosmo,
        clusterz: ClusterRedshift,
        clusterm: ClusterMass,
        lnM_obs_lower: typing.Sequence[float],
        lnM_obs_upper: typing.Sequence[float],
        lnM_obs_params: typing.Optional[typing.Sequence[float]],
        z_obs_lower: typing.Sequence[float],
        z_obs_upper: typing.Sequence[float],
        z_obs_params: typing.Optional[typing.Sequence[float]] = None,
    ) -> float: ...
    def intp_bin_d2n_bias(
        self,
        cosmo: HICosmo,
        clusterz: ClusterRedshift,
        clusterm: ClusterMass,
        lnM_obs_lower: typing.Sequence[float],
        lnM_obs_upper: typing.Sequence[float],
        lnM_obs_params: typing.Optional[typing.Sequence[float]],
        z_obs_lower: typing.Sequence[float],
        z_obs_upper: typing.Sequence[float],
        z_obs_params: typing.Optional[typing.Sequence[float]] = None,
    ) -> float: ...
    def intp_d2n(
        self,
        cosmo: HICosmo,
        clusterz: ClusterRedshift,
        clusterm: ClusterMass,
        lnM: float,
        z: float,
    ) -> float: ...
    def intp_d2n_bias(
        self,
        cosmo: HICosmo,
        clusterz: ClusterRedshift,
        clusterm: ClusterMass,
        lnM_obs: typing.Sequence[float],
        lnM_obs_params: typing.Optional[typing.Sequence[float]],
        z_obs: typing.Sequence[float],
        z_obs_params: typing.Optional[typing.Sequence[float]] = None,
    ) -> float: ...
    def lnM_p_d2n(
        self,
        cosmo: HICosmo,
        clusterz: ClusterRedshift,
        clusterm: ClusterMass,
        lnM_obs: typing.Sequence[float],
        lnM_obs_params: typing.Optional[typing.Sequence[float]],
        z: float,
    ) -> float: ...
    def mean_bias(
        self, cosmo: HICosmo, clusterz: ClusterRedshift, clusterm: ClusterMass
    ) -> float: ...
    def n(
        self, cosmo: HICosmo, clusterz: ClusterRedshift, clusterm: ClusterMass
    ) -> float: ...
    @classmethod
    def new(
        cls, mfp: HaloMassFunction, mbiasf: typing.Optional[HaloBias] = None
    ) -> ClusterAbundance: ...
    @classmethod
    def nodist_new(
        cls, mfp: HaloMassFunction, mbiasf: typing.Optional[HaloBias] = None
    ) -> ClusterAbundance: ...
    def prepare(
        self, cosmo: HICosmo, clusterz: ClusterRedshift, clusterm: ClusterMass
    ) -> None: ...
    def prepare_if_needed(
        self, cosmo: HICosmo, clusterz: ClusterRedshift, clusterm: ClusterMass
    ) -> None: ...
    def prepare_inv_dNdlnM_z(self, cosmo: HICosmo, lnMi: float, z: float) -> None: ...
    def prepare_inv_dNdz(self, cosmo: HICosmo, lnMi: float) -> None: ...
    def ref(self) -> ClusterAbundance: ...
    def set_area(self, area: float) -> None: ...
    def true_n(
        self, cosmo: HICosmo, clusterz: ClusterRedshift, clusterm: ClusterMass
    ) -> float: ...
    def z_p_d2n(
        self,
        cosmo: HICosmo,
        clusterz: ClusterRedshift,
        clusterm: ClusterMass,
        lnM: float,
        z_obs: typing.Sequence[float],
        z_obs_params: typing.Sequence[float],
    ) -> float: ...
    def z_p_lnM_p_d2n(
        self,
        cosmo: HICosmo,
        clusterz: ClusterRedshift,
        clusterm: ClusterMass,
        lnM_obs: typing.Sequence[float],
        lnM_obs_params: typing.Sequence[float],
        z_obs: typing.Sequence[float],
        z_obs_params: typing.Sequence[float],
    ) -> float: ...

class ClusterAbundanceClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        ClusterAbundanceClass()
    """

    parent_class: GObject.ObjectClass = ...

class ClusterAbundanceDataBin(GObject.GPointer): ...
class ClusterAbundanceDataBinM(GObject.GPointer): ...
class ClusterAbundanceDataBinZ(GObject.GPointer): ...
class ClusterAbundanceDataP(GObject.GPointer): ...

class ClusterMass(NumCosmoMath.Model):
    r"""
    :Constructors:

    ::

        ClusterMass(**properties)

    Object NcClusterMass

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: NumCosmoMath.Model = ...
    priv: ClusterMassPrivate = ...
    def __init__(
        self,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def clear(clusterm: ClusterMass) -> None: ...
    def do_N_limits(self, cosmo: HICosmo) -> typing.Tuple[float, float]: ...
    def do_P(
        self,
        cosmo: HICosmo,
        lnM: float,
        z: float,
        lnM_obs: typing.Sequence[float],
        lnM_obs_params: typing.Optional[typing.Sequence[float]] = None,
    ) -> float: ...
    def do_P_bin_limits(
        self,
        cosmo: HICosmo,
        lnM_obs_lower: typing.Sequence[float],
        lnM_obs_upper: typing.Sequence[float],
        lnM_obs_params: typing.Sequence[float],
    ) -> typing.Tuple[float, float]: ...
    def do_P_limits(
        self,
        cosmo: HICosmo,
        lnM_obs: typing.Sequence[float],
        lnM_obs_params: typing.Sequence[float],
    ) -> typing.Tuple[float, float]: ...
    def do_P_vec_z_lnMobs(
        self,
        cosmo: HICosmo,
        lnM: float,
        z: NumCosmoMath.Vector,
        lnM_obs: NumCosmoMath.Matrix,
        lnM_obs_params: NumCosmoMath.Matrix,
        res: typing.Sequence[float],
    ) -> None: ...
    def do_intP(self, cosmo: HICosmo, lnM: float, z: float) -> float: ...
    def do_intP_bin(
        self,
        cosmo: HICosmo,
        lnM: float,
        z: float,
        lnM_obs_lower: typing.Sequence[float],
        lnM_obs_upper: typing.Sequence[float],
        lnM_obs_params: typing.Optional[typing.Sequence[float]] = None,
    ) -> float: ...
    def do_resample(
        self,
        cosmo: HICosmo,
        lnM: float,
        z: float,
        lnM_obs: typing.Sequence[float],
        lnM_obs_params: typing.Sequence[float],
        rng: NumCosmoMath.RNG,
    ) -> bool: ...
    def do_volume(self) -> float: ...
    def free(self) -> None: ...
    @staticmethod
    def id() -> int: ...
    def intp(self, cosmo: HICosmo, lnM: float, z: float) -> float: ...
    def intp_bin(
        self,
        cosmo: HICosmo,
        lnM: float,
        z: float,
        lnM_obs_lower: typing.Sequence[float],
        lnM_obs_upper: typing.Sequence[float],
        lnM_obs_params: typing.Optional[typing.Sequence[float]] = None,
    ) -> float: ...
    @staticmethod
    def log_all_models() -> None: ...
    def n_limits(self, cosmo: HICosmo) -> typing.Tuple[float, float]: ...
    def obs_len(*args, **kwargs): ...  # FIXME Function
    def obs_params_len(*args, **kwargs): ...  # FIXME Function
    def p(
        self,
        cosmo: HICosmo,
        lnM: float,
        z: float,
        lnM_obs: typing.Sequence[float],
        lnM_obs_params: typing.Optional[typing.Sequence[float]] = None,
    ) -> float: ...
    def p_bin_limits(
        self,
        cosmo: HICosmo,
        lnM_obs_lower: typing.Sequence[float],
        lnM_obs_upper: typing.Sequence[float],
        lnM_obs_params: typing.Sequence[float],
    ) -> typing.Tuple[float, float]: ...
    def p_limits(
        self,
        cosmo: HICosmo,
        lnM_obs: typing.Sequence[float],
        lnM_obs_params: typing.Sequence[float],
    ) -> typing.Tuple[float, float]: ...
    def p_vec_z_lnMobs(
        self,
        cosmo: HICosmo,
        lnM: float,
        z: NumCosmoMath.Vector,
        lnM_obs: NumCosmoMath.Matrix,
        lnM_obs_params: NumCosmoMath.Matrix,
        res: typing.Sequence[float],
    ) -> None: ...
    def plcl_Msz_Ml_p_ndetone(
        self,
        lnMcut: float,
        z: float,
        Mpl: float,
        Mcl: float,
        sigma_pl: float,
        sigma_cl: float,
    ) -> float: ...
    def plcl_pdf(
        self,
        lnM_M0: float,
        w1: float,
        w2: float,
        Mobs: typing.Sequence[float],
        Mobs_params: typing.Sequence[float],
    ) -> float: ...
    def plcl_pdf_only_lognormal(
        self, lnM: float, lnMsz_M0: float, lnMl_M0: float
    ) -> float: ...
    def ref(self) -> ClusterMass: ...
    def resample(
        self,
        cosmo: HICosmo,
        lnM: float,
        z: float,
        lnM_obs: typing.Sequence[float],
        lnM_obs_params: typing.Sequence[float],
        rng: NumCosmoMath.RNG,
    ) -> bool: ...
    def volume(self) -> float: ...

class ClusterMassAscaso(ClusterMass):
    r"""
    :Constructors:

    ::

        ClusterMassAscaso(**properties)

    Object NcClusterMassAscaso

    Properties from NcClusterMassAscaso:
      M0 -> gdouble: M0
        Pivot mass
      z0 -> gdouble: z0
        Pivot redshift
      lnRichness-min -> gdouble: lnRichness-min
        Minimum LnRichness
      lnRichness-max -> gdouble: lnRichness-max
        Maximum LnRichness
      mup0 -> gdouble: mup0
        mu_p0
      mup1 -> gdouble: mup1
        mu_p1
      mup2 -> gdouble: mup2
        mu_p2
      sigmap0 -> gdouble: sigmap0
        \sigma_p0
      sigmap1 -> gdouble: sigmap1
        \sigma_p1
      sigmap2 -> gdouble: sigmap2
        \sigma_p2
      mup0-fit -> gboolean: mup0-fit
        mu_p0:fit
      mup1-fit -> gboolean: mup1-fit
        mu_p1:fit
      mup2-fit -> gboolean: mup2-fit
        mu_p2:fit
      sigmap0-fit -> gboolean: sigmap0-fit
        \sigma_p0:fit
      sigmap1-fit -> gboolean: sigmap1-fit
        \sigma_p1:fit
      sigmap2-fit -> gboolean: sigmap2-fit
        \sigma_p2:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        M0: float
        lnRichness_max: float
        lnRichness_min: float
        mup0: float
        mup0_fit: bool
        mup1: float
        mup1_fit: bool
        mup2: float
        mup2_fit: bool
        sigmap0: float
        sigmap0_fit: bool
        sigmap1: float
        sigmap1_fit: bool
        sigmap2: float
        sigmap2_fit: bool
        z0: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: ClusterMass = ...
    priv: ClusterMassAscasoPrivate = ...
    def __init__(
        self,
        M0: float = ...,
        lnRichness_max: float = ...,
        lnRichness_min: float = ...,
        mup0: float = ...,
        mup0_fit: bool = ...,
        mup1: float = ...,
        mup1_fit: bool = ...,
        mup2: float = ...,
        mup2_fit: bool = ...,
        sigmap0: float = ...,
        sigmap0_fit: bool = ...,
        sigmap1: float = ...,
        sigmap1_fit: bool = ...,
        sigmap2: float = ...,
        sigmap2_fit: bool = ...,
        z0: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...

class ClusterMassAscasoClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        ClusterMassAscasoClass()
    """

    parent_class: ClusterMassClass = ...

class ClusterMassAscasoPrivate(GObject.GPointer): ...

class ClusterMassBenson(ClusterMass):
    r"""
    :Constructors:

    ::

        ClusterMassBenson(**properties)

    Object NcClusterMassBenson

    Properties from NcClusterMassBenson:
      signif-obs-min -> gdouble: signif-obs-min
        Minimum obsevational significance
      signif-obs-max -> gdouble: signif-obs-max
        Maximum obsevational significance
      z0 -> gdouble: z0
        Reference redshift
      M0 -> gdouble: M0
        Reference mass
      Asz -> gdouble: Asz
        A_{SZ}
      Bsz -> gdouble: Bsz
        B_{SZ}
      Csz -> gdouble: Csz
        C_{SZ}
      Dsz -> gdouble: Dsz
        D_{SZ}
      Asz-fit -> gboolean: Asz-fit
        A_{SZ}:fit
      Bsz-fit -> gboolean: Bsz-fit
        B_{SZ}:fit
      Csz-fit -> gboolean: Csz-fit
        C_{SZ}:fit
      Dsz-fit -> gboolean: Dsz-fit
        D_{SZ}:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        Asz: float
        Asz_fit: bool
        Bsz: float
        Bsz_fit: bool
        Csz: float
        Csz_fit: bool
        Dsz: float
        Dsz_fit: bool
        M0: float
        signif_obs_max: float
        signif_obs_min: float
        z0: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: ClusterMass = ...
    signif_obs_min: float = ...
    signif_obs_max: float = ...
    z0: float = ...
    M0: float = ...
    def __init__(
        self,
        Asz: float = ...,
        Asz_fit: bool = ...,
        Bsz: float = ...,
        Bsz_fit: bool = ...,
        Csz: float = ...,
        Csz_fit: bool = ...,
        Dsz: float = ...,
        Dsz_fit: bool = ...,
        M0: float = ...,
        signif_obs_max: float = ...,
        signif_obs_min: float = ...,
        z0: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...

class ClusterMassBensonClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        ClusterMassBensonClass()
    """

    parent_class: ClusterMassClass = ...

class ClusterMassBensonXRay(ClusterMassBenson):
    r"""
    :Constructors:

    ::

        ClusterMassBensonXRay(**properties)

    Object NcClusterMassBensonXRay

    Properties from NcClusterMassBensonXRay:
      Yx-obs-min -> gdouble: Yx-obs-min
        Minimum obsevational Yx
      Yx-obs-max -> gdouble: Yx-obs-max
        Maximum obsevational Yx
      M0x -> gdouble: M0x
        X Ray Reference mass
      Y0 -> gdouble: Y0
        Reference Yx
      Ax -> gdouble: Ax
        A_{X}
      Bx -> gdouble: Bx
        B_{X}
      Cx -> gdouble: Cx
        C_{X}
      Dx -> gdouble: Dx
        D_{X}
      Ax-fit -> gboolean: Ax-fit
        A_{X}:fit
      Bx-fit -> gboolean: Bx-fit
        B_{X}:fit
      Cx-fit -> gboolean: Cx-fit
        C_{X}:fit
      Dx-fit -> gboolean: Dx-fit
        D_{X}:fit

    Properties from NcClusterMassBenson:
      signif-obs-min -> gdouble: signif-obs-min
        Minimum obsevational significance
      signif-obs-max -> gdouble: signif-obs-max
        Maximum obsevational significance
      z0 -> gdouble: z0
        Reference redshift
      M0 -> gdouble: M0
        Reference mass
      Asz -> gdouble: Asz
        A_{SZ}
      Bsz -> gdouble: Bsz
        B_{SZ}
      Csz -> gdouble: Csz
        C_{SZ}
      Dsz -> gdouble: Dsz
        D_{SZ}
      Asz-fit -> gboolean: Asz-fit
        A_{SZ}:fit
      Bsz-fit -> gboolean: Bsz-fit
        B_{SZ}:fit
      Csz-fit -> gboolean: Csz-fit
        C_{SZ}:fit
      Dsz-fit -> gboolean: Dsz-fit
        D_{SZ}:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        Ax: float
        Ax_fit: bool
        Bx: float
        Bx_fit: bool
        Cx: float
        Cx_fit: bool
        Dx: float
        Dx_fit: bool
        M0x: float
        Y0: float
        Yx_obs_max: float
        Yx_obs_min: float
        Asz: float
        Asz_fit: bool
        Bsz: float
        Bsz_fit: bool
        Csz: float
        Csz_fit: bool
        Dsz: float
        Dsz_fit: bool
        M0: float
        signif_obs_max: float
        signif_obs_min: float
        z0: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: ClusterMassBenson = ...
    Yx_obs_min: float = ...
    Yx_obs_max: float = ...
    M0x: float = ...
    Y0: float = ...
    def __init__(
        self,
        Ax: float = ...,
        Ax_fit: bool = ...,
        Bx: float = ...,
        Bx_fit: bool = ...,
        Cx: float = ...,
        Cx_fit: bool = ...,
        Dx: float = ...,
        Dx_fit: bool = ...,
        M0x: float = ...,
        Y0: float = ...,
        Yx_obs_max: float = ...,
        Yx_obs_min: float = ...,
        Asz: float = ...,
        Asz_fit: bool = ...,
        Bsz: float = ...,
        Bsz_fit: bool = ...,
        Csz: float = ...,
        Csz_fit: bool = ...,
        Dsz: float = ...,
        Dsz_fit: bool = ...,
        M0: float = ...,
        signif_obs_max: float = ...,
        signif_obs_min: float = ...,
        z0: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...

class ClusterMassBensonXRayClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        ClusterMassBensonXRayClass()
    """

    parent_class: ClusterMassBensonClass = ...

class ClusterMassClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        ClusterMassClass()
    """

    parent_class: NumCosmoMath.ModelClass = ...
    P: typing.Callable[
        [
            ClusterMass,
            HICosmo,
            float,
            float,
            typing.Sequence[float],
            typing.Optional[typing.Sequence[float]],
        ],
        float,
    ] = ...
    intP: typing.Callable[[ClusterMass, HICosmo, float, float], float] = ...
    intP_bin: typing.Callable[
        [
            ClusterMass,
            HICosmo,
            float,
            float,
            typing.Sequence[float],
            typing.Sequence[float],
            typing.Optional[typing.Sequence[float]],
        ],
        float,
    ] = ...
    resample: typing.Callable[
        [
            ClusterMass,
            HICosmo,
            float,
            float,
            typing.Sequence[float],
            typing.Sequence[float],
            NumCosmoMath.RNG,
        ],
        bool,
    ] = ...
    P_limits: typing.Callable[
        [ClusterMass, HICosmo, typing.Sequence[float], typing.Sequence[float]],
        typing.Tuple[float, float],
    ] = ...
    P_bin_limits: typing.Callable[
        [
            ClusterMass,
            HICosmo,
            typing.Sequence[float],
            typing.Sequence[float],
            typing.Sequence[float],
        ],
        typing.Tuple[float, float],
    ] = ...
    N_limits: typing.Callable[[ClusterMass, HICosmo], typing.Tuple[float, float]] = ...
    volume: typing.Callable[[ClusterMass], float] = ...
    P_vec_z_lnMobs: typing.Callable[
        [
            ClusterMass,
            HICosmo,
            float,
            NumCosmoMath.Vector,
            NumCosmoMath.Matrix,
            NumCosmoMath.Matrix,
            typing.Sequence[float],
        ],
        None,
    ] = ...
    _obs_len: int = ...
    _obs_params_len: int = ...
    def obs_len(self) -> int: ...
    def obs_params_len(self) -> int: ...

class ClusterMassLnnormal(ClusterMass):
    r"""
    :Constructors:

    ::

        ClusterMassLnnormal(**properties)

    Object NcClusterMassLnnormal

    Properties from NcClusterMassLnnormal:
      lnMobs-min -> gdouble: lnMobs-min
        Minimum LnMobs
      lnMobs-max -> gdouble: lnMobs-max
        Maximum LnMobs
      bias -> gdouble: bias
        bias
      sigma -> gdouble: sigma
        sigma
      bias-fit -> gboolean: bias-fit
        bias:fit
      sigma-fit -> gboolean: sigma-fit
        sigma:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        bias: float
        bias_fit: bool
        lnMobs_max: float
        lnMobs_min: float
        sigma: float
        sigma_fit: bool
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: ClusterMass = ...
    lnMobs_max: float = ...
    lnMobs_min: float = ...
    def __init__(
        self,
        bias: float = ...,
        bias_fit: bool = ...,
        lnMobs_max: float = ...,
        lnMobs_min: float = ...,
        sigma: float = ...,
        sigma_fit: bool = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...

class ClusterMassLnnormalClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        ClusterMassLnnormalClass()
    """

    parent_class: ClusterMassClass = ...

class ClusterMassNodist(ClusterMass):
    r"""
    :Constructors:

    ::

        ClusterMassNodist(**properties)

    Object NcClusterMassNodist

    Properties from NcClusterMassNodist:
      lnM-min -> gdouble: lnM-min
        Minimum mass
      lnM-max -> gdouble: lnM-max
        Maximum mass

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        lnM_max: float
        lnM_min: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: ClusterMass = ...
    priv: ClusterMassNodistPrivate = ...
    def __init__(
        self,
        lnM_max: float = ...,
        lnM_min: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...

class ClusterMassNodistClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        ClusterMassNodistClass()
    """

    parent_class: ClusterMassClass = ...

class ClusterMassNodistPrivate(GObject.GPointer): ...

class ClusterMassPlCL(ClusterMass):
    r"""
    :Constructors:

    ::

        ClusterMassPlCL(**properties)

    Object NcClusterMassPlCL

    Properties from NcClusterMassPlCL:
      M0 -> gdouble: M0
        Reference mass
      Asz -> gdouble: Asz
        \alpha_{SZ}
      Bsz -> gdouble: Bsz
        b_{SZ}
      sigma-sz -> gdouble: sigma-sz
        \sigma_{SZ}
      Al -> gdouble: Al
        \alpha_{L}
      Bl -> gdouble: Bl
        b_{L}
      sigma-l -> gdouble: sigma-l
        \sigma_{L}
      cor -> gdouble: cor
        \rho
      Asz-fit -> gboolean: Asz-fit
        \alpha_{SZ}:fit
      Bsz-fit -> gboolean: Bsz-fit
        b_{SZ}:fit
      sigma-sz-fit -> gboolean: sigma-sz-fit
        \sigma_{SZ}:fit
      Al-fit -> gboolean: Al-fit
        \alpha_{L}:fit
      Bl-fit -> gboolean: Bl-fit
        b_{L}:fit
      sigma-l-fit -> gboolean: sigma-l-fit
        \sigma_{L}:fit
      cor-fit -> gboolean: cor-fit
        \rho:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        Al: float
        Al_fit: bool
        Asz: float
        Asz_fit: bool
        Bl: float
        Bl_fit: bool
        Bsz: float
        Bsz_fit: bool
        M0: float
        cor: float
        cor_fit: bool
        sigma_l: float
        sigma_l_fit: bool
        sigma_sz: float
        sigma_sz_fit: bool
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: ClusterMass = ...
    M0: float = ...
    T: int = ...
    s: int = ...
    def __init__(
        self,
        Al: float = ...,
        Al_fit: bool = ...,
        Asz: float = ...,
        Asz_fit: bool = ...,
        Bl: float = ...,
        Bl_fit: bool = ...,
        Bsz: float = ...,
        Bsz_fit: bool = ...,
        M0: float = ...,
        cor: float = ...,
        cor_fit: bool = ...,
        sigma_l: float = ...,
        sigma_l_fit: bool = ...,
        sigma_sz: float = ...,
        sigma_sz_fit: bool = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def gsl_f(
        p: float,
        hx: float,
        n: int,
        mszl: ClusterMassPlCL,
        lnM: float,
        Mobs: typing.Sequence[float],
        Mobs_params: typing.Sequence[float],
    ) -> None: ...
    @staticmethod
    def peak_new_variables(
        N: float,
        lb: float,
        ub: float,
        mszl: ClusterMassPlCL,
        lnM: float,
        Mobs: typing.Sequence[float],
        Mobs_params: typing.Sequence[float],
    ) -> None: ...

class ClusterMassPlCLClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        ClusterMassPlCLClass()
    """

    parent_class: ClusterMassClass = ...

class ClusterMassPrivate(GObject.GPointer): ...

class ClusterMassVanderlinde(ClusterMass):
    r"""
    :Constructors:

    ::

        ClusterMassVanderlinde(**properties)

    Object NcClusterMassVanderlinde

    Properties from NcClusterMassVanderlinde:
      signif-obs-min -> gdouble: signif-obs-min
        Minimum observational significance
      signif-obs-max -> gdouble: signif-obs-max
        Maximum observational significance
      z0 -> gdouble: z0
        Reference redshift
      M0 -> gdouble: M0
        Reference mass
      Asz -> gdouble: Asz
        A_{SZ}
      Bsz -> gdouble: Bsz
        B_{SZ}
      Csz -> gdouble: Csz
        C_{SZ}
      Dsz -> gdouble: Dsz
        D_{SZ}
      Asz-fit -> gboolean: Asz-fit
        A_{SZ}:fit
      Bsz-fit -> gboolean: Bsz-fit
        B_{SZ}:fit
      Csz-fit -> gboolean: Csz-fit
        C_{SZ}:fit
      Dsz-fit -> gboolean: Dsz-fit
        D_{SZ}:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        Asz: float
        Asz_fit: bool
        Bsz: float
        Bsz_fit: bool
        Csz: float
        Csz_fit: bool
        Dsz: float
        Dsz_fit: bool
        M0: float
        signif_obs_max: float
        signif_obs_min: float
        z0: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: ClusterMass = ...
    signif_obs_min: float = ...
    signif_obs_max: float = ...
    z0: float = ...
    M0: float = ...
    def __init__(
        self,
        Asz: float = ...,
        Asz_fit: bool = ...,
        Bsz: float = ...,
        Bsz_fit: bool = ...,
        Csz: float = ...,
        Csz_fit: bool = ...,
        Dsz: float = ...,
        Dsz_fit: bool = ...,
        M0: float = ...,
        signif_obs_max: float = ...,
        signif_obs_min: float = ...,
        z0: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...

class ClusterMassVanderlindeClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        ClusterMassVanderlindeClass()
    """

    parent_class: ClusterMassClass = ...

class ClusterPhotozGauss(ClusterRedshift):
    r"""
    :Constructors:

    ::

        ClusterPhotozGauss(**properties)
        new() -> NumCosmo.ClusterRedshift

    Object NcClusterPhotozGauss

    Properties from NcClusterPhotozGauss:
      pz-min -> gdouble: pz-min
        Minimum photoz
      pz-max -> gdouble: pz-max
        Maximum photoz

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        pz_max: float
        pz_min: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: ClusterRedshift = ...
    pz_max: float = ...
    pz_min: float = ...
    def __init__(
        self,
        pz_max: float = ...,
        pz_min: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @classmethod
    def new(cls) -> ClusterPhotozGauss: ...

class ClusterPhotozGaussClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        ClusterPhotozGaussClass()
    """

    parent_class: ClusterRedshiftClass = ...

class ClusterPhotozGaussGlobal(ClusterRedshift):
    r"""
    :Constructors:

    ::

        ClusterPhotozGaussGlobal(**properties)
        new(pz_min:float, pz_max:float, z_bias:float, sigma0:float) -> NumCosmo.ClusterRedshift

    Object NcClusterPhotozGaussGlobal

    Properties from NcClusterPhotozGaussGlobal:
      pz-min -> gdouble: pz-min
        Minimum photoz
      pz-max -> gdouble: pz-max
        Maximum photoz
      z-bias -> gdouble: z-bias
        z-bias
      sigma0 -> gdouble: sigma0
        sigma0
      z-bias-fit -> gboolean: z-bias-fit
        z-bias:fit
      sigma0-fit -> gboolean: sigma0-fit
        sigma0:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        pz_max: float
        pz_min: float
        sigma0: float
        sigma0_fit: bool
        z_bias: float
        z_bias_fit: bool
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: ClusterRedshift = ...
    pz_min: float = ...
    pz_max: float = ...
    def __init__(
        self,
        pz_max: float = ...,
        pz_min: float = ...,
        sigma0: float = ...,
        sigma0_fit: bool = ...,
        z_bias: float = ...,
        z_bias_fit: bool = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    def get_sigma0(self) -> float: ...
    def get_z_bias(self) -> float: ...
    @classmethod
    def new(
        cls, pz_min: float, pz_max: float, z_bias: float, sigma0: float
    ) -> ClusterPhotozGaussGlobal: ...
    def set_sigma0(self, sigma0: float) -> None: ...
    def set_z_bias(self, z_bias: float) -> None: ...

class ClusterPhotozGaussGlobalClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        ClusterPhotozGaussGlobalClass()
    """

    parent_class: ClusterRedshiftClass = ...

class ClusterPseudoCounts(NumCosmoMath.Model):
    r"""
    :Constructors:

    ::

        ClusterPseudoCounts(**properties)
        new(nclusters:int) -> NumCosmo.ClusterPseudoCounts

    Object NcClusterPseudoCounts

    Properties from NcClusterPseudoCounts:
      number-clusters -> guint: number-clusters
        Number of clusters
      lnMCut -> gdouble: lnMCut
        \ln{M_{CUT}}
      sigma-Mcut -> gdouble: sigma-Mcut
        \sigma_{MCUT}
      zmin -> gdouble: zmin
        z_{min}
      Deltaz -> gdouble: Deltaz
        \delta{}z
      lnMCut-fit -> gboolean: lnMCut-fit
        \ln{M_{CUT}}:fit
      sigma-Mcut-fit -> gboolean: sigma-Mcut-fit
        \sigma_{MCUT}:fit
      zmin-fit -> gboolean: zmin-fit
        z_{min}:fit
      Deltaz-fit -> gboolean: Deltaz-fit
        \delta{}z:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        Deltaz: float
        Deltaz_fit: bool
        lnMCut: float
        lnMCut_fit: bool
        number_clusters: int
        sigma_Mcut: float
        sigma_Mcut_fit: bool
        zmin: float
        zmin_fit: bool
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: NumCosmoMath.Model = ...
    nclusters: int = ...
    T: int = ...
    s: int = ...
    workz: float = ...
    def __init__(
        self,
        Deltaz: float = ...,
        Deltaz_fit: bool = ...,
        lnMCut: float = ...,
        lnMCut_fit: bool = ...,
        number_clusters: int = ...,
        sigma_Mcut: float = ...,
        sigma_Mcut_fit: bool = ...,
        zmin: float = ...,
        zmin_fit: bool = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def clear(cpc: ClusterPseudoCounts) -> None: ...
    def free(self) -> None: ...
    @staticmethod
    def id() -> int: ...
    def mf_lognormal_integral(
        self,
        mfp: HaloMassFunction,
        clusterm: ClusterMass,
        cosmo: HICosmo,
        lnMsz: float,
        lnMl: float,
        z: float,
    ) -> float: ...
    def ndet(self, mfp: HaloMassFunction, cosmo: HICosmo) -> float: ...
    def ndet_no_z_integral(self, cosmo: HICosmo, z: float) -> float: ...
    @classmethod
    def new(cls, nclusters: int) -> ClusterPseudoCounts: ...
    def posterior_ndetone(
        self,
        mfp: HaloMassFunction,
        cosmo: HICosmo,
        clusterm: ClusterMass,
        z: float,
        Mpl: float,
        Mcl: float,
        sigma_pl: float,
        sigma_cl: float,
    ) -> float: ...
    def posterior_numerator(
        self,
        mfp: HaloMassFunction,
        clusterm: ClusterMass,
        cosmo: HICosmo,
        z: float,
        Mobs: typing.Sequence[float],
        Mobs_params: typing.Sequence[float],
    ) -> float: ...
    def posterior_numerator_plcl(
        self,
        mfp: HaloMassFunction,
        clusterm: ClusterMass,
        cosmo: HICosmo,
        z: float,
        Mpl: float,
        Mcl: float,
        sigma_pl: float,
        sigma_cl: float,
    ) -> float: ...
    def ref(self) -> ClusterPseudoCounts: ...
    def selection_function(self, lnM: float, z: float) -> float: ...
    def selection_function_lnMi(self, cosmo: HICosmo) -> float: ...

class ClusterPseudoCountsClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        ClusterPseudoCountsClass()
    """

    parent_class: NumCosmoMath.ModelClass = ...

class ClusterRedshift(NumCosmoMath.Model):
    r"""
    :Constructors:

    ::

        ClusterRedshift(**properties)

    Object NcClusterRedshift

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: NumCosmoMath.Model = ...
    priv: ClusterRedshiftPrivate = ...
    def __init__(
        self,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def clear(clusterz: ClusterRedshift) -> None: ...
    def do_N_limits(self, cosmo: HICosmo, z_lower: float, z_upper: float) -> None: ...
    def do_P(
        self, cosmo: HICosmo, lnM: float, z: float, z_obs: float, z_obs_params: float
    ) -> float: ...
    def do_P_bin_limits(
        self,
        cosmo: HICosmo,
        z_obs_lower: float,
        z_obs_upper: float,
        z_obs_params: float,
        z_lower: float,
        z_upper: float,
    ) -> None: ...
    def do_P_limits(
        self,
        cosmo: HICosmo,
        z_obs: float,
        z_obs_params: float,
        z_lower: float,
        z_upper: float,
    ) -> None: ...
    def do_intP(self, cosmo: HICosmo, lnM: float, z: float) -> float: ...
    def do_intP_bin(
        self,
        cosmo: HICosmo,
        lnM: float,
        z: float,
        z_obs_lower: float,
        z_obs_upper: float,
        z_obs_params: float,
    ) -> float: ...
    def do_resample(
        self, cosmo: HICosmo, lnM: float, z: float, rng: NumCosmoMath.RNG
    ) -> typing.Tuple[bool, float, float]: ...
    def do_volume(self) -> float: ...
    def free(self) -> None: ...
    @staticmethod
    def id() -> int: ...
    def intp(self, cosmo: HICosmo, lnM: float, z: float) -> float: ...
    def intp_bin(
        self,
        cosmo: HICosmo,
        lnM: float,
        z: float,
        z_obs_lower: typing.Sequence[float],
        z_obs_upper: typing.Sequence[float],
        z_obs_params: typing.Optional[typing.Sequence[float]] = None,
    ) -> float: ...
    @staticmethod
    def log_all_models() -> None: ...
    def n_limits(self, cosmo: HICosmo) -> typing.Tuple[float, float]: ...
    def obs_len(*args, **kwargs): ...  # FIXME Function
    def obs_params_len(*args, **kwargs): ...  # FIXME Function
    def p(
        self,
        cosmo: HICosmo,
        lnM: float,
        z: float,
        z_obs: typing.Sequence[float],
        z_obs_params: typing.Sequence[float],
    ) -> float: ...
    def p_bin_limits(
        self,
        cosmo: HICosmo,
        z_obs_lower: typing.Sequence[float],
        z_obs_upper: typing.Sequence[float],
        z_obs_params: typing.Sequence[float],
    ) -> typing.Tuple[float, float]: ...
    def p_limits(
        self,
        cosmo: HICosmo,
        z_obs: typing.Sequence[float],
        z_obs_params: typing.Sequence[float],
    ) -> typing.Tuple[float, float]: ...
    def ref(self) -> ClusterRedshift: ...
    def resample(
        self, cosmo: HICosmo, lnM: float, z: float, rng: NumCosmoMath.RNG
    ) -> typing.Tuple[bool, float, float]: ...
    def volume(self) -> float: ...

class ClusterRedshiftClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        ClusterRedshiftClass()
    """

    parent_class: NumCosmoMath.ModelClass = ...
    P: typing.Callable[
        [ClusterRedshift, HICosmo, float, float, float, float], float
    ] = ...
    intP: typing.Callable[[ClusterRedshift, HICosmo, float, float], float] = ...
    intP_bin: typing.Callable[
        [ClusterRedshift, HICosmo, float, float, float, float, float], float
    ] = ...
    resample: typing.Callable[
        [ClusterRedshift, HICosmo, float, float, NumCosmoMath.RNG],
        typing.Tuple[bool, float, float],
    ] = ...
    P_limits: typing.Callable[
        [ClusterRedshift, HICosmo, float, float, float, float], None
    ] = ...
    P_bin_limits: typing.Callable[
        [ClusterRedshift, HICosmo, float, float, float, float, float], None
    ] = ...
    N_limits: typing.Callable[[ClusterRedshift, HICosmo, float, float], None] = ...
    volume: typing.Callable[[ClusterRedshift], float] = ...
    _obs_len: int = ...
    _obs_params_len: int = ...
    def obs_len(self) -> int: ...
    def obs_params_len(self) -> int: ...

class ClusterRedshiftNodist(ClusterRedshift):
    r"""
    :Constructors:

    ::

        ClusterRedshiftNodist(**properties)

    Object NcClusterRedshiftNodist

    Properties from NcClusterRedshiftNodist:
      z-min -> gdouble: z-min
        Minimum z
      z-max -> gdouble: z-max
        Maximum z

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        z_max: float
        z_min: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: ClusterRedshift = ...
    priv: ClusterRedshiftNodistPrivate = ...
    def __init__(
        self,
        z_max: float = ...,
        z_min: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...

class ClusterRedshiftNodistClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        ClusterRedshiftNodistClass()
    """

    parent_class: ClusterRedshiftClass = ...

class ClusterRedshiftNodistPrivate(GObject.GPointer): ...
class ClusterRedshiftPrivate(GObject.GPointer): ...

class CorClusterCmbLensLimber(GObject.Object):
    r"""
    :Constructors:

    ::

        CorClusterCmbLensLimber(**properties)
        new() -> NumCosmo.CorClusterCmbLensLimber

    Object NcCorClusterCmbLensLimber

    Signals from GObject:
      notify (GParam)
    """

    parent_instance: GObject.Object = ...
    oneh_int_mass_spline: NumCosmoMath.Spline = ...
    @classmethod
    def new(cls) -> CorClusterCmbLensLimber: ...
    def oneh_int_mass(
        self,
        cad: ClusterAbundance,
        clusterm: ClusterMass,
        cosmo: HICosmo,
        dp: HaloDensityProfile,
        k: float,
        z: float,
        lnM_obs: typing.Sequence[float],
        lnM_obs_params: typing.Sequence[float],
    ) -> float: ...
    def oneh_term(
        self,
        cad: ClusterAbundance,
        cosmo: HICosmo,
        dist: Distance,
        dp: HaloDensityProfile,
        l: int,
        lnM_obs: typing.Sequence[float],
        lnM_obs_params: typing.Sequence[float],
        z_obs: typing.Sequence[float],
        z_obs_params: typing.Sequence[float],
    ) -> float: ...
    def twoh_int_mass1(
        self, cad: ClusterAbundance, clusterm: ClusterMass, cosmo: HICosmo, z: float
    ) -> float: ...
    def twoh_int_mass2(
        self,
        cad: ClusterAbundance,
        clusterm: ClusterMass,
        cosmo: HICosmo,
        dp: HaloDensityProfile,
        k: float,
        z: float,
    ) -> float: ...
    def twoh_int_mm(
        self,
        cad: ClusterAbundance,
        cosmo: HICosmo,
        dp: HaloDensityProfile,
        k: float,
        z: float,
    ) -> float: ...
    def twoh_term(
        self,
        cad: ClusterAbundance,
        cosmo: HICosmo,
        dist: Distance,
        dp: HaloDensityProfile,
        l: int,
        z_obs: typing.Sequence[float],
        z_obs_params: typing.Sequence[float],
    ) -> float: ...

class CorClusterCmbLensLimberClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        CorClusterCmbLensLimberClass()
    """

    parent_class: GObject.ObjectClass = ...

class DECont(NumCosmoMath.CSQ1D):
    r"""
    :Constructors:

    ::

        DECont(**properties)
        new(Omegaw:float, OmegaL:float, cs2:float, w:float) -> NumCosmo.DECont

    Object NcDECont

    Properties from NcDECont:
      Omegaw -> gdouble: Omegaw
        \Omega_w
      OmegaL -> gdouble: OmegaL
        \Omega_\Lambda
      cs2 -> gdouble: cs2
        c_s^2
      w -> gdouble: w
        w
      k -> gdouble: k
        k

    Properties from NcmCSQ1D:
      reltol -> gdouble: reltol
        Relative tolerance
      abstol -> gdouble: abstol
        Absolute tolerance tolerance
      ti -> gdouble: ti
        The initial time t_i
      tf -> gdouble: tf
        The final time t_f
      adiab-threshold -> gdouble: adiab-threshold
        The adiabatic threshold
      prop-threshold -> gdouble: prop-threshold
        The propagator threshold
      save-evol -> gboolean: save-evol
        Save the system evolution
      vacuum-type -> NcmCSQ1DInitialStateType: vacuum-type
        The vacuum type
      vacuum-reltol -> gdouble: vacuum-reltol
        The vacuum relative tolerance
      vacuum-max-time -> gdouble: vacuum-max-time
        The vacuum maximum time

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        OmegaL: float
        Omegaw: float
        cs2: float
        k: float
        w: float
        abstol: float
        adiab_threshold: float
        prop_threshold: float
        reltol: float
        save_evol: bool
        tf: float
        ti: float
        vacuum_max_time: float
        vacuum_reltol: float
        vacuum_type: NumCosmoMath.CSQ1DInitialStateType

    props: Props = ...
    def __init__(
        self,
        OmegaL: float = ...,
        Omegaw: float = ...,
        cs2: float = ...,
        k: float = ...,
        w: float = ...,
        abstol: float = ...,
        adiab_threshold: float = ...,
        prop_threshold: float = ...,
        reltol: float = ...,
        save_evol: bool = ...,
        tf: float = ...,
        ti: float = ...,
        vacuum_max_time: float = ...,
        vacuum_reltol: float = ...,
        vacuum_type: NumCosmoMath.CSQ1DInitialStateType = ...,
    ) -> None: ...
    @staticmethod
    def clear(dec: DECont) -> None: ...
    def free(self) -> None: ...
    def get_k(self) -> float: ...
    @classmethod
    def new(cls, Omegaw: float, OmegaL: float, cs2: float, w: float) -> DECont: ...
    def ref(self) -> DECont: ...
    def set_k(self, k: float) -> None: ...

class DEContClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DEContClass()
    """

    parent_class: NumCosmoMath.CSQ1DClass = ...

class DataBaoA(NumCosmoMath.DataGaussDiag):
    r"""
    :Constructors:

    ::

        DataBaoA(**properties)
        new_from_file(filename:str) -> NumCosmo.DataBaoA
        new_from_id(dist:NumCosmo.Distance, id:NumCosmo.DataBaoId) -> NumCosmo.DataBaoA

    Object NcDataBaoA

    Properties from NcDataBaoA:
      dist -> NcDistance: dist
        Distance object
      z -> NcmVector: z
        Data redshift

    Properties from NcmDataGaussDiag:
      n-points -> guint: n-points
        Data sample size
      w-mean -> gboolean: w-mean
        Whether to minimize analytically over the weighted mean
      mean -> NcmVector: mean
        Data mean
      sigma -> NcmVector: sigma
        Data standard deviation

    Properties from NcmData:
      name -> gchararray: name
        Data type name
      desc -> gchararray: desc
        Data description
      long-desc -> gchararray: long-desc
        Data detailed description
      init -> gboolean: init
        Data initialized state
      bootstrap -> NcmBootstrap: bootstrap
        Data bootstrap object

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        dist: Distance
        z: NumCosmoMath.Vector
        mean: NumCosmoMath.Vector
        n_points: int
        sigma: NumCosmoMath.Vector
        w_mean: bool
        bootstrap: NumCosmoMath.Bootstrap
        desc: str
        init: bool
        long_desc: str
        name: str

    props: Props = ...
    def __init__(
        self,
        dist: Distance = ...,
        z: NumCosmoMath.Vector = ...,
        mean: NumCosmoMath.Vector = ...,
        n_points: int = ...,
        sigma: NumCosmoMath.Vector = ...,
        w_mean: bool = ...,
        bootstrap: NumCosmoMath.Bootstrap = ...,
        desc: str = ...,
        init: bool = ...,
        long_desc: str = ...,
    ) -> None: ...
    @classmethod
    def new_from_file(cls, filename: str) -> DataBaoA: ...
    @classmethod
    def new_from_id(cls, dist: Distance, id: DataBaoId) -> DataBaoA: ...
    def set_dist(self, dist: Distance) -> None: ...

class DataBaoAClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DataBaoAClass()
    """

    parent_class: NumCosmoMath.DataGaussDiagClass = ...

class DataBaoDHrDAr(NumCosmoMath.DataGaussCov):
    r"""
    :Constructors:

    ::

        DataBaoDHrDAr(**properties)
        new_from_file(filename:str) -> NumCosmo.DataBaoDHrDAr
        new_from_id(dist:NumCosmo.Distance, id:NumCosmo.DataBaoId) -> NumCosmo.DataBaoDHrDAr

    Object NcDataBaoDHrDAr

    Properties from NcDataBaoDHrDAr:
      dist -> NcDistance: dist
        Distance object
      z -> NcmVector: z
        Data redshift

    Properties from NcmDataGaussCov:
      n-points -> guint: n-points
        Data sample size
      use-norma -> gboolean: use-norma
        Use the likelihood normalization to calculate -2lnL
      mean -> NcmVector: mean
        Data mean
      cov -> NcmMatrix: cov
        Data covariance

    Properties from NcmData:
      name -> gchararray: name
        Data type name
      desc -> gchararray: desc
        Data description
      long-desc -> gchararray: long-desc
        Data detailed description
      init -> gboolean: init
        Data initialized state
      bootstrap -> NcmBootstrap: bootstrap
        Data bootstrap object

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        dist: Distance
        z: NumCosmoMath.Vector
        cov: NumCosmoMath.Matrix
        mean: NumCosmoMath.Vector
        n_points: int
        use_norma: bool
        bootstrap: NumCosmoMath.Bootstrap
        desc: str
        init: bool
        long_desc: str
        name: str

    props: Props = ...
    def __init__(
        self,
        dist: Distance = ...,
        z: NumCosmoMath.Vector = ...,
        cov: NumCosmoMath.Matrix = ...,
        mean: NumCosmoMath.Vector = ...,
        n_points: int = ...,
        use_norma: bool = ...,
        bootstrap: NumCosmoMath.Bootstrap = ...,
        desc: str = ...,
        init: bool = ...,
        long_desc: str = ...,
    ) -> None: ...
    @classmethod
    def new_from_file(cls, filename: str) -> DataBaoDHrDAr: ...
    @classmethod
    def new_from_id(cls, dist: Distance, id: DataBaoId) -> DataBaoDHrDAr: ...
    def set_dist(self, dist: Distance) -> None: ...

class DataBaoDHrDArClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DataBaoDHrDArClass()
    """

    parent_class: NumCosmoMath.DataGaussCovClass = ...

class DataBaoDMrHr(NumCosmoMath.DataGaussCov):
    r"""
    :Constructors:

    ::

        DataBaoDMrHr(**properties)
        new_from_file(filename:str) -> NumCosmo.DataBaoDMrHr
        new_from_id(dist:NumCosmo.Distance, id:NumCosmo.DataBaoId) -> NumCosmo.DataBaoDMrHr

    Object NcDataBaoDMrHr

    Properties from NcDataBaoDMrHr:
      dist -> NcDistance: dist
        Distance object
      z -> NcmVector: z
        Data redshift
      rs-fiduc -> gdouble: rs-fiduc
        r_s fiducial

    Properties from NcmDataGaussCov:
      n-points -> guint: n-points
        Data sample size
      use-norma -> gboolean: use-norma
        Use the likelihood normalization to calculate -2lnL
      mean -> NcmVector: mean
        Data mean
      cov -> NcmMatrix: cov
        Data covariance

    Properties from NcmData:
      name -> gchararray: name
        Data type name
      desc -> gchararray: desc
        Data description
      long-desc -> gchararray: long-desc
        Data detailed description
      init -> gboolean: init
        Data initialized state
      bootstrap -> NcmBootstrap: bootstrap
        Data bootstrap object

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        dist: Distance
        rs_fiduc: float
        z: NumCosmoMath.Vector
        cov: NumCosmoMath.Matrix
        mean: NumCosmoMath.Vector
        n_points: int
        use_norma: bool
        bootstrap: NumCosmoMath.Bootstrap
        desc: str
        init: bool
        long_desc: str
        name: str

    props: Props = ...
    def __init__(
        self,
        dist: Distance = ...,
        rs_fiduc: float = ...,
        z: NumCosmoMath.Vector = ...,
        cov: NumCosmoMath.Matrix = ...,
        mean: NumCosmoMath.Vector = ...,
        n_points: int = ...,
        use_norma: bool = ...,
        bootstrap: NumCosmoMath.Bootstrap = ...,
        desc: str = ...,
        init: bool = ...,
        long_desc: str = ...,
    ) -> None: ...
    @classmethod
    def new_from_file(cls, filename: str) -> DataBaoDMrHr: ...
    @classmethod
    def new_from_id(cls, dist: Distance, id: DataBaoId) -> DataBaoDMrHr: ...
    def set_dist(self, dist: Distance) -> None: ...

class DataBaoDMrHrClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DataBaoDMrHrClass()
    """

    parent_class: NumCosmoMath.DataGaussCovClass = ...

class DataBaoDV(NumCosmoMath.DataGaussDiag):
    r"""
    :Constructors:

    ::

        DataBaoDV(**properties)
        new_from_file(filename:str) -> NumCosmo.DataBaoDV
        new_from_id(dist:NumCosmo.Distance, id:NumCosmo.DataBaoId) -> NumCosmo.DataBaoDV

    Object NcDataBaoDV

    Properties from NcDataBaoDV:
      dist -> NcDistance: dist
        Distance object
      z -> NcmVector: z
        Data redshift

    Properties from NcmDataGaussDiag:
      n-points -> guint: n-points
        Data sample size
      w-mean -> gboolean: w-mean
        Whether to minimize analytically over the weighted mean
      mean -> NcmVector: mean
        Data mean
      sigma -> NcmVector: sigma
        Data standard deviation

    Properties from NcmData:
      name -> gchararray: name
        Data type name
      desc -> gchararray: desc
        Data description
      long-desc -> gchararray: long-desc
        Data detailed description
      init -> gboolean: init
        Data initialized state
      bootstrap -> NcmBootstrap: bootstrap
        Data bootstrap object

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        dist: Distance
        z: NumCosmoMath.Vector
        mean: NumCosmoMath.Vector
        n_points: int
        sigma: NumCosmoMath.Vector
        w_mean: bool
        bootstrap: NumCosmoMath.Bootstrap
        desc: str
        init: bool
        long_desc: str
        name: str

    props: Props = ...
    def __init__(
        self,
        dist: Distance = ...,
        z: NumCosmoMath.Vector = ...,
        mean: NumCosmoMath.Vector = ...,
        n_points: int = ...,
        sigma: NumCosmoMath.Vector = ...,
        w_mean: bool = ...,
        bootstrap: NumCosmoMath.Bootstrap = ...,
        desc: str = ...,
        init: bool = ...,
        long_desc: str = ...,
    ) -> None: ...
    @classmethod
    def new_from_file(cls, filename: str) -> DataBaoDV: ...
    @classmethod
    def new_from_id(cls, dist: Distance, id: DataBaoId) -> DataBaoDV: ...
    def set_dist(self, dist: Distance) -> None: ...

class DataBaoDVClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DataBaoDVClass()
    """

    parent_class: NumCosmoMath.DataGaussDiagClass = ...

class DataBaoDVDV(NumCosmoMath.DataGaussDiag):
    r"""
    :Constructors:

    ::

        DataBaoDVDV(**properties)
        new_from_file(filename:str) -> NumCosmo.DataBaoDVDV
        new_from_id(dist:NumCosmo.Distance, id:NumCosmo.DataBaoId) -> NumCosmo.DataBaoDVDV

    Object NcDataBaoDVDV

    Properties from NcDataBaoDVDV:
      dist -> NcDistance: dist
        Distance object

    Properties from NcmDataGaussDiag:
      n-points -> guint: n-points
        Data sample size
      w-mean -> gboolean: w-mean
        Whether to minimize analytically over the weighted mean
      mean -> NcmVector: mean
        Data mean
      sigma -> NcmVector: sigma
        Data standard deviation

    Properties from NcmData:
      name -> gchararray: name
        Data type name
      desc -> gchararray: desc
        Data description
      long-desc -> gchararray: long-desc
        Data detailed description
      init -> gboolean: init
        Data initialized state
      bootstrap -> NcmBootstrap: bootstrap
        Data bootstrap object

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        dist: Distance
        mean: NumCosmoMath.Vector
        n_points: int
        sigma: NumCosmoMath.Vector
        w_mean: bool
        bootstrap: NumCosmoMath.Bootstrap
        desc: str
        init: bool
        long_desc: str
        name: str

    props: Props = ...
    def __init__(
        self,
        dist: Distance = ...,
        mean: NumCosmoMath.Vector = ...,
        n_points: int = ...,
        sigma: NumCosmoMath.Vector = ...,
        w_mean: bool = ...,
        bootstrap: NumCosmoMath.Bootstrap = ...,
        desc: str = ...,
        init: bool = ...,
        long_desc: str = ...,
    ) -> None: ...
    @classmethod
    def new_from_file(cls, filename: str) -> DataBaoDVDV: ...
    @classmethod
    def new_from_id(cls, dist: Distance, id: DataBaoId) -> DataBaoDVDV: ...
    def set_dist(self, dist: Distance) -> None: ...

class DataBaoDVDVClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DataBaoDVDVClass()
    """

    parent_class: NumCosmoMath.DataGaussDiagClass = ...

class DataBaoDtrDHr(NumCosmoMath.DataGaussCov):
    r"""
    :Constructors:

    ::

        DataBaoDtrDHr(**properties)
        new_from_file(filename:str) -> NumCosmo.DataBaoDtrDHr
        new_from_id(dist:NumCosmo.Distance, id:NumCosmo.DataBaoId) -> NumCosmo.DataBaoDtrDHr

    Object NcDataBaoDtrDHr

    Properties from NcDataBaoDtrDHr:
      dist -> NcDistance: dist
        Distance object
      z -> NcmVector: z
        Data redshift

    Properties from NcmDataGaussCov:
      n-points -> guint: n-points
        Data sample size
      use-norma -> gboolean: use-norma
        Use the likelihood normalization to calculate -2lnL
      mean -> NcmVector: mean
        Data mean
      cov -> NcmMatrix: cov
        Data covariance

    Properties from NcmData:
      name -> gchararray: name
        Data type name
      desc -> gchararray: desc
        Data description
      long-desc -> gchararray: long-desc
        Data detailed description
      init -> gboolean: init
        Data initialized state
      bootstrap -> NcmBootstrap: bootstrap
        Data bootstrap object

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        dist: Distance
        z: NumCosmoMath.Vector
        cov: NumCosmoMath.Matrix
        mean: NumCosmoMath.Vector
        n_points: int
        use_norma: bool
        bootstrap: NumCosmoMath.Bootstrap
        desc: str
        init: bool
        long_desc: str
        name: str

    props: Props = ...
    def __init__(
        self,
        dist: Distance = ...,
        z: NumCosmoMath.Vector = ...,
        cov: NumCosmoMath.Matrix = ...,
        mean: NumCosmoMath.Vector = ...,
        n_points: int = ...,
        use_norma: bool = ...,
        bootstrap: NumCosmoMath.Bootstrap = ...,
        desc: str = ...,
        init: bool = ...,
        long_desc: str = ...,
    ) -> None: ...
    @classmethod
    def new_from_file(cls, filename: str) -> DataBaoDtrDHr: ...
    @classmethod
    def new_from_id(cls, dist: Distance, id: DataBaoId) -> DataBaoDtrDHr: ...
    def set_dist(self, dist: Distance) -> None: ...

class DataBaoDtrDHrClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DataBaoDtrDHrClass()
    """

    parent_class: NumCosmoMath.DataGaussCovClass = ...

class DataBaoDvrDtDh(NumCosmoMath.DataGaussCov):
    r"""
    :Constructors:

    ::

        DataBaoDvrDtDh(**properties)
        new_from_file(filename:str) -> NumCosmo.DataBaoDvrDtDh
        new_from_id(dist:NumCosmo.Distance, id:NumCosmo.DataBaoId) -> NumCosmo.DataBaoDvrDtDh

    Object NcDataBaoDvrDtDh

    Properties from NcDataBaoDvrDtDh:
      dist -> NcDistance: dist
        Distance object
      z -> NcmVector: z
        Data redshift

    Properties from NcmDataGaussCov:
      n-points -> guint: n-points
        Data sample size
      use-norma -> gboolean: use-norma
        Use the likelihood normalization to calculate -2lnL
      mean -> NcmVector: mean
        Data mean
      cov -> NcmMatrix: cov
        Data covariance

    Properties from NcmData:
      name -> gchararray: name
        Data type name
      desc -> gchararray: desc
        Data description
      long-desc -> gchararray: long-desc
        Data detailed description
      init -> gboolean: init
        Data initialized state
      bootstrap -> NcmBootstrap: bootstrap
        Data bootstrap object

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        dist: Distance
        z: NumCosmoMath.Vector
        cov: NumCosmoMath.Matrix
        mean: NumCosmoMath.Vector
        n_points: int
        use_norma: bool
        bootstrap: NumCosmoMath.Bootstrap
        desc: str
        init: bool
        long_desc: str
        name: str

    props: Props = ...
    def __init__(
        self,
        dist: Distance = ...,
        z: NumCosmoMath.Vector = ...,
        cov: NumCosmoMath.Matrix = ...,
        mean: NumCosmoMath.Vector = ...,
        n_points: int = ...,
        use_norma: bool = ...,
        bootstrap: NumCosmoMath.Bootstrap = ...,
        desc: str = ...,
        init: bool = ...,
        long_desc: str = ...,
    ) -> None: ...
    @classmethod
    def new_from_file(cls, filename: str) -> DataBaoDvrDtDh: ...
    @classmethod
    def new_from_id(cls, dist: Distance, id: DataBaoId) -> DataBaoDvrDtDh: ...
    def set_dist(self, dist: Distance) -> None: ...

class DataBaoDvrDtDhClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DataBaoDvrDtDhClass()
    """

    parent_class: NumCosmoMath.DataGaussCovClass = ...

class DataBaoEmpiricalFit(NumCosmoMath.DataDist1d):
    r"""
    :Constructors:

    ::

        DataBaoEmpiricalFit(**properties)
        new(m2lnp:NumCosmoMath.Spline, Dv_fiduc:float, rs_fiduc:float, z:float) -> NumCosmo.DataBaoEmpiricalFit
        new_from_file(filename:str) -> NumCosmo.DataBaoEmpiricalFit
        new_from_id(dist:NumCosmo.Distance, id:NumCosmo.DataBaoId) -> NumCosmo.DataBaoEmpiricalFit

    Object NcDataBaoEmpiricalFit

    Properties from NcDataBaoEmpiricalFit:
      Dv-fiduc -> gdouble: Dv-fiduc
        Dv fiducial
      rs-fiduc -> gdouble: rs-fiduc
        r_s fiducial
      z -> gdouble: z
        Redshift
      m2lnp -> NcmSpline: m2lnp
        Empirical m2lnp
      dist -> NcDistance: dist
        Distance object

    Properties from NcmDataDist1d:
      n-points -> guint: n-points
        Data sample size
      vector -> NcmVector: vector
        Data vector

    Properties from NcmData:
      name -> gchararray: name
        Data type name
      desc -> gchararray: desc
        Data description
      long-desc -> gchararray: long-desc
        Data detailed description
      init -> gboolean: init
        Data initialized state
      bootstrap -> NcmBootstrap: bootstrap
        Data bootstrap object

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        Dv_fiduc: float
        dist: Distance
        m2lnp: NumCosmoMath.Spline
        rs_fiduc: float
        z: float
        n_points: int
        vector: NumCosmoMath.Vector
        bootstrap: NumCosmoMath.Bootstrap
        desc: str
        init: bool
        long_desc: str
        name: str

    props: Props = ...
    def __init__(
        self,
        Dv_fiduc: float = ...,
        dist: Distance = ...,
        m2lnp: NumCosmoMath.Spline = ...,
        rs_fiduc: float = ...,
        z: float = ...,
        n_points: int = ...,
        vector: NumCosmoMath.Vector = ...,
        bootstrap: NumCosmoMath.Bootstrap = ...,
        desc: str = ...,
        init: bool = ...,
        long_desc: str = ...,
    ) -> None: ...
    def get_alpha(self, mset: NumCosmoMath.MSet) -> float: ...
    def get_mode(self) -> float: ...
    @classmethod
    def new(
        cls, m2lnp: NumCosmoMath.Spline, Dv_fiduc: float, rs_fiduc: float, z: float
    ) -> DataBaoEmpiricalFit: ...
    @classmethod
    def new_from_file(cls, filename: str) -> DataBaoEmpiricalFit: ...
    @classmethod
    def new_from_id(cls, dist: Distance, id: DataBaoId) -> DataBaoEmpiricalFit: ...
    def set_dist(self, dist: Distance) -> None: ...

class DataBaoEmpiricalFit2d(NumCosmoMath.DataDist2d):
    r"""
    :Constructors:

    ::

        DataBaoEmpiricalFit2d(**properties)
        new(m2lnp:NumCosmoMath.Spline2d, Dh_rd_fiduc:float, Dt_rd_fiduc:float, z:float) -> NumCosmo.DataBaoEmpiricalFit2d
        new_from_file(filename:str) -> NumCosmo.DataBaoEmpiricalFit2d
        new_from_id(dist:NumCosmo.Distance, id:NumCosmo.DataBaoId) -> NumCosmo.DataBaoEmpiricalFit2d

    Object NcDataBaoEmpiricalFit2d

    Properties from NcDataBaoEmpiricalFit2d:
      Dh-rd-fiduc -> gdouble: Dh-rd-fiduc
        Dh/rd fiducial
      Dt-rd-fiduc -> gdouble: Dt-rd-fiduc
        Dt/rd fiducial
      z -> gdouble: z
        Redshift
      m2lnp -> NcmSpline2d: m2lnp
        Empirical m2lnp
      dist -> NcDistance: dist
        Distance object

    Properties from NcmDataDist2d:
      n-points -> guint: n-points
        Data sample size
      matrix -> NcmMatrix: matrix
        Data matrix

    Properties from NcmData:
      name -> gchararray: name
        Data type name
      desc -> gchararray: desc
        Data description
      long-desc -> gchararray: long-desc
        Data detailed description
      init -> gboolean: init
        Data initialized state
      bootstrap -> NcmBootstrap: bootstrap
        Data bootstrap object

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        Dh_rd_fiduc: float
        Dt_rd_fiduc: float
        dist: Distance
        m2lnp: NumCosmoMath.Spline2d
        z: float
        matrix: NumCosmoMath.Matrix
        n_points: int
        bootstrap: NumCosmoMath.Bootstrap
        desc: str
        init: bool
        long_desc: str
        name: str

    props: Props = ...
    def __init__(
        self,
        Dh_rd_fiduc: float = ...,
        Dt_rd_fiduc: float = ...,
        dist: Distance = ...,
        m2lnp: NumCosmoMath.Spline2d = ...,
        z: float = ...,
        matrix: NumCosmoMath.Matrix = ...,
        n_points: int = ...,
        bootstrap: NumCosmoMath.Bootstrap = ...,
        desc: str = ...,
        init: bool = ...,
        long_desc: str = ...,
    ) -> None: ...
    def get_alpha_parallel(self, mset: NumCosmoMath.MSet) -> float: ...
    def get_alpha_perpendicular(self, mset: NumCosmoMath.MSet) -> float: ...
    def get_mode(self) -> float: ...
    @classmethod
    def new(
        cls,
        m2lnp: NumCosmoMath.Spline2d,
        Dh_rd_fiduc: float,
        Dt_rd_fiduc: float,
        z: float,
    ) -> DataBaoEmpiricalFit2d: ...
    @classmethod
    def new_from_file(cls, filename: str) -> DataBaoEmpiricalFit2d: ...
    @classmethod
    def new_from_id(cls, dist: Distance, id: DataBaoId) -> DataBaoEmpiricalFit2d: ...
    def set_dist(self, dist: Distance) -> None: ...

class DataBaoEmpiricalFit2dClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DataBaoEmpiricalFit2dClass()
    """

    parent_class: NumCosmoMath.DataDist2dClass = ...

class DataBaoEmpiricalFitClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DataBaoEmpiricalFitClass()
    """

    parent_class: NumCosmoMath.DataDist1dClass = ...

class DataBaoRDV(NumCosmoMath.DataGauss):
    r"""
    :Constructors:

    ::

        DataBaoRDV(**properties)
        new_from_file(filename:str) -> NumCosmo.DataBaoRDV
        new_from_id(dist:NumCosmo.Distance, id:NumCosmo.DataBaoId) -> NumCosmo.DataBaoRDV

    Object NcDataBaoRDV

    Properties from NcDataBaoRDV:
      dist -> NcDistance: dist
        Distance object
      z -> NcmVector: z
        Data redshift
      is-rDV -> gboolean: is-rDV
        Whether the format is r/DV or DV/r

    Properties from NcmDataGauss:
      n-points -> guint: n-points
        Data sample size
      mean -> NcmVector: mean
        Data mean
      inv-cov -> NcmMatrix: inv-cov
        Data covariance inverse

    Properties from NcmData:
      name -> gchararray: name
        Data type name
      desc -> gchararray: desc
        Data description
      long-desc -> gchararray: long-desc
        Data detailed description
      init -> gboolean: init
        Data initialized state
      bootstrap -> NcmBootstrap: bootstrap
        Data bootstrap object

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        dist: Distance
        is_rDV: bool
        z: NumCosmoMath.Vector
        inv_cov: NumCosmoMath.Matrix
        mean: NumCosmoMath.Vector
        n_points: int
        bootstrap: NumCosmoMath.Bootstrap
        desc: str
        init: bool
        long_desc: str
        name: str

    props: Props = ...
    def __init__(
        self,
        dist: Distance = ...,
        is_rDV: bool = ...,
        z: NumCosmoMath.Vector = ...,
        inv_cov: NumCosmoMath.Matrix = ...,
        mean: NumCosmoMath.Vector = ...,
        n_points: int = ...,
        bootstrap: NumCosmoMath.Bootstrap = ...,
        desc: str = ...,
        init: bool = ...,
        long_desc: str = ...,
    ) -> None: ...
    @classmethod
    def new_from_file(cls, filename: str) -> DataBaoRDV: ...
    @classmethod
    def new_from_id(cls, dist: Distance, id: DataBaoId) -> DataBaoRDV: ...
    def set_dist(self, dist: Distance) -> None: ...

class DataBaoRDVClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DataBaoRDVClass()
    """

    parent_class: NumCosmoMath.DataGaussClass = ...

class DataCMBDistPriors(NumCosmoMath.DataGauss):
    r"""
    :Constructors:

    ::

        DataCMBDistPriors(**properties)
        new_empty(dist:NumCosmo.Distance) -> NumCosmo.DataCMBDistPriors
        new_from_file(filename:str) -> NumCosmo.DataCMBDistPriors
        new_from_id(dist:NumCosmo.Distance, id:NumCosmo.DataCMBId) -> NumCosmo.DataCMBDistPriors

    Object NcDataCMBDistPriors

    Properties from NcDataCMBDistPriors:
      dist -> NcDistance: dist
        Distance object

    Properties from NcmDataGauss:
      n-points -> guint: n-points
        Data sample size
      mean -> NcmVector: mean
        Data mean
      inv-cov -> NcmMatrix: inv-cov
        Data covariance inverse

    Properties from NcmData:
      name -> gchararray: name
        Data type name
      desc -> gchararray: desc
        Data description
      long-desc -> gchararray: long-desc
        Data detailed description
      init -> gboolean: init
        Data initialized state
      bootstrap -> NcmBootstrap: bootstrap
        Data bootstrap object

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        dist: Distance
        inv_cov: NumCosmoMath.Matrix
        mean: NumCosmoMath.Vector
        n_points: int
        bootstrap: NumCosmoMath.Bootstrap
        desc: str
        init: bool
        long_desc: str
        name: str

    props: Props = ...
    parent_instance: NumCosmoMath.DataGauss = ...
    dist: Distance = ...
    def __init__(
        self,
        dist: Distance = ...,
        inv_cov: NumCosmoMath.Matrix = ...,
        mean: NumCosmoMath.Vector = ...,
        n_points: int = ...,
        bootstrap: NumCosmoMath.Bootstrap = ...,
        desc: str = ...,
        init: bool = ...,
        long_desc: str = ...,
    ) -> None: ...
    @classmethod
    def new_empty(cls, dist: Distance) -> DataCMBDistPriors: ...
    @classmethod
    def new_from_file(cls, filename: str) -> DataCMBDistPriors: ...
    @classmethod
    def new_from_id(cls, dist: Distance, id: DataCMBId) -> DataCMBDistPriors: ...
    def set_dist(self, dist: Distance) -> None: ...

class DataCMBDistPriorsClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DataCMBDistPriorsClass()
    """

    parent_class: NumCosmoMath.DataGaussClass = ...

class DataCMBShiftParam(NumCosmoMath.DataGaussDiag):
    r"""
    :Constructors:

    ::

        DataCMBShiftParam(**properties)
        new_empty(dist:NumCosmo.Distance) -> NumCosmo.DataCMBShiftParam
        new_from_file(filename:str) -> NumCosmo.DataCMBShiftParam
        new_from_id(dist:NumCosmo.Distance, id:NumCosmo.DataCMBId) -> NumCosmo.DataCMBShiftParam

    Object NcDataCMBShiftParam

    Properties from NcDataCMBShiftParam:
      dist -> NcDistance: dist
        Distance object
      z -> NcmVector: z
        Data redshift

    Properties from NcmDataGaussDiag:
      n-points -> guint: n-points
        Data sample size
      w-mean -> gboolean: w-mean
        Whether to minimize analytically over the weighted mean
      mean -> NcmVector: mean
        Data mean
      sigma -> NcmVector: sigma
        Data standard deviation

    Properties from NcmData:
      name -> gchararray: name
        Data type name
      desc -> gchararray: desc
        Data description
      long-desc -> gchararray: long-desc
        Data detailed description
      init -> gboolean: init
        Data initialized state
      bootstrap -> NcmBootstrap: bootstrap
        Data bootstrap object

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        dist: Distance
        z: NumCosmoMath.Vector
        mean: NumCosmoMath.Vector
        n_points: int
        sigma: NumCosmoMath.Vector
        w_mean: bool
        bootstrap: NumCosmoMath.Bootstrap
        desc: str
        init: bool
        long_desc: str
        name: str

    props: Props = ...
    parent_instance: NumCosmoMath.DataGaussDiag = ...
    dist: Distance = ...
    x: NumCosmoMath.Vector = ...
    def __init__(
        self,
        dist: Distance = ...,
        z: NumCosmoMath.Vector = ...,
        mean: NumCosmoMath.Vector = ...,
        n_points: int = ...,
        sigma: NumCosmoMath.Vector = ...,
        w_mean: bool = ...,
        bootstrap: NumCosmoMath.Bootstrap = ...,
        desc: str = ...,
        init: bool = ...,
        long_desc: str = ...,
    ) -> None: ...
    @classmethod
    def new_empty(cls, dist: Distance) -> DataCMBShiftParam: ...
    @classmethod
    def new_from_file(cls, filename: str) -> DataCMBShiftParam: ...
    @classmethod
    def new_from_id(cls, dist: Distance, id: DataCMBId) -> DataCMBShiftParam: ...
    def set_dist(self, dist: Distance) -> None: ...

class DataCMBShiftParamClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DataCMBShiftParamClass()
    """

    parent_class: NumCosmoMath.DataGaussDiagClass = ...

class DataClusterNCount(NumCosmoMath.Data):
    r"""
    :Constructors:

    ::

        DataClusterNCount(**properties)
        new(cad:NumCosmo.ClusterAbundance, redshift_type:str, mass_type:str) -> NumCosmo.DataClusterNCount

    Object NcDataClusterNCount

    Properties from NcDataClusterNCount:
      cluster-abundance -> NcClusterAbundance: cluster-abundance
        Cluster abundance
      mass-type -> gchararray: mass-type
        Cluster mass proxy type
      redshift-type -> gchararray: redshift-type
        Cluster redshift proxy type
      lnM-true -> NcmVector: lnM-true
        Clusters true masses
      z-true -> NcmVector: z-true
        Clusters true redshifts
      z-obs -> NcmMatrix: z-obs
        Clusters redshift observables
      z-obs-params -> NcmMatrix: z-obs-params
        Clusters redshift observables parameters
      lnM-obs -> NcmMatrix: lnM-obs
        Clusters mass observables
      lnM-obs-params -> NcmMatrix: lnM-obs-params
        Clusters mass observables parameters
      area -> gdouble: area
        Cluster observation area
      use-true -> gboolean: use-true
        If the true data must be used
      binned -> gboolean: binned
        Whether use binned data
      z-obs-bins -> NcmObjArray: z-obs-bins
        Clusters redshifts bins
      lnM-obs-bins -> NcmObjArray: lnM-obs-bins
        Clusters mass bins
      bin-count -> NcmVector: bin-count
        Bin count
      fiducial -> gboolean: fiducial
        If it is fiducial data
      rng-seed -> guint64: rng-seed
        Random number generator seed
      rng-name -> gchararray: rng-name
        Random number generator name

    Properties from NcmData:
      name -> gchararray: name
        Data type name
      desc -> gchararray: desc
        Data description
      long-desc -> gchararray: long-desc
        Data detailed description
      init -> gboolean: init
        Data initialized state
      bootstrap -> NcmBootstrap: bootstrap
        Data bootstrap object

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        area: float
        bin_count: NumCosmoMath.Vector
        binned: bool
        cluster_abundance: ClusterAbundance
        fiducial: bool
        lnM_obs: NumCosmoMath.Matrix
        lnM_obs_bins: NumCosmoMath.ObjArray
        lnM_obs_params: NumCosmoMath.Matrix
        lnM_true: NumCosmoMath.Vector
        mass_type: str
        redshift_type: str
        rng_name: str
        rng_seed: int
        use_true: bool
        z_obs: NumCosmoMath.Matrix
        z_obs_bins: NumCosmoMath.ObjArray
        z_obs_params: NumCosmoMath.Matrix
        z_true: NumCosmoMath.Vector
        bootstrap: NumCosmoMath.Bootstrap
        desc: str
        init: bool
        long_desc: str
        name: str

    props: Props = ...
    parent_instance: NumCosmoMath.Data = ...
    priv: DataClusterNCountPrivate = ...
    def __init__(
        self,
        area: float = ...,
        bin_count: NumCosmoMath.Vector = ...,
        binned: bool = ...,
        cluster_abundance: ClusterAbundance = ...,
        fiducial: bool = ...,
        lnM_obs: NumCosmoMath.Matrix = ...,
        lnM_obs_bins: NumCosmoMath.ObjArray = ...,
        lnM_obs_params: NumCosmoMath.Matrix = ...,
        lnM_true: NumCosmoMath.Vector = ...,
        mass_type: str = ...,
        redshift_type: str = ...,
        rng_name: str = ...,
        rng_seed: int = ...,
        use_true: bool = ...,
        z_obs: NumCosmoMath.Matrix = ...,
        z_obs_bins: NumCosmoMath.ObjArray = ...,
        z_obs_params: NumCosmoMath.Matrix = ...,
        z_true: NumCosmoMath.Vector = ...,
        bootstrap: NumCosmoMath.Bootstrap = ...,
        desc: str = ...,
        init: bool = ...,
        long_desc: str = ...,
    ) -> None: ...
    def add_bin(
        self,
        lnM_obs_lb: NumCosmoMath.Vector,
        lnM_obs_ub: NumCosmoMath.Vector,
        z_obs_lb: NumCosmoMath.Vector,
        z_obs_ub: NumCosmoMath.Vector,
    ) -> None: ...
    def bin_data(self) -> None: ...
    def catalog_load(self, filename: str) -> None: ...
    def catalog_save(self, filename: str, overwrite: bool) -> None: ...
    @staticmethod
    def clear(ncount: DataClusterNCount) -> None: ...
    def del_bins(self) -> None: ...
    def free(self) -> None: ...
    def get_len(self) -> int: ...
    def get_lnM_obs(self) -> NumCosmoMath.Matrix: ...
    def get_lnM_obs_params(self) -> NumCosmoMath.Matrix: ...
    def get_lnM_true(self) -> NumCosmoMath.Vector: ...
    def get_z_obs(self) -> NumCosmoMath.Matrix: ...
    def get_z_obs_params(self) -> NumCosmoMath.Matrix: ...
    def get_z_true(self) -> NumCosmoMath.Vector: ...
    def has_lnM_true(self) -> bool: ...
    def has_z_true(self) -> bool: ...
    def init_from_sampling(
        self, mset: NumCosmoMath.MSet, area_survey: float, rng: NumCosmoMath.RNG
    ) -> None: ...
    def lnM_obs_len(self) -> int: ...
    def lnM_obs_params_len(self) -> int: ...
    @classmethod
    def new(
        cls, cad: ClusterAbundance, redshift_type: str, mass_type: str
    ) -> DataClusterNCount: ...
    def ref(self) -> DataClusterNCount: ...
    def set_bin_count(self, bin_count: NumCosmoMath.Vector) -> None: ...
    def set_binned(self, on: bool) -> None: ...
    def set_lnM_obs(self, m: NumCosmoMath.Matrix) -> None: ...
    def set_lnM_obs_bins(self, lnM_obs_bins: NumCosmoMath.ObjArray) -> None: ...
    def set_lnM_obs_params(self, m: NumCosmoMath.Matrix) -> None: ...
    def set_lnM_true(self, v: NumCosmoMath.Vector) -> None: ...
    def set_z_obs(self, m: NumCosmoMath.Matrix) -> None: ...
    def set_z_obs_bins(self, z_obs_bins: NumCosmoMath.ObjArray) -> None: ...
    def set_z_obs_params(self, m: NumCosmoMath.Matrix) -> None: ...
    def set_z_true(self, v: NumCosmoMath.Vector) -> None: ...
    def true_data(self, use_true_data: bool) -> None: ...
    def using_true_data(self) -> bool: ...
    def z_obs_len(self) -> int: ...
    def z_obs_params_len(self) -> int: ...

class DataClusterNCountClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DataClusterNCountClass()
    """

    parent_class: NumCosmoMath.DataClass = ...

class DataClusterNCountPrivate(GObject.GPointer): ...

class DataClusterNCountsGauss(NumCosmoMath.DataGaussCov):
    r"""
    :Constructors:

    ::

        DataClusterNCountsGauss(**properties)
        new(cad:NumCosmo.ClusterAbundance) -> NumCosmo.DataClusterNCountsGauss

    Object NcDataClusterNCountsGauss

    Properties from NcDataClusterNCountsGauss:
      z-obs -> NcmVector: z-obs
        Clusters redshift observables
      z-obs-params -> NcmMatrix: z-obs-params
        Clusters redshift observables parameters
      lnM-obs -> NcmVector: lnM-obs
        Clusters mass observables
      lnM-obs-params -> NcmMatrix: lnM-obs-params
        Clusters mass observables parameters
      cluster-abundance -> NcClusterAbundance: cluster-abundance
        Cluster abundance
      has-ssc -> gboolean: has-ssc
        Whether use super sample covariance
      s-matrix -> NcmMatrix: s-matrix
        Super sample covariance matrix
      resample-s-matrix -> NcmMatrix: resample-s-matrix
        Super sample covariance resample matrix
      fix-cov -> gboolean: fix-cov
        Whether the covariance matrix is fixed or not

    Properties from NcmDataGaussCov:
      n-points -> guint: n-points
        Data sample size
      use-norma -> gboolean: use-norma
        Use the likelihood normalization to calculate -2lnL
      mean -> NcmVector: mean
        Data mean
      cov -> NcmMatrix: cov
        Data covariance

    Properties from NcmData:
      name -> gchararray: name
        Data type name
      desc -> gchararray: desc
        Data description
      long-desc -> gchararray: long-desc
        Data detailed description
      init -> gboolean: init
        Data initialized state
      bootstrap -> NcmBootstrap: bootstrap
        Data bootstrap object

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        cluster_abundance: ClusterAbundance
        fix_cov: bool
        has_ssc: bool
        lnM_obs: NumCosmoMath.Vector
        lnM_obs_params: NumCosmoMath.Matrix
        resample_s_matrix: NumCosmoMath.Matrix
        s_matrix: NumCosmoMath.Matrix
        z_obs: NumCosmoMath.Vector
        z_obs_params: NumCosmoMath.Matrix
        cov: NumCosmoMath.Matrix
        mean: NumCosmoMath.Vector
        n_points: int
        use_norma: bool
        bootstrap: NumCosmoMath.Bootstrap
        desc: str
        init: bool
        long_desc: str
        name: str

    props: Props = ...
    parent_instance: NumCosmoMath.DataGaussCov = ...
    priv: DataClusterNCountsGaussPrivate = ...
    def __init__(
        self,
        cluster_abundance: ClusterAbundance = ...,
        fix_cov: bool = ...,
        has_ssc: bool = ...,
        lnM_obs: NumCosmoMath.Vector = ...,
        lnM_obs_params: NumCosmoMath.Matrix = ...,
        resample_s_matrix: NumCosmoMath.Matrix = ...,
        s_matrix: NumCosmoMath.Matrix = ...,
        z_obs: NumCosmoMath.Vector = ...,
        z_obs_params: NumCosmoMath.Matrix = ...,
        cov: NumCosmoMath.Matrix = ...,
        mean: NumCosmoMath.Vector = ...,
        n_points: int = ...,
        use_norma: bool = ...,
        bootstrap: NumCosmoMath.Bootstrap = ...,
        desc: str = ...,
        init: bool = ...,
        long_desc: str = ...,
    ) -> None: ...
    def free(self) -> None: ...
    def get_fix_cov(self) -> bool: ...
    def get_has_ssc(self) -> bool: ...
    def get_lnM_obs(self) -> NumCosmoMath.Vector: ...
    def get_lnM_obs_params(self) -> NumCosmoMath.Matrix: ...
    def get_resample_s_matrix(self) -> NumCosmoMath.Matrix: ...
    def get_s_matrix(self) -> NumCosmoMath.Matrix: ...
    def get_z_obs(self) -> NumCosmoMath.Vector: ...
    def get_z_obs_params(self) -> NumCosmoMath.Matrix: ...
    @classmethod
    def new(cls, cad: ClusterAbundance) -> DataClusterNCountsGauss: ...
    def set_fix_cov(self, on: bool) -> None: ...
    def set_has_ssc(self, on: bool) -> None: ...
    def set_lnM_obs(self, lnM_obs: NumCosmoMath.Vector) -> None: ...
    def set_lnM_obs_params(self, lnM_obs_params: NumCosmoMath.Matrix) -> None: ...
    def set_resample_s_matrix(self, s_matrix: NumCosmoMath.Matrix) -> None: ...
    def set_s_matrix(self, s_matrix: NumCosmoMath.Matrix) -> None: ...
    def set_z_obs(self, z_obs: NumCosmoMath.Vector) -> None: ...
    def set_z_obs_params(self, z_obs_params: NumCosmoMath.Matrix) -> None: ...

class DataClusterNCountsGaussClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DataClusterNCountsGaussClass()
    """

    parent_class: NumCosmoMath.DataGaussCovClass = ...

class DataClusterNCountsGaussPrivate(GObject.GPointer): ...

class DataClusterPseudoCounts(NumCosmoMath.Data):
    r"""
    :Constructors:

    ::

        DataClusterPseudoCounts(**properties)
        new(cad:NumCosmo.ClusterAbundance) -> NumCosmo.DataClusterPseudoCounts
        new_from_file(filename:str) -> NumCosmo.DataClusterPseudoCounts

    Object NcDataClusterPseudoCounts

    Properties from NcDataClusterPseudoCounts:
      cluster-abundance -> NcClusterAbundance: cluster-abundance
        Cluster abundance
      np -> guint: np
        Number of clusters
      obs -> NcmMatrix: obs
        Cluster observables
      true-data -> NcmMatrix: true-data
        Cluster (halo) true data (redshift and mass)
      M-z-flat-prior -> gboolean: M-z-flat-prior
        Flat priors for halo mass and selection functions.

    Properties from NcmData:
      name -> gchararray: name
        Data type name
      desc -> gchararray: desc
        Data description
      long-desc -> gchararray: long-desc
        Data detailed description
      init -> gboolean: init
        Data initialized state
      bootstrap -> NcmBootstrap: bootstrap
        Data bootstrap object

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        M_z_flat_prior: bool
        cluster_abundance: ClusterAbundance
        np: int
        obs: NumCosmoMath.Matrix
        true_data: NumCosmoMath.Matrix
        bootstrap: NumCosmoMath.Bootstrap
        desc: str
        init: bool
        long_desc: str
        name: str

    props: Props = ...
    parent_instance: NumCosmoMath.Data = ...
    cad: ClusterAbundance = ...
    obs: NumCosmoMath.Matrix = ...
    true_data: NumCosmoMath.Matrix = ...
    np: int = ...
    M_Z_FlatPrior: bool = ...
    rnd_name: str = ...
    def __init__(
        self,
        M_z_flat_prior: bool = ...,
        cluster_abundance: ClusterAbundance = ...,
        np: int = ...,
        obs: NumCosmoMath.Matrix = ...,
        true_data: NumCosmoMath.Matrix = ...,
        bootstrap: NumCosmoMath.Bootstrap = ...,
        desc: str = ...,
        init: bool = ...,
        long_desc: str = ...,
    ) -> None: ...
    @staticmethod
    def clear(dcpc: DataClusterPseudoCounts) -> None: ...
    def free(self) -> None: ...
    def get_nclusters(self) -> int: ...
    def init_from_sampling(
        self, mset: NumCosmoMath.MSet, rng: NumCosmoMath.RNG, np: int
    ) -> None: ...
    @classmethod
    def new(cls, cad: ClusterAbundance) -> DataClusterPseudoCounts: ...
    @classmethod
    def new_from_file(cls, filename: str) -> DataClusterPseudoCounts: ...
    def ref(self) -> DataClusterPseudoCounts: ...
    def set_cad(self, cad: ClusterAbundance) -> None: ...
    def set_nclusters(self, np: int) -> None: ...
    def set_obs(self, m: NumCosmoMath.Matrix) -> None: ...
    def set_true_data(self, m: NumCosmoMath.Matrix) -> None: ...

class DataClusterPseudoCountsClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DataClusterPseudoCountsClass()
    """

    parent_class: NumCosmoMath.DataClass = ...

class DataClusterWL(NumCosmoMath.Data):
    r"""
    :Constructors:

    ::

        DataClusterWL(**properties)
        new() -> NumCosmo.DataClusterWL

    Object NcDataClusterWL

    Properties from NcDataClusterWL:
      obs -> NcGalaxyWLObs: obs
        Galaxy weak lensing observables
      r-min -> gdouble: r-min
        Minimum radius of the weak lensing observables
      r-max -> gdouble: r-max
        Maximum radius of the weak lensing observables
      prec -> gdouble: prec
        Precision for integral
      len -> guint: len
        Number of galaxies
      resample-flag -> NcDataClusterWLResampleFlag: resample-flag
        Resample flag
      enable-parallel -> gboolean: enable-parallel
        Enable parallelization

    Properties from NcmData:
      name -> gchararray: name
        Data type name
      desc -> gchararray: desc
        Data description
      long-desc -> gchararray: long-desc
        Data detailed description
      init -> gboolean: init
        Data initialized state
      bootstrap -> NcmBootstrap: bootstrap
        Data bootstrap object

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        enable_parallel: bool
        len: int
        obs: GalaxyWLObs
        prec: float
        r_max: float
        r_min: float
        resample_flag: DataClusterWLResampleFlag
        bootstrap: NumCosmoMath.Bootstrap
        desc: str
        init: bool
        long_desc: str
        name: str

    props: Props = ...
    def __init__(
        self,
        enable_parallel: bool = ...,
        len: int = ...,
        obs: GalaxyWLObs = ...,
        prec: float = ...,
        r_max: float = ...,
        r_min: float = ...,
        resample_flag: DataClusterWLResampleFlag = ...,
        bootstrap: NumCosmoMath.Bootstrap = ...,
        desc: str = ...,
        init: bool = ...,
        long_desc: str = ...,
    ) -> None: ...
    @staticmethod
    def clear(dcwl: DataClusterWL) -> None: ...
    def free(self) -> None: ...
    def get_resample_flag(self) -> DataClusterWLResampleFlag: ...
    @classmethod
    def new(cls) -> DataClusterWL: ...
    def peek_data_array(self) -> NumCosmoMath.ObjArray: ...
    def peek_obs(self) -> GalaxyWLObs: ...
    def ref(self) -> DataClusterWL: ...
    def set_cut(self, r_min: float, r_max: float) -> None: ...
    def set_obs(self, obs: GalaxyWLObs) -> None: ...
    def set_prec(self, prec: float) -> None: ...
    def set_resample_flag(self, resample_flag: DataClusterWLResampleFlag) -> None: ...

class DataClusterWLClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DataClusterWLClass()
    """

    parent_class: NumCosmoMath.DataClass = ...

class DataClusterWLPrivate(GObject.GPointer): ...

class DataDistMu(NumCosmoMath.DataGaussDiag):
    r"""
    :Constructors:

    ::

        DataDistMu(**properties)
        new_empty(dist:NumCosmo.Distance) -> NumCosmo.DataDistMu
        new_from_file(filename:str) -> NumCosmo.DataDistMu
        new_from_id(dist:NumCosmo.Distance, id:NumCosmo.DataSNIAId) -> NumCosmo.DataDistMu

    Object NcDataDistMu

    Properties from NcDataDistMu:
      dist -> NcDistance: dist
        Distance object
      z -> NcmVector: z
        Data redshift

    Properties from NcmDataGaussDiag:
      n-points -> guint: n-points
        Data sample size
      w-mean -> gboolean: w-mean
        Whether to minimize analytically over the weighted mean
      mean -> NcmVector: mean
        Data mean
      sigma -> NcmVector: sigma
        Data standard deviation

    Properties from NcmData:
      name -> gchararray: name
        Data type name
      desc -> gchararray: desc
        Data description
      long-desc -> gchararray: long-desc
        Data detailed description
      init -> gboolean: init
        Data initialized state
      bootstrap -> NcmBootstrap: bootstrap
        Data bootstrap object

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        dist: Distance
        z: NumCosmoMath.Vector
        mean: NumCosmoMath.Vector
        n_points: int
        sigma: NumCosmoMath.Vector
        w_mean: bool
        bootstrap: NumCosmoMath.Bootstrap
        desc: str
        init: bool
        long_desc: str
        name: str

    props: Props = ...
    parent_instance: NumCosmoMath.DataGaussDiag = ...
    dist: Distance = ...
    x: NumCosmoMath.Vector = ...
    def __init__(
        self,
        dist: Distance = ...,
        z: NumCosmoMath.Vector = ...,
        mean: NumCosmoMath.Vector = ...,
        n_points: int = ...,
        sigma: NumCosmoMath.Vector = ...,
        w_mean: bool = ...,
        bootstrap: NumCosmoMath.Bootstrap = ...,
        desc: str = ...,
        init: bool = ...,
        long_desc: str = ...,
    ) -> None: ...
    @classmethod
    def new_empty(cls, dist: Distance) -> DataDistMu: ...
    @classmethod
    def new_from_file(cls, filename: str) -> DataDistMu: ...
    @classmethod
    def new_from_id(cls, dist: Distance, id: DataSNIAId) -> DataDistMu: ...
    def set_dist(self, dist: Distance) -> None: ...

class DataDistMuClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DataDistMuClass()
    """

    parent_class: NumCosmoMath.DataGaussDiagClass = ...

class DataHubble(NumCosmoMath.DataGaussDiag):
    r"""
    :Constructors:

    ::

        DataHubble(**properties)
        new_empty() -> NumCosmo.DataHubble
        new_from_file(filename:str) -> NumCosmo.DataHubble
        new_from_id(id:NumCosmo.DataHubbleId) -> NumCosmo.DataHubble

    Object NcDataHubble

    Properties from NcDataHubble:
      z -> NcmVector: z
        Data redshifts

    Properties from NcmDataGaussDiag:
      n-points -> guint: n-points
        Data sample size
      w-mean -> gboolean: w-mean
        Whether to minimize analytically over the weighted mean
      mean -> NcmVector: mean
        Data mean
      sigma -> NcmVector: sigma
        Data standard deviation

    Properties from NcmData:
      name -> gchararray: name
        Data type name
      desc -> gchararray: desc
        Data description
      long-desc -> gchararray: long-desc
        Data detailed description
      init -> gboolean: init
        Data initialized state
      bootstrap -> NcmBootstrap: bootstrap
        Data bootstrap object

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        z: NumCosmoMath.Vector
        mean: NumCosmoMath.Vector
        n_points: int
        sigma: NumCosmoMath.Vector
        w_mean: bool
        bootstrap: NumCosmoMath.Bootstrap
        desc: str
        init: bool
        long_desc: str
        name: str

    props: Props = ...
    def __init__(
        self,
        z: NumCosmoMath.Vector = ...,
        mean: NumCosmoMath.Vector = ...,
        n_points: int = ...,
        sigma: NumCosmoMath.Vector = ...,
        w_mean: bool = ...,
        bootstrap: NumCosmoMath.Bootstrap = ...,
        desc: str = ...,
        init: bool = ...,
        long_desc: str = ...,
    ) -> None: ...
    @classmethod
    def new_empty(cls) -> DataHubble: ...
    @classmethod
    def new_from_file(cls, filename: str) -> DataHubble: ...
    @classmethod
    def new_from_id(cls, id: DataHubbleId) -> DataHubble: ...
    def set_sample(self, id: DataHubbleId) -> None: ...

class DataHubbleBao(NumCosmoMath.DataGaussDiag):
    r"""
    :Constructors:

    ::

        DataHubbleBao(**properties)
        new(dist:NumCosmo.Distance, id:NumCosmo.DataHubbleBaoId) -> NumCosmoMath.Data

    Object NcDataHubbleBao

    Properties from NcDataHubbleBao:
      dist -> NcDistance: dist
        Distance object
      z -> NcmVector: z
        Data redshifts

    Properties from NcmDataGaussDiag:
      n-points -> guint: n-points
        Data sample size
      w-mean -> gboolean: w-mean
        Whether to minimize analytically over the weighted mean
      mean -> NcmVector: mean
        Data mean
      sigma -> NcmVector: sigma
        Data standard deviation

    Properties from NcmData:
      name -> gchararray: name
        Data type name
      desc -> gchararray: desc
        Data description
      long-desc -> gchararray: long-desc
        Data detailed description
      init -> gboolean: init
        Data initialized state
      bootstrap -> NcmBootstrap: bootstrap
        Data bootstrap object

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        dist: Distance
        z: NumCosmoMath.Vector
        mean: NumCosmoMath.Vector
        n_points: int
        sigma: NumCosmoMath.Vector
        w_mean: bool
        bootstrap: NumCosmoMath.Bootstrap
        desc: str
        init: bool
        long_desc: str
        name: str

    props: Props = ...
    parent_instance: NumCosmoMath.DataGaussDiag = ...
    dist: Distance = ...
    x: NumCosmoMath.Vector = ...
    def __init__(
        self,
        dist: Distance = ...,
        z: NumCosmoMath.Vector = ...,
        mean: NumCosmoMath.Vector = ...,
        n_points: int = ...,
        sigma: NumCosmoMath.Vector = ...,
        w_mean: bool = ...,
        bootstrap: NumCosmoMath.Bootstrap = ...,
        desc: str = ...,
        init: bool = ...,
        long_desc: str = ...,
    ) -> None: ...
    @classmethod
    def new(cls, dist: Distance, id: DataHubbleBaoId) -> DataHubbleBao: ...
    def set_sample(self, id: DataHubbleBaoId) -> None: ...

class DataHubbleBaoClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DataHubbleBaoClass()
    """

    parent_class: NumCosmoMath.DataGaussDiagClass = ...

class DataHubbleClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DataHubbleClass()
    """

    parent_class: NumCosmoMath.DataGaussDiagClass = ...

class DataPlanckLKL(NumCosmoMath.Data):
    r"""
    :Constructors:

    ::

        DataPlanckLKL(**properties)
        full_new(filename:str, pb:NumCosmo.HIPertBoltzmann) -> NumCosmo.DataPlanckLKL
        full_new_id(id:NumCosmo.DataPlanckLKLType, pb:NumCosmo.HIPertBoltzmann) -> NumCosmo.DataPlanckLKL
        new(filename:str) -> NumCosmo.DataPlanckLKL

    Object NcDataPlanckLKL

    Properties from NcDataPlanckLKL:
      data-file -> gchararray: data-file
        Data file
      hipert-boltzmann -> NcHIPertBoltzmann: hipert-boltzmann
        NcHIPertBoltzmann object
      is-lensing -> gboolean: is-lensing
        Whether the likelihood has lensing
      nparams -> guint: nparams
        Number of expected params
      checksum -> gchararray: checksum
        Params names checksum

    Properties from NcmData:
      name -> gchararray: name
        Data type name
      desc -> gchararray: desc
        Data description
      long-desc -> gchararray: long-desc
        Data detailed description
      init -> gboolean: init
        Data initialized state
      bootstrap -> NcmBootstrap: bootstrap
        Data bootstrap object

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        checksum: str
        data_file: str
        hipert_boltzmann: HIPertBoltzmann
        is_lensing: bool
        nparams: int
        bootstrap: NumCosmoMath.Bootstrap
        desc: str
        init: bool
        long_desc: str
        name: str

    props: Props = ...
    def __init__(
        self,
        data_file: str = ...,
        hipert_boltzmann: HIPertBoltzmann = ...,
        bootstrap: NumCosmoMath.Bootstrap = ...,
        desc: str = ...,
        init: bool = ...,
        long_desc: str = ...,
    ) -> None: ...
    @staticmethod
    def download_baseline(dir: str) -> None: ...
    @classmethod
    def full_new(cls, filename: str, pb: HIPertBoltzmann) -> DataPlanckLKL: ...
    @classmethod
    def full_new_id(
        cls, id: DataPlanckLKLType, pb: HIPertBoltzmann
    ) -> DataPlanckLKL: ...
    def get_param_name(self, i: int) -> str: ...
    def get_param_names(self) -> list[str]: ...
    @classmethod
    def new(cls, filename: str) -> DataPlanckLKL: ...
    def set_hipert_boltzmann(self, pb: HIPertBoltzmann) -> None: ...

class DataPlanckLKLClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DataPlanckLKLClass()
    """

    parent_class: NumCosmoMath.DataClass = ...

class DataSNIACov(NumCosmoMath.DataGaussCov):
    r"""
    :Constructors:

    ::

        DataSNIACov(**properties)
        new(use_norma:bool, cat_version:int) -> NumCosmo.DataSNIACov
        new_from_cat_id(id:NumCosmo.DataSNIAId, use_norma:bool) -> NumCosmo.DataSNIACov
        new_full(filename:str, use_norma:bool) -> NumCosmo.DataSNIACov

    Object NcDataSNIACov

    Properties from NcDataSNIACov:
      cat-version -> guint: cat-version
        Catalog version
      catalog-file -> gchararray: catalog-file
        Catalog file
      magnitude-cut -> gdouble: magnitude-cut
        Threshold where to change absolute magnitude
      z-hd -> NcmVector: z-hd
        Data CMB redshifts (peculiar velocity corrected)
      z-cmb -> NcmVector: z-cmb
        Data cmb redshifts
      z-He -> NcmVector: z-He
        Data He redshifts
      sigma-z -> NcmVector: sigma-z
        Redshifts standard deviation
      magnitudes -> NcmVector: magnitudes
        Magnitudes
      magnitude-b-corrected -> NcmVector: magnitude-b-corrected
        Magnitude B corrected
      ceph-dist -> NcmVector: ceph-dist
        Cepheid distance
      width -> NcmVector: width
        Width
      colour -> NcmVector: colour
        Colour
      thirdpar -> NcmVector: thirdpar
        Thirdpar
      sigma-thirdpar -> NcmVector: sigma-thirdpar
        Thirdpar standard deviation
      absmag-set -> GVariant: absmag-set
        Absolute magnitude set
      is-calib -> GVariant: is-calib
        Whether the SNIa is a calibrator
      used-in-sh0es -> GVariant: used-in-sh0es
        Whether the SNIa was used in SH0ES
      cov-full -> NcmMatrix: cov-full
        Full covariance matrix
      has-complete-cov -> gboolean: has-complete-cov
        Whether the covariance matrix is complete
      cov-mbc-mbc -> NcmMatrix: cov-mbc-mbc
        Covariance matrix for mag b corr

    Properties from NcmDataGaussCov:
      n-points -> guint: n-points
        Data sample size
      use-norma -> gboolean: use-norma
        Use the likelihood normalization to calculate -2lnL
      mean -> NcmVector: mean
        Data mean
      cov -> NcmMatrix: cov
        Data covariance

    Properties from NcmData:
      name -> gchararray: name
        Data type name
      desc -> gchararray: desc
        Data description
      long-desc -> gchararray: long-desc
        Data detailed description
      init -> gboolean: init
        Data initialized state
      bootstrap -> NcmBootstrap: bootstrap
        Data bootstrap object

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        absmag_set: GLib.Variant
        cat_version: int
        catalog_file: str
        ceph_dist: NumCosmoMath.Vector
        colour: NumCosmoMath.Vector
        cov_full: NumCosmoMath.Matrix
        cov_mbc_mbc: NumCosmoMath.Matrix
        has_complete_cov: bool
        is_calib: GLib.Variant
        magnitude_b_corrected: NumCosmoMath.Vector
        magnitude_cut: float
        magnitudes: NumCosmoMath.Vector
        sigma_thirdpar: NumCosmoMath.Vector
        sigma_z: NumCosmoMath.Vector
        thirdpar: NumCosmoMath.Vector
        used_in_sh0es: GLib.Variant
        width: NumCosmoMath.Vector
        z_He: NumCosmoMath.Vector
        z_cmb: NumCosmoMath.Vector
        z_hd: NumCosmoMath.Vector
        cov: NumCosmoMath.Matrix
        mean: NumCosmoMath.Vector
        n_points: int
        use_norma: bool
        bootstrap: NumCosmoMath.Bootstrap
        desc: str
        init: bool
        long_desc: str
        name: str

    props: Props = ...
    def __init__(
        self,
        absmag_set: GLib.Variant = ...,
        cat_version: int = ...,
        catalog_file: str = ...,
        ceph_dist: NumCosmoMath.Vector = ...,
        colour: NumCosmoMath.Vector = ...,
        cov_full: NumCosmoMath.Matrix = ...,
        cov_mbc_mbc: NumCosmoMath.Matrix = ...,
        has_complete_cov: bool = ...,
        is_calib: GLib.Variant = ...,
        magnitude_b_corrected: NumCosmoMath.Vector = ...,
        magnitude_cut: float = ...,
        magnitudes: NumCosmoMath.Vector = ...,
        sigma_thirdpar: NumCosmoMath.Vector = ...,
        sigma_z: NumCosmoMath.Vector = ...,
        thirdpar: NumCosmoMath.Vector = ...,
        used_in_sh0es: GLib.Variant = ...,
        width: NumCosmoMath.Vector = ...,
        z_He: NumCosmoMath.Vector = ...,
        z_cmb: NumCosmoMath.Vector = ...,
        z_hd: NumCosmoMath.Vector = ...,
        cov: NumCosmoMath.Matrix = ...,
        mean: NumCosmoMath.Vector = ...,
        n_points: int = ...,
        use_norma: bool = ...,
        bootstrap: NumCosmoMath.Bootstrap = ...,
        desc: str = ...,
        init: bool = ...,
        long_desc: str = ...,
    ) -> None: ...
    def apply_filter_sh0es_z(self, z_min: float, use_calib: bool) -> DataSNIACov: ...
    def estimate_width_colour(self, mset: NumCosmoMath.MSet) -> float: ...
    @staticmethod
    def get_catalog(id: str) -> str: ...
    @staticmethod
    def get_catalog_by_id(id: DataSNIAId) -> str: ...
    @staticmethod
    def get_catalog_id(id: str) -> DataSNIAId: ...
    def get_estimated_colour(self, mset: NumCosmoMath.MSet) -> NumCosmoMath.Vector: ...
    def get_estimated_mag(self, mset: NumCosmoMath.MSet) -> NumCosmoMath.Vector: ...
    def get_estimated_width(self, mset: NumCosmoMath.MSet) -> NumCosmoMath.Vector: ...
    @staticmethod
    def get_fits(filename: str, check_size: bool) -> str: ...
    def get_mag_cut(self) -> float: ...
    def load_txt(self, filename: str) -> None: ...
    @classmethod
    def new(cls, use_norma: bool, cat_version: int) -> DataSNIACov: ...
    @classmethod
    def new_from_cat_id(cls, id: DataSNIAId, use_norma: bool) -> DataSNIACov: ...
    @classmethod
    def new_full(cls, filename: str, use_norma: bool) -> DataSNIACov: ...
    def peek_abs_mag_set(self) -> list[int]: ...
    def peek_ceph_dist(self) -> NumCosmoMath.Vector: ...
    def peek_colour(self) -> NumCosmoMath.Vector: ...
    def peek_cov_full(self) -> NumCosmoMath.Matrix: ...
    def peek_cov_mbc_mbc(self) -> NumCosmoMath.Matrix: ...
    def peek_cov_packed(self) -> NumCosmoMath.Vector: ...
    def peek_dataset(self) -> list[int]: ...
    def peek_is_calib(self) -> list[int]: ...
    def peek_mag(self) -> NumCosmoMath.Vector: ...
    def peek_sigma_z(self) -> NumCosmoMath.Vector: ...
    def peek_thirdpar(self) -> NumCosmoMath.Vector: ...
    def peek_used_in_sh0es(self) -> list[int]: ...
    def peek_width(self) -> NumCosmoMath.Vector: ...
    def peek_z_cmb(self) -> NumCosmoMath.Vector: ...
    def peek_z_hd(self) -> NumCosmoMath.Vector: ...
    def peek_z_he(self) -> NumCosmoMath.Vector: ...
    def save(self, filename: str, overwrite: bool) -> None: ...
    def set_abs_mag_set(self, abs_mag_set: typing.Sequence[int]) -> None: ...
    def set_ceph_dist(self, ceph_dist: NumCosmoMath.Vector) -> None: ...
    def set_colour(self, colour: NumCosmoMath.Vector) -> None: ...
    def set_cov_full(self, cov_full: NumCosmoMath.Matrix) -> None: ...
    def set_cov_mbc_mbc(self, cov_mbc_mbc: NumCosmoMath.Matrix) -> None: ...
    def set_is_calib(self, is_calib: typing.Sequence[int]) -> None: ...
    def set_mag(self, mag: NumCosmoMath.Vector) -> None: ...
    def set_mag_b_corr(self, mag_b_corr: NumCosmoMath.Vector) -> None: ...
    def set_mag_cut(self, mag_cut: float) -> None: ...
    def set_sigma_thirdpar(self, sigma_thirdpar: NumCosmoMath.Vector) -> None: ...
    def set_sigma_z(self, sigma_z: NumCosmoMath.Vector) -> None: ...
    def set_thirdpar(self, thirdpar: NumCosmoMath.Vector) -> None: ...
    def set_used_in_sh0es(self, used_in_sh0es: typing.Sequence[int]) -> None: ...
    def set_width(self, width: NumCosmoMath.Vector) -> None: ...
    def set_z_cmb(self, z_cmb: NumCosmoMath.Vector) -> None: ...
    def set_z_hd(self, z_hd: NumCosmoMath.Vector) -> None: ...
    def set_z_he(self, z_he: NumCosmoMath.Vector) -> None: ...
    def sigma_int_len(self) -> int: ...
    def snia_len(self) -> int: ...

class DataSNIACovClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DataSNIACovClass()
    """

    parent_class: NumCosmoMath.DataGaussCovClass = ...

class DataXcor(NumCosmoMath.DataGaussCov):
    r"""
    :Constructors:

    ::

        DataXcor(**properties)
        new_full(nobs:int, xc:NumCosmo.Xcor, use_norma:bool) -> NumCosmo.DataXcor

    Object NcDataXcor

    Properties from NcDataXcor:
      nobs -> guint: nobs
        Number of observables
      xcab-oa -> NcmObjArray: xcab-oa
        NcXcorAB array
      X1 -> NcmMatrix: X1
        X matrix
      X2 -> NcmMatrix: X2
        X matrix
      xc -> NcXcor: xc
        Xcor object to compute theoretical spectra

    Properties from NcmDataGaussCov:
      n-points -> guint: n-points
        Data sample size
      use-norma -> gboolean: use-norma
        Use the likelihood normalization to calculate -2lnL
      mean -> NcmVector: mean
        Data mean
      cov -> NcmMatrix: cov
        Data covariance

    Properties from NcmData:
      name -> gchararray: name
        Data type name
      desc -> gchararray: desc
        Data description
      long-desc -> gchararray: long-desc
        Data detailed description
      init -> gboolean: init
        Data initialized state
      bootstrap -> NcmBootstrap: bootstrap
        Data bootstrap object

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        X1: NumCosmoMath.Matrix
        X2: NumCosmoMath.Matrix
        nobs: int
        xc: Xcor
        xcab_oa: NumCosmoMath.ObjArray
        cov: NumCosmoMath.Matrix
        mean: NumCosmoMath.Vector
        n_points: int
        use_norma: bool
        bootstrap: NumCosmoMath.Bootstrap
        desc: str
        init: bool
        long_desc: str
        name: str

    props: Props = ...
    def __init__(
        self,
        X1: NumCosmoMath.Matrix = ...,
        X2: NumCosmoMath.Matrix = ...,
        nobs: int = ...,
        xc: Xcor = ...,
        xcab_oa: NumCosmoMath.ObjArray = ...,
        cov: NumCosmoMath.Matrix = ...,
        mean: NumCosmoMath.Vector = ...,
        n_points: int = ...,
        use_norma: bool = ...,
        bootstrap: NumCosmoMath.Bootstrap = ...,
        desc: str = ...,
        init: bool = ...,
        long_desc: str = ...,
    ) -> None: ...
    def cov_func_abcd(
        self, cov: NumCosmoMath.Matrix, a: int, b: int, c: int, d: int
    ) -> None: ...
    def get_cl_obs(self, vp: NumCosmoMath.Vector, a: int, b: int) -> None: ...
    def mean_func_ab(self, vp: NumCosmoMath.Vector, a: int, b: int) -> None: ...
    @classmethod
    def new_full(cls, nobs: int, xc: Xcor, use_norma: bool) -> DataXcor: ...
    def set_2(
        self,
        a: int,
        b: int,
        ell_th_cut_off: int,
        ell_lik_min: int,
        ell_lik_max: int,
        clobs_filename: str,
        mixing_filename: str,
        mixing_filelength: int,
    ) -> None: ...
    def set_3(self) -> None: ...
    def set_4(
        self,
        a: int,
        b: int,
        c: int,
        d: int,
        X1_filename: str,
        X2_filename: str,
        X_filelength: int,
    ) -> None: ...
    def set_5(self) -> None: ...
    def set_AB(self, xcab: XcorAB) -> None: ...

class DataXcorClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DataXcorClass()
    """

    parent_class: NumCosmoMath.DataGaussCovClass = ...

class Distance(GObject.Object):
    r"""
    :Constructors:

    ::

        Distance(**properties)
        new(zf:float) -> NumCosmo.Distance

    Object NcDistance

    Properties from NcDistance:
      zf -> gdouble: zf
        Final cached redshift
      recomb -> NcRecomb: recomb
        Recombination object
      compute-inv-comoving -> gboolean: compute-inv-comoving
        Whether to compute the inverse comoving function

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        compute_inv_comoving: bool
        recomb: Recomb
        zf: float

    props: Props = ...
    parent_instance: GObject.Object = ...
    comoving_distance_spline: NumCosmoMath.OdeSpline = ...
    inv_comoving_dist: NumCosmoMath.Spline = ...
    comoving_distance_cache: NumCosmoMath.FunctionCache = ...
    comoving_infinity: NumCosmoMath.FunctionCache = ...
    time_cache: NumCosmoMath.FunctionCache = ...
    lookback_time_cache: NumCosmoMath.FunctionCache = ...
    conformal_time_cache: NumCosmoMath.FunctionCache = ...
    sound_horizon_cache: NumCosmoMath.FunctionCache = ...
    ctrl: NumCosmoMath.ModelCtrl = ...
    zf: float = ...
    use_cache: bool = ...
    cpu_inv_comoving: bool = ...
    max_comoving: float = ...
    recomb: Recomb = ...
    cmethod: DistanceComovingMethod = ...
    def __init__(
        self, compute_inv_comoving: bool = ..., recomb: Recomb = ..., zf: float = ...
    ) -> None: ...
    def DA_r(self, cosmo: HICosmo, z: float) -> float: ...
    def DH_r(self, cosmo: HICosmo, z: float) -> float: ...
    def Dt_r(self, cosmo: HICosmo, z: float) -> float: ...
    def acoustic_scale(self, cosmo: HICosmo) -> float: ...
    def angular_diameter(self, cosmo: HICosmo, z: float) -> float: ...
    def angular_diameter_array(
        self, cosmo: HICosmo, z: typing.Sequence[float]
    ) -> list[float]: ...
    def angular_diameter_curvature_scale(self, cosmo: HICosmo) -> float: ...
    def angular_diameter_vector(
        self, cosmo: HICosmo, z: NumCosmoMath.Vector, DA: NumCosmoMath.Vector
    ) -> None: ...
    def angular_diameter_z1_z2(self, cosmo: HICosmo, z1: float, z2: float) -> float: ...
    def bao_A_scale(self, cosmo: HICosmo, z: float) -> float: ...
    def bao_r_Dv(self, cosmo: HICosmo, z: float) -> float: ...
    @staticmethod
    def clear(dist: Distance) -> None: ...
    def comoving(self, cosmo: HICosmo, z: float) -> float: ...
    def comoving_array(
        self, cosmo: HICosmo, z: typing.Sequence[float]
    ) -> list[float]: ...
    def comoving_lss(self, cosmo: HICosmo) -> float: ...
    def comoving_vector(
        self, cosmo: HICosmo, z: NumCosmoMath.Vector, Dc: NumCosmoMath.Vector
    ) -> None: ...
    def comoving_volume_element(self, cosmo: HICosmo, z: float) -> float: ...
    def comoving_z1_z2(self, cosmo: HICosmo, z1: float, z2: float) -> float: ...
    def comoving_z_to_infinity(self, cosmo: HICosmo, z: float) -> float: ...
    def compute_inv_comoving(self, cpu_inv_xi: bool) -> None: ...
    def conformal_lookback_time(self, cosmo: HICosmo, z: float) -> float: ...
    def conformal_time(self, cosmo: HICosmo, z: float) -> float: ...
    def cosmic_time(self, cosmo: HICosmo, z: float) -> float: ...
    def decoupling_redshift(self, cosmo: HICosmo) -> float: ...
    def dilation_scale(self, cosmo: HICosmo, z: float) -> float: ...
    def dmodulus(self, cosmo: HICosmo, z: float) -> float: ...
    def dmodulus_array(
        self, cosmo: HICosmo, z: typing.Sequence[float]
    ) -> list[float]: ...
    def dmodulus_hef(self, cosmo: HICosmo, z_he: float, z_cmb: float) -> float: ...
    def dmodulus_vector(
        self, cosmo: HICosmo, z: NumCosmoMath.Vector, dmod: NumCosmoMath.Vector
    ) -> None: ...
    def drag_redshift(self, cosmo: HICosmo) -> float: ...
    def dsound_horizon_dz(self, cosmo: HICosmo, z: float) -> float: ...
    def dtransverse_dz(self, cosmo: HICosmo, z: float) -> float: ...
    def free(self) -> None: ...
    def hubble(self, cosmo: HICosmo) -> float: ...
    def inv_comoving(self, cosmo: HICosmo, xi: float) -> float: ...
    def lookback_time(self, cosmo: HICosmo, z: float) -> float: ...
    def luminosity(self, cosmo: HICosmo, z: float) -> float: ...
    def luminosity_array(
        self, cosmo: HICosmo, z: typing.Sequence[float]
    ) -> list[float]: ...
    def luminosity_hef(self, cosmo: HICosmo, z_he: float, z_cmb: float) -> float: ...
    def luminosity_vector(
        self, cosmo: HICosmo, z: NumCosmoMath.Vector, Dl: NumCosmoMath.Vector
    ) -> None: ...
    @classmethod
    def new(cls, zf: float) -> Distance: ...
    def prepare(self, cosmo: HICosmo) -> None: ...
    def prepare_if_needed(self, cosmo: HICosmo) -> None: ...
    def r_zd(self, cosmo: HICosmo) -> float: ...
    def r_zd_Mpc(self, cosmo: HICosmo) -> float: ...
    def ref(self) -> Distance: ...
    def require_zf(self, zf: float) -> None: ...
    def set_recomb(self, recomb: Recomb) -> None: ...
    def shift_parameter(self, cosmo: HICosmo, z: float) -> float: ...
    def shift_parameter_lss(self, cosmo: HICosmo) -> float: ...
    def sigma_critical(self, cosmo: HICosmo, zs: float, zl: float) -> float: ...
    def sigma_critical_infinity(self, cosmo: HICosmo, zl: float) -> float: ...
    def sound_horizon(self, cosmo: HICosmo, z: float) -> float: ...
    def theta100CMB(self, cosmo: HICosmo) -> float: ...
    def transverse(self, cosmo: HICosmo, z: float) -> float: ...
    def transverse_array(
        self, cosmo: HICosmo, z: typing.Sequence[float]
    ) -> list[float]: ...
    def transverse_vector(
        self, cosmo: HICosmo, z: NumCosmoMath.Vector, Dt: NumCosmoMath.Vector
    ) -> None: ...
    def transverse_z1_z2(self, cosmo: HICosmo, z1: float, z2: float) -> float: ...
    def transverse_z_to_infinity(self, cosmo: HICosmo, z: float) -> float: ...

class DistanceClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DistanceClass()
    """

    parent_class: GObject.ObjectClass = ...

class DistanceFunc(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DistanceFunc()
    """

    name: str = ...
    desc: str = ...
    f: typing.Callable[[Distance, HICosmo], float] = ...
    impl: HICosmoImpl = ...

class DistanceFuncZ(GObject.GPointer):
    r"""
    :Constructors:

    ::

        DistanceFuncZ()
    """

    name: str = ...
    desc: str = ...
    f: typing.Callable[[Distance, HICosmo, float], float] = ...
    impl: HICosmoImpl = ...

class GalaxyAcf(GObject.Object):
    r"""
    :Constructors:

    ::

        GalaxyAcf(**properties)
        new(gf:NumCosmo.GrowthFunc, dist:NumCosmo.Distance, tf:NumCosmo.TransferFunc) -> NumCosmo.GalaxyAcf

    Object NcGalaxyAcf

    Signals from GObject:
      notify (GParam)
    """

    parent_instance: GObject.Object = ...
    gf: GrowthFunc = ...
    dist: Distance = ...
    tf: TransferFunc = ...
    s: NumCosmoMath.Spline = ...
    b: float = ...
    @classmethod
    def new(cls, gf: GrowthFunc, dist: Distance, tf: TransferFunc) -> GalaxyAcf: ...
    def prepare_psi(self, cosmo: HICosmo, l: int) -> None: ...
    def psi(self, cosmo: HICosmo, k: float, l: int) -> float: ...

class GalaxyAcfClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        GalaxyAcfClass()
    """

    parent_class: GObject.ObjectClass = ...

class GalaxySDObsRedshift(NumCosmoMath.Model):
    r"""
    :Constructors:

    ::

        GalaxySDObsRedshift(**properties)

    Object NcGalaxySDObsRedshift

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: NumCosmoMath.Model = ...
    def __init__(
        self,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def clear(gsdor: GalaxySDObsRedshift) -> None: ...
    def do_data_init(self, data: GalaxySDObsRedshiftData) -> None: ...
    def do_gen(self, data: GalaxySDObsRedshiftData, rng: NumCosmoMath.RNG) -> None: ...
    def do_gen1(self, data: GalaxySDObsRedshiftData, rng: NumCosmoMath.RNG) -> bool: ...
    def do_get_lim(
        self, data: GalaxySDObsRedshiftData
    ) -> typing.Tuple[float, float]: ...
    def do_integ(self) -> GalaxySDObsRedshiftIntegrand: ...
    def do_prepare(self, data: GalaxySDObsRedshiftData) -> None: ...
    def free(self) -> None: ...
    def gen(self, data: GalaxySDObsRedshiftData, rng: NumCosmoMath.RNG) -> None: ...
    def gen1(self, data: GalaxySDObsRedshiftData, rng: NumCosmoMath.RNG) -> bool: ...
    def get_lim(self, data: GalaxySDObsRedshiftData) -> typing.Tuple[float, float]: ...
    @staticmethod
    def id() -> int: ...
    def integ(self) -> GalaxySDObsRedshiftIntegrand: ...
    def prepare(self, data: GalaxySDObsRedshiftData) -> None: ...
    def ref(self) -> GalaxySDObsRedshift: ...

class GalaxySDObsRedshiftClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        GalaxySDObsRedshiftClass()
    """

    parent_class: NumCosmoMath.ModelClass = ...
    gen: typing.Callable[
        [GalaxySDObsRedshift, GalaxySDObsRedshiftData, NumCosmoMath.RNG], None
    ] = ...
    gen1: typing.Callable[
        [GalaxySDObsRedshift, GalaxySDObsRedshiftData, NumCosmoMath.RNG], bool
    ] = ...
    integ: typing.Callable[[GalaxySDObsRedshift], GalaxySDObsRedshiftIntegrand] = ...
    prepare: typing.Callable[[GalaxySDObsRedshift, GalaxySDObsRedshiftData], None] = ...
    get_lim: typing.Callable[
        [GalaxySDObsRedshift, GalaxySDObsRedshiftData], typing.Tuple[float, float]
    ] = ...
    data_init: typing.Callable[[GalaxySDObsRedshift, GalaxySDObsRedshiftData], None] = (
        ...
    )
    padding: list[None] = ...

class GalaxySDObsRedshiftData(GObject.GBoxed):
    r"""
    :Constructors:

    ::

        GalaxySDObsRedshiftData()
        new(gsdor:NumCosmo.GalaxySDObsRedshift) -> NumCosmo.GalaxySDObsRedshiftData
    """

    z: float = ...
    ldata: None = ...
    ldata_destroy: typing.Callable[[None], None] = ...
    ldata_read_row: typing.Callable[
        [GalaxySDObsRedshiftData, GalaxyWLObs, int], None
    ] = ...
    ldata_write_row: typing.Callable[
        [GalaxySDObsRedshiftData, GalaxyWLObs, int], None
    ] = ...
    ldata_required_columns: None = ...
    ref_count: int = ...
    @classmethod
    def new(cls, gsdor: GalaxySDObsRedshift) -> GalaxySDObsRedshiftData: ...
    def read_row(self, obs: GalaxyWLObs, i: int) -> None: ...
    def ref(self) -> GalaxySDObsRedshiftData: ...
    def required_columns(self) -> list[str]: ...
    def unref(self) -> None: ...
    def write_row(self, obs: GalaxyWLObs, i: int) -> None: ...

class GalaxySDObsRedshiftGauss(GalaxySDObsRedshift):
    r"""
    :Constructors:

    ::

        GalaxySDObsRedshiftGauss(**properties)
        new(sdz:NumCosmo.GalaxySDTrueRedshift, zp_min:float, zp_max:float) -> NumCosmo.GalaxySDObsRedshiftGauss

    Object NcGalaxySDObsRedshiftGauss

    Properties from NcGalaxySDObsRedshiftGauss:
      lim -> NcmDTuple2: lim
        Galaxy sample photometric redshift distribution limits
      use-true-z -> gboolean: use-true-z
        Use the true redshift distribution

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        lim: NumCosmoMath.DTuple2
        use_true_z: bool
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    def __init__(
        self,
        lim: NumCosmoMath.DTuple2 = ...,
        use_true_z: bool = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def clear(gsdorgauss: GalaxySDObsRedshiftGauss) -> None: ...
    def data_get(
        self, data: GalaxySDObsRedshiftData, zp: float, sigma0: float, sigma_z: float
    ) -> None: ...
    def data_set(
        self, data: GalaxySDObsRedshiftData, zp: float, sigma0: float, sigma_z: float
    ) -> None: ...
    def free(self) -> None: ...
    def gen(
        self,
        mset: NumCosmoMath.MSet,
        data: GalaxySDObsRedshiftData,
        sigma0: float,
        rng: NumCosmoMath.RNG,
    ) -> None: ...
    def gen1(
        self,
        mset: NumCosmoMath.MSet,
        data: GalaxySDObsRedshiftData,
        sigma0: float,
        rng: NumCosmoMath.RNG,
    ) -> bool: ...
    def get_lim(self) -> typing.Tuple[float, float]: ...
    def get_use_true_z(self) -> bool: ...
    @classmethod
    def new(
        cls, sdz: GalaxySDTrueRedshift, zp_min: float, zp_max: float
    ) -> GalaxySDObsRedshiftGauss: ...
    def ref(self) -> GalaxySDObsRedshiftGauss: ...
    def set_lim(self, zp_min: float, zp_max: float) -> None: ...
    def set_use_true_z(self, use_true_z: bool) -> None: ...

class GalaxySDObsRedshiftGaussClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        GalaxySDObsRedshiftGaussClass()
    """

    parent_class: GalaxySDObsRedshiftClass = ...

class GalaxySDObsRedshiftIntegrand(GObject.GBoxed):
    r"""
    :Constructors:

    ::

        GalaxySDObsRedshiftIntegrand()
        new(func:NumCosmo.GalaxySDObsRedshiftIntegrandFunc, callback_data_free:NumCosmo.GalaxySDObsRedshiftIntegrandFreeData, callback_data_copy:NumCosmo.GalaxySDObsRedshiftIntegrandCopyData, callback_data_prepare:NumCosmo.GalaxySDObsRedshiftIntegrandPrepareData, callback_data=None) -> NumCosmo.GalaxySDObsRedshiftIntegrand
    """

    func: typing.Callable[[None, float, GalaxySDObsRedshiftData], float] = ...
    callback_data_free: typing.Callable[[None], None] = ...
    callback_data_copy: typing.Callable[[None], None] = ...
    callback_data_prepare: typing.Callable[[None, NumCosmoMath.MSet], None] = ...
    callback_data: None = ...
    def copy(self) -> GalaxySDObsRedshiftIntegrand: ...
    def eval(self, z: float, data: GalaxySDObsRedshiftData) -> float: ...
    def free(self) -> None: ...
    @classmethod
    def new(
        cls,
        func: typing.Callable[..., float],
        callback_data_free: typing.Callable[..., None],
        callback_data_copy: typing.Callable[..., None],
        callback_data_prepare: typing.Callable[..., None],
        *callback_data: typing.Any,
    ) -> GalaxySDObsRedshiftIntegrand: ...
    def prepare(self, mset: NumCosmoMath.MSet) -> None: ...

class GalaxySDObsRedshiftPz(GalaxySDObsRedshift):
    r"""
    :Constructors:

    ::

        GalaxySDObsRedshiftPz(**properties)
        new() -> NumCosmo.GalaxySDObsRedshiftPz

    Object NcGalaxySDObsRedshiftPz

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    def __init__(
        self,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def clear(gsdorpz: GalaxySDObsRedshiftPz) -> None: ...
    def data_get(self, data: GalaxySDObsRedshiftData) -> NumCosmoMath.Spline: ...
    def data_set(
        self, data: GalaxySDObsRedshiftData, spline: NumCosmoMath.Spline
    ) -> None: ...
    def free(self) -> None: ...
    def gen(
        self,
        mset: NumCosmoMath.MSet,
        data: GalaxySDObsRedshiftData,
        rng: NumCosmoMath.RNG,
    ) -> None: ...
    @classmethod
    def new(cls) -> GalaxySDObsRedshiftPz: ...
    def prepare(self, data: GalaxySDObsRedshiftData) -> None: ...
    def ref(self) -> GalaxySDObsRedshiftPz: ...

class GalaxySDObsRedshiftPzClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        GalaxySDObsRedshiftPzClass()
    """

    parent_class: GalaxySDObsRedshiftClass = ...

class GalaxySDObsRedshiftSpec(GalaxySDObsRedshift):
    r"""
    :Constructors:

    ::

        GalaxySDObsRedshiftSpec(**properties)
        new(sdz:NumCosmo.GalaxySDTrueRedshift, z_min:float, z_max:float) -> NumCosmo.GalaxySDObsRedshiftSpec

    Object NcGalaxySDObsRedshiftSpec

    Properties from NcGalaxySDObsRedshiftSpec:
      lim -> NcmDTuple2: lim
        Galaxy sample photometric redshift distribution limits

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        lim: NumCosmoMath.DTuple2
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    def __init__(
        self,
        lim: NumCosmoMath.DTuple2 = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def clear(gsdorspec: GalaxySDObsRedshiftSpec) -> None: ...
    def free(self) -> None: ...
    def gen(
        self,
        mset: NumCosmoMath.MSet,
        data: GalaxySDObsRedshiftData,
        rng: NumCosmoMath.RNG,
    ) -> None: ...
    def gen1(
        self,
        mset: NumCosmoMath.MSet,
        data: GalaxySDObsRedshiftData,
        rng: NumCosmoMath.RNG,
    ) -> bool: ...
    def get_lim(self) -> typing.Tuple[float, float]: ...
    @classmethod
    def new(
        cls, sdz: GalaxySDTrueRedshift, z_min: float, z_max: float
    ) -> GalaxySDObsRedshiftSpec: ...
    def ref(self) -> GalaxySDObsRedshiftSpec: ...
    def set_lim(self, z_min: float, z_max: float) -> None: ...

class GalaxySDObsRedshiftSpecClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        GalaxySDObsRedshiftSpecClass()
    """

    parent_class: GalaxySDObsRedshiftClass = ...

class GalaxySDPosition(NumCosmoMath.Model):
    r"""
    :Constructors:

    ::

        GalaxySDPosition(**properties)

    Object NcGalaxySDPosition

    Properties from NcGalaxySDPosition:
      ra-lim -> NcmDTuple2: ra-lim
        Galaxy sample right ascension distribution limits
      dec-lim -> NcmDTuple2: dec-lim
        Galaxy sample declination distribution limits

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        dec_lim: NumCosmoMath.DTuple2
        ra_lim: NumCosmoMath.DTuple2
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: NumCosmoMath.Model = ...
    def __init__(
        self,
        dec_lim: NumCosmoMath.DTuple2 = ...,
        ra_lim: NumCosmoMath.DTuple2 = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def clear(gsdp: GalaxySDPosition) -> None: ...
    def do_data_init(
        self, sdz_data: GalaxySDObsRedshiftData, data: GalaxySDPositionData
    ) -> None: ...
    def do_gen(self, data: GalaxySDPositionData, rng: NumCosmoMath.RNG) -> None: ...
    def do_get_dec_lim(self) -> typing.Tuple[bool, float, float]: ...
    def do_get_ra_lim(self) -> typing.Tuple[bool, float, float]: ...
    def do_integ(self) -> GalaxySDPositionIntegrand: ...
    def do_set_dec_lim(self, dec_min: float, dec_max: float) -> bool: ...
    def do_set_ra_lim(self, ra_min: float, ra_max: float) -> bool: ...
    def free(self) -> None: ...
    def gen(self, data: GalaxySDPositionData, rng: NumCosmoMath.RNG) -> None: ...
    def get_dec_lim(self) -> typing.Tuple[bool, float, float]: ...
    def get_ra_lim(self) -> typing.Tuple[bool, float, float]: ...
    @staticmethod
    def id() -> int: ...
    def integ(self) -> GalaxySDPositionIntegrand: ...
    def ref(self) -> GalaxySDPosition: ...
    def set_dec_lim(self, dec_min: float, dec_max: float) -> bool: ...
    def set_ra_lim(self, ra_min: float, ra_max: float) -> bool: ...

class GalaxySDPositionClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        GalaxySDPositionClass()
    """

    parent_class: NumCosmoMath.ModelClass = ...
    gen: typing.Callable[
        [GalaxySDPosition, GalaxySDPositionData, NumCosmoMath.RNG], None
    ] = ...
    integ: typing.Callable[[GalaxySDPosition], GalaxySDPositionIntegrand] = ...
    set_ra_lim: typing.Callable[[GalaxySDPosition, float, float], bool] = ...
    set_dec_lim: typing.Callable[[GalaxySDPosition, float, float], bool] = ...
    get_ra_lim: typing.Callable[
        [GalaxySDPosition], typing.Tuple[bool, float, float]
    ] = ...
    get_dec_lim: typing.Callable[
        [GalaxySDPosition], typing.Tuple[bool, float, float]
    ] = ...
    data_init: typing.Callable[
        [GalaxySDPosition, GalaxySDObsRedshiftData, GalaxySDPositionData], None
    ] = ...
    padding: list[None] = ...

class GalaxySDPositionData(GObject.GBoxed):
    r"""
    :Constructors:

    ::

        GalaxySDPositionData()
        new(gsdp:NumCosmo.GalaxySDPosition, sdz_data:NumCosmo.GalaxySDObsRedshiftData) -> NumCosmo.GalaxySDPositionData
    """

    sdz_data: GalaxySDObsRedshiftData = ...
    ra: float = ...
    dec: float = ...
    ldata: None = ...
    ldata_destroy: typing.Callable[[None], None] = ...
    ldata_read_row: typing.Callable[[GalaxySDPositionData, GalaxyWLObs, int], None] = (
        ...
    )
    ldata_write_row: typing.Callable[[GalaxySDPositionData, GalaxyWLObs, int], None] = (
        ...
    )
    ldata_required_columns: None = ...
    ref_count: int = ...
    @classmethod
    def new(
        cls, gsdp: GalaxySDPosition, sdz_data: GalaxySDObsRedshiftData
    ) -> GalaxySDPositionData: ...
    def read_row(self, obs: GalaxyWLObs, i: int) -> None: ...
    def ref(self) -> GalaxySDPositionData: ...
    def required_columns(self) -> list[str]: ...
    def unref(self) -> None: ...
    def write_row(self, obs: GalaxyWLObs, i: int) -> None: ...

class GalaxySDPositionFlat(GalaxySDPosition):
    r"""
    :Constructors:

    ::

        GalaxySDPositionFlat(**properties)
        new(ra_min:float, ra_max:float, dec_min:float, dec_max:float) -> NumCosmo.GalaxySDPositionFlat

    Object NcGalaxySDPositionFlat

    Properties from NcGalaxySDPosition:
      ra-lim -> NcmDTuple2: ra-lim
        Galaxy sample right ascension distribution limits
      dec-lim -> NcmDTuple2: dec-lim
        Galaxy sample declination distribution limits

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        dec_lim: NumCosmoMath.DTuple2
        ra_lim: NumCosmoMath.DTuple2
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    def __init__(
        self,
        dec_lim: NumCosmoMath.DTuple2 = ...,
        ra_lim: NumCosmoMath.DTuple2 = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def clear(gsdpflat: GalaxySDPositionFlat) -> None: ...
    def free(self) -> None: ...
    def gen(
        self, mset: NumCosmoMath.MSet, data: GalaxySDPositionData, rng: NumCosmoMath.RNG
    ) -> None: ...
    @classmethod
    def new(
        cls, ra_min: float, ra_max: float, dec_min: float, dec_max: float
    ) -> GalaxySDPositionFlat: ...
    def ref(self) -> GalaxySDPositionFlat: ...

class GalaxySDPositionFlatClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        GalaxySDPositionFlatClass()
    """

    parent_class: GalaxySDPositionClass = ...

class GalaxySDPositionIntegrand(GObject.GBoxed):
    r"""
    :Constructors:

    ::

        GalaxySDPositionIntegrand()
        new(func:NumCosmo.GalaxySDPositionIntegrandFunc, callback_data_free:NumCosmo.GalaxySDPositionIntegrandFreeData, callback_data_copy:NumCosmo.GalaxySDPositionIntegrandCopyData, callback_data_prepare:NumCosmo.GalaxySDPositionIntegrandPrepareData, callback_data=None) -> NumCosmo.GalaxySDPositionIntegrand
    """

    func: typing.Callable[[None, GalaxySDPositionData], float] = ...
    callback_data_free: typing.Callable[[None], None] = ...
    callback_data_copy: typing.Callable[[None], None] = ...
    callback_data_prepare: typing.Callable[[None, NumCosmoMath.MSet], None] = ...
    callback_data: None = ...
    def copy(self) -> GalaxySDPositionIntegrand: ...
    def eval(self, data: GalaxySDPositionData) -> float: ...
    def free(self) -> None: ...
    @classmethod
    def new(
        cls,
        func: typing.Callable[..., float],
        callback_data_free: typing.Callable[..., None],
        callback_data_copy: typing.Callable[..., None],
        callback_data_prepare: typing.Callable[..., None],
        *callback_data: typing.Any,
    ) -> GalaxySDPositionIntegrand: ...
    def prepare(self, mset: NumCosmoMath.MSet) -> None: ...

class GalaxySDShape(NumCosmoMath.Model):
    r"""
    :Constructors:

    ::

        GalaxySDShape(**properties)

    Object NcGalaxySDShape

    Properties from NcGalaxySDShape:
      ellip-conv -> NcGalaxyWLObsEllipConv: Ellipticity convention
        Weak lensing observables ellipticity convention

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        ellip_conv: GalaxyWLObsEllipConv
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: NumCosmoMath.Model = ...
    def __init__(
        self,
        ellip_conv: GalaxyWLObsEllipConv = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    def apply_shear(
        self,
        g: NumCosmoMath.Complex,
        E: NumCosmoMath.Complex,
        E_obs: NumCosmoMath.Complex,
    ) -> None: ...
    def apply_shear_inv(
        self,
        g: NumCosmoMath.Complex,
        E_obs: NumCosmoMath.Complex,
        E: NumCosmoMath.Complex,
    ) -> None: ...
    @staticmethod
    def clear(gsds: GalaxySDShape) -> None: ...
    def direct_estimate(
        self, mset: NumCosmoMath.MSet, data_array: typing.Sequence[GalaxySDShapeData]
    ) -> typing.Tuple[float, float, float, float, float]: ...
    def do_data_init(
        self, sdpos_data: GalaxySDPositionData, data: GalaxySDShapeData
    ) -> None: ...
    def do_direct_estimate(
        self, mset: NumCosmoMath.MSet, data_array: typing.Sequence[GalaxySDShapeData]
    ) -> typing.Tuple[float, float, float, float, float]: ...
    def do_gen(
        self, mset: NumCosmoMath.MSet, data: GalaxySDShapeData, rng: NumCosmoMath.RNG
    ) -> None: ...
    def do_integ(self) -> GalaxySDShapeIntegrand: ...
    def do_prepare_data_array(
        self, mset: NumCosmoMath.MSet, data_array: typing.Sequence[GalaxySDShapeData]
    ) -> bool: ...
    def free(self) -> None: ...
    def gen(
        self, mset: NumCosmoMath.MSet, data: GalaxySDShapeData, rng: NumCosmoMath.RNG
    ) -> None: ...
    def get_ellip_conv(self) -> GalaxyWLObsEllipConv: ...
    @staticmethod
    def id() -> int: ...
    def integ(self) -> GalaxySDShapeIntegrand: ...
    def lndet_jac(
        self, g: NumCosmoMath.Complex, E_obs: NumCosmoMath.Complex
    ) -> float: ...
    def prepare_data_array(
        self, mset: NumCosmoMath.MSet, data_array: typing.Sequence[GalaxySDShapeData]
    ) -> bool: ...
    def ref(self) -> GalaxySDShape: ...

class GalaxySDShapeClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        GalaxySDShapeClass()
    """

    parent_class: NumCosmoMath.ModelClass = ...
    gen: typing.Callable[
        [GalaxySDShape, NumCosmoMath.MSet, GalaxySDShapeData, NumCosmoMath.RNG], None
    ] = ...
    integ: typing.Callable[[GalaxySDShape], GalaxySDShapeIntegrand] = ...
    prepare_data_array: typing.Callable[
        [GalaxySDShape, NumCosmoMath.MSet, typing.Sequence[GalaxySDShapeData]], bool
    ] = ...
    data_init: typing.Callable[
        [GalaxySDShape, GalaxySDPositionData, GalaxySDShapeData], None
    ] = ...
    direct_estimate: typing.Callable[
        [GalaxySDShape, NumCosmoMath.MSet, typing.Sequence[GalaxySDShapeData]],
        typing.Tuple[float, float, float, float, float],
    ] = ...
    padding: list[None] = ...

class GalaxySDShapeData(GObject.GBoxed):
    r"""
    :Constructors:

    ::

        GalaxySDShapeData()
        new(gsds:NumCosmo.GalaxySDShape, sdpos_data:NumCosmo.GalaxySDPositionData) -> NumCosmo.GalaxySDShapeData
    """

    sdpos_data: GalaxySDPositionData = ...
    coord: GalaxyWLObsCoord = ...
    epsilon_int_1: float = ...
    epsilon_int_2: float = ...
    ldata: None = ...
    ldata_destroy: typing.Callable[[None], None] = ...
    ldata_read_row: typing.Callable[[GalaxySDShapeData, GalaxyWLObs, int], None] = ...
    ldata_write_row: typing.Callable[[GalaxySDShapeData, GalaxyWLObs, int], None] = ...
    ldata_required_columns: None = ...
    ldata_get_radius: typing.Callable[[GalaxySDShapeData], float] = ...
    ref_count: int = ...
    def get_radius(self) -> float: ...
    @classmethod
    def new(
        cls, gsds: GalaxySDShape, sdpos_data: GalaxySDPositionData
    ) -> GalaxySDShapeData: ...
    def read_row(self, obs: GalaxyWLObs, i: int) -> None: ...
    def ref(self) -> GalaxySDShapeData: ...
    def required_columns(self) -> list[str]: ...
    def unref(self) -> None: ...
    def write_row(self, obs: GalaxyWLObs, i: int) -> None: ...

class GalaxySDShapeGauss(GalaxySDShape):
    r"""
    :Constructors:

    ::

        GalaxySDShapeGauss(**properties)
        new(ellip_conv:NumCosmo.GalaxyWLObsEllipConv) -> NumCosmo.GalaxySDShapeGauss

    Object NcGalaxySDShapeGauss

    Properties from NcGalaxySDShapeGauss:
      sigma -> gdouble: sigma
        \sigma
      sigma-fit -> gboolean: sigma-fit
        \sigma:fit

    Properties from NcGalaxySDShape:
      ellip-conv -> NcGalaxyWLObsEllipConv: Ellipticity convention
        Weak lensing observables ellipticity convention

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        sigma: float
        sigma_fit: bool
        ellip_conv: GalaxyWLObsEllipConv
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    def __init__(
        self,
        sigma: float = ...,
        sigma_fit: bool = ...,
        ellip_conv: GalaxyWLObsEllipConv = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def clear(gsdsgauss: GalaxySDShapeGauss) -> None: ...
    def data_get(
        self, data: GalaxySDShapeData
    ) -> typing.Tuple[float, float, float]: ...
    def data_set(
        self,
        data: GalaxySDShapeData,
        epsilon_obs_1: float,
        epsilon_obs_2: float,
        std_noise: float,
    ) -> None: ...
    def free(self) -> None: ...
    def gen(
        self,
        mset: NumCosmoMath.MSet,
        data: GalaxySDShapeData,
        std_noise: float,
        coord: GalaxyWLObsCoord,
        rng: NumCosmoMath.RNG,
    ) -> None: ...
    @classmethod
    def new(cls, ellip_conv: GalaxyWLObsEllipConv) -> GalaxySDShapeGauss: ...
    def ref(self) -> GalaxySDShapeGauss: ...
    @staticmethod
    def sigma_from_std_shape(std_shape: float) -> float: ...
    @staticmethod
    def std_shape_from_sigma(sigma: float) -> float: ...

class GalaxySDShapeGaussClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        GalaxySDShapeGaussClass()
    """

    parent_class: GalaxySDShapeClass = ...

class GalaxySDShapeGaussHSC(GalaxySDShape):
    r"""
    :Constructors:

    ::

        GalaxySDShapeGaussHSC(**properties)
        new(ellip_conv:NumCosmo.GalaxyWLObsEllipConv) -> NumCosmo.GalaxySDShapeGaussHSC

    Object NcGalaxySDShapeGaussHSC

    Properties from NcGalaxySDShape:
      ellip-conv -> NcGalaxyWLObsEllipConv: Ellipticity convention
        Weak lensing observables ellipticity convention

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        ellip_conv: GalaxyWLObsEllipConv
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    def __init__(
        self,
        ellip_conv: GalaxyWLObsEllipConv = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def clear(gsdshsc: GalaxySDShapeGaussHSC) -> None: ...
    def data_get(
        self, data: GalaxySDShapeData
    ) -> typing.Tuple[float, float, float, float, float, float, float]: ...
    def data_set(
        self,
        data: GalaxySDShapeData,
        epsilon_obs_1: float,
        epsilon_obs_2: float,
        std_shape: float,
        std_noise: float,
        c1: float,
        c2: float,
        m: float,
    ) -> None: ...
    def free(self) -> None: ...
    def gen(
        self,
        mset: NumCosmoMath.MSet,
        data: GalaxySDShapeData,
        std_shape: float,
        std_noise: float,
        c1: float,
        c2: float,
        m: float,
        coord: GalaxyWLObsCoord,
        rng: NumCosmoMath.RNG,
    ) -> None: ...
    @classmethod
    def new(cls, ellip_conv: GalaxyWLObsEllipConv) -> GalaxySDShapeGaussHSC: ...
    def ref(self) -> GalaxySDShapeGaussHSC: ...

class GalaxySDShapeGaussHSCClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        GalaxySDShapeGaussHSCClass()
    """

    parent_class: GalaxySDShapeClass = ...

class GalaxySDShapeIntegrand(GObject.GBoxed):
    r"""
    :Constructors:

    ::

        GalaxySDShapeIntegrand()
        new(func:NumCosmo.GalaxySDShapeIntegrandFunc, callback_data_free:NumCosmo.GalaxySDShapeIntegrandFreeData, callback_data_copy:NumCosmo.GalaxySDShapeIntegrandCopyData, callback_data_prepare:NumCosmo.GalaxySDShapeIntegrandPrepareData, callback_data=None) -> NumCosmo.GalaxySDShapeIntegrand
    """

    func: typing.Callable[[None, float, GalaxySDShapeData], float] = ...
    callback_data_free: typing.Callable[[None], None] = ...
    callback_data_copy: typing.Callable[[None], None] = ...
    callback_data_prepare: typing.Callable[[None, NumCosmoMath.MSet], None] = ...
    callback_data: None = ...
    def copy(self) -> GalaxySDShapeIntegrand: ...
    def eval(self, z: float, data: GalaxySDShapeData) -> float: ...
    def free(self) -> None: ...
    @classmethod
    def new(
        cls,
        func: typing.Callable[..., float],
        callback_data_free: typing.Callable[..., None],
        callback_data_copy: typing.Callable[..., None],
        callback_data_prepare: typing.Callable[..., None],
        *callback_data: typing.Any,
    ) -> GalaxySDShapeIntegrand: ...
    def prepare(self, mset: NumCosmoMath.MSet) -> None: ...

class GalaxySDTrueRedshift(NumCosmoMath.Model):
    r"""
    :Constructors:

    ::

        GalaxySDTrueRedshift(**properties)

    Object NcGalaxySDTrueRedshift

    Properties from NcGalaxySDTrueRedshift:
      lim -> NcmDTuple2: lim
        Galaxy sample redshift distribution limits

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        lim: NumCosmoMath.DTuple2
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: NumCosmoMath.Model = ...
    def __init__(
        self,
        lim: NumCosmoMath.DTuple2 = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def clear(gsdtr: GalaxySDTrueRedshift) -> None: ...
    def do_gen(self, rng: NumCosmoMath.RNG) -> float: ...
    def do_get_lim(self, z_min: float, z_max: float) -> None: ...
    def do_integ(self, z: float) -> float: ...
    def do_set_lim(self, z_min: float, z_max: float) -> None: ...
    def free(self) -> None: ...
    def gen(self, rng: NumCosmoMath.RNG) -> float: ...
    def get_lim(self, z_min: float, z_max: float) -> None: ...
    @staticmethod
    def id() -> int: ...
    def integ(self, z: float) -> float: ...
    def ref(self) -> GalaxySDTrueRedshift: ...
    def set_lim(self, z_min: float, z_max: float) -> None: ...

class GalaxySDTrueRedshiftClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        GalaxySDTrueRedshiftClass()
    """

    parent_class: NumCosmoMath.ModelClass = ...
    gen: typing.Callable[[GalaxySDTrueRedshift, NumCosmoMath.RNG], float] = ...
    integ: typing.Callable[[GalaxySDTrueRedshift, float], float] = ...
    set_lim: typing.Callable[[GalaxySDTrueRedshift, float, float], None] = ...
    get_lim: typing.Callable[[GalaxySDTrueRedshift, float, float], None] = ...
    padding: list[None] = ...

class GalaxySDTrueRedshiftLSSTSRD(GalaxySDTrueRedshift):
    r"""
    :Constructors:

    ::

        GalaxySDTrueRedshiftLSSTSRD(**properties)
        new() -> NumCosmo.GalaxySDTrueRedshiftLSSTSRD
        new_y10() -> NumCosmo.GalaxySDTrueRedshiftLSSTSRD

    Object NcGalaxySDTrueRedshiftLSSTSRD

    Properties from NcGalaxySDTrueRedshiftLSSTSRD:
      alpha -> gdouble: alpha
        \alpha
      beta -> gdouble: beta
        \beta
      z0 -> gdouble: z0
        z_0
      alpha-fit -> gboolean: alpha-fit
        \alpha:fit
      beta-fit -> gboolean: beta-fit
        \beta:fit
      z0-fit -> gboolean: z0-fit
        z_0:fit

    Properties from NcGalaxySDTrueRedshift:
      lim -> NcmDTuple2: lim
        Galaxy sample redshift distribution limits

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        alpha: float
        alpha_fit: bool
        beta: float
        beta_fit: bool
        z0: float
        z0_fit: bool
        lim: NumCosmoMath.DTuple2
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    def __init__(
        self,
        alpha: float = ...,
        alpha_fit: bool = ...,
        beta: float = ...,
        beta_fit: bool = ...,
        z0: float = ...,
        z0_fit: bool = ...,
        lim: NumCosmoMath.DTuple2 = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def clear(gsdtrlsst: GalaxySDTrueRedshiftLSSTSRD) -> None: ...
    def free(self) -> None: ...
    @classmethod
    def new(cls) -> GalaxySDTrueRedshiftLSSTSRD: ...
    @classmethod
    def new_y10(cls) -> GalaxySDTrueRedshiftLSSTSRD: ...
    def ref(self) -> GalaxySDTrueRedshiftLSSTSRD: ...

class GalaxySDTrueRedshiftLSSTSRDClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        GalaxySDTrueRedshiftLSSTSRDClass()
    """

    parent_class: GalaxySDTrueRedshiftClass = ...

class GalaxySelfunc(GObject.Object):
    r"""
    :Constructors:

    ::

        GalaxySelfunc(**properties)
        new(nshells:int) -> NumCosmo.GalaxySelfunc

    Object NcGalaxySelfunc

    Properties from NcGalaxySelfunc:
      nshells -> guint: nshells
        Galaxy survey number of redshift shells
      shell-splines -> NcmObjArray: shell-splines
        Galaxy survey shell splines

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        nshells: int
        shell_splines: NumCosmoMath.ObjArray

    props: Props = ...
    parent_instance: GObject.Object = ...
    priv: GalaxySelfuncPrivate = ...
    def __init__(
        self, nshells: int = ..., shell_splines: NumCosmoMath.ObjArray = ...
    ) -> None: ...
    @staticmethod
    def clear(gsf: GalaxySelfunc) -> None: ...
    def eval(self, shell: int, z: float) -> float: ...
    def free(self) -> None: ...
    def get_nshells(self) -> int: ...
    def get_shell_splines(self) -> NumCosmoMath.ObjArray: ...
    def get_zmax(self, shell: int) -> float: ...
    def get_zmean(self, shell: int) -> float: ...
    def get_zmin(self, shell: int) -> float: ...
    def load_from_txts(
        self, prefix: str, suffix: typing.Optional[str] = None
    ) -> None: ...
    @classmethod
    def new(cls, nshells: int) -> GalaxySelfunc: ...
    def ref(self) -> GalaxySelfunc: ...
    def set_nshells(self, nshells: int) -> None: ...
    def set_shell_splines(self, dNdz_a: NumCosmoMath.ObjArray) -> None: ...

class GalaxySelfuncClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        GalaxySelfuncClass()
    """

    parent_class: GObject.ObjectClass = ...

class GalaxySelfuncPrivate(GObject.GPointer): ...

class GalaxyWLObs(GObject.Object):
    r"""
    :Constructors:

    ::

        GalaxyWLObs(**properties)
        new(ellip_conv:NumCosmo.GalaxyWLObsEllipConv, coord:NumCosmo.GalaxyWLObsCoord, nrows:int, col_names:list) -> NumCosmo.GalaxyWLObs

    Object NcGalaxyWLObs

    Properties from NcGalaxyWLObs:
      data -> NcmMatrix: Data
        Weak lensing observation data matrix
      pz -> NcmObjDictInt: P(z)
        P(z) splines
      columns -> GStrv: Colunms
        Data columns names
      coord -> NcGalaxyWLObsCoord: Coordinate system
        Coordinate system used to store the data
      ellip-conv -> NcGalaxyWLObsEllipConv: Ellipticity convention
        Weak lensing observables ellipticity convention
      len -> guint: Length
        Number of data rows

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        columns: list[str]
        coord: GalaxyWLObsCoord
        data: NumCosmoMath.Matrix
        ellip_conv: GalaxyWLObsEllipConv
        len: int
        pz: NumCosmoMath.ObjDictInt

    props: Props = ...
    def __init__(
        self,
        columns: typing.Sequence[str] = ...,
        coord: GalaxyWLObsCoord = ...,
        data: NumCosmoMath.Matrix = ...,
        ellip_conv: GalaxyWLObsEllipConv = ...,
        len: int = ...,
        pz: NumCosmoMath.ObjDictInt = ...,
    ) -> None: ...
    @staticmethod
    def clear(obs: GalaxyWLObs) -> None: ...
    def free(self) -> None: ...
    def get(self, col: str, i: int) -> float: ...
    def get_coord(self) -> GalaxyWLObsCoord: ...
    def get_ellip_conv(self) -> GalaxyWLObsEllipConv: ...
    def get_index(self, col: str, i: int) -> bool: ...
    def len(self) -> int: ...
    @classmethod
    def new(
        cls,
        ellip_conv: GalaxyWLObsEllipConv,
        coord: GalaxyWLObsCoord,
        nrows: int,
        col_names: typing.Sequence[str],
    ) -> GalaxyWLObs: ...
    def peek_columns(self) -> list[str]: ...
    def peek_pz(self, i: int) -> NumCosmoMath.Spline: ...
    def ref(self) -> GalaxyWLObs: ...
    def set(self, col: str, i: int, val: float) -> None: ...
    def set_coord(self, coord: GalaxyWLObsCoord) -> None: ...
    def set_pz(self, i: int, pz: NumCosmoMath.Spline) -> None: ...

class GalaxyWLObsClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        GalaxyWLObsClass()
    """

    parent_class: GObject.ObjectClass = ...

class GalaxyWLObsPrivate(GObject.GPointer): ...

class GrowthFunc(GObject.Object):
    r"""
    :Constructors:

    ::

        GrowthFunc(**properties)
        new() -> NumCosmo.GrowthFunc

    Object NcGrowthFunc

    Properties from NcGrowthFunc:
      reltol -> gdouble: reltol
        Relative tolerance
      abstol -> gdouble: abstol
        Absolute tolerance tolerance
      x-i -> gdouble: x-i
        Initial value for $x_i$

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        abstol: float
        reltol: float
        x_i: float

    props: Props = ...
    parent_instance: GObject.Object = ...
    priv: GrowthFuncPrivate = ...
    s: NumCosmoMath.Spline = ...
    Da0: float = ...
    def __init__(
        self, abstol: float = ..., reltol: float = ..., x_i: float = ...
    ) -> None: ...
    @staticmethod
    def clear(gf: GrowthFunc) -> None: ...
    def eval(self, cosmo: HICosmo, z: float) -> float: ...
    def eval_both(self, cosmo: HICosmo, z: float) -> typing.Tuple[float, float]: ...
    def eval_deriv(self, cosmo: HICosmo, z: float) -> float: ...
    def free(self) -> None: ...
    def get_abstol(self) -> float: ...
    def get_dust_norma_Da0(self) -> float: ...
    def get_reltol(self) -> float: ...
    def get_x_i(self) -> float: ...
    @classmethod
    def new(cls) -> GrowthFunc: ...
    def prepare(self, cosmo: HICosmo) -> None: ...
    def prepare_if_needed(self, cosmo: HICosmo) -> None: ...
    def ref(self) -> GrowthFunc: ...
    def set_abstol(self, abstol: float) -> None: ...
    def set_reltol(self, reltol: float) -> None: ...
    def set_x_i(self, x_i: float) -> None: ...

class GrowthFuncClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        GrowthFuncClass()
    """

    parent_class: GObject.ObjectClass = ...

class GrowthFuncPrivate(GObject.GPointer): ...

class HICosmo(NumCosmoMath.Model):
    r"""
    :Constructors:

    ::

        HICosmo(**properties)

    Object NcHICosmo

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: NumCosmoMath.Model = ...
    is_eternal: bool = ...
    prim: HIPrim = ...
    reion: HIReion = ...
    T: int = ...
    s: int = ...
    Tmin: int = ...
    smin: int = ...
    def __init__(
        self,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    def Dc(self, z: float) -> float: ...
    def E(self, z: float) -> float: ...
    def E2(self, z: float) -> float: ...
    def E2Omega_b(self, z: float) -> float: ...
    def E2Omega_c(self, z: float) -> float: ...
    def E2Omega_g(self, z: float) -> float: ...
    def E2Omega_k(self, z: float) -> float: ...
    def E2Omega_m(self, z: float) -> float: ...
    def E2Omega_mnu(self, z: float) -> float: ...
    def E2Omega_mnu_n(self, n: int, z: float) -> float: ...
    def E2Omega_nu(self, z: float) -> float: ...
    def E2Omega_r(self, z: float) -> float: ...
    def E2Omega_t(self, z: float) -> float: ...
    def E2Press_mnu(self, z: float) -> float: ...
    def E2Press_mnu_n(self, n: int, z: float) -> float: ...
    def Em2(self, z: float) -> float: ...
    def H(self, z: float) -> float: ...
    def H0(self) -> float: ...
    def H_number_density(self) -> float: ...
    def He_number_density(self) -> float: ...
    def MassNuInfo(self, nu_i: int) -> typing.Tuple[float, float, float, float]: ...
    def NMassNu(self) -> int: ...
    def Neff(self) -> float: ...
    def Omega_b0(self) -> float: ...
    def Omega_b0h2(self) -> float: ...
    def Omega_c0(self) -> float: ...
    def Omega_c0h2(self) -> float: ...
    def Omega_g0(self) -> float: ...
    def Omega_g0h2(self) -> float: ...
    def Omega_k0(self) -> float: ...
    def Omega_m0(self) -> float: ...
    def Omega_m0h2(self) -> float: ...
    def Omega_mnu0(self) -> float: ...
    def Omega_mnu0_n(self, n: int) -> float: ...
    def Omega_mnu0h2(self) -> float: ...
    def Omega_nu0(self) -> float: ...
    def Omega_nu0h2(self) -> float: ...
    def Omega_r0(self) -> float: ...
    def Omega_r0h2(self) -> float: ...
    def Omega_t0(self) -> float: ...
    def Press_mnu0(self) -> float: ...
    def Press_mnu0_n(self, n: int) -> float: ...
    def RH_Mpc(self) -> float: ...
    def RH_planck(self) -> float: ...
    def T_gamma0(self) -> float: ...
    def XHe(self) -> float: ...
    def Yp_1H(self) -> float: ...
    def Yp_4He(self) -> float: ...
    def abs_alpha(self, x: float) -> float: ...
    def as_drag(self) -> float: ...
    def baryon_density(self) -> float: ...
    def bgp_cs2(self, z: float) -> float: ...
    @staticmethod
    def clear(cosmo: HICosmo) -> None: ...
    def crit_density(self) -> float: ...
    def d2E2_dz2(self, z: float) -> float: ...
    def dE2_dz(self, z: float) -> float: ...
    def dH_dz(self, z: float) -> float: ...
    def dec(self, z: float) -> float: ...
    def dec_min(self, z_max: float) -> typing.Tuple[float, float]: ...
    def do_Dc(self, z: float) -> float: ...
    def do_E2(self, z: float) -> float: ...
    def do_E2Omega_b(self, z: float) -> float: ...
    def do_E2Omega_c(self, z: float) -> float: ...
    def do_E2Omega_g(self, z: float) -> float: ...
    def do_E2Omega_m(self, z: float) -> float: ...
    def do_E2Omega_mnu(self, z: float) -> float: ...
    def do_E2Omega_mnu_n(self, n: int, z: float) -> float: ...
    def do_E2Omega_nu(self, z: float) -> float: ...
    def do_E2Omega_r(self, z: float) -> float: ...
    def do_E2Omega_t(self, z: float) -> float: ...
    def do_E2Press_mnu(self, z: float) -> float: ...
    def do_E2Press_mnu_n(self, n: int, z: float) -> float: ...
    def do_H0(self) -> float: ...
    def do_MassNuInfo(self, nu_i: int) -> typing.Tuple[float, float, float, float]: ...
    def do_NMassNu(self) -> int: ...
    def do_Omega_b0(self) -> float: ...
    def do_Omega_c0(self) -> float: ...
    def do_Omega_g0(self) -> float: ...
    def do_Omega_m0(self) -> float: ...
    def do_Omega_mnu0(self) -> float: ...
    def do_Omega_mnu0_n(self, n: int) -> float: ...
    def do_Omega_nu0(self) -> float: ...
    def do_Omega_r0(self) -> float: ...
    def do_Omega_t0(self) -> float: ...
    def do_Press_mnu0(self) -> float: ...
    def do_Press_mnu0_n(self, n: int) -> float: ...
    def do_T_gamma0(self) -> float: ...
    def do_Yp_4He(self) -> float: ...
    def do_as_drag(self) -> float: ...
    def do_bgp_cs2(self, z: float) -> float: ...
    def do_d2E2_dz2(self, z: float) -> float: ...
    def do_dE2_dz(self, z: float) -> float: ...
    def do_get_bg_var(self, t: float, bg_var: HIPertBGVar) -> None: ...
    def do_xb(self) -> float: ...
    def do_z_lss(self) -> float: ...
    def free(self) -> None: ...
    def get_bg_var(self, t: float, bg_var: HIPertBGVar) -> None: ...
    def h(self) -> float: ...
    def h2(self) -> float: ...
    @staticmethod
    def id() -> int: ...
    def j(self, z: float) -> float: ...
    def kinetic_w(self, z: float) -> float: ...
    @staticmethod
    def log_all_models(parent: typing.Type[typing.Any]) -> None: ...
    def mqE2(self, z: float) -> float: ...
    def mqE2_max(self, z_max: float) -> typing.Tuple[float, float]: ...
    def nec(self, z: float) -> float: ...
    def peek_prim(self) -> HIPrim: ...
    def peek_reion(self) -> HIReion: ...
    @staticmethod
    def priors_stub() -> None: ...
    def q(self, z: float) -> float: ...
    def q_min(self, z_max: float) -> typing.Tuple[float, float]: ...
    def qp(self, z: float) -> float: ...
    def ref(self) -> HICosmo: ...
    def sigma8(self, psf: NumCosmoMath.PowspecFilter) -> float: ...
    def wec(self, z: float) -> float: ...
    def x_alpha(self, alpha: float) -> float: ...
    def xb(self) -> float: ...
    def z_lss(self) -> float: ...
    def zt(self, z_max: float) -> float: ...

class HICosmoClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoClass()
    """

    parent_class: NumCosmoMath.ModelClass = ...
    H0: typing.Callable[[HICosmo], float] = ...
    Omega_b0: typing.Callable[[HICosmo], float] = ...
    Omega_c0: typing.Callable[[HICosmo], float] = ...
    Omega_g0: typing.Callable[[HICosmo], float] = ...
    Omega_nu0: typing.Callable[[HICosmo], float] = ...
    Omega_mnu0: typing.Callable[[HICosmo], float] = ...
    Press_mnu0: typing.Callable[[HICosmo], float] = ...
    Omega_m0: typing.Callable[[HICosmo], float] = ...
    Omega_r0: typing.Callable[[HICosmo], float] = ...
    Omega_t0: typing.Callable[[HICosmo], float] = ...
    T_gamma0: typing.Callable[[HICosmo], float] = ...
    Yp_4He: typing.Callable[[HICosmo], float] = ...
    z_lss: typing.Callable[[HICosmo], float] = ...
    as_drag: typing.Callable[[HICosmo], float] = ...
    xb: typing.Callable[[HICosmo], float] = ...
    Omega_mnu0_n: typing.Callable[[HICosmo, int], float] = ...
    Press_mnu0_n: typing.Callable[[HICosmo, int], float] = ...
    E2Omega_b: typing.Callable[[HICosmo, float], float] = ...
    E2Omega_c: typing.Callable[[HICosmo, float], float] = ...
    E2Omega_g: typing.Callable[[HICosmo, float], float] = ...
    E2Omega_nu: typing.Callable[[HICosmo, float], float] = ...
    E2Omega_mnu: typing.Callable[[HICosmo, float], float] = ...
    E2Press_mnu: typing.Callable[[HICosmo, float], float] = ...
    E2Omega_m: typing.Callable[[HICosmo, float], float] = ...
    E2Omega_r: typing.Callable[[HICosmo, float], float] = ...
    E2Omega_t: typing.Callable[[HICosmo, float], float] = ...
    E2: typing.Callable[[HICosmo, float], float] = ...
    dE2_dz: typing.Callable[[HICosmo, float], float] = ...
    d2E2_dz2: typing.Callable[[HICosmo, float], float] = ...
    bgp_cs2: typing.Callable[[HICosmo, float], float] = ...
    Dc: typing.Callable[[HICosmo, float], float] = ...
    E2Omega_mnu_n: typing.Callable[[HICosmo, int, float], float] = ...
    E2Press_mnu_n: typing.Callable[[HICosmo, int, float], float] = ...
    NMassNu: typing.Callable[[HICosmo], int] = ...
    MassNuInfo: typing.Callable[
        [HICosmo, int], typing.Tuple[float, float, float, float]
    ] = ...
    get_bg_var: typing.Callable[[HICosmo, float, HIPertBGVar], None] = ...

class HICosmoDE(HICosmo):
    r"""
    :Constructors:

    ::

        HICosmoDE(**properties)

    Object NcHICosmoDE

    Properties from NcHICosmoDE:
      CCL-comp -> gboolean: CCL-comp
        Whether to use CCL compatible mode
      H0 -> gdouble: H0
        H_0
      Omegac -> gdouble: Omegac
        \Omega_{c0}
      Omegax -> gdouble: Omegax
        \Omega_{x0}
      Tgamma0 -> gdouble: Tgamma0
        T_{\gamma0}
      Yp -> gdouble: Yp
        Y_p
      ENnu -> gdouble: ENnu
        N_\nu
      Omegab -> gdouble: Omegab
        \Omega_{b0}
      massnu -> NcmVector: massnu
        m_\nu
      Tnu -> NcmVector: Tnu
        T_{\nu0}
      munu -> NcmVector: munu
        \mu_{\nu}
      gnu -> NcmVector: gnu
        g_{\nu}
      massnu-length -> guint: massnu-length
        m_\nu:length
      Tnu-length -> guint: Tnu-length
        T_{\nu0}:length
      munu-length -> guint: munu-length
        \mu_{\nu}:length
      gnu-length -> guint: gnu-length
        g_{\nu}:length
      H0-fit -> gboolean: H0-fit
        H_0:fit
      Omegac-fit -> gboolean: Omegac-fit
        \Omega_{c0}:fit
      Omegax-fit -> gboolean: Omegax-fit
        \Omega_{x0}:fit
      Tgamma0-fit -> gboolean: Tgamma0-fit
        T_{\gamma0}:fit
      Yp-fit -> gboolean: Yp-fit
        Y_p:fit
      ENnu-fit -> gboolean: ENnu-fit
        N_\nu:fit
      Omegab-fit -> gboolean: Omegab-fit
        \Omega_{b0}:fit
      massnu-fit -> GVariant: massnu-fit
        m_\nu:fit
      Tnu-fit -> GVariant: Tnu-fit
        T_{\nu0}:fit
      munu-fit -> GVariant: munu-fit
        \mu_{\nu}:fit
      gnu-fit -> GVariant: gnu-fit
        g_{\nu}:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        CCL_comp: bool
        ENnu: float
        ENnu_fit: bool
        H0: float
        H0_fit: bool
        Omegab: float
        Omegab_fit: bool
        Omegac: float
        Omegac_fit: bool
        Omegax: float
        Omegax_fit: bool
        Tgamma0: float
        Tgamma0_fit: bool
        Tnu: NumCosmoMath.Vector
        Tnu_fit: GLib.Variant
        Tnu_length: int
        Yp: float
        Yp_fit: bool
        gnu: NumCosmoMath.Vector
        gnu_fit: GLib.Variant
        gnu_length: int
        massnu: NumCosmoMath.Vector
        massnu_fit: GLib.Variant
        massnu_length: int
        munu: NumCosmoMath.Vector
        munu_fit: GLib.Variant
        munu_length: int
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: HICosmo = ...
    priv: HICosmoDEPrivate = ...
    def __init__(
        self,
        CCL_comp: bool = ...,
        ENnu: float = ...,
        ENnu_fit: bool = ...,
        H0: float = ...,
        H0_fit: bool = ...,
        Omegab: float = ...,
        Omegab_fit: bool = ...,
        Omegac: float = ...,
        Omegac_fit: bool = ...,
        Omegax: float = ...,
        Omegax_fit: bool = ...,
        Tgamma0: float = ...,
        Tgamma0_fit: bool = ...,
        Tnu: NumCosmoMath.Vector = ...,
        Tnu_fit: GLib.Variant = ...,
        Tnu_length: int = ...,
        Yp: float = ...,
        Yp_fit: bool = ...,
        gnu: NumCosmoMath.Vector = ...,
        gnu_fit: GLib.Variant = ...,
        gnu_length: int = ...,
        massnu: NumCosmoMath.Vector = ...,
        massnu_fit: GLib.Variant = ...,
        massnu_length: int = ...,
        munu: NumCosmoMath.Vector = ...,
        munu_fit: GLib.Variant = ...,
        munu_length: int = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    def E2Omega_de(self, z: float) -> float: ...
    def E2Omega_de_onepw(self, z: float) -> float: ...
    def cmb_params(self) -> None: ...
    def d2E2Omega_de_dz2(self, z: float) -> float: ...
    def dE2Omega_de_dz(self, z: float) -> float: ...
    def do_E2Omega_de(self, z: float) -> float: ...
    def do_d2E2Omega_de_dz2(self, z: float) -> float: ...
    def do_dE2Omega_de_dz(self, z: float) -> float: ...
    def do_w_de(self, z: float) -> float: ...
    @staticmethod
    def new_add_bbn(lh: NumCosmoMath.Likelihood) -> None: ...
    def omega_x2omega_k(self) -> None: ...
    def set_wmap5_params(self) -> None: ...
    def w_de(self, z: float) -> float: ...

class HICosmoDEClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoDEClass()
    """

    parent_class: HICosmoClass = ...
    E2Omega_de: typing.Callable[[HICosmoDE, float], float] = ...
    dE2Omega_de_dz: typing.Callable[[HICosmoDE, float], float] = ...
    d2E2Omega_de_dz2: typing.Callable[[HICosmoDE, float], float] = ...
    w_de: typing.Callable[[HICosmoDE, float], float] = ...

class HICosmoDECpl(HICosmoDE):
    r"""
    :Constructors:

    ::

        HICosmoDECpl(**properties)
        new() -> NumCosmo.HICosmoDECpl

    Object NcHICosmoDECpl

    Properties from NcHICosmoDECpl:
      w0 -> gdouble: w0
        w_0
      w1 -> gdouble: w1
        w_1
      w0-fit -> gboolean: w0-fit
        w_0:fit
      w1-fit -> gboolean: w1-fit
        w_1:fit

    Properties from NcHICosmoDE:
      CCL-comp -> gboolean: CCL-comp
        Whether to use CCL compatible mode
      H0 -> gdouble: H0
        H_0
      Omegac -> gdouble: Omegac
        \Omega_{c0}
      Omegax -> gdouble: Omegax
        \Omega_{x0}
      Tgamma0 -> gdouble: Tgamma0
        T_{\gamma0}
      Yp -> gdouble: Yp
        Y_p
      ENnu -> gdouble: ENnu
        N_\nu
      Omegab -> gdouble: Omegab
        \Omega_{b0}
      massnu -> NcmVector: massnu
        m_\nu
      Tnu -> NcmVector: Tnu
        T_{\nu0}
      munu -> NcmVector: munu
        \mu_{\nu}
      gnu -> NcmVector: gnu
        g_{\nu}
      massnu-length -> guint: massnu-length
        m_\nu:length
      Tnu-length -> guint: Tnu-length
        T_{\nu0}:length
      munu-length -> guint: munu-length
        \mu_{\nu}:length
      gnu-length -> guint: gnu-length
        g_{\nu}:length
      H0-fit -> gboolean: H0-fit
        H_0:fit
      Omegac-fit -> gboolean: Omegac-fit
        \Omega_{c0}:fit
      Omegax-fit -> gboolean: Omegax-fit
        \Omega_{x0}:fit
      Tgamma0-fit -> gboolean: Tgamma0-fit
        T_{\gamma0}:fit
      Yp-fit -> gboolean: Yp-fit
        Y_p:fit
      ENnu-fit -> gboolean: ENnu-fit
        N_\nu:fit
      Omegab-fit -> gboolean: Omegab-fit
        \Omega_{b0}:fit
      massnu-fit -> GVariant: massnu-fit
        m_\nu:fit
      Tnu-fit -> GVariant: Tnu-fit
        T_{\nu0}:fit
      munu-fit -> GVariant: munu-fit
        \mu_{\nu}:fit
      gnu-fit -> GVariant: gnu-fit
        g_{\nu}:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        w0: float
        w0_fit: bool
        w1: float
        w1_fit: bool
        CCL_comp: bool
        ENnu: float
        ENnu_fit: bool
        H0: float
        H0_fit: bool
        Omegab: float
        Omegab_fit: bool
        Omegac: float
        Omegac_fit: bool
        Omegax: float
        Omegax_fit: bool
        Tgamma0: float
        Tgamma0_fit: bool
        Tnu: NumCosmoMath.Vector
        Tnu_fit: GLib.Variant
        Tnu_length: int
        Yp: float
        Yp_fit: bool
        gnu: NumCosmoMath.Vector
        gnu_fit: GLib.Variant
        gnu_length: int
        massnu: NumCosmoMath.Vector
        massnu_fit: GLib.Variant
        massnu_length: int
        munu: NumCosmoMath.Vector
        munu_fit: GLib.Variant
        munu_length: int
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: HICosmoDE = ...
    def __init__(
        self,
        w0: float = ...,
        w0_fit: bool = ...,
        w1: float = ...,
        w1_fit: bool = ...,
        CCL_comp: bool = ...,
        ENnu: float = ...,
        ENnu_fit: bool = ...,
        H0: float = ...,
        H0_fit: bool = ...,
        Omegab: float = ...,
        Omegab_fit: bool = ...,
        Omegac: float = ...,
        Omegac_fit: bool = ...,
        Omegax: float = ...,
        Omegax_fit: bool = ...,
        Tgamma0: float = ...,
        Tgamma0_fit: bool = ...,
        Tnu: NumCosmoMath.Vector = ...,
        Tnu_fit: GLib.Variant = ...,
        Tnu_length: int = ...,
        Yp: float = ...,
        Yp_fit: bool = ...,
        gnu: NumCosmoMath.Vector = ...,
        gnu_fit: GLib.Variant = ...,
        gnu_length: int = ...,
        massnu: NumCosmoMath.Vector = ...,
        massnu_fit: GLib.Variant = ...,
        massnu_length: int = ...,
        munu: NumCosmoMath.Vector = ...,
        munu_fit: GLib.Variant = ...,
        munu_length: int = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @classmethod
    def new(cls) -> HICosmoDECpl: ...

class HICosmoDECplClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoDECplClass()
    """

    parent_class: HICosmoDEClass = ...

class HICosmoDEJbp(HICosmoDE):
    r"""
    :Constructors:

    ::

        HICosmoDEJbp(**properties)
        new() -> NumCosmo.HICosmoDEJbp

    Object NcHICosmoDEJbp

    Properties from NcHICosmoDEJbp:
      w0 -> gdouble: w0
        w_0
      w1 -> gdouble: w1
        w_1
      w0-fit -> gboolean: w0-fit
        w_0:fit
      w1-fit -> gboolean: w1-fit
        w_1:fit

    Properties from NcHICosmoDE:
      CCL-comp -> gboolean: CCL-comp
        Whether to use CCL compatible mode
      H0 -> gdouble: H0
        H_0
      Omegac -> gdouble: Omegac
        \Omega_{c0}
      Omegax -> gdouble: Omegax
        \Omega_{x0}
      Tgamma0 -> gdouble: Tgamma0
        T_{\gamma0}
      Yp -> gdouble: Yp
        Y_p
      ENnu -> gdouble: ENnu
        N_\nu
      Omegab -> gdouble: Omegab
        \Omega_{b0}
      massnu -> NcmVector: massnu
        m_\nu
      Tnu -> NcmVector: Tnu
        T_{\nu0}
      munu -> NcmVector: munu
        \mu_{\nu}
      gnu -> NcmVector: gnu
        g_{\nu}
      massnu-length -> guint: massnu-length
        m_\nu:length
      Tnu-length -> guint: Tnu-length
        T_{\nu0}:length
      munu-length -> guint: munu-length
        \mu_{\nu}:length
      gnu-length -> guint: gnu-length
        g_{\nu}:length
      H0-fit -> gboolean: H0-fit
        H_0:fit
      Omegac-fit -> gboolean: Omegac-fit
        \Omega_{c0}:fit
      Omegax-fit -> gboolean: Omegax-fit
        \Omega_{x0}:fit
      Tgamma0-fit -> gboolean: Tgamma0-fit
        T_{\gamma0}:fit
      Yp-fit -> gboolean: Yp-fit
        Y_p:fit
      ENnu-fit -> gboolean: ENnu-fit
        N_\nu:fit
      Omegab-fit -> gboolean: Omegab-fit
        \Omega_{b0}:fit
      massnu-fit -> GVariant: massnu-fit
        m_\nu:fit
      Tnu-fit -> GVariant: Tnu-fit
        T_{\nu0}:fit
      munu-fit -> GVariant: munu-fit
        \mu_{\nu}:fit
      gnu-fit -> GVariant: gnu-fit
        g_{\nu}:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        w0: float
        w0_fit: bool
        w1: float
        w1_fit: bool
        CCL_comp: bool
        ENnu: float
        ENnu_fit: bool
        H0: float
        H0_fit: bool
        Omegab: float
        Omegab_fit: bool
        Omegac: float
        Omegac_fit: bool
        Omegax: float
        Omegax_fit: bool
        Tgamma0: float
        Tgamma0_fit: bool
        Tnu: NumCosmoMath.Vector
        Tnu_fit: GLib.Variant
        Tnu_length: int
        Yp: float
        Yp_fit: bool
        gnu: NumCosmoMath.Vector
        gnu_fit: GLib.Variant
        gnu_length: int
        massnu: NumCosmoMath.Vector
        massnu_fit: GLib.Variant
        massnu_length: int
        munu: NumCosmoMath.Vector
        munu_fit: GLib.Variant
        munu_length: int
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: HICosmoDE = ...
    def __init__(
        self,
        w0: float = ...,
        w0_fit: bool = ...,
        w1: float = ...,
        w1_fit: bool = ...,
        CCL_comp: bool = ...,
        ENnu: float = ...,
        ENnu_fit: bool = ...,
        H0: float = ...,
        H0_fit: bool = ...,
        Omegab: float = ...,
        Omegab_fit: bool = ...,
        Omegac: float = ...,
        Omegac_fit: bool = ...,
        Omegax: float = ...,
        Omegax_fit: bool = ...,
        Tgamma0: float = ...,
        Tgamma0_fit: bool = ...,
        Tnu: NumCosmoMath.Vector = ...,
        Tnu_fit: GLib.Variant = ...,
        Tnu_length: int = ...,
        Yp: float = ...,
        Yp_fit: bool = ...,
        gnu: NumCosmoMath.Vector = ...,
        gnu_fit: GLib.Variant = ...,
        gnu_length: int = ...,
        massnu: NumCosmoMath.Vector = ...,
        massnu_fit: GLib.Variant = ...,
        massnu_length: int = ...,
        munu: NumCosmoMath.Vector = ...,
        munu_fit: GLib.Variant = ...,
        munu_length: int = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @classmethod
    def new(cls) -> HICosmoDEJbp: ...

class HICosmoDEJbpClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoDEJbpClass()
    """

    parent_class: HICosmoDEClass = ...

class HICosmoDEPrivate(GObject.GPointer): ...

class HICosmoDEReparamCMB(NumCosmoMath.Reparam):
    r"""
    :Constructors:

    ::

        HICosmoDEReparamCMB(**properties)
        new(length:int) -> NumCosmo.HICosmoDEReparamCMB

    Object NcHICosmoDEReparamCMB

    Properties from NcmReparam:
      length -> guint: length
        System's length
      params-desc -> NcmObjDictInt: params-desc
        News parameter descriptions
      compat-type -> gchararray: compat-type
        Compatible type

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        compat_type: str
        length: int
        params_desc: NumCosmoMath.ObjDictInt

    props: Props = ...
    parent_instance: NumCosmoMath.Reparam = ...
    def __init__(
        self,
        compat_type: str = ...,
        length: int = ...,
        params_desc: NumCosmoMath.ObjDictInt = ...,
    ) -> None: ...
    @classmethod
    def new(cls, length: int) -> HICosmoDEReparamCMB: ...

class HICosmoDEReparamCMBClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoDEReparamCMBClass()
    """

    parent_class: NumCosmoMath.ReparamClass = ...

class HICosmoDEReparamOk(NumCosmoMath.Reparam):
    r"""
    :Constructors:

    ::

        HICosmoDEReparamOk(**properties)
        new(length:int) -> NumCosmo.HICosmoDEReparamOk

    Object NcHICosmoDEReparamOk

    Properties from NcmReparam:
      length -> guint: length
        System's length
      params-desc -> NcmObjDictInt: params-desc
        News parameter descriptions
      compat-type -> gchararray: compat-type
        Compatible type

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        compat_type: str
        length: int
        params_desc: NumCosmoMath.ObjDictInt

    props: Props = ...
    parent_instance: NumCosmoMath.Reparam = ...
    def __init__(
        self,
        compat_type: str = ...,
        length: int = ...,
        params_desc: NumCosmoMath.ObjDictInt = ...,
    ) -> None: ...
    @classmethod
    def new(cls, length: int) -> HICosmoDEReparamOk: ...

class HICosmoDEReparamOkClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoDEReparamOkClass()
    """

    parent_class: NumCosmoMath.ReparamClass = ...

class HICosmoDEWSpline(HICosmoDE):
    r"""
    :Constructors:

    ::

        HICosmoDEWSpline(**properties)
        new(nknots:int, z_f:float) -> NumCosmo.HICosmoDEWSpline

    Object NcHICosmoDEWSpline

    Properties from NcHICosmoDEWSpline:
      z1 -> gdouble: z1
        second redshift knot
      zf -> gdouble: zf
        final redshift
      w -> NcmVector: w
        w
      w-length -> guint: w-length
        w:length
      w-fit -> GVariant: w-fit
        w:fit

    Properties from NcHICosmoDE:
      CCL-comp -> gboolean: CCL-comp
        Whether to use CCL compatible mode
      H0 -> gdouble: H0
        H_0
      Omegac -> gdouble: Omegac
        \Omega_{c0}
      Omegax -> gdouble: Omegax
        \Omega_{x0}
      Tgamma0 -> gdouble: Tgamma0
        T_{\gamma0}
      Yp -> gdouble: Yp
        Y_p
      ENnu -> gdouble: ENnu
        N_\nu
      Omegab -> gdouble: Omegab
        \Omega_{b0}
      massnu -> NcmVector: massnu
        m_\nu
      Tnu -> NcmVector: Tnu
        T_{\nu0}
      munu -> NcmVector: munu
        \mu_{\nu}
      gnu -> NcmVector: gnu
        g_{\nu}
      massnu-length -> guint: massnu-length
        m_\nu:length
      Tnu-length -> guint: Tnu-length
        T_{\nu0}:length
      munu-length -> guint: munu-length
        \mu_{\nu}:length
      gnu-length -> guint: gnu-length
        g_{\nu}:length
      H0-fit -> gboolean: H0-fit
        H_0:fit
      Omegac-fit -> gboolean: Omegac-fit
        \Omega_{c0}:fit
      Omegax-fit -> gboolean: Omegax-fit
        \Omega_{x0}:fit
      Tgamma0-fit -> gboolean: Tgamma0-fit
        T_{\gamma0}:fit
      Yp-fit -> gboolean: Yp-fit
        Y_p:fit
      ENnu-fit -> gboolean: ENnu-fit
        N_\nu:fit
      Omegab-fit -> gboolean: Omegab-fit
        \Omega_{b0}:fit
      massnu-fit -> GVariant: massnu-fit
        m_\nu:fit
      Tnu-fit -> GVariant: Tnu-fit
        T_{\nu0}:fit
      munu-fit -> GVariant: munu-fit
        \mu_{\nu}:fit
      gnu-fit -> GVariant: gnu-fit
        g_{\nu}:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        w: NumCosmoMath.Vector
        w_fit: GLib.Variant
        w_length: int
        z1: float
        zf: float
        CCL_comp: bool
        ENnu: float
        ENnu_fit: bool
        H0: float
        H0_fit: bool
        Omegab: float
        Omegab_fit: bool
        Omegac: float
        Omegac_fit: bool
        Omegax: float
        Omegax_fit: bool
        Tgamma0: float
        Tgamma0_fit: bool
        Tnu: NumCosmoMath.Vector
        Tnu_fit: GLib.Variant
        Tnu_length: int
        Yp: float
        Yp_fit: bool
        gnu: NumCosmoMath.Vector
        gnu_fit: GLib.Variant
        gnu_length: int
        massnu: NumCosmoMath.Vector
        massnu_fit: GLib.Variant
        massnu_length: int
        munu: NumCosmoMath.Vector
        munu_fit: GLib.Variant
        munu_length: int
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: HICosmoDE = ...
    priv: HICosmoDEWSplinePrivate = ...
    def __init__(
        self,
        w: NumCosmoMath.Vector = ...,
        w_fit: GLib.Variant = ...,
        w_length: int = ...,
        z1: float = ...,
        zf: float = ...,
        CCL_comp: bool = ...,
        ENnu: float = ...,
        ENnu_fit: bool = ...,
        H0: float = ...,
        H0_fit: bool = ...,
        Omegab: float = ...,
        Omegab_fit: bool = ...,
        Omegac: float = ...,
        Omegac_fit: bool = ...,
        Omegax: float = ...,
        Omegax_fit: bool = ...,
        Tgamma0: float = ...,
        Tgamma0_fit: bool = ...,
        Tnu: NumCosmoMath.Vector = ...,
        Tnu_fit: GLib.Variant = ...,
        Tnu_length: int = ...,
        Yp: float = ...,
        Yp_fit: bool = ...,
        gnu: NumCosmoMath.Vector = ...,
        gnu_fit: GLib.Variant = ...,
        gnu_length: int = ...,
        massnu: NumCosmoMath.Vector = ...,
        massnu_fit: GLib.Variant = ...,
        massnu_length: int = ...,
        munu: NumCosmoMath.Vector = ...,
        munu_fit: GLib.Variant = ...,
        munu_length: int = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    def get_alpha(self) -> NumCosmoMath.Vector: ...
    @classmethod
    def new(cls, nknots: int, z_f: float) -> HICosmoDEWSpline: ...

class HICosmoDEWSplineClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoDEWSplineClass()
    """

    parent_class: HICosmoDEClass = ...

class HICosmoDEWSplinePrivate(GObject.GPointer): ...

class HICosmoDEXcdm(HICosmoDE):
    r"""
    :Constructors:

    ::

        HICosmoDEXcdm(**properties)
        new() -> NumCosmo.HICosmoDEXcdm

    Object NcHICosmoDEXcdm

    Properties from NcHICosmoDEXcdm:
      w -> gdouble: w
        w
      w-fit -> gboolean: w-fit
        w:fit

    Properties from NcHICosmoDE:
      CCL-comp -> gboolean: CCL-comp
        Whether to use CCL compatible mode
      H0 -> gdouble: H0
        H_0
      Omegac -> gdouble: Omegac
        \Omega_{c0}
      Omegax -> gdouble: Omegax
        \Omega_{x0}
      Tgamma0 -> gdouble: Tgamma0
        T_{\gamma0}
      Yp -> gdouble: Yp
        Y_p
      ENnu -> gdouble: ENnu
        N_\nu
      Omegab -> gdouble: Omegab
        \Omega_{b0}
      massnu -> NcmVector: massnu
        m_\nu
      Tnu -> NcmVector: Tnu
        T_{\nu0}
      munu -> NcmVector: munu
        \mu_{\nu}
      gnu -> NcmVector: gnu
        g_{\nu}
      massnu-length -> guint: massnu-length
        m_\nu:length
      Tnu-length -> guint: Tnu-length
        T_{\nu0}:length
      munu-length -> guint: munu-length
        \mu_{\nu}:length
      gnu-length -> guint: gnu-length
        g_{\nu}:length
      H0-fit -> gboolean: H0-fit
        H_0:fit
      Omegac-fit -> gboolean: Omegac-fit
        \Omega_{c0}:fit
      Omegax-fit -> gboolean: Omegax-fit
        \Omega_{x0}:fit
      Tgamma0-fit -> gboolean: Tgamma0-fit
        T_{\gamma0}:fit
      Yp-fit -> gboolean: Yp-fit
        Y_p:fit
      ENnu-fit -> gboolean: ENnu-fit
        N_\nu:fit
      Omegab-fit -> gboolean: Omegab-fit
        \Omega_{b0}:fit
      massnu-fit -> GVariant: massnu-fit
        m_\nu:fit
      Tnu-fit -> GVariant: Tnu-fit
        T_{\nu0}:fit
      munu-fit -> GVariant: munu-fit
        \mu_{\nu}:fit
      gnu-fit -> GVariant: gnu-fit
        g_{\nu}:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        w: float
        w_fit: bool
        CCL_comp: bool
        ENnu: float
        ENnu_fit: bool
        H0: float
        H0_fit: bool
        Omegab: float
        Omegab_fit: bool
        Omegac: float
        Omegac_fit: bool
        Omegax: float
        Omegax_fit: bool
        Tgamma0: float
        Tgamma0_fit: bool
        Tnu: NumCosmoMath.Vector
        Tnu_fit: GLib.Variant
        Tnu_length: int
        Yp: float
        Yp_fit: bool
        gnu: NumCosmoMath.Vector
        gnu_fit: GLib.Variant
        gnu_length: int
        massnu: NumCosmoMath.Vector
        massnu_fit: GLib.Variant
        massnu_length: int
        munu: NumCosmoMath.Vector
        munu_fit: GLib.Variant
        munu_length: int
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: HICosmoDE = ...
    def __init__(
        self,
        w: float = ...,
        w_fit: bool = ...,
        CCL_comp: bool = ...,
        ENnu: float = ...,
        ENnu_fit: bool = ...,
        H0: float = ...,
        H0_fit: bool = ...,
        Omegab: float = ...,
        Omegab_fit: bool = ...,
        Omegac: float = ...,
        Omegac_fit: bool = ...,
        Omegax: float = ...,
        Omegax_fit: bool = ...,
        Tgamma0: float = ...,
        Tgamma0_fit: bool = ...,
        Tnu: NumCosmoMath.Vector = ...,
        Tnu_fit: GLib.Variant = ...,
        Tnu_length: int = ...,
        Yp: float = ...,
        Yp_fit: bool = ...,
        gnu: NumCosmoMath.Vector = ...,
        gnu_fit: GLib.Variant = ...,
        gnu_length: int = ...,
        massnu: NumCosmoMath.Vector = ...,
        massnu_fit: GLib.Variant = ...,
        massnu_length: int = ...,
        munu: NumCosmoMath.Vector = ...,
        munu_fit: GLib.Variant = ...,
        munu_length: int = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @classmethod
    def new(cls) -> HICosmoDEXcdm: ...

class HICosmoDEXcdmClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoDEXcdmClass()
    """

    parent_class: HICosmoDEClass = ...

class HICosmoFunc(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoFunc()
    """

    name: str = ...
    desc: str = ...
    f: typing.Callable[[HICosmo], float] = ...
    impl: HICosmoImpl = ...

class HICosmoFuncZ(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoFuncZ()
    """

    name: str = ...
    desc: str = ...
    f: typing.Callable[[HICosmo, float], float] = ...
    impl: HICosmoImpl = ...

class HICosmoGCG(HICosmo):
    r"""
    :Constructors:

    ::

        HICosmoGCG(**properties)

    Object NcHICosmoGCG

    Properties from NcHICosmoGCG:
      H0 -> gdouble: H0
        H_0
      Omegac -> gdouble: Omegac
        \Omega_{c0}
      Omegax -> gdouble: Omegax
        \Omega_{x0}
      Tgamma0 -> gdouble: Tgamma0
        T_{\gamma0}
      Yp -> gdouble: Yp
        Y_p
      ENnu -> gdouble: ENnu
        N_\nu
      Omegab -> gdouble: Omegab
        \Omega_{b0}
      gamma -> gdouble: gamma
        \gamma
      massnu -> NcmVector: massnu
        m_\nu
      Tnu -> NcmVector: Tnu
        T_{\nu0}
      munu -> NcmVector: munu
        \mu_{\nu}
      gnu -> NcmVector: gnu
        g_{\nu}
      massnu-length -> guint: massnu-length
        m_\nu:length
      Tnu-length -> guint: Tnu-length
        T_{\nu0}:length
      munu-length -> guint: munu-length
        \mu_{\nu}:length
      gnu-length -> guint: gnu-length
        g_{\nu}:length
      H0-fit -> gboolean: H0-fit
        H_0:fit
      Omegac-fit -> gboolean: Omegac-fit
        \Omega_{c0}:fit
      Omegax-fit -> gboolean: Omegax-fit
        \Omega_{x0}:fit
      Tgamma0-fit -> gboolean: Tgamma0-fit
        T_{\gamma0}:fit
      Yp-fit -> gboolean: Yp-fit
        Y_p:fit
      ENnu-fit -> gboolean: ENnu-fit
        N_\nu:fit
      Omegab-fit -> gboolean: Omegab-fit
        \Omega_{b0}:fit
      gamma-fit -> gboolean: gamma-fit
        \gamma:fit
      massnu-fit -> GVariant: massnu-fit
        m_\nu:fit
      Tnu-fit -> GVariant: Tnu-fit
        T_{\nu0}:fit
      munu-fit -> GVariant: munu-fit
        \mu_{\nu}:fit
      gnu-fit -> GVariant: gnu-fit
        g_{\nu}:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        ENnu: float
        ENnu_fit: bool
        H0: float
        H0_fit: bool
        Omegab: float
        Omegab_fit: bool
        Omegac: float
        Omegac_fit: bool
        Omegax: float
        Omegax_fit: bool
        Tgamma0: float
        Tgamma0_fit: bool
        Tnu: NumCosmoMath.Vector
        Tnu_fit: GLib.Variant
        Tnu_length: int
        Yp: float
        Yp_fit: bool
        gamma: float
        gamma_fit: bool
        gnu: NumCosmoMath.Vector
        gnu_fit: GLib.Variant
        gnu_length: int
        massnu: NumCosmoMath.Vector
        massnu_fit: GLib.Variant
        massnu_length: int
        munu: NumCosmoMath.Vector
        munu_fit: GLib.Variant
        munu_length: int
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: HICosmo = ...
    priv: HICosmoGCGPrivate = ...
    def __init__(
        self,
        ENnu: float = ...,
        ENnu_fit: bool = ...,
        H0: float = ...,
        H0_fit: bool = ...,
        Omegab: float = ...,
        Omegab_fit: bool = ...,
        Omegac: float = ...,
        Omegac_fit: bool = ...,
        Omegax: float = ...,
        Omegax_fit: bool = ...,
        Tgamma0: float = ...,
        Tgamma0_fit: bool = ...,
        Tnu: NumCosmoMath.Vector = ...,
        Tnu_fit: GLib.Variant = ...,
        Tnu_length: int = ...,
        Yp: float = ...,
        Yp_fit: bool = ...,
        gamma: float = ...,
        gamma_fit: bool = ...,
        gnu: NumCosmoMath.Vector = ...,
        gnu_fit: GLib.Variant = ...,
        gnu_length: int = ...,
        massnu: NumCosmoMath.Vector = ...,
        massnu_fit: GLib.Variant = ...,
        massnu_length: int = ...,
        munu: NumCosmoMath.Vector = ...,
        munu_fit: GLib.Variant = ...,
        munu_length: int = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    def cmb_params(self) -> None: ...
    def omega_x2omega_k(self) -> None: ...

class HICosmoGCGClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoGCGClass()
    """

    parent_class: HICosmoClass = ...

class HICosmoGCGPrivate(GObject.GPointer): ...

class HICosmoGCGReparamCMB(NumCosmoMath.Reparam):
    r"""
    :Constructors:

    ::

        HICosmoGCGReparamCMB(**properties)
        new(length:int) -> NumCosmo.HICosmoGCGReparamCMB

    Object NcHICosmoGCGReparamCMB

    Properties from NcmReparam:
      length -> guint: length
        System's length
      params-desc -> NcmObjDictInt: params-desc
        News parameter descriptions
      compat-type -> gchararray: compat-type
        Compatible type

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        compat_type: str
        length: int
        params_desc: NumCosmoMath.ObjDictInt

    props: Props = ...
    parent_instance: NumCosmoMath.Reparam = ...
    def __init__(
        self,
        compat_type: str = ...,
        length: int = ...,
        params_desc: NumCosmoMath.ObjDictInt = ...,
    ) -> None: ...
    @classmethod
    def new(cls, length: int) -> HICosmoGCGReparamCMB: ...

class HICosmoGCGReparamCMBClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoGCGReparamCMBClass()
    """

    parent_class: NumCosmoMath.ReparamClass = ...

class HICosmoGCGReparamOk(NumCosmoMath.Reparam):
    r"""
    :Constructors:

    ::

        HICosmoGCGReparamOk(**properties)
        new(length:int) -> NumCosmo.HICosmoGCGReparamOk

    Object NcHICosmoGCGReparamOk

    Properties from NcmReparam:
      length -> guint: length
        System's length
      params-desc -> NcmObjDictInt: params-desc
        News parameter descriptions
      compat-type -> gchararray: compat-type
        Compatible type

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        compat_type: str
        length: int
        params_desc: NumCosmoMath.ObjDictInt

    props: Props = ...
    parent_instance: NumCosmoMath.Reparam = ...
    def __init__(
        self,
        compat_type: str = ...,
        length: int = ...,
        params_desc: NumCosmoMath.ObjDictInt = ...,
    ) -> None: ...
    @classmethod
    def new(cls, length: int) -> HICosmoGCGReparamOk: ...

class HICosmoGCGReparamOkClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoGCGReparamOkClass()
    """

    parent_class: NumCosmoMath.ReparamClass = ...

class HICosmoIDEM2(HICosmo):
    r"""
    :Constructors:

    ::

        HICosmoIDEM2(**properties)

    Object NcHICosmoIDEM2

    Properties from NcHICosmoIDEM2:
      H0 -> gdouble: H0
        H_0
      Omegac -> gdouble: Omegac
        \Omega_{c0}
      Omegax -> gdouble: Omegax
        \Omega_{x0}
      Tgamma0 -> gdouble: Tgamma0
        T_{\gamma0}
      Yp -> gdouble: Yp
        Y_p
      ENnu -> gdouble: ENnu
        N_\nu
      Omegab -> gdouble: Omegab
        \Omega_{b0}
      gamma -> gdouble: gamma
        \gamma
      massnu -> NcmVector: massnu
        m_\nu
      Tnu -> NcmVector: Tnu
        T_{\nu0}
      munu -> NcmVector: munu
        \mu_{\nu}
      gnu -> NcmVector: gnu
        g_{\nu}
      massnu-length -> guint: massnu-length
        m_\nu:length
      Tnu-length -> guint: Tnu-length
        T_{\nu0}:length
      munu-length -> guint: munu-length
        \mu_{\nu}:length
      gnu-length -> guint: gnu-length
        g_{\nu}:length
      H0-fit -> gboolean: H0-fit
        H_0:fit
      Omegac-fit -> gboolean: Omegac-fit
        \Omega_{c0}:fit
      Omegax-fit -> gboolean: Omegax-fit
        \Omega_{x0}:fit
      Tgamma0-fit -> gboolean: Tgamma0-fit
        T_{\gamma0}:fit
      Yp-fit -> gboolean: Yp-fit
        Y_p:fit
      ENnu-fit -> gboolean: ENnu-fit
        N_\nu:fit
      Omegab-fit -> gboolean: Omegab-fit
        \Omega_{b0}:fit
      gamma-fit -> gboolean: gamma-fit
        \gamma:fit
      massnu-fit -> GVariant: massnu-fit
        m_\nu:fit
      Tnu-fit -> GVariant: Tnu-fit
        T_{\nu0}:fit
      munu-fit -> GVariant: munu-fit
        \mu_{\nu}:fit
      gnu-fit -> GVariant: gnu-fit
        g_{\nu}:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        ENnu: float
        ENnu_fit: bool
        H0: float
        H0_fit: bool
        Omegab: float
        Omegab_fit: bool
        Omegac: float
        Omegac_fit: bool
        Omegax: float
        Omegax_fit: bool
        Tgamma0: float
        Tgamma0_fit: bool
        Tnu: NumCosmoMath.Vector
        Tnu_fit: GLib.Variant
        Tnu_length: int
        Yp: float
        Yp_fit: bool
        gamma: float
        gamma_fit: bool
        gnu: NumCosmoMath.Vector
        gnu_fit: GLib.Variant
        gnu_length: int
        massnu: NumCosmoMath.Vector
        massnu_fit: GLib.Variant
        massnu_length: int
        munu: NumCosmoMath.Vector
        munu_fit: GLib.Variant
        munu_length: int
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: HICosmo = ...
    priv: HICosmoIDEM2Private = ...
    def __init__(
        self,
        ENnu: float = ...,
        ENnu_fit: bool = ...,
        H0: float = ...,
        H0_fit: bool = ...,
        Omegab: float = ...,
        Omegab_fit: bool = ...,
        Omegac: float = ...,
        Omegac_fit: bool = ...,
        Omegax: float = ...,
        Omegax_fit: bool = ...,
        Tgamma0: float = ...,
        Tgamma0_fit: bool = ...,
        Tnu: NumCosmoMath.Vector = ...,
        Tnu_fit: GLib.Variant = ...,
        Tnu_length: int = ...,
        Yp: float = ...,
        Yp_fit: bool = ...,
        gamma: float = ...,
        gamma_fit: bool = ...,
        gnu: NumCosmoMath.Vector = ...,
        gnu_fit: GLib.Variant = ...,
        gnu_length: int = ...,
        massnu: NumCosmoMath.Vector = ...,
        massnu_fit: GLib.Variant = ...,
        massnu_length: int = ...,
        munu: NumCosmoMath.Vector = ...,
        munu_fit: GLib.Variant = ...,
        munu_length: int = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    def cmb_params(self) -> None: ...
    def omega_x2omega_k(self) -> None: ...

class HICosmoIDEM2Class(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoIDEM2Class()
    """

    parent_class: HICosmoClass = ...

class HICosmoIDEM2Private(GObject.GPointer): ...

class HICosmoIDEM2ReparamCMB(NumCosmoMath.Reparam):
    r"""
    :Constructors:

    ::

        HICosmoIDEM2ReparamCMB(**properties)
        new(length:int) -> NumCosmo.HICosmoIDEM2ReparamCMB

    Object NcHICosmoIDEM2ReparamCMB

    Properties from NcmReparam:
      length -> guint: length
        System's length
      params-desc -> NcmObjDictInt: params-desc
        News parameter descriptions
      compat-type -> gchararray: compat-type
        Compatible type

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        compat_type: str
        length: int
        params_desc: NumCosmoMath.ObjDictInt

    props: Props = ...
    parent_instance: NumCosmoMath.Reparam = ...
    def __init__(
        self,
        compat_type: str = ...,
        length: int = ...,
        params_desc: NumCosmoMath.ObjDictInt = ...,
    ) -> None: ...
    @classmethod
    def new(cls, length: int) -> HICosmoIDEM2ReparamCMB: ...

class HICosmoIDEM2ReparamCMBClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoIDEM2ReparamCMBClass()
    """

    parent_class: NumCosmoMath.ReparamClass = ...

class HICosmoIDEM2ReparamOk(NumCosmoMath.Reparam):
    r"""
    :Constructors:

    ::

        HICosmoIDEM2ReparamOk(**properties)
        new(length:int) -> NumCosmo.HICosmoIDEM2ReparamOk

    Object NcHICosmoIDEM2ReparamOk

    Properties from NcmReparam:
      length -> guint: length
        System's length
      params-desc -> NcmObjDictInt: params-desc
        News parameter descriptions
      compat-type -> gchararray: compat-type
        Compatible type

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        compat_type: str
        length: int
        params_desc: NumCosmoMath.ObjDictInt

    props: Props = ...
    parent_instance: NumCosmoMath.Reparam = ...
    def __init__(
        self,
        compat_type: str = ...,
        length: int = ...,
        params_desc: NumCosmoMath.ObjDictInt = ...,
    ) -> None: ...
    @classmethod
    def new(cls, length: int) -> HICosmoIDEM2ReparamOk: ...

class HICosmoIDEM2ReparamOkClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoIDEM2ReparamOkClass()
    """

    parent_class: NumCosmoMath.ReparamClass = ...

class HICosmoLCDM(HICosmo):
    r"""
    :Constructors:

    ::

        HICosmoLCDM(**properties)
        new() -> NumCosmo.HICosmoLCDM

    Object NcHICosmoLCDM

    Properties from NcHICosmoLCDM:
      H0 -> gdouble: H0
        H_0
      Omegac -> gdouble: Omegac
        \Omega_{c0}
      Omegax -> gdouble: Omegax
        \Omega_{x0}
      Tgamma0 -> gdouble: Tgamma0
        T_{\gamma0}
      Yp -> gdouble: Yp
        Y_p
      ENnu -> gdouble: ENnu
        N_\nu
      Omegab -> gdouble: Omegab
        \Omega_{b0}
      H0-fit -> gboolean: H0-fit
        H_0:fit
      Omegac-fit -> gboolean: Omegac-fit
        \Omega_{c0}:fit
      Omegax-fit -> gboolean: Omegax-fit
        \Omega_{x0}:fit
      Tgamma0-fit -> gboolean: Tgamma0-fit
        T_{\gamma0}:fit
      Yp-fit -> gboolean: Yp-fit
        Y_p:fit
      ENnu-fit -> gboolean: ENnu-fit
        N_\nu:fit
      Omegab-fit -> gboolean: Omegab-fit
        \Omega_{b0}:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        ENnu: float
        ENnu_fit: bool
        H0: float
        H0_fit: bool
        Omegab: float
        Omegab_fit: bool
        Omegac: float
        Omegac_fit: bool
        Omegax: float
        Omegax_fit: bool
        Tgamma0: float
        Tgamma0_fit: bool
        Yp: float
        Yp_fit: bool
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: HICosmo = ...
    def __init__(
        self,
        ENnu: float = ...,
        ENnu_fit: bool = ...,
        H0: float = ...,
        H0_fit: bool = ...,
        Omegab: float = ...,
        Omegab_fit: bool = ...,
        Omegac: float = ...,
        Omegac_fit: bool = ...,
        Omegax: float = ...,
        Omegax_fit: bool = ...,
        Tgamma0: float = ...,
        Tgamma0_fit: bool = ...,
        Yp: float = ...,
        Yp_fit: bool = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @classmethod
    def new(cls) -> HICosmoLCDM: ...

class HICosmoLCDMClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoLCDMClass()
    """

    parent_class: HICosmoClass = ...

class HICosmoQConst(HICosmo):
    r"""
    :Constructors:

    ::

        HICosmoQConst(**properties)
        new() -> NumCosmo.HICosmoQConst

    Object NcHICosmoQConst

    Properties from NcHICosmoQConst:
      H0 -> gdouble: H0
        H_0
      Omegat -> gdouble: Omegat
        \Omega_{t0}
      Dc -> gdouble: Dc
        D_c
      E -> gdouble: E
        E
      q -> gdouble: q
        q
      zs -> gdouble: zs
        z_\star
      H0-fit -> gboolean: H0-fit
        H_0:fit
      Omegat-fit -> gboolean: Omegat-fit
        \Omega_{t0}:fit
      Dc-fit -> gboolean: Dc-fit
        D_c:fit
      E-fit -> gboolean: E-fit
        E:fit
      q-fit -> gboolean: q-fit
        q:fit
      zs-fit -> gboolean: zs-fit
        z_\star:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        Dc: float
        Dc_fit: bool
        E: float
        E_fit: bool
        H0: float
        H0_fit: bool
        Omegat: float
        Omegat_fit: bool
        q: float
        q_fit: bool
        zs: float
        zs_fit: bool
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: HICosmo = ...
    def __init__(
        self,
        Dc: float = ...,
        Dc_fit: bool = ...,
        E: float = ...,
        E_fit: bool = ...,
        H0: float = ...,
        H0_fit: bool = ...,
        Omegat: float = ...,
        Omegat_fit: bool = ...,
        q: float = ...,
        q_fit: bool = ...,
        zs: float = ...,
        zs_fit: bool = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @classmethod
    def new(cls) -> HICosmoQConst: ...

class HICosmoQConstClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoQConstClass()
    """

    parent_class: HICosmoClass = ...

class HICosmoQGRW(HICosmo, HIPertIGW, HIPertITwoFluids):
    r"""
    :Constructors:

    ::

        HICosmoQGRW(**properties)
        new() -> NumCosmo.HICosmoQGRW

    Object NcHICosmoQGRW

    Properties from NcHICosmoQGRW:
      H0 -> gdouble: H0
        H_0
      Omegar -> gdouble: Omegar
        \Omega_{r0}
      Omegaw -> gdouble: Omegaw
        \Omega_{w0}
      w -> gdouble: w
        w
      xb -> gdouble: xb
        x_b
      H0-fit -> gboolean: H0-fit
        H_0:fit
      Omegar-fit -> gboolean: Omegar-fit
        \Omega_{r0}:fit
      Omegaw-fit -> gboolean: Omegaw-fit
        \Omega_{w0}:fit
      w-fit -> gboolean: w-fit
        w:fit
      xb-fit -> gboolean: xb-fit
        x_b:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        H0: float
        H0_fit: bool
        Omegar: float
        Omegar_fit: bool
        Omegaw: float
        Omegaw_fit: bool
        w: float
        w_fit: bool
        xb: float
        xb_fit: bool
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: HICosmo = ...
    eom_two_fluids: HIPertITwoFluidsEOM = ...
    tv_two_fluids: HIPertITwoFluidsTV = ...
    def __init__(
        self,
        H0: float = ...,
        H0_fit: bool = ...,
        Omegar: float = ...,
        Omegar_fit: bool = ...,
        Omegaw: float = ...,
        Omegaw_fit: bool = ...,
        w: float = ...,
        w_fit: bool = ...,
        xb: float = ...,
        xb_fit: bool = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @classmethod
    def new(cls) -> HICosmoQGRW: ...

class HICosmoQGRWClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoQGRWClass()
    """

    parent_class: HICosmoClass = ...

class HICosmoQGW(HICosmo, HIPertIAdiab):
    r"""
    :Constructors:

    ::

        HICosmoQGW(**properties)
        new() -> NumCosmo.HICosmoQGW

    Object NcHICosmoQGW

    Properties from NcHICosmoQGW:
      H0 -> gdouble: H0
        H_0
      Omegaw -> gdouble: Omegaw
        \Omega_{w0}
      w -> gdouble: w
        w
      xb -> gdouble: xb
        x_b
      H0-fit -> gboolean: H0-fit
        H_0:fit
      Omegaw-fit -> gboolean: Omegaw-fit
        \Omega_{w0}:fit
      w-fit -> gboolean: w-fit
        w:fit
      xb-fit -> gboolean: xb-fit
        x_b:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        H0: float
        H0_fit: bool
        Omegaw: float
        Omegaw_fit: bool
        w: float
        w_fit: bool
        xb: float
        xb_fit: bool
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: HICosmo = ...
    priv: HICosmoQGWPrivate = ...
    def __init__(
        self,
        H0: float = ...,
        H0_fit: bool = ...,
        Omegaw: float = ...,
        Omegaw_fit: bool = ...,
        w: float = ...,
        w_fit: bool = ...,
        xb: float = ...,
        xb_fit: bool = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @classmethod
    def new(cls) -> HICosmoQGW: ...

class HICosmoQGWClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoQGWClass()
    """

    parent_class: HICosmoClass = ...

class HICosmoQGWPrivate(GObject.GPointer): ...

class HICosmoQLinear(HICosmo):
    r"""
    :Constructors:

    ::

        HICosmoQLinear(**properties)
        new() -> NumCosmo.HICosmoQLinear

    Object NcHICosmoQLinear

    Properties from NcHICosmoQLinear:
      H0 -> gdouble: H0
        H_0
      Omegat -> gdouble: Omegat
        \Omega_{t0}
      Dc -> gdouble: Dc
        D_c
      E -> gdouble: E
        E
      q -> gdouble: q
        q
      qp -> gdouble: qp
        q^\prime
      zs -> gdouble: zs
        z_\star
      H0-fit -> gboolean: H0-fit
        H_0:fit
      Omegat-fit -> gboolean: Omegat-fit
        \Omega_{t0}:fit
      Dc-fit -> gboolean: Dc-fit
        D_c:fit
      E-fit -> gboolean: E-fit
        E:fit
      q-fit -> gboolean: q-fit
        q:fit
      qp-fit -> gboolean: qp-fit
        q^\prime:fit
      zs-fit -> gboolean: zs-fit
        z_\star:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        Dc: float
        Dc_fit: bool
        E: float
        E_fit: bool
        H0: float
        H0_fit: bool
        Omegat: float
        Omegat_fit: bool
        q: float
        q_fit: bool
        qp: float
        qp_fit: bool
        zs: float
        zs_fit: bool
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: HICosmo = ...
    def __init__(
        self,
        Dc: float = ...,
        Dc_fit: bool = ...,
        E: float = ...,
        E_fit: bool = ...,
        H0: float = ...,
        H0_fit: bool = ...,
        Omegat: float = ...,
        Omegat_fit: bool = ...,
        q: float = ...,
        q_fit: bool = ...,
        qp: float = ...,
        qp_fit: bool = ...,
        zs: float = ...,
        zs_fit: bool = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def dE(z2: float, z1: float, q: float, qp: float) -> float: ...
    @classmethod
    def new(cls) -> HICosmoQLinear: ...

class HICosmoQLinearClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoQLinearClass()
    """

    parent_class: HICosmoClass = ...

class HICosmoQRBF(HICosmo):
    r"""
    :Constructors:

    ::

        HICosmoQRBF(**properties)
        new(np:int, z_f:float) -> NumCosmo.HICosmoQRBF

    Object NcHICosmoQRBF

    Properties from NcHICosmoQRBF:
      zf -> gdouble: zf
        final redshift
      H0 -> gdouble: H0
        H_0
      Omegat -> gdouble: Omegat
        Omega_t0
      asdrag -> gdouble: asdrag
        A_s
      hr -> gdouble: hr
        h_r
      xi -> NcmVector: xi
        x_i
      ci -> NcmVector: ci
        c_i
      xi-length -> guint: xi-length
        x_i:length
      ci-length -> guint: ci-length
        c_i:length
      H0-fit -> gboolean: H0-fit
        H_0:fit
      Omegat-fit -> gboolean: Omegat-fit
        Omega_t0:fit
      asdrag-fit -> gboolean: asdrag-fit
        A_s:fit
      hr-fit -> gboolean: hr-fit
        h_r:fit
      xi-fit -> GVariant: xi-fit
        x_i:fit
      ci-fit -> GVariant: ci-fit
        c_i:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        H0: float
        H0_fit: bool
        Omegat: float
        Omegat_fit: bool
        asdrag: float
        asdrag_fit: bool
        ci: NumCosmoMath.Vector
        ci_fit: GLib.Variant
        ci_length: int
        hr: float
        hr_fit: bool
        xi: NumCosmoMath.Vector
        xi_fit: GLib.Variant
        xi_length: int
        zf: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: HICosmo = ...
    priv: HICosmoQRBFPrivate = ...
    def __init__(
        self,
        H0: float = ...,
        H0_fit: bool = ...,
        Omegat: float = ...,
        Omegat_fit: bool = ...,
        asdrag: float = ...,
        asdrag_fit: bool = ...,
        ci: NumCosmoMath.Vector = ...,
        ci_fit: GLib.Variant = ...,
        ci_length: int = ...,
        hr: float = ...,
        hr_fit: bool = ...,
        xi: NumCosmoMath.Vector = ...,
        xi_fit: GLib.Variant = ...,
        xi_length: int = ...,
        zf: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @classmethod
    def new(cls, np: int, z_f: float) -> HICosmoQRBF: ...
    def q_roughness(self) -> float: ...
    def set_z_f(self, z_f: float) -> None: ...

class HICosmoQRBFClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoQRBFClass()
    """

    parent_class: HICosmoClass = ...

class HICosmoQRBFPrivate(GObject.GPointer): ...

class HICosmoQRBFRprior(NumCosmoMath.Prior):
    r"""
    :Constructors:

    ::

        HICosmoQRBFRprior(**properties)
        new(lambda_:float) -> NumCosmo.HICosmoQRBFRprior

    Object NcHICosmoQRBFRprior

    Properties from NcHICosmoQRBFRprior:
      lambda -> gdouble: lambda
        \lambda

    Properties from NcmMSetFunc:
      nvariables -> guint: nvariables
        Number of variables
      dimension -> guint: dimension
        Function dimension
      eval-x -> NcmVector: eval-x
        Evaluation point x

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        dimension: int
        eval_x: NumCosmoMath.Vector
        nvariables: int

    props: Props = ...
    parent_instance: NumCosmoMath.Prior = ...
    priv: HICosmoQRBFRpriorPrivate = ...
    def __init__(
        self,
        dimension: int = ...,
        eval_x: NumCosmoMath.Vector = ...,
        nvariables: int = ...,
    ) -> None: ...
    @classmethod
    def new(cls, lambda_: float) -> HICosmoQRBFRprior: ...

class HICosmoQRBFRpriorClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoQRBFRpriorClass()
    """

    parent_class: NumCosmoMath.PriorClass = ...

class HICosmoQRBFRpriorPrivate(GObject.GPointer): ...

class HICosmoQSpline(HICosmo):
    r"""
    :Constructors:

    ::

        HICosmoQSpline(**properties)
        new(s:NumCosmoMath.Spline, np:int, z_f:float) -> NumCosmo.HICosmoQSpline

    Object NcHICosmoQSpline

    Properties from NcHICosmoQSpline:
      spline -> NcmSpline: spline
        Spline object
      zf -> gdouble: zf
        final redshift
      H0 -> gdouble: H0
        H_0
      Omegat -> gdouble: Omegat
        \Omega_{t0}
      asdrag -> gdouble: asdrag
        A_s
      qparam -> NcmVector: qparam
        q
      qparam-length -> guint: qparam-length
        q:length
      H0-fit -> gboolean: H0-fit
        H_0:fit
      Omegat-fit -> gboolean: Omegat-fit
        \Omega_{t0}:fit
      asdrag-fit -> gboolean: asdrag-fit
        A_s:fit
      qparam-fit -> GVariant: qparam-fit
        q:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        H0: float
        H0_fit: bool
        Omegat: float
        Omegat_fit: bool
        asdrag: float
        asdrag_fit: bool
        qparam: NumCosmoMath.Vector
        qparam_fit: GLib.Variant
        qparam_length: int
        spline: NumCosmoMath.Spline
        zf: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: HICosmo = ...
    nknots: int = ...
    size: int = ...
    z_f: float = ...
    q_z: NumCosmoMath.Spline = ...
    E2_z: NumCosmoMath.OdeSpline = ...
    def __init__(
        self,
        H0: float = ...,
        H0_fit: bool = ...,
        Omegat: float = ...,
        Omegat_fit: bool = ...,
        asdrag: float = ...,
        asdrag_fit: bool = ...,
        qparam: NumCosmoMath.Vector = ...,
        qparam_fit: GLib.Variant = ...,
        qparam_length: int = ...,
        spline: NumCosmoMath.Spline = ...,
        zf: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    def add_continuity_priors(
        self, lh: NumCosmoMath.Likelihood, sigma: float, abstol: float
    ) -> HICosmoQSplineContPrior: ...
    @classmethod
    def new(cls, s: NumCosmoMath.Spline, np: int, z_f: float) -> HICosmoQSpline: ...

class HICosmoQSplineClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoQSplineClass()
    """

    parent_class: HICosmoClass = ...

class HICosmoQSplineContPrior(NumCosmoMath.Model):
    r"""
    :Constructors:

    ::

        HICosmoQSplineContPrior(**properties)
        new(npriors:int) -> NumCosmo.HICosmoQSplineContPrior

    Object NcHICosmoQSplineContPrior

    Properties from NcHICosmoQSplineContPrior:
      abstol -> gdouble: abstol
        abstol
      lnsigma -> NcmVector: lnsigma
        lnsigma
      lnsigma-length -> guint: lnsigma-length
        lnsigma:length
      abstol-fit -> gboolean: abstol-fit
        abstol:fit
      lnsigma-fit -> GVariant: lnsigma-fit
        lnsigma:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        abstol: float
        abstol_fit: bool
        lnsigma: NumCosmoMath.Vector
        lnsigma_fit: GLib.Variant
        lnsigma_length: int
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: NumCosmoMath.Model = ...
    def __init__(
        self,
        abstol: float = ...,
        abstol_fit: bool = ...,
        lnsigma: NumCosmoMath.Vector = ...,
        lnsigma_fit: GLib.Variant = ...,
        lnsigma_length: int = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    def free(self) -> None: ...
    def get_abstol(self) -> float: ...
    def get_lnsigma(self, i: int) -> float: ...
    @staticmethod
    def id() -> int: ...
    @classmethod
    def new(cls, npriors: int) -> HICosmoQSplineContPrior: ...
    def ref(self) -> HICosmoQSplineContPrior: ...
    def set_abstol(self, abstol: float) -> None: ...
    def set_all_lnsigma(self, ln_sigma: float) -> None: ...
    def set_lnsigma(self, i: int, ln_sigma: float) -> None: ...

class HICosmoQSplineContPriorClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoQSplineContPriorClass()
    """

    parent_class: NumCosmoMath.ModelClass = ...

class HICosmoVexp(HICosmo, HIPertIAdiab, HIPertIEM, HIPertIGW):
    r"""
    :Constructors:

    ::

        HICosmoVexp(**properties)
        new() -> NumCosmo.HICosmoVexp

    Object NcHICosmoVexp

    Properties from NcHICosmoVexp:
      glue-de -> gboolean: glue-de
        Whether to glue to a DE phase
      set-xb-max -> gboolean: set-xb-max
        Whether to use max xb allowed by the matching
      em-coupling -> NcHICosmoVexpEMCoupling: em-coupling
        Electromagnetic coupling
      H0 -> gdouble: H0
        H_0
      Omegac -> gdouble: Omegac
        \Omega_{c0}
      OmegaL -> gdouble: OmegaL
        \Omega_{\Lambda0}
      sigmaphi -> gdouble: sigmaphi
        \sigma_{\phi}
      dphi -> gdouble: dphi
        d_\phi
      alphab -> gdouble: alphab
        \alpha_b
      xb -> gdouble: xb
        x_b
      alphaem -> gdouble: alphaem
        \alpha_{\mathrm{em}}
      betaem -> gdouble: betaem
        \beta_\mathrm{em}
      H0-fit -> gboolean: H0-fit
        H_0:fit
      Omegac-fit -> gboolean: Omegac-fit
        \Omega_{c0}:fit
      OmegaL-fit -> gboolean: OmegaL-fit
        \Omega_{\Lambda0}:fit
      sigmaphi-fit -> gboolean: sigmaphi-fit
        \sigma_{\phi}:fit
      dphi-fit -> gboolean: dphi-fit
        d_\phi:fit
      alphab-fit -> gboolean: alphab-fit
        \alpha_b:fit
      xb-fit -> gboolean: xb-fit
        x_b:fit
      alphaem-fit -> gboolean: alphaem-fit
        \alpha_{\mathrm{em}}:fit
      betaem-fit -> gboolean: betaem-fit
        \beta_\mathrm{em}:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        H0: float
        H0_fit: bool
        OmegaL: float
        OmegaL_fit: bool
        Omegac: float
        Omegac_fit: bool
        alphab: float
        alphab_fit: bool
        alphaem: float
        alphaem_fit: bool
        betaem: float
        betaem_fit: bool
        dphi: float
        dphi_fit: bool
        em_coupling: HICosmoVexpEMCoupling
        glue_de: bool
        set_xb_max: bool
        sigmaphi: float
        sigmaphi_fit: bool
        xb: float
        xb_fit: bool
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: HICosmo = ...
    def __init__(
        self,
        H0: float = ...,
        H0_fit: bool = ...,
        OmegaL: float = ...,
        OmegaL_fit: bool = ...,
        Omegac: float = ...,
        Omegac_fit: bool = ...,
        alphab: float = ...,
        alphab_fit: bool = ...,
        alphaem: float = ...,
        alphaem_fit: bool = ...,
        betaem: float = ...,
        betaem_fit: bool = ...,
        dphi: float = ...,
        dphi_fit: bool = ...,
        em_coupling: HICosmoVexpEMCoupling = ...,
        glue_de: bool = ...,
        set_xb_max: bool = ...,
        sigmaphi: float = ...,
        sigmaphi_fit: bool = ...,
        xb: float = ...,
        xb_fit: bool = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    def E_tau(self, tau: float) -> float: ...
    def Ricci_scale(self, tau: float) -> float: ...
    def alpha(self, tau: float) -> float: ...
    def alpha_0c(self) -> float: ...
    def alpha_0e(self) -> float: ...
    def get_em_coupling(self) -> HICosmoVexpEMCoupling: ...
    @classmethod
    def new(cls) -> HICosmoVexp: ...
    def phi(self, tau: float) -> float: ...
    def set_em_coupling(self, coupling: HICosmoVexpEMCoupling) -> None: ...
    def tau_max(self) -> float: ...
    def tau_min(self) -> float: ...
    def tau_qt_c(self) -> float: ...
    def tau_qt_e(self) -> float: ...
    def tau_xc(self, xc: float) -> float: ...
    def tau_xe(self, xe: float) -> float: ...
    def x_y(self, tau: float) -> typing.Tuple[float, float]: ...
    def xbc(self) -> float: ...
    def xbe(self) -> float: ...
    def xc_tau(self, tau: float) -> float: ...
    def xe_tau(self, tau: float) -> float: ...

class HICosmoVexpClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HICosmoVexpClass()
    """

    parent_class: HICosmoClass = ...

class HICosmoVexpPrivate(GObject.GPointer): ...

class HIPert(GObject.Object):
    r"""
    :Constructors:

    ::

        HIPert(**properties)

    Object NcHIPert

    Properties from NcHIPert:
      k -> gdouble: k
        Mode k
      sys-size -> guint: sys-size
        System size
      reltol -> gdouble: reltol
        Relative tolerance
      abstol -> gdouble: abstol
        Absolute tolerance
      alphai -> gdouble: alphai
        Initial time

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        abstol: float
        alphai: float
        k: float
        reltol: float
        sys_size: int

    props: Props = ...
    parent_instance: GObject.Object = ...
    priv: HIPertPrivate = ...
    def __init__(
        self,
        abstol: float = ...,
        alphai: float = ...,
        k: float = ...,
        reltol: float = ...,
        sys_size: int = ...,
    ) -> None: ...
    def do_set_abstol(self, abstol: float) -> None: ...
    def do_set_mode_k(self, k: float) -> None: ...
    def do_set_reltol(self, reltol: float) -> None: ...
    def get_abstol(self) -> float: ...
    def get_mode_k(self) -> float: ...
    def get_reltol(self) -> float: ...
    def prepared(self) -> bool: ...
    def reset_solver(self) -> None: ...
    def set_abstol(self, abstol: float) -> None: ...
    def set_mode_k(self, k: float) -> None: ...
    def set_prepared(self, prepared: bool) -> None: ...
    def set_reltol(self, reltol: float) -> None: ...
    def set_stiff_solver(self, stiff: bool) -> None: ...
    def set_sys_size(self, sys_size: int) -> None: ...

class HIPertAdiab(NumCosmoMath.CSQ1D):
    r"""
    :Constructors:

    ::

        HIPertAdiab(**properties)
        new() -> NumCosmo.HIPertAdiab

    Object NcHIPertAdiab

    Properties from NcHIPertAdiab:
      k -> gdouble: k
        Wave number

    Properties from NcmCSQ1D:
      reltol -> gdouble: reltol
        Relative tolerance
      abstol -> gdouble: abstol
        Absolute tolerance tolerance
      ti -> gdouble: ti
        The initial time t_i
      tf -> gdouble: tf
        The final time t_f
      adiab-threshold -> gdouble: adiab-threshold
        The adiabatic threshold
      prop-threshold -> gdouble: prop-threshold
        The propagator threshold
      save-evol -> gboolean: save-evol
        Save the system evolution
      vacuum-type -> NcmCSQ1DInitialStateType: vacuum-type
        The vacuum type
      vacuum-reltol -> gdouble: vacuum-reltol
        The vacuum relative tolerance
      vacuum-max-time -> gdouble: vacuum-max-time
        The vacuum maximum time

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        k: float
        abstol: float
        adiab_threshold: float
        prop_threshold: float
        reltol: float
        save_evol: bool
        tf: float
        ti: float
        vacuum_max_time: float
        vacuum_reltol: float
        vacuum_type: NumCosmoMath.CSQ1DInitialStateType

    props: Props = ...
    def __init__(
        self,
        k: float = ...,
        abstol: float = ...,
        adiab_threshold: float = ...,
        prop_threshold: float = ...,
        reltol: float = ...,
        save_evol: bool = ...,
        tf: float = ...,
        ti: float = ...,
        vacuum_max_time: float = ...,
        vacuum_reltol: float = ...,
        vacuum_type: NumCosmoMath.CSQ1DInitialStateType = ...,
    ) -> None: ...
    @staticmethod
    def clear(pa: HIPertAdiab) -> None: ...
    def eval_cosmic_time(self, model: NumCosmoMath.Model, tau: float) -> float: ...
    def eval_delta_critial(self, model: NumCosmoMath.Model, tau: float) -> float: ...
    def eval_powspec_Psi(
        self, model: NumCosmoMath.Model
    ) -> NumCosmoMath.PowspecSpline2d: ...
    def eval_powspec_Psi_at(self, model: NumCosmoMath.Model, tau: float) -> float: ...
    def eval_powspec_drho(
        self, model: NumCosmoMath.Model
    ) -> NumCosmoMath.PowspecSpline2d: ...
    def eval_powspec_drho_at(self, model: NumCosmoMath.Model, tau: float) -> float: ...
    def eval_powspec_zeta(
        self, model: NumCosmoMath.Model
    ) -> NumCosmoMath.PowspecSpline2d: ...
    def eval_powspec_zeta_at(self, model: NumCosmoMath.Model, tau: float) -> float: ...
    def free(self) -> None: ...
    def get_k(self) -> float: ...
    @classmethod
    def new(cls) -> HIPertAdiab: ...
    def prepare_spectrum(
        self,
        model: NumCosmoMath.Model,
        k_array: typing.Sequence[float],
        tau_array: typing.Sequence[float],
    ) -> None: ...
    def ref(self) -> HIPertAdiab: ...
    def set_k(self, k: float) -> None: ...

class HIPertAdiabClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPertAdiabClass()
    """

    parent_class: NumCosmoMath.CSQ1DClass = ...

class HIPertBGVar(GObject.Object):
    r"""
    :Constructors:

    ::

        HIPertBGVar(**properties)
        new() -> NumCosmo.HIPertBGVar
        new_full(dist:NumCosmo.Distance, recomb:NumCosmo.Recomb, a:NumCosmo.Scalefactor) -> NumCosmo.HIPertBGVar

    Object NcHIPertBGVar

    Properties from NcHIPertBGVar:
      distance -> NcDistance: distance
        Distance object
      recomb -> NcRecomb: recomb
        Recombination object
      scalefactor -> NcScalefactor: scalefactor
        Scalefactor object
      zf -> gdouble: zf
        Maximum redshift

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        distance: Distance
        recomb: typing.Optional[Recomb]
        scalefactor: typing.Optional[Scalefactor]
        zf: float

    props: Props = ...
    parent_instance: GObject.Object = ...
    priv: HIPertBGVarPrivate = ...
    cstructs: list[None] = ...
    recomb: Recomb = ...
    dist: Distance = ...
    a: Scalefactor = ...
    t: float = ...
    eta: float = ...
    k: float = ...
    x: float = ...
    E: float = ...
    def __init__(
        self,
        distance: Distance = ...,
        recomb: Recomb = ...,
        scalefactor: Scalefactor = ...,
        zf: float = ...,
    ) -> None: ...
    def activate_id_array(self, ids: typing.Sequence[int]) -> None: ...
    @staticmethod
    def clear(bg_var: HIPertBGVar) -> None: ...
    def free(self) -> None: ...
    def get_dist(self) -> typing.Optional[Distance]: ...
    def get_recomb(self) -> typing.Optional[Recomb]: ...
    def get_scalefactor(self) -> typing.Optional[Scalefactor]: ...
    def get_zf(self) -> float: ...
    def len(self) -> int: ...
    @classmethod
    def new(cls) -> HIPertBGVar: ...
    @classmethod
    def new_full(
        cls, dist: Distance, recomb: Recomb, a: Scalefactor
    ) -> HIPertBGVar: ...
    def peek_dist(self) -> typing.Optional[Distance]: ...
    def peek_recomb(self) -> typing.Optional[Recomb]: ...
    def peek_scalefactor(self) -> typing.Optional[Scalefactor]: ...
    def prepare(self, cosmo: HICosmo) -> None: ...
    def prepare_if_needed(self, cosmo: HICosmo) -> None: ...
    def ref(self) -> HIPertBGVar: ...
    def set_dist(self, dist: Distance) -> None: ...
    def set_recomb(self, recomb: Recomb) -> None: ...
    def set_scalefactor(self, a: Scalefactor) -> None: ...
    def set_zf(self, zf: float) -> None: ...

class HIPertBGVarClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPertBGVarClass()
    """

    parent_class: GObject.ObjectClass = ...
    bg_var_id_len: int = ...
    ns_table: dict[None, None] = ...
    bg_var_desc_array: list[None] = ...

class HIPertBGVarDesc(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPertBGVarDesc()
    """

    init: bool = ...
    ns: str = ...
    desc: str = ...
    long_desc: str = ...
    cstruct_size: int = ...

class HIPertBGVarPrivate(GObject.GPointer): ...

class HIPertBGVarYDY(GObject.GBoxed):
    r"""
    :Constructors:

    ::

        HIPertBGVarYDY()
        new() -> NumCosmo.HIPertBGVarYDY
    """

    y: float = ...
    dy: float = ...
    start_index: int = ...
    perm: int = ...
    perm_inv: int = ...
    def dup(self) -> HIPertBGVarYDY: ...
    def free(self) -> None: ...
    def get_dy_i(self, i: int) -> float: ...
    def get_y_i(self, i: int) -> float: ...
    @classmethod
    def new(cls) -> HIPertBGVarYDY: ...
    def set_dy_i(self, i: int, dy_i: float) -> None: ...

class HIPertBoltzmann(HIPert):
    r"""
    :Constructors:

    ::

        HIPertBoltzmann(**properties)

    Object NcHIPertBoltzmann

    Properties from NcHIPertBoltzmann:
      recomb -> NcRecomb: recomb
        Recombination object
      target-Cls -> NcDataCMBDataType: target-Cls
        Which Cls must be calculated
      calc-transfer -> gboolean: calc-transfer
        Whether to calculate the matter transfer function
      use-lensed-Cls -> gboolean: use-lensed-Cls
        Whether use the lensed corrected Cls
      use-tensor -> gboolean: use-tensor
        Whether use tensor contribution
      PHIPHI-l-max -> guint: PHIPHI-l-max
        Last multipole in the PHIPHI correlation
      TT-l-max -> guint: TT-l-max
        Last multipole in the TT correlation
      EE-l-max -> guint: EE-l-max
        Last multipole in the EE correlation
      BB-l-max -> guint: BB-l-max
        Last multipole in the BB correlation
      TE-l-max -> guint: TE-l-max
        Last multipole in the TE correlation
      TB-l-max -> guint: TB-l-max
        Last multipole in the TB correlation
      EB-l-max -> guint: EB-l-max
        Last multipole in the EB correlation

    Properties from NcHIPert:
      k -> gdouble: k
        Mode k
      sys-size -> guint: sys-size
        System size
      reltol -> gdouble: reltol
        Relative tolerance
      abstol -> gdouble: abstol
        Absolute tolerance
      alphai -> gdouble: alphai
        Initial time

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        BB_l_max: int
        EB_l_max: int
        EE_l_max: int
        PHIPHI_l_max: int
        TB_l_max: int
        TE_l_max: int
        TT_l_max: int
        calc_transfer: bool
        recomb: Recomb
        target_Cls: DataCMBDataType
        use_lensed_Cls: bool
        use_tensor: bool
        abstol: float
        alphai: float
        k: float
        reltol: float
        sys_size: int

    props: Props = ...
    parent_instance: HIPert = ...
    recomb: Recomb = ...
    cosmo: HICosmo = ...
    a: Scalefactor = ...
    eta0: float = ...
    lambdai: float = ...
    lambdaf: float = ...
    lambda_opt_cutoff: float = ...
    lambda_rec: float = ...
    lambda_rec_10m2_max: list[float] = ...
    lambda_: float = ...
    target_Cls: DataCMBDataType = ...
    calc_transfer: bool = ...
    use_lensed_Cls: bool = ...
    use_tensor: bool = ...
    PHIPHI_lmax: int = ...
    TT_lmax: int = ...
    EE_lmax: int = ...
    BB_lmax: int = ...
    TE_lmax: int = ...
    TB_lmax: int = ...
    EB_lmax: int = ...
    tight_coupling: bool = ...
    ctrl_cosmo: NumCosmoMath.ModelCtrl = ...
    ctrl_prim: NumCosmoMath.ModelCtrl = ...
    def __init__(
        self,
        BB_l_max: int = ...,
        EB_l_max: int = ...,
        EE_l_max: int = ...,
        PHIPHI_l_max: int = ...,
        TB_l_max: int = ...,
        TE_l_max: int = ...,
        TT_l_max: int = ...,
        calc_transfer: bool = ...,
        recomb: Recomb = ...,
        target_Cls: DataCMBDataType = ...,
        use_lensed_Cls: bool = ...,
        use_tensor: bool = ...,
        abstol: float = ...,
        alphai: float = ...,
        k: float = ...,
        reltol: float = ...,
        sys_size: int = ...,
    ) -> None: ...
    def append_target_Cls(self, tCls: DataCMBDataType) -> None: ...
    @staticmethod
    def clear(pb: HIPertBoltzmann) -> None: ...
    def do_evol(self, g: float) -> None: ...
    def do_evol_step(self, g: float) -> None: ...
    def do_get(self, n: int) -> float: ...
    def do_get_BB_Cls(self, Cls: NumCosmoMath.Vector) -> None: ...
    def do_get_EB_Cls(self, Cls: NumCosmoMath.Vector) -> None: ...
    def do_get_EE_Cls(self, Cls: NumCosmoMath.Vector) -> None: ...
    def do_get_PHIPHI_Cls(self, Cls: NumCosmoMath.Vector) -> None: ...
    def do_get_TB_Cls(self, Cls: NumCosmoMath.Vector) -> None: ...
    def do_get_TE_Cls(self, Cls: NumCosmoMath.Vector) -> None: ...
    def do_get_TT_Cls(self, Cls: NumCosmoMath.Vector) -> None: ...
    def do_get_b0(self) -> float: ...
    def do_get_b1(self) -> float: ...
    def do_get_c0(self) -> float: ...
    def do_get_c1(self) -> float: ...
    def do_get_los_theta(self, n: int) -> float: ...
    def do_get_phi(self) -> float: ...
    def do_get_sources(self, S0: float, S1: float, S2: float) -> None: ...
    def do_get_theta(self, n: int) -> float: ...
    def do_get_theta_p(self, n: int) -> float: ...
    def do_get_z(self) -> float: ...
    def do_init(self, cosmo: HICosmo) -> None: ...
    def do_prepare(self, cosmo: HICosmo) -> None: ...
    def do_prepare_if_needed(self, cosmo: HICosmo) -> None: ...
    def do_print_all(self) -> None: ...
    def do_print_stats(self) -> None: ...
    def do_reset(self) -> None: ...
    def do_set_opts(self) -> None: ...
    def free(self) -> None: ...
    def get_BB_Cls(self, Cls: NumCosmoMath.Vector) -> None: ...
    def get_BB_lmax(self) -> int: ...
    def get_EB_Cls(self, Cls: NumCosmoMath.Vector) -> None: ...
    def get_EB_lmax(self) -> int: ...
    def get_EE_Cls(self, Cls: NumCosmoMath.Vector) -> None: ...
    def get_EE_lmax(self) -> int: ...
    def get_PHIPHI_Cls(self, Cls: NumCosmoMath.Vector) -> None: ...
    def get_PHIPHI_lmax(self) -> int: ...
    def get_TB_Cls(self, Cls: NumCosmoMath.Vector) -> None: ...
    def get_TB_lmax(self) -> int: ...
    def get_TE_Cls(self, Cls: NumCosmoMath.Vector) -> None: ...
    def get_TE_lmax(self) -> int: ...
    def get_TT_Cls(self, Cls: NumCosmoMath.Vector) -> None: ...
    def get_TT_lmax(self) -> int: ...
    def get_calc_transfer(self) -> bool: ...
    def get_target_Cls(self) -> DataCMBDataType: ...
    def lensed_Cls(self) -> bool: ...
    def prepare(self, cosmo: HICosmo) -> None: ...
    def prepare_if_needed(self, cosmo: HICosmo) -> None: ...
    def ref(self) -> HIPertBoltzmann: ...
    def set_BB_lmax(self, lmax: int) -> None: ...
    def set_EB_lmax(self, lmax: int) -> None: ...
    def set_EE_lmax(self, lmax: int) -> None: ...
    def set_PHIPHI_lmax(self, lmax: int) -> None: ...
    def set_TB_lmax(self, lmax: int) -> None: ...
    def set_TE_lmax(self, lmax: int) -> None: ...
    def set_TT_lmax(self, lmax: int) -> None: ...
    def set_calc_transfer(self, calc_transfer: bool) -> None: ...
    def set_lensed_Cls(self, use_lensed_Cls: bool) -> None: ...
    def set_recomb(self, recomb: Recomb) -> None: ...
    def set_target_Cls(self, tCls: DataCMBDataType) -> None: ...
    def set_tensor(self, use_tensor: bool) -> None: ...
    def tensor(self) -> bool: ...

class HIPertBoltzmannCBE(HIPertBoltzmann):
    r"""
    :Constructors:

    ::

        HIPertBoltzmannCBE(**properties)
        full_new(cbe:NumCosmo.CBE) -> NumCosmo.HIPertBoltzmannCBE
        new() -> NumCosmo.HIPertBoltzmannCBE

    Object NcHIPertBoltzmannCBE

    Properties from NcHIPertBoltzmannCBE:
      cbe -> NcCBE: cbe
        CLASS backend object

    Properties from NcHIPertBoltzmann:
      recomb -> NcRecomb: recomb
        Recombination object
      target-Cls -> NcDataCMBDataType: target-Cls
        Which Cls must be calculated
      calc-transfer -> gboolean: calc-transfer
        Whether to calculate the matter transfer function
      use-lensed-Cls -> gboolean: use-lensed-Cls
        Whether use the lensed corrected Cls
      use-tensor -> gboolean: use-tensor
        Whether use tensor contribution
      PHIPHI-l-max -> guint: PHIPHI-l-max
        Last multipole in the PHIPHI correlation
      TT-l-max -> guint: TT-l-max
        Last multipole in the TT correlation
      EE-l-max -> guint: EE-l-max
        Last multipole in the EE correlation
      BB-l-max -> guint: BB-l-max
        Last multipole in the BB correlation
      TE-l-max -> guint: TE-l-max
        Last multipole in the TE correlation
      TB-l-max -> guint: TB-l-max
        Last multipole in the TB correlation
      EB-l-max -> guint: EB-l-max
        Last multipole in the EB correlation

    Properties from NcHIPert:
      k -> gdouble: k
        Mode k
      sys-size -> guint: sys-size
        System size
      reltol -> gdouble: reltol
        Relative tolerance
      abstol -> gdouble: abstol
        Absolute tolerance
      alphai -> gdouble: alphai
        Initial time

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        cbe: CBE
        BB_l_max: int
        EB_l_max: int
        EE_l_max: int
        PHIPHI_l_max: int
        TB_l_max: int
        TE_l_max: int
        TT_l_max: int
        calc_transfer: bool
        recomb: Recomb
        target_Cls: DataCMBDataType
        use_lensed_Cls: bool
        use_tensor: bool
        abstol: float
        alphai: float
        k: float
        reltol: float
        sys_size: int

    props: Props = ...
    parent_instance: HIPertBoltzmann = ...
    cbe: CBE = ...
    PHIPHI_Cls: NumCosmoMath.Vector = ...
    TT_Cls: NumCosmoMath.Vector = ...
    EE_Cls: NumCosmoMath.Vector = ...
    BB_Cls: NumCosmoMath.Vector = ...
    TE_Cls: NumCosmoMath.Vector = ...
    TB_Cls: NumCosmoMath.Vector = ...
    EB_Cls: NumCosmoMath.Vector = ...
    def __init__(
        self,
        cbe: CBE = ...,
        BB_l_max: int = ...,
        EB_l_max: int = ...,
        EE_l_max: int = ...,
        PHIPHI_l_max: int = ...,
        TB_l_max: int = ...,
        TE_l_max: int = ...,
        TT_l_max: int = ...,
        calc_transfer: bool = ...,
        recomb: Recomb = ...,
        target_Cls: DataCMBDataType = ...,
        use_lensed_Cls: bool = ...,
        use_tensor: bool = ...,
        abstol: float = ...,
        alphai: float = ...,
        k: float = ...,
        reltol: float = ...,
        sys_size: int = ...,
    ) -> None: ...
    @staticmethod
    def clear(boltzmann_cbe: HIPertBoltzmannCBE) -> None: ...
    def free(self) -> None: ...
    @classmethod
    def full_new(cls, cbe: CBE) -> HIPertBoltzmannCBE: ...
    @classmethod
    def new(cls) -> HIPertBoltzmannCBE: ...
    def peek_cbe(self) -> CBE: ...
    def ref(self) -> HIPertBoltzmannCBE: ...

class HIPertBoltzmannCBEClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPertBoltzmannCBEClass()
    """

    parent_class: HIPertBoltzmannClass = ...

class HIPertBoltzmannClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPertBoltzmannClass()
    """

    parent_class: HIPertClass = ...
    init: typing.Callable[[HIPertBoltzmann, HICosmo], None] = ...
    set_opts: typing.Callable[[HIPertBoltzmann], None] = ...
    reset: typing.Callable[[HIPertBoltzmann], None] = ...
    evol_step: typing.Callable[[HIPertBoltzmann, float], None] = ...
    evol: typing.Callable[[HIPertBoltzmann, float], None] = ...
    prepare: typing.Callable[[HIPertBoltzmann, HICosmo], None] = ...
    prepare_if_needed: typing.Callable[[HIPertBoltzmann, HICosmo], None] = ...
    get_sources: typing.Callable[[HIPertBoltzmann, float, float, float], None] = ...
    print_stats: typing.Callable[[HIPertBoltzmann], None] = ...
    get_z: typing.Callable[[HIPertBoltzmann], float] = ...
    get_phi: typing.Callable[[HIPertBoltzmann], float] = ...
    get_c0: typing.Callable[[HIPertBoltzmann], float] = ...
    get_b0: typing.Callable[[HIPertBoltzmann], float] = ...
    get_c1: typing.Callable[[HIPertBoltzmann], float] = ...
    get_b1: typing.Callable[[HIPertBoltzmann], float] = ...
    get: typing.Callable[[HIPertBoltzmann, int], float] = ...
    get_theta: typing.Callable[[HIPertBoltzmann, int], float] = ...
    get_theta_p: typing.Callable[[HIPertBoltzmann, int], float] = ...
    get_los_theta: typing.Callable[[HIPertBoltzmann, int], float] = ...
    get_PHIPHI_Cls: typing.Callable[[HIPertBoltzmann, NumCosmoMath.Vector], None] = ...
    get_TT_Cls: typing.Callable[[HIPertBoltzmann, NumCosmoMath.Vector], None] = ...
    get_EE_Cls: typing.Callable[[HIPertBoltzmann, NumCosmoMath.Vector], None] = ...
    get_BB_Cls: typing.Callable[[HIPertBoltzmann, NumCosmoMath.Vector], None] = ...
    get_TE_Cls: typing.Callable[[HIPertBoltzmann, NumCosmoMath.Vector], None] = ...
    get_TB_Cls: typing.Callable[[HIPertBoltzmann, NumCosmoMath.Vector], None] = ...
    get_EB_Cls: typing.Callable[[HIPertBoltzmann, NumCosmoMath.Vector], None] = ...
    print_all: typing.Callable[[HIPertBoltzmann], None] = ...
    data: None = ...

class HIPertBoltzmannStd(HIPertBoltzmann):
    r"""
    :Constructors:

    ::

        HIPertBoltzmannStd(**properties)
        new(recomb:NumCosmo.Recomb, lmax:int) -> NumCosmo.HIPertBoltzmannStd

    Object NcHIPertBoltzmannStd

    Properties from NcHIPertBoltzmannStd:
      l-maxa -> guint: l-maxa
        Last multipole

    Properties from NcHIPertBoltzmann:
      recomb -> NcRecomb: recomb
        Recombination object
      target-Cls -> NcDataCMBDataType: target-Cls
        Which Cls must be calculated
      calc-transfer -> gboolean: calc-transfer
        Whether to calculate the matter transfer function
      use-lensed-Cls -> gboolean: use-lensed-Cls
        Whether use the lensed corrected Cls
      use-tensor -> gboolean: use-tensor
        Whether use tensor contribution
      PHIPHI-l-max -> guint: PHIPHI-l-max
        Last multipole in the PHIPHI correlation
      TT-l-max -> guint: TT-l-max
        Last multipole in the TT correlation
      EE-l-max -> guint: EE-l-max
        Last multipole in the EE correlation
      BB-l-max -> guint: BB-l-max
        Last multipole in the BB correlation
      TE-l-max -> guint: TE-l-max
        Last multipole in the TE correlation
      TB-l-max -> guint: TB-l-max
        Last multipole in the TB correlation
      EB-l-max -> guint: EB-l-max
        Last multipole in the EB correlation

    Properties from NcHIPert:
      k -> gdouble: k
        Mode k
      sys-size -> guint: sys-size
        System size
      reltol -> gdouble: reltol
        Relative tolerance
      abstol -> gdouble: abstol
        Absolute tolerance
      alphai -> gdouble: alphai
        Initial time

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        l_maxa: int
        BB_l_max: int
        EB_l_max: int
        EE_l_max: int
        PHIPHI_l_max: int
        TB_l_max: int
        TE_l_max: int
        TT_l_max: int
        calc_transfer: bool
        recomb: Recomb
        target_Cls: DataCMBDataType
        use_lensed_Cls: bool
        use_tensor: bool
        abstol: float
        alphai: float
        k: float
        reltol: float
        sys_size: int

    props: Props = ...
    parent_instance: HIPertBoltzmann = ...
    def __init__(
        self,
        l_maxa: int = ...,
        BB_l_max: int = ...,
        EB_l_max: int = ...,
        EE_l_max: int = ...,
        PHIPHI_l_max: int = ...,
        TB_l_max: int = ...,
        TE_l_max: int = ...,
        TT_l_max: int = ...,
        calc_transfer: bool = ...,
        recomb: Recomb = ...,
        target_Cls: DataCMBDataType = ...,
        use_lensed_Cls: bool = ...,
        use_tensor: bool = ...,
        abstol: float = ...,
        alphai: float = ...,
        k: float = ...,
        reltol: float = ...,
        sys_size: int = ...,
    ) -> None: ...
    @classmethod
    def new(cls, recomb: Recomb, lmax: int) -> HIPertBoltzmannStd: ...

class HIPertBoltzmannStdClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPertBoltzmannStdClass()
    """

    parent_class: HIPertBoltzmannClass = ...

class HIPertClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPertClass()
    """

    parent_class: GObject.ObjectClass = ...
    set_mode_k: typing.Callable[[HIPert, float], None] = ...
    set_reltol: typing.Callable[[HIPert, float], None] = ...
    set_abstol: typing.Callable[[HIPert, float], None] = ...

class HIPertComp(GObject.Object):
    r"""
    :Constructors:

    ::

        HIPertComp(**properties)

    Object NcHIPertComp

    Properties from NcHIPertComp:
      gauge -> NcHIPertGravGauge: gauge
        gauge

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        gauge: HIPertGravGauge

    props: Props = ...
    parent_instance: GObject.Object = ...
    priv: HIPertCompPrivate = ...
    def __init__(self, gauge: HIPertGravGauge = ...) -> None: ...
    @staticmethod
    def clear(comp: HIPertComp) -> None: ...
    def do_get_T_scalar(
        self, bg_var: HIPertBGVar, ydy: HIPertBGVarYDY, T_scalar: HIPertGravTScalar
    ) -> None: ...
    def do_get_T_scalar_info(self) -> HIPertGravTScalarInfo: ...
    def do_get_T_tensor(
        self, bg_var: HIPertBGVar, ydy: HIPertBGVarYDY, T_tensor: HIPertGravTTensor
    ) -> None: ...
    def do_get_T_vector(
        self, bg_var: HIPertBGVar, ydy: HIPertBGVarYDY, T_vector: HIPertGravTVector
    ) -> None: ...
    def do_get_deps(self, vindex: int) -> list[int]: ...
    def do_get_dy_scalar(
        self,
        bg_var: HIPertBGVar,
        ydy: HIPertBGVarYDY,
        T_scalar: HIPertGravTScalar,
        G_scalar: HIPertGravScalar,
    ) -> None: ...
    def do_get_gauge(self) -> HIPertGravGauge: ...
    def do_ndyn_var(self) -> int: ...
    def do_set_gauge(self, gauge: HIPertGravGauge) -> None: ...
    def free(self) -> None: ...
    def get_T_scalar(
        self, bg_var: HIPertBGVar, ydy: HIPertBGVarYDY, T_scalar: HIPertGravTScalar
    ) -> None: ...
    def get_T_scalar_info(self) -> HIPertGravTScalarInfo: ...
    def get_deps(self, vindex: int) -> list[int]: ...
    def get_dy_scalar(
        self,
        bg_var: HIPertBGVar,
        ydy: HIPertBGVarYDY,
        T_scalar: HIPertGravTScalar,
        G_scalar: HIPertGravScalar,
    ) -> None: ...
    def get_gauge(self) -> HIPertGravGauge: ...
    def get_id(self) -> int: ...
    def ndyn_var(self) -> int: ...
    def ref(self) -> HIPertComp: ...
    def set_gauge(self, gauge: HIPertGravGauge) -> None: ...

class HIPertCompClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPertCompClass()
    """

    parent_class: GObject.ObjectClass = ...
    ndyn_var: typing.Callable[[HIPertComp], int] = ...
    get_deps: typing.Callable[[HIPertComp, int], list[int]] = ...
    set_gauge: typing.Callable[[HIPertComp, HIPertGravGauge], None] = ...
    get_gauge: typing.Callable[[HIPertComp], HIPertGravGauge] = ...
    get_T_scalar_info: typing.Callable[[HIPertComp], HIPertGravTScalarInfo] = ...
    get_T_scalar: typing.Callable[
        [HIPertComp, HIPertBGVar, HIPertBGVarYDY, HIPertGravTScalar], None
    ] = ...
    get_T_vector: typing.Callable[
        [HIPertComp, HIPertBGVar, HIPertBGVarYDY, HIPertGravTVector], None
    ] = ...
    get_T_tensor: typing.Callable[
        [HIPertComp, HIPertBGVar, HIPertBGVarYDY, HIPertGravTTensor], None
    ] = ...
    get_dy_scalar: typing.Callable[
        [HIPertComp, HIPertBGVar, HIPertBGVarYDY, HIPertGravTScalar, HIPertGravScalar],
        None,
    ] = ...

class HIPertCompPB(HIPertComp):
    r"""
    :Constructors:

    ::

        HIPertCompPB(**properties)
        new() -> NumCosmo.HIPertCompPB

    Object NcHIPertCompPB

    Properties from NcHIPertCompPB:
      l-max -> guint: l-max
        l_max

    Properties from NcHIPertComp:
      gauge -> NcHIPertGravGauge: gauge
        gauge

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        l_max: int
        gauge: HIPertGravGauge

    props: Props = ...
    parent_instance: HIPertComp = ...
    priv: HIPertCompPBPrivate = ...
    def __init__(self, l_max: int = ..., gauge: HIPertGravGauge = ...) -> None: ...
    @staticmethod
    def clear(pb: HIPertCompPB) -> None: ...
    def free(self) -> None: ...
    def get_lmax(self) -> int: ...
    @staticmethod
    def id() -> int: ...
    @classmethod
    def new(cls) -> HIPertCompPB: ...
    def ref(self) -> HIPertCompPB: ...
    def set_lmax(self, lmax: int) -> None: ...

class HIPertCompPBClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPertCompPBClass()
    """

    parent_class: HIPertCompClass = ...

class HIPertCompPBPrivate(GObject.GPointer): ...
class HIPertCompPrivate(GObject.GPointer): ...

class HIPertEM(NumCosmoMath.CSQ1D):
    r"""
    :Constructors:

    ::

        HIPertEM(**properties)
        new() -> NumCosmo.HIPertEM

    Object NcHIPertEM

    Properties from NcHIPertEM:
      k -> gdouble: k
        Wave number

    Properties from NcmCSQ1D:
      reltol -> gdouble: reltol
        Relative tolerance
      abstol -> gdouble: abstol
        Absolute tolerance tolerance
      ti -> gdouble: ti
        The initial time t_i
      tf -> gdouble: tf
        The final time t_f
      adiab-threshold -> gdouble: adiab-threshold
        The adiabatic threshold
      prop-threshold -> gdouble: prop-threshold
        The propagator threshold
      save-evol -> gboolean: save-evol
        Save the system evolution
      vacuum-type -> NcmCSQ1DInitialStateType: vacuum-type
        The vacuum type
      vacuum-reltol -> gdouble: vacuum-reltol
        The vacuum relative tolerance
      vacuum-max-time -> gdouble: vacuum-max-time
        The vacuum maximum time

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        k: float
        abstol: float
        adiab_threshold: float
        prop_threshold: float
        reltol: float
        save_evol: bool
        tf: float
        ti: float
        vacuum_max_time: float
        vacuum_reltol: float
        vacuum_type: NumCosmoMath.CSQ1DInitialStateType

    props: Props = ...
    def __init__(
        self,
        k: float = ...,
        abstol: float = ...,
        adiab_threshold: float = ...,
        prop_threshold: float = ...,
        reltol: float = ...,
        save_evol: bool = ...,
        tf: float = ...,
        ti: float = ...,
        vacuum_max_time: float = ...,
        vacuum_reltol: float = ...,
        vacuum_type: NumCosmoMath.CSQ1DInitialStateType = ...,
    ) -> None: ...
    @staticmethod
    def clear(pem: HIPertEM) -> None: ...
    def eval_PE_PB(
        self, model: NumCosmoMath.Model, tau: float
    ) -> typing.Tuple[float, float]: ...
    def free(self) -> None: ...
    def get_k(self) -> float: ...
    @classmethod
    def new(cls) -> HIPertEM: ...
    def prepare_spectrum(
        self,
        model: NumCosmoMath.Model,
        k_array: typing.Sequence[float],
        tau_array: typing.Sequence[float],
    ) -> typing.Tuple[NumCosmoMath.PowspecSpline2d, NumCosmoMath.PowspecSpline2d]: ...
    def ref(self) -> HIPertEM: ...
    def set_k(self, k: float) -> None: ...

class HIPertEMClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPertEMClass()
    """

    parent_class: NumCosmoMath.CSQ1DClass = ...

class HIPertFirstOrder(HIPertBoltzmann):
    r"""
    :Constructors:

    ::

        HIPertFirstOrder(**properties)
        new() -> NumCosmo.HIPertFirstOrder
        new_full(dist:NumCosmo.Distance, recomb:NumCosmo.Recomb, a:NumCosmo.Scalefactor) -> NumCosmo.HIPertFirstOrder

    Object NcHIPertFirstOrder

    Properties from NcHIPertFirstOrder:
      gauge -> NcHIPertGravGauge: gauge
        Gauge
      grav -> NcHIPertGrav: grav
        Gravitation object
      comp-array -> NcmObjArray: comp-array
        Components array
      distance -> NcDistance: distance
        Distance object
      recomb -> NcRecomb: recomb
        Recombination object
      scalefactor -> NcScalefactor: scalefactor
        Scale factor object
      reltol -> gdouble: reltol
        Relative tolerance
      abstol -> gdouble: abstol
        Absolute tolerance tolerance
      integ -> NcHIPertFirstOrderInteg: integ
        ODE integrator

    Properties from NcHIPertBoltzmann:
      recomb -> NcRecomb: recomb
        Recombination object
      target-Cls -> NcDataCMBDataType: target-Cls
        Which Cls must be calculated
      calc-transfer -> gboolean: calc-transfer
        Whether to calculate the matter transfer function
      use-lensed-Cls -> gboolean: use-lensed-Cls
        Whether use the lensed corrected Cls
      use-tensor -> gboolean: use-tensor
        Whether use tensor contribution
      PHIPHI-l-max -> guint: PHIPHI-l-max
        Last multipole in the PHIPHI correlation
      TT-l-max -> guint: TT-l-max
        Last multipole in the TT correlation
      EE-l-max -> guint: EE-l-max
        Last multipole in the EE correlation
      BB-l-max -> guint: BB-l-max
        Last multipole in the BB correlation
      TE-l-max -> guint: TE-l-max
        Last multipole in the TE correlation
      TB-l-max -> guint: TB-l-max
        Last multipole in the TB correlation
      EB-l-max -> guint: EB-l-max
        Last multipole in the EB correlation

    Properties from NcHIPert:
      k -> gdouble: k
        Mode k
      sys-size -> guint: sys-size
        System size
      reltol -> gdouble: reltol
        Relative tolerance
      abstol -> gdouble: abstol
        Absolute tolerance
      alphai -> gdouble: alphai
        Initial time

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        abstol: float
        comp_array: NumCosmoMath.ObjArray
        distance: Distance
        gauge: HIPertGravGauge
        grav: typing.Optional[HIPertGrav]
        integ: HIPertFirstOrderInteg
        recomb: Recomb
        reltol: float
        scalefactor: Scalefactor
        BB_l_max: int
        EB_l_max: int
        EE_l_max: int
        PHIPHI_l_max: int
        TB_l_max: int
        TE_l_max: int
        TT_l_max: int
        calc_transfer: bool
        target_Cls: DataCMBDataType
        use_lensed_Cls: bool
        use_tensor: bool
        alphai: float
        k: float
        sys_size: int

    props: Props = ...
    parent_instance: HIPertBoltzmann = ...
    priv: HIPertFirstOrderPrivate = ...
    def __init__(
        self,
        abstol: float = ...,
        comp_array: NumCosmoMath.ObjArray = ...,
        distance: Distance = ...,
        gauge: HIPertGravGauge = ...,
        grav: HIPertGrav = ...,
        integ: HIPertFirstOrderInteg = ...,
        recomb: Recomb = ...,
        reltol: float = ...,
        scalefactor: Scalefactor = ...,
        BB_l_max: int = ...,
        EB_l_max: int = ...,
        EE_l_max: int = ...,
        PHIPHI_l_max: int = ...,
        TB_l_max: int = ...,
        TE_l_max: int = ...,
        TT_l_max: int = ...,
        calc_transfer: bool = ...,
        target_Cls: DataCMBDataType = ...,
        use_lensed_Cls: bool = ...,
        use_tensor: bool = ...,
        alphai: float = ...,
        k: float = ...,
        sys_size: int = ...,
    ) -> None: ...
    def add_comp(self, comp: HIPertComp) -> None: ...
    @staticmethod
    def clear(fo: HIPertFirstOrder) -> None: ...
    def free(self) -> None: ...
    def get_abstol(self) -> float: ...
    def get_gauge(self) -> HIPertGravGauge: ...
    def get_grav(self) -> typing.Optional[HIPertGrav]: ...
    def get_integ(self) -> HIPertFirstOrderInteg: ...
    def get_reltol(self) -> float: ...
    @classmethod
    def new(cls) -> HIPertFirstOrder: ...
    @classmethod
    def new_full(
        cls, dist: Distance, recomb: Recomb, a: Scalefactor
    ) -> HIPertFirstOrder: ...
    def peek_grav(self) -> typing.Optional[HIPertGrav]: ...
    def prepare(self, cosmo: HICosmo) -> None: ...
    def ref(self) -> HIPertFirstOrder: ...
    def set_abstol(self, abstol: float) -> None: ...
    def set_gauge(self, gauge: HIPertGravGauge) -> None: ...
    def set_grav(self, grav: HIPertGrav) -> None: ...
    def set_integ(self, integ: HIPertFirstOrderInteg) -> None: ...
    def set_reltol(self, reltol: float) -> None: ...

class HIPertFirstOrderClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPertFirstOrderClass()
    """

    parent_class: HIPertBoltzmannClass = ...

class HIPertFirstOrderPrivate(GObject.GPointer): ...

class HIPertGW(NumCosmoMath.CSQ1D):
    r"""
    :Constructors:

    ::

        HIPertGW(**properties)
        new() -> NumCosmo.HIPertGW

    Object NcHIPertGW

    Properties from NcHIPertGW:
      k -> gdouble: k
        Wave number

    Properties from NcmCSQ1D:
      reltol -> gdouble: reltol
        Relative tolerance
      abstol -> gdouble: abstol
        Absolute tolerance tolerance
      ti -> gdouble: ti
        The initial time t_i
      tf -> gdouble: tf
        The final time t_f
      adiab-threshold -> gdouble: adiab-threshold
        The adiabatic threshold
      prop-threshold -> gdouble: prop-threshold
        The propagator threshold
      save-evol -> gboolean: save-evol
        Save the system evolution
      vacuum-type -> NcmCSQ1DInitialStateType: vacuum-type
        The vacuum type
      vacuum-reltol -> gdouble: vacuum-reltol
        The vacuum relative tolerance
      vacuum-max-time -> gdouble: vacuum-max-time
        The vacuum maximum time

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        k: float
        abstol: float
        adiab_threshold: float
        prop_threshold: float
        reltol: float
        save_evol: bool
        tf: float
        ti: float
        vacuum_max_time: float
        vacuum_reltol: float
        vacuum_type: NumCosmoMath.CSQ1DInitialStateType

    props: Props = ...
    def __init__(
        self,
        k: float = ...,
        abstol: float = ...,
        adiab_threshold: float = ...,
        prop_threshold: float = ...,
        reltol: float = ...,
        save_evol: bool = ...,
        tf: float = ...,
        ti: float = ...,
        vacuum_max_time: float = ...,
        vacuum_reltol: float = ...,
        vacuum_type: NumCosmoMath.CSQ1DInitialStateType = ...,
    ) -> None: ...
    @staticmethod
    def clear(pgw: HIPertGW) -> None: ...
    def free(self) -> None: ...
    def get_k(self) -> float: ...
    @classmethod
    def new(cls) -> HIPertGW: ...
    def ref(self) -> HIPertGW: ...
    def set_k(self, k: float) -> None: ...

class HIPertGWClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPertGWClass()
    """

    parent_class: NumCosmoMath.CSQ1DClass = ...

class HIPertGrav(GObject.Object):
    r"""
    :Constructors:

    ::

        HIPertGrav(**properties)

    Object NcHIPertGrav

    Properties from NcHIPertGrav:
      gauge -> NcHIPertGravGauge: gauge
        gauge

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        gauge: HIPertGravGauge

    props: Props = ...
    parent_instance: GObject.Object = ...
    priv: HIPertGravPrivate = ...
    def __init__(self, gauge: HIPertGravGauge = ...) -> None: ...
    @staticmethod
    def clear(grav: HIPertGrav) -> None: ...
    def do_get_G_scalar(
        self,
        bg_var: HIPertBGVar,
        ydy: HIPertBGVarYDY,
        T_scalar: HIPertGravTScalar,
        G_scalar: HIPertGravScalar,
    ) -> None: ...
    def do_get_G_scalar_info(self) -> HIPertGravInfo: ...
    def do_get_deps(self, vindex: int) -> list[int]: ...
    def do_get_dy_scalar(
        self,
        bg_var: HIPertBGVar,
        ydy: HIPertBGVarYDY,
        T_scalar: HIPertGravTScalar,
        G_scalar: HIPertGravScalar,
    ) -> None: ...
    def do_get_gauge(self) -> HIPertGravGauge: ...
    def do_ndyn_var(self) -> int: ...
    def do_set_gauge(self, gauge: HIPertGravGauge) -> None: ...
    def free(self) -> None: ...
    def get_G_scalar(
        self,
        bg_var: HIPertBGVar,
        ydy: HIPertBGVarYDY,
        T_scalar: HIPertGravTScalar,
        G_scalar: HIPertGravScalar,
    ) -> None: ...
    def get_G_scalar_info(self) -> HIPertGravInfo: ...
    def get_deps(self, vindex: int) -> list[int]: ...
    def get_dy_scalar(
        self,
        bg_var: HIPertBGVar,
        ydy: HIPertBGVarYDY,
        T_scalar: HIPertGravTScalar,
        G_scalar: HIPertGravScalar,
    ) -> None: ...
    def get_gauge(self) -> HIPertGravGauge: ...
    def get_id(self) -> int: ...
    def ndyn_var(self) -> int: ...
    def ref(self) -> HIPertGrav: ...
    def set_gauge(self, gauge: HIPertGravGauge) -> None: ...

class HIPertGravClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPertGravClass()
    """

    parent_class: GObject.ObjectClass = ...
    ndyn_var: typing.Callable[[HIPertGrav], int] = ...
    get_deps: typing.Callable[[HIPertGrav, int], list[int]] = ...
    set_gauge: typing.Callable[[HIPertGrav, HIPertGravGauge], None] = ...
    get_gauge: typing.Callable[[HIPertGrav], HIPertGravGauge] = ...
    get_G_scalar_info: typing.Callable[[HIPertGrav], HIPertGravInfo] = ...
    get_G_scalar: typing.Callable[
        [HIPertGrav, HIPertBGVar, HIPertBGVarYDY, HIPertGravTScalar, HIPertGravScalar],
        None,
    ] = ...
    get_dy_scalar: typing.Callable[
        [HIPertGrav, HIPertBGVar, HIPertBGVarYDY, HIPertGravTScalar, HIPertGravScalar],
        None,
    ] = ...

class HIPertGravEinstein(HIPertGrav):
    r"""
    :Constructors:

    ::

        HIPertGravEinstein(**properties)
        new() -> NumCosmo.HIPertGravEinstein

    Object NcHIPertGravEinstein

    Properties from NcHIPertGravEinstein:
      nhoc -> gint: nhoc
        nhoc

    Properties from NcHIPertGrav:
      gauge -> NcHIPertGravGauge: gauge
        gauge

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        nhoc: int
        gauge: HIPertGravGauge

    props: Props = ...
    parent_instance: HIPertGrav = ...
    priv: HIPertGravEinsteinPrivate = ...
    def __init__(self, nhoc: int = ..., gauge: HIPertGravGauge = ...) -> None: ...
    @staticmethod
    def clear(gr: HIPertGravEinstein) -> None: ...
    def free(self) -> None: ...
    @staticmethod
    def id() -> int: ...
    @classmethod
    def new(cls) -> HIPertGravEinstein: ...
    def ref(self) -> HIPertGravEinstein: ...

class HIPertGravEinsteinClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPertGravEinsteinClass()
    """

    parent_class: HIPertGravClass = ...

class HIPertGravEinsteinPrivate(GObject.GPointer): ...

class HIPertGravInfo(GObject.GBoxed):
    r"""
    :Constructors:

    ::

        HIPertGravInfo()
        new() -> NumCosmo.HIPertGravInfo
    """

    phi_deps: list[int] = ...
    dsigma_deps: list[int] = ...
    psi_deps: list[int] = ...
    dotpsi_deps: list[int] = ...
    def dup(self) -> HIPertGravInfo: ...
    def free(self) -> None: ...
    def get_dotpsi_deps(self) -> list[int]: ...
    def get_dsigma_deps(self) -> list[int]: ...
    def get_phi_deps(self) -> list[int]: ...
    def get_psi_deps(self) -> list[int]: ...
    @classmethod
    def new(cls) -> HIPertGravInfo: ...
    def set_dotpsi_deps(self, dotpsi_deps: typing.Sequence[int]) -> None: ...
    def set_dsigma_deps(self, dsigma_deps: typing.Sequence[int]) -> None: ...
    def set_phi_deps(self, phi_deps: typing.Sequence[int]) -> None: ...
    def set_psi_deps(self, psi_deps: typing.Sequence[int]) -> None: ...
    def set_zero(self) -> None: ...

class HIPertGravPrivate(GObject.GPointer): ...

class HIPertGravScalar(GObject.GBoxed):
    r"""
    :Constructors:

    ::

        HIPertGravScalar()
        new() -> NumCosmo.HIPertGravScalar
    """

    phi: float = ...
    dsigma: float = ...
    psi: float = ...
    def dup(self) -> HIPertGravScalar: ...
    def free(self) -> None: ...
    @classmethod
    def new(cls) -> HIPertGravScalar: ...
    def set_zero(self) -> None: ...

class HIPertGravTScalar(GObject.GBoxed):
    r"""
    :Constructors:

    ::

        HIPertGravTScalar()
        new() -> NumCosmo.HIPertGravTScalar
    """

    drho_m_Aphi: float = ...
    A: float = ...
    rhopp_v: float = ...
    dp: float = ...
    Pi: float = ...
    def add(self, Ts1: HIPertGravTScalar, Ts2: HIPertGravTScalar) -> None: ...
    def dup(self) -> HIPertGravTScalar: ...
    def free(self) -> None: ...
    @classmethod
    def new(cls) -> HIPertGravTScalar: ...
    def set_zero(self) -> None: ...

class HIPertGravTScalarInfo(GObject.GBoxed):
    r"""
    :Constructors:

    ::

        HIPertGravTScalarInfo()
        new() -> NumCosmo.HIPertGravTScalarInfo
    """

    drho_deps: list[None] = ...
    rhoppv_deps: list[None] = ...
    dp_deps: list[None] = ...
    dPi_deps: list[None] = ...
    def append(self, Tsinfo1: HIPertGravTScalarInfo) -> None: ...
    def dup(self) -> HIPertGravTScalarInfo: ...
    def free(self) -> None: ...
    @classmethod
    def new(cls) -> HIPertGravTScalarInfo: ...
    def set_zero(self) -> None: ...

class HIPertGravTTensor(GObject.GBoxed):
    r"""
    :Constructors:

    ::

        HIPertGravTTensor()
        new() -> NumCosmo.HIPertGravTTensor
    """

    a: float = ...
    def add(self, Tt1: HIPertGravTTensor, Tt2: HIPertGravTTensor) -> None: ...
    def dup(self) -> HIPertGravTTensor: ...
    def free(self) -> None: ...
    @classmethod
    def new(cls) -> HIPertGravTTensor: ...
    def set_zero(self) -> None: ...

class HIPertGravTVector(GObject.GBoxed):
    r"""
    :Constructors:

    ::

        HIPertGravTVector()
        new() -> NumCosmo.HIPertGravTVector
    """

    a: float = ...
    def add(self, Tv1: HIPertGravTVector, Tv2: HIPertGravTVector) -> None: ...
    def dup(self) -> HIPertGravTVector: ...
    def free(self) -> None: ...
    @classmethod
    def new(cls) -> HIPertGravTVector: ...
    def set_zero(self) -> None: ...

class HIPertGravTensor(GObject.GBoxed):
    r"""
    :Constructors:

    ::

        HIPertGravTensor()
        new() -> NumCosmo.HIPertGravTensor
    """

    h: list[float] = ...
    def dup(self) -> HIPertGravTensor: ...
    def free(self) -> None: ...
    @classmethod
    def new(cls) -> HIPertGravTensor: ...
    def set_zero(self) -> None: ...

class HIPertGravVector(GObject.GBoxed):
    r"""
    :Constructors:

    ::

        HIPertGravVector()
        new() -> NumCosmo.HIPertGravVector
    """

    dsigma: list[float] = ...
    def dup(self) -> HIPertGravVector: ...
    def free(self) -> None: ...
    @classmethod
    def new(cls) -> HIPertGravVector: ...
    def set_zero(self) -> None: ...

class HIPertIAdiab(GObject.GInterface):
    r"""
    Interface NcHIPertIAdiab

    Signals from GObject:
      notify (GParam)
    """

    def eval_F1(self, tau: float, k: float) -> float: ...
    def eval_hubble(self, tau: float) -> float: ...
    def eval_lapse(self, tau: float) -> float: ...
    def eval_m(self, tau: float, k: float) -> float: ...
    def eval_nu(self, tau: float, k: float) -> float: ...
    def eval_p2Psi(self, tau: float, k: float) -> float: ...
    def eval_p2drho(self, tau: float, k: float) -> float: ...
    def eval_tau_hubble(self, k: float) -> float: ...
    def eval_tau_jeans(self, k: float) -> float: ...
    def eval_unit(self) -> float: ...
    def eval_x(self, tau: float) -> float: ...
    def eval_xi(self, tau: float, k: float) -> float: ...

class HIPertIAdiabInterface(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPertIAdiabInterface()
    """

    parent: GObject.TypeInterface = ...
    eval_xi: typing.Callable[[HIPertIAdiab, float, float], float] = ...
    eval_F1: typing.Callable[[HIPertIAdiab, float, float], float] = ...
    eval_nu: typing.Callable[[HIPertIAdiab, float, float], float] = ...
    eval_m: typing.Callable[[HIPertIAdiab, float, float], float] = ...
    eval_unit: typing.Callable[[HIPertIAdiab], float] = ...
    eval_x: typing.Callable[[HIPertIAdiab, float], float] = ...
    eval_p2Psi: typing.Callable[[HIPertIAdiab, float, float], float] = ...
    eval_p2drho: typing.Callable[[HIPertIAdiab, float, float], float] = ...
    eval_lapse: typing.Callable[[HIPertIAdiab, float], float] = ...
    eval_tau_hubble: typing.Callable[[HIPertIAdiab, float], float] = ...
    eval_tau_jeans: typing.Callable[[HIPertIAdiab, float], float] = ...
    eval_hubble: typing.Callable[[HIPertIAdiab, float], float] = ...
    padding: list[None] = ...

class HIPertIEM(GObject.GInterface):
    r"""
    Interface NcHIPertIEM

    Signals from GObject:
      notify (GParam)
    """

    def eval_F1(self, tau: float, k: float) -> float: ...
    def eval_m(self, tau: float, k: float) -> float: ...
    def eval_nu(self, tau: float, k: float) -> float: ...
    def eval_unit(self) -> float: ...
    def eval_x(self, tau: float) -> float: ...
    def eval_xi(self, tau: float, k: float) -> float: ...

class HIPertIEMInterface(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPertIEMInterface()
    """

    parent: GObject.TypeInterface = ...
    eval_xi: typing.Callable[[HIPertIEM, float, float], float] = ...
    eval_F1: typing.Callable[[HIPertIEM, float, float], float] = ...
    eval_nu: typing.Callable[[HIPertIEM, float, float], float] = ...
    eval_m: typing.Callable[[HIPertIEM, float, float], float] = ...
    eval_unit: typing.Callable[[HIPertIEM], float] = ...
    eval_x: typing.Callable[[HIPertIEM, float], float] = ...
    eval_lapse: typing.Callable[[HIPertIEM, float], float] = ...
    padding: list[None] = ...

class HIPertIGW(GObject.GInterface):
    r"""
    Interface NcHIPertIGW

    Signals from GObject:
      notify (GParam)
    """

    def eval_F1(self, tau: float, k: float) -> float: ...
    def eval_m(self, tau: float, k: float) -> float: ...
    def eval_nu(self, tau: float, k: float) -> float: ...
    def eval_unit(self) -> float: ...
    def eval_x(self, tau: float) -> float: ...
    def eval_xi(self, tau: float, k: float) -> float: ...

class HIPertIGWInterface(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPertIGWInterface()
    """

    parent: GObject.TypeInterface = ...
    eval_xi: typing.Callable[[HIPertIGW, float, float], float] = ...
    eval_F1: typing.Callable[[HIPertIGW, float, float], float] = ...
    eval_nu: typing.Callable[[HIPertIGW, float, float], float] = ...
    eval_m: typing.Callable[[HIPertIGW, float, float], float] = ...
    eval_unit: typing.Callable[[HIPertIGW], float] = ...
    eval_x: typing.Callable[[HIPertIGW, float], float] = ...
    padding: list[None] = ...

class HIPertITwoFluids(GObject.GInterface):
    r"""
    Interface NcHIPertITwoFluids

    Signals from GObject:
      notify (GParam)
    """

    def eom_eval(self, alpha: float, k: float) -> HIPertITwoFluidsEOM: ...
    def tv_eval(self, alpha: float, k: float) -> HIPertITwoFluidsTV: ...

class HIPertITwoFluidsEOM(GObject.GBoxed):
    r"""
    :Constructors:

    ::

        HIPertITwoFluidsEOM()
    """

    skey: int = ...
    alpha: float = ...
    k: float = ...
    nu1: float = ...
    nu2: float = ...
    gammabar11: float = ...
    gammabar22: float = ...
    gammabar12: float = ...
    taubar: float = ...
    m_zeta: float = ...
    m_s: float = ...
    mnu2_zeta: float = ...
    mnu2_s: float = ...
    y: float = ...
    sin2phi: float = ...
    cos2phi: float = ...
    cs2: float = ...
    cm2: float = ...
    def dup(self) -> HIPertITwoFluidsEOM: ...
    def free(self) -> None: ...

class HIPertITwoFluidsInterface(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPertITwoFluidsInterface()
    """

    parent: GObject.TypeInterface = ...
    eom: typing.Callable[[HIPertITwoFluids, float, float], HIPertITwoFluidsEOM] = ...
    tv: typing.Callable[[HIPertITwoFluids, float, float], HIPertITwoFluidsTV] = ...

class HIPertITwoFluidsTV(GObject.GBoxed):
    r"""
    :Constructors:

    ::

        HIPertITwoFluidsTV()
    """

    skey: int = ...
    alpha: float = ...
    k: float = ...
    zeta: list[float] = ...
    s: list[float] = ...
    Pzeta: list[float] = ...
    Ps: list[float] = ...
    def dup(self) -> HIPertITwoFluidsTV: ...
    def free(self) -> None: ...

class HIPertPrivate(GObject.GPointer): ...

class HIPertTwoFluids(HIPert):
    r"""
    :Constructors:

    ::

        HIPertTwoFluids(**properties)
        new() -> NumCosmo.HIPertTwoFluids

    Object NcHIPertTwoFluids

    Properties from NcHIPert:
      k -> gdouble: k
        Mode k
      sys-size -> guint: sys-size
        System size
      reltol -> gdouble: reltol
        Relative tolerance
      abstol -> gdouble: abstol
        Absolute tolerance
      alphai -> gdouble: alphai
        Initial time

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        abstol: float
        alphai: float
        k: float
        reltol: float
        sys_size: int

    props: Props = ...
    parent_instance: HIPert = ...
    priv: HIPertTwoFluidsPrivate = ...
    def __init__(
        self,
        abstol: float = ...,
        alphai: float = ...,
        k: float = ...,
        reltol: float = ...,
        sys_size: int = ...,
    ) -> None: ...
    @staticmethod
    def clear(ptf: HIPertTwoFluids) -> None: ...
    def compute_zeta_spectrum(
        self,
        cosmo: HICosmo,
        mode: int,
        alpha_i: float,
        alpha: float,
        ki: float,
        kf: float,
        nnodes: int,
    ) -> NumCosmoMath.Spline: ...
    def eom(self, cosmo: HICosmo, alpha: float) -> HIPertITwoFluidsEOM: ...
    def evolve(self, cosmo: HICosmo, alphaf: float) -> None: ...
    def evolve_array(
        self, cosmo: HICosmo, alphaf: float, step_reltol: float, step_abstol: float
    ) -> NumCosmoMath.Matrix: ...
    def evolve_mode1sub(self, cosmo: HICosmo, alphaf: float) -> None: ...
    def free(self) -> None: ...
    def get_cross_time(
        self, cosmo: HICosmo, cross: HIPertTwoFluidsCross, alpha_i: float, prec: float
    ) -> float: ...
    def get_init_cond_QP(
        self,
        cosmo: HICosmo,
        alpha: float,
        main_mode: int,
        beta_R: float,
        init_cond: NumCosmoMath.Vector,
    ) -> None: ...
    def get_init_cond_zetaS(
        self,
        cosmo: HICosmo,
        alpha: float,
        main_mode: int,
        beta_R: float,
        init_cond: NumCosmoMath.Vector,
    ) -> None: ...
    def get_state_mod(self) -> float: ...
    @classmethod
    def new(cls) -> HIPertTwoFluids: ...
    def peek_state(
        self, cosmo: HICosmo
    ) -> typing.Tuple[NumCosmoMath.Vector, float]: ...
    def ref(self) -> HIPertTwoFluids: ...
    def set_init_cond(
        self,
        cosmo: HICosmo,
        alpha: float,
        main_mode: int,
        useQP: bool,
        init_cond: NumCosmoMath.Vector,
    ) -> None: ...
    def set_init_cond_mode1sub(
        self, cosmo: HICosmo, alpha: float, init_cond: NumCosmoMath.Vector
    ) -> None: ...
    def to_zeta_s(
        self, cosmo: HICosmo, alpha: float, state: NumCosmoMath.Vector
    ) -> None: ...

class HIPertTwoFluidsClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPertTwoFluidsClass()
    """

    parent_class: HIPertClass = ...

class HIPertTwoFluidsPrivate(GObject.GPointer): ...

class HIPertWKB(HIPert):
    r"""
    :Constructors:

    ::

        HIPertWKB(**properties)

    Object NcHIPertWKB

    Properties from NcHIPert:
      k -> gdouble: k
        Mode k
      sys-size -> guint: sys-size
        System size
      reltol -> gdouble: reltol
        Relative tolerance
      abstol -> gdouble: abstol
        Absolute tolerance
      alphai -> gdouble: alphai
        Initial time

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        abstol: float
        alphai: float
        k: float
        reltol: float
        sys_size: int

    props: Props = ...
    parent_instance: HIPert = ...
    nuA: NumCosmoMath.Spline = ...
    lnF: NumCosmoMath.Spline = ...
    dlnF: NumCosmoMath.Spline = ...
    alpha_phase: float = ...
    cur_phase: float = ...
    alpha_i: float = ...
    alpha_f: float = ...
    alpha_p: float = ...
    def __init__(
        self,
        abstol: float = ...,
        alphai: float = ...,
        k: float = ...,
        reltol: float = ...,
        sys_size: int = ...,
    ) -> None: ...
    @staticmethod
    def clear(wkb: HIPertWKB) -> None: ...
    def do_get_dVnu2(
        self, model: NumCosmoMath.Model, alpha: float, k: float
    ) -> float: ...
    def do_get_m(self, model: NumCosmoMath.Model, alpha: float, k: float) -> float: ...
    def do_get_mnu_dmnu(
        self, model: NumCosmoMath.Model, alpha: float, k: float
    ) -> typing.Tuple[float, float]: ...
    def do_get_nu2(
        self, model: NumCosmoMath.Model, alpha: float, k: float
    ) -> float: ...
    def do_get_nu_V(
        self, model: NumCosmoMath.Model, alpha: float, k: float
    ) -> typing.Tuple[float, float]: ...
    def free(self) -> None: ...
    def get_dVnu2(self, model: NumCosmoMath.Model, alpha: float, k: float) -> float: ...
    def get_m(self, model: NumCosmoMath.Model, alpha: float, k: float) -> float: ...
    def get_mnu_dmnu(
        self, model: NumCosmoMath.Model, alpha: float, k: float
    ) -> typing.Tuple[float, float]: ...
    def get_nu2(self, model: NumCosmoMath.Model, alpha: float, k: float) -> float: ...
    def get_nu_V(
        self, model: NumCosmoMath.Model, alpha: float, k: float
    ) -> typing.Tuple[float, float]: ...
    def maxtime(
        self, model: NumCosmoMath.Model, alpha0: float, alpha1: float
    ) -> float: ...
    def maxtime_prec(
        self, model: NumCosmoMath.Model, cmp: HIPertWKBCmp, alpha0: float, alpha1: float
    ) -> float: ...
    def phase(self, model: NumCosmoMath.Model, alpha: float) -> float: ...
    def prepare(self, model: NumCosmoMath.Model) -> None: ...
    def q(
        self, model: NumCosmoMath.Model, alpha: float
    ) -> typing.Tuple[float, float]: ...
    def q_p(
        self, model: NumCosmoMath.Model, alpha: float
    ) -> typing.Tuple[float, float, float, float]: ...
    def ref(self) -> HIPertWKB: ...
    def set_interval(self, alpha_i: float, alpha_f: float) -> None: ...

class HIPertWKBClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPertWKBClass()
    """

    parent_class: HIPertClass = ...
    get_nu_V: typing.Callable[
        [HIPertWKB, NumCosmoMath.Model, float, float], typing.Tuple[float, float]
    ] = ...
    get_mnu_dmnu: typing.Callable[
        [HIPertWKB, NumCosmoMath.Model, float, float], typing.Tuple[float, float]
    ] = ...
    get_m: typing.Callable[[HIPertWKB, NumCosmoMath.Model, float, float], float] = ...
    get_nu2: typing.Callable[[HIPertWKB, NumCosmoMath.Model, float, float], float] = ...
    get_dVnu2: typing.Callable[[HIPertWKB, NumCosmoMath.Model, float, float], float] = (
        ...
    )

class HIPertWKBQgrwZeta(HIPertWKB):
    r"""
    :Constructors:

    ::

        HIPertWKBQgrwZeta(**properties)

    Object NcHIPertWKBQgrwZeta

    Properties from NcHIPert:
      k -> gdouble: k
        Mode k
      sys-size -> guint: sys-size
        System size
      reltol -> gdouble: reltol
        Relative tolerance
      abstol -> gdouble: abstol
        Absolute tolerance
      alphai -> gdouble: alphai
        Initial time

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        abstol: float
        alphai: float
        k: float
        reltol: float
        sys_size: int

    props: Props = ...
    def __init__(
        self,
        abstol: float = ...,
        alphai: float = ...,
        k: float = ...,
        reltol: float = ...,
        sys_size: int = ...,
    ) -> None: ...

class HIPrim(NumCosmoMath.Model):
    r"""
    :Constructors:

    ::

        HIPrim(**properties)

    Object NcHIPrim

    Properties from NcHIPrim:
      k-pivot -> gdouble: k-pivot
        Pivotal value of k

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        k_pivot: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: NumCosmoMath.Model = ...
    k_pivot: float = ...
    lnk_pivot: float = ...
    def __init__(
        self,
        k_pivot: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    def SA_Ampl(self) -> float: ...
    def SA_powspec_k(self, k: float) -> float: ...
    def T_Ampl(self) -> float: ...
    def T_SA_ratio(self) -> float: ...
    def T_powspec_k(self, k: float) -> float: ...
    @staticmethod
    def clear(prim: HIPrim) -> None: ...
    def do_lnSA_powspec_lnk(self, lnk: float) -> float: ...
    def do_lnT_powspec_lnk(self, lnk: float) -> float: ...
    def do_testee(self, x: float) -> float: ...
    def free(self) -> None: ...
    def get_k_pivot(self) -> float: ...
    def get_lnk_pivot(self) -> float: ...
    @staticmethod
    def id() -> int: ...
    def lnSA_powspec_lnk(self, lnk: float) -> float: ...
    def lnT_powspec_lnk(self, lnk: float) -> float: ...
    @staticmethod
    def log_all_models(parent: typing.Type[typing.Any]) -> None: ...
    def ref(self) -> HIPrim: ...
    def set_k_pivot(self, k_pivot: float) -> None: ...

class HIPrimAtan(HIPrim):
    r"""
    :Constructors:

    ::

        HIPrimAtan(**properties)
        new() -> NumCosmo.HIPrimAtan

    Object NcHIPrimAtan

    Properties from NcHIPrimAtan:
      ln10e10ASA -> gdouble: ln10e10ASA
        \log(10^{10}A_{\mathrm{SA}})
      n-SA -> gdouble: n-SA
        n_{\mathrm{SA}}
      lnkc -> gdouble: lnkc
        \ln(k_\mathrm{c})
      c2 -> gdouble: c2
        c_2
      c3 -> gdouble: c3
        c_3
      lambda -> gdouble: lambda
        \lambda
      T-SA-ratio -> gdouble: T-SA-ratio
        A_T/A_{\mathrm{SA}}
      n-T -> gdouble: n-T
        n_{\mathrm{T}}
      ln10e10ASA-fit -> gboolean: ln10e10ASA-fit
        \log(10^{10}A_{\mathrm{SA}}):fit
      n-SA-fit -> gboolean: n-SA-fit
        n_{\mathrm{SA}}:fit
      lnkc-fit -> gboolean: lnkc-fit
        \ln(k_\mathrm{c}):fit
      c2-fit -> gboolean: c2-fit
        c_2:fit
      c3-fit -> gboolean: c3-fit
        c_3:fit
      lambda-fit -> gboolean: lambda-fit
        \lambda:fit
      T-SA-ratio-fit -> gboolean: T-SA-ratio-fit
        A_T/A_{\mathrm{SA}}:fit
      n-T-fit -> gboolean: n-T-fit
        n_{\mathrm{T}}:fit

    Properties from NcHIPrim:
      k-pivot -> gdouble: k-pivot
        Pivotal value of k

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        T_SA_ratio: float
        T_SA_ratio_fit: bool
        c2: float
        c2_fit: bool
        c3: float
        c3_fit: bool
        lambda_fit: bool
        ln10e10ASA: float
        ln10e10ASA_fit: bool
        lnkc: float
        lnkc_fit: bool
        n_SA: float
        n_SA_fit: bool
        n_T: float
        n_T_fit: bool
        k_pivot: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: HIPrim = ...
    def __init__(
        self,
        T_SA_ratio: float = ...,
        T_SA_ratio_fit: bool = ...,
        c2: float = ...,
        c2_fit: bool = ...,
        c3: float = ...,
        c3_fit: bool = ...,
        lambda_fit: bool = ...,
        ln10e10ASA: float = ...,
        ln10e10ASA_fit: bool = ...,
        lnkc: float = ...,
        lnkc_fit: bool = ...,
        n_SA: float = ...,
        n_SA_fit: bool = ...,
        n_T: float = ...,
        n_T_fit: bool = ...,
        k_pivot: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @classmethod
    def new(cls) -> HIPrimAtan: ...

class HIPrimAtanClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPrimAtanClass()
    """

    parent_class: HIPrimClass = ...

class HIPrimBPL(HIPrim):
    r"""
    :Constructors:

    ::

        HIPrimBPL(**properties)
        new() -> NumCosmo.HIPrimBPL

    Object NcHIPrimBPL

    Properties from NcHIPrimBPL:
      ln10e10ASA -> gdouble: ln10e10ASA
        \log(10^{10}A_{\mathrm{SA}})
      n-SA -> gdouble: n-SA
        n_{\mathrm{SA}}
      delta -> gdouble: delta
        \delta
      lnkb -> gdouble: lnkb
        \ln(k_\mathrm{b})
      T-SA-ratio -> gdouble: T-SA-ratio
        A_T/A_{\mathrm{SA}}
      n-T -> gdouble: n-T
        n_{\mathrm{T}}
      ln10e10ASA-fit -> gboolean: ln10e10ASA-fit
        \log(10^{10}A_{\mathrm{SA}}):fit
      n-SA-fit -> gboolean: n-SA-fit
        n_{\mathrm{SA}}:fit
      delta-fit -> gboolean: delta-fit
        \delta:fit
      lnkb-fit -> gboolean: lnkb-fit
        \ln(k_\mathrm{b}):fit
      T-SA-ratio-fit -> gboolean: T-SA-ratio-fit
        A_T/A_{\mathrm{SA}}:fit
      n-T-fit -> gboolean: n-T-fit
        n_{\mathrm{T}}:fit

    Properties from NcHIPrim:
      k-pivot -> gdouble: k-pivot
        Pivotal value of k

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        T_SA_ratio: float
        T_SA_ratio_fit: bool
        delta: float
        delta_fit: bool
        ln10e10ASA: float
        ln10e10ASA_fit: bool
        lnkb: float
        lnkb_fit: bool
        n_SA: float
        n_SA_fit: bool
        n_T: float
        n_T_fit: bool
        k_pivot: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: HIPrim = ...
    def __init__(
        self,
        T_SA_ratio: float = ...,
        T_SA_ratio_fit: bool = ...,
        delta: float = ...,
        delta_fit: bool = ...,
        ln10e10ASA: float = ...,
        ln10e10ASA_fit: bool = ...,
        lnkb: float = ...,
        lnkb_fit: bool = ...,
        n_SA: float = ...,
        n_SA_fit: bool = ...,
        n_T: float = ...,
        n_T_fit: bool = ...,
        k_pivot: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @classmethod
    def new(cls) -> HIPrimBPL: ...

class HIPrimBPLClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPrimBPLClass()
    """

    parent_class: HIPrimClass = ...

class HIPrimClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPrimClass()
    """

    parent_class: NumCosmoMath.ModelClass = ...
    lnSA_powspec_lnk: typing.Callable[[HIPrim, float], float] = ...
    lnT_powspec_lnk: typing.Callable[[HIPrim, float], float] = ...
    testee: typing.Callable[[HIPrim, float], float] = ...

class HIPrimExpc(HIPrim):
    r"""
    :Constructors:

    ::

        HIPrimExpc(**properties)
        new() -> NumCosmo.HIPrimExpc

    Object NcHIPrimExpc

    Properties from NcHIPrimExpc:
      ln10e10ASA -> gdouble: ln10e10ASA
        \log(10^{10}A_{\mathrm{SA}})
      n-SA -> gdouble: n-SA
        n_{\mathrm{SA}}
      lambdac -> gdouble: lambdac
        \lambda_\mathrm{c}
      lnkc -> gdouble: lnkc
        \ln(k_\mathrm{c})
      c -> gdouble: c
        c
      T-SA-ratio -> gdouble: T-SA-ratio
        A_T/A_{\mathrm{SA}}
      n-T -> gdouble: n-T
        n_{\mathrm{T}}
      ln10e10ASA-fit -> gboolean: ln10e10ASA-fit
        \log(10^{10}A_{\mathrm{SA}}):fit
      n-SA-fit -> gboolean: n-SA-fit
        n_{\mathrm{SA}}:fit
      lambdac-fit -> gboolean: lambdac-fit
        \lambda_\mathrm{c}:fit
      lnkc-fit -> gboolean: lnkc-fit
        \ln(k_\mathrm{c}):fit
      c-fit -> gboolean: c-fit
        c:fit
      T-SA-ratio-fit -> gboolean: T-SA-ratio-fit
        A_T/A_{\mathrm{SA}}:fit
      n-T-fit -> gboolean: n-T-fit
        n_{\mathrm{T}}:fit

    Properties from NcHIPrim:
      k-pivot -> gdouble: k-pivot
        Pivotal value of k

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        T_SA_ratio: float
        T_SA_ratio_fit: bool
        c: float
        c_fit: bool
        lambdac: float
        lambdac_fit: bool
        ln10e10ASA: float
        ln10e10ASA_fit: bool
        lnkc: float
        lnkc_fit: bool
        n_SA: float
        n_SA_fit: bool
        n_T: float
        n_T_fit: bool
        k_pivot: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: HIPrim = ...
    def __init__(
        self,
        T_SA_ratio: float = ...,
        T_SA_ratio_fit: bool = ...,
        c: float = ...,
        c_fit: bool = ...,
        lambdac: float = ...,
        lambdac_fit: bool = ...,
        ln10e10ASA: float = ...,
        ln10e10ASA_fit: bool = ...,
        lnkc: float = ...,
        lnkc_fit: bool = ...,
        n_SA: float = ...,
        n_SA_fit: bool = ...,
        n_T: float = ...,
        n_T_fit: bool = ...,
        k_pivot: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @classmethod
    def new(cls) -> HIPrimExpc: ...

class HIPrimExpcClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPrimExpcClass()
    """

    parent_class: HIPrimClass = ...

class HIPrimPowerLaw(HIPrim):
    r"""
    :Constructors:

    ::

        HIPrimPowerLaw(**properties)
        new() -> NumCosmo.HIPrimPowerLaw

    Object NcHIPrimPowerLaw

    Properties from NcHIPrimPowerLaw:
      ln10e10ASA -> gdouble: ln10e10ASA
        \log(10^{10}A_{\mathrm{SA}})
      T-SA-ratio -> gdouble: T-SA-ratio
        A_T/A_{\mathrm{SA}}
      n-SA -> gdouble: n-SA
        n_{\mathrm{SA}}
      n-T -> gdouble: n-T
        n_{\mathrm{T}}
      ln10e10ASA-fit -> gboolean: ln10e10ASA-fit
        \log(10^{10}A_{\mathrm{SA}}):fit
      T-SA-ratio-fit -> gboolean: T-SA-ratio-fit
        A_T/A_{\mathrm{SA}}:fit
      n-SA-fit -> gboolean: n-SA-fit
        n_{\mathrm{SA}}:fit
      n-T-fit -> gboolean: n-T-fit
        n_{\mathrm{T}}:fit

    Properties from NcHIPrim:
      k-pivot -> gdouble: k-pivot
        Pivotal value of k

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        T_SA_ratio: float
        T_SA_ratio_fit: bool
        ln10e10ASA: float
        ln10e10ASA_fit: bool
        n_SA: float
        n_SA_fit: bool
        n_T: float
        n_T_fit: bool
        k_pivot: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: HIPrim = ...
    def __init__(
        self,
        T_SA_ratio: float = ...,
        T_SA_ratio_fit: bool = ...,
        ln10e10ASA: float = ...,
        ln10e10ASA_fit: bool = ...,
        n_SA: float = ...,
        n_SA_fit: bool = ...,
        n_T: float = ...,
        n_T_fit: bool = ...,
        k_pivot: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @classmethod
    def new(cls) -> HIPrimPowerLaw: ...

class HIPrimPowerLawClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPrimPowerLawClass()
    """

    parent_class: HIPrimClass = ...

class HIPrimSBPL(HIPrim):
    r"""
    :Constructors:

    ::

        HIPrimSBPL(**properties)
        new() -> NumCosmo.HIPrimSBPL

    Object NcHIPrimSBPL

    Properties from NcHIPrimSBPL:
      ln10e10ASA -> gdouble: ln10e10ASA
        \log(10^{10}A_{SA})
      n-SA -> gdouble: n-SA
        n_{\mathrm{SA}}
      delta -> gdouble: delta
        \delta
      RA -> gdouble: RA
        R_\mathrm{A}
      lnkb -> gdouble: lnkb
        \ln(k_\mathrm{b})
      lambda -> gdouble: lambda
        \lambda
      T-SA-ratio -> gdouble: T-SA-ratio
        A_T/A_{\mathrm{SA}}
      n-T -> gdouble: n-T
        n_{\mathrm{T}}
      ln10e10ASA-fit -> gboolean: ln10e10ASA-fit
        \log(10^{10}A_{SA}):fit
      n-SA-fit -> gboolean: n-SA-fit
        n_{\mathrm{SA}}:fit
      delta-fit -> gboolean: delta-fit
        \delta:fit
      RA-fit -> gboolean: RA-fit
        R_\mathrm{A}:fit
      lnkb-fit -> gboolean: lnkb-fit
        \ln(k_\mathrm{b}):fit
      lambda-fit -> gboolean: lambda-fit
        \lambda:fit
      T-SA-ratio-fit -> gboolean: T-SA-ratio-fit
        A_T/A_{\mathrm{SA}}:fit
      n-T-fit -> gboolean: n-T-fit
        n_{\mathrm{T}}:fit

    Properties from NcHIPrim:
      k-pivot -> gdouble: k-pivot
        Pivotal value of k

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        RA: float
        RA_fit: bool
        T_SA_ratio: float
        T_SA_ratio_fit: bool
        delta: float
        delta_fit: bool
        lambda_fit: bool
        ln10e10ASA: float
        ln10e10ASA_fit: bool
        lnkb: float
        lnkb_fit: bool
        n_SA: float
        n_SA_fit: bool
        n_T: float
        n_T_fit: bool
        k_pivot: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: HIPrim = ...
    def __init__(
        self,
        RA: float = ...,
        RA_fit: bool = ...,
        T_SA_ratio: float = ...,
        T_SA_ratio_fit: bool = ...,
        delta: float = ...,
        delta_fit: bool = ...,
        lambda_fit: bool = ...,
        ln10e10ASA: float = ...,
        ln10e10ASA_fit: bool = ...,
        lnkb: float = ...,
        lnkb_fit: bool = ...,
        n_SA: float = ...,
        n_SA_fit: bool = ...,
        n_T: float = ...,
        n_T_fit: bool = ...,
        k_pivot: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @classmethod
    def new(cls) -> HIPrimSBPL: ...

class HIPrimSBPLClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPrimSBPLClass()
    """

    parent_class: HIPrimClass = ...

class HIPrimTwoFluids(HIPrim):
    r"""
    :Constructors:

    ::

        HIPrimTwoFluids(**properties)
        new() -> NumCosmo.HIPrimTwoFluids

    Object NcHIPrimTwoFluids

    Properties from NcHIPrimTwoFluids:
      lnk-lnw-spline -> NcmSpline2d: lnk-lnw-spline
        Spline for the primordial adiabatic scalar power spectrum as a function of ln(k) and ln(w)
      use-default-calib -> gboolean: use-default-calib
        Use default calibration
      ln10e10ASA -> gdouble: ln10e10ASA
        \log(10^{10}A_{\mathrm{SA}})
      T-SA-ratio -> gdouble: T-SA-ratio
        A_T/A_{\mathrm{SA}}
      lnk0 -> gdouble: lnk0
        \ln(k_0)
      lnw -> gdouble: lnw
        \ln(w)
      n-T -> gdouble: n-T
        n_{\mathrm{T}}
      ln10e10ASA-fit -> gboolean: ln10e10ASA-fit
        \log(10^{10}A_{\mathrm{SA}}):fit
      T-SA-ratio-fit -> gboolean: T-SA-ratio-fit
        A_T/A_{\mathrm{SA}}:fit
      lnk0-fit -> gboolean: lnk0-fit
        \ln(k_0):fit
      lnw-fit -> gboolean: lnw-fit
        \ln(w):fit
      n-T-fit -> gboolean: n-T-fit
        n_{\mathrm{T}}:fit

    Properties from NcHIPrim:
      k-pivot -> gdouble: k-pivot
        Pivotal value of k

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        T_SA_ratio: float
        T_SA_ratio_fit: bool
        ln10e10ASA: float
        ln10e10ASA_fit: bool
        lnk_lnw_spline: NumCosmoMath.Spline2d
        lnk0: float
        lnk0_fit: bool
        lnw: float
        lnw_fit: bool
        n_T: float
        n_T_fit: bool
        use_default_calib: bool
        k_pivot: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: HIPrim = ...
    def __init__(
        self,
        T_SA_ratio: float = ...,
        T_SA_ratio_fit: bool = ...,
        ln10e10ASA: float = ...,
        ln10e10ASA_fit: bool = ...,
        lnk_lnw_spline: NumCosmoMath.Spline2d = ...,
        lnk0: float = ...,
        lnk0_fit: bool = ...,
        lnw: float = ...,
        lnw_fit: bool = ...,
        n_T: float = ...,
        n_T_fit: bool = ...,
        use_default_calib: bool = ...,
        k_pivot: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    def get_use_default_calib(self) -> bool: ...
    @classmethod
    def new(cls) -> HIPrimTwoFluids: ...
    def peek_lnk_lnw_spline(self) -> NumCosmoMath.Spline2d: ...
    def set_lnk_lnw_spline(
        self, lnSA_powspec_lnk_lnw: NumCosmoMath.Spline2d
    ) -> None: ...
    def set_use_default_calib(self, use_default_calib: bool) -> None: ...

class HIPrimTwoFluidsClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIPrimTwoFluidsClass()
    """

    parent_class: HIPrimClass = ...

class HIQG1D(GObject.Object):
    r"""
    :Constructors:

    ::

        HIQG1D(**properties)
        new() -> NumCosmo.HIQG1D
        new_full(nknots:int, lambda_:float) -> NumCosmo.HIQG1D

    Object NcHIQG1D

    Properties from NcHIQG1D:
      lambda -> gdouble: lambda
        \lambda
      abstol -> gdouble: abstol
        absolute tolerance
      reltol -> gdouble: reltol
        relative tolerance
      nknots -> guint: nknots
        n_k
      noboundary -> gboolean: noboundary
        no boundary condition at x_f

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        abstol: float
        nknots: int
        noboundary: bool
        reltol: float

    props: Props = ...
    parent_instance: GObject.Object = ...
    priv: HIQG1DPrivate = ...
    def __init__(
        self,
        abstol: float = ...,
        nknots: int = ...,
        noboundary: bool = ...,
        reltol: float = ...,
    ) -> None: ...
    def Bohm(self, i: int) -> float: ...
    def Bohm_p(self, i: int) -> float: ...
    def Hbasis(self, x: float, y: float, h: float, a: float) -> float: ...
    def Sbasis_x3(
        self, x: float, y1: float, y2: float, h: float, a: float
    ) -> float: ...
    def basis(self, x: float, y: float, h: float, a: float) -> float: ...
    @staticmethod
    def clear(qg1d: HIQG1D) -> None: ...
    def eval_dS(self, x: float) -> float: ...
    def eval_ev(self, i: int, x: float) -> float: ...
    def eval_psi(self, x: float) -> list[float]: ...
    def eval_psi0(self, x: float) -> list[float]: ...
    def evol(self, t: float) -> None: ...
    def expect_d(self) -> float: ...
    def expect_p(self) -> float: ...
    def free(self) -> None: ...
    def get_acs_a(self) -> float: ...
    def get_basis_a(self) -> float: ...
    def get_lambda(self) -> float: ...
    def get_mu(self) -> float: ...
    def get_nknots(self) -> int: ...
    def get_nu(self) -> float: ...
    def int_rho_0_inf(self) -> float: ...
    def int_x2rho_0_inf(self) -> float: ...
    def int_xrho_0_inf(self) -> float: ...
    def nBohm(self) -> int: ...
    @classmethod
    def new(cls) -> HIQG1D: ...
    @classmethod
    def new_full(cls, nknots: int, lambda_: float) -> HIQG1D: ...
    def peek_knots(self) -> NumCosmoMath.Vector: ...
    def prepare(self) -> None: ...
    def ref(self) -> HIQG1D: ...
    def set_init_cond(
        self,
        psi0_lnRS: typing.Callable[..., list[float]],
        xi: float,
        xf: float,
        *psi_data: typing.Any,
    ) -> None: ...
    def set_init_cond_exp(self, qm_exp: HIQG1DExp, xi: float, xf: float) -> None: ...
    def set_init_cond_gauss(
        self, qm_gauss: HIQG1DGauss, xi: float, xf: float
    ) -> None: ...
    def set_init_cond_sq(self, qm_sq: HIQG1DSQ, xi: float, xf: float) -> None: ...
    def set_nknots(self, nknots: int) -> None: ...

class HIQG1DClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIQG1DClass()
    """

    parent_class: GObject.ObjectClass = ...

class HIQG1DExp(GObject.GBoxed):
    r"""
    :Constructors:

    ::

        HIQG1DExp()
        new(n:float, V:float, pV:float) -> NumCosmo.HIQG1DExp
    """

    n: float = ...
    V: float = ...
    pV: float = ...
    lnNorm: float = ...
    def dup(self) -> HIQG1DExp: ...
    def eval(self, x: float) -> list[float]: ...
    def eval_lnRS(self, x: float) -> list[float]: ...
    def free(self) -> None: ...
    @classmethod
    def new(cls, n: float, V: float, pV: float) -> HIQG1DExp: ...

class HIQG1DGauss(GObject.GBoxed):
    r"""
    :Constructors:

    ::

        HIQG1DGauss()
        new(mean:float, alpha:float, sigma:float, Hi:float) -> NumCosmo.HIQG1DGauss
    """

    mean: float = ...
    alpha: float = ...
    sigma: float = ...
    Hi: float = ...
    lnNorm: float = ...
    def dup(self) -> HIQG1DGauss: ...
    def eval(self, x: float) -> list[float]: ...
    def eval_hermit(self, x: float) -> list[float]: ...
    def eval_lnRS(self, x: float) -> list[float]: ...
    def free(self) -> None: ...
    @classmethod
    def new(cls, mean: float, alpha: float, sigma: float, Hi: float) -> HIQG1DGauss: ...

class HIQG1DPrivate(GObject.GPointer): ...

class HIQG1DSQ(GObject.GBoxed):
    r"""
    :Constructors:

    ::

        HIQG1DSQ()
        new(mu:float, V:float, pV:float) -> NumCosmo.HIQG1DSQ
    """

    mu: float = ...
    V: float = ...
    pV: float = ...
    lnNorm: float = ...
    def dup(self) -> HIQG1DSQ: ...
    def eval(self, x: float) -> list[float]: ...
    def eval_lnRS(self, x: float) -> list[float]: ...
    def free(self) -> None: ...
    @classmethod
    def new(cls, mu: float, V: float, pV: float) -> HIQG1DSQ: ...

class HIReion(NumCosmoMath.Model):
    r"""
    :Constructors:

    ::

        HIReion(**properties)

    Object NcHIReion

    Properties from NcHIReion:
      prec -> gdouble: prec
        Precision for reionization calculations

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        prec: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: NumCosmoMath.Model = ...
    prec: float = ...
    def __init__(
        self,
        prec: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def clear(reion: HIReion) -> None: ...
    def do_get_Xe(self, cosmo: HICosmo, lambda_: float, Xe_recomb: float) -> float: ...
    def do_get_init_x(self, cosmo: HICosmo) -> float: ...
    def do_get_tau(self, cosmo: HICosmo) -> float: ...
    def free(self) -> None: ...
    def get_Xe(self, cosmo: HICosmo, lambda_: float, Xe_recomb: float) -> float: ...
    def get_init_x(self, cosmo: HICosmo) -> float: ...
    def get_tau(self, cosmo: HICosmo) -> float: ...
    @staticmethod
    def id() -> int: ...
    def ref(self) -> HIReion: ...

class HIReionCamb(HIReion):
    r"""
    :Constructors:

    ::

        HIReionCamb(**properties)
        new() -> NumCosmo.HIReionCamb

    Object NcHIReionCamb

    Properties from NcHIReionCamb:
      HII-HeII-reion-delta -> gdouble: HII-HeII-reion-delta
        Window size for HII and HeII reionization
      HeIII-reion-delta -> gdouble: HeIII-reion-delta
        Window size for HeIII reionization
      HII-HeII-reion-exponent -> gdouble: HII-HeII-reion-exponent
        Exponent for HII and HeII reionization transition
      HeII-reionized -> gboolean: HeII-reionized
        Whether HeIII is reionized
      z-re -> gdouble: z-re
        z_\mathrm{re}
      z-He-re -> gdouble: z-He-re
        z^\mathrm{He}_\mathrm{re}
      z-re-fit -> gboolean: z-re-fit
        z_\mathrm{re}:fit
      z-He-re-fit -> gboolean: z-He-re-fit
        z^\mathrm{He}_\mathrm{re}:fit

    Properties from NcHIReion:
      prec -> gdouble: prec
        Precision for reionization calculations

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        HII_HeII_reion_delta: float
        HII_HeII_reion_exponent: float
        HeII_reionized: bool
        HeIII_reion_delta: float
        z_He_re: float
        z_He_re_fit: bool
        z_re: float
        z_re_fit: bool
        prec: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: HIReion = ...
    HII_HeII_reion_delta: float = ...
    HeIII_reion_delta: float = ...
    HII_HeII_reion_expo: float = ...
    HII_HeII_reion_delta_eff: float = ...
    HII_HeII_reion_x_pow_expo: float = ...
    HEII_reionized: bool = ...
    fsol: int = ...
    tau_ctrl: NumCosmoMath.ModelCtrl = ...
    def __init__(
        self,
        HII_HeII_reion_delta: float = ...,
        HII_HeII_reion_exponent: float = ...,
        HeII_reionized: bool = ...,
        HeIII_reion_delta: float = ...,
        z_He_re: float = ...,
        z_He_re_fit: bool = ...,
        z_re: float = ...,
        z_re_fit: bool = ...,
        prec: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    def calc_z_from_tau(self, cosmo: HICosmo, tau: float) -> float: ...
    @classmethod
    def new(cls) -> HIReionCamb: ...
    def set_z_from_tau(self, cosmo: HICosmo, tau: float) -> None: ...
    def z_to_tau(self, cosmo: HICosmo) -> None: ...

class HIReionCambClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIReionCambClass()
    """

    parent_class: HIReionClass = ...

class HIReionCambReparamTau(NumCosmoMath.Reparam):
    r"""
    :Constructors:

    ::

        HIReionCambReparamTau(**properties)
        new(length:int, cosmo:NumCosmo.HICosmo) -> NumCosmo.HIReionCambReparamTau

    Object NcHIReionCambReparamTau

    Properties from NcHIReionCambReparamTau:
      cosmo -> NcHICosmo: cosmo
        Cosmological model used to transform tau <=> z

    Properties from NcmReparam:
      length -> guint: length
        System's length
      params-desc -> NcmObjDictInt: params-desc
        News parameter descriptions
      compat-type -> gchararray: compat-type
        Compatible type

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        cosmo: HICosmo
        compat_type: str
        length: int
        params_desc: NumCosmoMath.ObjDictInt

    props: Props = ...
    parent_instance: NumCosmoMath.Reparam = ...
    ctrl: NumCosmoMath.ModelCtrl = ...
    def __init__(
        self,
        cosmo: HICosmo = ...,
        compat_type: str = ...,
        length: int = ...,
        params_desc: NumCosmoMath.ObjDictInt = ...,
    ) -> None: ...
    @classmethod
    def new(cls, length: int, cosmo: HICosmo) -> HIReionCambReparamTau: ...

class HIReionCambReparamTauClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIReionCambReparamTauClass()
    """

    parent_class: NumCosmoMath.ReparamClass = ...

class HIReionClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HIReionClass()
    """

    parent_class: NumCosmoMath.ModelClass = ...
    get_init_x: typing.Callable[[HIReion, HICosmo], float] = ...
    get_Xe: typing.Callable[[HIReion, HICosmo, float, float], float] = ...
    get_tau: typing.Callable[[HIReion, HICosmo], float] = ...

class HaloBias(GObject.Object):
    r"""
    :Constructors:

    ::

        HaloBias(**properties)

    Object NcHaloBias

    Properties from NcHaloBias:
      mass-function -> NcHaloMassFunction: mass-function
        Mass Function.

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        mass_function: HaloMassFunction

    props: Props = ...
    parent_instance: GObject.Object = ...
    mfp: HaloMassFunction = ...
    def __init__(self, mass_function: HaloMassFunction = ...) -> None: ...
    @staticmethod
    def clear(bias: HaloBias) -> None: ...
    def do_eval(self, cosmo: HICosmo, sigma: float, z: float) -> float: ...
    def eval(self, cosmo: HICosmo, sigma: float, z: float) -> float: ...
    def free(self) -> None: ...
    def integrand(self, cosmo: HICosmo, lnM: float, z: float) -> float: ...

class HaloBiasClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HaloBiasClass()
    """

    parent_class: GObject.ObjectClass = ...
    eval: typing.Callable[[HaloBias, HICosmo, float, float], float] = ...

class HaloBiasPS(HaloBias):
    r"""
    :Constructors:

    ::

        HaloBiasPS(**properties)
        new(mfp:NumCosmo.HaloMassFunction) -> NumCosmo.HaloBiasPS
        new_full(mfp:NumCosmo.HaloMassFunction, delta_c:float) -> NumCosmo.HaloBiasPS

    Object NcHaloBiasPS

    Properties from NcHaloBiasPS:
      critical-delta -> gdouble: critical-delta
        Critical delta

    Properties from NcHaloBias:
      mass-function -> NcHaloMassFunction: mass-function
        Mass Function.

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        critical_delta: float
        mass_function: HaloMassFunction

    props: Props = ...
    parent_instance: HaloBias = ...
    delta_c: float = ...
    def __init__(
        self, critical_delta: float = ..., mass_function: HaloMassFunction = ...
    ) -> None: ...
    def get_delta_c(self) -> float: ...
    @classmethod
    def new(cls, mfp: HaloMassFunction) -> HaloBiasPS: ...
    @classmethod
    def new_full(cls, mfp: HaloMassFunction, delta_c: float) -> HaloBiasPS: ...
    def set_delta_c(self, delta_c: float) -> None: ...

class HaloBiasPSClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HaloBiasPSClass()
    """

    parent_class: HaloBiasClass = ...

class HaloBiasSTEllip(HaloBias):
    r"""
    :Constructors:

    ::

        HaloBiasSTEllip(**properties)
        new(mfp:NumCosmo.HaloMassFunction) -> NumCosmo.HaloBiasSTEllip
        new_full(mfp:NumCosmo.HaloMassFunction, delta_c:float, a:float, b:float, c:float) -> NumCosmo.HaloBiasSTEllip

    Object NcHaloBiasSTEllip

    Properties from NcHaloBiasSTEllip:
      critical-delta -> gdouble: critical-delta
        Critical delta
      a -> gdouble: a
        a
      b -> gdouble: b
        b
      c -> gdouble: c
        c

    Properties from NcHaloBias:
      mass-function -> NcHaloMassFunction: mass-function
        Mass Function.

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        a: float
        b: float
        c: float
        critical_delta: float
        mass_function: HaloMassFunction

    props: Props = ...
    parent_instance: HaloBias = ...
    delta_c: float = ...
    a: float = ...
    b: float = ...
    c: float = ...
    def __init__(
        self,
        a: float = ...,
        b: float = ...,
        c: float = ...,
        critical_delta: float = ...,
        mass_function: HaloMassFunction = ...,
    ) -> None: ...
    def get_a(self) -> float: ...
    def get_b(self) -> float: ...
    def get_c(self) -> float: ...
    def get_delta_c(self) -> float: ...
    @classmethod
    def new(cls, mfp: HaloMassFunction) -> HaloBiasSTEllip: ...
    @classmethod
    def new_full(
        cls, mfp: HaloMassFunction, delta_c: float, a: float, b: float, c: float
    ) -> HaloBiasSTEllip: ...
    def set_a(self, a: float) -> None: ...
    def set_b(self, b: float) -> None: ...
    def set_c(self, c: float) -> None: ...
    def set_delta_c(self, delta_c: float) -> None: ...

class HaloBiasSTEllipClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HaloBiasSTEllipClass()
    """

    parent_class: HaloBiasClass = ...

class HaloBiasSTSpher(HaloBias):
    r"""
    :Constructors:

    ::

        HaloBiasSTSpher(**properties)
        new(mfp:NumCosmo.HaloMassFunction) -> NumCosmo.HaloBiasSTSpher
        new_full(mfp:NumCosmo.HaloMassFunction, delta_c:float, a:float, p:float) -> NumCosmo.HaloBiasSTSpher

    Object NcHaloBiasSTSpher

    Properties from NcHaloBiasSTSpher:
      critical-delta -> gdouble: critical-delta
        Critical delta
      a -> gdouble: a
        a
      p -> gdouble: p
        p

    Properties from NcHaloBias:
      mass-function -> NcHaloMassFunction: mass-function
        Mass Function.

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        a: float
        critical_delta: float
        p: float
        mass_function: HaloMassFunction

    props: Props = ...
    parent_instance: HaloBias = ...
    delta_c: float = ...
    a: float = ...
    p: float = ...
    def __init__(
        self,
        a: float = ...,
        critical_delta: float = ...,
        p: float = ...,
        mass_function: HaloMassFunction = ...,
    ) -> None: ...
    def get_a(self) -> float: ...
    def get_delta_c(self) -> float: ...
    def get_p(self) -> float: ...
    @classmethod
    def new(cls, mfp: HaloMassFunction) -> HaloBiasSTSpher: ...
    @classmethod
    def new_full(
        cls, mfp: HaloMassFunction, delta_c: float, a: float, p: float
    ) -> HaloBiasSTSpher: ...
    def set_a(self, a: float) -> None: ...
    def set_delta_c(self, delta_c: float) -> None: ...
    def set_p(self, p: float) -> None: ...

class HaloBiasSTSpherClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HaloBiasSTSpherClass()
    """

    parent_class: HaloBiasClass = ...

class HaloBiasTinker(HaloBias):
    r"""
    :Constructors:

    ::

        HaloBiasTinker(**properties)
        new(mfp:NumCosmo.HaloMassFunction) -> NumCosmo.HaloBiasTinker
        new_full(mfp:NumCosmo.HaloMassFunction, delta_c:float, B:float, b:float, c:float) -> NumCosmo.HaloBiasTinker

    Object NcHaloBiasTinker

    Properties from NcHaloBiasTinker:
      critical-delta -> gdouble: critical-delta
        Critical delta
      B -> gdouble: B
        B
      b -> gdouble: b
        b
      c -> gdouble: c
        c

    Properties from NcHaloBias:
      mass-function -> NcHaloMassFunction: mass-function
        Mass Function.

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        B: float
        b: float
        c: float
        critical_delta: float
        mass_function: HaloMassFunction

    props: Props = ...
    parent_instance: HaloBias = ...
    delta_c: float = ...
    B: float = ...
    b: float = ...
    c: float = ...
    def __init__(
        self,
        B: float = ...,
        b: float = ...,
        c: float = ...,
        critical_delta: float = ...,
        mass_function: HaloMassFunction = ...,
    ) -> None: ...
    def get_B(self) -> float: ...
    def get_b(self) -> float: ...
    def get_c(self) -> float: ...
    def get_delta_c(self) -> float: ...
    @classmethod
    def new(cls, mfp: HaloMassFunction) -> HaloBiasTinker: ...
    @classmethod
    def new_full(
        cls, mfp: HaloMassFunction, delta_c: float, B: float, b: float, c: float
    ) -> HaloBiasTinker: ...
    def set_B(self, B: float) -> None: ...
    def set_b(self, b: float) -> None: ...
    def set_c(self, c: float) -> None: ...
    def set_delta_c(self, delta_c: float) -> None: ...

class HaloBiasTinkerClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HaloBiasTinkerClass()
    """

    parent_class: HaloBiasClass = ...

class HaloCMDuffy08(HaloMassSummary):
    r"""
    :Constructors:

    ::

        HaloCMDuffy08(**properties)
        new(mdef:NumCosmo.HaloMassSummaryMassDef, Delta:float) -> NumCosmo.HaloCMDuffy08

    Object NcHaloCMDuffy08

    Properties from NcHaloCMDuffy08:
      log10MDelta -> gdouble: log10MDelta
        \log_{10}(M_{\Delta})
      log10MDelta-fit -> gboolean: log10MDelta-fit
        \log_{10}(M_{\Delta}):fit

    Properties from NcHaloMassSummary:
      mass-def -> NcHaloMassSummaryMassDef: mass-def
        Mass definition
      Delta -> gdouble: Delta
        Overdensity constant

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        log10MDelta: float
        log10MDelta_fit: bool
        Delta: float
        mass_def: HaloMassSummaryMassDef
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    def __init__(
        self,
        log10MDelta: float = ...,
        log10MDelta_fit: bool = ...,
        Delta: float = ...,
        mass_def: HaloMassSummaryMassDef = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def clear(hcmd: HaloCMDuffy08) -> None: ...
    def free(self) -> None: ...
    @classmethod
    def new(cls, mdef: HaloMassSummaryMassDef, Delta: float) -> HaloCMDuffy08: ...
    def ref(self) -> HaloCMDuffy08: ...

class HaloCMDuffy08Class(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HaloCMDuffy08Class()
    """

    parent_class: HaloMassSummaryClass = ...

class HaloCMKlypin11(HaloMassSummary):
    r"""
    :Constructors:

    ::

        HaloCMKlypin11(**properties)
        new(mdef:NumCosmo.HaloMassSummaryMassDef, Delta:float) -> NumCosmo.HaloCMKlypin11

    Object NcHaloCMKlypin11

    Properties from NcHaloCMKlypin11:
      log10MDelta -> gdouble: log10MDelta
        \log_{10}(M_{\Delta})
      log10MDelta-fit -> gboolean: log10MDelta-fit
        \log_{10}(M_{\Delta}):fit

    Properties from NcHaloMassSummary:
      mass-def -> NcHaloMassSummaryMassDef: mass-def
        Mass definition
      Delta -> gdouble: Delta
        Overdensity constant

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        log10MDelta: float
        log10MDelta_fit: bool
        Delta: float
        mass_def: HaloMassSummaryMassDef
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    def __init__(
        self,
        log10MDelta: float = ...,
        log10MDelta_fit: bool = ...,
        Delta: float = ...,
        mass_def: HaloMassSummaryMassDef = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def clear(hcmk: HaloCMKlypin11) -> None: ...
    def free(self) -> None: ...
    @classmethod
    def new(cls, mdef: HaloMassSummaryMassDef, Delta: float) -> HaloCMKlypin11: ...
    def ref(self) -> HaloCMKlypin11: ...

class HaloCMKlypin11Class(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HaloCMKlypin11Class()
    """

    parent_class: HaloMassSummaryClass = ...

class HaloCMParam(HaloMassSummary):
    r"""
    :Constructors:

    ::

        HaloCMParam(**properties)
        new(mdef:NumCosmo.HaloMassSummaryMassDef, Delta:float) -> NumCosmo.HaloCMParam

    Object NcHaloCMParam

    Properties from NcHaloCMParam:
      cDelta -> gdouble: cDelta
        c_{\Delta}
      log10MDelta -> gdouble: log10MDelta
        \log_{10}(M_{\Delta})
      cDelta-fit -> gboolean: cDelta-fit
        c_{\Delta}:fit
      log10MDelta-fit -> gboolean: log10MDelta-fit
        \log_{10}(M_{\Delta}):fit

    Properties from NcHaloMassSummary:
      mass-def -> NcHaloMassSummaryMassDef: mass-def
        Mass definition
      Delta -> gdouble: Delta
        Overdensity constant

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        cDelta: float
        cDelta_fit: bool
        log10MDelta: float
        log10MDelta_fit: bool
        Delta: float
        mass_def: HaloMassSummaryMassDef
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    def __init__(
        self,
        cDelta: float = ...,
        cDelta_fit: bool = ...,
        log10MDelta: float = ...,
        log10MDelta_fit: bool = ...,
        Delta: float = ...,
        mass_def: HaloMassSummaryMassDef = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def clear(hcmp: HaloCMParam) -> None: ...
    def free(self) -> None: ...
    @classmethod
    def new(cls, mdef: HaloMassSummaryMassDef, Delta: float) -> HaloCMParam: ...
    def ref(self) -> HaloCMParam: ...

class HaloCMParamClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HaloCMParamClass()
    """

    parent_class: HaloMassSummaryClass = ...

class HaloDensityProfile(NumCosmoMath.Model):
    r"""
    :Constructors:

    ::

        HaloDensityProfile(**properties)

    Object NcHaloDensityProfile

    Properties from NcHaloDensityProfile:
      reltol -> gdouble: reltol
        Relative tolerance
      lnXi -> gdouble: lnXi
        Computation interval lower limit
      lnXf -> gdouble: lnXf
        Computation interval upper limit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        lnXf: float
        lnXi: float
        reltol: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: NumCosmoMath.Model = ...
    def __init__(
        self,
        lnXf: float = ...,
        lnXi: float = ...,
        reltol: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def clear(dp: HaloDensityProfile) -> None: ...
    def do_eval_dl_2d_density(self, X: float) -> float: ...
    def do_eval_dl_cyl_mass(self, X: float) -> float: ...
    def do_eval_dl_density(self, x: float) -> float: ...
    def do_eval_dl_spher_mass(self, x: float) -> float: ...
    def eval_2d_density(self, cosmo: HICosmo, R: float, z: float) -> float: ...
    def eval_2d_density_array(
        self,
        cosmo: HICosmo,
        R: typing.Sequence[float],
        fin: float,
        fout: float,
        z: float,
    ) -> list[float]: ...
    def eval_cyl_mass(self, cosmo: HICosmo, R: float, z: float) -> float: ...
    def eval_cyl_mass_array(
        self,
        cosmo: HICosmo,
        R: typing.Sequence[float],
        fin: float,
        fout: float,
        z: float,
    ) -> list[float]: ...
    def eval_density(self, cosmo: HICosmo, r: float, z: float) -> float: ...
    def eval_density_array(
        self,
        cosmo: HICosmo,
        r: typing.Sequence[float],
        fin: float,
        fout: float,
        z: float,
    ) -> list[float]: ...
    def eval_dl_2d_density(self, X: float) -> float: ...
    def eval_dl_cyl_mass(self, X: float) -> float: ...
    def eval_dl_density(self, x: float) -> float: ...
    def eval_dl_spher_mass(self, x: float) -> float: ...
    def eval_numint_dl_2d_density(self, X: float) -> float: ...
    def eval_numint_dl_cyl_mass(self, X: float) -> float: ...
    def eval_numint_dl_spher_mass(self, X: float) -> float: ...
    def eval_spher_mass(self, cosmo: HICosmo, r: float, z: float) -> float: ...
    def eval_spher_mass_delta(self, cosmo: HICosmo, z: float) -> float: ...
    def free(self) -> None: ...
    def get_lnXf(self) -> float: ...
    def get_lnXi(self) -> float: ...
    def get_numint_splines(
        self,
    ) -> typing.Tuple[
        NumCosmoMath.Spline, NumCosmoMath.Spline, NumCosmoMath.Spline
    ]: ...
    def get_phys_limts(
        self, cosmo: HICosmo, z: float
    ) -> typing.Tuple[float, float]: ...
    def get_reltol(self) -> float: ...
    @staticmethod
    def id() -> int: ...
    def peek_mass_summary(self) -> HaloMassSummary: ...
    def r_s(self, cosmo: HICosmo, z: float) -> float: ...
    def r_s_rho_s(self, cosmo: HICosmo, z: float) -> typing.Tuple[float, float]: ...
    def ref(self) -> HaloDensityProfile: ...
    def rho_s(self, cosmo: HICosmo, z: float) -> float: ...
    def set_lnXf(self, lnXf: float) -> None: ...
    def set_lnXi(self, lnXi: float) -> None: ...
    def set_reltol(self, reltol: float) -> None: ...

class HaloDensityProfileClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HaloDensityProfileClass()
    """

    parent_class: NumCosmoMath.ModelClass = ...
    eval_dl_density: typing.Callable[[HaloDensityProfile, float], float] = ...
    eval_dl_spher_mass: typing.Callable[[HaloDensityProfile, float], float] = ...
    eval_dl_2d_density: typing.Callable[[HaloDensityProfile, float], float] = ...
    eval_dl_cyl_mass: typing.Callable[[HaloDensityProfile, float], float] = ...

class HaloDensityProfileDK14(HaloDensityProfile):
    r"""
    :Constructors:

    ::

        HaloDensityProfileDK14(**properties)
        new() -> NumCosmo.HaloDensityProfile

    Object NcHaloDensityProfileDK14

    Properties from NcHaloDensityProfileDK14:
      Delta -> gdouble: Delta
        Overdensity constant
      rt -> gdouble: rt
        r_{t}
      beta -> gdouble: beta
        \beta
      rt-fit -> gboolean: rt-fit
        r_{t}:fit
      beta-fit -> gboolean: beta-fit
        \beta:fit

    Properties from NcHaloDensityProfile:
      reltol -> gdouble: reltol
        Relative tolerance
      lnXi -> gdouble: lnXi
        Computation interval lower limit
      lnXf -> gdouble: lnXf
        Computation interval upper limit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        Delta: float
        beta: float
        beta_fit: bool
        rt: float
        rt_fit: bool
        lnXf: float
        lnXi: float
        reltol: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: HaloDensityProfile = ...
    Delta: float = ...
    r_Delta: float = ...
    def __init__(
        self,
        Delta: float = ...,
        beta: float = ...,
        beta_fit: bool = ...,
        rt: float = ...,
        rt_fit: bool = ...,
        lnXf: float = ...,
        lnXi: float = ...,
        reltol: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @classmethod
    def new(cls) -> HaloDensityProfileDK14: ...

class HaloDensityProfileDK14Class(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HaloDensityProfileDK14Class()
    """

    parent_class: HaloDensityProfileClass = ...

class HaloDensityProfileEinasto(HaloDensityProfile):
    r"""
    :Constructors:

    ::

        HaloDensityProfileEinasto(**properties)
        new(hms:NumCosmo.HaloMassSummary) -> NumCosmo.HaloDensityProfileEinasto

    Object NcHaloDensityProfileEinasto

    Properties from NcHaloDensityProfileEinasto:
      alpha -> gdouble: alpha
        \alpha
      alpha-fit -> gboolean: alpha-fit
        \alpha:fit

    Properties from NcHaloDensityProfile:
      reltol -> gdouble: reltol
        Relative tolerance
      lnXi -> gdouble: lnXi
        Computation interval lower limit
      lnXf -> gdouble: lnXf
        Computation interval upper limit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        alpha: float
        alpha_fit: bool
        lnXf: float
        lnXi: float
        reltol: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    def __init__(
        self,
        alpha: float = ...,
        alpha_fit: bool = ...,
        lnXf: float = ...,
        lnXi: float = ...,
        reltol: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @classmethod
    def new(cls, hms: HaloMassSummary) -> HaloDensityProfileEinasto: ...

class HaloDensityProfileEinastoClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HaloDensityProfileEinastoClass()
    """

    parent_class: HaloDensityProfileClass = ...

class HaloDensityProfileHernquist(HaloDensityProfile):
    r"""
    :Constructors:

    ::

        HaloDensityProfileHernquist(**properties)
        new(hms:NumCosmo.HaloMassSummary) -> NumCosmo.HaloDensityProfileHernquist

    Object NcHaloDensityProfileHernquist

    Properties from NcHaloDensityProfile:
      reltol -> gdouble: reltol
        Relative tolerance
      lnXi -> gdouble: lnXi
        Computation interval lower limit
      lnXf -> gdouble: lnXf
        Computation interval upper limit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        lnXf: float
        lnXi: float
        reltol: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    def __init__(
        self,
        lnXf: float = ...,
        lnXi: float = ...,
        reltol: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @classmethod
    def new(cls, hms: HaloMassSummary) -> HaloDensityProfileHernquist: ...

class HaloDensityProfileHernquistClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HaloDensityProfileHernquistClass()
    """

    parent_class: HaloDensityProfileClass = ...

class HaloDensityProfileNFW(HaloDensityProfile):
    r"""
    :Constructors:

    ::

        HaloDensityProfileNFW(**properties)
        new(hms:NumCosmo.HaloMassSummary) -> NumCosmo.HaloDensityProfileNFW

    Object NcHaloDensityProfileNFW

    Properties from NcHaloDensityProfile:
      reltol -> gdouble: reltol
        Relative tolerance
      lnXi -> gdouble: lnXi
        Computation interval lower limit
      lnXf -> gdouble: lnXf
        Computation interval upper limit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        lnXf: float
        lnXi: float
        reltol: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    def __init__(
        self,
        lnXf: float = ...,
        lnXi: float = ...,
        reltol: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @classmethod
    def new(cls, hms: HaloMassSummary) -> HaloDensityProfileNFW: ...
    def set_ni(*args, **kwargs): ...  # FIXME Function

class HaloDensityProfileNFWClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HaloDensityProfileNFWClass()
    """

    parent_class: HaloDensityProfileClass = ...
    @staticmethod
    def set_ni(num: bool) -> None: ...

class HaloMassFunction(GObject.Object):
    r"""
    :Constructors:

    ::

        HaloMassFunction(**properties)
        new(dist:NumCosmo.Distance, psf:NumCosmoMath.PowspecFilter, mulf:NumCosmo.MultiplicityFunc) -> NumCosmo.HaloMassFunction

    Object NcHaloMassFunction

    Properties from NcHaloMassFunction:
      distance -> NcDistance: distance
        Distance
      powerspectrum-filtered -> NcmPowspecFilter: powerspectrum-filtered
        Filtered power-spectrum
      multiplicity -> NcMultiplicityFunc: multiplicity
        Multiplicity function
      area -> gdouble: area
        Angular area in steradian
      prec -> gdouble: prec
        Precision
      lnMi -> gdouble: lnMi
        Lower mass
      lnMf -> gdouble: lnMf
        Upper mass
      zi -> gdouble: zi
        Lower redshift
      zf -> gdouble: zf
        Upper redshift
      mf-lb -> gdouble: mf-lb
        Upper redshift

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        area: float
        distance: Distance
        lnMf: float
        lnMi: float
        mf_lb: float
        multiplicity: MultiplicityFunc
        powerspectrum_filtered: NumCosmoMath.PowspecFilter
        prec: float
        zf: float
        zi: float

    props: Props = ...
    parent_instance: GObject.Object = ...
    priv: HaloMassFunctionPrivate = ...
    d2NdzdlnM: NumCosmoMath.Spline2d = ...
    mf_lb: float = ...
    def __init__(
        self,
        area: float = ...,
        distance: Distance = ...,
        lnMf: float = ...,
        lnMi: float = ...,
        mf_lb: float = ...,
        multiplicity: MultiplicityFunc = ...,
        powerspectrum_filtered: NumCosmoMath.PowspecFilter = ...,
        prec: float = ...,
        zf: float = ...,
        zi: float = ...,
    ) -> None: ...
    @staticmethod
    def clear(mfp: HaloMassFunction) -> None: ...
    def d2n_dzdlnM(self, cosmo: HICosmo, lnM: float, z: float) -> float: ...
    def dn_dlnM(self, cosmo: HICosmo, lnM: float, z: float) -> float: ...
    def dn_dlnR(self, cosmo: HICosmo, lnR: float, z: float) -> float: ...
    def dn_dz(
        self, cosmo: HICosmo, lnMl: float, lnMu: float, z: float, spline: bool
    ) -> float: ...
    def dv_dzdomega(self, cosmo: HICosmo, z: float) -> float: ...
    def free(self) -> None: ...
    def lnM_to_lnR(self, cosmo: HICosmo, lnM: float) -> float: ...
    def lnR_to_lnM(self, cosmo: HICosmo, lnR: float) -> float: ...
    def n(
        self,
        cosmo: HICosmo,
        lnMl: float,
        lnMu: float,
        zl: float,
        zu: float,
        spline: HaloMassFunctionSplineOptimize,
    ) -> float: ...
    @classmethod
    def new(
        cls, dist: Distance, psf: NumCosmoMath.PowspecFilter, mulf: MultiplicityFunc
    ) -> HaloMassFunction: ...
    def peek_multiplicity_function(self) -> MultiplicityFunc: ...
    def peek_psf(self) -> NumCosmoMath.PowspecFilter: ...
    def peek_survey_area(self) -> float: ...
    def prepare(self, cosmo: HICosmo) -> None: ...
    def prepare_if_needed(self, cosmo: HICosmo) -> None: ...
    def set_area(self, area: float) -> None: ...
    def set_area_sd(self, area_sd: float) -> None: ...
    def set_eval_limits(
        self, cosmo: HICosmo, lnMi: float, lnMf: float, zi: float, zf: float
    ) -> None: ...
    def set_prec(self, prec: float) -> None: ...
    def sigma_lnM(self, cosmo: HICosmo, lnM: float, z: float) -> float: ...
    def sigma_lnR(self, cosmo: HICosmo, lnR: float, z: float) -> float: ...

class HaloMassFunctionClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HaloMassFunctionClass()
    """

    parent_class: GObject.ObjectClass = ...

class HaloMassFunctionPrivate(GObject.GPointer): ...

class HaloMassSummary(NumCosmoMath.Model):
    r"""
    :Constructors:

    ::

        HaloMassSummary(**properties)

    Object NcHaloMassSummary

    Properties from NcHaloMassSummary:
      mass-def -> NcHaloMassSummaryMassDef: mass-def
        Mass definition
      Delta -> gdouble: Delta
        Overdensity constant

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        Delta: float
        mass_def: HaloMassSummaryMassDef
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: NumCosmoMath.Model = ...
    def __init__(
        self,
        Delta: float = ...,
        mass_def: HaloMassSummaryMassDef = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    def Delta(self, cosmo: HICosmo, z: float) -> float: ...
    def Delta_rho_bg(self, cosmo: HICosmo, z: float) -> float: ...
    @staticmethod
    def clear(hms: HaloMassSummary) -> None: ...
    def concentration(self, cosmo: HICosmo, z: float) -> float: ...
    def do_concentration(self, cosmo: HICosmo, z: float) -> float: ...
    def do_mass(self) -> float: ...
    def do_set_Delta(self, Delta: float) -> None: ...
    def do_set_mdef(self, mdef: HaloMassSummaryMassDef) -> None: ...
    def free(self) -> None: ...
    def get_Delta(self) -> float: ...
    def get_mdef(self) -> HaloMassSummaryMassDef: ...
    @staticmethod
    def id() -> int: ...
    def mass(self) -> float: ...
    def ref(self) -> HaloMassSummary: ...
    def rho_bg(self, cosmo: HICosmo, z: float) -> float: ...

class HaloMassSummaryClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HaloMassSummaryClass()
    """

    parent_class: NumCosmoMath.ModelClass = ...
    mass: typing.Callable[[HaloMassSummary], float] = ...
    concentration: typing.Callable[[HaloMassSummary, HICosmo, float], float] = ...
    set_Delta: typing.Callable[[HaloMassSummary, float], None] = ...
    set_mdef: typing.Callable[[HaloMassSummary, HaloMassSummaryMassDef], None] = ...
    padding: list[None] = ...

class HaloPosition(NumCosmoMath.Model):
    r"""
    :Constructors:

    ::

        HaloPosition(**properties)
        new(dist:NumCosmo.Distance) -> NumCosmo.HaloPosition

    Object NcHaloPosition

    Properties from NcHaloPosition:
      dist -> NcDistance: Distance
        A NcDistance object
      ra -> gdouble: ra
        \mathrm{RA}
      dec -> gdouble: dec
        \mathrm{Dec}
      z -> gdouble: z
        z
      ra-fit -> gboolean: ra-fit
        \mathrm{RA}:fit
      dec-fit -> gboolean: dec-fit
        \mathrm{Dec}:fit
      z-fit -> gboolean: z-fit
        z:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        dec: float
        dec_fit: bool
        dist: Distance
        ra: float
        ra_fit: bool
        z: float
        z_fit: bool
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    def __init__(
        self,
        dec: float = ...,
        dec_fit: bool = ...,
        dist: Distance = ...,
        ra: float = ...,
        ra_fit: bool = ...,
        z: float = ...,
        z_fit: bool = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def clear(hp: HaloPosition) -> None: ...
    def free(self) -> None: ...
    def get_ra_dec(self) -> typing.Tuple[float, float]: ...
    def get_redshift(self) -> float: ...
    @staticmethod
    def id() -> int: ...
    @classmethod
    def new(cls, dist: Distance) -> HaloPosition: ...
    def polar_angles(self, ra: float, dec: float) -> typing.Tuple[float, float]: ...
    def prepare(self, cosmo: HICosmo) -> None: ...
    def prepare_if_needed(self, cosmo: HICosmo) -> None: ...
    def projected_radius(self, cosmo: HICosmo, theta: float) -> float: ...
    def projected_radius_from_ra_dec(
        self, cosmo: HICosmo, ra: float, dec: float
    ) -> float: ...
    def ref(self) -> HaloPosition: ...

class HaloPositionClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        HaloPositionClass()
    """

    parent_class: NumCosmoMath.ModelClass = ...

class MultiplicityFunc(GObject.Object):
    r"""
    :Constructors:

    ::

        MultiplicityFunc(**properties)

    Object NcMultiplicityFunc

    Properties from NcMultiplicityFunc:
      mass-def -> NcMultiplicityFuncMassDef: mass-def
        Mass definition
      Delta -> gdouble: Delta
        Delta

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        Delta: float
        mass_def: MultiplicityFuncMassDef

    props: Props = ...
    parent_instance: GObject.Object = ...
    priv: MultiplicityFuncPrivate = ...
    def __init__(
        self, Delta: float = ..., mass_def: MultiplicityFuncMassDef = ...
    ) -> None: ...
    @staticmethod
    def clear(mulf: MultiplicityFunc) -> None: ...
    def correction_factor(
        self, cosmo: HICosmo, sigma: float, z: float, lnM: float
    ) -> float: ...
    def do_correction_factor(
        self, cosmo: HICosmo, sigma: float, z: float, lnM: float
    ) -> float: ...
    def do_eval(self, cosmo: HICosmo, sigma: float, z: float) -> float: ...
    def do_get_Delta(self) -> float: ...
    def do_get_matter_Delta(self, cosmo: HICosmo, z: float) -> float: ...
    def do_get_mdef(self) -> MultiplicityFuncMassDef: ...
    def do_has_correction_factor(self) -> bool: ...
    def do_set_Delta(self, Delta: float) -> None: ...
    def do_set_mdef(self, mdef: MultiplicityFuncMassDef) -> None: ...
    def eval(self, cosmo: HICosmo, sigma: float, z: float) -> float: ...
    def free(self) -> None: ...
    def get_Delta(self) -> float: ...
    def get_matter_Delta(self, cosmo: HICosmo, z: float) -> float: ...
    def get_mdef(self) -> MultiplicityFuncMassDef: ...
    def has_correction_factor(self) -> bool: ...
    def set_Delta(self, Delta: float) -> None: ...
    def set_mdef(self, mdef: MultiplicityFuncMassDef) -> None: ...

class MultiplicityFuncBocquet(MultiplicityFunc):
    r"""
    :Constructors:

    ::

        MultiplicityFuncBocquet(**properties)
        new() -> NumCosmo.MultiplicityFuncBocquet
        new_full(mdef:NumCosmo.MultiplicityFuncMassDef, sim:NumCosmo.MultiplicityFuncBocquetSim, Delta:float) -> NumCosmo.MultiplicityFuncBocquet

    Object NcMultiplicityFuncBocquet

    Properties from NcMultiplicityFuncBocquet:
      sim -> NcMultiplicityFuncBocquetSim: sim
        Simulation type

    Properties from NcMultiplicityFunc:
      mass-def -> NcMultiplicityFuncMassDef: mass-def
        Mass definition
      Delta -> gdouble: Delta
        Delta

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        sim: MultiplicityFuncBocquetSim
        Delta: float
        mass_def: MultiplicityFuncMassDef

    props: Props = ...
    parent_instance: MultiplicityFunc = ...
    priv: MultiplicityFuncBocquetPrivate = ...
    def __init__(
        self,
        sim: MultiplicityFuncBocquetSim = ...,
        Delta: float = ...,
        mass_def: MultiplicityFuncMassDef = ...,
    ) -> None: ...
    @staticmethod
    def clear(mb: MultiplicityFuncBocquet) -> None: ...
    def free(self) -> None: ...
    def get_sim(self) -> MultiplicityFuncBocquetSim: ...
    @classmethod
    def new(cls) -> MultiplicityFuncBocquet: ...
    @classmethod
    def new_full(
        cls,
        mdef: MultiplicityFuncMassDef,
        sim: MultiplicityFuncBocquetSim,
        Delta: float,
    ) -> MultiplicityFuncBocquet: ...
    def ref(self) -> MultiplicityFuncBocquet: ...
    def set_sim(self, sim: MultiplicityFuncBocquetSim) -> None: ...

class MultiplicityFuncBocquetClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        MultiplicityFuncBocquetClass()
    """

    parent_class: MultiplicityFuncClass = ...

class MultiplicityFuncBocquetPrivate(GObject.GPointer): ...

class MultiplicityFuncClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        MultiplicityFuncClass()
    """

    parent_class: GObject.ObjectClass = ...
    set_mdef: typing.Callable[[MultiplicityFunc, MultiplicityFuncMassDef], None] = ...
    set_Delta: typing.Callable[[MultiplicityFunc, float], None] = ...
    get_Delta: typing.Callable[[MultiplicityFunc], float] = ...
    get_matter_Delta: typing.Callable[[MultiplicityFunc, HICosmo, float], float] = ...
    get_mdef: typing.Callable[[MultiplicityFunc], MultiplicityFuncMassDef] = ...
    eval: typing.Callable[[MultiplicityFunc, HICosmo, float, float], float] = ...
    has_correction_factor: typing.Callable[[MultiplicityFunc], bool] = ...
    correction_factor: typing.Callable[
        [MultiplicityFunc, HICosmo, float, float, float], float
    ] = ...

class MultiplicityFuncCrocce(MultiplicityFunc):
    r"""
    :Constructors:

    ::

        MultiplicityFuncCrocce(**properties)
        new() -> NumCosmo.MultiplicityFuncCrocce

    Object NcMultiplicityFuncCrocce

    Properties from NcMultiplicityFunc:
      mass-def -> NcMultiplicityFuncMassDef: mass-def
        Mass definition
      Delta -> gdouble: Delta
        Delta

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        Delta: float
        mass_def: MultiplicityFuncMassDef

    props: Props = ...
    parent_instance: MultiplicityFunc = ...
    priv: MultiplicityFuncCroccePrivate = ...
    def __init__(
        self, Delta: float = ..., mass_def: MultiplicityFuncMassDef = ...
    ) -> None: ...
    @staticmethod
    def clear(mc: MultiplicityFuncCrocce) -> None: ...
    def free(self) -> None: ...
    @classmethod
    def new(cls) -> MultiplicityFuncCrocce: ...
    def ref(self) -> MultiplicityFuncCrocce: ...

class MultiplicityFuncCrocceClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        MultiplicityFuncCrocceClass()
    """

    parent_class: MultiplicityFuncClass = ...

class MultiplicityFuncCroccePrivate(GObject.GPointer): ...

class MultiplicityFuncDespali(MultiplicityFunc):
    r"""
    :Constructors:

    ::

        MultiplicityFuncDespali(**properties)
        new() -> NumCosmo.MultiplicityFuncDespali
        new_full(mdef:NumCosmo.MultiplicityFuncMassDef, Delta:float) -> NumCosmo.MultiplicityFuncDespali

    Object NcMultiplicityFuncDespali

    Properties from NcMultiplicityFuncDespali:
      E0 -> gboolean: E0
        Whether the halo finder uses elliptical overdensity
      CMF -> gboolean: CMF
        Whether the use of the cluster mass function

    Properties from NcMultiplicityFunc:
      mass-def -> NcMultiplicityFuncMassDef: mass-def
        Mass definition
      Delta -> gdouble: Delta
        Delta

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        CMF: bool
        E0: bool
        Delta: float
        mass_def: MultiplicityFuncMassDef

    props: Props = ...
    parent_instance: MultiplicityFunc = ...
    def __init__(
        self,
        CMF: bool = ...,
        E0: bool = ...,
        Delta: float = ...,
        mass_def: MultiplicityFuncMassDef = ...,
    ) -> None: ...
    @staticmethod
    def clear(md: MultiplicityFuncDespali) -> None: ...
    def delta_c(self, cosmo: HICosmo, z: float) -> float: ...
    def delta_vir(self, cosmo: HICosmo, z: float) -> float: ...
    def free(self) -> None: ...
    def get_cmf(self) -> bool: ...
    def get_eo(self) -> bool: ...
    @classmethod
    def new(cls) -> MultiplicityFuncDespali: ...
    @classmethod
    def new_full(
        cls, mdef: MultiplicityFuncMassDef, Delta: float
    ) -> MultiplicityFuncDespali: ...
    def ref(self) -> MultiplicityFuncDespali: ...
    def set_cmf(self, on: bool) -> None: ...
    def set_eo(self, on: bool) -> None: ...

class MultiplicityFuncDespaliClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        MultiplicityFuncDespaliClass()
    """

    parent_class: MultiplicityFuncClass = ...

class MultiplicityFuncDespaliPrivate(GObject.GPointer): ...

class MultiplicityFuncJenkins(MultiplicityFunc):
    r"""
    :Constructors:

    ::

        MultiplicityFuncJenkins(**properties)
        new() -> NumCosmo.MultiplicityFuncJenkins

    Object NcMultiplicityFuncJenkins

    Properties from NcMultiplicityFunc:
      mass-def -> NcMultiplicityFuncMassDef: mass-def
        Mass definition
      Delta -> gdouble: Delta
        Delta

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        Delta: float
        mass_def: MultiplicityFuncMassDef

    props: Props = ...
    parent_instance: MultiplicityFunc = ...
    priv: MultiplicityFuncJenkinsPrivate = ...
    def __init__(
        self, Delta: float = ..., mass_def: MultiplicityFuncMassDef = ...
    ) -> None: ...
    @staticmethod
    def clear(mj: MultiplicityFuncJenkins) -> None: ...
    def free(self) -> None: ...
    @classmethod
    def new(cls) -> MultiplicityFuncJenkins: ...
    def ref(self) -> MultiplicityFuncJenkins: ...

class MultiplicityFuncJenkinsClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        MultiplicityFuncJenkinsClass()
    """

    parent_class: MultiplicityFuncClass = ...

class MultiplicityFuncJenkinsPrivate(GObject.GPointer): ...

class MultiplicityFuncPS(MultiplicityFunc):
    r"""
    :Constructors:

    ::

        MultiplicityFuncPS(**properties)
        new() -> NumCosmo.MultiplicityFuncPS

    Object NcMultiplicityFuncPS

    Properties from NcMultiplicityFuncPS:
      critical-delta -> gdouble: critical-delta
        Critical delta

    Properties from NcMultiplicityFunc:
      mass-def -> NcMultiplicityFuncMassDef: mass-def
        Mass definition
      Delta -> gdouble: Delta
        Delta

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        critical_delta: float
        Delta: float
        mass_def: MultiplicityFuncMassDef

    props: Props = ...
    parent_instance: MultiplicityFunc = ...
    priv: MultiplicityFuncPSPrivate = ...
    def __init__(
        self,
        critical_delta: float = ...,
        Delta: float = ...,
        mass_def: MultiplicityFuncMassDef = ...,
    ) -> None: ...
    @staticmethod
    def clear(mps: MultiplicityFuncPS) -> None: ...
    def free(self) -> None: ...
    def get_delta_c(self) -> float: ...
    @classmethod
    def new(cls) -> MultiplicityFuncPS: ...
    def ref(self) -> MultiplicityFuncPS: ...
    def set_delta_c(self, delta_c: float) -> None: ...

class MultiplicityFuncPSClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        MultiplicityFuncPSClass()
    """

    parent_class: MultiplicityFuncClass = ...

class MultiplicityFuncPSPrivate(GObject.GPointer): ...
class MultiplicityFuncPrivate(GObject.GPointer): ...

class MultiplicityFuncST(MultiplicityFunc):
    r"""
    :Constructors:

    ::

        MultiplicityFuncST(**properties)
        new() -> NumCosmo.MultiplicityFuncST

    Object NcMultiplicityFuncST

    Properties from NcMultiplicityFuncST:
      A -> gdouble: A
        A
      b -> gdouble: b
        b
      p -> gdouble: p
        p
      critical-delta -> gdouble: critical-delta
        Critical delta

    Properties from NcMultiplicityFunc:
      mass-def -> NcMultiplicityFuncMassDef: mass-def
        Mass definition
      Delta -> gdouble: Delta
        Delta

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        A: float
        b: float
        critical_delta: float
        p: float
        Delta: float
        mass_def: MultiplicityFuncMassDef

    props: Props = ...
    parent_instance: MultiplicityFunc = ...
    priv: MultiplicityFuncSTPrivate = ...
    def __init__(
        self,
        A: float = ...,
        b: float = ...,
        critical_delta: float = ...,
        p: float = ...,
        Delta: float = ...,
        mass_def: MultiplicityFuncMassDef = ...,
    ) -> None: ...
    @staticmethod
    def clear(mst: MultiplicityFuncST) -> None: ...
    def free(self) -> None: ...
    def get_A(self) -> float: ...
    def get_b(self) -> float: ...
    def get_delta_c(self) -> float: ...
    def get_p(self) -> float: ...
    @classmethod
    def new(cls) -> MultiplicityFuncST: ...
    def ref(self) -> MultiplicityFuncST: ...
    def set_A(self, A: float) -> None: ...
    def set_b(self, b: float) -> None: ...
    def set_delta_c(self, delta_c: float) -> None: ...
    def set_p(self, p: float) -> None: ...

class MultiplicityFuncSTClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        MultiplicityFuncSTClass()
    """

    parent_class: MultiplicityFuncClass = ...

class MultiplicityFuncSTPrivate(GObject.GPointer): ...

class MultiplicityFuncTinker(MultiplicityFunc):
    r"""
    :Constructors:

    ::

        MultiplicityFuncTinker(**properties)
        new() -> NumCosmo.MultiplicityFuncTinker
        new_full(mdef:NumCosmo.MultiplicityFuncMassDef, Delta:float) -> NumCosmo.MultiplicityFuncTinker

    Object NcMultiplicityFuncTinker

    Properties from NcMultiplicityFuncTinker:
      linear-interp -> gboolean: linear-interp
        Use linear interpolation

    Properties from NcMultiplicityFunc:
      mass-def -> NcMultiplicityFuncMassDef: mass-def
        Mass definition
      Delta -> gdouble: Delta
        Delta

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        linear_interp: bool
        Delta: float
        mass_def: MultiplicityFuncMassDef

    props: Props = ...
    parent_instance: MultiplicityFunc = ...
    priv: MultiplicityFuncTinkerPrivate = ...
    def __init__(
        self,
        linear_interp: bool = ...,
        Delta: float = ...,
        mass_def: MultiplicityFuncMassDef = ...,
    ) -> None: ...
    @staticmethod
    def clear(mt: MultiplicityFuncTinker) -> None: ...
    def free(self) -> None: ...
    @classmethod
    def new(cls) -> MultiplicityFuncTinker: ...
    @classmethod
    def new_full(
        cls, mdef: MultiplicityFuncMassDef, Delta: float
    ) -> MultiplicityFuncTinker: ...
    def ref(self) -> MultiplicityFuncTinker: ...
    def set_linear_interp(self, lin_interp: bool) -> None: ...

class MultiplicityFuncTinkerClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        MultiplicityFuncTinkerClass()
    """

    parent_class: MultiplicityFuncClass = ...

class MultiplicityFuncTinkerMeanNormalized(MultiplicityFunc):
    r"""
    :Constructors:

    ::

        MultiplicityFuncTinkerMeanNormalized(**properties)
        new() -> NumCosmo.MultiplicityFuncTinkerMeanNormalized

    Object NcMultiplicityFuncTinkerMeanNormalized

    Properties from NcMultiplicityFunc:
      mass-def -> NcMultiplicityFuncMassDef: mass-def
        Mass definition
      Delta -> gdouble: Delta
        Delta

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        Delta: float
        mass_def: MultiplicityFuncMassDef

    props: Props = ...
    parent_instance: MultiplicityFunc = ...
    priv: MultiplicityFuncTinkerMeanNormalizedPrivate = ...
    def __init__(
        self, Delta: float = ..., mass_def: MultiplicityFuncMassDef = ...
    ) -> None: ...
    @staticmethod
    def clear(mt10: MultiplicityFuncTinkerMeanNormalized) -> None: ...
    def free(self) -> None: ...
    @classmethod
    def new(cls) -> MultiplicityFuncTinkerMeanNormalized: ...
    def ref(self) -> MultiplicityFuncTinkerMeanNormalized: ...

class MultiplicityFuncTinkerMeanNormalizedClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        MultiplicityFuncTinkerMeanNormalizedClass()
    """

    parent_class: MultiplicityFuncClass = ...

class MultiplicityFuncTinkerMeanNormalizedPrivate(GObject.GPointer): ...
class MultiplicityFuncTinkerPrivate(GObject.GPointer): ...

class MultiplicityFuncWarren(MultiplicityFunc):
    r"""
    :Constructors:

    ::

        MultiplicityFuncWarren(**properties)
        new() -> NumCosmo.MultiplicityFuncWarren

    Object NcMultiplicityFuncWarren

    Properties from NcMultiplicityFunc:
      mass-def -> NcMultiplicityFuncMassDef: mass-def
        Mass definition
      Delta -> gdouble: Delta
        Delta

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        Delta: float
        mass_def: MultiplicityFuncMassDef

    props: Props = ...
    parent_instance: MultiplicityFunc = ...
    priv: MultiplicityFuncWarrenPrivate = ...
    def __init__(
        self, Delta: float = ..., mass_def: MultiplicityFuncMassDef = ...
    ) -> None: ...
    @staticmethod
    def clear(mw: MultiplicityFuncWarren) -> None: ...
    def free(self) -> None: ...
    @classmethod
    def new(cls) -> MultiplicityFuncWarren: ...
    def ref(self) -> MultiplicityFuncWarren: ...

class MultiplicityFuncWarrenClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        MultiplicityFuncWarrenClass()
    """

    parent_class: MultiplicityFuncClass = ...

class MultiplicityFuncWarrenPrivate(GObject.GPointer): ...

class MultiplicityFuncWatson(MultiplicityFunc):
    r"""
    :Constructors:

    ::

        MultiplicityFuncWatson(**properties)
        new() -> NumCosmo.MultiplicityFuncWatson

    Object NcMultiplicityFuncWatson

    Properties from NcMultiplicityFunc:
      mass-def -> NcMultiplicityFuncMassDef: mass-def
        Mass definition
      Delta -> gdouble: Delta
        Delta

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        Delta: float
        mass_def: MultiplicityFuncMassDef

    props: Props = ...
    parent_instance: MultiplicityFunc = ...
    priv: MultiplicityFuncWatsonPrivate = ...
    def __init__(
        self, Delta: float = ..., mass_def: MultiplicityFuncMassDef = ...
    ) -> None: ...
    @staticmethod
    def clear(mwat: MultiplicityFuncWatson) -> None: ...
    def free(self) -> None: ...
    @classmethod
    def new(cls) -> MultiplicityFuncWatson: ...
    def ref(self) -> MultiplicityFuncWatson: ...

class MultiplicityFuncWatsonClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        MultiplicityFuncWatsonClass()
    """

    parent_class: MultiplicityFuncClass = ...

class MultiplicityFuncWatsonPrivate(GObject.GPointer): ...

class PlanckFI(NumCosmoMath.Model):
    r"""
    :Constructors:

    ::

        PlanckFI(**properties)

    Object NcPlanckFI

    Properties from NcPlanckFI:
      version -> guint: version
        Planck compatible version

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        version: int
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: NumCosmoMath.Model = ...
    version: int = ...
    def __init__(
        self,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def clear(pfi: PlanckFI) -> None: ...
    def free(self) -> None: ...
    @staticmethod
    def id() -> int: ...
    @staticmethod
    def log_all_models() -> None: ...
    def ref(self) -> PlanckFI: ...

class PlanckFIClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        PlanckFIClass()
    """

    parent_class: NumCosmoMath.ModelClass = ...

class PlanckFICorTT(PlanckFI):
    r"""
    :Constructors:

    ::

        PlanckFICorTT(**properties)

    Object NcPlanckFICorTT

    Properties from NcPlanckFICorTT:
      A-cib-217 -> gdouble: A-cib-217
        A^{\mathrm{CIB}}_{217}
      cib-index -> gdouble: cib-index
        n^{\mathrm{CIB}}
      xi-sz-cib -> gdouble: xi-sz-cib
        \xi^{\mathrm{tSZ}\times \mathrm{CIB}}
      A-sz -> gdouble: A-sz
        A^{\mathrm{tSZ}}
      ps-A-100-100 -> gdouble: ps-A-100-100
        A^{\mathrm{PS}}_{100}
      ps-A-143-143 -> gdouble: ps-A-143-143
        A^{\mathrm{PS}}_{143}
      ps-A-143-217 -> gdouble: ps-A-143-217
        A^{\mathrm{PS}}_{143\times 217}
      ps-A-217-217 -> gdouble: ps-A-217-217
        A^{\mathrm{PS}}_{217}
      ksz-norm -> gdouble: ksz-norm
        A^{\mathrm{kSZ}}
      gal545-A-100 -> gdouble: gal545-A-100
        A^{\mathrm{dust}TT}_{100}
      gal545-A-143 -> gdouble: gal545-A-143
        A^{\mathrm{dust}TT}_{143}
      gal545-A-143-217 -> gdouble: gal545-A-143-217
        A^{\mathrm{dust}TT}_{143 \times 217}
      gal545-A-217 -> gdouble: gal545-A-217
        A^{\mathrm{dust}TT}_{217}
      A-sbpx-100-100-TT -> gdouble: A-sbpx-100-100-TT
        A^{\mathrm{sbpx}TT}_{100 \times 100}
      A-sbpx-143-143-TT -> gdouble: A-sbpx-143-143-TT
        A^{\mathrm{sbpx}TT}_{143 \times 143}
      A-sbpx-143-217-TT -> gdouble: A-sbpx-143-217-TT
        A^{\mathrm{sbpx}TT}_{143 \times 217}
      A-sbpx-217-217-TT -> gdouble: A-sbpx-217-217-TT
        A^{\mathrm{sbpx}TT}_{217 \times 217}
      calib-100T -> gdouble: calib-100T
        c_{100}
      calib-217T -> gdouble: calib-217T
        c_{217}
      A-planck -> gdouble: A-planck
        y_{\mathrm{cal}}
      A-cib-217-fit -> gboolean: A-cib-217-fit
        A^{\mathrm{CIB}}_{217}:fit
      cib-index-fit -> gboolean: cib-index-fit
        n^{\mathrm{CIB}}:fit
      xi-sz-cib-fit -> gboolean: xi-sz-cib-fit
        \xi^{\mathrm{tSZ}\times \mathrm{CIB}}:fit
      A-sz-fit -> gboolean: A-sz-fit
        A^{\mathrm{tSZ}}:fit
      ps-A-100-100-fit -> gboolean: ps-A-100-100-fit
        A^{\mathrm{PS}}_{100}:fit
      ps-A-143-143-fit -> gboolean: ps-A-143-143-fit
        A^{\mathrm{PS}}_{143}:fit
      ps-A-143-217-fit -> gboolean: ps-A-143-217-fit
        A^{\mathrm{PS}}_{143\times 217}:fit
      ps-A-217-217-fit -> gboolean: ps-A-217-217-fit
        A^{\mathrm{PS}}_{217}:fit
      ksz-norm-fit -> gboolean: ksz-norm-fit
        A^{\mathrm{kSZ}}:fit
      gal545-A-100-fit -> gboolean: gal545-A-100-fit
        A^{\mathrm{dust}TT}_{100}:fit
      gal545-A-143-fit -> gboolean: gal545-A-143-fit
        A^{\mathrm{dust}TT}_{143}:fit
      gal545-A-143-217-fit -> gboolean: gal545-A-143-217-fit
        A^{\mathrm{dust}TT}_{143 \times 217}:fit
      gal545-A-217-fit -> gboolean: gal545-A-217-fit
        A^{\mathrm{dust}TT}_{217}:fit
      A-sbpx-100-100-TT-fit -> gboolean: A-sbpx-100-100-TT-fit
        A^{\mathrm{sbpx}TT}_{100 \times 100}:fit
      A-sbpx-143-143-TT-fit -> gboolean: A-sbpx-143-143-TT-fit
        A^{\mathrm{sbpx}TT}_{143 \times 143}:fit
      A-sbpx-143-217-TT-fit -> gboolean: A-sbpx-143-217-TT-fit
        A^{\mathrm{sbpx}TT}_{143 \times 217}:fit
      A-sbpx-217-217-TT-fit -> gboolean: A-sbpx-217-217-TT-fit
        A^{\mathrm{sbpx}TT}_{217 \times 217}:fit
      calib-100T-fit -> gboolean: calib-100T-fit
        c_{100}:fit
      calib-217T-fit -> gboolean: calib-217T-fit
        c_{217}:fit
      A-planck-fit -> gboolean: A-planck-fit
        y_{\mathrm{cal}}:fit

    Properties from NcPlanckFI:
      version -> guint: version
        Planck compatible version

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        A_cib_217: float
        A_cib_217_fit: bool
        A_planck: float
        A_planck_fit: bool
        A_sbpx_100_100_TT: float
        A_sbpx_100_100_TT_fit: bool
        A_sbpx_143_143_TT: float
        A_sbpx_143_143_TT_fit: bool
        A_sbpx_143_217_TT: float
        A_sbpx_143_217_TT_fit: bool
        A_sbpx_217_217_TT: float
        A_sbpx_217_217_TT_fit: bool
        A_sz: float
        A_sz_fit: bool
        calib_100T: float
        calib_100T_fit: bool
        calib_217T: float
        calib_217T_fit: bool
        cib_index: float
        cib_index_fit: bool
        gal545_A_100: float
        gal545_A_100_fit: bool
        gal545_A_143: float
        gal545_A_143_217: float
        gal545_A_143_217_fit: bool
        gal545_A_143_fit: bool
        gal545_A_217: float
        gal545_A_217_fit: bool
        ksz_norm: float
        ksz_norm_fit: bool
        ps_A_100_100: float
        ps_A_100_100_fit: bool
        ps_A_143_143: float
        ps_A_143_143_fit: bool
        ps_A_143_217: float
        ps_A_143_217_fit: bool
        ps_A_217_217: float
        ps_A_217_217_fit: bool
        xi_sz_cib: float
        xi_sz_cib_fit: bool
        version: int
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: PlanckFI = ...
    def __init__(
        self,
        A_cib_217: float = ...,
        A_cib_217_fit: bool = ...,
        A_planck: float = ...,
        A_planck_fit: bool = ...,
        A_sbpx_100_100_TT: float = ...,
        A_sbpx_100_100_TT_fit: bool = ...,
        A_sbpx_143_143_TT: float = ...,
        A_sbpx_143_143_TT_fit: bool = ...,
        A_sbpx_143_217_TT: float = ...,
        A_sbpx_143_217_TT_fit: bool = ...,
        A_sbpx_217_217_TT: float = ...,
        A_sbpx_217_217_TT_fit: bool = ...,
        A_sz: float = ...,
        A_sz_fit: bool = ...,
        calib_100T: float = ...,
        calib_100T_fit: bool = ...,
        calib_217T: float = ...,
        calib_217T_fit: bool = ...,
        cib_index: float = ...,
        cib_index_fit: bool = ...,
        gal545_A_100: float = ...,
        gal545_A_100_fit: bool = ...,
        gal545_A_143: float = ...,
        gal545_A_143_217: float = ...,
        gal545_A_143_217_fit: bool = ...,
        gal545_A_143_fit: bool = ...,
        gal545_A_217: float = ...,
        gal545_A_217_fit: bool = ...,
        ksz_norm: float = ...,
        ksz_norm_fit: bool = ...,
        ps_A_100_100: float = ...,
        ps_A_100_100_fit: bool = ...,
        ps_A_143_143: float = ...,
        ps_A_143_143_fit: bool = ...,
        ps_A_143_217: float = ...,
        ps_A_143_217_fit: bool = ...,
        ps_A_217_217: float = ...,
        ps_A_217_217_fit: bool = ...,
        xi_sz_cib: float = ...,
        xi_sz_cib_fit: bool = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def add_all_default18_priors(lh: NumCosmoMath.Likelihood) -> None: ...
    @staticmethod
    def add_all_default_priors(lh: NumCosmoMath.Likelihood) -> None: ...
    @staticmethod
    def add_calib_priors(
        lh: NumCosmoMath.Likelihood,
        mean: NumCosmoMath.Vector,
        sigma: NumCosmoMath.Vector,
    ) -> None: ...
    @staticmethod
    def add_default18_calib_priors(lh: NumCosmoMath.Likelihood) -> None: ...
    @staticmethod
    def add_default18_gal_priors(lh: NumCosmoMath.Likelihood) -> None: ...
    @staticmethod
    def add_default18_sz_prior(lh: NumCosmoMath.Likelihood) -> None: ...
    @staticmethod
    def add_default_calib_priors(lh: NumCosmoMath.Likelihood) -> None: ...
    @staticmethod
    def add_default_gal_priors(lh: NumCosmoMath.Likelihood) -> None: ...
    @staticmethod
    def add_default_sz_prior(lh: NumCosmoMath.Likelihood) -> None: ...
    @staticmethod
    def add_gal_priors(
        lh: NumCosmoMath.Likelihood,
        mean: NumCosmoMath.Vector,
        sigma: NumCosmoMath.Vector,
    ) -> None: ...
    @staticmethod
    def add_sz_prior(
        lh: NumCosmoMath.Likelihood, f_tSZ: float, mean: float, sigma: float
    ) -> None: ...

class PlanckFICorTTClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        PlanckFICorTTClass()
    """

    parent_class: PlanckFIClass = ...

class PlanckFICorTTTEEE(PlanckFICorTT):
    r"""
    :Constructors:

    ::

        PlanckFICorTTTEEE(**properties)

    Object NcPlanckFICorTTTEEE

    Properties from NcPlanckFICorTTTEEE:
      galf-EE-A-100 -> gdouble: galf-EE-A-100
        A^{\mathrm{dust}EE}_{100}
      galf-EE-A-100-143 -> gdouble: galf-EE-A-100-143
        A^{\mathrm{dust}EE}_{100 \times 143}
      galf-EE-A-100-217 -> gdouble: galf-EE-A-100-217
        A^{\mathrm{dust}EE}_{100 \times 217}
      galf-EE-A-143 -> gdouble: galf-EE-A-143
        A^{\mathrm{dust}EE}_{143}
      galf-EE-A-143-217 -> gdouble: galf-EE-A-143-217
        A^{\mathrm{dust}EE}_{143 \times 217}
      galf-EE-A-217 -> gdouble: galf-EE-A-217
        A^{\mathrm{dust}EE}_{217}
      galf-EE-index -> gdouble: galf-EE-index
        n^{\mathrm{dust}EE}
      galf-TE-A-100 -> gdouble: galf-TE-A-100
        A^{\mathrm{dust}TE}_{100}
      galf-TE-A-100-143 -> gdouble: galf-TE-A-100-143
        A^{\mathrm{dust}TE}_{100 \times 143}
      galf-TE-A-100-217 -> gdouble: galf-TE-A-100-217
        A^{\mathrm{dust}TE}_{100 \times 217}
      galf-TE-A-143 -> gdouble: galf-TE-A-143
        A^{\mathrm{dust}TE}_{143}
      galf-TE-A-143-217 -> gdouble: galf-TE-A-143-217
        A^{\mathrm{dust}TE}_{143 \times 217}
      galf-TE-A-217 -> gdouble: galf-TE-A-217
        A^{\mathrm{dust}TE}_{217}
      galf-TE-index -> gdouble: galf-TE-index
        n^{\mathrm{dust}TE}
      A-cnoise-e2e-100-100-EE -> gdouble: A-cnoise-e2e-100-100-EE
        A^{\mathrm{cnoise}EE}_{100 \times 100}
      A-cnoise-e2e-143-143-EE -> gdouble: A-cnoise-e2e-143-143-EE
        A^{\mathrm{cnoise}EE}_{143 \times 143}
      A-cnoise-e2e-217-217-EE -> gdouble: A-cnoise-e2e-217-217-EE
        A^{\mathrm{cnoise}EE}_{217 \times 217}
      bleak-epsilon-0-0T-0E -> gdouble: bleak-epsilon-0-0T-0E
        \epsilon^{\mathrm{bleak}TE}_{0, 100}
      bleak-epsilon-1-0T-0E -> gdouble: bleak-epsilon-1-0T-0E
        \epsilon^{\mathrm{bleak}TE}_{1, 100}
      bleak-epsilon-2-0T-0E -> gdouble: bleak-epsilon-2-0T-0E
        \epsilon^{\mathrm{bleak}TE}_{2, 100}
      bleak-epsilon-3-0T-0E -> gdouble: bleak-epsilon-3-0T-0E
        \epsilon^{\mathrm{bleak}TE}_{3, 100}
      bleak-epsilon-4-0T-0E -> gdouble: bleak-epsilon-4-0T-0E
        \epsilon^{\mathrm{bleak}TE}_{4, 100}
      bleak-epsilon-0-0T-1E -> gdouble: bleak-epsilon-0-0T-1E
        \epsilon^{\mathrm{bleak}TE}_{0, 100 \times 143}
      bleak-epsilon-1-0T-1E -> gdouble: bleak-epsilon-1-0T-1E
        \epsilon^{\mathrm{bleak}TE}_{1, 100 \times 143}
      bleak-epsilon-2-0T-1E -> gdouble: bleak-epsilon-2-0T-1E
        \epsilon^{\mathrm{bleak}TE}_{2, 100 \times 143}
      bleak-epsilon-3-0T-1E -> gdouble: bleak-epsilon-3-0T-1E
        \epsilon^{\mathrm{bleak}TE}_{3, 100 \times 143}
      bleak-epsilon-4-0T-1E -> gdouble: bleak-epsilon-4-0T-1E
        \epsilon^{\mathrm{bleak}TE}_{4, 100 \times 143}
      bleak-epsilon-0-0T-2E -> gdouble: bleak-epsilon-0-0T-2E
        \epsilon^{\mathrm{bleak}TE}_{0, 100 \times 217}
      bleak-epsilon-1-0T-2E -> gdouble: bleak-epsilon-1-0T-2E
        \epsilon^{\mathrm{bleak}TE}_{1, 100 \times 217}
      bleak-epsilon-2-0T-2E -> gdouble: bleak-epsilon-2-0T-2E
        \epsilon^{\mathrm{bleak}TE}_{2, 100 \times 217}
      bleak-epsilon-3-0T-2E -> gdouble: bleak-epsilon-3-0T-2E
        \epsilon^{\mathrm{bleak}TE}_{3, 100 \times 217}
      bleak-epsilon-4-0T-2E -> gdouble: bleak-epsilon-4-0T-2E
        \epsilon^{\mathrm{bleak}TE}_{4, 100 \times 217}
      bleak-epsilon-0-1T-1E -> gdouble: bleak-epsilon-0-1T-1E
        \epsilon^{\mathrm{bleak}TE}_{0, 143}
      bleak-epsilon-1-1T-1E -> gdouble: bleak-epsilon-1-1T-1E
        \epsilon^{\mathrm{bleak}TE}_{1, 143}
      bleak-epsilon-2-1T-1E -> gdouble: bleak-epsilon-2-1T-1E
        \epsilon^{\mathrm{bleak}TE}_{2, 143}
      bleak-epsilon-3-1T-1E -> gdouble: bleak-epsilon-3-1T-1E
        \epsilon^{\mathrm{bleak}TE}_{3, 143}
      bleak-epsilon-4-1T-1E -> gdouble: bleak-epsilon-4-1T-1E
        \epsilon^{\mathrm{bleak}TE}_{4, 143}
      bleak-epsilon-0-1T-2E -> gdouble: bleak-epsilon-0-1T-2E
        \epsilon^{\mathrm{bleak}TE}_{0, 143 \times 217}
      bleak-epsilon-1-1T-2E -> gdouble: bleak-epsilon-1-1T-2E
        \epsilon^{\mathrm{bleak}TE}_{1, 143 \times 217}
      bleak-epsilon-2-1T-2E -> gdouble: bleak-epsilon-2-1T-2E
        \epsilon^{\mathrm{bleak}TE}_{2, 143 \times 217}
      bleak-epsilon-3-1T-2E -> gdouble: bleak-epsilon-3-1T-2E
        \epsilon^{\mathrm{bleak}TE}_{3, 143 \times 217}
      bleak-epsilon-4-1T-2E -> gdouble: bleak-epsilon-4-1T-2E
        \epsilon^{\mathrm{bleak}TE}_{4, 143 \times 217}
      bleak-epsilon-0-2T-2E -> gdouble: bleak-epsilon-0-2T-2E
        \epsilon^{\mathrm{bleak}TE}_{0, 217}
      bleak-epsilon-1-2T-2E -> gdouble: bleak-epsilon-1-2T-2E
        \epsilon^{\mathrm{bleak}TE}_{1, 217}
      bleak-epsilon-2-2T-2E -> gdouble: bleak-epsilon-2-2T-2E
        \epsilon^{\mathrm{bleak}TE}_{2, 217}
      bleak-epsilon-3-2T-2E -> gdouble: bleak-epsilon-3-2T-2E
        \epsilon^{\mathrm{bleak}TE}_{3, 217}
      bleak-epsilon-4-2T-2E -> gdouble: bleak-epsilon-4-2T-2E
        \epsilon^{\mathrm{bleak}TE}_{4, 217}
      bleak-epsilon-0-0E-0E -> gdouble: bleak-epsilon-0-0E-0E
        \epsilon^{\mathrm{bleak}EE}_{0, 100}
      bleak-epsilon-1-0E-0E -> gdouble: bleak-epsilon-1-0E-0E
        \epsilon^{\mathrm{bleak}EE}_{1, 100}
      bleak-epsilon-2-0E-0E -> gdouble: bleak-epsilon-2-0E-0E
        \epsilon^{\mathrm{bleak}EE}_{2, 100}
      bleak-epsilon-3-0E-0E -> gdouble: bleak-epsilon-3-0E-0E
        \epsilon^{\mathrm{bleak}EE}_{3, 100}
      bleak-epsilon-4-0E-0E -> gdouble: bleak-epsilon-4-0E-0E
        \epsilon^{\mathrm{bleak}EE}_{4, 100}
      bleak-epsilon-0-0E-1E -> gdouble: bleak-epsilon-0-0E-1E
        \epsilon^{\mathrm{bleak}EE}_{0, 100 \times 143}
      bleak-epsilon-1-0E-1E -> gdouble: bleak-epsilon-1-0E-1E
        \epsilon^{\mathrm{bleak}EE}_{1, 100 \times 143}
      bleak-epsilon-2-0E-1E -> gdouble: bleak-epsilon-2-0E-1E
        \epsilon^{\mathrm{bleak}EE}_{2, 100 \times 143}
      bleak-epsilon-3-0E-1E -> gdouble: bleak-epsilon-3-0E-1E
        \epsilon^{\mathrm{bleak}EE}_{3, 100 \times 143}
      bleak-epsilon-4-0E-1E -> gdouble: bleak-epsilon-4-0E-1E
        \epsilon^{\mathrm{bleak}EE}_{4, 100 \times 143}
      bleak-epsilon-0-0E-2E -> gdouble: bleak-epsilon-0-0E-2E
        \epsilon^{\mathrm{bleak}EE}_{0, 100 \times 217}
      bleak-epsilon-1-0E-2E -> gdouble: bleak-epsilon-1-0E-2E
        \epsilon^{\mathrm{bleak}EE}_{1, 100 \times 217}
      bleak-epsilon-2-0E-2E -> gdouble: bleak-epsilon-2-0E-2E
        \epsilon^{\mathrm{bleak}EE}_{2, 100 \times 217}
      bleak-epsilon-3-0E-2E -> gdouble: bleak-epsilon-3-0E-2E
        \epsilon^{\mathrm{bleak}EE}_{3, 100 \times 217}
      bleak-epsilon-4-0E-2E -> gdouble: bleak-epsilon-4-0E-2E
        \epsilon^{\mathrm{bleak}EE}_{4, 100 \times 217}
      bleak-epsilon-0-1E-1E -> gdouble: bleak-epsilon-0-1E-1E
        \epsilon^{\mathrm{bleak}EE}_{0, 143}
      bleak-epsilon-1-1E-1E -> gdouble: bleak-epsilon-1-1E-1E
        \epsilon^{\mathrm{bleak}EE}_{1, 143}
      bleak-epsilon-2-1E-1E -> gdouble: bleak-epsilon-2-1E-1E
        \epsilon^{\mathrm{bleak}EE}_{2, 143}
      bleak-epsilon-3-1E-1E -> gdouble: bleak-epsilon-3-1E-1E
        \epsilon^{\mathrm{bleak}EE}_{3, 143}
      bleak-epsilon-4-1E-1E -> gdouble: bleak-epsilon-4-1E-1E
        \epsilon^{\mathrm{bleak}EE}_{4, 143}
      bleak-epsilon-0-1E-2E -> gdouble: bleak-epsilon-0-1E-2E
        \epsilon^{\mathrm{bleak}EE}_{0, 143 \times 217}
      bleak-epsilon-1-1E-2E -> gdouble: bleak-epsilon-1-1E-2E
        \epsilon^{\mathrm{bleak}EE}_{1, 143 \times 217}
      bleak-epsilon-2-1E-2E -> gdouble: bleak-epsilon-2-1E-2E
        \epsilon^{\mathrm{bleak}EE}_{2, 143 \times 217}
      bleak-epsilon-3-1E-2E -> gdouble: bleak-epsilon-3-1E-2E
        \epsilon^{\mathrm{bleak}EE}_{3, 143 \times 217}
      bleak-epsilon-4-1E-2E -> gdouble: bleak-epsilon-4-1E-2E
        \epsilon^{\mathrm{bleak}EE}_{4, 143 \times 217}
      bleak-epsilon-0-2E-2E -> gdouble: bleak-epsilon-0-2E-2E
        \epsilon^{\mathrm{bleak}EE}_{0, 217}
      bleak-epsilon-1-2E-2E -> gdouble: bleak-epsilon-1-2E-2E
        \epsilon^{\mathrm{bleak}EE}_{1, 217}
      bleak-epsilon-2-2E-2E -> gdouble: bleak-epsilon-2-2E-2E
        \epsilon^{\mathrm{bleak}EE}_{2, 217}
      bleak-epsilon-3-2E-2E -> gdouble: bleak-epsilon-3-2E-2E
        \epsilon^{\mathrm{bleak}EE}_{3, 217}
      bleak-epsilon-4-2E-2E -> gdouble: bleak-epsilon-4-2E-2E
        \epsilon^{\mathrm{bleak}EE}_{4, 217}
      A-sbpx-100-100-EE -> gdouble: A-sbpx-100-100-EE
        A^{\mathrm{sbpx}EE}_{100 \times 100}
      A-sbpx-100-143-EE -> gdouble: A-sbpx-100-143-EE
        A^{\mathrm{sbpx}EE}_{100 \times 143}
      A-sbpx-100-217-EE -> gdouble: A-sbpx-100-217-EE
        A^{\mathrm{sbpx}EE}_{100 \times 217}
      A-sbpx-143-143-EE -> gdouble: A-sbpx-143-143-EE
        A^{\mathrm{sbpx}EE}_{143 \times 143}
      A-sbpx-143-217-EE -> gdouble: A-sbpx-143-217-EE
        A^{\mathrm{sbpx}EE}_{143 \times 217}
      A-sbpx-217-217-EE -> gdouble: A-sbpx-217-217-EE
        A^{\mathrm{sbpx}EE}_{217 \times 217}
      calib-100P -> gdouble: calib-100P
        c_{100P}
      calib-143P -> gdouble: calib-143P
        c_{143P}
      calib-217P -> gdouble: calib-217P
        c_{217P}
      A-pol -> gdouble: A-pol
        A_{\mathrm{pol}}
      galf-EE-A-100-fit -> gboolean: galf-EE-A-100-fit
        A^{\mathrm{dust}EE}_{100}:fit
      galf-EE-A-100-143-fit -> gboolean: galf-EE-A-100-143-fit
        A^{\mathrm{dust}EE}_{100 \times 143}:fit
      galf-EE-A-100-217-fit -> gboolean: galf-EE-A-100-217-fit
        A^{\mathrm{dust}EE}_{100 \times 217}:fit
      galf-EE-A-143-fit -> gboolean: galf-EE-A-143-fit
        A^{\mathrm{dust}EE}_{143}:fit
      galf-EE-A-143-217-fit -> gboolean: galf-EE-A-143-217-fit
        A^{\mathrm{dust}EE}_{143 \times 217}:fit
      galf-EE-A-217-fit -> gboolean: galf-EE-A-217-fit
        A^{\mathrm{dust}EE}_{217}:fit
      galf-EE-index-fit -> gboolean: galf-EE-index-fit
        n^{\mathrm{dust}EE}:fit
      galf-TE-A-100-fit -> gboolean: galf-TE-A-100-fit
        A^{\mathrm{dust}TE}_{100}:fit
      galf-TE-A-100-143-fit -> gboolean: galf-TE-A-100-143-fit
        A^{\mathrm{dust}TE}_{100 \times 143}:fit
      galf-TE-A-100-217-fit -> gboolean: galf-TE-A-100-217-fit
        A^{\mathrm{dust}TE}_{100 \times 217}:fit
      galf-TE-A-143-fit -> gboolean: galf-TE-A-143-fit
        A^{\mathrm{dust}TE}_{143}:fit
      galf-TE-A-143-217-fit -> gboolean: galf-TE-A-143-217-fit
        A^{\mathrm{dust}TE}_{143 \times 217}:fit
      galf-TE-A-217-fit -> gboolean: galf-TE-A-217-fit
        A^{\mathrm{dust}TE}_{217}:fit
      galf-TE-index-fit -> gboolean: galf-TE-index-fit
        n^{\mathrm{dust}TE}:fit
      A-cnoise-e2e-100-100-EE-fit -> gboolean: A-cnoise-e2e-100-100-EE-fit
        A^{\mathrm{cnoise}EE}_{100 \times 100}:fit
      A-cnoise-e2e-143-143-EE-fit -> gboolean: A-cnoise-e2e-143-143-EE-fit
        A^{\mathrm{cnoise}EE}_{143 \times 143}:fit
      A-cnoise-e2e-217-217-EE-fit -> gboolean: A-cnoise-e2e-217-217-EE-fit
        A^{\mathrm{cnoise}EE}_{217 \times 217}:fit
      bleak-epsilon-0-0T-0E-fit -> gboolean: bleak-epsilon-0-0T-0E-fit
        \epsilon^{\mathrm{bleak}TE}_{0, 100}:fit
      bleak-epsilon-1-0T-0E-fit -> gboolean: bleak-epsilon-1-0T-0E-fit
        \epsilon^{\mathrm{bleak}TE}_{1, 100}:fit
      bleak-epsilon-2-0T-0E-fit -> gboolean: bleak-epsilon-2-0T-0E-fit
        \epsilon^{\mathrm{bleak}TE}_{2, 100}:fit
      bleak-epsilon-3-0T-0E-fit -> gboolean: bleak-epsilon-3-0T-0E-fit
        \epsilon^{\mathrm{bleak}TE}_{3, 100}:fit
      bleak-epsilon-4-0T-0E-fit -> gboolean: bleak-epsilon-4-0T-0E-fit
        \epsilon^{\mathrm{bleak}TE}_{4, 100}:fit
      bleak-epsilon-0-0T-1E-fit -> gboolean: bleak-epsilon-0-0T-1E-fit
        \epsilon^{\mathrm{bleak}TE}_{0, 100 \times 143}:fit
      bleak-epsilon-1-0T-1E-fit -> gboolean: bleak-epsilon-1-0T-1E-fit
        \epsilon^{\mathrm{bleak}TE}_{1, 100 \times 143}:fit
      bleak-epsilon-2-0T-1E-fit -> gboolean: bleak-epsilon-2-0T-1E-fit
        \epsilon^{\mathrm{bleak}TE}_{2, 100 \times 143}:fit
      bleak-epsilon-3-0T-1E-fit -> gboolean: bleak-epsilon-3-0T-1E-fit
        \epsilon^{\mathrm{bleak}TE}_{3, 100 \times 143}:fit
      bleak-epsilon-4-0T-1E-fit -> gboolean: bleak-epsilon-4-0T-1E-fit
        \epsilon^{\mathrm{bleak}TE}_{4, 100 \times 143}:fit
      bleak-epsilon-0-0T-2E-fit -> gboolean: bleak-epsilon-0-0T-2E-fit
        \epsilon^{\mathrm{bleak}TE}_{0, 100 \times 217}:fit
      bleak-epsilon-1-0T-2E-fit -> gboolean: bleak-epsilon-1-0T-2E-fit
        \epsilon^{\mathrm{bleak}TE}_{1, 100 \times 217}:fit
      bleak-epsilon-2-0T-2E-fit -> gboolean: bleak-epsilon-2-0T-2E-fit
        \epsilon^{\mathrm{bleak}TE}_{2, 100 \times 217}:fit
      bleak-epsilon-3-0T-2E-fit -> gboolean: bleak-epsilon-3-0T-2E-fit
        \epsilon^{\mathrm{bleak}TE}_{3, 100 \times 217}:fit
      bleak-epsilon-4-0T-2E-fit -> gboolean: bleak-epsilon-4-0T-2E-fit
        \epsilon^{\mathrm{bleak}TE}_{4, 100 \times 217}:fit
      bleak-epsilon-0-1T-1E-fit -> gboolean: bleak-epsilon-0-1T-1E-fit
        \epsilon^{\mathrm{bleak}TE}_{0, 143}:fit
      bleak-epsilon-1-1T-1E-fit -> gboolean: bleak-epsilon-1-1T-1E-fit
        \epsilon^{\mathrm{bleak}TE}_{1, 143}:fit
      bleak-epsilon-2-1T-1E-fit -> gboolean: bleak-epsilon-2-1T-1E-fit
        \epsilon^{\mathrm{bleak}TE}_{2, 143}:fit
      bleak-epsilon-3-1T-1E-fit -> gboolean: bleak-epsilon-3-1T-1E-fit
        \epsilon^{\mathrm{bleak}TE}_{3, 143}:fit
      bleak-epsilon-4-1T-1E-fit -> gboolean: bleak-epsilon-4-1T-1E-fit
        \epsilon^{\mathrm{bleak}TE}_{4, 143}:fit
      bleak-epsilon-0-1T-2E-fit -> gboolean: bleak-epsilon-0-1T-2E-fit
        \epsilon^{\mathrm{bleak}TE}_{0, 143 \times 217}:fit
      bleak-epsilon-1-1T-2E-fit -> gboolean: bleak-epsilon-1-1T-2E-fit
        \epsilon^{\mathrm{bleak}TE}_{1, 143 \times 217}:fit
      bleak-epsilon-2-1T-2E-fit -> gboolean: bleak-epsilon-2-1T-2E-fit
        \epsilon^{\mathrm{bleak}TE}_{2, 143 \times 217}:fit
      bleak-epsilon-3-1T-2E-fit -> gboolean: bleak-epsilon-3-1T-2E-fit
        \epsilon^{\mathrm{bleak}TE}_{3, 143 \times 217}:fit
      bleak-epsilon-4-1T-2E-fit -> gboolean: bleak-epsilon-4-1T-2E-fit
        \epsilon^{\mathrm{bleak}TE}_{4, 143 \times 217}:fit
      bleak-epsilon-0-2T-2E-fit -> gboolean: bleak-epsilon-0-2T-2E-fit
        \epsilon^{\mathrm{bleak}TE}_{0, 217}:fit
      bleak-epsilon-1-2T-2E-fit -> gboolean: bleak-epsilon-1-2T-2E-fit
        \epsilon^{\mathrm{bleak}TE}_{1, 217}:fit
      bleak-epsilon-2-2T-2E-fit -> gboolean: bleak-epsilon-2-2T-2E-fit
        \epsilon^{\mathrm{bleak}TE}_{2, 217}:fit
      bleak-epsilon-3-2T-2E-fit -> gboolean: bleak-epsilon-3-2T-2E-fit
        \epsilon^{\mathrm{bleak}TE}_{3, 217}:fit
      bleak-epsilon-4-2T-2E-fit -> gboolean: bleak-epsilon-4-2T-2E-fit
        \epsilon^{\mathrm{bleak}TE}_{4, 217}:fit
      bleak-epsilon-0-0E-0E-fit -> gboolean: bleak-epsilon-0-0E-0E-fit
        \epsilon^{\mathrm{bleak}EE}_{0, 100}:fit
      bleak-epsilon-1-0E-0E-fit -> gboolean: bleak-epsilon-1-0E-0E-fit
        \epsilon^{\mathrm{bleak}EE}_{1, 100}:fit
      bleak-epsilon-2-0E-0E-fit -> gboolean: bleak-epsilon-2-0E-0E-fit
        \epsilon^{\mathrm{bleak}EE}_{2, 100}:fit
      bleak-epsilon-3-0E-0E-fit -> gboolean: bleak-epsilon-3-0E-0E-fit
        \epsilon^{\mathrm{bleak}EE}_{3, 100}:fit
      bleak-epsilon-4-0E-0E-fit -> gboolean: bleak-epsilon-4-0E-0E-fit
        \epsilon^{\mathrm{bleak}EE}_{4, 100}:fit
      bleak-epsilon-0-0E-1E-fit -> gboolean: bleak-epsilon-0-0E-1E-fit
        \epsilon^{\mathrm{bleak}EE}_{0, 100 \times 143}:fit
      bleak-epsilon-1-0E-1E-fit -> gboolean: bleak-epsilon-1-0E-1E-fit
        \epsilon^{\mathrm{bleak}EE}_{1, 100 \times 143}:fit
      bleak-epsilon-2-0E-1E-fit -> gboolean: bleak-epsilon-2-0E-1E-fit
        \epsilon^{\mathrm{bleak}EE}_{2, 100 \times 143}:fit
      bleak-epsilon-3-0E-1E-fit -> gboolean: bleak-epsilon-3-0E-1E-fit
        \epsilon^{\mathrm{bleak}EE}_{3, 100 \times 143}:fit
      bleak-epsilon-4-0E-1E-fit -> gboolean: bleak-epsilon-4-0E-1E-fit
        \epsilon^{\mathrm{bleak}EE}_{4, 100 \times 143}:fit
      bleak-epsilon-0-0E-2E-fit -> gboolean: bleak-epsilon-0-0E-2E-fit
        \epsilon^{\mathrm{bleak}EE}_{0, 100 \times 217}:fit
      bleak-epsilon-1-0E-2E-fit -> gboolean: bleak-epsilon-1-0E-2E-fit
        \epsilon^{\mathrm{bleak}EE}_{1, 100 \times 217}:fit
      bleak-epsilon-2-0E-2E-fit -> gboolean: bleak-epsilon-2-0E-2E-fit
        \epsilon^{\mathrm{bleak}EE}_{2, 100 \times 217}:fit
      bleak-epsilon-3-0E-2E-fit -> gboolean: bleak-epsilon-3-0E-2E-fit
        \epsilon^{\mathrm{bleak}EE}_{3, 100 \times 217}:fit
      bleak-epsilon-4-0E-2E-fit -> gboolean: bleak-epsilon-4-0E-2E-fit
        \epsilon^{\mathrm{bleak}EE}_{4, 100 \times 217}:fit
      bleak-epsilon-0-1E-1E-fit -> gboolean: bleak-epsilon-0-1E-1E-fit
        \epsilon^{\mathrm{bleak}EE}_{0, 143}:fit
      bleak-epsilon-1-1E-1E-fit -> gboolean: bleak-epsilon-1-1E-1E-fit
        \epsilon^{\mathrm{bleak}EE}_{1, 143}:fit
      bleak-epsilon-2-1E-1E-fit -> gboolean: bleak-epsilon-2-1E-1E-fit
        \epsilon^{\mathrm{bleak}EE}_{2, 143}:fit
      bleak-epsilon-3-1E-1E-fit -> gboolean: bleak-epsilon-3-1E-1E-fit
        \epsilon^{\mathrm{bleak}EE}_{3, 143}:fit
      bleak-epsilon-4-1E-1E-fit -> gboolean: bleak-epsilon-4-1E-1E-fit
        \epsilon^{\mathrm{bleak}EE}_{4, 143}:fit
      bleak-epsilon-0-1E-2E-fit -> gboolean: bleak-epsilon-0-1E-2E-fit
        \epsilon^{\mathrm{bleak}EE}_{0, 143 \times 217}:fit
      bleak-epsilon-1-1E-2E-fit -> gboolean: bleak-epsilon-1-1E-2E-fit
        \epsilon^{\mathrm{bleak}EE}_{1, 143 \times 217}:fit
      bleak-epsilon-2-1E-2E-fit -> gboolean: bleak-epsilon-2-1E-2E-fit
        \epsilon^{\mathrm{bleak}EE}_{2, 143 \times 217}:fit
      bleak-epsilon-3-1E-2E-fit -> gboolean: bleak-epsilon-3-1E-2E-fit
        \epsilon^{\mathrm{bleak}EE}_{3, 143 \times 217}:fit
      bleak-epsilon-4-1E-2E-fit -> gboolean: bleak-epsilon-4-1E-2E-fit
        \epsilon^{\mathrm{bleak}EE}_{4, 143 \times 217}:fit
      bleak-epsilon-0-2E-2E-fit -> gboolean: bleak-epsilon-0-2E-2E-fit
        \epsilon^{\mathrm{bleak}EE}_{0, 217}:fit
      bleak-epsilon-1-2E-2E-fit -> gboolean: bleak-epsilon-1-2E-2E-fit
        \epsilon^{\mathrm{bleak}EE}_{1, 217}:fit
      bleak-epsilon-2-2E-2E-fit -> gboolean: bleak-epsilon-2-2E-2E-fit
        \epsilon^{\mathrm{bleak}EE}_{2, 217}:fit
      bleak-epsilon-3-2E-2E-fit -> gboolean: bleak-epsilon-3-2E-2E-fit
        \epsilon^{\mathrm{bleak}EE}_{3, 217}:fit
      bleak-epsilon-4-2E-2E-fit -> gboolean: bleak-epsilon-4-2E-2E-fit
        \epsilon^{\mathrm{bleak}EE}_{4, 217}:fit
      A-sbpx-100-100-EE-fit -> gboolean: A-sbpx-100-100-EE-fit
        A^{\mathrm{sbpx}EE}_{100 \times 100}:fit
      A-sbpx-100-143-EE-fit -> gboolean: A-sbpx-100-143-EE-fit
        A^{\mathrm{sbpx}EE}_{100 \times 143}:fit
      A-sbpx-100-217-EE-fit -> gboolean: A-sbpx-100-217-EE-fit
        A^{\mathrm{sbpx}EE}_{100 \times 217}:fit
      A-sbpx-143-143-EE-fit -> gboolean: A-sbpx-143-143-EE-fit
        A^{\mathrm{sbpx}EE}_{143 \times 143}:fit
      A-sbpx-143-217-EE-fit -> gboolean: A-sbpx-143-217-EE-fit
        A^{\mathrm{sbpx}EE}_{143 \times 217}:fit
      A-sbpx-217-217-EE-fit -> gboolean: A-sbpx-217-217-EE-fit
        A^{\mathrm{sbpx}EE}_{217 \times 217}:fit
      calib-100P-fit -> gboolean: calib-100P-fit
        c_{100P}:fit
      calib-143P-fit -> gboolean: calib-143P-fit
        c_{143P}:fit
      calib-217P-fit -> gboolean: calib-217P-fit
        c_{217P}:fit
      A-pol-fit -> gboolean: A-pol-fit
        A_{\mathrm{pol}}:fit

    Properties from NcPlanckFICorTT:
      A-cib-217 -> gdouble: A-cib-217
        A^{\mathrm{CIB}}_{217}
      cib-index -> gdouble: cib-index
        n^{\mathrm{CIB}}
      xi-sz-cib -> gdouble: xi-sz-cib
        \xi^{\mathrm{tSZ}\times \mathrm{CIB}}
      A-sz -> gdouble: A-sz
        A^{\mathrm{tSZ}}
      ps-A-100-100 -> gdouble: ps-A-100-100
        A^{\mathrm{PS}}_{100}
      ps-A-143-143 -> gdouble: ps-A-143-143
        A^{\mathrm{PS}}_{143}
      ps-A-143-217 -> gdouble: ps-A-143-217
        A^{\mathrm{PS}}_{143\times 217}
      ps-A-217-217 -> gdouble: ps-A-217-217
        A^{\mathrm{PS}}_{217}
      ksz-norm -> gdouble: ksz-norm
        A^{\mathrm{kSZ}}
      gal545-A-100 -> gdouble: gal545-A-100
        A^{\mathrm{dust}TT}_{100}
      gal545-A-143 -> gdouble: gal545-A-143
        A^{\mathrm{dust}TT}_{143}
      gal545-A-143-217 -> gdouble: gal545-A-143-217
        A^{\mathrm{dust}TT}_{143 \times 217}
      gal545-A-217 -> gdouble: gal545-A-217
        A^{\mathrm{dust}TT}_{217}
      A-sbpx-100-100-TT -> gdouble: A-sbpx-100-100-TT
        A^{\mathrm{sbpx}TT}_{100 \times 100}
      A-sbpx-143-143-TT -> gdouble: A-sbpx-143-143-TT
        A^{\mathrm{sbpx}TT}_{143 \times 143}
      A-sbpx-143-217-TT -> gdouble: A-sbpx-143-217-TT
        A^{\mathrm{sbpx}TT}_{143 \times 217}
      A-sbpx-217-217-TT -> gdouble: A-sbpx-217-217-TT
        A^{\mathrm{sbpx}TT}_{217 \times 217}
      calib-100T -> gdouble: calib-100T
        c_{100}
      calib-217T -> gdouble: calib-217T
        c_{217}
      A-planck -> gdouble: A-planck
        y_{\mathrm{cal}}
      A-cib-217-fit -> gboolean: A-cib-217-fit
        A^{\mathrm{CIB}}_{217}:fit
      cib-index-fit -> gboolean: cib-index-fit
        n^{\mathrm{CIB}}:fit
      xi-sz-cib-fit -> gboolean: xi-sz-cib-fit
        \xi^{\mathrm{tSZ}\times \mathrm{CIB}}:fit
      A-sz-fit -> gboolean: A-sz-fit
        A^{\mathrm{tSZ}}:fit
      ps-A-100-100-fit -> gboolean: ps-A-100-100-fit
        A^{\mathrm{PS}}_{100}:fit
      ps-A-143-143-fit -> gboolean: ps-A-143-143-fit
        A^{\mathrm{PS}}_{143}:fit
      ps-A-143-217-fit -> gboolean: ps-A-143-217-fit
        A^{\mathrm{PS}}_{143\times 217}:fit
      ps-A-217-217-fit -> gboolean: ps-A-217-217-fit
        A^{\mathrm{PS}}_{217}:fit
      ksz-norm-fit -> gboolean: ksz-norm-fit
        A^{\mathrm{kSZ}}:fit
      gal545-A-100-fit -> gboolean: gal545-A-100-fit
        A^{\mathrm{dust}TT}_{100}:fit
      gal545-A-143-fit -> gboolean: gal545-A-143-fit
        A^{\mathrm{dust}TT}_{143}:fit
      gal545-A-143-217-fit -> gboolean: gal545-A-143-217-fit
        A^{\mathrm{dust}TT}_{143 \times 217}:fit
      gal545-A-217-fit -> gboolean: gal545-A-217-fit
        A^{\mathrm{dust}TT}_{217}:fit
      A-sbpx-100-100-TT-fit -> gboolean: A-sbpx-100-100-TT-fit
        A^{\mathrm{sbpx}TT}_{100 \times 100}:fit
      A-sbpx-143-143-TT-fit -> gboolean: A-sbpx-143-143-TT-fit
        A^{\mathrm{sbpx}TT}_{143 \times 143}:fit
      A-sbpx-143-217-TT-fit -> gboolean: A-sbpx-143-217-TT-fit
        A^{\mathrm{sbpx}TT}_{143 \times 217}:fit
      A-sbpx-217-217-TT-fit -> gboolean: A-sbpx-217-217-TT-fit
        A^{\mathrm{sbpx}TT}_{217 \times 217}:fit
      calib-100T-fit -> gboolean: calib-100T-fit
        c_{100}:fit
      calib-217T-fit -> gboolean: calib-217T-fit
        c_{217}:fit
      A-planck-fit -> gboolean: A-planck-fit
        y_{\mathrm{cal}}:fit

    Properties from NcPlanckFI:
      version -> guint: version
        Planck compatible version

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        A_cnoise_e2e_100_100_EE: float
        A_cnoise_e2e_100_100_EE_fit: bool
        A_cnoise_e2e_143_143_EE: float
        A_cnoise_e2e_143_143_EE_fit: bool
        A_cnoise_e2e_217_217_EE: float
        A_cnoise_e2e_217_217_EE_fit: bool
        A_pol: float
        A_pol_fit: bool
        A_sbpx_100_100_EE: float
        A_sbpx_100_100_EE_fit: bool
        A_sbpx_100_143_EE: float
        A_sbpx_100_143_EE_fit: bool
        A_sbpx_100_217_EE: float
        A_sbpx_100_217_EE_fit: bool
        A_sbpx_143_143_EE: float
        A_sbpx_143_143_EE_fit: bool
        A_sbpx_143_217_EE: float
        A_sbpx_143_217_EE_fit: bool
        A_sbpx_217_217_EE: float
        A_sbpx_217_217_EE_fit: bool
        bleak_epsilon_0_0E_0E: float
        bleak_epsilon_0_0E_0E_fit: bool
        bleak_epsilon_0_0E_1E: float
        bleak_epsilon_0_0E_1E_fit: bool
        bleak_epsilon_0_0E_2E: float
        bleak_epsilon_0_0E_2E_fit: bool
        bleak_epsilon_0_0T_0E: float
        bleak_epsilon_0_0T_0E_fit: bool
        bleak_epsilon_0_0T_1E: float
        bleak_epsilon_0_0T_1E_fit: bool
        bleak_epsilon_0_0T_2E: float
        bleak_epsilon_0_0T_2E_fit: bool
        bleak_epsilon_0_1E_1E: float
        bleak_epsilon_0_1E_1E_fit: bool
        bleak_epsilon_0_1E_2E: float
        bleak_epsilon_0_1E_2E_fit: bool
        bleak_epsilon_0_1T_1E: float
        bleak_epsilon_0_1T_1E_fit: bool
        bleak_epsilon_0_1T_2E: float
        bleak_epsilon_0_1T_2E_fit: bool
        bleak_epsilon_0_2E_2E: float
        bleak_epsilon_0_2E_2E_fit: bool
        bleak_epsilon_0_2T_2E: float
        bleak_epsilon_0_2T_2E_fit: bool
        bleak_epsilon_1_0E_0E: float
        bleak_epsilon_1_0E_0E_fit: bool
        bleak_epsilon_1_0E_1E: float
        bleak_epsilon_1_0E_1E_fit: bool
        bleak_epsilon_1_0E_2E: float
        bleak_epsilon_1_0E_2E_fit: bool
        bleak_epsilon_1_0T_0E: float
        bleak_epsilon_1_0T_0E_fit: bool
        bleak_epsilon_1_0T_1E: float
        bleak_epsilon_1_0T_1E_fit: bool
        bleak_epsilon_1_0T_2E: float
        bleak_epsilon_1_0T_2E_fit: bool
        bleak_epsilon_1_1E_1E: float
        bleak_epsilon_1_1E_1E_fit: bool
        bleak_epsilon_1_1E_2E: float
        bleak_epsilon_1_1E_2E_fit: bool
        bleak_epsilon_1_1T_1E: float
        bleak_epsilon_1_1T_1E_fit: bool
        bleak_epsilon_1_1T_2E: float
        bleak_epsilon_1_1T_2E_fit: bool
        bleak_epsilon_1_2E_2E: float
        bleak_epsilon_1_2E_2E_fit: bool
        bleak_epsilon_1_2T_2E: float
        bleak_epsilon_1_2T_2E_fit: bool
        bleak_epsilon_2_0E_0E: float
        bleak_epsilon_2_0E_0E_fit: bool
        bleak_epsilon_2_0E_1E: float
        bleak_epsilon_2_0E_1E_fit: bool
        bleak_epsilon_2_0E_2E: float
        bleak_epsilon_2_0E_2E_fit: bool
        bleak_epsilon_2_0T_0E: float
        bleak_epsilon_2_0T_0E_fit: bool
        bleak_epsilon_2_0T_1E: float
        bleak_epsilon_2_0T_1E_fit: bool
        bleak_epsilon_2_0T_2E: float
        bleak_epsilon_2_0T_2E_fit: bool
        bleak_epsilon_2_1E_1E: float
        bleak_epsilon_2_1E_1E_fit: bool
        bleak_epsilon_2_1E_2E: float
        bleak_epsilon_2_1E_2E_fit: bool
        bleak_epsilon_2_1T_1E: float
        bleak_epsilon_2_1T_1E_fit: bool
        bleak_epsilon_2_1T_2E: float
        bleak_epsilon_2_1T_2E_fit: bool
        bleak_epsilon_2_2E_2E: float
        bleak_epsilon_2_2E_2E_fit: bool
        bleak_epsilon_2_2T_2E: float
        bleak_epsilon_2_2T_2E_fit: bool
        bleak_epsilon_3_0E_0E: float
        bleak_epsilon_3_0E_0E_fit: bool
        bleak_epsilon_3_0E_1E: float
        bleak_epsilon_3_0E_1E_fit: bool
        bleak_epsilon_3_0E_2E: float
        bleak_epsilon_3_0E_2E_fit: bool
        bleak_epsilon_3_0T_0E: float
        bleak_epsilon_3_0T_0E_fit: bool
        bleak_epsilon_3_0T_1E: float
        bleak_epsilon_3_0T_1E_fit: bool
        bleak_epsilon_3_0T_2E: float
        bleak_epsilon_3_0T_2E_fit: bool
        bleak_epsilon_3_1E_1E: float
        bleak_epsilon_3_1E_1E_fit: bool
        bleak_epsilon_3_1E_2E: float
        bleak_epsilon_3_1E_2E_fit: bool
        bleak_epsilon_3_1T_1E: float
        bleak_epsilon_3_1T_1E_fit: bool
        bleak_epsilon_3_1T_2E: float
        bleak_epsilon_3_1T_2E_fit: bool
        bleak_epsilon_3_2E_2E: float
        bleak_epsilon_3_2E_2E_fit: bool
        bleak_epsilon_3_2T_2E: float
        bleak_epsilon_3_2T_2E_fit: bool
        bleak_epsilon_4_0E_0E: float
        bleak_epsilon_4_0E_0E_fit: bool
        bleak_epsilon_4_0E_1E: float
        bleak_epsilon_4_0E_1E_fit: bool
        bleak_epsilon_4_0E_2E: float
        bleak_epsilon_4_0E_2E_fit: bool
        bleak_epsilon_4_0T_0E: float
        bleak_epsilon_4_0T_0E_fit: bool
        bleak_epsilon_4_0T_1E: float
        bleak_epsilon_4_0T_1E_fit: bool
        bleak_epsilon_4_0T_2E: float
        bleak_epsilon_4_0T_2E_fit: bool
        bleak_epsilon_4_1E_1E: float
        bleak_epsilon_4_1E_1E_fit: bool
        bleak_epsilon_4_1E_2E: float
        bleak_epsilon_4_1E_2E_fit: bool
        bleak_epsilon_4_1T_1E: float
        bleak_epsilon_4_1T_1E_fit: bool
        bleak_epsilon_4_1T_2E: float
        bleak_epsilon_4_1T_2E_fit: bool
        bleak_epsilon_4_2E_2E: float
        bleak_epsilon_4_2E_2E_fit: bool
        bleak_epsilon_4_2T_2E: float
        bleak_epsilon_4_2T_2E_fit: bool
        calib_100P: float
        calib_100P_fit: bool
        calib_143P: float
        calib_143P_fit: bool
        calib_217P: float
        calib_217P_fit: bool
        galf_EE_A_100: float
        galf_EE_A_100_143: float
        galf_EE_A_100_143_fit: bool
        galf_EE_A_100_217: float
        galf_EE_A_100_217_fit: bool
        galf_EE_A_100_fit: bool
        galf_EE_A_143: float
        galf_EE_A_143_217: float
        galf_EE_A_143_217_fit: bool
        galf_EE_A_143_fit: bool
        galf_EE_A_217: float
        galf_EE_A_217_fit: bool
        galf_EE_index: float
        galf_EE_index_fit: bool
        galf_TE_A_100: float
        galf_TE_A_100_143: float
        galf_TE_A_100_143_fit: bool
        galf_TE_A_100_217: float
        galf_TE_A_100_217_fit: bool
        galf_TE_A_100_fit: bool
        galf_TE_A_143: float
        galf_TE_A_143_217: float
        galf_TE_A_143_217_fit: bool
        galf_TE_A_143_fit: bool
        galf_TE_A_217: float
        galf_TE_A_217_fit: bool
        galf_TE_index: float
        galf_TE_index_fit: bool
        A_cib_217: float
        A_cib_217_fit: bool
        A_planck: float
        A_planck_fit: bool
        A_sbpx_100_100_TT: float
        A_sbpx_100_100_TT_fit: bool
        A_sbpx_143_143_TT: float
        A_sbpx_143_143_TT_fit: bool
        A_sbpx_143_217_TT: float
        A_sbpx_143_217_TT_fit: bool
        A_sbpx_217_217_TT: float
        A_sbpx_217_217_TT_fit: bool
        A_sz: float
        A_sz_fit: bool
        calib_100T: float
        calib_100T_fit: bool
        calib_217T: float
        calib_217T_fit: bool
        cib_index: float
        cib_index_fit: bool
        gal545_A_100: float
        gal545_A_100_fit: bool
        gal545_A_143: float
        gal545_A_143_217: float
        gal545_A_143_217_fit: bool
        gal545_A_143_fit: bool
        gal545_A_217: float
        gal545_A_217_fit: bool
        ksz_norm: float
        ksz_norm_fit: bool
        ps_A_100_100: float
        ps_A_100_100_fit: bool
        ps_A_143_143: float
        ps_A_143_143_fit: bool
        ps_A_143_217: float
        ps_A_143_217_fit: bool
        ps_A_217_217: float
        ps_A_217_217_fit: bool
        xi_sz_cib: float
        xi_sz_cib_fit: bool
        version: int
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: PlanckFICorTT = ...
    def __init__(
        self,
        A_cnoise_e2e_100_100_EE: float = ...,
        A_cnoise_e2e_100_100_EE_fit: bool = ...,
        A_cnoise_e2e_143_143_EE: float = ...,
        A_cnoise_e2e_143_143_EE_fit: bool = ...,
        A_cnoise_e2e_217_217_EE: float = ...,
        A_cnoise_e2e_217_217_EE_fit: bool = ...,
        A_pol: float = ...,
        A_pol_fit: bool = ...,
        A_sbpx_100_100_EE: float = ...,
        A_sbpx_100_100_EE_fit: bool = ...,
        A_sbpx_100_143_EE: float = ...,
        A_sbpx_100_143_EE_fit: bool = ...,
        A_sbpx_100_217_EE: float = ...,
        A_sbpx_100_217_EE_fit: bool = ...,
        A_sbpx_143_143_EE: float = ...,
        A_sbpx_143_143_EE_fit: bool = ...,
        A_sbpx_143_217_EE: float = ...,
        A_sbpx_143_217_EE_fit: bool = ...,
        A_sbpx_217_217_EE: float = ...,
        A_sbpx_217_217_EE_fit: bool = ...,
        bleak_epsilon_0_0E_0E: float = ...,
        bleak_epsilon_0_0E_0E_fit: bool = ...,
        bleak_epsilon_0_0E_1E: float = ...,
        bleak_epsilon_0_0E_1E_fit: bool = ...,
        bleak_epsilon_0_0E_2E: float = ...,
        bleak_epsilon_0_0E_2E_fit: bool = ...,
        bleak_epsilon_0_0T_0E: float = ...,
        bleak_epsilon_0_0T_0E_fit: bool = ...,
        bleak_epsilon_0_0T_1E: float = ...,
        bleak_epsilon_0_0T_1E_fit: bool = ...,
        bleak_epsilon_0_0T_2E: float = ...,
        bleak_epsilon_0_0T_2E_fit: bool = ...,
        bleak_epsilon_0_1E_1E: float = ...,
        bleak_epsilon_0_1E_1E_fit: bool = ...,
        bleak_epsilon_0_1E_2E: float = ...,
        bleak_epsilon_0_1E_2E_fit: bool = ...,
        bleak_epsilon_0_1T_1E: float = ...,
        bleak_epsilon_0_1T_1E_fit: bool = ...,
        bleak_epsilon_0_1T_2E: float = ...,
        bleak_epsilon_0_1T_2E_fit: bool = ...,
        bleak_epsilon_0_2E_2E: float = ...,
        bleak_epsilon_0_2E_2E_fit: bool = ...,
        bleak_epsilon_0_2T_2E: float = ...,
        bleak_epsilon_0_2T_2E_fit: bool = ...,
        bleak_epsilon_1_0E_0E: float = ...,
        bleak_epsilon_1_0E_0E_fit: bool = ...,
        bleak_epsilon_1_0E_1E: float = ...,
        bleak_epsilon_1_0E_1E_fit: bool = ...,
        bleak_epsilon_1_0E_2E: float = ...,
        bleak_epsilon_1_0E_2E_fit: bool = ...,
        bleak_epsilon_1_0T_0E: float = ...,
        bleak_epsilon_1_0T_0E_fit: bool = ...,
        bleak_epsilon_1_0T_1E: float = ...,
        bleak_epsilon_1_0T_1E_fit: bool = ...,
        bleak_epsilon_1_0T_2E: float = ...,
        bleak_epsilon_1_0T_2E_fit: bool = ...,
        bleak_epsilon_1_1E_1E: float = ...,
        bleak_epsilon_1_1E_1E_fit: bool = ...,
        bleak_epsilon_1_1E_2E: float = ...,
        bleak_epsilon_1_1E_2E_fit: bool = ...,
        bleak_epsilon_1_1T_1E: float = ...,
        bleak_epsilon_1_1T_1E_fit: bool = ...,
        bleak_epsilon_1_1T_2E: float = ...,
        bleak_epsilon_1_1T_2E_fit: bool = ...,
        bleak_epsilon_1_2E_2E: float = ...,
        bleak_epsilon_1_2E_2E_fit: bool = ...,
        bleak_epsilon_1_2T_2E: float = ...,
        bleak_epsilon_1_2T_2E_fit: bool = ...,
        bleak_epsilon_2_0E_0E: float = ...,
        bleak_epsilon_2_0E_0E_fit: bool = ...,
        bleak_epsilon_2_0E_1E: float = ...,
        bleak_epsilon_2_0E_1E_fit: bool = ...,
        bleak_epsilon_2_0E_2E: float = ...,
        bleak_epsilon_2_0E_2E_fit: bool = ...,
        bleak_epsilon_2_0T_0E: float = ...,
        bleak_epsilon_2_0T_0E_fit: bool = ...,
        bleak_epsilon_2_0T_1E: float = ...,
        bleak_epsilon_2_0T_1E_fit: bool = ...,
        bleak_epsilon_2_0T_2E: float = ...,
        bleak_epsilon_2_0T_2E_fit: bool = ...,
        bleak_epsilon_2_1E_1E: float = ...,
        bleak_epsilon_2_1E_1E_fit: bool = ...,
        bleak_epsilon_2_1E_2E: float = ...,
        bleak_epsilon_2_1E_2E_fit: bool = ...,
        bleak_epsilon_2_1T_1E: float = ...,
        bleak_epsilon_2_1T_1E_fit: bool = ...,
        bleak_epsilon_2_1T_2E: float = ...,
        bleak_epsilon_2_1T_2E_fit: bool = ...,
        bleak_epsilon_2_2E_2E: float = ...,
        bleak_epsilon_2_2E_2E_fit: bool = ...,
        bleak_epsilon_2_2T_2E: float = ...,
        bleak_epsilon_2_2T_2E_fit: bool = ...,
        bleak_epsilon_3_0E_0E: float = ...,
        bleak_epsilon_3_0E_0E_fit: bool = ...,
        bleak_epsilon_3_0E_1E: float = ...,
        bleak_epsilon_3_0E_1E_fit: bool = ...,
        bleak_epsilon_3_0E_2E: float = ...,
        bleak_epsilon_3_0E_2E_fit: bool = ...,
        bleak_epsilon_3_0T_0E: float = ...,
        bleak_epsilon_3_0T_0E_fit: bool = ...,
        bleak_epsilon_3_0T_1E: float = ...,
        bleak_epsilon_3_0T_1E_fit: bool = ...,
        bleak_epsilon_3_0T_2E: float = ...,
        bleak_epsilon_3_0T_2E_fit: bool = ...,
        bleak_epsilon_3_1E_1E: float = ...,
        bleak_epsilon_3_1E_1E_fit: bool = ...,
        bleak_epsilon_3_1E_2E: float = ...,
        bleak_epsilon_3_1E_2E_fit: bool = ...,
        bleak_epsilon_3_1T_1E: float = ...,
        bleak_epsilon_3_1T_1E_fit: bool = ...,
        bleak_epsilon_3_1T_2E: float = ...,
        bleak_epsilon_3_1T_2E_fit: bool = ...,
        bleak_epsilon_3_2E_2E: float = ...,
        bleak_epsilon_3_2E_2E_fit: bool = ...,
        bleak_epsilon_3_2T_2E: float = ...,
        bleak_epsilon_3_2T_2E_fit: bool = ...,
        bleak_epsilon_4_0E_0E: float = ...,
        bleak_epsilon_4_0E_0E_fit: bool = ...,
        bleak_epsilon_4_0E_1E: float = ...,
        bleak_epsilon_4_0E_1E_fit: bool = ...,
        bleak_epsilon_4_0E_2E: float = ...,
        bleak_epsilon_4_0E_2E_fit: bool = ...,
        bleak_epsilon_4_0T_0E: float = ...,
        bleak_epsilon_4_0T_0E_fit: bool = ...,
        bleak_epsilon_4_0T_1E: float = ...,
        bleak_epsilon_4_0T_1E_fit: bool = ...,
        bleak_epsilon_4_0T_2E: float = ...,
        bleak_epsilon_4_0T_2E_fit: bool = ...,
        bleak_epsilon_4_1E_1E: float = ...,
        bleak_epsilon_4_1E_1E_fit: bool = ...,
        bleak_epsilon_4_1E_2E: float = ...,
        bleak_epsilon_4_1E_2E_fit: bool = ...,
        bleak_epsilon_4_1T_1E: float = ...,
        bleak_epsilon_4_1T_1E_fit: bool = ...,
        bleak_epsilon_4_1T_2E: float = ...,
        bleak_epsilon_4_1T_2E_fit: bool = ...,
        bleak_epsilon_4_2E_2E: float = ...,
        bleak_epsilon_4_2E_2E_fit: bool = ...,
        bleak_epsilon_4_2T_2E: float = ...,
        bleak_epsilon_4_2T_2E_fit: bool = ...,
        calib_100P: float = ...,
        calib_100P_fit: bool = ...,
        calib_143P: float = ...,
        calib_143P_fit: bool = ...,
        calib_217P: float = ...,
        calib_217P_fit: bool = ...,
        galf_EE_A_100: float = ...,
        galf_EE_A_100_143: float = ...,
        galf_EE_A_100_143_fit: bool = ...,
        galf_EE_A_100_217: float = ...,
        galf_EE_A_100_217_fit: bool = ...,
        galf_EE_A_100_fit: bool = ...,
        galf_EE_A_143: float = ...,
        galf_EE_A_143_217: float = ...,
        galf_EE_A_143_217_fit: bool = ...,
        galf_EE_A_143_fit: bool = ...,
        galf_EE_A_217: float = ...,
        galf_EE_A_217_fit: bool = ...,
        galf_EE_index: float = ...,
        galf_EE_index_fit: bool = ...,
        galf_TE_A_100: float = ...,
        galf_TE_A_100_143: float = ...,
        galf_TE_A_100_143_fit: bool = ...,
        galf_TE_A_100_217: float = ...,
        galf_TE_A_100_217_fit: bool = ...,
        galf_TE_A_100_fit: bool = ...,
        galf_TE_A_143: float = ...,
        galf_TE_A_143_217: float = ...,
        galf_TE_A_143_217_fit: bool = ...,
        galf_TE_A_143_fit: bool = ...,
        galf_TE_A_217: float = ...,
        galf_TE_A_217_fit: bool = ...,
        galf_TE_index: float = ...,
        galf_TE_index_fit: bool = ...,
        A_cib_217: float = ...,
        A_cib_217_fit: bool = ...,
        A_planck: float = ...,
        A_planck_fit: bool = ...,
        A_sbpx_100_100_TT: float = ...,
        A_sbpx_100_100_TT_fit: bool = ...,
        A_sbpx_143_143_TT: float = ...,
        A_sbpx_143_143_TT_fit: bool = ...,
        A_sbpx_143_217_TT: float = ...,
        A_sbpx_143_217_TT_fit: bool = ...,
        A_sbpx_217_217_TT: float = ...,
        A_sbpx_217_217_TT_fit: bool = ...,
        A_sz: float = ...,
        A_sz_fit: bool = ...,
        calib_100T: float = ...,
        calib_100T_fit: bool = ...,
        calib_217T: float = ...,
        calib_217T_fit: bool = ...,
        cib_index: float = ...,
        cib_index_fit: bool = ...,
        gal545_A_100: float = ...,
        gal545_A_100_fit: bool = ...,
        gal545_A_143: float = ...,
        gal545_A_143_217: float = ...,
        gal545_A_143_217_fit: bool = ...,
        gal545_A_143_fit: bool = ...,
        gal545_A_217: float = ...,
        gal545_A_217_fit: bool = ...,
        ksz_norm: float = ...,
        ksz_norm_fit: bool = ...,
        ps_A_100_100: float = ...,
        ps_A_100_100_fit: bool = ...,
        ps_A_143_143: float = ...,
        ps_A_143_143_fit: bool = ...,
        ps_A_143_217: float = ...,
        ps_A_143_217_fit: bool = ...,
        ps_A_217_217: float = ...,
        ps_A_217_217_fit: bool = ...,
        xi_sz_cib: float = ...,
        xi_sz_cib_fit: bool = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def add_all_default18_priors(lh: NumCosmoMath.Likelihood) -> None: ...
    @staticmethod
    def add_all_default_priors(lh: NumCosmoMath.Likelihood) -> None: ...
    @staticmethod
    def add_calib_priors(
        lh: NumCosmoMath.Likelihood,
        mean: NumCosmoMath.Vector,
        sigma: NumCosmoMath.Vector,
    ) -> None: ...
    @staticmethod
    def add_default18_galf_priors(lh: NumCosmoMath.Likelihood) -> None: ...
    @staticmethod
    def add_default_galf_priors(lh: NumCosmoMath.Likelihood) -> None: ...
    @staticmethod
    def add_galf_priors(
        lh: NumCosmoMath.Likelihood,
        mean: NumCosmoMath.Vector,
        sigma: NumCosmoMath.Vector,
    ) -> None: ...
    @staticmethod
    def add_sz_prior(
        lh: NumCosmoMath.Likelihood, f_tSZ: float, mean: float, sigma: float
    ) -> None: ...

class PlanckFICorTTTEEEClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        PlanckFICorTTTEEEClass()
    """

    parent_class: PlanckFICorTTClass = ...

class PowspecML(NumCosmoMath.Powspec):
    r"""
    :Constructors:

    ::

        PowspecML(**properties)

    Object NcPowspecML

    Properties from NcPowspecML:
      zi -> gdouble: zi
        Initial redshift
      zf -> gdouble: zf
        Final redshift
      kmin -> gdouble: kmin
        Minimum mode value
      kmax -> gdouble: kmax
        Maximum mode value

    Properties from NcmPowspec:
      zi -> gdouble: zi
        Initial time
      zf -> gdouble: zf
        Final time
      kmin -> gdouble: kmin
        Minimum mode value
      kmax -> gdouble: kmax
        Maximum mode value
      reltol -> gdouble: reltol
        Relative tolerance on the interpolation error

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        kmax: float
        kmin: float
        zf: float
        zi: float
        reltol: float

    props: Props = ...
    parent_instance: NumCosmoMath.Powspec = ...
    def __init__(
        self,
        kmax: float = ...,
        kmin: float = ...,
        zf: float = ...,
        zi: float = ...,
        reltol: float = ...,
    ) -> None: ...
    @staticmethod
    def clear(ps_ml: PowspecML) -> None: ...
    def free(self) -> None: ...
    def ref(self) -> PowspecML: ...

class PowspecMLCBE(PowspecML):
    r"""
    :Constructors:

    ::

        PowspecMLCBE(**properties)
        new() -> NumCosmo.PowspecMLCBE
        new_full(cbe:NumCosmo.CBE) -> NumCosmo.PowspecMLCBE

    Object NcPowspecMLCBE

    Properties from NcPowspecMLCBE:
      cbe -> NcCBE: cbe
        Class backend object
      intern-k-min -> gdouble: intern-k-min
        Class minimum mode k
      intern-k-max -> gdouble: intern-k-max
        Class maximum mode k

    Properties from NcPowspecML:
      zi -> gdouble: zi
        Initial redshift
      zf -> gdouble: zf
        Final redshift
      kmin -> gdouble: kmin
        Minimum mode value
      kmax -> gdouble: kmax
        Maximum mode value

    Properties from NcmPowspec:
      zi -> gdouble: zi
        Initial time
      zf -> gdouble: zf
        Final time
      kmin -> gdouble: kmin
        Minimum mode value
      kmax -> gdouble: kmax
        Maximum mode value
      reltol -> gdouble: reltol
        Relative tolerance on the interpolation error

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        cbe: CBE
        intern_k_max: float
        intern_k_min: float
        kmax: float
        kmin: float
        zf: float
        zi: float
        reltol: float

    props: Props = ...
    parent_instance: PowspecML = ...
    priv: PowspecMLCBEPrivate = ...
    def __init__(
        self,
        cbe: CBE = ...,
        intern_k_max: float = ...,
        intern_k_min: float = ...,
        kmax: float = ...,
        kmin: float = ...,
        zf: float = ...,
        zi: float = ...,
        reltol: float = ...,
    ) -> None: ...
    def get_intern_k_max(self) -> float: ...
    def get_intern_k_min(self) -> float: ...
    @classmethod
    def new(cls) -> PowspecMLCBE: ...
    @classmethod
    def new_full(cls, cbe: CBE) -> PowspecMLCBE: ...
    def peek_cbe(self) -> CBE: ...
    def set_cbe(self, cbe: CBE) -> None: ...
    def set_intern_k_max(self, k_max: float) -> None: ...
    def set_intern_k_min(self, k_min: float) -> None: ...

class PowspecMLCBEClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        PowspecMLCBEClass()
    """

    parent_class: PowspecMLClass = ...

class PowspecMLCBEPrivate(GObject.GPointer): ...

class PowspecMLClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        PowspecMLClass()
    """

    parent_class: NumCosmoMath.PowspecClass = ...
    padding: list[None] = ...

class PowspecMLSpline(PowspecML):
    r"""
    :Constructors:

    ::

        PowspecMLSpline(**properties)
        new(Pk:NumCosmoMath.Spline) -> NumCosmo.PowspecMLSpline

    Object NcPowspecMLSpline

    Properties from NcPowspecMLSpline:
      spline -> NcmSpline: spline
        Spline representing Pk at z=0

    Properties from NcPowspecML:
      zi -> gdouble: zi
        Initial redshift
      zf -> gdouble: zf
        Final redshift
      kmin -> gdouble: kmin
        Minimum mode value
      kmax -> gdouble: kmax
        Maximum mode value

    Properties from NcmPowspec:
      zi -> gdouble: zi
        Initial time
      zf -> gdouble: zf
        Final time
      kmin -> gdouble: kmin
        Minimum mode value
      kmax -> gdouble: kmax
        Maximum mode value
      reltol -> gdouble: reltol
        Relative tolerance on the interpolation error

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        spline: NumCosmoMath.Spline
        kmax: float
        kmin: float
        zf: float
        zi: float
        reltol: float

    props: Props = ...
    def __init__(
        self,
        spline: NumCosmoMath.Spline = ...,
        kmax: float = ...,
        kmin: float = ...,
        zf: float = ...,
        zi: float = ...,
        reltol: float = ...,
    ) -> None: ...
    @classmethod
    def new(cls, Pk: NumCosmoMath.Spline) -> PowspecMLSpline: ...
    def peek_spline(self) -> NumCosmoMath.Spline: ...
    def set_spline(self, Pk: NumCosmoMath.Spline) -> None: ...

class PowspecMLSplineClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        PowspecMLSplineClass()
    """

    parent_class: PowspecMLClass = ...

class PowspecMLTransfer(PowspecML):
    r"""
    :Constructors:

    ::

        PowspecMLTransfer(**properties)
        new(tf:NumCosmo.TransferFunc) -> NumCosmo.PowspecMLTransfer

    Object NcPowspecMLTransfer

    Properties from NcPowspecMLTransfer:
      transfer -> NcTransferFunc: transfer
        Transfer function
      growth -> NcGrowthFunc: growth
        Growth function

    Properties from NcPowspecML:
      zi -> gdouble: zi
        Initial redshift
      zf -> gdouble: zf
        Final redshift
      kmin -> gdouble: kmin
        Minimum mode value
      kmax -> gdouble: kmax
        Maximum mode value

    Properties from NcmPowspec:
      zi -> gdouble: zi
        Initial time
      zf -> gdouble: zf
        Final time
      kmin -> gdouble: kmin
        Minimum mode value
      kmax -> gdouble: kmax
        Maximum mode value
      reltol -> gdouble: reltol
        Relative tolerance on the interpolation error

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        growth: GrowthFunc
        transfer: TransferFunc
        kmax: float
        kmin: float
        zf: float
        zi: float
        reltol: float

    props: Props = ...
    def __init__(
        self,
        growth: GrowthFunc = ...,
        transfer: TransferFunc = ...,
        kmax: float = ...,
        kmin: float = ...,
        zf: float = ...,
        zi: float = ...,
        reltol: float = ...,
    ) -> None: ...
    @classmethod
    def new(cls, tf: TransferFunc) -> PowspecMLTransfer: ...
    def peek_gf(self) -> GrowthFunc: ...
    def peek_tf(self) -> TransferFunc: ...
    def set_gf(self, gf: GrowthFunc) -> None: ...
    def set_tf(self, tf: TransferFunc) -> None: ...

class PowspecMLTransferClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        PowspecMLTransferClass()
    """

    parent_class: PowspecMLClass = ...

class PowspecMNL(NumCosmoMath.Powspec):
    r"""
    :Constructors:

    ::

        PowspecMNL(**properties)

    Object NcPowspecMNL

    Properties from NcmPowspec:
      zi -> gdouble: zi
        Initial time
      zf -> gdouble: zf
        Final time
      kmin -> gdouble: kmin
        Minimum mode value
      kmax -> gdouble: kmax
        Maximum mode value
      reltol -> gdouble: reltol
        Relative tolerance on the interpolation error

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        kmax: float
        kmin: float
        reltol: float
        zf: float
        zi: float

    props: Props = ...
    parent_instance: NumCosmoMath.Powspec = ...
    def __init__(
        self,
        kmax: float = ...,
        kmin: float = ...,
        reltol: float = ...,
        zf: float = ...,
        zi: float = ...,
    ) -> None: ...
    @staticmethod
    def clear(ps_mnl: PowspecMNL) -> None: ...
    def free(self) -> None: ...
    def ref(self) -> PowspecMNL: ...

class PowspecMNLClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        PowspecMNLClass()
    """

    parent_class: NumCosmoMath.PowspecClass = ...

class PowspecMNLHaloFit(PowspecMNL):
    r"""
    :Constructors:

    ::

        PowspecMNLHaloFit(**properties)
        new(psml:NumCosmo.PowspecML, zmaxnl:float, reltol:float) -> NumCosmo.PowspecMNLHaloFit

    Object NcPowspecMNLHaloFit

    Properties from NcPowspecMNLHaloFit:
      power-spec -> NcPowspecML: power-spec
        Linear power spectrum.
      zmaxnl -> gdouble: zmaxnl
        Max redshift for halofit correction
      reltol -> gdouble: reltol
        Relative tolerance (precision) for halofit computations
      use-pkequal -> gboolean: use-pkequal
        Whether to use PKEqual

    Properties from NcmPowspec:
      zi -> gdouble: zi
        Initial time
      zf -> gdouble: zf
        Final time
      kmin -> gdouble: kmin
        Minimum mode value
      kmax -> gdouble: kmax
        Maximum mode value
      reltol -> gdouble: reltol
        Relative tolerance on the interpolation error

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        power_spec: PowspecML
        reltol: float
        use_pkequal: bool
        zmaxnl: float
        kmax: float
        kmin: float
        zf: float
        zi: float

    props: Props = ...
    parent_instance: PowspecMNL = ...
    priv: PowspecMNLHaloFitPrivate = ...
    def __init__(
        self,
        power_spec: PowspecML = ...,
        reltol: float = ...,
        use_pkequal: bool = ...,
        zmaxnl: float = ...,
        kmax: float = ...,
        kmin: float = ...,
        zf: float = ...,
        zi: float = ...,
    ) -> None: ...
    @classmethod
    def new(
        cls, psml: PowspecML, zmaxnl: float, reltol: float
    ) -> PowspecMNLHaloFit: ...
    def pkequal(self, on: bool) -> None: ...
    def set_kbounds_from_ml(self) -> None: ...

class PowspecMNLHaloFitClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        PowspecMNLHaloFitClass()
    """

    parent_class: PowspecMNLClass = ...

class PowspecMNLHaloFitPrivate(GObject.GPointer): ...

class PriorQSplineCont(NumCosmoMath.Prior):
    r"""
    :Constructors:

    ::

        PriorQSplineCont(**properties)
        new() -> NumCosmo.PriorQSplineCont

    Object NcPriorQSplineCont

    Properties from NcmMSetFunc:
      nvariables -> guint: nvariables
        Number of variables
      dimension -> guint: dimension
        Function dimension
      eval-x -> NcmVector: eval-x
        Evaluation point x

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        dimension: int
        eval_x: NumCosmoMath.Vector
        nvariables: int

    props: Props = ...
    parent_instance: NumCosmoMath.Prior = ...
    def __init__(
        self,
        dimension: int = ...,
        eval_x: NumCosmoMath.Vector = ...,
        nvariables: int = ...,
    ) -> None: ...
    @classmethod
    def new(cls) -> PriorQSplineCont: ...

class PriorQSplineContClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        PriorQSplineContClass()
    """

    parent_class: NumCosmoMath.PriorClass = ...

class Recomb(GObject.Object):
    r"""
    :Constructors:

    ::

        Recomb(**properties)

    Object NcRecomb

    Properties from NcRecomb:
      zi -> gdouble: zi
        Initial redshift for recombination calculations
      init-frac -> gdouble: init-frac
        Initial fraction to start numerical integration
      prec -> gdouble: prec
        Precision for recombination calculations

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        init_frac: float
        prec: float
        zi: float

    props: Props = ...
    parent_instance: GObject.Object = ...
    zi: float = ...
    lambdai: float = ...
    lambdaf: float = ...
    prec: float = ...
    init_frac: float = ...
    fmin: int = ...
    fsol: int = ...
    dtau_dlambda_s: NumCosmoMath.Spline = ...
    tau_s: NumCosmoMath.Spline = ...
    tau_ode_s: NumCosmoMath.OdeSpline = ...
    tau_drag_ode_s: NumCosmoMath.OdeSpline = ...
    ctrl_cosmo: NumCosmoMath.ModelCtrl = ...
    ctrl_reion: NumCosmoMath.ModelCtrl = ...
    v_tau_max_z: float = ...
    v_tau_max_lambda: float = ...
    tau_z: float = ...
    tau_lambda: float = ...
    tau_drag_z: float = ...
    tau_drag_lambda: float = ...
    tau_cutoff_z: float = ...
    tau_cutoff_lambda: float = ...
    def __init__(
        self, init_frac: float = ..., prec: float = ..., zi: float = ...
    ) -> None: ...
    @staticmethod
    def HI_ion_saha(cosmo: HICosmo, x: float) -> float: ...
    @staticmethod
    def HeII_ion_saha(cosmo: HICosmo, x: float) -> float: ...
    @staticmethod
    def HeII_ion_saha_x(cosmo: HICosmo, f: float) -> float: ...
    @staticmethod
    def HeII_ion_saha_x_by_HeIII_He(cosmo: HICosmo, f: float) -> float: ...
    @staticmethod
    def HeI_ion_saha(cosmo: HICosmo, x: float) -> float: ...
    @staticmethod
    def He_fully_ionized_Xe(cosmo: HICosmo, x: float) -> float: ...
    @staticmethod
    def He_fully_ionized_dtau_dlambda(cosmo: HICosmo, lambda_: float) -> float: ...
    def XHII(self, cosmo: HICosmo, lambda_: float) -> float: ...
    def XHeII(self, cosmo: HICosmo, lambda_: float) -> float: ...
    def Xe(self, cosmo: HICosmo, lambda_: float) -> float: ...
    @staticmethod
    def clear(recomb: Recomb) -> None: ...
    def d2tau_dlambda2(self, cosmo: HICosmo, lambda_: float) -> float: ...
    def d2v_tau_dlambda2(self, cosmo: HICosmo, lambda_: float) -> float: ...
    def d3tau_dlambda3(self, cosmo: HICosmo, lambda_: float) -> float: ...
    def do_XHII(self, cosmo: HICosmo, lambda_: float) -> float: ...
    def do_XHeII(self, cosmo: HICosmo, lambda_: float) -> float: ...
    def do_Xe(self, cosmo: HICosmo, lambda_: float) -> float: ...
    def do_prepare(self, cosmo: HICosmo) -> None: ...
    def dtau_dlambda(self, cosmo: HICosmo, lambda_: float) -> float: ...
    @staticmethod
    def dtau_dlambda_Xe(cosmo: HICosmo, lambda_: float) -> float: ...
    def dtau_dx(self, cosmo: HICosmo, lambda_: float) -> float: ...
    def dv_tau_dlambda(self, cosmo: HICosmo, lambda_: float) -> float: ...
    def equilibrium_XHI(self, cosmo: HICosmo, x: float) -> float: ...
    def equilibrium_XHII(self, cosmo: HICosmo, x: float) -> float: ...
    def equilibrium_XHeI(self, cosmo: HICosmo, x: float) -> float: ...
    def equilibrium_XHeII(self, cosmo: HICosmo, x: float) -> float: ...
    def equilibrium_XHeIII(self, cosmo: HICosmo, x: float) -> float: ...
    def equilibrium_Xe(self, cosmo: HICosmo, x: float) -> float: ...
    def free(self) -> None: ...
    def get_tau_cutoff_lambda(self, cosmo: HICosmo) -> float: ...
    def get_tau_cutoff_z(self, cosmo: HICosmo) -> float: ...
    def get_tau_drag_lambda(self, cosmo: HICosmo) -> float: ...
    def get_tau_drag_z(self, cosmo: HICosmo) -> float: ...
    def get_tau_lambda(self, cosmo: HICosmo) -> float: ...
    def get_tau_z(self, cosmo: HICosmo) -> float: ...
    def get_v_tau_max_lambda(self, cosmo: HICosmo) -> float: ...
    def get_v_tau_max_z(self, cosmo: HICosmo) -> float: ...
    def get_zi(self) -> float: ...
    def log_v_tau(self, cosmo: HICosmo, lambda_: float) -> float: ...
    def prepare(self, cosmo: HICosmo) -> None: ...
    def prepare_if_needed(self, cosmo: HICosmo) -> None: ...
    def ref(self) -> Recomb: ...
    def require_zi(self, zi: float) -> None: ...
    def set_zi(self, zi: float) -> None: ...
    def tau(self, cosmo: HICosmo, lambda_: float) -> float: ...
    def tau_drag(self, cosmo: HICosmo, lambda_: float) -> float: ...
    def tau_lambda0_lambda1(
        self, cosmo: HICosmo, lambda0: float, lambda1: float
    ) -> float: ...
    def v_tau(self, cosmo: HICosmo, lambda_: float) -> float: ...
    def v_tau_lambda_features(
        self, cosmo: HICosmo, logref: float
    ) -> typing.Tuple[float, float, float]: ...

class RecombCBE(Recomb):
    r"""
    :Constructors:

    ::

        RecombCBE(**properties)
        full_new(cbe:NumCosmo.CBE) -> NumCosmo.RecombCBE
        new() -> NumCosmo.RecombCBE

    Object NcRecombCBE

    Properties from NcRecombCBE:
      cbe -> NcCBE: cbe
        Class backend

    Properties from NcRecomb:
      zi -> gdouble: zi
        Initial redshift for recombination calculations
      init-frac -> gdouble: init-frac
        Initial fraction to start numerical integration
      prec -> gdouble: prec
        Precision for recombination calculations

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        cbe: CBE
        init_frac: float
        prec: float
        zi: float

    props: Props = ...
    parent_instance: Recomb = ...
    cbe: CBE = ...
    Xe_s: NumCosmoMath.Spline = ...
    def __init__(
        self, cbe: CBE = ..., init_frac: float = ..., prec: float = ..., zi: float = ...
    ) -> None: ...
    @staticmethod
    def clear(recomb_cbe: RecombCBE) -> None: ...
    def free(self) -> None: ...
    @classmethod
    def full_new(cls, cbe: CBE) -> RecombCBE: ...
    @classmethod
    def new(cls) -> RecombCBE: ...
    def peek_cbe(self) -> CBE: ...
    def ref(self) -> RecombCBE: ...
    def set_cbe(self, cbe: CBE) -> None: ...

class RecombCBEClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        RecombCBEClass()
    """

    parent_class: RecombClass = ...

class RecombClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        RecombClass()
    """

    parent_class: GObject.ObjectClass = ...
    prepare: typing.Callable[[Recomb, HICosmo], None] = ...
    Xe: typing.Callable[[Recomb, HICosmo, float], float] = ...
    XHII: typing.Callable[[Recomb, HICosmo, float], float] = ...
    XHeII: typing.Callable[[Recomb, HICosmo, float], float] = ...

class RecombSeager(Recomb):
    r"""
    :Constructors:

    ::

        RecombSeager(**properties)
        new() -> NumCosmo.RecombSeager
        new_full(init_frac:float, zi:float, prec:float) -> NumCosmo.RecombSeager

    Object NcRecombSeager

    Properties from NcRecombSeager:
      options -> NcRecombSeagerOpt: options
        Integration options

    Properties from NcRecomb:
      zi -> gdouble: zi
        Initial redshift for recombination calculations
      init-frac -> gdouble: init-frac
        Initial fraction to start numerical integration
      prec -> gdouble: prec
        Precision for recombination calculations

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        options: RecombSeagerOpt
        init_frac: float
        prec: float
        zi: float

    props: Props = ...
    parent_instance: Recomb = ...
    priv: RecombSeagerPrivate = ...
    def __init__(
        self,
        options: RecombSeagerOpt = ...,
        init_frac: float = ...,
        prec: float = ...,
        zi: float = ...,
    ) -> None: ...
    @staticmethod
    def clear(recomb_seager: RecombSeager) -> None: ...
    def free(self) -> None: ...
    def get_options(self) -> RecombSeagerOpt: ...
    def hummer_HeI_case_B(self, cosmo: HICosmo, Tm: float) -> float: ...
    def hummer_HeI_case_B_dTm(self, cosmo: HICosmo, Tm: float) -> float: ...
    def hummer_HeI_case_B_trip(self, cosmo: HICosmo, Tm: float) -> float: ...
    def hummer_HeI_case_B_trip_dTm(self, cosmo: HICosmo, Tm: float) -> float: ...
    @classmethod
    def new(cls) -> RecombSeager: ...
    @classmethod
    def new_full(cls, init_frac: float, zi: float, prec: float) -> RecombSeager: ...
    def pequignot_HI_case_B(self, cosmo: HICosmo, Tm: float) -> float: ...
    def pequignot_HI_case_B_dTm(self, cosmo: HICosmo, Tm: float) -> float: ...
    def ref(self) -> RecombSeager: ...
    def set_options(self, opts: RecombSeagerOpt) -> None: ...
    def set_switch(self, H_switch: int, He_switch: int) -> None: ...

class RecombSeagerClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        RecombSeagerClass()
    """

    parent_class: RecombClass = ...

class RecombSeagerPrivate(GObject.GPointer): ...

class ReducedShearCalib(NumCosmoMath.Model):
    r"""
    :Constructors:

    ::

        ReducedShearCalib(**properties)

    Object NcReducedShearCalib

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: NumCosmoMath.Model = ...
    priv: ReducedShearCalibPrivate = ...
    def __init__(
        self,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def clear(rs_calib: ReducedShearCalib) -> None: ...
    def do_eval(self, g_th: float, psf_size: float, gal_size: float) -> float: ...
    def eval(self, g_th: float, psf_size: float, gal_size: float) -> float: ...
    def free(self) -> None: ...
    @staticmethod
    def id() -> int: ...
    def ref(self) -> ReducedShearCalib: ...

class ReducedShearCalibClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        ReducedShearCalibClass()
    """

    parent_class: NumCosmoMath.ModelClass = ...
    eval: typing.Callable[[ReducedShearCalib, float, float, float], float] = ...

class ReducedShearCalibPrivate(GObject.GPointer): ...

class ReducedShearCalibWtg(ReducedShearCalib):
    r"""
    :Constructors:

    ::

        ReducedShearCalibWtg(**properties)
        new() -> NumCosmo.ReducedShearCalibWtg

    Object NcReducedShearCalibWtg

    Properties from NcReducedShearCalibWtg:
      mslope -> gdouble: mslope
        m_s
      mb -> gdouble: mb
        m_b
      c -> gdouble: c
        c
      xp -> gdouble: xp
        x_p
      mslope-fit -> gboolean: mslope-fit
        m_s:fit
      mb-fit -> gboolean: mb-fit
        m_b:fit
      c-fit -> gboolean: c-fit
        c:fit
      xp-fit -> gboolean: xp-fit
        x_p:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        c: float
        c_fit: bool
        mb: float
        mb_fit: bool
        mslope: float
        mslope_fit: bool
        xp: float
        xp_fit: bool
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: ReducedShearCalib = ...
    priv: ReducedShearCalibWtgPrivate = ...
    def __init__(
        self,
        c: float = ...,
        c_fit: bool = ...,
        mb: float = ...,
        mb_fit: bool = ...,
        mslope: float = ...,
        mslope_fit: bool = ...,
        xp: float = ...,
        xp_fit: bool = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def clear(rs_wtg: ReducedShearCalibWtg) -> None: ...
    def free(self) -> None: ...
    @classmethod
    def new(cls) -> ReducedShearCalibWtg: ...
    def ref(self) -> ReducedShearCalibWtg: ...

class ReducedShearCalibWtgClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        ReducedShearCalibWtgClass()
    """

    parent_class: ReducedShearCalibClass = ...

class ReducedShearCalibWtgPrivate(GObject.GPointer): ...

class ReducedShearClusterMass(NumCosmoMath.Model):
    r"""
    :Constructors:

    ::

        ReducedShearClusterMass(**properties)
        new() -> NumCosmo.ReducedShearClusterMass

    Object NcReducedShearClusterMass

    Properties from NcReducedShearClusterMass:
      R -> gdouble: R
        Distance from the center of the lens
      number-z-bins -> guint: number-z-bins
        Number of redshift bins
      a -> gdouble: a
        a
      b -> gdouble: b
        b
      c -> gdouble: c
        c
      xp -> gdouble: xp
        xp
      sigma -> gdouble: sigma
        \sigma
      Gamma -> gdouble: Gamma
        \Gamma
      a-fit -> gboolean: a-fit
        a:fit
      b-fit -> gboolean: b-fit
        b:fit
      c-fit -> gboolean: c-fit
        c:fit
      xp-fit -> gboolean: xp-fit
        xp:fit
      sigma-fit -> gboolean: sigma-fit
        \sigma:fit
      Gamma-fit -> gboolean: Gamma-fit
        \Gamma:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        Gamma: float
        Gamma_fit: bool
        R: float
        a: float
        a_fit: bool
        b: float
        b_fit: bool
        c: float
        c_fit: bool
        number_z_bins: int
        sigma: float
        sigma_fit: bool
        xp: float
        xp_fit: bool
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: NumCosmoMath.Model = ...
    R_Mpc: float = ...
    nzbins: float = ...
    T: int = ...
    s: int = ...
    workz: float = ...
    def __init__(
        self,
        Gamma: float = ...,
        Gamma_fit: bool = ...,
        R: float = ...,
        a: float = ...,
        a_fit: bool = ...,
        b: float = ...,
        b_fit: bool = ...,
        c: float = ...,
        c_fit: bool = ...,
        number_z_bins: int = ...,
        sigma: float = ...,
        sigma_fit: bool = ...,
        xp: float = ...,
        xp_fit: bool = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    def P_z_gth_gobs(
        self, cosmo: HICosmo, z: float, g_th: float, g_obs: float
    ) -> float: ...
    @staticmethod
    def clear(rscm: ReducedShearClusterMass) -> None: ...
    def free(self) -> None: ...
    @staticmethod
    def id() -> int: ...
    @classmethod
    def new(cls) -> ReducedShearClusterMass: ...
    def posterior_no_shear_calibration(
        self, cosmo: HICosmo, z: float, g_obs: float
    ) -> float: ...
    def ref(self) -> ReducedShearClusterMass: ...

class ReducedShearClusterMassClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        ReducedShearClusterMassClass()
    """

    parent_class: NumCosmoMath.ModelClass = ...

class SNIADistCov(NumCosmoMath.Model):
    r"""
    :Constructors:

    ::

        SNIADistCov(**properties)
        new(dist:NumCosmo.Distance, sigma_int_len:int) -> NumCosmo.SNIADistCov
        new_by_id(dist:NumCosmo.Distance, snia_id:NumCosmo.DataSNIAId) -> NumCosmo.SNIADistCov

    Object NcSNIADistCov

    Properties from NcSNIADistCov:
      dist -> NcDistance: dist
        Distance object
      empty-fac -> gboolean: empty-fac
        Empty universe approximation factor
      alpha -> gdouble: alpha
        \alpha
      beta -> gdouble: beta
        \beta
      M1 -> gdouble: M1
        \mathcal{M}_1
      M2 -> gdouble: M2
        \mathcal{M}_2
      lnsigma-pecz -> gdouble: lnsigma-pecz
        \ln(\sigma_{\mathrm{pecz}})
      lnsigma-lens -> gdouble: lnsigma-lens
        \ln(\sigma_{\mathrm{lens}})
      lnsigma-int -> NcmVector: lnsigma-int
        \ln(\sigma_{\mathrm{int}})
      mu -> NcmVector: mu
        \mu
      lnsigma-int-length -> guint: lnsigma-int-length
        \ln(\sigma_{\mathrm{int}}):length
      mu-length -> guint: mu-length
        \mu:length
      alpha-fit -> gboolean: alpha-fit
        \alpha:fit
      beta-fit -> gboolean: beta-fit
        \beta:fit
      M1-fit -> gboolean: M1-fit
        \mathcal{M}_1:fit
      M2-fit -> gboolean: M2-fit
        \mathcal{M}_2:fit
      lnsigma-pecz-fit -> gboolean: lnsigma-pecz-fit
        \ln(\sigma_{\mathrm{pecz}}):fit
      lnsigma-lens-fit -> gboolean: lnsigma-lens-fit
        \ln(\sigma_{\mathrm{lens}}):fit
      lnsigma-int-fit -> GVariant: lnsigma-int-fit
        \ln(\sigma_{\mathrm{int}}):fit
      mu-fit -> GVariant: mu-fit
        \mu:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        M1: float
        M1_fit: bool
        M2: float
        M2_fit: bool
        alpha: float
        alpha_fit: bool
        beta: float
        beta_fit: bool
        dist: Distance
        empty_fac: bool
        lnsigma_int: NumCosmoMath.Vector
        lnsigma_int_fit: GLib.Variant
        lnsigma_int_length: int
        lnsigma_lens: float
        lnsigma_lens_fit: bool
        lnsigma_pecz: float
        lnsigma_pecz_fit: bool
        mu: NumCosmoMath.Vector
        mu_fit: GLib.Variant
        mu_length: int
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: NumCosmoMath.Model = ...
    dist: Distance = ...
    var_int: list[None] = ...
    empty_fac: bool = ...
    cov_cpu: None = ...
    alpha_cpu: float = ...
    beta_cpu: float = ...
    lnsigma_pecz_cpu: float = ...
    lnsigma_lens_cpu: float = ...
    def __init__(
        self,
        M1: float = ...,
        M1_fit: bool = ...,
        M2: float = ...,
        M2_fit: bool = ...,
        alpha: float = ...,
        alpha_fit: bool = ...,
        beta: float = ...,
        beta_fit: bool = ...,
        dist: Distance = ...,
        empty_fac: bool = ...,
        lnsigma_int: NumCosmoMath.Vector = ...,
        lnsigma_int_fit: GLib.Variant = ...,
        lnsigma_int_length: int = ...,
        lnsigma_lens: float = ...,
        lnsigma_lens_fit: bool = ...,
        lnsigma_pecz: float = ...,
        lnsigma_pecz_fit: bool = ...,
        mu: NumCosmoMath.Vector = ...,
        mu_fit: GLib.Variant = ...,
        mu_length: int = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    def alpha_beta(self) -> typing.Tuple[float, float]: ...
    def calc(self, snia_cov: DataSNIACov, cov: NumCosmoMath.Matrix) -> bool: ...
    @staticmethod
    def clear(dcov: SNIADistCov) -> None: ...
    def extra_var(self, snia_cov: DataSNIACov, i: int) -> float: ...
    def free(self) -> None: ...
    @staticmethod
    def id() -> int: ...
    def mag(
        self,
        cosmo: HICosmo,
        snia_cov: DataSNIACov,
        i: int,
        width_th: float,
        colour_th: float,
    ) -> float: ...
    def mag_to_width_colour(
        self,
        cosmo: HICosmo,
        snia_cov: DataSNIACov,
        obs: NumCosmoMath.Vector,
        X: NumCosmoMath.Matrix,
        colmajor: bool,
    ) -> None: ...
    def mean(
        self, cosmo: HICosmo, snia_cov: DataSNIACov, y: NumCosmoMath.Vector
    ) -> None: ...
    def mean_V2(
        self, cosmo: HICosmo, snia_cov: DataSNIACov, y: NumCosmoMath.Vector
    ) -> None: ...
    @classmethod
    def new(cls, dist: Distance, sigma_int_len: int) -> SNIADistCov: ...
    @classmethod
    def new_by_id(cls, dist: Distance, snia_id: DataSNIAId) -> SNIADistCov: ...
    def prepare(self, mset: NumCosmoMath.MSet) -> None: ...
    def prepare_if_needed(self, mset: NumCosmoMath.MSet) -> None: ...
    def ref(self) -> SNIADistCov: ...
    def set_dist(self, dist: Distance) -> None: ...
    def set_empty_fac(self, enable: bool) -> None: ...

class SNIADistCovClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        SNIADistCovClass()
    """

    parent_class: NumCosmoMath.ModelClass = ...

class Scalefactor(GObject.Object):
    r"""
    :Constructors:

    ::

        Scalefactor(**properties)
        new(zf:float, dist:NumCosmo.Distance) -> NumCosmo.Scalefactor

    Object NcScalefactor

    Properties from NcScalefactor:
      zf -> gdouble: zf
        Initial redshift
      a0 -> gdouble: a0
        Scale factor today a_0
      a0-conformal-normal -> gboolean: a0-conformal-normal
        Scale factor today a_0 from normalized curvature radius
      dist -> NcDistance: dist
        Distance object
      reltol -> gdouble: reltol
        Relative tolerance
      abstol -> gdouble: abstol
        Absolute tolerance

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        a0: float
        a0_conformal_normal: bool
        abstol: float
        dist: Distance
        reltol: float
        zf: float

    props: Props = ...
    parent_instance: GObject.Object = ...
    priv: ScalefactorPrivate = ...
    def __init__(
        self,
        a0: float = ...,
        a0_conformal_normal: bool = ...,
        abstol: float = ...,
        dist: Distance = ...,
        reltol: float = ...,
        zf: float = ...,
    ) -> None: ...
    @staticmethod
    def clear(a: Scalefactor) -> None: ...
    def eval_a_eta(self, eta: float) -> float: ...
    def eval_eta_t(self, t: float) -> float: ...
    def eval_eta_x(self, x: float) -> float: ...
    def eval_eta_z(self, z: float) -> float: ...
    def eval_t_eta(self, eta: float) -> float: ...
    def eval_z_eta(self, eta: float) -> float: ...
    def free(self) -> None: ...
    def get_a0(self) -> float: ...
    def get_abstol(self) -> float: ...
    def get_reltol(self) -> float: ...
    def get_zf(self) -> float: ...
    @classmethod
    def new(cls, zf: float, dist: Distance) -> Scalefactor: ...
    def prepare(self, cosmo: HICosmo) -> None: ...
    def prepare_if_needed(self, cosmo: HICosmo) -> None: ...
    def ref(self) -> Scalefactor: ...
    def require_zf(self, zf: float) -> None: ...
    def set_a0(self, a0: float) -> None: ...
    def set_a0_conformal_normal(self, enable: bool) -> None: ...
    def set_abstol(self, abstol: float) -> None: ...
    def set_reltol(self, reltol: float) -> None: ...
    def set_zf(self, zf: float) -> None: ...

class ScalefactorClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        ScalefactorClass()
    """

    parent_class: GObject.ObjectClass = ...

class ScalefactorPrivate(GObject.GPointer): ...

class TransferFunc(GObject.Object):
    r"""
    :Constructors:

    ::

        TransferFunc(**properties)

    Object NcTransferFunc

    Signals from GObject:
      notify (GParam)
    """

    parent_instance: GObject.Object = ...
    ctrl_cosmo: NumCosmoMath.ModelCtrl = ...
    ctrl_reion: NumCosmoMath.ModelCtrl = ...
    @staticmethod
    def clear(tf: TransferFunc) -> None: ...
    def do_calc(self, k: float) -> float: ...
    def do_prepare(self, cosmo: HICosmo) -> None: ...
    def eval(self, cosmo: HICosmo, kh: float) -> float: ...
    def free(self) -> None: ...
    def prepare(self, cosmo: HICosmo) -> None: ...
    def prepare_if_needed(self, cosmo: HICosmo) -> None: ...
    def ref(self) -> TransferFunc: ...

class TransferFuncBBKS(TransferFunc):
    r"""
    :Constructors:

    ::

        TransferFuncBBKS(**properties)
        new() -> NumCosmo.TransferFunc

    Object NcTransferFuncBBKS

    Properties from NcTransferFuncBBKS:
      type -> NcTransferFuncBBKSType: type
        BBKS variant type

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        type: TransferFuncBBKSType

    props: Props = ...
    parent_instance: TransferFunc = ...
    priv: TransferFuncBBKSPrivate = ...
    def __init__(self, type: TransferFuncBBKSType = ...) -> None: ...
    @classmethod
    def new(cls) -> TransferFuncBBKS: ...
    def set_type(self, bbks_type: TransferFuncBBKSType) -> None: ...

class TransferFuncBBKSClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        TransferFuncBBKSClass()
    """

    parent_class: TransferFuncClass = ...

class TransferFuncBBKSPrivate(GObject.GPointer): ...

class TransferFuncCAMB(TransferFunc):
    r"""
    :Constructors:

    ::

        TransferFuncCAMB(**properties)
        new() -> NumCosmo.TransferFunc

    Object NcTransferFuncCAMB

    Signals from GObject:
      notify (GParam)
    """

    parent_instance: TransferFunc = ...
    T_spline: NumCosmoMath.Spline = ...
    init: bool = ...
    @classmethod
    def new(cls) -> TransferFuncCAMB: ...

class TransferFuncCAMBClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        TransferFuncCAMBClass()
    """

    parent_class: TransferFuncClass = ...

class TransferFuncClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        TransferFuncClass()
    """

    parent_class: GObject.ObjectClass = ...
    alloc: typing.Callable[[], None] = ...
    prepare: typing.Callable[[TransferFunc, HICosmo], None] = ...
    calc: typing.Callable[[TransferFunc, float], float] = ...

class TransferFuncEH(TransferFunc):
    r"""
    :Constructors:

    ::

        TransferFuncEH(**properties)
        new() -> NumCosmo.TransferFunc

    Object NcTransferFuncEH

    Properties from NcTransferFuncEH:
      CCL-comp -> gboolean: CCL-comp
        Whether to use CCL compatible mode

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        CCL_comp: bool

    props: Props = ...
    parent_instance: TransferFunc = ...
    priv: TransferFuncEHPrivate = ...
    def __init__(self, CCL_comp: bool = ...) -> None: ...
    @classmethod
    def new(cls) -> TransferFuncEH: ...
    def set_CCL_comp(self, CCL_comp: bool) -> None: ...

class TransferFuncEHClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        TransferFuncEHClass()
    """

    parent_class: TransferFuncClass = ...

class TransferFuncEHPrivate(GObject.GPointer): ...

class WLSurfaceMassDensity(NumCosmoMath.Model):
    r"""
    :Constructors:

    ::

        WLSurfaceMassDensity(**properties)
        new(dist:NumCosmo.Distance) -> NumCosmo.WLSurfaceMassDensity

    Object NcWLSurfaceMassDensity

    Properties from NcWLSurfaceMassDensity:
      distance -> NcDistance: distance
        Distance
      pcc -> gdouble: pcc
        p_{cc}
      Roff -> gdouble: Roff
        R_{off}
      pcc-fit -> gboolean: pcc-fit
        p_{cc}:fit
      Roff-fit -> gboolean: Roff-fit
        R_{off}:fit

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        Roff: float
        Roff_fit: bool
        distance: Distance
        pcc: float
        pcc_fit: bool
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: NumCosmoMath.Model = ...
    dist: Distance = ...
    ctrl_cosmo: NumCosmoMath.ModelCtrl = ...
    ctrl_dp: NumCosmoMath.ModelCtrl = ...
    def __init__(
        self,
        Roff: float = ...,
        Roff_fit: bool = ...,
        distance: Distance = ...,
        pcc: float = ...,
        pcc_fit: bool = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @staticmethod
    def clear(smd: WLSurfaceMassDensity) -> None: ...
    def convergence(
        self,
        dp: HaloDensityProfile,
        cosmo: HICosmo,
        R: float,
        zs: float,
        zl: float,
        zc: float,
    ) -> float: ...
    def convergence_infinity(
        self, dp: HaloDensityProfile, cosmo: HICosmo, R: float, zl: float, zc: float
    ) -> float: ...
    def free(self) -> None: ...
    @staticmethod
    def id() -> int: ...
    def magnification(
        self,
        dp: HaloDensityProfile,
        cosmo: HICosmo,
        R: float,
        zs: float,
        zl: float,
        zc: float,
    ) -> float: ...
    @classmethod
    def new(cls, dist: Distance) -> WLSurfaceMassDensity: ...
    def prepare(self, cosmo: HICosmo) -> None: ...
    def prepare_if_needed(self, cosmo: HICosmo) -> None: ...
    def reduced_shear(
        self,
        dp: HaloDensityProfile,
        cosmo: HICosmo,
        R: float,
        zs: float,
        zl: float,
        zc: float,
    ) -> float: ...
    def reduced_shear_array(
        self,
        dp: HaloDensityProfile,
        cosmo: HICosmo,
        R: typing.Sequence[float],
        fin: float,
        fout: float,
        zs: typing.Sequence[float],
        zl: float,
        zc: float,
    ) -> list[float]: ...
    def reduced_shear_array_equal(
        self,
        dp: HaloDensityProfile,
        cosmo: HICosmo,
        R: typing.Sequence[float],
        fin: float,
        fout: float,
        zs: typing.Sequence[float],
        zl: float,
        zc: float,
    ) -> list[float]: ...
    def reduced_shear_infinity(
        self,
        dp: HaloDensityProfile,
        cosmo: HICosmo,
        R: float,
        zs: float,
        zl: float,
        zc: float,
    ) -> float: ...
    def reduced_shear_optzs(
        self,
        dp: HaloDensityProfile,
        cosmo: HICosmo,
        zs: float,
        zl: float,
        optzs: WLSurfaceMassDensityOptzs,
    ) -> float: ...
    def ref(self) -> WLSurfaceMassDensity: ...
    def shear(
        self,
        dp: HaloDensityProfile,
        cosmo: HICosmo,
        R: float,
        zs: float,
        zl: float,
        zc: float,
    ) -> float: ...
    def shear_infinity(
        self, dp: HaloDensityProfile, cosmo: HICosmo, R: float, zl: float, zc: float
    ) -> float: ...
    def sigma(
        self, dp: HaloDensityProfile, cosmo: HICosmo, R: float, zc: float
    ) -> float: ...
    def sigma_array(
        self,
        dp: HaloDensityProfile,
        cosmo: HICosmo,
        R: typing.Sequence[float],
        fin: float,
        fout: float,
        zc: float,
    ) -> list[float]: ...
    def sigma_critical(
        self, cosmo: HICosmo, zs: float, zl: float, zc: float
    ) -> float: ...
    def sigma_critical_infinity(
        self, cosmo: HICosmo, zl: float, zc: float
    ) -> float: ...
    def sigma_excess(
        self, dp: HaloDensityProfile, cosmo: HICosmo, R: float, zc: float
    ) -> float: ...
    def sigma_excess_array(
        self,
        dp: HaloDensityProfile,
        cosmo: HICosmo,
        R: typing.Sequence[float],
        fin: float,
        fout: float,
        zc: float,
    ) -> list[float]: ...
    def sigma_mean(
        self, dp: HaloDensityProfile, cosmo: HICosmo, R: float, zc: float
    ) -> float: ...

class WLSurfaceMassDensityClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        WLSurfaceMassDensityClass()
    """

    parent_class: NumCosmoMath.ModelClass = ...

class WLSurfaceMassDensityOptzs(GObject.GPointer):
    r"""
    :Constructors:

    ::

        WLSurfaceMassDensityOptzs()
    """

    k: int = ...
    sqrt_Omega_k0: float = ...
    dl: float = ...
    sc_Dls_Ds: float = ...
    sigma: float = ...
    mean_sigma: float = ...

class Window(GObject.Object):
    r"""
    :Constructors:

    ::

        Window(**properties)

    Object NcWindow

    Signals from GObject:
      notify (GParam)
    """

    parent_instance: GObject.Object = ...
    @staticmethod
    def clear(wf: Window) -> None: ...
    def deriv_fourier(self, k: float, R: float) -> float: ...
    def do_deriv_fourier(self, k: float, R: float) -> float: ...
    def do_eval_fourier(self, k: float, R: float) -> float: ...
    def do_eval_real(self, r: float, R: float) -> float: ...
    def eval_fourier(self, k: float, R: float) -> float: ...
    def eval_realspace(self, r: float, R: float) -> float: ...
    def free(self) -> None: ...
    def volume(self) -> float: ...

class WindowClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        WindowClass()
    """

    parent_class: GObject.ObjectClass = ...
    volume: float = ...
    eval_fourier: typing.Callable[[Window, float, float], float] = ...
    deriv_fourier: typing.Callable[[Window, float, float], float] = ...
    eval_real: typing.Callable[[Window, float, float], float] = ...

class WindowGaussian(Window):
    r"""
    :Constructors:

    ::

        WindowGaussian(**properties)
        new() -> NumCosmo.Window

    Object NcWindowGaussian

    Signals from GObject:
      notify (GParam)
    """

    parent_instance: Window = ...
    @classmethod
    def new(cls) -> WindowGaussian: ...

class WindowGaussianClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        WindowGaussianClass()
    """

    parent_class: WindowClass = ...

class WindowTophat(Window):
    r"""
    :Constructors:

    ::

        WindowTophat(**properties)
        new() -> NumCosmo.Window

    Object NcWindowTophat

    Signals from GObject:
      notify (GParam)
    """

    parent_instance: Window = ...
    @classmethod
    def new(cls) -> WindowTophat: ...

class WindowTophatClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        WindowTophatClass()
    """

    parent_class: WindowClass = ...

class Xcor(GObject.Object):
    r"""
    :Constructors:

    ::

        Xcor(**properties)
        new(dist:NumCosmo.Distance, ps:NumCosmoMath.Powspec, meth:NumCosmo.XcorLimberMethod) -> NumCosmo.Xcor

    Object NcXcor

    Properties from NcXcor:
      distance -> NcDistance: distance
        Distance.
      power-spec -> NcmPowspec: power-spec
        Matter power spectrum.
      meth -> NcXcorLimberMethod: meth
        Method.
      reltol -> gdouble: reltol
        Relative tolerance.

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        distance: Distance
        meth: XcorLimberMethod
        power_spec: NumCosmoMath.Powspec
        reltol: float

    props: Props = ...
    def __init__(
        self,
        distance: Distance = ...,
        meth: XcorLimberMethod = ...,
        power_spec: NumCosmoMath.Powspec = ...,
        reltol: float = ...,
    ) -> None: ...
    @staticmethod
    def clear(xc: Xcor) -> None: ...
    def free(self) -> None: ...
    def get_reltol(self) -> float: ...
    def limber(
        self,
        xclk1: XcorLimberKernel,
        xclk2: XcorLimberKernel,
        cosmo: HICosmo,
        lmin: int,
        lmax: int,
        vp: NumCosmoMath.Vector,
    ) -> None: ...
    @classmethod
    def new(
        cls, dist: Distance, ps: NumCosmoMath.Powspec, meth: XcorLimberMethod
    ) -> Xcor: ...
    def prepare(self, cosmo: HICosmo) -> None: ...
    def ref(self) -> Xcor: ...
    def set_reltol(self, reltol: float) -> None: ...

class XcorAB(GObject.Object):
    r"""
    :Constructors:

    ::

        XcorAB(**properties)
        new(a:int, b:int, ell_th_cut_off:int, ell_lik_min:int, ell_lik_max:int, clobs_filename:str=None, mixing_filename:str=None, mixing_filelength:int) -> NumCosmo.XcorAB

    Object NcXcorAB

    Properties from NcXcorAB:
      a -> guint: a
        a
      b -> guint: b
        b
      ell-th-cut-off -> guint: ell-th-cut-off
        ell_th_cut_off
      ell-lik-min -> guint: ell-lik-min
        ell_lik_min
      ell-lik-max -> guint: ell-lik-max
        ell_lik_max
      mixing -> NcmMatrix: mixing
        mixing
      cl-th -> NcmMatrix: cl-th
        cl_th
      cl-obs -> NcmVector: cl-obs
        cl_obs

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        a: int
        b: int
        cl_obs: NumCosmoMath.Vector
        cl_th: NumCosmoMath.Matrix
        ell_lik_max: int
        ell_lik_min: int
        ell_th_cut_off: int
        mixing: NumCosmoMath.Matrix

    props: Props = ...
    parent_instance: GObject.Object = ...
    a: int = ...
    b: int = ...
    ell_th_cut_off: int = ...
    ell_lik_min: int = ...
    ell_lik_max: int = ...
    nell_lik: int = ...
    mixing: NumCosmoMath.Matrix = ...
    cl_th: NumCosmoMath.Matrix = ...
    cl_obs: NumCosmoMath.Vector = ...
    def __init__(
        self,
        a: int = ...,
        b: int = ...,
        cl_obs: NumCosmoMath.Vector = ...,
        cl_th: NumCosmoMath.Matrix = ...,
        ell_lik_max: int = ...,
        ell_lik_min: int = ...,
        ell_th_cut_off: int = ...,
        mixing: NumCosmoMath.Matrix = ...,
    ) -> None: ...
    @staticmethod
    def clear(xcab: XcorAB) -> None: ...
    def free(self) -> None: ...
    @classmethod
    def new(
        cls,
        a: int,
        b: int,
        ell_th_cut_off: int,
        ell_lik_min: int,
        ell_lik_max: int,
        clobs_filename: typing.Optional[str],
        mixing_filename: typing.Optional[str],
        mixing_filelength: int,
    ) -> XcorAB: ...
    def ref(self) -> XcorAB: ...

class XcorABClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        XcorABClass()
    """

    parent_class: GObject.ObjectClass = ...

class XcorClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        XcorClass()
    """

    parent_class: GObject.ObjectClass = ...

class XcorKinetic(GObject.GBoxed):
    r"""
    :Constructors:

    ::

        XcorKinetic()
    """

    xi_z: float = ...
    E_z: float = ...
    def copy(self) -> XcorKinetic: ...
    def free(self) -> None: ...

class XcorLimberKernel(NumCosmoMath.Model):
    r"""
    :Constructors:

    ::

        XcorLimberKernel(**properties)

    Object NcXcorLimberKernel

    Properties from NcXcorLimberKernel:
      zmin -> gdouble: zmin
        Minimum redshift
      zmax -> gdouble: zmax
        Maximum redshift

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        zmax: float
        zmin: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    parent_instance: NumCosmoMath.Model = ...
    def __init__(
        self,
        zmax: float = ...,
        zmin: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    def add_noise(
        self, vp1: NumCosmoMath.Vector, vp2: NumCosmoMath.Vector, lmin: int
    ) -> None: ...
    @staticmethod
    def clear(xclk: XcorLimberKernel) -> None: ...
    def do_add_noise(
        self, vp1: NumCosmoMath.Vector, vp2: NumCosmoMath.Vector, lmin: int
    ) -> None: ...
    def do_eval(self, cosmo: HICosmo, z: float, xck: XcorKinetic, l: int) -> float: ...
    def do_obs_len(self) -> int: ...
    def do_obs_params_len(self) -> int: ...
    def do_prepare(self, cosmo: HICosmo) -> None: ...
    def eval(self, cosmo: HICosmo, z: float, xck: XcorKinetic, l: int) -> float: ...
    def eval_full(self, cosmo: HICosmo, z: float, dist: Distance, l: int) -> float: ...
    def free(self) -> None: ...
    def get_const_factor(self) -> float: ...
    def get_z_range(self) -> typing.Tuple[float, float, float]: ...
    @staticmethod
    def id() -> int: ...
    @staticmethod
    def log_all_models() -> None: ...
    def obs_len(self) -> int: ...
    def obs_params_len(self) -> int: ...
    def prepare(self, cosmo: HICosmo) -> None: ...
    def ref(self) -> XcorLimberKernel: ...
    def set_const_factor(self, cf: float) -> None: ...
    def set_z_range(self, zmin: float, zmax: float, zmid: float) -> None: ...

class XcorLimberKernelCMBISW(XcorLimberKernel):
    r"""
    :Constructors:

    ::

        XcorLimberKernelCMBISW(**properties)
        new(dist:NumCosmo.Distance, ps:NumCosmoMath.Powspec, recomb:NumCosmo.Recomb, Nl:NumCosmoMath.Vector) -> NumCosmo.XcorLimberKernelCMBISW

    Object NcXcorLimberKernelCMBISW

    Properties from NcXcorLimberKernelCMBISW:
      dist -> NcDistance: dist
        Distance object
      ps -> NcmPowspec: ps
        Power Spectrum object
      recomb -> NcRecomb: recomb
        Recombination object
      Nl -> NcmVector: Nl
        Noise spectrum

    Properties from NcXcorLimberKernel:
      zmin -> gdouble: zmin
        Minimum redshift
      zmax -> gdouble: zmax
        Maximum redshift

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        Nl: NumCosmoMath.Vector
        dist: Distance
        ps: NumCosmoMath.Powspec
        recomb: Recomb
        zmax: float
        zmin: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    def __init__(
        self,
        Nl: NumCosmoMath.Vector = ...,
        dist: Distance = ...,
        ps: NumCosmoMath.Powspec = ...,
        recomb: Recomb = ...,
        zmax: float = ...,
        zmin: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @classmethod
    def new(
        cls,
        dist: Distance,
        ps: NumCosmoMath.Powspec,
        recomb: Recomb,
        Nl: NumCosmoMath.Vector,
    ) -> XcorLimberKernelCMBISW: ...

class XcorLimberKernelCMBISWClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        XcorLimberKernelCMBISWClass()
    """

    parent_class: XcorLimberKernelClass = ...

class XcorLimberKernelCMBLensing(XcorLimberKernel):
    r"""
    :Constructors:

    ::

        XcorLimberKernelCMBLensing(**properties)
        new(dist:NumCosmo.Distance, recomb:NumCosmo.Recomb, Nl:NumCosmoMath.Vector) -> NumCosmo.XcorLimberKernelCMBLensing

    Object NcXcorLimberKernelCMBLensing

    Properties from NcXcorLimberKernelCMBLensing:
      dist -> NcDistance: dist
        Distance object
      recomb -> NcRecomb: recomb
        Recombination object
      Nl -> NcmVector: Nl
        Noise spectrum

    Properties from NcXcorLimberKernel:
      zmin -> gdouble: zmin
        Minimum redshift
      zmax -> gdouble: zmax
        Maximum redshift

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        Nl: NumCosmoMath.Vector
        dist: Distance
        recomb: Recomb
        zmax: float
        zmin: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    def __init__(
        self,
        Nl: NumCosmoMath.Vector = ...,
        dist: Distance = ...,
        recomb: Recomb = ...,
        zmax: float = ...,
        zmin: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @classmethod
    def new(
        cls, dist: Distance, recomb: Recomb, Nl: NumCosmoMath.Vector
    ) -> XcorLimberKernelCMBLensing: ...

class XcorLimberKernelCMBLensingClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        XcorLimberKernelCMBLensingClass()
    """

    parent_class: XcorLimberKernelClass = ...

class XcorLimberKernelClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        XcorLimberKernelClass()
    """

    parent_class: NumCosmoMath.ModelClass = ...
    eval: typing.Callable[
        [XcorLimberKernel, HICosmo, float, XcorKinetic, int], float
    ] = ...
    prepare: typing.Callable[[XcorLimberKernel, HICosmo], None] = ...
    add_noise: typing.Callable[
        [XcorLimberKernel, NumCosmoMath.Vector, NumCosmoMath.Vector, int], None
    ] = ...
    obs_len: typing.Callable[[XcorLimberKernel], int] = ...
    obs_params_len: typing.Callable[[XcorLimberKernel], int] = ...

class XcorLimberKernelGal(XcorLimberKernel):
    r"""
    :Constructors:

    ::

        XcorLimberKernelGal(**properties)
        new(zmin:float, zmax:float, np:int, nbarm1:float, dn_dz:NumCosmoMath.Spline, dist:NumCosmo.Distance, domagbias:bool) -> NumCosmo.XcorLimberKernelGal

    Object NcXcorLimberKernelGal

    Properties from NcXcorLimberKernelGal:
      dndz -> NcmSpline: dndz
        Galaxy redshift distribution
      bias -> NcmSpline: bias
        Bias spline object
      domagbias -> gboolean: domagbias
        Do magnification bias
      nbarm1 -> gdouble: nbarm1
        One over nbar (galaxy angular density)
      dist -> NcDistance: dist
        Distance object
      mag-bias -> gdouble: mag-bias
        mag_bias
      noise-bias -> gdouble: noise-bias
        noise_bias
      bparam -> NcmVector: bparam
        bparam
      bparam-length -> guint: bparam-length
        bparam:length
      mag-bias-fit -> gboolean: mag-bias-fit
        mag_bias:fit
      noise-bias-fit -> gboolean: noise-bias-fit
        noise_bias:fit
      bparam-fit -> GVariant: bparam-fit
        bparam:fit

    Properties from NcXcorLimberKernel:
      zmin -> gdouble: zmin
        Minimum redshift
      zmax -> gdouble: zmax
        Maximum redshift

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        bias: NumCosmoMath.Spline
        bparam: NumCosmoMath.Vector
        bparam_fit: GLib.Variant
        bparam_length: int
        dist: Distance
        dndz: NumCosmoMath.Spline
        domagbias: bool
        mag_bias: float
        mag_bias_fit: bool
        nbarm1: float
        noise_bias: float
        noise_bias_fit: bool
        zmax: float
        zmin: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    def __init__(
        self,
        bias: NumCosmoMath.Spline = ...,
        bparam: NumCosmoMath.Vector = ...,
        bparam_fit: GLib.Variant = ...,
        bparam_length: int = ...,
        dist: Distance = ...,
        dndz: NumCosmoMath.Spline = ...,
        domagbias: bool = ...,
        mag_bias: float = ...,
        mag_bias_fit: bool = ...,
        nbarm1: float = ...,
        noise_bias: float = ...,
        noise_bias_fit: bool = ...,
        zmax: float = ...,
        zmin: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    def get_bias(self) -> typing.Tuple[float, float, float]: ...
    def get_fast_update(self) -> bool: ...
    @classmethod
    def new(
        cls,
        zmin: float,
        zmax: float,
        np: int,
        nbarm1: float,
        dn_dz: NumCosmoMath.Spline,
        dist: Distance,
        domagbias: bool,
    ) -> XcorLimberKernelGal: ...
    def set_bias_old(self, bias_old: float, noise_bias_old: float) -> None: ...
    def set_fast_update(self, fast_update: bool) -> None: ...

class XcorLimberKernelGalClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        XcorLimberKernelGalClass()
    """

    parent_class: XcorLimberKernelClass = ...

class XcorLimberKernelWeakLensing(XcorLimberKernel):
    r"""
    :Constructors:

    ::

        XcorLimberKernelWeakLensing(**properties)
        new(zmin:float, zmax:float, dn_dz:NumCosmoMath.Spline, nbar:float, intr_shear:float, dist:NumCosmo.Distance) -> NumCosmo.XcorLimberKernelWeakLensing

    Object NcXcorLimberKernelWeakLensing

    Properties from NcXcorLimberKernelWeakLensing:
      dndz -> NcmSpline: dndz
        Source redshift distribution
      nbar -> gdouble: nbar
        nbar (galaxy angular density)
      intr-shear -> gdouble: intr-shear
        Intrinsic galaxy shear
      dist -> NcDistance: dist
        Distance object
      reltol -> gdouble: reltol
        Relative tolerance
      abstol -> gdouble: abstol
        Absolute tolerance tolerance

    Properties from NcXcorLimberKernel:
      zmin -> gdouble: zmin
        Minimum redshift
      zmax -> gdouble: zmax
        Maximum redshift

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        abstol: float
        dist: Distance
        dndz: NumCosmoMath.Spline
        intr_shear: float
        nbar: float
        reltol: float
        zmax: float
        zmin: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    def __init__(
        self,
        abstol: float = ...,
        dist: Distance = ...,
        dndz: NumCosmoMath.Spline = ...,
        intr_shear: float = ...,
        nbar: float = ...,
        reltol: float = ...,
        zmax: float = ...,
        zmin: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @classmethod
    def new(
        cls,
        zmin: float,
        zmax: float,
        dn_dz: NumCosmoMath.Spline,
        nbar: float,
        intr_shear: float,
        dist: Distance,
    ) -> XcorLimberKernelWeakLensing: ...

class XcorLimberKernelWeakLensingClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        XcorLimberKernelWeakLensingClass()
    """

    parent_class: XcorLimberKernelClass = ...

class XcorLimberKerneltSZ(XcorLimberKernel):
    r"""
    :Constructors:

    ::

        XcorLimberKerneltSZ(**properties)
        new(zmax:float) -> NumCosmo.XcorLimberKerneltSZ

    Object NcXcorLimberKerneltSZ

    Properties from NcXcorLimberKerneltSZ:
      noise -> gdouble: noise
        Constant noise level

    Properties from NcXcorLimberKernel:
      zmin -> gdouble: zmin
        Minimum redshift
      zmax -> gdouble: zmax
        Maximum redshift

    Properties from NcmModel:
      name -> gchararray: name
        Model's name
      nick -> gchararray: nick
        Model's nick
      scalar-params-len -> guint: scalar-params-len
        Number of scalar parameters
      vector-params-len -> guint: vector-params-len
        Number of vector parameters
      implementation -> guint64: implementation
        Bitwise specification of functions implementation
      sparam-array -> NcmObjDictInt: sparam-array
        NcmModel array of NcmSParam
      params-types -> GArray: params-types
        Parameters' types
      reparam -> NcmReparam: reparam
        Model reparametrization
      submodel-array -> NcmObjArray: submodel-array
        NcmModel array of submodels

    Signals from GObject:
      notify (GParam)
    """

    class Props:
        noise: float
        zmax: float
        zmin: float
        implementation: int
        name: str
        nick: str
        params_types: list[None]
        reparam: NumCosmoMath.Reparam
        scalar_params_len: int
        sparam_array: NumCosmoMath.ObjDictInt
        submodel_array: NumCosmoMath.ObjArray
        vector_params_len: int

    props: Props = ...
    def __init__(
        self,
        noise: float = ...,
        zmax: float = ...,
        zmin: float = ...,
        reparam: NumCosmoMath.Reparam = ...,
        sparam_array: NumCosmoMath.ObjDictInt = ...,
        submodel_array: NumCosmoMath.ObjArray = ...,
    ) -> None: ...
    @classmethod
    def new(cls, zmax: float) -> XcorLimberKerneltSZ: ...

class XcorLimberKerneltSZClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        XcorLimberKerneltSZClass()
    """

    parent_class: XcorLimberKernelClass = ...

class _HaloPositionClass(GObject.GPointer):
    r"""
    :Constructors:

    ::

        _HaloPositionClass()
    """

    parent_class: NumCosmoMath.ModelClass = ...
    padding: list[None] = ...

class DataCMBDataType(GObject.GFlags):
    ALL: DataCMBDataType = ...
    BB: DataCMBDataType = ...
    EB: DataCMBDataType = ...
    EE: DataCMBDataType = ...
    PHIPHI: DataCMBDataType = ...
    TB: DataCMBDataType = ...
    TE: DataCMBDataType = ...
    TT: DataCMBDataType = ...

class DataClusterWLResampleFlag(GObject.GFlags):
    ALL: DataClusterWLResampleFlag = ...
    POSITION: DataClusterWLResampleFlag = ...
    REDSHIFT: DataClusterWLResampleFlag = ...
    SHAPE: DataClusterWLResampleFlag = ...

class HICosmoDEImpl(GObject.GFlags):
    D2E2OMEGA_DE_DZ2: HICosmoDEImpl = ...
    DE2OMEGA_DE_DZ: HICosmoDEImpl = ...
    E2OMEGA_DE: HICosmoDEImpl = ...
    W_DE: HICosmoDEImpl = ...

class HICosmoImpl(GObject.GFlags):
    AS_DRAG: HICosmoImpl = ...
    BGP_CS2: HICosmoImpl = ...
    D2E2_DZ2: HICosmoImpl = ...
    DC: HICosmoImpl = ...
    DE2_DZ: HICosmoImpl = ...
    E2: HICosmoImpl = ...
    E2OMEGA_B: HICosmoImpl = ...
    E2OMEGA_C: HICosmoImpl = ...
    E2OMEGA_G: HICosmoImpl = ...
    E2OMEGA_M: HICosmoImpl = ...
    E2OMEGA_MNU: HICosmoImpl = ...
    E2OMEGA_MNU_N: HICosmoImpl = ...
    E2OMEGA_NU: HICosmoImpl = ...
    E2OMEGA_R: HICosmoImpl = ...
    E2OMEGA_T: HICosmoImpl = ...
    E2PRESS_MNU: HICosmoImpl = ...
    E2PRESS_MNU_N: HICosmoImpl = ...
    GET_BG_VAR: HICosmoImpl = ...
    H0: HICosmoImpl = ...
    MASSNUINFO: HICosmoImpl = ...
    NMASSNU: HICosmoImpl = ...
    OMEGA_B0: HICosmoImpl = ...
    OMEGA_C0: HICosmoImpl = ...
    OMEGA_G0: HICosmoImpl = ...
    OMEGA_M0: HICosmoImpl = ...
    OMEGA_MNU0: HICosmoImpl = ...
    OMEGA_MNU0_N: HICosmoImpl = ...
    OMEGA_NU0: HICosmoImpl = ...
    OMEGA_R0: HICosmoImpl = ...
    OMEGA_T0: HICosmoImpl = ...
    PRESS_MNU0: HICosmoImpl = ...
    PRESS_MNU0_N: HICosmoImpl = ...
    T_GAMMA0: HICosmoImpl = ...
    XB: HICosmoImpl = ...
    YP_4HE: HICosmoImpl = ...
    Z_LSS: HICosmoImpl = ...

class HIPrimImpl(GObject.GFlags):
    LNSA_POWSPEC_LNK: HIPrimImpl = ...
    LNT_POWSPEC_LNK: HIPrimImpl = ...

class RecombSeagerOpt(GObject.GFlags):
    ALL: RecombSeagerOpt = ...
    HEII_SOBOLEV_1P1: RecombSeagerOpt = ...
    HEII_SOBOLEV_1P1_CO: RecombSeagerOpt = ...
    HEII_SOBOLEV_3P012: RecombSeagerOpt = ...
    HEII_SOBOLEV_3P012_CO: RecombSeagerOpt = ...
    HII_FUDGE: RecombSeagerOpt = ...
    HII_FUDGE_GAUSS_COR: RecombSeagerOpt = ...

class ClusterMassAscasoSParams(GObject.GEnum):
    MU_P0: ClusterMassAscasoSParams = ...
    MU_P1: ClusterMassAscasoSParams = ...
    MU_P2: ClusterMassAscasoSParams = ...
    SIGMA_P0: ClusterMassAscasoSParams = ...
    SIGMA_P1: ClusterMassAscasoSParams = ...
    SIGMA_P2: ClusterMassAscasoSParams = ...

class ClusterMassBensonSParams(GObject.GEnum):
    A_SZ: ClusterMassBensonSParams = ...
    B_SZ: ClusterMassBensonSParams = ...
    C_SZ: ClusterMassBensonSParams = ...
    D_SZ: ClusterMassBensonSParams = ...

class ClusterMassBensonXRaySParams(GObject.GEnum):
    A_X: ClusterMassBensonXRaySParams = ...
    B_X: ClusterMassBensonXRaySParams = ...
    C_X: ClusterMassBensonXRaySParams = ...
    D_X: ClusterMassBensonXRaySParams = ...

class ClusterMassImpl(GObject.GEnum):
    INTP: ClusterMassImpl = ...
    N_LIMITS: ClusterMassImpl = ...
    P: ClusterMassImpl = ...
    P_LIMITS: ClusterMassImpl = ...
    RESAMPLE: ClusterMassImpl = ...

class ClusterMassLnnormalSParams(GObject.GEnum):
    BIAS: ClusterMassLnnormalSParams = ...
    SIGMA: ClusterMassLnnormalSParams = ...

class ClusterMassPlCLSParams(GObject.GEnum):
    A_L: ClusterMassPlCLSParams = ...
    A_SZ: ClusterMassPlCLSParams = ...
    B_L: ClusterMassPlCLSParams = ...
    B_SZ: ClusterMassPlCLSParams = ...
    COR: ClusterMassPlCLSParams = ...
    SD_L: ClusterMassPlCLSParams = ...
    SD_SZ: ClusterMassPlCLSParams = ...

class ClusterMassVanderlindeSParams(GObject.GEnum):
    A_SZ: ClusterMassVanderlindeSParams = ...
    B_SZ: ClusterMassVanderlindeSParams = ...
    C_SZ: ClusterMassVanderlindeSParams = ...
    D_SZ: ClusterMassVanderlindeSParams = ...

class ClusterPhotozGaussGlobalSParams(GObject.GEnum):
    SIGMA0: ClusterPhotozGaussGlobalSParams = ...
    Z_BIAS: ClusterPhotozGaussGlobalSParams = ...

class ClusterPseudoCountsSParams(GObject.GEnum):
    DELTAZ: ClusterPseudoCountsSParams = ...
    LNMCUT: ClusterPseudoCountsSParams = ...
    SD_MCUT: ClusterPseudoCountsSParams = ...
    ZMIN: ClusterPseudoCountsSParams = ...

class ClusterRedshiftImpl(GObject.GEnum):
    INTP: ClusterRedshiftImpl = ...
    N_LIMTS: ClusterRedshiftImpl = ...
    P: ClusterRedshiftImpl = ...
    P_LIMITS: ClusterRedshiftImpl = ...
    RESAMPLE: ClusterRedshiftImpl = ...

class DataBaoId(GObject.GEnum):
    A_EISENSTEIN2005: DataBaoId = ...
    DHR_DAR_SDSS_DR11_2015: DataBaoId = ...
    DHR_DAR_SDSS_DR11_2015_LYAF_AUTO_CROSS: DataBaoId = ...
    DMR_HR_SDSS_DR12_2016: DataBaoId = ...
    DTR_DHR_DESI_DR1_LYA_2025: DataBaoId = ...
    DTR_DHR_SDSS_DR12_2016_DR16_COMPATIBLE: DataBaoId = ...
    DTR_DHR_SDSS_DR16_LRG_2021: DataBaoId = ...
    DTR_DHR_SDSS_DR16_QSO_2021: DataBaoId = ...
    DVDV_PERCIVAL2007: DataBaoId = ...
    DVDV_PERCIVAL2010: DataBaoId = ...
    DVR_DTDH_DESI_DR1_2024: DataBaoId = ...
    DVR_DTDH_DESI_DR2_2025: DataBaoId = ...
    DV_EISENSTEIN2005: DataBaoId = ...
    EMPIRICAL_FIT_1D_SDSS_DR16_ELG_2021: DataBaoId = ...
    EMPIRICAL_FIT_2D_BAUTISTA2017: DataBaoId = ...
    EMPIRICAL_FIT_2D_SDSS_DR16_LYAUTO_2021: DataBaoId = ...
    EMPIRICAL_FIT_2D_SDSS_DR16_LYXQSO_2021: DataBaoId = ...
    EMPIRICAL_FIT_ROSS2015: DataBaoId = ...
    RDV_ANDERSON2012: DataBaoId = ...
    RDV_BEUTLER2011: DataBaoId = ...
    RDV_BLAKE2012: DataBaoId = ...
    RDV_BOSS_QSO_ATA2017: DataBaoId = ...
    RDV_DESI_DR1_BGS_QSO_2024: DataBaoId = ...
    RDV_DESI_DR2_BGS_2025: DataBaoId = ...
    RDV_KAZIN2014: DataBaoId = ...
    RDV_PADMANABHAN2012: DataBaoId = ...
    RDV_PERCIVAL2007: DataBaoId = ...
    RDV_PERCIVAL2010: DataBaoId = ...

class DataCMBId(GObject.GEnum):
    DIST_PRIORS_WMAP5: DataCMBId = ...
    DIST_PRIORS_WMAP7: DataCMBId = ...
    DIST_PRIORS_WMAP9: DataCMBId = ...
    SHIFT_PARAM_WMAP3: DataCMBId = ...
    SHIFT_PARAM_WMAP5: DataCMBId = ...
    SHIFT_PARAM_WMAP7: DataCMBId = ...

class DataClusterAbundanceId(GObject.GEnum):
    FIT: DataClusterAbundanceId = ...
    SAMPLING: DataClusterAbundanceId = ...
    TXT: DataClusterAbundanceId = ...

class DataClusterPseudoCountsObs(GObject.GEnum):
    MCL: DataClusterPseudoCountsObs = ...
    MPL: DataClusterPseudoCountsObs = ...
    SD_MCL: DataClusterPseudoCountsObs = ...
    SD_MPL: DataClusterPseudoCountsObs = ...
    Z: DataClusterPseudoCountsObs = ...

class DataHubbleBaoId(GObject.GEnum):
    BUSCA2013: DataHubbleBaoId = ...

class DataHubbleId(GObject.GEnum):
    BORGHI2022: DataHubbleId = ...
    BUSCA2013_BAO_WMAP: DataHubbleId = ...
    CABRE: DataHubbleId = ...
    GOMEZ_VALENT_COMP2018: DataHubbleId = ...
    JIAO2023: DataHubbleId = ...
    JIMENEZ2023: DataHubbleId = ...
    MORESCO2012_BC03: DataHubbleId = ...
    MORESCO2012_MASTRO: DataHubbleId = ...
    MORESCO2015: DataHubbleId = ...
    MORESCO2016_DR9_BC03: DataHubbleId = ...
    MORESCO2016_DR9_MASTRO: DataHubbleId = ...
    RATSIMBAZAFY2017: DataHubbleId = ...
    RIESS2008_HST: DataHubbleId = ...
    RIESS2016_HST_WFC3: DataHubbleId = ...
    RIESS2018: DataHubbleId = ...
    SIMON2005: DataHubbleId = ...
    STERN2009: DataHubbleId = ...
    TOMASETTI2023: DataHubbleId = ...
    ZHANG2012: DataHubbleId = ...

class DataPlanckLKLType(GObject.GEnum):
    BASELINE_18_HIGHL_TT: DataPlanckLKLType = ...
    BASELINE_18_HIGHL_TTTEEE: DataPlanckLKLType = ...
    BASELINE_18_HIGHL_TTTEEE_LITE: DataPlanckLKLType = ...
    BASELINE_18_HIGHL_TT_LITE: DataPlanckLKLType = ...
    BASELINE_18_LENSING: DataPlanckLKLType = ...
    BASELINE_18_LENSING_CMB_MARGINALIZED: DataPlanckLKLType = ...
    BASELINE_18_LOWL_BB: DataPlanckLKLType = ...
    BASELINE_18_LOWL_EE: DataPlanckLKLType = ...
    BASELINE_18_LOWL_EEBB: DataPlanckLKLType = ...
    BASELINE_18_LOWL_TT: DataPlanckLKLType = ...
    LENGTH: DataPlanckLKLType = ...

class DataSNIACovError(GObject.GEnum):
    ID_NOT_FOUND: DataSNIACovError = ...
    INVALID_ID: DataSNIACovError = ...
    INVALID_SAMPLE: DataSNIACovError = ...
    @staticmethod
    def quark() -> int: ...

class DataSNIACovOrder(GObject.GEnum):
    COLOUR_COLOUR: DataSNIACovOrder = ...
    MAG_COLOUR: DataSNIACovOrder = ...
    MAG_MAG: DataSNIACovOrder = ...
    MAG_WIDTH: DataSNIACovOrder = ...
    WIDTH_COLOUR: DataSNIACovOrder = ...
    WIDTH_WIDTH: DataSNIACovOrder = ...

class DataSNIAId(GObject.GEnum):
    COV_JLA_SNLS3_SDSS_SYS_STAT: DataSNIAId = ...
    COV_JLA_SNLS3_SDSS_SYS_STAT_CMPL: DataSNIAId = ...
    COV_PANTHEON: DataSNIAId = ...
    COV_PANTHEON_PLUS_SH0ES_STAT: DataSNIAId = ...
    COV_PANTHEON_PLUS_SH0ES_SYS_STAT: DataSNIAId = ...
    COV_SNLS3_STAT_ONLY: DataSNIAId = ...
    COV_SNLS3_SYS_STAT: DataSNIAId = ...
    SIMPLE_CFA3: DataSNIAId = ...
    SIMPLE_ESSENCE: DataSNIAId = ...
    SIMPLE_GOLD_157: DataSNIAId = ...
    SIMPLE_GOLD_182: DataSNIAId = ...
    SIMPLE_GOLD_182_FULL: DataSNIAId = ...
    SIMPLE_LEGACY: DataSNIAId = ...
    SIMPLE_SDSS_EMILLE: DataSNIAId = ...
    SIMPLE_UNION: DataSNIAId = ...
    SIMPLE_UNION2: DataSNIAId = ...
    SIMPLE_UNION2_1: DataSNIAId = ...

class DistanceComovingMethod(GObject.GEnum):
    FROM_MODEL: DistanceComovingMethod = ...
    INT_E: DistanceComovingMethod = ...

class GalaxySDShapeGaussParams(GObject.GEnum):
    SIGMA: GalaxySDShapeGaussParams = ...

class GalaxySDTrueRedshiftLSSTSRDSParams(GObject.GEnum):
    ALPHA: GalaxySDTrueRedshiftLSSTSRDSParams = ...
    BETA: GalaxySDTrueRedshiftLSSTSRDSParams = ...
    Z0: GalaxySDTrueRedshiftLSSTSRDSParams = ...

class GalaxyWLObsCoord(GObject.GEnum):
    CELESTIAL: GalaxyWLObsCoord = ...
    EUCLIDEAN: GalaxyWLObsCoord = ...

class GalaxyWLObsEllipConv(GObject.GEnum):
    TRACE: GalaxyWLObsEllipConv = ...
    TRACE_DET: GalaxyWLObsEllipConv = ...

class HICosmoDECplSParams(GObject.GEnum):
    W0: HICosmoDECplSParams = ...
    W1: HICosmoDECplSParams = ...

class HICosmoDEJbpSParams(GObject.GEnum):
    W0: HICosmoDEJbpSParams = ...
    W1: HICosmoDEJbpSParams = ...

class HICosmoDESParams(GObject.GEnum):
    ENNU: HICosmoDESParams = ...
    H0: HICosmoDESParams = ...
    HE_YP: HICosmoDESParams = ...
    OMEGA_B: HICosmoDESParams = ...
    OMEGA_C: HICosmoDESParams = ...
    OMEGA_X: HICosmoDESParams = ...
    T_GAMMA0: HICosmoDESParams = ...

class HICosmoDEVParams(GObject.GEnum):
    G: HICosmoDEVParams = ...
    M: HICosmoDEVParams = ...
    MU: HICosmoDEVParams = ...
    T: HICosmoDEVParams = ...

class HICosmoDEWSplineVParams(GObject.GEnum):
    W: HICosmoDEWSplineVParams = ...

class HICosmoDEXCDMSParams(GObject.GEnum):
    W: HICosmoDEXCDMSParams = ...

class HICosmoGCGSParams(GObject.GEnum):
    ENNU: HICosmoGCGSParams = ...
    GAMMA: HICosmoGCGSParams = ...
    H0: HICosmoGCGSParams = ...
    HE_YP: HICosmoGCGSParams = ...
    OMEGA_B: HICosmoGCGSParams = ...
    OMEGA_C: HICosmoGCGSParams = ...
    OMEGA_X: HICosmoGCGSParams = ...
    T_GAMMA0: HICosmoGCGSParams = ...

class HICosmoGCGVParams(GObject.GEnum):
    G: HICosmoGCGVParams = ...
    M: HICosmoGCGVParams = ...
    MU: HICosmoGCGVParams = ...
    T: HICosmoGCGVParams = ...

class HICosmoIDEM2SParams(GObject.GEnum):
    ENNU: HICosmoIDEM2SParams = ...
    GAMMA: HICosmoIDEM2SParams = ...
    H0: HICosmoIDEM2SParams = ...
    HE_YP: HICosmoIDEM2SParams = ...
    OMEGA_B: HICosmoIDEM2SParams = ...
    OMEGA_C: HICosmoIDEM2SParams = ...
    OMEGA_X: HICosmoIDEM2SParams = ...
    T_GAMMA0: HICosmoIDEM2SParams = ...

class HICosmoIDEM2VParams(GObject.GEnum):
    G: HICosmoIDEM2VParams = ...
    M: HICosmoIDEM2VParams = ...
    MU: HICosmoIDEM2VParams = ...
    T: HICosmoIDEM2VParams = ...

class HICosmoQConstSParams(GObject.GEnum):
    CD: HICosmoQConstSParams = ...
    E: HICosmoQConstSParams = ...
    H0: HICosmoQConstSParams = ...
    OMEGA_T: HICosmoQConstSParams = ...
    Q: HICosmoQConstSParams = ...
    Z1: HICosmoQConstSParams = ...

class HICosmoQGRWSParams(GObject.GEnum):
    H0: HICosmoQGRWSParams = ...
    OMEGA_R: HICosmoQGRWSParams = ...
    OMEGA_W: HICosmoQGRWSParams = ...
    W: HICosmoQGRWSParams = ...
    X_B: HICosmoQGRWSParams = ...

class HICosmoQGWSParams(GObject.GEnum):
    H0: HICosmoQGWSParams = ...
    OMEGA_W: HICosmoQGWSParams = ...
    W: HICosmoQGWSParams = ...
    X_B: HICosmoQGWSParams = ...

class HICosmoQLinearSParams(GObject.GEnum):
    CD: HICosmoQLinearSParams = ...
    E: HICosmoQLinearSParams = ...
    H0: HICosmoQLinearSParams = ...
    OMEGA_T: HICosmoQLinearSParams = ...
    Q: HICosmoQLinearSParams = ...
    QP: HICosmoQLinearSParams = ...
    Z1: HICosmoQLinearSParams = ...

class HICosmoQRBFSParams(GObject.GEnum):
    AS_DRAG: HICosmoQRBFSParams = ...
    H0: HICosmoQRBFSParams = ...
    OMEGA_T: HICosmoQRBFSParams = ...
    RBF_H: HICosmoQRBFSParams = ...

class HICosmoQRBFVParams(GObject.GEnum):
    CENTERS: HICosmoQRBFVParams = ...
    COEFFS: HICosmoQRBFVParams = ...

class HICosmoQSplineSParams(GObject.GEnum):
    AS_DRAG: HICosmoQSplineSParams = ...
    H0: HICosmoQSplineSParams = ...
    OMEGA_T: HICosmoQSplineSParams = ...

class HICosmoQSplineVParams(GObject.GEnum):
    Q: HICosmoQSplineVParams = ...

class HICosmoVexpEMCoupling(GObject.GEnum):
    CAUCHY: HICosmoVexpEMCoupling = ...
    GAUSS: HICosmoVexpEMCoupling = ...
    INVALID: HICosmoVexpEMCoupling = ...
    NONE: HICosmoVexpEMCoupling = ...

class HICosmoVexpSParams(GObject.GEnum):
    ALPHA_B: HICosmoVexpSParams = ...
    D_PHI: HICosmoVexpSParams = ...
    EM_ALPHA: HICosmoVexpSParams = ...
    EM_BETA: HICosmoVexpSParams = ...
    H0: HICosmoVexpSParams = ...
    OMEGA_C: HICosmoVexpSParams = ...
    OMEGA_L: HICosmoVexpSParams = ...
    SIGMA_PHI: HICosmoVexpSParams = ...
    X_B: HICosmoVexpSParams = ...

class HIPertAdiabVars(GObject.GEnum):
    IM_PZETA: HIPertAdiabVars = ...
    IM_ZETA: HIPertAdiabVars = ...
    RE_PZETA: HIPertAdiabVars = ...
    RE_ZETA: HIPertAdiabVars = ...

class HIPertBoltzmannVars(GObject.GEnum):
    B0: HIPertBoltzmannVars = ...
    B1: HIPertBoltzmannVars = ...
    C0: HIPertBoltzmannVars = ...
    C1: HIPertBoltzmannVars = ...
    PHI: HIPertBoltzmannVars = ...
    THETA0: HIPertBoltzmannVars = ...
    THETA1: HIPertBoltzmannVars = ...
    THETA2: HIPertBoltzmannVars = ...
    THETA_P0: HIPertBoltzmannVars = ...
    THETA_P1: HIPertBoltzmannVars = ...
    THETA_P2: HIPertBoltzmannVars = ...

class HIPertCompPBVar(GObject.GEnum):
    DELTA_B: HIPertCompPBVar = ...
    DELTA_G: HIPertCompPBVar = ...
    F_G3: HIPertCompPBVar = ...
    THETA_G: HIPertCompPBVar = ...
    V_B: HIPertCompPBVar = ...
    V_G: HIPertCompPBVar = ...

class HIPertEMVars(GObject.GEnum):
    IM_H: HIPertEMVars = ...
    IM_PH: HIPertEMVars = ...
    RE_H: HIPertEMVars = ...
    RE_PH: HIPertEMVars = ...

class HIPertFirstOrderInteg(GObject.GEnum):
    ARKODE: HIPertFirstOrderInteg = ...
    CVODE: HIPertFirstOrderInteg = ...

class HIPertGWVars(GObject.GEnum):
    IM_H: HIPertGWVars = ...
    IM_PH: HIPertGWVars = ...
    RE_H: HIPertGWVars = ...
    RE_PH: HIPertGWVars = ...

class HIPertGravGauge(GObject.GEnum):
    CONST_CURV: HIPertGravGauge = ...
    CONST_EXP: HIPertGravGauge = ...
    NEWTONIAN: HIPertGravGauge = ...
    SYNCHRONOUS: HIPertGravGauge = ...

class HIPertGravSElem(GObject.GEnum):
    DOTPSI: HIPertGravSElem = ...
    DP: HIPertGravSElem = ...
    DPI: HIPertGravSElem = ...
    DRHO: HIPertGravSElem = ...
    DSIGMA: HIPertGravSElem = ...
    PHI: HIPertGravSElem = ...
    PSI: HIPertGravSElem = ...
    RHOPPV: HIPertGravSElem = ...

class HIPertITwoFluidsVars(GObject.GEnum):
    PS_I: HIPertITwoFluidsVars = ...
    PS_R: HIPertITwoFluidsVars = ...
    PZETA_I: HIPertITwoFluidsVars = ...
    PZETA_R: HIPertITwoFluidsVars = ...
    S_I: HIPertITwoFluidsVars = ...
    S_R: HIPertITwoFluidsVars = ...
    ZETA_I: HIPertITwoFluidsVars = ...
    ZETA_R: HIPertITwoFluidsVars = ...

class HIPertTwoFluidsCross(GObject.GEnum):
    MODE1MAIN: HIPertTwoFluidsCross = ...
    MODE1SUB: HIPertTwoFluidsCross = ...
    MODE2MAIN: HIPertTwoFluidsCross = ...
    MODE2SUB: HIPertTwoFluidsCross = ...

class HIPertWKBCmp(GObject.GEnum):
    ALPHA2: HIPertWKBCmp = ...
    POTENTIAL: HIPertWKBCmp = ...

class HIPertWKBVars(GObject.GEnum):
    IM_P: HIPertWKBVars = ...
    IM_Q: HIPertWKBVars = ...
    RE_P: HIPertWKBVars = ...
    RE_Q: HIPertWKBVars = ...

class HIPrimAtanSParams(GObject.GEnum):
    C2: HIPrimAtanSParams = ...
    C3: HIPrimAtanSParams = ...
    LAMBDA: HIPrimAtanSParams = ...
    LN10E10ASA: HIPrimAtanSParams = ...
    LNKC: HIPrimAtanSParams = ...
    N_SA: HIPrimAtanSParams = ...
    N_T: HIPrimAtanSParams = ...
    T_SA_RATIO: HIPrimAtanSParams = ...

class HIPrimBPLSParams(GObject.GEnum):
    DELTA: HIPrimBPLSParams = ...
    LN10E10ASA: HIPrimBPLSParams = ...
    LNKB: HIPrimBPLSParams = ...
    N_SA: HIPrimBPLSParams = ...
    N_T: HIPrimBPLSParams = ...
    T_SA_RATIO: HIPrimBPLSParams = ...

class HIPrimExpcSParams(GObject.GEnum):
    C: HIPrimExpcSParams = ...
    LAMBDAC: HIPrimExpcSParams = ...
    LN10E10ASA: HIPrimExpcSParams = ...
    LNKC: HIPrimExpcSParams = ...
    N_SA: HIPrimExpcSParams = ...
    N_T: HIPrimExpcSParams = ...
    T_SA_RATIO: HIPrimExpcSParams = ...

class HIPrimPowerLawSParams(GObject.GEnum):
    LN10E10ASA: HIPrimPowerLawSParams = ...
    N_SA: HIPrimPowerLawSParams = ...
    N_T: HIPrimPowerLawSParams = ...
    T_SA_RATIO: HIPrimPowerLawSParams = ...

class HIPrimSBPLSParams(GObject.GEnum):
    DELTA: HIPrimSBPLSParams = ...
    LAMBDA: HIPrimSBPLSParams = ...
    LN10E10ASA: HIPrimSBPLSParams = ...
    LNKB: HIPrimSBPLSParams = ...
    N_SA: HIPrimSBPLSParams = ...
    N_T: HIPrimSBPLSParams = ...
    RA: HIPrimSBPLSParams = ...
    T_SA_RATIO: HIPrimSBPLSParams = ...

class HIPrimTwoFluidsSParams(GObject.GEnum):
    LN10E10ASA: HIPrimTwoFluidsSParams = ...
    LNK0: HIPrimTwoFluidsSParams = ...
    LNW: HIPrimTwoFluidsSParams = ...
    N_T: HIPrimTwoFluidsSParams = ...
    T_SA_RATIO: HIPrimTwoFluidsSParams = ...

class HIReionCambSParams(GObject.GEnum):
    HEIII_Z: HIReionCambSParams = ...
    HII_HEII_Z: HIReionCambSParams = ...

class HaloCMDuffy08SParams(GObject.GEnum):
    DELTA: HaloCMDuffy08SParams = ...

class HaloCMKlypin11SParams(GObject.GEnum):
    DELTA: HaloCMKlypin11SParams = ...

class HaloCMParamSParams(GObject.GEnum):
    C_DELTA: HaloCMParamSParams = ...
    LOG10M_DELTA: HaloCMParamSParams = ...

class HaloDensityProfileDK14MethodParams(GObject.GEnum):
    DIRECT_RHOSRS: HaloDensityProfileDK14MethodParams = ...
    MC2RHOSRS: HaloDensityProfileDK14MethodParams = ...

class HaloDensityProfileDK14Params(GObject.GEnum):
    BE: HaloDensityProfileDK14Params = ...
    BETA: HaloDensityProfileDK14Params = ...
    GAMMA: HaloDensityProfileDK14Params = ...
    RT: HaloDensityProfileDK14Params = ...
    SE: HaloDensityProfileDK14Params = ...

class HaloDensityProfileEinastoParams(GObject.GEnum):
    ALPHA: HaloDensityProfileEinastoParams = ...

class HaloMassFunctionSplineOptimize(GObject.GEnum):
    LNM: HaloMassFunctionSplineOptimize = ...
    NONE: HaloMassFunctionSplineOptimize = ...
    Z: HaloMassFunctionSplineOptimize = ...

class HaloMassSummaryMassDef(GObject.GEnum):
    CRITICAL: HaloMassSummaryMassDef = ...
    MEAN: HaloMassSummaryMassDef = ...
    VIRIAL: HaloMassSummaryMassDef = ...

class HaloPositionSParams(GObject.GEnum):
    DEC: HaloPositionSParams = ...
    RA: HaloPositionSParams = ...
    Z: HaloPositionSParams = ...

class MultiplicityFuncBocquetSim(GObject.GEnum):
    DM: MultiplicityFuncBocquetSim = ...
    HYDRO: MultiplicityFuncBocquetSim = ...

class MultiplicityFuncMassDef(GObject.GEnum):
    CRITICAL: MultiplicityFuncMassDef = ...
    FOF: MultiplicityFuncMassDef = ...
    MEAN: MultiplicityFuncMassDef = ...
    VIRIAL: MultiplicityFuncMassDef = ...

class PlanckFICorTTSParams(GObject.GEnum):
    A_CIB_217: PlanckFICorTTSParams = ...
    A_PLANCK: PlanckFICorTTSParams = ...
    A_SBPX_100_100_TT: PlanckFICorTTSParams = ...
    A_SBPX_143_143_TT: PlanckFICorTTSParams = ...
    A_SBPX_143_217_TT: PlanckFICorTTSParams = ...
    A_SBPX_217_217_TT: PlanckFICorTTSParams = ...
    A_SZ: PlanckFICorTTSParams = ...
    CALIB_100T: PlanckFICorTTSParams = ...
    CALIB_217T: PlanckFICorTTSParams = ...
    CIB_INDEX: PlanckFICorTTSParams = ...
    GAL545_A_100: PlanckFICorTTSParams = ...
    GAL545_A_143: PlanckFICorTTSParams = ...
    GAL545_A_143_217: PlanckFICorTTSParams = ...
    GAL545_A_217: PlanckFICorTTSParams = ...
    KSZ_NORM: PlanckFICorTTSParams = ...
    PS_A_100_100: PlanckFICorTTSParams = ...
    PS_A_143_143: PlanckFICorTTSParams = ...
    PS_A_143_217: PlanckFICorTTSParams = ...
    PS_A_217_217: PlanckFICorTTSParams = ...
    XI_SZ_CIB: PlanckFICorTTSParams = ...

class PlanckFICorTTTEEESParams(GObject.GEnum):
    A_CNOISE_E2E_100_100_EE: PlanckFICorTTTEEESParams = ...
    A_CNOISE_E2E_143_143_EE: PlanckFICorTTTEEESParams = ...
    A_CNOISE_E2E_217_217_EE: PlanckFICorTTTEEESParams = ...
    A_POL: PlanckFICorTTTEEESParams = ...
    A_SBPX_100_100_EE: PlanckFICorTTTEEESParams = ...
    A_SBPX_100_143_EE: PlanckFICorTTTEEESParams = ...
    A_SBPX_100_217_EE: PlanckFICorTTTEEESParams = ...
    A_SBPX_143_143_EE: PlanckFICorTTTEEESParams = ...
    A_SBPX_143_217_EE: PlanckFICorTTTEEESParams = ...
    A_SBPX_217_217_EE: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_0_0E_0E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_0_0E_1E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_0_0E_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_0_0T_0E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_0_0T_1E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_0_0T_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_0_1E_1E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_0_1E_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_0_1T_1E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_0_1T_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_0_2E_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_0_2T_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_1_0E_0E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_1_0E_1E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_1_0E_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_1_0T_0E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_1_0T_1E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_1_0T_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_1_1E_1E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_1_1E_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_1_1T_1E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_1_1T_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_1_2E_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_1_2T_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_2_0E_0E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_2_0E_1E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_2_0E_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_2_0T_0E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_2_0T_1E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_2_0T_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_2_1E_1E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_2_1E_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_2_1T_1E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_2_1T_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_2_2E_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_2_2T_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_3_0E_0E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_3_0E_1E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_3_0E_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_3_0T_0E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_3_0T_1E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_3_0T_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_3_1E_1E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_3_1E_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_3_1T_1E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_3_1T_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_3_2E_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_3_2T_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_4_0E_0E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_4_0E_1E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_4_0E_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_4_0T_0E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_4_0T_1E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_4_0T_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_4_1E_1E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_4_1E_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_4_1T_1E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_4_1T_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_4_2E_2E: PlanckFICorTTTEEESParams = ...
    BLEAK_EPSILON_4_2T_2E: PlanckFICorTTTEEESParams = ...
    CALIB_100P: PlanckFICorTTTEEESParams = ...
    CALIB_143P: PlanckFICorTTTEEESParams = ...
    CALIB_217P: PlanckFICorTTTEEESParams = ...
    GALF_EE_A_100: PlanckFICorTTTEEESParams = ...
    GALF_EE_A_100_143: PlanckFICorTTTEEESParams = ...
    GALF_EE_A_100_217: PlanckFICorTTTEEESParams = ...
    GALF_EE_A_143: PlanckFICorTTTEEESParams = ...
    GALF_EE_A_143_217: PlanckFICorTTTEEESParams = ...
    GALF_EE_A_217: PlanckFICorTTTEEESParams = ...
    GALF_EE_INDEX: PlanckFICorTTTEEESParams = ...
    GALF_TE_A_100: PlanckFICorTTTEEESParams = ...
    GALF_TE_A_100_143: PlanckFICorTTTEEESParams = ...
    GALF_TE_A_100_217: PlanckFICorTTTEEESParams = ...
    GALF_TE_A_143: PlanckFICorTTTEEESParams = ...
    GALF_TE_A_143_217: PlanckFICorTTTEEESParams = ...
    GALF_TE_A_217: PlanckFICorTTTEEESParams = ...
    GALF_TE_INDEX: PlanckFICorTTTEEESParams = ...

class ReducedShearCalibWtgSParams(GObject.GEnum):
    C: ReducedShearCalibWtgSParams = ...
    MB: ReducedShearCalibWtgSParams = ...
    MSLOPE: ReducedShearCalibWtgSParams = ...
    SIZE_RATIO: ReducedShearCalibWtgSParams = ...

class ReducedShearClusterMassParams(GObject.GEnum):
    A: ReducedShearClusterMassParams = ...
    B: ReducedShearClusterMassParams = ...
    C: ReducedShearClusterMassParams = ...
    VGAMMA: ReducedShearClusterMassParams = ...
    VSIGMA: ReducedShearClusterMassParams = ...
    XP: ReducedShearClusterMassParams = ...

class SNIADistCovSParams(GObject.GEnum):
    ALPHA: SNIADistCovSParams = ...
    BETA: SNIADistCovSParams = ...
    LNSIGMA_LENS: SNIADistCovSParams = ...
    LNSIGMA_PECZ: SNIADistCovSParams = ...
    M1: SNIADistCovSParams = ...
    M2: SNIADistCovSParams = ...

class SNIADistCovVParams(GObject.GEnum):
    LNSIGMA_INT: SNIADistCovVParams = ...
    MU: SNIADistCovVParams = ...

class TransferFuncBBKSType(GObject.GEnum):
    BARYONS: TransferFuncBBKSType = ...
    CCL: TransferFuncBBKSType = ...
    NOBARYONS: TransferFuncBBKSType = ...

class WLSurfaceMassDensityParams(GObject.GEnum):
    PCC: WLSurfaceMassDensityParams = ...
    ROFF: WLSurfaceMassDensityParams = ...

class XcorLimberKernelCMBLensingSParams(GObject.GEnum):
    LEN: XcorLimberKernelCMBLensingSParams = ...

class XcorLimberKernelGalSParams(GObject.GEnum):
    MAG_BIAS: XcorLimberKernelGalSParams = ...
    NOISE_BIAS: XcorLimberKernelGalSParams = ...

class XcorLimberKernelGalVParams(GObject.GEnum):
    BIAS: XcorLimberKernelGalVParams = ...

class XcorLimberKernelImpl(GObject.GEnum):
    ADD_NOISE: XcorLimberKernelImpl = ...
    EVAL: XcorLimberKernelImpl = ...
    PREPARE: XcorLimberKernelImpl = ...

class XcorLimberKernelWeakLensingSParams(GObject.GEnum):
    LEN: XcorLimberKernelWeakLensingSParams = ...

class XcorLimberKernelWeakLensingVParams(GObject.GEnum):
    LEN: XcorLimberKernelWeakLensingVParams = ...

class XcorLimberKerneltSZSParams(GObject.GEnum):
    LEN: XcorLimberKerneltSZSParams = ...

class XcorLimberMethod(GObject.GEnum):
    CUBATURE: XcorLimberMethod = ...
    GSL: XcorLimberMethod = ...
