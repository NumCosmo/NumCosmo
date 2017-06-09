
/* Generated data (by glib-mkenums) */


#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

/**
 * SECTION:ncm_enum_types     
 * @title: NcmEnumTypes                           
 * @short_description: Automaticaly generated enum types from NumCosmoMath library.
 *   
 */
/* enumerations from "math/function_cache.h" */
GType
ncm_function_cache_search_type_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NC_FUNCTION_CACHE_SEARCH_BOTH, "NC_FUNCTION_CACHE_SEARCH_BOTH", "both" },
        { NC_FUNCTION_CACHE_SEARCH_GT, "NC_FUNCTION_CACHE_SEARCH_GT", "gt" },
        { NC_FUNCTION_CACHE_SEARCH_LT, "NC_FUNCTION_CACHE_SEARCH_LT", "lt" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmFunctionCacheSearchType"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
/* enumerations from "math/grid_one.h" */
GType
ncm_grid_nodes_end_points_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_GRID_NODES_START, "NCM_GRID_NODES_START", "start" },
        { NCM_GRID_NODES_END, "NCM_GRID_NODES_END", "end" },
        { NCM_GRID_NODES_BOTH, "NCM_GRID_NODES_BOTH", "both" },
        { NCM_GRID_NODES_NONE, "NCM_GRID_NODES_NONE", "none" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmGridNodesEndPoints"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
/* enumerations from "math/ncm_data_poisson.h" */
GType
ncm_data_poisson_type_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_DATA_POISSON_INT, "NCM_DATA_POISSON_INT", "int" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmDataPoissonType"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
/* enumerations from "math/ncm_dataset.h" */
GType
ncm_dataset_bstrap_type_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_DATASET_BSTRAP_DISABLE, "NCM_DATASET_BSTRAP_DISABLE", "disable" },
        { NCM_DATASET_BSTRAP_PARTIAL, "NCM_DATASET_BSTRAP_PARTIAL", "partial" },
        { NCM_DATASET_BSTRAP_TOTAL, "NCM_DATASET_BSTRAP_TOTAL", "total" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmDatasetBStrapType"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
/* enumerations from "math/ncm_fit.h" */
GType
ncm_fit_type_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_FIT_TYPE_GSL_LS, "NCM_FIT_TYPE_GSL_LS", "gsl-ls" },
        { NCM_FIT_TYPE_GSL_MM, "NCM_FIT_TYPE_GSL_MM", "gsl-mm" },
        { NCM_FIT_TYPE_GSL_MMS, "NCM_FIT_TYPE_GSL_MMS", "gsl-mms" },
        { NCM_FIT_TYPE_LEVMAR, "NCM_FIT_TYPE_LEVMAR", "levmar" },
        { NCM_FIT_TYPE_NLOPT, "NCM_FIT_TYPE_NLOPT", "nlopt" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmFitType"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
GType
ncm_fit_grad_type_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_FIT_GRAD_ANALYTICAL, "NCM_FIT_GRAD_ANALYTICAL", "analytical" },
        { NCM_FIT_GRAD_NUMDIFF_FORWARD, "NCM_FIT_GRAD_NUMDIFF_FORWARD", "numdiff-forward" },
        { NCM_FIT_GRAD_NUMDIFF_CENTRAL, "NCM_FIT_GRAD_NUMDIFF_CENTRAL", "numdiff-central" },
        { NCM_FIT_GRAD_NUMDIFF_ACCURATE, "NCM_FIT_GRAD_NUMDIFF_ACCURATE", "numdiff-accurate" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmFitGradType"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
GType
ncm_fit_run_msgs_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_FIT_RUN_MSGS_NONE, "NCM_FIT_RUN_MSGS_NONE", "none" },
        { NCM_FIT_RUN_MSGS_SIMPLE, "NCM_FIT_RUN_MSGS_SIMPLE", "simple" },
        { NCM_FIT_RUN_MSGS_FULL, "NCM_FIT_RUN_MSGS_FULL", "full" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmFitRunMsgs"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
/* enumerations from "math/ncm_fit_gsl_mm.h" */
GType
ncm_fit_gslmm_algos_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_FIT_GSL_MM_CONJUGATE_FR, "NCM_FIT_GSL_MM_CONJUGATE_FR", "conjugate-fr" },
        { NCM_FIT_GSL_MM_CONJUGATE_PR, "NCM_FIT_GSL_MM_CONJUGATE_PR", "conjugate-pr" },
        { NCM_FIT_GSL_MM_VECTOR_BFGS, "NCM_FIT_GSL_MM_VECTOR_BFGS", "vector-bfgs" },
        { NCM_FIT_GSL_MM_VECTOR_BFGS2, "NCM_FIT_GSL_MM_VECTOR_BFGS2", "vector-bfgs2" },
        { NCM_FIT_GSL_MM_STEEPEST_DESCENT, "NCM_FIT_GSL_MM_STEEPEST_DESCENT", "steepest-descent" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmFitGSLMMAlgos"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
/* enumerations from "math/ncm_fit_gsl_mms.h" */
GType
ncm_fit_gslmms_algos_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_FIT_GSL_MMS_NMSIMPLEX2, "NCM_FIT_GSL_MMS_NMSIMPLEX2", "nmsimplex2" },
        { NCM_FIT_GSL_MMS_NMSIMPLEX, "NCM_FIT_GSL_MMS_NMSIMPLEX", "nmsimplex" },
        { NCM_FIT_GSL_MMS_NMSIMPLEX2RAND, "NCM_FIT_GSL_MMS_NMSIMPLEX2RAND", "nmsimplex2rand" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmFitGSLMMSAlgos"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
/* enumerations from "math/ncm_fit_levmar.h" */
GType
ncm_fit_levmar_algos_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_FIT_LEVMAR_DER, "NCM_FIT_LEVMAR_DER", "der" },
        { NCM_FIT_LEVMAR_DIF, "NCM_FIT_LEVMAR_DIF", "dif" },
        { NCM_FIT_LEVMAR_BC_DER, "NCM_FIT_LEVMAR_BC_DER", "bc-der" },
        { NCM_FIT_LEVMAR_BC_DIF, "NCM_FIT_LEVMAR_BC_DIF", "bc-dif" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmFitLevmarAlgos"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
/* enumerations from "math/ncm_fit_mc.h" */
GType
ncm_fit_mc_resample_type_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_FIT_MC_RESAMPLE_FROM_MODEL, "NCM_FIT_MC_RESAMPLE_FROM_MODEL", "from-model" },
        { NCM_FIT_MC_RESAMPLE_BOOTSTRAP_NOMIX, "NCM_FIT_MC_RESAMPLE_BOOTSTRAP_NOMIX", "bootstrap-nomix" },
        { NCM_FIT_MC_RESAMPLE_BOOTSTRAP_MIX, "NCM_FIT_MC_RESAMPLE_BOOTSTRAP_MIX", "bootstrap-mix" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmFitMCResampleType"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
/* enumerations from "math/ncm_hoaa.h" */
GType
ncm_hoaa_opt_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_HOAA_OPT_FULL, "NCM_HOAA_OPT_FULL", "full" },
        { NCM_HOAA_OPT_V_ONLY, "NCM_HOAA_OPT_V_ONLY", "v-only" },
        { NCM_HOAA_OPT_DLNMNU_ONLY, "NCM_HOAA_OPT_DLNMNU_ONLY", "dlnmnu-only" },
        { NCM_HOAA_OPT_INVALID, "NCM_HOAA_OPT_INVALID", "invalid" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmHOAAOpt"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
GType
ncm_hoaa_sing_type_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_HOAA_SING_TYPE_ZERO, "NCM_HOAA_SING_TYPE_ZERO", "zero" },
        { NCM_HOAA_SING_TYPE_INF, "NCM_HOAA_SING_TYPE_INF", "inf" },
        { NCM_HOAA_SING_TYPE_INVALID, "NCM_HOAA_SING_TYPE_INVALID", "invalid" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmHOAASingType"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
GType
ncm_hoaa_var_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_HOAA_VAR_TS, "NCM_HOAA_VAR_TS", "ts" },
        { NCM_HOAA_VAR_TC, "NCM_HOAA_VAR_TC", "tc" },
        { NCM_HOAA_VAR_EPSILON, "NCM_HOAA_VAR_EPSILON", "epsilon" },
        { NCM_HOAA_VAR_GAMMA, "NCM_HOAA_VAR_GAMMA", "gamma" },
        { NCM_HOAA_VAR_SYS_SIZE, "NCM_HOAA_VAR_SYS_SIZE", "sys-size" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmHOAAVar"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
/* enumerations from "math/ncm_lh_ratio1d.h" */
GType
ncm_lh_ratio1d_root_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_LH_RATIO1D_ROOT_BRACKET, "NCM_LH_RATIO1D_ROOT_BRACKET", "bracket" },
        { NCM_LH_RATIO1D_ROOT_NUMDIFF, "NCM_LH_RATIO1D_ROOT_NUMDIFF", "numdiff" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmLHRatio1dRoot"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
/* enumerations from "math/ncm_lh_ratio2d.h" */
GType
ncm_lh_ratio2d_root_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_LH_RATIO2D_ROOT_BRACKET, "NCM_LH_RATIO2D_ROOT_BRACKET", "bracket" },
        { NCM_LH_RATIO2D_ROOT_NUMDIFF, "NCM_LH_RATIO2D_ROOT_NUMDIFF", "numdiff" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmLHRatio2dRoot"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
/* enumerations from "math/ncm_matrix.h" */
GType
ncm_matrix_internal_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_MATRIX_SLICE, "NCM_MATRIX_SLICE", "slice" },
        { NCM_MATRIX_GSL_MATRIX, "NCM_MATRIX_GSL_MATRIX", "gsl-matrix" },
        { NCM_MATRIX_MALLOC, "NCM_MATRIX_MALLOC", "malloc" },
        { NCM_MATRIX_GARRAY, "NCM_MATRIX_GARRAY", "garray" },
        { NCM_MATRIX_DERIVED, "NCM_MATRIX_DERIVED", "derived" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmMatrixInternal"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
/* enumerations from "math/ncm_mset_catalog.h" */
GType
ncm_mset_catalog_sync_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_MSET_CATALOG_SYNC_DISABLE, "NCM_MSET_CATALOG_SYNC_DISABLE", "disable" },
        { NCM_MSET_CATALOG_SYNC_AUTO, "NCM_MSET_CATALOG_SYNC_AUTO", "auto" },
        { NCM_MSET_CATALOG_SYNC_TIMED, "NCM_MSET_CATALOG_SYNC_TIMED", "timed" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmMSetCatalogSync"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
GType
ncm_mset_catalog_trim_type_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GFlagsValue values[] = {
        { NCM_MSET_CATALOG_TRIM_TYPE_ESS, "NCM_MSET_CATALOG_TRIM_TYPE_ESS", "ess" },
        { NCM_MSET_CATALOG_TRIM_TYPE_HEIDEL, "NCM_MSET_CATALOG_TRIM_TYPE_HEIDEL", "heidel" },
        { NCM_MSET_CATALOG_TRIM_TYPE_ALL, "NCM_MSET_CATALOG_TRIM_TYPE_ALL", "all" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_flags_register_static (g_intern_static_string ("NcmMSetCatalogTrimType"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
GType
ncm_mset_catalog_tau_method_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_MSET_CATALOG_TAU_METHOD_ACOR, "NCM_MSET_CATALOG_TAU_METHOD_ACOR", "acor" },
        { NCM_MSET_CATALOG_TAU_METHOD_AR_MODEL, "NCM_MSET_CATALOG_TAU_METHOD_AR_MODEL", "ar-model" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmMSetCatalogTauMethod"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
/* enumerations from "math/ncm_powspec_filter.h" */
GType
ncm_powspec_filter_type_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_POWSPEC_FILTER_TYPE_TOPHAT, "NCM_POWSPEC_FILTER_TYPE_TOPHAT", "tophat" },
        { NCM_POWSPEC_FILTER_TYPE_GAUSS, "NCM_POWSPEC_FILTER_TYPE_GAUSS", "gauss" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmPowspecFilterType"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
/* enumerations from "math/ncm_serialize.h" */
GType
ncm_serialize_opt_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GFlagsValue values[] = {
        { NCM_SERIALIZE_OPT_NONE, "NCM_SERIALIZE_OPT_NONE", "none" },
        { NCM_SERIALIZE_OPT_AUTOSAVE_SER, "NCM_SERIALIZE_OPT_AUTOSAVE_SER", "autosave-ser" },
        { NCM_SERIALIZE_OPT_AUTONAME_SER, "NCM_SERIALIZE_OPT_AUTONAME_SER", "autoname-ser" },
        { NCM_SERIALIZE_OPT_CLEAN_DUP, "NCM_SERIALIZE_OPT_CLEAN_DUP", "clean-dup" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_flags_register_static (g_intern_static_string ("NcmSerializeOpt"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
/* enumerations from "math/ncm_sparam.h" */
GType
ncm_param_type_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_PARAM_TYPE_FREE, "NCM_PARAM_TYPE_FREE", "free" },
        { NCM_PARAM_TYPE_FIXED, "NCM_PARAM_TYPE_FIXED", "fixed" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmParamType"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
/* enumerations from "math/ncm_sphere_map.h" */
GType
ncm_sphere_map_order_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NC_SPHERE_MAP_ORDER_NEST, "NC_SPHERE_MAP_ORDER_NEST", "nest" },
        { NC_SPHERE_MAP_ORDER_RING, "NC_SPHERE_MAP_ORDER_RING", "ring" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmSphereMapOrder"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
GType
ncm_sphere_map_type_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GFlagsValue values[] = {
        { NC_SPHERE_MAP_TYPE_TEMPERATURE, "NC_SPHERE_MAP_TYPE_TEMPERATURE", "temperature" },
        { NC_SPHERE_MAP_TYPE_Q_POLARIZATION, "NC_SPHERE_MAP_TYPE_Q_POLARIZATION", "q-polarization" },
        { NC_SPHERE_MAP_TYPE_U_POLARISATION, "NC_SPHERE_MAP_TYPE_U_POLARISATION", "u-polarisation" },
        { NC_SPHERE_MAP_TYPE_SPUR_SIGNAL, "NC_SPHERE_MAP_TYPE_SPUR_SIGNAL", "spur-signal" },
        { NC_SPHERE_MAP_TYPE_N_OBS, "NC_SPHERE_MAP_TYPE_N_OBS", "n-obs" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_flags_register_static (g_intern_static_string ("NcmSphereMapType"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
/* enumerations from "math/ncm_sphere_map_pix.h" */
GType
ncm_sphere_map_pix_order_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_SPHERE_MAP_PIX_ORDER_NEST, "NCM_SPHERE_MAP_PIX_ORDER_NEST", "nest" },
        { NCM_SPHERE_MAP_PIX_ORDER_RING, "NCM_SPHERE_MAP_PIX_ORDER_RING", "ring" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmSphereMapPixOrder"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
GType
ncm_sphere_map_pix_coord_sys_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_SPHERE_MAP_PIX_COORD_SYS_GALACTIC, "NCM_SPHERE_MAP_PIX_COORD_SYS_GALACTIC", "galactic" },
        { NCM_SPHERE_MAP_PIX_COORD_SYS_ECLIPTIC, "NCM_SPHERE_MAP_PIX_COORD_SYS_ECLIPTIC", "ecliptic" },
        { NCM_SPHERE_MAP_PIX_COORD_SYS_CELESTIAL, "NCM_SPHERE_MAP_PIX_COORD_SYS_CELESTIAL", "celestial" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmSphereMapPixCoordSys"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
/* enumerations from "math/ncm_spline_func.h" */
GType
ncm_spline_func_type_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_SPLINE_FUNCTION_4POINTS, "NCM_SPLINE_FUNCTION_4POINTS", "4points" },
        { NCM_SPLINE_FUNCTION_2x2POINTS, "NCM_SPLINE_FUNCTION_2x2POINTS", "2x2points" },
        { NCM_SPLINE_FUNCTION_SPLINE, "NCM_SPLINE_FUNCTION_SPLINE", "spline" },
        { NCM_SPLINE_FUNCTION_SPLINE_LNKNOT, "NCM_SPLINE_FUNCTION_SPLINE_LNKNOT", "spline-lnknot" },
        { NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT, "NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT", "spline-sinhknot" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmSplineFuncType"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
/* enumerations from "math/ncm_spline_gsl.h" */
GType
ncm_spline_gsl_type_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_SPLINE_GSL_LINEAR, "NCM_SPLINE_GSL_LINEAR", "linear" },
        { NCM_SPLINE_GSL_POLYNOMIAL, "NCM_SPLINE_GSL_POLYNOMIAL", "polynomial" },
        { NCM_SPLINE_GSL_CSPLINE, "NCM_SPLINE_GSL_CSPLINE", "cspline" },
        { NCM_SPLINE_GSL_CSPLINE_PERIODIC, "NCM_SPLINE_GSL_CSPLINE_PERIODIC", "cspline-periodic" },
        { NCM_SPLINE_GSL_AKIMA, "NCM_SPLINE_GSL_AKIMA", "akima" },
        { NCM_SPLINE_GSL_AKIMA_PERIODIC, "NCM_SPLINE_GSL_AKIMA_PERIODIC", "akima-periodic" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmSplineGslType"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
/* enumerations from "math/ncm_stats_dist1d_epdf.h" */
GType
ncm_stats_dist1d_epdf_bw_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_STATS_DIST1D_EPDF_BW_FIXED, "NCM_STATS_DIST1D_EPDF_BW_FIXED", "fixed" },
        { NCM_STATS_DIST1D_EPDF_BW_RoT, "NCM_STATS_DIST1D_EPDF_BW_RoT", "rot" },
        { NCM_STATS_DIST1D_EPDF_BW_AUTO, "NCM_STATS_DIST1D_EPDF_BW_AUTO", "auto" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmStatsDist1dEPDFBw"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
/* enumerations from "math/ncm_stats_vec.h" */
GType
ncm_stats_vec_type_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_STATS_VEC_MEAN, "NCM_STATS_VEC_MEAN", "mean" },
        { NCM_STATS_VEC_VAR, "NCM_STATS_VEC_VAR", "var" },
        { NCM_STATS_VEC_COV, "NCM_STATS_VEC_COV", "cov" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmStatsVecType"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
GType
ncm_stats_vec_ar_type_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_STATS_VEC_AR_NONE, "NCM_STATS_VEC_AR_NONE", "none" },
        { NCM_STATS_VEC_AR_FPE, "NCM_STATS_VEC_AR_FPE", "fpe" },
        { NCM_STATS_VEC_AR_AIC, "NCM_STATS_VEC_AR_AIC", "aic" },
        { NCM_STATS_VEC_AR_AICC, "NCM_STATS_VEC_AR_AICC", "aicc" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmStatsVecARType"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
/* enumerations from "math/ncm_vector.h" */
GType
ncm_vector_internal_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NCM_VECTOR_SLICE, "NCM_VECTOR_SLICE", "slice" },
        { NCM_VECTOR_GSL_VECTOR, "NCM_VECTOR_GSL_VECTOR", "gsl-vector" },
        { NCM_VECTOR_MALLOC, "NCM_VECTOR_MALLOC", "malloc" },
        { NCM_VECTOR_ARRAY, "NCM_VECTOR_ARRAY", "array" },
        { NCM_VECTOR_DERIVED, "NCM_VECTOR_DERIVED", "derived" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmVectorInternal"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 

/* Generated data ends here */

