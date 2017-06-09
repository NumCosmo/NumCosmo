
/* Generated data (by glib-mkenums) */


#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

/**
 * SECTION:ncm_fit_nlopt_enum
 * @title: NcmFitNLOptEnum
 * @short_description: Automaticaly imported enum from NLOpt library.
 *
 */
/* enumerations from "ncm_fit_nlopt_enum_meta.h" */
GType
ncm_fit_nlopt_algorithm_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NLOPT_GN_DIRECT, "NLOPT_GN_DIRECT", "gn-direct" },
        { NLOPT_GN_DIRECT_L, "NLOPT_GN_DIRECT_L", "gn-direct-l" },
        { NLOPT_GN_DIRECT_L_RAND, "NLOPT_GN_DIRECT_L_RAND", "gn-direct-l-rand" },
        { NLOPT_GN_DIRECT_NOSCAL, "NLOPT_GN_DIRECT_NOSCAL", "gn-direct-noscal" },
        { NLOPT_GN_DIRECT_L_NOSCAL, "NLOPT_GN_DIRECT_L_NOSCAL", "gn-direct-l-noscal" },
        { NLOPT_GN_DIRECT_L_RAND_NOSCAL, "NLOPT_GN_DIRECT_L_RAND_NOSCAL", "gn-direct-l-rand-noscal" },
        { NLOPT_GN_ORIG_DIRECT, "NLOPT_GN_ORIG_DIRECT", "gn-orig-direct" },
        { NLOPT_GN_ORIG_DIRECT_L, "NLOPT_GN_ORIG_DIRECT_L", "gn-orig-direct-l" },
        { NLOPT_GD_STOGO, "NLOPT_GD_STOGO", "gd-stogo" },
        { NLOPT_GD_STOGO_RAND, "NLOPT_GD_STOGO_RAND", "gd-stogo-rand" },
        { NLOPT_LD_LBFGS_NOCEDAL, "NLOPT_LD_LBFGS_NOCEDAL", "ld-lbfgs-nocedal" },
        { NLOPT_LD_LBFGS, "NLOPT_LD_LBFGS", "ld-lbfgs" },
        { NLOPT_LN_PRAXIS, "NLOPT_LN_PRAXIS", "ln-praxis" },
        { NLOPT_LD_VAR1, "NLOPT_LD_VAR1", "ld-var1" },
        { NLOPT_LD_VAR2, "NLOPT_LD_VAR2", "ld-var2" },
        { NLOPT_LD_TNEWTON, "NLOPT_LD_TNEWTON", "ld-tnewton" },
        { NLOPT_LD_TNEWTON_RESTART, "NLOPT_LD_TNEWTON_RESTART", "ld-tnewton-restart" },
        { NLOPT_LD_TNEWTON_PRECOND, "NLOPT_LD_TNEWTON_PRECOND", "ld-tnewton-precond" },
        { NLOPT_LD_TNEWTON_PRECOND_RESTART, "NLOPT_LD_TNEWTON_PRECOND_RESTART", "ld-tnewton-precond-restart" },
        { NLOPT_GN_CRS2_LM, "NLOPT_GN_CRS2_LM", "gn-crs2-lm" },
        { NLOPT_GN_MLSL, "NLOPT_GN_MLSL", "gn-mlsl" },
        { NLOPT_GD_MLSL, "NLOPT_GD_MLSL", "gd-mlsl" },
        { NLOPT_GN_MLSL_LDS, "NLOPT_GN_MLSL_LDS", "gn-mlsl-lds" },
        { NLOPT_GD_MLSL_LDS, "NLOPT_GD_MLSL_LDS", "gd-mlsl-lds" },
        { NLOPT_LD_MMA, "NLOPT_LD_MMA", "ld-mma" },
        { NLOPT_LN_COBYLA, "NLOPT_LN_COBYLA", "ln-cobyla" },
        { NLOPT_LN_NEWUOA, "NLOPT_LN_NEWUOA", "ln-newuoa" },
        { NLOPT_LN_NEWUOA_BOUND, "NLOPT_LN_NEWUOA_BOUND", "ln-newuoa-bound" },
        { NLOPT_LN_NELDERMEAD, "NLOPT_LN_NELDERMEAD", "ln-neldermead" },
        { NLOPT_LN_SBPLX, "NLOPT_LN_SBPLX", "ln-sbplx" },
        { NLOPT_LN_AUGLAG, "NLOPT_LN_AUGLAG", "ln-auglag" },
        { NLOPT_LD_AUGLAG, "NLOPT_LD_AUGLAG", "ld-auglag" },
        { NLOPT_LN_AUGLAG_EQ, "NLOPT_LN_AUGLAG_EQ", "ln-auglag-eq" },
        { NLOPT_LD_AUGLAG_EQ, "NLOPT_LD_AUGLAG_EQ", "ld-auglag-eq" },
        { NLOPT_LN_BOBYQA, "NLOPT_LN_BOBYQA", "ln-bobyqa" },
        { NLOPT_GN_ISRES, "NLOPT_GN_ISRES", "gn-isres" },
        { NLOPT_AUGLAG, "NLOPT_AUGLAG", "auglag" },
        { NLOPT_AUGLAG_EQ, "NLOPT_AUGLAG_EQ", "auglag-eq" },
        { NLOPT_G_MLSL, "NLOPT_G_MLSL", "g-mlsl" },
        { NLOPT_G_MLSL_LDS, "NLOPT_G_MLSL_LDS", "g-mlsl-lds" },
        { NLOPT_LD_SLSQP, "NLOPT_LD_SLSQP", "ld-slsqp" },
        { NLOPT_LD_CCSAQ, "NLOPT_LD_CCSAQ", "ld-ccsaq" },
        { NLOPT_GN_ESCH, "NLOPT_GN_ESCH", "gn-esch" },
        { NLOPT_NUM_ALGORITHMS, "NLOPT_NUM_ALGORITHMS", "num-algorithms" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmFitNloptAlgorithm"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 
GType
ncm_fit_nlopt_result_get_type (void)
{
  static volatile gsize g_define_type_id__volatile = 0;

  if (g_once_init_enter (&g_define_type_id__volatile))
    {
      static const GEnumValue values[] = {
        { NLOPT_FAILURE, "NLOPT_FAILURE", "failure" },
        { NLOPT_INVALID_ARGS, "NLOPT_INVALID_ARGS", "invalid-args" },
        { NLOPT_OUT_OF_MEMORY, "NLOPT_OUT_OF_MEMORY", "out-of-memory" },
        { NLOPT_ROUNDOFF_LIMITED, "NLOPT_ROUNDOFF_LIMITED", "roundoff-limited" },
        { NLOPT_FORCED_STOP, "NLOPT_FORCED_STOP", "forced-stop" },
        { NLOPT_SUCCESS, "NLOPT_SUCCESS", "success" },
        { NLOPT_STOPVAL_REACHED, "NLOPT_STOPVAL_REACHED", "stopval-reached" },
        { NLOPT_FTOL_REACHED, "NLOPT_FTOL_REACHED", "ftol-reached" },
        { NLOPT_XTOL_REACHED, "NLOPT_XTOL_REACHED", "xtol-reached" },
        { NLOPT_MAXEVAL_REACHED, "NLOPT_MAXEVAL_REACHED", "maxeval-reached" },
        { NLOPT_MAXTIME_REACHED, "NLOPT_MAXTIME_REACHED", "maxtime-reached" },
        { 0, NULL, NULL }
      };
      GType g_define_type_id =
        g_enum_register_static (g_intern_static_string ("NcmFitNloptResult"), values);
      g_once_init_leave (&g_define_type_id__volatile, g_define_type_id);
    }

  return g_define_type_id__volatile;
}
 

/* Generated data ends here */

