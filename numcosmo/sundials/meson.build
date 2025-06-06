sundials_arkode_sources = files(
    'include/arkode/arkode.h',
    'include/arkode/arkode_arkstep.h',
    'include/arkode/arkode_bandpre.h',
    'include/arkode/arkode_bbdpre.h',
    'include/arkode/arkode_butcher.h',
    'include/arkode/arkode_butcher_dirk.h',
    'include/arkode/arkode_butcher_erk.h',
    'include/arkode/arkode_erkstep.h',
    'include/arkode/arkode_forcingstep.h',
    'include/arkode/arkode_ls.h',
    'include/arkode/arkode_lsrkstep.h',
    'include/arkode/arkode_mristep.h',
    'include/arkode/arkode_splittingstep.h',
    'include/arkode/arkode_sprk.h',
    'include/arkode/arkode_sprkstep.h',
    'src/arkode/arkode.c',
    'src/arkode/arkode_adapt.c',
    'src/arkode/arkode_adapt_impl.h',
    'src/arkode/arkode_arkstep.c',
    'src/arkode/arkode_arkstep_impl.h',
    'src/arkode/arkode_arkstep_io.c',
    'src/arkode/arkode_arkstep_nls.c',
    'src/arkode/arkode_bandpre.c',
    'src/arkode/arkode_bandpre_impl.h',
    'src/arkode/arkode_bbdpre.c',
    'src/arkode/arkode_bbdpre_impl.h',
    'src/arkode/arkode_butcher.c',
    'src/arkode/arkode_butcher_dirk.c',
    'src/arkode/arkode_butcher_erk.c',
    'src/arkode/arkode_erkstep.c',
    'src/arkode/arkode_erkstep_impl.h',
    'src/arkode/arkode_erkstep_io.c',
    'src/arkode/arkode_forcingstep.c',
    'src/arkode/arkode_forcingstep_impl.h',
    'src/arkode/arkode_impl.h',
    'src/arkode/arkode_interp.c',
    'src/arkode/arkode_interp_impl.h',
    'src/arkode/arkode_io.c',
    'src/arkode/arkode_ls.c',
    'src/arkode/arkode_ls_impl.h',
    'src/arkode/arkode_lsrkstep.c',
    'src/arkode/arkode_lsrkstep_impl.h',
    'src/arkode/arkode_lsrkstep_io.c',
    'src/arkode/arkode_mri_tables.c',
    'src/arkode/arkode_mri_tables_impl.h',
    'src/arkode/arkode_mristep.c',
    'src/arkode/arkode_mristep_controller.c',
    'src/arkode/arkode_mristep_impl.h',
    'src/arkode/arkode_mristep_io.c',
    'src/arkode/arkode_mristep_nls.c',
    'src/arkode/arkode_relaxation.c',
    'src/arkode/arkode_relaxation_impl.h',
    'src/arkode/arkode_root.c',
    'src/arkode/arkode_root_impl.h',
    'src/arkode/arkode_splittingstep.c',
    'src/arkode/arkode_splittingstep_coefficients.c',
    'src/arkode/arkode_splittingstep_impl.h',
    'src/arkode/arkode_sprk.c',
    'src/arkode/arkode_sprkstep.c',
    'src/arkode/arkode_sprkstep_impl.h',
    'src/arkode/arkode_sprkstep_io.c',
    'src/arkode/arkode_sunstepper.c',
    'src/arkode/arkode_types_impl.h',
    'src/arkode/arkode_user_controller.c',
    'src/arkode/arkode_user_controller.h',
)

sunadaptcontroller_sources = files(
    'include/sunadaptcontroller/sunadaptcontroller_imexgus.h',
    'include/sunadaptcontroller/sunadaptcontroller_mrihtol.h',
    'include/sunadaptcontroller/sunadaptcontroller_soderlind.h',
    'src/sunadaptcontroller/imexgus/sunadaptcontroller_imexgus.c',
    'src/sunadaptcontroller/mrihtol/sunadaptcontroller_mrihtol.c',
    'src/sunadaptcontroller/soderlind/sunadaptcontroller_soderlind.c',
)

sundials_cvode_sources = files(
    'include/cvode/cvode.h',
    'include/cvode/cvode_bandpre.h',
    'include/cvode/cvode_bbdpre.h',
    'include/cvode/cvode_diag.h',
    'include/cvode/cvode_ls.h',
    'include/cvode/cvode_proj.h',
    'src/cvode/cvode.c',
    'src/cvode/cvode_bandpre.c',
    'src/cvode/cvode_bandpre_impl.h',
    'src/cvode/cvode_bbdpre.c',
    'src/cvode/cvode_bbdpre_impl.h',
    'src/cvode/cvode_diag.c',
    'src/cvode/cvode_diag_impl.h',
    'src/cvode/cvode_fused_stubs.c',
    'src/cvode/cvode_impl.h',
    'src/cvode/cvode_io.c',
    'src/cvode/cvode_ls.c',
    'src/cvode/cvode_ls_impl.h',
    'src/cvode/cvode_nls.c',
    'src/cvode/cvode_proj.c',
    'src/cvode/cvode_proj_impl.h',
)

sundials_nvector_sources = files(
    'include/nvector/nvector_manyvector.h',
    'include/nvector/nvector_pthreads.h',
    'include/nvector/nvector_serial.h',
    'src/nvector/manyvector/nvector_manyvector.c',
    'src/nvector/pthreads/nvector_pthreads.c',
    'src/nvector/serial/nvector_serial.c',
)

if omp_c.found()
    sundials_nvector_sources += files(
        'include/nvector/nvector_openmp.h',
        'src/nvector/openmp/nvector_openmp.c',
    )
endif

sundials_sunlinsol_sources = files(
    'include/sunlinsol/sunlinsol_band.h',
    'include/sunlinsol/sunlinsol_dense.h',
    'include/sunlinsol/sunlinsol_lapackband.h',
    'include/sunlinsol/sunlinsol_lapackdense.h',
    'include/sunlinsol/sunlinsol_pcg.h',
    'include/sunlinsol/sunlinsol_spbcgs.h',
    'include/sunlinsol/sunlinsol_spfgmr.h',
    'include/sunlinsol/sunlinsol_spgmr.h',
    'include/sunlinsol/sunlinsol_sptfqmr.h',
    'src/sunlinsol/band/sunlinsol_band.c',
    'src/sunlinsol/dense/sunlinsol_dense.c',
    'src/sunlinsol/lapackband/sunlinsol_lapackband.c',
    'src/sunlinsol/lapackdense/sunlinsol_lapackdense.c',
    'src/sunlinsol/pcg/sunlinsol_pcg.c',
    'src/sunlinsol/spbcgs/sunlinsol_spbcgs.c',
    'src/sunlinsol/spfgmr/sunlinsol_spfgmr.c',
    'src/sunlinsol/spgmr/sunlinsol_spgmr.c',
    'src/sunlinsol/sptfqmr/sunlinsol_sptfqmr.c',
)

sundials_sunmatrix_sources = files(
    'include/sunmatrix/sunmatrix_band.h',
    'include/sunmatrix/sunmatrix_dense.h',
    'include/sunmatrix/sunmatrix_sparse.h',
    'src/sunmatrix/band/sunmatrix_band.c',
    'src/sunmatrix/dense/sunmatrix_dense.c',
    'src/sunmatrix/sparse/sunmatrix_sparse.c',
)

sundials_sunnonlinsol_sources = files(
    'include/sunnonlinsol/sunnonlinsol_fixedpoint.h',
    'include/sunnonlinsol/sunnonlinsol_newton.h',
    'src/sunnonlinsol/fixedpoint/sunnonlinsol_fixedpoint.c',
    'src/sunnonlinsol/newton/sunnonlinsol_newton.c',
)

sundials_sources = files(
    'include/sundials/sundials_adaptcontroller.h',
    'include/sundials/sundials_band.h',
    'include/sundials/sundials_config.h',
    'include/sundials/sundials_context.h',
    'include/sundials/sundials_core.h',
    'include/sundials/sundials_dense.h',
    'include/sundials/sundials_direct.h',
    'include/sundials/sundials_errors.h',
    'include/sundials/sundials_export.h',
    'include/sundials/sundials_futils.h',
    'include/sundials/sundials_iterative.h',
    'include/sundials/sundials_linearsolver.h',
    'include/sundials/sundials_logger.h',
    'include/sundials/sundials_math.h',
    'include/sundials/sundials_matrix.h',
    'include/sundials/sundials_memory.h',
    'include/sundials/sundials_nonlinearsolver.h',
    'include/sundials/sundials_nvector.h',
    'include/sundials/sundials_nvector_senswrapper.h',
    'include/sundials/sundials_profiler.h',
    'include/sundials/sundials_stepper.h',
    'include/sundials/sundials_types.h',
    'include/sundials/sundials_types_deprecated.h',
    'include/sundials/sundials_version.h',
    'src/sundials/sundials_adaptcontroller.c',
    'src/sundials/sundials_adiak_metadata.h',
    'src/sundials/sundials_band.c',
    'src/sundials/sundials_context.c',
    'src/sundials/sundials_debug.h',
    'src/sundials/sundials_dense.c',
    'src/sundials/sundials_direct.c',
    'src/sundials/sundials_errors.c',
    'src/sundials/sundials_futils.c',
    'src/sundials/sundials_hashmap.c',
    'src/sundials/sundials_hashmap_impl.h',
    'src/sundials/sundials_iterative.c',
    'src/sundials/sundials_iterative_impl.h',
    'src/sundials/sundials_linearsolver.c',
    'src/sundials/sundials_logger.c',
    'src/sundials/sundials_logger_impl.h',
    'src/sundials/sundials_macros.h',
    'src/sundials/sundials_math.c',
    'src/sundials/sundials_matrix.c',
    'src/sundials/sundials_memory.c',
    'src/sundials/sundials_nonlinearsolver.c',
    'src/sundials/sundials_nvector.c',
    'src/sundials/sundials_nvector_senswrapper.c',
    'src/sundials/sundials_profiler.c',
    'src/sundials/sundials_stepper.c',
    'src/sundials/sundials_stepper_impl.h',
    'src/sundials/sundials_sycl.h',
    'src/sundials/sundials_utils.h',
    'src/sundials/sundials_version.c',
    'src/sundials/sundials_version.c',
)

inc = include_directories('include')
src_sundials = include_directories('src/sundials')

sundials_c_args = ['-DG_LOG_DOMAIN="NUMCOSMO"', '-DHAVE_CONFIG_H']

libsundials = static_library(
    'sundials',
    sources: [
        sundials_arkode_sources,
        sundials_cvode_sources,
        sundials_nvector_sources,
        sundials_sunlinsol_sources,
        sundials_sunmatrix_sources,
        sundials_sunnonlinsol_sources,
        sunadaptcontroller_sources,
        sundials_sources,
    ],
    dependencies: [omp_c],
    include_directories: [configinc, inc, src_sundials],
    install: false,
    gnu_symbol_visibility: 'hidden',
    c_args: sundials_c_args,
)