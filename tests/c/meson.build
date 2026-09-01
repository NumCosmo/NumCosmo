# Math tests (ncm_*)
math_tests = [
    {
        'name': 'ncm_cfg',
        'sources': ['ncm/core/test_ncm_cfg.c'],
    },
    {
        'name': 'ncm_generic',
        'sources': ['ncm/core/test_ncm_generic.c'],
    },
    {
        'name': 'ncm_sparam',
        'sources': ['ncm/model/test_ncm_sparam.c'],
    },
    {
        'name': 'ncm_diff',
        'sources': ['ncm/integration/test_ncm_diff.c'],
    },
    {
        'name': 'ncm_dtuple',
        'sources': ['ncm/core/test_ncm_dtuple.c'],
    },
    {
        'name': 'ncm_catalog',
        'sources': ['ncm/data/test_ncm_catalog.c'],
    },
    {
        'name': 'ncm_sky_footprint',
        'sources': ['ncm/sphere/test_ncm_sky_footprint.c'],
    },
    {
        'name': 'ncm_fftlog',
        'sources': ['ncm/fftlog/test_ncm_fftlog.c'],
        'priority': 2,
        'timeout': 0,
    },
    {
        'name': 'ncm_vector',
        'sources': ['ncm/algebra/test_ncm_vector.c'],
    },
    {
        'name': 'ncm_poly_roots',
        'sources': ['ncm/algebra/test_ncm_poly_roots.c'],
    },
    {
        'name': 'ncm_laurent_series',
        'sources': ['ncm/algebra/test_ncm_laurent_series.c'],
    },
    {
        'name': 'ncm_matrix',
        'sources': ['ncm/algebra/test_ncm_matrix.c'],
    },
    {
        'name': 'ncm_stats_vec',
        'sources': ['ncm/stats/test_ncm_stats_vec.c'],
    },
    {
        'name': 'ncm_stats_dist1d_epdf',
        'sources': ['ncm/stats/test_ncm_stats_dist1d_epdf.c'],
    },
    {
        'name': 'ncm_stats_dist',
        'sources': ['ncm/stats/test_ncm_stats_dist.c'],
        'timeout': 0,
        'suite': ['c-statistical', 'omp'],
    },
    {
        'name': 'ncm_stats_dist_kernel',
        'sources': ['ncm/stats/test_ncm_stats_dist_kernel.c'],
    },
    {
        'name': 'ncm_spline',
        'sources': ['ncm/spline/test_ncm_spline.c'],
        'priority': 3,
    },
    {
        'name': 'ncm_spline2d',
        'sources': ['ncm/spline/test_ncm_spline2d.c'],
    },
    {
        'name': 'ncm_spline_vec',
        'sources': ['ncm/spline/test_ncm_spline_vec.c'],
    },
    # Bit-exact agreement between serial and concurrent evaluation of one shared
    # prepared NcmSplineBSpline: its GSL workspace is per-call scratch, so this is
    # the regression test for evaluating it from parallel kernel integrands.
    {
        'name': 'ncm_spline_bspline_threads',
        'sources': ['ncm/spline/test_ncm_spline_bspline_threads.c'],
        'suite': ['c', 'omp'],
    },
    {
        'name': 'ncm_ode_spline',
        'sources': ['ncm/spline/test_ncm_ode_spline.c'],
    },
    {
        'name': 'ncm_integral1d',
        'sources': ['ncm/integration/test_ncm_integral1d.c'],
    },
    {
        'name': 'ncm_integral_nd',
        'sources': ['ncm/integration/test_ncm_integral_nd.c'],
    },
    {
        'name': 'ncm_integral_fixed',
        'sources': ['ncm/integration/test_ncm_integral_fixed.c'],
    },
    {
        'name': 'ncm_sf_sbessel',
        'sources': ['ncm/specfunc/test_ncm_sf_sbessel.c'],
    },
    {
        'name': 'ncm_mpsf_0F1',
        'sources': ['ncm/specfunc/test_ncm_mpsf_0F1.c'],
    },
    {
        'name': 'ncm_mpsf_trig_int',
        'sources': ['ncm/specfunc/test_ncm_mpsf_trig_int.c'],
    },
    {
        'name': 'ncm_model',
        'sources': [
            'ncm/model/test_ncm_model.c',
            'ncm/model/ncm_model_test.c',
            'ncm/model/ncm_model_test.h',
        ],
    },
    {
        'name': 'ncm_model_ctrl',
        'sources': ['ncm/model/test_ncm_model_ctrl.c'],
    },
    {
        'name': 'ncm_mset',
        'sources': ['ncm/model/test_ncm_mset.c'],
    },
    {
        'name': 'ncm_mset_catalog',
        'sources': ['ncm/fit/test_ncm_mset_catalog.c'],
        'priority': 1,
        'timeout': 0,
    },
    {
        'name': 'ncm_obj_array',
        'sources': ['ncm/core/test_ncm_obj_array.c'],
    },
    {
        'name': 'ncm_serialize',
        'sources': ['ncm/core/test_ncm_serialize.c'],
    },
    # One executable per NcmFit algorithm. The checks are identical for all of
    # them and live in test_ncm_fit_common.c, which every one of these links; the
    # per-algorithm source only names its algorithm. Splitting this far keeps a
    # failure attributable without reading the whole log, lets each carry a
    # timeout its own size, and gives `meson test --slice` units it can balance
    # by count -- as one binary these were 177 s, which no amount of
    # --num-processes could divide.
    {
        'name': 'ncm_fit_nlopt_neldermead',
        'sources': ['ncm/fit/test_ncm_fit_nlopt_neldermead.c', 'ncm/fit/test_ncm_fit_common.c'],
        'priority': 6,
        'timeout': 0,
    },
    {
        'name': 'ncm_fit_nlopt_slsqp',
        'sources': ['ncm/fit/test_ncm_fit_nlopt_slsqp.c', 'ncm/fit/test_ncm_fit_common.c'],
        'priority': 6,
        'timeout': 0,
    },
    {
        'name': 'ncm_fit_gsl_ls',
        'sources': ['ncm/fit/test_ncm_fit_gsl_ls.c', 'ncm/fit/test_ncm_fit_common.c'],
        'priority': 6,
        'timeout': 0,
    },
    {
        'name': 'ncm_fit_gsl_mm_conjugate_fr',
        'sources': ['ncm/fit/test_ncm_fit_gsl_mm_conjugate_fr.c', 'ncm/fit/test_ncm_fit_common.c'],
        'priority': 6,
        'timeout': 0,
    },
    {
        'name': 'ncm_fit_gsl_mm_conjugate_pr',
        'sources': ['ncm/fit/test_ncm_fit_gsl_mm_conjugate_pr.c', 'ncm/fit/test_ncm_fit_common.c'],
        'priority': 6,
        'timeout': 0,
    },
    {
        'name': 'ncm_fit_gsl_mm_vector_bfgs',
        'sources': ['ncm/fit/test_ncm_fit_gsl_mm_vector_bfgs.c', 'ncm/fit/test_ncm_fit_common.c'],
        'priority': 6,
        'timeout': 0,
    },
    {
        'name': 'ncm_fit_gsl_mm_vector_bfgs2',
        'sources': ['ncm/fit/test_ncm_fit_gsl_mm_vector_bfgs2.c', 'ncm/fit/test_ncm_fit_common.c'],
        'priority': 6,
        'timeout': 0,
    },
    {
        'name': 'ncm_fit_gsl_mm_steepest_descent',
        'sources': ['ncm/fit/test_ncm_fit_gsl_mm_steepest_descent.c', 'ncm/fit/test_ncm_fit_common.c'],
        'priority': 6,
        'timeout': 0,
    },
    {
        'name': 'ncm_fit_gsl_nmsimplex',
        'sources': ['ncm/fit/test_ncm_fit_gsl_nmsimplex.c', 'ncm/fit/test_ncm_fit_common.c'],
        'priority': 6,
        'timeout': 0,
    },
    {
        'name': 'ncm_fit_gsl_nmsimplex2',
        'sources': ['ncm/fit/test_ncm_fit_gsl_nmsimplex2.c', 'ncm/fit/test_ncm_fit_common.c'],
        'priority': 6,
        'timeout': 0,
    },
    {
        'name': 'ncm_fit_gsl_nmsimplex2rand',
        'sources': ['ncm/fit/test_ncm_fit_gsl_nmsimplex2rand.c', 'ncm/fit/test_ncm_fit_common.c'],
        'priority': 6,
        'timeout': 0,
    },
    {
        'name': 'ncm_fit_levmar_der',
        'sources': ['ncm/fit/test_ncm_fit_levmar_der.c', 'ncm/fit/test_ncm_fit_common.c'],
        'priority': 6,
        'timeout': 0,
    },
    {
        'name': 'ncm_fit_levmar_dif',
        'sources': ['ncm/fit/test_ncm_fit_levmar_dif.c', 'ncm/fit/test_ncm_fit_common.c'],
        'priority': 6,
        'timeout': 0,
    },
    {
        'name': 'ncm_fit_levmar_bc_der',
        'sources': ['ncm/fit/test_ncm_fit_levmar_bc_der.c', 'ncm/fit/test_ncm_fit_common.c'],
        'priority': 6,
        'timeout': 0,
    },
    {
        'name': 'ncm_fit_levmar_bc_dif',
        'sources': ['ncm/fit/test_ncm_fit_levmar_bc_dif.c', 'ncm/fit/test_ncm_fit_common.c'],
        'priority': 6,
        'timeout': 0,
    },
    {
        'name': 'ncm_fit_state',
        'sources': ['ncm/fit/test_ncm_fit_state.c'],
    },
    {
        'name': 'ncm_mset_trans_kern_gauss',
        'sources': ['ncm/fit/test_ncm_mset_trans_kern_gauss.c'],
    },
    # Split by claim, not by check: driving a chain is mechanics a short run covers, while
    # requiring the sample covariance to reproduce the data covariance converges as
    # 1/sqrt(n) and is chased by loops that double the chain length up to fifteen times.
    # The mechanics run in the instrumented lane at a fixed short length; the recovery
    # assertions, with their retries, run in the statistical tier.
    {
        'name': 'ncm_fit_esmcmc',
        'sources': [
            'ncm/fit/test_ncm_fit_esmcmc.c',
            'ncm/fit/test_ncm_fit_esmcmc_common.c',
        ],
        'timeout': 0,
        'suite': ['c', 'omp'],
    },
    {
        'name': 'ncm_fit_esmcmc_recovery',
        'sources': [
            'ncm/fit/test_ncm_fit_esmcmc_recovery.c',
            'ncm/fit/test_ncm_fit_esmcmc_common.c',
        ],
        'timeout': 0,
        'suite': ['c-statistical', 'omp'],
    },
    {
        'name': 'ncm_func_eval',
        'sources': ['ncm/core/test_ncm_func_eval.c'],
    },
    {
        'name': 'ncm_sf_spherical_harmonics',
        'sources': ['ncm/specfunc/test_ncm_sf_spherical_harmonics.c'],
    },
    {
        'name': 'ncm_sphere_map',
        'sources': ['ncm/sphere/test_ncm_sphere_map.c'],
        'timeout': 0,
    },
    {
        'name': 'ncm_util',
        'sources': ['ncm/core/test_ncm_util.c'],
    },
    {
        'name': 'ncm_quaternion',
        'sources': ['ncm/algebra/test_ncm_quaternion.c'],
    },
]

# Model tests (NcHICosmo, recombination, distance, CLASS backend)
model_tests = [
    {
        'name': 'nc_hicosmo_de',
        'sources': ['nc/background/test_nc_hicosmo_de.c'],
    },
    {
        'name': 'nc_hicosmo_qgw',
        'sources': ['nc/background/test_nc_hicosmo_qgw.c'],
    },
    {
        'name': 'nc_recomb',
        'sources': ['nc/recomb/test_nc_recomb.c'],
    },
    # Split by cost: Cls and calc_ps run CLASS deep enough to take 12.3 s and 4.4 s
    # per model, where every other check is 0.4 s or less. Measured per line those
    # two are 91% of this file's runtime for about 2% of the lines only they
    # reach, so they run in the acceptance tier and the cheap checks stay in the
    # instrumented lane. Split by check rather than by sharing a prepared state,
    # which would couple the subtests: calc_ps mutates the object with
    # nc_cbe_set_calc_transfer().
    {
        'name': 'nc_cbe',
        'sources': ['nc/cmb/test_nc_cbe.c', 'nc/cmb/test_nc_cbe_common.c'],
        'priority': 4,
        'timeout': 0,
    },
    {
        'name': 'nc_cbe_class',
        'sources': ['nc/cmb/test_nc_cbe_class.c', 'nc/cmb/test_nc_cbe_common.c'],
        'suite': ['c-acceptance'],
        'priority': 4,
        'timeout': 0,
    },
    {
        'name': 'nc_hiqg_1d',
        'sources': ['nc/quantum/test_nc_hiqg_1d.c'],
    },
    {
        'name': 'nc_distance',
        'sources': ['nc/background/test_nc_distance.c'],
    },
]

# Data tests
data_tests = [
    {
        'name': 'ncm_data_gauss_cov',
        'sources': [
            'ncm/data/test_ncm_data_gauss_cov.c',
            'ncm/data/ncm_data_gauss_cov_test.c',
            'ncm/data/ncm_data_gauss_cov_test.h',
        ],
        'timeout': 0,
    },
    {
        'name': 'nc_data_bao_dtr_dhr',
        'sources': ['nc/data/test_nc_data_bao_dtr_dhr.c'],
    },
    {
        'name': 'nc_data_bao_rdv',
        'sources': ['nc/data/test_nc_data_bao_rdv.c'],
    },
    {
        'name': 'nc_data_bao_dvdv',
        'sources': ['nc/data/test_nc_data_bao_dvdv.c'],
    },
    {
        'name': 'nc_data_bao_dvr_dtdh',
        'sources': ['nc/data/test_nc_data_bao_dvr_dtdh.c'],
    },
    {
        'name': 'nc_data_cluster_ncount',
        'sources': ['nc/data/test_nc_data_cluster_ncount.c'],
    },
    {
        'name': 'nc_data_cluster_ncounts_gauss',
        'sources': ['nc/data/test_nc_data_cluster_ncounts_gauss.c'],
    },
    {
        'name': 'nc_data_hubble',
        'sources': ['nc/data/test_nc_data_hubble.c'],
    },
    # Truth-table validation of NcDataClusterWLFactor's integration methods:
    # regenerates a curated sample from a fixed seed (nothing stored) and
    # checks every method against a recomputed high-node self-converged
    # reference. Curated matrix across shape schemes (VarAdd/SeriesLensed/
    # Quad) and populations (Gauss/Beta), plus a direct SeriesLensed-vs-Quad
    # cross-check. Deterministic and fast.
    {
        'name': 'nc_data_cluster_wl_factor_truth',
        'sources': ['nc/data/test_nc_data_cluster_wl_factor_truth.c'],
        'suite': ['c', 'omp'],
    },
]

# Large-scale structure tests (halos, clusters, bias, transfer, power spectrum, window)
lss_tests = [
    {
        'name': 'nc_halo_catalog',
        'sources': ['nc/lss/halo/test_nc_halo_catalog.c'],
    },
    {
        'name': 'nc_halo_catalog_generator',
        'sources': ['nc/lss/halo/test_nc_halo_catalog_generator.c'],
    },
    {
        'name': 'nc_halo_catalog_member_generator',
        'sources': ['nc/lss/halo/test_nc_halo_catalog_member_generator.c'],
    },
    {
        'name': 'nc_halo_bias',
        'sources': ['nc/lss/halo/test_nc_halo_bias.c'],
    },
    {
        'name': 'nc_cluster_pseudo_counts',
        'sources': ['nc/lss/cluster/test_nc_cluster_pseudo_counts.c'],
    },
    {
        'name': 'nc_cluster_abundance',
        'sources': ['nc/lss/cluster/test_nc_cluster_abundance.c'],
    },
    {
        'name': 'nc_halo_density_profile',
        'sources': ['nc/lss/halo/test_nc_halo_density_profile.c'],
    },
    {
        'name': 'nc_halo_position',
        'sources': ['nc/lss/halo/test_nc_halo_position.c'],
    },
    {
        'name': 'nc_window',
        'sources': ['nc/powspec/test_nc_window.c'],
    },
    {
        'name': 'nc_transfer_func',
        'sources': ['nc/powspec/test_nc_transfer_func.c'],
    },
    # nc_powspec is built twice from one source (via -D, see the C file). Measured: the
    # CLASS-backed cbe half is fast (~4s, splined once) -> fast-lane unit; the analytic
    # transfer-function half (EH/BBKS + halofit + corr3d) is the slow one (>2 min) -> heavy
    # acceptance tier (coverage-only, unbounded timeout).
    {
        'name': 'nc_powspec_cbe',
        'sources': ['nc/powspec/test_nc_powspec.c'],
        'c_args': ['-DPOWSPEC_SPLIT_CBE'],
    },
    {
        'name': 'nc_powspec_transfer',
        'sources': ['nc/powspec/test_nc_powspec.c'],
        'c_args': ['-DPOWSPEC_SPLIT_TRANSFER'],
        'timeout': 0,
        'priority': 5,
        'suite': ['c-acceptance'],
    },
    {
        'name': 'nc_galaxy_acf',
        'sources': ['nc/lss/galaxy/test_nc_galaxy_acf.c'],
    },
    {
        'name': 'nc_wl_surface_mass_density',
        'sources': ['nc/lss/wl/test_nc_wl_surface_mass_density.c'],
    },
    {
        'name': 'nc_wl_ellipticity',
        'sources': ['nc/lss/wl/test_nc_wl_ellipticity.c'],
    },
    {
        'name': 'nc_wl_ellipticity_series',
        'sources': ['nc/lss/wl/test_nc_wl_ellipticity_series.c'],
    },
]

# Galaxy and weak lensing tests
galaxy_tests = [
    {
        'name': 'nc_galaxy_hod',
        'sources': ['nc/lss/galaxy/test_nc_galaxy_hod.c'],
    },
    {
        'name': 'nc_galaxy_wl_obs',
        'sources': ['nc/lss/galaxy/test_nc_galaxy_wl_obs.c'],
    },
    {
        'name': 'nc_galaxy_shape_intrinsic_mode',
        'sources': ['nc/lss/galaxy/test_nc_galaxy_shape_intrinsic_mode.c'],
    },
    {
        'name': 'nc_galaxy_shape_pop_series',
        'sources': ['nc/lss/galaxy/test_nc_galaxy_shape_pop_series.c'],
    },
    {
        'name': 'nc_galaxy_position_factor',
        'sources': ['nc/lss/galaxy/test_nc_galaxy_position_factor.c'],
    },
]

# Cross-correlation tests
xcor_tests = [
    {
        'name': 'nc_xcor',
        'sources': ['nc/xcor/test_nc_xcor.c'],
    },
    {
        'name': 'nc_xcor_kernel',
        'sources': ['nc/xcor/test_nc_xcor_kernel.c'],
    },
    {
        'name': 'nc_xcor_kquad',
        'sources': ['nc/xcor/test_nc_xcor_kquad.c'],
    },
    {
        'name': 'nc_xcor_kernel_analytic',
        'sources': ['nc/xcor/test_nc_xcor_kernel_analytic.c'],
    },
    {
        'name': 'nc_xcor_kernel_integrand',
        'sources': ['nc/xcor/test_nc_xcor_kernel_integrand.c'],
    },
    {
        'name': 'nc_xcor_solver',
        'sources': ['nc/xcor/test_nc_xcor_solver.c'],
    },
    {
        'name': 'nc_xcor_ssc_sij',
        'sources': ['nc/xcor/test_nc_xcor_ssc_sij.c'],
    },
]

# Combine all test arrays
all_c_tests = math_tests + model_tests + data_tests + lss_tests + galaxy_tests + xcor_tests

# Build and register all tests
tests_c_args = ['-DG_LOG_DOMAIN="NUMCOSMO"', '-DHAVE_CONFIG_H']

foreach c_test : all_c_tests
    extra_c_args = c_test.has_key('c_args') ? c_test['c_args'] : []
    executable_test = executable(
        c_test['name'],
        c_test['sources'],
        c_args: tests_c_args + extra_c_args,
        dependencies: [libnumcosmo_dep, gsl_dep],
        include_directories: [configinc, numcosmoinc],
    )
    if c_test.has_key('timeout')
        timeout = c_test['timeout']
    else
        timeout = 60
    endif
    if c_test.has_key('priority')
        priority = c_test['priority']
    else
        priority = 0
    endif
    if c_test.has_key('suite')
        suite = c_test['suite']
    else
        suite = ['c']
    endif
    # OMP-path tests run alone with all cores; the rest run concurrently, pinned to
    # a single thread (see tests/meson.build).
    test_is_omp = suite.contains('omp')
    test_env = [
        'G_TEST_SRCDIR=@0@'.format(meson.current_source_dir()),
        'G_TEST_BUILDDIR=@0@'.format(meson.current_build_dir()),
    ]
    foreach k, v : (test_is_omp ? omp_vars_threaded : omp_vars_serial)
        test_env += '@0@=@1@'.format(k, v)
    endforeach
    if test_is_omp
        priority = 100
    endif
    test(
        c_test['name'],
        test_is_omp ? omp_thread_wrapper : executable_test,
        args: test_is_omp ? [executable_test.full_path()]
        + test_seed_args : test_seed_args,
        depends: test_is_omp ? [executable_test] : [],
        timeout: timeout,
        env: test_env,
        protocol: 'tap',
        priority: priority,
        suite: suite,
        is_parallel: not test_is_omp,
    )
endforeach

# MPI tests
if mpi_c_dep.found()
    mpiexec = find_program('mpiexec', required: false)
    if mpiexec.found()

        c_tests_mpi = [
            {
                'name': 'ncm_fit_esmcmc_mpi',
                'sources': ['ncm/fit/test_ncm_fit_esmcmc_mpi.c'],
                'timeout': 0,
            },
        ]
        foreach c_test_mpi : c_tests_mpi
            executable_test = executable(
                c_test_mpi['name'],
                c_test_mpi['sources'],
                c_args: tests_c_args,
                dependencies: [libnumcosmo_dep],
                include_directories: [configinc, numcosmoinc],
            )
            if c_test_mpi.has_key('timeout')
                timeout = c_test_mpi['timeout']
            else
                timeout = 60
            endif
            if c_test_mpi.has_key('suite')
                suite = c_test_mpi['suite']
            else
                suite = ['c', 'mpi']
            endif

            mpiexec_args = (
                [
                    '-n', '2',
                    executable_test.full_path(),
                ]
                + test_seed_args
            )
            mpi_test_env = [
                'G_TEST_SRCDIR=@0@'.format(meson.current_source_dir()),
                'G_TEST_BUILDDIR=@0@'.format(meson.current_build_dir()),
            ]
            foreach k, v : omp_vars_serial
                mpi_test_env += '@0@=@1@'.format(k, v)
            endforeach
            test(
                c_test_mpi['name'],
                mpiexec,
                args: mpiexec_args,
                depends: [executable_test],
                timeout: timeout,
                env: mpi_test_env,
                protocol: 'tap',
                suite: suite,
            )
        endforeach
    endif
endif