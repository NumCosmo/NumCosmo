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
    {
        'name': 'ncm_fit',
        'sources': ['ncm/fit/test_ncm_fit.c'],
        'priority': 6,
        'timeout': 0,
    },
    {
        'name': 'ncm_fit_state',
        'sources': ['ncm/fit/test_ncm_fit_state.c'],
    },
    {
        'name': 'ncm_fit_esmcmc',
        'sources': ['ncm/fit/test_ncm_fit_esmcmc.c'],
        'timeout': 0,
        'suite': ['c-acceptance', 'omp'],
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
    {
        'name': 'nc_cbe',
        'sources': ['nc/cmb/test_nc_cbe.c'],
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
    # The cluster_wl test is one combinatorial 24x6 executable; it is built twice from the
    # same source (via -D, see the C file) so meson --slice can balance CI wall time. The
    # cheap property checks are a normal unit test; the heavy MC/fit half is acceptance.
    {
        'name': 'nc_data_cluster_wl_prop',
        'sources': ['nc/data/test_nc_data_cluster_wl.c'],
        'c_args': ['-DCLUSTER_WL_SPLIT_CHEAP'],
        'timeout': 0,
        # m2lnP (the OpenMP-parallel likelihood eval) is in the cheap half, so this exec
        # also exercises the OpenMP path and stays in the omp lane.
        'suite': ['c', 'omp'],
    },
    {
        'name': 'nc_data_cluster_wl',
        'sources': ['nc/data/test_nc_data_cluster_wl.c'],
        'c_args': ['-DCLUSTER_WL_SPLIT_EXPENSIVE'],
        'timeout': 0,
        'suite': ['c-acceptance', 'omp'],
    },
    # Truth-table validation of the integration methods: regenerates a curated
    # sample from a fixed seed (nothing stored) and checks every method against a
    # recomputed FIXED-320x7 reference. Deterministic and fast.
    {
        'name': 'nc_data_cluster_wl_truth',
        'sources': ['nc/data/test_nc_data_cluster_wl_truth.c'],
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
        'name': 'nc_galaxy_sd_position',
        'sources': ['nc/lss/galaxy/test_nc_galaxy_sd_position.c'],
    },
    {
        'name': 'nc_galaxy_sd_true_redshift',
        'sources': ['nc/lss/galaxy/test_nc_galaxy_sd_true_redshift.c'],
    },
    {
        'name': 'nc_galaxy_sd_obs_redshift',
        'sources': ['nc/lss/galaxy/test_nc_galaxy_sd_obs_redshift.c'],
    },
    {
        'name': 'nc_galaxy_sd_shape',
        'sources': ['nc/lss/galaxy/test_nc_galaxy_sd_shape.c'],
    },
]

# Combine all test arrays
all_c_tests = math_tests + model_tests + data_tests + lss_tests + galaxy_tests

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
        executable_test,
        args: test_seed_args,
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