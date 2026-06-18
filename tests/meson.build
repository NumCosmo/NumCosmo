# OpenMP / BLAS thread policy for the test lanes (see TESTING.md).
#
# Tests run concurrently (is_parallel: true) must pin every threading backend to a
# single thread, otherwise N concurrent processes x N threads each = N^2
# oversubscription, which on the mixed OpenBLAS stack can deadlock. Tests that run
# alone (is_parallel: false) instead get OMP_NUM_THREADS = online cores, so their
# `#pragma omp parallel` paths are actually exercised inside the normal `meson test`
# run -- no separate OMP lane or external env vars required. This is selected per
# test below from the `omp` suite tag (C) / the non-xdist pytest lane (Python).
getconf_prog = find_program('getconf', required: false)
if getconf_prog.found()
    n_omp_threads = run_command(getconf_prog, '_NPROCESSORS_ONLN', check: true).stdout().strip()
else
    n_omp_threads = '4'
endif

omp_vars_serial = {
    'OMP_NUM_THREADS': '1',
    'OMP_MAX_ACTIVE_LEVELS': '1',
    'OMP_THREAD_LIMIT': '1',
    'OPENBLAS_NUM_THREADS': '1',
    'BLIS_NUM_THREADS': '1',
    'MKL_NUM_THREADS': '1',
}
omp_vars_threaded = {
    'OMP_NUM_THREADS': n_omp_threads,
    'OMP_MAX_ACTIVE_LEVELS': '1',
    'OMP_THREAD_LIMIT': n_omp_threads,
}

# Reproducibility (see TESTING.md). The default (push-CI) lane pins a fixed g_test seed so
# every run is bit-identical; the flaky lane (-Dflaky_tests=true, run weekly with
# `meson test --repeat=N`) uses a fresh random seed each run to stress the seeded/statistical
# assertions across many draws. C tests parse `--seed`; Python tests seed their own RNGs.
if get_option('flaky_tests')
    test_seed_args = []
else
    test_seed_args = ['--seed=R02S0123456789abcdef0123456789abcdef']
endif

# C tests are now in a separate subdirectory
subdir('c')

# Python tests are now in a separate subdirectory
subdir('python')