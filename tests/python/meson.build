if enable_gir
    common_env = {
        'GI_TYPELIB_PATH': numcosmo_build_dir / 'numcosmo',
        'LD_LIBRARY_PATH': numcosmo_build_dir / 'numcosmo',
        'DYLD_LIBRARY_PATH': numcosmo_build_dir / 'numcosmo',
        'PYTHONPATH': numcosmo_source_dir,
    }

    env = environment(common_env)

    env_sphere_map = environment(common_env)
    env_sphere_map.set('MALLOC_PERTURB_', '0')

    pytest = find_program('pytest', required: true)

    pytest_base_args = ['-ra', '--tb=short']
    pytest_protocol = 'exitcode'
    # We are using xdist for parallel test execution, so we set this to false to avoid
    # conflicts with pytest's internal parallelization mechanisms.
    pytest_is_parallel = false

    if get_option('b_coverage')
        python_module.find_installation(
            'python3',
            modules: ['pytest_cov'],
            required: true,
        )
        pytest_base_args += ['--cov', 'numcosmo_py']
        pytest_base_args += ['--cov-report', 'xml']
        pytest_base_args += ['--cov-append']
    endif

    py_tap = python_module.find_installation(
        'python3',
        modules: ['pytest_tap'],
        required: false,
    )
    if py_tap.found() and get_option('pytest_tap')
        r = run_command(
            python,
            '-c', 'import importlib.metadata; import sys; sys.stdout.write(importlib.metadata.version("pytest-tap"))',
            check: true,
        )
        message('pytest-tap version: @0@'.format(r.stdout().strip()))

        if r.stdout().version_compare('>= 3.4')
            pytest_base_args += ['--tap']
        else
            pytest_base_args += ['--tap-stream']
        endif
        pytest_protocol = 'tap'
    endif

    # If xdist is available, use it for parallel test execution
    py_xdist = python_module.find_installation(
        'python3',
        modules: ['xdist'],
        required: false,
    )
    if py_xdist.found()
        pytest_parallel_args = pytest_base_args + ['-n', 'auto']
    else
        pytest_parallel_args = pytest_base_args
    endif

    # Full Python test suite (excluding specialized markers)
    test(
        'pytest',
        pytest,
        args: pytest_parallel_args + [meson.current_source_dir()],
        depends: [numcosmo_gir, numcosmo_math_gir],
        env: env,
        timeout: 0,
        protocol: pytest_protocol,
        is_parallel: pytest_is_parallel,
        suite: ['python'],
    )

    # Powspec subset
    test(
        'pytest-powspec',
        pytest,
        args: pytest_parallel_args
        + ['--run-powspec', '-m', 'powspec', meson.current_source_dir()],
        depends: [numcosmo_gir, numcosmo_math_gir],
        env: env,
        timeout: 0,
        protocol: pytest_protocol,
        is_parallel: pytest_is_parallel,
        suite: ['py-powspec'],
    )

    # Xcor subset
    test(
        'pytest-xcor',
        pytest,
        args: pytest_parallel_args
        + ['--run-xcor', '-m', 'xcor', meson.current_source_dir()],
        depends: [numcosmo_gir, numcosmo_math_gir],
        env: env,
        timeout: 0,
        protocol: pytest_protocol,
        is_parallel: pytest_is_parallel,
        suite: ['py-xcor'],
    )

    # Sphere_map subset
    test(
        'pytest-sphere_map',
        pytest,
        args: pytest_parallel_args
        + ['--run-sphere-map', '-m', 'sphere_map', meson.current_source_dir()],
        depends: [numcosmo_gir, numcosmo_math_gir],
        env: env_sphere_map,
        timeout: 0,
        protocol: pytest_protocol,
        is_parallel: pytest_is_parallel,
        suite: ['python'],
    )

    # MPI tests (if available)
    if mpi_c_dep.found()
        mpiexec = find_program('mpiexec', required: false)
        if mpiexec.found()
            test(
                'pytest-mpi',
                mpiexec,
                args: ['-n', '2', pytest.full_path()]
                + pytest_base_args
                + ['--run-mpi', '-m', 'mpi', meson.current_source_dir()],
                depends: [numcosmo_gir, numcosmo_math_gir],
                env: env,
                timeout: 0,
                is_parallel: pytest_is_parallel,
                suite: ['python', 'mpi'],
            )
        endif
    endif
endif