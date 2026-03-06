if enable_gir
    env = environment()
    env.prepend(
        'GI_TYPELIB_PATH',
        numcosmo_build_dir / 'numcosmo',
    )

    env.prepend(
        'LD_LIBRARY_PATH',
        numcosmo_build_dir / 'numcosmo',
    )

    env.prepend(
        'DYLD_LIBRARY_PATH',
        numcosmo_build_dir / 'numcosmo',
    )

    env.prepend(
        'PYTHONPATH',
        numcosmo_source_dir,
    )

    pytest = find_program('pytest', required: true)

    pytest_base_args = ['-ra', '--tb=short']
    pytest_protocol = 'exitcode'
    pytest_is_parallel = true

    if get_option('b_coverage')
        python_module.find_installation(
            'python3',
            modules: ['pytest_cov'],
            required: true,
        )
        pytest_base_args += ['--cov', 'numcosmo_py']
        pytest_base_args += ['--cov-report', 'xml']
        pytest_base_args += ['--cov-append']
        pytest_is_parallel = false
    endif

    py_tap = python_module.find_installation(
        'python3',
        modules: ['pytest_tap'],
        required: false,
    )
    if py_tap.found() and get_option('pytest_tap')
        r = run_command(
            python,
            '-c',
            'import importlib.metadata; import sys; sys.stdout.write(importlib.metadata.version("pytest-tap"))',
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

    # Full Python test suite
    test(
        'pytest',
        pytest,
        args: pytest_base_args + [meson.current_source_dir()],
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
        args: pytest_base_args + ['-m', 'powspec', meson.current_source_dir()],
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
        args: pytest_base_args + ['-m', 'xcor', meson.current_source_dir()],
        depends: [numcosmo_gir, numcosmo_math_gir],
        env: env,
        timeout: 0,
        protocol: pytest_protocol,
        is_parallel: pytest_is_parallel,
        suite: ['py-xcor'],
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
                + ['-m', 'mpi', meson.current_source_dir()],
                depends: [numcosmo_gir, numcosmo_math_gir],
                env: env,
                timeout: 0,
                suite: ['python', 'mpi'],
            )
        endif
    endif
endif
