# Offline reference generators.
#
# These are not tests and are never installed. They produce the certified values
# the test suite compares against, using arbitrary-precision arithmetic that is
# not a NumCosmo dependency.
#
# arb_dep/flint_dep come from the top-level meson.build, which already probes
# for FLINT (setting HAVE_FLINT) and copes with arb.h living either at the top
# level or under flint/. Without them the tools are simply not built, and the
# suites that would use their output fall back to the closed-form checks
# written in the tests themselves.
mathlib_dep = cc.find_library('m', required: false)

if arb_dep.found()
    # Build-checked only, not registered as a test. The comparison driver
    # (compare_powspec_analytic_arb.py) compiles its own copy through
    # pkg-config so it can run against an installed NumCosmo; this target
    # exists so a compile break is caught wherever FLINT is available.
    # FLINT's own headers trip -Wstringop-overread and -Warray-bounds under GCC 16,
    # which reads the interior pointers of an acb_t as if they sized the whole struct.
    # The diagnostics are about flint/acb.h, not about anything here, so they are
    # filtered rather than read past.
    arb_tool_args = cc.get_supported_arguments(
        '-Wno-stringop-overread',
        '-Wno-array-bounds',
    )

    foreach tool : ['ncm_powspec_analytic_arb', 'nc_xcor_kernel_analytic_arb', 'nc_xcor_kquad_arb']
        executable(
            tool,
            tool + '.c',
            c_args: arb_tool_args,
            dependencies: [arb_dep, flint_dep, mathlib_dep],
            install: false,
            build_by_default: true,
        )
    endforeach
else
    message('FLINT not found: reference generators in tests/tools are not built.')
endif
