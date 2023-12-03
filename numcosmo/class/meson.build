
class_sources = files(
    'hyrec/helium.c',
    'hyrec/helium.h',
    'hyrec/history.c',
    'hyrec/history.h',
    'hyrec/hydrogen.c',
    'hyrec/hydrogen.h',
    'hyrec/hyrec.h',
    'hyrec/hyrectools.c',
    'hyrec/hyrectools.h',
    'include/arrays.h',
    'include/background.h',
    'include/class.h',
    'include/common.h',
    'include/dei_rkck.h',
    'include/evolver_ndf15.h',
    'include/evolver_rkck.h',
    'include/growTable.h',
    'include/hyperspherical.h',
    'include/input.h',
    'include/lensing.h',
    'include/nonlinear.h',
    'include/output.h',
    'include/parser.h',
    'include/perturbations.h',
    'include/primordial.h',
    'include/quadrature.h',
    'include/sparse.h',
    'include/spectra.h',
    'include/svnversion.h',
    'include/thermodynamics.h',
    'include/transfer.h',
    'source/background.c',
    'source/input.c',
    'source/lensing.c',
    'source/nonlinear.c',
    'source/perturbations.c',
    'source/primordial.c',
    'source/spectra.c',
    'source/thermodynamics.c',
    'source/transfer.c',
    'tools/arrays.c',
    'tools/common.c',
    'tools/dei_rkck.c',
    'tools/evolver_ndf15.c',
    'tools/evolver_rkck.c',
    'tools/growTable.c',
    'tools/hermite3_interpolation_csource.h',
    'tools/hermite4_interpolation_csource.h',
    'tools/hermite6_interpolation_csource.h',
    'tools/hyperspherical.c',
    'tools/parser.c',
    'tools/quadrature.c',
    'tools/sparse.c',
)

class_c_args = [
    '-DG_LOG_DOMAIN="NUMCOSMO"',
    '-DHAVE_CONFIG_H',
    '-D__CLASSDIR__="@0@/class_data"'.format(numcosmo_pkgdatadir),
]

disable_warnings = ['-Wno-sign-compare', '-Wno-implicit-fallthrough']

foreach flag : disable_warnings
    if cc.has_argument(flag)
        class_c_args += [flag]
    endif
endforeach

inc = include_directories('include')

libclass = static_library(
    'class',
    sources: class_sources,
    dependencies: [omp_c, glib_dep],
    include_directories: [configinc, numcosmoinc, inc],
    install: false,
    gnu_symbol_visibility: 'hidden',
    c_args: class_c_args,
)
