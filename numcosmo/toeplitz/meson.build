

toeplitz_sources = files(
    'solvers/dcommonsub.c',
    'solvers/dcommonsub.h',
    'solvers/dlev.c',
    'solvers/dlevcl.c',
    'solvers/dlevin.c',
    'solvers/dopt_k.c',
    'solvers/dopt_k.h',
    'solvers/drefine.c',
    'solvers/drefine.h',
    'solvers/dschur.c',
    'solvers/dschurin.c',
    'solvers/dschurinup.c',
    'solvers/dschurlev.c',
    'solvers/dschurlevin.c',
    'solvers/dschurup.c',
    'solvers/fblas.h',
    'solvers/ftypes.h',
    'solvers/polytool.c',
    'solvers/polytool.h',
    'solvers/toeplitz.h',
    'solvers/toetool.c',
    'solvers/toetool.h',
)

toeplitz_c_args = ['-DG_LOG_DOMAIN="NUMCOSMO"', '-DHAVE_CONFIG_H', '-Wno-sign-compare']

libtoeplitz = static_library(
    'toeplitz',
    sources: toeplitz_sources,
    dependencies: [omp_c, glib_dep],
    include_directories: [configinc, numcosmoinc],
    install: false,
    gnu_symbol_visibility: 'hidden',
    c_args: toeplitz_c_args,
)
