

levmar_sources = files(
    'Axb.c',
    'compiler.h',
    'levmar.h',
    'lm.h',
    'misc.c',
    'misc.h',
    'levmar.h',
    'lmbc.c',
    'lmblec.c',
    'lmbleic.c',
    'lm.c',
    'lmlec.c'
)

levmar_c_args = ['-DG_LOG_DOMAIN="NUMCOSMO"', '-DHAVE_CONFIG_H', '-Wno-implicit-fallthrough']

liblevmar = static_library(
    'levmar',
    sources: levmar_sources,
    dependencies: [omp_c],
    include_directories: [configinc, numcosmoinc],
    install: false,
    gnu_symbol_visibility: 'hidden',
    c_args: levmar_c_args,
)
