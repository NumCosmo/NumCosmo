

levmar_sources = files(
    'Axb.c',
    'compiler.h',
    'levmar.h',
    'levmar.h',
    'lm.c',
    'lm.h',
    'lmbc.c',
    'lmblec.c',
    'lmbleic.c',
    'lmlec.c',
    'misc.c',
    'misc.h',
)

levmar_c_args = ['-DG_LOG_DOMAIN="NUMCOSMO"', '-DHAVE_CONFIG_H']

disable_warnings = ['-Wno-implicit-fallthrough']

foreach flag : disable_warnings
    if cc.has_argument(flag)
        levmar_c_args += flag
    endif
endforeach

liblevmar = static_library(
    'levmar',
    sources: levmar_sources,
    dependencies: [omp_c],
    include_directories: [configinc, numcosmoinc],
    install: false,
    gnu_symbol_visibility: 'hidden',
    c_args: levmar_c_args,
)
