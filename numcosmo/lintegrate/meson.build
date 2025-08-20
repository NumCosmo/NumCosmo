lintegrate_c_args = [
    '-DG_LOG_DOMAIN="NUMCOSMO"',
    '-DHAVE_CONFIG_H',
]

disable_warnings = []

foreach flag : disable_warnings
    if cc.has_argument(flag)
        lintegrate_c_args += flag
    endif
endforeach

lintegrate_sources = files(
    'lintegrate.h',
    'lintegrate_cquad.c',
    'lintegrate_qag.c',
    'lintegrate_qng.c',
    'qng.h',
)

liblintegrate = static_library(
    'libintegrate',
    sources: lintegrate_sources,
    dependencies: [gsl_dep],
    install: false,
    gnu_symbol_visibility: 'hidden',
    c_args: lintegrate_c_args,
)