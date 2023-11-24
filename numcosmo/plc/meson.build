

plc_sources = files(
    'cldf/cldf.c',
    'cldf/cldf.h',
    'minipmc/distribution.c',
    'minipmc/errorlist.c',
    'minipmc/errorlist.h',
    'minipmc/io.c',
    'minipmc/io.h',
    'minipmc/maths_base.h',
    'minipmc/pmc.h',
    'clik_dic.c',
    'clik_dic.h',
    'lapack_clik.h',
)

plc_inc = include_directories(['.', 'cldf', 'minipmc', 'plik'])

plc_c_args = [
    '-DG_LOG_DOMAIN="NUMCOSMO"',
    '-DHAVE_CONFIG_H',
    '-DHAS_LAPACK',
    '-DLAPACK_CLIK',
    '-DNOHEALPIX',
    '-DCLIK_LENSING',
    '-DCAMSPEC_V1',
    '-DCLIKSVNVERSION="0031e3604f0c"',
]

if numcosmo_conf.has('HAVE_RTLD_DEFAULT')
    plc_c_args += '-DHAS_RTLD_DEFAULT'
endif

disable_warnings = ['-Wno-implicit-fallthrough', '-Wno-sign-compare']

foreach flag : disable_warnings
    if cc.has_argument(flag)
        plc_c_args += flag
    endif
endforeach

libplc = static_library(
    'plc',
    sources: plc_sources,
    dependencies: [omp_c],
    include_directories: [configinc, numcosmoinc, plc_inc],
    install: false,
    gnu_symbol_visibility: 'hidden',
    c_args: plc_c_args,
)
