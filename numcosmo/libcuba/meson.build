


common_inc = include_directories('common')
cuba_c_args = [
    '-DG_LOG_DOMAIN="NUMCOSMO"',
    '-DHAVE_CONFIG_H',
    '-DNOUNDERSCORE',
]
if not cc_supports_c99_vla
    cuba_c_args += ['-DNDIM=@0@'.format(get_option('libcuba-max-dim'))]
    cuba_c_args += ['-DNCOMP=@0@'.format(get_option('libcuba-max-comp'))]
endif

disable_warnings = ['-Wno-sign-compare', '-Wno-clobbered']

foreach flag : disable_warnings
    if cc.has_argument(flag)
        cuba_c_args += flag
    endif
endforeach



cuhre_sources = files(
    'cuhre/Cuhre.c',
    'cuhre/decl.h',
)

cuhre_inc = include_directories('cuhre')

libcuba_cuhre = static_library(
    'cuba_cuhre',
    sources: cuhre_sources,
    dependencies: [],
    include_directories: [configinc, cuhre_inc, common_inc],
    install: false,
    gnu_symbol_visibility: 'hidden',
    c_args: cuba_c_args,
)

divonne_sources = files(
    'divonne/Divonne.c',
    'divonne/decl.h',
)

divonne_inc = include_directories('divonne')

libcuba_divonne = static_library(
    'cuba_divonne',
    sources: divonne_sources,
    dependencies: [],
    include_directories: [configinc, divonne_inc, common_inc],
    install: false,
    gnu_symbol_visibility: 'hidden',
    c_args: cuba_c_args,
)

suave_sources = files(
    'suave/Suave.c',
    'suave/decl.h',
)

suave_inc = include_directories('suave')

libcuba_suave = static_library(
    'cuba_suave',
    sources: suave_sources,
    dependencies: [],
    include_directories: [configinc, suave_inc, common_inc],
    install: false,
    gnu_symbol_visibility: 'hidden',
    c_args: cuba_c_args,
)

vegas_sources = files(
    'vegas/Vegas.c',
    'vegas/decl.h',
)

vegas_inc = include_directories('vegas')

libcuba_vegas = static_library(
    'cuba_vegas',
    sources: vegas_sources,
    dependencies: [],
    include_directories: [configinc, vegas_inc, common_inc],
    install: false,
    gnu_symbol_visibility: 'hidden',
    c_args: cuba_c_args,
)

libcuba_sources = files(
    'common/Data.c',
    'common/Fork.c',
    'common/Global.c',
    'cuba.h',
)

libcuba_inc = include_directories('.')

libcuba = static_library(
    'cuba',
    sources: libcuba_sources,
    dependencies: [],
    include_directories: [configinc, libcuba_inc],
    link_with: [
        libcuba_cuhre,
        libcuba_divonne,
        libcuba_suave,
        libcuba_vegas,
    ],
    install: false,
    gnu_symbol_visibility: 'hidden',
    c_args: cuba_c_args,
)
