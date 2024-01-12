

misc_sources = files(
    'Faddeeva.c',
    'Faddeeva.h',
    'LowRankQP.c',
    'clencurt.h',
    'converged.h',
    'cqpminimizer.c',
    'cubature.h',
    'gsl_cqp.h',
    'hcubature.c',
    'initial_point.c',
    'initial_point.h',
    'kdtree.c',
    'kdtree.h',
    'libqp.h',
    'libqp_gsmo.c',
    'libqp_splx.c',
    'mg_pdip.c',
    'nnls.c',
    'pcubature.c',
    'rb_knn_list.c',
    'rb_knn_list.h',
    'rcm.c',
    'rcm.h',
    'vwrapper.h',
)

misc_c_args = ['-DG_LOG_DOMAIN="NUMCOSMO"', '-DHAVE_CONFIG_H', '-Wno-sign-compare']

libmisc = static_library(
    'misc',
    sources: misc_sources,
    dependencies: [omp_c, gsl_dep, glib_dep],
    include_directories: [configinc, numcosmoinc],
    install: false,
    gnu_symbol_visibility: 'hidden',
    c_args: misc_c_args,
)
