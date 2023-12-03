

executable(
    'cmb_maps',
    'cmb_maps.c',
    install: true,
    include_directories: configinc,
    dependencies: [libnumcosmo_dep],
)

executable(
    'gobj_itest',
    'gobj_itest.c',
    install: true,
    include_directories: configinc,
    dependencies: [libnumcosmo_dep],
)

executable(
    'mcat_analyze',
    'mcat_analyze.c',
    install: true,
    include_directories: configinc,
    dependencies: [libnumcosmo_dep],
)

executable(
    'mcat_join',
    'mcat_join.c',
    install: true,
    include_directories: configinc,
    dependencies: [libnumcosmo_dep],
)

executable(
    'mcat_trim_oob',
    'mcat_trim_oob.c',
    install: true,
    include_directories: configinc,
    dependencies: [libnumcosmo_dep],
)

executable(
    'mset_gen',
    'mset_gen.c',
    install: true,
    include_directories: configinc,
    dependencies: [libnumcosmo_dep],
)

python_scripts = [
    'mcat_calibrate_apes',
    'mcat_plot_corner',
]

foreach script : python_scripts
    install_data(
        script,
        install_dir: numcosmo_bindir,
        install_mode: 'rwxr-xr-x',
    )
endforeach

install_data(
    'mcat_plot_corner',
    install_dir: numcosmo_bindir,
    install_mode: 'rwxr-xr-x',
)

executable(
    'darkenergy',
    [
        'darkenergy.c',
        'data_cluster.c',
        'data_cluster.h',
        'de_options.c',
        'de_options.h',
        'savedata.c',
        'savedata.h',
    ],
    install: true,
    include_directories: configinc,
    dependencies: [libnumcosmo_dep],
)
