src_site_files = files(
    '_definitions.qmd',
    '_functions.qmd',
    'about.qmd',
    'container.css',
    'custom.scss',
    'index.qmd',
    'install.qmd',
    'references.bib',
)

fs = import('fs')
site_files = []
subdir('reference')
subdir('images')
subdir('tutorials')
subdir('examples')
subdir('benchmarks')

quarto_conf = configuration_data()
quarto_conf.set('version', meson.project_version())

quarto_yml = configure_file(
    input: '_quarto.yml.in',
    output: '_quarto.yml',
    configuration: quarto_conf,
    install: false,
)

quarto = find_program('quarto', required: get_option('documentation'), native: true)

if get_option('documentation')
    foreach file : src_site_files
        site_files += fs.copyfile(file)
    endforeach

    env = environment()
    env.prepend(
        'GI_TYPELIB_PATH',
        numcosmo_build_dir / 'numcosmo',
    )

    env.prepend(
        'LD_LIBRARY_PATH',
        numcosmo_build_dir / 'numcosmo',
    )

    env.prepend(
        'DYLD_LIBRARY_PATH',
        numcosmo_build_dir / 'numcosmo',
    )

    env.prepend(
        'PYTHONPATH',
        numcosmo_source_dir,
    )

    custom_target(
        'numcosmo-site',
        input: quarto_yml,
        output: 'numcosmo-site',
        command: [
            quarto,
            'render',
            meson.current_build_dir(),
            '--output-dir', 'numcosmo-site',
            '--use-freezer',
        ],
        build_by_default: true,
        env: env,
        install: true,
        install_dir: docs_dir,
        depends: [numcosmo_math_doc, numcosmo_doc, site_files],
    )

endif