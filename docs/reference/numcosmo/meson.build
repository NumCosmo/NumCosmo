expand_content_md_files = []

if get_option('documentation')
    numcosmo_toml = configure_file(
        input: 'numcosmo.toml.in',
        output: 'numcosmo.toml',
        configuration: toml_conf,
        install: true,
        install_dir: docs_dir / 'numcosmo',
    )

    custom_target(
        'numcosmo-doc',
        input: [numcosmo_toml, numcosmo_gir[0]],
        output: 'numcosmo',
        command: [
            gidocgen,
            'generate',
            gidocgen_common_args,
            '--add-include-path=@0@'.format(meson.current_build_dir() / '../../../numcosmo'),
            '--config=@INPUT0@',
            '--output-dir=@OUTPUT@',
            '--content-dir=@0@'.format(meson.current_source_dir()),
            '@INPUT1@',
        ],
        depend_files: [expand_content_md_files],
        build_by_default: true,
        install: true,
        install_dir: docs_dir,
    )

    test(
        'doc-check-gdk',
        gidocgen,
        args: [
            'check',
            '--config', numcosmo_toml,
            '--add-include-path=@0@'.format(meson.current_build_dir() / '../../../gtk'),
            numcosmo_gir[0],
        ],
        depends: numcosmo_gir[0],
        suite: ['docs'],
    )
endif