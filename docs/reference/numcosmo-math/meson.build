expand_content_md_files = []

if get_option('documentation')
    numcosmo_math_toml = configure_file(
        input: 'numcosmo-math.toml.in',
        output: 'numcosmo-math.toml',
        configuration: toml_conf,
        install: true,
        install_dir: docs_dir / 'numcosmo-math',
    )

    custom_target(
        'numcosmo-math-doc',
        input: [numcosmo_math_toml, numcosmo_math_gir[0]],
        output: 'numcosmo-math',
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
        'doc-check-numcosmo-math',
        gidocgen,
        args: [
            'check',
            '--config', numcosmo_math_toml,
            '--add-include-path=@0@'.format(meson.current_build_dir() / '../../../gtk'),
            numcosmo_math_gir[0],
        ],
        depends: numcosmo_math_gir[0],
        suite: ['docs'],
    )
endif