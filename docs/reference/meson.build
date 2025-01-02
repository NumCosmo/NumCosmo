toml_conf = configuration_data()
toml_conf.set('version', meson.project_version())

gidocgen = find_program('gi-docgen', required: get_option('documentation'), native: true)
gidocgen_common_args = ['--quiet', '--no-namespace-dir', '--fatal-warnings']

docs_dir = numcosmo_datadir / 'doc'

if get_option('documentation') and not enable_gir
  error('API reference requires introspection.')
endif

numcosmo_math_expand_content_md_files = ['overview.md']
numcosmo_expand_content_md_files = ['overview.md']

if get_option('documentation')
  numcosmo_math_toml = configure_file(
    input: 'numcosmo-math.toml.in',
    output: 'numcosmo-math.toml',
    configuration: toml_conf,
    install: true,
    install_dir: docs_dir / 'numcosmo-math',
  )
  numcosmo_toml = configure_file(
    input: 'numcosmo.toml.in',
    output: 'numcosmo.toml',
    configuration: toml_conf,
    install: true,
    install_dir: docs_dir / 'numcosmo',
  )

  numcosmo_math_doc = custom_target(
    'numcosmo-math-doc',
    input: [numcosmo_math_toml, numcosmo_math_gir[0]],
    output: 'numcosmo-math',
    command: [
      gidocgen,
      'generate',
      gidocgen_common_args,
      '--add-include-path=@0@'.format(meson.current_build_dir() / '../../numcosmo'),
      '--add-include-path=@0@'.format(girdir),
      '--config=@INPUT0@',
      '--output-dir=@OUTPUT@',
      '--content-dir=@0@'.format(meson.current_source_dir()),
      '@INPUT1@',
    ],
    depend_files: [numcosmo_math_expand_content_md_files],
    build_by_default: true,
    install: true,
    install_dir: docs_dir,
  )
  numcosmo_doc = custom_target(
    'numcosmo-doc',
    input: [numcosmo_toml, numcosmo_gir[0]],
    output: 'numcosmo',
    command: [
      gidocgen,
      'generate',
      gidocgen_common_args,
      '--add-include-path=@0@'.format(meson.current_build_dir() / '../../numcosmo'),
      '--add-include-path=@0@'.format(girdir),
      '--config=@INPUT0@',
      '--output-dir=@OUTPUT@',
      '--content-dir=@0@'.format(meson.current_source_dir()),
      '@INPUT1@',
    ],
    depend_files: [numcosmo_expand_content_md_files],
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
      '--add-include-path=@0@'.format(meson.current_build_dir() / '../../numcosmo'),
      numcosmo_math_gir[0],
    ],
    depends: numcosmo_math_gir[0],
    suite: ['docs'],
  )
  test(
    'doc-check-gdk',
    gidocgen,
    args: [
      'check',
      '--config', numcosmo_toml,
      '--add-include-path=@0@'.format(meson.current_build_dir() / '../../numcosmo'),
      numcosmo_gir[0],
    ],
    depends: numcosmo_gir[0],
    suite: ['docs'],
  )

endif