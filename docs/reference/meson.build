toml_conf = configuration_data()
toml_conf.set('version', meson.project_version())

gidocgen = find_program('gi-docgen', required: get_option('documentation'), native: true)
gidocgen_common_args = ['--quiet', '--no-namespace-dir', '--fatal-warnings']

docs_dir = numcosmo_datadir / 'doc'

if get_option('documentation') and not enable_gir
  error('API reference requires introspection.')
endif

subdir('numcosmo-math')
subdir('numcosmo')