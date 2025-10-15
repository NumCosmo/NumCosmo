src_site_notes_files = files('_adiabatic.qmd', '_numcosmo_py.qmd')

if get_option('documentation')
    foreach file : src_site_notes_files
        site_files += fs.copyfile(file)
    endforeach
endif