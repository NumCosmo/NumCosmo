src_site_theory_ncm_integration_files = files(
    'diff.qmd',
)

if get_option('documentation')
    foreach file : src_site_theory_ncm_integration_files
        site_files += fs.copyfile(file)
    endforeach
endif
