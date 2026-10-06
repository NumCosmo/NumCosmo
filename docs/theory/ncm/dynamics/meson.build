src_site_theory_ncm_dynamics_files = files(
    'csq1d.qmd',
)

if get_option('documentation')
    foreach file : src_site_theory_ncm_dynamics_files
        site_files += fs.copyfile(file)
    endforeach
endif
