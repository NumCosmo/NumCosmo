src_site_theory_ncm_core_files = files(
    'poisson_lognormal.qmd',
)

if get_option('documentation')
    foreach file : src_site_theory_ncm_core_files
        site_files += fs.copyfile(file)
    endforeach
endif
