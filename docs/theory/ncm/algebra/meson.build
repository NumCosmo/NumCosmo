src_site_theory_ncm_algebra_files = files(
    'spectral.qmd',
)

if get_option('documentation')
    foreach file : src_site_theory_ncm_algebra_files
        site_files += fs.copyfile(file)
    endforeach
endif
