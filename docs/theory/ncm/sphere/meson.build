src_site_theory_ncm_sphere_files = files(
    'sphere_xi.qmd',
)

if get_option('documentation')
    foreach file : src_site_theory_ncm_sphere_files
        site_files += fs.copyfile(file)
    endforeach
endif
