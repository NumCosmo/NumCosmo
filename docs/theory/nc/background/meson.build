src_site_theory_nc_background_files = files(
    'distances.qmd',
    'hcosmo_bi.qmd',
)

if get_option('documentation')
    foreach file : src_site_theory_nc_background_files
        site_files += fs.copyfile(file)
    endforeach
endif
