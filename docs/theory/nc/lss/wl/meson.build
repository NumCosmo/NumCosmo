src_site_theory_nc_lss_wl_files = files(
    'wl_surface_mass_density.qmd',
)

if get_option('documentation')
    foreach file : src_site_theory_nc_lss_wl_files
        site_files += fs.copyfile(file)
    endforeach
endif
