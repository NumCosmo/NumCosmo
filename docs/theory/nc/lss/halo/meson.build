src_site_theory_nc_lss_halo_files = files(
    'castro_hmf.qmd',
    'halo_density_profile.qmd',
)

if get_option('documentation')
    foreach file : src_site_theory_nc_lss_halo_files
        site_files += fs.copyfile(file)
    endforeach
endif
