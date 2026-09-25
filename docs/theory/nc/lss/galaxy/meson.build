src_site_theory_nc_lss_galaxy_files = files(
    'galaxy_wl_framework.qmd',
    'wl_ellipticity.qmd',
    'wl_shape_marginalization_fixed_quad.qmd',
    'wl_shape_marginalization_series.qmd',
)

if get_option('documentation')
    foreach file : src_site_theory_nc_lss_galaxy_files
        site_files += fs.copyfile(file)
    endforeach
endif
