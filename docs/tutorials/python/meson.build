src_site_tutorials_python_files = files(
    'bounce.qmd',
    'bounce_spectra.qmd',
    'cluster_wl_simul.qmd',
    'galaxy_wl_mock.qmd',
    'halo_catalog_mock.qmd',
    'sky_match.qmd',
    'sky_match_id.qmd',
)

if get_option('documentation')
    foreach file : src_site_tutorials_python_files
        site_files += fs.copyfile(file)
    endforeach
endif