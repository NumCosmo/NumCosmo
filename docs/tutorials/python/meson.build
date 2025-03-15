src_site_tutorials_python_files = files(
    'bouncing.qmd',
    'cluster_wl_simul.qmd',
    'sky_match.qmd',
)

if get_option('documentation')
    foreach file : src_site_tutorials_python_files
        site_files += fs.copyfile(file)
    endforeach
endif