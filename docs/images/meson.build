src_site_images_files = files(
    'favicon_numcosmo.png',
    'numcosmo.png',
    'numcosmo_logo.png',
    'spline_func_knots_evolution.png',
)

if get_option('documentation')
    foreach file : src_site_images_files
        site_files += fs.copyfile(file)
    endforeach
endif