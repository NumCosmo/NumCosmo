src_site_tutorials_c_files = files(
    'curve_tutorial.qmd',
    'gobject.qmd',
    'nc_curve.c',
    'nc_curve.h',
    'nc_object_name.c',
    'nc_object_name.h',
)

if get_option('documentation')
    foreach file : src_site_tutorials_c_files
        site_files += fs.copyfile(file)
    endforeach
endif