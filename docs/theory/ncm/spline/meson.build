src_site_theory_ncm_spline_files = files(
    'spline_func.qmd',
)

if get_option('documentation')
    foreach file : src_site_theory_ncm_spline_files
        site_files += fs.copyfile(file)
    endforeach
endif
