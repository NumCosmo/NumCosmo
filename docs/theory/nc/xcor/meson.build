src_site_theory_nc_xcor_files = files(
    'ssc.qmd',
)

if get_option('documentation')
    foreach file : src_site_theory_nc_xcor_files
        site_files += fs.copyfile(file)
    endforeach
endif
