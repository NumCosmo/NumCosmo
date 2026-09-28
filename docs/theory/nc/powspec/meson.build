src_site_theory_nc_powspec_files = files(
    'transfer_func_eh.qmd',
)

if get_option('documentation')
    foreach file : src_site_theory_nc_powspec_files
        site_files += fs.copyfile(file)
    endforeach
endif
