src_site_theory_files = files(
    'index.qmd',
)

if get_option('documentation')
    foreach file : src_site_theory_files
        site_files += fs.copyfile(file)
    endforeach
endif

subdir('nc')
subdir('ncm')
