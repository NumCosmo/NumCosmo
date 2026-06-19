src_site_theory_files = files(
    'csq1d.qmd',
    'distances.qmd',
)

if get_option('documentation')
    foreach file : src_site_theory_files
        site_files += fs.copyfile(file)
    endforeach
endif