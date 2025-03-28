src_site_examples_data_analysis_files = files(
    'snia.qmd',
)

if get_option('documentation')
    foreach file : src_site_examples_data_analysis_files
        site_files += fs.copyfile(file)
    endforeach
endif