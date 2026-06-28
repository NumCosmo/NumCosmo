src_site_examples_files = files(
    'index.qmd',
)

if get_option('documentation')
    foreach file : src_site_examples_files
        site_files += fs.copyfile(file)
    endforeach
endif

subdir('data_analysis')
subdir('intro')