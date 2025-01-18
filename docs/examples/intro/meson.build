src_site_examples_intro_files = files(
    'simple.qmd',
)

if get_option('documentation')
    foreach file : src_site_examples_intro_files
        site_files += fs.copyfile(file)
    endforeach
endif