src_site_examples_intro_files = files(
    'hiprim.qmd',
    'hiprim_tensor.qmd',
    'recombination.qmd',
    'simple.qmd',
)

if get_option('documentation')
    foreach file : src_site_examples_intro_files
        site_files += fs.copyfile(file)
    endforeach
endif