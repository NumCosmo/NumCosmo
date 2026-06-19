src_site_tutorials_files = files(
    'index.qmd',
)

if get_option('documentation')
    foreach file : src_site_tutorials_files
        site_files += fs.copyfile(file)
    endforeach
endif

subdir('c')
subdir('python')
