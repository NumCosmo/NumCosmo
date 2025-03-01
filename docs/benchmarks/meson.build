src_site_benchmark_files = files(
    'ccl_background.qmd',
)

if get_option('documentation')
    foreach file : src_site_benchmark_files
        site_files += fs.copyfile(file)
    endforeach
endif