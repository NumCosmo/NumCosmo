src_site_benchmark_files = files(
    '_setup_models.qmd',
    'ccl_background.qmd',
    'ccl_power_spectrum.qmd',
)

if get_option('documentation')
    foreach file : src_site_benchmark_files
        site_files += fs.copyfile(file)
    endforeach
endif