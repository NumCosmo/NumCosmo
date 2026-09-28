src_site_theory_ncm_fftlog_files = files(
    'fftlog.qmd',
)

if get_option('documentation')
    foreach file : src_site_theory_ncm_fftlog_files
        site_files += fs.copyfile(file)
    endforeach
endif
