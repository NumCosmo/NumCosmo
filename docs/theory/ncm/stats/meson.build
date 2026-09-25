src_site_theory_ncm_stats_files = files(
    'autocorrelation.qmd',
    'stats_dist.qmd',
)

if get_option('documentation')
    foreach file : src_site_theory_ncm_stats_files
        site_files += fs.copyfile(file)
    endforeach
endif
