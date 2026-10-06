src_site_theory_nc_lss_cluster_files = files(
    'cluster_richness_projection.qmd',
)

if get_option('documentation')
    foreach file : src_site_theory_nc_lss_cluster_files
        site_files += fs.copyfile(file)
    endforeach
endif
