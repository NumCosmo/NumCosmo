src_site_theory_files = files(
    'index.qmd',
    'csq1d.qmd',
    'distances.qmd',
    'recombination.qmd',
    'transfer_func_eh.qmd',
    'wl_surface_mass_density.qmd',
    'halo_density_profile.qmd',
    'fftlog.qmd',
)

if get_option('documentation')
    foreach file : src_site_theory_files
        site_files += fs.copyfile(file)
    endforeach
endif