src_site_theory_ncm_specfunc_files = files(
    'sbessel_ode_solver.qmd',
    'sbessel_projection.qmd',
)

if get_option('documentation')
    foreach file : src_site_theory_ncm_specfunc_files
        site_files += fs.copyfile(file)
    endforeach
endif
