if build_gtk_doc

    manual_head_xml = configure_file(
        input: 'manual_head.xml.in',
        output: 'manual_head.xml',
        configuration: numcosmo_build_cfg,
    )

    ignore_h_files = [
        'class',
        'class/hyrec',
        'class/include',
        'class/source',
        'class/tools',
        'clencurt.h',
        'common.c',
        'config_extra.h',
        'config.h',
        'converged.h',
        'cubature.h',
        'Faddeeva.h',
        'gsl_cqp.h',
        'gsl_rstat.h',
        'hcubature.c',
        'initial_point.h',
        'kdtree.c',
        'kdtree.h',
        'levmar',
        'libcuba',
        'libqp_gsmo.c',
        'libqp_splx.c',
        'libqp.h',
        'nc_hipert_private.h',
        'ncm_fit_nlopt_enum_meta.h',
        'ncm_flapack.h',
        'ncm_gsl_blas_types.h',
        'ncm_sphere_map_block.c',
        'ncm_stats_dist_kde_private.h',
        'ncm_stats_dist_kernel_private.h',
        'ncm_stats_dist_private.h',
        'ncm_stats_dist_vkde_private.h',
        'nnls.c',
        'pcubature.c',
        'plc',
        'rb_knn_list.c',
        'rb_knn_list.h',
        'rcm.c',
        'rcm.h',
        'rquantile.c',
        'stddecl.h',
        'sundials',
        'toeplitz',
        'vwrapper.h',
    ]

    if not nlopt_dep.found()
        ignore_h_files += ['ncm_fit_nlopt_enum.h', 'ncm_fit_nlopt.h']
    endif
    ignore_files_str = '@0@'.format(' '.join(ignore_h_files))

    docfiles = gnome.gtkdoc(
        'numcosmo',
        mode: 'sgml',
        main_sgml: 'numcosmo-docs.sgml',
        gobject_typesfile: 'numcosmo.types',
        dependencies: [libnumcosmo_dep],
        content_files: ['introduction.xml', 'refs.xml', manual_head_xml, 'references.xml'],
        expand_content_files: [manual_head_xml],
        src_dir: ['numcosmo'],
        scan_args: [
            '--rebuild-types',
            '--rebuild-sections',
            '--ignore-decorators=NCM_INLINE',
        ],
        mkdb_args: [
            '--sgml-mode',
            '--output-format=xml',
            '--name-space=nc',
            '--name-space=ncm',
            '--ignore-files',
            ignore_files_str,
        ],
        html_assets: [
            'container.css',
            'images/spline_func_knots_evolution.png',
            'images/vkde.png',
            'images/kde.png',
            'images/apes.png',
        ],
        install: true,
        ignore_headers: ignore_h_files,
        include_directories: [numcosmo_includedir],
        check: false,
    )

    python_bf4 = python_module.find_installation(
        'python3',
        modules: ['bs4'],
        required: false,
    )

    if not python_bf4.found()
        message('WARNING: python3-bs4 not found, skipping adding mathjax to html docs')
    else
        addmathjax = find_program('addmathjax.py', required: true)
        custom_target(
            'numcosmo-doc-mathjax',
            input: docfiles,
            output: 'numcosmo-mathjax-update',
            command: [addmathjax, '@OUTDIR@' / 'html'],
            build_by_default: true,
            build_always_stale: true,
            depends: docfiles,
        )
    endif

endif
