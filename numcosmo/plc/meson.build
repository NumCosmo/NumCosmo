

clik_sources = files(
    'cldf/cldf.c',
    'cldf/cldf.h',
    'minipmc/distribution.c',
    'minipmc/errorlist.c',
    'minipmc/errorlist.h',
    'minipmc/io.c',
    'minipmc/io.h',
    'minipmc/maths_base.h',
    'minipmc/pmc.h',
    'clik.c',
    'clik.h',
    'clik_dic.c',
    'clik_dic.h',
    'clik_helper.c',
    'clik_helper.h',
    'clik_lensing.h',
    'lapack_clik.h',
    'lklbs.c',
    'lklbs.h',
    'lowly_common.c',
    'lowly_common.h',
)

clik_act_spt_lkl_sources = files(
    'actspt/ACT_equa_likelihood.F90',
    'actspt/ACT_south_likelihood.F90',
    'actspt/Foregrounds_loading.F90',
    'actspt/Highell_likelihood.F90',
    'actspt/Highell_options.F90',
    'actspt/Highell_subroutines.F90',
    'actspt/SPT_keisler_likelihood.F90',
    'actspt/SPT_reichardt_likelihood.F90',
    'actspt/clik_actspt.c',
    'actspt/clik_actspt_f90.f90',
)

clik_camspec_lkl_sources = files(
    'camspec/CAMspec.f90',
    'camspec/clik_CAMspec.c',
    'camspec/clik_CAMspec_f90.F90',
    'camspec/temp_like.f90',
    'camspec/temp_like_v3.f90',
)

clik_gibbs_lkl_sources = files(
    'gibbs/clik_gibbs.c',
    'gibbs/clik_gibbs_f90.f90',
    'gibbs/comm_br_mod.f90',
    'gibbs/comm_gauss_br_mod.f90',
    'gibbs/comm_gauss_br_mod_v3.f90',
    'gibbs/comm_lowl_mod_dist.f90',
)

clik_lowlike_lkl_sources = files(
    'lowlike/Planck_likelihood.F90',
    'lowlike/Planck_options.F90',
    'lowlike/Planck_teeebb_pixlike.F90',
    'lowlike/br_mod_dist.f90',
    'lowlike/clik_lowlike.c',
    'lowlike/clik_lowlike_f90.f90',
    'lowlike/healpix_types.f90',
    'lowlike/read_archive_map.f90',
    'lowlike/read_fits.f90',
)

clik_bflike_lkl_sources = files(
    'bflike/bflike.f90',
    'bflike/bflike_QUonly.f90',
    'bflike/bflike_smw.f90',
    'bflike/clik_bflike.c',
    'bflike/clik_bflike_f90.f90',
    'bflike/fitstools_smw.f90',
    'bflike/long_intrinsic_smw.f90',
)

clik_pliklite_lkl_sources = files(
    'cmbonly/clik_cmbonly.c',
    'cmbonly/clik_cmbonly_f90.f90',
    'cmbonly/plik_cmbonly.f90',
)

clik_plik_lkl_sources = files(
    'plik/clik_hfipack.c',
    'plik/clik_parametric.c',
    'plik/clik_parametric.h',
    'plik/clik_parametric_addon.c',
    'plik/clik_parametric_addon.h',
    'plik/component_plugin/rel2015/corrnoise.c',
    'plik/component_plugin/rel2015/fg2015.c',
    'plik/component_plugin/rel2015/leakage.c',
    'plik/smica.c',
    'plik/smica.h',
)

clik_lens_lkl_sources = files(
    'lenslike/plenslike/plenslike.h',
    'lenslike/plenslike/plenslike_dat_full.c',
    'lenslike/plenslike/plenslike_dat_mono.c',
    'lenslike/plenslike/plenslike_dat_qecl.c',
    'lenslike/plenslike/plenslike_dat_quad.c',
    'lenslike/plenslike/qest.c',
    'lenslike/plenslike/wignerd.c',
    'clik_lensing.c',
)

clik_simall_sources = files('simall/clik_simall.c')


plc_inc = include_directories(['.', 'cldf', 'minipmc', 'plik'])

plc_preprocessor_args = [
    '-DG_LOG_DOMAIN="NUMCOSMO"',
    '-DHAVE_CONFIG_H',
    '-DHAS_LAPACK',
    '-DLAPACK_CLIK',
    '-DNOHEALPIX',
    '-DCLIK_LENSING',
    '-DCAMSPEC_V1',
    '-DCLIKSVNVERSION="0031e3604f0c"',
]

plc_c_args = [] + plc_preprocessor_args

if numcosmo_conf.has('HAVE_RTLD_DEFAULT')
    plc_c_args += '-DHAS_RTLD_DEFAULT'
endif

disable_warnings = ['-Wno-implicit-fallthrough', '-Wno-sign-compare', '-Wno-maybe-uninitialized']

foreach flag : disable_warnings
    if cc.has_argument(flag)
        plc_c_args += flag
    endif
endforeach

plc_fc_args = [] + plc_preprocessor_args
freeform_args = [
    '-ffree-form',
    '-FR',
    '-free',
    '-qfree',
    '-Mfree',
    '-Mfreeform',
    '-freeform',
    '-f free',
    '-8',
    '+source=free',
    '-nfix',
    '--nfix',
    '-Free',
]

free_form_code = '''

  program freeform
       ! FIXME: how to best confuse non-freeform compilers?
       print *, 'Hello ', &
          'world.'
       end
'''

if (
    not (
        fc.compiles(
            free_form_code,
            name: 'freeform',
        )
    )
)
    foreach flag : freeform_args
        if fc.has_argument(flag)
            if fc.compiles(
                free_form_code,
                args: flag,
                name: 'freeform',
            )
                plc_fc_args += flag
                break
            endif
        endif
    endforeach
endif

line_length_args = [
    '-ffree-line-length-none',
    '-ffixed-line-length-none',
    '-ffree-line-length-huge',
    '-ffree-line-length-unlimited',
    '-ffixed-line-length-unlimited',
    '-qfixed=unlimited',
    '-Mextend',
    '-unlimited',
    '-extend_source',
    '-Wunlimited',
    '-W',
    '+extend_source',
    '+es',
    '-wide',
    '--wide',
    '-w',
    '-e',
    '-f',
    '-Wf,-f',
    '-xline',
]

line_length_code = '''
subroutine longer_than_132(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20)
end subroutine
'''

if not fc.compiles(line_length_code, name: 'line_length')
    foreach flag : line_length_args
        if fc.has_argument(flag)
            if fc.compiles(
                line_length_code,
                args: flag,
                name: 'line_length',
            )
                plc_fc_args += flag
                break
            endif
        endif
    endforeach
endif

disable_warnings_fortran = [
    '-Wno-tabs',
    '-Wno-unused-variable',
    '-Wno-unused-dummy-argument',
    '-Wno-unused-function',
    '-Wno-conversion',
    '-Wno-compare-reals',
    '-Wno-do-subscript',
    '-Wno-maybe-uninitialized',
    '-fautomatic',
    '-frecursive',
]

foreach flag : disable_warnings_fortran
    if fc.has_argument(flag)
        plc_fc_args += flag
    endif
endforeach


libplc = static_library(
    'plc',
    sources: [
        clik_sources,
        clik_act_spt_lkl_sources,
        clik_camspec_lkl_sources,
        clik_gibbs_lkl_sources,
        clik_lowlike_lkl_sources,
        clik_bflike_lkl_sources,
        clik_pliklite_lkl_sources,
        clik_plik_lkl_sources,
        clik_lens_lkl_sources,
        clik_simall_sources,
    ],
    dependencies: [omp_c],
    include_directories: [configinc, numcosmoinc, plc_inc],
    install: false,
    gnu_symbol_visibility: 'hidden',
    c_args: plc_c_args,
    fortran_args: plc_fc_args,
)
