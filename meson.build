project(
  'numcosmo',
  ['c', 'fortran'],
  version: '0.21.2',
  license: 'GPL-3.0-or-later',
  meson_version: '>= 1.2.0',
  default_options: [
    'buildtype=debugoptimized',
    'warning_level=2',
    'c_std=gnu99',
    'b_lto=true',
    'b_lto_threads=2',
  ],
)

cc = meson.get_compiler('c')
c_standards = {}

foreach std : ['90', '99', '11', '17']
  arg = (cc.get_id() == 'msvc' ? '/std:' : '-std=') + 'c' + std
  if cc.has_argument(arg)
    c_standards += {std: arg}
  endif
endforeach

fc = meson.get_compiler('fortran')
fortran_standards = {}

foreach std : ['legacy', 'f77', 'f90', 'f95', 'f2003', 'f2008', 'f2018']
  arg = (fc.get_id() == 'msvc' ? '/std:' : '-std=') + std
  if fc.has_argument(arg)
    fortran_standards += {std: arg}
  endif
endforeach

global_warn_suppress = ['-Wno-unused', '-Wno-pedantic', '-Wno-unused-parameter']

foreach warn : global_warn_suppress
  if cc.has_argument(warn)
    add_global_arguments(warn, language: 'c')
  endif
endforeach

cc_can_run = meson.can_run_host_binaries()

if cc.get_id() == 'msvc'
  # Set the input and exec encoding to utf-8, like is the default with GCC
  add_project_arguments(cc.get_supported_arguments(['/utf-8']), language: 'c')
  # Disable SAFESEH with MSVC for plugins and libs that use external deps that
  # are built with MinGW
  noseh_link_args = ['/SAFESEH:NO']
else
  noseh_link_args = []
  # -mms-bitfields vs -fnative-struct ?
endif

host_system = host_machine.system()

numcosmo_version = meson.project_version()
numcosmo_api_version = '1.0'
version_arr = numcosmo_version.split('.')
major_version = version_arr[0].to_int()
minor_version = version_arr[1].to_int()
micro_version = version_arr[2].to_int()
numcosmo_fullname = '@0@-@1@'.format('numcosmo', numcosmo_version)

interface_age = minor_version.is_odd() ? 0 : micro_version
binary_age = 100 * minor_version + micro_version

soversion = 0
# Maintain compatibility with previous libtool versioning
# current = minor * 100 + micro
current = binary_age - interface_age
library_version = '@0@.@1@.@2@'.format(soversion, current, interface_age)
darwin_versions = [current + 1, '@0@.@1@'.format(current + 1, interface_age)]

configinc = include_directories('.')
numcosmoinc = include_directories('numcosmo')

numcosmo_prefix = get_option('prefix')
numcosmo_bindir = join_paths(numcosmo_prefix, get_option('bindir'))
numcosmo_libdir = join_paths(numcosmo_prefix, get_option('libdir'))
numcosmo_libexecdir = join_paths(numcosmo_prefix, get_option('libexecdir'))
numcosmo_datadir = join_paths(numcosmo_prefix, get_option('datadir'))
numcosmo_pkgdatadir = join_paths(numcosmo_datadir, numcosmo_fullname)
numcosmo_includedir = join_paths(numcosmo_prefix, get_option('includedir'), 'numcosmo')

numcosmo_pkgconfigreldir = join_paths(numcosmo_libdir, 'pkgconfig')

# When building numcosmo and gobject-introspection with subprojects, gobject-introspection
# requires to know the path of the sources and the build directory for the subproject.
# We provide it here with a variable.
numcosmo_source_dir = meson.current_source_dir()
numcosmo_build_dir = meson.current_build_dir()

# Don’t build the tests unless we can run them (either natively, in an exe wrapper, or by installing them for later use)
build_tests = get_option('tests') and meson.can_run_host_binaries()

common_test_env = [
  'G_DEBUG=gc-friendly',
  'G_ENABLE_DIAGNOSTIC=1',
  'MALLOC_CHECK_=2',
]

# Note: this may cause the tests output not to be printed when running in
# verbose mode, see https://github.com/mesonbuild/meson/issues/11185
# Can be changed it to 'exitcode' if required during development.
test_protocol = 'tap'
test_timeout = 30
test_timeout_slow = 90

add_test_setup(
  'default',
  is_default: true,
  exclude_suites: ['flaky', 'failing'],
  env: common_test_env,
  timeout_multiplier: 2,
)

# dummy/empty dependency() object to declare fallbacks and simpler dependencies
not_found = dependency('', required: false)

#######################################################################################
# Configuration begins:
#
# numcosmo_conf is used to generate config.h and its used by numcosmo.
#
# numcosmo_build_cfg is used to generate numcosmo/build_cfg.h and its used by codes
# that use numcosmo.
#######################################################################################
numcosmo_conf = configuration_data()
numcosmo_build_cfg = configuration_data()

# accumulated list of defines as we check for them, so we can easily
# use them later in test programs (autoconf does this automatically)
numcosmo_conf_prefix = ''

numcosmo_conf.set('NUMCOSMO_MAJOR_VERSION', major_version)
numcosmo_conf.set('NUMCOSMO_MINOR_VERSION', minor_version)
numcosmo_conf.set('NUMCOSMO_MICRO_VERSION', micro_version)
numcosmo_conf.set('NUMCOSMO_INTERFACE_AGE', interface_age)
numcosmo_conf.set('NUMCOSMO_BINARY_AGE', binary_age)
numcosmo_conf.set_quoted('PACKAGE_BUGREPORT', 'https://github.com/NumCosmo/NumCosmo/issues')
numcosmo_conf.set_quoted('PACKAGE_NAME', 'numcosmo')
numcosmo_conf.set_quoted('PACKAGE_STRING', 'numcosmo @0@'.format(meson.project_version()))
numcosmo_conf.set_quoted('PACKAGE_TARNAME', 'numcosmo')
numcosmo_conf.set_quoted('PACKAGE_URL', 'https://github.com/NumCosmo/NumCosmo')
numcosmo_conf.set_quoted('PACKAGE_VERSION', meson.project_version())
numcosmo_conf.set_quoted('PACKAGE_DATA_DIR', numcosmo_pkgdatadir)
numcosmo_conf.set_quoted('PACKAGE_SOURCE_DIR', numcosmo_source_dir)

numcosmo_build_cfg.set('NUMCOSMO_MAJOR_VERSION', major_version)
numcosmo_build_cfg.set('NUMCOSMO_MINOR_VERSION', minor_version)
numcosmo_build_cfg.set('NUMCOSMO_MICRO_VERSION', micro_version)
numcosmo_build_cfg.set('NUMCOSMO_VERSION', numcosmo_version)

numcosmo_build_cfg.set('have_inline_support', '#define NUMCOSMO_HAVE_INLINE 1')
numcosmo_build_cfg.set('have_inline_macro', '#define NCM_INLINE static inline')
numcosmo_build_cfg.set('numcosmo_source_dir', numcosmo_source_dir)
numcosmo_build_cfg.set('numcosmo_build_dir', numcosmo_build_dir)
numcosmo_conf.set('HAVE_INLINE', 1)

# used by the .rc.in files
numcosmo_build_cfg.set('LT_CURRENT_MINUS_AGE', soversion)

numcosmo_build_shared = false
numcosmo_build_static = false
if get_option('default_library') == 'both'
  numcosmo_build_static = true
  numcosmo_build_shared = true
elif get_option('default_library') == 'static'
  numcosmo_build_static = true
elif get_option('default_library') == 'shared'
  numcosmo_build_shared = true
endif

numcosmo_build_both = numcosmo_build_static and numcosmo_build_shared
numcosmo_build_static_only = numcosmo_build_static and not numcosmo_build_shared
numcosmo_build_shared_only = numcosmo_build_shared and not numcosmo_build_static

if (
  numcosmo_build_shared
  and numcosmo_build_static
  and (host_system == 'windows'
  or host_system == 'cygwin')
)
  error('On Windows default_library must be "shared" or "static" but not "both"')
endif

if numcosmo_build_static_only
  numcosmo_build_cfg.set('NUMCOSMO_STATIC_COMPILATION', '1')
  numcosmo_build_cfg.set('GOBJECT_STATIC_COMPILATION', '1')
  numcosmo_build_cfg.set('GIO_STATIC_COMPILATION', '1')
  numcosmo_build_cfg.set('GMODULE_STATIC_COMPILATION', '1')
  numcosmo_build_cfg.set('G_INTL_STATIC_COMPILATION', '1')
  numcosmo_build_cfg.set('FFI_STATIC_BUILD', '1')
endif

# Use debug/optimization flags to determine whether to enable debug or disable
# cast checks. We have a non-production (debug) build if debug is true and if
# optimization is 0 or g; otherwise, we have a production build.
numcosmo_debug_cflags = []
numcosmo_debug = get_option('numcosmo_debug')
if (
  numcosmo_debug.enabled()
  or (numcosmo_debug.auto()
  and get_option('debug')
  and get_option('optimization') in ['0', 'g'])
)
  numcosmo_debug_cflags += ['-DG_ENABLE_DEBUG']
  message('Enabling various debug infrastructure')
else
  numcosmo_debug_cflags += ['-DG_DISABLE_CAST_CHECKS']
  message('Disabling cast checks')
endif

if not get_option('numcosmo_assert')
  numcosmo_debug_cflags += ['-DG_DISABLE_ASSERT']
  message('Disabling GLib asserts')
endif

if not get_option('numcosmo_checks')
  numcosmo_debug_cflags += ['-DG_DISABLE_CHECKS']
  message('Disabling GLib checks')
endif

gtkdoc_dep = dependency('gtk-doc', required: get_option('gtkdoc'))
build_gtk_doc = gtkdoc_dep.found()

add_project_arguments(numcosmo_debug_cflags, language: 'c')

#######################################################################################
# Checking headers:
#######################################################################################
message('Checking for headers availability')

headers = [
  'alloca.h',
  'dlfcn.h',
  'execinfo.h',
  'float.h',
  'math.h',
  'stddef.h',
  'stdint.h',
  'stdio.h',
  'stdlib.h',
]

foreach h : headers
  if cc.has_header(h)
    define = 'HAVE_' + h.underscorify().to_upper()
    numcosmo_conf.set(define, 1)
  endif
endforeach

if numcosmo_conf.has('HAVE_DLFCN_H')
  if cc.has_header_symbol('dlfcn.h', 'RTLD_DEFAULT', required: false)
    numcosmo_conf.set('HAVE_RTLD_DEFAULT', 1)
  endif
endif

#######################################################################################
# Checking functions:
#######################################################################################
message('Checking for functions availability')

functions = [
  'alloca',
  'backtrace_symbols',
  'backtrace',
  'cos',
  'erf',
  'erf',
  'erfc',
  'exp10',
  'finite',
  'fma',
  'fork',
  'getloadavg',
  'isfinite',
  'lgamma_r',
  'lgammaf_r',
  'powl',
  'powl',
  'shmget',
  'sin',
  'sincos',
]

foreach f : functions
  if cc.has_function(f)
    define = 'HAVE_' + f.underscorify().to_upper()
    numcosmo_conf.set(define, 1)
  endif
endforeach

if cc.has_header_symbol('math.h', 'isfinite', required: false)
  numcosmo_conf.set('HAVE_DECL_ISFINITE', 1)
endif

#######################################################################################
# Compiler features:
#######################################################################################
message('Checking for compiler features')

if cc.compiles(
  '''static int x;
                  char some_func (void) {
                    char a[++x];
                    a[sizeof a - 1] = 0;
                    int N;
                    return a[0];
                  }''',
  name: 'C99 variable-size arrays',
)
  cc_supports_c99_vla = true
endif

#######################################################################################
# Types:
#######################################################################################
message('Checking for types availability')

if cc.compiles(
  '''/* The Stardent Vistra knows sizeof(long double), but does not support it.  */
                  long double foo = 0.0;
                  /* On Ultrix 4.3 cc, long double is 4 and double is 8.  */
                  int array [2*(sizeof(long double) >= sizeof(double)) - 1];''',
  name: 'long double',
)
  numcosmo_conf.set('HAVE_LONG_DOUBLE', 1)
endif

#######################################################################################
# Fortran mangling:
#######################################################################################

numcosmo_conf.set('F77_FUNC(name,NAME)', 'name ## _')

#######################################################################################
# Checking for basic libraries:
#######################################################################################

dl_dep = dependency('dl', method: 'auto')

#######################################################################################
# Dependencies versions:
#######################################################################################

glib_req = '>= 2.54.0'
gsl_req = '>= 2.4'
fftw_req = '>= 3.1.2'
cfitsio_req = '>= 3.25'
flint_req = '>= 3.0'

#######################################################################################
# Checking for MPI:
#
# MPI is not required, but if it is found, we will use it.
#######################################################################################
message('Checking for MPI')

mpi_c_dep = dependency(
  'mpi',
  language: 'c',
  required: false,
)

mpi_fortran_dep = dependency(
  'mpi',
  language: 'fortran',
  required: false,
)

# MPICH is not found by the mpi dependency, so we need to check for it manually
if not mpi_c_dep.found()
  mpi_c_dep = dependency('mpich', required: false)
endif

# MPICH is not found by the mpi dependency, so we need to check for it manually
if not mpi_fortran_dep.found()
  mpi_fortran_dep = dependency('mpichf90', required: false)
endif

if mpi_c_dep.found()
  numcosmo_conf.set('HAVE_MPI', true)
  numcosmo_build_cfg.set('have_mpi_support', '#define NUMCOSMO_HAVE_MPI 1')
  add_global_arguments(
    '-DUSE_NCM_MPI',
    language: 'c',
  )
else
  numcosmo_build_cfg.set('have_mpi_support', '')
endif

if get_option('debug-mpi') and mpi_c_dep.found()
  add_global_arguments(
    '-DNCM_MPI_DEBUG',
    language: 'c',
  )
endif

#######################################################################################
# Checking for OpenMP:
#
# OpenMP is not required, but if it is found, we will use it.
#######################################################################################
message('Checking for OpenMP')

omp_c = dependency(
  'openmp',
  language: 'c',
  required: false,
)
omp_fortran = dependency(
  'openmp',
  language: 'fortran',
  required: false,
)

#######################################################################################
# Checking for GObject instrospection:
#
# GObject introspection is not required, but if it is found, we will use it.
#
# Also check for Vala, which is required for the Vala bindings.
#######################################################################################
message('Checking for GObject introspection')

gir_scanner = find_program('g-ir-scanner', required: get_option('introspection'))

enable_gir = false
if gir_scanner.found()
  if get_option('introspection').auto() and meson.can_run_host_binaries()
    enable_gir = true
  elif get_option('introspection').enabled()
    enable_gir = true
  endif
endif

if get_option('introspection').enabled() and not meson.can_run_host_binaries()
  message(
    'To build GObject introspection files in a cross build, you need to set the appropriate GI_CROSS_LAUNCHER variable.',
  )
endif

gir_args = [
  '--quiet',
]

have_vala = add_languages(
  'vala',
  required: get_option('vala'),
  native: false,
)

#######################################################################################
# Checking for Glib:
#######################################################################################
message('Checking for Glib, GObject and GIO')

glib_dep = dependency('glib-2.0', version: glib_req)
gobject_dep = dependency('gobject-2.0', version: glib_req)
gio_dep = dependency('gio-2.0', version: glib_req)

#######################################################################################
# Checking for GSL:
#######################################################################################
message('Checking for GSL')

gsl_dep = dependency('gsl', version: gsl_req)

if not get_option('gsl_range_check')
  numcosmo_conf.set('GSL_RANGE_CHECK_OFF', true)
endif

#######################################################################################
# Checking for GMP and MPFR:
#######################################################################################
message('Checking for GMP and MPFR')

gmp_dep = dependency('gmp')
mpfr_dep = dependency('mpfr')

#######################################################################################
# Checking for FFTW (optional):
#######################################################################################
message('Checking for FFTW')

fftw3_dep = dependency(
  'fftw3',
  version: fftw_req,
  required: false,
)

fftw3f_dep = dependency(
  'fftw3f',
  version: fftw_req,
  required: false,
)

if fftw3_dep.found()
  numcosmo_conf.set('HAVE_FFTW3', true)
endif

if fftw3f_dep.found()
  numcosmo_conf.set('HAVE_FFTW3F', true)
endif

#######################################################################################
# Checking for cfitsio (optional):
#######################################################################################
message('Checking for cfitsio')

cfitsio_dep = dependency(
  'cfitsio',
  version: cfitsio_req,
  required: false,
)

if cfitsio_dep.found()
  numcosmo_conf.set('HAVE_CFITSIO', true)
endif

#######################################################################################
# Checking for libfyaml
#######################################################################################
message('Checking for libfyaml')

libfyaml_dep = dependency(
  'libfyaml',
  required: false,
)

if libfyaml_dep.found()
  numcosmo_conf.set('HAVE_LIBFYAML', true)
endif

#######################################################################################
# Using internal libcuba:
#######################################################################################

numcosmo_conf.set('HAVE_LIBCUBA', 1)
numcosmo_conf.set('HAVE_LIBCUBA_4_0', 1)

#######################################################################################
# Checking for NLopt (optional):
#######################################################################################
message('Checking for NLopt')

nlopt_dep = dependency(
  'nlopt',
  required: false,
)

if nlopt_dep.found()
  numcosmo_conf.set('HAVE_NLOPT', true)
  numcosmo_build_cfg.set('have_nlopt_support', '#define NUMCOSMO_HAVE_NLOPT 1')
endif

#######################################################################################
# Checking for Flint (optional):
#######################################################################################
message('Checking for Flint')

# Currently the flint pkg-config file is broken, so we need to check for the
# library manually. When this is fixed, we can replace 'flint_does_not_work' with
# 'flint' in the dependency() call below.
flint_dep = dependency(
  'flint_does_not_work',
  version: flint_req,
  required: false,
)

if flint_dep.found()
  arb_dep = flint_dep
else
  flint_dep = cc.find_library(
    'flint',
    required: false,
  )

  if flint_dep.found()
    if cc.has_function('acb_set', dependencies: [flint_dep])
      arb_dep = flint_dep
    else
      arb_lib_names = [
        'arb',
        'flint-arb',
      ]
      foreach arb_lib_name : arb_lib_names
        arb_dep = cc.find_library(
          arb_lib_name,
          required: false,
        )
        if arb_dep.found()
          if cc.has_function('acb_set', dependencies: [arb_dep, flint_dep])
            break
          endif
        endif
      endforeach
    endif
  else
    arb_dep = flint_dep
  endif
endif

if arb_dep.found()
  numcosmo_conf.set('HAVE_FLINT', true)
  flint_headers = [
    'arb.h',
    'acb.h',
  ]
  foreach h : flint_headers
    if cc.has_header(
      h,
      dependencies: [arb_dep, flint_dep],
      required: false,
    )
      define = 'HAVE_' + h.underscorify().to_upper()
      numcosmo_conf.set(define, 1)
    elif cc.has_header(
      'flint' / h,
      dependencies: [arb_dep, flint_dep],
      required: false,
    )
      define = 'HAVE_' + h.underscorify().to_upper()
      define_flint = 'HAVE_FLINT_' + h.underscorify().to_upper()
      numcosmo_conf.set(define, 1)
      numcosmo_conf.set(define_flint, 1)
    endif
  endforeach
endif

#######################################################################################
# Checking for BLAS/LAPACK:
#######################################################################################
message('Checking for BLAS/LAPACK')

blas_name = get_option('blas')
lapack_name = get_option('lapack')

# Trying to find the blas library, if 'auto' is specified, try to find the first
# available blas library in the list of blas libraries specified in blas-order
if blas_name == 'auto'
  foreach _name : get_option('blas-order')
    blas_dep = dependency(
      _name,
      modules: ['cblas', 'interface: lp64'],
      required: false,
    )
    if blas_dep.found()
      blas_name = _name
      break
    endif
  endforeach
else
  blas_dep = dependency(
    blas_name,
    modules: ['cblas', 'interface: lp64'],
    required: false,
  )
endif

# If we still haven't found a blas library, try to find the blas library
# specified by the user
if not blas_dep.found()
  blas_dep = cc.find_library('blas', required: false)
  blas_dep_name = 'generic'
else
  blas_dep_name = blas_dep.name()
endif
# If non of the above worked, we rely on the blas from gsl

# Check if the blas library has a cblas header and if so, check if
# it defines CBLAS_ORDER, if not, we need to define it ourselves
if blas_name == 'mkl'
  if cc.has_header('mkl_cblas.h', required: true)
    numcosmo_conf.set('HAVE_MKL_CBLAS_H', true)
    cblas_header = 'mkl_cblas.h'
  endif
else
  if cc.has_header('cblas.h', required: false)
    numcosmo_conf.set('HAVE_CBLAS_H', true)
    cblas_header = 'cblas.h'
  else
    numcosmo_conf.set('HAVE_GSL_CBLAS_H', true)
    cblas_header = 'gsl/gsl_cblas.h'
  endif
endif
if not cc.has_header_symbol(cblas_header, 'CBLAS_ORDER', required: false)
  numcosmo_conf.set('BLAS_NOT_TYPEDEFED', true)
endif

# If we have either mkl or accelerate, we don't need lapack since they
# provide lapack routines
if 'mkl' in blas_dep_name or blas_dep_name == 'accelerate'
  lapack_dep = blas_dep
else
  if lapack_name == 'auto'
    foreach _name : get_option('lapack-order')
      lapack_dep = dependency(
        _name,
        modules: ['lapack', 'interface: lp64'],
        required: false,
      )
      if lapack_dep.found()
        break
      endif
    endforeach
  else
    lapack_dep = dependency(
      lapack_name,
      modules: ['lapack', 'interface: lp64'],
      required: false,
    )
  endif
endif

if not lapack_dep.found()
  lapack_dep = cc.find_library('lapack', required: false)
endif

have_lapack = lapack_dep.found()
if have_lapack
  numcosmo_conf.set('HAVE_LAPACK', true)
  lapack_functions = [
    'dgeev_',
    'dgeevx_',
    'dgelqf_',
    'dgels_',
    'dgelsd_',
    'dgelsy_',
    'dgeqlf_',
    'dgeqrf_',
    'dgerqf_',
    'dggglm_',
    'dposv_',
    'dpotrf_',
    'dpotri_',
    'dpotrs_',
    'dpstrf_',
    'dptsv_',
    'dsyevd_',
    'dsyevr_',
    'dsysv_',
    'dsysvx_',
    'dsysvxx_',
    'dsytrf_',
    'dsytri_',
    'dsytrs_',
    'dtrsv_',
  ]
  foreach f : lapack_functions
    if cc.has_function(f, dependencies: [lapack_dep, blas_dep])
      define = 'HAVE_' + f.underscorify().to_upper()
      numcosmo_conf.set(define, 1)
    endif
  endforeach
endif

#######################################################################################
# Checking for Python:
#######################################################################################

python_module = import('python')
python = python_module.find_installation('python3', modules: ['pytest', 'gi'])

python_pip = python_module.find_installation(
  'python3',
  modules: ['pip'],
  required: false,
)

# Installing the python module:
if get_option('numcosmo_py')
  if python_pip.found()
    pip_cmd = [
      python.full_path(),
      '-m',
      'pip',
      'install',
      numcosmo_source_dir,
      '-t',
      python.get_install_dir(),
      '--no-dependencies',
    ]
    meson.add_install_script(pip_cmd, skip_if_destdir: true)
  else
    error('Python pip not found, cannot install the NumCosmoPy python module')
  endif
endif

#######################################################################################
# Creating the configuration files:
#######################################################################################

configure_file(
  input: 'config.h.meson',
  output: 'config.h',
  configuration: numcosmo_conf,
)

configure_file(
  output: 'config_deps.h',
  configuration: numcosmo_conf,
)

configure_file(
  input: 'numcosmo_export.sh.in',
  output: 'numcosmo_export.sh',
  configuration: numcosmo_build_cfg,
)

#######################################################################################
# Installing data and examples files:
#######################################################################################

install_subdir('data', install_dir: numcosmo_pkgdatadir)
install_subdir('examples', install_dir: numcosmo_pkgdatadir)

#######################################################################################
# Including the subdirectories:
#######################################################################################

pkg = import('pkgconfig')

subdir('numcosmo')
subdir('tools')
subdir('tests')
subdir('docs')

building_options = {
  'Build type': get_option('buildtype'),
  'Optimization level': get_option('optimization'),
  'Debug': get_option('debug'),
  'C-Compiler': cc.get_id(),
  'C-Compiler version': cc.version(),
  'Fortran-Compiler': fc.get_id(),
  'Fortran-Compiler version': fc.version(),
  'Build tests': build_tests,
  'Build documentation': build_gtk_doc,
  'Build Python module': get_option('numcosmo_py'),
  'Build Vala library': have_vala,
  'Build GObject introspection': enable_gir,
  'Build MPI support': mpi_c_dep,
  'Build OpenMP support': omp_c,
  'Build prefix': numcosmo_prefix,
  'Build bindir': numcosmo_bindir,
  'Build libdir': numcosmo_libdir,
  'Build libexecdir': numcosmo_libexecdir,
  'Build datadir': numcosmo_datadir,
  'Build includedir': numcosmo_includedir,
  'Build pkgdatadir': numcosmo_pkgdatadir,
  'Build pkgconfigdir': numcosmo_pkgconfigreldir,
}

summary(
  building_options,
  section: 'Building options',
)

dependencies_summary = {
  'blas': blas_dep,
  'cfitsio': cfitsio_dep,
  'fftw3': fftw3_dep,
  'fftw3f': fftw3f_dep,
  'flint': flint_dep,
  'glib-gio': gio_dep,
  'glib-gobject': gobject_dep,
  'glib': glib_dep,
  'gmp': gmp_dep,
  'gsl': gsl_dep,
  'lapack': lapack_dep,
  'libfyaml': libfyaml_dep,
  'mpfr': mpfr_dep,
  'mpi': mpi_c_dep,
  'nlopt': nlopt_dep,
  'openmp': omp_c,
  'python': python,
}

summary(dependencies_summary, section: 'Dependencies')

