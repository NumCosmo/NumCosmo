
# Vendored subset of libfyaml v0.9.6 (https://github.com/pantoniou/libfyaml, MIT
# license, see LICENSE in this directory).
#
# Only the sources that make up libfyaml itself are vendored -- the upstream CLI
# tool, reflection/libclang backend, python bindings, and test harness are not
# needed by NumCosmo and are not included. Within src/blake3/, only the portable
# hashing backend is built (no SSE2/SSE41/AVX2/AVX512/NEON variants, which would
# also need matching hand-written assembly): blake3_backend.c gracefully
# registers only the portable hasher when the TARGET_HAS_* macros are undefined,
# and content-hash throughput doesn't matter for parsing config-sized YAML
# documents.
#
# config.h in this directory is a hand-written replacement for the header
# upstream's own configure/CMake would generate, defining only the small set of
# HAVE_*/TARGET_HAS_* macros the vendored sources actually branch on.

libfyaml_thread_dep = dependency('threads')

libfyaml_sources = files(
    'src/lib/fy-accel.c',
    'src/lib/fy-atom.c',
    'src/lib/fy-composer.c',
    'src/lib/fy-composer-diag.c',
    'src/lib/fy-diag.c',
    'src/lib/fy-docbuilder.c',
    'src/lib/fy-docbuilder-diag.c',
    'src/lib/fy-doc.c',
    'src/lib/fy-doc-diag.c',
    'src/lib/fy-docstate.c',
    'src/lib/fy-dump.c',
    'src/lib/fy-emit.c',
    'src/lib/fy-event.c',
    'src/lib/fy-input.c',
    'src/lib/fy-input-diag.c',
    'src/lib/fy-parse.c',
    'src/lib/fy-parse-diag.c',
    'src/lib/fy-path.c',
    'src/lib/fy-token.c',
    'src/lib/fy-types.c',
    'src/lib/fy-walk.c',
    'src/util/fy-blob.c',
    'src/util/fy-ctype.c',
    'src/util/fy-utf8.c',
    'src/util/fy-utils.c',
    'src/xxhash/xxhash.c',
    'src/thread/fy-thread.c',
    'src/allocator/fy-allocator.c',
    'src/allocator/fy-allocator-linear.c',
    'src/allocator/fy-allocator-malloc.c',
    'src/allocator/fy-allocator-mremap.c',
    'src/allocator/fy-allocator-dedup.c',
    'src/allocator/fy-allocator-auto.c',
    'src/blake3/blake3_portable.c',
    'src/blake3/blake3.c',
    'src/blake3/blake3_backend.c',
    'src/blake3/blake3_be_cpusimd.c',
    'src/blake3/blake3_host_state.c',
    'src/blake3/fy-blake3.c',
)

libfyaml_inc_local = include_directories(
    '.',
    'include',
    'src/lib',
    'src/util',
    'src/xxhash',
    'src/thread',
    'src/allocator',
    'src/blake3',
)

# blake3.c is upstream's generic hasher driver, built once per SIMD backend with
# these two defines set by the build system (CMake builds it once per variant
# object library). We only build the portable backend, so these are fixed.
libfyaml_c_args = [
    '-DG_LOG_DOMAIN="NUMCOSMO"',
    '-DHAVE_CONFIG_H',
    '-DHASHER_SUFFIX=portable',
    '-DSIMD_DEGREE=1',
]

disable_warnings = [
    '-Wno-unused-but-set-variable',
    '-Wno-unused-but-set-parameter',
    '-Wno-type-limits',
    '-Wno-array-bounds',
]

foreach flag : disable_warnings
  if cc.has_argument(flag)
    libfyaml_c_args += flag
  endif
endforeach

libfyaml = static_library(
    'fyaml',
    sources: libfyaml_sources,
    dependencies: [libfyaml_thread_dep],
    include_directories: libfyaml_inc_local,
    install: false,
    gnu_symbol_visibility: 'hidden',
    c_args: libfyaml_c_args,
    # Upstream requires C17 (CMAKE_C_STANDARD 17 in their CMakeLists.txt): its
    # atomics header and use of DBL_DECIMAL_DIG depend on C11-or-later features
    # that NumCosmo's project-wide '-std=gnu99' default doesn't provide.
    override_options: ['c_std=gnu17'],
)
