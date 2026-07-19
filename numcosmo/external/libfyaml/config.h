/*
 * config.h - hand-written configuration header for the vendored libfyaml subset
 *
 * NumCosmo only vendors the core parser/emitter library (see meson.build in this
 * directory), not the upstream CLI tool, reflection/libclang backend, python
 * bindings, or test harness. This header only defines the macros that the
 * vendored sources actually branch on, using values that are safe and portable
 * across the Linux/macOS toolchains NumCosmo supports.
 */
#ifndef NUMCOSMO_LIBFYAML_CONFIG_H
#define NUMCOSMO_LIBFYAML_CONFIG_H

#define VERSION "0.9.6"

/* alloca.h is available on Linux and macOS. */
#define HAVE_ALLOCA_H 1

/*
 * byteswap.h is glibc/Linux-specific and absent on macOS. Skip it and rely on
 * the __builtin_bswap* fallback below instead, which is portable to every
 * GCC/Clang toolchain this project already requires.
 */
/* #undef HAVE_BYTESWAP_H */
#define HAVE___BUILTIN_BSWAP16 1
#define HAVE___BUILTIN_BSWAP32 1
#define HAVE___BUILTIN_BSWAP64 1

/*
 * qsort_r() has an incompatible argument order between glibc and BSD/macOS.
 * Rather than probe and branch on it, keep it disabled: the fallback path
 * (plain qsort() plus a static sort context) is fully correct for NumCosmo's
 * single-threaded (de)serialization use.
 */
#define HAVE_QSORT_R 0

/*
 * mremap() is a Linux-only syscall (absent on macOS/BSD). The dedup allocator
 * falls back to a portable alloc+copy path when this is 0.
 */
#define HAVE_MREMAP 0

/* No libclang/LLVM reflection backend is vendored or linked. */
#define HAVE_LIBCLANG 0
/* #undef HAVE_LLVM_C_CORE_H */

/* Portable BLAKE3 backend only -- see meson.build for the rationale. */
/* #undef TARGET_HAS_SSE2 */
/* #undef TARGET_HAS_SSE41 */
/* #undef TARGET_HAS_AVX2 */
/* #undef TARGET_HAS_AVX512 */
/* #undef TARGET_HAS_NEON */

#endif /* NUMCOSMO_LIBFYAML_CONFIG_H */
