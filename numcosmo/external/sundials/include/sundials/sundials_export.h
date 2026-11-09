
#ifndef SUNDIALS_EXPORT_H
#define SUNDIALS_EXPORT_H

#ifdef SUNDIALS_STATIC_DEFINE
#  define SUNDIALS_EXPORT
#  define SUNDIALS_NO_EXPORT
#else
#  ifndef SUNDIALS_EXPORT
#    ifdef sundials_core_EXPORTS
        /* We are building this library */
#      define SUNDIALS_EXPORT 
#    else
        /* We are using this library */
#      define SUNDIALS_EXPORT 
#    endif
#  endif

#  ifndef SUNDIALS_NO_EXPORT
#    define SUNDIALS_NO_EXPORT 
#  endif
#endif

#ifndef SUNDIALS_DEPRECATED_ATTRIBUTE
#  define SUNDIALS_DEPRECATED_ATTRIBUTE __attribute__ ((__deprecated__))
#endif

#ifndef SUNDIALS_DEPRECATED_ATTRIBUTE_EXPORT
#  define SUNDIALS_DEPRECATED_ATTRIBUTE_EXPORT SUNDIALS_EXPORT SUNDIALS_DEPRECATED_ATTRIBUTE
#endif

#ifndef SUNDIALS_DEPRECATED_ATTRIBUTE_NO_EXPORT
#  define SUNDIALS_DEPRECATED_ATTRIBUTE_NO_EXPORT SUNDIALS_NO_EXPORT SUNDIALS_DEPRECATED_ATTRIBUTE
#endif

/* NOLINTNEXTLINE(readability-avoid-unconditional-preprocessor-if) */
#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef SUNDIALS_NO_DEPRECATED
#    define SUNDIALS_NO_DEPRECATED
#  endif
#endif

#endif /* SUNDIALS_EXPORT_H */
