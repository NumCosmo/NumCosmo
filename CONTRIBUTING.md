# Contributing to NumCosmo

Thank you for your interest in contributing to NumCosmo! We welcome contributions from the community to help improve and expand the library. Before getting started, please review the guidelines and requirements outlined below.

## Adding New Files

When adding new files to NumCosmo, please follow these steps:

1. Add the new files to the appropriate sections in the `numcosmo/Makefile.am` file:
   - If it is a generic math or statistics code, include the file in `ncm_sources` and `ncm_headers`.
   - If it is cosmology-related, include the file in `nc_sources` and `nc_headers`.

2. Ensure the new files are also included in `numcosmo-math.h` or `numcosmo.h` as necessary.

3. Include the new files in the appropriate section of `docs/numcosmo-docs.sgml`. If no suitable section exists, create a new one.

4. Add the header inclusion to `numcosmo/math/ncm_cfg.h`. All new objects should be registered using `ncm_cfg_register_obj (NCM_TYPE_NEW_OBJECT);`.

## Unit Testing

All new code must have at least minimal unit testing. Follow these steps to include unit tests for your new code:

1. Create unit tests for the new code in the `tests` directory.

2. Include the new files in `tests/Makefile.am` following the provided pattern:

```
test_nc_new_test_SOURCES =  \
        test_nc_new_test.c

test_programs =  \
        test_ncm_cfg    \
        test_ncm_vector \
...
        test_nc_new_test \
...

test_nc_new_test_LDADD = \
        $(top_builddir)/numcosmo/libnumcosmo.la \
        $(GLIB_LIBS) \
        $(GSL_LIBS) \
        $(COVLIBS)

```

3. Ensure the unit tests link against the necessary libraries and dependencies.

## Code Formatting

All C files, including headers, must be formatted using uncrustify with the provided configuration file `numcosmo_uncrustify.cfg`.

## Submitting Contributions

When submitting a contribution to NumCosmo, please follow these steps:

1. Fork the NumCosmo repository to your GitHub account. If you in the development team you can skip this.

2. Create a new branch for your contribution based on the latest `master` branch.

3. Make your changes and ensure all new code adheres to the project guidelines, including file inclusion, unit testing, and code formatting.

4. Commit your changes and push them to your forked repository.

5. Submit a pull request to the main NumCosmo repository. Provide a clear description of your changes and any relevant details.

6. Your contribution will be reviewed by the maintainers, and they may provide feedback or request additional changes.

7. Once approved, your contribution will be merged into the main repository.

## Code of Conduct

Please note that by contributing to NumCosmo, you are expected to abide by the [Code of Conduct](CODE_OF_CONDUCT.md). Ensure your interactions and contributions align with the principles of respect and inclusivity.

We appreciate your interest in contributing to NumCosmo and look forward to your contributions!

