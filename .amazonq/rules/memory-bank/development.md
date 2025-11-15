# NumCosmo Development Guide

## Development Setup

### Prerequisites
- C compiler (GCC or Clang)
- Meson (≥1.3.0) and Ninja
- Python (≥3.9)
- GLib/GObject development files
- GSL (GNU Scientific Library)
- Git for version control

### Initial Setup
```bash
# Clone repository
git clone https://github.com/NumCosmo/NumCosmo.git
cd NumCosmo

# Configure build
meson setup builddir

# Compile
meson compile -C builddir

# Run tests
meson test -C builddir

# Install (optional)
meson install -C builddir
```

### Development Environment
```bash
# Create conda environment
conda env create -f devel_environment.yml
conda activate numcosmo-dev

# Or use pip
pip install -e ".[dev,test]"
```

## Code Style

### C Code Style
- Follow GNU coding style with modifications
- Configuration in `.clang-format`
- 2-space indentation
- 120 character line limit
- Use `clang-format` for automatic formatting

### Formatting C Code
```bash
# Format single file
clang-format -i file.c

# Format all C files
find numcosmo -name "*.c" -o -name "*.h" | xargs clang-format -i
```

### Python Code Style
- Follow PEP 8 style guide
- Configuration in `.flake8` and `.pylintrc`
- 4-space indentation
- 88 character line limit (Black style)

### Linting Python Code
```bash
# Check style
flake8 numcosmo_py/

# Lint code
pylint numcosmo_py/

# Type checking
mypy numcosmo_py/
```

## Build Configurations

### Debug Build
```bash
meson setup --buildtype=debug builddir
```
- Enables debug symbols
- Disables optimizations
- Enables GLib debug checks

### Release Build
```bash
meson setup --buildtype=release Optimized
```
- Maximum optimizations
- Disables debug checks
- Smaller binary size

### Custom Build Options
```bash
# Enable MPI
meson setup -Dmpi=enabled builddir

# Enable optional dependencies
meson setup -Dfftw=enabled -Dcfitsio=enabled -Dnlopt=enabled builddir

# Build documentation
meson setup -Ddocumentation=true builddir

# Enable coverage
meson setup -Db_coverage=true builddir
```

## Adding New Features

### Adding a New C Class

1. Create header file `numcosmo/nc_new_class.h`:
```c
#ifndef _NC_NEW_CLASS_H_
#define _NC_NEW_CLASS_H_

#include <glib-object.h>

G_BEGIN_DECLS

#define NC_TYPE_NEW_CLASS (nc_new_class_get_type ())
G_DECLARE_FINAL_TYPE (NcNewClass, nc_new_class, NC, NEW_CLASS, GObject)

NcNewClass *nc_new_class_new (void);
void nc_new_class_method (NcNewClass *self);

G_END_DECLS

#endif /* _NC_NEW_CLASS_H_ */
```

2. Create implementation file `numcosmo/nc_new_class.c`:
```c
#include "nc_new_class.h"

struct _NcNewClass
{
  GObject parent_instance;
  /* Private fields */
};

G_DEFINE_TYPE (NcNewClass, nc_new_class, G_TYPE_OBJECT)

static void
nc_new_class_init (NcNewClass *self)
{
  /* Initialize instance */
}

static void
nc_new_class_class_init (NcNewClassClass *klass)
{
  /* Initialize class */
}

NcNewClass *
nc_new_class_new (void)
{
  return g_object_new (NC_TYPE_NEW_CLASS, NULL);
}
```

3. Add to `numcosmo/meson.build`:
```python
numcosmo_sources += files('nc_new_class.c')
numcosmo_headers += files('nc_new_class.h')
```

4. Add to `numcosmo/numcosmo.h`:
```c
#include <numcosmo/nc_new_class.h>
```

### Adding a New Python Module

1. Create file `numcosmo_py/new_module.py`:
```python
"""New module description."""

from numcosmo_py import Nc, Ncm

def new_function():
    """Function description."""
    pass
```

2. Add to `numcosmo_py/__init__.py`:
```python
from . import new_module
```

3. Add tests in `tests/test_py_new_module.py`

## Documentation

### C Documentation (GTK-Doc)
```c
/**
 * nc_new_class_method:
 * @self: a #NcNewClass
 * @param: parameter description
 *
 * Method description.
 *
 * Returns: return value description
 */
```

### Python Documentation (Sphinx/NumPy style)
```python
def new_function(param):
    """Brief description.
    
    Detailed description.
    
    Parameters
    ----------
    param : type
        Parameter description
        
    Returns
    -------
    type
        Return value description
    """
```

### Building Documentation
```bash
# Configure with documentation
meson setup -Ddocumentation=true builddir

# Build documentation
meson compile -C builddir docs

# Build ReadTheDocs documentation
cd docs && ./build_rtd.sh
```

## Version Control

### Branch Strategy
- `main` - Stable release branch
- `develop` - Development branch
- Feature branches: `feature/description`
- Bug fix branches: `fix/description`

### Commit Messages
- Use conventional commits format
- Format: `type(scope): description`
- Types: feat, fix, docs, style, refactor, test, chore
- Example: `feat(distance): add new comoving distance method`

### Pull Request Process
1. Create feature branch from `develop`
2. Implement changes with tests
3. Ensure all tests pass
4. Update documentation
5. Submit pull request to `develop`
6. Address review comments
7. Merge after approval

## Debugging

### Debug Macros
```c
// Enable debug output
#define NCM_DEBUG 1

// Debug print
NCM_DEBUG_PRINT ("Debug message: %f\n", value);
```

### GDB Debugging
```bash
# Run with gdb
gdb --args builddir/program args

# Common gdb commands
(gdb) break nc_function
(gdb) run
(gdb) print variable
(gdb) backtrace
```

### Memory Debugging
```bash
# Valgrind memory check
valgrind --leak-check=full builddir/program

# With suppressions
valgrind --suppressions=numcosmo-valgrind.supp builddir/program
```

## Performance Optimization

### Profiling
```bash
# Compile with profiling
meson setup --buildtype=debugoptimized builddir

# Run with gprof
gprof builddir/program gmon.out > analysis.txt

# Run with perf
perf record builddir/program
perf report
```

### Optimization Flags
- `-O3` - Maximum optimization
- `-march=native` - CPU-specific optimizations
- `-funroll-loops` - Loop unrolling
- `-ftree-vectorize` - Auto-vectorization

## Release Process

### Version Bumping
1. Update version in `meson.build`
2. Update `ChangeLog.md`
3. Tag release: `git tag -a v0.X.Y -m "Release v0.X.Y"`
4. Push tag: `git push origin v0.X.Y`

### Release Checklist
- [ ] All tests pass
- [ ] Documentation updated
- [ ] ChangeLog updated
- [ ] Version bumped
- [ ] Tag created
- [ ] GitHub release created
- [ ] PyPI package uploaded (if applicable)

## Troubleshooting

### Common Build Issues
- Missing dependencies: Install required packages
- Meson version: Ensure Meson ≥1.3.0
- Python version: Ensure Python ≥3.9
- BLAS/LAPACK: Check library detection

### Common Runtime Issues
- GObject warnings: Check object lifecycle
- Memory leaks: Run valgrind
- Segmentation faults: Use gdb for backtrace
- Python import errors: Check PYTHONPATH and GI_TYPELIB_PATH
