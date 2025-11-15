# NumCosmo Technology Stack

## Programming Languages

### Primary Languages
- **C (gnu99)** - Core library implementation
- **Fortran** - Legacy code integration and numerical routines
- **Python (≥3.9)** - High-level interface and scripting

### Language Standards
- C standard: gnu99
- Fortran standards: legacy, f77, f90, f95, f2003, f2008, f2018
- Python versions: 3.9, 3.10, 3.11, 3.12, 3.13

## Build System

### Build Tools
- **Meson (≥1.3.0)** - Primary build system
- **Ninja** - Build backend
- **pkg-config** - Dependency management

### Build Configuration
```bash
# Default build (debug optimized)
meson setup builddir
meson compile -C builddir

# Optimized build
meson setup --buildtype=release Optimized
meson compile -C Optimized

# With MPI support
meson setup -Dmpi=enabled builddir

# With all optional features
meson setup -Dfftw=enabled -Dcfitsio=enabled -Dnlopt=enabled builddir
```

## Core Dependencies

### Required Dependencies
- **GLib (≥2.54.0)** - Core utilities, object system, I/O
- **GObject (≥2.54.0)** - Object-oriented framework
- **GIO (≥2.54.0)** - I/O and networking
- **GSL (≥2.4)** - GNU Scientific Library for numerical computations
- **GMP** - Arbitrary precision arithmetic
- **MPFR** - Multiple-precision floating-point arithmetic

### Optional Dependencies
- **FFTW3 (≥3.1.2)** - Fast Fourier Transform library
- **cfitsio (≥3.25)** - FITS file I/O
- **NLopt** - Nonlinear optimization library
- **Flint (≥3.0)** - Fast Library for Number Theory
- **libfyaml** - YAML parsing library

### Linear Algebra
- **BLAS** - Basic Linear Algebra Subprograms
  - Supported implementations: MKL, OpenBLAS, Accelerate, generic
  - Auto-detection with fallback order
- **LAPACK** - Linear Algebra Package
  - Required for advanced matrix operations

### Parallel Computing
- **MPI** - Message Passing Interface for distributed computing
  - Supports MPICH and generic MPI implementations
- **OpenMP** - Shared memory parallelization
  - Available for both C and Fortran code

## Python Environment

### Python Dependencies
- **pytest** - Testing framework
- **pytest-lazy-fixtures** - Fixture utilities for pytest
- **gi (PyGObject)** - GObject Introspection bindings

### Optional Python Packages
- **firecrown** - Likelihood framework
- **pyccl** - Core Cosmology Library Python interface
- **pydantic** - Data validation
- **sacc** - Summary statistic data format
- **rich** - Terminal formatting
- **tabulate** - Table formatting
- **typer** - CLI framework

### Development Tools
- **mypy** - Static type checking
- **pylint** - Code linting
- **flake8** - Style guide enforcement

## Integrated Libraries

### Bundled Components
- **libcuba (4.0)** - Multidimensional integration
- **SUNDIALS ARKODE** - ODE/DAE solver suite
- **levmar** - Levenberg-Marquardt optimization
- **toeplitz** - Toeplitz matrix solver
- **lintegrate** - Integration routines
- **CLASS** - Cosmic Linear Anisotropy Solving System

## Development Commands

### Testing
```bash
# Run all tests
meson test -C builddir

# Run specific test
meson test -C builddir test_name

# Run Python tests
pytest tests/

# Run with coverage
meson setup -Db_coverage=true builddir
meson test -C builddir
ninja coverage -C builddir
```

### Code Quality
```bash
# Format C code
clang-format -i file.c

# Check Python style
flake8 numcosmo_py/
pylint numcosmo_py/
mypy numcosmo_py/
```

### Documentation
```bash
# Build documentation
meson setup -Ddocumentation=true builddir
meson compile -C builddir docs

# Generate API documentation
cd docs && ./build_rtd.sh
```

### Installation
```bash
# Install library
meson install -C builddir

# Install Python package
pip install -e .

# Install with optional dependencies
pip install -e ".[test,dev]"
```

## Compiler Features

### Optimization Flags
- `-funroll-loops` - Loop unrolling
- `-fomit-frame-pointer` - Frame pointer optimization
- `-ftree-vectorize` - Auto-vectorization
- `-fno-stack-protector` - Disable stack protection for performance
- `-fno-math-errno` - Disable errno for math functions
- `-march=native` - CPU-specific optimizations (optional)

### Debug Options
- `G_ENABLE_DEBUG` - Enable GLib debug infrastructure
- `G_DISABLE_CAST_CHECKS` - Disable type cast checks (release builds)
- `G_DISABLE_ASSERT` - Disable assertions (optional)
- `G_DISABLE_CHECKS` - Disable GLib checks (optional)
- `GSL_RANGE_CHECK_OFF` - Disable GSL range checking (optional)

## Platform Support
- **Linux** - Primary development platform
- **macOS** - Supported via Accelerate framework
- **Unix/POSIX** - General Unix compatibility
- **Windows/Cygwin** - Limited support (static or shared builds only)

## Version Information
- Current version: 0.26.0
- API version: 1.0
- License: GPL-3.0-or-later
- Library versioning: libtool-compatible
