# NumCosmo Project Structure

## Root Directory Organization

### Core Source Code
- **numcosmo/** - Main C library source code
  - **math/** - Mathematical utilities and numerical methods
  - **model/** - Cosmological model implementations
  - **lss/** - Large Scale Structure calculations
  - **perturbations/** - Cosmological perturbation theory
  - **galaxy/** - Galaxy-related calculations
  - **xcor/** - Cross-correlation analysis
  - **data/** - Data handling and likelihood implementations
  - **class/** - CLASS (Cosmic Linear Anisotropy Solving System) integration
  - **sundials/** - Integrated SUNDIALS ODE solver library
  - **libcuba/** - Cuba integration library
  - **levmar/** - Levenberg-Marquardt optimization
  - **toeplitz/** - Toeplitz matrix operations
  - **lintegrate/** - Integration routines
  - **plc/** - Planck likelihood code

### Python Interface
- **numcosmo_py/** - Python package and bindings
  - **app/** - Application-level Python tools
  - **ccl/** - CCL (Core Cosmology Library) compatibility
  - **datasets/** - Dataset handling utilities
  - **experiments/** - Experimental analysis tools
  - **external/** - External library interfaces
  - **interpolation/** - Interpolation utilities
  - **plotting/** - Visualization tools
  - **sampling/** - Sampling algorithms

### Testing and Examples
- **tests/** - Comprehensive test suite (C and Python)
- **examples/** - Example programs demonstrating library usage
- **notebooks/** - Jupyter notebooks for tutorials and analysis

### Documentation
- **docs/** - Documentation source files
  - **reference/** - API reference documentation
  - **tutorials/** - Tutorial documents
  - **examples/** - Example documentation
  - **benchmarks/** - Performance benchmarks
  - **notes/** - Development notes

### Data and Configuration
- **data/** - Observational data files and serialized objects
  - BAO measurements from various surveys
  - Hubble parameter measurements
  - SNIa datasets
  - CMB priors
  - Cluster data
- **tools/** - Command-line tools and utilities

### Build Directories
- **builddir/** - Default build directory (Meson)
- **BDocs/** - Documentation build directory
- **Optimized/** - Optimized build directory

## Architectural Patterns

### Object-Oriented C with GObject
- Uses GObject type system for object-oriented programming in C
- Reference counting for memory management
- Properties and signals for object interaction
- Introspection support for language bindings

### Modular Design
- Clear separation between mathematical utilities (ncm_*) and cosmology-specific code (nc_*)
- Plugin architecture for models and data types
- Extensible framework for adding new models and likelihoods

### Namespace Convention
- **NCM** prefix - NumCosmo Math (general mathematical utilities)
- **NC** prefix - NumCosmo (cosmology-specific implementations)
- Hierarchical naming: nc_[module]_[class]_[method]

### Core Components

#### Mathematical Foundation (NCM)
- NcmModel - Base class for all models
- NcmMSet - Model set management
- NcmData - Base class for data/likelihood
- NcmDataset - Collection of data objects
- NcmFit - Fitting framework
- NcmMatrix, NcmVector - Linear algebra
- NcmSpline, NcmSpline2d - Interpolation
- NcmIntegral - Numerical integration
- NcmStatsDist - Statistical distributions

#### Cosmology Framework (NC)
- NcHICosmo - Base class for homogeneous and isotropic cosmologies
- NcDistance - Cosmological distance calculations
- NcPowspec - Power spectrum interface
- NcTransferFunc - Transfer function calculations
- NcHaloDensityProfile - Halo profiles
- NcClusterAbundance - Cluster abundance calculations
- NcXcor - Cross-correlation framework

### Build System
- Meson build system with Ninja backend
- Modular configuration with feature detection
- Support for multiple BLAS/LAPACK implementations
- Optional dependencies (MPI, OpenMP, FFTW, cfitsio, NLopt, Flint)
- Python integration via pip installation
