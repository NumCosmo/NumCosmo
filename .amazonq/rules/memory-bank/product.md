# NumCosmo Product Overview

## Purpose
NumCosmo is a comprehensive numerical cosmology library designed for cosmological calculations and statistical model analysis. It provides researchers and developers with powerful tools to compute cosmological observables, perform statistical inference, and analyze data in cosmology and astrophysics.

## Key Features

### Cosmological Calculations
- Distance calculations and cosmological evolution
- Power spectrum computations (linear and non-linear)
- Transfer functions and recombination physics
- Halo mass functions and density profiles
- Weak lensing surface mass density
- Baryon Acoustic Oscillations (BAO) analysis
- Cosmic Microwave Background (CMB) calculations
- Supernova Type Ia (SNIa) distance modulus

### Statistical Analysis
- Advanced fitting algorithms (MCMC, ESMCMC, Monte Carlo)
- Multiple optimization backends (GSL, NLopt)
- Statistical distributions and kernel density estimation
- Likelihood analysis for various cosmological probes
- Model parameter estimation and constraints
- Catalog management for MCMC chains

### Data Analysis
- Support for multiple observational datasets (BAO, Hubble, SNIa, CMB, clusters)
- Covariance matrix handling
- Data serialization and deserialization
- FITS file support for astronomical data

### Mathematical Tools
- Numerical integration (1D and N-dimensional)
- Spline interpolation (1D and 2D)
- Matrix and vector operations
- Special functions (spherical Bessel, spherical harmonics)
- FFT and FFTLog transformations
- ODE solvers (via integrated SUNDIALS)

## Target Users

### Primary Users
- Cosmologists and astrophysicists conducting research
- Scientists analyzing observational cosmological data
- Researchers developing new cosmological models
- Graduate students learning computational cosmology

### Use Cases
- Constraining cosmological parameters from observational data
- Computing theoretical predictions for cosmological observables
- Performing Bayesian inference on cosmological models
- Comparing different cosmological models
- Generating mock catalogs for survey simulations
- Cross-correlation analysis between different cosmological probes

## Technology Stack
- Core library written in C with GObject framework
- Python bindings via GObject Introspection
- Fortran components for legacy code integration
- Support for parallel computing (MPI, OpenMP)
- Integration with scientific libraries (GSL, FFTW, LAPACK/BLAS)
