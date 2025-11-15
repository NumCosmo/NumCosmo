# NumCosmo Core Modules

## Mathematical Foundation (NCM)

### Model Framework
- **NcmModel** - Base class for all parametric models
- **NcmMSet** - Model set container for multiple models
- **NcmMSetFunc** - Functions of model sets
- **NcmModelCtrl** - Model state control and caching

### Data and Likelihood
- **NcmData** - Base class for data/likelihood implementations
- **NcmDataset** - Collection of NcmData objects
- **NcmDataGauss** - Gaussian likelihood base class
- **NcmDataGaussCov** - Gaussian with covariance matrix
- **NcmDataGaussDiag** - Gaussian with diagonal covariance
- **NcmDataPoisson** - Poisson likelihood

### Fitting and Optimization
- **NcmFit** - Main fitting framework
- **NcmFitState** - Fit state management
- **NcmLikelihood** - Likelihood computation
- **NcmPrior** - Prior distributions
- **NcmFitNLOpt** - NLopt optimization backend
- **NcmFitGSL** - GSL optimization backend

### Statistical Sampling
- **NcmFitMC** - Monte Carlo sampling
- **NcmFitESMCMC** - Ensemble Sampler MCMC
- **NcmMSetCatalog** - MCMC chain catalog management
- **NcmStatsDist** - Statistical distributions
- **NcmStatsDistKernel** - Kernel density estimation
- **NcmStatsVec** - Statistical vector operations

### Linear Algebra
- **NcmMatrix** - Matrix operations
- **NcmVector** - Vector operations
- **NcmLapack** - LAPACK interface
- **NcmBlas** - BLAS interface

### Interpolation
- **NcmSpline** - 1D spline interpolation
- **NcmSpline2d** - 2D spline interpolation
- **NcmOdeSpline** - ODE-based spline
- **NcmSplineCubic** - Cubic spline
- **NcmSplineFunc** - Function-based spline

### Integration
- **NcmIntegral1d** - 1D numerical integration
- **NcmIntegralND** - N-dimensional integration
- **NcmFftlog** - FFTLog transformations
- **NcmFftlogTophatwin2** - FFTLog with tophat window

### Special Functions
- **NcmSFSBessel** - Spherical Bessel functions
- **NcmSFSphericalHarmonics** - Spherical harmonics
- **NcmMPSF** - Multiple precision special functions

### Utilities
- **NcmSerialize** - Object serialization
- **NcmCfg** - Configuration management
- **NcmDiff** - Numerical differentiation
- **NcmFunc** - Function evaluation framework
- **NcmRNG** - Random number generation
- **NcmTimer** - Timing utilities

## Cosmology Framework (NC)

### Background Cosmology
- **NcHICosmo** - Base class for homogeneous isotropic cosmologies
- **NcHICosmoDe** - Dark energy cosmologies
- **NcHICosmoDeXcdm** - XCDM parametrization
- **NcHICosmoDeWSpline** - w(z) spline parametrization
- **NcHICosmoQGRW** - Quantum gravity cosmology
- **NcDistance** - Cosmological distance calculations
- **NcScalefactor** - Scale factor evolution

### Primordial Perturbations
- **NcHIPrim** - Base class for primordial power spectra
- **NcHIPrimPowerLaw** - Power-law primordial spectrum
- **NcHIPrimAtan** - Arctan primordial spectrum
- **NcHIPrimTwoFluids** - Two-fluid primordial perturbations

### Recombination and Reionization
- **NcRecomb** - Base class for recombination
- **NcRecombSeager** - Seager recombination
- **NcRecombCBE** - CLASS-based recombination
- **NcHIReion** - Base class for reionization
- **NcHIReionCamb** - CAMB-style reionization

### Power Spectra
- **NcPowspecML** - Matter linear power spectrum
- **NcPowspecMLTransfer** - Transfer function-based power spectrum
- **NcPowspecMLCBE** - CLASS-based power spectrum
- **NcPowspecMNL** - Matter non-linear power spectrum
- **NcPowspecMNLHaloFit** - HaloFit non-linear corrections

### Transfer Functions
- **NcTransferFunc** - Base class for transfer functions
- **NcTransferFuncEH** - Eisenstein-Hu transfer function
- **NcTransferFuncBBKS** - BBKS transfer function
- **NcTransferFuncCAMB** - CAMB transfer function

### Large Scale Structure
- **NcHaloDensityProfile** - Halo density profiles (NFW, Einasto, etc.)
- **NcHaloMassFunction** - Halo mass function
- **NcHaloBias** - Halo bias models
- **NcMultiplicityFunc** - Multiplicity function
- **NcWLSurfaceMassDensity** - Weak lensing surface mass density
- **NcWindow** - Window functions

### Cluster Cosmology
- **NcClusterMass** - Cluster mass observables
- **NcClusterRedshift** - Cluster redshift observables
- **NcClusterAbundance** - Cluster abundance calculations
- **NcClusterPseudoCounts** - Pseudo-counts for clusters
- **NcDataClusterNCount** - Cluster number counts data
- **NcDataClusterWL** - Cluster weak lensing data

### Galaxy Surveys
- **NcGalaxySDPosition** - Galaxy position distribution
- **NcGalaxySDShapeGauss** - Galaxy shape distribution
- **NcGalaxySDObsRedshift** - Observed redshift distribution
- **NcGalaxySDTrueRedshift** - True redshift distribution
- **NcGalaxyWL** - Galaxy weak lensing
- **NcGalaxyACF** - Galaxy angular correlation function

### Cross-Correlation
- **NcXcor** - Cross-correlation framework
- **NcXcorLimberKernel** - Limber approximation kernels
- **NcXcorLimberKernelGal** - Galaxy kernel
- **NcXcorLimberKernelCMBLensing** - CMB lensing kernel
- **NcXcorLimberKernelWeakLensing** - Weak lensing kernel

### Observational Data
- **NcDataBao** - Baryon Acoustic Oscillations data
- **NcDataHubble** - Hubble parameter measurements
- **NcDataSNIa** - Type Ia Supernovae data
- **NcDataCMB** - CMB distance priors
- **NcPlanckFI** - Planck Fisher information

### CLASS Integration
- **NcCBE** - CLASS backend engine
- **NcCBEPrecision** - CLASS precision parameters
- **NcPowspecMLCBE** - CLASS power spectrum
- **NcRecombCBE** - CLASS recombination
