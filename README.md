NumCosmo
========

NumCosmo is a numerical cosmology library. It contains a comprehensive set of
tools for calculating cosmology observables and to analyze statistical models.

Description
-----------

The library is written in C, but since it uses the [GObject](https://wiki.gnome.org/action/show/Projects/GObjectIntrospection) 
framework, it is developed in a object oriented fashion. Additionally, it has automatic
bindings for every language which supports GObject introspection (Perl,
Python, etc. For a complete list see https://wiki.gnome.org/Projects/GObjectIntrospection/Users ).

The available observables objects are:
  - Type Ia Supernovae.
  - Baryon Acoustic Oscillations.
  - Cosmic Microwave Background (shift parameter and distance priors, full
  analysis is on the way).
  - Cluster Number counts.
  - Hubble data H(z).

Currently it has the following statistical tools:
  - Monte Carlo (NcmFitMC) -- resampling and fitting.
  - Monte Carlo Bootstrap (NcmFitMCBS) -- resampling/bootstraping and fitting.
  - Markov Chain Monte Carlo -- MCMC with the Metropolis-Hastings sampler, it supports
    general samplers through transtion kernel object NcmMSetTransKern.
  - Ensemble Sampler MCMC -- Ensemble Markov Chain Monte Carlo consists in
    every point of the MCMC chain being a emsemble of points in the
    parameter space. It implements an affine invariant move method (stretch move).

All methods above generate a catalog using the NcmMSetCatalog which provide
a unified way to analyze the results. Besides, the use of a catalog
provides the support for restarting the algorithms from a previous crash
or to extend the precision.

Links:
------

The Savannah project page can be found [here](https://savannah.nongnu.org/projects/numcosmo/).
The NumCosmo homepage is [here](http://www.nongnu.org/numcosmo/), it contains the documentation and some examples.

Building from repository:
-------------------------

To build from the git repository first run ./autogen.sh to build the
configure script. Note that this requires the autotools developer 
enviroment (autoconf, automake, intltool, libtool, ...).

To build from a release package (or after running ./autogen.sh), run: 
./configure (--help to see options)
make
make install

For a generic installation instructions, see INSTALL.

Requirements:
-------------

  - Glib >= 2.28.0
    Data structures, threads, portability, memory allocation, etc.
    * http://www.gtk.org/
  - GSL  >= 1.15
    Several computational tools.
    * http://www.gnu.org/software/gsl/
  - GMP  >= 4.3.2
    Big integers library.
    * http://gmplib.org/
  - MPFR >= 2.4.2
    Multiple precision float library.
    * http://www.mpfr.org/
  - Sundials >= 2.4.0
    ODE solver library. 
    * https://computation.llnl.gov/casc/sundials/main.html

Optional packages:
------------------

NumCosmo, provides it own internal version of Levmar and Cuba for a full version 
of the original libraries visit their website linked below. It is, therefore,
*not* necessary to install them.

  - FFTW3 >= 3.1.2
    Needed to build the spherical harmonic decomposition of CMB data.
    * http://www.fftw.org/
  - Any optimized BLAS library (ATLAS, OpenBLAS, MKL, etc)
    Improve speed in linear algebra calculations.
    * http://math-atlas.sourceforge.net/
    * http://www.openblas.net/
    * https://software.intel.com/en-us/intel-mkl
  - Lapack
    Linear Algebra PACKage    
    * http://www.netlib.org/lapack/
  - Extra (besides gsl's) minimization packages
    - Levmar
      Least squares minimization library.
      * http://www.ics.forth.gr/~lourakis/levmar/
    - Cuba 
      A library for multidimensional numerical integration
      * http://www.feynarts.de/cuba/  
    - NLOpt
      Several general purpose minimization algorithms.
      * http://ab-initio.mit.edu/wiki/index.php/NLopt
  - CFITSIO
    C library used to manipulate fits files. It will build the support to
    read CMB maps and astronomical data in general in fits format.
    * http://heasarc.nasa.gov/fitsio/fitsio.html
  - gtk-doc
    GTK-Doc is used to generate API documentation from comments added to 
    C code, needed only to generate new releases.
