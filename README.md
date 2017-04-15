NumCosmo
========

NumCosmo is a numerical cosmology library. It contains a comprehensive set of
tools for calculating cosmology observables and to analyze statistical models.

[![Build Status](https://travis-ci.org/NumCosmo/NumCosmo.svg?branch=master)](https://travis-ci.org/NumCosmo/NumCosmo)

Description
-----------

The library is written in C, but since it uses the [GObject](https://wiki.gnome.org/action/show/Projects/GObjectIntrospection) 
framework, it is developed in a object oriented fashion. Additionally, it has automatic
bindings for every language which supports GObject introspection (Perl,
Python, etc. For a complete list see https://wiki.gnome.org/Projects/GObjectIntrospection/Users ).

The available observables objects are:
  - Type Ia Supernovae
  - Baryon Acoustic Oscillations
  - Cosmic Microwave Background
  - Cluster number counts
  - Cluster pseudo number counts
  - Hubble data H(z)

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

To build from a release package, for which the configure script is ready, run: 
  - ./configure (--help to see options)
  - make
  - make install ([optional step](http://www.nongnu.org/numcosmo/manual/compiling.html))

To build from the git repository, run:
  - ./autogen.sh 
    - The configure script is built at this point. 
      Note that this requires the autotools developer enviroment (latest version): 
      - autoconf 
        * http://ftp.gnu.org/gnu/autoconf/ 
      - automake
        * http://ftp.gnu.org/gnu/automake/
      - libtool
        * http://ftp.gnu.org/gnu/libtool/
      And also:
      - gtk-doc
        * http://www.gtk.org/gtk-doc/
      - gobject-introspection
        * https://wiki.gnome.org/action/show/Projects/GObjectIntrospection?action=show&redirect=GObjectIntrospection
  - ./cofigure (--help to see options)
  - make 
  - make install ([optional step](http://www.nongnu.org/numcosmo/manual/compiling.html))

For a generic installation instructions, see INSTALL.

The requirements below can be found on most Linux distribution, see [here](http://www.nongnu.org/numcosmo/manual/compiling.html)
for a list of packages names for some distributions. 

A pre-compiled version of NumCosmo can be found [here](https://build.opensuse.org/package/show/home:vitenti/numcosmo)

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
  - GObject-introspection
    Middleware layer between C libraries (using GObject) and language bindings. 
    *https://wiki.gnome.org/action/show/Projects/GObjectIntrospection?action=show&redirect=GObjectIntrospection
    This package is optional but highly recommended.
    To use NumCosmo from Python, for example, you also need PyGObject.
    *https://wiki.gnome.org/action/show/Projects/PyGObject?action=show&redirect=PyGObject
    Note that pygobject3 refers to the PyGObject version (not Python's version).

Optional packages:
------------------

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
  - ARB
    C library for arbitrary-precision interval arithmetic.
    * http://fredrikj.net/arb/ 
