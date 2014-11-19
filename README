NumCosmo
========

This is the readme file of the numcosmo library.

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

  - SQLite3 >= 3.6.10
    Standalone database used to store several observation data (SN Ia, BAO
    Dv, etc), some data are already included in numcosmo.
    * http://www.sqlite.org/
  - FFTW3   >= 3.1.2
    Needed to build the spherical harmonic decomposition of CMB data.
    * http://www.fftw.org/
  - Atlas (cblas) any version
    Improve speed in linear algebra calculations.
    * http://math-atlas.sourceforge.net/
  - Extra (besides gsl's) minimization packages
    - Levmar
      Least squares minimization library.
      * http://www.ics.forth.gr/~lourakis/levmar/
    - NLOpt
      Several general purpose minimization algorithms.
      * http://ab-initio.mit.edu/wiki/index.php/NLopt
  - CFITSIO
    C library used to manipulate fits files. It will build the support to
    read CMB maps and astronomical data in general in fits format.
    * http://heasarc.nasa.gov/fitsio/fitsio.html
  - gtk-doc
    GTK-Doc is used to generate API documentation from comments added to 
    C code.
