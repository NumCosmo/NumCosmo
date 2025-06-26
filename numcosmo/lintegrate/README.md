# lintegrate

A numerical integration library for when you want/need to work with the natural logarithm of the function requiring integration.

This library provides three numerical integration functions, heavily based on [GSL](https://www.gnu.org/software/gsl/) functions, to integrate a function when only its natural logarithm is given, and return the natural logarithm of that integral. The three functions:

 * `lintegrate_qag`
 * `lintegrate_qng`
 * `lintegrate_cquad`

 are equivalents of the GSL functions:

 * [`gsl_integration_qag`](https://www.gnu.org/software/gsl/doc/html/integration.html#qag-adaptive-integration)
 * [`gsl_integration_qng`](https://www.gnu.org/software/gsl/doc/html/integration.html#qng-non-adaptive-gauss-kronrod-integration)
 * [`gsl_integration_cquad`](https://www.gnu.org/software/gsl/doc/html/integration.html#cquad-doubly-adaptive-integration)

respectively. These can be useful when, e.g., you can calculate the natural logarithm of a Gaussian likelihood function (in cases where the exponentiation of the Gaussian function would lead to zeros or infinities) and you want to numerically find the integral of the Gaussian function itself.

The functions `lintegrate_qag`, `lintegrate_qng`, and `lintegrate_cquad`, all have wrappers functions (with `_split` appended to their names) that allow the user to specify a set of intervals that the integrals will be split into when performing the calculation. The intervals could, for example, be spaced evenly in log-space, for cases where the integral function has a very pronounced peak as it approaches zero.

The full API documentation and examples can be found [here](https://lintegrate.readthedocs.io/).

## Example

An [example](example/example.c) of the use the functions is:

```C
/* example using lintegrate functionality */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include <lintegrate.h>

/* create function for integration */
double lintegrand(double x, void *params);

struct intparams {
  double mu;
  double sig;
};

double lintegrand(double x, void *params){
  struct intparams * p = (struct intparams *)params;
  double mu = p->mu;
  double sig = p->sig;

  return -0.5*(mu-x)*(mu-x)/(sig*sig);
}

double integrand(double x, void *params){
  struct intparams * p = (struct intparams *)params;
  double mu = p->mu;
  double sig = p->sig;

  return exp(-0.5*(mu-x)*(mu-x)/(sig*sig));
}

int main( int argv, char **argc ){
  gsl_function F;
  struct intparams params;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc (100);
  gsl_integration_cquad_workspace *cw = gsl_integration_cquad_workspace_alloc(50);
  double qaganswer = 0., qnganswer = 0., cquadanswer = 0., answer = 0.;
  double abserr = 0.;
  size_t neval = 0;

  double minlim = -6.; /* minimum for integration range */
  double maxlim = 6.;  /* maximum for integration range */

  double abstol = 1e-10; /* absolute tolerance */
  double reltol = 1e-10; /* relative tolerance */

  params.mu = 0.;
  params.sig = 1.;

  F.function = &lintegrand;
  F.params = &params;

  /* integrate log of function using QAG */
  lintegration_qag(&F, minlim, maxlim, abstol, reltol, 100, GSL_INTEG_GAUSS31, w, &qaganswer, &abserr);

  /* integrate log of function using QNG */
  lintegration_qng(&F, minlim, maxlim, abstol, reltol, &qnganswer, &abserr, &neval);

  /* integrate log of function using CQUAD */
  lintegration_cquad(&F, minlim, maxlim, abstol, reltol, cw, &cquadanswer, &abserr, &neval);

  /* integrate function using GSL QAG */
  F.function = &integrand;
  gsl_integration_qag(&F, minlim, maxlim, abstol, reltol, 100, GSL_INTEG_GAUSS31, w, &answer, &abserr);

  gsl_integration_workspace_free(w);
  gsl_integration_cquad_workspace_free(cw);

  fprintf(stdout, "Answer \"lintegrate QAG\" = %.8lf\n", qaganswer);
  fprintf(stdout, "Answer \"lintegrate QNG\" = %.8lf\n", qnganswer);
  fprintf(stdout, "Answer \"lintegrate CQUAD\" = %.8lf\n", cquadanswer);
  fprintf(stdout, "Answer \"gsl_integrate_qag\" = %.8lf\n", log(answer));
  fprintf(stdout, "Analytical answer = %.8lf\n", log(sqrt(2.*M_PI)));

  return 0;
}
```

## Requirements

* [GSL](https://www.gnu.org/software/gsl/) - on Debian/Ubuntu (16.04) install with e.g. `sudo apt-get install libgsl-dev`

## Installation

The library can be built using [scons](http://scons.org) by just typing `sudo scons` in the base directory. To install
the library system-wide (in `/usr/local/lib` by default) run:
```
sudo scons
sudo scons install
```

A Python module containing wrappers to the functions can be built and installed from source for the user by running, e.g.:
```bash
pip install .
```
from within the repository directory.

The Python module can also be installed from [PyPI](https://pypi.org/project/lintegrate/) using pip with:
```bash
pip install lintegrate
```

or in a Conda environment with:
```bash
conda install -c conda-forge lintegrate
```

## Python

If the Python module has been installed it has the following functions:
 * `lqng` - a wrapper to `lintegration_qng`
 * `lqag` - a wrapper to `lintegration_qag`
 * `lcquad` - a wrapper to `lintegration_cquad`
 * `logtrapz` - using the trapezium rule for integration on a grid of values

The `lqng`, `lqag`, and `lcquad` functions are used in a similar way to the scipy [`quad`](https://docs.scipy.org/doc/scipy-0.19.0/reference/generated/scipy.integrate.quad.html) function.

An example of their use would be:

```python
from lintegrate import lqag, lqng, lcquad, logtrapz
import numpy as np

# define the log of the function to be integrated
def integrand(x, args):
    mu, sig = args # unpack extra arguments
    return -0.5*((x-mu)/sig)**2

# set integration limits
xmin = -6.
xmax = 6.

# set additional arguments
mu = 0.
sig = 1.

resqag = lqag(integrand, xmin, xmax, args=(mu, sig))
resqng = lqng(integrand, xmin, xmax, args=(mu, sig))
rescquad = lcquad(integrand, xmin, xmax, args=(mu, sig))
restrapz = logtrapz(integrand, np.linspace(xmin, xmax, 100), args=(mu, sig))
```

## R

In [R](https://www.r-project.org/) one can use the [**reticulate**](https://github.com/rstudio/reticulate) package to call the functions in `lintegrate`.
The above example would be:
```R
library(reticulate)
py_install("lintegrate", pip = TRUE) ## run once to make sure lintegrate is installed and visible to reticulate.
lint <- import("lintegrate", convert = FALSE)
integrand <- function(x, args){
  mu = args[1]
  sig = args[2]
  return(-.5 * ((x-mu)/sig)^2 )
} 
integrand <- Vectorize(integrand)
mu <- 0
sig <- 1
mmin <- -10
mmax <- 10
lint$lqag(py_func(integrand), r_to_py(mmin), r_to_py(mmax), c(mu, sig))
```

## Citation

If using `lintegrate` in your research, I would be grateful if you cite the associated [JOSS paper](https://joss.theoj.org/papers/10.21105/joss.04231) for the software. The following BibTeX citation can be used:

```bibtex
@article{Pitkin2022,
  doi = {10.21105/joss.04231},
  url = {https://doi.org/10.21105/joss.04231},
  year = {2022},
  publisher = {The Open Journal},
  volume = {7},
  number = {73},
  pages = {4231},
  author = {Matthew Pitkin},
  title = {lintegrate: A C/Python numerical integration library for working in log-space},
  journal = {Journal of Open Source Software}
}
```

You may also want to cite the [GSL](https://www.gnu.org/software/gsl/) reference "_M. Galassi et al, GNU Scientific Library Reference Manual (3rd Ed.), ISBN 0954612078_" and the URL http://www.gnu.org/software/gsl/.


[![DOI](https://joss.theoj.org/papers/10.21105/joss.04231/status.svg)](https://doi.org/10.21105/joss.04231)
[![Build Status](https://github.com/mattpitkin/lintegrate/workflows/build/badge.svg)](https://github.com/mattpitkin/lintegrate/actions?query=workflow%3Abuild)
[![PyPI version](https://badge.fury.io/py/lintegrate.svg)](https://badge.fury.io/py/lintegrate)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/lintegrate/badges/version.svg)](https://anaconda.org/conda-forge/lintegrate)
[![Documentation Status](https://readthedocs.org/projects/lintegrate/badge/?version=latest)](http://lintegrate.readthedocs.io/en/latest/?badge=latest)

&copy; 2017 Matthew Pitkin
