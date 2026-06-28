# NumCosmo

NumCosmo is a C library for cosmological calculations and statistical analysis. It is written in C using GObject, so every class is usable from Python (and other GObject-Introspection languages) through the `numcosmo_py` package without writing binding code.

[![Build Status](https://github.com/NumCosmo/NumCosmo/workflows/Build%20and%20Check/badge.svg)](https://github.com/NumCosmo/NumCosmo/actions) [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![codecov](https://codecov.io/gh/NumCosmo/NumCosmo/graph/badge.svg?token=FZ3PX0PKWG)](https://codecov.io/gh/NumCosmo/NumCosmo)

Visit [NumCosmo's website](https://numcosmo.readthedocs.io/en/latest/) for more information.

## What it provides

- **Cosmological observables** — background and distances ($\Lambda$CDM, wCDM, kinematic reconstructions), matter power spectra (linear and Halofit), large-scale structure, galaxy clusters, weak lensing, SNIa, BAO, CMB, and cross-correlations.

- **Statistical framework** — best-fit estimation, Fisher forecasts, and MCMC samplers (ESMCMC and the APES sampler), parallelized with MPI and OpenMP.

- **API documentation** — [NumCosmo](https://numcosmo.readthedocs.io/en/latest/reference/numcosmo/) (cosmology) and [NumCosmoMath](https://numcosmo.readthedocs.io/en/latest/reference/numcosmo-math/) (math foundation).

## Get Started

To get started with NumCosmo, follow these steps:

1. **Installation:** Clone the repository and install the necessary dependencies as outlined in the [installation guide](https://numcosmo.readthedocs.io/en/latest/install.html).

2. **Examples and Tutorials:** Work through the [examples and tutorials](https://numcosmo.readthedocs.io/en/latest/) on the website to see how to drive NumCosmo for cosmological calculations and statistical analysis.

3. **Contribute:** NumCosmo is an open-source project, and we welcome contributions from the community.

## Citation

If you find NumCosmo useful for your research or project, please consider citing it:
```
@Misc{Vitenti2012c,
  author        = {S. D. P. Vitenti and M. Penna-Lima},
  title         = {Numerical Cosmology -- {NumCosmo}},
  year          = {2014},
  howpublished  = {Astrophysics Source Code Library},
  eprint        = {1408.013},
  url           = {https://github.com/NumCosmo/NumCosmo},
  adsurl        = {http://adsabs.harvard.edu/abs/2014ascl.soft08013D},
  archiveprefix = {ascl},
  primaryclass  = {ascl}
}
```

## License

NumCosmo is released under the GNU General Public License v3.0. See the [license file](https://github.com/NumCosmo/NumCosmo/blob/master/COPYING) for more details.

