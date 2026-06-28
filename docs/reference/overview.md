Title: NumCosmo Overview

# NumCosmo API Overview

This is the generated API reference for NumCosmo. The library is split into two
GObject namespaces; pick the one matching what you need:

1. **NumCosmoMath** — [API Reference](../numcosmo-math/index.html)
   : The math foundation, independent of cosmology: vectors and matrices,
     splines, integration and differentiation, special functions, the
     model/MSet parameter system, the fitting and MCMC tools (including the APES
     sampler), FFTLog, and HEALPix. Class names are prefixed `Ncm`.

2. **NumCosmo** — [API Reference](../numcosmo/index.html)
   : The cosmology layer built on top: background models and distances,
     perturbations, the matter power spectrum, large-scale structure, galaxy
     clusters, weak lensing, and the observational likelihoods (SNIa, BAO, CMB,
     cross-correlations). Class names are prefixed `Nc`.

## How the two relate

`Nc` cosmology models are concrete realizations of the generic `Ncm` model and
fitting machinery. A cosmological model is an `NcmModel` with a defined
parameter set; it plugs into the same `NcmFit`, `NcmFitMC`, and `NcmFitESMCMC`
tools used for any statistical model. This means the prediction code and the
constraint code share one object.

## Where to start

- Conceptual and physics background, tutorials, and worked examples live on the
  [NumCosmo website](https://numcosmo.readthedocs.io/en/latest/).
- For the mathematical derivations behind specific classes, follow the
  *theoretical background* links in the class documentation (for example,
  `NcmCSQ1D` links to the CSQ1D Formalism page).
- Downstream code should include only the umbrella headers `numcosmo.h` and
  `numcosmo-math.h`, never individual subdirectory headers.
