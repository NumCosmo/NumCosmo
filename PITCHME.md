# NumCosmo 

Numerical Cosmology Library

---

## NumCosmo Python example

In this example we build a model and a data object in Python and the use all NumCosmo
tools to analyze the data.

### Model object

The objective of the model is to describe a simple function $f(x) \in \mathbb{R}$ of some 
independent variable $x \in \mathbb{R}$. Below we implement a model where $f(x) = a e^{\alpha x}$.

+++?code=examples/pydata_simple/py_sline_model.py&lang=python&title=SLine model
@[5-14](Importing NumCosmo and NumCosmoMath using GI.)
@[20-21](Creating a new model using the *NcmModelBuilder* object.)
@[28-29, 36-37](Adding parameters to the model.)
@[45-48](All necessary steps to bring out the new object from C to Python)
@[53](Defines a implementation of the class built above)
@[58](If we need a property that is not a model parameter that's how we define it.)
@[63-64](In the constructor we need to chain-up the father's constructor.)
@[69-70](Defines the only method of this object, the function $f(x) = e^{x \alpha}a.$)
@[75](As the last step we register the new object in the type system.)

+++

### Data object

The data object will describe the observation error distribution. For simplicity we describe 
the variable $ \delta f_i = f^\mathrm{obs}_i - f(x_i) $, with a multivariate Gaussian distribution with zero 
mean and covariance `C` randomly generated. The standard deviation of $ \delta f_i $ will 
be draw from a uniform distribution between $[0.5, 2]$.

+++?code=examples/pydata_simple/py_sline_gauss.py&lang=python&title=SLine data

@[7-16](Importing NumCosmo and NumCosmoMath using GI.)
@[18](Imports our model.)
@[23](Define a new data class, child of *NcmDataGaussCov*.)
@[27](Our distribution have an one dimensional vector of independent variables $x$ which we define here.)
@[32-40](In the initialization process, we first chain-up to the father's initialization.)
@[52-53, 58-59](Here we implement the necessary virtual functions of *NcmDataGaussCov*.)
@[68-69, 76-77](In our simple likelihood it is not necessary to have a begin or a prepare function, usually needed when precalculations are necessary.)
@[84-92](Here we inform the *NcmDataGaussCov* object how to compute $f(x_i)$.)
@[94-99](We define a method to create a random covariance matrix, with the correlation factor 15 and the standard deviations ranging from 0.5 to 2.)
@[105](As the last step we register the new object in the type system.)

+++

### Generating data randomly

We can generate a sample for our likelihood using the resample method already 
implemented in the *NcmDataGaussCov* object, in the script below we explain how
to do it

+++?code=examples/pydata_simple/example_create_data.py&lang=python&title=SLine data generate

@[9-22](Loads everything.)
@[28](In a executable script using NumCosmo must call Ncm.cfg_init before any other NumCosmo functions.)
@[33](Creates a new pseudo-random number generator with the default algorithm and seed 123.)
@[39-41](Here we define our fiducial model that will be used to generate the data.)
@[53-55](The data object is created with length = $50$, the covariance will be a $50\times50$ matrix. Our independent vector is just a uniform grid between $[0, 10]$.)
@[61-64](We need a model set *NcmMSet* to transport our model and the free parameters map.)
@[73-76](Finally, we create a serialization object, use the resample method to create the data and save it to disk.)

+++

### Fitting the generated data

In the next script we exemplify how to fit and compute the Fisher matrix, both observed and expected.
The observed fisher matrix is simple proportional to $\partial_i\partial_j(-2\ln L)$
and must be computed at the maximum of `L` to be a good estimate of the expected 
$\langle\partial_i\partial_j(-2\ln L)\rangle$.

+++?code=examples/pydata_simple/example_fit.py&lang=python&title=SLine model fit

@[9-22](Loads everything.)
@[28](In a executable script using NumCosmo must call Ncm.cfg_init before any other NumCosmo functions.)
@[34-36](Here we define model and the initial guess point.)
@[42-45](We need a model set *NcmMSet* to transport our model and the free parameters map.)
@[53, 58](We create a new serialization object and load the data from the datafile.)
@[63-65, 72-73](The data set *NcmDataset* contains our data and is used to build the likelihood *NcmLikelihood* and the fit *NcmFit* objects.)
@[79, 84](We then run the minimization process and log the results.)
@[89-90, 95-96](Now, at the best-fit we compute the Fisher matrix both ways.)

+++

### Results from fit

```bash
#----------------------------------------------------------------------------------
# Model fitting. Interating using:
#  - solver:            NLOpt:ln-neldermead
#  - differentiation:   Numerical differentiantion (forward)
#................
#  Minimum found with precision: |df|/f =  1.00000e-08 and |dx| =  1.00000e-05
#  Elapsed time: 00 days, 00:00:00.0134530
#  iteration            [000073]
#  function evaluations [000074]
#  gradient evaluations [000000]
#  degrees of freedom   [000048]
#  m2lnL     =     48.2024604329839 (      48.20246 )
#  Fit parameters:
#     1.00686253235027     0.47436330251462
#----------------------------------------------------------------------------------
# Data used:
#   - py_sline_gauss+PySLineGauss
#----------------------------------------------------------------------------------
# Model[00000]:
#   - NcPySLineModel : A simple python example model
#----------------------------------------------------------------------------------
# Model parameters
#   - alpha[00]:  1.00686253235027    [FREE]
#   -     a[01]:  0.47436330251462    [FREE]
# Observed Fisher Matrix
#----------------------------------------------------------------------------------
# NcmMSet parameters covariance matrix
#                                                 -------------------------------
# alpha[00000:00] =  1.007       +/-  0.01927     |  1           | -0.8932      |
#     a[00000:01] =  0.4744      +/-  0.05232     | -0.8932      |  1           |
#                                                 -------------------------------
# Expected Fisher matrix
#----------------------------------------------------------------------------------
# NcmMSet parameters covariance matrix
#                                                 -------------------------------
# alpha[00000:00] =  1.007       +/-  0.01839     |  1           | -0.8844      |
#     a[00000:01] =  0.4744      +/-  0.05041     | -0.8844      |  1           |
#                                                 -------------------------------
```
@[3-4](The algorithm used for the minimization process and differentiation (if needed).)
@[8](Number of steps to find the minimum.)
@[12, 14](The best-fit value for $-2\ln L$ and the parameters at the minimum.)
@[17, 20](The data and models which were used in the run.)
@[23-24](Here we show all paremeters, fixed or not.)
@[29-30](The observed Fisher matrix.)
@[36-37](The expected Fisher matrix.)

+++

### Monte Carlo study

We now perform a Monte Carlo study: we first get the best-fit and then use it as a
fiducial model. From the fiducial model we compute several realizations of the data
and for each one we compute the best fit.

+++?code=examples/pydata_simple/example_mc.py&lang=python&title=SLine model MC

@[104-106](The script is the same for the fit up to this point, here we create a *NcmFitMC* object.)
@[123-125, 131-132](Then we just need to run it, and print out the results.)

+++

### Monte Carlo results

```bash
# NcmMSetCatalog: Current mean:   48.059       1.0068       0.47597    
# NcmMSetCatalog: Current msd:    0.092678     0.00017523   0.00047581 
# NcmMSetCatalog: Current sd:     9.7263       0.01839      0.049936   
# NcmMSetCatalog: Current var:    94.601       0.0003382    0.0024936  
# NcmMSetCatalog: Current tau:    1            1            1          
# Task:NcmFitMC, completed: 11014 of 11014, elapsed time: 00:11:42.5975
# Task:NcmFitMC, mean time: 00:00:00.0638 +/- 00:00:00.0006
# Task:NcmFitMC, time left: 00:00:00.0000 +/- 00:00:00.0000
# Task:NcmFitMC, current time:        Fri Jun 22 2018, 09:34:48
# Task:NcmFitMC, estimated to end at: Fri Jun 22 2018, 09:34:48 +/- 00:00:00.0000
# NcmFitMC: Largest relative error 1.000000e-03 attained: 9.996788e-04
#----------------------------------------------------------------------------------
# NcmMSet parameters covariance matrix
#                                                 -------------------------------
# alpha[00000:00] =  1.007       +/-  0.01839     |  1           | -0.886       |
#     a[00000:01] =  0.476       +/-  0.04994     | -0.886       |  1           |
#                                                 -------------------------------
```
@[1](The means for $-2\ln L$, $\alpha$ and $a$.)
@[2](The standard deviation on the means.)
@[3](The standard deviation on the parameters.)
@[15-16](The covariance matrix from the sample of best-fits. Remember that the fiducial model had $\alpha = 1.007$ and $a = 0.4744$.)

+++

### Monte Carlo results

![MC](examples/pydata_simple/example_mc_out.png)

---

# Questions?

---



















