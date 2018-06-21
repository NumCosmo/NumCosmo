# NumCosmo 

Numerical Cosmology Library

---

### NumCosmo Python example

In this example we build a model and a data object in Python and the use all NumCosmo
tools to analyze the data.

## Model object

The objective of the model is to describe a simple function $f(x) \in \mathbb{R}$ of some 
independent variable $x \in \mathbb{R}$. Below we implement a model where $f(x) = a e^{\alpha x}$.

+++?code=examples/pydata_simple/py_sline_model.py&lang=python&title=SLine model
@[5-14](Importing NumCosmo and NumCosmoMath using GI.)
@[20-21](Creating a new model using the NcmModelBuilder object.)
@[28-29, 36-37](Adding parameters to the model.)
@[45-48](All necessary steps to bring out the new object from C to Python)
@[53](Defines a implementation of the class built above)
@[58](If we need a property that is not a model parameter that's how we define it.)
@[63-64](In the constructor we need to chain-up the father's constructor.)
@[69-70](Defines the only method of this object, the function $f(x) = e^{x m}b.$)
@[75](As the last step we register the new object in the type system.)

+++

## Data object

The data object will describe the observation error distribution. For simplicity we describe 
the variable $\delta f_i = f^o_i - f(x_i),$ with a multivariate Gaussian distribution with zero 
mean and covariance $C_{ij}$ randomly generated the standard deviation of $\delta_ f_i$ will 
be draw from a uniform distribution between $[0.5, 2]$.

---



