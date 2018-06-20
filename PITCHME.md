# NumCosmo 

Numerical Cosmology Library

---

### NumCosmo Python example

---?code=examples/pydata_simple/py_sline_model.py&lang=python&title=SLine model

@[5-14](Importing NumCosmo and NumCosmoMath using GI.)
@[20-21](Creating a new model using the NcmModelBuilder object.)
@[28-29, 36-37](Adding parameters to the model.)
@[45-48](All necessary steps to bring out the new object from C to Python)
@[53](Defines a implementation of the class built above)
@[58](If we need a property that is not a model parameter that's how we define it.)
@[63-64](In the constructor we need to chain-up the father's constructor.)
@[69-70](Defines the only method of this object, the function $$f(x) = e^{x m}b.$$)
@[75](As the last step we register the new object in the type system.)
---


