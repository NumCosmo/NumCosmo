
#!/usr/bin/env python 
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import os.path

try:
    import gi
    gi.require_version('NumCosmo', '1.0')  
    gi.require_version('NumCosmoMath', '1.0')
except:
    pass
    
from gi.repository import GObject
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm
__name___ = "NcContext"

Ncm.cfg_init ()
Ncm.cfg_set_log_handler (lambda msg: sys.stdout.write (msg) and sys.stdout.flush ())

dim = 5
np.random.seed(seed=123)
p = np.random.random_sample((dim,))
print (p)

rng = Ncm.RNG.seeded_new (None, 123)

fmodel = Ncm.ModelMVND.new (dim)
fdata = Ncm.DataGaussCovMVND.new_full (dim, 0.1, 0.4, 10.0, -1.0, 1.0, rng)

fdata.props.use_norma = True

mset = Ncm.MSet.new_array ([fmodel])
mset.param_set_all_ftype (Ncm.ParamType.FREE)
mset.prepare_fparam_map ()

mset.fparams_set_array (p)
mset.pretty_log ()

y_a = []
x_a = []
xn_a = []

#interp = Ncm.StatsDistNdKDEStudentt.new (len (p), Ncm.StatsDistNdCV.NONE, 3.0)
interp = Ncm.StatsDistNdVBKStudentt.new (len (p), Ncm.StatsDistNdCV.NONE, 3.0)
nps = 500

for a in range (nps):
    v, N = fdata.gen (mset, None, None, rng)
    interp.add_obs (v)
    y_a.append (fdata.m2lnL_val (mset))
    x_a.append (v.dup_array ())

for a in range (nps):
    v, N = fdata.gen (mset, None, None, rng)
    xn_a.append (v.dup_array ())

y_a = np.array (y_a)
x_a = np.array (x_a)
xn_a = np.array (xn_a)

fdata.y.set_array (p)

fdata.m2lnL_val (mset)
interp.set_cv_type (Ncm.StatsDistNdCV.NONE)
interp.set_split_frac (1.0)
interp.prepare_interp(Ncm.Vector.new_array (y_a))
       
yi_a = []
ya_a = []
yi_b = []

for x in xn_a:
    mset.fparams_set_array (x)
    fdata.prepare (mset)
    yi_a.append (interp.eval_m2lnp (Ncm.Vector.new_array (x)))
    ya_a.append (fdata.m2lnL_val (mset))

yi_a = np.array (yi_a)
ya_a = np.array (ya_a)

m = -np.mean (yi_a - ya_a)
diff = (yi_a + m) / ya_a - 1.0