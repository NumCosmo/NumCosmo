#!/usr/bin/env python

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

import math
import argparse
import numpy as np

from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

from scipy.stats import chi2
from chainconsumer import ChainConsumer

parser = argparse.ArgumentParser(description='Process mset catalogs')

parser.add_argument('-C', '--catalog', metavar='file.fits', 
                    help='catalog fits file', 
                    action='append', required = True)

parser.add_argument('-B', '--burnin', metavar='N', 
                    help='catalog burnin', type=int, 
                    action='append')

parser.add_argument('--kde',
                    help='whether to use kde interpolation',  
                    action='store_true')

parser.add_argument('--col', type=int, nargs="*",
                    help='Columns to include')

parser.add_argument('--truth', type=float, nargs="*",
                    help='Columns to include')

parser.add_argument('--sigma', type=int, nargs="+", default = [1, 2],
                    help='Sigmas to compute the confidence regions')

parser.add_argument('--out', default = "corner.pdf",
                    help='Output filename')

parser.add_argument('--mode', choices=['corner', 'walks'], default='corner')

args = parser.parse_args()

Ncm.cfg_init ()

c = ChainConsumer()

for cat in args.catalog:

    bin = 0
    if args.burnin and (len (args.burnin) > 0):
        bin = args.burnin.pop (0)
        
    print (f"# Adding {cat} with burnin {bin}")

    mcat = Ncm.MSetCatalog.new_from_file_ro (cat, bin)
    nwalkers = mcat.nchains ()
    
    m2lnL = mcat.get_m2lnp_var ()

    rows = np.array ([mcat.peek_row (i).dup_array () for i in range (mcat.len ())])
    params = ["$" + mcat.col_symb (i) + "$" for i in range (mcat.ncols ())]
        
    posterior = -0.5 * rows[:,m2lnL]
    
    rows   = np.delete (rows,   m2lnL, 1)
    params = np.delete (params, m2lnL, 0)

    if args.col:
        assert max (args.col) < mcat.ncols ()
        indices = np.array (args.col)

        rows   = rows[:,indices]
        params = params[indices]

    #c.add_chain(rows, posterior = posterior, parameters=list(params), plot_point = True, name = cat.replace ("_", "-"))
    c.add_chain(rows, posterior = posterior, parameters=list(params), name = cat.replace ("_", "-"))

c.configure (kde = args.kde, label_font_size=8, sigma2d=False, sigmas = args.sigma, spacing = 0.0, tick_font_size=8)

plot_args = {}

if args.truth is not None:
    plot_args['truth'] = args.truth

if args.mode == "corner":
    fig = c.plotter.plot(**plot_args)
elif args.mode == "walks":
    fig = c.plotter.plot_walks(**plot_args, convolve=100)
else:
    assert False

fig.savefig(args.out, bbox_inches="tight")

