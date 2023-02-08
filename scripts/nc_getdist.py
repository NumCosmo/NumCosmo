import gi

gi.require_version('NumCosmo', '1.0')
gi.require_version('NumCosmoMath', '1.0')

from gi.repository import NumCosmoMath as Ncm

import numpy as np

from getdist import plots, MCSamples
import getdist

def mcat_to_mcsamples(mcat: Ncm.MSetCatalog, name: str) -> MCSamples:
    nchains = mcat.nchains()
    rows = np.array ([mcat.peek_row (i).dup_array () for i in range (0, mcat.len ())])
    params = [mcat.col_symb (i) for i in range (mcat.ncols ())]
    m2lnL = mcat.get_m2lnp_var ()
    posterior = 0.5 * rows[:,m2lnL]
    rows   = np.delete (rows,   m2lnL, 1)
    params = np.delete (params, m2lnL, 0)
    
    split_chains = [rows[n::nchains] for n in range(nchains)]
    split_posterior = [posterior[n::nchains] for n in range(nchains)]
    
    mcsample = MCSamples(samples=split_chains, loglikes=split_posterior, names = params, labels = params, label = name)

    return mcsample, rows, posterior
