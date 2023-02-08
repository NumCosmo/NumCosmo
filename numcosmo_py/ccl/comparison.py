try:
    import gi

    gi.require_version("NumCosmo", "1.0")
    gi.require_version("NumCosmoMath", "1.0")
except:
    pass

from gi.repository import GObject
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

import pyccl
import math
import numpy as np
import pylab as plt


def compare_ccl_nc_func(
    x,
    y_ccl,
    y_nc,
    x_name="x",
    y_name="func",
    subplots_pars={"figsize": (12, 6)},
    xscale="linear",
    yscale="log",
):

    ccl_name, nc_name = "%s^\\mathrm{ccl}" % y_name, "%s^\\mathrm{nc}" % y_name

    x = np.array(x)
    y_ccl = np.array(y_ccl)
    y_nc = np.array(y_nc)
    diff = np.zeros_like(y_ccl)

    non_zind = np.where(y_ccl != 0.0)[0]
    zind = np.where(y_ccl == 0.0)[0]
    diff[non_zind] = y_nc[non_zind] / y_ccl[non_zind] - 1.0
    diff[zind] = y_nc[zind] - y_ccl[zind]
    print(
        "[%10s]: rel diff min: %e\trel diff max: %e"
        % (y_name, min(abs(diff)), max(abs(diff)))
    )

    fig, axs = plt.subplots(2, sharex=True, **subplots_pars)
    fig.subplots_adjust(hspace=0)

    axs[0].plot(x, y_ccl, label="ccl", lw=3)
    axs[0].plot(x, y_nc, label="nc")
    axs[1].plot(x, np.abs(diff), c="r")
    axs[1].set_xscale(xscale)
    axs[1].set_yscale(yscale)

    axs[0].legend()
    axs[0].set_ylabel("$%s$" % y_name)
    axs[1].set_xlabel("$%s$" % x_name)
    axs[1].set_ylabel("$%s/%s-1$" % (nc_name, ccl_name))

    return fig, axs
