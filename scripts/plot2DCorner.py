import pylab
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pl
import matplotlib.patches
import matplotlib.colors as mcolors
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import plot2Ddist

import scipy.stats

def plotCorner (xs, bins=50, labels=None, label_size=9, range=None, fig=None, max_n_ticks=5, 
                plotscatter=True, use_math_text=False, thin=1, plotcontours=True,
                labelcontours=True, contourKDEthin=1, contourNGrid=100, 
                contourFractions=[0.6827, 0.9545, 0.9973],
                contourstyle={}, **styleArgs):

    """
    This script was created based on corner.py and plot2Ddist.py.
    
    Parameters
        ----------
        xs : array_like (nsamples, ndim)
            The samples. This should be a 1- or 2-dimensional array. For a 1-D
            array this results in a simple histogram. For a 2-D array, the zeroth
            axis is the list of samples and the next axis are the dimensions of
            the space.
            
        labels : iterable (ndim,) (optional)
            A list of names for the dimensions. If a ``xs`` is a
            ``pandas.DataFrame``, labels will default to column names.
            
       range : iterable (ndim,) (optional)
         A list where each element is either a length 2 tuple containing
         lower and upper bounds or a float in range (0., 1.)
         giving the fraction of samples to include in bounds, e.g.,
         [(0.,10.), (1.,5), 0.999, etc.].
         If a fraction, the bounds are chosen to be equal-tailed.
         
    """
    quantiles = []
    title_kwargs = dict()
    label_kwargs = dict()

    # Try filling in labels from pandas.DataFrame columns.
    if labels is None:
        try:
            labels = xs.columns
        except AttributeError:
            pass
    
    # Deal with 1D sample lists.
    xs = np.atleast_1d(xs)
    if len(xs.shape) == 1:
        xs = np.atleast_2d(xs)
    else:
        assert len(xs.shape) == 2, "The input sample array must be 1- or 2-D."
        xs = xs.T
    assert xs.shape[0] <= xs.shape[1], "I don't believe that you want more " \
                                       "dimensions than samples!"
                                       
    if range is None:
        if "extents" in hist2d_kwargs:
            logging.warn("Deprecated keyword argument 'extents'. "
                         "Use 'range' instead.")
            range = hist2d_kwargs.pop("extents")
        else:
            range = [[x.min(), x.max()] for x in xs]
            # Check for parameters that never change.
            m = np.array([e[0] == e[1] for e in range], dtype=bool)
            if np.any(m):
                raise ValueError(("It looks like the parameter(s) in "
                                  "column(s) {0} have no dynamic range. "
                                  "Please provide a `range` argument.")
                                 .format(", ".join(map(
                                     "{0}".format, np.arange(len(m))[m]))))

    else:
        # If any of the extents are percentiles, convert them to ranges.
        # Also make sure it's a normal list.
        range = list(range)
        for i, _ in enumerate(range):
            try:
                emin, emax = range[i]
            except TypeError:
                q = [0.5 - 0.5*range[i], 0.5 + 0.5*range[i]]
                range[i] = quantile(xs[i], q) #, weights=weights)

    if len(range) != xs.shape[0]:
        raise ValueError("Dimension mismatch between samples and range")

    # Parse the bin specifications.
    try:
        bins = [int(bins) for _ in range]
    except TypeError:
        if len(bins) != len(range):
            raise ValueError("Dimension mismatch between bins and range")                                                   

    # Some magic numbers for pretty axis layout.
    K = len(xs)
    factor = 2.0           # size of one side of one panel
    lbdim = 0.5 * factor   # size of left/bottom margin
    trdim = 0.2 * factor   # size of top/right margin
    whspace = 0.05         # w/hspace size
    plotdim = factor * K + factor * (K - 1.) * whspace
    dim = lbdim + plotdim + trdim

    # Create a new figure if one wasn't provided.
    if fig is None:
        fig, axes = pl.subplots(K, K, figsize=(dim, dim))
    else:
        try:
            axes = np.array(fig.axes).reshape((K, K))
        except:
            raise ValueError("Provided figure has {0} axes, but data has "
                             "dimensions K={1}".format(len(fig.axes), K))

    # Format the figure.
    lb = lbdim / dim
    tr = (lbdim + plotdim) / dim
    fig.subplots_adjust(left=lb, bottom=lb, right=tr, top=tr,
                        wspace=whspace, hspace=whspace)

    # Set up the default histogram keywords.
    hist_kwargs = dict()
    hist_kwargs["color"] = hist_kwargs.get("color", 'k')
    hist_kwargs["histtype"] = hist_kwargs.get("histtype", "step")
    
    trimto = len(xs[0])

    for i, x in enumerate(xs):
        # Deal with masked arrays.
        if hasattr(x, "compressed"):
            x = x.compressed()

        if np.shape(xs)[0] == 1:
            ax = axes
        else:
            ax = axes[i, i]
        # Plot the histograms.
        n, _, _ = ax.hist(x, bins=bins[i], range=range[i], **hist_kwargs) # weights=weights

        # Set up the axes.
        ax.set_xlim(range[i])
        ax.set_ylim(0, 1.1 * np.max(n))
        ax.set_yticklabels([])
        ax.xaxis.set_major_locator(MaxNLocator(max_n_ticks, prune="lower"))

        if i < K - 1:
            ax.set_xticklabels([])
        else:
            [l.set_rotation(45) for l in ax.get_xticklabels()]
            if labels is not None:
                ax.set_xlabel(labels[i], fontsize = label_size) #**label_kwargs
                ax.xaxis.set_label_coords(0.5, -0.3)

            # use MathText for axes ticks
            ax.xaxis.set_major_formatter(
                ScalarFormatter(useMathText=use_math_text))
                
        for j, y in enumerate(xs):
            if np.shape(xs)[0] == 1:
                ax = axes
            else:
                ax = axes[i, j]
            if j > i:
                ax.set_frame_on(False)
                ax.set_xticks([])
                ax.set_yticks([])
                continue
            elif j == i:
                continue
        
            xlim = [min(range[i]), max(range[i])]
            ylim = [min(range[j]), max(range[j])]

            # Plot 2D scatter of variables.
            if plotscatter:
                c = mcolors.ColorConverter().to_rgb
                rvb = plot2Ddist.make_colormap([c('#FFFFFF'), 0.01, c('#C8C8C8'), 0.20, c('#969696'), 0.40, c('#646464'), 0.60, c('#323232'), 0.80, c('#000000')])
                ax.hist2d (y, x, bins=80, cmap=rvb, range = [ylim, xlim])
                print ("Plotscatter: (%d, %d) done." % (i, j))

            if plotcontours:
                xkde = x[-trimto::contourKDEthin]
                ykde = y[-trimto::contourKDEthin]
                # Inspired by Abraham Flaxman's https://gist.github.com/626689
                style = {'linewidths':2.0, 'alpha':0.75, 'colors':'r',
                         #'cmap':matplotlib.cm.Greys,
                         'zorder':10}
                style.update(styleArgs)
                style.update(contourstyle)
                if 'color' in style:
                    style['colors'] = style['color']
                gkde = scipy.stats.gaussian_kde([ykde,xkde])
                #xgrid, ygrid = numpy.mgrid[xlim[0]:xlim[1]:contourNGrid * 1j,
                #                           ylim[0]:ylim[1]:contourNGrid * 1j]
                ygrid, xgrid = np.mgrid[ylim[0]:ylim[1]:contourNGrid * 1j,
                                           xlim[0]:xlim[1]:contourNGrid * 1j]
                zvals = np.array(gkde.evaluate([ygrid.flatten(),
                                                   xgrid.flatten()])
                                    ).reshape(ygrid.shape)
                contours = plot2Ddist.contour_enclosing(y, x, contourFractions, 
                                             ygrid, xgrid, zvals, 
                                             ax, **style)      
                print ("Plot countour: (%d, %d) done." % (i, j))

            #if plotcontours and labelcontours:
                #plot2Ddist.frac_label_contours(y, x, contours)

            ax.xaxis.set_major_locator(MaxNLocator(max_n_ticks, prune="lower"))
            ax.yaxis.set_major_locator(MaxNLocator(max_n_ticks, prune="lower"))                

            if i < K - 1:
                ax.set_xticklabels([])
            else:
                [l.set_rotation(45) for l in ax.get_xticklabels()]
                if labels is not None:
                    ax.set_xlabel(labels[j], fontsize = label_size) #**label_kwargs
                    ax.xaxis.set_label_coords(0.5, -0.3)

                # use MathText for axes ticks
                ax.xaxis.set_major_formatter(
                    ScalarFormatter(useMathText=use_math_text))

            if j > 0:
                ax.set_yticklabels([])
            else:
                [l.set_rotation(45) for l in ax.get_yticklabels()]
                if labels is not None:
                    ax.set_ylabel(labels[i], fontsize = label_size) #**label_kwargs
                    ax.yaxis.set_label_coords(-0.3, 0.5)

               # use MathText for axes ticks
                ax.yaxis.set_major_formatter(
                    ScalarFormatter(useMathText=use_math_text))
                                                    
    return fig
