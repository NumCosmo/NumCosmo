#!/usr/bin/env python

from astropy.io import fits
from plot2Ddist import plot2Ddist
import matplotlib.pyplot as pl
import numpy as np
import sys
import os.path
from typing import List, Dict, Any

sys.argv.pop(0)

files: List[Dict[str, Any]] = []
col_name_att: Dict[str, Any] = {}

while len(sys.argv) > 1:

    fitsfile = sys.argv.pop(0)
    padding_cols = int(sys.argv.pop(0))
    burning_in = int(sys.argv.pop(0))
    (prefix, ext) = os.path.splitext(fitsfile)

    print("# Opening: %s." % fitsfile)
    hdulist = fits.open(fitsfile)
    cols = hdulist[1].columns
    ncols = len(cols)
    col_names = []
    col_symbols = []
    col_data = []
    print(
        "# Number of columns %d padding columns %d burning in %d."
        % (ncols, padding_cols, burning_in)
    )

    for i in range(ncols):
        col = hdulist[1].data.field(i)
        col = col[burning_in:]
        col_data.append(col)
        col_name = cols.names[i]
        col_names.append(col_name)
        col_symbol = None

        if i < padding_cols:
            col_symbol = col_name
        else:
            col_symbol = r"$%s$" % (
                hdulist[1].header["FSYMB%d" % (i - padding_cols + 1)]
                )
        col_symbols.append(col_symbol)

        col_min = np.percentile(col, 0.0035)  # min (col)
        col_max = np.percentile(col, 99.9965)  # max (col)
        print(col_min, col_max)
        # interval = col_max - col_min
        # col_min = col_min - interval * 0.05
        # col_max = col_max + interval * 0.05

        if col_name in col_name_att:
            cur_min = col_name_att[col_name]["min"]
            cur_max = col_name_att[col_name]["max"]
            col_name_att[col_name] = {
                "min": min(col_min, cur_min),
                "max": max(col_max, cur_max),
            }
        else:
            col_name_att[col_name] = {"min": col_min, "max": col_max}

    files.append(
        {
            "fitsfile": fitsfile,
            "pad": padding_cols,
            "prefix": prefix,
            "ncols": ncols,
            "cols": col_data,
            "col_names": col_names,
            "col_symbols": col_symbols,
        }
    )
    hdulist.close()
    del hdulist
    del fitsfile
    del padding_cols
    del prefix
    del ext

for file in files:

    scatterstyle = {"color": "r", "alpha": 0.5}
    styleargs = {"color": "k", "scatterstyle": scatterstyle}
    padding_cols = file["pad"]
    ncols = file["ncols"]

    for i in range(padding_cols, ncols):
        for j in range(i + 1, ncols):
            x = file["cols"][i]
            y = file["cols"][j]
            xlim = [
                col_name_att[file["col_names"][i]]["min"],
                col_name_att[file["col_names"][i]]["max"],
            ]
            ylim = [
                col_name_att[file["col_names"][j]]["min"],
                col_name_att[file["col_names"][j]]["max"],
            ]

            labelx = file["col_symbols"][i]
            labely = file["col_symbols"][j]

            plot2Ddist(
                [x, y],
                plothists=True,
                labels=[labelx, labely],
                xlim=xlim,
                ylim=ylim,
                **styleargs
            )

            figfile = "%s_%s_%s.pdf" % (
                file["prefix"],
                file["col_names"][i],
                file["col_names"][j],
            )
            print("#  Saving figure: %s." % figfile)
            pl.savefig(figfile)

    pl.close()
