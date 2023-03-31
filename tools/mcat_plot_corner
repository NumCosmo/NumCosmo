#
# plot_corner.py
#
# Wed Feb 8 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# plot_corner.py
# Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
#
# numcosmo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# numcosmo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""NumCosmoPy corner plot utilities."""

import argparse
import numpy as np

from chainconsumer import ChainConsumer
import gi

gi.require_version("NumCosmo", "1.0")
gi.require_version("NumCosmoMath", "1.0")

# pylint:disable-next=wrong-import-position
from gi.repository import NumCosmoMath as Ncm

parser = argparse.ArgumentParser(description="Process mset catalogs")

parser.add_argument(
    "-C",
    "--catalog",
    metavar="file.fits",
    help="catalog fits file",
    action="append",
    required=True,
)

parser.add_argument(
    "-B", "--burnin", metavar="N", help="catalog burnin", type=int, action="append"
)

parser.add_argument(
    "--kde", help="whether to use kde interpolation", action="store_true"
)

parser.add_argument("--col", type=int, nargs="*", help="Columns to include")

parser.add_argument("--truth", type=float, nargs="*", help="Columns to include")

parser.add_argument(
    "--sigma",
    type=int,
    nargs="+",
    default=[1, 2],
    help="Sigmas to compute the confidence regions",
)

parser.add_argument("--out", default="corner.pdf", help="Output filename")

parser.add_argument("--mode", choices=["corner", "walks"], default="corner")

args = parser.parse_args()

Ncm.cfg_init()

c = ChainConsumer()

for cat in args.catalog:

    burnin = 0
    if args.burnin and (len(args.burnin) > 0):
        burnin = args.burnin.pop(0)

    print(f"# Adding {cat} with burnin {burnin}")

    mcat = Ncm.MSetCatalog.new_from_file_ro(cat, burnin)
    nwalkers = mcat.nchains()

    m2lnL = mcat.get_m2lnp_var()

    rows = np.array([mcat.peek_row(i).dup_array() for i in range(mcat.len())])
    params = ["$" + mcat.col_symb(i) + "$" for i in range(mcat.ncols())]

    posterior = -0.5 * rows[:, m2lnL]

    rows = np.delete(rows, m2lnL, 1)
    params = list(np.delete(params, m2lnL, 0))
    print(params)

    if args.col:
        assert max(args.col) < mcat.ncols()
        indices = np.array(args.col)

        rows = rows[:, indices]
        params = params[indices]

    c.add_chain(
        rows, posterior=posterior, parameters=list(params), name=cat.replace("_", "-")
    )

c.configure(
    kde=args.kde,
    label_font_size=8,
    sigma2d=False,
    sigmas=args.sigma,
    spacing=0.0,
    tick_font_size=8,
)

plot_args = {}

if args.truth is not None:
    plot_args["truth"] = args.truth

if args.mode == "corner":
    fig = c.plotter.plot(**plot_args)
elif args.mode == "walks":
    fig = c.plotter.plot_walks(**plot_args, convolve=100)
else:
    assert False

fig.savefig(args.out, bbox_inches="tight")
