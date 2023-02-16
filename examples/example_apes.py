#!/usr/bin/env python
#
# example_apes.py
#
# Wed Feb 8 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_apes.py
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

"""Example of using the APES MCMC sampler on test posteriors."""

import typer
from numcosmo_py import Ncm
from numcosmo_py.experiments.rosenbrock import run_rosenbrock_mcmc
from numcosmo_py.experiments.gaussmix2d import run_gaussmix2d_mcmc
from numcosmo_py.experiments.funnel import run_funnel_mcmc


Ncm.cfg_init()
app = typer.Typer()

app.command()(run_rosenbrock_mcmc)
app.command()(run_gaussmix2d_mcmc)
app.command()(run_funnel_mcmc)

if __name__ == "__main__":
    app()
