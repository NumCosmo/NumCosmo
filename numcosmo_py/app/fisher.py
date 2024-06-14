#
# fisher.py
#
# Wed Feb 14 18:59:34 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# fisher.py
# Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""NumCosmo APP subcommand to compute the Fisher matrix and related quantities."""


import dataclasses
from pathlib import Path
from typing import Optional, Annotated, cast

import typer

from .. import Ncm
from ..sampling import FisherType
from .loading import LoadExperiment
from .run_fit import RunCommonOptions


@dataclasses.dataclass(kw_only=True)
class ComputeTheoryVector(LoadExperiment):
    """Compute theory vectory for a given experiment."""

    def __post_init__(self) -> None:
        """Compute theory vector for a given experiment."""
        super().__post_init__()

        dset: Ncm.Dataset = self.likelihood.peek_dataset()
        if not dset.has_mean_vector():
            raise RuntimeError("mean vector computation not supported by this dataset.")

        theory_vector = Ncm.Vector.new(dset.get_n())
        dset.mean_vector(self.mset, theory_vector)
        if self.output is not None:
            self.output_dict.add("theory-vector", theory_vector)
        else:
            self.console.print(theory_vector.dup_array())

        self.end_experiment()


@dataclasses.dataclass(kw_only=True)
class RunFisher(RunCommonOptions):
    """Compute the Fisher matrix of the model to the data."""

    fisher_type: Annotated[
        FisherType,
        typer.Option(
            help="Type of Fisher matrix to compute.",
        ),
    ] = FisherType.OBSERVED

    def __post_init__(self) -> None:
        """Compute the Fisher matrix of the model to the data."""
        super().__post_init__()
        if self.fisher_type == FisherType.OBSERVED:
            self.fit.obs_fisher()
            self.fit.log_covar()
        elif self.fisher_type == FisherType.EXPECTED:
            self.fit.fisher()
            self.fit.log_covar()
        else:
            raise RuntimeError(f"Invalid Fisher type {self.fisher_type}.")

        if self.output is not None:
            self.output_dict.add("model-set", self.fit.peek_mset())
            self.output_dict.add("covariance", self.fit.get_covar())

        self.end_experiment()


@dataclasses.dataclass(kw_only=True)
class RunFisherBias(RunCommonOptions):
    """Computes the Fisher matrix of the model to the data and the bias."""

    theory_vector: Annotated[
        Optional[Path],
        typer.Option(
            help="Path to the theory vector file to compute the bias relative to."
        ),
    ] = None

    def __post_init__(self) -> None:
        """Compute the Fisher matrix of the model to the data and the bias."""
        super().__post_init__()

        if self.product_file and self.theory_vector is not None:
            raise RuntimeError(
                "The theory vector option is incompatible with the product file option."
            )

        if self.theory_vector is None:
            if self.product_file:
                theory_vector: Ncm.Vector = cast(
                    Ncm.Vector, self.output_dict.get("theory-vector")
                )
                if theory_vector is None:
                    raise RuntimeError(
                        "No theory vector found in the product file, "
                        "cannot compute bias. Use the --theory-vector option or "
                        "a product file with a theory vector."
                    )
            else:
                raise RuntimeError("No theory vector file given.")
        else:
            ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
            theory_vector_dict = ser.dict_str_from_yaml_file(
                self.theory_vector.absolute().as_posix()
            )
            theory_vector = cast(Ncm.Vector, theory_vector_dict.get("theory-vector"))
            if not isinstance(theory_vector, Ncm.Vector):
                raise RuntimeError("Invalid theory vector file.")

        dset = self.likelihood.peek_dataset()
        if not dset.has_mean_vector():
            raise RuntimeError(
                "mean vector computation not supported by this dataset, "
                "cannot compute bias."
            )

        if theory_vector.len() != dset.get_n():
            raise RuntimeError(
                "Theory vector and dataset have different sizes, "
                "cannot compute bias."
            )

        delta_theta = self.fit.fisher_bias(theory_vector)

        if self.output is not None:
            self.output_dict.add("covariance", self.fit.get_covar())
            self.output_dict.add("delta-theta", delta_theta)

        self.console.print(delta_theta.dup_array())

        self.end_experiment()
