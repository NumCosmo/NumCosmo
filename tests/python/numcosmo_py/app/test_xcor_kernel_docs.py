#!/usr/bin/env python
#
# test_xcor_kernel_docs.py
#
# Thu Sep 10 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_xcor_kernel_docs.py
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Tests for the kernel documentation the xcor CLI generates."""

import pytest
from typer.testing import CliRunner

import matplotlib

matplotlib.use("Agg")

# flake8: noqa: E402
# pylint: disable=wrong-import-position

from pydantic import BaseModel

from numcosmo_py import Ncm
from numcosmo_py.app import app
from numcosmo_py.app.xcor.common import XcorKernelCommon
from numcosmo_py.app.xcor.kernels import (
    KERNEL_CONFIG_REGISTRY,
    _field_descriptions,
    kernel_parameter_summary,
)

pytestmark = pytest.mark.app
runner = CliRunner()

Ncm.cfg_init()

REGISTRY = list(KERNEL_CONFIG_REGISTRY.items())
NAMES = [name for name, _ in REGISTRY]


def test_kernel_list_summarises_every_type() -> None:
    """The summary names every registered kernel and points at the detail view."""
    result = runner.invoke(app, ["xcor", "kernel", "list"])

    assert result.exit_code == 0, result.output
    for name in NAMES:
        assert name in result.output
    assert "numcosmo xcor kernel list <kernel_type>" in result.output


@pytest.mark.parametrize("kernel_name", NAMES)
def test_kernel_list_summary_names_every_parameter(kernel_name: str) -> None:
    """No parameter is missing from the summary.

    The summary used to be a hand-written string per kernel and had drifted:
    bessel_deriv was absent from five of the seven radial shapes, so there was
    no way to find out from the command line that they take it.
    """
    result = runner.invoke(app, ["xcor", "kernel", "list"])

    assert result.exit_code == 0, result.output
    summary = kernel_parameter_summary(KERNEL_CONFIG_REGISTRY[kernel_name])
    for field in KERNEL_CONFIG_REGISTRY[kernel_name].model_fields:
        assert f"{field}=" in summary


@pytest.mark.parametrize(("kernel_name", "config_class"), REGISTRY)
def test_kernel_list_documents_one_type(
    kernel_name: str, config_class: type[BaseModel]
) -> None:
    """Each kernel documents its own parameters, defaults and library object."""
    result = runner.invoke(app, ["xcor", "kernel", "list", kernel_name])

    assert result.exit_code == 0, result.output
    assert kernel_name in result.output
    assert config_class.nc_type.__name__ in result.output  # type: ignore[attr-defined]
    for field in config_class.model_fields:
        assert field in result.output


@pytest.mark.parametrize(("kernel_name", "config_class"), REGISTRY)
def test_every_parameter_is_documented(
    kernel_name: str, config_class: type[BaseModel]
) -> None:
    """Every field carries a description, which the class docstring supplies.

    A field added without an ``:ivar:`` entry reaches the command line as a
    blank Description cell, which is what this guards against.
    """
    documented = _field_descriptions(config_class)

    for field in config_class.model_fields:
        assert documented.get(field), f"{kernel_name}.{field} has no description"


def test_kernel_list_rejects_an_unknown_type() -> None:
    """An unknown name is refused, with the known ones named."""
    result = runner.invoke(app, ["xcor", "kernel", "list", "radial-banana"])

    assert result.exit_code != 0
    assert "radial-gauss" in result.output


def test_documented_library_object_is_the_one_built() -> None:
    """Each config's nc_type is what the command actually constructs.

    The documentation is generated from nc_type, so it is only as good as that
    pointing at the right class -- which nothing but this checks for the six
    kernels whose builders name their class themselves.
    """
    command = XcorKernelCommon(kernel=NAMES)

    assert len(command.kernels) == len(REGISTRY)
    for (label, kernel_obj), (kernel_name, config_class) in zip(
        command.kernels, REGISTRY
    ):
        nc_type = config_class.nc_type  # type: ignore[attr-defined]
        assert isinstance(kernel_obj, nc_type), kernel_name
        assert label.startswith(config_class.label)  # type: ignore[attr-defined]
