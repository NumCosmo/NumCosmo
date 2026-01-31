#!/usr/bin/env python
# Simple test for NcmSBesselOdeSolver

import pytest
import numpy as np
from numpy.testing import assert_allclose
from numcosmo_py import Ncm


class TestSBesselOdeSolver:
    """Tests for NcmSBesselOdeSolver."""

    @pytest.fixture
    def solver(self) -> Ncm.SBesselOdeSolver:
        """Create an ODE solver with l=1."""
        return Ncm.SBesselOdeSolver.new(1, -1.0, 1.0)

    def test_creation(self, solver: Ncm.SBesselOdeSolver):
        """Test solver creation."""
        assert solver is not None
        assert solver.get_l() == 1
        a, b = solver.get_interval()
        assert a == -1.0
        assert b == 1.0

    def test_solve(self, solver: Ncm.SBesselOdeSolver):
        """Test solving ODE."""
        N = 20
        rhs = Ncm.Vector.new(N)
        for i in range(N):
            rhs.set(i, 1.0 / (i + 1) ** 2)

        solution = solver.solve(rhs)
        assert solution is not None
        assert solution.len() > 0

    def test_get_operator_matrix(self, solver: Ncm.SBesselOdeSolver):
        """Test matrix extraction."""
        mat = solver.get_operator_matrix(10)
        assert mat is not None
        assert mat.nrows() == 10
        assert mat.ncols() > 0
