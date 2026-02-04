#!/usr/bin/env python3
"""
Plot spherical Bessel ODE solutions for Gaussian and Rational RHS functions.

This script reads parameters from the truth tables and plots the solution
for a specific l and k value.
"""

import json
import gzip
from pathlib import Path
from typing import Tuple
import numpy as np
import matplotlib.pyplot as plt
import typer
from typing_extensions import Annotated
from numcosmo_py import Ncm

# Initialize NumCosmo
Ncm.cfg_init()

app = typer.Typer()


def load_truth_table(func_type: str) -> dict:
    """Load truth table from data directory.

    :param func_type: Either "gaussian" or "rational"
    :return: Truth table data
    """
    # Fix the filename prefix
    prefix = "gauss" if func_type == "gaussian" else "rational"
    filename = f"{prefix}_jl_500.json.gz"
    truth_table_path = Path(Ncm.cfg_get_data_filename(f"truth_tables/{filename}", True))

    with gzip.open(truth_table_path, "rt") as f:
        return json.load(f)


def compute_rhs_values(
    func_type: str, x_physical: np.ndarray, k: float, center: float, std: float
) -> np.ndarray:
    """Compute RHS function values.

    :param func_type: Either "gaussian" or "rational"
    :param x_physical: Physical domain x values
    :param k: Wavenumber
    :param center: Function center parameter
    :param std: Function width parameter
    :return: RHS function values
    """
    if func_type == "gaussian":
        inv_std2 = 1.0 / (std * std)
        return np.array(
            [np.exp(-((x / k - center) ** 2) * inv_std2 / 2.0) for x in x_physical]
        )
    else:  # rational
        inv_std = 1.0 / std
        return np.array(
            [
                ((x / k) ** 2) / ((1.0 + ((x / k - center) * inv_std) ** 2) ** 3)
                for x in x_physical
            ]
        )


def solve_ode(
    func_type: str,
    ell: int,
    k_index: int,
    N: int,
    truth_table: dict,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Ncm.Vector, dict]:
    """Solve ODE and compute solution and derivative values.

    :param func_type: Either "gaussian" or "rational"
    :param ell: Angular momentum quantum number
    :param k_index: Index into k values array
    :param N: Number of Chebyshev nodes
    :param truth_table: Truth table data
    :return: Tuple (x_physical, sol_vals, sol_deriv_vals, rhs_vals, solution, params)
    """
    # Extract parameters
    center = truth_table["center"]
    std = truth_table["std"]
    lb = truth_table["lower-bound"]
    ub = truth_table["upper-bound"]
    kvals = truth_table["kvals"]
    k = kvals[k_index]

    # Create solver
    solver = Ncm.SBesselOdeSolver.new(ell, lb, ub)

    # Get RHS
    if func_type == "gaussian":
        T_rhs_func = solver.get_gaussian_rhs(center, std, k, N)
    else:  # rational
        T_rhs_func = solver.get_rational_rhs(center, std, k, N)

    rhs_func = T_rhs_func.dup()
    solver.chebT_to_gegenbauer_lambda2(T_rhs_func, rhs_func)

    rhs = Ncm.Vector.new(rhs_func.len() + 2)
    rhs.set(0, 0.0)  # Boundary condition at lb
    rhs.set(1, 0.0)  # Boundary condition at ub
    rhs.memcpy2(rhs_func, 2, 0, rhs_func.len())

    # Update solver interval
    solver.set_interval(lb * k, ub * k)

    # Solve
    solution = solver.solve(rhs)

    # Evaluate solution
    xi = np.linspace(-1, 1, 5000)
    a_scaled, b_scaled = lb * k, ub * k
    mid = 0.5 * (a_scaled + b_scaled)
    half_h = 0.5 * (b_scaled - a_scaled)
    x_physical = mid + half_h * xi

    sol_vals = np.array(
        [Ncm.SBesselOdeSolver.chebyshev_eval(solution, xi_val) for xi_val in xi]
    )

    # Evaluate derivative (in computational domain)
    sol_deriv_xi = np.array(
        [Ncm.SBesselOdeSolver.chebyshev_deriv(solution, xi_val) for xi_val in xi]
    )
    # Transform to physical domain: dy/dx = (dy/dxi) / half_h
    sol_deriv_vals = sol_deriv_xi / half_h

    # Compute boundary terms for Green's identity
    y_prime_a = Ncm.SBesselOdeSolver.chebyshev_deriv(solution, -1.0) / half_h
    y_prime_b = Ncm.SBesselOdeSolver.chebyshev_deriv(solution, 1.0) / half_h
    j_l_a = Ncm.sf_sbessel(ell, a_scaled)
    j_l_b = Ncm.sf_sbessel(ell, b_scaled)
    boundary_term_a = a_scaled * a_scaled * j_l_a * y_prime_a
    boundary_term_b = b_scaled * b_scaled * j_l_b * y_prime_b
    green_identity = boundary_term_b - boundary_term_a

    # Compute RHS function values
    rhs_vals = compute_rhs_values(func_type, x_physical, k, center, std)

    params = {
        "k": k,
        "center": center,
        "std": std,
        "lb": lb,
        "ub": ub,
        "boundary_term_a": boundary_term_a,
        "boundary_term_b": boundary_term_b,
        "green_identity": green_identity,
    }
    return x_physical, sol_vals, sol_deriv_vals, rhs_vals, solution, params


def add_plots_to_axes(
    ax1,
    ax2,
    ax3,
    x_physical: np.ndarray,
    sol_vals: np.ndarray,
    sol_deriv_vals: np.ndarray,
    rhs_vals: np.ndarray,
    func_type: str,
    ell: int,
    k: float,
):
    """Add solution, derivative, and RHS plots to given axes.

    :param ax1: Axes for solution plot
    :param ax2: Axes for derivative plot
    :param ax3: Axes for RHS plot
    :param x_physical: Physical domain x values
    :param sol_vals: Solution values
    :param sol_deriv_vals: Solution derivative values
    :param rhs_vals: RHS function values
    :param func_type: Either "gaussian" or "rational"
    :param ell: Angular momentum quantum number
    :param k: Wavenumber
    """
    # Plot solution
    ax1.plot(x_physical, np.abs(sol_vals), linewidth=2)
    ax1.set_xlabel("x (physical domain)", fontsize=11)
    ax1.set_ylabel("Solution y(x)", fontsize=11)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_title(
        f"{func_type.capitalize()} - Solution (l={ell}, k={k:.3f})",
        fontsize=12,
        fontweight="bold",
    )
    ax1.grid(True, alpha=0.3)

    # Plot derivative
    ax2.plot(x_physical, np.abs(sol_deriv_vals), "g-", linewidth=2)
    ax2.set_xlabel("x (physical domain)", fontsize=11)
    ax2.set_ylabel("dy/dx", fontsize=11)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_title(
        f"{func_type.capitalize()} - Derivative", fontsize=12, fontweight="bold"
    )
    ax2.grid(True, alpha=0.3)

    # Plot RHS
    ax3.plot(x_physical, rhs_vals, "r-", linewidth=2)
    ax3.set_xlabel("x (physical domain)", fontsize=11)
    ax3.set_ylabel("f(x)", fontsize=11)
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.set_title(
        f"{func_type.capitalize()} - RHS Function", fontsize=12, fontweight="bold"
    )
    ax3.grid(True, alpha=0.3, which="both")


def plot_solution(func_type: str, ell: int, k_index: int = 0, N: int = 128):
    """Plot the ODE solution for a given function type.

    :param func_type: Either "gaussian" or "rational"
    :param ell: Angular momentum quantum number
    :param k_index: Index into the k values array from truth table
    :param N: Number of Chebyshev nodes for the solution
    """
    # Load truth table
    truth_table = load_truth_table(func_type)

    # Solve ODE
    x_physical, sol_vals, sol_deriv_vals, rhs_vals, _, params = solve_ode(
        func_type, ell, k_index, N, truth_table
    )

    print(f"Plotting {func_type} solution:")
    print(f"  l      = {ell}")
    print(f"  k      = {params['k']}")
    print(f"  center = {params['center']}")
    print(f"  std    = {params['std']}")
    print(f"  bounds = [{params['lb']}, {params['ub']}]")
    print(f"  N      = {N}")
    print("\nGreen's identity boundary terms:")
    print(f"  a²j_l(a)y'(a) = {params['boundary_term_a']:.10e}")
    print(f"  b²j_l(b)y'(b) = {params['boundary_term_b']:.10e}")
    print(f"  Difference    = {params['green_identity']:.10e}")

    # Create figure
    _, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))

    # Add plots
    add_plots_to_axes(
        ax1,
        ax2,
        ax3,
        x_physical,
        sol_vals,
        sol_deriv_vals,
        rhs_vals,
        func_type,
        ell,
        params["k"],
    )

    plt.tight_layout()

    # Save and show
    output_file = f"sbessel_solution_{func_type}_l{ell}_k{k_index}.png"
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    print(f"\nPlot saved to: {output_file}")

    plt.show()


def compare_both_functions(ell: int, k_index: int = 0, N: int = 128):
    """Plot both Gaussian and Rational solutions side by side.

    :param ell: Angular momentum quantum number
    :param k_index: Index into the k values array from truth table
    :param N: Number of Chebyshev nodes for the solution
    """
    _, axes = plt.subplots(2, 3, figsize=(18, 10))

    for i, func_type in enumerate(["gaussian", "rational"]):
        # Load truth table
        truth_table = load_truth_table(func_type)

        # Solve ODE
        x_physical, sol_vals, sol_deriv_vals, rhs_vals, _, params = solve_ode(
            func_type, ell, k_index, N, truth_table
        )

        # Add plots to the appropriate row
        add_plots_to_axes(
            axes[i, 0],
            axes[i, 1],
            axes[i, 2],
            x_physical,
            sol_vals,
            sol_deriv_vals,
            rhs_vals,
            func_type,
            ell,
            params["k"],
        )

    plt.tight_layout()

    # Save figure
    output_file = f"sbessel_comparison_l{ell}_k{k_index}.png"
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    print(f"\nComparison plot saved to: {output_file}")

    plt.show()


@app.command()
def main(
    func_type: Annotated[
        str,
        typer.Option(help="Function type to plot"),
    ] = "both",
    l: Annotated[
        int,
        typer.Option(help="Angular momentum quantum number"),
    ] = 10,
    k_index: Annotated[
        int,
        typer.Option(help="Index into k values array"),
    ] = 0,
    n: Annotated[
        int,
        typer.Option(help="Number of Chebyshev nodes"),
    ] = 128,
):
    """Plot spherical Bessel ODE solutions from truth tables."""
    if func_type not in ["gaussian", "rational", "both"]:
        typer.echo("Error: func_type must be 'gaussian', 'rational', or 'both'")
        raise typer.Exit(1)

    if func_type == "both":
        compare_both_functions(l, k_index, n)
    else:
        plot_solution(func_type, l, k_index, n)


if __name__ == "__main__":
    app()
