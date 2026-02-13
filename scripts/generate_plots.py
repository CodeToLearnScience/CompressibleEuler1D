#!/usr/bin/env python3
"""
Generate plots and error tables from existing results.
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Add scripts to path
sys.path.insert(0, str(Path(__file__).parent))
from validate import read_numerical_csv, read_analytical_dat, compute_errors, interpolate_to_grid

WORKSPACE = Path(__file__).parent.parent
DATA_DIR = WORKSPACE / "data"
RESULTS_DIR = WORKSPACE / "results"
PLOTS_DIR = WORKSPACE / "plots"

FLUX_SCHEMES = ["llf", "hll", "hllc"]
ORDERS = [1, 2]


def compute_all_errors(cases: list) -> dict:
    """Compute errors for specified test cases."""
    errors = {}
    
    for flux in FLUX_SCHEMES:
        errors[flux] = {}
        for order in ORDERS:
            errors[flux][order] = {}
            for case in cases:
                result_file = RESULTS_DIR / f"test_case{case}_{flux}_order{order}.csv"
                ref_file = DATA_DIR / f"analytical_ref_test_case{case}.dat"
                
                if not result_file.exists() or not ref_file.exists():
                    continue
                
                try:
                    num_sol = read_numerical_csv(result_file)
                    ana_sol = read_analytical_dat(ref_file)
                    
                    if len(num_sol.x) != len(ana_sol.x):
                        ana_sol = interpolate_to_grid(ana_sol, num_sol.x)
                    
                    dx = num_sol.x[1] - num_sol.x[0]
                    
                    errors[flux][order][case] = {
                        'rho': compute_errors(num_sol.rho, ana_sol.rho, dx),
                        'u': compute_errors(num_sol.u, ana_sol.u, dx),
                        'p': compute_errors(num_sol.p, ana_sol.p, dx),
                    }
                except Exception as e:
                    print(f"Error for {flux}/order{order}/case{case}: {e}")
    
    return errors


def generate_sod_comparison(output_dir: Path):
    """Generate Sod shock tube (case 1) comparison plot."""
    case = 1
    ref_file = DATA_DIR / f"analytical_ref_test_case{case}.dat"
    ana_sol = read_analytical_dat(ref_file)
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle("Sod Shock Tube (Test Case 1): Flux Scheme Comparison", fontsize=14)
    
    colors = {'llf': 'C0', 'hll': 'C1', 'hllc': 'C2'}
    
    variables = [('rho', 'Density'), ('u', 'Velocity'), ('p', 'Pressure')]
    
    # First row: 1st order
    for ax, (var, label) in zip(axes[0], variables):
        ax.plot(ana_sol.x, getattr(ana_sol, var), 'k-', lw=2, label='Analytical')
        for flux in FLUX_SCHEMES:
            f = RESULTS_DIR / f"test_case{case}_{flux}_order1.csv"
            if f.exists():
                sol = read_numerical_csv(f)
                ax.plot(sol.x, getattr(sol, var), '--', color=colors[flux], lw=1.2, label=flux.upper())
        ax.set_xlabel('x')
        ax.set_ylabel(label)
        ax.set_title(f'{label} (1st Order)')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    
    # Second row: 2nd order
    for ax, (var, label) in zip(axes[1], variables):
        ax.plot(ana_sol.x, getattr(ana_sol, var), 'k-', lw=2, label='Analytical')
        for flux in FLUX_SCHEMES:
            f = RESULTS_DIR / f"test_case{case}_{flux}_order2.csv"
            if f.exists():
                sol = read_numerical_csv(f)
                ax.plot(sol.x, getattr(sol, var), '-', color=colors[flux], lw=1.2, label=flux.upper())
        ax.set_xlabel('x')
        ax.set_ylabel(label)
        ax.set_title(f'{label} (2nd Order)')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    out = output_dir / "sod_shock_tube_comparison.png"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"Saved: {out}")


def generate_order_comparison(case: int, output_dir: Path):
    """Compare 1st vs 2nd order for all fluxes."""
    ref_file = DATA_DIR / f"analytical_ref_test_case{case}.dat"
    if not ref_file.exists():
        return
    ana_sol = read_analytical_dat(ref_file)
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(f"Test Case {case}: 1st Order vs 2nd Order (Density)", fontsize=14)
    
    for ax, flux in zip(axes, FLUX_SCHEMES):
        ax.plot(ana_sol.x, ana_sol.rho, 'k-', lw=2, label='Analytical')
        
        for order in ORDERS:
            f = RESULTS_DIR / f"test_case{case}_{flux}_order{order}.csv"
            if f.exists():
                sol = read_numerical_csv(f)
                style = '--' if order == 1 else '-'
                ax.plot(sol.x, sol.rho, style, lw=1.2, label=f'Order {order}')
        
        ax.set_xlabel('x')
        ax.set_ylabel('Density')
        ax.set_title(f'{flux.upper()} Flux')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    out = output_dir / f"case{case}_order_comparison.png"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"Saved: {out}")


def generate_flux_comparison(case: int, order: int, output_dir: Path):
    """Compare all fluxes for given order."""
    ref_file = DATA_DIR / f"analytical_ref_test_case{case}.dat"
    if not ref_file.exists():
        return
    ana_sol = read_analytical_dat(ref_file)
    
    colors = {'llf': 'C0', 'hll': 'C1', 'hllc': 'C2'}
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(f"Test Case {case}: Flux Comparison (Order {order})", fontsize=14)
    
    variables = [('rho', 'Density'), ('u', 'Velocity'), ('p', 'Pressure')]
    
    for ax, (var, label) in zip(axes, variables):
        ax.plot(ana_sol.x, getattr(ana_sol, var), 'k-', lw=2, label='Analytical')
        
        for flux in FLUX_SCHEMES:
            f = RESULTS_DIR / f"test_case{case}_{flux}_order{order}.csv"
            if f.exists():
                sol = read_numerical_csv(f)
                ax.plot(sol.x, getattr(sol, var), '-', color=colors[flux], lw=1.2, label=flux.upper())
        
        ax.set_xlabel('x')
        ax.set_ylabel(label)
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    out = output_dir / f"case{case}_flux_comparison_order{order}.png"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"Saved: {out}")


def format_error_table(errors: dict, cases: list) -> str:
    """Generate markdown error table."""
    lines = []
    
    for flux in FLUX_SCHEMES:
        lines.append(f"\n### {flux.upper()} Flux\n")
        lines.append("| Case | Order | ρ L1 | ρ L2 | ρ L∞ | u L1 | p L1 |")
        lines.append("|------|-------|------|------|------|------|------|")
        
        for case in cases:
            for order in ORDERS:
                if flux in errors and order in errors[flux] and case in errors[flux][order]:
                    e = errors[flux][order][case]
                    lines.append(
                        f"| {case} | {order} | "
                        f"{e['rho'].L1:.2e} | {e['rho'].L2:.2e} | {e['rho'].Linf:.2e} | "
                        f"{e['u'].L1:.2e} | {e['p'].L1:.2e} |"
                    )
    
    return "\n".join(lines)


def main():
    PLOTS_DIR.mkdir(exist_ok=True)
    
    cases = list(range(1, 13))  # All 12 test cases
    
    print("Generating plots...")
    
    # Main Sod shock tube comparison
    generate_sod_comparison(PLOTS_DIR)
    
    # Order comparison for each case
    for case in cases:
        generate_order_comparison(case, PLOTS_DIR)
    
    # Flux comparison for each case/order
    for case in cases:
        for order in ORDERS:
            generate_flux_comparison(case, order, PLOTS_DIR)
    
    print("\nComputing errors...")
    errors = compute_all_errors(list(range(1, 13)))
    
    print("\nGenerating error table...")
    table = format_error_table(errors, cases)
    print(table)
    
    # Save tables
    (RESULTS_DIR / "error_tables.md").write_text(table)
    print(f"\nSaved: {RESULTS_DIR / 'error_tables.md'}")


if __name__ == "__main__":
    main()
