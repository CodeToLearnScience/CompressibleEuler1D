#!/usr/bin/env python3
"""
Run all flux/order combinations and generate comparison plots/tables.
"""

import subprocess
import sys
from pathlib import Path
from dataclasses import dataclass
import shutil

import numpy as np
import matplotlib.pyplot as plt

# Add scripts to path
sys.path.insert(0, str(Path(__file__).parent))
from validate import read_numerical_csv, read_analytical_dat, compute_errors, interpolate_to_grid

WORKSPACE = Path(__file__).parent.parent
DATA_DIR = WORKSPACE / "data"
RESULTS_DIR = WORKSPACE / "results"
PLOTS_DIR = WORKSPACE / "plots"
BUILD_DIR = WORKSPACE / "build"

FLUX_SCHEMES = ["llf", "rusanov", "hll", "hllc"]
ORDERS = [1, 2]
TEST_CASES = list(range(1, 13))


def create_modified_toml(base_toml: Path, flux: str, order: int) -> Path:
    """Create a modified TOML file with specified flux and order."""
    content = base_toml.read_text()
    
    # Replace flux
    import re
    content = re.sub(r'flux\s*=\s*"[^"]+"', f'flux = "{flux}"', content)
    content = re.sub(r'order\s*=\s*\d+', f'order = {order}', content)
    
    output_path = RESULTS_DIR / f"{base_toml.stem}_{flux}_order{order}.toml"
    output_path.write_text(content)
    return output_path


def run_simulation(toml_path: Path, output_name: str) -> bool:
    """Run simulation with given TOML config."""
    exe = BUILD_DIR / "euler1d"
    if not exe.exists():
        print(f"ERROR: Executable not found: {exe}")
        return False
    
    result = subprocess.run(
        [str(exe), str(toml_path)],
        capture_output=True,
        text=True,
        cwd=WORKSPACE
    )
    
    # Move output files
    for ext in [".csv", ".vtk"]:
        # Output is based on test_name in toml which is test_caseN
        import re
        match = re.search(r'test_case(\d+)', toml_path.stem)
        if match:
            case_num = match.group(1)
            src = WORKSPACE / f"test_case{case_num}{ext}"
            if src.exists():
                dst = RESULTS_DIR / f"{output_name}{ext}"
                shutil.move(str(src), str(dst))
    
    return result.returncode == 0


def compute_all_errors():
    """Compute errors for all result files."""
    errors = {}
    
    for flux in FLUX_SCHEMES:
        errors[flux] = {}
        for order in ORDERS:
            errors[flux][order] = {}
            for case in TEST_CASES:
                result_file = RESULTS_DIR / f"test_case{case}_{flux}_order{order}.csv"
                ref_file = DATA_DIR / f"analytical_ref_test_case{case}.dat"
                
                if not result_file.exists():
                    continue
                if not ref_file.exists():
                    continue
                
                try:
                    num_sol = read_numerical_csv(result_file)
                    ana_sol = read_analytical_dat(ref_file)
                    
                    # Interpolate if needed
                    if len(num_sol.x) != len(ana_sol.x):
                        ana_sol = interpolate_to_grid(ana_sol, num_sol.x)
                    
                    dx = num_sol.x[1] - num_sol.x[0]
                    
                    errors[flux][order][case] = {
                        'rho': compute_errors(num_sol.rho, ana_sol.rho, dx),
                        'u': compute_errors(num_sol.u, ana_sol.u, dx),
                        'p': compute_errors(num_sol.p, ana_sol.p, dx),
                    }
                except Exception as e:
                    print(f"Error computing for {flux}/order{order}/case{case}: {e}")
    
    return errors


def generate_comparison_plot(case: int, output_dir: Path):
    """Generate comparison plot for a single test case showing all flux/order combos."""
    ref_file = DATA_DIR / f"analytical_ref_test_case{case}.dat"
    if not ref_file.exists():
        return
    
    ana_sol = read_analytical_dat(ref_file)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"Test Case {case}: Density Comparison", fontsize=14)
    
    colors = {'llf': 'blue', 'rusanov': 'green', 'hll': 'orange', 'hllc': 'red'}
    linestyles = {1: '--', 2: '-'}
    
    for idx, (ax, flux) in enumerate(zip(axes.flat, FLUX_SCHEMES)):
        ax.plot(ana_sol.x, ana_sol.rho, 'k-', linewidth=2, label='Analytical')
        
        for order in ORDERS:
            result_file = RESULTS_DIR / f"test_case{case}_{flux}_order{order}.csv"
            if result_file.exists():
                num_sol = read_numerical_csv(result_file)
                ax.plot(num_sol.x, num_sol.rho, 
                       linestyle=linestyles[order], 
                       color=colors[flux],
                       linewidth=1.2,
                       alpha=0.8 if order == 1 else 1.0,
                       label=f'{flux.upper()} Order {order}')
        
        ax.set_xlabel('x')
        ax.set_ylabel('Density')
        ax.set_title(f'{flux.upper()} Flux')
        ax.legend(loc='best', fontsize=8)
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = output_dir / f"comparison_case{case}_density.png"
    plt.savefig(output_file, dpi=150)
    plt.close()
    print(f"Saved: {output_file}")


def generate_order_comparison_plot(case: int, flux: str, output_dir: Path):
    """Generate plot comparing 1st vs 2nd order for specific flux."""
    ref_file = DATA_DIR / f"analytical_ref_test_case{case}.dat"
    if not ref_file.exists():
        return
    
    ana_sol = read_analytical_dat(ref_file)
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f"Test Case {case}: {flux.upper()} Flux - 1st vs 2nd Order", fontsize=14)
    
    variables = [
        ('rho', 'Density'),
        ('u', 'Velocity'),
        ('p', 'Pressure'),
    ]
    
    for ax, (var, label) in zip(axes.flat[:3], variables):
        ana_data = getattr(ana_sol, var)
        ax.plot(ana_sol.x, ana_data, 'k-', linewidth=2, label='Analytical')
        
        for order in ORDERS:
            result_file = RESULTS_DIR / f"test_case{case}_{flux}_order{order}.csv"
            if result_file.exists():
                num_sol = read_numerical_csv(result_file)
                num_data = getattr(num_sol, var)
                style = '--' if order == 1 else '-'
                ax.plot(num_sol.x, num_data, style, linewidth=1.2, 
                       label=f'Order {order}')
        
        ax.set_xlabel('x')
        ax.set_ylabel(label)
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    # Error plot
    ax = axes[1, 1]
    for order in ORDERS:
        result_file = RESULTS_DIR / f"test_case{case}_{flux}_order{order}.csv"
        if result_file.exists():
            num_sol = read_numerical_csv(result_file)
            ana_interp = np.interp(num_sol.x, ana_sol.x, ana_sol.rho)
            error = np.abs(num_sol.rho - ana_interp)
            ax.semilogy(num_sol.x, error + 1e-16, label=f'Order {order}')
    
    ax.set_xlabel('x')
    ax.set_ylabel('|ρ_num - ρ_ana|')
    ax.set_title('Density Error')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = output_dir / f"case{case}_{flux}_order_comparison.png"
    plt.savefig(output_file, dpi=150)
    plt.close()
    print(f"Saved: {output_file}")


def generate_flux_comparison_plot(case: int, order: int, output_dir: Path):
    """Generate plot comparing all fluxes for given order."""
    ref_file = DATA_DIR / f"analytical_ref_test_case{case}.dat"
    if not ref_file.exists():
        return
    
    ana_sol = read_analytical_dat(ref_file)
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f"Test Case {case}: All Flux Schemes - Order {order}", fontsize=14)
    
    colors = {'llf': 'blue', 'rusanov': 'green', 'hll': 'orange', 'hllc': 'red'}
    
    variables = [
        ('rho', 'Density'),
        ('u', 'Velocity'),
        ('p', 'Pressure'),
    ]
    
    for ax, (var, label) in zip(axes.flat[:3], variables):
        ana_data = getattr(ana_sol, var)
        ax.plot(ana_sol.x, ana_data, 'k-', linewidth=2, label='Analytical')
        
        for flux in FLUX_SCHEMES:
            result_file = RESULTS_DIR / f"test_case{case}_{flux}_order{order}.csv"
            if result_file.exists():
                num_sol = read_numerical_csv(result_file)
                num_data = getattr(num_sol, var)
                ax.plot(num_sol.x, num_data, '-', color=colors[flux], 
                       linewidth=1.2, label=flux.upper())
        
        ax.set_xlabel('x')
        ax.set_ylabel(label)
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    # Error comparison
    ax = axes[1, 1]
    for flux in FLUX_SCHEMES:
        result_file = RESULTS_DIR / f"test_case{case}_{flux}_order{order}.csv"
        if result_file.exists():
            num_sol = read_numerical_csv(result_file)
            ana_interp = np.interp(num_sol.x, ana_sol.x, ana_sol.rho)
            error = np.abs(num_sol.rho - ana_interp)
            ax.semilogy(num_sol.x, error + 1e-16, color=colors[flux], label=flux.upper())
    
    ax.set_xlabel('x')
    ax.set_ylabel('|ρ_num - ρ_ana|')
    ax.set_title('Density Error')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = output_dir / f"case{case}_flux_comparison_order{order}.png"
    plt.savefig(output_file, dpi=150)
    plt.close()
    print(f"Saved: {output_file}")


def format_error_table_markdown(errors: dict, cases: list = None) -> str:
    """Format errors as markdown table."""
    if cases is None:
        cases = TEST_CASES
    
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
    print("=" * 70)
    print("Running all flux/order combinations")
    print("=" * 70)
    
    # Run simulations
    for flux in FLUX_SCHEMES:
        for order in ORDERS:
            print(f"\n--- {flux.upper()} Order {order} ---")
            for case in TEST_CASES:
                base_toml = DATA_DIR / f"test_case{case}.toml"
                if not base_toml.exists():
                    continue
                
                toml_path = create_modified_toml(base_toml, flux, order)
                output_name = f"test_case{case}_{flux}_order{order}"
                
                print(f"  Running test_case{case}...", end=" ", flush=True)
                if run_simulation(toml_path, output_name):
                    print("OK")
                else:
                    print("FAILED")
    
    # Compute errors
    print("\n" + "=" * 70)
    print("Computing errors")
    print("=" * 70)
    errors = compute_all_errors()
    
    # Generate plots
    print("\n" + "=" * 70)
    print("Generating plots")
    print("=" * 70)
    
    # Key test cases for documentation
    key_cases = [1, 2, 3, 4, 5]
    
    for case in key_cases:
        generate_comparison_plot(case, PLOTS_DIR)
        for flux in FLUX_SCHEMES:
            generate_order_comparison_plot(case, flux, PLOTS_DIR)
        for order in ORDERS:
            generate_flux_comparison_plot(case, order, PLOTS_DIR)
    
    # Generate error table
    print("\n" + "=" * 70)
    print("Error Table (Markdown)")
    print("=" * 70)
    table = format_error_table_markdown(errors, key_cases)
    print(table)
    
    # Save error table
    (RESULTS_DIR / "error_tables.md").write_text(table)
    print(f"\nSaved error tables to: {RESULTS_DIR / 'error_tables.md'}")


if __name__ == "__main__":
    main()
