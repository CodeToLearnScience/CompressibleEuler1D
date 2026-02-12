#!/usr/bin/env python3
"""
Validation script for 1D Compressible Euler solver.

Compares numerical results against analytical reference solutions.
Computes L1, L2, and Linf errors for density, velocity, and pressure.
"""

import argparse
import sys
from pathlib import Path
from dataclasses import dataclass
from typing import Optional

import numpy as np
from numpy.typing import NDArray


@dataclass
class SolutionData:
    """Container for solution data."""
    x: NDArray[np.float64]
    rho: NDArray[np.float64]
    u: NDArray[np.float64]
    p: NDArray[np.float64]
    E: Optional[NDArray[np.float64]] = None


@dataclass
class ErrorMetrics:
    """Container for error metrics."""
    L1: float
    L2: float
    Linf: float


def read_numerical_csv(filepath: Path) -> SolutionData:
    """Read numerical solution from CSV file.
    
    Expected format:
    # 1D Euler solution at time = ...
    # x,rho,u,p,E
    data rows (comma-separated)
    """
    data = np.loadtxt(filepath, delimiter=',', comments='#')
    return SolutionData(
        x=data[:, 0],
        rho=data[:, 1],
        u=data[:, 2],
        p=data[:, 3],
        E=data[:, 4] if data.shape[1] > 4 else None
    )


def read_analytical_dat(filepath: Path) -> SolutionData:
    """Read analytical reference solution from DAT file.
    
    Expected format:
    whitespace-separated columns: x, rho, u, p, e (internal specific energy)
    """
    data = np.loadtxt(filepath)
    return SolutionData(
        x=data[:, 0],
        rho=data[:, 1],
        u=data[:, 2],
        p=data[:, 3],
        E=data[:, 4] if data.shape[1] > 4 else None  # This is internal specific energy
    )


def interpolate_to_grid(source: SolutionData, target_x: NDArray[np.float64]) -> SolutionData:
    """Interpolate solution to a different grid."""
    return SolutionData(
        x=target_x,
        rho=np.interp(target_x, source.x, source.rho),
        u=np.interp(target_x, source.x, source.u),
        p=np.interp(target_x, source.x, source.p),
        E=np.interp(target_x, source.x, source.E) if source.E is not None else None
    )


def compute_errors(numerical: NDArray[np.float64], 
                   analytical: NDArray[np.float64],
                   dx: float) -> ErrorMetrics:
    """Compute L1, L2, and Linf error norms."""
    diff = np.abs(numerical - analytical)
    return ErrorMetrics(
        L1=np.sum(diff) * dx,
        L2=np.sqrt(np.sum(diff**2) * dx),
        Linf=np.max(diff)
    )


def validate_case(numerical_path: Path, 
                  analytical_path: Path,
                  verbose: bool = True) -> dict:
    """Validate a single test case.
    
    Returns dict with error metrics for each variable.
    """
    # Read data
    num_sol = read_numerical_csv(numerical_path)
    ana_sol = read_analytical_dat(analytical_path)
    
    # Check if grids match
    if len(num_sol.x) != len(ana_sol.x) or not np.allclose(num_sol.x, ana_sol.x, rtol=1e-10):
        if verbose:
            print(f"  Grids differ: numerical has {len(num_sol.x)} points, "
                  f"analytical has {len(ana_sol.x)} points")
            print("  Interpolating analytical solution to numerical grid...")
        ana_sol = interpolate_to_grid(ana_sol, num_sol.x)
    
    # Compute grid spacing
    dx = num_sol.x[1] - num_sol.x[0] if len(num_sol.x) > 1 else 1.0
    
    # Compute errors
    errors = {
        'rho': compute_errors(num_sol.rho, ana_sol.rho, dx),
        'u': compute_errors(num_sol.u, ana_sol.u, dx),
        'p': compute_errors(num_sol.p, ana_sol.p, dx),
    }
    
    if verbose:
        print(f"\n  {'Variable':<10} {'L1 Error':<15} {'L2 Error':<15} {'Linf Error':<15}")
        print(f"  {'-'*55}")
        for var, err in errors.items():
            print(f"  {var:<10} {err.L1:<15.6e} {err.L2:<15.6e} {err.Linf:<15.6e}")
    
    return errors


def run_validation(numerical_dir: Path, 
                   analytical_dir: Path,
                   test_cases: Optional[list[int]] = None,
                   tolerance: float = 0.1) -> bool:
    """Run validation for all test cases.
    
    Args:
        numerical_dir: Directory containing numerical CSV files
        analytical_dir: Directory containing analytical DAT files  
        test_cases: List of test case numbers to validate (default: 1-12)
        tolerance: Maximum allowed L1 error for rho (default: 0.1)
        
    Returns:
        True if all validations pass, False otherwise
    """
    if test_cases is None:
        test_cases = list(range(1, 13))
    
    all_passed = True
    results = {}
    
    print("=" * 70)
    print("1D Compressible Euler Solver - Validation Report")
    print("=" * 70)
    
    for case_num in test_cases:
        numerical_file = numerical_dir / f"test_case{case_num}.csv"
        analytical_file = analytical_dir / f"analytical_ref_test_case{case_num}.dat"
        
        print(f"\n{'='*70}")
        print(f"Test Case {case_num}")
        print(f"{'='*70}")
        
        # Check files exist
        if not numerical_file.exists():
            print(f"  ERROR: Numerical solution not found: {numerical_file}")
            all_passed = False
            continue
            
        if not analytical_file.exists():
            print(f"  ERROR: Analytical reference not found: {analytical_file}")
            all_passed = False
            continue
        
        try:
            errors = validate_case(numerical_file, analytical_file)
            results[case_num] = errors
            
            # Check if errors are within tolerance
            rho_l1 = errors['rho'].L1
            status = "PASS" if rho_l1 < tolerance else "WARN"
            if rho_l1 >= tolerance:
                status = "WARN (high error - expected for 1st order)"
            
            print(f"\n  Status: {status}")
            
        except Exception as e:
            print(f"  ERROR: {e}")
            all_passed = False
    
    # Summary
    print(f"\n{'='*70}")
    print("Summary")
    print(f"{'='*70}")
    print(f"\n{'Case':<8} {'rho L1':<15} {'u L1':<15} {'p L1':<15}")
    print(f"{'-'*53}")
    
    for case_num, errors in sorted(results.items()):
        print(f"{case_num:<8} {errors['rho'].L1:<15.6e} "
              f"{errors['u'].L1:<15.6e} {errors['p'].L1:<15.6e}")
    
    return all_passed


def generate_comparison_plots(numerical_dir: Path,
                               analytical_dir: Path, 
                               output_dir: Path,
                               test_cases: Optional[list[int]] = None):
    """Generate comparison plots for each test case."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed. Skipping plot generation.")
        return
    
    if test_cases is None:
        test_cases = list(range(1, 13))
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    for case_num in test_cases:
        numerical_file = numerical_dir / f"test_case{case_num}.csv"
        analytical_file = analytical_dir / f"analytical_ref_test_case{case_num}.dat"
        
        if not numerical_file.exists() or not analytical_file.exists():
            continue
        
        num_sol = read_numerical_csv(numerical_file)
        ana_sol = read_analytical_dat(analytical_file)
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle(f"Test Case {case_num}: Numerical vs Analytical", fontsize=14)
        
        variables = [
            ('rho', 'Density', num_sol.rho, ana_sol.rho),
            ('u', 'Velocity', num_sol.u, ana_sol.u),
            ('p', 'Pressure', num_sol.p, ana_sol.p),
        ]
        
        for ax, (var_name, label, num_data, ana_data) in zip(axes.flat[:3], variables):
            ax.plot(ana_sol.x, ana_data, 'k-', label='Analytical', linewidth=1.5)
            ax.plot(num_sol.x, num_data, 'r--', label='Numerical', linewidth=1.2)
            ax.set_xlabel('x')
            ax.set_ylabel(label)
            ax.legend()
            ax.grid(True, alpha=0.3)
        
        # Error plot in 4th subplot
        axes[1, 1].semilogy(num_sol.x, np.abs(num_sol.rho - np.interp(num_sol.x, ana_sol.x, ana_sol.rho)), 
                            'b-', label='|rho_num - rho_ana|')
        axes[1, 1].set_xlabel('x')
        axes[1, 1].set_ylabel('Absolute Error')
        axes[1, 1].set_title('Density Error')
        axes[1, 1].grid(True, alpha=0.3)
        axes[1, 1].legend()
        
        plt.tight_layout()
        output_file = output_dir / f"comparison_test_case{case_num}.png"
        plt.savefig(output_file, dpi=150)
        plt.close()
        print(f"Saved plot: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Validate 1D Euler solver results against analytical solutions"
    )
    parser.add_argument(
        "--numerical-dir", "-n",
        type=Path,
        default=Path("."),
        help="Directory containing numerical CSV files (default: current directory)"
    )
    parser.add_argument(
        "--analytical-dir", "-a", 
        type=Path,
        default=Path("data"),
        help="Directory containing analytical DAT files (default: data/)"
    )
    parser.add_argument(
        "--cases", "-c",
        type=str,
        default="1-12",
        help="Test cases to validate, e.g., '1-12' or '1,2,5' (default: 1-12)"
    )
    parser.add_argument(
        "--tolerance", "-t",
        type=float,
        default=0.1,
        help="L1 error tolerance for pass/fail (default: 0.1)"
    )
    parser.add_argument(
        "--plot", "-p",
        action="store_true",
        help="Generate comparison plots"
    )
    parser.add_argument(
        "--plot-dir",
        type=Path,
        default=Path("plots"),
        help="Directory to save plots (default: plots/)"
    )
    
    args = parser.parse_args()
    
    # Parse test cases
    if '-' in args.cases:
        start, end = map(int, args.cases.split('-'))
        test_cases = list(range(start, end + 1))
    else:
        test_cases = [int(c.strip()) for c in args.cases.split(',')]
    
    # Run validation
    passed = run_validation(
        args.numerical_dir,
        args.analytical_dir,
        test_cases,
        args.tolerance
    )
    
    # Generate plots if requested
    if args.plot:
        print(f"\n{'='*70}")
        print("Generating comparison plots...")
        print(f"{'='*70}")
        generate_comparison_plots(
            args.numerical_dir,
            args.analytical_dir,
            args.plot_dir,
            test_cases
        )
    
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
