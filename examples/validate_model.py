"""Generate a machine-readable validation report; exit nonzero if any check fails."""

import argparse
import json
import platform
import warnings
from math import exp, sqrt
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import brentq
from scipy.special import zeta

from carrier_stats import (
    SILICON,
    ModelScopeWarning,
    boltzmann_complete_solution,
    fermi_half,
    solve_equilibrium,
)


def validation_report() -> dict:
    """Exercise a declared grid and compare independent analytical expressions."""
    checks = []

    def add_check(name: str, measured: float, maximum: float) -> None:
        checks.append(
            {
                "name": name,
                "measured_error": float(measured),
                "maximum_error": maximum,
                "passed": bool(measured <= maximum),
            }
        )

    fd_exact = (1 - 1 / sqrt(2)) * zeta(1.5, 1)
    add_check(
        "F_half(0): relative error against zeta identity", abs(fermi_half(0) / fd_exact - 1), 2e-12
    )
    convergence = max(
        abs(
            fermi_half(eta, relative_tolerance=1e-7) / fermi_half(eta, relative_tolerance=1e-12)
            - 1
        )
        for eta in [-30, -1, 0, 1, 10, 100, 10000]
    )
    add_check("F_half: relative change when quadrature tolerance tightened", convergence, 1e-8)
    neutrality_errors, analytical_errors = [], []
    temperatures = [40, 50, 100, 300, 500, 1000]
    doping_pairs = [
        (0, 0),
        (1e12, 0),
        (1e15, 0),
        (0, 1e16),
        (1e16, 9e15),
        (1e15, 1e15),
        (1e17, 1e17),
        (1e20, 0),
    ]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelScopeWarning)
        for temperature_K in temperatures:
            for donors, acceptors in doping_pairs:
                for statistics in ["boltzmann", "fermi-dirac"]:
                    for ionization in ["complete", "incomplete"]:
                        result = solve_equilibrium(
                            temperature_K,
                            donors,
                            acceptors,
                            statistics=statistics,
                            ionization=ionization,
                        )
                        neutrality_errors.append(abs(result.neutrality_relative_residual))
                        if statistics == "boltzmann" and ionization == "complete":
                            expected_n, expected_p = boltzmann_complete_solution(
                                temperature_K, donors, acceptors
                            )
                            analytical_errors.extend(
                                [
                                    abs(result.electron_cm3 / expected_n - 1),
                                    abs(result.hole_cm3 / expected_p - 1),
                                ]
                            )
        room = solve_equilibrium(300, 1e15)
        cold = solve_equilibrium(40, 1e15)
    add_check("Equilibrium: maximum reduced charge-balance residual", max(neutrality_errors), 2e-9)
    add_check(
        "Boltzmann/complete: maximum relative error against quadratic",
        max(analytical_errors),
        2e-9,
    )
    return {
        "report_type": "NUMERICAL VERIFICATION; NOT EXPERIMENTAL VALIDATION",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "equilibrium_cases": len(neutrality_errors),
        "temperature_grid_K": temperatures,
        "doping_pairs_cm3": doping_pairs,
        "scope": "Includes formal stress tests outside the physical model's intended scope.",
        "checks": checks,
        "all_passed": all(check["passed"] for check in checks),
        "reference_model_values": {
            "bandgap_300K_eV": SILICON.bandgap_eV(300),
            "effective_dos_300K_cm3": SILICON.density_of_states_cm3(300),
            "intrinsic_boltzmann_300K_cm3": SILICON.intrinsic_boltzmann_cm3(300),
            "donor_1e15_300K_electron_cm3": room.electron_cm3,
            "donor_1e15_300K_fermi_level_eV": room.fermi_level_eV,
            "donor_1e15_40K_ionized_fraction": cold.ionized_donor_cm3 / 1e15,
            "fixed_fermi_energy_1percent_boltzmann_error_eta": brentq(
                lambda eta: exp(eta) / fermi_half(eta) - 1.01, -10, 0
            ),
            "fixed_fermi_energy_5percent_boltzmann_error_eta": brentq(
                lambda eta: exp(eta) / fermi_half(eta) - 1.05, -10, 0
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("outputs/validation.json"))
    args = parser.parse_args()
    report = validation_report()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    if not report["all_passed"]:
        raise SystemExit("Scientific verification failed; inspect the report.")


if __name__ == "__main__":
    main()
