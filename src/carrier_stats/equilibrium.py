"""Charge-neutral equilibrium with selectable statistics and dopant ionization."""

from dataclasses import dataclass
from math import exp, expm1, hypot, isfinite, log
from typing import Literal
from warnings import warn

import numpy as np
from scipy.optimize import brentq

from .material import BOLTZMANN_EV_K, SILICON, SiliconParameters, validate_doping
from .statistics import log_fermi_half

Statistics = Literal["fermi-dirac", "boltzmann"]
Ionization = Literal["incomplete", "complete"]
NEGATIVE_INFINITY = float("-inf")


class ModelScopeWarning(UserWarning):
    """The numerical result extends beyond the conservative demonstration scope."""


@dataclass(frozen=True)
class EquilibriumResult:
    """Uniform bulk equilibrium; energies are referenced to E_v=0 eV.

    neutrality_relative_residual uses the reduced charge balance documented
    in docs/model.md. It remains sensitive under exact compensation.
    """

    temperature_K: float
    donor_cm3: float
    acceptor_cm3: float
    statistics: Statistics
    ionization: Ionization
    bandgap_eV: float
    fermi_level_eV: float
    electron_cm3: float
    hole_cm3: float
    ionized_donor_cm3: float
    ionized_acceptor_cm3: float
    neutral_donor_cm3: float
    neutral_acceptor_cm3: float
    neutrality_relative_residual: float
    iterations: int
    scope_notes: tuple[str, ...]


def _scaled_difference(log_left: float, log_right: float, log_scale: float) -> float:
    """Compute (left-right)/scale without subtracting almost equal exponentials."""
    if log_left == log_right:
        return 0.0
    if log_left > log_right:
        return exp(log_left - log_scale) * -expm1(log_right - log_left)
    return -exp(log_right - log_scale) * -expm1(log_left - log_right)


def solve_equilibrium(
    temperature_K: float,
    donor_cm3: float = 0.0,
    acceptor_cm3: float = 0.0,
    *,
    statistics: Statistics = "fermi-dirac",
    ionization: Ionization = "incomplete",
    material: SiliconParameters = SILICON,
    energy_tolerance_eV: float = 1e-12,
    quadrature_relative_tolerance: float = 1e-10,
) -> EquilibriumResult:
    """Solve n+N_A^- = p+N_D^+ for one Fermi level, using Brent's method.

    Temperatures: 40..1000 K. Active donor/acceptor densities: 0..1e20 cm^-3.
    Results above total doping 1e17 cm^-3 are marked as formal ideal-model
    studies, not calibrated high-doping silicon predictions. This deliberately
    conservative threshold is a project policy, not a universal phase boundary.
    No position, contacts, voltage, illumination, or transport is represented.
    """
    validate_doping(donor_cm3, acceptor_cm3)
    if statistics not in ("fermi-dirac", "boltzmann"):
        raise ValueError("statistics must be 'fermi-dirac' or 'boltzmann'")
    if ionization not in ("incomplete", "complete"):
        raise ValueError("ionization must be 'incomplete' or 'complete'")
    if not isfinite(energy_tolerance_eV) or not 1e-14 <= energy_tolerance_eV <= 1e-6:
        raise ValueError("energy_tolerance_eV must be between 1e-14 and 1e-6")
    if (
        not isfinite(quadrature_relative_tolerance)
        or not 1e-12 <= quadrature_relative_tolerance <= 1e-4
    ):
        raise ValueError("quadrature_relative_tolerance must be between 1e-12 and 1e-4")

    gap_eV = material.bandgap_eV(temperature_K)
    thermal_eV = BOLTZMANN_EV_K * temperature_K
    conduction_dos, valence_dos = material.density_of_states_cm3(temperature_K)
    log_conduction_dos, log_valence_dos = log(conduction_dos), log(valence_dos)
    log_donors = log(donor_cm3) if donor_cm3 else NEGATIVE_INFINITY
    log_acceptors = log(acceptor_cm3) if acceptor_cm3 else NEGATIVE_INFINITY
    net_doping = donor_cm3 - acceptor_cm3
    log_net_doping = log(abs(net_doping)) if net_doping else NEGATIVE_INFINITY

    def state_logs(fermi_eV: float) -> tuple[float, float, float, float]:
        electron_eta = (fermi_eV - gap_eV) / thermal_eV
        hole_eta = -fermi_eV / thermal_eV
        if statistics == "fermi-dirac":
            electron_log = log_conduction_dos + log_fermi_half(
                electron_eta, relative_tolerance=quadrature_relative_tolerance
            )
            hole_log = log_valence_dos + log_fermi_half(
                hole_eta, relative_tolerance=quadrature_relative_tolerance
            )
        else:
            electron_log = log_conduction_dos + electron_eta
            hole_log = log_valence_dos + hole_eta
        if ionization == "complete":
            return electron_log, hole_log, NEGATIVE_INFINITY, NEGATIVE_INFINITY
        donor_argument = (
            log(material.donor_degeneracy)
            + (fermi_eV - gap_eV + material.donor_binding_eV) / thermal_eV
        )
        acceptor_argument = (
            log(material.acceptor_degeneracy)
            + (material.acceptor_binding_eV - fermi_eV) / thermal_eV
        )
        # log(sigmoid(argument)) retains tiny neutral fractions under compensation.
        neutral_donor_log = log_donors - float(np.logaddexp(0, -donor_argument))
        neutral_acceptor_log = log_acceptors - float(np.logaddexp(0, -acceptor_argument))
        return electron_log, hole_log, neutral_donor_log, neutral_acceptor_log

    def residual(fermi_eV: float) -> float:
        electron_log, hole_log, donor_neutral_log, acceptor_neutral_log = state_logs(fermi_eV)
        log_scale = max(
            electron_log, hole_log, donor_neutral_log, acceptor_neutral_log, log_net_doping
        )
        net_scaled = np.sign(net_doping) * exp(log_net_doping - log_scale) if net_doping else 0.0
        return float(
            _scaled_difference(electron_log, hole_log, log_scale)
            + _scaled_difference(donor_neutral_log, acceptor_neutral_log, log_scale)
            - net_scaled
        )

    # Bracketing limits are numerical bounds, not spatial boundary conditions.
    fermi_eV, solver = brentq(
        residual, -1.0, gap_eV + 1.0, xtol=energy_tolerance_eV, rtol=1e-14, full_output=True
    )
    electron_log, hole_log, donor_neutral_log, acceptor_neutral_log = state_logs(fermi_eV)
    electron_cm3, hole_cm3 = exp(electron_log), exp(hole_log)
    donor_neutral_cm3, acceptor_neutral_cm3 = exp(donor_neutral_log), exp(acceptor_neutral_log)
    if ionization == "complete":
        donor_ion_cm3, acceptor_ion_cm3 = donor_cm3, acceptor_cm3
    else:
        donor_argument = (
            log(material.donor_degeneracy)
            + (fermi_eV - gap_eV + material.donor_binding_eV) / thermal_eV
        )
        acceptor_argument = (
            log(material.acceptor_degeneracy)
            + (material.acceptor_binding_eV - fermi_eV) / thermal_eV
        )
        donor_ion_cm3 = exp(log_donors - float(np.logaddexp(0, donor_argument)))
        acceptor_ion_cm3 = exp(log_acceptors - float(np.logaddexp(0, acceptor_argument)))
    notes = []
    if donor_cm3 + acceptor_cm3 > 1e17:
        notes.append(
            "Total doping exceeds 1e17 cm^-3: formal ideal-model study only; "
            "dopant interactions and band-gap narrowing are omitted."
        )
    if temperature_K > 500:
        notes.append(
            "Above 500 K, fixed effective masses limit quantitative accuracy; "
            "use this result for qualitative temperature trends."
        )
    if notes:
        warn(" ".join(notes), ModelScopeWarning, stacklevel=2)
    return EquilibriumResult(
        temperature_K,
        donor_cm3,
        acceptor_cm3,
        statistics,
        ionization,
        gap_eV,
        fermi_eV,
        electron_cm3,
        hole_cm3,
        donor_ion_cm3,
        acceptor_ion_cm3,
        donor_neutral_cm3,
        acceptor_neutral_cm3,
        residual(fermi_eV),
        solver.iterations,
        tuple(notes),
    )


def boltzmann_complete_solution(
    temperature_K: float,
    donor_cm3: float = 0.0,
    acceptor_cm3: float = 0.0,
    *,
    material: SiliconParameters = SILICON,
) -> tuple[float, float]:
    """Analytical (electron, hole) densities in cm^-3 with complete ionization.

    Solve n-p=N_D-N_A and np=n_i^2. Calculate the majority density first
    to avoid cancellation in the minority-carrier quadratic root.
    """
    validate_doping(donor_cm3, acceptor_cm3)
    intrinsic_cm3 = material.intrinsic_boltzmann_cm3(temperature_K)
    net_doping = donor_cm3 - acceptor_cm3
    majority_cm3 = 0.5 * (abs(net_doping) + hypot(net_doping, 2 * intrinsic_cm3))
    minority_cm3 = (intrinsic_cm3 / majority_cm3) * intrinsic_cm3
    if net_doping >= 0:
        return majority_cm3, minority_cm3
    return minority_cm3, majority_cm3
