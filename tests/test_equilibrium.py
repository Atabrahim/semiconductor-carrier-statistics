"""Test conservation, analytical cases, compensation, and model-scope handling."""

import json
import subprocess
import sys
import warnings
from dataclasses import replace
from math import exp, sqrt

import pytest

from carrier_stats import (
    BOLTZMANN_EV_K,
    SILICON,
    ModelScopeWarning,
    boltzmann_complete_solution,
    solve_equilibrium,
)


@pytest.mark.parametrize("temperature_K", [40, 100, 300, 500])
@pytest.mark.parametrize("donors,acceptors", [(0, 0), (1e15, 0), (0, 1e16), (1e17, 1e17)])
def test_numerical_boltzmann_matches_stable_analytical_solution(temperature_K, donors, acceptors):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelScopeWarning)
        result = solve_equilibrium(
            temperature_K, donors, acceptors, statistics="boltzmann", ionization="complete"
        )
    expected_n, expected_p = boltzmann_complete_solution(temperature_K, donors, acceptors)
    assert result.electron_cm3 == pytest.approx(expected_n, rel=2e-9)
    assert result.hole_cm3 == pytest.approx(expected_p, rel=2e-9)
    assert abs(result.neutrality_relative_residual) < 2e-9


@pytest.mark.parametrize("statistics", ["fermi-dirac", "boltzmann"])
@pytest.mark.parametrize("ionization", ["incomplete", "complete"])
@pytest.mark.parametrize("temperature_K", [40, 100, 300, 500])
@pytest.mark.parametrize("donors,acceptors", [(0, 0), (1e15, 0), (0, 1e15), (1e16, 9e15)])
def test_charge_and_dopant_conservation(statistics, ionization, temperature_K, donors, acceptors):
    result = solve_equilibrium(
        temperature_K, donors, acceptors, statistics=statistics, ionization=ionization
    )
    positive = result.hole_cm3 + result.ionized_donor_cm3
    negative = result.electron_cm3 + result.ionized_acceptor_cm3
    assert positive == pytest.approx(negative, rel=2e-9, abs=0)
    assert abs(result.neutrality_relative_residual) < 2e-9
    assert result.ionized_donor_cm3 + result.neutral_donor_cm3 == pytest.approx(donors, rel=1e-12)
    assert result.ionized_acceptor_cm3 + result.neutral_acceptor_cm3 == pytest.approx(
        acceptors, rel=1e-12
    )
    assert result.electron_cm3 > 0 and result.hole_cm3 > 0


@pytest.mark.parametrize("donors,acceptors", [(1e15, 0), (0, 1e15), (1e16, 8e15)])
def test_mass_action_in_boltzmann_statistics_with_incomplete_ionization(donors, acceptors):
    result = solve_equilibrium(300, donors, acceptors, statistics="boltzmann")
    assert result.electron_cm3 * result.hole_cm3 == pytest.approx(
        SILICON.intrinsic_boltzmann_cm3(300) ** 2, rel=2e-12
    )


def test_donor_freeze_out_against_analytical_quadratic():
    temperature_K, donor_cm3 = 50, 1e15
    conduction_dos, _ = SILICON.density_of_states_cm3(temperature_K)
    coefficient_cm3 = (
        SILICON.donor_degeneracy
        / conduction_dos
        * exp(SILICON.donor_binding_eV / (BOLTZMANN_EV_K * temperature_K))
    )
    # Neglect holes: a*n^2+n-N_D=0; stable positive root, not the solver's residual.
    expected = 2 * donor_cm3 / (1 + sqrt(1 + 4 * coefficient_cm3 * donor_cm3))
    result = solve_equilibrium(temperature_K, donor_cm3, statistics="boltzmann")
    assert result.electron_cm3 == pytest.approx(expected, rel=2e-10)
    assert result.electron_cm3 < 0.2 * donor_cm3


def test_acceptor_freeze_out_against_analytical_quadratic():
    temperature_K, acceptor_cm3 = 50, 1e15
    _, valence_dos = SILICON.density_of_states_cm3(temperature_K)
    coefficient_cm3 = (
        SILICON.acceptor_degeneracy
        / valence_dos
        * exp(SILICON.acceptor_binding_eV / (BOLTZMANN_EV_K * temperature_K))
    )
    expected = 2 * acceptor_cm3 / (1 + sqrt(1 + 4 * coefficient_cm3 * acceptor_cm3))
    result = solve_equilibrium(temperature_K, acceptor_cm3=acceptor_cm3, statistics="boltzmann")
    assert result.hole_cm3 == pytest.approx(expected, rel=2e-10)


def test_freeze_out_extrinsic_and_intrinsic_limits():
    cold = solve_equilibrium(40, 1e15)
    room = solve_equilibrium(300, 1e15)
    with pytest.warns(ModelScopeWarning, match="fixed effective masses"):
        hot = solve_equilibrium(1000, 1e15)
    assert cold.electron_cm3 < 0.05 * 1e15
    assert room.electron_cm3 == pytest.approx(1e15, rel=0.001)
    assert hot.electron_cm3 / hot.hole_cm3 == pytest.approx(1, rel=0.05)


def test_exact_compensation_with_incomplete_ionization_remains_resolved():
    result = solve_equilibrium(40, 1e15, 1e15)
    assert abs(result.neutrality_relative_residual) < 1e-8
    assert 0.45 < result.fermi_level_eV < 0.7
    assert result.electron_cm3 > 0 and result.hole_cm3 > 0


def test_statistics_difference_at_fixed_doping_is_mostly_fermi_level_shift():
    with pytest.warns(ModelScopeWarning):
        fd = solve_equilibrium(300, 1e20, ionization="complete")
    with pytest.warns(ModelScopeWarning):
        mb = solve_equilibrium(300, 1e20, statistics="boltzmann", ionization="complete")
    assert fd.electron_cm3 == pytest.approx(mb.electron_cm3, rel=1e-9)
    assert fd.fermi_level_eV > mb.fermi_level_eV + 0.01
    assert fd.hole_cm3 < mb.hole_cm3
    assert fd.scope_notes


@pytest.mark.parametrize("temperature_K", [39, 1001, 0, -1, float("nan"), float("inf")])
def test_invalid_temperature(temperature_K):
    with pytest.raises(ValueError, match="temperature_K"):
        solve_equilibrium(temperature_K)


@pytest.mark.parametrize(
    "key,value",
    [
        ("donor_cm3", -1),
        ("acceptor_cm3", -1),
        ("donor_cm3", 1e21),
        ("donor_cm3", float("nan")),
        ("acceptor_cm3", float("inf")),
        ("statistics", "wrong"),
        ("ionization", "wrong"),
        ("energy_tolerance_eV", 0),
        ("quadrature_relative_tolerance", 0),
    ],
)
def test_invalid_solver_arguments(key, value):
    with pytest.raises(ValueError):
        solve_equilibrium(300, **{key: value})


def test_gap_units_and_density_of_states_scaling():
    assert SILICON.bandgap_eV(300) == pytest.approx(1.1240192307692307, rel=1e-12)
    conduction_300, valence_300 = SILICON.density_of_states_cm3(300)
    conduction_150, valence_150 = SILICON.density_of_states_cm3(150)
    assert 2.7e19 < conduction_300 < 2.9e19
    assert 1.0e19 < valence_300 < 1.1e19
    assert conduction_300 / conduction_150 == pytest.approx(2**1.5)
    assert valence_300 / valence_150 == pytest.approx(2**1.5)
    assert SILICON.bandgap_eV(500) < SILICON.bandgap_eV(100)


def test_parameter_validation_and_binding_energy_sensitivity():
    with pytest.raises(ValueError):
        replace(SILICON, donor_binding_eV=-0.1)
    with pytest.raises(ValueError):
        replace(SILICON, donor_binding_eV=2)
    nominal = solve_equilibrium(50, 1e15)
    stronger_binding = solve_equilibrium(
        50, 1e15, material=replace(SILICON, donor_binding_eV=0.046)
    )
    assert stronger_binding.electron_cm3 < nominal.electron_cm3


def test_cli_valid_json_and_invalid_input():
    good = subprocess.run(
        [sys.executable, "-m", "carrier_stats", "--temperature", "300", "--donors", "1e15"],
        capture_output=True,
        text=True,
        check=True,
    )
    payload = json.loads(good.stdout)
    assert payload["result_type"] == "NUMERICAL RESULT"
    assert 9.9e14 < payload["electron_cm3"] < 1.01e15
    bad = subprocess.run(
        [sys.executable, "-m", "carrier_stats", "--temperature", "-1"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert bad.returncode == 2
    assert "temperature_K" in bad.stderr and "Traceback" not in bad.stderr
