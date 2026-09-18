"""Independent representations and asymptotic limits of the Fermi integral."""

from math import exp, pi, sqrt

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import expit, zeta

from carrier_stats import fermi_half, log_fermi_half


def test_exact_zero_argument_identity():
    expected = (1 - 2**-0.5) * zeta(1.5, 1)
    assert fermi_half(0) == pytest.approx(expected, rel=2e-12)


@pytest.mark.parametrize("eta", [-12, -2, 0, 2, 10, 100])
def test_independent_energy_variable_quadrature(eta):
    reference, _ = quad(
        lambda energy: sqrt(energy) * expit(eta - energy),
        0,
        max(100, eta + 100),
        epsabs=1e-13,
        epsrel=1e-11,
        limit=300,
    )
    assert fermi_half(eta) == pytest.approx(2 / sqrt(pi) * reference, rel=2e-10)


@pytest.mark.parametrize("eta", [-20, -100, -1000])
def test_boltzmann_tail_without_underflow(eta):
    assert log_fermi_half(eta) == pytest.approx(eta, abs=1e-9)


@pytest.mark.parametrize("eta", [30, 100, 10000])
def test_sommerfeld_limit(eta):
    expected = 4 / (3 * sqrt(pi)) * eta**1.5 * (1 + pi**2 / (8 * eta**2))
    assert fermi_half(eta) == pytest.approx(expected, rel=3e-6)


def test_fermi_integral_is_positive_increasing_and_below_boltzmann():
    grid = np.linspace(-15, 10, 45)
    values = np.array([fermi_half(eta) for eta in grid])
    assert np.all(values > 0)
    assert np.all(np.diff(values) > 0)
    assert np.all(values < np.exp(grid))


def test_quadrature_convergence():
    for eta in [-30, -0.01, 0.01, 5, 50]:
        loose = fermi_half(eta, relative_tolerance=1e-7)
        tight = fermi_half(eta, relative_tolerance=1e-12)
        assert loose == pytest.approx(tight, rel=1e-8)


def test_negative_tail_underflow_is_explicit():
    assert fermi_half(-1000) == 0
    assert log_fermi_half(-1000) == pytest.approx(-1000)
    assert fermi_half(-10) == pytest.approx(exp(-10), rel=2e-5)


@pytest.mark.parametrize("eta", [float("nan"), float("inf"), -10001, 10001])
def test_reject_invalid_reduced_energy(eta):
    with pytest.raises(ValueError, match="reduced_fermi_energy"):
        fermi_half(eta)


@pytest.mark.parametrize("tolerance", [0, 1e-16, 0.01, float("nan")])
def test_reject_invalid_quadrature_tolerance(tolerance):
    with pytest.raises(ValueError, match="relative_tolerance"):
        fermi_half(0, relative_tolerance=tolerance)
