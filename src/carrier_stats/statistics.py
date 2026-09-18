"""Normalized Fermi-Dirac integral of order 1/2 (NIST DLMF 25.12.14)."""

from math import exp, isfinite, log, pi, sqrt

from scipy.integrate import quad
from scipy.special import expit


def log_fermi_half(reduced_fermi_energy: float, *, relative_tolerance: float = 1e-10) -> float:
    """Return ln(F_1/2(eta)) for dimensionless eta, using adaptive quadrature.

    F_1/2(eta) = 2/sqrt(pi) integral_0^inf sqrt(t)/(1+exp(t-eta)) dt.
    Substitution t=u^2 smooths the endpoint. For eta<=0, exp(eta) is
    factored out before integration to preserve relative accuracy in the tail.
    For eta>0, the upper bound t=eta+50 leaves a negligible exponential tail.
    The supported eta interval is [-10000, 10000]. No interpolation table or
    fitted approximation is used. Tolerance is numerical, not physical uncertainty.
    """
    eta = reduced_fermi_energy
    if not isfinite(eta) or abs(eta) > 10000:
        raise ValueError("reduced_fermi_energy must be finite and within [-10000, 10000]")
    if not isfinite(relative_tolerance) or not 1e-12 <= relative_tolerance <= 1e-4:
        raise ValueError("relative_tolerance must be between 1e-12 and 1e-4")
    options = {"epsabs": 1e-13, "epsrel": relative_tolerance, "limit": 150}
    normalization = 4 / sqrt(pi)
    if eta <= 0:
        integral, _ = quad(
            lambda u: u * u * exp(-u * u) / (1 + exp(eta - u * u)),
            0,
            float("inf"),
            **options,
        )
        return eta + log(normalization * integral)
    integrand = lambda u: u * u * expit(eta - u * u)  # noqa: E731
    edge = sqrt(eta)
    # Resolve both sides of the narrow Fermi edge even when eta is very large.
    lower_edge = sqrt(max(eta - 50, 0))
    bulk, _ = quad(integrand, 0, lower_edge, **options)
    lower, _ = quad(integrand, lower_edge, edge, **options)
    upper, _ = quad(integrand, edge, sqrt(eta + 50), **options)
    return log(normalization * (bulk + lower + upper))


def fermi_half(reduced_fermi_energy: float, *, relative_tolerance: float = 1e-10) -> float:
    """Return normalized F_1/2(eta); underflow may return zero for very negative eta.

    Use log_fermi_half when representing extremely dilute carrier populations.
    """
    return exp(log_fermi_half(reduced_fermi_energy, relative_tolerance=relative_tolerance))
