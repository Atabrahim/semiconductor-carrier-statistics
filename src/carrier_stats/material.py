"""Silicon parameters, with provenance and model scope in docs/model.md."""

from dataclasses import dataclass
from math import exp, isfinite, pi, sqrt

from scipy.constants import Boltzmann, electron_mass, elementary_charge, h

BOLTZMANN_EV_K = Boltzmann / elementary_charge
MIN_TEMPERATURE_K = 40.0
MAX_TEMPERATURE_K = 1000.0


def validate_temperature(temperature_K: float) -> None:
    """Restrict calculations to the documented demonstration range, in kelvin."""
    if not isfinite(temperature_K) or not MIN_TEMPERATURE_K <= temperature_K <= MAX_TEMPERATURE_K:
        raise ValueError("temperature_K must be finite and between 40 and 1000 K")


def validate_doping(donor_cm3: float, acceptor_cm3: float) -> None:
    """Validate active substitutional dopant densities, in cm^-3."""
    for name, density in (("donor_cm3", donor_cm3), ("acceptor_cm3", acceptor_cm3)):
        if not isfinite(density) or not 0 <= density <= 1e20:
            raise ValueError(f"{name} must be finite and between 0 and 1e20 cm^-3")


@dataclass(frozen=True)
class SiliconParameters:
    """Fixed-mass bulk Si model; energy in eV, relative masses in units of m_e.

    Gap and masses: Palankovski (2000), Tables 3.13 and 3.16.
    Donor: isolated phosphorus, Smith et al. (2016), 45.59 meV.
    Acceptor: isolated boron, Meng et al. (2023), 45 meV.
    Overrides support controlled sensitivity studies, not arbitrary materials.
    """

    bandgap_0_eV: float = 1.1695
    varshni_alpha_eV_K: float = 4.73e-4
    varshni_beta_K: float = 636.0
    electron_transverse_mass_ratio: float = 0.19
    electron_longitudinal_mass_ratio: float = 0.98
    heavy_hole_mass_ratio: float = 0.49
    light_hole_mass_ratio: float = 0.16
    conduction_valleys: int = 6
    donor_binding_eV: float = 0.04559
    acceptor_binding_eV: float = 0.045
    donor_degeneracy: float = 2.0
    acceptor_degeneracy: float = 4.0

    def __post_init__(self) -> None:
        for name, value in vars(self).items():
            if not isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        if self.conduction_valleys != 6:
            raise ValueError("This silicon model requires six conduction valleys")
        smallest_gap = self.bandgap_eV(MAX_TEMPERATURE_K)
        if max(self.donor_binding_eV, self.acceptor_binding_eV) >= smallest_gap:
            raise ValueError("Dopant binding energies must lie within the band gap")

    def bandgap_eV(self, temperature_K: float) -> float:
        """Return the Varshni band gap in eV; this is a semi-empirical model."""
        validate_temperature(temperature_K)
        return self.bandgap_0_eV - (
            self.varshni_alpha_eV_K * temperature_K**2 / (temperature_K + self.varshni_beta_K)
        )

    def density_of_states_cm3(self, temperature_K: float) -> tuple[float, float]:
        """Return (N_c, N_v) in cm^-3, including spin and valley multiplicity.

        SI constants yield m^-3 before the explicit factor of 1e-6.
        Effective masses are held constant: temperature-dependent DOS masses,
        valence-band warping, and the split-off band are not modeled.
        """
        validate_temperature(temperature_K)
        electron_dos_mass = (
            self.electron_transverse_mass_ratio**2 * self.electron_longitudinal_mass_ratio
        ) ** (1 / 3)
        hole_dos_mass = (self.heavy_hole_mass_ratio**1.5 + self.light_hole_mass_ratio**1.5) ** (
            2 / 3
        )
        prefactor_cm3 = (
            2 * (2 * pi * electron_mass * Boltzmann * temperature_K / h**2) ** 1.5 / 1e6
        )
        return (
            prefactor_cm3 * self.conduction_valleys * electron_dos_mass**1.5,
            prefactor_cm3 * hole_dos_mass**1.5,
        )

    def intrinsic_boltzmann_cm3(self, temperature_K: float) -> float:
        """Return analytical Boltzmann n_i in cm^-3 for this parameter set."""
        conduction_dos_cm3, valence_dos_cm3 = self.density_of_states_cm3(temperature_K)
        return sqrt(conduction_dos_cm3 * valence_dos_cm3) * exp(
            -self.bandgap_eV(temperature_K) / (2 * BOLTZMANN_EV_K * temperature_K)
        )


SILICON = SiliconParameters()
