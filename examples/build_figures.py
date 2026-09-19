"""Reproduce the scientific figures and their numerical CSV data.

Run from the repository root: python examples/build_figures.py --output .
The default output directory is outputs/. No experimental data are used.
"""

import argparse
import csv
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from carrier_stats import SILICON, ModelScopeWarning, fermi_half, solve_equilibrium

COLORS = {"electron": "#126E82", "hole": "#C05526", "reference": "#606B77", "accent": "#7347A0"}


def write_csv(path: Path, rows: list[dict]) -> None:
    """Export explicitly labeled model results; field names carry physical units."""
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save_figure(figure, directory: Path, stem: str) -> None:
    """Save a README preview and an editable vector figure."""
    figure.savefig(directory / f"{stem}.png", dpi=180, bbox_inches="tight")
    figure.savefig(directory / f"{stem}.svg", bbox_inches="tight", metadata={"Date": None})
    plt.close(figure)


def carrier_regimes(figures: Path, data: Path) -> None:
    """Show where donor ionization and intrinsic carriers determine equilibrium."""
    temperature_grid = np.geomspace(40, 1000, 180)
    incomplete = [solve_equilibrium(t, 1e15) for t in temperature_grid]
    complete = [solve_equilibrium(t, 1e15, ionization="complete") for t in temperature_grid]
    intrinsic = [SILICON.intrinsic_boltzmann_cm3(t) for t in temperature_grid]
    figure, axes = plt.subplots(1, 3, figsize=(14, 4.6), layout="constrained")
    figure.suptitle(r"When is one free electron per donor a good approximation?", fontsize=16)
    axes[0].loglog(
        temperature_grid,
        [r.electron_cm3 for r in incomplete],
        color=COLORS["electron"],
        label="Electrons: incomplete ionization",
    )
    axes[0].loglog(
        temperature_grid,
        [r.electron_cm3 for r in complete],
        "--",
        color=COLORS["hole"],
        label="Electrons: complete ionization",
    )
    axes[0].loglog(
        temperature_grid,
        intrinsic,
        ":",
        color=COLORS["reference"],
        label=r"Analytical $n_i$: Boltzmann",
    )
    axes[0].set(
        ylim=(1e12, 2e18),
        ylabel=r"Carrier concentration (cm$^{-3}$)",
        title="A  |  Freeze-out to intrinsic conduction",
    )
    axes[0].legend(loc="upper left", fontsize=8)
    axes[1].semilogx(
        temperature_grid,
        [r.ionized_donor_cm3 / 1e15 for r in incomplete],
        color=COLORS["electron"],
    )
    axes[1].axhline(0.99, color=COLORS["reference"], ls=":", label="99% ionized")
    axes[1].set(
        ylim=(0, 1.05),
        ylabel=r"Ionized fraction $N_D^+/N_D$",
        title="B  |  Thermal occupation of donors",
    )
    axes[1].legend(fontsize=8)
    axes[2].semilogx(
        temperature_grid,
        [r.bandgap_eV for r in incomplete],
        color=COLORS["reference"],
        label=r"Conduction edge $E_C$",
    )
    axes[2].semilogx(
        temperature_grid,
        [r.bandgap_eV - SILICON.donor_binding_eV for r in incomplete],
        color=COLORS["reference"],
        ls=":",
        label=r"Donor level $E_D$",
    )
    axes[2].semilogx(
        temperature_grid,
        [r.fermi_level_eV for r in incomplete],
        color=COLORS["accent"],
        label=r"Fermi level $E_F$",
    )
    axes[2].set(
        ylabel=r"Energy relative to $E_V$ (eV)",
        ylim=(0, 1.25),
        title="C  |  Equilibrium chemical potential",
    )
    axes[2].legend(fontsize=8, loc="lower left")
    for axis in axes:
        axis.set(xlabel="Temperature (K)", xlim=(40, 1000))
        axis.set_xticks([40, 100, 300, 1000], labels=["40", "100", "300", "1000"])
        axis.axvspan(500, 1000, color="#EBCB8B", alpha=0.24, zorder=0)
    figure.supxlabel(
        "NUMERICAL MODEL · Si:P, active donors = 10¹⁵ cm⁻³, no acceptors · "
        "Shading: fixed-mass extrapolation above 500 K",
        fontsize=9,
    )
    save_figure(figure, figures, "carrier_regimes")
    write_csv(
        data / "carrier_regimes.csv",
        [
            {
                "result_type": "NUMERICAL RESULT",
                "temperature_K": r.temperature_K,
                "electron_incomplete_cm3": r.electron_cm3,
                "hole_incomplete_cm3": r.hole_cm3,
                "electron_complete_cm3": full.electron_cm3,
                "intrinsic_boltzmann_analytical_cm3": ni,
                "ionized_donor_fraction": r.ionized_donor_cm3 / 1e15,
                "fermi_level_eV": r.fermi_level_eV,
                "bandgap_eV": r.bandgap_eV,
                "scope_notes": " | ".join(r.scope_notes),
            }
            for r, full, ni in zip(incomplete, complete, intrinsic, strict=True)
        ],
    )


def approximation_limits(figures: Path, data: Path) -> None:
    """Separate ionization error at fixed doping from statistics error at fixed eta."""
    temperatures = np.linspace(40, 400, 36)
    donors = np.geomspace(1e12, 1e17, 32)
    errors = np.empty((len(donors), len(temperatures)))
    rows = []
    for row, donor_cm3 in enumerate(donors):
        for column, temperature_K in enumerate(temperatures):
            partial = solve_equilibrium(temperature_K, donor_cm3)
            full = solve_equilibrium(temperature_K, donor_cm3, ionization="complete")
            error_percent = 100 * (full.electron_cm3 / partial.electron_cm3 - 1)
            errors[row, column] = error_percent
            rows.append(
                {
                    "result_type": "NUMERICAL RESULT",
                    "temperature_K": temperature_K,
                    "donor_cm3": donor_cm3,
                    "complete_ionization_error_percent": error_percent,
                }
            )
    figure, axes = plt.subplots(1, 2, figsize=(12, 4.8), layout="constrained")
    figure.suptitle("Two approximations, two different error questions", fontsize=16)
    mesh = axes[0].pcolormesh(
        temperatures,
        donors,
        errors,
        shading="auto",
        norm=LogNorm(vmin=0.001, vmax=10000),
        cmap="magma",
    )
    axes[0].set(
        yscale="log",
        xlabel="Temperature (K)",
        ylabel=r"Active donors (cm$^{-3}$)",
        title="A  |  Complete-ionization error at fixed doping",
    )
    # Set the log scale before choosing label positions in display coordinates.
    contours = axes[0].contour(
        temperatures, donors, errors, levels=[1, 5, 50], colors="white", linewidths=1
    )
    axes[0].clabel(contours, fmt=lambda level: f"{level:g}%", fontsize=9)
    figure.colorbar(mesh, ax=axes[0], label="Electron-density overestimate (%)", extend="both")
    reduced_energies = np.linspace(-10, 4, 240)
    statistical_errors = np.array(
        [100 * (np.exp(eta) / fermi_half(eta) - 1) for eta in reduced_energies]
    )
    axes[1].semilogy(reduced_energies, statistical_errors, color=COLORS["electron"])
    for level, color in [(1, COLORS["reference"]), (5, COLORS["hole"])]:
        axes[1].axhline(level, ls="--", color=color, label=f"{level}% error")
    axes[1].set(
        xlabel=r"Reduced Fermi energy $\eta=(E_F-E_C)/(k_B T)$",
        ylabel=r"$100\,[\exp(\eta)/F_{1/2}(\eta)-1]$ (%)",
        title="B  |  Boltzmann error at fixed Fermi energy",
        xlim=(-10, 4),
    )
    axes[1].legend(fontsize=9)
    figure.supxlabel(
        "MODEL COMPARISONS · A: Fermi–Dirac carriers in both models, isolated Si:P donors · "
        "B: ideal parabolic band",
        fontsize=9,
    )
    save_figure(figure, figures, "approximation_limits")
    write_csv(data / "ionization_errors.csv", rows)
    write_csv(
        data / "statistics_errors.csv",
        [
            {
                "result_type": "NUMERICAL RESULT",
                "reduced_fermi_energy": eta,
                "boltzmann_density_error_percent": error,
            }
            for eta, error in zip(reduced_energies, statistical_errors, strict=True)
        ],
    )


def degeneracy_shift(figures: Path, data: Path) -> None:
    """At fixed complete ionization, neutrality fixes n but statistics shift E_F."""
    donors = np.geomspace(1e15, 1e20, 110)
    fd = [solve_equilibrium(300, density, ionization="complete") for density in donors]
    mb = [
        solve_equilibrium(300, density, ionization="complete", statistics="boltzmann")
        for density in donors
    ]
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.7), layout="constrained")
    figure.suptitle("Degeneracy changes the Fermi level before the majority density", fontsize=15)
    axes[0].semilogx(
        donors,
        [r.fermi_level_eV - r.bandgap_eV for r in fd],
        color=COLORS["electron"],
        label="Fermi–Dirac",
    )
    axes[0].semilogx(
        donors,
        [r.fermi_level_eV - r.bandgap_eV for r in mb],
        "--",
        color=COLORS["hole"],
        label="Boltzmann",
    )
    axes[0].axhline(0, color=COLORS["reference"], lw=0.8)
    axes[0].set(
        xlabel=r"Active donors (cm$^{-3}$)",
        ylabel=r"$E_F-E_C$ (eV)",
        title="A  |  Energy required to accommodate carriers",
    )
    axes[0].legend()
    axes[1].semilogx(
        donors,
        [a.electron_cm3 / b.electron_cm3 for a, b in zip(fd, mb, strict=True)],
        color=COLORS["electron"],
        label=r"Electron ratio $n_{FD}/n_{MB}$",
    )
    axes[1].semilogx(
        donors,
        [a.hole_cm3 / b.hole_cm3 for a, b in zip(fd, mb, strict=True)],
        color=COLORS["hole"],
        label=r"Hole ratio $p_{FD}/p_{MB}$",
    )
    axes[1].set(
        xlabel=r"Active donors (cm$^{-3}$)",
        ylabel="Density ratio (dimensionless)",
        ylim=(0, 1.08),
        title="B  |  Neutrality constrains the majority carrier",
    )
    axes[1].legend(loc="lower left")
    for axis in axes:
        axis.axvspan(1e17, 1e20, color="#EBCB8B", alpha=0.24, zorder=0)
        axis.set_xlim(1e15, 1e20)
    figure.supxlabel(
        "FORMAL IDEAL-MODEL COMPARISON · 300 K, complete ionization · "
        "Shading: beyond conservative doping scope; no band-gap narrowing",
        fontsize=9,
    )
    save_figure(figure, figures, "degeneracy_shift")
    write_csv(
        data / "degeneracy_shift.csv",
        [
            {
                "result_type": "NUMERICAL RESULT",
                "donor_cm3": density,
                "fermi_dirac_fermi_level_eV": a.fermi_level_eV,
                "boltzmann_fermi_level_eV": b.fermi_level_eV,
                "electron_ratio_fd_mb": a.electron_cm3 / b.electron_cm3,
                "hole_ratio_fd_mb": a.hole_cm3 / b.hole_cm3,
                "scope_notes": " | ".join(a.scope_notes),
            }
            for density, a, b in zip(donors, fd, mb, strict=True)
        ],
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    figures, data = args.output / "figures", args.output / "data"
    figures.mkdir(parents=True, exist_ok=True)
    data.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.15,
            "lines.linewidth": 2,
            "svg.hashsalt": "carrier-statistics-v0.1",
        }
    )
    # Formal extrapolations are labeled in figures and CSV metadata; avoid repeated warnings.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelScopeWarning)
        carrier_regimes(figures, data)
        approximation_limits(figures, data)
        degeneracy_shift(figures, data)
    print(f"Saved three model figures (PNG/SVG) and four CSV tables to {args.output.resolve()}")
    print("High-temperature/high-doping formal extrapolations are marked. No experimental data.")


if __name__ == "__main__":
    main()
