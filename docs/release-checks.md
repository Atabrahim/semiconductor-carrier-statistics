# v0.1.0 release checks

Release review: 19 September 2026. These checks verify the implementation and
distribution of the documented idealized silicon model. They do not establish
experimental accuracy.

## Distribution and clean installation

- `python -m build` produced an sdist and a platform-independent Python wheel.
  The wheel was built from the sdist, with no build warnings after updating the
  MIT license metadata to an SPDX expression.
- The source archive contains the model documentation, examples, tests, generated
  data and PNG/SVG figures. The wheel contains the reusable package, CLI entry
  point, dependency metadata and license.
- A fresh virtual environment installed the wheel with its `plots` and `dev`
  extras. This environment had no editable installation or system-site packages.
- Outside the checkout, every exported API symbol imported successfully;
  representative material, integral and equilibrium calculations passed.
  Both the installed `carrier-stats` command and `python -m carrier_stats`
  produced valid JSON. The imported module path was inside `site-packages`.
- `python -m pip check` reported no broken requirements.

Clean-install environment: Python 3.12.14, NumPy 2.5.3, SciPy 1.18.1,
Matplotlib 3.11.2, pytest 9.1.1 and Ruff 0.16.8. This is a second dependency
combination beyond the original [reference verification](validation.md).

## Automated and documented checks

- Full suite: **130 passed** from the wheel installation, with warnings treated
  as errors except explicitly expected physical-scope warnings.
- Both `ruff check` and `ruff format --check` passed.
- All README Bash commands and the Python API snippet executed successfully
  from an extracted source archive. The documented editable setup was checked
  in a separate environment, preserving the clean wheel installation.
- The validation report passed all checks across **192 equilibrium cases**.
  Its maximum reduced charge residual was **5.96 × 10⁻¹¹** and its maximum
  relative discrepancy from the analytical Boltzmann quadratic was
  **1.23 × 10⁻¹¹**. Both agree with the original report to the shown precision.
- Relative links in the README and documentation resolve to existing files.
- Tracked files exclude environments, caches, credentials, temporary outputs
  and distribution build directories.

[GitHub Actions](https://github.com/Atabrahim/semiconductor-carrier-statistics/actions/workflows/ci.yml)
builds and installs the wheel, exercises the API/CLI outside the checkout,
and runs the tests, Ruff and verification report. Its matrix covers Linux
on Python 3.11–3.13 and Windows on Python 3.12. The Windows job also executes
the PowerShell virtual-environment activation command documented in the README.
Python 3.12 jobs execute the figure-generation example. Consult the run for
the release commit or tag for its actual status.

## Visual and physical review

All three PNG figures and all three independently rendered SVG exports were
inspected for readable text, units, legends, scales and clipping.

| Figure | Interpretation checked |
|---|---|
| Carrier regimes | Freeze-out, the extrinsic plateau and intrinsic excitation; K, cm⁻³ and eV labels; analytical intrinsic-density curve distinguished from numerical solutions; high-temperature extrapolation shaded. |
| Approximation limits | Ionization error at fixed doping is separate from statistical error at fixed Fermi energy; logarithmic concentration/error scales; 1%, 5% and 50% contours readable. |
| Degeneracy shift | Charge neutrality keeps majority density near the donor density under complete ionization, while Fermi level and minority density change; density ratios are dimensionless; high-doping ideal-model extrapolation shaded. |

The initial contour labels were positioned before the logarithmic axis scale
was applied, clipping one label. Applying the scale first corrected the layout.
No equations, material parameters or numerical data were adjusted during this
release review. A LaTeX spacing escape that GitHub displayed as an exclamation
mark was also removed from the rendered documentation.

## Scope of completion

The release is a reproducible equilibrium-statistics study with numerical
verification, not an experimentally calibrated silicon model or a transport
simulator. Temperature-dependent effective masses, dopant interactions,
band-gap narrowing and comparison with licensed measurements remain future
work. The [learning notes](learning.md) explain the concepts a maintainer should
be able to defend before claiming independent technical ownership.
