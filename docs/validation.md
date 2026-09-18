# Verification and scientific limits

The tests establish that the software solves its documented model. They do not establish experimental accuracy for real silicon. Published material parameters are **REFERENCE DATA**; plots and solver outputs are **NUMERICAL RESULTS**; closed-form benchmarks are **ANALYTICAL RESULTS**. No **EXPERIMENTAL DATA** are included.

## Reproduce the checks

From a checkout with the development dependencies installed:

```bash
python -m pytest -q
python examples/validate_model.py
python examples/build_figures.py
python -m ruff check .
python -m ruff format --check .
python -m build
```

The scripts write fresh output to ignored `outputs/`. The committed [validation.json](validation.json) records a reference run, including Python, NumPy and SciPy versions. Its report exits unsuccessfully if any reported error exceeds its tolerance. It is not a substitute for running the full test suite.

## Benchmarks

| Check | Independent expectation | Failure it can detect |
|---|---|---|
| Normalized Fermi integral at zero | $F_{1/2}(0)=(1-2^{-1/2})\zeta(3/2)$ | Wrong normalization or quadrature |
| Independent integration variable | Integrate in energy rather than its square root | Substitution/Jacobian errors |
| Dilute limit | $\log F_{1/2}(\eta)\to\eta$ | Tail underflow or sign error |
| Degenerate limit | Sommerfeld expansion through its first correction | Unresolved narrow Fermi edge |
| Quadrature convergence | Tightened tolerance changes the result negligibly | Poorly converged numerical integral |
| Complete-ionization Boltzmann solution | Independent stable quadratic | Incorrect equilibrium equation or root |
| Freeze-out for each dopant type | Independent single-dopant quadratic with negligible minority carriers | Binding-energy/degeneracy/sign errors |
| Charge and dopant balance | Positive/negative charge equality and neutral + ionized = total | Inconsistent occupation accounting |
| Exact compensation at low temperature | Tiny intrinsic populations still determine the root | False convergence caused by cancellation |
| Mass action in Boltzmann statistics | $np=n_i^2$, even for incomplete ionization | Inconsistent electron/hole statistics |
| Dimensional/scaling checks | Known gap value, DOS magnitude and $T^{3/2}$ dependence | eV/J or m⁻³/cm⁻³ mistakes |
| Interface behavior | Reject invalid inputs; CLI emits parseable JSON | Silent invalid states or unusable output |

Here $\zeta$ is the Riemann zeta function; $\eta$ is dimensionless reduced Fermi energy; $n,p,n_i$ are concentrations in cm⁻³. Other symbols follow [model.md](model.md).

The 130 pytest cases include parameterized combinations rather than 130 unrelated physical claims. Warnings are treated as errors, except explicitly expected model-scope warnings.

## Reference run

Environment: Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0. The separate report covers 192 combinations: six temperatures, eight donor/acceptor pairs, two statistics options and two ionization options.

| Metric | Measured error | Acceptance limit |
|---|---:|---:|
| Relative error in $F_{1/2}(0)$ | 1.11 × 10⁻¹⁶ | 2 × 10⁻¹² |
| Maximum relative quadrature change | 1.13 × 10⁻¹⁴ | 1 × 10⁻⁸ |
| Maximum reduced charge residual | 5.96 × 10⁻¹¹ | 2 × 10⁻⁹ |
| Maximum relative error against the Boltzmann quadratic | 1.23 × 10⁻¹¹ | 2 × 10⁻⁹ |

Residual normalization uses the reduced charge-balance scale described in the model document; it is not simply divided by the total chemical doping. Stress tests beyond the intended physical scope are included to test numerical behavior. Passing them does not validate the missing high-doping physics.

## A defect the tests caught

An early quadrature implementation split only at the Fermi edge. At reduced energy 10,000, the integrator missed the thin region just below that edge, producing roughly 10⁻⁴ relative error. Splitting at both sides of the edge fixed the defect; the extreme-degeneracy test remains as a regression check.

An initial qualitative freeze-out assertion also used an unsupported 2% threshold at 40 K. The independently checked model predicts 3.46365% for the documented donor density. The test now uses a 5% qualitative bound, while the quantitative freeze-out test uses an independently derived analytical quadratic. No material parameters or plotted data were changed to force agreement.

## What remains unvalidated

- Absolute agreement with temperature-dependent Hall measurements.
- Temperature-dependent masses, impurity-band formation and band-gap narrowing.
- Uncertainty in physical parameters; numerical tolerances are not experimental error bars.
- Cryogenic kinetics, spatial inhomogeneity, surfaces and transport.

A useful experimental comparison would need known active doping, compensation, temperature calibration, the Hall factor, units and a redistribution license. No unrelated dataset is used to suggest that comparison has already happened.
