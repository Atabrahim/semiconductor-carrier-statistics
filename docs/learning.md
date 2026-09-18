# Understand and defend the project

## What the software answers

For uniform silicon at a chosen temperature and active doping, how many mobile electrons and holes exist, where is the Fermi level, and which common approximations are accurate enough? This project studies equilibrium populations; it does not simulate a working diode or a fabrication process.

The implementation was developed with AI assistance. A portfolio claim should distinguish directing an assisted implementation from independently deriving, coding or reviewing it. The exercises below provide concrete ways to develop and demonstrate that ownership.

## Development milestones

| Milestone | What was built | What to understand |
|---|---|---|
| Physical model | Silicon gap, DOS and isolated dopant parameters with sources | Why a material parameter has conditions and a model range |
| Statistics | Normalized Fermi integral and Boltzmann comparison | State occupancy, integral normalization and the dilute limit |
| Equilibrium | One bracketed root with complete/incomplete ionization | Neutrality couples carrier populations and dopant occupation |
| Verification | Analytical checks, limiting cases and conservation tests | Passing code is weaker evidence than an independent benchmark |
| Presentation | Reproducible figures, tables, CLI and documentation | Show the scientific question, numerical result and model limitation together |

## A 30-second explanation

“This project calculates equilibrium carrier concentrations and the Fermi level in silicon. It compares Fermi–Dirac and Boltzmann statistics and includes incomplete dopant ionization, so it captures freeze-out within an isolated-dopant model. It uses numerical integration and a charge-neutrality root solver, with analytical checks and automated tests. The main result is an explanation of when simpler assumptions fail, including why degeneracy can shift the Fermi level even when doping fixes the majority-carrier density.”

When asked about authorship, explain the AI assistance and identify the parts you personally checked or extended. Do not claim unaided implementation or experimental validation.

## A two-minute explanation

Start with the question: donor atoms do not always supply one mobile electron, particularly at low temperature. Describe the inputs—temperature and active donor/acceptor concentrations—and outputs—electron/hole densities, dopant charge states and Fermi level.

Explain that the density of states tells us how many states are available and Fermi–Dirac statistics tell us how occupied they are. Dopants also have occupancy probabilities. All populations depend on the Fermi level, which is selected by charge neutrality.

Describe the numerical work: an adaptive integral and Brent's bracketed root finder. The computation uses eV and cm⁻³ consistently and uses logarithms to retain very small populations. Exact compensation required special care to avoid subtracting large, nearly equal charge terms.

Finish with evidence and limits: 130 automated tests plus a 192-case verification report, including an analytical quadratic and Fermi-integral identities. Fixed masses and isolated dopants limit quantitative predictions; high-doping calculations are explicitly formal comparisons.

## Technical explanation

1. Derive the three-dimensional DOS dependence proportional to the square root of energy above a parabolic band edge. Explain why anisotropic silicon conduction valleys need a DOS mass and valley multiplicity.
2. Integrate DOS × occupancy to obtain the normalized order-1/2 Fermi integral. Identify the thermal-energy substitution and its units.
3. Derive donor occupation from empty/occupied states, retaining the spin multiplicity. Explain the opposite sign for acceptors and the degeneracy convention.
4. Write charge neutrality. Explain why its unscaled imbalance is monotonic in Fermi energy, making a bracketed method appropriate.
5. Explain the substitution that smooths quadrature, logarithmic treatment of dilute tails, and compensation residual. `expm1(z)` accurately computes exp(z)−1 when z is small; `logaddexp(a,b)` computes log(exp(a)+exp(b)) without overflow.
6. Derive the complete-ionization Boltzmann quadratic and calculate the minority population from mass action to avoid cancellation.
7. Interpret two comparisons separately: statistics at fixed Fermi energy versus equilibrium at fixed chemical doping. Under complete ionization, the latter constrains the majority population while allowing the Fermi level and minority population to change.
8. Separate numerical error, parameter uncertainty and model inadequacy. A smaller solver tolerance only addresses the first.

## Vocabulary

| Term | Simple meaning and scientific meaning | Engineering use and example |
|---|---|---|
| Fermi level | Chemical potential controlling state occupancy; equilibrium occupancy is 1/2 at this energy | Determine band populations; it moves toward the conduction band in donor-doped silicon |
| Density of states | Number of available states per energy and volume | Integrate available states with occupancy; different effective masses change carrier populations |
| Effective mass | Band curvature expressed as a mass parameter | DOS masses count states; transport masses describe acceleration and need not be equal |
| Freeze-out | Many carriers remain bound to dopants at low temperature | Predict reduced free-carrier density; at 40 K this model ionizes about 3.46% of 10¹⁵ cm⁻³ phosphorus donors |
| Compensation | Donors and acceptors are both present and partially cancel electrically | Distinguish chemical doping from net free-carrier supply |
| Degenerate semiconductor | Fermi occupation cannot be replaced accurately by a dilute exponential | Select statistics near a band edge; Boltzmann density errors grow as reduced Fermi energy increases |
| Charge neutrality | Positive and negative bulk charge balance | Determine a uniform material's equilibrium Fermi level; this is not a substitute for Poisson's equation in a junction |
| Active dopant | Substitutional dopant capable of electronic ionization | Distinguish active density from implanted dose and from the fraction currently ionized |
| Residual | Remaining mismatch in an equation after a numerical solve | Diagnose convergence; a small residual does not establish that the physical model is sufficient |

## Interview questions and answer direction

1. **Why is intrinsic Fermi energy not necessarily at mid-gap?** Conduction and valence DOS differ.
2. **Why does the carrier density not equal doping at low temperature?** Dopant occupation, including its binding energy and multiplicity.
3. **Does incomplete ionization invalidate mass action?** Not by itself under equilibrium Boltzmann statistics; degeneracy changes the product relation.
4. **Why use Brent's method rather than Newton's method?** A reliable bracket and no derivative requirement; discuss its cost and convergence tradeoff.
5. **What fails in naive exact compensation?** Large charged-dopant populations can hide the small carrier imbalance in floating-point arithmetic.
6. **What caused the high-degeneracy integration defect?** An unresolved narrow region below the Fermi edge; show the regression test.
7. **Why do FD and Boltzmann majority densities almost agree at fixed high donor density?** Neutrality and assumed complete ionization constrain them.
8. **Can this predict a wafer's resistivity?** No mobility or transport model is included; Hall measurements also require the Hall factor.
9. **Why is the 300 K intrinsic density different from another reference?** Parameterization and fixed-mass assumptions; never adjust a plot to hide the difference.
10. **How would you extend this to a PN junction?** Introduce spatial electrostatics, boundary conditions and eventually transport; bulk neutrality no longer holds everywhere.

## Small ownership exercises

- Derive the minority-carrier quadratic on paper and explain the numerically stable alternative.
- Add an example for 50% compensation and interpret its low-temperature behavior. Keep model logic in the package.
- Use `dataclasses.replace` to perturb the donor binding energy by ±1 meV. Label it a sensitivity study, not a confidence interval.
- Add an independently reasoned test for that extension, run the suite, and commit it as `Add donor binding-energy sensitivity study`.

Use a small feature branch when exploring; inspect `git diff`, commit one coherent change, and publish after each substantial milestone. Do not claim a pushed change until the remote file or commit has been verified.

## Portfolio wording

**GitHub description:** Python models of equilibrium silicon carrier statistics, dopant ionization and approximation limits, with analytical verification and scientific figures.

**Accurate CV bullet:** Built an AI-assisted Python semiconductor equilibrium package using Fermi–Dirac integration and charge-neutrality root finding, with incomplete dopant ionization, analytical verification and 130 automated tests.

**LinkedIn/project description:** A reproducible Python study of carrier populations in silicon, comparing statistical and dopant-ionization assumptions across temperature and doping. Includes analytical benchmarks, numerical verification and documented limits; developed with AI assistance.

## Hiring signal and remaining gap

The project supplies concrete discussion material for scientific-software, computational-physics and device-modeling student roles: physical assumptions, dimensional reasoning, quadrature, nonlinear solves, tests and reproducibility. It does not demonstrate cleanroom practice, fabrication expertise, commercial TCAD proficiency or experimental measurement skills. Its value in an interview depends on being able to defend and modify the work.

Project 2—PN-Junction Electrostatics and Diode Diagnostics—will add spatial potential, electric field, junction boundary conditions and device-level analytical comparisons. That is a new capability beyond uniform equilibrium statistics.
