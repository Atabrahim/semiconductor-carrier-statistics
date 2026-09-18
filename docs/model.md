# Scientific model and implementation

## Problem and scope

Given temperature and active substitutional donor/acceptor concentrations, find the equilibrium Fermi level and the populations of mobile carriers and charged dopants in uniform bulk silicon. This is a zero-dimensional equilibrium calculation. There is no spatial boundary-value problem: no device length, electrodes, applied voltage, surface charge, or electric field is supplied.

The crystal is unstrained. Electrons, holes and dopants share a temperature and equilibrium chemical potential. A continuous bulk density of states is assumed; finite-volume discreteness, confinement and nonequilibrium quasi-Fermi levels are excluded. The root interval is an algorithmic bracket, not a physical boundary condition.

The model combines analytical state counting, numerical integration and root finding, and semi-empirical material parameters. No parameters are fitted to project results.

## Units and parameter provenance

| Parameter | Value | Units / condition | Provenance |
|---|---:|---|---|
| Zero-temperature gap parameter | 1.1695 | eV; silicon X-valley gap model | Palankovski Table 3.13 |
| Varshni coefficient α | 4.73 × 10⁻⁴ | eV K⁻¹ | Same table |
| Varshni coefficient β | 636 | K | Same table |
| Electron transverse mass / free electron mass | 0.19 | Dimensionless; fixed band-edge approximation | Palankovski Table 3.16 |
| Electron longitudinal mass / free electron mass | 0.98 | Dimensionless; same approximation | Same table |
| Heavy-hole mass / free electron mass | 0.49 | Dimensionless | Same table |
| Light-hole mass / free electron mass | 0.16 | Dimensionless | Same table |
| Equivalent conduction valleys | 6 | Dimensionless; unstrained bulk Si | Palankovski Table 3.19 |
| Phosphorus donor binding energy ΔE_D | 0.04559 | eV; isolated bulk donor | Accepted value reported by Smith et al. |
| Boron acceptor binding energy ΔE_A | 0.045 | eV; isolated bulk acceptor approximation | Meng et al., ionization-energy discussion |
| Donor occupation degeneracy g_D | 2 | Neutral donor spin multiplicity | Single-level partition function below |
| Acceptor occupation degeneracy g_A | 4 | Unstrained bulk acceptor ground state | Mol et al., introduction |

The source masses are held constant **as a deliberate approximation**; the temperature-dependent mass models on the same source page are not implemented. Binding energies are held fixed relative to the band edges. These are reference inputs, not an experimental dataset or a declaration of precision over the entire input range.

SI constants $k_B$, $h$, $m_0$, and $q$ come from `scipy.constants`: Boltzmann's constant (J K⁻¹), Planck's constant (J s), free-electron mass (kg), and elementary charge (C). For eV-based energies use $k_{B,\mathrm{eV}}=k_B/q$. Densities in m⁻³ are multiplied by $10^{-6}$ to obtain cm⁻³. The validation report records the SciPy version.

## Temperature-dependent gap

$$E_g(T)=E_g(0)-\frac{\alpha T^2}{T+\beta}.$$

$T$ is temperature in K; $E_g$ and $E_g(0)$ are in eV; $\alpha$ is in eV K⁻¹; $\beta$ is in K. This is the semi-empirical Varshni model. We set the valence-band reference $E_V=0$ eV and conduction edge $E_C=E_g(T)$. Temperature comparisons use that reference, not a calculated absolute vacuum alignment.

## Effective densities of states

The single-valley electron DOS mass and combined hole DOS mass are

$$m_{n,\mathrm{DOS}}=(m_t^2m_l)^{1/3},\qquad
m_{p,\mathrm{DOS}}=(m_{hh}^{3/2}+m_{lh}^{3/2})^{2/3}.$$

$m_t,m_l$ are transverse and longitudinal electron masses; $m_{hh},m_{lh}$ are heavy- and light-hole masses. All masses in these physical equations are in kg. Code stores dimensionless ratios and multiplies by $m_0$ in the DOS prefactor.

$$N_C=2M_C\left(\frac{2\pi m_{n,\mathrm{DOS}}k_BT}{h^2}\right)^{3/2},\qquad
N_V=2\left(\frac{2\pi m_{p,\mathrm{DOS}}k_BT}{h^2}\right)^{3/2}.$$

$M_C=6$ is conduction-valley multiplicity; the factor 2 counts spin. $\pi$ is the circle constant. These expressions give m⁻³ before unit conversion. Constant masses imply $N_C,N_V\propto T^{3/2}$. The split-off band and nonparabolicity are omitted. DOS masses are not transport masses.

## Carrier populations

$$F_{1/2}(\eta)=\frac{2}{\sqrt\pi}\int_0^\infty
\frac{\sqrt t}{1+\exp(t-\eta)}\,dt.$$

Here $t$ is a **dimensionless energy integration variable**, not time; $\eta$ is reduced Fermi energy; $F_{1/2}$ is dimensionless. Normalization follows NIST DLMF 25.12.14. Mixing normalized and unnormalized integral conventions causes a factor error.

$$\eta_n=\frac{E_F-E_C}{k_{B,\mathrm{eV}}T},\qquad
\eta_p=\frac{E_V-E_F}{k_{B,\mathrm{eV}}T},\qquad
n=N_CF_{1/2}(\eta_n),\quad p=N_VF_{1/2}(\eta_p).$$

$E_F$ is the chemical potential in eV. $n,p$ are electron/hole concentrations in cm⁻³ when $N_C,N_V$ use cm⁻³. Holes are unoccupied valence states; their reduced-energy sign differs from that for electrons.

Boltzmann statistics replace $F_{1/2}(\eta)$ by $\exp(\eta)$ when occupation is dilute. In this limit,

$$n_i=\sqrt{N_CN_V}\exp\!\left(-\frac{E_g}{2k_{B,\mathrm{eV}}T}\right),\qquad np=n_i^2.$$

$n_i$ is the analytical intrinsic carrier concentration for the **Boltzmann model with this parameter set**. The product relation does not generally hold under degeneracy. Incomplete dopant ionization alone does not invalidate Boltzmann mass action.

## Dopant occupation and neutrality

$$E_D=E_C-\Delta E_D,\qquad E_A=E_V+\Delta E_A,$$

$$N_D^+=\frac{N_D}{1+g_D\exp[(E_F-E_D)/(k_{B,\mathrm{eV}}T)]},$$

$$N_A^-=\frac{N_A}{1+g_A\exp[(E_A-E_F)/(k_{B,\mathrm{eV}}T)]}.$$

$E_D,E_A$ are dopant levels in eV; $\Delta E_D,\Delta E_A$ are positive binding energies in eV. $N_D,N_A$ are active dopant concentrations in cm⁻³; superscripts label charged populations. $g_D,g_A$ are dimensionless multiplicity ratios **in the equations as written**. Other conventions can use reciprocal factors.

A neutral donor has a bound electron with two spin states; the ionized empty donor has one state, giving the donor factor 2. The bulk acceptor neutral ground state has four bound-hole states. Strain or an interface can split these states, so the factor is not universal for nanostructures.

Complete ionization separately sets $N_D^+=N_D$ and $N_A^-=N_A$. It is not implied by selecting Fermi–Dirac statistics.

$$n+N_A^-=p+N_D^+.$$

This is local bulk charge neutrality after canceling the common elementary-charge factor. It determines $E_F$. Increasing $E_F$ increases negative charge and reduces positive charge; the unscaled imbalance is monotonic, yielding a unique root for this model.

## Stable implementation

1. Compute $E_g,N_C,N_V$ once per input condition.
2. For each trial $E_F$, evaluate carrier and neutral-dopant log densities.
3. Evaluate a scaled charge residual and solve its sign change with `scipy.optimize.brentq`.
4. Return populations, solver iterations, residual and scope notes.

In the integral, substitute $t=u^2$, with dimensionless $u$, to smooth the endpoint. For negative $\eta$, factor $\exp(\eta)$ outside quadrature and retain its logarithm. For positive $\eta$, split both sides of the Fermi edge at $t=\max(\eta-50,0),\eta,\eta+50$. The omitted upper tail is exponentially smaller than the requested tolerances over the supported interval. The cutoff 50 is measured in thermal-energy units; it is not a material parameter.

The integral supports $-10^4\le\eta\le10^4$. `fermi_half` may underflow far into the negative tail; `log_fermi_half` retains the log value and is used by the solver. The quadrature relative tolerance is $10^{-10}$; absolute tolerance on scaled integrals is $10^{-13}$. The root absolute energy tolerance is $10^{-12}$ eV.

Define neutral populations $N_D^0=N_D-N_D^+$, $N_A^0=N_A-N_A^-$ and net doping $\Delta N=N_D-N_A$. To avoid subtracting large charged populations during compensation, evaluate

$$R=n-p+N_D^0-N_A^0-\Delta N,\qquad
S=\max(n,p,N_D^0,N_A^0,|\Delta N|).$$

$R$ and positive scale $S$ have units cm⁻³. The returned residual is dimensionless $R/S$. Paired differences use `expm1`; occupation factors use `logaddexp`. At exact complete-ionization compensation, the neutral populations and $\Delta N$ vanish, so even tiny intrinsic carriers determine the root. Normalizing by total chemical doping could misleadingly accept a broad range of Fermi energies.

## Independent analytical check

For Boltzmann carriers and complete ionization, solve $n-p=\Delta N$ and $np=n_i^2$. The discriminant is **$\Delta N^2+4n_i^2$**, with a plus sign. Calculate

$$c_{\mathrm{major}}=\frac{|\Delta N|+\sqrt{\Delta N^2+4n_i^2}}2,
\qquad c_{\mathrm{minor}}=\frac{n_i^2}{c_{\mathrm{major}}}.$$

$c_{\mathrm{major}},c_{\mathrm{minor}}$ are densities in cm⁻³. They are $n,p$ for positive net doping, and $p,n$ for negative net doping. This stable expression avoids subtracting almost equal values and is implemented separately from the numerical solver.

## Interpretation limits

Inputs outside 40–1000 K or 0–10²⁰ cm⁻³ per dopant species are rejected. Within that numerical domain, total doping above 10¹⁷ cm⁻³ or temperatures above 500 K produce `ModelScopeWarning` and result metadata. These conservative thresholds are project policies, not measured transition points. Results inside them still lack an experimental accuracy guarantee.

Band-gap narrowing, impurity bands, screening, dopant clustering, density-dependent binding energies, temperature-dependent masses and transport are absent. Low-temperature equilibration may be slow in a real sample. A converged solution cannot repair missing physics. High-doping comparisons deliberately isolate the mathematical statistics effect rather than validate either model for those conditions.

## Sources

- [NIST DLMF 25.12.14](https://dlmf.nist.gov/25.12): normalized integral.
- Palankovski (2000), [band-gap parameters](https://www.iue.tuwien.ac.at/phd/palankovski/node37.html), [masses](https://www.iue.tuwien.ac.at/phd/palankovski/node40.html), [DOS and multiplicity](https://www.iue.tuwien.ac.at/phd/palankovski/node41.html).
- Smith et al. (2016), [arXiv:1612.00569](https://arxiv.org/abs/1612.00569): accepted phosphorus binding energy 45.59 meV; their calculated 41 meV is **not** used here.
- Meng et al. (2023), [arXiv:2303.01564](https://arxiv.org/abs/2303.01564): boron ionization energy $I_p=45$ meV, used only as an input parameter.
- Mol et al. (2015), [arXiv:1501.05669](https://arxiv.org/abs/1501.05669): fourfold bulk acceptor degeneracy and its lifting near interfaces.
- [SciPy constants](https://docs.scipy.org/doc/scipy/reference/constants.html): constants from the installed SciPy release.

References checked on 18 September 2026. No source code, figures or measurement datasets were copied from these works.
