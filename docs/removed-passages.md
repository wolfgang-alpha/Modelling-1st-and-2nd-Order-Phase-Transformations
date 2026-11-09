# Passages removed from earlier versions of the paper

This file preserves text that was written for earlier drafts of

> W. Flachberger et al., *A novel methodology for modelling first and second order phase transformations* (arXiv:2411.16430)

but that did not make it into the current version (v3, 2026). The passages are
reproduced as written at the time (including the original notation: mole fraction
$x$, Lagrangian $\mathcal{L}$, Onsager coefficient $L$, vacancy generation rate
$\phi$); only the LaTeX was converted to Markdown. In v3 the mole fraction is
called $\phi$, the vacancy generation rate $\psi$, and the Lagrangian is the
Rayleighian $\mathcal{R}=\dot{\mathcal F}+\tfrac12\mathcal P$.

Sources and their dates:

| tag | source | date |
|---|---|---|
| **[Jun 2024]** | notebook `CH-limiting-cases.ipynb` (first draft) | June 2024 |
| **[Jul 2024]** | *Limiting cases of the Cahn-Hilliard equation – thermodynamic aspects, variational methods and applications* (tex export of the notebook) | July 2024 |
| **[Nov 2024]** | journal submission *A Novel Methodology for Modelling First and Second Order Phase Transformations* | November 2024 |
| **[arXiv v2]** | arXiv:2411.16430v2 | 2024/2025 |

Equation numbers in the text refer to the respective old version and are kept only
where the sentence needs them.

**Reinstated in v3 (November 2026).** After this file was compiled, four of the passages were
put back into the paper in condensed form: the Cahn–Hilliard quotation on $\kappa \le 0$ and
the "unstable region" (Section 2 below, now in Section 1.2 of the paper), the paragraph on mole
fraction, phase fraction and order parameter (Section 3.1, now in the Motivation), the remark on
affinity versus chemical potential (Section 3.5, now in Section 2.1, with the diffusion potential
of Larché and Cahn), and the facet (jump-penalty) form of the gradient energy for the
Raviart–Thomas/DG pair (Section 6, now equation (62) in Section 4.2). The remaining passages
stay out of the paper and are preserved here only.

---

## 1. Original abstracts and keywords

### 1.1 First draft [Jun 2024]

A new methodology for the mathematical description and modeling of first and second
order phase transformations is introduced. It will be shown that the methodology can
be related to certain limiting cases of the Cahn-Hilliard equation, namely the cases
of having (i) a convex molar free energy function and (ii) a convex molar free energy
function and no regularization. Case (ii) is widely regarded as plainly unstable.
However, by changing the variational approach and by additionally solving for rate
dependent variables, a stabilized method is obtained that can handle the lacking
regularization. While the specific numerical method that will be used to solve the
equations, namely a mixed finite element approach, is not entirely new (since a
related method has been employed to stabilize solutions of the Laplace equation, see
the FEniCS mixed-Poisson demo), its utilization in the field of diffusion and more
importantly diffusional phase transformations is also a novelty. Since the special
case (ii) leads to a first order phase transformation and therefore to what has by
some authors been loosely referred to as "sharp interface" model
[Flachberger2023][Svoboda2017][Svoboda2013], we must emphasize the difference to
sharp interface models which employ a different variational approach, clearly
separating the balance laws for bulk and interface (such as derived in
[Dreyer2015]).

*Keywords:* Linear nonequilibrium thermodynamics, Phase Field method, mixed Finite
Element method, Theorem of minimal entropy production, Thermodynamic Extremal
Principle

### 1.2 Journal submission [Nov 2024]

This paper introduces a novel methodology for the mathematical modelling of first
and second order phase transformations. It will be shown that this methodology can
be related to certain limiting cases of the Cahn-Hilliard equation, specifically the
cases of having (i) a convex molar free energy function and (ii) a convex molar free
energy function with no regularization. The latter case is commonly regarded as
unstable; however, by modifying the variational approach and solving for
rate-dependent variables, we obtain a stabilized method capable of handling the
missing regularization.

While the specific numerical method used to solve the equations (a mixed finite
element approach) has been previously employed in related contexts (e.g., to
stabilize solutions of the Laplace equation), its application to diffusion and
diffusional phase transformations is novel. Moreover, our approach differs from
traditional "sharp interface" models, which typically employ a different variational
approach and clearly separate balance laws for bulk and interface.

### 1.3 arXiv v2

There are various methods for modeling phase transformations in materials science,
including general classes of phase-field methods [biner2017programming] and reactive
diffusion methodologies [svoboda2013new], which most importantly differ in their
treatment of interface energy. These methodologies appear mutually exclusive since
the respective numerical schemes only allow for their primary use case. To address
this issue, a novel methodology for modeling phase transformations in multi-phase,
multi-component systems, with particular emphasis on applications in materials
science and the study of substitutional alloys is introduced. The fundamental role
of interface energy in the evolution of a material's morphology will be studied by
example of binary and ternary systems. Allowing full control over the interface
energy quantity enables more detailed investigations and bridges the gaps between
known methods. We prove the thermodynamic consistency of the derived method and
discuss several use cases, such as vacancy-mediated diffusion. Furthermore a scheme
for relating Onsager and Diffusion coefficients is proposed, which allows us to
study the intricate coupling that is observed in multi-component systems. We hope
to contribute to the development of new mathematical tools for modeling complex
phase transformations in materials science.

---

## 2. The Cahn–Hilliard ansatz and the convex hull (introduction) [Nov 2024]

In their famous publication from 1958 [cahn1958free], Cahn and Hilliard introduced a
new way to formulate the Helmholtz energy of a binary alloy. Their starting point was
that the total Helmholtz energy of a volume $V$ of nonuniform composition can be
expressed as a functional of a Taylor series expansion of the mole fraction $x$,
where $\Omega$ denotes the molar volume of the alloy. For a stress-free system with
equal partial molar volumes of the components, it can be considered constant. For a
cubic crystal or isotropic medium, the functional can be simplified to include only
two summands: the molar free energy $f_0(x)$ of a uniform composition system, and a
"gradient energy" term $\kappa |\nabla x|^2$. The derivation involves neglecting
higher-order terms and applying the divergence theorem.

For $f_0(x)$, a double-well potential is typically employed, although this choice is
often without justification in the literature. However, the assumption of an
interface energy is implicit, as it implies that the mixed state ($x=0.5$) has a
higher free energy than a linear mixture of the individual phases in their
equilibrium states ($x_{\alpha}=0.25$, $x_{\beta}=0.75$). If no interface energy is
intended, the double-well potential can be replaced by its convex hull curve. This
represents the first limiting case (i). Using the convex hull curve instead of the
double-well potential introduces an intriguing change in the behavior of a system's
evolution, which has not been explored previously. This may be due to the fact that
finding an analytical expression for the convex hull can be more challenging than
working with the double-well. By employing the convex free energy, the system's
evolution will no longer focus on minimizing the surface area between the two
phases. Instead, it will solely aim to minimize the total free energy of the system
through phase growth, subject to boundary conditions and mass conservation. Another
intriguing property of using the convex hull curve is that it eliminates the need for
regularization, which typically requires $\kappa > 0$ for stability. Although Cahn
and Hilliard referred to the term $\kappa |\nabla x|^2$ as "gradient energy", it has
become common to call it a regularization or homogenization term due to its
smoothing effect on the interface. This methodology is also known as "Laplacian
smoothing" in fields like image processing, where similar techniques are used for
denoising and smoothing images. Interestingly, in the classical Cahn-Hilliard model,
which employs a typical double-well potential, this regularization term is necessary
for stability, meaning $\kappa > 0$. However, when using the convex counterpart of
the double-well, curiously, the homogenization term becomes unnecessary (i.e.,
$\kappa = 0$) for stability, provided an appropriate numerical scheme is employed.
This corresponds to limiting case (ii). As expected, the interface between the
phases will exhibit a very sharp profile in this case. Furthermore, we will
demonstrate that the overall quality of the system's evolution remains unaffected by
regularization, i.e., both cases (i) and (ii) evolve accordingly, differing only in
their representation, with (i) featuring a smooth interface and (ii) displaying a
sharp interface. Notably, the evolution will be significantly different from the
classical Cahn-Hilliard model that incorporates interface energy. Interestingly,
Cahn and Hilliard themselves noted that if $\kappa$ is set to zero, the interface
becomes sharp [cahn1958free]:

> *We can see no reason why $\kappa$ should change sign or vanish in the unstable
> region, but we are unable to prove that it could not happen. If for some system
> there were a range over which $\kappa \leq 0$, then there would be a corresponding
> discontinuity in the interface profile and the treatment would have to be modified
> accordingly.*

The term "unstable region" refers to the unstable region of the Van der Waals
equation below the critical temperature $T_c$, where a phase cannot exist in a
homogeneous equilibrium state (where particles must decide whether to form a gas or
a liquid, and hence an interface). The fact that Cahn and Hilliard hint at the Van
der Waals equation shows the generality and applicability of their concept. For now,
we simply want to motivate the special case of $\kappa = 0$ by noting that it can be
beneficial to have full control over the interface width when performing
multiphysics simulations.

---

## 3. The double-well potential: implications, interface energy, limitations [Nov 2024, Jul 2024]

### 3.1 Implications [Nov 2024]

The classical approach to introducing the Cahn-Hilliard equation involves the
example of spinodal decomposition in a two-phase binary solution (such as oil
droplets floating on water). Therefore, a simple assumption for the double-well
potential is made:

$$\Phi(x) = \tfrac{1}{4}x^2(1-x)^2 .$$

This is the usual form in which the double-well potential is presented in the
literature and it implies that the components have to be entirely separated from
each other (either $x=0$ or $x=1$) to reach equilibrium. From a physical point of
view this is a specific case but it is also convenient since the notions of mole
fraction, phase fraction and order parameter become identical. While the free energy
$\Phi$ is, due to its simplifying implications, beneficial for introducing the
concept of the Cahn-Hilliard equation or Phase Field Methods (PFM) in general, it
nevertheless caused some confusion between the PFM community (focused on
mathematics) and the reactive diffusion community (focused on materials science). It
led materials scientists to believe that there is no distinction between mole
fraction, phase fraction and order function in the Cahn-Hilliard equation, since
these concepts become blurred for this free energy. The misunderstanding can be
resolved by investigation of the proposed limiting cases of the Cahn-Hilliard
equation and their relation to reactive diffusion. This is because reactive
diffusion *is* the limiting case (ii) discussed in this paper. All models of Svoboda
[svoboda2006diffusion][svoboda2013new][svoboda2017incorporation], for example, can
be re-implemented to correspond formally to the Cahn-Hilliard equation. For
completeness, it needs to be mentioned that it is possible, even in materials
science, to employ the double-well $\Phi$ solely to introduce an order parameter to
describe different phases (see for example [miehe2014formulation]). Moreover, $\Phi$
is also helpful in modelling more interdisciplinary applications outside of
thermodynamics, for example in medical applications such as tumor growth
[acosta2024analysis]. A numerical problem that can occur, due to the fact that the
minima of $\Phi$ lie at the boundary of the permissible domain $[0,1]$, is the
phenomenon that during the iterations of the nonlinear solver the increment might
drift into the nonphysical realm of the domain (i.e. below zero or above one). This
obviously reduces the quality of a numerical solution. There exists however a new
numerical approach by Acosta-Soba et al. which eliminates this problem
[acosta2023upwind].

### 3.2 The effect of interface energy [Nov 2024]

To comprehend the significance of neglecting the interface energy as in the proposed
limiting cases (i) and (ii), it is essential to outline the thermodynamic situation
encountered in the simple case of the double-well potential $\Phi$. The striving for
total separation of the components as imposed is rather unusual in materials
science. In an alloy for example, it would be expected that the stationary (minimum)
points of the double-well are located within the $[0,1]$ domain, for the mole
fraction should obtain values corresponding to the mixing ratios of the solid
solution crystals to be described. This is also what Cahn and Hilliard primarily
intended to model in their famous publication and therefore gave no explicit
equation for the double-well but introduced it purely graphically. However, the case
of $\Phi$ is also of physical relevance, e.g. for modelling gas bubbles in a liquid
or for modelling the hydrophobic effect that is observed when mixing water and oil.
The physical example of mixing water and oil, or, more generally water and
hydrocarbons, is also a good starting point for explaining the influence of the
interface energy on the system evolution. Most hydrocarbons have a negligible mixing
enthalpy and negative, high-magnitude mixing entropy in water. This gives rise to
the hydrophobic effect [silverstein1998real]. The polar $H_2O$ molecules and
nonpolar hydrocarbons arrange themselves in a specific way at the interface, leading
to a reduction in microstates per macrostate, i.e., less entropy at the interface.
Therefore, the system strives to minimize free energy by maximizing entropy, and it
maximizes entropy by minimizing the area of contact (i.e., the interface). A key
result from the calculus of variations [walker2015shapes] reveals that globally
minimizing a surface area requires locally minimizing its curvature, resulting in
typical spherical shapes (or at least shapes with little curvature) that the
Cahn-Hilliard equation predicts. If the interface energy is neglected on the other
hand, the effect of spinodal decomposition can not be observed, and there will be no
tendency for the system to minimize its interface area. In conclusion, it can be
said that if the double-well $\Phi$ is used, we effectively impose the following
conditions on the system of interest:

1. The system strives for total separation of the components.
2. The relative height of the "bump" in the double-well potential (i.e. the
   interface energy) corresponds to the free energy of mixing of the components.
3. The standard free energies of the unmixed individual components are equal to one
   another.

While the standard Cahn-Hilliard equation may appear restrictive in its
application, we argue that this perception is reversed when considering the
equation's underlying structure. By modifying the free energy function, the CH
equation can be transformed into a highly versatile tool for modeling a wide range
of phase transformation phenomena. This flexibility arises from the ability to
tailor the double-well potential to specific problems, allowing the equation to
capture diverse characteristics. The double-well potential should therefore be
regarded as a placeholder that showcases what can be achieved with minimal effort.
The real challenge lies in selecting an appropriate free energy function to model
specific phenomena. A problem arises when attempting to assign names to these
modified equations. Is the CH equation still the CH equation if we replace the
double-well? Is it still a Phase Field Method? Replacing the double-well with a
quadratic function, for instance, would result in a simple diffusion equation,
albeit homogenized. Our goal is not to invent new nomenclature, but rather to
demonstrate how to adapt and extend the CH equation to tackle various real-world
phenomena.

It's worth noting that this idea is not novel; numerous publications and books have
already explored similar goals [emmerich2003diffuse][kim2016basic]. The novelty of
our publication lies in explicitly addressing the special case where no interface
energy exists (i.e., no "bump" and thus $\Delta G^{mix} \ngtr 0$). This seemingly
trivial consideration can, in fact, be a source of error for many applications.
*[Jul 2024]*

### 3.3 Thermodynamics of the double well: mixing, enthalpy and entropy [Jul 2024]

While it is a well-established fact that the overall entropy of an isolated system
must always increase over time, the CH equation initially appears to contradict
thermodynamics by creating a new, ordered arrangement from an initially fully mixed
solution. At first glance, this is surprising since entropy is loosely associated
with disorder, as exemplified by the irreversible mixing of two gases. However,
Kondepudi and Prigogine point out that in nonequilibrium thermodynamics, the
increase of total entropy can be accompanied by both creation of order and disorder
on different scales or in different quantities within a system (e.g., spreading heat
and temperature while achieving order in the arrangement of matter). For examples,
see [kondepudi2014modern].

The direction of a chemical reaction (either reactants → products or products →
reactants) is, for an equilibrium thermodynamic system at constant pressure and
temperature, determined by the temperature as well as the enthalpy and entropy
changes associated with the reaction,

$$\Delta G = \Delta H - T\,\Delta S .$$

A spontaneous chemical reaction is therefore either enthalpically, entropically, or
both enthalpically and entropically favorable, but not neither. The change in
enthalpy represents the actual change in energy of the matter, while the change in
entropy reflects the change in probability for the matter to be in a certain state.
However, this concept is not limited solely to chemical reactions, where molecules
exchange atoms to form different molecules. It also applies to the arrangement of
molecules in solutions or atoms in alloys.

Due to the double-well potential, the Cahn-Hilliard problem makes clear assumptions
about the overall Gibbs energy change encountered when mixing or separating the two
components, namely $\Delta G^{mix}$. Since the evolution tends towards phase
separation, it is evident that $\Delta G^{mix}>0$ holds. As a result, the mixed
solution can be viewed as the "reactants" and the unmixed phases as the "products"
of a chemical reaction, which explains the fourth-order form of the double-well
potential.

The double-well potential has three stationary points: two minima ($x=0$ and $x=1$),
which correspond to the Gibbs energies of the individual, unmixed phases; and one
extremum ($x=\tfrac12$), which corresponds to the mixed state. The difference in
Gibbs energies between the individual and mixed states must be equal to the Gibbs
energy change due to mixing:

$$f_0(\tfrac{1}{2}) - f_0(0) = \Delta G^{mix} .$$

This particular choice of the double-well potential also implies that the Gibbs
energies of the individual components are equal to one another.

A crucial aspect that is often overlooked is the influence of temperature changes
due to phase separation. While it may be acceptable to neglect this influence in
many cases, we will explore some implications arising from the Cahn-Hilliard
equation. As the Cahn-Hilliard equation describes the isothermal and linear regime
of nonequilibrium thermodynamics, we can only accurately describe systems where no
heat is released during phase separation, meaning that $\Delta H^{mix} = 0$ (an
isenthalpic phase separation). This limitation arises from the fact that a heat
release would introduce inhomogeneities in both temperature distribution and Gibbs
energies of the components. In such cases, the double-well potential would need to
be temperature-dependent itself.

There are certain hydrocarbons that exhibit the unique feature of having zero
mixing enthalpy when dissolved in water, with hexane ($C_6H_{14}$) being one
example. However, most hydrocarbons have a negligible mixing enthalpy
($\Delta H^{mix} \approx 0$) and negative, high-magnitude mixing entropy in water,
which gives rise to the hydrophobic effect [silverstein1998real].

In conclusion, if the double-well potential is used, we effectively impose the
following conditions on the system of interest:

1. $\Delta G^{mix} > 0$ (ideally, $\Delta H^{mix} = 0$ and $\Delta S^{mix} < 0$);
2. the "bump" in the double-well potential has height $\Delta G^{mix}$ relative to
   the local minima;
3. the standard Gibbs energies of the unmixed individual components are equal to one
   another;
4. if the two phases are fluids, we assume equal densities and neglect possible
   volume forces and effects due to fluid flow.

These conditions are obviously quite specific and do not resemble the problem
setting of any phase transformation phenomenon in physics. A scientist who employs
the CH equation to model a kind of phase transformation where few or none of the
above conditions are satisfied might be drawn to the conclusion that the CH is just
a phenomenological description that can at best deliver qualitative assessment of a
phase transformation. The authors' opinion on the matter is that the Cahn-Hilliard
equation is in fact a very accurate description of a specific case and that it
should not merely be understood as the description of this case, but as a way of
introducing a whole set of ideas to a scientist. *[Jun 2024]*

### 3.4 Limitations and thermodynamic consistency (Landau theory) [Nov 2024]

The Cahn-Hilliard model is often linked to Landau's theory of first- and
second-order phase transformations near the critical point (of the Van der Waals
equation). In this context, the CH model corresponds to modeling behavior below the
critical point. The special case (ii) corresponds specifically to the behavior at
the critical point, where two phases can coexist but without any driving force for
separation. One can envision a thermodynamic state where air doesn't form bubbles in
water but rather a non-gaseous and non-aqueous solution – indeed a very critical
phenomenon. However, Landau's theory has been criticized for its inability to
reproduce experimental measurement data, leading to its replacement by the
Renormalization Group (RG) theory. Despite this, Landau's theory still retains
relevance in the field of FEM simulation and continuum modeling. As a
macroscopic/thermodynamic theory, it provides a fundamental framework for
understanding phase transitions at the mesoscopic scale. In contrast, RG theory
combines aspects of both thermodynamics and atomistics, making it challenging to
reconcile with continuum modeling.

Additionally it is worth noting that the main point of the critique of Landau's
theory is concerned with interfaces of gas and liquid phases and not with
solid-solid interfaces [kondepudi2014modern].

The CH model is, however, a thermodynamically consistent model. Therefore, it does
not deserve to be merely called phenomenological. In fact, one could argue that this
criticism can be extended further: if we consider thermodynamics itself as
phenomenological (a subject of ongoing debate), then what does it mean to label a
specific model within thermodynamics as phenomenological? By recognizing the CH
model's thermodynamic consistency, we can avoid unnecessary classification and focus
on its fundamental role in understanding phase transitions.

### 3.5 On the affinity versus the chemical potential [Nov 2024]

The debate surrounding thermodynamic consistency in Phase Field Models (PFMs) has
been ongoing, with one key issue being the common mistake of interpreting $\mu$, the
functional derivative of the free energy with respect to the mole fraction, as the
chemical potential. In fact, within the Cahn-Hilliard equation, $\mu$ is an affinity
rather than a chemical potential. This misconception can lead to difficulties in
relating PFMs to real-world data and parameters. To clarify this point and avoid
further confusion, we adopt standard notation from materials science, where $\mu$
indeed represents the chemical potential, and use $\bar{\mu}$ to denote the
parameter of the CH equation, which is often mistakenly referred to as the chemical
potential.

---

## 4. Remarks on the Cahn–Hilliard equation in concentration form [Nov 2024, arXiv v2]

Equations (32) and (35) [the flux and affinity derived from the Lagrangian] are
reminiscent of equations (4) and (5) [the classical Onsager flux law and the
definition of the chemical potential], however they are not the same. As soon as the
constraint $\sum_i \underline{j}_i = 0$ is imposed on the classical equations the
system becomes overdetermined and may be reformulated to satisfy (32) and (35).
Considered on their own, both representations are however thermodynamically
consistent. It also must be noted that the Onsager coefficients $L_{ik}$ must be
defined differently for the two cases. However, $L_{ik}$ must always be a symmetric
and positive semi-definite matrix, as postulated by Lars Onsager. It is therefore
proven that reactive diffusion models and the Cahn-Hilliard model are inherently the
same, just distinguishable by their definitions of the molar free energy $f_m$.

By choosing $\Omega = 1$, $\tfrac{L}{T}=M$ and
$f_m = f_0(c) + \tfrac{\kappa}{2} |\nabla c|^2$ the evolution equation can also be
presented in terms of the concentration $c$, which is the usual way the
Cahn-Hilliard equation is presented in the literature. This equation also hints at
the problem that is encountered when a convex molar free energy is employed; the
function $f_m$ has to be at least $C^3$ continuous if the equation is to be solved
directly. This gives rise to the assumption that the relation between reactive
diffusion and Phase Field Methods was often overlooked because the formalisms seem
incompatible.

*[arXiv v2 version of the same remark:]* The differential operators in the
Cahn-Hilliard equation demand at least $C^3$ continuity of the function $f_m(x)$.
While this continuity requirement can be lowered to $C^2$ continuity (this is also
necessary since the convex hull curve $f_0^{**}$ is only $C^2$-continuous) by
deriving variational formulations, the problem remains that if $x$ itself is
discontinuous, $f_m(x)$ can also become discontinuous, i.e. $C^0$-continuous. In
practice, case (ii) should therefore be used with caution, and regularization
should generally be applied to restore higher-order continuity of the mole fraction
$x$.

---

## 5. Variational principles in nonequilibrium thermodynamics (history of the TEP) [Nov 2024]

The Lagrangian can also be derived using the Thermodynamic Extremal Principle (TEP)
[svoboda2023note]. This approach is particularly useful for systems defined by a
total free energy functional $\mathcal{F}$. The TEP is formulated as follows:

$$\delta\big(\dot{\mathcal{F}} + \tfrac{T}{2}\,\mathcal{P}\big) = 0 ,$$

where $\mathcal{P}$ represents the rate of total entropy production. For a detailed
derivation, see references such as [emmerich2003diffuse] or
[svoboda2005application]. The existence of the TEP was initially postulated by
Kelvin (the first draft said: by Lord Rayleigh), while Lars Onsager was the first to
formulate the relations of linear thermodynamics in this manner. However, it was
Prigogine who presented a general formulation and demonstrated the validity of the
principle [kondepudi2014modern]. Prigogine referred to this principle as the
"Theorem of Minimum Entropy Production" because he found that the rate of entropy
production reaches a minimum value at the nonequilibrium stationary state. As noted
in [svoboda2023note], this title can however be misleading since the rate of
entropy production actually reaches its constrained maximum during the evolution of
thermodynamic systems that transition from a general nonequilibrium state to either
equilibrium or a stationary nonequilibrium state. The TEP was formulated in
[svoboda2023note] with respect to the fluxes $\underline{j}$ to demonstrate that the
evolution of a general nonequilibrium state corresponds to maximum dissipation:

$$\delta_{\underline{j}} \int_V \Big( \underline{j} \cdot \nabla \bar{\mu} + \frac{T}{2L} | \underline{j} |^2 \Big) \, dV = 0
\quad\Longrightarrow\quad
\underline{j} = - \frac{L}{T} \nabla \bar{\mu} .$$

Emmerich [emmerich2003diffuse] formulated the TEP in terms of the affinity
$\bar{\mu}$:

$$\delta_{\bar{\mu}} \int_V \Big( \bar{\mu}\, \dot{c} + \frac{L}{2T} \nabla \bar{\mu} \cdot \nabla \bar{\mu} \Big) \, dV = 0 .$$

Here the first term of the integrand, $\bar{\mu}\dot{c}$, corresponds to the rate of
free energy and the second term, $\frac{L}{2T}\nabla\bar\mu\cdot\nabla\bar\mu$, to
the rate of entropy production. Variation with respect to the affinity $\bar{\mu}$
yields the well known evolution equation

$$\dot{c} + \nabla \cdot \Big(- \frac{L}{T} \nabla \bar{\mu}\Big) = 0 ,$$

with the common variational form for usual continuous Galerkin type numerical
schemes:

$$\int_V \dot{c}\, v \, dV = \int_V L\, \nabla v \cdot \nabla \bar{\mu}(c) \, dV .$$

While a simple approach like the one described above might seem appealing for
solving the Cahn-Hilliard (CH) equation, it may not always be favourable to employ.
A crucial consideration is that the variational form of this evolution equation will
only be stable for continuous Galerkin methods if $\bar{\mu}$ is at least $C^1$
continuous [flachberger2024numerical]. This can be a problem when dealing with
relaxed potentials, such as the convex hull curve of a double well potential. In
such cases, it is necessary to solve for the flux directly to avoid the
discontinuity of the solution. Furthermore, coupling diffusion phenomena with
mechanical models can be challenging since it is not immediately clear which part
of the functional corresponds to dissipative and energy terms. To overcome these
difficulties, we introduced the Lagrangian, providing a straightforward
computational methodology for solving the CH equation.

Here, the affinities $\bar{\mu}_k$ were used as Lagrangian multipliers to
incorporate the mass conservation constraint. There are also other examples where
physical quantities emerge as Lagrangian multipliers in fundamental variational
problems, see for example [flachberger2021investigating].

*[arXiv v2, footnote:]* In this context, "Lagrangian" refers to the general
framework of the calculus of variations, distinct from the Lagrangian of classical
mechanics, which is defined in terms of kinetic and potential energy.

---

## 6. The original mixed formulation with $R^1/D^0$ elements and the jump penalty [Nov 2024]

> In v3 the time discretization is imposed after the variation and the
> RT/DG pair appears as a cross-check (Section 4.2). The November 2024 version
> discretized the Lagrangian in time first and used RT/DG as the primary element pair.

For simulation purposes the Lagrangian can even be adapted into a discrete time
setting. This is particularly useful when performing time-dependent simulations in a
backward time discretization scheme. In this case, the variational formulation for
Ritz/Galerkin-type methods can be easily obtained by computing the first variation
of the functional with respect to all independent variables
$(\{x\}, \underline{j}, \bar{\mu})$. Furthermore, newer finite element tools like
FEniCSx or COMSOL even include automated differentiation functions, which allows for
rapid implementation of models and simpler coupling of many phenomena. For
simplicity, we will consider a binary system:

$$\int_V \mathcal{L}(x, \underline{j}, \bar{\mu}) \, dV :=
\int_V \Big( \frac{1}{\Omega} \frac{f_m(x)-f_m(x_t)}{\Delta t} + \frac{T}{2L} |\underline{j}|^2
+ \bar{\mu} \Big( \frac{1}{\Omega} \frac{x - x_t}{\Delta t} + \nabla \cdot \underline{j}\Big) \Big) \, dV .$$

*[arXiv v2 remark:]* Interestingly, the difference between total and partial time
derivatives in the discrete setting is negligible: the total time derivative of the
molar free energy is implemented as a difference quotient, just as the partial time
derivative of the mole fraction $x$.

Computing the variational derivative $\delta \int_V \mathcal{L}\, dV = 0$ yields

$$\frac{\partial}{\partial \varepsilon} \int_V \Big( \frac{1}{\Omega} \frac{f_m(x+\varepsilon \hat{x})-f_m(x_t)}{\Delta t}
+ \frac{T}{2L} | \underline{j}+\varepsilon \underline{\hat{j}} |^2
+ (\bar{\mu} + \varepsilon \hat{\bar{\mu}}) \Big( \frac{1}{\Omega} \frac{x+\varepsilon \hat{x} - x_t}{\Delta t}
+ \nabla \cdot (\underline{j}+\varepsilon \underline{\hat{j}}) \Big) \Big) \, dV \; \Big|_{\varepsilon=0} = 0 ,$$

$$\int_V \Big( \frac{1}{\Omega} \frac{1}{\Delta t} \frac{\partial f_m(x)}{\partial x} \hat{x}
+ \frac{T}{L} \underline{j} \cdot \underline{\hat{j}}
+ \hat{\bar{\mu}} \Big( \frac{1}{\Omega} \frac{x - x_t}{\Delta t} + \nabla \cdot \underline{j} \Big)
+ \bar{\mu} \Big( \frac{1}{\Omega} \frac{\hat{x}}{\Delta t} + \nabla \cdot \underline{\hat{j}} \Big) \Big) \, dV = 0
\quad \forall \; \hat{x},\underline{\hat{j}},\hat{\bar{\mu}} .$$

If the test functions are varied individually the individual variational
formulations are obtained:

$$0 = \int_V \Big( \frac{\partial f_m(x)}{\partial x} + \bar{\mu} \Big) \hat{x} \, dV \quad \forall \; \hat{x} \in D^0 ,$$

$$0 = \int_V \Big( \frac{T}{L} \underline{j} \cdot \underline{\hat{j}} + \bar{\mu} \, \nabla \cdot \underline{\hat{j}} \Big) \, dV \quad \forall \; \underline{\hat{j}} \in R^1 ,$$

$$0 = \int_V \Big( \frac{1}{\Omega} \frac{x - x_t}{\Delta t} + \nabla \cdot \underline{j} \Big) \hat{\bar{\mu}} \, dV \quad \forall \; \hat{\bar{\mu}} \in D^0 .$$

While the first form is straightforward, as it simply ensures that $\bar{\mu}$
represents the chemical affinity, the second and third are notable for their ability
to stabilize solutions of the Laplace equation. These forms were first introduced in
[brezzi1985two] and utilize a special mixed function space for their solutions. To
this end, we employ a discontinuous Galerkin finite element space of degree zero for
both the mole fraction and the affinity. The flux is typically represented by
interpolation functions of the $H(\mathrm{div})$ family, which are vector-valued and
rotation-free finite elements that are conforming to the divergence operator. This
ensures that the variational form strongly enforces the properties of the divergence
operator in the weak form, thereby conserving mass for transport-type equations.
This approximation is particularly beneficial when accurate results are required for
multi-component systems. The rotation-free property also eliminates unnecessary
degrees of freedom, as fluxes arise from the gradient of a potential field that is
inherently rotation-free. We denote the space of $H(\mathrm{div})$ conforming
elements of degree one as $R^1$ and the discontinuous Galerkin space of degree zero
as $D^0$. Due to this choice of function space a regularization term as in the CH
equation can be included by means of the jump operator $[\cdot]$. Altering the
first variational form accordingly yields:

$$0 = \int_V \Big( \frac{\partial f_m(x)}{\partial x} + \bar{\mu} \Big) \hat{x} \, dV
+ \int_{\partial \Omega_e} \kappa\, [x]\,[\hat{x}] \, dS \quad \forall \; \hat{x} \in D^0 .$$

Here $\partial \Omega_e$ denotes the boundaries of individual finite element cells.
An innovative aspect of this variational form lies in its treatment of the affinity,
which is not differentiated. By solving for the flux explicitly the second order
elliptic PDE is split into two first order PDEs. In [flachberger2024numerical], a
detailed analysis reveals why traditional variational forms for reactive diffusion
problems fail due to the continuity properties of the affinity, being just $C^1$
continuous if a convex hull is employed. By avoiding differentiation of the
affinity, our approach circumvents these issues and is able to solve the discussed
limiting cases. Note that while this mixed finite element method is due to its
properties tailored for this problem it is also possible to use standard
(continuous) Galerkin elements for all degrees of freedom.

*[Jun 2024 wording:]* The ingenious part of this variational form is that the
affinity is not differentiated, which is good since it is only $C^1$ continuous
with respect to the mole fraction (and the mole fraction can itself be discontinuous
over the domain if a sharp interface is assumed).

---

## 7. Variable molar volumes and Appendix B [Nov 2024]

> Removed already in arXiv v2. v3 assumes a constant molar volume throughout.

### 7.1 Variable molar volumes

So far, the molar volume was always considered constant, as in the Cahn-Hilliard
equation. The presented method is however also particularly useful for modeling the
behaviour of substitutional alloys with different partial molar volumes of the
components or different molar volumes of the phases. This is of special interest,
given that many alloys exhibit multiple crystal structures with potentially
different molar volumes, making it essential to account for these differences in
order to accurately capture their behavior. Therefore it is necessary to define the
dependence of the molar volume $\Omega$ on the composition in the Lagrangian.
Computing the variational derivative of the Lagrangian will then obviously yield
weak forms that are more complex than the three forms above. This is however no
problem as long as $\Omega(x)$ is at least $C^2$ continuous. Additionally the step
of applying the variational derivative can, as already pointed out, be automated.
Curiously, the variational forms above can also be employed if the molar volume
depends on the phase and a convex molar free energy $f_0^{**}(x)$ is used. The proof
for this is however mathematically involved and not relevant for all readers and
was therefore moved to Appendix B.

### 7.2 Appendix B: Special properties of the limiting cases for different molar volumes of the phases

Consider the Lagrangian of a binary system:

$$\mathcal{L}(x, \underline{j}, \bar{\mu}) := \frac{d}{dt} \Big(\frac{f_m}{\Omega}\Big)
+ \frac{T}{2L} |\underline{j}|^2
+ \bar{\mu} \Big( \frac{d}{dt} \Big(\frac{x}{\Omega}\Big) + \nabla \cdot \underline{j}\Big) .$$

This time a more general case is studied where not only the molar free energy $f_m$
but also the molar volume $\Omega$ are explicitly dependent on the mole fraction
$x$: $f_m := f_m(x)$, $\Omega := \Omega(x)$. Application of the Euler-Lagrange
equation with respect to $x$ yields the affinity $\bar{\mu}$:

$$\frac{\delta \mathcal{L}}{\delta x} := \frac{d}{dt} \Big( \frac{f_m' \Omega - f_m \Omega'}{\Omega^2} \Big)
- \bar{\mu} \frac{d}{dt} \Big( \frac{\Omega - x \Omega'}{\Omega^2} \Big) = 0
\quad\Longrightarrow\quad
\bar{\mu} = \frac{d \Big( \frac{f_m' \Omega - f_m \Omega'}{\Omega^2} \Big)}{d \Big( \frac{\Omega - x \Omega'}{\Omega^2} \Big)} .$$

As expected, the dependencies are now more involved. However, if $f_m$ is a linear
function of $x$, which is the case for limiting cases (i) and (ii) at least in the
interface region, we can express $f_m$ as $f_m = f_m' x$, where $f_m'$ is of course
a constant. This simplifies the expression:

$$\bar{\mu} = \frac{d \Big( f_m' \frac{ \Omega - x \Omega'}{\Omega^2} \Big)}{d \Big( \frac{\Omega - x \Omega'}{\Omega^2} \Big)} = f_m' .$$

The condition is therefore satisfied at the interface. For the bulk area another
consideration can be made: if $\Omega(x)$ is constructed by typical functions for
the phase fraction, that is functions that are constant over each specific phase and
have a smooth transition at the interface, we can make the following simplification
for non-interface, i.e. bulk, areas ($\Omega' = 0$):

$$\bar{\mu} = \frac{d \Big( \frac{f_m' \Omega}{\Omega^2} \Big)}{d \Big( \frac{\Omega}{\Omega^2} \Big)}
= \frac{d \big( f_m' / \Omega \big)}{d \big( 1 / \Omega \big)} = f_m' .$$

Therefore it is obvious that the variational forms of the constant-volume case can
be employed for the special cases (i) and (ii) if phases with different molar
volumes are to be modeled. For the general case, the full expression for $\bar\mu$
has to be employed. The variational framework however also allows the variational
derivatives to be automated for the general case.

---

## 8. Variable diffusion coefficients (long version) [Nov 2024, arXiv v2]

> v3 keeps a one-paragraph version without equation.

It is commonly observed that the components of an alloy exhibit different diffusion
coefficients in various phases. To account for this phenomenon, the diffusion
coefficients $D_i$ or, more generally, the Onsager coefficients $L_{ij}$ can be
chosen to be phase-dependent and expressed as a function of the mole fraction $x$.
However, at the fundamental level of the variational problem defined by the
Lagrangian, it is crucial to maintain constant values for the Onsager coefficients
when applying variations to the Lagrangian. This is due to the definition of the
Onsager coefficients, which, as is typical in thermodynamics, requires maintaining
other quantities as fixed when evaluating partial derivatives. An updated scaled
Lagrangian for the binary system could therefore be realized as:

$$\mathcal{L}(x, \underline{j}, \bar{\mu}) := \frac{f_m(x)-f_m(x_t)}{\Delta t}
+ \frac{k}{2\tilde{D}(x_t)} |\underline{j}|^2
+ \bar{\mu} \Big( \frac{x - x_t}{\Delta t} + \nabla \cdot \underline{j}\Big) .$$

Note that the inter-diffusion coefficient $\tilde{D}$ is made dependent on the mole
fraction of the "previous" time step $x_t$. This avoids variation of the dissipative
term with respect to $x$ when applying automated derivation of the variational
forms. Alternatively, if the variational form is directly implemented, it can still
be required that $\tilde{D}:=\tilde{D}(x)$, rather than deriving them automatically
from the Lagrangian. *[Nov 2024 wording:]* The Onsager coefficient $L$ can also be
chosen phase dependent and therefore as a function of the mole fraction $x$.
However, the value of the Onsager coefficient should be held constant when applying
the variation of the Lagrangian – otherwise one receives invalid evolution
equations.

---

## 9. The ternary model: original presentation [Nov 2024]

> The ternary section survives in v3 in rewritten form (Rayleighian with a
> rate-dependent constraint functional). This is the original presentation with
> the time-discretized constraint functional and the dissipation written with the
> coefficients $A_i$ of Svoboda's model.

The dissipative part of the problem is defined, according to
[svoboda2006diffusion], by the following functional:

$$\mathcal{P}[\underline{j}_0, \underline{j}_1, \phi] := \int_V \Big( \frac{1}{A_0} |\underline{j}_0|^2
+ \frac{1}{A_1} |\underline{j}_1|^2 + \frac{1}{A_2} |-\underline{j}_0-\underline{j}_1|^2
+ \frac{1}{A_{\phi}} \phi^2 \Big) \, dV .$$

The coefficients $A_i$ are in general dependent on the composition and the
diffusion coefficients. A definition can be found in [svoboda2006diffusion]. Note
that the flux of component 2 was substituted as
$\underline{j}_2 := - \underline{j}_0 - \underline{j}_1$. This not only constrains
the system and enforces the vacancy mechanism but also causes the resulting Onsager
coefficient matrix of the system to have non-diagonal entries. Therefore the fluxes
become coupled and a driving force of one component may also cause fluxes of other
components, as postulated by Onsager's reciprocal relations. The function $\phi$ is
characteristic for this very model and represents the rate of generation and
annihilation of vacancies. Since vacancies are a non-conserved component they can
emerge and vanish in the bulk so as to facilitate diffusion. Their evolution in time
is therefore also not determined by a pure conservation equation but rather a
transport-type equation. The evolution equations are considered in the constraint
functional $\mathcal{C}$:

$$\mathcal{C} := \int_V \Big( \bar{\mu}_0 \Big( \frac{1}{\Omega} \frac{x_0 - x_0^{t}}{\Delta t}
+ \nabla \cdot \underline{j}_0 - \frac{\phi}{\Omega} (1-x_0) \Big)
+ \bar{\mu}_1 \Big( \frac{1}{\Omega} \frac{x_1 - x_1^{t}}{\Delta t}
+ \nabla \cdot \underline{j}_1 + \frac{\phi}{\Omega} x_1 \Big) \Big) \, dV .$$

The final problem is defined by its Lagrangian
$\int_V \mathcal{L}\, dV := \mathcal{F} + \mathcal{P} + \mathcal{C}$ and
$\delta \int_V \mathcal{L}\, dV = 0$.

---

## 10. Results narrative of the November 2024 version (unscaled parameters)

> The simulations of the Nov 2024 version used $L = T = \Omega = 1$ and
> $\kappa = 5$ (1D), $\kappa = 100$ (interface-width study), $\kappa = 1$ (2D).
> v3 uses the scaled parameters $k = 2$ J/mol, $D = 2$ mm²/s, $\kappa = 2c$ that
> match the public notebooks.

To explore the capabilities of our method, we employed the open-source finite
element solver FEniCSx, which is continuously developed by the FEniCS community. The
software relies on several key modules, including dolfinx [BarattaEtal2023], basix
[ScroggsEtal2022][BasixJoss], and ufl [AlnaesEtal2014]. FEniCSx can be easily
integrated with Python 3 and was used in conjunction with matplotlib and pyvista to
generate all graphics. We start by exploring the one-dimensional case, where we
contrast solutions for the classical Cahn-Hilliard model and its special cases (i)
and (ii), all initialized with a common condition, namely a linearly rising
distribution of the mole fraction. As anticipated, the classical Cahn-Hilliard model
exhibits phase separation, characterized by a smooth transition between the phases.
In contrast, special case (ii) features a sharp, discontinuous interface due to the
absence of homogenization terms and the presence of a convex hull curve.
Interestingly, the solution is unstable as long as the interface is not infinitely
thin. By incorporating an added regularization, special case (i) yields a smooth
phase transition, similar to that observed for the classical Cahn-Hilliard model,
but with a wider phase transition region. The disparity in the interface thickness
of special case (i) and the Cahn-Hilliard equation can be attributed to the convex
free energy not imposing an energetic penalty on the interface. To prove this, we
have performed the calculation again for many different interface energies: the
interface energy unarguably impacts the width of the interface. In the classical
Cahn-Hilliard model, the "bump" in the double-well potential counteracts the
regularization term, resulting in a smooth yet locally restricted interface. The
interplay between these mechanisms drives the interface surface area towards
minimization in the multidimensional case.

In the two-dimensional study we consider a square domain of unit length, initialized
with a random pattern for the mole fractions. The binary system is composed of two
components that are entirely separated from each other at this stage. The blue
region represents 100 % component 2 and phase α, while the red region corresponds
to 100 % component 1 and phase β. For the classical Cahn-Hilliard model the mole
fraction approaches its equilibrium value for each present phase respectively, and
the model exhibits a tendency to minimize contact area between the two phases,
resulting in typical spherical shapes. In contrast, limiting cases (i) and (ii)
exhibit distinct patterns: there are regions where the interface exhibits higher
curvature, indicating that the model does not prioritize contact area minimization.
Nevertheless, the mole fraction reaches equilibrium values for each phase
respectively, as in the Cahn-Hilliard model. For limiting case (i), the overall
evolution of the structure corresponds to the evolution of limiting case (ii), with
the notable exception being the smooth representation of the interface. It is
important to note that the calculation time is significantly lower in limiting case
(ii), which is due to the fact that the sharp representation of the interface also
causes the matrices of the nonlinear solver to be more sparse.

Finally, we consider the ternary system with vacancies. The same initial condition
of the previous example is applied for component $x_1$. The initial condition for
the vacancies will be uniformly set to its equilibrium site fraction
$x_0 = x_0^{eq}$. Note that the dissipation term of the vacancy flux is, unlike
suggested in [svoboda2006diffusion], neglected (implying that
$A_0 \rightarrow \infty$) because it is assumed that moving a vacancy by itself
does not account for any dissipation. However, since diffusion in substitutional
alloys is primarily enabled by atoms switching sites with vacancies, we nevertheless
expect dissipation which is accounted for by the other terms, containing $A_1$ and
$A_2$. These coefficients are chosen as $A_1=2$, $A_2=1$ to trigger the so-called
Kirkendall effect. Due to this choice the component $x_1$ will obviously diffuse two
times faster than component $x_2$ ($= 1 - x_0 - x_1$). Thereby, the need to
counterbalance this disparity by a vacancy flux is created. The flux of vacancies
$\underline{j}_0$ will create an excess in vacancies in the area of the inclusions
and a lack of vacancies in the bulk regions, as expected. The overall evolution of
$x_1$ is not drastically changed, but rather slowed down by the involved vacancy
mechanism.

Original figure captions:

- *Figure 2 (1D Cahn-Hilliard):* the molar free energy $f_0(x)$ is employed and the
  homogenization parameter was set to $\kappa = 5$; all other parameters were set
  to one ($L = T = \Omega = 1$).
- *Figure 3 (1D case (ii)):* $f_0^{**}(x)$, $\kappa = 0$, $L = T = \Omega = 1$.
  Note that the solution is just well behaved when the interface is infinitely
  thin.
- *Figure 6 (interface-energy study):* distributions of the mole fraction after
  0.2 s for different choices of $\Delta F_{int}$; $\kappa = 100$,
  $L = T = \Omega = 1$.
- *Figure 8 (2D Cahn-Hilliard):* mole fraction after 0.1 s, $\kappa = 1$,
  $L = T = \Omega = 1$; the simulation took approximately 48 minutes (single
  threaded).
- *Figure 11 (ternary):* site fraction of component $x_1$ after 0.1 s;
  $x_0^{eq}=0.001$, $k=100$, $A_0 \rightarrow \infty$, $A_1=2$,
  $A_2=A_{\phi}=\Omega=1$; the simulation took approximately 22 minutes (single
  threaded).
- *Initial condition (Jul 2024):* the grey lines indicate the regions where
  individual phases are stable. The lowest region $0 \leq x \leq 0.25$ corresponds
  to phase α while the region $0.75 \leq x \leq 1$ corresponds to phase β. The
  middle $0.25 \leq x \leq 0.75$ marks the interface region where both phases are
  stable at different fractions.

---

## 11. CRediT statement of the November 2024 submission

W. Flachberger: Conceptualization, Methodology, Formal Analysis, Writing – original
draft. T. Antretter, J. Svoboda, S. Leitner, M. Petersmann: Supervision, Project
administration, Funding acquisition, Writing – review & editing.

---

## 12. References cited only in the removed passages

- Acosta-Soba, D., Guillén-González, F., Rodríguez-Galván, J. R.: An upwind DG
  scheme preserving the maximum principle for the convective Cahn-Hilliard model.
  Numerical Algorithms 92(3), 1589–1619 (2023).
- Acosta-Soba, D.: Analysis and numerical simulation of tumor growth models. PhD
  thesis, Universidad de Cádiz (2024).
- Dreyer, W., Guhlke, C.: Sharp limit of the viscous Cahn-Hilliard equation and
  thermodynamic consistency. Continuum Mech. Thermodyn. 29, 913–934 (2017).
- Flachberger, W.: Investigating the possibility to solve the
  Hamilton-Jacobi-Bellman equation by the finite element method to enable feedback
  control of nonlinear dynamical systems (2021).
- Kim, J. et al.: Basic principles and practical applications of the Cahn-Hilliard
  equation. Mathematical Problems in Engineering 2016, 9532608 (2016).
- Kondepudi, D., Prigogine, I.: Modern Thermodynamics: From Heat Engines to
  Dissipative Structures. Wiley (2014).
- Manning, J. R.: Correlation factors for diffusion in nondilute alloys. Phys. Rev.
  B 4(4), 1111 (1971).
- Miehe, C., Mauthe, S., Ulmer, H.: Formulation and numerical exploitation of mixed
  variational principles for coupled problems of Cahn-Hilliard-type and standard
  diffusion in elastic solids. Int. J. Numer. Meth. Engng 99(10), 737–762 (2014).
- Nagaraja, S. G., Flachberger, W., Antretter, T.: Variational formulation of
  Cahn-Hilliard-type diffusion coupled with crystal plasticity. Int. J. Plasticity
  167, 103652 (2023).
- Silverstein, T. P.: The real reason why oil and water don't mix. J. Chem. Educ.
  75(1), 116 (1998).
- Walker, S. W.: The Shapes of Things: A Practical Guide to Differential Geometry
  and the Shape Derivative. SIAM (2015).
