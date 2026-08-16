# Shallow Water Equations using Roe's Riemann Solver

> **Course Project : AS5420 — Intro to CFD**

This project implements an approximate Riemann-solver-based finite-volume method for solving the **Shallow Water Equations (SWEs)** in one and two dimensions. The implementation is based on Roe's approximate Riemann solver and is developed to study the numerical treatment of discontinuous free-surface flows, including dam-break problems and flow over uneven bottom topography.

The project reproduces and extends the methodology presented in *"Approximation of Shallow Water Equations by Roe's Riemann Solver"* by Ambrosi (1995), with numerical experiments covering both 1D and 2D configurations.

---

## Overview

The inviscid shallow water equations form a hyperbolic system used to model a wide range of free-surface flows, including open-channel flows, flood propagation, and dam-break problems.

The main objective of this project is to implement a numerical solver capable of capturing rapidly varying and discontinuous solutions while maintaining numerical stability and avoiding non-physical oscillations.

The implementation includes:

- Finite-volume discretization of the shallow water equations
- Roe's approximate Riemann solver
- Roe-averaged velocity and wave celerity
- Treatment of geometric source terms arising from uneven bottom topography
- CFL-based adaptive timestep restriction
- Wet and dry dam-break configurations
- Extension from one to two spatial dimensions using dimensional splitting
- Visualization and export of numerical solutions

---

## Numerical Method

The solver follows a finite-volume formulation in which the solution is evolved using numerical fluxes across cell interfaces.

The main components of the numerical method are:

### Roe Approximate Riemann Solver

The Riemann problem at each interface is approximated using Roe-averaged flow variables. This provides a characteristic-based flux formulation suitable for capturing discontinuities such as hydraulic bores and dam-break fronts.

### Flux and Wave Propagation

The eigenvalues of the Roe-averaged system determine the propagation speeds of the characteristic waves. The numerical flux is constructed from the corresponding wave strengths and eigenvectors.

### Source Terms

For non-flat bottom topography, geometric source terms are included in the shallow water equations. A still-water test over an uneven bottom is used to investigate the resulting numerical balance and truncation error.

### CFL Stability Condition

The timestep is restricted according to the CFL condition,

$$
\Delta t
\leq
C
\frac{\Delta x}
{\max(|v|+c)},
$$

where $C$ is the Courant number and

$$
c=\sqrt{gh}
$$

is the local wave celerity.

The simulations use a Courant number close to 0.9.

---

## Test Cases

The solver is tested on four representative configurations.

### 1. Still Water over an Uneven Bottom

A still-water configuration over a sloping bottom is used to examine the treatment of geometric source terms.

The analytical solution corresponds to a stationary free surface, making this a useful test of the numerical balance between fluxes and source terms.

---

### 2. 1D Dam Break on a Dry Bottom

A dam-break problem with an initially dry downstream region is used to examine the behavior of the solver near a wet-dry interface.

The simulation tests the ability of the numerical scheme to maintain monotonicity and avoid non-physical negative water depths as the advancing front approaches the dry region.

---

### 3. 1D Dam Break in a Wet Basin

A second dam-break configuration is considered with a non-zero downstream water depth.

This avoids the dry-bed singularity and provides a standard test for shock propagation and numerical oscillation control.

---

### 4. 2D Partial Dam Break

The 1D solver is extended to two dimensions using **dimensional splitting**, in which the numerical solution is advanced successively in the $x$ and $y$ directions.

A partial breach in a rectangular reservoir is simulated to demonstrate two-dimensional shock propagation and wave reflection from the reservoir boundaries.

---

## Results

The numerical experiments demonstrate the ability of the solver to:

- Capture dam-break fronts and hydraulic bores
- Maintain monotonicity near strong discontinuities
- Handle wet-dry interfaces with a small depth threshold
- Simulate wave propagation in two spatial dimensions
- Represent the effects of uneven bottom topography
- Maintain numerical stability under a CFL restriction

For the 1D wet-bottom dam-break case, the computed bore is captured sharply within a small number of computational cells, while the 2D simulation demonstrates the propagation and reflection of the resulting wave field.

│
├── requirements.txt
├── .gitignore
└── README.md
