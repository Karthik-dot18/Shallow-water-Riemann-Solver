# Shallow Water Riemann Solver
Course Project : AS5420 - Introduction to CFD
A finite-volume numerical solver for the one- and two-dimensional shallow water equations using Roe's approximate Riemann solver.

## Features

- 1D shallow water equations
- 2D shallow water equations
- Roe approximate Riemann solver
- Reflective wall boundary conditions
- Dimensional splitting for 2D simulations
- Uneven-bottom still-water test
- 1D dam break on a dry bed
- 1D dam break in a wet basin
- 2D partial dam-break simulation
- Numerical results exported to CSV
- Visualization using Matplotlib

## Test Cases

### 1. Still Water over an Uneven Bottom

Tests the treatment of the geometric source term and preservation of a steady free surface.

### 2. 1D Dam Break on a Dry Bed

Tests the solver near a wet-dry interface and verifies monotonic propagation of the dam-break front.

### 3. 1D Dam Break in a Wet Basin

Tests shock propagation in a non-zero downstream water depth.

### 4. 2D Partial Dam Break

A two-dimensional reservoir simulation using dimensional splitting.

## Numerical Method

The solver is based on a finite-volume formulation of the shallow water equations.

The numerical flux is evaluated using Roe-averaged quantities and an approximate Riemann solver. The 2D solver uses dimensional splitting to apply the 1D solver successively in the coordinate directions.
