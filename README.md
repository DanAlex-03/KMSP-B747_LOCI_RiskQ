# LOC-I Risk Modeling using Subset Simulation and Energy-Based Metrics

This repository contains MATLAB code and supporting datasets developed for a thesis project that quantifies the probability of **Loss of Control In-Flight (LOC-I)** during the final approach phase of a **Boeing 747**, using a combination of physical energy modeling, flight data analysis, and rare-event simulation via **Subset Simulation**.

## Project Overview

The model focuses on the last **1000 ft AGL** before touchdown at **KMSP Runway 30R**, using Quick Access Recorder (QAR) data from over 4200 flights. It estimates the failure probability that the aircraft’s total energy drops below a defined threshold, indicating potential LOC-I.

Key components:
- Energy-rate modeling based on aerodynamic, thrust, and gravity force components.
- Data-driven statistical fitting of parameters across all flights.
- Definition of failure based on cumulative energy deficit (below threshold).
- Use of **Subset Simulation with MCMC sampling** to estimate rare-event probabilities.

The simulation includes both:
- An **overall failure probability** (integrated over entire approach),
- A **time-dependent per-second failure probability**, revealing how risk evolves during descent.

## Requirements

- MATLAB R2021b or newer
- [Subset Simulation Toolbox (TUM FSD)](https://www.fs.tum.de/en/aircraft/projects/software-tools/)
- Curve Fitting Toolbox
- Distribution Fitter App
- Parallel Computing Toolbox

## Key Features

- Energy model integrating control inputs: AoA, N1, pitch, bank, elevator, fuel mass, airspeed, altitude, and crosswind
- Dynamic per-second failure probability estimation
- Ranking of input parameters using both manual and toolbox-assisted sensitivity analysis
- Export-ready plots and Excel outputs for visualization

## Example Results

- Estimated LOC-I probability over full approach: **~3.57 × 10⁻⁸**
- The time-dependent failure probability was fit using a power-law model:
  P_F(t) ≈ 10^(a·t^b + c), with a goodness of fit R² = 0.9913.
- Sensitivity analysis reveals **angle of attack**, **crosswind**, and **N1** as most influential

## License and Permissions

This repository is **not licensed for reuse or redistribution**. Please contact the author for permission before using any part of the code or data.

## Acknowledgments

Developed as part of an undergraduate thesis (FTMD ITB, Student ID: 13621021) titled *Risk Quantification of Aerodynamic Stall and In-Flight Loss of Control (LOC-I) Events During Aircraft Approach*. Simulation methodology based on the TUM FSD Subset Simulation Toolbox and supervisor Dr. -Ing. Ir. Javensius Sembiring.
