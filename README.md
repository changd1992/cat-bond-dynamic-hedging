# Cat Bond Dynamic Hedging

## Overview

This repository contains MATLAB code for simulating tropical cyclone (TC) trajectories, evaluating catastrophe-linked securities, and implementing dynamic hedging strategies.

The framework simulates stochastic TC evolution, calculates payoffs, and evaluates hedging performance.

Two hedging approaches are implemented:

- **Delta hedging:** A dynamic hedging strategy based on the sensitivity of the CAT bond's value to changes in underlying state variables.
- **Minimum Variance Ratio (MVR) hedging:** A hedging strategy designed to minimize the variance of the hedged position.

## MATLAB Files

| File | Description |
|---|---|
| `initialize_parameters.m` | Initializes model parameters used in the simulations. |
| `simulate_tc_path.m` | Simulates a single stochastic tropical cyclone path. |
| `simulate_tc_ensemble.m` | Simulates an ensemble of TC paths, calculates payoffs, and produces aggregated simulation results. |
| `simulate_tc_ensemble_hedge.m` | Extends the ensemble simulation to calculate delta and MVR hedges, along with payoffs and aggregated results. |
| `simulate.m` | Demonstrates the simulation and hedging framework using a simple example. |

## Usage

Run `simulate.m` in MATLAB to execute the example implementation.

The remaining functions can also be used individually for TC path simulation, ensemble analysis, and hedging calculations.

## Development Status

The code is intended for research purposes and may be updated as the simulation and hedging framework is further developed.
