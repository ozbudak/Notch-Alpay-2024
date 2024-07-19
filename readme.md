# Segmentation Clock Model Simulation

## Overview

This repository contains a set of MATLAB scripts for simulating the segmentation clock model using delay differential equations and is designed to perform optimization tasks using genetic algorithms. These scripts are intended to run on MATLAB R2024a and leverage the capabilities of the Optimization Toolbox and Parallel Processing Toolbox.

## Files

- `checkPeriod.m` - Calculates the period of oscillations in a given solution.
- `checkSusOsc.m` - Checks for sustained oscillations in a given solution.
- `dCmodel.m` - Simulates the segmentation clock model with given parameters and perturbations for two cells.
- `Figure_singleset.m` - Generates figures for single data sets.
- `FigureGenerator.m` - Generates average value figures for all the parameter sets.
- `getscoredC.m` - Evaluates the score of a parameter set based on the presence of sustained oscillations, period, and other criteria.
- `optimizedC.m` - Performs optimization using a genetic algorithm to find the best parameter set.

## Prerequisites

1. **MATLAB R2024a**: Ensure you have MATLAB R2024a installed.
2. **Optimization Toolbox**: Required for running genetic algorithms.
3. **Parallel Processing Toolbox**: Needed for executing computations across multiple processors or on a computer cluster.

## Usage

1. **Simulation**:
   To simulate the model with a given set of parameters and perturbations, use the `dCmodel` function:
   ```matlab
   [t, mh1Matrix, mh7Matrix] = dCmodel(parameters, pb1, pb2);
   ```

2. **Scoring**:
   To evaluate the score of a parameter set, use the `getscoredC` function:
   ```matlab
   score = getscoredC(parameters, pb1, pb2);
   ```

3. **Optimization**:
   To find the optimal parameters, run the `optimizedC.m` script:
   ```matlab
   run optimizedC.m
   ```

## Notes

- Ensure that all the files are in the same directory.
- The parameter ranges and perturbations can be adjusted within the scripts as needed.

