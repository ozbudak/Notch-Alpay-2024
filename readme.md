# Segmentation Clock Model - DeltaC Positive Feedback

## Overview

This repository contains a set of MATLAB scripts for simulating the segmentation clock model using delay differential equations and is designed to perform optimization tasks using genetic algorithms. These scripts are intended to run on MATLAB R2024a and leverage the capabilities of the Optimization Toolbox and Parallel Processing Toolbox.

## Files

- `checkPeriod.m` - Calculates the period of oscillations in a given solution.
- `checkSusOsc.m` - Checks for sustained oscillations in a given solution.
- `syncBreak.m` - Determines the time at which synchrony between two cells breaks.
- `dCmodelPF.m` - Simulates the Positive Feedback model with given parameters and perturbations for two cells.
- `dCmodelNPF.m` - Simulates the No Positive Feedback model with given parameters and perturbations for two cells.
- `Figure_singleset.m` - Generates synchronization break time figures for single data sets on different mutants.
- `FigureGenerator.m` - Generates average synchronization break time and mRNA level figures for all parameter sets.
- `getscoredC.m` - Evaluates the score of a parameter set based on the presence of sustained oscillations, period, and other criteria in mutants for the Positive Feedback model.
- `getscoredC2cond.m` - Evaluates the score of a parameter set based on the presence of sustained oscillations and period in Wild Type embryos for the Positive Feedback model.
- `getscoredC3cond.m` - Evaluates the score of a parameter set based on the presence of sustained oscillations, period, and synchronization score in Wild Type embryos for the Positive Feedback model.
- `getscoredCNPF2cond.m` - Evaluates the score of a parameter set based on the presence of sustained oscillations and period in Wild Type embryos for the No Positive Feedback model.
- `getscoredCNPF3cond.m` - Evaluates the score of a parameter set based on the presence of sustained oscillations, period, and synchronization score in Wild Type embryos for the No Positive Feedback model.
- `optimizedC.m` - Performs optimization using a genetic algorithm to find the best parameter set.
- `optimizedC2conditions.m` - Performs parameter search to compare two models based on oscillation conditions and generates figures of synchronization scores for parameter sets.
- `optimizedC3conditions.m` - Performs parameter search to compare two models based on oscillation conditions and synchronization score threshold, and generates figures of synchronization scores for parameter sets.

## Prerequisites

1. **MATLAB R2024a**: Ensure you have MATLAB R2024a installed.
2. **Optimization Toolbox**: Required for running genetic algorithms.
3. **Parallel Processing Toolbox**: Needed for executing computations across multiple processors or on a computer cluster.

## Usage

1. **Simulation**:
   To simulate the model with a given set of parameters and perturbations, use the `dCmodelPF` function:
   ```matlab
   [t, mh1Matrix, mh7Matrix] = dCmodelPF(parameters, pb1, pb2);
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

4. **Generate Figures for the Parameter Sets**:
   Put all the sets `.mat` files in the same folder and customize the code in `Figure_singleset.m` and `FigureGenerator.m` based on the filenames. Then, run these scripts in MATLAB.

5. **Compare Models with and without Positive Feedback**:
   Run the `optimizedC2conditions.m` script:
   ```matlab
   run optimizedC2conditions.m
   ```

## Notes

- Ensure that all the files are in the same directory.
- The parameter ranges and perturbations can be adjusted within the scripts as needed.

