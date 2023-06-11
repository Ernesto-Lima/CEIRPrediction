# Predicting Response to Combination Evofosfamide and Immunotherapy

This repository contains the code and data files associated with the research paper titled "Predicting response to combination evofosfamide and immunotherapy under hypoxic conditions in murine models of colon cancer." The goal of the study is to develop a mathematical model that captures the interaction between evofosfamide, immunotherapy, and the hypoxic landscape of the tumor in the treatment of tumors.

## Code Structure

The repository is organized as follows:

- `calibration/`: Contains the C++ code implementing the mathematical model given by Eqs (1) and (2) from our paper. The model is solved using a fourth-order Runge-Kutta method. Includes Python files for post-processing the results. It contains scripts to analyze the data obtained from the model simulations.

- `sensitivity_analysis/`: Includes Python files for performing sensitivity analysis.

- `data/`: Contains the data files used for calibration and analysis. The data obtained from mice experiments with colon adenocarcinoma cells, as well as other relevant data, are included in this directory.

## Calibration Framework

The calibration framework employed in this study is implemented in C++ using the QUESO (Quantification of Uncertainty for Estimation, Simulation, and Optimization) library. The model parameters are estimated using a parallel, adaptive, multilevel Markov Chain Monte Carlo (MCMC) sampling method.

## Getting Started

To use the code and reproduce the results:

1. Clone this repository to your local machine.

2. Follow the instructions in the respective directories (`calibration/` and `sensitivity_analysis/`) to compile and run the code.

3. Refer to the README files within each directory for more detailed instructions on how to use the code, preprocess data, and interpret the results.

## Citation

If you use the code or findings from this work, please consider citing the original paper.

## License

GNU General Public License, Version 3, 29 June 2007

## Contact

For any questions or inquiries regarding the code or the research, please contact: ernesto.lima AT utexas.edu.
