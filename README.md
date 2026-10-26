# Supplementary materials

This repository contains the verification code, numerical data, and output logs accompanying the paper [arXiv:2610.32600](https://arxiv.org/abs/2610.32600).

The materials support the computer-assisted proof of the existence of new time-periodic solutions of the conformal cubic wave equation on the Einstein cylinder.

## Contents

- `Proof.m`: the _Mathematica_ script used to verify the prescribed norm bounds and the inequalities establishing existence. The same script is used for all three solutions, with the appropriate data set and parameters selected before execution.
- `Distance.m`: a separate script that verifies separation of the three solutions.
- `Data/`: the data sets for the three approximate solutions.
- `Logs/`: archived verification output.

## Reproducing the verification

The verification was performed using _Mathematica_ 15.0.0. No external packages are required, and all mathematical verification computations use exact rational arithmetic.

Run the scripts from the repository root, following the configuration and execution instructions in the paper. Run `Proof.m` separately for each data set, preferably in a fresh kernel. The script `Distance.m` loads all three data sets.

## Licensing

- **Code:** `Proof.m` and `Distance.m` are licensed under the [MIT License](LICENSE).
- **Data and logs:** the contents of `Data/` and `Logs/` are licensed under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](LICENSE-DATA).

## Citation

If you use the code or data in your research, please cite the accompanying paper linked above.
