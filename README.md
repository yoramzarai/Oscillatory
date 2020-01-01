# TP, TN and Oscillatory Matrices
This repository contains Python and Matlab functions, and examples related to TP, TN, and oscillatory (OSC) matrices, the SEB factorization of I-TN (i.e. invertible-TN) matrices, and the associated planar networks and triangle diagrams.


## Functions
- TP_TN_OSC_funcs.py contains several functions related to the number of sign variations in a vector, multiplicative compound matrices, SEB factorization, planar networks and other TP/TN/OSC matrices related functions.

## Notebooks
- eb_factorization_example.ipynb contains several exmples of SEB factorization and asscoiated planar networks.

- SEB_paths.ipynb uses the triangle graph (derived from the planar networks of an SEB factorization) to deduce the exponent of an oscillatory matrix $A$, and to compute the families of vertex-disjoint paths of each corner minor of $A^r$, where $r$ is the exponent of $A$.

- analyze_osc.ipynb analyzes the exponent of an oscillatory matrices.

- MC.ipynb contains functions and example of computing the multiplicative compound of a matrix.

## Matlab

- mat_call_py_example.m is a Matlab example using the Python functions in TP_TN_OSC_funcs.py.

- compute_sign_variations.m contains Matlab function for computing the number of sign variations in a vector.

