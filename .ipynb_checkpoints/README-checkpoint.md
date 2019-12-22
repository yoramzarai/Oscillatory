# TP, TN and Oscillatory Matrices
This repository contains Python and Matlab functions and examples related to TP, TN, and oscillatory (OSC) matrices, the EB factorization of I-TN (i.e. invertible-TN) matrices, and the associated planar networks and triangle diagrams.


## Functions
- TP_TN_OSC_funcs.py contains functions related to EB factorization, planar networks and TP/TN/OSC matrices.

## Notebooks
- eb_factorization_example.ipynb contains several exmples of EB factorization and asscoiated planar networks.

- SEB_paths.ipynb uses the triangle graphs (derived from the planar networks of an SEB factorization) to deduce the exponent of an oscillatory matrix $A$, and to compute the families of vertex-disjoint paths of each corner minor of $A^r$, where $r$ is the exponent of $A$.

- analyze_osc.ipynb analyzes the exponent of an oscillatory matrices.

- MC.ipynb contains functions and example of computing the multiplicative compound of a matrix. Note that there functions there are included in TP_TN_OSC_funcs.py

## Matlab

- mat_call_py_example.m is a Matlab example using the Python functions in the files above.

- compute_sign_variations.m contains Matlab function for computing the number of sign variations in a vector.

