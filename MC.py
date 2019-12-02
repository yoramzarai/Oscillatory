
import numpy as np
from itertools import combinations

def matrix_minor(A, r_indxs, c_indxs):
    '''This function returns the minor of a matrix A with rows indexed by r_indxs and 
    columns by c_indxs. r_indxs and c_indxs are lists (or 1D numpy arrays), and these indexes 
    start from 0 (which is the first row/column index).'''
    return np.linalg.det(A[np.ix_(r_indxs, c_indxs)]) if len(r_indxs)==len(c_indxs) else None

def compute_MC_matrix( A, p ):
    '''This function computes the p'th order multiplicative
    compound matrix of the given matrix A. It returns the MC
    matrix and the lexicography order (with 0 as the first index)'''
    x = np.arange(np.minimum(*A.shape), dtype=int)  # 0, .., n-1, where n:=min(matrix dimensions)
    lp = np.array(list(combinations(x, p))) # lexicography order of the p inxedes in x 
    lp_len = len(lp)
    Q = np.array([matrix_minor(A, lp[r], lp[c]) for r in range(lp_len) for c in range(lp_len)]).reshape(lp_len, lp_len)
    return Q, lp
