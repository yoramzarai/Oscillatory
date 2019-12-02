#!/opt/local/bin/python3
'''The method here to determine the entries in A^{(p)} that positive follows the proof of Thm. 5.2 in Pinkus book
"Totally Positive Matrices".
''' 
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt
import time
import argparse

# Command-line parsing information
def parse_in_args():
    ''' Defines input arguments '''
    # place description
    parser = argparse.ArgumentParser(description='Computes the positive entries of a multiplicative compound (MC) of an oscillatory matrix.', add_help=False)

    # required arguments
    group = parser.add_argument_group('Required arguments')
    group.add_argument('-n', dest='n', help='Square matrix dimension', type=int, metavar='<n>', required='True')
    group.add_argument('-p', dest='p', help='MC order', type=int, metavar='<p>', required='True')


    # optional arguments
    group_opt = parser.add_argument_group('Optional arguments')
    group_opt.add_argument('-v', help='Verbose (display evolution of (A^(p))^i).', action='store_true')
    #group_opt.add_argument('-D', help='Enables debug prints.', action='store_true')
    group_opt.add_argument('-h', '--help', help='Displays usage information and exits.', action='help')

    return parser.parse_args(), parser.print_help

def is_in_set(i,j):
    '''This function checks if the following two conditions hold
    for all k = 1, .., len(i):
    1. |i_k-j_k|<=1.
    2. max{i_k, j_k}<min{i_{k+1},j_{k+1}}.
    See more details in the proof of Thm. 5.2. in Pinkus book.'''
    
    # checking condition 1
    if not np.all(np.absolute(i-j)<=1):
        return False
    
    # checking condition 2
    mx = np.maximum(i, j)
    mn = np.minimum(i, j)
    if not np.all(mx[:-1]<mn[1:]):
        return False
        
    return True


def compute_MC_mask(n, p):
    '''This function computes the mask matrix of the multiplicative compound matrix
    of a TN, nonsingular and irreducible matrix A (size nxn), 
    i.e. it rerurns a matrix size (n over k)x(n over k), where the (i,j) entry is 1 [0] 
    if the corresponding minor is >0 [>=0].'''
    
    x = np.arange(1, n+1, dtype=int)
    lp = np.array(list(combinations(x, p))) # lexicography order of the p inxedes in x 
    lp_len = len(lp)
    y = np.zeros((lp_len, lp_len), dtype=int)

    for i in range(lp_len):
        for j in range(lp_len):
            y[i,j] = is_in_set(lp[i], lp[j])
    
    return y, lp

# =======================================================================================================================================

# Main function
# =============
def main():
    ''' Main body '''
    args, _ = parse_in_args()  # parse, validate and return input arguments

    n, p = (args.n, args.p)  # n = the square matrix size, p = MC order
    if p>n:
        print("n must be larger than p !!")
        return

    y, lp = compute_MC_mask(n, p)

    # check if y^(n-1)>>0
    yn1 = np.linalg.matrix_power(y, n-1)
    if np.all(yn1>0):
        print("(A^({}))^(n-1)>>0".format(p))
    else:
        print("Not all entries of (A^({}))^(n-1) are positive !!".format(p))

    print("Mask of A^({}) (A size {}x{}) (Blue means >0, white means >=0)".format(p,n,n))
    lp_len = len(lp)

    fig = plt.figure(figsize=(8,8)); ax = fig.add_subplot(1,1,1)
    h = ax.matshow(y, cmap=plt.cm.Blues)
    ax.set_xticks(np.arange(-.5, lp_len, 1))
    ax.set_yticks(np.arange(-.5, lp_len, 1))
    ax.set_xticklabels(np.arange(0, lp_len+1))
    ax.set_yticklabels(np.arange(0, lp_len+1))
    ax.set_title("A^({}), with A size {}x{} (blue: >0, white: >=0)".format(p, n, n))
    plt.grid()
    #plt.show(block=False)
    plt.show()
    
    # printing the indexes corresponding to the figure indexes
    print("(A^({}))^(n-1) row and column indexes:".format(p))
    for i, r in zip(range(1,lp_len+1), lp):
        print("\t{} =>{}".format(i, r))
        
    print("\n(A^({}))^(n-1)=\n{}".format(p, yn1))

    # A^(i) evolution
    if args.v:
        fig = plt.figure(figsize=(8,8)); ax = fig.add_subplot(1,1,1)
        for i in range(1,n):
            ax.cla()
            yi = np.linalg.matrix_power(y, i)
            #h = ax.matshow(yi>0, cmap=plt.cm.Blues)
            h = ax.imshow(yi>0, cmap=plt.cm.Blues)
            ax.set_xticks(np.arange(-.5, lp_len, 1))
            ax.set_yticks(np.arange(-.5, lp_len, 1))
            ax.set_xticklabels(np.arange(0, lp_len+1))
            ax.set_yticklabels(np.arange(0, lp_len+1))
            plt.grid()
            plt.title("(A^({}))^{}, with A size {}x{}".format(p,i, n, n))
            plt.pause(2)

    
# ================================================================================================
if __name__ == "__main__":
    main()
else:
    print(__name__, 'has been imported.')
