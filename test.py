import numpy as np
import random
from scipy.linalg import hadamard
from matplotlib import pyplot as plt

#Consts
SCALE = 2000

def random_bernoulli(n,N):
    return (2*np.random.randint(2, size=(n,N))-1)/np.sqrt(n)
 
def random_gaussian(n,N):
    gauss_mat = np.random.normal(size=(n,N))
    for i in xrange(N):
        gauss_mat[:,i] = gauss_mat[:,i]/np.linalg.norm(gauss_mat[:,i])
    return gauss_mat

def fjlt_derive(n,N):
    S = np.zeros(n,N)
    for i in xrange(n):
        S[i,np.random.randint(0,N)] = 1
    S = S * np.sqrt(N)/np.sqrt(n)
    T = random_bernoulli(n,n)
    
    H = hadamard(N)
    return np.matrix(T) * np.matrix(S) * np.matrix(H)

def convert_codes_to_vector(codes,q,t):
    vector = np.zeros( (q*t,1) )
    for i in xrange(len(codes)):
        vector[q*i+codes[i]] = 1
    return vector
    
def error_code_mat(q,t,N):
    #generate random codes
    mat = np.random.randomint(q, size=(t,N))
    new_mat = np.zeros( (q*t,N) )
    for i in xrange(N):
        new_mat[:, i] = convert_codes_to_vector(mat[:,i])
    new_mat = new_mat / np.sqrt(t)
    return new_mat
    
 
def coherence(matrix):
    row, col = matrix.shape
    return max([ np.dot(matrix[:,i],matrix[:,j]) for i in xrange(col) for j in xrange(i+1, col)  ])

def error_info(errors):
    print "Average Error", np.average(errors)
    print "Max Error", max(errors)
    print "Min Error", min(errors)
    print "Error Std", np.std(errors)
    plt.hist(errors, bins = len(set(errors)))
    plt.show()
    
def test_matrix(mat, k):
    m,n = mat.shape
    num_unif_trials = 10
    print "Coherence", coherence(mat)
    
    # Generate random vectors of length n that are k sparse uniform random values
    errors = []
    for i in xrange(num_unif_trials):
        rand = SCALE * np.random.random(k) - (SCALE/2)
        locs = random.sample(xrange(n), k)
        vec = np.matrix(np.zeros(n)).T
        for i in xrange(k):
            vec[locs[i]][0] = rand[i]
        # Test
        error =  abs(1 - np.linalg.norm(mat * vec, ord='fro')**2 / np.linalg.norm(vec, ord='fro')**2)
        errors.append(error)
    error_info(errors)

if __name__ == "__main__":
    print coherence(np.eye(10))