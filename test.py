import numpy as np
import random
from scipy.linalg import hadamard
from matplotlib import pyplot as plt

#Consts
SCALE = 2000


# k = O(delta n/log(2N/n) )
# 
def random_bernoulli(n,N):
    return (2*np.random.randint(2, size=(n,N))-1)/np.sqrt(n)
 
def random_gaussian(n,N):
    gauss_mat = np.random.normal(loc=0,scale=1.0/np.sqrt(n),size=(n,N))
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
    #plt.hist(errors, bins = len(set(errors)))
    #plt.show()

def generate_unif_vec(k,N):
    rand = 2*np.random.random(k) -1
    locs = random.sample(xrange(N), k)
    vec = np.matrix(np.zeros(N)).T
    for i in xrange(k):
        vec[locs[i]][0] = rand[i]
    return vec

def test_matrix(mat, k):
    m,n = mat.shape
    print k
    num_unif_trials = 100*k

    # Generate random vectors of length n that are k sparse uniform random values
    errors = []
    for i in xrange(num_unif_trials):
        vec = generate_unif_vec(k, n)
        error =  abs(1 - np.linalg.norm(mat * vec, ord='fro')**2 / np.linalg.norm(vec, ord='fro')**2)
        errors.append(error)
    return max(errors)

def gradient_descent(mat, k):
    MAX_TRIAL = 10000
    m,n = mat.shape
    alpha = 0.3
    # pick k random coordinates and random numbs
    rand = 2*np.random.random(k) -1
    locs = random.sample(xrange(n), k)
    vec = np.matrix(np.zeros(n)).T
    for i in xrange(k):
        vec[locs[i]][0] = rand[i]
    prev_vec = np.copy(vec)
    curr_error = np.linalg.norm(mat * vec )/np.linalg.norm(vec)
    prev_error = 1
    go_pos = curr_error > 1
    print curr_error
    while (go_pos and curr_error - prev_error > 0) or (not go_pos and prev_error - curr_error > 0) and MAX_TRIAL > 0:
        gradient = np.matrix(np.zeros(n)).T
        #find the gradient
        for index in xrange(k):
            deriv = 0
            for i in xrange(m):
                for j in xrange(n):
                    deriv += 2*mat[i][locs[index]]*mat[i][j]*vec[j]
            gradient[locs[index]][0] = deriv - 2 * vec[locs[index]]
        #move in direction of gradient
        prev_vec = np.copy(vec)
        if go_pos:
             vec = vec + alpha * gradient
        else:
            vec = vec - alpha*gradient
        MAX_TRIAL -= 1
        prev_error = curr_error
        curr_error = np.linalg.norm(mat * vec )/np.linalg.norm(vec)
        print curr_error, prev_error

    return prev_vec

if __name__ == "__main__":
    """
    N = 1000
    ks = np.linspace(N/10, N/2, 10)
    unif_error = []
    for k in ks:
        for n in xrange(int(k), N, int((N-k)/10)):
            print (k,n,N)
            #unif_error.append( (k,n,N,test_matrix(random_bernoulli(n, N), int(k)) ) )
    print unif_error
    """
    gradient_descent(random_bernoulli(25,100), 2)