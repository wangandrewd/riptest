import numpy as np
import random
import Image
#import pywt
from scipy.linalg import hadamard
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import math
from matplotlib import rc

# k = O(delta n/log(2N/n) )
# 
def random_bernoulli(n,N):
    return (2*np.random.randint(2, size=(n,N))-1)/np.sqrt(n)
 
def random_gaussian(n,N):
    gauss_mat = np.random.normal(loc=0,scale=1.0/np.sqrt(n),size=(n,N))
    for i in xrange(N):
        gauss_mat[:,i] = gauss_mat[:,i]/np.linalg.norm(gauss_mat[:,i])
    return gauss_mat

def fjlt_derive(n,N):
    S = np.zeros((n,N))
    for i in xrange(n):
        S[i,np.random.randint(0,N)] = 1
    S = S * np.sqrt(N)/np.sqrt(n)
    T = random_bernoulli(n,n)
    
    H = hadamard(N)
    return np.matrix(T) * np.matrix(S) * np.matrix(H)

def evaluate(a, j, q):
    return sum((a[i] * (j ** i)) % q for i in range(len(a))) % q

def reed_solomon_code(q,d):
    mat = np.matrix(np.zeros((q * q, q ** (d + 1))))
    a = [0] * (d + 1)
    for i in range(q ** (d + 1)):
        for j in range(q):
            val = evaluate(a, j, q)
            mat[j * q + val,i] = 1 / np.sqrt(q)

        curr_ind = 0
        a[0] += 1
        if a[0] == q:
            a[0] = 0
        while a[curr_ind] == 0 and curr_ind < d:
            curr_ind += 1
            a[curr_ind] += 1
            if a[curr_ind] == q:
                a[curr_ind] = 0
    print mat

def error_correcting_code(k, N, epsilon):
    alpha = epsilon / k
    q_guess = int((math.log(n) / alpha) / (math.log(math.log(n)) + math.log(1 / alpha)))
    d_guess = int(alpha * q_guess)
    #print q_guess
    #print d_guess
    while q_guess ** (d_guess + 1) < N:
        if alpha * q_guess >= d_guess + 1:
            d_guess += 1
        else:
            q_guess += 1
    last = 0
    while q_guess ** (d_guess + 1) >= N:
        if alpha * (q_guess - 1) >= d_guess:
            q_guess -= 1
            last = 1
        else:
            d_guess -= 1
            last = 2
    if last == 1:
        q_guess += 1
    elif last == 2:
        d_guess += 1
    else:
        print "SHIIIIIET"
    #print q_guess
    #print d_guess
    #print q_guess ** (d_guess + 1)
    mat = reed_solomon_code(q_guess, d_guess)
    return mat
 
def coherence(matrix):
    row, col = matrix.shape
    return max([ np.dot(matrix[:,i],matrix[:,j]) for i in xrange(col) for j in xrange(i+1, col)  ])

def error_info(errors, reps, N, n, k, epsilon):
    print "Average Error", np.average(errors)
    print "Max Error", max(errors)
    print "Min Error", min(errors)
    print "Error Std", np.std(errors)
    str_format = "N=%d, n=%d, k=%d, epsilon=%f"
    fig, ax = plt.subplots()
    ax.set_yscale('log', basey=10)
    ax.hist(errors, bins = reps/10)
    ax.set_xlabel('Values')
    ax.set_ylabel('Frequency')
    ax.set_title(str_format % (N,n,k,epsilon) )
    plt.show()

def generate_unif_vec(k,N):
    rand = 2*np.random.random(k) -1
    locs = random.sample(xrange(N), k)
    vec = np.matrix(np.zeros(N)).T
    for i in xrange(k):
        vec[locs[i]][0] = rand[i]
    return vec, locs

def test_matrix(mat, k, reps):
    m,n = mat.shape
    num_unif_trials = reps

    # Generate random vectors of length n that are k sparse uniform random values
    errors = []
    for i in xrange(num_unif_trials):
        vec, locs = generate_unif_vec(k, n)
        error =  abs(1 - np.linalg.norm(mat * vec, ord='fro') / np.linalg.norm(vec, ord='fro'))
        errors.append(error)
    return errors

def gradient_descent(mat, k): #vec, locs):

    MAX_TRIAL = 1000
    thresh = .000005
    m,n = mat.shape
    vec, locs = generate_unif_vec(k, n)
    alpha = 0.3
    prev_vec = np.copy(vec)
    curr_error = np.linalg.norm(mat * vec )/np.linalg.norm(vec)-1
    prev_error = 1
    go_pos = curr_error > 0
    #print curr_error
    while (go_pos and curr_error - prev_error > thresh) or (not go_pos and prev_error - curr_error > thresh) and MAX_TRIAL > 0:
        gradient = np.matrix(np.zeros(n)).T
        #find the gradient
        for index in xrange(k):
            gprimex = 0
            for i in xrange(m):
                assert(locs[index] < n)
                gprimex += 2*mat[i][locs[index]]*np.dot(np.array(mat[i]),np.array(vec))
            hprimex = 2*vec[locs[index]]
            gx = np.linalg.norm(mat*vec)**2
            hx = np.linalg.norm(vec)**2
            gradient[locs[index]][0] = (gprimex*hx - hprimex*gx)/(hx*hx)
        #move in direction of gradient
        prev_vec = np.copy(vec)
        if go_pos:
             vec = vec + alpha * gradient
        else:
            vec = vec - alpha*gradient
        MAX_TRIAL -= 1
        prev_error = curr_error
        curr_error = np.linalg.norm(mat * vec )/np.linalg.norm(vec)-1
        #print curr_error, prev_error

    return (prev_vec, curr_error)

def run_matrix_tests(mat, k, reps):
    errors = test_matrix(mat, k, reps)
    #return error
    #for i in xrange(reps/100):
    #    errors.append(gradient_descent(mat, k))
    return errors
    #vec, grad_error = gradient_descent(mat, k, worst_error_vec, worst_locs)
    #print grad_error
    #return error
    #return max(error, abs(grad_error))

def run_test_suite(n, N, k, trials, reps, mat_gen):
    results = []
    for i in xrange(trials):
        mat = mat_gen(n, N)
        result = run_matrix_tests(mat, k, reps)
        results.append(result)
    results.sort()
    return results

if __name__ == "__main__":
    N = 8241
    n = 435
    k = 3264
    epsilon = np.sqrt(2)-1
    reps = 20000
    error_info(run_matrix_tests(random_gaussian(n,N), k, reps), reps, N, n, k, epsilon)
