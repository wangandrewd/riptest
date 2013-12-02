import numpy as np
import random
import Image
#import pywt
from scipy.linalg import hadamard
from matplotlib import pyplot as plt
import math

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

def test_matrix(mat, k, reps):
    m,n = mat.shape
    #print k
    num_unif_trials = reps
    #num_unif_trials = 100*k

    # Generate random vectors of length n that are k sparse uniform random values
    errors = []
    for i in xrange(num_unif_trials):
        vec = generate_unif_vec(k, n)
        error =  abs(1 - np.linalg.norm(mat * vec, ord='fro')**2 / np.linalg.norm(vec, ord='fro')**2)
        errors.append(error)
    errors.sort()
    #print errors
    return max(errors)

def gradient_descent(mat, k):
    MAX_TRIAL = 10000
    thresh = .000005
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
    while (go_pos and curr_error - prev_error > thresh) or (not go_pos and prev_error - curr_error > thresh) and MAX_TRIAL > 0:
        gradient = np.matrix(np.zeros(n)).T
        #find the gradient
        for index in xrange(k):
            gprimex = 0
            for i in xrange(m):
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
        curr_error = np.linalg.norm(mat * vec )/np.linalg.norm(vec)
        print curr_error, prev_error

    return (prev_vec, curr_error)

def run_matrix_tests(mat, k, reps):
    error = test_matrix(mat, k, reps)
    return error

    #errors = []
    #errors_b = []
    #for i in xrange(trials):
    #    _, error = gradient_descent(mat, k)
    #    errors.append(error)
     
    #print error
    #errors.sort()
    #errors_b.sort()
    #print errors
    #print errors_b
    #print max(errors_b)

def run_test_suite(n, N, k, trials, reps):
    results = []
    for i in xrange(trials):
        print i
        mat = random_bernoulli(n, N)
        result = run_matrix_tests(mat, k, reps)
        results.append(result)
    results.sort()
    #print results
    return results

def find_n(k, N, trials, reps, epsilon, min_good, max_good):
    guesses = []
    guess_n = int(1 / (epsilon ** 2) * k * math.log(N / (1.0 * k)))
    increasing = 0
    round = 1
    lower_bound = -1
    upper_bound = -1
    while True:
        if lower_bound > 0 and guess_n <= lower_bound:
            round *= 2
            guess_n = int(1.5 ** (1.0 / round) * guess_n)
            continue
        if upper_bound > 0 and guess_n >= upper_bound:
            round *= 2
            guess_n = int(guess_n / (1.5 ** (1.0 / round)))
            continue

        guesses.append(guess_n)
        results = run_test_suite(guess_n, N, k, trials, reps)
        print guess_n
        print results
        if results[min_good - 1] > epsilon:
            if increasing == -1:
                round *= 2
            lower_bound = guess_n
            guess_n = int(1.5 ** (1.0 / round) * guess_n)
            increasing = 1
        elif results[max_good - 1] < epsilon:
            if increasing == 1:
                round *= 2
            upper_bound = guess_n
            guess_n = int(guess_n / (1.5 ** (1.0 / round)))
            increasing = -1
        else:
            break


def haar_decomp(img_file_name):
    img = Image.open(img_file_name).convert('L')
    x,y = img.size
    coeffs = pywt.wavedec2(img, 'haar', level=pywt.dwt_max_level(max(x,y), pywt.Wavelet('haar')))
    params = coeffs[0]
    final_arr = []
    for i in xrange(1,len(coeffs)):
	for j in coeffs[i]:
		final_arr = np.concatenate((final_arr, np.ndarray.flatten(j)))
    return params, final_arr

def harr_recomp(params, final_arr):
    return 0
    

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
    find_n(2, 100, 10, 20000, .5, 5, 8)
    #run_test_suite(25, 100, 2, 10, 50000)
    #haar_decomp("natural.jpg")

