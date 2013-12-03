import numpy as np
import random
import Image
import pywt
from scipy.linalg import hadamard
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import math

# k = O(delta n/log(2N/n) )
# 
def random_bernoulli(n,N):
    return (2*np.random.randint(2, size=(n,N))-1)/np.sqrt(n)
 
def random_gaussian(n,N):
    gauss_mat = np.random.normal(loc=0,scale=1.0/np.sqrt(n),size=(n,N))
    return gauss_mat

def fjlt_derive(n,N):
    S = np.zeros((n,N))
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
    #print curr_error
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
        #print curr_error, prev_error

    return (prev_vec, curr_error)

def run_matrix_tests(mat, k, reps):
    error = test_matrix(mat, k, reps)
    #return error

    #errors = []
    #errors_b = []
    #for i in xrange(reps/100):
    #    _, error = gradient_descent(mat, k)
    #    errors_b.append(error)
     
    #print error
    #errors.sort()
    #errors_b.sort()
    #print errors
    #print errors_b
    #print max(errors_b)
    return error
    #return max(error, max(errors_b))

def run_test_suite(n, N, k, trials, reps, mat_gen):
    results = []
    for i in xrange(trials):
        #print i
        mat = mat_gen(n, N)
        result = run_matrix_tests(mat, k, reps)
        results.append(result)
    results.sort()
    #print results
    return results

def find_n(k, N, trials, reps, epsilon, min_good, max_good, mat_gen):
    guesses = []
    guess_n = int(1 / (epsilon ** 2) * k * math.log(N / (1.0 * k)))
    increasing = 0
    round_fac = 1
    lower_bound = -1
    upper_bound = -1
    int_min_good = int(min_good * trials) - 1
    int_max_good = int(max_good * trials) - 1

    while True:
        #print round_fac
        #if lower_bound > 0 and guess_n <= lower_bound:
        #    round_fac *= 2
        #    guess_n = int(1.5 ** (1.0 / round_fac) * guess_n)
        #    continue
        #if upper_bound > 0 and guess_n >= upper_bound:
        #    round_fac *= 2
        #    guess_n = int(guess_n / (1.5 ** (1.0 / round_fac)))
        #    continue

        guesses.append(guess_n)
        results = run_test_suite(guess_n, N, k, trials, reps, mat_gen)
        print guess_n
        print results
        if results[int_min_good] > epsilon:
            if increasing == -1:
                round_fac *= 2
            lower_bound = guess_n
            guess_n = int(1.5 ** (1.0 / round_fac) * guess_n)
            increasing = 1
        elif results[int_max_good] < epsilon:
            if increasing == 1:
                round_fac *= 2
            upper_bound = guess_n
            guess_n = int(guess_n / (1.5 ** (1.0 / round_fac)))
            increasing = -1
        else:
            return guess_n


def haar_decomp(img_file_name):
    img = Image.open(img_file_name).convert('L')
    x,y = img.size
    coeffs = pywt.wavedec2(img, 'haar', level=pywt.dwt_max_level(max(x,y), pywt.Wavelet('haar')))
    params = coeffs[0]
    final_arr = []
    for i in xrange(1,len(coeffs)):
	for j in coeffs[i]:
		final_arr = np.concatenate((final_arr, np.ndarray.flatten(j)))
    return params, coeffs[1:], final_arr

def reverse_flatten(arr_format, values, index):
    if not hasattr(arr_format, '__iter__'):
        return index

    for i in xrange(len(arr_format)):
        if hasattr(arr_format[i], '__iter__'):
            index = reverse_flatten(arr_format[i], values, index)
        else:
            arr_format[i] = values[index]
            index += 1

    return index

def haar_recomp(params, rest_coeff, final_arr):
    reverse_flatten(rest_coeff, final_arr, 0)
    final_tuple = [ params]
    for i in rest_coeff:
        final_tuple.append(i)
    final_tuple = tuple(final_tuple)
    return pywt.waverec2(final_tuple, 'haar')


#format of result_dict
#result_dict -> matrix type -> N,k,epsilon -> n
def write_csv(file_name, result_dict):
    f = open(file_name, 'w')
    print_str = "%d, %d, %f,"
    for matrix_type in result_dict.keys():
        f.write(matrix_type)
        f.write("\n")
        for nkepsilon in result_dict[matrix_type].keys():
            print nkepsilon
            f.write( print_str % nkepsilon )
            f.write( str(result_dict[matrix_type][nkepsilon]))
            f.write("\n")
    f.close()

if __name__ == "__main__":
    result_dict = dict(bern=dict(), gauss=dict())
    Ns = [1000]#, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000]
    epsilons = [0.1] #, 0.25, 0.5]
    for N in Ns:
        for epsilon in epsilons:
            values_of_k = 5
            value_of_k = 2
            while 1.0/ epsilon ** 2 * value_of_k * math.log(N * 1.0 / value_of_k) < N:
                value_of_k  = int(value_of_k * 1.5)
            value_of_k = int(value_of_k / 1.5)
            k_set = [int((value_of_k / 2) ** (i / (1.0 * values_of_k)) * 2) for i in range(values_of_k)]
            k_set = list(set(k_set))
            k_set.sort()
            print k_set
            for k in k_set:
                k = int(k)
                print "K"
                print k
                print "Bernoulli"

                result_dict['bern'][(N,k,epsilon)] =   find_n(k, 100, 10, 20000, .5, .5, .8, random_bernoulli)
                print "Gaussian"
                result_dict['gauss'][(N,k,epsilon)] = find_n(k, 100, 10, 20000, .5, .5, .8, random_gaussian)
                #print "FJLT"
                #find_n(k, 100, 10, 20000, .5, .5, .8, fjlt_derive)
                #find_n(k, 100, 10, 20000, .5, 5, 8)
    write_csv("THUNDERBEAR.txt", result_dict)

