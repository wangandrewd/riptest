import math
import numpy as np

m = 7
epsilon = 1.0/(402*m)

def generate_primes_until(n):
    reverse_prime_dict = dict()
    for i in xrange(2, n):
        if i not in reverse_prime_dict:
            j = 2
            while i*j in reverse_prime_dict:
                j += 1
            reverse_prime_dict[i*j] = i
        else:
            prime = reverse_prime_dict[i]
            factor = i/prime
            reverse_prime_dict.pop(i, None)
            j = 1
            while (factor+j)*prime in reverse_prime_dict:
                j += 1
            reverse_prime_dict[(factor+j)*prime] = prime
    return reverse_prime_dict

def get_prime(k, N, epsilon):
    assert(k**(2-epsilon) <= N)
    assert(k**(2+epsilon) >= N)
    # find a prime in the interval [k**(2-epsilon), 2k**(2-epsilon) ]
    reverse_prime_dict = generate_primes_until(int(math.ceil(math.sqrt(2*k**(2-epsilon)))))
    for i in xrange(int(math.floor(k**(2-epsilon))), int(math.ceil(2*k**(2-epsilon)))):
        if i not in reverse_prime_dict:
            return i
    raise Exception("bug code")

def generate_fancy_A(p, alpha):
    return xrange(1,int(math.floor(p**alpha))+1)

def perm_gen(k, arr):
    if k==0:
        yield []
    else:
        for subarr in perm_gen(k-1, arr):
            for i in arr:
                yield [i] + subarr 
    
def generate_fancy_B(r, M):
    tally = 0
    for num_set in perm_gen(int(r), range(int(M))):
        for i in xrange(int(r)):
            tally += num_set[i]*(2*M)**(i)
        yield tally
    
def bourgain(N, k=None):
    if k is None:
        k = math.sqrt(N)
    p = get_prime(k,N,epsilon)
    alpha = 1.0/(2*m)
    fancy_A = generate_fancy_A(p, alpha)
    beta = 1/(2.01 * m)
    r = math.floor(beta * math.log(p)/math.log(2))
    M = math.floor(2**(2.01*m)-1)
    fancy_B = generate_fancy_B(r, M)
    
    bourgain_mat = np.zeros((p,N), dtype=np.complex128)
    
    N_counter = 0
    for a in fancy_A:
        for b in fancy_B:
            
            if N_counter >= N:
                return bourgain_mat
            
            x_range = np.linspace(1,p,p)
            bourgain_mat[N_counter][0:p] = 1.0/math.sqrt(p) * np.exp( 2*np.pi/p*(a*x_range*x_range + b*x_range)*1j)
            N_counter += 1
    return bourgain_mat
    
print bourgain(2**14)