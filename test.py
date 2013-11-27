import numpy as np
import random

#Consts
SCALE = 2000



def error_info(errors):
    print "Average Error", np.average(errors)
    print "Max Error", max(errors)
    print "Min Error", min(errors)
    print "Error Std", np.std(errors)

def test_matrix(mat, k):
    m,n = mat.shape
    num_unif_trials = 10
    
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
    test_matrix(np.eye(10), 2)