import numpy as np


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

reed_solomon_code(3, 1)
