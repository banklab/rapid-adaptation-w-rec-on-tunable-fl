import numpy as np
from numpy.linalg import inv

import sys

if len(sys.argv) != 2:
    print("Error: incorrect arguments. Command should be 'python3 SSWM.py filename'")

filename = sys.argv[1]
data = np.loadtxt(filename)

L = len(data[0]) - 1
n = 2**L

if L > 12:
    print("Error: landscape too large to analyze with this program.")

seqs = data[:, :L].astype(int)
fitness = data[:, L]

# list all neighbors in the sequence space
neighbors = np.zeros((n, L), dtype=int)
for i in range(n):
    idx = 0
    for j in range(n):
        if np.abs(seqs[i] - seqs[j]).sum() == 1:
            neighbors[i, idx] = j
            idx += 1

# Find the maxima
M = [[] for _ in seqs]
maxima, non_maxima = [], []
for i in range(n):
    f = fitness[i]
    for j in neighbors[i]:
        if f < fitness[j]:
            M[i].append(j)

    M[i] = np.array(M[i])
    if M[i].size > 0:
        non_maxima.append(i)
    else:
        maxima.append(i)

maxima = np.array(maxima)
non_maxima = np.array(non_maxima)

# calculate the transition matrix P
P = np.zeros((n, n))
for i in range(n):
    f = fitness[i]
    if i in maxima:
        P[i, i] = 1.
    else:
        for j in M[i]:
            P[i, j] = 1. - np.exp(-2.*(fitness[j] - f))

for i in non_maxima:
    P[i] /= P[i].sum()

# Obtain the submatrices
Q = np.delete(np.delete(P, maxima, axis=0), maxima,     axis=1)
R = np.delete(np.delete(P, maxima, axis=0), non_maxima, axis=1)
# Ik = np.delete(np.delete(P, non_maxima, axis=0), non_maxima, axis=1)
# O = np.delete(np.delete(P, non_maxima, axis=0), maxima, axis=1)

k = len(maxima)

print('maxima:')
for i in range(k):
    print(seqs[maxima[i]])
print()


I = np.identity(n - k)
N = inv(I - Q)

expected_t = np.dot(N, np.ones(n - k))
variance_t = np.dot(2*N - I, expected_t) - expected_t*expected_t

print('expected time to fixation:')
for i in range(n - k):
    print(seqs[non_maxima[i]], expected_t[i], variance_t[i])
print()

endstate = np.dot(N, R)
print('endstate probability:')
for j in range(k):
    for i in range(n - k):
        print(seqs[non_maxima[i]], ' -> ', seqs[maxima[j]], endstate[i, j])
    print()
