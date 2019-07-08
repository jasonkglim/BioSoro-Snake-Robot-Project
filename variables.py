import numpy as np
import math

# robot parameters
n = 3

# link parameters
link_mass = 1
link_length = 1
lG = 0.3 # distance to center of gravity of link
masses = np.repeat(link_mass, n)
lengths = np.repeat(link_length, n)
Glengths = np.repeat(lG, n)

L = np.sum(lengths)
mass_total = np.sum(masses)


  
# initialize all matrices and vectors
M = np.asmatrix(np.ones(shape=(n, n)))

M_0 = np.asmatrix(np.zeros(shape=(n, 2)))

m_0 =  mass_total * np.asmatrix(np.diag(np.ones(2)))

m = np.asmatrix(np.zeros(shape=(2, n)))

ff_vector = np.asmatrix(np.ones(shape=(2, 1)))

f0_vector = ff_vector

D = np.asmatrix(np.zeros(shape=(n, n-1)))
for i in range(n):
  for j in range(n-1):
    if i == j:
      D[i, j] = -1
    elif i == j+1:
      D[i, j] = 1
      
e = np.asmatrix(np.ones(shape=(n, 1)))


