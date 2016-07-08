import random
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse.linalg as spla
from visual import*
'''
    Simple evolvement for a wave function.
    phi =E ** (-0.01 * x ** 2)
    dphi_dt = 0 * x
    d2phi_dx2 = (E ** (-0.01 * x ** 2)) * (4 * 0.01 * x ** 2 - 2)
'''
j = 750    # number of points on X-axis.
i = 50     # number of iteration steps. (possibly not necessary)

lp = -50 # left limit
rp = 50   # right limit
DT = 0.01    # discrete time step
DX = (rp - lp) / j  # discrete spacial step
A = - DT # coeffi of the upper part of the matrix
B = - DT / (DX ** 2 )     # coeffi of the lower part of the matrix
x = np.linspace(lp, rp, j)     # initialize x axis
#phi = 0.1 * sin(2 * np.pi * x/ (rp - lp))        # initial condition
phi = exp( - 0.1 * x **2)
print(phi)
#for i in range(len(phi)):
#    if i < j * 0.25 or i > j * 0.75:
#        phi[i] = 0
dphi_dt = 0 * x
b = append(phi, dphi_dt)               # The entire right hand side
b_t = np.empty([2 * j, 1])
for m in range(len(b)):
    b_t[m][0] = b[m]
b = b_t
#print(b)
#x_n_1 # the solution of the next time step

'''
    Construct matrixes
'''
m_n_1 = np.identity(2 * j)            # initialize left matrix
m_n = np.identity(2 * j)          # initialize right matrix
for m in range(len(m_n_1)):
    for n in range(len(m_n_1[0])):
        if n - j == m:
            m_n_1[m][n] = A
        elif (m - j) == n:
            if n == 0:
                m_n_1[m][n] = 2 * B
                m_n_1[m][n + 1] = B
            elif n == j - 1:
                m_n_1[m][n] = 2 * B
                m_n_1[m][n - 1] = B
            else:
                m_n_1[m][n] = 2 * B
                m_n_1[m][n + 1] = B
                m_n_1[m][n - 1] = B
print('The left matrix should be:\n', m_n_1)
#print('The right matrix should be:\n', m_n)
#print('Test element of the matrix.', m_n_1[0][j])
#print(r)
#print('The unknown matrix is', sol)
#print('aM is', np.dot(m_n_1, sol))
#print('aM should be', r)


print('Check point 1')
'''
    Plot the result
'''
while i >= 0:
    print('Check point 2')
    b = np.linalg.solve(m_n_1, b)
    i -= 1
b_p = np.linspace(0, 1, 2 * j)
for m in range(len(b)):
    b_p[m] = b[m][0]
print(b_p)
p, = plt.plot(x, b_p[:j])
plt.xlim([lp, rp])
plt.ylim([-10, 10])
plt.show()
