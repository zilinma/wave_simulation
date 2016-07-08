#!/usr/remote/python-3.2/bin/python
'''
    Numerical solution of wave equation using Crank Nicolson method.
    author@Zilinma-bucknell
'''
from scipy import*
from numpy import*
from scipy.sparse.linalg import factorized
from scipy.sparse import lil_matrix
from scipy.sparse import csc_matrix
from scipy.sparse import linalg
from scipy.sparse import identity
from scipy.sparse.linalg import spsolve
from scipy.sparse import diags
import matplotlib.pyplot as plt

lep = -10         # left limit
rep = 10         # right limit
np = 1999         # number of points
dim = 2*np      # total dimension of the matrix

dx = (rep - lep) / (np + 1)
dt = 0.01
d = dt/ (2 * dx ** 2)
x = linspace(lep + dx, rep - dx, np)

# initialize phi /psi
sigma = 1
phi = exp(-((x)/sigma)**2)
psi = 0 * x
v0 = append(phi, psi).flatten()

h1 = identity(np, float).toarray()

d1 = 2 * ones(np, float)
d2 = -ones(np - 1, float)
diagonals = [d1, d2, d2]
h2 = diags(diagonals, [0, -1, 1]).toarray()

z = zeros((np, np), float)

top = concatenate((z, -dt * h1 /2), axis = 1)
bottom = concatenate((d * h2, z), axis = 1)
h = concatenate((top, bottom), axis = 0)

h = csc_matrix(h)
lhs = identity(dim, float, format = 'csc') + h
rhs = identity(dim, float, format = 'csc') - h

solve = factorized(lhs)
nrhs = rhs.dot(v0)
v_old = v0

p, = plt.plot(x, v0[:np])
plt.ylim([-1,1])
plt.xlim([lep, rep])
plt.ion()
while True:
    p.set_ydata(v_old[:np])
    v0 = solve(rhs.dot(v_old))
    v_old = v0
    plt.draw()
    plt.show()

