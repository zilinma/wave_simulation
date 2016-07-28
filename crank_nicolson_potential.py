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
import pickle

lep = -50         # left limit
rep = 50         # right limit
np = 40000         # number of points
dim = 2 * np      # total dimension of the matrix
total_run_time = 50
dx = (rep - lep) / (np + 1)
dt = 0.0001
v = 1           # wave velocity / not implemented
d = dt/ (2 * dx ** 2)
x = linspace(lep + dx, rep - dx, np)
rho = 1 * x         # mass density function
A = 1000                  # amplitude of the wave function
# initialize phi
sigma = 0.1
start_point = -49.5

if start_point > rep or start_point < lep:
    raise RuntimeError

phi = exp(-((x - start_point)/sigma)**2)
psi = 200 * exp(-((x - start_point)/sigma) **2) * (x + start_point)
gamma = dt / 2
potential = 0 * x
#potential = [random.choice(range(10)) * gamma for x in x]
#exp(-((x + start_point)/sigma)**2) * 5
#potential *= gamma
print(potential)
for l in range(20):
    start_point = random.uniform(-1, 1)
    potential += exp(-((x + start_point)/0.05)**2) * random.choice([-1, 1]) * A

potential *= gamma

'''
print(len(potential))
for i in range(len(potential)):
    if i > 0.75 * np:
        potential[i]  = 2 * gamma
print(len(potential))
'''
if False:
    potential = 0 * x # initialize potential
v0 = append(phi, psi).flatten()

h1 = identity(np, float).toarray()

d1 = 2 * ones(np, float)
d2 = -ones(np - 1, float)
diagonals = [d1, d2, d2]
h2 = diags(diagonals, [0, -1, 1]).toarray()
potential_diag = diags(potential, 0).toarray()
z = zeros((np, np), float)

top = concatenate((z, -dt * h1 /2), axis = 1)
bottom = concatenate((d * h2, z), axis = 1)
h = concatenate((top, bottom), axis = 0)

top = concatenate((z, z), axis = 1)
bottom = concatenate((potential_diag, z), axis = 1)
p = concatenate((top, bottom), axis = 0)

h = csc_matrix(h)
p = csc_matrix(p)
lhs = identity(dim, float, format = 'csc') + h - p
rhs = identity(dim, float, format = 'csc') - h + p

solve = factorized(lhs)
nrhs = rhs.dot(v0)
v_old = v0

'''
p, = plt.plot(x, v0[:np])

plt.ylim([-1,1])
plt.xlim([lep, rep])
plt.ion()
T = 0

while True:
    p.set_ydata(v_old[:np])
    v0 = solve(rhs.dot(v_old))
    v_old = v0

    plt.draw()
    plt.show()

    T += dt
    print(T)
'''



for i in range(total_run_time / dt):
    v0 = solve(rhs.dot(v_old))
    v_old = v0

#p, = plt.plot(x, v_old[:np])
#potential *= 1/ (gamma * A)
#p, = plt.plot(x, potential)
#plt.xlim([lep, rep])
#plt.ylim([-1, 1])
#plt.xlabel(r'$x$ (t = 20000)', fontsize = 30)
#plt.ylabel(r'$\phi$', fontsize = 30)
#plt.tick_params(labelsize = 20)
#plt.show()

'''
    dump v_old
'''
pickle.dump(v_old, open('v.p', 'wb'))
pickle.dump(x, open('x.p', 'wb'))
