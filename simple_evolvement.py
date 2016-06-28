'''
    Simple evolvement for a function.
    Printed values are: 1. The energy of the system, which should remain constant. 2. The heigh of the central point. 
    Should be 0 after the wave is restored. 
'''
import random
import numpy as np
import matplotlib.pyplot as plt
from visual import *
# constance
DT = 0.1
E = 2.71828
plt.ion()

# innitialize x
'''
100, 500 most workable
'''
lep = -100.0
rep = 100.0
np = 500
x = linspace(lep, rep, np)
DX = (rep - lep)/ (np - 1)

'''
    Initial conditions
'''
'''
initialize phi
Good initial condition 1 (DT = 0.1, np = 1000/500/250(wierd), np/lp = +-1000)
phi = E ** (-0.01 * x ** 2)
semi-good initial condition
phi = sin(x)
'''
phi = E ** (-0.01 * x ** 2)

''' parabola is fun

phi = -x ** 2 + 100
'''
'''
phi = 0 * x
for i in range(len(phi)):
    if i < 0.75 * np and i > 25 * np:
        phi[i] = 1
    else:
        phi[i] = 0
'''

'''
initialize dphi/dt
Good initial condition 1
dphi_dt = 0 * x
'''
dphi_dt = 0 * x

'''
define d2phi_dx2
Good initial condition 1
d2phi_dx2 = (E ** (-0.01 * x ** 2)) * (4 * 0.01 * x ** 2 - 2)
semi-good initial condition
d2phi_dx2 = -sin(x)
'''
d2phi_dx2 = (E ** (-0.01 * x ** 2)) * (4 * 0.01 * x ** 2 - 2)
''' parabola is fun
d2phi_dx2 = -2 * (dphi_dt + 1)
'''

'''
    initialize potential
'''
v = 0 * x + 0.0005 # potential function
p = len(v)
for i in range(len(v)):
    if i < p * 3/4 and i > p/4:
        v[i] = 0
# iterate

p, = plt.plot(x, phi)
plt.ylim([-1, 1])
plt.xlim([-100, 100])

while True:
    rate(1000)
    '''
        Energy ploting
    '''
    p.set_ydata(phi)


    dphi_dx = 0 * x
    for i in range(len(phi)):
        '''
            Update phi (possibly not accurate)
        '''
        '''
            RK4
        '''

        phi[i] = phi[i] + (dphi_dt[(i + 1) % np] + dphi_dt[i - 1]) * DT / 2  # Lax for phi
        dphi_dt[i] = ((d2phi_dx2[(i + 1) % np] + d2phi_dx2[i - 1]) / 2) * DT + \
                     (dphi_dt[i - 1] + dphi_dt[(i + 1) % np]) / 2 - v[i] * phi[i]
                         #dphi_dt[i]  # Lax for phi_t
        '''
            Update d2phi_dx2
        '''

            #d2phi_dx2[i] = (phi[i + 2] - 2 * phi[i + 1] + phi[i]) / (DX ** 2)
        d2phi_dx2[i] = (((phi[(i + 1)%np] + phi[(i + 3)%np])/2 - 2 * (phi[i] + phi[(i + 2)%np])/2 + phi[i]) / (DX ** 2))
        '''
            Calculate dphi_dx
        '''
        dphi_dx[i] = (phi[i] - phi[(i + 1)%np]) / DX

    '''
        Calculate energy
    '''
    energy = 0
    for i in range(len(phi)):
        energy += DX * ((1/2) * (dphi_dt[i]) ** 2)
        energy += DX * (1 + (dphi_dx[i]) ** 2) ** (1/2)

    print(energy)
    print(phi[(np)/2])

    '''
        Graph phi
    '''
    plt.draw()
