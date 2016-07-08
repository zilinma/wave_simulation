'''
    Simple evolvement for a function.
'''
import random
import numpy as np
import matplotlib.pyplot as plt
from visual import *
# constance
DT = 0.05
E = 2.71828
plt.ion()

# innitialize x
lep =  -100.0
rep =  100.0
np = 501
x = linspace(lep, rep, np)
print(x)
DX = (rep - lep)/ (np)
a = DT/DX
'''
    Initial conditions
'''
u = 0.1 * exp(- 0.01 * x ** 2)
print(u)
#0.5 * sin(x)
s = 0 * x       # du/dt
r = 0.1 * exp( - 0.01 * x **2) * ( - 2 * 0.01 * x)
#0.1 * cos(x)
      # du/dx
v = 0 * x + 0.0005 # potential function

# iterate
p, = plt.plot(x, u)
plt.ylim([-1, 1])
plt.xlim([-100, 100])

while True:
    rate(1000)
    p.set_ydata(u)

    for i in range(len(u)):
        '''
            Update u
        '''
        u[i] +=  DX * (r[(i + 1) % np] - r[(i)])

        '''
            update s
        '''
        s[i] += a * ((s[(i + 1) % np] - s[(i - 1)]) / 2 + a * (r[(i + 1) % np] - 2 * r[i] + r[(i - 1)]) / 2)
        '''
            update r
        '''
        r[i] += a * ((r[(i + 1) % np] - r[(i - 1)]) / 2 + a * (s[(i + 1) % np] - 2 * s[i] + s[(i - 1)]) / 2)

    '''
        Calculate energy
    '''
    energy = 0
    for i in range(len(u)):
        energy += DX * ((1/2) * (s[i]) ** 2)
        energy += DX * (1 + (r[i]) ** 2) ** (1/2)

    print(energy)

    '''
        Graph u
    '''
    plt.draw()

