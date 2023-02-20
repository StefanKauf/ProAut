"""
Date: 19.12.2022
Author: Kaufmann Stefan

RNL Aufgabe 3.1
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from Plot import phasenplot



def model(x,t,k1=0.1,k2=0.5):
    """ Nonlinear System Model
        Params
         --------
        x:             steady states as [x1,x2]
        t:             time as int
        k1:            stiffnes joint 1
        k2:            stiffnes joint 2
                              

        Returns
        --------
        dx:       chance of the state as a vektor [dx1, dx2]     
                
    """
    dx =[0,0]  
    
    dx[0] = x[1]
    dx[1] = -k1*x[1]-k2*np.sin(x[0])    
    return dx

# set the initial conditions
x0=[0,1]

# define the discretization points
timePoints=np.linspace(0,100,300)

#define the constants 
k1=0.1
k2=0.5

solutionOde=odeint(model,x0,timePoints, args=(k1,k2))

plt.plot(timePoints, solutionOde[:, 0], 'b', label='x1')
plt.plot(timePoints, solutionOde[:, 1], 'g', label='x2')
plt.legend(loc='best')
plt.xlabel('time')
plt.ylabel('x1(t), x2(t)')
plt.grid()
plt.savefig('simulation.png')
#plt.show()


phasenplot(model,limit=[-10,10])









    
 


