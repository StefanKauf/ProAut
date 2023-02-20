"""
Date: 19.12.2022
Author: Kaufmann Stefan

RNL Aufgabe 3.1
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from Plot import phasenplot



def model(x,t):
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
    u = eingang(x,t)
    
    dx[0] = x[0]*x[1]
    dx[1] = x[0] + u
    return dx


def eingang(x,t,k=[1,1]):
    """ System Input
        Params
         --------
        x:             steady states as [x1,x2]
        t:             time as int
        k:             Controler Gain k = [k1,k2]   mit k1,k2 > 0       
                              

        Returns
        --------
        u:              System input   
                
    """
    k1 = k[0]
    k2 = k[1]
    x1 = x[0]
    x2 = x[1]

    #u = -k2*(x2+k1)-x1-x1**2  # Version 1
    u = -k2*(x2+k1*(x1**2)) - x1 -x1**2 -2*k1*x1**2*(x2+k1*x1**2)

    return u 

# set the initial conditions
x0=[2,2]

# define the discretization points
t_start = 0
t_stop = 10
dt = 1/1000

t_sim=np.linspace(t_start, t_stop, int((t_stop - t_start) / dt + 1))



solutionOde=odeint(model,x0,t_sim)

plt.plot(t_sim, solutionOde[:, 0], 'b', label='x1')
plt.plot(t_sim, solutionOde[:, 1], 'g', label='x2')
plt.legend(loc='best')
plt.xlabel('time')
plt.ylabel('x1(t), x2(t)')
plt.grid()
plt.savefig('simulation.png')
plt.show()


#phasenplot(model,limit=[0,3],N=10,T=10)









    
 


