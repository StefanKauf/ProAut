"""
Date: 20.12.2022
Author: Kaufmann Stefan

Hilfreiche Funktionen
"""


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import control as ctr




## ********************* Phase Plot   ************************
def phasenplot(model,limit = [-5,5],N=5,T=50):
    """ Nonlinear System Model
        Params
         --------
        model          second order ode dx = [x1,x2]     
        limit          Size limitation
        N              Initial Points N
        T              Simulation Time               

        Returns
        --------
                    
            
                
    """
    from control.phaseplot import phase_plot

  
    plt.figure()
    plt.clf()
    plt.axis([limit[0], limit[1], limit[0], limit[1]])
    plt.title('Phase plot')
    plt.xlabel('$x_1$')
    plt.ylabel('$x_2$')




    X0 = np.array([0,0])
    for i in np.linspace(limit[0],limit[1],N):
        for j in np.linspace(limit[0],limit[1],N):
            X0 = np.append(X0,np.array([i,j]))
    x0 = np.reshape(X0,(-1,2))

    phase_plot(model, X0=x0, T=np.linspace(0, T, T*20), lingrid=0)
    plt.show()
