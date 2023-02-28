import numpy as np


######## ***************************************  
## Nummerische Differenttion
######## ***************************************  

def Forward_diff(f,x):   
    # Vorwärtsdifferenzen 
    h_opt = (np.finfo(float).eps)**(1/2)     
    h = np.zeros(len(x))
    gradient = np.zeros(len(x))
    for i in range(len(x)):
        h   = h*0
        h[i] += h_opt 
        gradient[i] = (f(x+h)-f(x))/h_opt        


    return gradient


######## ***************************************  
## Nummerische Integration
######## ***************************************  

def runge_kutta_k4(f,x,u):
    #RK4 integration with zero-order hold on u
    h = np.ones(len(x))*0.05      # time step (20Hz) 

    f1 = f(x, u)
    f2 = f(x + 0.5*h*f1, u)
    f3 = f(x + 0.5*h*f2, u)
    f4 = f(x + h*f3, u)

    return x + (h/6.0)*(f1 + 2*f2 + 2*f3 + f4)

######## ***************************************  
## Optimierungsverfahren
######## ***************************************  

def quasiNewton(x,x_alt,df, B_alt):
    # Approzimierung der Hesse Matrix
    epsilon = 1e-10
    d = x -x_alt
    y = df(x)-df(x_alt)

    rho = 1/(y*d+epsilon)
    # aktualisieren der Inversen
    B = (1-rho*d*y)*B_alt*(1-rho*y*d) + rho*d**2
            
    # Regularisierung
    beta = 1
    while B <= 0:
            B += beta   
    s = -B*df(x)
                   
        


    # Schrittweite bestimmen mittels Powell-Wolfe Verfahren
    alpha = 1
    c = 0.5
    b = 1
    k = 0
    while f(x +alpha*s) > f(x) + b*alpha*df(x) and k<10:
        alpha = alpha*c
        k+=1
    #print(s*alpha)
    #alpha = 1
    return x + alpha*s , B

def quasiNewton_3D(x,x_alt,df, B_alt):
    # Approzimierung der Hesse Matrix
    epsilon = 1e-30
    d = np.array([x -x_alt])
    y = np.array([df(x)-df(x_alt)])

    rho = 1/(y@d.T +epsilon)
    # aktualisieren der Inversen
    B = (np.eye(len(x))-rho*d.T@y)*B_alt*(np.eye(len(x))-rho*y.T@d) + rho*d.T@d
            
    # Regularisierung
    beta = np.eye(len(x))*0.01
    while not np.all(np.linalg.eigvals(B)):
            B += beta       
    s = -B@df(x)

 
    # Schrittweitenstrategie
    
    # Amerigo
    '''
    b = 1
    alpha = 1
    c= 0.5 
    while f(x +alpha*s) >= f(x) + b*alpha*np.linalg.norm(df(x)):
        alpha = alpha*c        
    '''
    # Wolf Powel, vereinfacht
    
    alpha = 1
   
    b= 0.5
    while f(x+alpha*s)>f(x):
        alpha = alpha*b        

    return x + alpha*s , B


def Linesearch(x0,f,df,ddf,N,e_x = 1e-5, e_f = 1e-15):
    """ Animation
        Params
         --------
        x0:            Startwerte
        e_x:           Abbruchkriterium für die Schrittweite
        e_f:           Abbruchkriterium für den Funktionswert
        f:             Funktion
        df:            Ableitung der Funktion
        Returns
        --------     

        x:              Minimum     
        k:              steps
        X:        Step size                                       
    """

    x = x0  
    x_alt = x0    
    k = 0
    B0 = 1 #np.eye(np.size(x))  
    x = x0 - 1e-3*Forward_diff(f,x0)
    X = np.zeros(N)
    while k < N:       
              
        
        #x = NewtonStep(x_alt,df,ddf)    
        #x = bactracking_reg_newton(x_alt,df,ddf)  
        #x,B0 = quasiNewton(x,x_alt,df,B0)
        x,B0 = quasiNewton_3D(x,x_alt,df,B0)       
        
        #ax.scatter(x[0],x[1],f(x), c='r', marker='o', s = 50)        
                    
        X[k] = np.linalg.norm(x_alt-x)
        k += 1 
        # Abbruchkriterium
        if np.linalg.norm(x_alt-x) < e_x or np.linalg.norm(df(x_alt)-df(x)) < e_f:            
            return x,k,X

        x_alt = x 
        
        
            

    return x,k,np.linalg.norm(x_alt-x),X