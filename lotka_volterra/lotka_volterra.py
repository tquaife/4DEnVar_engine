"""
Some code taken from here:

https://scientific-python.readthedocs.io/en/latest/notebooks_rst/3_Ordinary_Differential_Equations/02_Examples/Lotka_Volterra_model.html
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.optimize import minimize
from copy import copy

def lv_derivative(X, t, alpha, beta, delta, gamma):
    x, y = X
    dotx = x * (alpha - beta * y)
    doty = y * (-delta + gamma * x)
    return np.array([dotx, doty])

def lv_forward_run(X0, alpha, beta, delta, gamma, dt, tmax):
    
    t=get_time_steps(dt, tmax)
    res = integrate.odeint(lv_derivative, X0, t, args = (alpha, beta, delta, gamma))
    x1, x2 = res.T
    
    return [x1,x2]

def get_time_steps(dt, tmax):
    nt = int(np.ceil(tmax/float(dt)))
    t = np.linspace(0.,tmax, nt)
    return t
        
def gen_obs(X, nobs, dt, tmax):        
        
    t=get_time_steps(dt, tmax)
    x1_obs=np.zeros(nobs)        
    x2_obs=np.zeros(nobs)        
    t_obs=np.zeros(nobs)        
    obs_intrvl=len(t)/nobs    
     
    x1,x2=X 
    
    for i in range(nobs):
    
        ti=int((i*obs_intrvl)+(obs_intrvl/2.)) 
        t_obs[i]=t[ti]
        t_obs[i]=ti
        x1_obs[i]=x1[ti]
        x2_obs[i]=x2[ti]
                
    return x1_obs, x2_obs, t_obs            



#def J(X, Y, Xb, hX, Binv, Rinv):
def J(X, alpha, beta, delta, gamma, Y, Xb, Binv, Rinv, dt, tmax, nobs):
#def J(X, args):

    x1_0=X[0]
    x2_0=X[1]
    #alpha=X[2]
    #beta=X[3]
    #gamma=X[4]
    #delta=X[5]

    #Y=args[0]
    #Xb=args[1] 
    #Binv=args[2] 
    #Rinv=args[3]
    #dt=args[4]
    #tmax=args[5]
    #nobs=args[6]

    [hx1_full_run,hx2_full_run]=lv_forward_run([x1_0,x2_0], alpha, beta, delta, gamma, dt, tmax )
    hx1,hx2,t_obs=gen_obs([hx1_full_run,hx2_full_run], nobs, dt, tmax) 
    hX=np.array([hx1,hx2]).flatten()
   
    #print(np.shape(Y))
    #print(np.shape(hX),np.shape(hx1))
    #print(np.shape(Rinv))

    J1=np.matmul((Xb-X).T, np.matmul(Binv,(Xb-X)))
    J2=np.matmul((Y-hX).T, np.matmul(Rinv,(Y-hX)))

    return 0.5*(J1+J2)


def J_test(X, Y, Xb, Binv, Rinv):

    J1=np.matmul((Xb-X).T, np.matmul(Binv,(Xb-X)))
    J2=np.matmul((Y-X).T, np.matmul(Rinv,(Y-X)))

    return 0.5*(J1+J2)

def J_grad_test(X, Y, Xb, Binv, Rinv):

    J1=np.matmul(Binv,(Xb-X))
    J2=np.matmul(Rinv,(Y-X))

    return J1-J2

                
def basic_example():
    """This is the odeint example from:     
    https://scientific-python.readthedocs.io/
    """
    
    alpha = 1. #mortality rate due to predators
    beta = 1.
    delta = 1.
    gamma = 1.
    x0 = 4.
    y0 = 2.

    Nt = 1000
    tmax = 30.
    t = np.linspace(0.,tmax, Nt)
    X0 = [x0, y0]
    res = integrate.odeint(lv_derivative, X0, t, args = (alpha, beta, delta, gamma))
    x, y = res.T
       
    plt.figure()
    plt.grid()
    plt.title("odeint method")
    plt.plot(t, x, '-b', label = 'Deer')
    plt.plot(t, y, '-r', label = "Wolves")
    plt.xlabel('Time t, [days]')
    plt.ylabel('Population')
    plt.legend()

    plt.show()

if __name__=="__main__":    
    
    alpha = 1. #mortality rate due to predators
    beta = 2.
    delta = 1.
    gamma = 1.
    x1_0 = 4.
    x2_0 = 2
    
    dt=0.1
    tmax=100
    nobs=30
    
    
    Xtrue=np.array([x1_0, x2_0])    
    Xb=np.array([x1_0, x2_0])
    X0=np.array([x1_0+0.5, x2_0])
    
    [y1_full_run,y2_full_run]=lv_forward_run(Xtrue, alpha, beta, delta, gamma, dt, tmax )
    y1,y2,t_obs=gen_obs([y1_full_run,y2_full_run], nobs, dt, tmax)        
    Y=np.array([y1,y2]).flatten()
    Binv=np.eye(len(Xb))*0.01
    Rinv=np.eye(nobs*2)
     
    args=(alpha, beta, delta, gamma, Y,Xb,Binv,Rinv,dt,tmax,nobs)    
    print(J(X0,alpha, beta, delta, gamma, Y, Xb, Binv, Rinv, dt, tmax, nobs))

    res=minimize(J, X0, args=args,method='L-BFGS-B',jac='3-point')    
    print(res.x)
    print(J(res.x, alpha, beta, delta, gamma, Y, Xb, Binv, Rinv, dt, tmax, nobs))
    [hxa1_full_run,hxa2_full_run]=lv_forward_run(res.x, alpha, beta, delta, gamma, dt, tmax )

    [hxb1_full_run,hxb2_full_run]=lv_forward_run(X0, alpha, beta, delta, gamma, dt, tmax )



    if True:
        plt.figure()
        plt.grid()
        plt.plot(y1_full_run, '-b', label = 'x1')
        plt.plot(y2_full_run, '-r', label = "x2")
        plt.plot(t_obs,y1,'ob',label="y1")
        plt.plot(t_obs,y2,'or',label="y2")
        plt.plot(hxb1_full_run, '--b', label = 'xb1', alpha=0.3)
        plt.plot(hxb2_full_run, '--r', label = "xb2", alpha=0.3)
        plt.plot(hxa1_full_run, '-.b', label = 'xa1')
        plt.plot(hxa2_full_run, '-.r', label = "xa2")

        plt.xlabel('time')
        plt.ylabel('state')
        plt.legend()
        plt.show()


