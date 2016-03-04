import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os,sys,time
import re


class doit():
    def __init__(self):
        return

def get_inputs(g1):
    g1.progress = 0
    g1.xi = -2.0
    g1.xf = 2.0
    g1.dx = 0.1
    g1.dt = 0.1
    g1.cfl = 0.5 
    g1.tsteps = 50
    g1.method = 2 
    return

def define_domain(g1):
    g1.msg = 'Defining domain'
    g1.progress = g1.progress + 1
    g1.jmax = int(abs(g1.xf - g1.xi)/g1.dx) + 1  
    g1.x = np.linspace(g1.xi,g1.xf,g1.jmax)
    g1.u = np.zeros((g1.jmax,g1.tsteps+1))
    return

def define_initial_condn(g1):
    g1.msg = 'Initial condition'
    # g1.u[:,0] = np.sin(2*g1.x)
    for j in range(g1.jmax):
        if g1.x[j] < 0.0:
            g1.u[j,0]=0.0
        else:
            g1.u[j,0]=1.0
    return

def initialize(g1):
    g1.msg = '...'
    g1.progress = g1.progress + 1
    return

def take_step(g1):
    t = g1.step*g1.dt
    n = g1.step - 1
    c = g1.cfl*g1.dx/g1.dt
    if(g1.method == 1):
        # exact solution
        for j in range(1,g1.jmax-1):
            g1.u[j,n+1] = np.sin(2*(g1.x[j] + c*t))
    if(g1.method == 2):
        # upwind method
        for j in range(1,g1.jmax-1):
            g1.u[j,n+1] =g1.u[j,n] + g1.cfl*(g1.u[j+1,n] - g1.u[j,n])
    return

def apply_bc(g1,xind,bctype):
    g1.msg = 'Applying BC'
    tind = g1.step
    t = g1.step*g1.dt
    c = g1.cfl*g1.dx/g1.dt
    if(bctype == 1):
        # free boundary condition
        g1.u[xind,tind]=np.sin(2*(g1.x[xind] + c*t))
    if(bctype == 2):
        # fixed boundary condition
        g1.u[xind,tind]=g1.u[xind,0]
    return

def plot_results(g1):
    plt.ion()
    plt.figure()
    for tind in range(g1.tsteps+1):
        plt.xlim(-2,2)
        plt.ylim(-1.5,1.5)
        plt.grid()
        plt.plot(g1.x,g1.u[:,tind])
        plt.draw()
        plt.clf()
    return

if __name__ == '__main__':
    g1 = doit()
    get_inputs(g1)
    define_domain(g1)
    define_initial_condn(g1)
    for n in range(g1.tsteps):
        g1.step = n+1
        take_step(g1)
        apply_bc(g1,0,2)
        apply_bc(g1,-1,2)
    plot_results(g1)
