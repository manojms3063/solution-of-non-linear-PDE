import numpy as np
from scipy.linalg import solve

# v1 = 1  # Initial condition for y
# w1= 1  # Initial condition for y'
# t0 = 0      # Initial time
# tf = 10     # Final time
Nx = 10
Nt = 10    # Number of time steps
# tol=1e-6
# n =100  #max_iter

nu = 0.1  # viscosity parameter
L = 10.0  # length of spatial domain
T = 2.0   # total time
dx = L/(Nx - 1)
dt = T/Nt 
# t = np.linspace(t0, tf, Nx)
# x = np.linspace()
# print(t)
# dt = t[1] - t[0]
# dx = 
# print(dt)
# Initial guess
v = np.zeros(Nx)
w = np.zeros(Nx)
v[0] = 0
w[0] = 0
x = np.linspace(0, L, Nx)

def ic(x):
    return -0.2 * np.sin(np.pi * x / L)

for i in range(Nx):
    v[i] = ic(x[i])
    w[i] = -0.2*np.pi/L*np.cos(np.pi*x[i]/L)

for i in range(Nx-1):
    v_guess = v[i]
    w_guess = w[i]
    
    # for k in range(n):
    F = np.array([v[i+1] - dt*w[i] - v[i],v[i+1]-v[i]*(1-dt*(w[i+1]-w[i])/dx)])
    J = np.array([[1, -dt],[-1, v[i]*(dt/dx)]])
    J=np.linalg.inv(J)
    delta = np.dot(J, -F)
    v_guess += delta[0]
    w_guess += delta[1]
        # if np.linalg.norm(delta) < tol:
            # break

    # v[i+1] = v_guess
    # w[i+1] = w_guess
    m = [v_guess,w_guess]
    m = m+delta
    print(m)
    # print(delta)
# # print("t values:", t)
# print(" v values:", v)
# print(" w values:", w)


