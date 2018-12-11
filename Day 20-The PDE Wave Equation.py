#Tarik Dahnoun
#Worked with Alex Towe
#Wave equation PDE
# November 28, 2017

import numpy as np
import pylab as py

# Initial Values
mu = 0.01
T = 25.0
c = np.sqrt(T/mu)
print 'verify c speed',c
xa = -5.0
xb = 5.0
J = 100
N = J
sigma = 1.4 
A = 5.3
x = np.linspace(xa, xb, J+1)
dx = x[1]-x[0]

# CFL_condition
#order unity is .25
nuu=.25
dt = nuu*dx/c


y = np.zeros([N+1, J+1])
y1 = np.zeros([N+1, J+1])
y2 = np.zeros([N+1, J+1])
v = np.zeros([N+1, J+1])
v2 = np.zeros([N+1, J+1])
x2 = np.zeros(J+1)
t = np.zeros(N+1)
yt = np.zeros(J+1)
vt = np.zeros(J+1)

# Function for initial amplitude
def f(x):
    return A*np.exp((-x**2)/(sigma**2))

# Function for initial velocity 
def g(x):
    return (2*A*c*x /sigma**2)*np.exp((-x**2)/sigma**2)
    
def PDE(y, v, N, J):    
    for i in range(0, N):
        # First and last amplitudes and velocities for each timestep don't change
        y[i+1,0] = y[i,0]
        y[i+1,J] = y[i,J]
        v[i+1,0] = 0
        v[i+1,J] = 0
        
        # Find amplitudes and velocites for each half-timestep
        for j in range(1, J):
            yt[j]=y[i,j]+0.5*dt*v[i,j]
            vt[j]=v[i,j]+0.5*dt*(c**2)*((y[i,j+1]-2*y[i,j]+y[i,j-1])/dx**2)        
        
        # Use half-timesteps to find amplitudes and velocites for each timestep
        for k in range(1, J):
            y[i+1,k]=y[i,k]+dt*vt[k]
            v[i+1,k]=v[i,k]+dt*(c**2)*((yt[k+1]-2*yt[k]+yt[k-1])/dx**2)
    return y

def plot(num, x, y): 
    py.figure(num)
    for z in range(100):
        py.plot(x,y[z])
    py.title('Various Time Steps')
    py.xlabel('x')
    py.ylabel('y(x)')
    py.legend()

# Find values for initial timestep
for i in range(0, J+1):
    t[i]=i*dt
    x2[i]=xa+i*dx
    y[0,i]=f(x2[i])
    
    # Initial velocities for part 2
    v[0,i]=0
    
    # Initial velocites for part 3
    v2[0,i]=g(x2[i])

# print v2[0,0]

# For part 1
y1 = PDE(y, v, N, J)
plot(1, x2, y1)

# For part 2
y2 = PDE(y, v2, N, J)
plot(2, x2, y2)

py.show()