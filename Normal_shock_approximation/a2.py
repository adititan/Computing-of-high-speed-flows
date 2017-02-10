import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc,rcParams
import re

M1 = 10.0
n = 100
alpha = 0.5
CFL = 0.2

gamma = 1.4
Rsp   = 287                            # sp. gas constant R
Rspnd = 1/gamma                        # non-dimensional sp. gas constant

dp1   = 1.01325e5                      # Dimensional pressure
drho1 = 1.225                           # Dimensional density
dT1   = 288
a1  = np.sqrt(gamma*Rsp*dT1) # Dimensional temperature and freestream speed of sound 
u1   = M1                              # velocity non-dimensionalized by a1
rho1 = drho1/drho1                     # density  non-dimensionalized by its freestream value 1.225
T1   = dT1/dT1                         # temperature non-dimensionalized by its freestream value 288
p1   = dp1/(drho1*a1*a1) 

u   = []
rho = []
p   = []
T   = []
rhoe = []

M2 = np.sqrt((M1**2*(gamma-1)+2)/(2*gamma*M1**2-(gamma-1)))
p2 = p1*((1+gamma*M1**2)/(1+gamma*M2**2))
T2 = T1*((1+(gamma-1)/2*M1**2)*((2*gamma/(gamma-1)*M1**2-1)))/(M1**2*(2*gamma/(gamma-1)+(gamma-1)/2))
dT2 = T2*dT1
a2  = np.sqrt(gamma*Rsp*dT2)
u2 = M2*(a2/a1)
rho2 = u1/u2

L     = 35                 # Domain = [-5,30] shock being at x = 0
dx    = 35/(n-1)
dt    = CFL*dx/(abs(u1+1)) # Computing the time-step based on largest eigen value

x = []
for i in range(n):    
    x.append(-5 +(i)*dx)     
    if x[i]<0:           # free-stream conditions before shock   
        rho.append(rho1)
        u.append(u1)
        p.append(p1)
        T.append(T1)
    else:                 # conditions after shock   
        rho.append(rho2)
        u.append(u2)
        p.append(p2)
        T.append(T2)     
        rhoe.append(p[i]/(gamma-1)+0.5*rho[i]*u[i]**2)

#Exact Plots
plt.close("all")
plt.figure(figsize=(12.0,10.0))
plt.subplot(221)
plt.plot(x,rho)
plt.xlabel('x')
plt.ylabel(r'$\rho$')

plt.subplot(222)
plt.plot(x,u)
plt.xlabel('x')
plt.ylabel('u')

plt.subplot(223)
plt.plot(x,p)
plt.xlabel('x')
plt.ylabel('p')

plt.subplot(224)
plt.plot(x,T)
plt.xlabel('x')
plt.ylabel('T')

plt.savefig('exact.png',bbox_inches='tight')

#Alpha Variation
alpha = [0.9,0.6,0.4,0.3,0.2,0.15]

Rho = np.zeros((len(alpha),100))
U = np.zeros((len(alpha),100))
P = np.zeros((len(alpha),100))
T = np.zeros((len(alpha),100))
Err_a = []

plt.close("all")
plt.figure(figsize=(12.0,10.0))
for i,item in enumerate(alpha):
    f = open('dat_files/data_10_100_'+str(item)+'_0.2_.dat', "r")
    lines = f.readlines()[1:]
    sum_ = 0
    for j,line in enumerate(lines):
        line = re.sub(' +',' ',line)
        name = line.strip().split(" ")
        Rho[i][j]=(float(name[1]))
        U[i][j]=(float(name[2]))
        P[i][j]=(float(name[3]))
        T[i][j]=(float(name[4]))
        sum_ = sum_+(U[i][j]-u[j])**2
    Err_a.append(np.sqrt(sum_)/100)
    f.close()
    
    plt.subplot(221)
    plt.xlim([-5,10])
    plt.plot(x,Rho[i],label=r'$\alpha$ ='+str(item))
    plt.xlabel('x')
    plt.ylabel(r'$\rho$')
    plt.legend(loc = 'center right',shadow = True )
    
    plt.subplot(222)
    plt.xlim([-5,10])
    plt.plot(x,U[i],label=r'$\alpha$ ='+str(item))
    plt.xlabel('x')
    plt.ylabel('u')
    plt.legend(loc = 'upper right',shadow = True )
    
    plt.subplot(223)
    plt.xlim([-5,10])
    plt.plot(x,P[i],label=r'$\alpha$ ='+str(item))
    plt.xlabel('x')
    plt.ylabel('P')
    plt.legend(loc = 'center right',shadow = True )
    
    plt.subplot(224)
    plt.xlim([-5,10])
    plt.plot(x,T[i],label=r'$\alpha$ ='+str(item))
    plt.xlabel('x')
    plt.ylabel('T')
    plt.legend(loc = 'center right',shadow = True )
plt.savefig('alpha.png',bbox_inches='tight')

plt.close("all")
plt.scatter(alpha,Err_a)
plt.xlabel(r'$\alpha$')
plt.ylabel('Error')
plt.savefig('err_alpha.png',bbox_inches='tight')


#CFL Variation
cfl = [0.9,0.7,0.5,0.3,0.2,0.1]

Rho = np.zeros((len(cfl),100))
U = np.zeros((len(cfl),100))
P = np.zeros((len(cfl),100))
T = np.zeros((len(cfl),100))
Err_cfl = []

plt.close("all")
plt.figure(figsize=(12.0,10.0))
for i,item in enumerate(cfl):
    f = open('dat_files/data_10_100_0.6_'+str(item)+'_.dat', "r")
    lines = f.readlines()[1:]
    sum_ = 0
    for j,line in enumerate(lines):
        line = re.sub(' +',' ',line)
        name = line.strip().split(" ")
        Rho[i][j]=(float(name[1]))
        U[i][j]=(float(name[2]))
        P[i][j]=(float(name[3]))
        T[i][j]=(float(name[4]))
        sum_ = sum_+(U[i][j]-u[j])**2
    Err_cfl.append(np.sqrt(sum_)/100)
    f.close()
    
    plt.subplot(221)
    plt.xlim([-5,10])
    plt.plot(x,Rho[i],label='CFL = '+str(item))
    plt.xlabel('x')
    plt.ylabel(r'$\rho$')
    plt.legend(loc = 'center right',shadow = True )
    
    plt.subplot(222)
    plt.xlim([-5,10])
    plt.plot(x,U[i],label='CFL = '+str(item))
    plt.xlabel('x')
    plt.ylabel('u')
    plt.legend(loc = 'upper right',shadow = True )
    
    plt.subplot(223)
    plt.xlim([-5,10])
    plt.plot(x,P[i],label='CFL = '+str(item))
    plt.xlabel('x')
    plt.ylabel('P')
    plt.legend(loc = 'center right',shadow = True )
    
    plt.subplot(224)
    plt.xlim([-5,10])
    plt.plot(x,T[i],label='CFL = '+str(item))
    plt.xlabel('x')
    plt.ylabel('T')
    plt.legend(loc = 'center right',shadow = True )
plt.savefig('cfl.png',bbox_inches='tight')

plt.close("all")
plt.scatter(cfl,Err_cfl)
plt.xlabel('CFL')
plt.ylabel('Error')
plt.savefig('err_cfl.png',bbox_inches='tight')

#Grid Variation
N = [100,200,400,600,800,1000]

Rho = [[0.0 for i in range(0)] for j in range(len(N))]
U = [[0.0 for i in range(0)] for j in range(len(N))]
P = [[0.0 for i in range(0)] for j in range(len(N))]
T = [[0.0 for i in range(0)] for j in range(len(N))]
x = [[0.0 for i in range(0)] for j in range(len(N))]
Err_N = []

plt.close("all")
plt.figure(figsize=(12.0,10.0))
for i,item in enumerate(N):
    f = open('dat_files/data_10_'+str(item)+'_0.5_0.2_.dat', "r")
    lines = f.readlines()[1:]
    sum_ = 0
    for j,line in enumerate(lines):
        line = re.sub(' +',' ',line)
        name = line.strip().split(" ")
        x[i].append(float(name[0]))
        Rho[i].append(float(name[1]))
        U[i].append(float(name[2]))
        P[i].append(float(name[3]))
        T[i].append(float(name[4]))
        if float(name[0])<0:
            sum_ = sum_+(float(name[2])-u1)**2
        else:
            sum_ = sum_+(float(name[2])-u2)**2
    Err_N.append(np.sqrt(sum_)/item)
    f.close()
    
    plt.subplot(221)
    plt.xlim([-5,10])
    plt.plot(x[i],Rho[i],label='N ='+str(item))
    plt.xlabel('x')
    plt.ylabel(r'$\rho$')
    plt.legend(loc = 'center right',shadow = True )
    
    plt.subplot(222)
    plt.xlim([-5,10])
    plt.plot(x[i],U[i],label='N ='+str(item))
    plt.xlabel('x')
    plt.ylabel('u')
    plt.legend(loc = 'upper right',shadow = True )
    
    plt.subplot(223)
    plt.xlim([-5,10])
    plt.plot(x[i],P[i],label='N ='+str(item))
    plt.xlabel('x')
    plt.ylabel('P')
    plt.legend(loc = 'center right',shadow = True )
    
    plt.subplot(224)
    plt.xlim([-5,10])
    plt.plot(x[i],T[i],label='N ='+str(item))
    plt.xlabel('x')
    plt.ylabel('T')
    plt.legend(loc = 'center right',shadow = True )
plt.savefig('grid.png',bbox_inches='tight')

plt.close("all")
plt.scatter(N,Err_N)
plt.xlabel('Grid Size')
plt.ylabel('Error')
plt.savefig('err_grid.png',bbox_inches='tight')
