import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc,rcParams
import re
M1 = 10.0
alpha = 0.5
n = 800
CFL = 0.2

gamma = 1.4
Rsp   = 287.0                            # sp. gas constant R
Rspnd = 1.0/gamma                        # non-dimensional sp. gas constant

dp1   = 1.01325e5                      # Dimensional pressure
drho1 = 1.225                           # Dimensional density
dT1   = 288
a1  = np.sqrt(gamma*Rsp*dT1) # Dimensional temperature and freestream speed of sound 
u1   = M1*a1                              # velocity non-dimensionalized by a1
rho1 = drho1                     # density  non-dimensionalized by its freestream value 1.225
T1   = dT1                         # temperature non-dimensionalized by its freestream value 288
p1   = dp1 

u   = []
rho = []
p   = []
T   = []
rhoe = []

M2 = np.sqrt((M1**2*(gamma-1)+2)/(2*gamma*M1**2-(gamma-1)))
p2 = p1*((1+gamma*M1**2)/(1+gamma*M2**2))
T2 = T1*((1+(gamma-1)/2*M1**2)*((2*gamma/(gamma-1)*M1**2-1)))/(M1**2*(2*gamma/(gamma-1)+(gamma-1)/2))
dT2 = T2
a2  = np.sqrt(gamma*Rsp*dT2)
u2 = M2*(a2)
rho2 = u1/u2*rho1

L     = 8                 # Domain = [0,8] shock being at x = 0.5
dx    = 8.0/(n-1)
dt    = CFL*dx/(abs(u1+1)) # Computing the time-step based on largest eigen value

x = []


for i in range(n):    
    x.append((i)*dx)     
    if x[i]<0.5:           # free-stream conditions before shock   
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


x_exact = x
T_exact = T
U_exact = u
P_exact = p
Rho_exact = rho

N = [200,400,800]


# steger waming data extraction
def Steger_warming(filename):
    f = open(filename, "r")
    lines = f.readlines()[7:]
    A = []
    for line in lines: 
        name = line.strip().split(" ")
        if len(name) == 2:
            break
        for item in name:
            A.append(float(item))
    f.close()
    return A
Rho_sw = [[0.0 for i in range(0)] for j in range(len(N))]
U_sw = [[0.0 for i in range(0)] for j in range(len(N))]
P_sw = [[0.0 for i in range(0)] for j in range(len(N))]
T_sw = [[0.0 for i in range(0)] for j in range(len(N))]
Err_sw = [[0.0 for i in range(0)] for j in range(len(N))]

for i,item in enumerate(N):
    Rho_sw[i] = Steger_warming("r_" + str(item)+".dat")[:item]
    P_sw[i] = Steger_warming("p_" + str(item)+".dat")[:item]
    U_sw[i] = Steger_warming("u_" + str(item)+".dat")[:item]
    T_sw[i] = Steger_warming("t_" + str(item)+".dat")[:item]
    x = np.linspace(0,8,item)
    sum_ = 0
    for j in range(item):
        if x[j]<0.5:
            Err_sw[i].append(abs(U_sw[i][j]-u1)/u1)
        else:
            Err_sw[i].append(abs(U_sw[i][j]-u2)/u1)
    
# Lax-Fredrich Data extraction
Rho = [[0.0 for i in range(0)] for j in range(len(N))]
U = [[0.0 for i in range(0)] for j in range(len(N))]
P = [[0.0 for i in range(0)] for j in range(len(N))]
T = [[0.0 for i in range(0)] for j in range(len(N))]
x = [[0.0 for i in range(0)] for j in range(len(N))]
Err_Lax = [[0.0 for i in range(0)] for j in range(len(N))]


for i,item in enumerate(N):
    f = open('data_10_'+str(item)+'_0.5_0.2_.dat', "r")
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
        if float(name[0])<0.5:
            Err_Lax[i].append(abs(float(name[2])-u1)/u1)
        else:
            Err_Lax[i].append(abs(float(name[2])-u2)/u1)
    f.close()
    plt.figure()
    plt.plot(x[i],Err_sw[i],label="Steger Warming")
    plt.plot(x[i],Err_Lax[i],label="Lax Friedrich")
    plt.xlabel('x')
    plt.ylabel('Relative Error in "u"')
    plt.legend()
    plt.xlim([0.25,0.75])
    plt.savefig('Err_'+str(item)+'.png',bbox_inches='tight')

plt.close("all")
for i in range(3):
    plt.figure(i,figsize=(12.0,10.0))
    plt.subplot(221)
    plt.plot(x[i],U_sw[i],label = 'Steger Warming')
    plt.plot(x[i],U[i],label='Lax Friedrich')
    plt.plot(x_exact,U_exact,label='Exact')
    plt.xlim([0.25,1.25])
    plt.xlabel('x')
    plt.ylabel('u')
    plt.legend(loc = 'center right',shadow = True )
    
    plt.subplot(222)
    plt.plot(x[i],T_sw[i],label = 'Steger Warming')
    plt.plot(x[i],T[i],label='Lax Friedrich')
    plt.plot(x_exact,T_exact,label='Exact')
    plt.xlim([0.25,1.25])
    plt.xlabel('x')
    plt.ylabel('Temperature')
    plt.legend(loc = 'center right',shadow = True )
    
    plt.subplot(223)
    plt.plot(x[i],P_sw[i],label = 'Steger Warming')
    plt.plot(x[i],P[i],label='Lax Friedrich')
    plt.plot(x_exact,P_exact,label='Exact')
    plt.xlim([0.25,1.25])
    plt.xlabel('x')
    plt.ylabel('Pressure')
    plt.legend(loc = 'center right',shadow = True )
    
    plt.subplot(224)
    plt.plot(x[i],Rho_sw[i],label = 'Steger Warming')
    plt.plot(x[i],Rho[i],label='Lax Friedrich')
    plt.plot(x_exact,Rho_exact,label='Exact')
    plt.xlim([0.25,1.25])
    plt.xlabel('x')
    plt.ylabel('Density')
    plt.legend(loc = 'center right',shadow = True )

plt.figure(0)
plt.savefig("N_200.png",bbox_inches='tight')
plt.figure(1)
plt.title('Solution for 400 Grid points')
plt.savefig("N_400.png",bbox_inches='tight')
plt.figure(2)
plt.title('Solution for 800 Grid points')
plt.savefig("N_800.png",bbox_inches='tight')

Err_Lax_sum = []
Err_sw_sum = []
for i,item in enumerate(N):
    Err_Lax_sum.append(sum(Err_Lax[i])/item)
    Err_sw_sum.append(sum(Err_sw[i])/item)
plt.close("all")
plt.figure(figsize=(12.0,10.0))
plt.plot(N,Err_sw_sum,label='Steger Warming')
plt.plot(N,Err_Lax_sum,label='Lax-Friedrich')
plt.title("Error per unit length using different schemes")
plt.xlabel('Number of grid points')
plt.ylabel('Error')
plt.legend(shadow=True)
plt.savefig('Err_per_length.png',bbox_inches='tight')
