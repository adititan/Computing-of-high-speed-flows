import numpy as np
import matplotlib.pyplot as plt

plt.figure()
niter1,residue1,junk1,junk2 = np.loadtxt("200x50/resid.dat", unpack=True, skiprows=1)
plt.semilogy(niter1,residue1,label = "200x50")

niter4,residue4,junk1,junk2 = np.loadtxt("resid_0.1_0.7.dat", unpack=True, skiprows=1)
plt.plot(niter4,residue4,label = "200x50, xmax=0.1, ymax=0.7 ")

niter2,residue2,junk1,junk2 = np.loadtxt("200x100/resid.dat", unpack=True, skiprows=1)
plt.semilogy(niter2,residue2,label = "200x100")

niter3,residue3,junk1,junk2 = np.loadtxt("05x04/resid.dat", unpack=True, skiprows=1)
plt.semilogy(niter3,residue3,label = "200x50, xmax=0.5, ymax=0.4")

plt.xlabel("Number of Iterations")
plt.ylabel("Residue")
plt.legend(loc="best")
plt.savefig('residue.png')

plt.figure()
niter1,residue1,junk1,junk2 = np.loadtxt("001/resid.dat", unpack=True, skiprows=1)
plt.semilogy(niter1,residue1,label = "CFL = 0.01")


niter3,residue3,junk1,junk2 = np.loadtxt("005/resid.dat", unpack=True, skiprows=1)
plt.semilogy(niter3,residue3,label = "CFL = 0.05")

niter4,residue4,junk1,junk2 = np.loadtxt("01/resid.dat", unpack=True, skiprows=1)
plt.plot(niter4,residue4,label = "CFL = 0.1 ")

niter2,residue2,junk1,junk2 = np.loadtxt("02/resid.dat", unpack=True, skiprows=1)
plt.semilogy(niter2,residue2,label = "CFl = 0.2")


niter5,residue5,junk1,junk2 = np.loadtxt("05/resid.dat", unpack=True, skiprows=1)
plt.semilogy(niter5,residue5,label = "CFL = 0.5")

plt.xlabel("Number of Iterations")
plt.ylabel("Residue")
plt.legend(loc="best")
plt.savefig('residue1.png')


plt.figure()
niter6,residue6,junk1,junk2 = np.loadtxt("best_solution/resid.dat", unpack=True, skiprows=1)
plt.semilogy(niter6,residue6,label = "CFl = Optimal")

plt.xlabel("Number of Iterations")
plt.ylabel("Residue")
plt.legend(loc="best")
plt.savefig('best.png')
