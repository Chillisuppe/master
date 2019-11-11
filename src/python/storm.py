"""

  Plotting routines for Holland's storm model
  Nicole Beisiegel, August 2012

"""
import matplotlib.pyplot as plt
import numpy as np
import math
import os

SFDIR = os.getcwd()

#Determine parameters
rho= 1.15
e = math.exp(1)
pc = np.array([950])
pn = np.array([1005])
A  = np.array([30.05])
B  = np.array([0.25, 0.75,1,1.25,2.25])

#Parameters for plotting
N=100

def profile(r,pc,pn,A,B):
    erg = pc+(pn-pc)*math.exp(-(A/math.pow(r,B)))
    return erg

def normalize(p,pc,pn):
    erg = (p-pc)/(pn-pc)
    return erg

def wind(A,B,pc, pn, r, rho):
    erg=math.sqrt(A*B*(pn-pc)*math.exp(-(A/math.pow(r,B)))/(rho*math.pow(r,B)))
    return erg

#The ratio of radius_of_max_pressure/radius_of_max_winds approaches 1 for increasing B

def radius_of_max_winds(A,B):
    erg=math.pow(A,1/B)
    return erg

def radius_of_max_pressure(A,B):
    c=(A*B)/(B+1)
    erg = math.pow(c,1/B)
    return erg

def max_wind_speed(B,pn,pc,rho,e):
    C = math.sqrt(B/(rho*e))
    erg = C*math.sqrt(pn-pc)
    return erg

#Assemble vector of radii
r=np.linspace(1,160,N,endpoint=True)
rr=r
rmax = r.max()
nop  = np.zeros(N)
p    = np.zeros(N)
vc   = np.zeros(N)

plt.ioff()
fig = plt.figure(1)

ax1 = fig.add_subplot(111) # Aufteilung Gebiet 2x2, 1 Quadrant
plt.xlim(0,rmax)
plt.xlabel("radial distance from centre [r]")
plt.ylabel("normalized pressure")
plt.title("Profile of a hurricane")

fig=plt.figure(2)

ax2 = fig.add_subplot(111)
plt.xlim(0,rmax)
plt.xlabel("radial distance from centre [r]")
plt.ylabel("sea level pressure [mbar]")
plt.title("")

fig=plt.figure(3)

ax3 = fig.add_subplot(111)
plt.xlim(0,rmax)
plt.xlabel("radial distance from centre [r]")
plt.ylabel("wind speed [m$s^{-1}$")
plt.title("")

M=pc.size
sA=A.size
sB=B.size


for j in range(M):
    for k in range(sA):
        for l in range(sB):
            for i in range(N):
#Determine pressure at radius
                r[i]=r[i]*10
                p[i]  = profile(r[i],pc[j],pn[j],A[k],B[l])
                nop[i] = normalize(p[i],pc[j],pn[j])
                vc[i] = wind(A[k],B[l],pc[j], pn[j], r[i], rho)
                r[i]=r[i]/10
            ax1.plot(r,nop, label="A="+str(A[k])+"; B="+str(B[l]))
            ax2.plot(r,p, label="A="+str(A[k])+"; B="+str(B[l]))
            ax3.plot(r,vc, label="A="+str(A[k])+"; B="+str(B[l]))
ax1.legend()
ax2.legend()
#ax1.legend(loc='upper center', bbox_to_anchor=(1.15,0.55))

plt.show()





