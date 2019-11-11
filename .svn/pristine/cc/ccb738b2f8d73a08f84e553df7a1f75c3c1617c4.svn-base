"""

   Script for comparing basically 1D problems using vtu file output
   Nicole Beisiegel
   May, 2013

"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math as m

# a) Read through every line and search for a) coordinates and dof numbers
# b) find ssh value acc to dof number
# c) interpolate for the adaptive grid

#----------------------------------------------------------------------------
"""
  FUNCTION reading 
  INPUT  file (in vtu format from from StormFlash2d)
  OUTPUT numdofs - number of degrees of freedom (dofs)
         coo     - coordinates and fluid height of dofs
"""

def reading(file):
    fo=open(file,"r")
    contents=fo.readlines()
    idx=0
    coo=np.zeros(0)
    for line in open(file):
        li = line.strip()
        idx=idx+1
        if li.startswith("<"):
            if (li.find("Piece NumberOfPoints",0,30) !=-1):
                numdofs= int(str(li[24:31]).strip())
            if idx==5:
                for i in range(idx,numdofs+idx):
                    tmp=np.array((str(contents[i]).strip()).split('      '))
                    tmp2 = tmp.astype(np.float)
                    if (i == idx):
                        coo = tmp2
                    else:
                        coo=np.vstack([coo,tmp2])
                break
    return numdofs, coo

#-----------------------------------------------------------------------------    
"""
  FUNCTION reading_lineplot
  INPUT   file (.txt file --> VALIDATION SF2D)
  
"""

def reading_lineplot(file):
    fo=open(file,"r")
    cont=fo.readlines()
    numdofs=int(cont[0])
##    idx=[]
    ssh=np.zeros(numdofs)
    bat=np.zeros(numdofs)
    coo=np.zeros(numdofs)
    cooy=np.zeros(numdofs)
    for i in range(1,numdofs):
        ssh[i]=cont[i+1]
        bat[i]=cont[numdofs+i+1]
        coo[i]=cont[numdofs*2+i+1]
        cooy[i]=cont[numdofs*3+i+1]
    order = coo.argsort()
    coo = np.take(coo, order, 0)
    ssh = np.take(ssh, order,0)
    bat = np.take(bat, order,0)
    cooy = np.take(cooy, order, 0)
    idx = getlinedofs(numdofs,cooy)
 ##   for i in range(1,numdofs):
 ##       if (float(cooy[i]) == 0.0):
 ##           idx.append(i)
    #now sort the entries
    return ssh, bat, coo, idx

#---------------------------------------------------------------
"""
   FUNCTION getlinedofs 
   Computes the indices of dofs that lie on the y-axis (y=0)
"""

def getlinedofs(numdofs, coo):
    idx=[]
    for i in range(numdofs):
        if (coo[i,1]== 0.0):
            idx.append(i)
    return idx

#----------------------------------------------------------------------
"""
FUNCTIONS bathymetry*
for computing bathymetry for Balzano tests analytically
"""

def bathymetry1(alpha,x):
    bathy=alpha*x
    return bathy

def bathymetry2(alpha,x,c1,c2):
    if ((x <= float(3600.0)) | (x >= float(6000.0))):
        bathy= x*alpha
    elif (x >= float(3600.0)) & (x <= float(4800.0)):
        bathy= -x*alpha+c1
    else:
        bathy= float(x/920.0) -c2
    return bathy

#----------------------------------------------------------------------
"""
FUNCTION getbathy
Determines the bathymetry for all coordinates coo
"""

def getbathy(coo):
    m=np.size(coo)
    bat=np.zeros(m)
    alpha=float(5.0/13800.0)
    c1 = float(60.0/23.0)
    c2 = float(100.0/23.0)
    for i in range(m):
#        bat[i]=bathymetry1(alpha,coo[i])
       bat[i]=bathymetry2(alpha,coo[i],c1,c2)
    return bat

#-----------------------------------------------------------------------    
"""
FUNCTION getheight
Function that adds fluid height and bathymetry in the numpy-sense
"""

def getheight(ssh, bathy):
    height=np.add(ssh, bathy)
    return height

#------------------------------------------------------------------------
"""
FUNCTION interpolate
INPUT  coo_fine - points of fine uniform mesh
coo_adap - points of adaptive mesh
ssh_adap - fluid height of adaptive numerical solution

OUTPUT ssh_new  - interpolated fluid height from adaptive numerical 
solution onto points from fine uniform mesh
deltax   - grid spacing of uniform fine mesh

"""

def interpolate(coo_fine, coo_adap, ssh_adap):
    m = np.size(coo_fine)
    ssh_new=np.zeros(m)
    deltax= coo_fine[2]-coo_fine[1]
    j=0
    k=0
#just works for sfc order and linear basis functions
    for i in range(m):
        if (coo_fine[i] == coo_adap[j]):
            ssh_new[i] = ssh_adap[j]
            j=j+1
            if (k >=1):
                a = ((ssh_new[i]-ssh_new[i-k-1])/(coo_fine[i]-coo_fine[i-k-1]))
                c = ssh_new[i] - a*coo_fine[i]
                for l in range(0,k):
                    ssh_new[i-k+l]=a*coo_fine[i-k+l]+c
            k=0                    
        else:
            k=k+1
            ssh_new[i] = 0.0
    return ssh_new, deltax
#--------------------------------------------------------------------------    
"""
    FUNCTION sort

    INPUT  coo - X-coordinates of mesh points
    OUTPUT idx - Indices of mesh points after elimination of multiple points 
                 with the same x-coordinates
"""

def sort(coo):
    m=np.size(coo)
    idx=[]
    idx.append(0)
    for i in range(1,m):
        if coo[i] != coo[i-1]:
            idx.append(i)
    #print 'reduce', m, 'to', np.size(idx)
    return idx
#----------------------------------------------------------------------------

# Plotting script

# Determine indices of vtu files to be read in
indices=range(0,140000,5000)                 

#Open file for writing integrals
f=open('workfile', 'w')         
    
# For every vtu file..
for indi in indices:
#... determine path of files for fine and adaptive simulation
    file_fine= "/home/nicole/Beachtwo/StormFlash2d00"+str(indi).zfill(6)+".vtu"
    file_adap="/home/nicole/Beachtwo_84/StormFlash2d00"+str(indi).zfill(6)+".vtu"

# Determine path of lineplot file (--> Run SF2D with VALIDATION)
    file_linp="/home/nicole/lineplot 5.txt"

# Determine coordinates, fluid height and number of dofs for fine and adaptive simulation
    (numdofs_fine, coo_fine) = reading(file_fine)
    (numdofs_adap, coo_adap) = reading(file_adap)

#    (ssh_lin, bat_lin, coo_lin, idx_lin) = reading_lineplot(file_linp)
#    siz= np.size(bat_lin)
#    hei_lin=np.zeros(siz)
#    for i in range(siz):
#        hei_lin[i] = ssh_lin[i]+ bat_lin[i]

# order coordinates and fluid height according to x-coordinates
    order = coo_fine[:, 0].argsort()
    coo_fine  = np.take(coo_fine, order, 0)
    order2=coo_adap[:,0].argsort()
    coo_adap = np.take(coo_adap, order2,0)

# Get indices of dofs that lie on x-axis
    idx_fine = getlinedofs(numdofs_fine, coo_fine)
    idx_adap = getlinedofs(numdofs_adap, coo_adap)

# Interpolate fluid height of adaptive simulation onto dofs of uniform run
    (ssh_new, deltax)=interpolate(coo_fine[idx_fine,0], coo_adap[idx_adap,0], coo_adap[idx_adap,2])

# Get bathymetry on adaptive and uniform mesh
    bat_fine = getbathy(coo_fine[idx_fine,0])
    bat_adap = getbathy(coo_adap[idx_adap,0])

# Determine total height = bathy + ssh for fine, interpolated and adaptive mesh
    hei_fine = getheight(coo_fine[idx_fine,2],bat_fine)
    hei_new  = getheight(ssh_new, bat_fine)
    hei_adap = getheight(coo_adap[idx_adap,2], bat_adap)

# Determine min and max of xaxis
    xmax=max(coo_fine[idx_fine,0])
    xmin=min(coo_fine[idx_fine,0])

    plt.ioff()

# Set parameters for plotting
    mpl.rcParams.update({'font.size': 18})
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 22}
#    plt.rc('font', **font)

# Initialize first figure: 
# Plot of Uniform vs Adaptive Simulation (fluid height)
    fig=plt.figure(1)

# USER: Please define x limits
    ax=fig.add_subplot(111)
#    plt.xlim(xmin+(xmax-xmin)/2,xmax)
    plt.xlim(2000,8000)
#    plt.xlim(xmin,xmax)
    plt.xlabel("x")
    plt.ylabel("$\phi$")
    plt.title("Fluid Height: Uniform Vs. Adaptive Simulation")
    
    ax.plot(coo_fine[idx_fine,0], coo_fine[idx_fine,2], 'k')
    ax.plot(coo_adap[idx_adap,0], coo_adap[idx_adap,2],'k*')
    ax.plot(coo_fine[idx_fine,0], ssh_new,'k-')
    # ax.plot(coo_fine[idx_fine,0], hei_fine,'k')
    # ax.plot(coo_fine[idx_fine,0], bat_fine,'k', linewidth=3)
    # ax.plot(coo_lin[idx_lin], hei_lin[idx_lin],'k', linewidth=3)
    # ax.plot(coo_lin[idx_lin], bat_lin[idx_lin],'k', linewidth=3)
    
# Initialize second figure:
# Plot of Uniform vs Adaptive Simulation (bathymetry & height = b+ssh)
    fig=plt.figure(2)

    ax=fig.add_subplot(111)
  #  plt.xlim(xmin+(xmax-xmin)/2,xmax)
    plt.xlim(xmin, xmax)
    plt.xlabel("x")
    plt.ylabel("$\phi$")
    plt.title("Plot Uniform Vs. Adaptive Simulation")
    
    ax.plot(coo_fine[idx_fine,0], bat_fine, 'k', linewidth=3)
    ax.plot(coo_fine[idx_fine,0], hei_fine, 'k-')
    ax.plot(coo_adap[idx_adap,0], hei_adap, 'k')


# Check for double entries
    iin = sort(coo_fine[idx_fine,0])
# Fluid height on fine mesh
    coof = coo_fine[idx_fine,2]

# Now compute the integrals
    domainwidth = float(2760)
    vec=np.zeros((np.size(idx_fine)))

# Difference fine - crs(interpolated) fluid height
    ssh_diff=np.subtract(coof, ssh_new)

    sum1 =np.sum(np.square(coof))
    sum2 =np.sum(np.square(ssh_diff))

# print np.size(iin)
    ssh2 = np.subtract(coof[iin], ssh_new[iin])
    ssh3 = coof[iin]
    
# Compute mean values of fluid height for every interval [x_i - x_{i-1}]
    kk=np.size(ssh2)-1
    meanssh2=np.zeros(kk)
    meanssh3=np.zeros(kk)
    for i in range(kk):
        meanssh2[i]=float(0.5)*(m.pow(ssh2[i],2)+m.pow(ssh2[i+1],2))
        meanssh3[i]=float(0.5)*(m.pow(ssh3[i],2)+m.pow(ssh3[i+1],2))

# L2 Norm of uniform fine and diff solution
    int_fine = sum1 #*deltax*domainwidth
    int_diff = sum2 

# Max Norm of uniform fine solution and the difference between numerical solutions
    max_fine = max(np.absolute(coof))
    max_diff = max(np.absolute(ssh_diff))

#L2 Norm of quadratic ..
    int_fine2 = np.sum(meanssh3)
    int_diff2 = np.sum(meanssh2)

# Max Norm of quadratic ..
    max_fine2 = max(np.absolute(meanssh3))
    max_diff2 = max(np.absolute(meanssh2))
    
#Compute relative errors for the L2 norm
    error = int_diff/int_fine
    error2 = int_diff2/int_fine2

#Compute relative error in max Norm
    errormax = max_diff/max_fine
    errormax2 = max_diff2/max_fine2

    #for wetting and drying
    if (int_fine <= 10E-10):
        error  = int_diff
    if (int_fine2 <= 10E-10):
        error2 = int_diff2
    if (max_fine <=10E-10):
        errormax = max_diff
    if (max_fine2 <= 10E-10):
        errormax2 = max_diff2

    f.write('Timestep No.'+str(indi)+'\n')
    f.write('$L_2$ Norm of adaptive and fine solution'+str(int_diff2)+','+str(int_fine2)+'\n')
    f.write('Relative Errors in $L_2$-Norm:' + str(error)+','+str(error2)+'\n')
    f.write('Relative Errors in Max Norm:' +str(errormax)+','+str(errormax2)+'\n')

f.close()


# plot the timeline of plots
indic=range(0,140000,5000)

# Third figure: Plot of the uniform solution
fig=plt.figure(3)
#plt.xlim(xmin+3*(xmax-xmin)/4,xmax)
#plt.ylim(2,5.5)
#plt.xlim(xmin, xmax)
plt.xlim(2000,8000)
plt.xlabel("x")
plt.ylabel("$\phi$")
#ax.yaxis.label.set_size(40)
mpl.rcParams.update({'font.size': 18})
#mpl.axes.titlesize
#mpl.axes.labelsize

plt.title("Plot Uniform Simulation")
for indx in indic:
#    file="/home/nicole/Beachtwo/StormFlash2d00"+str(indx).zfill(6)+".vtu"
    file="/home/nicole/Beachtwo/StormFlash2d00"+str(indx).zfill(6)+".vtu"
    (numdofs, coo) = reading(file)
    order = coo[:, 0].argsort()
    coo   = np.take(coo, order,0)
    idx = getlinedofs(numdofs, coo)   
    bat = getbathy(coo[idx,0])
    hei = getheight(coo[idx,2],bat)
    plt.plot(coo[idx,0],bat,'k',linewidth=3)
    plt.plot(coo[idx,0],hei,'k',linewidth=1)

plt.show()
