import ConfigParser
config = ConfigParser.RawConfigParser()
config.read('pylocal.cfg')
rootdir = config.get('DEFAULT', 'rootdir')
datdir  = config.get('DEFAULT', 'datdir')
figdir  = config.get('DEFAULT', 'figdir')

import os
import site
site.addsitedir(os.path.join(rootdir,'src/python'))

import numpy as np
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from matplotlib import cm

import pyugrid
import SF2D_ugrid_utils as sfug

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

testsuite = 'testsuite20190826-120210/vortexstrong'

ugrid0 = pyugrid.UGrid.from_ncfile(os.path.join(rootdir,datdir,testsuite,'StormFlash2d00001000.nc'), load_data=True)
ugrid1 = pyugrid.UGrid.from_ncfile(os.path.join(rootdir,datdir,testsuite,'StormFlash2d00003000.nc'), load_data=True)

#vnames = ['bathy', 'depth']
vnames = ['depth', 'm_x']

ugridnew = sfug.getcommonmesh_remesh(ugrid0, ugrid1, vnames)

triang0 = tri.Triangulation(ugrid0.nodes[:, 0], ugrid0.nodes[:, 1], triangles=ugrid0.faces)
triang1 = tri.Triangulation(ugrid1.nodes[:, 0], ugrid1.nodes[:, 1], triangles=ugrid1.faces)
triangnew = tri.Triangulation(ugridnew.nodes[:, 0], ugridnew.nodes[:, 1], triangles=ugridnew.faces)

d0    = ugrid0.data['depth'].data[0,:]
d1    = ugrid1.data['depth'].data[0,:]
d0new = ugridnew.data['depth_0'].data[0,:]
d1new = ugridnew.data['depth_1'].data[0,:]

levels = np.arange(0.97, 1.001, 0.001)
cmap   = cm.get_cmap(name='jet', lut=None)

fig = plt.figure(1)
plt.clf()

plt.subplot(221)
plt.gca().set_aspect('equal')
cs   = plt.tricontourf(triang0, d0, levels=levels, cmap=cmap, extend='both')
plt.triplot(triang0, lw=0.2, color='black')
cbar = plt.colorbar(cs, shrink=0.65, extend='both')
plt.title('1st variable on original 1st mesh')
plt.xlabel('$x$')
plt.ylabel('$y$')

plt.subplot(222)
plt.gca().set_aspect('equal')
cs   = plt.tricontourf(triang1, d1, levels=levels, cmap=cmap, extend='both')
plt.triplot(triang1, lw=0.2, color='black')
cbar = plt.colorbar(cs, shrink=0.65, extend='both')
plt.title('2nd variable on original 2nd mesh')
plt.xlabel('$x$')
plt.ylabel('$y$')

plt.subplot(223)
plt.gca().set_aspect('equal')
cs   = plt.tricontourf(triangnew, d0new, levels=levels, cmap=cmap, extend='both')
plt.triplot(triangnew, lw=0.2, color='black')
cbar = plt.colorbar(cs, shrink=0.65, extend='both')
plt.title('1st variable on refined common mesh')
plt.xlabel('$x$')
plt.ylabel('$y$')

plt.subplot(224)
plt.gca().set_aspect('equal')
cs   = plt.tricontourf(triangnew, d1new, levels=levels, cmap=cmap, extend='both')
plt.triplot(triangnew, lw=0.2, color='black')
cbar = plt.colorbar(cs, shrink=0.65, extend='both')
plt.title('2nd variable on refined common mesh')
plt.xlabel('$x$')
plt.ylabel('$y$')


plt.show()
