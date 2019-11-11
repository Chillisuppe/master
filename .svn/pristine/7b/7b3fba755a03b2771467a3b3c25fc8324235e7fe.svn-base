#-------------------------------------------------------------------------------
# This is a script to compare two runs of the testsuite.
#
# Usage: has to be called as a script with two or three arguments.
# The first two arguments have to be names of testsuite directories
# and need to be present in the compile directory.
# The third argument is an optional tolerance (default 1E-12).

#-------------------------------------------------------------------------------
# import modules and define several utilities

import ConfigParser
config = ConfigParser.RawConfigParser()
config.read('pylocal.cfg')
rootdir = config.get('DEFAULT', 'rootdir')
datdir  = config.get('DEFAULT', 'datdir')

import sys
import os
import site
site.addsitedir(os.path.join(rootdir,'src/python'))

import numpy as np
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import MaxNLocator

import pyugrid
import SF2D_ugrid_utils as sfug

#-------------------------------------------------------------------------------

def error(testsuite1, testsuite2, tol, case):

  # get type of equations from case string
  if (case[:7] == 'eulerPT'):
    # Euler equations with pot. temperature
    unknowns = ['density', 'm_x', 'm_y', 'pot_temp']
  elif (case[:7] == 'eulerTE'):
    # Euler equations with total energy
    unknowns = ['density', 'm_x', 'm_y', 'energy']
  elif (case[:6] == 'linadv'):
    # linear advection equation
    unknowns = ['tracer']
  else:
    # shallow water equations
    unknowns =  ['depth', 'm_x', 'm_y']

  # define files for comparison (should be last netcdf file)
  ncfile1 = [fn for fn in sorted(os.listdir(os.path.join(rootdir,datdir,testsuite1,case)))
             if fn.endswith('0.nc')][-1]
  ncfile2 = [fn for fn in sorted(os.listdir(os.path.join(rootdir,datdir,testsuite2,case)))
             if fn.endswith('0.nc')][-1]

  ugridA = pyugrid.UGrid.from_ncfile(os.path.join(rootdir,datdir,testsuite1,case,ncfile1), load_data=True)
  ugridB = pyugrid.UGrid.from_ncfile(os.path.join(rootdir,datdir,testsuite2,case,ncfile2), load_data=True)

  txtfile.write('\n'+case+':\n')
  print '{0:s}:'.format(case)

  # Compare the names of the NetCDF files
  if (ncfile1 != ncfile2):
    txtfile.write('The test failed. Final NetCDF files have different names.\n')
    print '  Test failed. The NetCDF files have different names.'
    return 1

  # Compare model times
  if (ugridA.data[unknowns[0]].time[0] != ugridB.data[unknowns[0]].time[0]):
    txtfile.write('The test failed. Model times for the final output are different.\n')
    print '  Test failed. The model times for the final output are different.'
    return 1

  # Compare number of cells
  nodesA = len(ugridA.nodes[:,0])
  nodesB = len(ugridB.nodes[:,0])
  if (nodesA != nodesB):
    # plot data and report failure
    triangA = tri.Triangulation(ugridA.nodes[:, 0], ugridA.nodes[:, 1], triangles=ugridA.faces)
    triangB = tri.Triangulation(ugridB.nodes[:, 0], ugridB.nodes[:, 1], triangles=ugridB.faces)
    fig = plt.figure()
    iunknown = 0
    for unknown in unknowns:
      dA  = ugridA.data[unknown]
      dB  = ugridB.data[unknown]
      if (np.any(np.isnan(dB))):
        levels = MaxNLocator(nbins=15).tick_values(np.min(dA.data[0,:]), np.max(dA.data[0,:]))
        txtfile.write('NaNs found in {0:s} field!\n'.format(unknown))
        print '  NaNs found in {0:s} field!'.format(unknown)
      else:
        levels = MaxNLocator(nbins=15).tick_values(np.minimum(np.min(dA.data[0,:]), np.min(dB.data[0,:])),
                                                   np.maximum(np.max(dA.data[0,:]), np.max(dB.data[0,:])))
      plt.subplot(len(unknowns),2,iunknown*2+1)
      plt.gca().set_aspect('equal')
      cmap = cm.get_cmap(name='jet', lut=None)
      cs   = plt.tricontourf(triangA, dA.data[0,:], levels=levels, cmap=cmap, extend='both')
      cbar = plt.colorbar(cs, shrink=0.65, extend='both')
      plt.xlabel('$x$')
      plt.ylabel('$y$')
      plt.title(unknown+' '+testsuite1)

      plt.subplot(len(unknowns),2,iunknown*2+2)
      plt.gca().set_aspect('equal')
      cmap = cm.get_cmap(name='jet', lut=None)
      cs   = plt.tricontourf(triangB, dB.data[0,:], levels=levels, cmap=cmap, extend='both')
      cbar = plt.colorbar(cs, shrink=0.65, extend='both')
      plt.xlabel('$x$')
      plt.ylabel('$y$')
      plt.title(unknown+' '+testsuite2)

      iunknown = iunknown+1

    plt.suptitle('Testcase: {0:s}, $t={1:5.2f}$ (different number of nodes!)'.format(case, ugridA.data[unknowns[0]].time[0]))
    fig.set_size_inches(3*fig.get_size_inches())
    fig.savefig(os.path.join(rootdir,datdir,case+'-error.pdf'))
    plt.close()

    txtfile.write('The test failed. Meshes have different number of nodes. Data has been plotted.\n')
    print '  Test failed. Meshes have different number of nodes. Data has been plotted.'

    return 1

  qug = sfug.QuadUGRID(1, 2)
  sqr = lambda x: x**2

  # compute area to get normalized L2 error
  faces = ugridA.faces
  x = ugridA.nodes[:,0]
  y = ugridA.nodes[:,1]

  area = 0.0
  ctri = np.zeros((3,2))
  for i in range(len(faces[:,0])):
    ctri[:,0] = x[faces[i,:]]
    ctri[:,1] = y[faces[i,:]]
    area     += 0.5*np.abs((ctri[0,0]-ctri[2,0])*(ctri[1,1]-ctri[0,1]) -
                        (ctri[0,0]-ctri[1,0])*(ctri[2,1]-ctri[0,1]))

  L2err = {}
  Lierr = {}
  for unknown in unknowns:
    # compute error functions
    dA  = ugridA.data[unknown]
    dB  = ugridB.data[unknown]
    err = pyugrid.UVar(unicode('err_'+unknown), location='node',
                       data=dA.data-dB.data, time=0.0)
    ugridA.add_data(err)

    # calculate norms of errors
    L2err[unknown] = np.sqrt(qug(ugridA, 'err_'+unknown, sqr))/area
    Lierr[unknown] = np.max(np.abs(ugridA.data['err_'+unknown].data[0]))

    txtfile.write('L2   error {0:s}: {1:+6.4e}\n'.format(unknown, L2err[unknown]))
    txtfile.write('Linf error {0:s}: {1:+6.4e}\n'.format(unknown, Lierr[unknown]))


  if (max(L2err.values()) < tol and max(Lierr.values()) < tol):
    txtfile.write('The test passed.\n')
    print '  Test passed.'
    return 0
  else:
    # plot data and report failure
    triang = tri.Triangulation(ugridA.nodes[:, 0], ugridA.nodes[:, 1], triangles=ugridA.faces)
    fig = plt.figure()
    iunknown = 0
    for unknown in unknowns:
      dA  = ugridA.data[unknown]
      dB  = ugridB.data[unknown]
      if (np.any(np.isnan(dB))):
        levels = MaxNLocator(nbins=15).tick_values(np.min(dA.data[0,:]), np.max(dA.data[0,:]))
        txtfile.write('NaNs found in {0:s} field!\n'.format(unknown))
        print '  NaNs found in {0:s} field!'.format(unknown)
      else:
        levels = MaxNLocator(nbins=15).tick_values(np.minimum(np.min(dA.data[0,:]), np.min(dB.data[0,:])),
                                                   np.maximum(np.max(dA.data[0,:]), np.max(dB.data[0,:])))
      plt.subplot(len(unknowns),3,iunknown*3+1)
      plt.gca().set_aspect('equal')
      cmap = cm.get_cmap(name='jet', lut=None)
      cs   = plt.tricontourf(triang, dA.data[0,:], levels=levels, cmap=cmap, extend='both')
      cbar = plt.colorbar(cs, shrink=0.65, extend='both')
      plt.xlabel('$x$')
      plt.ylabel('$y$')
      plt.title(unknown+' '+testsuite1)

      plt.subplot(len(unknowns),3,iunknown*3+2)
      plt.gca().set_aspect('equal')
      cmap = cm.get_cmap(name='jet', lut=None)
      cs   = plt.tricontourf(triang, dB.data[0,:], levels=levels, cmap=cmap, extend='both')
      cbar = plt.colorbar(cs, shrink=0.65, extend='both')
      plt.xlabel('$x$')
      plt.ylabel('$y$')
      plt.title(unknown+' '+testsuite2)

      plt.subplot(len(unknowns),3,iunknown*3+3)
      plt.gca().set_aspect('equal')
      cmap = cm.get_cmap(name='jet', lut=None)
      cs   = plt.tricontourf(triang, dA.data[0,:]-dB.data[0,:], cmap=cmap, extend='both')
      cbar = plt.colorbar(cs, shrink=0.65, extend='both')
      plt.xlabel('$x$')
      plt.ylabel('$y$')
      plt.title(unknown+' difference')

      iunknown = iunknown+1

    plt.suptitle('Testcase: {0:s}, $t={1:5.2f}$'.format(case, ugridA.data[unknowns[0]].time[0]))
    fig.set_size_inches(3*fig.get_size_inches())
    fig.savefig(os.path.join(rootdir,datdir,case+'-error.pdf'))
    plt.close()

    txtfile.write('The test failed. Data has been plotted.\n')
    print '  Test failed. Data has been plotted.'

    return 1

#-------------------------------------------------------------------------------
# this is the main routine
def main(testsuite1, testsuite2, tol):

  global txtfile
  txtfile = open(os.path.join(rootdir,datdir,'testsuite_comparison.txt'),'w')
  txtfile.write('Testsuite runs {0:s} and {1:s} are compared.\n'.format(testsuite1, testsuite2))
  txtfile.write('The tolerance is {0:6.4e}\n'.format(tol))

  tests1 = set(os.listdir(os.path.join(rootdir,datdir,testsuite1)))
  tests2 = set(os.listdir(os.path.join(rootdir,datdir,testsuite2)))

  difft1t2 = ', '.join(str(tst) for tst in tests1.difference(tests2))
  difft2t1 = ', '.join(str(tst) for tst in tests2.difference(tests1))
  txtfile.write('Tests only present in {0:s}: '.format(testsuite1) + difft1t2 + '\n')
  txtfile.write('Tests only present in {0:s}: '.format(testsuite2) + difft2t1 + '\n')

  print 'Tests only present in {0:s}: '.format(testsuite1), difft1t2
  print 'Tests only present in {0:s}: '.format(testsuite2), difft2t1
  print

  nfailed = 0
  for test in tests1.intersection(tests2):
    nfailed = nfailed + error(testsuite1, testsuite2, tol, test)

  txtfile.write('\nNumber of tests which failed: {0:2d}\n'.format(nfailed))
  txtfile.close()

  print 'Number of tests which failed: {0:2d}\n'.format(nfailed)
  return

#-------------------------------------------------------------------------------
# run script only when it is actually called as a script
if __name__ == '__main__':
# check whether console input is valid and set default if not
  if len(sys.argv) < 3:
    raise ValueError('Input has to include two testsuite directory names.')

# read testsuite directory names from console input
  testsuite1 = sys.argv[1]
  testsuite2 = sys.argv[2]
  tol = 1E-12
  if len(sys.argv) == 4:
    tol = np.float(sys.argv[3])

# run main routine
  main(testsuite1, testsuite2, tol)
