#*******************************************************************************
# This is a script to create plots for the testsuite documentation.
#
# Usage: has to be called as a script with one argument.
# The argument has to be a name of testsuite directory and needs to be present
# in the compile directory.
# This script creates the documentation plots for all testsuite cases.

#*******************************************************************************
# import modules and define several utilities

import ConfigParser
config = ConfigParser.RawConfigParser()
config.read('pylocal.cfg')
rootdir = config.get('DEFAULT', 'rootdir')
datdir  = config.get('DEFAULT', 'datdir')
figdir  = config.get('DEFAULT', 'figdir')

import sys
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

prop_cycle = plt.rcParams['axes.prop_cycle']
col = prop_cycle.by_key()['color']

#*******************************************************************************
def shockplot(testsuite, case):
# check if directory exists
  (destination) = os.path.isdir(os.path.join(rootdir,datdir,testsuite,case))
  if destination == False:
    print 'Directory for ' + case + ' does not exist. Skipping this testcase.'
    return 1
  else:
    print 'Plotting testcase {0:s}'.format(case)

  # define files and saves
  if case == 'shockquad':
    files = ['StormFlash2d00000000.nc', 'StormFlash2d00000200.nc', 'StormFlash2d00000400.nc']
  else:
    files = ['StormFlash2d00000000.nc', 'StormFlash2d00000050.nc', 'StormFlash2d00000100.nc']
  saves = [f.replace('StormFlash2d',case) for f in files]
  saves = [f.replace('.nc','.pdf') for f in saves]

  # load data
  for cnt in xrange(len(files)):
    ugrid = pyugrid.UGrid.from_ncfile(os.path.join(rootdir,datdir,testsuite,case,files[cnt]), load_data=True)

    # cross section along y = 1.0
    vnames = ['bathy', 'depth', 'm_x', 'm_y', 'exact_depth', 'exact_m_x', 'exact_m_y']

    xv, yv, varpts = sfug.get_crosssection([0,1], [9.99999,1], ugrid, vnames)

    fig = plt.figure()
    # plot fluid depth
    plt.subplot(211)
    plt.plot(xv, varpts['bathy']+varpts['exact_depth'], linewidth=2.0, label='exact')
    plt.plot(xv, varpts['bathy']+varpts['depth']      , linewidth=1.0, label='numerical')
    plt.xlim([0, 10])
    if case == 'shockdry':
      plt.ylim([-0.05, 0.3])
    else:
      plt.ylim([0.7, 3.3])
    plt.grid(True)
    plt.ylabel('$h+b$')
    plt.legend(loc=1)

    # x-momentum
    plt.subplot(212)
    plt.plot(xv, varpts['exact_m_x'], linewidth=2.0, label='exact')
    plt.plot(xv, varpts['m_x']      , linewidth=1.0, label='numerical')
    plt.xlim([0, 10])
    if case == 'shockdry':
      plt.ylim([-0.01, 0.13])
    else:
      plt.ylim([-0.2, 4.5])
    plt.grid(True)
    plt.xlabel('$x$')
    plt.ylabel('$hu$')
    #plt.legend(loc=1)

    #plt.suptitle(case)
    fig.subplots_adjust(left=0.15, right=0.96, bottom=0.1, top=0.97)
    fig.set_size_inches(4,4.2)
    fig.savefig(os.path.join(rootdir,figdir,saves[cnt]))
    plt.close()

  return

#*******************************************************************************
def eulersimplewave(testsuite, case):
# check if directory exists
  if case == 'limited':
    casePT = 'eulerPTsimplewaveshu'
    caseTE = 'eulerTEsimplewaveshu'
    tag = 'shu'
  elif case == 'non-limited':
    casePT = 'eulerPTsimplewave'
    caseTE = 'eulerTEsimplewave'
    tag = ''
  else:
    print 'Wrong case input for simplewave.'
    return 1

  (destination1) = os.path.isdir(os.path.join(rootdir,datdir,testsuite,casePT))
  (destination2) = os.path.isdir(os.path.join(rootdir,datdir,testsuite,caseTE))

  case = 'eulersimplewave'+tag

  if destination1 == False:
    print 'Directory for eulerPTsimplewave' + tag + ' does not exist. Skipping this testcase.'
    return 1
  elif destination2 == False:
    print 'Directory for eulerTEsimplewave' + tag + ' does not exist. Skipping this testcase.'
    return 1
  else:
    print 'Plotting testcase {0:s}'.format(case)

  # define files and saves
  files = ['StormFlash2d00000000.nc', 'StormFlash2d00000400.nc', 'StormFlash2d00001000.nc']
  saves = [f.replace('StormFlash2d',case) for f in files]
  saves = [f.replace('.nc','.pdf') for f in saves]

  # load data
  for cnt in xrange(len(files)):
    ugridPT = pyugrid.UGrid.from_ncfile(os.path.join(rootdir,datdir,testsuite,casePT,files[cnt]), load_data=True)
    ugridTE = pyugrid.UGrid.from_ncfile(os.path.join(rootdir,datdir,testsuite,caseTE,files[cnt]), load_data=True)

    # cross section along y = 1.0
    vnames = ['density', 'm_x']

    xvPT, yvPT, varptsPT = sfug.get_crosssection([0,0.5], [1,0.5], ugridPT, vnames)
    xvTE, yvTE, varptsTE = sfug.get_crosssection([0,0.5], [1,0.5], ugridTE, vnames)

    fig = plt.figure()
    # density
    plt.subplot(211)
    plt.plot(xvTE, varptsTE['density'], linewidth=2.0, color='r', label='TE')
    plt.plot(xvPT, varptsPT['density'], linewidth=1.0, color='b', label='PT',linestyle='-.')
    plt.xlim([0, 1])
    plt.ylim([0.85, 1.15])
    plt.ylabel(r'$\rho$ [kg m$^{-3}$]')
    plt.grid(True)
    plt.legend()

    # x-momentum
    plt.subplot(212)
    plt.plot(xvTE, varptsTE['m_x'], linewidth=2.0, color='r', label='TE')
    plt.plot(xvPT, varptsPT['m_x'], linewidth=1.0, color='b', label='PT',linestyle='-.')
    plt.xlim([0, 1])
    plt.ylim([-75, 75])
    plt.xlabel('$x$')
    plt.ylabel(r'$\rho u$ [kg m$^{-2}$ s$^{-1}$]')
    plt.grid(True)
    plt.legend()

    #plt.suptitle(case)
    fig.subplots_adjust(left=0.15, right=0.96, bottom=0.1, top=0.97)
    fig.set_size_inches(4,4.2)
    fig.savefig(os.path.join(rootdir,figdir,saves[cnt]))
    plt.close()

  return

#*******************************************************************************
def beachplot(testsuite, case):
# check if directory exists
  (destination) = os.path.isdir(os.path.join(rootdir,datdir,testsuite,case))
  if destination == False:
    print 'Directory for ' + case + ' does not exist. Skipping this testcase.'
    return 1
  else:
    print 'Plotting testcase {0:s}'.format(case)

# define files and saves
  files = ['StormFlash2d00000000.nc', 'StormFlash2d00000400.nc',
           'StormFlash2d00000800.nc', 'StormFlash2d00001200.nc',
           'StormFlash2d00001600.nc', 'StormFlash2d00002000.nc']
  saves = [f.replace('StormFlash2d',case) for f in files]
  saves = [f.replace('.nc','.pdf') for f in saves]

# load data
  for cnt in xrange(len(files)):
    ugrid = pyugrid.UGrid.from_ncfile(os.path.join(rootdir,datdir,testsuite,case,files[cnt]), load_data=True)

# cross section along y = 1.0
    vnames = ['bathy', 'depth', 'm_x', 'm_y']
    xv, yv, varpts = sfug.get_crosssection([0,1], [14000,1], ugrid, vnames)

    fig = plt.figure()
    plt.plot(xv, varpts['bathy'],                 linewidth=2.0, label='$b$')
    plt.plot(xv, varpts['bathy']+varpts['depth'], linewidth=1.0, label='$h+b$')
    plt.xlim([0, 14000])
    plt.ylim([0, 5.5])
    plt.xlabel('$x$')
    plt.ylabel('$z$')
    plt.legend(loc=4)

    fig.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.97)
    fig.set_size_inches(4,2.2)
    fig.savefig(os.path.join(rootdir,figdir,saves[cnt]))
    plt.close()

  return

#*******************************************************************************
def basinplot(testsuite, case):
# check if directory exists
  (destination) = os.path.isdir(os.path.join(rootdir,datdir,testsuite,case))
  if destination == False:
    print 'Directory for ' + case + ' does not exist. Skipping this testcase.'
    return 1
  else:
    print 'Plotting testcase {0:s}'.format(case)

  # define files and saves
  files = ['StormFlash2d00000000.nc', 'StormFlash2d00000500.nc', 'StormFlash2d00001000.nc']
  saves = [f.replace('StormFlash2d',case) for f in files]
  saves = [f.replace('.nc','.pdf') for f in saves]

  # load data
  for cnt in xrange(len(files)):
    ugrid = pyugrid.UGrid.from_ncfile(os.path.join(rootdir,datdir,testsuite,case,files[cnt]), load_data=True)

  # cross section along y = 0.0
    vnames = ['bathy', 'depth', 'm_x', 'm_y', 'exact_depth', 'exact_m_x', 'exact_m_y']
    xv, yv, varpts = sfug.get_crosssection([-4000,0], [4000,0], ugrid, vnames)

    fig = plt.figure()
    # fluid depth
    plt.subplot(211)
    plt.plot(xv, varpts['bathy'],                       color='0.6', linestyle='dashed', linewidth=1.0)
    plt.plot(xv, varpts['bathy']+varpts['exact_depth'], color=col[0], linewidth=2.0, label='exact')
    plt.plot(xv, varpts['bathy']+varpts['depth']      , color=col[1], linewidth=1.0, label='numerical')
    plt.xlim([-4000, 4000])
    plt.ylim([0, 2.5])
    #plt.xlabel('$x$')
    plt.ylabel('$z$')
    plt.legend(loc=1)

    # x-momentum
    plt.subplot(212)
    plt.plot(xv, varpts['exact_m_x'], color=col[0], linewidth=2.0, label='exact')
    plt.plot(xv, varpts['m_x']      , color=col[1], linewidth=1.0, label='numerical')
    plt.xlim([-4000, 4000])
    plt.xlabel('$x$')
    plt.ylabel('$hu$')
    #plt.legend(loc=1)

    #plt.suptitle(case)
    fig.subplots_adjust(left=0.15, right=0.96, bottom=0.1, top=0.97)
    fig.set_size_inches(4,4.2)
    fig.savefig(os.path.join(rootdir,figdir,saves[cnt]))
    plt.close()

  return

#*******************************************************************************
def lakeplot(testsuite, case):
# check if directory exists
  (destination) = os.path.isdir(os.path.join(rootdir,datdir,testsuite,case))
  if destination == False:
    print 'Directory for ' + case + ' does not exist. Skipping this testcase.'
    return 1
  else:
    print 'Plotting testcase {0:s}'.format(case)

  # define files and saves
  files = 'StormFlash2d00000000.nc'
  saves = files.replace('StormFlash2d',case)
  saves = saves.replace('.nc','.pdf')

  # load data
  ugrid = pyugrid.UGrid.from_ncfile(os.path.join(rootdir,datdir,testsuite,'wellbalanced',files), load_data=True)

  # cross section along y = 1.0
  vnames = ['bathy', 'depth', 'm_x', 'm_y', 'exact_depth', 'exact_m_x', 'exact_m_y']
  xv, yv, varpts = sfug.get_crosssection([0,15], [30,15], ugrid, vnames)

  fig = plt.figure()
  plt.plot(xv, varpts['bathy'],                 color = 'b', linewidth=2.0, label='bathy')
  plt.plot(xv, varpts['bathy']+varpts['depth'], color = 'r', linewidth=1.0, label='total height')
  plt.xlim([0, 30])
  plt.ylim([0, 12])
  plt.xlabel('$x$')
  plt.ylabel('$h+b$ or $b$')
  plt.legend()

  fig.subplots_adjust(left=0.15, right=0.96, bottom=0.17, top=0.97)
  fig.set_size_inches(4,2.2)
  fig.savefig(os.path.join(rootdir,figdir,saves))
  plt.close()

  return

#*******************************************************************************
def windplot(testsuite, case):
# check if directory exists
  (destination) = os.path.isdir(os.path.join(rootdir,datdir,testsuite,case))
  if destination == False:
    print 'Directory for ' + case + ' does not exist. Skipping this testcase.'
    return 1
  else:
    print 'Plotting testcase {0:s}'.format(case)

  # define files and saves
  files = ['StormFlash2d00004000.nc', 'StormFlash2d00008000.nc', 'StormFlash2d00012000.nc', 'StormFlash2d00016000.nc']
  saves = [f.replace('StormFlash2d',case) for f in files]
  saves = [f.replace('.nc','.pdf') for f in saves]

  for cnt in xrange(len(files)):
    ugrid = pyugrid.UGrid.from_ncfile(os.path.join(rootdir,datdir,testsuite,'wind',files[cnt]), load_data=True)
    # load depth and momentum
    d  = ugrid.data['depth'].data[0,:]
    mx = ugrid.data['m_x'].data[0,:]
    my = ugrid.data['m_y'].data[0,:]

    # calculate velocities
    u = mx/d
    v = my/d
    z = np.zeros(d.shape)

    x = ugrid.nodes[:,0]
    y = ugrid.nodes[:,1]

    triang = tri.Triangulation(ugrid.nodes[:, 0], ugrid.nodes[:, 1], triangles=ugrid.faces)

    fig = plt.figure()
    plt.gca().set_aspect('equal')
    #plt.triplot(triang, lw=0.2, color='black')

    levels = np.arange(0.0, 2.7e-8, 0.2e-8)
    cmap   = cm.get_cmap(name='jet', lut=None)
    cs     = plt.tricontourf(triang, np.sqrt(u**2+v**2), levels=levels, cmap=cmap, extend='both')
    plt.quiver(ugrid.nodes[:,0],ugrid.nodes[:,1],u,v,scale=1e-8, scale_units='xy')
    cbar = plt.colorbar(cs, shrink=0.65, extend='both')
    plt.xlabel('$x$')
    plt.ylabel('$y$')

    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.13, top=0.95)
    fig.set_size_inches(4.5,3.5)
    fig.savefig(os.path.join(rootdir,figdir,saves[cnt]))
    plt.close()

  return


#*******************************************************************************
def humpplot(testsuite, case):
# check if directory exists
  (destination) = os.path.isdir(os.path.join(rootdir,datdir,testsuite,case))
  if destination == False:
    print 'Directory for ' + case + ' does not exist. Skipping this testcase.'
    return 1
  else:
    print 'Plotting testcase {0:s}'.format(case)

  # define files and saves
  files = ['StormFlash2d00000000.nc', 'StormFlash2d00000400.nc']
  saves = [f.replace('StormFlash2d',case) for f in files]
  saves = [f.replace('.nc','.pdf') for f in saves]

  for cnt in xrange(len(files)):
    ugrid = pyugrid.UGrid.from_ncfile(os.path.join(rootdir,datdir,testsuite,'hump',files[cnt]), load_data=True)

    d  = ugrid.data['depth'].data[0,:]

    triang = tri.Triangulation(ugrid.nodes[:, 0], ugrid.nodes[:, 1], triangles=ugrid.faces)

    fig = plt.figure()
    plt.gca().set_aspect('equal')

    levels = np.arange(0.5, 0.61, 0.01)

    cmap   = cm.get_cmap(name='jet', lut=None)
    cs     = plt.tricontourf(triang, d, levels=levels, cmap=cmap, extend='both')
    plt.triplot(triang, lw=0.2, color='black')
    cbar = plt.colorbar(cs, shrink=0.65, extend='both')
    plt.xlabel('$x$')
    plt.ylabel('$y$')

    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.13, top=0.95)
    fig.set_size_inches(4.5,3.5)
    fig.savefig(os.path.join(rootdir,figdir,saves[cnt]))
    plt.close()

  return

#*******************************************************************************
def vortexplot(testsuite, case):
# check if directory exists
  (destination) = os.path.isdir(os.path.join(rootdir,datdir,testsuite,case))
  if destination == False:
    print 'Directory for ' + case + ' does not exist. Skipping this testcase.'
    return 1
  else:
    print 'Plotting testcase {0:s}'.format(case)

  # define files and saves
  files = ['StormFlash2d00000000.nc', 'StormFlash2d00004000.nc']

  saves = [f.replace('StormFlash2d',case) for f in files]
  saves = [f.replace('.nc','.pdf') for f in saves]

  for cnt in xrange(len(files)):
    ugrid = pyugrid.UGrid.from_ncfile(os.path.join(rootdir,datdir,testsuite,case,files[cnt]), load_data=True)

    d  = ugrid.data['depth'].data[0,:]

    triang = tri.Triangulation(ugrid.nodes[:, 0], ugrid.nodes[:, 1], triangles=ugrid.faces)

    fig = plt.figure()
    plt.gca().set_aspect('equal')

    levels = np.arange(0.97, 1.001, 0.001)

    cmap = cm.get_cmap(name='jet', lut=None)
    cs   = plt.tricontourf(triang, d, levels=levels, cmap=cmap, extend='both')
    plt.triplot(triang, lw=0.2, color='black')
    cbar = plt.colorbar(cs, shrink=0.65, extend='both')
    plt.xlabel('$x$')
    plt.ylabel('$y$')

    fig.subplots_adjust(left=0.1, right=1.05, bottom=0.12, top=0.97)
    fig.set_size_inches(7.5,3.5)
    fig.savefig(os.path.join(rootdir,figdir,saves[cnt]))
    plt.close()

  return

#*******************************************************************************
def eulervortexplot(testsuite, case):
# check if directory exists
  (destination) = os.path.isdir(os.path.join(rootdir,datdir,testsuite,case))
  if destination == False:
    print 'Directory for ' + case + ' does not exist. Skipping this testcase.'
    return 1
  else:
    print 'Plotting testcase {0:s}'.format(case)

  # define files and saves
  files = ['StormFlash2d00000000.nc', 'StormFlash2d00001000.nc']

  saves = [f.replace('StormFlash2d',case) for f in files]
  saves = [f.replace('.nc','.pdf') for f in saves]

  for cnt in xrange(len(files)):
    ugrid = pyugrid.UGrid.from_ncfile(os.path.join(rootdir,datdir,testsuite,case,files[cnt]), load_data=True)

    if case == 'eulerPTvortex':
      d = ugrid.data['pot_temp'].data[0,:]
    else:
      d = ugrid.data['energy'].data[0,:]

    triang = tri.Triangulation(ugrid.nodes[:, 0], ugrid.nodes[:, 1], triangles=ugrid.faces)

    fig = plt.figure()
    plt.gca().set_aspect('equal')

    if case[0:7] == 'eulerPT':
      levels = np.linspace(2.75E-3, 3.50E-3, 51)
    elif (case[0:7] == 'eulerTE'):
      levels = np.linspace(2.4,4.2,51)

    cmap = cm.get_cmap(name='jet', lut=None)
    cs   = plt.tricontourf(triang, d, levels=levels, cmap=cmap, extend='both')
    cbar = plt.colorbar(cs, shrink=0.65, extend='both')
    plt.xlabel('$x$')
    plt.ylabel('$y$')

    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.13, top=0.95)
    fig.set_size_inches(4.5,3.5)
    fig.savefig(os.path.join(rootdir,figdir,saves[cnt]))
    plt.close()

  return

#*******************************************************************************
def linadvsines(testsuite, case):
# check if directory exists
  (destination) = os.path.isdir(os.path.join(rootdir,datdir,testsuite,case))
  if destination == False:
    print 'Directory for ' + case + ' does not exist. Skipping this testcase.'
    return 1
  else:
    print 'Plotting testcase {0:s}'.format(case)

  # define files and saves
  files = ['StormFlash2d00000000.nc', 'StormFlash2d00001000.nc']

  saves = [f.replace('StormFlash2d',case) for f in files]
  saves = [f.replace('.nc','.pdf') for f in saves]

  for cnt in xrange(len(files)):
    ugrid = pyugrid.UGrid.from_ncfile(os.path.join(rootdir,datdir,testsuite,case,files[cnt]), load_data=True)

    d = ugrid.data['tracer'].data[0,:]

    triang = tri.Triangulation(ugrid.nodes[:, 0], ugrid.nodes[:, 1], triangles=ugrid.faces)

    fig = plt.figure()
    plt.gca().set_aspect('equal')

    levels = np.linspace(1.,3.,51)

    cmap = cm.get_cmap(name='jet', lut=None)
    cs   = plt.tricontourf(triang, d, levels=levels, cmap=cmap, extend='both')
    cbar = plt.colorbar(cs, shrink=0.65, extend='both')
    plt.xlabel('$x$')
    plt.ylabel('$y$')

    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.13, top=0.95)
    fig.set_size_inches(4.5,3.5)
    fig.savefig(os.path.join(rootdir,figdir,saves[cnt]))
    plt.close()

  return

#*******************************************************************************
# this is the main routine
def main(testsuite):
  # first create figure directory
  try:
    (destination) = os.makedirs(os.path.join(rootdir,figdir), 0755)
  except OSError:
    print 'Skipping creation of {0:s} because it exists already.'.format(figdir)

  print 'Plotting figures to figure directory. Please wait.'
  shockplot( testsuite, 'shockdry')
  shockplot( testsuite, 'shocknone')
  shockplot( testsuite, 'shockquad')
  shockplot( testsuite, 'shockshu')
  shockplot( testsuite, 'shockgiraldo')
  shockplot( testsuite, 'shockBJSV')
  shockplot( testsuite, 'shockKuSV')
  windplot(  testsuite, 'wind')
  beachplot( testsuite, 'beach')
  basinplot( testsuite, 'freebasinBJSV')
  basinplot( testsuite, 'freebasinKuSV')
  lakeplot(  testsuite, 'wellbalanced')
  humpplot(  testsuite, 'hump')
  vortexplot(testsuite, 'vortexstrong')
  vortexplot(testsuite, 'vortexweak')
# euler test cases
  eulervortexplot(testsuite, 'eulerPTvortex')
  eulervortexplot(testsuite, 'eulerTEvortex')
  eulersimplewave( testsuite, 'limited')
  eulersimplewave( testsuite, 'non-limited')
# linear advection
  linadvsines(testsuite, 'linadvsines')

  print 'Plotting has finished.'

  return

#*******************************************************************************
# run script only when it is actually called as a script
if __name__ == '__main__':
  # check whether console input is valid and set default if not
  if len(sys.argv) < 2:
    raise ValueError('Input has to include a testsuite directory name.')

  # read testsuite directory names from console input
  testsuite = sys.argv[1]

  # run main routine
  main(testsuite)






