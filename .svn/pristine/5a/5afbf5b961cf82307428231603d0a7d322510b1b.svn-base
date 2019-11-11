# This is a script to create plots for the testsuite documentation.
#
# Usage: has to be called as a script with one argument.
# The argument has to be a name of testsuite directory and needs to be present in the compile directory.
# This script creates plots for all testsuite cases.

#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
def shockplot(testsuite, case):
  # define files and saves
  files = []
  saves = []
  for ncfile in sorted(os.listdir(os.path.join(rootdir,datdir,testsuite,case))):
    if ncfile.endswith('.nc'):
      files.append(ncfile)
  saves = [f.replace('StormFlash2d',case) for f in files]
  saves = [f.replace('.nc','.pdf') for f in saves]

  # load data
  sys.stdout.write('t = ')
  for cnt in xrange(len(files)):
    ugrid = pyugrid.UGrid.from_ncfile(os.path.join(rootdir,datdir,testsuite,case,files[cnt]), load_data=True)

    sys.stdout.write('{0:5.2f}, '.format(ugrid.data['depth'].time[0]))
    sys.stdout.flush()

    # cross section along y = 1.0
    vnames = ['bathy', 'depth', 'm_x', 'm_y', 'exact_depth', 'exact_m_x', 'exact_m_y']
    xv, yv, varpts = sfug.get_crosssection([0,1], [9.99999,1], ugrid, vnames)

    # fluid depth plot
    fig = plt.figure(1)
    plt.clf()
    plt.subplot(311)
    plt.plot(xv, varpts['bathy']+varpts['exact_depth'], color=col[0], linewidth=2.0, label='exact')
    plt.plot(xv, varpts['bathy']+varpts['depth']      , color=col[1], linewidth=1.0, label='numerical')
    plt.xlim([0, 10])
    if case == 'shockdry':
      plt.ylim([-0.05, 0.3])
    else:
      plt.ylim([0.5, 3.5])
    plt.ylabel('$h+b$')
    #plt.xlabel('$x$')
    plt.title('height')
    plt.legend()

    # x-momentum
    plt.subplot(312)
    plt.plot(xv, varpts['exact_m_x'], color=col[0], linewidth=2.0, label='exact')
    plt.plot(xv, varpts['m_x']      , color=col[1], linewidth=1.0, label='numerical')
    plt.xlim([0, 10])
    #plt.xlabel('$x$')
    plt.ylabel('$hu$')
    plt.title('$x$-momentum')
    #plt.legend()

    # y-momentum
    plt.subplot(313)
    plt.plot(xv, varpts['exact_m_y'], color=col[0], linewidth=2.0, label='exact')
    plt.plot(xv, varpts['m_y']      , color=col[1], linewidth=1.0, label='numerical')
    plt.xlim([0, 10])
    plt.xlabel('$x$')
    plt.ylabel('$hv$')
    plt.title('$y$-momentum')
    #plt.legend()

    plt.suptitle('Testcase: {0:s}, $t={1:5.2f}$'.format(case, ugrid.data['depth'].time[0]))
    #plt.show()

    fig.set_size_inches(21/2.54,29.7/2.54)
    fig.savefig(os.path.join(rootdir,figdir,saves[cnt]))
    
  print
  plt.close()

  return

#-------------------------------------------------------------------------------
def beachplot(testsuite, case):
  # define files and saves
  files = []
  saves = []
  for ncfile in sorted(os.listdir(os.path.join(rootdir,datdir,testsuite,case))):
    if ncfile.endswith('.nc'):
      files.append(ncfile)
  saves = [f.replace('StormFlash2d',case) for f in files]
  saves = [f.replace('.nc','.pdf') for f in saves]

  # load data
  sys.stdout.write('t = ')
  for cnt in xrange(len(files)):
    ugrid = pyugrid.UGrid.from_ncfile(os.path.join(rootdir,datdir,testsuite,case,files[cnt]), load_data=True)

    sys.stdout.write('{0:5.1f}, '.format(ugrid.data['depth'].time[0]))
    sys.stdout.flush()

    # cross section along y = 1.0
    vnames = ['bathy', 'depth', 'm_x', 'm_y', 'exact_depth', 'exact_m_x', 'exact_m_y']
    xv, yv, varpts = sfug.get_crosssection([0,1], [14000,1], ugrid, vnames)

    fig = plt.figure(1)
    plt.clf()
    plt.plot(xv, varpts['bathy'],                 color='0.6' , linewidth=1.0,
             linestyle='dashed', label='bathymetry')
    plt.plot(xv, varpts['bathy']+varpts['depth'], color=col[1], linewidth=1.0, label='total height')
    plt.xlim([0, 14000])
    plt.ylim([0, 5.5])
    plt.xlabel('$x$')
    plt.ylabel('$z$')
    plt.legend()

    plt.suptitle('Testcase: {0:s}, $t={1:5.2f}$'.format(case, ugrid.data['depth'].time[0]))
    plt.show()

    #fig.subplots_adjust(left=0.2, right=0.96, bottom=0.17, top=0.97)
    #fig.set_size_inches(4.5,2.5)
    fig.savefig(os.path.join(rootdir,figdir,saves[cnt]))

  print
  plt.close()

  return

#-------------------------------------------------------------------------------
def basinplot(testsuite, case):
  # define files and saves
  files = []
  saves = []
  for ncfile in sorted(os.listdir(os.path.join(rootdir,datdir,testsuite,case))):
    if ncfile.endswith('.nc'):
      files.append(ncfile)
  saves = [f.replace('StormFlash2d',case) for f in files]
  saves = [f.replace('.nc','.pdf') for f in saves]

  # load data
  sys.stdout.write('t = ')
  for cnt in xrange(len(files)):
    ugrid = pyugrid.UGrid.from_ncfile(os.path.join(rootdir,datdir,testsuite,case,files[cnt]), load_data=True)

    sys.stdout.write('{0:5.1f}, '.format(ugrid.data['depth'].time[0]))
    sys.stdout.flush()

    # cross section along y = 1.0
    vnames = ['bathy', 'depth', 'm_x', 'm_y', 'exact_depth', 'exact_m_x', 'exact_m_y']
    xv, yv, varpts = sfug.get_crosssection([-4000,0], [4000,0], ugrid, vnames)

    fig = plt.figure(1)
    plt.clf()
    plt.subplot(311)
    plt.plot(xv, varpts['bathy'],                       color='0.6',  linewidth=1.0,
             linestyle='dashed', label='bathymetry')
    plt.plot(xv, varpts['bathy']+varpts['exact_depth'], color=col[0], linewidth=2.0, label='exact')
    plt.plot(xv, varpts['bathy']+varpts['depth']      , color=col[1], linewidth=1.0, label='numerical')
    plt.xlim([-4000, 4000])
    plt.ylim([0, 2.5])
    #plt.xlabel('$x$')
    plt.ylabel('$z$')
    plt.legend()

    # x-momentum
    plt.subplot(312)
    plt.plot(xv, varpts['exact_m_x'], color=col[0], linewidth=2.0, label='exact')
    plt.plot(xv, varpts['m_x']      , color=col[1], linewidth=1.0, label='numerical')
    plt.xlim([-4000, 4000])
    #plt.xlabel('$x$')
    plt.ylabel('$hu$')
    plt.title('$x$-momentum')
    #plt.legend()

    # y-momentum
    plt.subplot(313)
    plt.plot(xv, varpts['exact_m_y'], color=col[0], linewidth=2.0, label='exact')
    plt.plot(xv, varpts['m_y']      , color=col[1], linewidth=1.0, label='numerical')
    plt.xlim([-4000, 4000])
    plt.xlabel('$x$')
    plt.ylabel('$hv$')
    plt.title('$y$-momentum')
    #plt.legend()

    plt.suptitle('Testcase: {0:s}, $t={1:5.2f}$'.format(case, ugrid.data['depth'].time[0]))
    plt.show()

    fig.set_size_inches(21/2.54,29.7/2.54)
    fig.savefig(os.path.join(rootdir,figdir,saves[cnt]))

  print
  plt.close()

  return

#-------------------------------------------------------------------------------

def wellbplot(testsuite, case):
  # define files and saves
  files = []
  saves = []
  for ncfile in sorted(os.listdir(os.path.join(rootdir,datdir,testsuite,case))):
    if ncfile.endswith('.nc'):
      files.append(ncfile)
  saves = [f.replace('StormFlash2d',case) for f in files]
  saves = [f.replace('.nc','.pdf') for f in saves]

  # load data
  sys.stdout.write('t = ')
  for cnt in xrange(len(files)):
    ugrid = pyugrid.UGrid.from_ncfile(os.path.join(rootdir,datdir,testsuite,case,files[cnt]), load_data=True)

    sys.stdout.write('{0:5.1f}, '.format(ugrid.data['depth'].time[0]))
    sys.stdout.flush()

    # cross section along y = 1.0
    vnames = ['bathy', 'depth', 'm_x', 'm_y', 'exact_depth', 'exact_m_x', 'exact_m_y']
    xv, yv, varpts = sfug.get_crosssection([0,15], [30,15], ugrid, vnames)

    fig = plt.figure(1)
    plt.clf()
    plt.plot(xv, varpts['bathy'],                 color='0.6',  linewidth=1.0,
             linestyle='dashed', label='bathymetry')
    plt.plot(xv, varpts['bathy']+varpts['depth'], color=col[1], linewidth=1.0, label='total height')
    plt.xlim([0, 30])
    plt.ylim([0, 1.2])
    plt.xlabel('$x$')
    plt.ylabel('$h+b$ or $b$')
    plt.legend()

    plt.suptitle('Testcase: {0:s}, $t={1:5.2f}$'.format(case, ugrid.data['depth'].time[0]))
    plt.show()

    fig.savefig(os.path.join(rootdir,figdir,saves[cnt]))

  print
  plt.close()

  return

#-------------------------------------------------------------------------------

def windplot(testsuite, case):
  # define files and saves
  files = []
  saves = []
  for ncfile in sorted(os.listdir(os.path.join(rootdir,datdir,testsuite,case))):
    if ncfile.endswith('.nc'):
      files.append(ncfile)
  saves = [f.replace('StormFlash2d',case) for f in files]
  saves = [f.replace('.nc','.pdf') for f in saves]

  sys.stdout.write('t = ')
  for cnt in xrange(len(files)):
    ugrid = pyugrid.UGrid.from_ncfile(os.path.join(rootdir,datdir,testsuite,'wind',files[cnt]), load_data=True)

    sys.stdout.write('{0:5.1f}, '.format(ugrid.data['depth'].time[0]))
    sys.stdout.flush()

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

    fig = plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('equal')
    #plt.triplot(triang, lw=0.2, color='black')

    levels = np.arange(0.0, 2.7e-8, 0.2e-8)
    cmap   = cm.get_cmap(name='jet', lut=None)
    cs     = plt.tricontourf(triang, np.sqrt(u**2+v**2), levels=levels, cmap=cmap, extend='both')
    plt.quiver(ugrid.nodes[:,0],ugrid.nodes[:,1],u,v,scale=1e-8, scale_units='xy')
    cbar = plt.colorbar(cs, shrink=0.65, extend='both')
    plt.xlabel("$x$")
    plt.ylabel("$y$")

    plt.suptitle('Testcase: {0:s}, $t={1:5.2f}$'.format(case, ugrid.data['depth'].time[0]))
    plt.show()

    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.13, top=0.95)
    fig.set_size_inches(4.5,3.5)
    fig.savefig(os.path.join(rootdir,figdir,saves[cnt]))

  print
  plt.close()

  return

#-------------------------------------------------------------------------------

def humpplot(testsuite, case):
  # define files and saves
  files = []
  saves = []
  for ncfile in sorted(os.listdir(os.path.join(rootdir,datdir,testsuite,case))):
    if ncfile.endswith('.nc'):
      files.append(ncfile)
  saves = [f.replace('StormFlash2d',case) for f in files]
  saves = [f.replace('.nc','.pdf') for f in saves]

  sys.stdout.write('t = ')
  for cnt in xrange(len(files)):
    ugrid = pyugrid.UGrid.from_ncfile(os.path.join(rootdir,datdir,testsuite,case,files[cnt]), load_data=True)

    sys.stdout.write('{0:5.2f}, '.format(ugrid.data['depth'].time[0]))
    sys.stdout.flush()

    d  = ugrid.data['depth'].data[0,:]

    triang = tri.Triangulation(ugrid.nodes[:, 0], ugrid.nodes[:, 1], triangles=ugrid.faces)

    fig = plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('equal')
    plt.triplot(triang, lw=0.2, color='black')

    levels = np.arange(0.5, 0.61, 0.01)

    cmap   = cm.get_cmap(name='jet', lut=None)
    plt.tricontour(triang, d, levels=levels,
                colors=['0.25', '0.5', '0.5', '0.5', '0.5'],
                linewidths=[1.0, 0.5, 0.5, 0.5, 0.5])
    cbar = plt.colorbar(cs, shrink=0.65, extend='both', ticks=levels[::5])
    plt.xlabel('$x$')
    plt.ylabel('$y$')

    plt.suptitle('Testcase: {0:s}, $t={1:5.2f}$'.format(case, ugrid.data['depth'].time[0]))
    plt.show()

    #fig.subplots_adjust(left=0.2, right=0.96, bottom=0.17, top=0.97)
    #fig.set_size_inches(4.5,4)
    fig.savefig(os.path.join(rootdir,figdir,saves[cnt]))

  print
  plt.close()

  return

#-------------------------------------------------------------------------------

def vortexplot(testsuite, case):
  # define files and saves
  files = []
  saves = []
  for ncfile in sorted(os.listdir(os.path.join(rootdir,datdir,testsuite,case))):
    if ncfile.endswith('.nc'):
      files.append(ncfile)
  saves = [f.replace('StormFlash2d',case) for f in files]
  saves = [f.replace('.nc','.pdf') for f in saves]

  sys.stdout.write('t = ')
  for cnt in xrange(len(files)):
    ugrid = pyugrid.UGrid.from_ncfile(os.path.join(rootdir,datdir,testsuite,case,files[cnt]), load_data=True)

    sys.stdout.write('{0:5.2f}, '.format(ugrid.data['depth'].time[0]))
    sys.stdout.flush()

    d  = ugrid.data['depth'].data[0,:]
    mx = ugrid.data['m_x'].data[0,:]
    my = ugrid.data['m_y'].data[0,:]

    triang = tri.Triangulation(ugrid.nodes[:, 0], ugrid.nodes[:, 1], triangles=ugrid.faces)

    fig = plt.figure(1)
    plt.clf()
    plt.subplot(311)
    plt.gca().set_aspect('equal')
    plt.triplot(triang, lw=0.2, color='black')

    cmap   = cm.get_cmap(name='jet', lut=None)
    cs     = plt.tricontourf(triang, d, cmap=cmap, extend='both')
    cbar = plt.colorbar(cs, shrink=0.65, extend='both')
    plt.xlabel('$x$')
    plt.ylabel('$y$')

    plt.subplot(312)
    plt.gca().set_aspect('equal')
    plt.triplot(triang, lw=0.2, color='black')

    cmap   = cm.get_cmap(name='jet', lut=None)
    cs     = plt.tricontourf(triang, mx, cmap=cmap, extend='both')
    cbar = plt.colorbar(cs, shrink=0.65, extend='both')
    plt.xlabel('$x$')
    plt.ylabel('$y$')

    plt.subplot(313)
    plt.gca().set_aspect('equal')
    plt.triplot(triang, lw=0.2, color='black')

    cmap   = cm.get_cmap(name='jet', lut=None)
    cs     = plt.tricontourf(triang, my, cmap=cmap, extend='both')
    cbar = plt.colorbar(cs, shrink=0.65, extend='both')
    plt.xlabel('$x$')
    plt.ylabel('$y$')

    plt.suptitle('Testcase: {0:s}, $t={1:5.2f}$'.format(case, ugrid.data['depth'].time[0]))
    plt.show()

    fig.set_size_inches(21/2.54,29.7/2.54)
    fig.savefig(os.path.join(rootdir,figdir,saves[cnt]))

  print
  plt.close()

  return

#-------------------------------------------------------------------------------
# this is the main routine
def main(testsuite, icase):
  # check if testsuite directory exists
  if (not os.path.isdir(os.path.join(rootdir,datdir,testsuite))):
    print 'Testsuite directory does not exist!'
    return

  # get testcase to be plotted and check if it exists
  tests = sorted(os.listdir(os.path.join(rootdir,datdir,testsuite)))
  if icase == None:
    print 'Existing testcases: ' + ', '.join(map(str, tests))
    icase = raw_input('Choose testcase (default = all):\n')

  if (icase != 'all' and icase != ''):
    # check if testcase directory exists
    if (not icase in tests):
      print 'Testcase {0:s} does not exist! Skipping this testcase.'.format(icase)
      return
    else:
      tests = [icase]

  # create figure directory
  try:
    (destination) = os.makedirs(os.path.join(rootdir,figdir), 0755)
  except OSError:
    print 'Skipping creation of {0:s} because it exists already.'.format(os.path.join(rootdir,figdir))
  print 'Plotting all figures to {0:s} directory. Please wait.'.format(os.path.join(rootdir,figdir))

  for case in tests:
    print 'Plotting testcase {0:s}'.format(case)

    if (case[0:5] == 'shock'):
      shockplot(testsuite,  case)
    elif (case == 'wind'):
      windplot(testsuite,   case)
    elif (case == 'beach'):
      beachplot(testsuite,  case)
    elif (case[0:5] == 'freeb'):
      basinplot(testsuite,  case)
    elif (case[0:4] == 'well'):
      wellbplot(testsuite,   case)
    elif (case == 'hump'):
      humpplot(testsuite,   case)
    elif (case[0:4] == 'vort'):
      vortexplot(testsuite, case)

  print 'Plotting has finished.'
  return

#-------------------------------------------------------------------------------
# run script only when it is actually called as a script
if __name__ == '__main__':
  # check whether console input is valid and set default if not
  if (len(sys.argv) < 2 or len(sys.argv) > 3):
    print 'Input has to include a testsuite directory name and optionally a testcase name.'
  else:
    if len(sys.argv) == 2:
      case = None
    else:
      case = sys.argv[2]

    main(sys.argv[1], case)
