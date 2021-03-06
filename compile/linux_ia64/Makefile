##################################################################
#   STORMFLASH                                                   #
#   Adaptive STORM surge model based on the                      #
#   FLexible Adaptive Semi-Lagrangian Hack                       #
#   written in Fortran 90                                        #
##################################################################

# MACHINE
MACHINE = linux_ia64

# SET MAKETHING CORRESPONDING TO MACHINE
MAKETHING = DGM

# OPTIMIZATION SETTINGS [debug|norm|opt]
MODE := opt

# OpenMP parallelization
OMP := no

# LIBRARY SETTINGS [yes|no]
# Usage: $(MAKE) USE_NETCDF=yes
USE_NETCDF := yes

# SET MAIN DIRECTORY PATH
# !! This has to be alterd by user !!
ROOTDIR = $(HOME)/rechnen/amatos

# SET atlas/blas DIRECTORY PATH
# !! This has to be alterd by user !!
BLASLIB = -L/usr/lib64 -lblas

# SET LAPACK DIRECTORY PATH/BLAS INCLUDED IN INTEL MKL
# !! This has to be alterd by user !!
LAPACKLIB = -L/usr/lib64 -llapack -lpthread -lc #-lmkl_ia32 -lguide -lptrhead
LAPACKINC =

# SET NETCDF PATHS
NETCDFLIB =
NETCDFINC =
NCUGRIDLIB =
ifeq ($(strip $(USE_NETCDF)), yes)
  NETCDFLIB = -L/home/system/programming/intel-netcdf/lib64 -lnetcdff -lnetcdf
  NETCDFINC = -I/home/system/programming/intel-netcdf/include
  NCUGRIDLIB = -lncugrid
endif

# SET MORE DIRECTORY PATHS
LIBDIR  = $(ROOTDIR)/amatos2d/lib/$(MACHINE)
INCDIR  = $(ROOTDIR)/amatos2d/include/$(MACHINE)
MODDIR  = $(ROOTDIR)/amatos2d/include/$(MACHINE)
MAINDIR = $(ROOTDIR)/StormFlash2d

SRCDIR  = $(MAINDIR)/src/flash
SYSDIR  = $(MAINDIR)/src/system/ifortran
TIMDIR  = $(MAINDIR)/src/timing
DATDIR  = $(MAINDIR)/data
BUILDIR = $(MAINDIR)/compile/$(MACHINE)

# SET PATH TO GRID_API.mod
AMATOS = $(MODDIR)/grid_api.mod

# SET THE ENDING OF MODULE FILES
MODEND = mod

#----------------------------------------------------------------#
# FLAGS FOR LINUX / Intel Fortran Compiler                       #
#----------------------------------------------------------------#
F90     = ifort
cc      = icc
CC      = icc
LOADER  = ifort -i-dynamic

# -------------- next for debugging -----------------------------#
ifeq ($(strip $(MODE)),debug)
  FFLAGS  = -g -fpic
  CFLAGS  = -g
  LDFLAGS = -g -fpic
endif
# -------------- next for normal compilation --------------------#
ifeq ($(strip $(MODE)),norm)
  FFLAGS  =
  CFLAGS  =
  LDFLAGS = -fpic
endif
# -------------- next with aggresive optimization (Pentium M) ---#
ifeq ($(strip $(MODE)),opt)
  FFLAGS  = -O3 -fpic
  CFLAGS  = -O3
  LDFLAGS = -O3 -fpic
endif

#------------- Use OpenMP ---------------------------------------#
ifeq ($(strip $(OMP)),yes)
  LIBS    += -openmp
  ##-lgomp -lpthread
  FFLAGS  += -openmp
  LDFLAGS += -openmp
  CFLAGS  += -openmp
endif

# --------------------- include NetCDF or not -------------------#
ifeq ($(strip $(USE_NETCDF)),yes)
  FFLAGS += -DUSE_NETCDF
endif

# -------- next flag for using dummy graphics library -----------#
LIBS = -L$(LIBDIR) -lamatos $(NCUGRIDLIB) $(NETCDFLIB) $(LAPACKLIB) \
       $(BLASLIB) $(CCLIB)
INCS = $(LAPACKINC) $(NETCDFINC) -I$(INCDIR) -I$(MODDIR)

#----------------------------------------------------------------#
# common stuff                                                   #
#----------------------------------------------------------------#

include $(MAINDIR)/compile/Makefiles/Makefile.common
-include $(MAINDIR)/compile/Makefiles/dependency.mk
-include $(MAINDIR)/compile/Makefiles/Makefile.Testsuite
