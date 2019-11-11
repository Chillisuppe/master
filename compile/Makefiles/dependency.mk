FLASH_parameters.o:
MISC_system.o:
MISC_utils.o:
MISC_bathy.o:
MISC_diag.o: FLASH_parameters.o DG_initial.o MISC_timing.o DG_utils.o
MISC_quad.o: DG_utils.o MISC_utils.o

IO_utils.o: FLASH_parameters.o MISC_timing.o MISC_system.o DG_initial.o IO_equation.o DG_timestepping.o
IO_equation.o: FLASH_parameters.o
IO_paraplot_dg.o: FLASH_parameters.o DG_equation.o

DG_errorestimate.o: FLASH_parameters.o DG_equation.o

DG_initial.o: FLASH_parameters.o DG_equation.o IO_equation.o MISC_quad.o MISC_bathy.o
DG_utils.o:
DG_limiter.o: FLASH_parameters.o DG_utils.o DG_limiter_utils.o DG_equation.o
DG_limiter_utils.o: FLASH_parameters.o
DG_boundary.o: FLASH_parameters.o DG_equation.o IO_equation.o DG_initial.o
DG_equation.o: FLASH_parameters.o IO_equation.o DG_utils.o
DG_riemann_solver.o: DG_equation.o
DG_flux.o: FLASH_parameters.o DG_boundary.o DG_riemann_solver.o DG_equation.o IO_equation.o
DG_timestepping.o: FLASH_parameters.o DG_limiter.o DG_flux.o DG_time_utils.o
DG_cfl.o: FLASH_parameters.o DG_equation.o IO_equation.o
ADV_dg.o: FLASH_parameters.o MISC_timing.o MISC_utils.o IO_paraplot_dg.o IO_utils.o DG_errorestimate.o DG_initial.o DG_timestepping.o DG_utils.o MISC_diag.o DG_cfl.o
TAM_main.o: FLASH_parameters.o IO_utils.o ADV_dg.o DG_equation.o

ifeq ($(strip $(USE_NETCDF)), yes)
  IO_netcdfplot.o: FLASH_parameters.o DG_equation.o
  ADV_dg.o: IO_netcdfplot.o
endif

DG_initial.o: DG_storm_holland.o
DG_storm_drag.o:
DG_storm_holland.o: DG_storm_drag.o
