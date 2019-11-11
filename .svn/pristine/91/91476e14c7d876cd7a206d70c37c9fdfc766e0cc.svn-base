!*******************************************************************************
!
! MODULE NAME:
!   CCI_postprocessing
! FUNCTION:
!   organize the different steps of postprocessing.
!   - compute track by weighted vorticity
!   - write down coordinates of track in temp file
!   - break run in case storm reaches its initial position
!   - compute energy, enstrophy and palinstrophy
! CONTAINS:
!-----------------------------------------------------------------
!
! NAME:
!   fvm_evolution
! FUNCTION:
!   this implements the actual evaluation in U(t+dt)= U(t) + dt* F(U,t)
! SYNTAX:
!   CALL fvm_initialize(grid, param)
! ON INPUT:
!   p_param: parameter data structure   TYPE (global_param)
! ON OUTPUT:
!   p_ghand: grid handling data structure   TYPE (grid_handle)
! CALLS:
!
! COMMENTS:
!*******************************************************************************
MODULE CCI_postprocessing

  USE FLASH_parameters
  USE IO_utils
  USE GRID_api
  USE ADV_update

  PRIVATE
  PUBLIC :: compute_vorticity, compute_track, compute_energy, &
            check_position, compute_palinstrophy, compute_enstrophy, postprocess

  CONTAINS

! Subroutine to compute track as a weighted mean of vorticity
!*******************************************************************************
  SUBROUTINE compute_vorticity(p_ghand, r_metrics_inv)

    IMPLICIT NONE

    TYPE (grid_handle),DIMENSION(:)                             :: p_ghand
  !------- For vorticity computation
    REAL (KIND = GRID_SR), DIMENSION(:,:,:), ALLOCATABLE        :: r_grad_u, r_grad_v
    REAL (KIND = GRID_SR), DIMENSION(:,:,:,:), ALLOCATABLE      :: r_laplace_u, r_laplace_v
    REAL (KIND = GRID_SR), DIMENSION(:,:,:), INTENT(in)         :: r_metrics_inv
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE, TARGET  :: r_deru, r_derv
    REAL (KIND = GRID_SR), DIMENSION(:,:),ALLOCATABLE           :: r_vort
    INTEGER (KIND = GRID_SI)                                    :: i_elmt, i_face, i_numelmt, &
                                                                   i_faceunknowns, i_numdofs, &
                                                                   i_alct, i_numedge, i_numnode
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), ALLOCATABLE       :: i_elementdofs
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE, TARGET  :: r_uv
    INTEGER (KIND = GRID_SI), DIMENSION(2)                      :: i_valvelo

    i_faceunknowns = GRID_femtypes%p_type(FEM_VELO)%sig%i_unknowns !instead of elementpoints
    i_numelmt      = p_ghand(i_timeplus)%i_enumfine
    i_numedge      = p_ghand(i_timeplus)%i_gnumfine
    i_numnode      = p_ghand(i_timeplus)%i_nnumber
    i_numdofs      = i_numelmt*i_faceunknowns

    !---- For diagnostics
    ALLOCATE(r_grad_u(i_numelmt,2,i_faceunknowns), r_grad_v(i_numelmt,2,i_faceunknowns), &
             r_laplace_u(i_numelmt, i_faceunknowns, 2, 2), r_laplace_v(i_numelmt, i_faceunknowns,2,2), &
             r_deru(2, i_numdofs), r_derv(2, i_numdofs), r_vort(1, i_numdofs), &
             i_elementdofs(i_faceunknowns, i_numelmt), r_uv(7,i_numdofs), STAT=i_alct)

    i_valvelo = (/ VAR2D_UVMNUS, VAR2D_UVMNUS+1/) !(/ VAR2D_UVZERO, VAR2D_UVZERO+1/)
    CALL grid_getinfo(p_ghand(i_timeplus), i_femtype=FEM_VELO, &
                      i_arraypoint=i_valvelo, i_elementdofs=i_elementdofs, r_dofvalues=r_uv)

    !------------ Compute vorticity etc.
    CALL vorticity(r_uv(1,:), r_uv(2,:), r_metrics_inv, r_grad_u, r_grad_v, &
                   r_laplace_u, r_laplace_v, p_ghand(i_timeplus), i_elementdofs)

    DO i_elmt = 1, i_numelmt
      DO i_face = 1, i_faceunknowns
        r_vort(1,(i_elmt-1)*i_faceunknowns+i_face) = r_grad_v(i_elmt,1,i_face) - r_grad_u(i_elmt,2,i_face)
        r_deru(1,(i_elmt-1)*i_faceunknowns+i_face) = r_grad_u(i_elmt,1,i_face)
        r_deru(2,(i_elmt-1)*i_faceunknowns+i_face) = r_grad_u(i_elmt,2,i_face)
        r_derv(1,(i_elmt-1)*i_faceunknowns+i_face) = r_grad_v(i_elmt,1,i_face)
        r_derv(2,(i_elmt-1)*i_faceunknowns+i_face) = r_grad_v(i_elmt,2,i_face)
      END DO
    END DO
    WRITE(*,*) maxval(r_uv), maxval(r_vort), minval(r_vort)


    CALL grid_putinfo(p_ghand(i_timeplus), l_finelevel=.TRUE., &
        i_arraypoint=(/VAR1D_VORT/), r_dofvalues=r_vort)

    CALL grid_putinfo(p_ghand(i_timeplus), l_finelevel=.TRUE., &
        i_arraypoint=(/VAR2D_DERU/), r_dofvalues=r_deru)

    CALL grid_putinfo(p_ghand(i_timeplus), l_finelevel=.TRUE., &
        i_arraypoint=(/VAR2D_DERV/), r_dofvalues=r_derv)

    DEALLOCATE(r_grad_u, r_grad_v, r_laplace_u, r_laplace_v, r_deru, r_derv, r_vort, r_uv)

    RETURN

  END SUBROUTINE compute_vorticity

! Subroutine to compute track as a weighted mean of vorticity
!*******************************************************************************
  SUBROUTINE compute_track(p_ghand, x_out, y_out)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(in)                              :: p_ghand
    REAL, INTENT(out)                                           :: x_out, y_out
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE, TARGET  :: r_vort, r_xy
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), ALLOCATABLE       :: i_elementdofs
    INTEGER                                                     :: i_numdofs, i_dofs, &
                                                                   i_numelmt, i_elmt, &
                                                                   i_faceunknowns, i_face
    REAL                                                        :: max_vorticity_dof, vorticity_bound
    INTEGER                                                     :: i_max_vorticity
    REAL                                                        :: x, y, vort, x_offset, y_offset, &
                                                                   x_tilde, y_tilde, weight, &
                                                                   weighted_upper, mesh_size, &
                                                                   weighted_upper_x, weighted_upper_y

    max_vorticity_dof = 0.0
    vorticity_bound   = 0.0002 !0.0006 !//0.0005 für 30_deg
    ! TODO viel zu kleine Vorticity-Werte!!!!!!!!1

    i_faceunknowns = GRID_femtypes%p_type(FEM_VELO)%sig%i_unknowns
    i_numelmt      = p_ghand%i_enumfine
    i_numdofs      = i_numelmt*i_faceunknowns

    ALLOCATE(r_vort(1, i_numdofs), r_xy(GRID_dimension,i_numdofs), &
             i_elementdofs(i_faceunknowns,i_numelmt))

    CALL grid_getinfo(p_ghand, i_femtype=FEM_VELO, i_arraypoint=(/ VAR1D_VORT /), &
                      r_dofvalues=r_vort)

    CALL grid_getinfo(p_ghand, i_femtype=FEM_ZETA, l_relative=.TRUE., &
                      l_finelevel=.TRUE., r_dofcoordinates=r_xy, &
                      i_elementdofs=i_elementdofs)

    ! Search for maximum vorticity on all DOFs
    DO  i_dofs = 1, i_numdofs
      IF (max_vorticity_dof <= ABS(r_vort(1, i_dofs))) THEN
        max_vorticity_dof = ABS(r_vort(1, i_dofs))
        i_max_vorticity   = i_dofs
      END IF
    END DO

    IF (max_vorticity_dof <= vorticity_bound*1.05) THEN
      vorticity_bound = 0.75 * vorticity_bound
      IF (max_vorticity_dof <= vorticity_bound*1.05) THEN
        vorticity_bound = 0.75 * vorticity_bound
      END IF
      IF (max_vorticity_dof <= vorticity_bound*1.05) THEN
        vorticity_bound = 0.75 * vorticity_bound
      END IF
      IF (max_vorticity_dof <= vorticity_bound*1.05) THEN
        vorticity_bound = 0.75 * vorticity_bound
      END IF
      IF (max_vorticity_dof <= vorticity_bound*1.05) THEN
        vorticity_bound = 0.75 * vorticity_bound
      END IF
    END IF

    x_offset         = 0.0
    y_offset         = 0.0
    weighted_upper   = 0.0
    weighted_upper_x = 0.0
    weighted_upper_y = 0.0

    DO i_elmt = 1, i_numelmt
      mesh_size = ABS(r_xy(1, i_elementdofs(1,i_elmt)) - r_xy(1, i_elementdofs(2,i_elmt)))
      DO i_face = 1, i_faceunknowns
        i_dofs = i_elementdofs(i_face,i_elmt)
        x = r_xy(1, i_dofs)
        y = r_xy(2, i_dofs)
        vort = abs(r_vort(1, i_dofs))

        IF ((vort >= vorticity_bound) .AND. ( x_offset == 0.) .AND. (y_offset == 0.)) THEN
          x_offset = x - 1500000.
          y_offset = y - 1000000.
        END IF

        x_tilde = x - x_offset
        y_tilde = y - y_offset

        IF (x_tilde < 0.)       x_tilde = x_tilde + 3333000.
        IF (x_tilde > 3333000.) x_tilde = x_tilde - 3333000.
        IF (y_tilde < 0.)       y_tilde = y_tilde + 2000000.
        IF (y_tilde > 2000000.) y_tilde = y_tilde - 2000000.


        IF (vort > vorticity_bound) THEN
          weight           = vort - vorticity_bound
          weight           = weight * mesh_size*mesh_size
          weighted_upper_x = weighted_upper_x + x_tilde*weight
          weighted_upper_y = weighted_upper_y + y_tilde*weight
          weighted_upper   = weighted_upper + weight
        END IF
      END DO
    END DO

    weighted_upper_x = weighted_upper_x/weighted_upper
    weighted_upper_y = weighted_upper_y/weighted_upper

    weighted_upper_x = weighted_upper_x + x_offset
    weighted_upper_y = weighted_upper_y + y_offset


    IF (weighted_upper_x < 0.) weighted_upper_x = weighted_upper_x + 3333000.
    IF (weighted_upper_x > 3333000.) weighted_upper_x = weighted_upper_x - 3333000.
    IF (weighted_upper_y < 0.) weighted_upper_y = weighted_upper_y + 2000000.
    IF (weighted_upper_y > 2000000.) weighted_upper_y = weighted_upper_y - 2000000.

    x_out = weighted_upper_x
    y_out = weighted_upper_y

    DEALLOCATE(r_vort, r_xy, i_elementdofs)

    RETURN

  END SUBROUTINE compute_track

!*******************************************************************************
  SUBROUTINE compute_energy(p_ghand, energy, energy_nd, r_element_vol)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(in)                              :: p_ghand
    REAL, INTENT(out)                                           :: energy, energy_nd
    REAL, DIMENSION(:), INTENT(in)                              :: r_element_vol

    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE, TARGET  :: r_vel,  r_xy
    INTEGER , DIMENSION(:,:), ALLOCATABLE, TARGET               :: i_elementdofs

    INTEGER                                                     :: i_dof, i_numdofs, i_faceunknowns, &
                                                                   i_elmt, i_face, i_numelmt
    REAL                                                        :: alpha, PI, velocity, u_drift, v_drift
    REAL                                                        :: r_energy_nd, r_energy
    INTEGER                                                     :: time
    REAL                                                        :: max_vol, min_vol

    i_faceunknowns = GRID_femtypes%p_type(FEM_VELO)%sig%i_unknowns !instead of elementpoints
    i_numelmt      = p_ghand%i_enumfine
    i_numdofs      = i_numelmt*i_faceunknowns
    max_vol = 0.0
    min_vol = 99999999999999999999999999999999999.0

    ALLOCATE(r_vel(2, i_numdofs), r_xy(GRID_dimension,i_numdofs), &
             i_elementdofs(i_faceunknowns, i_numelmt))

    CALL grid_getinfo(p_ghand, i_femtype=FEM_ZETA, l_finelevel=.TRUE., &
                      i_elementdofs=i_elementdofs)

    CALL grid_getinfo(p_ghand, i_femtype=FEM_VELO, &
                      i_arraypoint=(/ VAR2D_UVZERO, VAR2D_UVZERO+1 /), r_dofvalues=r_vel)
    energy    = 0.0
    energy_nd = 0.0
    velocity  = p_contr%phy%r_vdrift
    u_drift   = velocity



    DO i_elmt = 1, i_numelmt
      IF (r_element_vol(i_elmt) > max_vol) THEN
        max_vol = r_element_vol(i_elmt)
      ELSE IF (r_element_vol(i_elmt) < min_vol) THEN
        min_vol = r_element_vol(i_elmt)
      END IF
    END DO
    ! WRITE(*,*) min_vol, max_vol

    DO i_elmt=1,i_numelmt
      r_energy    = 0.0
      r_energy_nd = 0.0
      DO i_face=1,i_faceunknowns
        i_dof = i_elementdofs(i_face, i_elmt)
        r_energy    = r_energy + (r_vel(1,i_dof))**2 + (r_vel(2,i_dof))**2
        r_energy_nd = r_energy_nd + (r_vel(1,i_dof)-u_drift)**2 + (r_vel(2,i_dof))**2
      END DO
      r_energy    = r_energy*r_element_vol(i_elmt)/2000000/1600000/3
      r_energy_nd = r_energy_nd*r_element_vol(i_elmt)/2000000/1600000/3
      energy      = energy + r_energy
      energy_nd   = energy_nd + r_energy_nd
    END DO
    WRITE(*,*) energy_nd

    DEALLOCATE(r_vel, r_xy, i_elementdofs)

    RETURN

  END SUBROUTINE compute_energy

!*******************************************************************************
  SUBROUTINE check_position(x, count_steps, init_quitting)

    IMPLICIT NONE

    REAL, INTENT(in)                          ::  x
    INTEGER, INTENT(inout)                    ::  count_steps
    LOGICAL, INTENT(inout)                    ::  init_quitting
    REAL                                      ::  start_x

    start_x = p_contr%phy%r_xinit * (10**3)

    IF (.NOT. init_quitting) THEN
      IF (x > start_x + 400000.) THEN
        init_quitting = .TRUE.
      END IF
    END IF

    RETURN
  END SUBROUTINE check_position

!*******************************************************************************
  SUBROUTINE compute_palinstrophy(p_ghand, r_palin)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(in)                              :: p_ghand
    REAL, INTENT(out)                                           :: r_palin
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE, TARGET  :: r_derderu, r_derderv

    INTEGER                                                     :: i_dof, i_numdofs, &
                                                                   i_faceunknowns, i_numelmt

    i_faceunknowns = GRID_femtypes%p_type(FEM_VELO)%sig%i_unknowns !instead of elementpoints
    i_numelmt      = p_ghand%i_enumfine
    i_numdofs      = i_numelmt*i_faceunknowns

    ALLOCATE(r_derderu(3, i_numdofs), r_derderv(3, i_numdofs))

    r_palin = 0.

    CALL grid_getinfo(p_ghand, i_femtype=FEM_VELO, i_arraypoint= (/ VAR3D_DER2U /), r_dofvalues = r_derderu)
    CALL grid_getinfo(p_ghand, i_femtype=FEM_VELO, i_arraypoint= (/ VAR3D_DER2V /), r_dofvalues = r_derderv)


    DO i_dof = 1, i_numdofs
      r_palin = r_palin + 0.5*((r_derderv(1,i_dof)-r_derderu(2,i_dof))**2 + &
                               (r_derderv(2,i_dof)-r_derderu(3,i_dof))**2)
    END DO

    DEALLOCATE(r_derderu, r_derderv)

    RETURN

  END SUBROUTINE compute_palinstrophy

!*******************************************************************************
  SUBROUTINE compute_enstrophy(p_ghand, r_enst)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(in)                              :: p_ghand
    REAL, INTENT(out)                                           :: r_enst
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE, TARGET  :: r_deru, r_derv

    INTEGER                                                     :: i_dof, i_numdofs, &
                                                                   i_faceunknowns, i_numelmt

    i_faceunknowns =  GRID_femtypes%p_type(FEM_VELO)%sig%i_unknowns !instead of elementpoints
    i_numelmt      =  p_ghand%i_enumfine
    i_numdofs      =  i_numelmt*i_faceunknowns

    ALLOCATE(r_deru(2, i_numdofs), r_derv(2, i_numdofs))

    r_enst = 0.

    CALL grid_getinfo(p_ghand, i_femtype=FEM_VELO, i_arraypoint= (/ VAR2D_DERU /), r_dofvalues = r_deru)
    CALL grid_getinfo(p_ghand, i_femtype=FEM_VELO, i_arraypoint= (/ VAR2D_DERV /), r_dofvalues = r_derv)

    DO i_dof = 1, i_numdofs
      r_enst = r_enst + 0.5*(r_derv(1,i_dof)-r_deru(2,i_dof))**2
    END DO

    DEALLOCATE(r_deru, r_derv)

    RETURN

  END SUBROUTINE compute_enstrophy

!*******************************************************************************
  SUBROUTINE postprocess(p_ghand,x, y, r_element_vol)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(in)        :: p_ghand
    REAL                                  :: energy, energy_nd, palinst, &
                                             enstrophy,vorticity, max_vel, impuls
    INTEGER                               :: time
    REAL, INTENT(out)                     :: x, y
    REAL, DIMENSION(:), INTENT(in)        :: r_element_vol

    CALL max_velocity(p_ghand, max_vel)
    CALL compute_energy(p_ghand, energy, energy_nd, r_element_vol)
    CALL compute_palinstrophy(p_ghand, palinst)
    CALL compute_enstrophy(p_ghand, enstrophy)
    CALL vorticity_sum(p_ghand, vorticity)
    CALL compute_track(p_ghand, x, y)

    time = p_timestepinfo%r_modeltime

    IF (time == 1) THEN
      OPEN(unit=25, file='PP/e_kin.txt', action='write', status ='replace')
      WRITE(25,*) "       Zeit           x           y      energy   energy_nd ",&
      "   enstrophy    palinstr    vorticity   "
      CLOSE(25)

      OPEN(unit=26, file='PP/pp_cmg.txt', action='write', status ='replace')
      WRITE(26,*) "       Zeit           x           y    max_velocity"
      CLOSE(26)
    ENDIF

    OPEN(unit=25, file='PP/e_kin.txt', action='write', access = 'append')
    WRITE(25,*) time,",", x,",", y,",", energy,",", energy_nd,",", enstrophy,",", palinst,",", vorticity
    CLOSE(25)

    OPEN(unit=26, file='PP/pp_cmg.txt', action='write', access = 'append')
    WRITE(26,*) time,",", x,",", y,",", max_vel
    CLOSE(26)

    RETURN

  END SUBROUTINE postprocess

!*******************************************************************************
  SUBROUTINE max_velocity(p_ghand, max_vel)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(in)                              :: p_ghand
    REAL, INTENT(out)                                           :: max_vel
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE, TARGET  :: r_vel, r_xy

    INTEGER                                                     :: i_dof, i_numdofs, &
                                                                   i_faceunknowns, i_numelmt
  ! REAL                                                           ::  x, y
    INTEGER                                                     :: time

    max_vel = 0.
    i_faceunknowns = GRID_femtypes%p_type(FEM_VELO)%sig%i_unknowns !instead of elementpoints
    i_numelmt      = p_ghand%i_enumfine
    i_numdofs      = i_numelmt*i_faceunknowns

    ALLOCATE(r_vel(2, i_numdofs), r_xy(GRID_dimension,i_numdofs))

    CALL grid_getinfo(p_ghand, i_femtype=FEM_ZETA, l_finelevel=.TRUE., &
                      r_dofcoordinates=r_xy)

    CALL grid_getinfo(p_ghand, i_femtype=FEM_VELO, &
                      i_arraypoint=(/ VAR2D_UVZERO, VAR2D_UVZERO+1 /), r_dofvalues=r_vel)

    DO i_dof = 1, i_numdofs
      IF (sqrt(r_vel(1,i_dof)**2 + r_vel(2,i_dof)**2) > max_vel) THEN
          max_vel = sqrt(r_vel(1,i_dof)**2 + r_vel(2,i_dof)**2)
        !  x = r_xy(1, i_dof)
        !  y = r_xy(2, i_dof)
      END IF
    END DO

    DEALLOCATE(r_vel, r_xy)
  RETURN

  END SUBROUTINE max_velocity

!*******************************************************************************
  SUBROUTINE vorticity_sum(p_ghand, vorticity)

    IMPLICIT NONE

    TYPE (grid_handle), INTENT(in)                              :: p_ghand
    REAL, INTENT(inout)                                         :: vorticity
    REAL (KIND = GRID_SR), DIMENSION(:,:), ALLOCATABLE, TARGET  :: r_vort
    INTEGER (KIND = GRID_SI), DIMENSION(:,:), ALLOCATABLE       :: i_elementdofs
    INTEGER                                                     :: i_numdofs, i_dofs, &
                                                                   i_numelmt, i_elmt, &
                                                                   i_faceunknowns, i_face

    vorticity      = 0.0
    i_faceunknowns =  GRID_femtypes%p_type(FEM_VELO)%sig%i_unknowns !instead of elementpoints
    i_numelmt      =  p_ghand%i_enumfine
    i_numdofs      =  i_numelmt*i_faceunknowns

    ALLOCATE(r_vort(1, i_numdofs))

    CALL grid_getinfo(p_ghand, i_femtype=FEM_VELO, i_arraypoint=(/ VAR1D_VORT /), &
                      r_dofvalues=r_vort)

    DO i_dofs=1,i_numdofs
      vorticity = vorticity + r_vort(1, i_dofs)
    END DO

    DEALLOCATE(r_vort)

    RETURN

  END SUBROUTINE vorticity_sum

!*******************************************************************************
END MODULE CCI_postprocessing
