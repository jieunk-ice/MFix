!  ======================================================================
!  usr1.f  -  dynamic cyclone char-return coupling
!             (shared by gasifier_2d.mfx and gasifier_3d.mfx)
!  ----------------------------------------------------------------------
!  Closes the recirculation loop: each timestep this routine integrates
!  the char mass flux leaving the TOP outlet (elutriation), then sets the
!  char-return point-source rate to a fraction (the cyclone collection
!  efficiency) of that flux:
!
!      ps_massflow_s(RETURN_PS, 1) = ETA * (char mass leaving the top)
!
!  USR1 is called once per timestep, so the return tracks elutriation as
!  the bed evolves. Requires CALL_USR = .True. in the .mfx (already set).
!
!  >>> VERIFY (could not be tested without building) <<<
!    * RETURN_PS / OUTLET_BC indices below must match the .mfx.
!    * Face-area array AXZ = north (xz-plane, y-normal) face area, and
!      V_s(.,1) = north-face y-velocity of the cell. Confirm for your
!      MFiX version; swap AXZ/AYZ/AXY + the velocity component if the
!      outlet is on a different face.
!    * GLOBAL_ALL_SUM interface (DMP reduction) - adjust if your build
!      uses a different name/signature.
!  ======================================================================

      SUBROUTINE USR1

      USE bc
      USE compar
      USE fldvar,      only: ROP_s, X_s, V_s
      USE functions,   only: funijk, jm_of, fluid_at, is_on_mype_owns
      USE geometry,    only: axz
      USE mpi_utility, only: global_all_sum
      USE param1,      only: zero
      USE ps,          only: ps_massflow_s
      USE run,         only: time

      IMPLICIT NONE

! Configuration
!---------------------------------------------------------------------
      INTEGER, PARAMETER :: OUTLET_BC = 2     ! top pressure-outflow BC
      INTEGER, PARAMETER :: RETURN_PS = 3     ! char-return point source
      INTEGER, PARAMETER :: CHAR_PHASE = 1    ! solids phase carrying char
! Cyclone collection efficiency (fraction of elutriated char returned)
      DOUBLE PRECISION, PARAMETER :: ETA = 0.90d0
! Don't start returning until a little char has had time to reach the top
      DOUBLE PRECISION, PARAMETER :: T_START = 1.0d0

! Local variables
!---------------------------------------------------------------------
      INTEGER :: I, J, K, IJK, IJKM, lun
      DOUBLE PRECISION :: vface, flux_local, flux_global
! Throttle the CSV log to ~every LOG_DT seconds (not every timestep)
      DOUBLE PRECISION, PARAMETER :: LOG_DT = 0.05d0
      DOUBLE PRECISION, SAVE :: t_next = ZERO

!---------------------------------------------------------------------

      flux_local = ZERO

! Integrate upward char mass flux across the outlet plane. BC_J_s = BC_J_n
! for a planar boundary; the interior fluid cell is one row below (jm_of).
      J = BC_J_s(OUTLET_BC)
      DO K = BC_K_b(OUTLET_BC), BC_K_t(OUTLET_BC)
      DO I = BC_I_w(OUTLET_BC), BC_I_e(OUTLET_BC)
         IJK  = FUNIJK(I, J, K)
         IJKM = JM_OF(IJK)                       ! interior fluid cell
         IF (.NOT. FLUID_AT(IJKM)) CYCLE
         IF (.NOT. IS_ON_myPE_owns(I, J-1, K)) CYCLE   ! avoid double count
         vface = V_s(IJKM, CHAR_PHASE)           ! north-face y-velocity
         IF (vface <= ZERO) CYCLE                ! count outflow only
         flux_local = flux_local                                         &
            + ROP_s(IJKM,CHAR_PHASE) * X_s(IJKM,CHAR_PHASE,Char)         &
              * vface * AXZ(IJKM)
      ENDDO
      ENDDO

! Sum across MPI ranks (no-op in serial)
      CALL GLOBAL_ALL_SUM(flux_local, flux_global)

! Set the return rate (kg/s). Guard against negatives and the warm-up.
      IF (TIME >= T_START) THEN
         PS_MASSFLOW_S(RETURN_PS, CHAR_PHASE) = MAX(ZERO, ETA*flux_global)
      ELSE
         PS_MASSFLOW_S(RETURN_PS, CHAR_PHASE) = ZERO
      ENDIF

! Log the loop (I/O rank only, throttled) for postprocess.py
      IF (myPE == PE_IO .AND. TIME >= t_next) THEN
         OPEN(newunit=lun, file='recirc.csv', position='append', action='write')
         WRITE(lun,'(ES12.5,",",ES12.5,",",ES12.5)') TIME, flux_global,  &
               PS_MASSFLOW_S(RETURN_PS, CHAR_PHASE)
         CLOSE(lun)
         t_next = TIME + LOG_DT
      ENDIF

      RETURN
      END SUBROUTINE USR1
