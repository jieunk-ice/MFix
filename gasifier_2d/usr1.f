!  ======================================================================
!  usr1.f  -  dynamic cyclone char-return coupling + carbon-flux logging
!  ----------------------------------------------------------------------
!  Each timestep (CALL_USR = .True.) this routine:
!    1. integrates the char mass flux leaving the TOP outlet (elutriation),
!    2. sets the char-return point source to ETA x elutriation (cyclone),
!    3. integrates the char mass flux leaving the SIDE overflow drain,
!    4. logs time, elutriation, return, and overflow to recirc.csv.
!  Steps 3-4 let postprocess.py close a steady-state solid-carbon balance:
!    X_C = 1 - (overflow + (1-ETA)*elutriation) / (fresh char carbon fed).
!
!  >>> VERIFY (could not be tested without building) <<<
!    * OUTLET_BC / OVERFLOW_BC / RETURN_PS must match THIS deck. The
!      overflow BC index differs by case: 2D = 5, 3D = 7.
!    * Top outlet (y-normal): area AXZ, velocity V_s. Side overflow
!      (x-normal): area AYZ, velocity U_s. Swap if your faces differ.
!    * GLOBAL_ALL_SUM interface (DMP reduction).
!  ======================================================================

      SUBROUTINE USR1

      USE bc
      USE compar
      USE fldvar,      only: ROP_s, X_s, U_s, V_s
      USE functions,   only: funijk, im_of, jm_of, fluid_at,             &
                             is_on_mype_owns
      USE geometry,    only: axz, ayz
      USE mpi_utility, only: global_all_sum
      USE param1,      only: zero
      USE ps,          only: ps_massflow_s
      USE run,         only: time

      IMPLICIT NONE

! Configuration
!---------------------------------------------------------------------
      INTEGER, PARAMETER :: OUTLET_BC   = 2   ! top pressure-outflow BC
      INTEGER, PARAMETER :: OVERFLOW_BC = 5   ! side overflow drain (2D=5, 3D=7)
      INTEGER, PARAMETER :: RETURN_PS   = 3   ! char-return point source
      INTEGER, PARAMETER :: CHAR_PHASE  = 1   ! solids phase carrying char
      INTEGER, PARAMETER :: CHAR_SP     = 1   ! char species index in phase 1
      DOUBLE PRECISION, PARAMETER :: ETA = 0.90d0   ! cyclone efficiency
      DOUBLE PRECISION, PARAMETER :: T_START = 1.0d0
      DOUBLE PRECISION, PARAMETER :: LOG_DT  = 0.05d0

! Local variables
!---------------------------------------------------------------------
      INTEGER :: I, J, K, IJK, IJKM, lun
      DOUBLE PRECISION :: vface
      DOUBLE PRECISION :: elut_loc, elut_glb     ! top elutriation [kg/s]
      DOUBLE PRECISION :: ovfl_loc, ovfl_glb     ! side overflow   [kg/s]
      DOUBLE PRECISION, SAVE :: t_next = ZERO

!---------------------------------------------------------------------

      elut_loc = ZERO
      ovfl_loc = ZERO

! (1) Char mass flux out the TOP outlet (y-normal face: V_s, AXZ).
      J = BC_J_s(OUTLET_BC)
      DO K = BC_K_b(OUTLET_BC), BC_K_t(OUTLET_BC)
      DO I = BC_I_w(OUTLET_BC), BC_I_e(OUTLET_BC)
         IJK  = FUNIJK(I, J, K)
         IJKM = JM_OF(IJK)                       ! interior fluid cell below
         IF (.NOT. FLUID_AT(IJKM)) CYCLE
         IF (.NOT. IS_ON_myPE_owns(I, J-1, K)) CYCLE
         vface = V_s(IJKM, CHAR_PHASE)
         IF (vface <= ZERO) CYCLE                ! upward outflow only
         elut_loc = elut_loc                                             &
            + ROP_s(IJKM,CHAR_PHASE) * X_s(IJKM,CHAR_PHASE,CHAR_SP)      &
              * vface * AXZ(IJKM)
      ENDDO
      ENDDO

! (3) Char mass flux out the SIDE overflow (x-normal face: U_s, AYZ).
      I = BC_I_w(OVERFLOW_BC)
      DO K = BC_K_b(OVERFLOW_BC), BC_K_t(OVERFLOW_BC)
      DO J = BC_J_s(OVERFLOW_BC), BC_J_n(OVERFLOW_BC)
         IJK  = FUNIJK(I, J, K)
         IJKM = IM_OF(IJK)                        ! interior fluid cell to west
         IF (.NOT. FLUID_AT(IJKM)) CYCLE
         IF (.NOT. IS_ON_myPE_owns(I-1, J, K)) CYCLE
         vface = U_s(IJKM, CHAR_PHASE)
         IF (vface <= ZERO) CYCLE                ! eastward outflow only
         ovfl_loc = ovfl_loc                                             &
            + ROP_s(IJKM,CHAR_PHASE) * X_s(IJKM,CHAR_PHASE,CHAR_SP)      &
              * vface * AYZ(IJKM)
      ENDDO
      ENDDO

! Sum across MPI ranks (no-op in serial)
      CALL GLOBAL_ALL_SUM(elut_loc, elut_glb)
      CALL GLOBAL_ALL_SUM(ovfl_loc, ovfl_glb)

! (2) Set the return rate (kg/s). Guard against negatives and the warm-up.
      IF (TIME >= T_START) THEN
         PS_MASSFLOW_S(RETURN_PS, CHAR_PHASE) = MAX(ZERO, ETA*elut_glb)
      ELSE
         PS_MASSFLOW_S(RETURN_PS, CHAR_PHASE) = ZERO
      ENDIF

! (4) Log (I/O rank only, throttled) for postprocess.py
      IF (myPE == PE_IO .AND. TIME >= t_next) THEN
         OPEN(newunit=lun, file='recirc.csv', position='append', action='write')
         WRITE(lun,'(ES12.5,",",ES12.5,",",ES12.5,",",ES12.5)')          &
               TIME, elut_glb,                                           &
               PS_MASSFLOW_S(RETURN_PS, CHAR_PHASE), ovfl_glb
         CLOSE(lun)
         t_next = TIME + LOG_DT
      ENDIF

      RETURN
      END SUBROUTINE USR1
