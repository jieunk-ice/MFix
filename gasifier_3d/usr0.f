!  ======================================================================
!  usr0.f  -  one-time setup for the cyclone char-return diagnostics
!             (shared by gasifier_2d.mfx and gasifier_3d.mfx)
!  ----------------------------------------------------------------------
!  Called once at the start of the run (requires CALL_USR = .True.).
!  Creates/!truncates the recirculation log written by usr1.f so the
!  closed loop can be inspected and plotted (see postprocess.py).
!  Only the I/O rank writes, so the file is correct under MPI (DMP).
!  ======================================================================

      SUBROUTINE USR0

      USE compar, only: myPE, PE_IO

      IMPLICIT NONE

      INTEGER :: lun

      IF (myPE /= PE_IO) RETURN

      OPEN(newunit=lun, file='recirc.csv', status='replace', action='write')
      WRITE(lun,'(A)') 'Time,elutriation_char_kg_s,return_char_kg_s'
      CLOSE(lun)

      RETURN
      END SUBROUTINE USR0
