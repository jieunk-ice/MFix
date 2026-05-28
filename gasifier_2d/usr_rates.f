!  ======================================================================
!  usr_rates.f  -  reaction rates for the 2D char steam gasification case
!  ----------------------------------------------------------------------
!  MFiX calls this once per fluid cell (IJK) and expects RATES() filled
!  for every reaction declared in the @(RXNS) block of gasifier_2d.mfx.
!  The reaction names (Char_Gasification, Boudouard, Water_Gas_Shift) and
!  the species aliases (Char, H2O, CO, CO2, H2, ...) are auto-generated as
!  integer indices and are available here once the case is built.
!
!  >>> UNITS WARNING <<<
!  RATES() must be in your MFiX version's expected molar rate units
!  (SI builds: kmol . m^-3 . s^-1). Before trusting numbers, open the
!  bundled "silane_pyrolysis" tutorial's usr_rates.f and confirm the
!  convention, then rescale the pre-exponential factors A1..A3 below.
!
!  >>> KINETICS WARNING <<<
!  The pre-exponentials (A*) and activation energies (E*) here are
!  PLACEHOLDERS chosen only so the case runs. Replace them with kinetics
!  for YOUR char. Useful starting references:
!    - char-steam:  Barrio & Hustad (2001); Matsui et al. (1987)
!    - Boudouard:   Barrio et al. (2001)
!    - water-gas shift: Bustamante et al. (2004/2005)
!  ======================================================================

      SUBROUTINE USR_RATES(IJK, RATES)

      USE compar
      USE constant
      USE fldvar
      USE functions
      USE param1,    only: zero, small_number
      USE physprop
      USE run
      USE rxns
      USE toleranc

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJK
      DOUBLE PRECISION, DIMENSION(NO_OF_RXNS), INTENT(OUT) :: RATES

! Local variables
!---------------------------------------------------------------------
      DOUBLE PRECISION :: Ts, Tg          ! solids / gas temperature [K]
      DOUBLE PRECISION :: c_H2O, c_CO, c_CO2   ! gas molar conc [kmol/m3]
      DOUBLE PRECISION :: c_Char          ! carbon molar conc [kmol/m3]
      DOUBLE PRECISION :: k1, k2, k3      ! rate constants

! Universal gas constant consistent with kmol:  R = 8314.34 J/(kmol.K)
      DOUBLE PRECISION, PARAMETER :: Rg = 8314.34d0

! Placeholder Arrhenius parameters  (A = pre-exp, E = activation energy)
      DOUBLE PRECISION, PARAMETER :: A1 = 1.0d3,  E1 = 1.50d8  ! char + H2O
      DOUBLE PRECISION, PARAMETER :: A2 = 5.0d2,  E2 = 1.65d8  ! Boudouard
      DOUBLE PRECISION, PARAMETER :: A3 = 1.0d2,  E3 = 7.00d7  ! WGS

!---------------------------------------------------------------------

      RATES(:) = ZERO

! Temperatures (gasification rates are controlled by the solids temp)
      Tg = T_g(IJK)
      Ts = T_s(IJK,1)

! Gas-phase molar concentrations [kmol/m3] = rho_g * X_g / MW_g
      c_H2O = RO_g(IJK) * X_g(IJK,H2O) / MW_g(H2O)
      c_CO  = RO_g(IJK) * X_g(IJK,CO)  / MW_g(CO)
      c_CO2 = RO_g(IJK) * X_g(IJK,CO2) / MW_g(CO2)

! Solid carbon molar concentration [kmol/m3] = rho_bulk,s * X_char / MW_C
      c_Char = ROP_s(IJK,1) * X_s(IJK,1,Char) / MW_s(1,Char)

!---------------------------------------------------------------------
! (1) Char-steam gasification:  C + H2O --> CO + H2   (heterogeneous)
!     1st order in steam, tied to available carbon. Zero if no char/steam.
      IF (c_Char > SMALL_NUMBER .AND. c_H2O > SMALL_NUMBER) THEN
         k1 = A1 * EXP(-E1/(Rg*Ts))
         RATES(Char_Gasification) = k1 * c_H2O * c_Char
      ENDIF

!---------------------------------------------------------------------
! (2) Boudouard:  C + CO2 --> 2 CO   (heterogeneous)
      IF (c_Char > SMALL_NUMBER .AND. c_CO2 > SMALL_NUMBER) THEN
         k2 = A2 * EXP(-E2/(Rg*Ts))
         RATES(Boudouard) = k2 * c_CO2 * c_Char
      ENDIF

!---------------------------------------------------------------------
! (3) Water-gas shift:  CO + H2O --> CO2 + H2   (homogeneous gas phase)
!     Forward only here; add the reverse for true equilibrium behavior.
      IF (c_CO > SMALL_NUMBER .AND. c_H2O > SMALL_NUMBER) THEN
         k3 = A3 * EXP(-E3/(Rg*Tg))
         RATES(Water_Gas_Shift) = k3 * c_CO * c_H2O
      ENDIF

      RETURN
      END SUBROUTINE USR_RATES
