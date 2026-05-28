!  ======================================================================
!  usr_rates.f  -  reaction rates for the 2D char steam gasification case
!  ----------------------------------------------------------------------
!  MFiX calls this once per fluid cell (IJK) and expects RATES() filled
!  for every reaction declared in the @(RXNS) block of gasifier_2d.mfx:
!
!    1 Char_Gasification        Char + H2O --> CO + H2
!    2 Boudouard                Char + CO2 --> 2 CO
!    3 Methanation              Char + 2 H2 --> CH4
!    4 Water_Gas_Shift          CO + H2O --> CO2 + H2
!    5 Rev_Water_Gas_Shift      CO2 + H2 --> CO + H2O
!    6 Steam_Methane_Reforming  CH4 + H2O --> CO + 3 H2
!    7 Tar_Reforming            C6H6 + 6 H2O --> 6 CO + 9 H2
!
!  The reaction names above and the species aliases (Char, H2O, CO, CO2,
!  H2, CH4, Tar) are auto-generated as integer indices and are available
!  here once the case is built.
!
!  >>> UNITS WARNING <<<
!  RATES() must be in your MFiX version's expected molar rate units
!  (SI builds: kmol . m^-3 . s^-1). Confirm against the bundled
!  "silane_pyrolysis" tutorial's usr_rates.f, then rescale A1..A7.
!
!  >>> KINETICS WARNING <<<
!  The pre-exponentials (A*) and activation energies (E*) are PLACEHOLDERS
!  chosen only so the case runs. Replace with kinetics for YOUR char/tar:
!    - char-steam:  Barrio & Hustad (2001); Matsui et al. (1987)
!    - Boudouard:   Barrio et al. (2001)
!    - methanation / SMR: Jones & Lindstedt (1988); Xu & Froment (1989)
!    - water-gas shift: Bustamante et al. (2004/2005)
!    - tar (benzene) reforming: Jess (1996)
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
      DOUBLE PRECISION :: Ts, Tg                ! solids / gas temp [K]
      DOUBLE PRECISION :: c_H2O, c_CO, c_CO2    ! gas conc [kmol/m3]
      DOUBLE PRECISION :: c_H2, c_CH4, c_Tar
      DOUBLE PRECISION :: c_Char                ! solid carbon conc [kmol/m3]
      DOUBLE PRECISION :: k_wgs, Keq, wgs_net   ! water-gas-shift terms

! Universal gas constant consistent with kmol:  R = 8314.34 J/(kmol.K)
      DOUBLE PRECISION, PARAMETER :: Rg = 8314.34d0

! Placeholder Arrhenius parameters  (A = pre-exp, E = activation energy)
      DOUBLE PRECISION, PARAMETER :: A1 = 1.0d3,  E1 = 1.50d8  ! char + H2O
      DOUBLE PRECISION, PARAMETER :: A2 = 5.0d2,  E2 = 1.65d8  ! Boudouard
      DOUBLE PRECISION, PARAMETER :: A3 = 1.0d1,  E3 = 1.30d8  ! methanation
      DOUBLE PRECISION, PARAMETER :: A4 = 1.0d2,  E4 = 7.00d7  ! WGS (fwd)
      DOUBLE PRECISION, PARAMETER :: A6 = 5.0d1,  E6 = 1.25d8  ! SMR
      DOUBLE PRECISION, PARAMETER :: A7 = 1.0d2,  E7 = 1.00d8  ! tar reforming

!---------------------------------------------------------------------

      RATES(:) = ZERO

! Temperatures (char reactions follow the solids temperature)
      Tg = T_g(IJK)
      Ts = T_s(IJK,1)

! Gas-phase molar concentrations [kmol/m3] = rho_g * X_g / MW_g
      c_H2O = RO_g(IJK) * X_g(IJK,H2O) / MW_g(H2O)
      c_CO  = RO_g(IJK) * X_g(IJK,CO)  / MW_g(CO)
      c_CO2 = RO_g(IJK) * X_g(IJK,CO2) / MW_g(CO2)
      c_H2  = RO_g(IJK) * X_g(IJK,H2)  / MW_g(H2)
      c_CH4 = RO_g(IJK) * X_g(IJK,CH4) / MW_g(CH4)
      c_Tar = RO_g(IJK) * X_g(IJK,Tar) / MW_g(Tar)

! Solid carbon molar concentration [kmol/m3] = rho_bulk,s * X_char / MW_C
      c_Char = ROP_s(IJK,1) * X_s(IJK,1,Char) / MW_s(1,Char)

!---------------------------------------------------------------------
! (1) Char-steam gasification:  C + H2O --> CO + H2   (heterogeneous)
      IF (c_Char > SMALL_NUMBER .AND. c_H2O > SMALL_NUMBER) THEN
         RATES(Char_Gasification) = A1*EXP(-E1/(Rg*Ts)) * c_H2O * c_Char
      ENDIF

! (2) Boudouard:  C + CO2 --> 2 CO   (heterogeneous)
      IF (c_Char > SMALL_NUMBER .AND. c_CO2 > SMALL_NUMBER) THEN
         RATES(Boudouard) = A2*EXP(-E2/(Rg*Ts)) * c_CO2 * c_Char
      ENDIF

! (3) Methanation:  C + 2 H2 --> CH4   (heterogeneous)
      IF (c_Char > SMALL_NUMBER .AND. c_H2 > SMALL_NUMBER) THEN
         RATES(Methanation) = A3*EXP(-E3/(Rg*Ts)) * c_H2*c_H2 * c_Char
      ENDIF

!---------------------------------------------------------------------
! (4)/(5) Water-gas shift, reversible:  CO + H2O <--> CO2 + H2
!   Net rate uses the equilibrium constant Keq(T); the sign decides
!   whether the forward or the reverse reaction carries the rate.
!   Keq correlation (Moe, 1962):  Keq = exp(4577.8/T - 4.33)
      k_wgs = A4*EXP(-E4/(Rg*Tg))
      Keq   = EXP(4577.8d0/Tg - 4.33d0)
      wgs_net = k_wgs * (c_CO*c_H2O - c_CO2*c_H2/Keq)
      IF (wgs_net >= ZERO) THEN
         RATES(Water_Gas_Shift) = wgs_net
      ELSE
         RATES(Rev_Water_Gas_Shift) = -wgs_net
      ENDIF

!---------------------------------------------------------------------
! (6) Steam methane reforming:  CH4 + H2O --> CO + 3 H2   (homogeneous)
      IF (c_CH4 > SMALL_NUMBER .AND. c_H2O > SMALL_NUMBER) THEN
         RATES(Steam_Methane_Reforming) =                                &
            A6*EXP(-E6/(Rg*Tg)) * c_CH4 * c_H2O
      ENDIF

! (7) Tar (benzene) steam reforming:  C6H6 + 6 H2O --> 6 CO + 9 H2
      IF (c_Tar > SMALL_NUMBER .AND. c_H2O > SMALL_NUMBER) THEN
         RATES(Tar_Reforming) = A7*EXP(-E7/(Rg*Tg)) * c_Tar * c_H2O
      ENDIF

      RETURN
      END SUBROUTINE USR_RATES
