!  ======================================================================
!  usr_rates.f  -  reaction rates for the char steam gasification cases
!                  (shared by gasifier_2d.mfx and gasifier_3d.mfx)
!  ----------------------------------------------------------------------
!  Reactions (referenced by NAME, so @(RXNS) order is free):
!    Drying                   Moisture --> H2O
!    Pyrolysis                Biomass --> char + CO + CO2 + CH4 + H2 + H2O + Tar
!    Char_Gasification        Char + H2O --> CO + H2
!    Boudouard                Char + CO2 --> 2 CO
!    Methanation              Char + 2 H2 --> CH4
!    Water_Gas_Shift          CO + H2O --> CO2 + H2
!    Rev_Water_Gas_Shift      CO2 + H2 --> CO + H2O
!    Steam_Methane_Reforming  CH4 + H2O --> CO + 3 H2
!    Tar_Reforming            C6H6 + 6 H2O --> 6 CO + 9 H2
!    Char_Combustion          Char + O2 --> CO2          (autothermal)
!    CO_Combustion            CO + 0.5 O2 --> CO2        (autothermal)
!    H2_Combustion            H2 + 0.5 O2 --> H2O        (autothermal)
!    CH4_Combustion           CH4 + 2 O2 --> CO2 + 2 H2O (autothermal)
!
!  RATE FORMS (literature-based):
!    - Devolatilization (drying, pyrolysis) first order in the solid
!      species (Chan et al. 1985).
!    - Combustion mass-action (Westbrook & Dryer 1981; Smith 1982).
!    - Char-steam and Boudouard use Langmuir-Hinshelwood (LH) forms, which
!      capture the H2 / CO product inhibition seen in char gasification
!      (Barrio & Hustad 2001; compiled in Gomez-Barea & Leckner, Prog.
!      Energy Combust. Sci. 2010). The LH expression gives a specific
!      reactivity [1/s] on a carbon basis; multiplying by the local carbon
!      molar concentration gives the volumetric molar rate.
!    - Gas-phase reactions use mass-action kinetics (Arrhenius x conc.):
!        WGS (fwd) Bustamante et al. 2005; reverse via Keq(T) (Moe 1962)
!        SMR       Jones & Lindstedt 1988
!        Tar       Jess 1996 (benzene steam reforming)
!
!  >>> !!! RECALL / TODO: REPLACE WITH EXPERIMENTAL KINETICS !!! <<<
!    The A* and E* below are LITERATURE values (sources cited per line),
!    used as a stand-in so the model runs with physically reasonable
!    behaviour. They are NOT fitted to your char/tar. Before any
!    quantitative use, refit the pre-exponentials (and ideally the LH
!    inhibition constants) to your own gasification/TGA data and update
!    this block. Activation energies are the better-established part;
!    pre-exponentials vary by orders of magnitude between chars.
!
!  >>> UNITS <<<
!    RATES() must be in your MFiX build's molar-rate units (SI: kmol/m3/s).
!    Partial pressures here are in bar; concentrations in kmol/m3. Gas-phase
!    pre-exponentials below were converted from the papers' (mol,cm3,s) to
!    SI (kmol,m3,s) by x1e-3 for a bimolecular rate - re-derive if you
!    change the rate order. Confirm the overall convention against the
!    bundled silane_pyrolysis tutorial.
!  ======================================================================

      SUBROUTINE USR_RATES(IJK, RATES)

      USE compar
      USE constant
      USE fldvar
      USE functions
      USE param1,    only: zero, one, small_number
      USE physprop
      USE run
      USE rxns
      USE toleranc

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJK
      DOUBLE PRECISION, DIMENSION(NO_OF_RXNS), INTENT(OUT) :: RATES

! Local variables
!---------------------------------------------------------------------
      INTEGER :: N
      DOUBLE PRECISION :: Ts, Tg                ! solids / gas temp [K]
      DOUBLE PRECISION :: mw_mix                ! gas mixture MW [kg/kmol]
      DOUBLE PRECISION :: p_H2O, p_CO2, p_CO, p_H2   ! partial pressure [bar]
      DOUBLE PRECISION :: c_H2O, c_CO, c_CO2, c_H2   ! conc [kmol/m3]
      DOUBLE PRECISION :: c_CH4, c_Tar, c_O2
      DOUBLE PRECISION :: c_Char                ! solid carbon [kmol/m3]
      DOUBLE PRECISION :: c_Bio, c_Moist        ! biomass / moisture [kmol/m3]
      DOUBLE PRECISION :: r1, r2                ! LH reactivities [1/s]
      DOUBLE PRECISION :: k_wgs, Keq, wgs_net

! Universal gas constant consistent with kmol:  R = 8314.34 J/(kmol.K)
      DOUBLE PRECISION, PARAMETER :: Rg = 8314.34d0

! LITERATURE kinetics (see RECALL banner above). E in J/kmol; the LH
! reactivities r1,r2 come out in 1/s with p in bar; gas A's are in SI
! (kmol,m3,s). Inhibition terms use E=0 (adsorption ~ T-independent).

! --- Char-steam LH (Barrio & Hustad, Energy & Fuels 15 (2001) 1109):
!       r1 = k1*pH2O / (1 + k2*pH2 + k3*pH2O)      [1/s], p[bar]
!     E1 = 237 kJ/mol (birch char). Inhibition constants approximate -
!     refine from Barrio's tables for your char.
      DOUBLE PRECISION, PARAMETER :: A1=2.0d5,  E1=2.37d8  ! main [1/(bar.s)]
      DOUBLE PRECISION, PARAMETER :: A2=1.5d-2, E2=0.0d0   ! H2 inhibition [1/bar]
      DOUBLE PRECISION, PARAMETER :: A3=3.0d-2, E3=0.0d0   ! H2O term [1/bar]

! --- Boudouard LH (Barrio, Hustad et al., Prog. Thermochem. Biomass
!     Conversion 2001):  r2 = k4*pCO2 / (1 + k5*pCO + k6*pCO2)  [1/s]
!     E4 = 232 kJ/mol (birch char).
      DOUBLE PRECISION, PARAMETER :: A4=1.0d5,  E4=2.32d8  ! main [1/(bar.s)]
      DOUBLE PRECISION, PARAMETER :: A5=6.0d-2, E5=0.0d0   ! CO inhibition [1/bar]
      DOUBLE PRECISION, PARAMETER :: A6=2.0d-2, E6=0.0d0   ! CO2 term [1/bar]

! --- Methanation  C + 2H2 -> CH4  (slow; representative, E~150 kJ/mol;
!     e.g. Biba et al., Ind. Eng. Chem. Process Des. Dev. 17 (1978) 92):
      DOUBLE PRECISION, PARAMETER :: A7=1.0d1,  E7=1.50d8

! --- WGS forward (Bustamante et al., AIChE J. 51 (2005) 1440; homogeneous):
!     E = 288 kJ/mol; A = 2.34e10 cm3/mol/s -> 2.34e7 m3/kmol/s (x1e-3).
      DOUBLE PRECISION, PARAMETER :: A8=2.34d7, E8=2.88d8

! --- Steam methane reforming (Jones & Lindstedt, Combust. Flame 73 (1988)
!     233): E = 125 kJ/mol; A = 3e8 cm3/mol/s -> 3e5 m3/kmol/s (x1e-3).
      DOUBLE PRECISION, PARAMETER :: A9=3.0d5,  E9=1.25d8

! --- Tar (benzene) steam reforming (Jess, Chem. Eng. Process. 35 (1996)
!     487): E ~ 200 kJ/mol; A representative (refit to your tar).
      DOUBLE PRECISION, PARAMETER :: A10=1.0d4, E10=2.00d8

! --- Biomass drying  Moisture -> H2O  (Chan, Kelbon, Krieger, Fuel 64
!     (1985) 1505): k = 5.13e6 exp(-87.9 kJ/mol /RT) [1/s].
      DOUBLE PRECISION, PARAMETER :: A11=5.13d6, E11=8.79d7

! --- Biomass pyrolysis  Biomass -> char + volatiles (single step; Chan
!     et al. 1985 / Font et al.): A ~ 1e7 1/s, E ~ 140 kJ/mol.
      DOUBLE PRECISION, PARAMETER :: A12=1.0d7,  E12=1.40d8

! --- Combustion / partial oxidation (autothermal heat source). Gas A's in
!     SI (kmol,m3,s); refit/replace as for the other reactions.
!   Char + O2 -> CO2 (char combustion; Smith, 19th Symp. Combust. 1982):
      DOUBLE PRECISION, PARAMETER :: A13=5.0d2,  E13=1.40d8
!   CO + 0.5 O2 -> CO2 (Westbrook & Dryer, Combust. Sci. Tech. 1981):
      DOUBLE PRECISION, PARAMETER :: A14=1.0d10, E14=1.67d8
!   H2 + 0.5 O2 -> H2O (representative fast oxidation):
      DOUBLE PRECISION, PARAMETER :: A15=1.0d11, E15=1.10d8
!   CH4 + 2 O2 -> CO2 + 2 H2O (Westbrook & Dryer 1981): E ~ 202 kJ/mol:
      DOUBLE PRECISION, PARAMETER :: A16=2.0d8,  E16=2.02d8

!---------------------------------------------------------------------

      RATES(:) = ZERO

! Temperatures (char reactions follow the solids temperature)
      Tg = T_g(IJK)
      Ts = T_s(IJK,1)

! Gas mixture molecular weight (self-contained: no reliance on MW_MIX)
      mw_mix = ZERO
      DO N = 1, NMAX(0)
         mw_mix = mw_mix + X_g(IJK,N)/MW_g(N)
      ENDDO
      mw_mix = ONE / MAX(mw_mix, SMALL_NUMBER)

! Partial pressures [bar]:  p_i = y_i * P  ,  y_i = X_i * MW_mix / MW_i
      p_H2O = X_g(IJK,H2O)*mw_mix/MW_g(H2O) * P_g(IJK) / 1.0d5
      p_CO2 = X_g(IJK,CO2)*mw_mix/MW_g(CO2) * P_g(IJK) / 1.0d5
      p_CO  = X_g(IJK,CO) *mw_mix/MW_g(CO)  * P_g(IJK) / 1.0d5
      p_H2  = X_g(IJK,H2) *mw_mix/MW_g(H2)  * P_g(IJK) / 1.0d5

! Molar concentrations [kmol/m3] = rho_g * X_g / MW_g
      c_H2O = RO_g(IJK) * X_g(IJK,H2O) / MW_g(H2O)
      c_CO  = RO_g(IJK) * X_g(IJK,CO)  / MW_g(CO)
      c_CO2 = RO_g(IJK) * X_g(IJK,CO2) / MW_g(CO2)
      c_H2  = RO_g(IJK) * X_g(IJK,H2)  / MW_g(H2)
      c_CH4 = RO_g(IJK) * X_g(IJK,CH4) / MW_g(CH4)
      c_Tar = RO_g(IJK) * X_g(IJK,Tar) / MW_g(Tar)
      c_O2  = RO_g(IJK) * X_g(IJK,O2)  / MW_g(O2)

! Solid-phase molar concentrations [kmol/m3]
      c_Char  = ROP_s(IJK,1) * X_s(IJK,1,Char)     / MW_s(1,Char)
      c_Bio   = ROP_s(IJK,1) * X_s(IJK,1,Biomass)  / MW_s(1,Biomass)
      c_Moist = ROP_s(IJK,1) * X_s(IJK,1,Moisture) / MW_s(1,Moisture)

!---------------------------------------------------------------------
! Devolatilization (first order, solids temp)
!   Drying:    Moisture -> H2O
      IF (c_Moist > SMALL_NUMBER) THEN
         RATES(Drying) = A11*EXP(-E11/(Rg*Ts)) * c_Moist
      ENDIF
!   Pyrolysis: Biomass -> char + volatiles
      IF (c_Bio > SMALL_NUMBER) THEN
         RATES(Pyrolysis) = A12*EXP(-E12/(Rg*Ts)) * c_Bio
      ENDIF

!---------------------------------------------------------------------
! (1) Char-steam gasification  (Langmuir-Hinshelwood, solids temp)
      IF (c_Char > SMALL_NUMBER .AND. p_H2O > SMALL_NUMBER) THEN
         r1 = (A1*EXP(-E1/(Rg*Ts)) * p_H2O) /                            &
              (ONE + A2*EXP(-E2/(Rg*Ts))*p_H2                            &
                   + A3*EXP(-E3/(Rg*Ts))*p_H2O)
         RATES(Char_Gasification) = r1 * c_Char
      ENDIF

! (2) Boudouard  (Langmuir-Hinshelwood, solids temp)
      IF (c_Char > SMALL_NUMBER .AND. p_CO2 > SMALL_NUMBER) THEN
         r2 = (A4*EXP(-E4/(Rg*Ts)) * p_CO2) /                            &
              (ONE + A5*EXP(-E5/(Rg*Ts))*p_CO                            &
                   + A6*EXP(-E6/(Rg*Ts))*p_CO2)
         RATES(Boudouard) = r2 * c_Char
      ENDIF

! (3) Methanation  C + 2 H2 -> CH4
      IF (c_Char > SMALL_NUMBER .AND. c_H2 > SMALL_NUMBER) THEN
         RATES(Methanation) = A7*EXP(-E7/(Rg*Ts)) * c_H2*c_H2 * c_Char
      ENDIF

!---------------------------------------------------------------------
! (4)/(5) Water-gas shift, reversible:  CO + H2O <--> CO2 + H2
!   Net rate via Keq(T); its sign selects forward vs reverse reaction.
!   Keq (Moe 1962):  Keq = exp(4577.8/T - 4.33)
      k_wgs = A8*EXP(-E8/(Rg*Tg))
      Keq   = EXP(4577.8d0/Tg - 4.33d0)
      wgs_net = k_wgs * (c_CO*c_H2O - c_CO2*c_H2/Keq)
      IF (wgs_net >= ZERO) THEN
         RATES(Water_Gas_Shift) = wgs_net
      ELSE
         RATES(Rev_Water_Gas_Shift) = -wgs_net
      ENDIF

!---------------------------------------------------------------------
! (6) Steam methane reforming:  CH4 + H2O -> CO + 3 H2
      IF (c_CH4 > SMALL_NUMBER .AND. c_H2O > SMALL_NUMBER) THEN
         RATES(Steam_Methane_Reforming) =                               &
            A9*EXP(-E9/(Rg*Tg)) * c_CH4 * c_H2O
      ENDIF

! (7) Tar (benzene) steam reforming:  C6H6 + 6 H2O -> 6 CO + 9 H2
      IF (c_Tar > SMALL_NUMBER .AND. c_H2O > SMALL_NUMBER) THEN
         RATES(Tar_Reforming) = A10*EXP(-E10/(Rg*Tg)) * c_Tar * c_H2O
      ENDIF

!---------------------------------------------------------------------
! Combustion / partial oxidation (autothermal heat source). Only active
! where O2 is present (near the inlet); consumed fast.
      IF (c_O2 > SMALL_NUMBER) THEN
!   Char + O2 -> CO2   (heterogeneous, solids temp)
         IF (c_Char > SMALL_NUMBER)                                      &
            RATES(Char_Combustion) = A13*EXP(-E13/(Rg*Ts)) * c_O2 * c_Char
!   CO + 0.5 O2 -> CO2
         IF (c_CO > SMALL_NUMBER)                                        &
            RATES(CO_Combustion)  = A14*EXP(-E14/(Rg*Tg)) * c_CO  * c_O2
!   H2 + 0.5 O2 -> H2O
         IF (c_H2 > SMALL_NUMBER)                                        &
            RATES(H2_Combustion)  = A15*EXP(-E15/(Rg*Tg)) * c_H2  * c_O2
!   CH4 + 2 O2 -> CO2 + 2 H2O
         IF (c_CH4 > SMALL_NUMBER)                                       &
            RATES(CH4_Combustion) = A16*EXP(-E16/(Rg*Tg)) * c_CH4 * c_O2
      ENDIF

      RETURN
      END SUBROUTINE USR_RATES
