! This function calculates the differential cross-section for secondary electron production, with secondary electron energy Es, primary electron energy E and ionization potential I

  FUNCTION DIFCS_SEC(Es,E,I)
    USE PRECISION
    IMPLICIT NONE
    REAL(RP) :: E, Es, E0, I, DIFCS_SEC
    E0 = 13.0 ! shape parameter for N2
!    I = 15.6
    DIFCS_SEC = ((1. + (Es/E0))**(-2.1))/(E0*tanh((E-I)/(2.*E0)))
  END FUNCTION DIFCS_SEC
