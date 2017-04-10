  FUNCTION DIFCS_SEC(Es,E,I)
    USE PRECISION
    IMPLICIT NONE
    REAL(RP) :: E, Es, E0, I, DIFCS_SEC
    E0 = 13.0
!    I = 15.6
    DIFCS_SEC = ((1. + (Es/E0))**(-2.1))/(E0*tanh((E-I)/(2.*E0)))
  END FUNCTION DIFCS_SEC
