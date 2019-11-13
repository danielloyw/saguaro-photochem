  FUNCTION LOSS(E,NE,TE)
     USE PRECISION
     USE CONSTANTS
     IMPLICIT NONE
     REAL(RP) :: LOSS, E, NE, TE, Ee
     IF(E <= 0._RP) THEN
     WRITE(*,"(' LOSS: E = ',ES15.7)") E
        STOP
     END IF
     IF(Ne <= 0._RP) THEN
        WRITE(*,"(' LOSS: Ne = ',ES15.7)") Ne 
        STOP
     END IF
     Ee = 8.618E-5_RP*TE
     IF (Ee > E) THEN
        LOSS = zero
     ELSE 
     LOSS = (3.37e-12_RP/((E**0.94_RP)*(NE**0.03_RP)))*((E-Ee)/(E-0.53_RP*Ee))**2.36_RP
     END IF
!    IF (Ee > E) LOSS = zero !(3.37e-12/((E**0.94)*(NE**0.03)))*((Ee-E)/(Ee-0.53*E))**2.36  !zero
!   IF (Ee > E) write(*,*) 'LOSS = ',LOSS
  END FUNCTION LOSS
