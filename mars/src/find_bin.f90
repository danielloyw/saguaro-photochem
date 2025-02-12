! This function finds the index of the bin containing Etest in an ascending vector elctreV with bin size elctDeV.

  FUNCTION FIND_BIN(elctreV,elctDeV,Etest)
    USE PRECISION
    IMPLICIT NONE
    REAL(RP), INTENT(IN), DIMENSION(:) :: elctreV, elctDeV
    REAL(RP), INTENT(IN) :: Etest
    INTEGER :: FIND_BIN
    INTEGER :: i, nelb
    nelb = SIZE(elctreV)
    FIND_BIN = 1
    DO i = 1, nelb-1
       IF((elctreV(i)-0.5_RP*elctDeV(i) <= Etest) .and.   &
            (Etest <= elctreV(i+1)-0.5_RP*elctDeV(i+1)) ) THEN 
          FIND_BIN = i 
          EXIT
       ENDIF
    ENDDO
    IF(Etest > elctreV(nelb-1)+0.5_RP*elctDeV(nelb-1)) FIND_BIN=nelb
!    IF(FIND_BIN == 0) THEN
!       WRITE(*,*) ' ERROR In FIND_BIN, SET TO 1'
!       FIND_BIN = 1
!    END IF


!   DO ne = 1, nelb
!       egrid(ne) = elctreV(ne)-half*elctDeV(ne)
!   END DO
!   FIND_BIN = LOCATE(egrid,Etest)
!

  END FUNCTION FIND_BIN
