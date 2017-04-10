   SUBROUTINE INTRP(xtab,ytab,xval,yval)
     USE PRECISION
     IMPLICIT NONE
     INTERFACE
        FUNCTION LOCATE(xtab, x)
          USE PRECISION
          IMPLICIT NONE
          REAL(RP), DIMENSION(:), INTENT(IN) :: xtab
          REAL(RP), INTENT(IN) :: x
          INTEGER :: LOCATE
        END FUNCTION LOCATE
     END INTERFACE
     REAL(RP), INTENT(IN), DIMENSION(:) :: xtab, ytab, xval
     REAL(RP), INTENT(OUT), DIMENSION(:) :: yval
     INTEGER :: ntab, nval, nv, l1, l2
     ntab = SIZE(xtab)
     nval = SIZE(xval)
     DO nv = 1, nval
        IF (xval(nv) <= xtab(1)) THEN
            yval(nv) = ytab(1)
        ELSE IF (xval(nv) >= xtab(ntab)) THEN
            yval(nv) = ytab(ntab)
        ELSE
           l1 = LOCATE(xtab,xval(nv))
           l2 = l1 + 1
           yval(nv)=ytab(l1)+(xval(nv)-xtab(l1))*(ytab(l2)-ytab(l1))/   &
                (xtab(l2)-xtab(l1))
        END IF
     END DO
     RETURN
   END SUBROUTINE INTRP
