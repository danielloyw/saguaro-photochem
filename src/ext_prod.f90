SUBROUTINE EXT_PROD

  USE PRECISION
  USE CONSTANTS
  USE GLOBAL_VARIABLES
  USE SUBS, ONLY : FIND_NAME
  
  IMPLICIT NONE
  INTEGER :: nx, nm, nl
  REAL(RP) :: sm, pext

  ALLOCATE(prext(nlev,0:nsp)) 
  prext(:,:) = zero
  IF(next > 0) THEN
  DO nx = 1, next
     nm = FIND_NAME(name_ext(nx), name)
     IF((nm <= 1).or.(nm >= nsp)) THEN
        WRITE(*,*) nm,name_ext(nm),' EXT PRODUCTION: NAME NOT FOUND, EXITING ...'
        STOP
     ENDIF
     
     pext = fext(nx)/hext(nx)
     DO nl = 1, nlev
        prext(nl,nm) = pext*EXP(-(z(nl)-zext(nx))/hext(nx)-EXP(-(z(nl)-zext(nx))/hext(nx)))
     END DO
     sm = zero
     DO nl = 2, nlev
        sm = sm + half*(prext(nl,nm)*(rz(nl)/rz(1))**2                                              &
             + prext(nl-1,nm)*(rz(nl-1)/rz(1))**2)*(rz(nl)-rz(nl-1))
     END DO
     DO nl = 1, nlev
        prext(nl,nm) = (fext(nx)/sm)*prext(nl,nm)
     END DO
  END DO
  END IF

  RETURN
END SUBROUTINE EXT_PROD
