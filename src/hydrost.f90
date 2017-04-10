SUBROUTINE HYDROST

  USE PRECISION
  USE CONSTANTS
  USE GLOBAL_VARIABLES

  IMPLICIT NONE

  !
  !  .. Internal Variables
  !

  INTEGER :: nz,nm
  REAL(RP) :: alp, bet, fac
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: y, pnew,xnew
!  CHARACTER(len=1) :: iret

  ALLOCATE(y(nlev),pnew(nlev),xnew(nlev))

  !
  !  .. Recompute pressure
  !

  DO nz = 1, nlev
     y(nz) = RPLANET/rz(nz) - one
  END DO

  

  pnew(1) = rkb*tn(1)*den(1,0)
  DO nz = 2, nlev
     alp = mass(nz-1)/tn(nz-1)
     bet = (mass(nz)/tn(nz) - alp)
     fac = (GM*amu/rkb/RPLANET)*(y(nz)-y(nz-1))*(mass(nz)/tn(nz)+mass(nz-1)/tn(nz-1))/two
!     fac = (GM*amu/rkb/RPLANET)*(alp*(y(nz)-y(nz-1))+half*bet*(y(nz)+y(nz-1))) !GM/kT/R*dy*m
     pnew(nz) = pnew(nz-1)*EXP(fac)
  END DO

!  DO nz = 1, nlev
!     WRITE(*,"(F8.0,2X,F8.3,2(2X,ES11.3))") 1.E-5*z(nz),mass(nz),prs(nz),pnew(nz)
!  END DO
!  READ(*,"(A)") iret

  !
  !  .. Map mass, temperature and mole fractions to new pressure grid
  !

!  xnew = PMAP(mass, prs, pnew)
!  mass = xnew
!  xnew = PMAP(tn, prs, pnew)
!  tn = xnew
!  xnew = PMAP(te, prs, pnew)
!  te = xnew
!  DO nm = 1, nsp
!     xnew = PMAP(xmol(:,nm), prs, pnew)
!     xmol(1:nlev,nm) = xnew(1:nlev)
!  END DO

!  DO nz = 1, nlev
!     xmol(nz,1) = one-SUM(xmol(nz,2:nsp))
!  END DO

  prs(1:nlev) = pnew(1:nlev)
!  den(1:nlev,0) = prs(1:nlev)/rkb/tn(1:nlev)
!  DO nm = 1, nsp
!     den(1:nlev,nm) = xmol(1:nlev,nm)*den(1:nlev,0)
!  END DO

  DEALLOCATE(y,pnew,xnew)

  RETURN

CONTAINS

  FUNCTION PMAP(ytab,ptab,prs)
    USE PRECISION
    USE CONSTANTS
    USE SUBS, ONLY : LOCATE
    REAL(RP), INTENT(IN), DIMENSION(:) :: ytab, ptab, prs
    REAL(RP), DIMENSION(SIZE(prs)) :: PMAP
    REAL(RP), ALLOCATABLE, DIMENSION(:) :: plog
    REAL(RP) :: rr, pn
    INTEGER :: i1, i2, n
    
    ALLOCATE(plog(SIZE(ptab)))
    plog = LOG(ptab)
    DO n = 1, SIZE(prs)
       pn = LOG(prs(n))
       IF(prs(n) >= ptab(1)) THEN
          PMAP(n) = ytab(1)
       ELSE IF(prs(n) <= ptab(nlev)) THEN
          PMAP(n) = ytab(nlev)
       ELSE
          i1 = LOCATE(plog,pn)
          i2 = i1 + 1
          rr = (pn-plog(i1))/(plog(i2)-plog(i1))
          PMAP(n) = ytab(i1) + rr*(ytab(i2)-ytab(i1))
       END IF
    END DO
    DEALLOCATE(plog)

  RETURN
  END FUNCTION PMAP

END SUBROUTINE HYDROST
