SUBROUTINE HYDROST

  use types, only: wp => dp
  USE CONSTANTS
  USE GLOBAL_VARIABLES

  IMPLICIT NONE

  !
  !  .. Internal Variables
  !

  INTEGER :: nz,nm
  REAL(RP) :: alp, bet, fac
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: y, pnew, xnew
!  CHARACTER(len=1) :: iret

  ALLOCATE(y(n_z),pnew(n_z),xnew(n_z))

  !
  !  .. Recompute pressure
  !

  DO nz = 1, n_z
     y(nz) = RPLANET/rz(nz) - one
     prs(nz) = kB*tn(nz)*den(nz,0)
  END DO

  pnew(1) = kB*tn(1)*den(1,0)
  DO nz = 2, n_z
     fac = (GM*amu/kB/RPLANET)*(y(nz)-y(nz-1))*(mass(nz)/tn(nz)+mass(nz-1)/tn(nz-1))/two !GM/kT/R*dy*m
     pnew(nz) = pnew(nz-1)*EXP(fac)
  END DO

!  DO nz = 1, n_z
!     if (abs(pnew(nz)-prs(nz))/prs(nz) > 1.E-3_RP) THEN
!        WRITE(*,"(' HYDROST: nz, prs, pnew ',I6,2ES12.5)") nz, prs(nz), pnew(nz)
!        STOP
!     END IF
!  END DO

  !
  !  .. Map mass, temperature and mole fractions to new pressure grid
  !

  xnew = PMAP(mass, prs, pnew)
  mass = xnew
!  xnew = PMAP(tn, prs, pnew)
!  tn = xnew
!  xnew = PMAP(te, prs, pnew)
!  te = xnew
  DO nm = 1, n_sp
     xnew = PMAP(vmr(:,nm), prs, pnew)
     vmr(1:n_z,nm) = xnew(1:n_z)
  END DO

  DO nz = 1, n_z
     vmr(nz,1) = one-SUM(vmr(nz,2:n_sp))
  END DO


  prs(1:n_z) = pnew(1:n_z)
  den(1:n_z,0) = prs(1:n_z)/kB/tn(1:n_z)
  DO nm = 1, n_sp
     den(1:n_z,nm) = vmr(1:n_z,nm)*den(1:n_z,0)
  END DO
  

  DEALLOCATE(y,pnew,xnew)

  RETURN

CONTAINS

  FUNCTION PMAP(ytab,ptab,prs) ! interpolates ytab values at ptab pressures to new prs pressures
    use types, only: wp => dp
    USE CONSTANTS
    USE utils, ONLY : LOCATE
    REAL(RP), INTENT(IN), DIMENSION(:) :: ytab, ptab, prs
    REAL(RP), DIMENSION(SIZE(prs)) :: PMAP
    REAL(RP), ALLOCATABLE, DIMENSION(:) :: plog
    REAL(RP) :: rr, pn
    INTEGER :: i1, i2, n
    
    ALLOCATE(plog(SIZE(ptab)))
    plog = LOG(ptab)
    DO n = 1, SIZE(prs)
       pn = LOG(prs(n))
       if (prs(n) >= ptab(1)) THEN
          PMAP(n) = ytab(1)
       ELSE if (prs(n) <= ptab(n_z)) THEN
          PMAP(n) = ytab(n_z)
       ELSE
          i1 = LOCATE(pn,plog)
          i2 = i1 + 1
          rr = (pn-plog(i1))/(plog(i2)-plog(i1))
          PMAP(n) = ytab(i1) + rr*(ytab(i2)-ytab(i1))
       END IF
    END DO
    DEALLOCATE(plog)

  RETURN
  END FUNCTION PMAP

END SUBROUTINE HYDROST
