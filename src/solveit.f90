! Implements the tridiagonal matrix algorithm.
! c'(1) = c(1)/b(1)
! c'(i) = c(i)/(b(i) - a(i)*c'(i-1))
! d'(1) = d(1)/b(1)
! d'(i) = (d(i)-a(i)*d'(i-1))/(b(i) - a(i)*c'(i-1))
! x(i) = d'(i) - c'(i)*x(i+1)
! x(n) = d'(n)
! dum is the "denominator": b(1) when i=1, b(i) - a(i)*c'(i-1) when i>1
! del holds d' until it is converted to x at the last step
! gam(n) is c'(n-1)

SUBROUTINE SOLVEIT (iopt, a, b, c, d, del, ierr)

  USE PRECISION
  USE CONSTANTS
  USE F95_Precision, ONLY: WP => DP
  USE LAPACK95, ONLY:GETRF, GETRS

  IMPLICIT none

  !
  !  .. External Variables
  !

  INTEGER, INTENT(IN) :: iopt ! optimization
  REAL(RP), INTENT(IN), DIMENSION(:,:,:) :: a
  REAL(RP), INTENT(IN), DIMENSION(:,:,:) :: b
  REAL(RP), INTENT(IN), DIMENSION(:,:,:) :: c
  REAL(RP), INTENT(IN), DIMENSION(:,:) :: d
  REAL(RP), INTENT(OUT), DIMENSION(:,:) :: del
  INTEGER, INTENT(OUT) :: ierr

  !
  !  .. Local Variables
  !

  INTEGER :: nmol, nlev
  INTEGER, DIMENSION(SIZE(a,1)) :: indx
  INTEGER :: nl, j, i, m
  REAL(RP), DIMENSION(SIZE(a,1),SIZE(a,1),SIZE(a,3)) :: gam
  REAL(RP), DIMENSION(SIZE(a,1),SIZE(a,3)) :: dp, d_sav, delp
  REAL(RP), DIMENSION(SIZE(a,1),SIZE(a,1)) :: dum
  REAL(RP), DIMENSION(SIZE(a,1)) :: sav
  REAL(RP) :: sm
  INTEGER :: info
  
  !
  !  .. Initialize
  !

  nmol=SIZE(a,1)
  nlev=SIZE(a,3)

  ierr = 0
  d_sav = d
  dum = zero
  gam = zero
  del = zero

  ! .. solve b(1)*del(1)=d(1)

  nl = 1
  dum(1:nmol,1:nmol) = b(1:nmol,1:nmol,nl)
  sav(1:nmol) = d(1:nmol,nl)

  CALL GETRF ( dum(1:nmol,1:nmol), indx(1:nmol), info) ! A = P*L*U: A (replaced by LU), P
  IF(info /= 0) THEN
     WRITE(*,"(' SOLVEIT, ERROR IN GETRF: info = ',I3,' nl = ',I4)") info, nl
     ierr = nl
     RETURN
  END IF
  CALL GETRS ( dum(1:nmol,1:nmol), indx(1:nmol), sav(1:nmol), 'N', info ) ! A*X = B: A, P, B (replaced with X)
  IF(info /= 0) THEN
     WRITE(*,"(' SOLVEIT, ERROR IN GETRS: info = ',I3,' nl = ',I4)") info, nl
     ierr = nl
     RETURN
  END IF

  del(1:nmol,nl) = sav(1:nmol)

  DO nl = 2, nlev

     ! .. solve dum*gam=c                        

     DO j = 1, nmol
        gam(j,j,nl) = c(j,j,nl-1)
        sav(1:nmol) = gam(1:nmol,j,nl)
        CALL GETRS ( dum(1:nmol,1:nmol), indx(1:nmol), sav(1:nmol), 'N', info )
        IF(info /= 0) THEN
           WRITE(*,"(' SOLVEIT, ERROR IN GETRS: info = ',I3,' nl = ',I4)") info, nl
           ierr = nl
           RETURN
        END IF
        gam(1:nmol,j,nl) = sav(1:nmol)
     END DO
     
     ! .. update dum = b(nl)-a(nl)*gam(nl) for i>1
     
     DO j = 1, nmol
        DO i = 1, nmol
           dum(i,j) = b(i,j,nl) - a(i,i,nl) * gam(i,j,nl)
        END DO
     END DO

     ! .. solve dum*del(nl)=d(nl)-a(nl)*del(nl-1)

     DO i = 1, nmol
        del(i,nl) = d(i,nl) - a(i,i,nl) * del(i,nl-1)
     END DO
     
     CALL GETRF (dum(1:nmol,1:nmol), indx(1:nmol), info )
     IF(info /= 0) THEN
        WRITE(*,"(' SOLVEIT, ERROR IN GETRF: info = ',I3,' nl = ',I4)") info, nl
        ierr = nl
        RETURN
     END IF
     sav(1:nmol) = del(1:nmol,nl)
     CALL GETRS (dum(1:nmol,1:nmol), indx(1:nmol), sav(1:nmol), 'N', info )
     IF(info /= 0) THEN
        WRITE(*,"(' SOLVEIT, ERROR IN GETRS: info = ',I3,' nl = ',I4)") info, nl
        ierr = nl
        RETURN
     END IF
     del(1:nmol,nl) = sav(1:nmol)
     
  END DO

! Convert d' to x

  DO nl = nlev-1, 1, -1
     
     ! .. del(nl) = del(nl) - gam(nl+1)*del(nl+1)
     
     DO i = 1, nmol
        sm = zero
        DO m = 1, nmol
           sm = sm + gam(i,m,nl+1) * del(m,nl+1)
        END DO
        del(i,nl) = del(i,nl) - sm
     END DO
     
  END DO

  ! ############################################################################
  !
  !  .. Improve Solution, Create error vector
  !
  ! ############################################################################

  IF(iopt >= 1) THEN

  nl = 1
  DO i = 1, nmol
     sm = zero
     DO j = 1, nmol
        sm = sm + b(i,j,nl)*del(j,nl)+c(i,j,nl)*del(j,nl+1)
     END DO
     dp(i,nl)=sm
  END DO

  DO nl = 2, nlev-1
  DO i = 1, nmol
     sm = zero
     DO j = 1, nmol
        sm = sm + a(i,j,nl)*del(j,nl-1)+b(i,j,nl)*del(j,nl)+c(i,j,nl)*del(j,nl+1)
     END DO
     dp(i,nl)=sm
  END DO
  END DO

  nl = nlev
  DO i = 1, nmol
     sm = zero
     DO j = 1, nmol
        sm = sm + a(i,j,nl)*del(j,nl-1)+b(i,j,nl)*del(j,nl)
     END DO
     dp(i,nl)=sm
  END DO

  dp = dp - d_sav

  ! .. solve A*del=dp

  dum = zero
  gam = zero
  delp = zero

  nl = 1
  dum(1:nmol,1:nmol) = b(1:nmol,1:nmol,nl)
  sav(1:nmol) = dp(1:nmol,nl)

  CALL GETRF ( dum(1:nmol,1:nmol), indx(1:nmol))
  CALL GETRS ( dum(1:nmol,1:nmol), indx(1:nmol), sav(1:nmol) )
  delp(1:nmol,nl) = sav(1:nmol)

  !                                                                        
  !  ***** frontwards                                                      
  !                                                                        

  DO nl = 2, nlev

     ! .. solve dum*gam=c                        

     DO j = 1, nmol
        gam(j,j,nl) = c(j,j,nl-1)
     END DO
     
     DO j = 1, nmol
        sav(1:nmol) = gam(1:nmol,j,nl)
        CALL GETRS ( dum(1:nmol,1:nmol), indx(1:nmol), sav(1:nmol) )
        gam(1:nmol,j,nl) = sav(1:nmol)
     END DO
     
     ! .. dum = b(nl)-a(nl)*gam(nl)
     
     DO j = 1, nmol
        DO i = 1, nmol
           dum(i,j) = b(i,j,nl) - a(i,i,nl) * gam(i,j,nl)
        END DO
     END DO

     ! .. solve dum*del(nl)=d(nl)-a(nl)*del(nl-1)

     DO i = 1, nmol
        delp(i,nl) = dp(i,nl) - a(i,i,nl) * delp(i,nl-1)
     END DO
     
     CALL GETRF (dum(1:nmol,1:nmol), indx(1:nmol) )
     sav(1:nmol) = delp(1:nmol,nl)
     CALL GETRS (dum(1:nmol,1:nmol), indx(1:nmol), sav(1:nmol) )
     delp(1:nmol,nl) = sav(1:nmol)
     
  END DO
  
  !                                                                        
  !  ***** backwards                                                      
  !                                                                        
  
  DO nl = nlev-1, 1, -1
     
     ! .. del(nl) = del(nl) - gam(nl+1)*del(nl+1)
     
     DO i = 1, nmol
        sm = zero
        DO m = 1, nmol
           sm = sm + gam(i,m,nl+1) * delp(m,nl+1)
        END DO
        delp(i,nl) = delp(i,nl) - sm
     END DO
     
  END DO

  del = del - delp

  END IF

!     nl = 1
!     DO nx = 1, neq
!        sm = zero
!        DO nx1 = 1, neq
!           sm = sm + fb(nx,nx1,nl)*del(nx1,nl)+fc(nx,nx1,nl)*del(nx1,nl+1)
!        END DO
!        fd_tst(nx,nl)=sm
!     END DO

!     DO nl = 2, nlev-1
!     DO nx = 1, neq
!        sm = zero
!        DO nx1 = 1, neq
!           sm = sm + fa(nx,nx1,nl)*del(nx1,nl-1)+fb(nx,nx1,nl)*del(nx1,nl)+fc(nx,nx1,nl)*del(nx1,nl+1)
!        END DO
!        fd_tst(nx,nl)=sm
!     END DO
!     END DO

!     nl = nlev
!     DO nx = 1, neq
!        sm = zero
!        DO nx1 = 1, neq
!           sm = sm + fa(nx,nx1,nl)*del(nx1,nl-1)+fb(nx,nx1,nl)*del(nx1,nl)
!        END DO
!        fd_tst(nx,nl)=sm
!     END DO

!     imax = MAXLOC(ABS(fd_tst-fd_sav))
!     atst = ABS(fd_tst(imax(1),imax(2))-fd_sav(imax(1),imax(2)))
!     rtst = atst/ABS(fd_sav(imax(1),imax(2)))
!     IF(((rtst > 1.0E-10_RP).and.(ABS(fd_sav(imax(1),imax(2))) > 1.0E-30_RP))) THEN
!        DO nl = 1, nlev
!           WRITE(*,"(I3,3(2X,ES16.9))") nl,fd_tst(imax(1),nl),fd_sav(imax(1),nl),ABS(fd_tst(imax(1),nl)-fd_sav(imax(1),nl))
!        END DO
!        WRITE(*,"(' MATRIX INVERSION INACCURATE: ',3ES11.3,' NL=',I3,' SPECIES=',A12)") &
!             fd_tst(imax(1),imax(2))-fd_sav(imax(1),imax(2)),fd_sav(imax(1),imax(2)), &
!             fd_tst(imax(1),imax(2)),imax(2),name(locp(imax(1)))
!        STOP 
!     END IF     
    
       
  RETURN
END SUBROUTINE SOLVEIT
