SUBROUTINE SOLVEIT (iopt, a, b, c, g, del, ierr)

  USE PRECISION
  USE CONSTANTS
  USE F95_Precision, ONLY: WP => DP
  USE LAPACK95, ONLY:GETRF, GETRS

  IMPLICIT none

  !
  !  .. External Variables
  !

  INTEGER, INTENT(IN) :: iopt
  REAL(RP), INTENT(IN), DIMENSION(:,:,:) :: a
  REAL(RP), INTENT(IN), DIMENSION(:,:,:) :: b
  REAL(RP), INTENT(IN), DIMENSION(:,:,:) :: c
  REAL(RP), INTENT(IN), DIMENSION(:,:) :: g
  REAL(RP), INTENT(OUT), DIMENSION(:,:) :: del
  INTEGER, INTENT(OUT) :: ierr

  !
  !  .. Local Variables
  !

  INTEGER :: nmol, nlev
  INTEGER, DIMENSION(SIZE(a,1)) :: indx
  INTEGER :: k, j, i, m 
  REAL(RP), DIMENSION(SIZE(a,1),SIZE(a,1),SIZE(a,3)) :: gam
  REAL(RP), DIMENSION(SIZE(a,1),SIZE(a,3)) :: gp, g_sav, delp
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
  g_sav = g
  dum = zero
  gam = zero
  del = zero

  ! .. solve b*del=g                        

  k = 1
  dum(1:nmol,1:nmol) = b(1:nmol,1:nmol,k)
  sav(1:nmol) = g(1:nmol,k)

  CALL GETRF ( dum(1:nmol,1:nmol), indx(1:nmol), info)
  IF(info /= 0) THEN
     WRITE(*,"(' SOLVEIT, ERROR IN GETRF: info = ',I3,' nl = ',I4)") info, k
     ierr = k
     RETURN
  END IF
  CALL GETRS ( dum(1:nmol,1:nmol), indx(1:nmol), sav(1:nmol), 'N', info )
  IF(info /= 0) THEN
     WRITE(*,"(' SOLVEIT, ERROR IN GETRS: info = ',I3,' nl = ',I4)") info, k
     ierr = k
     RETURN
  END IF

  del(1:nmol,k) = sav(1:nmol)

  !                                                                        
  !  ***** frontwards                                                      
  !                                                                        

  DO k = 2, nlev

     ! .. solve dum*gam=c                        

     DO j = 1, nmol
        gam(j,j,k) = c(j,j,k-1)
     END DO
     
     DO j = 1, nmol
        sav(1:nmol) = gam(1:nmol,j,k)
        CALL GETRS ( dum(1:nmol,1:nmol), indx(1:nmol), sav(1:nmol), 'N', info )
        IF(info /= 0) THEN
           WRITE(*,"(' SOLVEIT, ERROR IN GETRS: info = ',I3,' nl = ',I4)") info, k
           ierr = k
           RETURN
        END IF
        gam(1:nmol,j,k) = sav(1:nmol)
     END DO
     
     ! .. dum = b(k)-a(k)*gam(k)               
     
     DO j = 1, nmol
        DO i = 1, nmol
           dum(i,j) = b(i,j,k) - a(i,i,k) * gam(i,j,k)
        END DO
     END DO

     ! .. solve dum*del(k)=g(k)-a(k)*del(k-1)   

     DO i = 1, nmol
        del(i,k) = g(i,k) - a(i,i,k) * del(i,k-1)
     END DO
     
     CALL GETRF (dum(1:nmol,1:nmol), indx(1:nmol), info )
     IF(info /= 0) THEN
        WRITE(*,"(' SOLVEIT, ERROR IN GETRF: info = ',I3,' nl = ',I4)") info, k
        ierr = k
        RETURN
     END IF
     sav(1:nmol) = del(1:nmol,k)
     CALL GETRS (dum(1:nmol,1:nmol), indx(1:nmol), sav(1:nmol), 'N', info )
     IF(info /= 0) THEN
        WRITE(*,"(' SOLVEIT, ERROR IN GETRS: info = ',I3,' nl = ',I4)") info, k
        ierr = k
        RETURN
     END IF
     del(1:nmol,k) = sav(1:nmol)
     
  END DO
  
  !                                                                        
  !  ***** backwards                                                      
  !                                                                        
  
  DO k = nlev-1, 1, -1
     
     ! .. del(k) = del(k) - gam(k+1)*del(k+1)   
     
     DO i = 1, nmol
        sm = zero
        DO m = 1, nmol
           sm = sm + gam(i,m,k+1) * del(m,k+1)
        END DO
        del(i,k) = del(i,k) - sm
     END DO
     
  END DO

  ! ############################################################################
  !
  !  .. Improve Solution, Create error vector
  !
  ! ############################################################################

  IF(iopt >= 1) THEN

  k = 1
  DO i = 1, nmol
     sm = zero
     DO j = 1, nmol
        sm = sm + b(i,j,k)*del(j,k)+c(i,j,k)*del(j,k+1)
     END DO
     gp(i,k)=sm
  END DO

  DO k = 2, nlev-1
  DO i = 1, nmol
     sm = zero
     DO j = 1, nmol
        sm = sm + a(i,j,k)*del(j,k-1)+b(i,j,k)*del(j,k)+c(i,j,k)*del(j,k+1)
     END DO
     gp(i,k)=sm
  END DO
  END DO

  k = nlev
  DO i = 1, nmol
     sm = zero
     DO j = 1, nmol
        sm = sm + a(i,j,k)*del(j,k-1)+b(i,j,k)*del(j,k)
     END DO
     gp(i,k)=sm
  END DO

  gp = gp - g_sav

  ! .. solve A*del=gp                        

  dum = zero
  gam = zero
  delp = zero

  k = 1
  dum(1:nmol,1:nmol) = b(1:nmol,1:nmol,k)
  sav(1:nmol) = gp(1:nmol,k)

  CALL GETRF ( dum(1:nmol,1:nmol), indx(1:nmol))
  CALL GETRS ( dum(1:nmol,1:nmol), indx(1:nmol), sav(1:nmol) )
  delp(1:nmol,k) = sav(1:nmol)

  !                                                                        
  !  ***** frontwards                                                      
  !                                                                        

  DO k = 2, nlev

     ! .. solve dum*gam=c                        

     DO j = 1, nmol
        gam(j,j,k) = c(j,j,k-1)
     END DO
     
     DO j = 1, nmol
        sav(1:nmol) = gam(1:nmol,j,k)
        CALL GETRS ( dum(1:nmol,1:nmol), indx(1:nmol), sav(1:nmol) )
        gam(1:nmol,j,k) = sav(1:nmol)
     END DO
     
     ! .. dum = b(k)-a(k)*gam(k)               
     
     DO j = 1, nmol
        DO i = 1, nmol
           dum(i,j) = b(i,j,k) - a(i,i,k) * gam(i,j,k)
        END DO
     END DO

     ! .. solve dum*del(k)=g(k)-a(k)*del(k-1)   

     DO i = 1, nmol
        delp(i,k) = gp(i,k) - a(i,i,k) * delp(i,k-1)
     END DO
     
     CALL GETRF (dum(1:nmol,1:nmol), indx(1:nmol) )
     sav(1:nmol) = delp(1:nmol,k)
     CALL GETRS (dum(1:nmol,1:nmol), indx(1:nmol), sav(1:nmol) )
     delp(1:nmol,k) = sav(1:nmol)
     
  END DO
  
  !                                                                        
  !  ***** backwards                                                      
  !                                                                        
  
  DO k = nlev-1, 1, -1
     
     ! .. del(k) = del(k) - gam(k+1)*del(k+1)   
     
     DO i = 1, nmol
        sm = zero
        DO m = 1, nmol
           sm = sm + gam(i,m,k+1) * delp(m,k+1)
        END DO
        delp(i,k) = delp(i,k) - sm
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
