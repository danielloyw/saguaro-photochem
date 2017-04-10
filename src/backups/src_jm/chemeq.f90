SUBROUTINE CHEMEQ 

  USE PRECISION
  USE CONSTANTS
  USE GLOBAL_VARIABLES
  USE F95_Precision, ONLY: WP => DP
  USE LAPACK95, ONLY:GESV
  
  !
  !  .. Local Variables
  !
  
  IMPLICIT none

  INTEGER :: itermax = 10
  REAL(RP), PARAMETER :: tol = 1.0E-15_RP
  REAL(RP), PARAMETER :: eps = 1.E-50_RP
  REAL(RP), PARAMETER :: ufac = 0.5_RP
  INTEGER :: neq, iter
  REAL(RP) :: tinv, test, delmax, abserr
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: rjac, fb
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: fd, del
  INTEGER, DIMENSION(1) :: imax
  INTEGER :: nm, nx, nr, nl, nm1, nm2, nm3, nm4, nm5, nx1, nx2, nx3, nx4, nx5, nln, nmn
  REAL(RP) :: rr, rta, r1, r2, d1, d2, sm, dnm, test_chem, tst1
  

  ! #################################################################################################
  ! #                                                                                               #
  ! #                                           MAIN BODY                                           #
  ! #                                                                                               #
  ! #################################################################################################

  !
  !  Array sizes, etc.
  !

  IF(lions) THEN
     neq = nchem + 1
  ELSE
     neq = nchem
  END IF

  ALLOCATE(rjac(0:neq,0:neq),fb(neq,neq),fd(neq),del(neq))

  tinv = one/tstep_chem

  IF (nchem == 0) THEN

     !
     !  .. Just calculate prod and loss rates
     !

     pr_chem = zero; ls_chem = zero
     DO nr = 1, nrct

        IF((itype(nr) == 1) .or. (itype(nr) == 7)) THEN

           nm1=irct(1,nr)
           nm3=irct(3,nr)
           nm4=irct(4,nr)
           nm5=irct(5,nr)
           DO nl = 1, nlev
              rr = rt(nr,nl); rta = rr*den(nl,nm1)
              rct(nr,nl) = rta
              ls_chem(nl,nm1) = ls_chem(nl,nm1) + rta
              pr_chem(nl,nm3) = pr_chem(nl,nm3) + rta
              pr_chem(nl,nm4) = pr_chem(nl,nm4) + rta
              pr_chem(nl,nm5) = pr_chem(nl,nm5) + rta              
           END DO

        ELSE

           nm1=irct(1,nr)
           nm2=irct(2,nr)
           nm3=irct(3,nr)
           nm4=irct(4,nr)
           nm5=irct(5,nr)
           DO nl = 1, nlev            
              rr = rt(nr,nl); d1 = den(nl,nm1); d2 = den(nl,nm2)
              r1 = rr*d1; r2 = rr*d2; rta = rr*d1*d2
              rct(nr,nl) = rta
              ls_chem(nl,nm1) = ls_chem(nl,nm1) + rta
              ls_chem(nl,nm2) = ls_chem(nl,nm2) + rta
              pr_chem(nl,nm3) = pr_chem(nl,nm3) + rta
              pr_chem(nl,nm4) = pr_chem(nl,nm4) + rta
              pr_chem(nl,nm5) = pr_chem(nl,nm5) + rta
           END DO

        END IF
     END DO

     DO nl = 1, nlev
        DO nm = 1, nsp
           pr(nl,nm) = prext(nl,nm) + pr_ph(nl,nm) + pr_pe(nl,nm) + pr_chem(nl,nm) 
           ls(nl,nm) = ls_ph(nl,nm) + ls_pe(nl,nm) + ls_chem(nl,nm) 
           dNdt(nl,nm)=(den(nl,nm)-den_old(nl,nm))*tinv 
        END DO
     END DO

     RETURN

  END IF

  !
  !  .. Enforce charge conservation for ions and electrons
  !
     
  DO nl = 1, nlev
     sm = zero
     DO nm = 1,nsp-1
        sm = sm + ichrg(nm)*den(nl,nm)
     END DO
     den(nl,nsp) = sm
  END DO

  !
  !  .. Loop through altitudes, calculate photochem equil at each level
  !

  DO nl = 1, nlev

     iter = 0; test = 10._RP*tol
     delmax = one; abserr = 10._RP*tol
     NEWTON: DO
     
        iter = iter + 1
        IF( (iter>itermax) .OR. (ABS(abserr)<tol) .or. (delmax<1.E-10_RP)) EXIT
          
        !  .. Chemical Reactions
     
        rjac=zero; del = zero; fb=zero; fd=zero
        pr_chem(nl,0:nsp) = zero
        ls_chem(nl,0:nsp) = zero
        DO nr = 1, nrct

        IF((itype(nr) == 1) .or. (itype(nr) == 7)) THEN

           nm1=irct(1,nr); nx1 = lopc(nm1)
           nm3=irct(3,nr); nx3 = lopc(nm3)
           nm4=irct(4,nr); nx4 = lopc(nm4)
           nm5=irct(5,nr); nx5 = lopc(nm5)
           rr = rt(nr,nl); d1 = den(nl,nm1)
           rta = rr*d1; rct(nr,nl) = rta
           ls_chem(nl,nm1) = ls_chem(nl,nm1) + rta
           pr_chem(nl,nm3) = pr_chem(nl,nm3) + rta
           pr_chem(nl,nm4) = pr_chem(nl,nm4) + rta
           pr_chem(nl,nm5) = pr_chem(nl,nm5) + rta
           rjac(nx3,nx1) = rjac(nx3,nx1) + rr
           rjac(nx4,nx1) = rjac(nx4,nx1) + rr
           rjac(nx5,nx1) = rjac(nx5,nx1) + rr
           rjac(nx1,nx1) = rjac(nx1,nx1) - rr

        ELSE

           nm1=irct(1,nr); nx1 = lopc(nm1)
           nm2=irct(2,nr); nx2 = lopc(nm2)
           nm3=irct(3,nr); nx3 = lopc(nm3)
           nm4=irct(4,nr); nx4 = lopc(nm4)
           nm5=irct(5,nr); nx5 = lopc(nm5)
           rr = rt(nr,nl); d1 = den(nl,nm1); d2 = den(nl,nm2)
           r1 = rr*d1; r2 = rr*d2; rta = rr*d1*d2
           rct(nr,nl) = rta
           
           ! .. Loss rates                                       
           
           ls_chem(nl,nm1) = ls_chem(nl,nm1) + rta
           ls_chem(nl,nm2) = ls_chem(nl,nm2) + rta

!           IF(nr == 2563) THEN
!              WRITE(*,"(' CHECK: ',A12,F9.0,2ES11.3)") name(nm1),1.E-5_RP*z(nl),rta,ls_chem(nl,nm1)
!           END IF
        
           ! .. Production rates                                 
        
           pr_chem(nl,nm3) = pr_chem(nl,nm3) + rta
           pr_chem(nl,nm4) = pr_chem(nl,nm4) + rta
           pr_chem(nl,nm5) = pr_chem(nl,nm5) + rta
           
           ! .. Terms in Jacobian
        
           rjac(nx3,nx1) = rjac(nx3,nx1) + r2
           rjac(nx3,nx2) = rjac(nx3,nx2) + r1
           rjac(nx4,nx1) = rjac(nx4,nx1) + r2
           rjac(nx4,nx2) = rjac(nx4,nx2) + r1
           rjac(nx5,nx1) = rjac(nx5,nx1) + r2
           rjac(nx5,nx2) = rjac(nx5,nx2) + r1
           rjac(nx1,nx2) = rjac(nx1,nx2) - r1
           rjac(nx1,nx1) = rjac(nx1,nx1) - r2
           rjac(nx2,nx2) = rjac(nx2,nx2) - r1
           rjac(nx2,nx1) = rjac(nx2,nx1) - r2
        END IF

        END DO

        DO nm = 1, nsp
           pr(nl,nm) = prext(nl,nm) + pr_ph(nl,nm) + pr_pe(nl,nm) + pr_chem(nl,nm) 
           ls(nl,nm) = ls_ph(nl,nm) + ls_pe(nl,nm) + ls_chem(nl,nm) 
           dNdt(nl,nm)=(den(nl,nm)-den_old(nl,nm))*tinv 
        END DO

        ! .. Calculate Jacobian 

        DO nx = 1, nchem
           nm = locp(nx)
           del(nx) = pr(nl,nm)-ls(nl,nm)-dNdt(nl,nm)
           fd(nx) = del(nx)
           fb(nx,nx) = tinv 
           DO nx1 = 1, nchem
              nm1 = locp(nx1)
              fb(nx,nx1) = fb(nx,nx1) - rjac(nx,nx1)
           END DO             
        END DO        
     
        !  .. Jacobian for electrons
     
        IF (lions) THEN
           sm = zero
           DO nm = neutrmax+1,nsp-1
              sm = sm + ichrg(nm)*den(nl,nm)
           END DO
           del(neq) = (sm-den(nl,nsp))*tinv
           fd(neq) = del(neq)
           DO nx = 1, nchem
              nm = locp(nx)
              fb(neq,nx) = -tinv*ichrg(nm)
           END DO
           fb(neq,neq) = tinv
        END IF
        
        !  .. Solve N' = N - (dF/dN)^-1 * F(N); fb*del = fd
        
        CALL GESV(fb,del)
           
        !  .. Update Densities
           
        DO nx = 1, nchem
           nm = locp(nx)
           IF (del(nx) < 0) THEN
              del(nx) = MAX(del(nx),-ufac*den(nl,nm))
           ELSE IF(del(nx) >= 0) THEN
              del(nx) = MIN(del(nx),ufac*den(nl,nm))
           END IF
           den(nl,nm) = MAX(den(nl,nm) + del(nx),eps)    
        END DO
        
        !  .. Update Electron Density
        
        sm = zero
        DO nm = neutrmax+1,nsp-1
           sm = sm + ichrg(nm)*den(nl,nm)
        END DO
        den(nl,nsp) = sm
           
        !  .. Locate Max Change
           
        imax = MAXLOC(ABS(fd(1:neq)))
        abserr = ABS(fd(imax(1)))
        imax = MAXLOC(ABS(del(1:neq)))
        delmax = ABS(del(imax(1)))

        !     DO nx = 1, neq
        !        nm = locp(nx)
        !        WRITE(*,"(I4,2X,A12,10(2X,ES11.3))") iter,name(nm),den(nm), den_old(nm), prz(nm), lsz(nm),dNdt(nm),fd(nx)
        !     END DO
!        WRITE(*,"(10X,'  CHEMEQ: NL = ',I3,2X,' ITER=',I3,2X,' F=',ES11.3,2X,' DelN=',ES11.3,  &
!         ' SPECIES = ',A12)") nl, iter, abserr, delmax,name(locp(imax(1)))
           
     END DO NEWTON
     
  END DO
  
  DEALLOCATE(rjac,fb,fd,del)
  
  RETURN
END SUBROUTINE CHEMEQ
