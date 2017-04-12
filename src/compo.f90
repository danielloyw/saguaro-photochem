SUBROUTINE COMPO 

  USE PRECISION
  USE CONSTANTS
  USE GLOBAL_VARIABLES
  USE SUBS_COMPO
  USE SUBS, ONLY : VAPOR


  ! #################################################################################################
  ! #                                                                                               #
  ! #                                           DECLARATIONS                                        #
  ! #                                                                                               #
  ! #################################################################################################


  IMPLICIT none
  INTEGER, PARAMETER :: itmin = 1
  INTEGER, PARAMETER :: itmax = 20
  REAL(RP), PARAMETER :: den_tol = 1.0E-06_RP
  REAL(RP), PARAMETER :: bal_tol = 1.0E-06_RP
  REAL(RP), PARAMETER :: eps = 1.E-50_RP 
  REAL(RP), PARAMETER :: ufac = 0.1_RP
  INTEGER :: iter, iopt, ierr
  REAL(RP) :: tinv, denerr, balerr
  REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: rjac, fa, fb, fc
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: fd, del, fdmax, den_lst
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: ldiffuse

  INTEGER :: nl, nm, nx, nr, nm1, nm2, nm3, nm4, nm5, nx1, nx2, nx3, nx4, nx5, np,  &
       imxdnx, imxdnl, imxbnx, imxbnl
  REAL(RP) :: rr, rta, r1, r2, pvap, pprs
  INTEGER :: nz
  CHARACTER(len=1) :: iret


  ! #################################################################################################
  ! #                                                                                               #
  ! #                                         INITIALIZATION                                        #
  ! #                                                                                               #
  ! #################################################################################################
  
  imxdnx = 1; imxdnl = 1; imxbnx = 1; imxbnl = 1

  IF (ndiff == 0) THEN

     !  Just calculate production and loss rates

     pr_chem = zero; ls_chem = zero
     DO nr = 1, nrct

        IF((itype(nr) == 1) .or. (itype(nr) == 7)) THEN

           nm1=irct(1,nr)
           nm3=irct(3,nr)
           nm4=irct(4,nr)
           nm5=irct(5,nr)
           DO nl = 1, nlev
              rct(nr,nl) = rt(nr,nl)*den(nl,nm1)
              ls_chem(nl,nm1) = ls_chem(nl,nm1) + rct(nr,nl)
              pr_chem(nl,nm3) = pr_chem(nl,nm3) + rct(nr,nl)
              pr_chem(nl,nm4) = pr_chem(nl,nm4) + rct(nr,nl)
              pr_chem(nl,nm5) = pr_chem(nl,nm5) + rct(nr,nl)
           END DO

        ELSE

           nm1=irct(1,nr)
           nm2=irct(2,nr)
           nm3=irct(3,nr)
           nm4=irct(4,nr)
           nm5=irct(5,nr)
           DO nl = 1, nlev
              r1 = rt(nr,nl)*den(nl,nm1); r2 = rt(nr,nl)*den(nl,nm2)
              rct(nr,nl) = rt(nr,nl)*den(nl,nm1)*den(nl,nm2)
              ls_chem(nl,nm1) = ls_chem(nl,nm1) + rct(nr,nl)
              ls_chem(nl,nm2) = ls_chem(nl,nm2) + rct(nr,nl)
              pr_chem(nl,nm3) = pr_chem(nl,nm3) + rct(nr,nl)
              pr_chem(nl,nm4) = pr_chem(nl,nm4) + rct(nr,nl)
              pr_chem(nl,nm5) = pr_chem(nl,nm5) + rct(nr,nl)
           END DO

        END IF
     END DO

     pr(:,:) = prext(:,:) + pr_ph(:,:) + pr_pe(:,:) + pr_chem(:,:)
     ls(:,:) = ls_ph(:,:) + ls_pe(:,:) + ls_chem(:,:)

     RETURN

  END IF

  
  iopt = 0

  ALLOCATE(rjac(0:ndiff,0:ndiff,nlev),fa(ndiff,ndiff,nlev),fb(ndiff,ndiff,nlev),                    & 
       fc(ndiff,ndiff,nlev),fd(ndiff,nlev),del(ndiff,nlev),fdmax(ndiff,nlev),den_lst(nlev,0:nsp),   &
       ldiffuse(ndiff,nlev))

  ! #################################################################################################
  ! #                                                                                               #
  ! #                      CALCULATE DENSITIES AT NEXT TIME STEP WITH NEWTON'S METHOD               #
  ! #                                                                                               #
  ! #################################################################################################
  
  tinv = one/tstep_diff
  iter = 0
  balerr = 10._RP*bal_tol
  denerr = 10._RP*den_tol
  NEWTON: DO
     
     iter = iter + 1
     
     IF( (iter>itmin) .and.( (iter>itmax) .or. (ABS(denerr)<den_tol) .or. (ABS(balerr)<bal_tol) )) EXIT
     
     rct=zero; del=zero; rcdn=zero; rjac=zero; fa=zero; fb=zero; fc=zero; fd=zero; flx=zero; div_flx=zero

     !    ---------------------------------------------------------------------------
     !          Photon Rates
     !    ---------------------------------------------------------------------------

     CALL SOLDEP1 

     pr_ph = zero; ls_ph = zero
     DO np = 1, nphrt
        nm1=irpt(1,np); nx1 = ldpc(nm1) 
        nm2=irpt(2,np); nx2 = ldpc(nm2) 
        nm3=irpt(3,np); nx3 = ldpc(nm3) 
        nm4=irpt(4,np); nx4 = ldpc(nm4) 
        nm5=irpt(5,np); nx5 = ldpc(nm5) 
        DO nl = 1, nlev
           rpt(np,nl) = rph(np,nl)*den(nl,nm1)
           ls_ph(nl,nm1) = ls_ph(nl,nm1) + rpt(np,nl)
           pr_ph(nl,nm2) = pr_ph(nl,nm2) + rpt(np,nl)
           pr_ph(nl,nm3) = pr_ph(nl,nm3) + rpt(np,nl)
           pr_ph(nl,nm4) = pr_ph(nl,nm4) + rpt(np,nl)
           pr_ph(nl,nm5) = pr_ph(nl,nm5) + rpt(np,nl)     
           rjac(nx1,nx1,nl) = rjac(nx1,nx1,nl) - rph(np,nl)
           rjac(nx2,nx1,nl) = rjac(nx2,nx1,nl) + rph(np,nl)
           rjac(nx3,nx1,nl) = rjac(nx3,nx1,nl) + rph(np,nl)
           rjac(nx4,nx1,nl) = rjac(nx4,nx1,nl) + rph(np,nl)
           rjac(nx5,nx1,nl) = rjac(nx5,nx1,nl) + rph(np,nl)

!     IF(rjac(nx1,nx1,nl) .ne. rjac(nx1,nx1,nl)) THEN
!        WRITE(*,"(' COMPO1: ',11I4,10(2X,ES11.3))") np, nx1,nx2,nx3,nx4,nx5,nm1,nm2,nm3,nm4,nm5,d1,d2,rr,r1,r2,rjac(nx1,nx1,nl)
!        WRITE(*,"(10ES11.3)") (den(nz,nm1),nz=1,nlev)
!        STOP
!     END IF

 
        END DO
     END DO

     !-------------------------------------------------------------------------
     !  Section 2:  Suprathermal Electron Production Rates
     !-------------------------------------------------------------------------

     pr_pe = zero; ls_pe = zero
     DO np = 1, nert
        nm1=iert(1,np)
        nm2=iert(2,np)
        nm3=iert(3,np)
        nm4=iert(4,np)
        nm5=iert(5,np)
        DO nl = nibot, nlev
           rr = eph(np,nl)
           rpe(np,nl) = eph(np,nl)*den(nl,nm1)
           ls_pe(nl,nm1) = ls_pe(nl,nm1) + rpe(np,nl)
           pr_pe(nl,nm2) = pr_pe(nl,nm2) + rpe(np,nl)
           pr_pe(nl,nm3) = pr_pe(nl,nm3) + rpe(np,nl)
           pr_pe(nl,nm4) = pr_pe(nl,nm4) + rpe(np,nl)
           pr_pe(nl,nm5) = pr_pe(nl,nm5) + rpe(np,nl) 
!           rjac(nx2,nx1,nl) = rjac(nx2,nx1,nl) + rr
!           rjac(nx3,nx1,nl) = rjac(nx3,nx1,nl) + rr
!           rjac(nx4,nx1,nl) = rjac(nx4,nx1,nl) + rr
!           rjac(nx5,nx1,nl) = rjac(nx5,nx1,nl) + rr
        END DO
     END DO

     !-------------------------------------------------------------------------
     ! Chemical Species
     !-------------------------------------------------------------------------

!     CALL CHEMEQ
          
     !-------------------------------------------------------------------------
     !  Section 3:  Chemical Reactions
     !-------------------------------------------------------------------------
     
     pr_chem = zero; ls_chem = zero
     DO nr = 1, nrct

        IF((itype(nr) == 1) .or. (itype(nr) == 7)) THEN

           nm1=irct(1,nr); nx1 = ldpc(nm1) 
           nm3=irct(3,nr); nx3 = ldpc(nm3) 
           nm4=irct(4,nr); nx4 = ldpc(nm4) 
           nm5=irct(5,nr); nx5 = ldpc(nm5) 
           DO nl = 1, nlev
           
              rct(nr,nl) = rt(nr,nl)*den(nl,nm1)
           
              ! .. Loss rates        
                               
              ls_chem(nl,nm1) = ls_chem(nl,nm1) + rct(nr,nl)

              ! .. Production rates                                 
           
              pr_chem(nl,nm3) = pr_chem(nl,nm3) + rct(nr,nl)
              pr_chem(nl,nm4) = pr_chem(nl,nm4) + rct(nr,nl)
              pr_chem(nl,nm5) = pr_chem(nl,nm5) + rct(nr,nl)
              
              ! .. Terms in Jacobian
           
              rjac(nx3,nx1,nl) = rjac(nx3,nx1,nl) + rt(nr,nl)
              rjac(nx4,nx1,nl) = rjac(nx4,nx1,nl) + rt(nr,nl)
              rjac(nx5,nx1,nl) = rjac(nx5,nx1,nl) + rt(nr,nl)
              rjac(nx1,nx1,nl) = rjac(nx1,nx1,nl) - rt(nr,nl)

!     IF(rjac(nx1,nx1,nl) .ne. rjac(nx1,nx1,nl)) THEN
!        WRITE(*,"(' COMPO2: ',11I4,10(2X,ES11.3))") nr, nx1,nx2,nx3,nx4,nx5,nm1,nm2,nm3,nm4,nm5,d1,d2,rr,r1,r2,rjac(nx1,nx1,nl)
!        STOP
!     END IF

           END DO

        ELSE

           nm1=irct(1,nr); nx1 = ldpc(nm1) 
           nm2=irct(2,nr); nx2 = ldpc(nm2) 
           nm3=irct(3,nr); nx3 = ldpc(nm3) 
           nm4=irct(4,nr); nx4 = ldpc(nm4) 
           nm5=irct(5,nr); nx5 = ldpc(nm5) 
           DO nl = 1, nlev
              
              r1 = rt(nr,nl)*den(nl,nm1); r2 = rt(nr,nl)*den(nl,nm2);
              rct(nr,nl) = rt(nr,nl)*den(nl,nm1)*den(nl,nm2)
              
              ! .. Loss rates        
              
              ls_chem(nl,nm1) = ls_chem(nl,nm1) + rct(nr,nl)
              ls_chem(nl,nm2) = ls_chem(nl,nm2) + rct(nr,nl)
              
              ! .. Production rates                                 
              
              pr_chem(nl,nm3) = pr_chem(nl,nm3) + rct(nr,nl)
              pr_chem(nl,nm4) = pr_chem(nl,nm4) + rct(nr,nl)
              pr_chem(nl,nm5) = pr_chem(nl,nm5) + rct(nr,nl)
              
              ! .. Terms in Jacobian
              
              rjac(nx3,nx1,nl) = rjac(nx3,nx1,nl) + r2
              rjac(nx3,nx2,nl) = rjac(nx3,nx2,nl) + r1
              rjac(nx4,nx1,nl) = rjac(nx4,nx1,nl) + r2
              rjac(nx4,nx2,nl) = rjac(nx4,nx2,nl) + r1
              rjac(nx5,nx1,nl) = rjac(nx5,nx1,nl) + r2
              rjac(nx5,nx2,nl) = rjac(nx5,nx2,nl) + r1
              rjac(nx1,nx2,nl) = rjac(nx1,nx2,nl) - r1
              rjac(nx1,nx1,nl) = rjac(nx1,nx1,nl) - r2
              rjac(nx2,nx2,nl) = rjac(nx2,nx2,nl) - r1
              rjac(nx2,nx1,nl) = rjac(nx2,nx1,nl) - r2

!     IF(rjac(nx1,nx1,nl) .ne. rjac(nx1,nx1,nl)) THEN
!        WRITE(*,"(' COMPO3: ',11I4,10(2X,ES11.3))") nr, nx1,nx2,nx3,nx4,nx5,nm1,nm2,nm3,nm4,nm5,d1,d2,rr,r1,r2,rjac(nx1,nx1,nl)
!        STOP
!     END IF

           END DO

        END IF

     END DO       

     !
     !  .. Add up production and loss rates
     !

     pr(:,:) = prext(:,:) + pr_ph(:,:) + pr_pe(:,:) + pr_chem(:,:)
     ls(:,:) = ls_ph(:,:) + ls_pe(:,:) + ls_chem(:,:)

     
     !-------------------------------------------------------------------------
     !  Section 3: Condensation
     !-------------------------------------------------------------------------
     
     DO nx = 1, ndiff
        nm = ldcp(nx)
        DO nl = 1, nlev
           pvap = VAPOR(name(nm),tn(nl))
           pprs = rkb*den(nl,nm)*tn(nl)
           IF(pprs > pvap) THEN
              rta = 1.E+12_RP*(pprs - pvap)
              rcdn(nl,nm) = rta
              rjac(nx,nx,nl) = rjac(nx,nx,nl) - 1.0E+012_RP*rkb*tn(nl)
              IF(rjac(nx,nx,nl) .ne. rjac(nx,nx,nl)) THEN !?
                 WRITE(*,"(' VAPOR PROBLEM: pvap = ',ES11.3)") rjac(nx,nx,nl)
                 STOP
              END IF
           END IF
        END DO
     END DO

     !-------------------------------------------------------------------------
     !  Section 4: Collect chemcial and diffusion terms, perpare for inversion
     !-------------------------------------------------------------------------
          
     DO nx = 1, ndiff
        nm = ldcp(nx)
        
        !  .. Interior
        
        DO nl = 2, nlev-1
           
           flx(nl,nm) = (alpha(nl,nx)-half*beta(nl,nx))*den(nl,nm)            & 
                      - (alpha(nl,nx)+half*beta(nl,nx))*den(nl+1,nm)

           div_flx(nl,nm) = a(nl,nx)*den(nl-1,nm)                             &
                          + b(nl,nx)*den(nl,nm)                               &
                          + c(nl,nx)*den(nl+1,nm)

           dNdt(nl,nm)=(den(nl,nm)-den_old(nl,nm))*tinv 
           
           fd(nx,nl)=pr(nl,nm)-rcdn(nl,nm)-ls(nl,nm)-div_flx(nl,nm)-dNdt(nl,nm)
           fdmax(nx,nl) = MAX(ABS(pr(nl,nm)),ABS(rcdn(nl,nm)),ABS(ls(nl,nm)),ABS(div_flx(nl,nm)),ABS(dNdt(nl,nm)))
           fa(nx,nx,nl) = a(nl,nx)
           fb(nx,nx,nl) = tinv + b(nl,nx)
           fc(nx,nx,nl) = c(nl,nx)
           DO nx1 = 1, ndiff
              fb(nx,nx1,nl) = fb(nx,nx1,nl) - rjac(nx,nx1,nl)
           END DO

        END DO
        
        !  .. Lower Boundary
        
        nl = 1
        IF (ibnd(nm,1) == 1) THEN                    !  .. Chemical Equilibrium
           flx(nl,nm) = zero
           div_flx(nl,nm) = zero
           dNdt(nl,nm) = (den(nl,nm)-den_old(nl,nm))*tinv
           fd(nx,nl) = pr(nl,nm)-rcdn(nl,nm)-ls(nl,nm)-dNdt(nl,nm)
           fdmax(nx,nl) = MAX(ABS(pr(nl,nm)),ABS(rcdn(nl,nm)),ABS(ls(nl,nm)),ABS(dNdt(nl,nm)))
           fb(nx,nx,nl) = tinv
           DO nx1 = 1, ndiff
              nm1 = ldcp(nx1)
              fb(nx,nx1,nl) = fb(nx,nx1,nl)-rjac(nx,nx1,nl)
           END DO
        ELSE IF (ibnd(nm,1) == 2) THEN               !  .. Fixed Velocity
           flx(nl,nm) = bval(nm,1)*half*(den(nl,nm)+den(nl+1,nm))
           div_flx(nl,nm) = b(nl,nx)*den(nl,nm)+c(nl,nx)*den(nl+1,nm)
           dNdt(nl,nm) = (den(nl,nm)-den_old(nl,nm))*tinv
           fd(nx,nl) = pr(nl,nm)-rcdn(nl,nm)-ls(nl,nm)-div_flx(nl,nm) - dNdt(nl,nm)
           fdmax(nx,nl) = MAX(ABS(pr(nl,nm)),ABS(rcdn(nl,nm)),ABS(ls(nl,nm)),ABS(div_flx(nl,nm)),ABS(dNdt(nl,nm)))
           fa(nx,nx,nl) = zero
           fb(nx,nx,nl) = tinv + b(nl,nx)
           fc(nx,nx,nl) = c(nl,nx)
           DO nx1 = 1, ndiff
              fb(nx,nx1,nl) = fb(nx,nx1,nl) - rjac(nx,nx1,nl)
           END DO

        ELSE IF (ibnd(nm,1) == 3) THEN                !  .. Fixed density
           flx(nl,nm) = (alpha(nl,nx)-half*beta(nl,nx))*den(nl,nm)            & 
                      - (alpha(nl,nx)+half*beta(nl,nx))*den(nl+1,nm)
           dNdt(nl,nm) = -tinv*(den(nl,nm)-bval(nm,1))
           fb(nx,nx,nl) = tinv
           fd(nx,nl) = dNdt(nl,nm)
           fdmax(nx,nl) = ABS(dNdt(nl,nm))
        ELSE
           WRITE(*,*) ' COMPO: ERROR IN LOWER B.C., EXITING ...'
           STOP
        END IF
        
        !  .. Top Boundary
        
        nl = nlev
        IF(ibnd(nm,2) == 1) THEN                     !  .. Chemical Equilibrium
           flx(nl,nm) = zero
           div_flx(nl,nm) = zero
           dNdt(nl,nm) = (den(nl,nm)-den_old(nl,nm))*tinv
           fd(nx,nl) =  pr(nl,nm)-rcdn(nl,nm)-ls(nl,nm) - dNdt(nl,nm)
           fdmax(nx,nl) = MAX(ABS(pr(nl,nm)),ABS(ls(nl,nm)),ABS(dNdt(nl,nm)))
           fb(nx,nx,nl) = tinv
           DO nx1 = 1, ndiff
              fb(nx,nx1,nl) = fb(nx,nx1,nl)-rjac(nx,nx1,nl)
           END DO
        ELSE IF(ibnd(nm,2) == 2) THEN                 !  .. Jeans Velocity
           flx(nl,nm) = bval(nm,2)*den(nl,nm)
           div_flx(nl,nm) = a(nl,nx)*den(nl-1,nm)+b(nl,nx)*den(nl,nm)
           dNdt(nl,nm) = (den(nl,nm)-den_old(nl,nm))*tinv
           fd(nx,nl) = -div_flx(nl,nm) - dNdt(nl,nm)
           fdmax(nx,nl) = MAX(ABS(div_flx(nl,nm)),ABS(dNdt(nl,nm)))
           fa(nx,nx,nl) = a(nl,nx)
           fb(nx,nx,nl) = tinv + b(nl,nx)
           DO nx1 = 1, ndiff
              fb(nx,nx1,nl) = fb(nx,nx1,nl)-rjac(nx,nx1,nl)
           END DO
        ELSE IF(ibnd(nm,2) == 3) THEN                 !  .. Specified Velocity
           flx(nl,nm) = bval(nm,2)*den(nl,nm)
           div_flx(nl,nm) = a(nl,nx)*den(nl-1,nm)+b(nl,nx)*den(nl,nm)
           dNdt(nl,nm) = (den(nl,nm)-den_old(nl,nm))*tinv
           fd(nx,nl) = -div_flx(nl,nm) - dNdt(nl,nm)
           fdmax(nx,nl) = MAX(ABS(div_flx(nl,nm)),ABS(dNdt(nl,nm)))
           fa(nx,nx,nl) = a(nl,nx)
           fb(nx,nx,nl) = tinv + b(nl,nx)
           DO nx1 = 1, ndiff
              fb(nx,nx1,nl) = fb(nx,nx1,nl)-rjac(nx,nx1,nl)
           END DO
        ELSE IF(ibnd(nm,2) == 4) THEN                 !  .. Specified Flux
           flx(nl,nm) = bval(nm,2)
           div_flx(nl,nm) = flx(nl,nm)/rm2(nl)/drp(nl-1)+a(nl,nx)*den(nl-1,nm)+b(nl,nx)*den(nl,nm)
           dNdt(nl,nm) = (den(nl,nm)-den_old(nl,nm))*tinv
           fd(nx,nl) = -div_flx(nl,nm) - dNdt(nl,nm)
           fdmax(nx,nl) = MAX(ABS(div_flx(nl,nm)),ABS(dNdt(nl,nm)))
           fa(nx,nx,nl) = a(nl,nx)
           fb(nx,nx,nl) = tinv + b(nl,nx)
           DO nx1 = 1, ndiff
              fb(nx,nx1,nl) = fb(nx,nx1,nl)-rjac(nx,nx1,nl)
           END DO
        ELSE
           WRITE(*,*) ' COMPO: ERROR IN UPPER B.C., EXITING ...'
           STOP
        END IF
        
     END DO

!     DO nx = 1, ndiff
!        DO nl = 1, nlev
!           IF(fb(nx,nx,nl) .ne. fb(nx,nx,nl)) THEN
!              nm = ldcp(nx)
!              WRITE(*,"(' COMPO:',3I5,20(2X,ES11.3))") nm,nx,nl,b(nl,nx),(rjac(nx,nx1,nl),nx1=1,ndiff),fb(nx,nx,nl)
!              DO nx1 = 1, ndiff
!                 IF(rjac(nx,nx1,nl) .ne. rjac(nx,nx1,nl)) THEN
!                    WRITE(*,"(4I5,2(2X,A12))") nx,nx1,ldcp(nx),ldcp(nx1),name(ldcp(nx)),name(ldcp(nx1))
!                 END IF
!              END DO
!              STOP
!           END IF
!        END DO
!     END DO
     
     !  .. Solve linear system of equations

     CALL SOLVEIT(iopt, fa, fb, fc, fd, del, ierr)

     IF(ierr /= 0) THEN
        WRITE(*,"(' ndiff = ',I5)") ndiff
        DO nx = 1, ndiff
           WRITE(*,"(2I5,2X,A12)") nx,ldcp(nx),name(ldcp(nx))
           DO nl = 1, nlev
              WRITE(*,"(I5,4(2X,ES11.3))") nl,a(nl,nx),b(nl,nx),c(nl,nx),div_flx(nl,nm)
           END DO
           READ(*,"(A)") iret
           DO nl = 1, nlev
              WRITE(*,"(I5,4(2X,ES11.3))") nl,fa(nx,nx,nl),fb(nx,nx,nl),fc(nx,nx,nl),fd(nx,nl)
           END DO
           READ(*,"(A)") iret
        END DO
        STOP
     END IF


!     DO nl = 1, nlev
!        WRITE(*,"(I4,10ES11.3)") nl, den(nl,nm),del(nx,nl),pr(nl,nm),ls(nl,nm),div_flx(nl,nm),dNdt(nl,nm)
!     END DO
!     READ(*,"(A)") iret

     !  .. Locate Max Change

     denerr = zero
     DO nl = 1, nlev
     DO nx = 1, ndiff
        nm = ldcp(nx)
        IF((ABS(del(nx,nl)/den(nl,nm)) > denerr) .and. (den(nl,nm) > 1.0_RP)) THEN
           denerr = ABS(del(nx,nl)/den(nl,nm))
           imxdnx = nx
           imxdnl = nl
        END IF
     END DO
     END DO

     !  .. Update Densities
     
     den_lst = den
     DO nx = 1, ndiff
        nm = ldcp(nx)
        DO  nl = 1, nlev
!           den(nl,nm) = den(nl,nm) + del(nx,nl)   
!           IF((ABS(del(nx,nl)/den(nl,nm)) > ufac) .and. (ABS(del(nx,nl)) > 1._RP)) THEN
!              WRITE(*,"(' WARNING: delN = ',ES11.3,' N = ',ES11.3,2X,' Species = ',A12,' Z = ',F10.1,' F = ',ES11.3)")    &
!                   del(nx,nl), den(nl,nm), name(nm), 1.E-5_RP*(rz(nl)-RPLANET), fd(nx,nl)
!              READ(*,"(A)") iret
!           END IF
           IF (del(nx,nl) < 0) THEN
              del(nx,nl) = MAX(del(nx,nl),-ufac*den(nl,nm))
           ELSE IF(del(nx,nl) >= 0) THEN
              del(nx,nl) = MIN(del(nx,nl),ufac*den(nl,nm))
           END IF
           den(nl,nm) = MAX(den(nl,nm) + del(nx,nl),eps)    
        END DO
     END DO

     !  .. Locate max imbalance

     ldiffuse(:,:) = .false.
     DO nl = 1, nlev
     DO nx = 1, ndiff
        nm = ldcp(nx)
        IF ( (ABS(pr(nl,nm)/div_flx(nl,nm)) < 1.E-3_RP)              &
            .and. (ABS(ls(nl,nm)/div_flx(nl,nm)) < 1.E-3_RP)         &
            .and. (ABS(rcdn(nl,nm)/div_flx(nl,nm)) < 1.E-3_RP)) THEN
           ldiffuse(nx,nl) = .true.
        END IF
     END DO
     END DO

     balerr = zero
     imxbnx = 1
     imxbnl = 1
     DO nl = 1, nlev
     DO nx = 1, ndiff
        nm = ldcp(nx)
        IF((ABS(fd(nx,nl)/fdmax(nx,nl)) > balerr) .and. (.not. ldiffuse(nx,nl)) .and. (fdmax(nx,nl) > 1.E-10_RP)) THEN
           balerr = ABS(fd(nx,nl)/fdmax(nx,nl))
           imxbnx = nx
           imxbnl = nl
        END IF
     END DO
     END DO
     
!     imax = MAXLOC(ABS(fd(1:ndiff,1:nlev)))
!     abserr = fd(imax(1),imax(2))
!     relerr = abserr/fdmax(imax(1),imax(2))
     
!     imax = MAXLOC(ABS(del(1:ndiff,1:nlev)))
!     delmax = del(imax(1),imax(2))
     
!     IF (lprnt) WRITE(*,900) iter, abserr, relerr, delmax, den(imax(2),locp(imax(1))),imax(2),name(locp(imax(1)))
!     IF (lprnt) WRITE(*,"(10X,'  COMPO  : ITER== ',I4,' delN/N = ',ES11.3,' N = ',ES11.3,' Mol = ',A12,' Z = ',F10.1, &
!                          5X,' delF/F = ',ES11.3, ' F = ',ES11.3,' Mol = ',A12,' Z = ',F10.1)")         &
!          iter, denerr, den(imxdnl,ldcp(imxdnx)), name(ldcp(imxdnx)), 1.E-5_RP*(rz(imxdnl)-RPLANET),    &
!                balerr, fd(imxbnx,imxbnl), name(ldcp(imxbnx)), 1.E-5_RP*(rz(imxbnl)-RPLANET)

     IF (lprnt) WRITE(*,"(10X,'  COMPO  : ITER== ',I4,2X,A12,' Z = ',F8.1,' N = ',ES11.3,' P = ',ES11.3,    &
          ' L = ',ES11.3,' C = ',ES11.3,' Div_Flx = ',ES11.3,' delN/N = ',ES11.3)")                         &
          iter, name(ldcp(imxdnx)), 1.E-5_RP*z(imxdnl), den(imxdnl,ldcp(imxdnx)), pr(imxdnl,ldcp(imxdnx)),  &
          ls(imxdnl,ldcp(imxdnx)), rcdn(imxdnl,ldcp(imxdnx)), div_flx(imxdnl,ldcp(imxdnx)), denerr

 !    IF(ABS(balerr)>ABS(balerr_lst)) THEN
 !       den=den_lst
 !    IF (lprnt) WRITE(*,"(' No Converge  ',' delN = ',ES11.3,' N = ',ES11.3,' Mol = ',A12,' Z = ',F10.1)") &
 !               denerr, den(imxdnl,locp(imxdnx)), name(locp(imxdnx)), 1.E-5_RP*(rz(imxdnl)-RPLANET),    &
 !               balerr, fd(imxbnx,imxbnl), name(locp(imxbnx)), 1.E-5_RP*(rz(imxdnl)-RPLANET)
!        IF (lprnt) WRITE(*,901) iter, abserr_lst, relerr, delmax, den(imax(2),locp(imax(1))),imax(2),name(locp(imax(1)))
!        EXIT
!     END IF

  END DO NEWTON
 
  DEALLOCATE(rjac,fa,fb,fc,fd,del,fdmax,den_lst,ldiffuse)

900 FORMAT(10X,'  COMPO: ITER=',I3,2X,' ABSERR=',ES11.3,2X,' RELERR = ',ES11.3,2X,' DelN=',ES11.3,   &
         ' N=',ES11.3,' NL=',I6,' SPECIES = ',A12)

901 FORMAT(10X,'  NO CONVERGE=',I3,2X,' ABSERR=',ES11.3,2X,' RELERR = ',ES11.3,2X,' DelN=',ES11.3,   &
         ' N=',ES11.3,' NL=',I6,' SPECIES = ',A12)
  
  RETURN

END SUBROUTINE COMPO
