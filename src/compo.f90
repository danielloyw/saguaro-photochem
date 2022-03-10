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
  REAL(RP), PARAMETER :: bal_tol = 1.0E-03_RP
  REAL(RP), PARAMETER :: eps = 1.E-50_RP
  REAL(RP), PARAMETER :: ufac = 0.1_RP ! damping factor
  INTEGER :: iter, iopt, ierr
  REAL(RP) :: tinv, denerr, balerr
  REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: rjac, fa, fb, fc
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: fd, del, fdmax, den_lst
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: ldiffuse

  INTEGER :: nl, nm, nx, nr, nm1, nm2, nm3, nm4, nm5, nx1, nx2, nx3, nx4, nx5, np,  &
       imxdnx, imxdnl, imxbnx, imxbnl
  REAL(RP) :: rr, r1, r2, pvap, pprs
!  REAL(RP) :: Hi, Ha, Hh
  INTEGER :: nz
  REAL(RP) :: den0, flx0, flx1
  CHARACTER(len=1) :: iret


  ! #################################################################################################
  ! #                                                                                               #
  ! #                                         INITIALIZATION                                        #
  ! #                                                                                               #
  ! #################################################################################################

  imxdnx = 1; imxdnl = 1; imxbnx = 1; imxbnl = 1

  IF (ndiff == 0) THEN

     !  Just calculate production and loss rates without updating densities

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

     del=zero; fa=zero; fb=zero; fc=zero; fd=zero; flx=zero; div_flx=zero

     !    ---------------------------------------------------------------------------
     !          Solar Energy Deposition
     !    ---------------------------------------------------------------------------

     IF(lcompo_sol) THEN
        CALL SOLDEP1   ! Solar Enenergy Deposition
!        CALL ELDEP1    ! Photoelectrons
     END IF

     !    ---------------------------------------------------------------------------
     !          Chemical Production and Loss Rates
     !    ---------------------------------------------------------------------------

     rjac = zero

     ! .. Photons
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
           rpe(np,nl) = eph(np,nl)*den(nl,nm1)
           ls_pe(nl,nm1) = ls_pe(nl,nm1) + rpe(np,nl)
           pr_pe(nl,nm2) = pr_pe(nl,nm2) + rpe(np,nl)
           pr_pe(nl,nm3) = pr_pe(nl,nm3) + rpe(np,nl)
           pr_pe(nl,nm4) = pr_pe(nl,nm4) + rpe(np,nl)
           pr_pe(nl,nm5) = pr_pe(nl,nm5) + rpe(np,nl)
!           rjac(nx2,nx1,nl) = rjac(nx2,nx1,nl) + eph(np,nl)
!           rjac(nx3,nx1,nl) = rjac(nx3,nx1,nl) + eph(np,nl)
!           rjac(nx4,nx1,nl) = rjac(nx4,nx1,nl) + eph(np,nl)
!           rjac(nx5,nx1,nl) = rjac(nx5,nx1,nl) + eph(np,nl)
        END DO
     END DO

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
              ls_chem(nl,nm1) = ls_chem(nl,nm1) + rct(nr,nl)
              pr_chem(nl,nm3) = pr_chem(nl,nm3) + rct(nr,nl)
              pr_chem(nl,nm4) = pr_chem(nl,nm4) + rct(nr,nl)
              pr_chem(nl,nm5) = pr_chem(nl,nm5) + rct(nr,nl)
              rjac(nx1,nx1,nl) = rjac(nx1,nx1,nl) - rt(nr,nl)
              rjac(nx3,nx1,nl) = rjac(nx3,nx1,nl) + rt(nr,nl)
              rjac(nx4,nx1,nl) = rjac(nx4,nx1,nl) + rt(nr,nl)
              rjac(nx5,nx1,nl) = rjac(nx5,nx1,nl) + rt(nr,nl)
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
              ls_chem(nl,nm1) = ls_chem(nl,nm1) + rct(nr,nl)
              ls_chem(nl,nm2) = ls_chem(nl,nm2) + rct(nr,nl)
              pr_chem(nl,nm3) = pr_chem(nl,nm3) + rct(nr,nl)
              pr_chem(nl,nm4) = pr_chem(nl,nm4) + rct(nr,nl)
              pr_chem(nl,nm5) = pr_chem(nl,nm5) + rct(nr,nl)
              rjac(nx1,nx2,nl) = rjac(nx1,nx2,nl) - r1
              rjac(nx1,nx1,nl) = rjac(nx1,nx1,nl) - r2
              rjac(nx2,nx2,nl) = rjac(nx2,nx2,nl) - r1
              rjac(nx2,nx1,nl) = rjac(nx2,nx1,nl) - r2
              rjac(nx3,nx1,nl) = rjac(nx3,nx1,nl) + r2
              rjac(nx3,nx2,nl) = rjac(nx3,nx2,nl) + r1
              rjac(nx4,nx1,nl) = rjac(nx4,nx1,nl) + r2
              rjac(nx4,nx2,nl) = rjac(nx4,nx2,nl) + r1
              rjac(nx5,nx1,nl) = rjac(nx5,nx1,nl) + r2
              rjac(nx5,nx2,nl) = rjac(nx5,nx2,nl) + r1
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

     rcdn = zero
     DO nx = 1, ndiff
        nm = ldcp(nx)
        DO nl = 1, nlev
           pvap = VAPOR(name(nm),tn(nl))
           pprs = rkb*den(nl,nm)*tn(nl)
           IF(pprs > pvap) THEN
              rcdn(nl,nm) = 1.E+12_RP*(pprs - pvap) ! what is 10^12?
              rjac(nx,nx,nl) = rjac(nx,nx,nl) - 1.0E+012_RP*rkb*tn(nl)
              IF(rjac(nx,nx,nl) .ne. rjac(nx,nx,nl)) THEN       !if it is NaN
                 WRITE(*,"(' VAPOR PROBLEM: pvap = ',ES11.3)") rjac(nx,nx,nl)
                 STOP
              END IF
           END IF
        END DO
     END DO

     IF( (iter>itmin) .and.( (iter>itmax) .or. (ABS(denerr)<den_tol) .or. (ABS(balerr)<bal_tol) )) EXIT

     !-------------------------------------------------------------------------
     !  Section 4: Collect chemcial and diffusion terms, perpare for inversion
     !-------------------------------------------------------------------------

     DO nx = 1, ndiff
        nm = ldcp(nx)

        !  .. Interior

        DO nl = 2, nlev-1

           flx(nl,nm) = (alpha(nl,nx)-beta(nl,nx))*den(nl,nm)            &
                      - (alpha(nl,nx)+beta(nl,nx))*den(nl+1,nm)
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
        ELSE IF (ibnd(nm,1) == 2) THEN               !  .. Specified Velocity
!           den0 = ((two*bval(nm,1)+alpha(nl-1,nx)+beta(nl-1,nx)-alpha(nl,nx)+beta(nl,nx))/(alpha(nl-1,nx)-beta(nl-1,nx)))*den(nl,nm) &
!                + ((alpha(nl,nx)+beta(nl,nx))/(alpha(nl-1,nx)-beta(nl-1,nx)))*den(nl+1,nm)
!           div_flx(nl,nm) = a(nl,nx)*den0                             &
!                          + b(nl,nx)*den(nl,nm)                               &
!                          + c(nl,nx)*den(nl+1,nm)
!           flx0 = (alpha(0,nx)-beta(0,nx))*den0-(alpha(0,nx)+beta(0,nx))*den(1,nm)
           flx(nl,nm) = bval(nm,1)*den(nl,nm)
           div_flx(nl,nm) = b(nl,nx)*den(nl,nm) + c(nl,nx)*den(nl+1,nm)
           dNdt(nl,nm) = (den(nl,nm)-den_old(nl,nm))*tinv
           fd(nx,nl) = pr(nl,nm)-rcdn(nl,nm)-ls(nl,nm)-div_flx(nl,nm) - dNdt(nl,nm)
           fdmax(nx,nl) = MAX(ABS(pr(nl,nm)),ABS(rcdn(nl,nm)),ABS(ls(nl,nm)),ABS(div_flx(nl,nm)),ABS(dNdt(nl,nm)))
           fa(nx,nx,nl) = zero
!           fb(nx,nx,nl) = tinv + b(nl,nx)
           fb(nx,nx,nl) = b(nl,nx)
           fc(nx,nx,nl) = c(nl,nx)
           DO nx1 = 1, ndiff
              fb(nx,nx1,nl) = fb(nx,nx1,nl) - rjac(nx,nx1,nl)
           END DO

        ELSE IF (ibnd(nm,1) == 3) THEN                !  .. Fixed mole fraction
           flx(nl,nm) = (alpha(nl,nx)-beta(nl,nx))*den(nl,nm)            &
                - (alpha(nl,nx)+beta(nl,nx))*den(nl+1,nm)
           div_flx(nl,nm) = b(nl,nx)*den(nl,nm) + c(nl,nx)*den(nl+1,nm)
           dNdt(nl,nm) = tinv*(den(nl,nm)-bval(nm,1)*den(nl,0))
           fd(nx,nl) = -dNdt(nl,nm)
           fb(nx,nx,nl) = tinv

           fdmax(nx,nl) = ABS(dNdt(nl,nm))

        ELSE IF (ibnd(nm,1) == 4) THEN                !  .. Specified flux
           flx(nl,nm) = bval(nm,1)
           div_flx(nl,nm) = flx(nl,nm)/rm2(nl)/drp(nl-1)+b(nl,nx)*den(nl,nm)+c(nl,nx)*den(nl+1,nm)
           dNdt(nl,nm) = (den(nl,nm)-den_old(nl,nm))*tinv
           fd(nx,nl) = -div_flx(nl,nm) - dNdt(nl,nm)
           fdmax(nx,nl) = MAX(ABS(div_flx(nl,nm)),ABS(dNdt(nl,nm)))
           fa(nx,nx,nl) = zero
           fb(nx,nx,nl) = b(nl,nx)
           fc(nx,nx,nl) = c(nl,nx)
           DO nx1 = 1, ndiff
              fb(nx,nx1,nl) = fb(nx,nx1,nl) - rjac(nx,nx1,nl)
           END DO

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
           fd(nx,nl)=pr(nl,nm)-rcdn(nl,nm)-ls(nl,nm)-div_flx(nl,nm)-dNdt(nl,nm)
           fdmax(nx,nl) = MAX(ABS(pr(nl,nm)),ABS(rcdn(nl,nm)),ABS(ls(nl,nm)),ABS(div_flx(nl,nm)),ABS(dNdt(nl,nm)))
           fa(nx,nx,nl) = a(nl,nx)
           fb(nx,nx,nl) = tinv + b(nl,nx)
           DO nx1 = 1, ndiff
              fb(nx,nx1,nl) = fb(nx,nx1,nl)-rjac(nx,nx1,nl)
           END DO
        ELSE IF(ibnd(nm,2) == 3) THEN                 !  .. Specified Velocity
           flx(nl,nm) = bval(nm,2)*den(nl,nm)
           div_flx(nl,nm) = a(nl,nx)*den(nl-1,nm)+b(nl,nx)*den(nl,nm)
           dNdt(nl,nm) = (den(nl,nm)-den_old(nl,nm))*tinv
           fd(nx,nl)=pr(nl,nm)-rcdn(nl,nm)-ls(nl,nm)-div_flx(nl,nm)-dNdt(nl,nm)
           fdmax(nx,nl) = MAX(ABS(pr(nl,nm)),ABS(rcdn(nl,nm)),ABS(ls(nl,nm)),ABS(div_flx(nl,nm)),ABS(dNdt(nl,nm)))
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
!        nm = ldcp(nx)
!        DO nl = 1, nlev
!           IF(isnan(fa(nx,nx,nl))) THEN
!              WRITE(*,"(' NaN found in fa: ',A12,3I6)") name(nm),nl, nx
!              STOP
!           END IF
!           IF(isnan(fc(nx,nx,nl))) THEN
!              WRITE(*,"(' NaN found in fc: ',A12,2I6)") name(nm), nl, nx


!              STOP
!           END IF
!           DO nx1 = 1, ndiff
!              nm1 = ldcp(nx1)
!              IF(isnan(fb(nx,nx1,nl))) THEN
!                 WRITE(*,"(' NaN found in fb: ',2A12,3I6)") name(nm),name(nm1), nx, nx1, nl
!                 STOP
!              END IF
!           END DO
!        END DO
!     END DO NEWTON

!     nx = 1
!     nm = ldcp(nx)
!     DO nl= nlev, 1, -1
!        WRITE(*,"(' COMPO: ',I6,10ES11.3)") nl, cm_to_km*z(nl), den(nl,nm), den(nl,nm)/den(nl,0), flx(nl,nm), div_flx(nl,nm)
!     END DO
!     READ(*,"(A)") iret

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

     DO nl = 1, nlev
        den(nl,0) = sum(den(nl,1:nsp))
     END DO

     !  .. Locate Max Change

     denerr = zero ! convergence by del->0
     DO nl = 1, nlev
     DO nx = 1, ndiff
        nm = ldcp(nx)
!        IF((ABS(del(nx,nl)/den(nl,nm)) > denerr) .and. (den(nl,nm) > 1.0_RP)) THEN
        IF((ABS(del(nx,nl)) > denerr) .and. (den(nl,nm) > 1.0E2_RP)) THEN
!           denerr = ABS(del(nx,nl)/den(nl,nm))
           denerr = ABS(del(nx,nl))
           imxdnx = nx
           imxdnl = nl
        END IF
     END DO
  END DO
!     nx = 1
!     nm = ldcp(nx)
!     DO nl=1, nlev
!        WRITE(*,"(' COMPO: ',I6,10ES11.3)") nl, cm_to_km*z(nl), pr(nl,nm), ls(nl,nm), div_flx(nl,nm), dNdt(nl,nm), fd(nx,nl), del(nx,nm), den_lst(nl,nm), den(nl,nm)
!     END DO
!     READ(*,"(A)") iret

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

     balerr = zero ! convergence by F->0
     imxbnx = 0
     imxbnl = 0
     DO nl = 1, nlev
     DO nx = 1, ndiff
        nm = ldcp(nx)
!        IF((ABS(fd(nx,nl)/fdmax(nx,nl)) > balerr) .and. (.not. ldiffuse(nx,nl)) .and. (fdmax(nx,nl) > 1.E-10_RP)) THEN
!           balerr = ABS(fd(nx,nl)/fdmax(nx,nl))
        IF(ABS(fd(nx,nl)) > balerr) THEN
           balerr = ABS(fd(nx,nl))
!           denerr = ABS(del(nx,nl))
           imxbnx = nx
           imxbnl = nl
        END IF
     END DO
     END DO

     IF(imxbnx == 0) THEN
        WRITE(*,"(' max balerr not found')")
     ELSE
        WRITE(*,"(10X,'  COMPO : ITER== ',I4,2X,A12,' Z = ',F8.1,' N = ',ES11.3,' P = ',ES11.3,    &
          ' L = ',ES11.3,' C = ',ES11.3,' Div_Flx = ',ES11.3,' dNdt = ',ES11.3,' bal = ',ES11.3,' delN = ',ES11.3)")                         &
          iter, name(ldcp(imxbnx)), 1.E-5_RP*z(imxbnl), den(imxbnl,ldcp(imxbnx)), pr(imxbnl,ldcp(imxbnx)),  &
          ls(imxbnl,ldcp(imxbnx)), rcdn(imxbnl,ldcp(imxbnx)), div_flx(imxbnl,ldcp(imxbnx)), dNdt(imxbnl,ldcp(imxbnx)),fd(imxbnx,imxbnl), del(imxbnx,imxbnl)
     END IF
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

 !    IF(ABS(balerr)>ABS(balerr_lst)) THEN
 !       den=den_lst
 !    IF (lprnt) WRITE(*,"(' No Converge  ',' delN = ',ES11.3,' N = ',ES11.3,' Mol = ',A12,' Z = ',F10.1)") &
 !               denerr, den(imxdnl,locp(imxdnx)), name(locp(imxdnx)), 1.E-5_RP*(rz(imxdnl)-RPLANET),    &
 !               balerr, fd(imxbnx,imxbnl), name(locp(imxbnx)), 1.E-5_RP*(rz(imxdnl)-RPLANET)
!        IF (lprnt) WRITE(*,901) iter, abserr_lst, relerr, delmax, den(imax(2),locp(imax(1))),imax(2),name(locp(imax(1)))
!        EXIT
!     END IF

  END DO NEWTON
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
        END DO
     END DO

     !  .. Photoelectrons

     pr_pe = zero; ls_pe = zero
     DO np = 1, nert
        nm1=iert(1,np)
        nm2=iert(2,np)
        nm3=iert(3,np)
        nm4=iert(4,np)
        nm5=iert(5,np)
        DO nl = nibot, nlev
           rpe(np,nl) = eph(np,nl)*den(nl,nm1)
           ls_pe(nl,nm1) = ls_pe(nl,nm1) + rpe(np,nl)
           pr_pe(nl,nm2) = pr_pe(nl,nm2) + rpe(np,nl)
           pr_pe(nl,nm3) = pr_pe(nl,nm3) + rpe(np,nl)
           pr_pe(nl,nm4) = pr_pe(nl,nm4) + rpe(np,nl)
           pr_pe(nl,nm5) = pr_pe(nl,nm5) + rpe(np,nl)
        END DO
     END DO

     !  .. Chemical Reactions

     pr_chem = zero; ls_chem = zero
     DO nr = 1, nrct

        IF((itype(nr) == 1) .or. (itype(nr) == 7)) THEN

           nm1=irct(1,nr); nx1 = ldpc(nm1)
           nm3=irct(3,nr); nx3 = ldpc(nm3)
           nm4=irct(4,nr); nx4 = ldpc(nm4)
           nm5=irct(5,nr); nx5 = ldpc(nm5)
           DO nl = 1, nlev
              rct(nr,nl) = rt(nr,nl)*den(nl,nm1)
              ls_chem(nl,nm1) = ls_chem(nl,nm1) + rct(nr,nl)
              pr_chem(nl,nm3) = pr_chem(nl,nm3) + rct(nr,nl)
              pr_chem(nl,nm4) = pr_chem(nl,nm4) + rct(nr,nl)
              pr_chem(nl,nm5) = pr_chem(nl,nm5) + rct(nr,nl)
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
              ls_chem(nl,nm1) = ls_chem(nl,nm1) + rct(nr,nl)
              ls_chem(nl,nm2) = ls_chem(nl,nm2) + rct(nr,nl)
              pr_chem(nl,nm3) = pr_chem(nl,nm3) + rct(nr,nl)
              pr_chem(nl,nm4) = pr_chem(nl,nm4) + rct(nr,nl)
              pr_chem(nl,nm5) = pr_chem(nl,nm5) + rct(nr,nl)
           END DO

        END IF

     END DO

     !  .. Add up production and loss rates

     pr(:,:) = prext(:,:) + pr_ph(:,:) + pr_pe(:,:) + pr_chem(:,:)
     ls(:,:) = ls_ph(:,:) + ls_pe(:,:) + ls_chem(:,:)


     DO nx = 1, ndiff
        nm = ldcp(nx)

        !  .. Interior

        DO nl = 2, nlev-1
           flx(nl,nm) = half*((alpha(nl,nx)-beta(nl,nx))*den(nl,nm)-(alpha(nl,nx)+beta(nl,nx))*den(nl+1,nm)       &
                      +       (alpha(nl-1,nx)-beta(nl-1,nx))*den(nl-1,nm)-(alpha(nl-1,nx)+beta(nl-1,nx))*den(nl,nm))
           div_flx(nl,nm) = a(nl,nx)*den(nl-1,nm) + b(nl,nx)*den(nl,nm) + c(nl,nx)*den(nl+1,nm)
           dNdt(nl,nm)=(den(nl,nm)-den_old(nl,nm))*tinv
        END DO

        !  .. Lower Boundary

        nl = 1
        IF (ibnd(nm,1) == 1) THEN                    !  .. Chemical Equilibrium
           flx(nl,nm) = zero
           div_flx(nl,nm) = zero
           dNdt(nl,nm) = (den(nl,nm)-den_old(nl,nm))*tinv
        ELSE IF (ibnd(nm,1) == 2) THEN               !  .. Fixed Velocity
           div_flx(nl,nm) = b(nl,nx)*den(nl,nm) + c(nl,nx)*den(nl+1,nm)
           flx(nl,nm) = bval(nm,1)*den(nl,nm)
           dNdt(nl,nm) = (den(nl,nm)-den_old(nl,nm))*tinv
        ELSE IF (ibnd(nm,1) == 3) THEN                !  .. Fixed mole fraction
!           flx(nl,nm) = (alpha(nl,nx)-beta(nl,nx))*den(nl,nm) - (alpha(nl,nx)+beta(nl,nx))*den(nl+1,nm)
!           pr(nl,nm) = zero
!           ls(nl,nm) = zero
!           div_flx(nl,nm) = zero
           den0 = (pr(nl,nm)-ls(nl,nm)-b(nl,nx)*den(nl,nm)-c(nl,nx)*den(nl+1,nm))/a(nl,nx)
           flx0 = (alpha(nl-1,nx)-beta(nl-1,nx))*den0 - (alpha(nl-1,nx)+beta(nl-1,nx))*den(nl,nm)
           flx1 = (alpha(nl,nx)-beta(nl,nx))*den(nl,nm) - (alpha(nl,nx)+beta(nl,nx))*den(nl+1,nm)
           flx(nl,nm) = half*(flx0+flx1)
           dNdt(nl,nm) = tinv*(den(nl,nm)-bval(nm,1)*den(nl,0))
           div_flx(nl,nm) = pr(nl,nm)-ls(nl,nm)-dNdt(nl,nm)
           WRITE(*,"(' COMPO: den0,flx0,flx1,flx,div_flx = ',10(2X,ES11.3))") den0, flx0, flx1, flx(nl,nm), div_flx(nl,nm)
        ELSE IF (ibnd(nm,1) == 4) THEN                !  .. Fixed flux
           flx(nl,nm) = bval(nm,1)
           div_flx(nl,nm) = flx(nl,nm)/rm2(nl)/drp(nl-1)+b(nl,nx)*den(nl,nm)+c(nl,nx)*den(nl+1,nm)
           dNdt(nl,nm) = (den(nl,nm)-den_old(nl,nm))*tinv
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
        ELSE IF(ibnd(nm,2) == 2) THEN                 !  .. Jeans Velocity
           flx(nl,nm) = bval(nm,2)*den(nl,nm)
           div_flx(nl,nm) = a(nl,nx)*den(nl-1,nm)+b(nl,nx)*den(nl,nm)
           dNdt(nl,nm) = (den(nl,nm)-den_old(nl,nm))*tinv
        ELSE IF(ibnd(nm,2) == 3) THEN                 !  .. Specified Velocity
           flx(nl,nm) = bval(nm,2)*den(nl,nm)
           div_flx(nl,nm) = a(nl,nx)*den(nl-1,nm)+b(nl,nx)*den(nl,nm)
           dNdt(nl,nm) = (den(nl,nm)-den_old(nl,nm))*tinv
        ELSE IF(ibnd(nm,2) == 4) THEN                 !  .. Specified Flux
           flx(nl,nm) = bval(nm,2)
           div_flx(nl,nm) = flx(nl,nm)/rm2(nl)/drp(nl-1)+a(nl,nx)*den(nl-1,nm)+b(nl,nx)*den(nl,nm)
           dNdt(nl,nm) = (den(nl,nm)-den_old(nl,nm))*tinv
        ELSE
           WRITE(*,*) ' COMPO: ERROR IN UPPER B.C., EXITING ...'
           STOP
        END IF

     END DO

!     DO nx = 1, ndiff
!        nm = ldcp(nx)
!        WRITE(*,"(' WRITE ',A)") name(nm)(1:LEN_TRIM(name(nm)))
!        DO nl = 1, nlev
!            WRITE(*,"(' COMPO: ',F10.1,20(1X,ES10.3))") cm_to_km*z(nl), den(nl,nm), den(nl,nm)/den(nl,0), flx(nl,nm),             &
!                prext(nl,nm), pr_ph(nl,nm), ls_ph(nl,nm), pr_pe(nl,nm), ls_pe(nl,nm),               &
!                pr_chem(nl,nm), ls_chem(nl,nm), pr(nl,nm), ls(nl,nm), rcdn(nl,nm), div_flx(nl,nm),  &
!                pr(nl,nm)-ls(nl,nm)-rcdn(nl,nm)-div_flx(nl,nm)
!         END DO
!        READ(*,"(A)") iret
!     END DO

!     DO nx = 1, ndiff
!        nm = ldcp(nx)
!        DO nl = 1, nlev-1
!           Hi = 1.E-5_RP*half*(den(nl+1,nm)+den(nl,nm))*drp(nl)/(den(nl+1,nm)-den(nl,nm))
!           Ha = 1.E-5_RP*half*(den(nl+1,0)+den(nl,0))*drp(nl)/(den(nl+1,0)-den(nl,0))
!           Hh = 1.E-5_RP*half*rkb*(tn(nl+1)*rz(nl+1)**2/GM/amu/mass(nl+1)+tn(nl)*rz(nl)**2/GM/amu/mass(nl))
!           WRITE(*,"(' COMPO: ',I5,20(2X,ES11.3))") nl, Ha, Hi, Hh
!           !DO nl = 1, nlev-1
!           !WRITE(*,"(I5,20(2X,ES11.3))") nl, (den(nl+1,nm)-den(nl,nm))/drp(nl), half*(den(nl+1,nm)+den(nl,nm))/(half*(ht(nl,0)+ht(nl+1,0))), &
!           !     (den(nl+1,nm)-den(nl,nm))/drp(nl) +  half*(den(nl+1,nm)+den(nl,nm))/(half*(ht(nl,0)+ht(nl+1,0))), &
!           !     -alpha(nl,nx)*(den(nl+1,nm)-den(nl,nm))-beta(nl,nx)*(den(nl+1,nm)+den(nl,nm))
!           !END DO
!        END DO
!        READ(*,"(A)") iret
!     END DO
     ! flx = -K*(dN/dr + N*(1/Ha + (1/T)*dT/dr) - D*(dN/dr + N*(1/H + (1/T)*dT/dr)

  DEALLOCATE(rjac,fa,fb,fc,fd,del,fdmax,den_lst,ldiffuse)

900 FORMAT(10X,'  COMPO: ITER=',I3,2X,' ABSERR=',ES11.3,2X,' RELERR = ',ES11.3,2X,' DelN=',ES11.3,   &
         ' N=',ES11.3,' NL=',I6,' SPECIES = ',A12)

901 FORMAT(10X,'  NO CONVERGE=',I3,2X,' ABSERR=',ES11.3,2X,' RELERR = ',ES11.3,2X,' DelN=',ES11.3,   &
         ' N=',ES11.3,' NL=',I6,' SPECIES = ',A12)

  RETURN

END SUBROUTINE COMPO
