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
  REAL(RP), PARAMETER :: tol = 1.0E-10_RP
  REAL(RP), PARAMETER :: eps = 1.E-50_RP ! machine accuracy
  REAL(RP), PARAMETER :: ufac = 0.5_RP
  INTEGER :: neq, iter
  REAL(RP) :: tinv, delmax, abserr
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: rjac, fb
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: fd, del
  INTEGER, DIMENSION(1) :: imax
  INTEGER :: nm, nx, nr, nl, nm1, nm2, nm3, nm4, nm5, nx1, nx2, nx3, nx4, nx5, np
  REAL(RP) :: r1, r2, sm


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
     !  .. Just calculate prod and loss rates without updating densities
     !
     !  .. Photon Reactions

     DO nl = 1, nlev

        pr_ph(nl,:) = zero; ls_ph(nl,:) = zero; rpt(:,nl) = zero
        DO np = 1, nphrt
           nm1=irpt(1,np); nx1 = lopc(nm1)
           nm2=irpt(2,np); nx2 = lopc(nm2)
           nm3=irpt(3,np); nx3 = lopc(nm3)
           nm4=irpt(4,np); nx4 = lopc(nm4)
           nm5=irpt(5,np); nx5 = lopc(nm5)
           rpt(np,nl) = rph(np,nl)*den(nl,nm1)
           ls_ph(nl,nm1) = ls_ph(nl,nm1) + rpt(np,nl)
           pr_ph(nl,nm2) = pr_ph(nl,nm2) + rpt(np,nl)
           pr_ph(nl,nm3) = pr_ph(nl,nm3) + rpt(np,nl)
           pr_ph(nl,nm4) = pr_ph(nl,nm4) + rpt(np,nl)
           pr_ph(nl,nm5) = pr_ph(nl,nm5) + rpt(np,nl)
        END DO

        !  .. Suprathermal electron production and loss

        pr_pe(nl,:) = zero; ls_pe(nl,:) = zero; rpe(:,nl) = zero
        DO np = 1, nert
           nm1=iert(1,np)
           nm2=iert(2,np)
           nm3=iert(3,np)
           nm4=iert(4,np)
           nm5=iert(5,np)
           rpe(np,nl) = eph(np,nl)*den(nl,nm1)
           ls_pe(nl,nm1) = ls_pe(nl,nm1) + rpe(np,nl)
           pr_pe(nl,nm2) = pr_pe(nl,nm2) + rpe(np,nl)
           pr_pe(nl,nm3) = pr_pe(nl,nm3) + rpe(np,nl)
           pr_pe(nl,nm4) = pr_pe(nl,nm4) + rpe(np,nl)
           pr_pe(nl,nm5) = pr_pe(nl,nm5) + rpe(np,nl)
        END DO

        pr_chem(nl,:) = zero; ls_chem(nl,:) = zero; rct(:,nl) = zero
        DO nr = 1, nrct

           IF((itype(nr) == 1) .or. (itype(nr) == 7)) THEN
              nm1=irct(1,nr)
              nm3=irct(3,nr)
              nm4=irct(4,nr)
              nm5=irct(5,nr)
              rct(nr,nl) = rt(nr,nl)*den(nl,nm1)
              ls_chem(nl,nm1) = ls_chem(nl,nm1) + rct(nr,nl)
              pr_chem(nl,nm3) = pr_chem(nl,nm3) + rct(nr,nl)
              pr_chem(nl,nm4) = pr_chem(nl,nm4) + rct(nr,nl)
              pr_chem(nl,nm5) = pr_chem(nl,nm5) + rct(nr,nl)

           ELSE
              nm1=irct(1,nr)
              nm2=irct(2,nr)
              nm3=irct(3,nr)
              nm4=irct(4,nr)
              nm5=irct(5,nr)
              r1 = rt(nr,nl)*den(nl,nm1); r2 = rt(nr,nl)*den(nl,nm2)
              rct(nr,nl) = rt(nr,nl)*den(nl,nm1)*den(nl,nm2)
              ls_chem(nl,nm1) = ls_chem(nl,nm1) + rct(nr,nl)
              ls_chem(nl,nm2) = ls_chem(nl,nm2) + rct(nr,nl)
              pr_chem(nl,nm3) = pr_chem(nl,nm3) + rct(nr,nl)
              pr_chem(nl,nm4) = pr_chem(nl,nm4) + rct(nr,nl)
              pr_chem(nl,nm5) = pr_chem(nl,nm5) + rct(nr,nl)
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

  ALT: DO nl = 1, nlev

     iter = 0
     delmax = one; ! maximum change from iteration
     abserr = 10._RP*tol ! how close F(N) is to 0
     NEWTON: DO

        iter = iter + 1
        IF( (iter>itermax) .OR. (ABS(abserr)<tol) .or. (delmax<1.E-10_RP)) EXIT

        !  .. Photon Reactions
        rjac=zero;
        pr_ph(nl,:) = zero; ls_ph(nl,:) = zero; rpt(:,nl) = zero
        DO np = 1, nphrt
           nm1=irpt(1,np); nx1 = lopc(nm1)
           nm2=irpt(2,np); nx2 = lopc(nm2)
           nm3=irpt(3,np); nx3 = lopc(nm3)
           nm4=irpt(4,np); nx4 = lopc(nm4)
           nm5=irpt(5,np); nx5 = lopc(nm5)
           rpt(np,nl) = rph(np,nl)*den(nl,nm1)
           ls_ph(nl,nm1) = ls_ph(nl,nm1) + rpt(np,nl)
           pr_ph(nl,nm2) = pr_ph(nl,nm2) + rpt(np,nl)
           pr_ph(nl,nm3) = pr_ph(nl,nm3) + rpt(np,nl)
           pr_ph(nl,nm4) = pr_ph(nl,nm4) + rpt(np,nl)
           pr_ph(nl,nm5) = pr_ph(nl,nm5) + rpt(np,nl)

           rjac(nx1,nx1) = rjac(nx1,nx1) - rph(np,nl)
           rjac(nx2,nx1) = rjac(nx2,nx1) + rph(np,nl)
           rjac(nx3,nx1) = rjac(nx3,nx1) + rph(np,nl)
           rjac(nx4,nx1) = rjac(nx4,nx1) + rph(np,nl)
           rjac(nx5,nx1) = rjac(nx5,nx1) + rph(np,nl)
        END DO

        !  .. Suprathermal electron production and loss

        pr_pe(nl,:) = zero; ls_pe(nl,:) = zero; rpe(:,nl) = zero
        DO np = 1, nert
           nm1=iert(1,np)
           nm2=iert(2,np)
           nm3=iert(3,np)
           nm4=iert(4,np)
           nm5=iert(5,np)
           rpe(np,nl) = eph(np,nl)*den(nl,nm1)
           ls_pe(nl,nm1) = ls_pe(nl,nm1) + rpe(np,nl)
           pr_pe(nl,nm2) = pr_pe(nl,nm2) + rpe(np,nl)
           pr_pe(nl,nm3) = pr_pe(nl,nm3) + rpe(np,nl)
           pr_pe(nl,nm4) = pr_pe(nl,nm4) + rpe(np,nl)
           pr_pe(nl,nm5) = pr_pe(nl,nm5) + rpe(np,nl)
        END DO

        !  .. Chemical Reactions

        pr_chem(nl,:) = zero; ls_chem(nl,:) = zero; rct(:,nl) = zero
        DO nr = 1, nrct
            IF((itype(nr) == 1) .or. (itype(nr) == 7)) THEN
               nm1=irct(1,nr); nx1 = lopc(nm1)
               nm3=irct(3,nr); nx3 = lopc(nm3)
               nm4=irct(4,nr); nx4 = lopc(nm4)
               nm5=irct(5,nr); nx5 = lopc(nm5)
               rct(nr,nl) = rt(nr,nl)*den(nl,nm1)
               ls_chem(nl,nm1) = ls_chem(nl,nm1) + rct(nr,nl)
               pr_chem(nl,nm3) = pr_chem(nl,nm3) + rct(nr,nl)
               pr_chem(nl,nm4) = pr_chem(nl,nm4) + rct(nr,nl)
               pr_chem(nl,nm5) = pr_chem(nl,nm5) + rct(nr,nl)

               rjac(nx1,nx1) = rjac(nx1,nx1) - rt(nr,nl)
               rjac(nx3,nx1) = rjac(nx3,nx1) + rt(nr,nl)
               rjac(nx4,nx1) = rjac(nx4,nx1) + rt(nr,nl)
               rjac(nx5,nx1) = rjac(nx5,nx1) + rt(nr,nl)

            ELSE
               nm1=irct(1,nr); nx1 = lopc(nm1)
               nm2=irct(2,nr); nx2 = lopc(nm2)
               nm3=irct(3,nr); nx3 = lopc(nm3)
               nm4=irct(4,nr); nx4 = lopc(nm4)
               nm5=irct(5,nr); nx5 = lopc(nm5)
               r1 = rt(nr,nl)*den(nl,nm1); r2 = rt(nr,nl)*den(nl,nm2)
               rct(nr,nl) = rt(nr,nl)*den(nl,nm1)*den(nl,nm2)
               ls_chem(nl,nm1) = ls_chem(nl,nm1) + rct(nr,nl)
               ls_chem(nl,nm2) = ls_chem(nl,nm2) + rct(nr,nl)
               pr_chem(nl,nm3) = pr_chem(nl,nm3) + rct(nr,nl)
               pr_chem(nl,nm4) = pr_chem(nl,nm4) + rct(nr,nl)
               pr_chem(nl,nm5) = pr_chem(nl,nm5) + rct(nr,nl)

               rjac(nx1,nx2) = rjac(nx1,nx2) - r1
               rjac(nx1,nx1) = rjac(nx1,nx1) - r2
               rjac(nx2,nx2) = rjac(nx2,nx2) - r1
               rjac(nx2,nx1) = rjac(nx2,nx1) - r2
               rjac(nx3,nx1) = rjac(nx3,nx1) + r2
               rjac(nx3,nx2) = rjac(nx3,nx2) + r1
               rjac(nx4,nx1) = rjac(nx4,nx1) + r2
               rjac(nx4,nx2) = rjac(nx4,nx2) + r1
               rjac(nx5,nx1) = rjac(nx5,nx1) + r2
               rjac(nx5,nx2) = rjac(nx5,nx2) + r1

            END IF

        END DO

        DO nm = 1, nsp
           pr(nl,nm) = prext(nl,nm) + pr_ph(nl,nm) + pr_pe(nl,nm) + pr_chem(nl,nm)
           ls(nl,nm) = ls_ph(nl,nm) + ls_pe(nl,nm) + ls_chem(nl,nm)
           dNdt(nl,nm)=(den(nl,nm)-den_old(nl,nm))*tinv
        END DO

        ! .. Set up equation

        del = zero; fb=zero; fd=zero
        DO nx = 1, nchem
           nm = locp(nx)
           fd(nx) = pr(nl,nm)-ls(nl,nm)-dNdt(nl,nm)
           del(nx) = fd(nx)
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

        CALL GESV(fb,del) ! solves fb*x = del (or fd), overwrites del with x

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

     END DO NEWTON

  END DO ALT

  DEALLOCATE(rjac,fb,fd,del)

  RETURN
END SUBROUTINE CHEMEQ