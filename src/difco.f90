SUBROUTINE DIFCO

  USE PRECISION
  USE CONSTANTS
  USE GLOBAL_VARIABLES

  !  .. Declarations

  REAL(RP), ALLOCATABLE, DIMENSION(:) :: ekp, grvp, tnp, dtndr, dtpdr, tpl, tip, tep, dtpldr
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: htp, dfp
  REAL(RP) :: dN2, dCO2, mN2, mCO2, rmN2, nuN2, bN2, rmCO2, nuCO2, bCO2, sigma, mu
  INTEGER :: nl, nx, nm
  CHARACTER(len=1) :: iret
  
  !  .. Allocate External Arrays

  ALLOCATE(df(nlev,nsp),alpha(nlev,ndiff),beta(nlev,ndiff),a(nlev,ndiff),b(nlev,ndiff),c(nlev,ndiff),    &
       alphax(nlev,ndiff),betax(nlev,ndiff))

  !  .. Allocate Local Arrays

  ALLOCATE(ekp(nlev),grvp(nlev),tnp(nlev),htp(nlev,0:nsp),dfp(nlev,nsp),dtndr(nlev),dtpdr(nlev),    &
       tpl(nlev),tip(nlev),tep(nlev),dtpldr(nlev))

  !
  !  .. Calculate Neutral Molecular diffusion coefficients
  !

  DO nm = 1, neutrmax
    DO nl = 1, nlev
         IF (dtype(nm) == 0) THEN                                                       ! hard sphere approximation
            sigma = 3.0E-15_RP
            mu = amu/(one/mmw(nm) + one/mmw(iCO2))
            df(nl,nm) = (three/eight*sqpi)*(one/sigma)*SQRT(rkb*tn(nl)/mu)/den(nl,0)
        ELSE IF (dtype(nm) == 1) THEN                                                  ! Mason & Morrero 1970 formula 136
            df(nl,nm) = ad(nm)*(tn(nl)**sd(nm))*EXP(-sd_2(nm)/tn(nl))/prs(nl)
        ELSE IF (dtype(nm) == 2) THEN
            df(nl,nm) = (ad(nm)*(tn(nl)**sd(nm))*EXP(-sd_2(nm)/tn(nl))            &     ! Mason & Morrero 1970 formula 135
                 *EXP(-sd_3(nm)/tn(nl)**2))/(prs(nl)*(LOG(phi(nm)/tn(nl)))**2)
        ELSE IF (dtype(nm) == 3) THEN                                                  !  Need to add Lennard-Jones approximation
            WRITE(*,"('Incorrect diffusion coefficient type, EXITING ... ')")
            STOP
        ELSE
            WRITE(*,"('Incorrect diffusion coefficient type, EXITING ... ')")
            STOP
        END IF
    END DO
  END DO

  !  .. Values at n+1/2 

  DO nl = 1, nlev-1
     tnp(nl) = half*(tn(nl+1)+tn(nl))
     grvp(nl) = half*(grv(nl+1)+grv(nl))
     dtndr(nl) = (tn(nl+1)-tn(nl))/drp(nl)
     htp(nl,0) = rkb*tnp(nl)/grvp(nl)/mass(nl)/amu
     ekp(nl) = half*(ek(nl+1)+ek(nl))

     DO nx = 1, ndiff
        nm = ldcp(nx)
        htp(nl,nm) = rkb*tnp(nl)/grvp(nl)/mmw(nm)/amu
        dfp(nl,nm) = half*(df(nl+1,nm)+df(nl,nm))
     END DO

  END DO

  !  .. Set boundary conditions

  DO nx = 1, ndiff
     nm = ldcp(nx)

     ! .. Bottom Boundary
   
     IF(ibnd(nm,1) == 1) THEN                           ! Chemical Eq.
        bval(nm,1) = zero
!     ELSE IF (ibnd(nm,1) == 2) THEN                     ! Max Velocity
!        IF(bval(nm,1) /= zero) THEN
!           bval(nm,1) = -ek(1)/ht(1,0)
!        ELSE
!           bval(nm,1) = zero
!        END IF
     END IF
     
     ! .. Top Boundary
     
     IF(ibnd(nm,2) == 1) THEN                           ! Chemical Eq.
        bval(nm,2) =    zero
     ELSE IF (ibnd(nm,2) == 2) THEN                     ! Jean's Velocity
           bval(nm,2) = bval(nm,2)*WJEANS (mmw(nm), rz(nlev), tn(nlev)) ! bval(nm,2) = 1: escaping at Jeans rate
     END IF
     
  END DO

  ! Set up diffusion

  DO nx = 1, ndiff
     nm = ldcp(nx)
     DO nl =1, nlev-1
        alpha(nl,nx) = (dfp(nl,nm) + ekp(nl))/drp(nl)
        beta(nl,nx) = half*(dfp(nl,nm)/htp(nl,nm) + ekp(nl)/htp(nl,0)            &
             + (dfp(nl,nm)+ekp(nl))*dtndr(nl)/tnp(nl))
     END DO
  END DO
  
  
  DO nx = 1, ndiff
     nm = ldcp(nx)
     
     !  .. Interior

     DO nl = 2, nlev-1                  
        a(nl,nx) = -(alpha(nl-1,nx)-beta(nl-1,nx))*rm2(nl)/dr(nl)
        b(nl,nx) = ((alpha(nl,nx)-beta(nl,nx))*rp2(nl)/dr(nl)         &
             + (alpha(nl-1,nx)+beta(nl-1,nx))*rm2(nl)/dr(nl))
        c(nl,nx) = -(alpha(nl,nx)+beta(nl,nx))*rp2(nl)/dr(nl)
     END DO
     
     !  .. Lower Boundary
     
     nl = 1
     IF (ibnd(nm,1) == 1) THEN                   !  .. Chemical Equilibrium
        b(nl,nx) = zero           
        c(nl,nx) = zero
     ELSE IF (ibnd(nm,1) == 2) THEN              !  .. Fixed Velocity        
        b(nl,nx) = -(-alpha(nl,nx)+beta(nl,nx)                              &
             +half*bval(nm,1))/drp(nl)
        c(nl,nx) = -(alpha(nl,nx)+beta(nl,nx)                               &
             +half*bval(nm,1))/drp(nl)
     ELSE IF (ibnd(nm,1) == 3) THEN              !  .. Fixed Mole Fraction
        b(nl,nx) = zero           
        c(nl,nx) = zero
     ELSE
        WRITE(*,*) ' DIFCO: ERROR IN LOWER B.C., EXITING ...'
        STOP
     END IF
     
     !  .. Top Boundary
     
     nl = nlev
     IF(ibnd(nm,2) == 1) THEN                     !  .. Chemical Equilibrium
        a(nl,nx) = zero
        b(nl,nx) = zero
     ELSE IF(ibnd(nm,2) == 2) THEN                !  .. Jeans Velocity
        a(nl,nx) = -(alpha(nl-1,nx)-beta(nl-1,nx)                           &
             -half*bval(nm,2))*rm2(nl)/drp(nl-1)
        b(nl,nx) = (alpha(nl-1,nx)+beta(nl-1,nx)                            &
             +half*bval(nm,2))*rm2(nl)/drp(nl-1)
     ELSE IF(ibnd(nm,2) == 3) THEN                !  .. Fixed Velocity
        a(nl,nx) = -(alpha(nl-1,nx)-beta(nl-1,nx)                           &
             -half*bval(nm,2))*rm2(nl)/drp(nl-1)
        b(nl,nx) = (alpha(nl-1,nx)+beta(nl-1,nx)                            &
             +half*bval(nm,2))*rm2(nl)/drp(nl-1)
     ELSE IF(ibnd(nm,2) == 4) THEN                !  .. Fixed Flux
        a(nl,nx) = -(alpha(nl-1,nx)                                         &
             -beta(nl-1,nx))*rm2(nl)/drp(nl-1)
        b(nl,nx) = (alpha(nl-1,nx)                                          &
             +beta(nl-1,nx))*rm2(nl)/drp(nl-1)
     ELSE
        WRITE(*,*) ' DIFCO: ERROR IN UPPER B.C., EXITING ...'
        STOP
     END IF
     
  END DO
  
  !
  !  .. Ions
  !

  !  .. Calculate plasma temperature

  tpl = ti + te

  !  .. Ion Molecular diffusion coefficients
  
  DO nl = 1, nlev
     dN2 = den(nl,iN2)
!     dCH4 = den(nl,iCH4) 
     dCO2 = den(nl,iCO2)    
     mN2 = mmw(iN2)
!     mCH4 = mmw(iCH4)
     mCO2 = mmw(iCO2)
     DO nm = neutrmax+1, nsp-1
        rmN2 = amu/(one/mN2 + one/mmw(nm))
        nuN2 = 3.33E-09_RP*SQRT(polN2/rmN2)*mN2*dN2/(mN2+mmw(nm))
        bN2 = rkb*ti(nl)/amu/mmw(nm)/nuN2
!        rmCH4 = amu/(one/mCH4 + one/mmw(nm))
!        nuCH4 = 3.33E-09_RP*SQRT(polCH4/rmCH4)*mCH4*dCH4/(mCH4+mmw(nm))
!        bCH4 = rkb*ti(nl)/amu/mmw(nm)/nuCH4
!        df(nl,nm) = one/(one/bN2+one/bCH4)
        rmCO2 = amu/(one/mCO2 + one/mmw(nm))
        nuCO2 = 3.33E-09_RP*SQRT(polCO2/rmCO2)*mCO2*dCO2/(mCO2+mmw(nm))
        bCO2 = rkb*ti(nl)/amu/mmw(nm)/nuCO2
        df(nl,nm) = one/(one/bN2+one/bCO2)
     END DO    
  END DO

  !  .. Diffusion Coeff at n+1/2

  DO nl = 1, nlev-1
     DO nm = neutrmax+1, nsp-1
        dfp(nl,nm) = half*(df(nl,nm)+df(nl+1,nm))
     END DO
  END DO

  !  .. Temperature at n+1/2

  DO nl = 1, nlev-1
     tip(nl) = half*(ti(nl+1)+ti(nl))
     tep(nl) = half*(te(nl+1)+te(nl))
  END DO

  !  .. Plasma temperature Derivatives at n+1/2

  DO nl = 1, nlev-1
     dtpldr(nl) = (tpl(nl+1)-tpl(nl))/drp(nl)
  END DO

  !  .. Scale Heights at n+1/2

  DO nl = 1, nlev-1
     DO nm = neutrmax+1, nsp-1
        htp(nl,nm) = rkb*tip(nl)/grvp(nl)/mmw(nm)/amu
     END DO
  END DO


!  DO nx = 1, ndiff
!     nm = lidf(nx)
!     DO nl = 2, nlev-1
!        alphax(nl,nx) = dfp(nl,nm)/drp(nl) 
!        betax(nl,nx)=dfp(nl,nm)/htp(nl,nm)                                &
!             + dfp(nl,nm)*dtpldr(nl)/tip(nl)                              
!     END DO
!  END DO
  
!  DO nx = 1, nidf
!     nm = lidf(nx)

     ! .. Bottom Boundary
   
!     IF(ibnd(nm,1) == 1) THEN                           ! Chemical Eq.
!        bval(nm,1) = zero
!     ELSE IF (ibnd(nm,1) == 2) THEN                     ! Max Velocity
!        IF(bval(nm,1) /= zero) THEN
!           bval(nm,1) = -ek(1)/ht(1,0)
!        ELSE
!           bval(nm,1) = zero
!        END IF
!     END IF
     
     ! .. Top Boundary
     
!     IF(ibnd(nm,2) == 1) THEN                           ! Chemical Eq.
!        bval(nm,2) =    zero
!     ELSE IF (ibnd(nm,2) == 2) THEN                     ! Jean's Velocity
!        bval(nm,2) = WJEANS (GM, half*mmw(nm), rz(nlev), tn(nlev))
!           bval(nm,2) = 1.E4_RP
!     END IF
     
!  END DO


  DEALLOCATE(ekp,grvp,tnp,htp,dfp,dtndr,dtpdr,tpl,tip,tep,dtpldr,alphax,betax)  

  RETURN

CONTAINS

  FUNCTION WJEANS (M, R, T)
    USE PRECISION
    IMPLICIT none
    REAL(RP) M, R, T, lam, u, ams, WJEANS
    REAL(RP), PARAMETER :: amu=1.66E-24, rkb=1.38E-16
    REAL(RP), PARAMETER :: one = 1., two = 2., pi = 3.1415926
    ams = M*amu
    lam = GM*ams/rkb/T/R
    u = SQRT(two*rkb*T/ams)
    WJEANS = (u/two/SQRT(pi)) * (one + lam) * EXP(-lam)
    RETURN
  END FUNCTION WJEANS

END SUBROUTINE DIFCO
