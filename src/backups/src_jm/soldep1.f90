SUBROUTINE SOLDEP1

  !
  !  .. Modules
  !

  USE PRECISION
  USE CONSTANTS
  USE GLOBAL_VARIABLES
  USE SUBS, ONLY : FIND_BIN

  !
  !  .. Declarations
  !

  IMPLICIT  none

  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: clm
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: trnA
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: trnB
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: trnC
  REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: prtA
  REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: prtB
  REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: prtC
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: sol1, sol2, sol3, sol4, sol5
  REAL(RP) :: sin_sza, sa, ca, xp, yv, htop, sm, fac, Eion, Ehv, Eel
  INTEGER :: nm, nz, nl, nw, na, nph, nb, ne

  !
  !  .. Initialize
  !
  
  ALLOCATE(clm(nlev,0:nsp),trnA(ncrsA,nlev),trnB(ncrsB,nlev),trnC(ncrsC,nlev),                     &
       prtA(ncrsA,nbrmaxA,nabsA),prtB(ncrsB,nbrmaxB,nabsB),prtC(ncrsC,nbrmaxC,nabsC))
 
  IF(lcrsA) THEN
  DO na = 1, nabsA
     DO nb = 1, nbrnchA(na)
        DO nw = 1, ncrsA
        prtA(nw,nb,na) = bratA(nw,nb,na)*xcrsA(nw,na)*fsolA(nw)
        END DO
     END DO
  END DO
  END IF

  IF(lcrsB) THEN
  DO na = 1, nabsB
     DO nb = 1, nbrnchB(na)
        DO nw = 1, ncrsB
        prtB(nw,nb,na) = bratB(nw,nb,na)*xcrsB(nw,na)*fsolB(nw)
        END DO
     END DO
  END DO
  END IF

  IF(lcrsC) THEN
  DO na = 1, nabsC
     DO nb = 1, nbrnchC(na)
        DO nw = 1, ncrsC
        prtC(nw,nb,na) = bratC(nw,nb,na)*xcrsC(nw,na)*fsolC(nw)
        END DO
     END DO
  END DO
  END IF
  
  !
  !  ..  Calculate Column Densities to Sun
  !

  sin_sza = SQRT(one-cos_sza*cos_sza)
  
  IF (illum == 1) THEN             !  .. dayside
     DO nm = 1, nsp-1
        DO nz = 1, nlev
           sa = rz(nz)*sin_sza/rz(nlev)    !  .. Above top layer
           ca = SQRT(one-sa*sa)
           IF( ht(nlev,nm) > 1.E-10_RP) THEN
              xp = rz(nlev)/ht(nlev,nm)
           ELSE
	      xp = zero
           END IF
           IF (xp < 100._RP) THEN
              yv = SQRT(half*xp) * ABS(ca)
              htop = ht(nlev,nm) * SQRT(pi*xp/two) * EXP(yv*yv) * ERFC(yv)
           ELSE
              htop = ht(nlev,nm)/cos_sza
           END IF
           sm = htop * den(nlev,nm)
           DO nl = nlev-1, nz, -1  
              sm = sm + half*(den(nl,nm)+den(nl+1,nm)) * ds(nl,nz)
           END DO
           clm(nz,nm) = sm
        END DO
     END DO

  ELSE IF (illum == 0) THEN       !  .. twilight

     DO nm = 1, nsp-1
        clm(1:nbot-1,nm) = 1.E33_RP
        DO nz = nbot, nlev
           sa = rz(nz)*sin_sza/rz(nlev)     ! .. Above top layer
           ca = SQRT(one-sa*sa)
           xp = rz(nlev)/ht(nlev,nm)
           IF (xp < 100._RP) THEN
              yv = SQRT(half*xp) * ABS(ca)
              htop = ht(nlev,nm) * SQRT(pi*xp/two) * EXP(yv*yv) * ERFC(yv)
           ELSE
              htop = ht(nlev,nm)/cos_sza
           END IF
           sm = htop * den(nlev,nm)
           DO nl = MAX(ntan(nz),1), nz-1           ! interior to rz
              sm = sm + (den(nl,nm)+den(nl+1,nm)) * ds(nl,nz)
           END DO
           DO nl = nz, nlev-1               ! exterior to rz
              sm = sm + half*(den(nl,nm)+den(nl+1,nm)) * ds(nl,nz)
           END DO
           clm(nz,nm) = sm
        END DO
     END DO

  ELSE IF (illum == -1) THEN      !  .. nightside

     DO nm = 1, nsp-1
        DO nz = 1, nlev
           clm(nz,nm) = 1.E33_RP
        END DO
     END DO

  END IF

  !
  !  .. Calculate gas absorption optical depth at each level and wavelength
  !
  
  IF(lcrsA) THEN
  DO nz = 1, nlev
     DO nw = 1, ncrsA
        sm = zero
        DO na = 1, nabsA
           nm = loabA(na)
           sm = sm + xcrsA(nw,na) * clm(nz,nm)
        END DO
        trnA(nw,nz) = EXP(-sm)
     END DO
  END DO
  END IF
  
  IF(lcrsB) THEN
  DO nz = 1, nlev
     DO nw = 1, ncrsB
        sm = zero
        DO na = 1, nabsB
           nm = loabB(na)
           sm = sm + xcrsB(nw,na) * clm(nz,nm)
        END DO
        trnB(nw,nz) = EXP(-sm)
     END DO
  END DO
  END IF
  
  IF(lcrsC) THEN
  DO nz = 1, nlev
     DO nw = 1, ncrsC
        sm = zero
        DO na = 1, nabsC
           nm = loabC(na)
           sm = sm + xcrsC(nw,na) * clm(nz,nm)
        END DO
        trnC(nw,nz) = EXP(-(sm + (tau_aer(nz,nw) + tau_ray(nz,nw))/cos_sza)) 
     END DO
  END DO
  END IF
  
  !
  !  .. Calculate absorption rate for each molecule at each level
  !     for each wavelength. 
  !

  nph = 0
  IF(lcrsA) THEN
  DO na = 1, nabsA
     DO nb = 1, nbrnchA(na)
        nph = nph + 1
           DO nz = 1, nlev
              sm = zero
              DO nw = 1, ncrsA
                 sm = sm + prtA(nw,nb,na)*trnA(nw,nz)
              END DO
              rph(nph,nz) = diurnal_average*sm      
           END DO
     END DO
  END DO
  ELSE
  nph = 0
  DO na = 1, nabsA
     DO nb = 1, nbrnchA(na)
        nph = nph + 1
     END DO
  END DO
  END IF
  
  IF(lcrsB) THEN
  DO na = 1, nabsB                                 
     DO nb = 1, nbrnchB(na)
        nph = nph + 1
           DO nz = 1, nlev
              sm = zero
              DO nw = 1, ncrsB
                 sm = sm + prtB(nw,nb,na)*trnB(nw,nz)
              END DO
              rph(nph,nz) = diurnal_average*sm     
           END DO
     END DO
  END DO
  ELSE
  DO na = 1, nabsB                                 
     DO nb = 1, nbrnchB(na)
        nph = nph + 1
     END DO
  END DO
  END IF

  IF(lcrsC) THEN  
  DO na = 1, nabsC                                 
     DO nb = 1, nbrnchC(na)
        nph = nph + 1
           DO nz = 1, nlev
              sm = zero
              DO nw = 1, ncrsC
                 sm = sm + prtC(nw,nb,na)*trnC(nw,nz)
              END DO
              rph(nph,nz) = diurnal_average*sm     
           END DO
     END DO
  END DO
  ELSE
  DO na = 1, nabsC                                 
     DO nb = 1, nbrnchC(na)
        nph = nph + 1
     END DO
  END DO
  END IF

  IF (illum >= 0.) THEN                                      !  .. Optically thin photolysis
  IF(lcrsJ) THEN     
     DO na = 1, nabsJ
        DO nb = 1, nbrnchJ(na)
           nph = nph + 1
           rph(nph,1:nlev) = diurnal_average*srateJ(nb,na)
        END DO
     END DO
  END IF
  ELSE 
     DO na = 1, nabsJ
        DO nb = 1, nbrnchJ(na)
           nph = nph + 1
           rph(nph,1:nlev) = zero
        END DO
     END DO
  END IF

  !
  !  .. Photoelectron production
  !

  nph = 0
  IF(lcrsA) THEN
  DO na = 1, nabsA
     nm = loabA(na)  
     DO nb = 1, nbrnchA(na)
        nph = nph + 1
        IF(charge_stateA(nb,na) > zero) THEN
           esrc(nph,:,:) = zero
           Eion = (1.24E4_RP/enrgIA(nb,na))                      ! Ionization limit (eV)
           DO nw = 1, ncrsA
              Ehv = (1.24E4_RP/wcrsA(nw))         ! Solar photon energy eV  (why do this everytime?)
              Eel = (Ehv - Eion)/charge_stateA(nb,na)
              IF (Eel > zero) THEN                     
                 ne = FIND_BIN(elctreV,elctDeV,Eel)
                 fac = (Eel/elctreV(ne)/elctDeV(ne))*prtA(nw,nb,na)*charge_stateA(nb,na)
                DO nz = nibot, nlev
                    esrc(nph,nz,ne)=esrc(nph,nz,ne)+fac*den(nz,nm)*trnA(nw,nz)
                 END DO
              END IF
           END DO
        END IF
     END DO
  END DO
  END IF

  IF(lcrsB) THEN
  DO na = 1, nabsB
     nm = loabB(na)  
     DO nb = 1, nbrnchB(na)
        nph = nph + 1
        IF(charge_stateB(nb,na) > zero) THEN
           esrc(nph,:,:) = zero
           Eion = (1.24E4_RP/enrgIB(nb,na))                      
           DO nw = 1, ncrsB
              Ehv = (1.24E4_RP/wcrsB(nw))         
              Eel = (Ehv - Eion)/charge_stateB(nb,na)           
              IF (Eel > zero) THEN                     
                 ne = FIND_BIN(elctreV,elctDeV,Eel)
                 fac = (Eel/elctreV(ne)/elctDeV(ne))*charge_stateB(nb,na)*prtB(nw,nb,na)
                 DO nz = nibot, nlev
                    esrc(nph,nz,ne)=esrc(nph,nz,ne)+fac*den(nz,nm)*trnB(nw,nz)
                 END DO
              END IF
           END DO
        END IF
     END DO
  END DO
  END IF

  IF(lcrsC) THEN
  DO na = 1, nabsC
     nm = loabC(na)  
     DO nb = 1, nbrnchC(na)
        nph = nph + 1
        IF(charge_stateC(nb,na) > zero) THEN
           esrc(nph,1:nlev,1:nelb) = zero
           Eion = (1.24E4_RP/enrgIC(nb,na))                      
           DO nw = 1, ncrsC
              Ehv = (1.24E4_RP/wcrsC(nw))         
              Eel = (Ehv - Eion)/charge_stateC(nb,na)           
              IF (Eel > zero) THEN                     
                 ne = FIND_BIN(elctreV,elctDeV,Eel)
                 fac = Eel/elctreV(ne)/elctDeV(ne)
!                 DO nz = nibot, nlev
!                    esrc(nph,nz,ne)=esrc(nph,nz,ne)+charge_stateC(nb,na)*fac*prtC(nw,nb,na)*den(nz,nm)*trnC(nw,nz)
!                 END DO
              END IF
           END DO
        END IF
     END DO
  END DO
  END IF

  esrc(:,:,:) = diurnal_average*esrc(:,:,:)

  Sel(:,:) = zero
  DO ne = 1, nelb
     DO nz = 1, nlev
        Sel(nz,ne) = SUM(esrc(1:nphrt,nz,ne))
     END DO
  END DO

  !
  !  .. Write output file  .. move this to compout
  !

  IF( lcrsA .and. lcrsB .and. lcrsC) THEN
  ALLOCATE(sol1(nlev),sol2(nlev),sol3(nlev),sol4(nlev),sol5(nlev))
  DO nz = 1, nlev
     sm = zero
     DO nw = 1, ncrsA
        sm = sm + fsolA(nw)*trnA(nw,nz)
     END DO
     sol1(nz) = sm
  END DO

  DO nz = 1, nlev
     sm = zero
     DO nw = 1, ncrsB
        sm = sm + fsolB(nw)*trnB(nw,nz)
     END DO
     sol2(nz) = sm
  END DO

  DO nz = 1, nlev
     sm = zero
     DO nw = 1, 45
        sm = sm + fsolC(nw)*trnC(nw,nz)
     END DO
     sol3(nz) = sm
  END DO

  DO nz = 1, nlev
     sm = zero
     DO nw = 46, 80
        sm = sm + fsolC(nw)*trnC(nw,nz)
     END DO
     sol4(nz) = sm
  END DO

  DO nz = 1, nlev
     sm = zero
     DO nw = 81, ncrsC
        sm = sm + fsolC(nw)*trnC(nw,nz)
     END DO
     sol5(nz) = sm
  END DO


  OPEN(unit=70,file='solar_flux.out',status='unknown')
     DO nz = 1, nlev
        WRITE(70,"(6ES11.3)") 1.E-5_RP*rz(nz)-2575._RP, sol1(nz),sol2(nz),sol3(nz),sol4(nz),sol5(nz)
     END DO
  CLOSE(unit=70)
  DEALLOCATE(sol1,sol2,sol3,sol4,sol5)
  END IF

  !
  !  .. All Done
  !

  DEALLOCATE(clm,trnA,trnB,trnC,prtA,prtB,prtC)
  
  RETURN
END SUBROUTINE SOLDEP1
