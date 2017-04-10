!###################################################################################################
!#                                                                                                 #
!# Calculation of electron mean intensity based on the photoelectron flux provided by soldep1.f90. #
!# The subroutine assumes no vertical transport and solves backwards for the mean intensity at     #
!# each energy bin assuming zero input for the highest energy bin.                                 #
!#                                                                                                 #
!###################################################################################################

SUBROUTINE ELDEP1

  USE PRECISION
  USE CONSTANTS
  USE GLOBAL_VARIABLES
  USE SUBS, ONLY : FIND_BIN, DIFCS_SEC

  !
  !  .. Declarations 
  !
 
  IMPLICIT NONE

  INTEGER :: nabs_el
  INTEGER :: i, is, j, k, m, n, ni, nm, npe, na, nb, nl, ne
  INTEGER, DIMENSION(2) :: isp
  REAL(RP), ALLOCATABLE, DIMENSION(:) :: Selz, src_plasma
  REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: Ael, Ael_ext, Ael_sec, Ael_ion, sum_cs_ion, sum_cs_ext
  REAL(RP) :: sum1,sum2, sum_tot, sm, rn, rnp1, Ann, An1n1, E1, E2, E2_min, E2_max, step, Ei,Eion,  &
       sm_ion, sm_sec, E3, E4, tei, eden, sigma_excite_new, sigma_dif_ion, sigma_dif_sec, pSz,  &
       cs_ion, cs_sec

  !
  !  .. Initialize
  ! 

  nabs_el = SIZE(ipath,1)
  eFLUX(:,:) = zero
  isp(1) = iN2
  isp(2) = iCH4

  ALLOCATE(Selz(nelb),src_plasma(nelb),Ael(nelb,nelb), Ael_ext(nelb,nelb), Ael_sec(nelb,nelb),      &
       Ael_ion(nelb,nelb),sum_cs_ion(nelb,nabs_el),sum_cs_ext(nelb,nabs_el))


  !
  !  .. Calculate electron flux at each energy and level
  !

  ALTI: DO i = nibot, nlev
     Ael(:,:) = zero
     Ael_ext(:,:) = zero
     Ael_sec(:,:) = zero
     Ael_ion(:,:) = zero
     Selz(:)=Sel(i,:)
     tei = te(i)
     eden = den(i,iELE)

     if( SUM(Selz(1:nelb)) > zero ) then 

     sum_cs_ext(:,:) = zero 
     sum_cs_ion(:,:) = zero


     ENERGY: DO m = nelb,2,-1    
!     WRITE(*,"(' Energy index = ',I6)") m 
     SPECIES: DO is = 1, nabs_el
!     WRITE(*,"(' Species index = ',I6)") is
        ! ==> Matrix elements for production at energy n from energy m


        ! m is the primary 
        ! n is the secondary
        ! nm is the degraded primary

        ! Secondary electron production
!        WRITE(*,"(' Secondary electron production')")
        DO j = 1, ipath(is,3)

           E1 = elctreV(m)-enrgE(is,ipath(is,1)+ipath(is,2)+j)
           if(E1 > zero) then 
              ni = FIND_BIN(elctreV,elctDeV,E1)
              if(ni.eq.m) stop
              if(elctreV(ni) > E1.and.ni /= 1) then
                 E4 = elctreV(ni)-half*elctDeV(ni)
                 ni =ni-1
              else
                 E4 = elctreV(ni)+half*elctDeV(ni)
              endif
              
              Eion = elctreV(m) - E4
              cs_ion = SUM_CS_REES_ion(m,is,j)
              cs_sec = SUM_CS_REES_sec(m,is,j)
             
              if(ni > 0) then  
                 DO nm = ni, 1, -1 
                    E3 = one 
                    if(nm.eq.ni) E3 = enrgE(is,ipath(is,1)+ipath(is,2)+j)/Eion
                    E2 = elctreV(nm) 
                    E2_min = elctreV(nm)-half*elctDeV(nm)
                    E2_max = elctreV(nm)+half*elctDeV(nm)
                    step = (E2_max-E2_min)/one
                    sm_ion = zero
                    sm_sec = zero
                    sm = zero
                    DO k=1,1
                       Ei = E2_min + (k-1)*step + half*step
                       sm_ion = sm_ion + step*DIFCS_SEC(elctreV(m)-Ei-Eion,elctreV(m),Eion) &
                            *eCS(m,is,ipath(is,1)+ipath(is,2)+j)/cs_ion
                       sm_sec = sm_sec + step*DIFCS_SEC(Ei,elctreV(m),Eion)*eCS(m,is,ipath(is,1) &
                            +ipath(is,2)+j)/cs_sec
                    END DO
                    sigma_dif_ion = sm_ion/elctDeV(nm)
                    sigma_dif_sec = sm_sec/elctDeV(nm)
                    Ael_ion(nm,m) = Ael_ion(nm,m) + den(i,isp(is))*elctDeV(m)*(sigma_dif_ion)*E3 
                    Ael_sec(nm,m) = Ael_sec(nm,m) + den(i,isp(is))*elctDeV(m)*(sigma_dif_sec)                     
                 END DO
              endif
              sum_cs_ion(m,is) = sum_cs_ion(m,is) + eCS(m,is,ipath(is,1)+ipath(is,2)+j) 
           endif
        END DO
        
        ! Excitation to discrete states
!        WRITE(*,"(' Excitation to discrete states')")
        DO j = 1, ipath(is,1)+ipath(is,2)
           sigma_excite_new = zero
           E1 = elctreV(m)-enrgE(is,j)
           nm = FIND_BIN(elctreV,elctDeV,E1)
           IF(E1.gt.0..and.nm.lt.m) then
              E1 = enrgE(is,j)/(elctreV(m)-elctreV(nm))
              E3=1.
              sigma_excite_new = eCS(m,is,j) !*E1
              Ael_ext(nm,m) = Ael_ext(nm,m) + den(i,isp(is))*sigma_excite_new*E1
           ENDIF
           IF(E1.gt.0..and.nm.eq.m) then
              E3 = enrgE(is,j)/(elctreV(nm)-elctreV(nm-1))
              sigma_excite_new = eCS(m,is,j) !/E3 
              Ael_ext(nm-1,m) = Ael_ext(nm-1,m) + den(i,isp(is))*sigma_excite_new/E3
           ENDIF
           sum_cs_ext(m,is) = sum_cs_ext(m,is) + sigma_excite_new !*E3
        END DO
        
        
     END DO SPECIES
     END DO ENERGY
  
     Ael = Ael_ext + Ael_ion + Ael_sec
  
!     WRITE(*,"(' Electron fulx at each energy level')")
     ! Electron flux at each energy for level i
     eFLUX(i,nelb) = Selz(nelb)/(den(i,iN2)*(sum_cs_ext(nelb,1)+sum_cs_ion(nelb,1)) &
          +den(i,iCH4)*(sum_cs_ext(nelb,2)+sum_cs_ion(nelb,2)))
     
     DO n = nelb-1, 1, -1
!        WRITE(*,"(' n = ',I6)") n
        sum1 = DOT_PRODUCT(Ael(n,n+1:nelb),eFLUX(i,n+1:nelb))
        sum2 = DOT_PRODUCT(Ael(n+1,n+2:nelb),eFLUX(i,n+2:nelb))
        sum_tot = half * (sum1 + sum2)
!        tei = te(i)
!        eden = den(i,iELE)
!        WRITE(*,"(' sum_tot = ',4ES11.3)") sum_tot , elctreV(n)
        rn = eden*LOSS(elctreV(n),eden,tei)
        rnp1 = eden*LOSS(elctreV(n+1),eden,tei)
!        WRITE(*,"(' rn, rnp1 = ',2ES11.3)") rn, rnp1 
        Ann = (rn / (elctreV(n+1)-elctreV(n))) + half*( den(i,iN2)* (sum_cs_ext(n,1) + sum_cs_ion(n,1) ) &
             + den(i,iCH4)* (sum_cs_ion(n,2)   + sum_cs_ext(n,2)) )
        An1n1 = (rnp1 / (elctreV(n+1)-elctreV(n))) - half*( den(i,iN2)* (sum_cs_ext(n+1,1) + sum_cs_ion(n+1,1)) & 
             + den(i,iCH4)* (sum_cs_ion(n+1,2) + sum_cs_ext(n+1,2)) )
!        WRITE(*,"(' Ann, An1n1 = ',2ES11.3)") Ann, An1n1
        IF (ABS(Ann) < 1.E-25_RP) THEN
           eFLUX(i,n) = zero
        ELSE
           eFLUX(i,n) = ( half*(Selz(n)+Selz(n+1)) + An1n1*eFLUX(i,n+1) + sum_tot) / Ann
        ENDIF
        if(eFLUX(i,n) <= zero ) eFLUX(i,n) = zero
!     WRITE(*,"(' Eflux = ',ES11.3)") eFLUX(i,n)     
     END DO

     endif
  
  END DO ALTI

   
  ! Secondary production rates
!  WRITE(*,"(' Secondary production rates')")
  DO k = nibot, nlev
     pSz = zero
     do i = 1, nabs_el
        do n = 1, nelb  !-1
           pSz = pSz + den(k,isp(i))*eCS_ion(n,i)*eflux(k,n)*elctDeV(n) !0.5*(f1+f2)*de
        enddo
     enddo
     pS(k) = pSz
  ENDDO

  
  npe = 0
  DO na = 1, nabs_el_thk
     DO nb = 1, nbrnch_el(na)
        npe = npe + 1
        DO nl = nibot, nlev
           sm = zero
           DO ne = 1, nelb
              sm = sm + eflux(nl,ne)*eCS(ne,na,ipath(na,1)+nb)*elctDeV(ne) 
           END DO
           eph(npe,nl) = sm
        END DO
     END DO
  END DO
  
  DO n = 1, nabs_el_thn ! why is this different?
     na = n + nabs_el_thk
     DO nb = 1, nbrnch_el(na)
        npe = npe + 1
        DO nl = nibot, nlev
           sm = zero
           DO ne = 1, nelb
              sm = sm + brat_el(ne,nb,na)*eFLUX(nl,ne)*crs_tot_inel(ne,na)*elctDeV(ne) 
           END DO
           eph(npe,nl) = sm
        END DO
     END DO
  END DO

  DEALLOCATE(Selz,src_plasma,Ael,Ael_ext,Ael_sec,Ael_ion,sum_cs_ion,sum_cs_ext)

  RETURN

  CONTAINS

  FUNCTION LOSS(E,NE,TE)
    USE PRECISION
    USE CONSTANTS
    IMPLICIT NONE
    REAL(RP) :: LOSS, E, NE, TE, Ee
!    IF(E <= 0._RP) THEN
!	WRITE(*,"(' LOSS: E = ',ES15.7)") E
!        STOP
!    END IF
!    IF(Ne <= 0._RP) THEN
!	WRITE(*,"(' LOSS: Ne = ',ES15.7)") Ne 
!        STOP
!    END IF
    Ee = 8.618E-5_RP*TE
    IF (Ee > E) THEN
	LOSS = zero
    ELSE 
	LOSS = (3.37e-12_RP/((E**0.94_RP)*(NE**0.03_RP)))*((E-Ee)/(E-0.53_RP*Ee))**2.36_RP
    END IF
!    IF (Ee > E) LOSS = zero !(3.37e-12/((E**0.94)*(NE**0.03)))*((Ee-E)/(Ee-0.53*E))**2.36  !zero
!   IF (Ee > E) write(*,*) 'LOSS = ',LOSS
  END FUNCTION LOSS

END SUBROUTINE eldep1
